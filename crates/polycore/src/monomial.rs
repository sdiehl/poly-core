use std::cmp::Ordering;
use std::fmt;
use std::ops::Mul;
use std::str::FromStr;
use std::sync::Arc;

use crate::field::{Field, pow};

/// An exponent vector in one shared allocation, with its total degree cached.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Monomial {
    deg: u32,
    exps: Arc<[u32]>,
}

impl Monomial {
    pub fn new(exps: impl Into<Arc<[u32]>>) -> Self {
        let exps = exps.into();
        Self {
            deg: exps.iter().sum(),
            exps,
        }
    }

    pub fn one(n: usize) -> Self {
        Self::new(vec![0; n])
    }

    pub fn var(i: usize, n: usize) -> Self {
        let mut e = vec![0; n];
        e[i] = 1;
        Self::new(e)
    }

    pub fn exps(&self) -> &[u32] {
        &self.exps
    }

    pub fn nvars(&self) -> usize {
        self.exps.len()
    }

    pub const fn degree(&self) -> u32 {
        self.deg
    }

    pub const fn is_one(&self) -> bool {
        self.deg == 0
    }

    pub fn divides(&self, o: &Self) -> bool {
        self.deg <= o.deg && self.exps.iter().zip(o.exps.iter()).all(|(a, b)| a <= b)
    }

    /// `self / o` when `o` divides `self`.
    pub fn quo(&self, o: &Self) -> Option<Self> {
        o.divides(self).then(|| self.zip(o, |a, b| a - b))
    }

    #[must_use]
    pub fn lcm(&self, o: &Self) -> Self {
        self.zip(o, Ord::max)
    }

    pub fn is_coprime(&self, o: &Self) -> bool {
        self.exps
            .iter()
            .zip(o.exps.iter())
            .all(|(a, b)| *a == 0 || *b == 0)
    }

    /// The value at `x`.
    pub fn eval<F: Field>(&self, x: &[F]) -> F {
        self.exps
            .iter()
            .zip(x)
            .filter(|t| *t.0 > 0)
            .fold(F::one(), |v, (&e, xi)| v * pow(xi, e.into()))
    }

    fn zip(&self, o: &Self, f: impl Fn(u32, u32) -> u32) -> Self {
        debug_assert_eq!(self.nvars(), o.nvars());
        Self::new(
            self.exps
                .iter()
                .zip(o.exps.iter())
                .map(|(a, b)| f(*a, *b))
                .collect::<Vec<_>>(),
        )
    }
}

impl Mul for &Monomial {
    type Output = Monomial;
    fn mul(self, o: Self) -> Monomial {
        self.zip(o, |a, b| a + b)
    }
}

/// `x0*x1^2`, or `1`.
impl fmt::Display for Monomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_one() {
            return f.write_str("1");
        }
        let mut sep = "";
        for (i, &e) in self.exps.iter().enumerate().filter(|t| *t.1 > 0) {
            write!(f, "{sep}x{i}")?;
            if e > 1 {
                write!(f, "^{e}")?;
            }
            sep = "*";
        }
        Ok(())
    }
}

/// A monomial order on exponent vectors, variable 0 most significant.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Order {
    Lex,
    GrLex,
    GRevLex,
    /// Weighted degree first, ties broken by the inner order.
    Weighted(Arc<[u32]>, Arc<Self>),
    /// Product order: each block of variables compared by its own order, left to right.
    Block(Arc<[(Self, usize)]>),
}

impl Order {
    pub fn weighted(weights: impl Into<Arc<[u32]>>, tie: Self) -> Self {
        Self::Weighted(weights.into(), Arc::new(tie))
    }

    pub fn block(blocks: impl Into<Arc<[(Self, usize)]>>) -> Self {
        Self::Block(blocks.into())
    }

    /// Graded reverse lex on the first `k` variables, then on the remaining `rest`.
    pub fn elimination(k: usize, rest: usize) -> Self {
        Self::block(vec![(Self::GRevLex, k), (Self::GRevLex, rest)])
    }

    /// Whether every monomial involving `x_0..x_{k-1}` exceeds every one free of them, so a Groebner
    /// basis meets `F[x_k..]` in a basis of the elimination ideal.
    pub fn eliminates(&self, k: usize) -> bool {
        match self {
            _ if k == 0 => true,
            Self::Lex => true,
            Self::GrLex | Self::GRevLex => false,
            Self::Weighted(w, _) => {
                w.len() >= k && w[..k].iter().all(|&x| x > 0) && w[k..].iter().all(|&x| x == 0)
            }
            Self::Block(blocks) => {
                let mut start = 0;
                for (order, size) in blocks.iter() {
                    if k < start + size {
                        return order.eliminates(k - start);
                    }
                    start += size;
                }
                true
            }
        }
    }

    pub fn compare(&self, a: &Monomial, b: &Monomial) -> Ordering {
        match self {
            Self::GrLex | Self::GRevLex if a.deg != b.deg => a.deg.cmp(&b.deg),
            Self::GrLex => a.exps.cmp(&b.exps),
            Self::GRevLex => b.exps.iter().rev().cmp(a.exps.iter().rev()),
            _ => self.compare_exps(&a.exps, &b.exps),
        }
    }

    pub fn compare_exps(&self, a: &[u32], b: &[u32]) -> Ordering {
        let deg = |e: &[u32]| e.iter().sum::<u32>();
        match self {
            Self::Lex => a.cmp(b),
            Self::GrLex => deg(a).cmp(&deg(b)).then_with(|| a.cmp(b)),
            Self::GRevLex => deg(a)
                .cmp(&deg(b))
                .then_with(|| b.iter().rev().cmp(a.iter().rev())),
            Self::Weighted(w, tie) => {
                let wdeg = |e: &[u32]| -> u64 {
                    e.iter()
                        .zip(w.iter())
                        .map(|(x, y)| u64::from(*x) * u64::from(*y))
                        .sum()
                };
                wdeg(a).cmp(&wdeg(b)).then_with(|| tie.compare_exps(a, b))
            }
            Self::Block(blocks) => {
                let mut start = 0;
                for (order, size) in blocks.iter() {
                    let end = (start + size).min(a.len());
                    match order.compare_exps(&a[start..end], &b[start..end]) {
                        Ordering::Equal => start = end,
                        other => return other,
                    }
                }
                Ordering::Equal
            }
        }
    }
}

impl FromStr for Order {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "lex" => Ok(Self::Lex),
            "grlex" => Ok(Self::GrLex),
            "grevlex" => Ok(Self::GRevLex),
            _ => Err(format!("unknown order {s}")),
        }
    }
}
