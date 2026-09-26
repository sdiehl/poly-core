use std::cmp::Ordering;
use std::ops::{Add, Mul, Neg, Sub};

use crate::field::{pow, Field};
use crate::monomial::{Monomial, Order};

pub type Term<F> = (Monomial, F);

/// A sparse multivariate polynomial: nonzero terms in strictly descending `order`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly<F> {
    pub terms: Vec<Term<F>>,
    pub nvars: usize,
    pub order: Order,
}

impl<F: Field> Poly<F> {
    /// Sorts, combines like terms and drops zeros.
    pub fn new(mut terms: Vec<Term<F>>, nvars: usize, order: Order) -> Self {
        terms.sort_by(|a, b| order.compare(&b.0, &a.0));
        let mut out: Vec<Term<F>> = Vec::with_capacity(terms.len());
        for (m, c) in terms {
            match out.last_mut() {
                Some(last) if last.0 == m => last.1 = last.1.clone() + c,
                _ => out.push((m, c)),
            }
        }
        out.retain(|t| !t.1.is_zero());
        Self {
            terms: out,
            nvars,
            order,
        }
    }

    pub const fn zero(nvars: usize, order: Order) -> Self {
        Self {
            terms: Vec::new(),
            nvars,
            order,
        }
    }

    pub fn constant(c: F, nvars: usize, order: Order) -> Self {
        Self::new(vec![(Monomial::one(nvars), c)], nvars, order)
    }

    pub fn var(i: usize, nvars: usize, order: Order) -> Self {
        Self::new(vec![(Monomial::var(i, nvars), F::one())], nvars, order)
    }

    fn with(&self, terms: Vec<Term<F>>) -> Self {
        Self {
            terms,
            nvars: self.nvars,
            order: self.order.clone(),
        }
    }

    pub const fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn is_constant(&self) -> bool {
        self.terms.iter().all(|t| t.0.is_one())
    }

    pub fn lt(&self) -> Option<&Term<F>> {
        self.terms.first()
    }

    pub fn lm(&self) -> Option<&Monomial> {
        self.lt().map(|t| &t.0)
    }

    pub fn lc(&self) -> Option<&F> {
        self.lt().map(|t| &t.1)
    }

    pub fn total_degree(&self) -> u32 {
        self.terms.iter().map(|t| t.0.degree()).max().unwrap_or(0)
    }

    pub fn degree(&self, k: usize) -> u32 {
        self.terms.iter().map(|t| t.0.exps()[k]).max().unwrap_or(0)
    }

    #[must_use]
    pub fn reorder(&self, order: Order) -> Self {
        Self::new(self.terms.clone(), self.nvars, order)
    }

    pub fn map<G: Field>(&self, f: impl Fn(&F) -> G) -> Poly<G> {
        let terms = self.terms.iter().map(|(m, c)| (m.clone(), f(c))).collect();
        Poly::new(terms, self.nvars, self.order.clone())
    }

    #[must_use]
    pub fn scale(&self, c: &F) -> Self {
        if c.is_zero() {
            return self.with(Vec::new());
        }
        self.with(
            self.terms
                .iter()
                .map(|(m, a)| (m.clone(), a.clone() * c.clone()))
                .collect(),
        )
    }

    #[must_use]
    pub fn mul_term(&self, c: &F, m: &Monomial) -> Self {
        let shifted = self
            .terms
            .iter()
            .map(|(n, a)| (n * m, a.clone() * c.clone()));
        self.with(shifted.filter(|t| !t.1.is_zero()).collect())
    }

    #[must_use]
    pub fn monic(&self) -> Self {
        self.lc()
            .and_then(Field::inverse)
            .map_or_else(|| self.clone(), |l| self.scale(&l))
    }

    /// `self - c * m * g` in one merge.
    #[must_use]
    pub fn sub_mul(&self, c: &F, m: &Monomial, g: &Self) -> Self {
        self.merge(&g.mul_term(&-c.clone(), m))
    }

    fn merge(&self, o: &Self) -> Self {
        debug_assert_eq!(self.nvars, o.nvars);
        let mut out = Vec::with_capacity(self.terms.len() + o.terms.len());
        let (mut a, mut b) = (self.terms.iter().peekable(), o.terms.iter().peekable());
        while let (Some(x), Some(y)) = (a.peek(), b.peek()) {
            match self.order.compare(&x.0, &y.0) {
                Ordering::Greater => out.push(a.next().cloned().expect("peeked")),
                Ordering::Less => out.push(b.next().cloned().expect("peeked")),
                Ordering::Equal => {
                    let c = x.1.clone() + y.1.clone();
                    if !c.is_zero() {
                        out.push((x.0.clone(), c));
                    }
                    a.next();
                    b.next();
                }
            }
        }
        out.extend(a.cloned());
        out.extend(b.cloned());
        self.with(out)
    }

    pub fn eval(&self, x: &[F]) -> F {
        self.terms.iter().fold(F::zero(), |acc, (m, c)| {
            let v = m
                .exps()
                .iter()
                .zip(x)
                .fold(c.clone(), |v, (&e, xi)| v * pow(xi, e.into()));
            acc + v
        })
    }

    /// Substitutes `a` for `x_k`.
    #[must_use]
    pub fn eval_var(&self, k: usize, a: &F) -> Self {
        let terms = self.terms.iter().map(|(m, c)| {
            let mut e = m.exps().to_vec();
            let v = c.clone() * pow(a, std::mem::take(&mut e[k]).into());
            (Monomial::new(e), v)
        });
        Self::new(terms.collect(), self.nvars, self.order.clone())
    }

    /// The S-polynomial of two nonzero polynomials.
    #[must_use]
    pub fn spoly(&self, o: &Self) -> Self {
        let ((m1, c1), (m2, c2)) = (self.lt().expect("nonzero"), o.lt().expect("nonzero"));
        let l = m1.lcm(m2);
        let a = self.mul_term(&(F::one() / c1.clone()), &l.quo(m1).expect("lcm"));
        a.sub_mul(&(F::one() / c2.clone()), &l.quo(m2).expect("lcm"), o)
    }

    /// Division by an ordered list: `self = sum q[i] * divs[i] + r`, no term of `r` divisible by a
    /// leading monomial, each step reducing by the first divisor that applies.
    pub fn divide(&self, divs: &[Self]) -> (Vec<Self>, Self) {
        let mut q: Vec<Vec<Term<F>>> = vec![Vec::new(); divs.len()];
        let mut rem = Vec::new();
        let mut work = self.clone();
        while let Some((m, c)) = work.lt().cloned() {
            let hit = divs.iter().enumerate().find_map(|(i, g)| {
                let (gm, gc) = g.lt()?;
                Some((i, m.quo(gm)?, c.clone() / gc.clone()))
            });
            match hit {
                Some((i, t, k)) => {
                    work = work.sub_mul(&k, &t, &divs[i]);
                    q[i].push((t, k));
                }
                None => {
                    rem.push(work.terms.remove(0));
                }
            }
        }
        let q = q
            .into_iter()
            .map(|t| Self::new(t, self.nvars, self.order.clone()))
            .collect();
        (q, self.with(rem))
    }

    /// The remainder of [`Poly::divide`].
    #[must_use]
    pub fn reduce(&self, divs: &[Self]) -> Self {
        self.divide(divs).1
    }
}

/// `sum h[i] * gens[i]`, the polynomial a certificate claims.
pub fn combination<F: Field>(h: &[Poly<F>], gens: &[Poly<F>]) -> Poly<F> {
    let zero = gens.first().map(|g| Poly::zero(g.nvars, g.order.clone()));
    h.iter().zip(gens).fold(
        zero.unwrap_or_else(|| Poly::zero(0, Order::Lex)),
        |acc, (a, g)| &acc + &(a * g),
    )
}

impl<F: Field> Add for &Poly<F> {
    type Output = Poly<F>;
    fn add(self, o: Self) -> Poly<F> {
        self.merge(o)
    }
}

impl<F: Field> Neg for &Poly<F> {
    type Output = Poly<F>;
    fn neg(self) -> Poly<F> {
        self.with(
            self.terms
                .iter()
                .map(|(m, c)| (m.clone(), -c.clone()))
                .collect(),
        )
    }
}

impl<F: Field> Sub for &Poly<F> {
    type Output = Poly<F>;
    fn sub(self, o: Self) -> Poly<F> {
        self.merge(&-o)
    }
}

impl<F: Field> Mul for &Poly<F> {
    type Output = Poly<F>;
    fn mul(self, o: Self) -> Poly<F> {
        let zero = self.with(Vec::new());
        o.terms
            .iter()
            .fold(zero, |acc, (m, c)| acc.merge(&self.mul_term(c, m)))
    }
}
