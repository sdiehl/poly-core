//! Incremental interpolation.
//!
//! Newton for polynomials, Thiele for rational functions, and the Berlekamp-Massey and transposed
//! Vandermonde steps of sparse interpolation. Each grows one point
//! at a time and reports when a new point was already predicted, so no degree bounds are needed.

use crate::field::Field;
use crate::uni::Uni;

/// Newton's divided differences at distinct abscissae.
#[derive(Clone, Debug)]
pub struct Newton<F> {
    xs: Vec<F>,
    cs: Vec<F>,
}

impl<F> Default for Newton<F> {
    fn default() -> Self {
        Self {
            xs: Vec::new(),
            cs: Vec::new(),
        }
    }
}

impl<F: Field> Newton<F> {
    /// Adds `(x, y)` at a new `x`, returning `false` when the interpolant already passed through it.
    pub fn add(&mut self, x: F, y: F) -> bool {
        let (mut v, mut w) = (F::zero(), F::one());
        for (xi, c) in self.xs.iter().zip(&self.cs) {
            v = v + c.clone() * w.clone();
            w = w * (x.clone() - xi.clone());
        }
        let c = (y - v) / w;
        let fresh = !c.is_zero();
        self.xs.push(x);
        self.cs.push(c);
        fresh
    }

    pub fn contains(&self, x: &F) -> bool {
        self.xs.contains(x)
    }

    pub const fn len(&self) -> usize {
        self.xs.len()
    }

    pub const fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }

    pub fn poly(&self) -> Uni<F> {
        self.xs
            .iter()
            .zip(&self.cs)
            .rev()
            .fold(Uni::zero(), |r, (xi, c)| {
                &(&r * &Uni::new(vec![-xi.clone(), F::one()])) + &Uni::constant(c.clone())
            })
    }
}

/// Thiele's continued fraction `c_0 + (t - t_0) / (c_1 + (t - t_1) / (c_2 + ...))`.
#[derive(Clone, Debug)]
pub struct Thiele<F> {
    ts: Vec<F>,
    cs: Vec<F>,
}

impl<F> Default for Thiele<F> {
    fn default() -> Self {
        Self {
            ts: Vec::new(),
            cs: Vec::new(),
        }
    }
}

impl<F: Field> Thiele<F> {
    /// Adds `(t, y)` by inverted differences, returning `false` if they break down there, as they
    /// do at a point the fraction already predicts.
    pub fn add(&mut self, t: F, y: F) -> bool {
        let mut v = y;
        for (ti, ci) in self.ts.iter().zip(&self.cs) {
            let Some(d) = (v - ci.clone()).inverse() else {
                return false;
            };
            v = (t.clone() - ti.clone()) * d;
        }
        self.ts.push(t);
        self.cs.push(v);
        true
    }

    pub fn contains(&self, t: &F) -> bool {
        self.ts.contains(t)
    }

    pub const fn len(&self) -> usize {
        self.ts.len()
    }

    pub const fn is_empty(&self) -> bool {
        self.ts.is_empty()
    }

    /// The fraction at `t`, or `None` at a pole and before the first point.
    pub fn eval(&self, t: &F) -> Option<F> {
        let mut pairs = self.ts.iter().zip(&self.cs).rev();
        let (_, last) = pairs.next()?;
        pairs.try_fold(last.clone(), |v, (ti, ci)| {
            Some(ci.clone() + (t.clone() - ti.clone()) * v.inverse()?)
        })
    }

    /// The fraction as `(num, den)` in lowest terms with `den` monic.
    pub fn rational(&self) -> (Uni<F>, Uni<F>) {
        let mut pairs = self.ts.iter().zip(&self.cs).rev();
        let Some((_, last)) = pairs.next() else {
            return (Uni::zero(), Uni::constant(F::one()));
        };
        let (num, den) = pairs.fold(
            (Uni::constant(last.clone()), Uni::constant(F::one())),
            |(num, den), (ti, ci)| {
                let shifted = &den * &Uni::new(vec![-ti.clone(), F::one()]);
                (&shifted + &num.scale(ci), num)
            },
        );
        let g = num.gcd(&den);
        let (num, den) = (&num / &g, &den / &g);
        let l = den.lc().inverse().expect("nonzero denominator");
        (num.scale(&l), den.scale(&l))
    }
}

/// `prod_l (z - v_l)`.
pub fn master<F: Field>(vals: &[F]) -> Uni<F> {
    vals.iter().fold(Uni::constant(F::one()), |acc, v| {
        &acc * &Uni::new(vec![-v.clone(), F::one()])
    })
}

/// Solves the transposed Vandermonde system `sum_l c_l v_l^(j+1) = w_j`.
///
/// For distinct nonzero `vals`, in quadratic time: the cofactor `master / (z - v_l)` annihilates every column but the
/// `l`-th. A monomial taking value `v` at a point contributes `v^j` at its `j`-th power.
pub fn solve<F: Field>(vals: &[F], master: &Uni<F>, w: &[F]) -> Vec<F> {
    vals.iter()
        .map(|v| {
            let q = master / &Uni::new(vec![-v.clone(), F::one()]);
            let num =
                q.0.iter()
                    .zip(w)
                    .fold(F::zero(), |acc, (a, b)| acc + a.clone() * b.clone());
            num / (q.eval(v) * v.clone())
        })
        .collect()
}

/// Berlekamp-Massey, one term at a time: the shortest linear recurrence generating a sequence.
#[derive(Clone, Debug)]
pub struct Massey<F> {
    seq: Vec<F>,
    c: Vec<F>,
    b: Vec<F>,
    l: usize,
    shift: usize,
    last: F,
    zeros: usize,
}

impl<F: Field> Default for Massey<F> {
    fn default() -> Self {
        Self {
            seq: Vec::new(),
            c: vec![F::one()],
            b: vec![F::one()],
            l: 0,
            shift: 1,
            last: F::one(),
            zeros: 0,
        }
    }
}

impl<F: Field> Massey<F> {
    pub fn push(&mut self, a: F) {
        self.seq.push(a);
        let k = self.seq.len() - 1;
        let d = self
            .c
            .iter()
            .zip(self.seq.iter().rev())
            .fold(F::zero(), |acc, (c, x)| acc + c.clone() * x.clone());
        if d.is_zero() {
            self.shift += 1;
            self.zeros += 1;
            return;
        }
        self.zeros = 0;
        let r = d.clone() / self.last.clone();
        let mut c = self.c.clone();
        c.resize(c.len().max(self.b.len() + self.shift), F::zero());
        for (i, x) in self.b.iter().enumerate() {
            c[i + self.shift] = c[i + self.shift].clone() - r.clone() * x.clone();
        }
        if 2 * self.l <= k {
            self.l = k + 1 - self.l;
            self.b = std::mem::replace(&mut self.c, c);
            self.last = d;
            self.shift = 1;
        } else {
            self.c = c;
            self.shift += 1;
        }
    }

    pub const fn len(&self) -> usize {
        self.seq.len()
    }

    pub const fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// The recurrence length.
    pub const fn complexity(&self) -> usize {
        self.l
    }

    /// Whether the recurrence survived the last `margin` terms and `2 l + margin` terms are in:
    /// early termination, correct with high probability over random points.
    pub const fn settled(&self, margin: usize) -> bool {
        self.zeros >= margin && self.seq.len() >= 2 * self.l + margin
    }

    /// `prod (z - r)` over the ratios `r`: the recurrence reversed.
    pub fn generator(&self) -> Uni<F> {
        let mut c = self.c.clone();
        c.resize(self.l + 1, F::zero());
        c.reverse();
        Uni::new(c)
    }
}
