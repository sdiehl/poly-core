use std::ops::{Add, Div, Mul, Neg, Rem, Sub};

use crate::field::{Field, nat, pow};
use crate::fp::Modular;
use crate::monomial::{Monomial, Order};
use crate::poly::Poly;
use crate::sample::Rng;

/// A dense univariate polynomial by ascending power, without trailing zeros.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Uni<F>(pub Vec<F>);

impl<F: Field> Uni<F> {
    pub fn new(mut v: Vec<F>) -> Self {
        while v.last().is_some_and(F::is_zero) {
            v.pop();
        }
        Self(v)
    }

    pub const fn zero() -> Self {
        Self(Vec::new())
    }

    pub fn constant(c: F) -> Self {
        Self::new(vec![c])
    }

    /// The variable itself.
    pub fn x() -> Self {
        Self(vec![F::zero(), F::one()])
    }

    pub const fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The degree, 0 for the zero polynomial.
    pub const fn deg(&self) -> usize {
        self.0.len().saturating_sub(1)
    }

    pub fn lc(&self) -> F {
        self.0.last().cloned().unwrap_or_else(F::zero)
    }

    #[must_use]
    pub fn scale(&self, c: &F) -> Self {
        Self::new(self.0.iter().map(|a| a.clone() * c.clone()).collect())
    }

    #[must_use]
    pub fn monic(&self) -> Self {
        self.lc()
            .inverse()
            .map_or_else(|| self.clone(), |l| self.scale(&l))
    }

    /// Multiplies by `x - a` in place.
    pub fn mul_linear(&mut self, a: &F) {
        if self.is_zero() {
            return;
        }
        self.0.push(F::zero());
        for i in (1..self.0.len()).rev() {
            self.0[i] = self.0[i - 1].clone() - a.clone() * self.0[i].clone();
        }
        self.0[0] = -(a.clone() * self.0[0].clone());
    }

    /// Synthetic division by `x - a`: the quotient and the remainder `self(a)`.
    pub fn deflate(&self, a: &F) -> (Self, F) {
        let Some((lead, rest)) = self.0.split_last() else {
            return (Self::zero(), F::zero());
        };
        let mut q = vec![F::zero(); rest.len()];
        let mut r = lead.clone();
        for (qi, c) in q.iter_mut().zip(rest).rev() {
            *qi = r.clone();
            r = c.clone() + a.clone() * r;
        }
        (Self(q), r)
    }

    pub fn divrem(&self, d: &Self) -> (Self, Self) {
        let l = d.lc().inverse().expect("division by the zero polynomial");
        let mut r = self.0.clone();
        if r.len() < d.0.len() {
            return (Self::zero(), self.clone());
        }
        let mut q = vec![F::zero(); r.len() - d.0.len() + 1];
        for i in (0..q.len()).rev() {
            let c = r[i + d.deg()].clone() * l.clone();
            for (j, y) in d.0.iter().enumerate() {
                r[i + j] = r[i + j].clone() - c.clone() * y.clone();
            }
            q[i] = c;
        }
        (Self::new(q), Self::new(r))
    }

    /// The quotient when `d` divides `self`.
    pub fn exact(&self, d: &Self) -> Option<Self> {
        let (q, r) = self.divrem(d);
        r.is_zero().then_some(q)
    }

    #[must_use]
    pub fn pow(&self, e: u32) -> Self {
        (0..u32::BITS - e.leading_zeros())
            .rev()
            .fold(Self::constant(F::one()), |acc, i| {
                let acc = &acc * &acc;
                if e >> i & 1 == 1 { &acc * self } else { acc }
            })
    }

    /// `self^e mod f`.
    #[must_use]
    pub fn powmod(&self, e: u64, f: &Self) -> Self {
        let base = self % f;
        (0..u64::BITS - e.leading_zeros())
            .rev()
            .fold(&Self::constant(F::one()) % f, |acc, i| {
                let acc = &(&acc * &acc) % f;
                if e >> i & 1 == 1 {
                    &(&acc * &base) % f
                } else {
                    acc
                }
            })
    }

    /// `self(g)`.
    #[must_use]
    pub fn compose(&self, g: &Self) -> Self {
        self.0.iter().rev().fold(Self::zero(), |acc, c| {
            &(&acc * g) + &Self::constant(c.clone())
        })
    }

    pub fn eval(&self, x: &F) -> F {
        self.0
            .iter()
            .rev()
            .fold(F::zero(), |acc, c| acc * x.clone() + c.clone())
    }

    #[must_use]
    pub fn derivative(&self) -> Self {
        let mut k = F::zero();
        let v = self.0.iter().skip(1).map(|a| {
            k = k.clone() + F::one();
            a.clone() * k.clone()
        });
        Self::new(v.collect())
    }

    /// The antiderivative with zero constant term. Panics in characteristic `p` at degree `p - 1`.
    #[must_use]
    pub fn integral(&self) -> Self {
        let v = self
            .0
            .iter()
            .enumerate()
            .map(|(i, a)| a.clone() / nat::<F>(i as u64 + 1));
        Self::new(std::iter::once(F::zero()).chain(v).collect())
    }

    /// The monic gcd.
    #[must_use]
    pub fn gcd(&self, o: &Self) -> Self {
        let (mut a, mut b) = (self.clone(), o.clone());
        while !b.is_zero() {
            let r = &a % &b;
            a = std::mem::replace(&mut b, r);
        }
        a.monic()
    }

    /// `(s, t, g)` with `s * self + t * o = g`, the gcd `g` monic.
    pub fn bezout(&self, o: &Self) -> (Self, Self, Self) {
        let (mut r0, mut r1) = (self.clone(), o.clone());
        let (mut s0, mut s1) = (Self::constant(F::one()), Self::zero());
        let (mut t0, mut t1) = (Self::zero(), Self::constant(F::one()));
        while !r1.is_zero() {
            let (q, r) = r0.divrem(&r1);
            r0 = std::mem::replace(&mut r1, r);
            let s = &s0 - &(&q * &s1);
            s0 = std::mem::replace(&mut s1, s);
            let t = &t0 - &(&q * &t1);
            t0 = std::mem::replace(&mut t1, t);
        }
        let l = r0.lc().inverse().unwrap_or_else(F::one);
        (s0.scale(&l), t0.scale(&l), r0.scale(&l))
    }

    /// The resultant by the Euclidean recurrence.
    pub fn resultant(&self, o: &Self) -> F {
        if self.is_zero() || o.is_zero() {
            return F::zero();
        }
        let (m, n) = (self.deg(), o.deg());
        if n == 0 {
            return pow(&o.lc(), m as u64);
        }
        let r = self % o;
        if r.is_zero() {
            return F::zero();
        }
        let sign = if m * n % 2 == 1 { -F::one() } else { F::one() };
        sign * pow(&o.lc(), (m - r.deg()) as u64) * o.resultant(&r)
    }

    /// Newton interpolation through `(xs[i], ys[i])` at distinct abscissae.
    pub fn interpolate(xs: &[F], ys: &[F]) -> Self {
        let mut out = Self::zero();
        let mut basis = Self::constant(F::one());
        for (x, y) in xs.iter().zip(ys) {
            let c = (y.clone() - out.eval(x)) / basis.eval(x);
            out = &out + &basis.scale(&c);
            basis.mul_linear(x);
        }
        out
    }

    pub fn map<G: Field>(&self, f: impl Fn(&F) -> G) -> Uni<G> {
        Uni::new(self.0.iter().map(f).collect())
    }

    /// `self` in variable `k` of an `n`-variate ring.
    pub fn to_poly(&self, k: usize, n: usize, order: Order) -> Poly<F> {
        let term = |(i, c): (usize, &F)| {
            let mut e = vec![0; n];
            e[k] = i as u32;
            (Monomial::new(e), c.clone())
        };
        Poly::new(self.0.iter().enumerate().map(term).collect(), n, order)
    }

    /// `p` as a polynomial in `x_k`, or `None` if another variable occurs.
    pub fn from_poly(p: &Poly<F>, k: usize) -> Option<Self> {
        let mut v = vec![F::zero(); p.degree(k) as usize + 1];
        for (m, c) in &p.terms {
            let others = m.exps().iter().enumerate().any(|(i, &e)| i != k && e != 0);
            if others {
                return None;
            }
            v[m.exps()[k] as usize] = c.clone();
        }
        Some(Self::new(v))
    }
}

impl<F: Modular> Uni<F> {
    /// The distinct roots in `GF(p)`, ascending by residue: the part of `gcd(f, x^p - x)` that
    /// splits into linear factors, then Cantor-Zassenhaus equal degree splitting.
    pub fn roots(&self) -> Vec<F> {
        let p = self.0.iter().map(Modular::modulus).max().unwrap_or(0);
        assert!(p != 0 || self.deg() == 0, "roots need a bound modulus");
        if self.deg() == 0 {
            return Vec::new();
        }
        let lift = |v: u64| F::from_residue(v, p);
        let mut out: Vec<F> = if p == 2 {
            (0..2)
                .map(lift)
                .filter(|a| self.eval(a).is_zero())
                .collect()
        } else {
            let f = self.monic();
            let x = Self::new(vec![lift(0), lift(1)]);
            let g = f.gcd(&(&x.powmod(p, &f) - &x));
            let mut out = Vec::new();
            split(&g, p, &mut Rng::new(p), &mut out);
            out
        };
        out.sort_by_key(|a| a.residue_mod(p));
        out
    }
}

fn split<F: Modular>(g: &Uni<F>, p: u64, rng: &mut Rng, out: &mut Vec<F>) {
    match g.deg() {
        0 => {}
        1 => out.push(-(g.0[0] / g.0[1])),
        _ => loop {
            let lin = Uni::new(vec![
                F::from_residue(rng.nonzero(p), p),
                F::from_residue(1, p),
            ]);
            let h = g.gcd(&(&lin.powmod((p - 1) / 2, g) - &Uni::constant(F::one())));
            if 0 < h.deg() && h.deg() < g.deg() {
                let rest = g / &h;
                split(&h, p, rng, out);
                split(&rest, p, rng, out);
                return;
            }
        },
    }
}

impl<F: Field> Add for &Uni<F> {
    type Output = Uni<F>;
    fn add(self, o: Self) -> Uni<F> {
        let at = |v: &[F], i: usize| v.get(i).cloned().unwrap_or_else(F::zero);
        let n = self.0.len().max(o.0.len());
        Uni::new((0..n).map(|i| at(&self.0, i) + at(&o.0, i)).collect())
    }
}

impl<F: Field> Neg for &Uni<F> {
    type Output = Uni<F>;
    fn neg(self) -> Uni<F> {
        Uni(self.0.iter().map(|c| -c.clone()).collect())
    }
}

impl<F: Field> Sub for &Uni<F> {
    type Output = Uni<F>;
    fn sub(self, o: Self) -> Uni<F> {
        self + &-o
    }
}

impl<F: Field> Mul for &Uni<F> {
    type Output = Uni<F>;
    fn mul(self, o: Self) -> Uni<F> {
        if self.is_zero() || o.is_zero() {
            return Uni::zero();
        }
        let mut v = vec![F::zero(); self.0.len() + o.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in o.0.iter().enumerate() {
                v[i + j] = v[i + j].clone() + a.clone() * b.clone();
            }
        }
        Uni::new(v)
    }
}

impl<F: Field> Div for &Uni<F> {
    type Output = Uni<F>;
    fn div(self, d: Self) -> Uni<F> {
        self.divrem(d).0
    }
}

impl<F: Field> Rem for &Uni<F> {
    type Output = Uni<F>;
    fn rem(self, d: Self) -> Uni<F> {
        self.divrem(d).1
    }
}
