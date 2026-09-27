//! Arithmetic in a number field `Q[a]/(m(a))`. An [`Alg`] is a polynomial in the generator
//! reduced modulo the minimal polynomial, so `Uni<Alg>` inherits division, gcd and Bezout.

use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::rc::Rc;

use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use polycore::Uni;

type Q = BigRational;

/// The field of a root of the monic irreducible `m`, with an interval `(lo, hi]` isolating it
/// when it is real.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NumberField {
    pub m: Uni<Q>,
    pub interval: Option<(Q, Q)>,
}

impl NumberField {
    /// The field of a root of an irreducible `m`, with no real approximation.
    pub fn new(m: &Uni<Q>) -> Self {
        Self {
            m: m.monic(),
            interval: None,
        }
    }
}

impl AsRef<Self> for NumberField {
    fn as_ref(&self) -> &Self {
        self
    }
}

/// An element of a number field `K`, which may carry more than its [`NumberField`].
///
/// The field is `None` only for the constants `Zero::zero` and `One::one` build before any field is known;
/// arithmetic takes the field from whichever side has one.
#[derive(Debug)]
pub struct Alg<K = NumberField> {
    pub c: Uni<Q>,
    k: Option<Rc<K>>,
}

impl<K> Clone for Alg<K> {
    fn clone(&self) -> Self {
        Self {
            c: self.c.clone(),
            k: self.k.clone(),
        }
    }
}

impl<K: AsRef<NumberField>> Alg<K> {
    pub fn new(c: Uni<Q>, k: &Rc<K>) -> Self {
        let m = &(**k).as_ref().m;
        let c = if c.deg() >= m.deg() { &c % m } else { c };
        Self {
            c,
            k: Some(k.clone()),
        }
    }

    pub fn generator(k: &Rc<K>) -> Self {
        Self::new(Uni::x(), k)
    }

    /// The field, unless this is a constant built without one.
    pub const fn field(&self) -> Option<&Rc<K>> {
        self.k.as_ref()
    }

    fn join(&self, o: &Self) -> Option<Rc<K>> {
        self.k.clone().or_else(|| o.k.clone())
    }

    fn times(&self, o: &Self) -> Self {
        let c = &self.c * &o.c;
        match self.join(o) {
            Some(k) => Self::new(c, &k),
            None => Self { c, k: None },
        }
    }

    fn inverse(&self) -> Self {
        match (&self.k, self.c.deg()) {
            (_, 0) => Self {
                c: Uni::constant(self.c.lc().recip()),
                k: self.k.clone(),
            },
            (Some(k), _) => Self::new(self.c.bezout(&(**k).as_ref().m).0, k),
            (None, _) => unreachable!("an element of no field"),
        }
    }

    pub fn rational(&self) -> Option<Q> {
        (self.c.deg() == 0).then(|| self.c.lc())
    }

    /// The trace down to Q, the sum of the conjugates. The field is passed because a constant may
    /// not know it.
    pub fn trace(&self, k: &K) -> Q {
        let p = k.as_ref().m.power_sums();
        self.c.0.iter().zip(&p).map(|(c, p)| c * p).sum()
    }

    /// The norm down to Q, the product of the conjugates.
    pub fn norm(&self, k: &K) -> Q {
        k.as_ref().m.resultant(&self.c)
    }

    /// Whether the element is positive in a real field: the interval is halved until the
    /// polynomial has no root in it, and the sign at the end is the sign at the generator.
    /// Nothing is positive in a complex field.
    pub fn positive(&self) -> bool {
        if self.c.is_zero() {
            return false;
        }
        let Some(k) = &self.k else {
            return self.c.lc().is_positive();
        };
        let k = (**k).as_ref();
        let Some((mut lo, mut hi)) = k.interval.clone() else {
            return false;
        };
        while self.c.sturm(&lo, &hi) > 0 {
            (lo, hi) = k.m.refine(&lo, &hi);
        }
        self.c.eval(&hi).is_positive()
    }
}

impl<K> From<Q> for Alg<K> {
    fn from(q: Q) -> Self {
        Self {
            c: Uni::constant(q),
            k: None,
        }
    }
}

impl<K> PartialEq for Alg<K> {
    fn eq(&self, o: &Self) -> bool {
        self.c == o.c
    }
}

impl<K: AsRef<NumberField>> Zero for Alg<K> {
    fn zero() -> Self {
        Self::from(Q::zero())
    }

    fn is_zero(&self) -> bool {
        self.c.is_zero()
    }
}

impl<K: AsRef<NumberField>> One for Alg<K> {
    fn one() -> Self {
        Self::from(Q::one())
    }
}

impl<K: AsRef<NumberField>> Add for Alg<K> {
    type Output = Self;

    fn add(self, o: Self) -> Self {
        Self {
            k: self.join(&o),
            c: &self.c + &o.c,
        }
    }
}

impl<K: AsRef<NumberField>> Sub for Alg<K> {
    type Output = Self;

    fn sub(self, o: Self) -> Self {
        Self {
            k: self.join(&o),
            c: &self.c - &o.c,
        }
    }
}

impl<K> Neg for Alg<K> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            c: -&self.c,
            k: self.k,
        }
    }
}

impl<K: AsRef<NumberField>> Mul for Alg<K> {
    type Output = Self;

    fn mul(self, o: Self) -> Self {
        self.times(&o)
    }
}

impl<K: AsRef<NumberField>> Div for Alg<K> {
    type Output = Self;

    fn div(self, o: Self) -> Self {
        self.times(&o.inverse())
    }
}
