//! Arithmetic in a simple extension `F[a]/(m(a))`: a number field over Q, or one step of a tower
//! over another. An [`Alg`] is a polynomial in the generator reduced modulo the minimal
//! polynomial, so `Uni<Alg>` inherits division, gcd and Bezout.

use std::ops::{Add, Div, Mul, Neg, Sub};
use std::rc::Rc;

use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use polycore::{Field, Uni};

type Q = BigRational;

/// The field of a root of the monic irreducible `m` over `F`, with an interval `(lo, hi]`
/// isolating it when it is real over Q.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NumberField<F = Q> {
    pub m: Uni<F>,
    pub interval: Option<(Q, Q)>,
}

impl<F: Field> NumberField<F> {
    /// The field of a root of an irreducible `m`, with no real approximation.
    pub fn new(m: &Uni<F>) -> Self {
        Self {
            m: m.monic(),
            interval: None,
        }
    }
}

impl<F> AsRef<Self> for NumberField<F> {
    fn as_ref(&self) -> &Self {
        self
    }
}

/// An element of a number field `K` over `F`, which may carry more than its [`NumberField`].
///
/// The field is `None` only for elements of `F` built before any field is known, such as
/// `Zero::zero` and `One::one`; arithmetic takes the field from whichever side has one.
#[derive(Debug)]
pub struct Alg<K = NumberField, F = Q> {
    pub c: Uni<F>,
    k: Option<Rc<K>>,
}

impl<K, F: Clone> Clone for Alg<K, F> {
    fn clone(&self) -> Self {
        Self {
            c: self.c.clone(),
            k: self.k.clone(),
        }
    }
}

impl<F: Field, K: AsRef<NumberField<F>>> Alg<K, F> {
    pub fn new(c: Uni<F>, k: &Rc<K>) -> Self {
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

    /// The field, unless this is an element of `F` built without one.
    pub const fn field(&self) -> Option<&Rc<K>> {
        self.k.as_ref()
    }

    /// The coordinate at the `i`-th power of the generator.
    pub fn coordinate(&self, i: usize) -> F {
        self.c.0.get(i).cloned().unwrap_or_else(F::zero)
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
                c: Uni::constant(F::one() / self.c.lc()),
                k: self.k.clone(),
            },
            (Some(k), _) => Self::new(self.c.bezout(&(**k).as_ref().m).0, k),
            (None, _) => unreachable!("an element of no field"),
        }
    }

    /// The element as one of `F`, when it is.
    pub fn rational(&self) -> Option<F> {
        (self.c.deg() == 0).then(|| self.c.lc())
    }

    /// The trace down to `F`, the sum of the conjugates. The field is passed because an element
    /// of `F` may not know it.
    pub fn trace(&self, k: &K) -> F {
        let p = k.as_ref().m.power_sums();
        self.c
            .0
            .iter()
            .zip(&p)
            .fold(F::zero(), |s, (c, p)| s + c.clone() * p.clone())
    }

    /// The norm down to `F`, the product of the conjugates.
    pub fn norm(&self, k: &K) -> F {
        k.as_ref().m.resultant(&self.c)
    }
}

impl<K: AsRef<NumberField>> Alg<K> {
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

impl<K, F: Field> From<F> for Alg<K, F> {
    fn from(a: F) -> Self {
        Self {
            c: Uni::constant(a),
            k: None,
        }
    }
}

impl<K, F: PartialEq> PartialEq for Alg<K, F> {
    fn eq(&self, o: &Self) -> bool {
        self.c == o.c
    }
}

impl<F: Field, K: AsRef<NumberField<F>>> Zero for Alg<K, F> {
    fn zero() -> Self {
        Self::from(F::zero())
    }

    fn is_zero(&self) -> bool {
        self.c.is_zero()
    }
}

impl<F: Field, K: AsRef<NumberField<F>>> One for Alg<K, F> {
    fn one() -> Self {
        Self::from(F::one())
    }
}

impl<F: Field, K: AsRef<NumberField<F>>> Add for Alg<K, F> {
    type Output = Self;

    fn add(self, o: Self) -> Self {
        Self {
            k: self.join(&o),
            c: &self.c + &o.c,
        }
    }
}

impl<F: Field, K: AsRef<NumberField<F>>> Sub for Alg<K, F> {
    type Output = Self;

    fn sub(self, o: Self) -> Self {
        Self {
            k: self.join(&o),
            c: &self.c - &o.c,
        }
    }
}

impl<F: Field, K> Neg for Alg<K, F> {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            c: -&self.c,
            k: self.k,
        }
    }
}

impl<F: Field, K: AsRef<NumberField<F>>> Mul for Alg<K, F> {
    type Output = Self;

    fn mul(self, o: Self) -> Self {
        self.times(&o)
    }
}

impl<F: Field, K: AsRef<NumberField<F>>> Div for Alg<K, F> {
    type Output = Self;

    fn div(self, o: Self) -> Self {
        self.times(&o.inverse())
    }
}
