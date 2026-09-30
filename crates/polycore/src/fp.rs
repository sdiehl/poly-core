use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::str::FromStr;

use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::field::Field;
use crate::{crt, modp};

/// A prime field of word-sized residues, for code that drops to bare `u64` arithmetic in
/// [`modp`] on its inner loops.
pub trait Modular: Field + Copy {
    /// The prime, or 0 for a constant not yet bound to one.
    fn modulus(&self) -> u64;
    /// The stored residue; a signed integer in two's complement while unbound.
    fn residue(&self) -> u64;
    fn from_residue(v: u64, p: u64) -> Self;

    /// The residue in `[0, p)`, binding an unbound constant to `p`.
    fn residue_mod(&self, p: u64) -> u64 {
        if self.modulus() == 0 {
            modp::from_signed(self.residue() as i64, p)
        } else {
            self.residue()
        }
    }
}

/// An element of `GF(p)` for a prime `p < 2^64` chosen at runtime.
///
/// `Fp::zero()` and `Fp::one()` carry no modulus (`p = 0`) and adopt the modulus of whatever they
/// meet, so generic code that starts from `F::zero()` or `F::one()` works unchanged.
#[derive(Clone, Copy, Debug)]
pub struct Fp {
    v: u64,
    p: u64,
}

impl Fp {
    /// `v` modulo `p`, or an unbound constant when `p = 0`.
    pub const fn new(v: u64, p: u64) -> Self {
        Self {
            v: if p == 0 { v } else { v % p },
            p,
        }
    }

    pub const fn from_i64(x: i64, p: u64) -> Self {
        Self {
            v: modp::from_signed(x, p),
            p,
        }
    }

    /// `q` modulo `p`, or `None` when `p` divides its denominator.
    pub fn from_rational(q: &BigRational, p: u64) -> Option<Self> {
        crt::reduce(q, p).map(|v| Self { v, p })
    }

    pub const fn value(self) -> u64 {
        self.v
    }

    /// The prime, or 0 for a constant not yet bound to one.
    pub const fn modulus(self) -> u64 {
        self.p
    }

    /// The inverse, or `None` for a non-unit: zero, a residue sharing a factor with a composite
    /// modulus, or an unbound constant other than `1` and `-1`.
    pub const fn try_inverse(self) -> Option<Self> {
        match (self.p, self.v) {
            (0, 1 | u64::MAX) => Some(self),
            (0, _) => None,
            (p, v) => match modp::try_inv(v, p) {
                Some(v) => Some(Self { v, p }),
                None => None,
            },
        }
    }

    #[must_use]
    pub const fn bind(self, p: u64) -> Self {
        if self.p == 0 && p != 0 {
            Self::from_i64(self.v as i64, p)
        } else {
            self
        }
    }

    fn unify(self, o: Self) -> (u64, u64, u64) {
        let p = self.p.max(o.p);
        (self.bind(p).v, o.bind(p).v, p)
    }

    #[inline]
    fn lift(self, o: Self, wrap: fn(u64, u64) -> u64, op: fn(u64, u64, u64) -> u64) -> Self {
        if self.p == o.p && self.p != 0 {
            return Self {
                v: op(self.v, o.v, self.p),
                p: self.p,
            };
        }
        self.lift_cold(o, wrap, op)
    }

    #[cold]
    fn lift_cold(self, o: Self, wrap: fn(u64, u64) -> u64, op: fn(u64, u64, u64) -> u64) -> Self {
        let (a, b, p) = self.unify(o);
        let v = if p == 0 { wrap(a, b) } else { op(a, b, p) };
        Self { v, p }
    }
}

impl PartialEq for Fp {
    fn eq(&self, o: &Self) -> bool {
        if self.p == o.p {
            return self.v == o.v;
        }
        let (a, b, _) = self.unify(*o);
        a == b
    }
}

impl Eq for Fp {}

/// Hashes the symmetric residue, so equal elements of one field hash alike and an unbound
/// constant hashes like its binding whenever it is below `p / 2` in magnitude.
impl Hash for Fp {
    fn hash<H: Hasher>(&self, h: &mut H) {
        let v = if self.p == 0 {
            i128::from(self.v as i64)
        } else {
            modp::symmetric(self.v, self.p)
        };
        v.hash(h);
    }
}

impl From<i32> for Fp {
    /// An unbound constant.
    fn from(x: i32) -> Self {
        Self {
            v: i64::from(x) as u64,
            p: 0,
        }
    }
}

impl Modular for Fp {
    fn modulus(&self) -> u64 {
        self.p
    }
    fn residue(&self) -> u64 {
        self.v
    }
    fn from_residue(v: u64, p: u64) -> Self {
        Self::new(v, p)
    }
}

impl Add for Fp {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        self.lift(o, u64::wrapping_add, modp::add)
    }
}

impl Sub for Fp {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self.lift(o, u64::wrapping_sub, modp::sub)
    }
}

impl Mul for Fp {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        self.lift(o, u64::wrapping_mul, modp::mul)
    }
}

impl Div for Fp {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        assert!(!o.is_zero(), "division by zero in GF(p)");
        let unit = |a: u64, b: u64| {
            assert!(b == 1 || b == u64::MAX, "division of unbound constants");
            a.wrapping_mul(b)
        };
        self.lift(o, unit, |a, b, p| modp::mul(a, modp::inv(b, p), p))
    }
}

impl Neg for Fp {
    type Output = Self;
    fn neg(self) -> Self {
        Self::zero() - self
    }
}

impl Zero for Fp {
    fn zero() -> Self {
        Self { v: 0, p: 0 }
    }
    fn is_zero(&self) -> bool {
        self.v == 0
    }
}

impl One for Fp {
    fn one() -> Self {
        Self { v: 1, p: 0 }
    }
}

impl fmt::Display for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.p == 0 {
            write!(f, "{}", self.v as i64)
        } else {
            write!(f, "{}", self.v)
        }
    }
}

/// An element of `GF(P)` for a prime `P` fixed at compile time.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Gf<const P: u64>(u64);

impl<const P: u64> Gf<P> {
    pub const fn new(v: u64) -> Self {
        Self(v % P)
    }

    pub const fn from_i64(x: i64) -> Self {
        Self(modp::from_signed(x, P))
    }

    pub const fn value(self) -> u64 {
        self.0
    }

    pub const fn modulus(self) -> u64 {
        P
    }
}

impl<const P: u64> Modular for Gf<P> {
    fn modulus(&self) -> u64 {
        P
    }
    fn residue(&self) -> u64 {
        self.0
    }
    fn from_residue(v: u64, _: u64) -> Self {
        Self::new(v)
    }
}

impl<const P: u64> Add for Gf<P> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self(modp::add(self.0, o.0, P))
    }
}

impl<const P: u64> Sub for Gf<P> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self(modp::sub(self.0, o.0, P))
    }
}

impl<const P: u64> Mul for Gf<P> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self(modp::mul(self.0, o.0, P))
    }
}

impl<const P: u64> Div for Gf<P> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        assert!(o.0 != 0, "division by zero in GF(p)");
        Self(modp::mul(self.0, modp::inv(o.0, P), P))
    }
}

impl<const P: u64> Neg for Gf<P> {
    type Output = Self;
    fn neg(self) -> Self {
        Self(modp::neg(self.0, P))
    }
}

impl<const P: u64> Zero for Gf<P> {
    fn zero() -> Self {
        Self(0)
    }
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl<const P: u64> One for Gf<P> {
    fn one() -> Self {
        Self(1 % P)
    }
}

impl<const P: u64> From<u64> for Gf<P> {
    fn from(v: u64) -> Self {
        Self::new(v)
    }
}

impl<const P: u64> From<u32> for Gf<P> {
    fn from(v: u32) -> Self {
        Self::new(v.into())
    }
}

impl<const P: u64> From<i32> for Gf<P> {
    fn from(x: i32) -> Self {
        Self::from_i64(x.into())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParseFpError;

impl fmt::Display for ParseFpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("expected an optionally signed integer")
    }
}

impl std::error::Error for ParseFpError {}

/// Reduces as it reads, so any number of digits is accepted.
impl<const P: u64> FromStr for Gf<P> {
    type Err = ParseFpError;
    fn from_str(s: &str) -> Result<Self, ParseFpError> {
        let s = s.trim();
        let (neg, digits) = s
            .strip_prefix('-')
            .map_or_else(|| (false, s.strip_prefix('+').unwrap_or(s)), |d| (true, d));
        if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
            return Err(ParseFpError);
        }
        let v = digits.bytes().fold(Self::new(0), |acc, b| {
            acc * Self::new(10) + Self::new(u64::from(b - b'0'))
        });
        Ok(if neg { -v } else { v })
    }
}

impl<const P: u64> fmt::Display for Gf<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
