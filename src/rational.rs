//! Content and primitive parts over Q: the primitive part has coprime integer coefficients and a
//! positive leading coefficient, and the content is the rational that scales it back.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::poly::Poly;
use crate::uni::Uni;

fn content<'a>(
    cs: impl Iterator<Item = &'a BigRational>,
    lead: Option<&BigRational>,
) -> BigRational {
    let (g, l) = cs.fold((BigInt::zero(), BigInt::one()), |(g, l), c| {
        (g.gcd(c.numer()), l.lcm(c.denom()))
    });
    let c = BigRational::new(g, l);
    match lead {
        Some(a) if a.is_negative() => -c,
        _ => c,
    }
}

impl Poly<BigRational> {
    /// Zero for the zero polynomial.
    pub fn content(&self) -> BigRational {
        content(self.terms.iter().map(|t| &t.1), self.lc())
    }

    #[must_use]
    pub fn primitive(&self) -> Self {
        let c = self.content();
        if c.is_zero() {
            self.clone()
        } else {
            self.scale(&c.recip())
        }
    }
}

impl Uni<BigRational> {
    /// Zero for the zero polynomial.
    pub fn content(&self) -> BigRational {
        content(self.0.iter(), self.0.last())
    }

    #[must_use]
    pub fn primitive(&self) -> Self {
        let c = self.content();
        if c.is_zero() {
            self.clone()
        } else {
            self.scale(&c.recip())
        }
    }
}
