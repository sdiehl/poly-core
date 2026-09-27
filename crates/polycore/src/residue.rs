//! Rational polynomials read modulo a prime, and residues read back as integers, either in
//! `[0, p)` or nearest zero.

use num_rational::BigRational;

use crate::fp::Fp;
use crate::modp;
use crate::poly::Poly;
use crate::uni::Uni;

fn residue(c: Fp) -> BigRational {
    BigRational::from_integer(c.value().into())
}

fn symmetric(c: Fp) -> BigRational {
    BigRational::from_integer(modp::symmetric(c.value(), c.modulus()).into())
}

impl Poly<BigRational> {
    /// `self` modulo `p`, or `None` when `p` divides a denominator.
    pub fn modp(&self, p: u64) -> Option<Poly<Fp>> {
        self.try_map(|c| Fp::from_rational(c, p))
    }
}

impl Uni<BigRational> {
    /// `self` modulo `p`, or `None` when `p` divides a denominator.
    pub fn modp(&self, p: u64) -> Option<Uni<Fp>> {
        let cs = self.0.iter().map(|c| Fp::from_rational(c, p));
        cs.collect::<Option<_>>().map(Uni::new)
    }
}

impl Poly<Fp> {
    /// The residues as the integers `0..p`.
    pub fn residues(&self) -> Poly<BigRational> {
        self.map(|c| residue(*c))
    }

    /// The residues as the integers in `(-p/2, p/2]`.
    pub fn symmetric(&self) -> Poly<BigRational> {
        self.map(|c| symmetric(*c))
    }
}

impl Uni<Fp> {
    /// The residues as the integers `0..p`.
    pub fn residues(&self) -> Uni<BigRational> {
        self.map(|c| residue(*c))
    }

    /// The residues as the integers in `(-p/2, p/2]`.
    pub fn symmetric(&self) -> Uni<BigRational> {
        self.map(|c| symmetric(*c))
    }
}
