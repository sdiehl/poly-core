//! Cached and geometric evaluation of sparse polynomials.

use crate::field::pow;
use crate::modp::mul;
use crate::{Field, Monomial, Poly, Uni};
use std::collections::BTreeMap;

/// Bare-residue powers `x[i]^j` for `j <= degrees[i]`.
/// Panics for mismatched dimensions or a modulus less than two.
pub fn power_table(x: &[u64], degrees: &[usize], p: u64) -> Vec<Vec<u64>> {
    assert!(p >= 2, "modulus must be at least two");
    assert_eq!(x.len(), degrees.len(), "point and degree dimensions differ");
    x.iter()
        .zip(degrees)
        .map(|(&x, &d)| {
            std::iter::successors(Some(1), |&v| Some(mul(v, x % p, p)))
                .take(d.checked_add(1).expect("power table too large"))
                .collect()
        })
        .collect()
}

/// Reusable powers at one point. Construction allocates `sum(degrees[i] + 1)` entries.
#[derive(Clone, Debug)]
pub struct PowerTable<F> {
    powers: Vec<Vec<F>>,
}

impl<F: Field> PowerTable<F> {
    pub fn new(point: &[F], degrees: &[u32]) -> Self {
        assert_eq!(
            point.len(),
            degrees.len(),
            "point and degree dimensions differ"
        );
        let powers = point
            .iter()
            .zip(degrees)
            .map(|(x, &d)| {
                std::iter::successors(Some(F::one()), |v| Some(v.clone() * x.clone()))
                    .take(
                        usize::try_from(d)
                            .expect("degree fits usize")
                            .checked_add(1)
                            .expect("power table too large"),
                    )
                    .collect()
            })
            .collect();
        Self { powers }
    }

    /// Evaluate a monomial within the cached degree box.
    pub fn monomial(&self, exponents: &[u32]) -> F {
        assert_eq!(
            exponents.len(),
            self.powers.len(),
            "monomial dimension differs"
        );
        exponents
            .iter()
            .zip(&self.powers)
            .fold(F::one(), |v, (&e, powers)| v * powers[e as usize].clone())
    }
}

impl<F: Field> Poly<F> {
    /// Degrees in every variable in a single pass through the support.
    pub fn degrees(&self) -> Vec<u32> {
        let mut degrees = vec![0; self.nvars];
        for (m, _) in &self.terms {
            for (d, &e) in degrees.iter_mut().zip(m.exps()) {
                *d = (*d).max(e);
            }
        }
        degrees
    }

    /// Evaluate using powers shared with other polynomials at the same point.
    pub fn eval_cached(&self, table: &PowerTable<F>) -> F {
        assert_eq!(self.nvars, table.powers.len(), "point dimension differs");
        self.terms.iter().fold(F::zero(), |sum, (m, c)| {
            sum + c.clone() * table.monomial(m.exps())
        })
    }

    /// Leave `x_k` free while using cached powers in the other variables.
    /// The cache need not contain powers above zero for variable `k`.
    pub fn eval_except_cached(&self, k: usize, table: &PowerTable<F>) -> Uni<F> {
        assert!(k < self.nvars, "variable out of bounds");
        assert_eq!(self.nvars, table.powers.len(), "point dimension differs");
        let mut coefficients = vec![F::zero(); self.degree(k) as usize + 1];
        for (m, c) in &self.terms {
            let value = m
                .exps()
                .iter()
                .enumerate()
                .filter(|&(i, _)| i != k)
                .fold(c.clone(), |v, (i, &e)| {
                    v * table.powers[i][e as usize].clone()
                });
            let coefficient = &mut coefficients[m.exps()[k] as usize];
            *coefficient = coefficient.clone() + value;
        }
        Uni::new(coefficients)
    }

    /// All leave-one-variable-free images in O(n T + sum(degrees)) operations.
    /// Prefix/suffix products also work when evaluation coordinates are zero.
    pub fn univariate_images(&self, point: &[F]) -> Vec<Uni<F>> {
        let degrees = self.degrees();
        let table = PowerTable::new(point, &degrees);
        let mut images: Vec<_> = degrees
            .iter()
            .map(|&d| vec![F::zero(); d as usize + 1])
            .collect();
        let mut prefix = vec![F::one(); self.nvars + 1];
        for (m, c) in &self.terms {
            for (i, &e) in m.exps().iter().enumerate() {
                prefix[i + 1] = prefix[i].clone() * table.powers[i][e as usize].clone();
            }
            let mut suffix = c.clone();
            for (i, &e) in m.exps().iter().enumerate().rev() {
                let slot = &mut images[i][e as usize];
                *slot = slot.clone() + prefix[i].clone() * suffix.clone();
                suffix = suffix * table.powers[i][e as usize].clone();
            }
        }
        images.into_iter().map(Uni::new).collect()
    }
}

/// Successive evaluations at `scale[i] * ratios[i]^j`, starting at j = 1,
/// retaining the first `keep` variables. Each sample updates every term once.
#[derive(Clone, Debug)]
pub struct GeometricEvaluator<F> {
    keep: usize,
    order: crate::Order,
    groups: BTreeMap<Vec<u32>, Vec<(F, F)>>,
}

impl<F: Field> GeometricEvaluator<F> {
    pub fn new(f: &Poly<F>, keep: usize, scale: &[F], ratios: &[F]) -> Self {
        assert!(keep <= f.nvars, "too many retained variables");
        assert_eq!(scale.len(), f.nvars - keep, "scale dimension differs");
        assert_eq!(ratios.len(), scale.len(), "ratio dimension differs");
        let mut groups: BTreeMap<Vec<u32>, Vec<(F, F)>> = BTreeMap::new();
        for (m, c) in &f.terms {
            let at = |point: &[F]| {
                m.exps()[keep..]
                    .iter()
                    .zip(point)
                    .fold(F::one(), |v, (&e, x)| v * pow(x, u64::from(e)))
            };
            groups
                .entry(m.exps()[..keep].to_vec())
                .or_default()
                .push((c.clone() * at(scale), at(ratios)));
        }
        Self {
            keep,
            order: f.order.clone(),
            groups,
        }
    }

    /// Next image, in the retained variables and the input monomial order.
    pub fn advance(&mut self) -> Poly<F> {
        let terms = self
            .groups
            .iter_mut()
            .map(|(e, terms)| {
                let c = terms.iter_mut().fold(F::zero(), |sum, (value, ratio)| {
                    *value = value.clone() * ratio.clone();
                    sum + value.clone()
                });
                (Monomial::new(e.clone()), c)
            })
            .collect();
        Poly::new(terms, self.keep, self.order.clone())
    }
}
