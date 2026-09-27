//! Real roots of polynomials over Q: counted by Sturm sequences and isolated by bisection of the
//! Cauchy bound.

use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::uni::Uni;

type Q = BigRational;

impl Uni<Q> {
    /// The number of distinct real roots in the half-open `(a, b]`.
    pub fn sturm(&self, a: &Q, b: &Q) -> usize {
        let p = self / &self.gcd(&self.derivative());
        let d = p.derivative();
        let mut chain = vec![p, d];
        while !chain[chain.len() - 1].is_zero() {
            let (f, g) = (&chain[chain.len() - 2], &chain[chain.len() - 1]);
            chain.push(-&(f % g));
        }
        let changes = |x: &Q| {
            let signs: Vec<bool> = chain
                .iter()
                .map(|f| f.eval(x))
                .filter(|v| !v.is_zero())
                .map(|v| v.is_positive())
                .collect();
            signs.windows(2).filter(|w| w[0] != w[1]).count()
        };
        changes(a).saturating_sub(changes(b))
    }

    /// Isolating intervals `(a, b]` for the real roots, ascending.
    pub fn isolate(&self) -> Vec<(Q, Q)> {
        if self.deg() == 0 {
            return Vec::new();
        }
        let lc = self.lc();
        let bound = self
            .0
            .iter()
            .fold(Q::one(), |m, c| m.max((c / &lc).abs() + Q::one()));
        let mut out = Vec::new();
        let mut todo = vec![(-bound.clone(), bound)];
        while let Some((lo, hi)) = todo.pop() {
            match self.sturm(&lo, &hi) {
                0 => {}
                1 => out.push((lo, hi)),
                _ => {
                    let mid = (&lo + &hi) / Q::from_integer(2.into());
                    todo.push((lo, mid.clone()));
                    todo.push((mid, hi));
                }
            }
        }
        out.sort();
        out
    }

    /// The half of an isolating interval `(lo, hi]` that keeps the root.
    pub fn refine(&self, lo: &Q, hi: &Q) -> (Q, Q) {
        let mid = (lo + hi) / Q::from_integer(2.into());
        if self.sturm(lo, &mid) == 1 {
            (lo.clone(), mid)
        } else {
            (mid, hi.clone())
        }
    }
}
