//! Deterministic sampling and black boxes: every point is named by a key, so reruns and snapshots
//! reproduce and callers that share a key share evaluations.

use rand_xoshiro::rand_core::{Rng as _, SeedableRng};
use rand_xoshiro::{SplitMix64, Xoshiro256PlusPlus};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// A hash of `key`, for deriving independent seeds from one.
pub fn hash(key: &[u64]) -> u64 {
    key.iter()
        .fold(0, |h, &k| SplitMix64::seed_from_u64(h ^ k).next_u64())
}

/// The nonzero residue modulo `p` named by `key`.
pub fn point(key: &[u64], p: u64) -> u64 {
    1 + hash(key) % (p - 1)
}

/// `xoshiro256++`, seeded through `SplitMix64`.
#[derive(Clone, Debug)]
pub struct Rng(Xoshiro256PlusPlus);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(Xoshiro256PlusPlus::seed_from_u64(seed))
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0.next_u64()
    }

    /// Uniform in `[1, p)`, up to a bias of `p / 2^64`.
    pub fn nonzero(&mut self, p: u64) -> u64 {
        1 + self.next_u64() % (p - 1)
    }
}

/// A polynomial or rational function known only through its values modulo `p`. `None` marks a
/// point where evaluation broke down and which should be avoided.
pub trait BlackBox {
    fn eval(&self, x: &[u64], p: u64) -> Option<u64>;

    /// Independent points at once, which an implementation may spread over threads.
    fn eval_many(&self, xs: &[Vec<u64>], p: u64) -> Vec<Option<u64>> {
        xs.iter().map(|x| self.eval(x, p)).collect()
    }
}

/// Closures are black boxes, evaluated in parallel with the `parallel` feature.
impl<F: Fn(&[u64], u64) -> Option<u64> + Sync> BlackBox for F {
    fn eval(&self, x: &[u64], p: u64) -> Option<u64> {
        self(x, p)
    }

    #[cfg(feature = "parallel")]
    fn eval_many(&self, xs: &[Vec<u64>], p: u64) -> Vec<Option<u64>> {
        xs.par_iter().map(|x| self(x, p)).collect()
    }
}
