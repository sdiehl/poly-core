//! Chinese remaindering and Wang's rational number reconstruction.

use std::collections::HashMap;
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::lehmer;
use crate::modp::{add, inv, mul, sub, try_inv};

/// `x mod p` in `[0, p)`.
pub fn residue(x: &BigInt, p: u64) -> u64 {
    x.mod_floor(&BigInt::from(p))
        .to_u64()
        .expect("residue below p")
}

/// `q` modulo `p`, or `None` when `p` divides its denominator.
pub fn reduce(q: &BigRational, p: u64) -> Option<u64> {
    let s = residue(q.denom(), p);
    (s != 0).then(|| mul(residue(q.numer(), p), inv(s, p), p))
}

/// Invalid input to a checked CRT update. Failed updates leave the accumulator unchanged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CrtError {
    LengthMismatch { expected: usize, actual: usize },
    InvalidModulus,
    NonCoprimeModuli,
    UnreducedResidue,
}

impl std::fmt::Display for CrtError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LengthMismatch { expected, actual } => {
                write!(
                    f,
                    "CRT vector length mismatch: expected {expected}, got {actual}"
                )
            }
            Self::InvalidModulus => f.write_str("CRT requires m > 0 and p > 1"),
            Self::NonCoprimeModuli => f.write_str("CRT moduli must be coprime"),
            Self::UnreducedResidue => f.write_str("CRT image residues must be less than p"),
        }
    }
}
impl std::error::Error for CrtError {}

/// Checked Garner step. Requires `m > 0`, `p > 1`, coprime moduli, equal vector
/// lengths, and image residues in `[0, p)`. Validates before changing either output.
pub fn try_garner(
    xs: &mut [BigInt],
    m: &mut BigInt,
    values: &[u64],
    p: u64,
) -> Result<(), CrtError> {
    if xs.len() != values.len() {
        return Err(CrtError::LengthMismatch {
            expected: xs.len(),
            actual: values.len(),
        });
    }
    if !m.is_positive() || p < 2 {
        return Err(CrtError::InvalidModulus);
    }
    if values.iter().any(|&v| v >= p) {
        return Err(CrtError::UnreducedResidue);
    }
    let minv = try_inv(residue(m, p), p).ok_or(CrtError::NonCoprimeModuli)?;
    for (x, &v) in xs.iter_mut().zip(values) {
        *x += &*m * mul(sub(v, residue(x, p), p), minv, p);
    }
    *m *= p;
    Ok(())
}

/// Garner's step: residues `xs` modulo `m` and `values` modulo a coprime `p`
/// become residues modulo `m * p`.
///
/// # Panics
/// Panics on invalid inputs, including mismatched vector lengths. See [`try_garner`]
/// for a non-panicking alternative and the input requirements.
pub fn garner(xs: &mut [BigInt], m: &mut BigInt, values: &[u64], p: u64) {
    try_garner(xs, m, values, p).expect("invalid Garner update");
}

/// Persistent CRT state for a fixed, caller-defined ordering of coefficients.
/// Prime selection, coefficient alignment, and candidate verification belong to the caller.
///
/// ```
/// use polycore::crt::{CrtAccumulator, WangContext};
/// let mut acc = CrtAccumulator::new(2);
/// acc.add(101, &[3, 100])?;
/// acc.add(103, &[3, 102])?;
/// let context = WangContext::new(acc.modulus()).unwrap();
/// let candidate = context.reconstruct_many(acc.residues()).unwrap();
/// assert_eq!(candidate[0], num_rational::BigRational::from_integer(3.into()));
/// assert_eq!(candidate[1], num_rational::BigRational::from_integer((-1).into()));
/// // Verify candidates independently before accepting them.
/// # Ok::<(), polycore::crt::CrtError>(())
/// ```
#[derive(Clone, Debug)]
pub struct CrtAccumulator {
    residues: Vec<BigInt>,
    modulus: BigInt,
    image_count: usize,
}

impl CrtAccumulator {
    /// Start a vector of `len` zeros modulo one, with no images accumulated.
    pub fn new(len: usize) -> Self {
        Self {
            residues: vec![BigInt::zero(); len],
            modulus: BigInt::one(),
            image_count: 0,
        }
    }

    /// Merge one image. On error, all accumulated state is preserved.
    pub fn add(&mut self, p: u64, values: &[u64]) -> Result<(), CrtError> {
        try_garner(&mut self.residues, &mut self.modulus, values, p)?;
        self.image_count += 1;
        Ok(())
    }

    pub fn residues(&self) -> &[BigInt] {
        &self.residues
    }
    pub const fn modulus(&self) -> &BigInt {
        &self.modulus
    }
    pub const fn image_count(&self) -> usize {
        self.image_count
    }

    /// Attempt reconstruction without consuming state. Success is only a candidate;
    /// the caller must verify it against independent images or exact arithmetic.
    pub fn reconstruct(&self) -> Option<Vec<BigRational>> {
        WangContext::new(&self.modulus)?.reconstruct_many(&self.residues)
    }
}

/// Persistent CRT state in Garner's mixed radix form, for primes below `2^32`.
///
/// Coefficient `i` is `d[0][i] + p[0] (d[1][i] + p[1] (d[2][i] + ...))`, one row of `u32`
/// digits per prime. Adding an image takes only machine word arithmetic, independent of how
/// many primes came before in all but a dot product, so long vectors accumulate much faster
/// than with [`CrtAccumulator`]. Big integers appear only when reconstructing.
///
/// ```
/// use polycore::crt::MixedRadixAccumulator;
/// let mut acc = MixedRadixAccumulator::new(2);
/// acc.add(101, &[3, 100])?;
/// acc.add(103, &[3, 102])?;
/// let candidate = acc.reconstruct().unwrap();
/// assert_eq!(candidate[0], num_rational::BigRational::from_integer(3.into()));
/// assert_eq!(candidate[1], num_rational::BigRational::from_integer((-1).into()));
/// # Ok::<(), polycore::crt::CrtError>(())
/// ```
#[derive(Clone, Debug)]
pub struct MixedRadixAccumulator {
    len: usize,
    primes: Vec<u64>,
    digits: Vec<Vec<u32>>,
}

impl MixedRadixAccumulator {
    // Coefficients per parallel task when adding images and reconstructing.
    const ADD_CHUNK: usize = 4096;
    const RECONSTRUCT_CHUNK: usize = 256;
    // Spread-out coefficients reconstructed first, so a premature attempt fails fast.
    const PROBES: usize = 64;

    /// Start a vector of `len` zeros modulo one, with no images accumulated.
    pub const fn new(len: usize) -> Self {
        Self {
            len,
            primes: Vec::new(),
            digits: Vec::new(),
        }
    }

    /// Merge one image modulo a prime `p < 2^32`. Moduli of `2^32` or more are an
    /// [`CrtError::InvalidModulus`]. On error, all accumulated state is preserved.
    pub fn add(&mut self, p: u64, values: &[u64]) -> Result<(), CrtError> {
        if values.len() != self.len {
            return Err(CrtError::LengthMismatch {
                expected: self.len,
                actual: values.len(),
            });
        }
        if !(2..1 << 32).contains(&p) {
            return Err(CrtError::InvalidModulus);
        }
        if values.iter().any(|&v| v >= p) {
            return Err(CrtError::UnreducedResidue);
        }
        // `x mod p` is the dot product of the digits with these weights.
        let mut weights = Vec::with_capacity(self.primes.len());
        let mut m = 1;
        for &q in &self.primes {
            weights.push(m as u32);
            m = mul(m, q % p, p);
        }
        let minv = try_inv(m, p).ok_or(CrtError::NonCoprimeModuli)?;
        let max = self.primes.iter().copied().fold(p, u64::max);
        let lazy = u128::from(max).pow(2) * (self.primes.len() as u128 + 1) < 1 << 64;
        let digits = &self.digits;
        let fill = |(c, out): (usize, &mut [u32])| {
            let span = c * Self::ADD_CHUNK..c * Self::ADD_CHUNK + out.len();
            let mut acc = vec![0u64; out.len()];
            for (row, &w) in digits.iter().zip(&weights) {
                let row = &row[span.clone()];
                if lazy {
                    for (a, &t) in acc.iter_mut().zip(row) {
                        *a += u64::from(t) * u64::from(w);
                    }
                } else {
                    for (a, &t) in acc.iter_mut().zip(row) {
                        *a = add(*a, mul(u64::from(t), u64::from(w), p), p);
                    }
                }
            }
            for ((o, &v), a) in out.iter_mut().zip(&values[span]).zip(acc) {
                *o = mul(sub(v, a % p, p), minv, p) as u32;
            }
        };
        let mut next = vec![0u32; self.len];
        #[cfg(feature = "parallel")]
        next.par_chunks_mut(Self::ADD_CHUNK)
            .enumerate()
            .for_each(fill);
        #[cfg(not(feature = "parallel"))]
        next.chunks_mut(Self::ADD_CHUNK).enumerate().for_each(fill);
        self.digits.push(next);
        self.primes.push(p);
        Ok(())
    }

    /// The primes accumulated so far, in order.
    pub fn primes(&self) -> &[u64] {
        &self.primes
    }

    pub const fn image_count(&self) -> usize {
        self.primes.len()
    }

    /// The product of the accumulated primes.
    pub fn modulus(&self) -> BigInt {
        self.primes.iter().product::<BigUint>().into()
    }

    /// Realign coefficients to a new support: entry `i` becomes the old entry
    /// `positions[i]`, or zero modulo every prime when it is `None`.
    ///
    /// # Panics
    /// Panics if a position is out of range.
    pub fn remap(&mut self, positions: &[Option<usize>]) {
        for row in &mut self.digits {
            *row = positions.iter().map(|i| i.map_or(0, |i| row[i])).collect();
        }
        self.len = positions.len();
    }

    /// Consecutive primes whose product fits a word, as one big radix digit each.
    fn radices(&self) -> Vec<(u64, std::ops::Range<usize>)> {
        let mut radices: Vec<(u64, std::ops::Range<usize>)> = Vec::new();
        for (j, &p) in self.primes.iter().enumerate() {
            match radices.last_mut() {
                Some((r, js)) if r.checked_mul(p).is_some() => {
                    *r *= p;
                    js.end = j + 1;
                }
                _ => radices.push((p, j..j + 1)),
            }
        }
        radices
    }

    fn value(&self, radices: &[(u64, std::ops::Range<usize>)], i: usize) -> BigInt {
        let mut x = BigUint::zero();
        for (r, js) in radices.iter().rev() {
            let d = js
                .clone()
                .rev()
                .fold(0, |d, j| d * self.primes[j] + u64::from(self.digits[j][i]));
            x *= *r;
            x += d;
        }
        x.into()
    }

    /// The accumulated residues in `[0, m)`.
    pub fn residues(&self) -> Vec<BigInt> {
        let radices = self.radices();
        (0..self.len).map(|i| self.value(&radices, i)).collect()
    }

    /// Attempt reconstruction without consuming state, equal to reconstructing
    /// [`Self::residues`] with a [`WangContext`]. Success is only a candidate; the caller
    /// must verify it against independent images or exact arithmetic.
    ///
    /// Neighbouring coefficients that share a denominator skip the half gcd, and a few
    /// spread-out coefficients are tried first so that a premature attempt fails fast.
    pub fn reconstruct(&self) -> Option<Vec<BigRational>> {
        if self.primes.is_empty() {
            return None;
        }
        let radices = self.radices();
        let modulus = self.modulus();
        let context = WangContext::new(&modulus)?;
        let n = self.len;
        let probes = Self::PROBES.min(n);
        for i in 0..probes {
            let index = i * (n - 1) / probes.saturating_sub(1).max(1);
            context.reconstruct(&self.value(&radices, index))?;
        }
        // A failure anywhere fails the attempt, so the other chunks stop early.
        let failed = AtomicBool::new(false);
        let run = |c: usize| -> Option<Vec<BigRational>> {
            // Once `d` holds a denominator shared with `x`, `x d` is a small integer. Both
            // parts are then within the Wang bound, so this is the fraction Wang would find.
            let mut d = BigInt::one();
            let start = c * Self::RECONSTRUCT_CHUNK;
            (start..n.min(start + Self::RECONSTRUCT_CHUNK))
                .map(|i| {
                    if failed.load(Ordering::Relaxed) {
                        return None;
                    }
                    let x = self.value(&radices, i);
                    let y = symmetric(&(&x * &d), &modulus);
                    if y.magnitude() <= context.bound.magnitude() {
                        return Some(lehmer::fraction(y, d.magnitude()));
                    }
                    let Some(q) = context.reconstruct(&x) else {
                        failed.store(true, Ordering::Relaxed);
                        return None;
                    };
                    let g = lehmer::gcd(d.magnitude(), q.denom().magnitude());
                    let lcm = &d / BigInt::from(g) * q.denom();
                    d = if lcm <= context.bound {
                        lcm
                    } else {
                        q.denom().clone()
                    };
                    Some(q)
                })
                .collect()
        };
        let chunks = 0..n.div_ceil(Self::RECONSTRUCT_CHUNK);
        #[cfg(feature = "parallel")]
        let parts = chunks.into_par_iter().map(run).collect::<Option<Vec<_>>>();
        #[cfg(not(feature = "parallel"))]
        let parts = chunks.map(run).collect::<Option<Vec<_>>>();
        parts.map(|v| v.concat())
    }
}

/// The integer in `[0, m)` with the given residues, and `m`, the product of the moduli.
pub fn crt(residues: &[(u64, u64)]) -> (BigInt, BigInt) {
    let (mut x, mut m) = ([BigInt::zero()], BigInt::one());
    for &(v, p) in residues {
        garner(&mut x, &mut m, &[v], p);
    }
    let [x] = x;
    (x, m)
}

/// The fraction `r / s` congruent to `x` modulo `m` with `|r|, |s| <= sqrt(m / 2)`, unique when it
/// exists, found by stopping the extended Euclidean algorithm halfway.
pub fn wang(x: &BigInt, m: &BigInt) -> Option<BigRational> {
    WangContext::new(m)?.reconstruct(x)
}

/// Immutable Wang reconstruction context sharing one bound across coefficients.
/// It can be borrowed by parallel workers without a parallel-runtime dependency.
#[derive(Debug)]
pub struct WangContext<'a> {
    modulus: &'a BigInt,
    bound: BigInt,
}

impl<'a> WangContext<'a> {
    /// Compute the reconstruction bound once. Returns `None` unless `modulus > 1`.
    pub fn new(modulus: &'a BigInt) -> Option<Self> {
        (modulus > &BigInt::one()).then(|| Self {
            modulus,
            bound: (modulus / 2u32).sqrt(),
        })
    }

    /// Reconstruct one coefficient using the shared modulus and bound, by Lehmer's half gcd.
    pub fn reconstruct(&self, x: &BigInt) -> Option<BigRational> {
        lehmer::wang(x, self.modulus, &self.bound)
    }

    /// Reconstruct a vector, stopping at the first coefficient that cannot be recovered.
    pub fn reconstruct_many(&self, xs: &[BigInt]) -> Option<Vec<BigRational>> {
        xs.iter().map(|x| self.reconstruct(x)).collect()
    }
}

/// `x mod m` in the symmetric range `(-m/2, m/2]`.
pub fn symmetric(x: &BigInt, m: &BigInt) -> BigInt {
    let r = x.mod_floor(m);
    if &r * 2u32 > *m { r - m } else { r }
}

/// Rationals from their images modulo successive primes, accepted once one further prime agrees
/// with the reconstruction. `image` may refuse an unlucky prime.
pub fn reconstruct(
    mut image: impl FnMut(u64) -> Option<Vec<u64>>,
    primes: impl IntoIterator<Item = u64>,
) -> Option<Vec<BigRational>> {
    reconstruct_keyed(|p| image(p).map(|v| ((), v)), primes).map(|((), q)| q)
}

/// [`reconstruct`] for images that also report a key.
///
/// The key, such as a degree or a support, is greatest at lucky primes: an image with a smaller key than the best seen is skipped as unlucky,
/// and a greater one discards what was accumulated so far.
///
/// # Panics
/// Images with the same key must have equal vector lengths and aligned coefficients.
/// Invalid CRT updates panic; see [`garner`].
pub fn reconstruct_keyed<K: Ord>(
    mut image: impl FnMut(u64) -> Option<(K, Vec<u64>)>,
    primes: impl IntoIterator<Item = u64>,
) -> Option<(K, Vec<BigRational>)> {
    let mut acc: Option<(K, Vec<BigInt>, BigInt)> = None;
    let mut candidate: Option<Vec<BigRational>> = None;
    for p in primes {
        let Some((key, values)) = image(p) else {
            continue;
        };
        match &mut acc {
            Some((best, ..)) if key < *best => continue,
            Some((best, xs, m)) if key == *best => {
                let agrees = |q: &Vec<BigRational>| {
                    q.len() == values.len()
                        && q.iter().zip(&values).all(|(r, &v)| reduce(r, p) == Some(v))
                };
                if candidate.as_ref().is_some_and(agrees) {
                    return acc.map(|(k, ..)| k).zip(candidate);
                }
                garner(xs, m, &values, p);
            }
            _ => acc = Some((key, values.iter().map(|&v| v.into()).collect(), p.into())),
        }
        candidate = acc
            .as_ref()
            .and_then(|(_, xs, m)| WangContext::new(m)?.reconstruct_many(xs));
    }
    None
}

/// Rationals from images modulo batches of primes, with majority voting over keys.
///
/// For problems whose unlucky primes give a key that is not ordered against the lucky one, such
/// as the leading monomials of a Groebner basis. Images are grouped by key, the most frequent group is reconstructed, and its candidate is
/// accepted once a later prime with that key agrees. `images` receives each batch of primes at
/// once so it can compute them in parallel, and may refuse a prime with `None`.
///
/// # Panics
/// The callback must return exactly one entry per requested prime, in order.
/// Images with the same key must have equal vector lengths and aligned coefficients.
/// Invalid CRT updates panic; see [`garner`].
pub fn reconstruct_voted<K: Clone + Eq + Hash>(
    mut images: impl FnMut(&[u64]) -> Vec<Option<(K, Vec<u64>)>>,
    primes: impl IntoIterator<Item = u64>,
) -> Option<(K, Vec<BigRational>)> {
    let mut primes = primes.into_iter();
    let mut groups: HashMap<K, (Vec<BigInt>, BigInt, usize)> = HashMap::new();
    let mut candidate: Option<(K, Vec<BigRational>)> = None;
    let mut batch = 1;
    loop {
        let ps: Vec<u64> = primes.by_ref().take(batch).collect();
        if ps.is_empty() {
            return None;
        }
        let batch_images = images(&ps);
        assert_eq!(
            batch_images.len(),
            ps.len(),
            "image callback must return one entry per prime"
        );
        for (&p, image) in ps.iter().zip(batch_images) {
            let Some((key, values)) = image else {
                continue;
            };
            if let Some((k, q)) = &candidate
                && *k == key
                && q.len() == values.len()
                && q.iter().zip(&values).all(|(r, &v)| reduce(r, p) == Some(v))
            {
                return candidate;
            }
            if let Some((xs, m, n)) = groups.get_mut(&key) {
                garner(xs, m, &values, p);
                *n += 1;
            } else {
                let xs = values.iter().map(|&v| v.into()).collect();
                groups.insert(key, (xs, p.into(), 1));
            }
        }
        candidate = groups
            .iter()
            .max_by_key(|(_, g)| g.2)
            .and_then(|(k, (xs, m, _))| {
                let context = WangContext::new(m)?;
                if let Some(x) = xs.last() {
                    context.reconstruct(x)?;
                }
                let q = context.reconstruct_many(xs)?;
                Some((k.clone(), q))
            });
        batch = (batch * 2).min(64);
    }
}
