//! Chinese remaindering and Wang's rational number reconstruction.

use std::collections::HashMap;
use std::hash::Hash;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::modp::{inv, mul, sub, try_inv};

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

    /// Reconstruct one coefficient using the shared modulus and bound.
    pub fn reconstruct(&self, x: &BigInt) -> Option<BigRational> {
        let (mut r0, mut r1) = (self.modulus.clone(), x.mod_floor(self.modulus));
        let (mut s0, mut s1) = (BigInt::zero(), BigInt::one());
        while r1 > self.bound {
            let q = &r0 / &r1;
            let r = &r0 - &q * &r1;
            r0 = std::mem::replace(&mut r1, r);
            let s = &s0 - &q * &s1;
            s0 = std::mem::replace(&mut s1, s);
        }
        (!s1.is_zero() && s1.abs() <= self.bound && r1.gcd(&s1).is_one())
            .then(|| BigRational::new(r1, s1))
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
