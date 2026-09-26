//! Chinese remaindering and Wang's rational number reconstruction.

use std::collections::HashMap;
use std::hash::Hash;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::modp::{inv, mul, sub};

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

/// Garner's step: residues `xs` modulo `m` and `values` modulo a prime `p` coprime to `m` become
/// residues modulo `m * p`.
pub fn garner(xs: &mut [BigInt], m: &mut BigInt, values: &[u64], p: u64) {
    let minv = inv(residue(m, p), p);
    for (x, &v) in xs.iter_mut().zip(values) {
        *x += &*m * mul(sub(v, residue(x, p), p), minv, p);
    }
    *m *= p;
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
    let bound = (m / 2u32).sqrt();
    let (mut r0, mut r1) = (m.clone(), x.mod_floor(m));
    let (mut s0, mut s1) = (BigInt::zero(), BigInt::one());
    while r1 > bound {
        let q = &r0 / &r1;
        let r = &r0 - &q * &r1;
        r0 = std::mem::replace(&mut r1, r);
        let s = &s0 - &q * &s1;
        s0 = std::mem::replace(&mut s1, s);
    }
    (!s1.is_zero() && s1.abs() <= bound && r1.gcd(&s1).is_one()).then(|| BigRational::new(r1, s1))
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
            .and_then(|(_, xs, m)| xs.iter().map(|x| wang(x, m)).collect());
    }
    None
}

/// Rationals from images modulo batches of primes, with majority voting over keys.
///
/// For problems whose unlucky primes give a key that is not ordered against the lucky one, such
/// as the leading monomials of a Groebner basis. Images are grouped by key, the most frequent group is reconstructed, and its candidate is
/// accepted once a later prime with that key agrees. `images` receives each batch of primes at
/// once so it can compute them in parallel, and may refuse a prime with `None`.
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
        for (&p, image) in ps.iter().zip(images(&ps)) {
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
                if let Some(x) = xs.last() {
                    wang(x, m)?;
                }
                let q = xs.iter().map(|x| wang(x, m)).collect::<Option<_>>()?;
                Some((k.clone(), q))
            });
        batch = (batch * 2).min(64);
    }
}
