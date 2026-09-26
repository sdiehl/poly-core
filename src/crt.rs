//! Chinese remaindering and Wang's rational number reconstruction.

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

/// Rationals from their images modulo successive primes, accepted once one further prime agrees
/// with the reconstruction. `image` may refuse an unlucky prime.
pub fn reconstruct(
    mut image: impl FnMut(u64) -> Option<Vec<u64>>,
    primes: impl IntoIterator<Item = u64>,
) -> Option<Vec<BigRational>> {
    let mut acc: Option<(Vec<BigInt>, BigInt)> = None;
    let mut candidate: Option<Vec<BigRational>> = None;
    for p in primes {
        let Some(values) = image(p) else { continue };
        let agrees = |q: &Vec<BigRational>| {
            q.len() == values.len() && q.iter().zip(&values).all(|(r, &v)| reduce(r, p) == Some(v))
        };
        if candidate.as_ref().is_some_and(agrees) {
            return candidate;
        }
        match &mut acc {
            Some((xs, m)) => garner(xs, m, &values, p),
            None => acc = Some((values.iter().map(|&v| v.into()).collect(), p.into())),
        }
        candidate = acc
            .as_ref()
            .and_then(|(xs, m)| xs.iter().map(|x| wang(x, m)).collect());
    }
    None
}
