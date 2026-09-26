//! Arithmetic on bare `u64` residues modulo a prime `p < 2^64`, for inner loops that keep the
//! modulus out of the element.

pub const fn add(a: u64, b: u64, p: u64) -> u64 {
    let (s, carry) = a.overflowing_add(b);
    if carry || s >= p {
        s.wrapping_sub(p)
    } else {
        s
    }
}

pub const fn sub(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a.wrapping_sub(b).wrapping_add(p)
    }
}

pub const fn neg(a: u64, p: u64) -> u64 {
    sub(0, a, p)
}

const TOP: u64 = 1 << 62;

/// For the primes just below `2^62` that [`Primes`] yields, `2^62 = c (mod p)` with `c` small, so
/// folding the high bits twice replaces the 128-bit division.
pub const fn mul(a: u64, b: u64, p: u64) -> u64 {
    let x = a as u128 * b as u128;
    let c = TOP.wrapping_sub(p);
    if c >= 1 << 20 {
        return (x % p as u128) as u64;
    }
    let y = (x >> 62) * c as u128 + (x as u64 & (TOP - 1)) as u128;
    let z = (y >> 62) as u64 * c + (y as u64 & (TOP - 1));
    if z >= p { z - p } else { z }
}

pub const fn pow(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    while e > 0 {
        if e & 1 == 1 {
            r = mul(r, a, p);
        }
        a = mul(a, a, p);
        e >>= 1;
    }
    r
}

/// The inverse of `a` modulo any `m`, or `None` unless `gcd(a, m) = 1`.
pub const fn try_inv(a: u64, m: u64) -> Option<u64> {
    let (mut r0, mut r1) = (m as i128, (a % m) as i128);
    let (mut s0, mut s1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (s0, s1) = (s1, s0 - q * s1);
    }
    if r0 == 1 {
        Some(s0.rem_euclid(m as i128) as u64)
    } else {
        None
    }
}

/// The inverse of a unit `a` modulo `p`.
pub const fn inv(a: u64, p: u64) -> u64 {
    match try_inv(a, p) {
        Some(x) => x,
        None => panic!("not invertible"),
    }
}

/// The residue as a signed integer in `(-p/2, p/2]`, which shows small integers as themselves.
pub const fn symmetric(a: u64, p: u64) -> i128 {
    if a > p / 2 {
        a as i128 - p as i128
    } else {
        a as i128
    }
}

/// Deterministic for every `u64`, by [`machine_prime`].
pub const fn is_prime(n: u64) -> bool {
    machine_prime::is_prime(n)
}

/// Primes descending from `2^62`, or from any start with [`Primes::below`].
#[derive(Clone, Debug)]
pub struct Primes(u64);

impl Primes {
    pub const fn new() -> Self {
        Self(TOP)
    }

    pub const fn below(n: u64) -> Self {
        Self(n)
    }
}

impl Default for Primes {
    fn default() -> Self {
        Self::new()
    }
}

impl Iterator for Primes {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        while self.0 > 2 {
            self.0 -= 1;
            if is_prime(self.0) {
                return Some(self.0);
            }
        }
        None
    }
}
