//! Power-of-two subgroups of word-sized prime fields.

use crate::modp::{inv, is_prime, mul, pow};
use crate::sample::Rng;

/// Descending primes below a bound with `v_2(p - 1) = bits`.
#[derive(Clone, Debug)]
pub struct SmoothPrimes {
    next: Option<u64>,
    step: u64,
}

impl SmoothPrimes {
    /// Primes below 2^62, or None for a bit count outside 1..=62.
    pub fn new(bits: u32) -> Option<Self> {
        Self::below(1u64 << 62, bits)
    }

    /// Primes strictly below `upper`. None for a bit count outside 1..=62.
    /// An empty interval yields an empty iterator.
    pub fn below(upper: u64, bits: u32) -> Option<Self> {
        if !(1..=62).contains(&bits) {
            return None;
        }
        let order = 1u64 << bits;
        let quotient = upper.saturating_sub(2) / order;
        let odd = if quotient.is_multiple_of(2) {
            quotient.checked_sub(1)
        } else {
            Some(quotient)
        };
        let next = odd.map(|c| c * order + 1).filter(|&p| p < upper);
        Some(Self {
            next,
            step: 2 * order,
        })
    }
}

impl Iterator for SmoothPrimes {
    type Item = u64;
    fn next(&mut self) -> Option<u64> {
        while let Some(p) = self.next {
            self.next = p.checked_sub(self.step).filter(|&q| q > 2);
            if is_prime(p) {
                return Some(p);
            }
        }
        None
    }
}

/// A generator of exact order `2^bits`, with binary Pohlig–Hellman logarithms.
#[derive(Clone, Copy, Debug)]
pub struct PowerOfTwoSubgroup {
    p: u64,
    root: u64,
    bits: u32,
}

impl PowerOfTwoSubgroup {
    /// Validate the prime, subgroup order and generator. Allows the trivial subgroup.
    pub const fn from_generator(p: u64, root: u64, bits: u32) -> Option<Self> {
        if bits >= 64 || !is_prime(p) {
            return None;
        }
        let order = 1u64 << bits;
        let root = root % p;
        if !(p - 1).is_multiple_of(order)
            || pow(root, order, p) != 1
            || (bits > 0 && pow(root, order / 2, p) == 1)
        {
            return None;
        }
        Some(Self { p, root, bits })
    }

    /// Sample a generator, making at most 32 attempts. None for invalid
    /// parameters or exhausted sampling; never returns a smaller-order generator.
    pub fn new(p: u64, bits: u32, rng: &mut Rng) -> Option<Self> {
        if bits >= 64 || !is_prime(p) || !(p - 1).is_multiple_of(1u64 << bits) {
            return None;
        }
        let order = 1u64 << bits;
        (0..32).find_map(|_| {
            let root = pow(rng.nonzero(p), (p - 1) / order, p);
            (bits == 0 || pow(root, order / 2, p) != 1).then_some(Self { p, root, bits })
        })
    }

    pub const fn modulus(self) -> u64 {
        self.p
    }
    pub const fn generator(self) -> u64 {
        self.root
    }
    pub const fn order(self) -> u64 {
        1u64 << self.bits
    }
    pub const fn pow(self, exponent: u64) -> u64 {
        pow(self.root, exponent, self.p)
    }

    /// The unique exponent in `[0, order)`, or None outside this subgroup.
    pub fn log(self, value: u64) -> Option<u64> {
        let mut value = value % self.p;
        let mut inverse = inv(self.root, self.p);
        let mut exponent = 0;
        for bit in 0..self.bits {
            match pow(value, 1u64 << (self.bits - bit - 1), self.p) {
                1 => {}
                minus_one if minus_one == self.p - 1 => {
                    exponent |= 1u64 << bit;
                    value = mul(value, inverse, self.p);
                }
                _ => return None,
            }
            inverse = mul(inverse, inverse, self.p);
        }
        (value == 1).then_some(exponent)
    }
}
