use polycore::modp::{self, MulBy};

// Knuth's MMIX linear congruential multiplier, for a reproducible residue stream.
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;

#[test]
fn mul_and_mul_by_match_wide_remainder() {
    let mut x = 1u64;
    for p in [
        3,
        32_003,
        (1 << 32) + 15,
        (1 << 62) - 57,
        (1 << 62) + 135,
        (1 << 63) - 25,
        (1 << 63) + 29,
        u64::MAX - 58,
    ] {
        assert!(modp::is_prime(p), "{p}");
        let edges = [0, 1, 2, p / 2, p - 2, p - 1];
        for _ in 0..2000 {
            x = x.wrapping_mul(LCG_MULTIPLIER).wrapping_add(1);
            let a = [x % p, edges[(x >> 60) as usize % edges.len()]];
            x = x.wrapping_mul(LCG_MULTIPLIER).wrapping_add(1);
            let b = [x % p, edges[(x >> 60) as usize % edges.len()]];
            for (a, b) in a.into_iter().zip(b) {
                let expected = u128::from(a) * u128::from(b) % u128::from(p);
                assert_eq!(
                    u128::from(modp::mul(a, b, p)),
                    expected,
                    "{a} * {b} mod {p}"
                );
                assert_eq!(
                    u128::from(MulBy::new(a, p).mul(b, p)),
                    expected,
                    "{a} * {b} mod {p}"
                );
            }
        }
    }
}
