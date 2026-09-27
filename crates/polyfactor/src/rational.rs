//! Factoring over Q: the factors modulo a prime are lifted by Hensel to a modulus beyond the
//! Mignotte bound on the coefficients of any factor, and the products of lifted factors that
//! divide the polynomial are its factors (Zassenhaus).

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Zero, pow};
use polycore::modp::is_prime;
use polycore::{Field, Fp, Uni, crt};

use crate::zp::berlekamp;

type Q = BigRational;

/// The most factors modulo one prime that recombination will search.
const MAX_MODULAR: usize = 20;

fn modp(f: &Uni<Q>, p: u64) -> Uni<Fp> {
    f.modp(p).expect("a coefficient prime to the modulus")
}

/// Integer coefficients folded into `(-m/2, m/2]`.
fn symmetric(f: &Uni<Q>, m: &BigInt) -> Uni<Q> {
    f.map(|c| Q::from_integer(crt::symmetric(&c.to_integer(), m)))
}

fn product<F: Field>(fs: impl Iterator<Item = Uni<F>>) -> Uni<F> {
    fs.fold(Uni::constant(F::one()), |acc, f| &acc * &f)
}

/// The Mignotte bound on the coefficients of any factor of a monic `f`.
fn mignotte(f: &Uni<Q>) -> BigInt {
    let norm: BigInt = f.0.iter().map(|c| c.to_integer().pow(2)).sum();
    pow(BigInt::from(2), f.deg()) * (norm.sqrt() + 1)
}

/// Factors of `f` modulo `p` lifted to factors modulo `m > 2 * bound`.
fn hensel(f: &Uni<Q>, gs: &[Uni<Fp>], p: u64, bound: &BigInt) -> (Vec<Uni<Q>>, BigInt) {
    let cofactors: Vec<Uni<Fp>> = gs
        .iter()
        .enumerate()
        .map(|(i, g)| {
            let others = product(
                gs.iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .map(|(_, h)| h.clone()),
            );
            &others.bezout(g).0 % g
        })
        .collect();
    let mut lifted: Vec<Uni<Q>> = gs.iter().map(Uni::residues).collect();
    let mut m = BigInt::from(p);
    while m <= bound * 2 {
        let scale = Q::from_integer(m.clone());
        let error = modp(
            &(f - &product(lifted.iter().cloned())).scale(&scale.recip()),
            p,
        );
        for ((g, g0), s) in lifted.iter_mut().zip(gs).zip(&cofactors) {
            let delta = &(&error * s) % g0;
            *g = &*g + &delta.residues().scale(&scale);
        }
        m *= p;
    }
    (lifted.iter().map(|g| symmetric(g, &m)).collect(), m)
}

/// The products of lifted factors that divide `f` over Z, by size.
fn zassenhaus(f: &Uni<Q>, lifted: Vec<Uni<Q>>, m: &BigInt) -> Vec<Uni<Q>> {
    let (mut f, mut pool, mut out) = (f.clone(), lifted, Vec::new());
    let mut size = 1;
    'search: while 2 * size <= pool.len() {
        let n = pool.len();
        for mask in (1..1usize << n).filter(|k| k.count_ones() as usize == size) {
            let h = symmetric(
                &product(
                    (0..n)
                        .filter(|i| mask >> i & 1 == 1)
                        .map(|i| pool[i].clone()),
                ),
                m,
            );
            let plausible = !h.0[0].is_zero() && (&f.0[0] / &h.0[0]).is_integer();
            let (q, r) = f.divrem(&h);
            if plausible && r.is_zero() {
                f = q;
                out.push(h);
                pool = (0..n)
                    .filter(|i| mask >> i & 1 == 0)
                    .map(|i| pool[i].clone())
                    .collect();
                continue 'search;
            }
        }
        size += 1;
    }
    if f.deg() > 0 {
        out.push(f);
    }
    out
}

/// The monic irreducible factors over Q of a squarefree `f`.
fn irreducible(f: &Uni<Q>) -> Vec<Uni<Q>> {
    if f.deg() <= 1 {
        return vec![f.monic()];
    }
    if f.0[0].is_zero() {
        let x = Uni::x();
        let mut out = irreducible(&(f / &x));
        out.push(x);
        return out;
    }
    let ints = f.primitive();
    let (a, n) = (ints.lc(), ints.deg());
    let monic = Uni::new(
        ints.0
            .iter()
            .enumerate()
            .map(|(i, c)| c * pow(a.clone(), n - i) / &a)
            .collect(),
    );
    let mut best: Option<(u64, Vec<Uni<Fp>>)> = None;
    let images = (2..).filter(|&p| is_prime(p)).map(|p| (p, modp(&monic, p)));
    for (p, fp) in images
        .filter(|(_, fp)| fp.gcd(&fp.derivative()).deg() == 0)
        .take(3)
    {
        let gs = berlekamp(&fp, p);
        let done = gs.len() == 1;
        if best.as_ref().is_none_or(|(_, b)| gs.len() < b.len()) {
            best = Some((p, gs));
        }
        if done {
            break;
        }
    }
    let Some((p, gs)) = best.filter(|(_, gs)| gs.len() > 1 && gs.len() <= MAX_MODULAR) else {
        return vec![f.monic()];
    };
    let (lifted, m) = hensel(&monic, &gs, p, &mignotte(&monic));
    zassenhaus(&monic, lifted, &m)
        .iter()
        .map(|h| {
            Uni::new(
                h.0.iter()
                    .enumerate()
                    .map(|(i, c)| c * pow(a.clone(), i))
                    .collect(),
            )
            .monic()
        })
        .collect()
}

/// The monic irreducible factors over Q with multiplicities.
pub fn factor(f: &Uni<Q>) -> Vec<(Uni<Q>, u32)> {
    f.squarefree()
        .into_iter()
        .flat_map(|(g, m)| irreducible(&g).into_iter().map(move |h| (h, m)))
        .collect()
}

/// Whether `f` is irreducible over Q, of positive degree.
pub fn is_irreducible(f: &Uni<Q>) -> bool {
    matches!(factor(f).as_slice(), [(_, 1)])
}
