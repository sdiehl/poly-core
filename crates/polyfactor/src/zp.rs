//! Berlekamp's algorithm over `GF(p)`: the Frobenius map `v -> v^p` is linear on `GF(p)[x]/(f)`,
//! its fixed points are the constants on each irreducible factor, and `gcd(f, v - s)` for a fixed
//! point `v` and a residue `s` splits `f`.

use num_traits::{One, Zero};
use polycore::dense::nullspace;
use polycore::{Fp, Uni};

/// The monic irreducible factors of a squarefree monic `f` over `GF(p)`.
pub(crate) fn berlekamp(f: &Uni<Fp>, p: u64) -> Vec<Uni<Fp>> {
    let n = f.deg();
    if n <= 1 {
        return vec![f.clone()];
    }
    let frobenius = Uni::x().powmod(p, f);
    let mut row = Uni::constant(Fp::one());
    let mut m = vec![Vec::new(); n];
    for i in 0..n {
        for (j, col) in m.iter_mut().enumerate() {
            let c = row.0.get(j).copied().unwrap_or_else(Fp::zero);
            col.push(if i == j { c - Fp::one() } else { c });
        }
        row = &(&row * &frobenius) % f;
    }
    let basis = nullspace(&m);
    let mut factors = vec![f.clone()];
    for v in basis.iter().map(|v| Uni::new(v.clone())) {
        for s in 0..p {
            if factors.len() == basis.len() {
                return factors;
            }
            let shift = &v - &Uni::constant(Fp::new(s, p));
            factors = factors
                .into_iter()
                .flat_map(|g| {
                    let h = g.gcd(&shift);
                    if h.deg() > 0 && h.deg() < g.deg() {
                        vec![&g / &h, h]
                    } else {
                        vec![g]
                    }
                })
                .collect();
        }
    }
    factors
}

/// Squarefree decomposition in characteristic `p`: a vanishing derivative means `f = g(x)^p`.
fn squarefree(f: &Uni<Fp>, p: u64) -> Vec<(Uni<Fp>, u32)> {
    let mut out = Vec::new();
    let mut c = f.gcd(&f.derivative());
    let mut w = f / &c;
    let mut i = 1;
    while w.deg() > 0 {
        let y = w.gcd(&c);
        let z = &w / &y;
        if z.deg() > 0 {
            out.push((z, i));
        }
        w = y;
        c = &c / &w;
        i += 1;
    }
    if c.deg() > 0 {
        let step = usize::try_from(p).unwrap_or(usize::MAX);
        let root = Uni::new(c.0.iter().step_by(step).copied().collect());
        let e = u32::try_from(p).unwrap_or(u32::MAX);
        out.extend(squarefree(&root, p).into_iter().map(|(g, m)| (g, m * e)));
    }
    out
}

/// The leading coefficient and the monic irreducible factors over `GF(p)` with multiplicities,
/// the modulus read from the coefficients.
pub fn factor_mod(f: &Uni<Fp>) -> (Fp, Vec<(Uni<Fp>, u32)>) {
    if f.is_zero() {
        return (Fp::zero(), Vec::new());
    }
    let p = f.0.iter().map(|c| c.modulus()).max().unwrap_or(0);
    let parts = squarefree(&f.monic(), p)
        .into_iter()
        .flat_map(|(g, m)| berlekamp(&g, p).into_iter().map(move |h| (h, m)))
        .collect();
    (f.lc(), parts)
}
