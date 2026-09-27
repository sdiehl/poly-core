//! Factoring over a number field by Trager's algorithm. For a squarefree `g` and a shift `s` that
//! makes the norm `N(x) = Res_y(m(y), g(x - s*y))` squarefree, the factors of `g` are
//! `gcd(N_j(x), g(x - s*a))` shifted back, over the irreducible factors `N_j` of `N` over Q.

use std::fmt::Debug;
use std::rc::Rc;

use num_rational::BigRational;
use num_traits::One;
use polycore::Uni;

use crate::field::{Alg, NumberField};
use crate::rational::factor;

type Q = BigRational;

/// `p(x + c)`.
fn shift<K: AsRef<NumberField> + Debug>(p: &Uni<Alg<K>>, c: &Alg<K>) -> Uni<Alg<K>> {
    p.compose(&Uni::new(vec![c.clone(), Alg::one()]))
}

/// The norm of `g`, the product of its conjugates, by interpolation at `deg g * deg m + 1`
/// points.
fn norm<K: AsRef<NumberField> + Debug>(g: &Uni<Alg<K>>, k: &K) -> Uni<Q> {
    let xs: Vec<Q> = (0..=g.deg() * k.as_ref().m.deg())
        .map(|i| Q::from_integer(i.into()))
        .collect();
    let ys: Vec<Q> = xs
        .iter()
        .map(|x| g.eval(&Alg::from(x.clone())).norm(k))
        .collect();
    Uni::interpolate(&xs, &ys)
}

/// The monic irreducible factors of a squarefree `g` over the field.
fn trager<K: AsRef<NumberField> + Debug>(g: &Uni<Alg<K>>, k: &Rc<K>) -> Vec<Uni<Alg<K>>> {
    if g.deg() <= 1 {
        return vec![g.monic()];
    }
    let a = Alg::generator(k);
    for s in 0i64.. {
        let sa = a.clone() * Alg::from(Q::from_integer(s.into()));
        let h = shift(g, &-sa.clone());
        let n = norm(&h, k);
        if n.gcd(&n.derivative()).deg() > 0 {
            continue;
        }
        let parts = factor(&n);
        if parts.len() == 1 {
            return vec![g.monic()];
        }
        return parts
            .iter()
            .map(|(p, _)| shift(&p.map(|c| Alg::from(c.clone())).gcd(&h), &sa))
            .collect();
    }
    unreachable!("some shift makes the norm squarefree")
}

/// The monic irreducible factors over the field `k` with multiplicities.
pub fn factor_over<K: AsRef<NumberField> + Debug>(
    f: &Uni<Alg<K>>,
    k: &Rc<K>,
) -> Vec<(Uni<Alg<K>>, u32)> {
    f.squarefree()
        .into_iter()
        .flat_map(|(g, m)| trager(&g, k).into_iter().map(move |h| (h, m)))
        .collect()
}
