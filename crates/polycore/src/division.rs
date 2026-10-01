//! Exact sparse division, using a packed lexicographic heap when exponents fit.

use crate::{Field, Monomial, Poly};
use std::collections::binary_heap::PeekMut;
use std::collections::{BTreeMap, BinaryHeap};

type Terms<F> = Vec<(Vec<u32>, F)>;

struct Work(Option<usize>);
impl Work {
    fn spend(&mut self, n: usize) -> Option<()> {
        if let Some(left) = &mut self.0 {
            *left = left.checked_sub(n)?;
        }
        Some(())
    }
}

fn canonical<F: Field>(n: usize, terms: &[(Vec<u32>, F)]) -> Option<Terms<F>> {
    let mut map = BTreeMap::new();
    for (e, c) in terms {
        if e.len() != n {
            return None;
        }
        let value = map.entry(e.clone()).or_insert_with(F::zero);
        *value = value.clone() + c.clone();
    }
    Some(
        map.into_iter()
            .rev()
            .filter(|(_, c)| !c.is_zero())
            .collect(),
    )
}

fn degrees<F>(n: usize, terms: &[(Vec<u32>, F)]) -> Vec<u32> {
    let mut d = vec![0; n];
    for (e, _) in terms {
        for (d, &e) in d.iter_mut().zip(e) {
            *d = (*d).max(e);
        }
    }
    d
}

fn coefficient_quotient<F: Field>(a: &F, b: &F) -> Option<F> {
    let q = a.clone() / b.clone();
    // Also supports integer coefficients: truncating division is never accepted.
    (q.clone() * b.clone() == *a).then_some(q)
}

/// Exact quotient of sparse exponent/coefficient lists, in descending lex order.
///
/// Input lists may be unsorted and contain duplicates or zeros. None means a
/// zero divisor, a dimension mismatch, a non-exact quotient, exponent overflow
/// or exhausted work. `budget` counts term reductions; None means unlimited.
/// Integer coefficients must divide exactly, even when `/` truncates.
pub fn exact_lex<F: Field>(
    n: usize,
    a: &[(Vec<u32>, F)],
    b: &[(Vec<u32>, F)],
    budget: Option<usize>,
) -> Option<Terms<F>> {
    let a = canonical(n, a)?;
    let b = canonical(n, b)?;
    if b.is_empty() {
        return None;
    }
    if a.is_empty() {
        return Some(Vec::new());
    }
    let (ad, bd) = (degrees(n, &a), degrees(n, &b));
    if bd.iter().zip(&ad).any(|(b, a)| b > a) {
        return None;
    }
    let mut work = Work(budget);
    if let Some(packing) = Packing::new(&ad) {
        packed(&a, &b, &ad, &bd, &packing, &mut work)
    } else {
        unpacked(&a, &b, &mut work)
    }
}

fn packed<F: Field>(
    a: &[(Vec<u32>, F)],
    b: &[(Vec<u32>, F)],
    ad: &[u32],
    bd: &[u32],
    pk: &Packing,
    work: &mut Work,
) -> Option<Terms<F>> {
    let bound = pk.pack(&ad.iter().zip(bd).map(|(a, b)| a - b).collect::<Vec<_>>());
    let bp: Vec<_> = b.iter().map(|(e, _)| pk.pack(e)).collect();
    let mut heap: BinaryHeap<(u64, usize, usize)> = BinaryHeap::new();
    let mut q: Vec<(u64, F)> = Vec::new();
    let mut input = a.iter().map(|(e, c)| (pk.pack(e), c)).peekable();
    loop {
        let m = match (input.peek().map(|t| t.0), heap.peek().map(|t| t.0)) {
            (None, None) => break,
            (a, b) => a.max(b).expect("one is some"),
        };
        work.spend(1)?;
        let mut c = input
            .next_if(|t| t.0 == m)
            .map_or_else(F::zero, |t| t.1.clone());
        while let Some((_, i, j)) = heap.peek_mut().filter(|t| t.0 == m).map(PeekMut::pop) {
            work.spend(1)?;
            c = c - q[i].1.clone() * b[j].1.clone();
            if j + 1 < bp.len() {
                heap.push((q[i].0 + bp[j + 1], i, j + 1));
            }
        }
        if c.is_zero() {
            continue;
        }
        let e = pk.sub(m, bp[0]).filter(|&e| pk.sub(bound, e).is_some())?;
        let c = coefficient_quotient(&c, &b[0].1)?;
        if bp.len() > 1 {
            heap.push((e + bp[1], q.len(), 1));
        }
        q.push((e, c));
    }
    Some(q.into_iter().map(|(e, c)| (pk.unpack(e), c)).collect())
}

fn unpacked<F: Field>(
    a: &[(Vec<u32>, F)],
    b: &[(Vec<u32>, F)],
    work: &mut Work,
) -> Option<Terms<F>> {
    let (lead, lc) = &b[0];
    let mut rem: BTreeMap<_, _> = a.iter().cloned().collect();
    let mut terms = Vec::new();
    while let Some((e, c)) = rem.pop_last() {
        work.spend(b.len())?;
        let q = e
            .iter()
            .zip(lead)
            .map(|(a, b)| a.checked_sub(*b))
            .collect::<Option<Vec<_>>>()?;
        let c = coefficient_quotient(&c, lc)?;
        for (e, v) in &b[1..] {
            let e = e
                .iter()
                .zip(&q)
                .map(|(e, q)| e.checked_add(*q))
                .collect::<Option<Vec<_>>>()?;
            let value = rem.get(&e).cloned().unwrap_or_else(F::zero) - c.clone() * v.clone();
            if value.is_zero() {
                rem.remove(&e);
            } else {
                rem.insert(e, value);
            }
        }
        terms.push((q, c));
    }
    Some(terms)
}

impl<F: Field> Poly<F> {
    /// Exact quotient, preserving this polynomial's monomial order.
    /// None for a zero divisor, incompatible variable counts or a remainder.
    pub fn exact(&self, divisor: &Self) -> Option<Self> {
        self.exact_impl(divisor, None)
    }

    /// Exact quotient with a limit on term reductions. None also signals
    /// exhausted work; use [`Self::exact`] when no work limit is wanted.
    pub fn exact_with_budget(&self, divisor: &Self, budget: usize) -> Option<Self> {
        self.exact_impl(divisor, Some(budget))
    }

    fn exact_impl(&self, divisor: &Self, budget: Option<usize>) -> Option<Self> {
        if self.nvars != divisor.nvars {
            return None;
        }
        let terms = |f: &Self| {
            f.terms
                .iter()
                .map(|(m, c)| (m.exps().to_vec(), c.clone()))
                .collect::<Vec<_>>()
        };
        let quotient = exact_lex(self.nvars, &terms(self), &terms(divisor), budget)?;
        Some(Self::new(
            quotient
                .into_iter()
                .map(|(e, c)| (Monomial::new(e), c))
                .collect(),
            self.nvars,
            self.order.clone(),
        ))
    }
}

/// One word, variable zero highest, with a guard bit above each exponent.
/// Integer comparison is lexicographic; overflow falls back to vector exponents.
struct Packing {
    shifts: Vec<u32>,
    guard: u64,
}

impl Packing {
    fn new(degrees: &[u32]) -> Option<Self> {
        let mut shifts = vec![0; degrees.len()];
        let (mut at, mut guard) = (0u32, 0u64);
        for (s, d) in shifts.iter_mut().zip(degrees).rev() {
            let width = 33 - d.leading_zeros();
            *s = at;
            at += width;
            if at > 64 {
                return None;
            }
            guard |= 1 << (at - 1);
        }
        Some(Self { shifts, guard })
    }

    fn pack(&self, e: &[u32]) -> u64 {
        e.iter()
            .zip(&self.shifts)
            .map(|(&x, s)| u64::from(x) << s)
            .sum()
    }

    fn unpack(&self, x: u64) -> Vec<u32> {
        let mut hi = 64;
        self.shifts
            .iter()
            .map(|&s| {
                let v = (x & (u64::MAX >> (64 - hi))) >> s;
                hi = s;
                v as u32
            })
            .collect()
    }

    fn sub(&self, a: u64, b: u64) -> Option<u64> {
        let r = (a | self.guard).wrapping_sub(b);
        (r & self.guard == self.guard).then_some(r & !self.guard)
    }
}
