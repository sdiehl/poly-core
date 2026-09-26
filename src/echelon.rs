use std::collections::{BTreeMap, BinaryHeap};

use crate::field::Field;

/// A sparse row as `(column, coefficient)` pairs, leading column first once reduced.
pub type SparseRow<F> = Vec<(usize, F)>;

/// Which end of a row leads: `Low` for F4 matrices whose columns ascend with descending
/// monomials, `High` for Laporta systems that eliminate the most complex unknown first.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Lead {
    Low,
    High,
}

/// An incremental row echelon form: monic rows keyed by their leading column, optionally recording
/// which pivots each row was reduced by, over a dense scratch row that grows on demand.
#[derive(Clone, Debug)]
pub struct Echelon<F> {
    lead: Lead,
    record: bool,
    rows: Vec<Option<SparseRow<F>>>,
    uses: Vec<Vec<usize>>,
    acc: Vec<F>,
    live: Vec<bool>,
}

impl<F: Field> Echelon<F> {
    pub const fn new(lead: Lead) -> Self {
        Self {
            lead,
            record: false,
            rows: Vec::new(),
            uses: Vec::new(),
            acc: Vec::new(),
            live: Vec::new(),
        }
    }

    /// Also keeps, for every pivot, the pivots its row was reduced by.
    pub const fn recording(lead: Lead) -> Self {
        let mut e = Self::new(lead);
        e.record = true;
        e
    }

    // An involution turning the leading column into the heap's maximum.
    const fn key(&self, c: usize) -> usize {
        match self.lead {
            Lead::High => c,
            Lead::Low => !c,
        }
    }

    fn grow(&mut self, c: usize) {
        if c >= self.acc.len() {
            self.rows.resize(c + 1, None);
            self.uses.resize(c + 1, Vec::new());
            self.acc.resize(c + 1, F::zero());
            self.live.resize(c + 1, false);
        }
    }

    pub fn rank(&self) -> usize {
        self.rows.iter().flatten().count()
    }

    /// Pivot columns, leading first.
    pub fn pivots(&self) -> Vec<usize> {
        let mut p: Vec<usize> = (0..self.rows.len())
            .filter(|&c| self.rows[c].is_some())
            .collect();
        p.sort_by_key(|&c| std::cmp::Reverse(self.key(c)));
        p
    }

    pub fn row(&self, c: usize) -> Option<&SparseRow<F>> {
        self.rows.get(c)?.as_ref()
    }

    /// The pivots row `c` was reduced by when inserted; empty unless recording.
    pub fn uses(&self, c: usize) -> &[usize] {
        self.uses.get(c).map_or(&[], Vec::as_slice)
    }

    /// Reduces `row` from its leading column down: by every pivot if `full`, otherwise only until
    /// the leading column has none. The remaining terms, leading first, and the pivots used.
    pub fn reduce(&mut self, row: &[(usize, F)], full: bool) -> (SparseRow<F>, Vec<usize>) {
        let mut heap = BinaryHeap::new();
        for (c, v) in row {
            self.grow(*c);
            self.acc[*c] = self.acc[*c].clone() + v.clone();
            if !self.live[*c] {
                self.live[*c] = true;
                heap.push(self.key(*c));
            }
        }
        let (mut out, mut used) = (SparseRow::new(), Vec::new());
        while let Some(k) = heap.pop() {
            let c = self.key(k);
            self.live[c] = false;
            let v = std::mem::replace(&mut self.acc[c], F::zero());
            match &self.rows[c] {
                _ if v.is_zero() => {}
                Some(pivot) if full || out.is_empty() => {
                    used.push(c);
                    for (j, x) in &pivot[1..] {
                        self.acc[*j] = self.acc[*j].clone() - v.clone() * x.clone();
                        if !self.live[*j] {
                            self.live[*j] = true;
                            heap.push(self.key(*j));
                        }
                    }
                }
                _ => out.push((c, v)),
            }
        }
        (out, used)
    }

    /// Adds `row` as a new monic pivot row and returns its column, or `None` if it reduced to zero.
    pub fn insert(&mut self, row: &[(usize, F)]) -> Option<usize> {
        let (mut r, used) = self.reduce(row, false);
        let c = r.first()?.0;
        let l = r[0].1.inverse()?;
        for t in &mut r {
            t.1 = t.1.clone() * l.clone();
        }
        self.rows[c] = Some(r);
        if self.record {
            self.uses[c] = used;
        }
        Some(c)
    }

    /// Column `t` as a combination of non-pivot columns (`u_t = sum a_j u_j`) and the pivots used.
    pub fn solve(&mut self, t: usize) -> (SparseRow<F>, Vec<usize>) {
        self.grow(t);
        let Some(row) = self.rows[t].clone() else {
            return (vec![(t, F::one())], Vec::new());
        };
        let (r, mut used) = self.reduce(&row[1..], true);
        used.push(t);
        (r.into_iter().map(|(j, a)| (j, -a)).collect(), used)
    }

    /// The reduced row echelon form, leading row first.
    pub fn rref(&mut self) -> Vec<SparseRow<F>> {
        let pivots = self.pivots();
        let mut out = Vec::with_capacity(pivots.len());
        for c in pivots {
            let row = self.rows[c].clone().expect("pivot");
            let mut r = vec![(c, F::one())];
            r.extend(self.reduce(&row[1..], true).0);
            out.push(r);
        }
        out
    }

    /// A basis of the null space among columns `0..ncols`, one vector per non-pivot column, each
    /// sorted by column.
    pub fn nullspace(&mut self, ncols: usize) -> Vec<SparseRow<F>> {
        self.grow(ncols.saturating_sub(1));
        let mut free: BTreeMap<usize, SparseRow<F>> = (0..ncols)
            .filter(|&j| self.rows[j].is_none())
            .map(|j| (j, vec![(j, F::one())]))
            .collect();
        for c in self.pivots() {
            for (j, a) in self.solve(c).0 {
                free.entry(j).or_default().push((c, a));
            }
        }
        free.into_values()
            .map(|mut v| {
                v.sort_by_key(|t| t.0);
                v
            })
            .collect()
    }
}
