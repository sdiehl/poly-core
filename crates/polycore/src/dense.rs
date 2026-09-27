//! Dense Gaussian elimination over any field, for the small square systems that sparse methods
//! bottom out in. Matrices are rows of equal length.

use std::fmt;

use crate::field::Field;

/// Why [`solve`] found no unique solution.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolveError {
    /// Some solution exists but it is not unique.
    Deficient,
    /// No solution exists.
    Inconsistent,
}

impl fmt::Display for SolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Deficient => "rank deficient",
            Self::Inconsistent => "inconsistent",
        })
    }
}

impl std::error::Error for SolveError {}

/// Reduces `a` in place to reduced row echelon form and returns the pivot columns.
pub fn rref<F: Field>(a: &mut [Vec<F>]) -> Vec<usize> {
    let ncols = a.first().map_or(0, Vec::len);
    let mut pivots = Vec::new();
    for col in 0..ncols {
        let r = pivots.len();
        let Some(i) = (r..a.len()).find(|&i| !a[i][col].is_zero()) else {
            continue;
        };
        a.swap(r, i);
        let l = a[r][col].inverse().expect("nonzero pivot");
        a[r].iter_mut().for_each(|x| *x = x.clone() * l.clone());
        let pivot = a[r].clone();
        for (k, row) in a.iter_mut().enumerate() {
            let c = row[col].clone();
            if k != r && !c.is_zero() {
                for (x, y) in row.iter_mut().zip(&pivot) {
                    *x = x.clone() - c.clone() * y.clone();
                }
            }
        }
        pivots.push(col);
        if pivots.len() == a.len() {
            break;
        }
    }
    pivots
}

/// A basis of `{x : a x = 0}`.
pub fn nullspace<F: Field>(a: &[Vec<F>]) -> Vec<Vec<F>> {
    let ncols = a.first().map_or(0, Vec::len);
    let mut m = a.to_vec();
    let pivots = rref(&mut m);
    (0..ncols)
        .filter(|c| !pivots.contains(c))
        .map(|free| {
            let mut x = vec![F::zero(); ncols];
            x[free] = F::one();
            for (row, &p) in m.iter().zip(&pivots) {
                x[p] = -row[free].clone();
            }
            x
        })
        .collect()
}

/// The unique `x` with `a x = b`.
pub fn solve<F: Field>(a: &[Vec<F>], b: &[F]) -> Result<Vec<F>, SolveError> {
    let n = a.first().map_or(0, Vec::len);
    let mut m: Vec<Vec<F>> = a
        .iter()
        .zip(b)
        .map(|(row, y)| row.iter().cloned().chain([y.clone()]).collect())
        .collect();
    let pivots = rref(&mut m);
    if pivots.last() == Some(&n) {
        return Err(SolveError::Inconsistent);
    }
    if pivots.len() < n {
        return Err(SolveError::Deficient);
    }
    Ok(m.into_iter().take(n).map(|row| row[n].clone()).collect())
}

/// The inverse of a square `a`, or `None` if it is singular.
pub fn invert<F: Field>(a: &[Vec<F>]) -> Option<Vec<Vec<F>>> {
    let n = a.len();
    if a.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut m: Vec<Vec<F>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let unit = (0..n).map(|j| if i == j { F::one() } else { F::zero() });
            row.iter().cloned().chain(unit).collect()
        })
        .collect();
    let pivots = rref(&mut m);
    (pivots.len() == n && pivots.last() < Some(&n))
        .then(|| m.into_iter().map(|row| row[n..].to_vec()).collect())
}
