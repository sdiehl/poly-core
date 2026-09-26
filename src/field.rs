use std::fmt::Debug;
use std::ops::{Div, Neg, Sub};

use num_traits::{One, Zero};

/// A field given by its operators. Every type with these operations is one, so `BigRational`,
/// [`crate::Fp`], [`crate::Gf`] and user extensions such as `Q(sqrt(m))` need no wrapper.
pub trait Field:
    Clone
    + PartialEq
    + Debug
    + Zero
    + One
    + Neg<Output = Self>
    + Sub<Output = Self>
    + Div<Output = Self>
{
    /// The multiplicative inverse, or `None` for zero.
    fn inverse(&self) -> Option<Self> {
        (!self.is_zero()).then(|| Self::one() / self.clone())
    }
}

impl<F> Field for F where
    F: Clone + PartialEq + Debug + Zero + One + Neg<Output = F> + Sub<Output = F> + Div<Output = F>
{
}

/// `x^e` by repeated squaring.
pub fn pow<F: Field>(x: &F, mut e: u64) -> F {
    let (mut base, mut acc) = (x.clone(), F::one());
    while e > 0 {
        if e & 1 == 1 {
            acc = acc * base.clone();
        }
        e >>= 1;
        if e > 0 {
            base = base.clone() * base;
        }
    }
    acc
}
