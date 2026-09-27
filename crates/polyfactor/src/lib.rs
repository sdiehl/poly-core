//! Factoring polynomials in one variable over `GF(p)`, over Q, and over number fields, on top of
//! [`polycore`].
#![allow(
    clippy::many_single_char_names,
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::redundant_pub_crate
)]

mod field;
mod rational;
mod trager;
mod zp;

pub use field::{Alg, NumberField};
pub use rational::{factor, is_irreducible};
pub use trager::factor_over;
pub use zp::factor_mod;
