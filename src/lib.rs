//! Fields, prime fields, sparse and dense polynomials, and sparse echelon forms: the shared
//! foundation beneath Groebner bases, sparse interpolation and linear system reduction.
#![allow(
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::module_name_repetitions,
    clippy::must_use_candidate,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

pub mod crt;
mod echelon;
mod field;
mod fp;
pub mod modp;
mod monomial;
mod parse;
mod poly;
mod uni;

pub use echelon::{Echelon, Lead, SparseRow};
pub use field::{pow, Field};
pub use fp::{Fp, Gf};
pub use modp::Primes;
pub use monomial::{Monomial, Order};
pub use parse::{ParseError, Ring};
pub use poly::{combination, Poly, Term};
pub use uni::Uni;
