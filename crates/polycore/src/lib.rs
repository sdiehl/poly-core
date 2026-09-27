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
pub mod dense;
mod echelon;
mod field;
mod fp;
pub mod interp;
pub mod modp;
mod monomial;
mod parse;
mod poly;
mod ratfunc;
mod rational;
mod real;
mod residue;
pub mod sample;
mod uni;

pub use echelon::{Echelon, Lead, SparseRow};
pub use field::{Field, nat, pow};
pub use fp::{Fp, Gf, Modular, ParseFpError};
pub use modp::Primes;
pub use monomial::{Monomial, Order};
pub use parse::{ParseError, Ring};
pub use poly::{Poly, Term, combination};
pub use ratfunc::RatFunc;
pub use uni::Uni;
