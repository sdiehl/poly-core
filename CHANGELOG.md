# Changelog

## 0.1.1 (2026-09-27)

- Split into a workspace of `polycore` and `polyfactor`.
- Add `polyfactor` for factoring over GF(p), Q, and number fields.
- Add `Uni::squarefree` and `Uni::power_sums`.
- Add Sturm counting, isolation and refinement on `Uni<Q>`.
- `Poly::modp` and `Uni::modp` reduce rational polynomials modulo a prime; `residues` and `symmetric` read them back.
- `Uni::mul_linear` multiplies by `x - a` in place.
- `Uni::deflate` performs synthetic division by `x - a`, returning the quotient and `self(a)`.
- `Fp` arithmetic and equality take a fast path when both operands share a modulus.
- Newton, Thiele, `master`, `solve`, and `Uni::interpolate` use in-place linear factors instead of full polynomial multiplication and division.

## 0.1.0 (2026-09-26)

- Initial release.
- Operator-based `Field` trait with a blanket impl.
- `Fp` runtime prime field and `Gf<P>` const prime field.
- Bare `u64` modular arithmetic and a word-sized `Primes` iterator.
- CRT, Garner steps, and Wang rational reconstruction.
- Persistent `crt::CrtAccumulator` with checked, non-consuming updates and reconstruction.
- `crt::WangContext` shares the reconstruction bound across scalar, batched, or caller-parallel reconstruction.
- `crt::try_garner` validates vector lengths, moduli, and image residues before mutation; `garner` rejects invalid inputs explicitly.
- Reconstruction helpers reuse Wang bounds and reject inconsistent image and batch lengths.
- `Monomial` and `Order` with lex, grlex, grevlex, weighted, block.
- Sparse multivariate `Poly<F>` with division certificates.
- Dense univariate `Uni<F>` with gcd, Bezout, resultant, interpolation.
- Sparse `Echelon<F>` with RREF, null space, pivot recording.
- `Ring` parser and printer for named variables, with validation, parameters, lists and LaTeX.
- Primality by `machine-prime`, `Modular` trait over `Fp` and `Gf<P>`.
- `Uni` powers, composition, integral, exact division and roots over `GF(p)`.
- `Poly` powers, derivatives, permutations, line and partial evaluation, coefficients in one variable.
- `RatFunc<F>`, content and primitive parts over Q, keyed reconstruction.
- `interp` (Newton, Thiele, Berlekamp-Massey, Vandermonde), `dense` linear algebra, `sample` and `BlackBox`.
