# Changelog

## 0.1.0 (2026-09-26)

- Initial release.
- Operator-based `Field` trait with a blanket impl.
- `Fp` runtime prime field and `Gf<P>` const prime field.
- Bare `u64` modular arithmetic and a word-sized `Primes` iterator.
- CRT, Garner steps, and Wang rational reconstruction.
- `Monomial` and `Order` with lex, grlex, grevlex, weighted, block.
- Sparse multivariate `Poly<F>` with division certificates.
- Dense univariate `Uni<F>` with gcd, Bezout, resultant, interpolation.
- Sparse `Echelon<F>` with RREF, null space, pivot recording.
- `Ring` parser and printer for named variables.
