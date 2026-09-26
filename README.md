# poly-core

A shared polynomial trait foundation for fields, word-sized prime fields, Chinese remaindering, sparse multivariate and dense univariate polynomials, and a sparse incremental echelon form.

- `Field`: any `Clone + PartialEq + Debug + Zero + One + Neg + Sub + Div` type, by a blanket impl.
- `Fp` (runtime prime below `2^64`), `Gf<P>` (const prime), `modp` bare `u64` arithmetic, `Primes` descending from `2^62`.
- `crt`: Garner steps, `crt`, Wang rational reconstruction, and `reconstruct` which lifts images until one more prime agrees.
- `Monomial`, `Order` (lex, grlex, grevlex, weighted, block), sparse `Poly<F>` with arithmetic, evaluation, S-polynomials, `divide` by a list and `combination` to check its certificate.
- `Uni<F>`: divrem, gcd, Bezout, resultant, derivative, Newton interpolation.
- `Echelon<F>` over `SparseRow<F>`: insert, reduce, solve, RREF, null space, low or high leading column, optional pivot recording.
- `Ring`: parse and print polynomials with named variables.

`Field` is operator-based rather than method-based (`add`, `multiply`, ...) so that `BigRational`, `Fp`, and extensions can be fields without wrappers or impls, generic code reads as the mathematics does, and nothing new is needed from `num-traits`. The cost is some `clone()` calls on big coefficients but whatever.

`Fp::zero()` and `Fp::one()` carry no modulus and adopt the one they meet, so nullary constructors still work at a runtime modulus.

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
