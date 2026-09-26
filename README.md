# poly-core

A shared polynomial trait foundation for fields, word-sized prime fields, Chinese remaindering, sparse multivariate and dense univariate polynomials, and a sparse incremental echelon form.

- [`Field`](src/field.rs): any `Clone + PartialEq + Debug + Zero + One + Neg + Sub + Div` type, by a blanket impl.
- [`Fp`](src/fp.rs) (runtime prime below `2^64`), [`Gf<P>`](src/fp.rs) (const prime), [`modp`](src/modp.rs) bare `u64` arithmetic, [`Primes`](src/modp.rs) descending from `2^62`.
- [`crt`](src/crt.rs): Garner steps, `crt`, Wang rational reconstruction, and `reconstruct` which lifts images until one more prime agrees.
- [`Monomial`](src/monomial.rs), [`Order`](src/monomial.rs) (lex, grlex, grevlex, weighted, block), sparse [`Poly<F>`](src/poly.rs) with arithmetic, evaluation, S-polynomials, `divide` by a list and `combination` to check its certificate.
- [`Uni<F>`](src/uni.rs): divrem, gcd, Bezout, resultant, derivative, Newton interpolation.
- [`Echelon<F>`](src/echelon.rs) over [`SparseRow<F>`](src/echelon.rs): insert, reduce, solve, RREF, null space, low or high leading column, optional pivot recording.
- [`Ring`](src/parse.rs): parse and print polynomials with named variables.

`Field` is operator-based rather than method-based (`add`, `multiply`, ...) so that `BigRational`, `Fp`, and extensions can be fields without wrappers or impls, generic code reads as the mathematics does, and nothing new is needed from `num-traits`. The cost is some `clone()` calls on big coefficients but whatever.

`Fp::zero()` and `Fp::one()` carry no modulus and adopt the one they meet, so nullary constructors still work at a runtime modulus.

**This library is for scientific computing uses not cryptography or anything timing attack sensitive.**

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
