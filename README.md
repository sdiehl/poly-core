# poly-core

A shared polynomial trait foundation for fields, word-sized prime fields, Chinese remaindering, sparse multivariate and dense univariate polynomials, and a sparse incremental echelon form.

- [`Field`](src/field.rs): blanket trait over the field operators.
- [`Fp`](src/fp.rs), [`Gf<P>`](src/fp.rs): prime fields below `2^64`.
- [`Poly<F>`](src/poly.rs), [`Monomial`](src/monomial.rs), [`Order`](src/monomial.rs): sparse multivariate polynomials.
- [`Uni<F>`](src/uni.rs): dense univariate polynomials.
- [`RatFunc<F>`](src/ratfunc.rs): univariate rational functions.
- [`Echelon<F>`](src/echelon.rs), [`dense`](src/dense.rs): sparse and dense linear algebra.
- [`crt`](src/crt.rs): Chinese remaindering and rational reconstruction.
- [`interp`](src/interp.rs): Newton, Thiele, Berlekamp-Massey.
- [`Ring`](src/parse.rs): parsing and printing.

`Field` is operator-based rather than method-based (`add`, `multiply`, ...) so that `BigRational`, `Fp`, and extensions can be fields without wrappers or impls, generic code reads as the mathematics does, and nothing new is needed from `num-traits`. The cost is some `clone()` calls on big coefficients but whatever.

`Fp::zero()` and `Fp::one()` carry no modulus and adopt the one they meet, so nullary constructors still work at a runtime modulus.

**This library is for scientific computing uses not cryptography or anything timing attack sensitive.**

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
