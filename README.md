# polycore

## polycore

A shared polynomial trait foundation for fields, word-sized prime fields, Chinese remaindering, sparse multivariate and dense univariate polynomials, and a sparse incremental echelon form.

- [`Field`](crates/polycore/src/field.rs): blanket trait over the field operators
- [`Modular`](crates/polycore/src/fp.rs): prime fields with access to bare residues
- [`Fp`](crates/polycore/src/fp.rs), [`Gf<P>`](crates/polycore/src/fp.rs): prime fields below $2^{64}$
- [`Poly<F>`](crates/polycore/src/poly.rs), [`Monomial`](crates/polycore/src/monomial.rs), [`Order`](crates/polycore/src/monomial.rs): sparse multivariate polynomials
- [`Uni<F>`](crates/polycore/src/uni.rs): dense univariate polynomials
- [`PrimeField`](crates/polycore/src/fast.rs): NTT multiplication, Newton division and half-GCD
- [`PowerTable<F>`](crates/polycore/src/evaluation.rs): cached monomial powers
- [`GeometricEvaluator<F>`](crates/polycore/src/evaluation.rs): successive geometric evaluations
- [`modp`](crates/polycore/src/modp.rs): word-sized modular arithmetic
- [`MulBy`](crates/polycore/src/modp.rs): Shoup multiplication by a fixed residue
- [`Primes`](crates/polycore/src/modp.rs): descending word-sized primes
- [`SmoothPrimes`](crates/polycore/src/subgroup.rs): primes with a prescribed power of two in p − 1
- [`PowerOfTwoSubgroup`](crates/polycore/src/subgroup.rs): subgroup generators and discrete logarithms
- [`RatFunc<F>`](crates/polycore/src/ratfunc.rs): univariate rational functions
- [`Echelon<F>`](crates/polycore/src/echelon.rs), [`dense`](crates/polycore/src/dense.rs): sparse and dense linear algebra
- [`crt`](crates/polycore/src/crt.rs): Chinese remaindering and rational reconstruction
- [`CrtAccumulator`](crates/polycore/src/crt.rs): incremental CRT over coefficient vectors
- [`MixedRadixAccumulator`](crates/polycore/src/crt.rs): CRT with coefficient support remapping
- [`WangContext`](crates/polycore/src/crt.rs): rational reconstruction with a shared bound
- [`lehmer`](crates/polycore/src/lehmer.rs): Lehmer gcd for big integers
- [`Newton<F>`](crates/polycore/src/interp.rs): incremental polynomial interpolation
- [`Thiele<F>`](crates/polycore/src/interp.rs): rational interpolation by continued fractions
- [`Massey<F>`](crates/polycore/src/interp.rs): incremental Berlekamp–Massey recurrence recovery
- [`Rng`](crates/polycore/src/sample.rs): deterministic sampling from keyed seeds
- [`BlackBox`](crates/polycore/src/sample.rs): evaluation of unknown functions over prime fields
- [`Ring`](crates/polycore/src/parse.rs): parsing and printing

`Field` is operator-based rather than method-based (`add`, `multiply`, ...) so that `BigRational`, `Fp`, and extensions can be fields without wrappers or impls, generic code reads as the mathematics does, and nothing new is needed from `num-traits`. The cost is some `clone()` calls on big coefficients but whatever.

`Fp::zero()` and `Fp::one()` carry no modulus and adopt the one they meet, so nullary constructors still work at a runtime modulus.

**This library is for scientific computing uses not cryptography or anything timing attack sensitive.**

## polyfactor

Factoring in one variable, on top of polycore.

- [`factor_mod`](crates/polyfactor/src/zp.rs): Berlekamp over $\mathrm{GF}(p)$
- [`factor`](crates/polyfactor/src/rational.rs): Hensel lifting and Zassenhaus over $\mathbb{Q}$
- [`Alg`, `NumberField`](crates/polyfactor/src/field.rs): number field arithmetic
- [`factor_over`](crates/polyfactor/src/trager.rs): Trager over number fields

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
