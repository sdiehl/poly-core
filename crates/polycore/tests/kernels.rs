use num_bigint::BigInt;
use polycore::evaluation::{GeometricEvaluator, PowerTable, power_table};
use polycore::fast::PrimeField;
use polycore::modp::{PowerOfTwoSubgroup, SmoothPrimes};
use polycore::sample::Rng;
use polycore::{Fp, Gf, Monomial, Order, Poly, Ring, Uni};

#[test]
fn fast_public_apis_normalize_and_match_classical_arithmetic() {
    for p in [2, 101, 1_000_000_007, 18_446_744_073_709_551_557] {
        let field = PrimeField::new(p).unwrap();
        let a = Uni::new((0..130).map(|i| Fp::new(i * 17 + 1, p)).collect());
        let b = Uni::new((0..67).map(|i| Fp::new(i * 11 + 1, p)).collect());
        assert_eq!(a.mul_fast(&b), &a * &b);
        assert_eq!(a.divrem_fast(&b), a.divrem(&b));
        assert_eq!(a.gcd_fast(&b), a.gcd(&b));
        if let Some(inverse) = b.inverse_mod(&a) {
            assert_eq!((&inverse * &b).divrem(&a).1, Uni::constant(Fp::new(1, p)));
        } else {
            assert!(a.gcd(&b).deg() > 0);
        }
        assert!(field.divrem(&[1, 0], &[p, 0]).is_none());
        assert!(field.inverse_mod(&[1], &[p, 0]).is_none());
        assert_eq!(field.mul(&[p, 1, 0], &[p, 1, 0]), vec![0, 0, 1]);
        assert_eq!(field.gcd(&[p, 0], &[p, 0]), Vec::<u64>::new());
        assert_eq!(field.scale(&[p, 1, 0], p), Vec::<u64>::new());
    }
    for p in [0, 1, 4, 100] {
        assert!(PrimeField::new(p).is_none());
    }
    let a = Uni::<Gf<101>>::new(vec![Gf::new(1), Gf::new(2)]);
    assert_eq!(a.mul_fast(&a), &a * &a);
}

#[test]
#[should_panic(expected = "incompatible prime fields")]
fn fast_univariate_rejects_mixed_fields() {
    let a = Uni::constant(Fp::new(1, 101));
    let b = Uni::constant(Fp::new(1, 103));
    let _ = a.gcd_fast(&b);
}

#[test]
fn subgroup_logs_validate_orders_membership_and_prime_ranges() {
    assert_eq!(
        SmoothPrimes::below(100, 3).unwrap().collect::<Vec<_>>(),
        vec![89, 73, 41]
    );
    for upper in [0, 1, 2, 3, 4, 5] {
        assert!(SmoothPrimes::below(upper, 3).unwrap().next().is_none());
    }
    assert!(SmoothPrimes::new(0).is_none());
    assert!(SmoothPrimes::new(63).is_none());
    let group = PowerOfTwoSubgroup::from_generator(17, 9, 3).unwrap();
    for e in 0..group.order() {
        assert_eq!(group.log(group.pow(e)), Some(e));
    }
    assert!(group.log(0).is_none());
    assert!(group.log(3).is_none());
    assert!(PowerOfTwoSubgroup::from_generator(17, 1, 3).is_none());
    assert!(PowerOfTwoSubgroup::from_generator(15, 2, 3).is_none());
    assert!(PowerOfTwoSubgroup::from_generator(17, 3, 64).is_none());
    assert!(PowerOfTwoSubgroup::new(17, 5, &mut Rng::new(1)).is_none());
    let trivial = PowerOfTwoSubgroup::new(2, 0, &mut Rng::new(1)).unwrap();
    assert_eq!(trivial.log(1), Some(0));
    assert_eq!(trivial.log(0), None);
}

#[test]
fn cached_images_match_direct_evaluation_including_zero_coordinates() {
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let f = ring
        .parse("3*x^4*y^2+5*x*z^3-y*z+2")
        .unwrap()
        .map(|q| Fp::from_rational(q, 101).unwrap());
    for coordinates in [[0, 0, 0], [0, 2, 3], [1, 0, 4], [2, 3, 5]] {
        let point: Vec<_> = coordinates.into_iter().map(|x| Fp::new(x, 101)).collect();
        let cache = PowerTable::new(&point, &f.degrees());
        assert_eq!(f.eval_cached(&cache), f.eval(&point));
        for (k, image) in f.univariate_images(&point).iter().enumerate() {
            assert_eq!(*image, f.eval_except(k, &point));
            assert_eq!(*image, f.eval_except_cached(k, &cache));
        }
    }
    assert_eq!(
        power_table(&[103, 0], &[3, 2], 101),
        vec![vec![1, 2, 4, 8], vec![1, 0, 0]]
    );
    assert!(
        Poly::<Fp>::zero(0, Order::Lex)
            .univariate_images(&[])
            .is_empty()
    );
}

#[test]
fn geometric_images_match_direct_substitution_for_every_retained_dimension() {
    let order = Order::block(vec![(Order::GRevLex, 2), (Order::Lex, 1)]);
    let ring = Ring::new(["x", "y", "z"], order);
    let f = ring
        .parse("x^3*y+2*x*z+7*y^2*z^3+1")
        .unwrap()
        .map(|q| Fp::from_rational(q, 101).unwrap());
    for keep in 0..=3 {
        let scale: Vec<_> = (keep..3).map(|i| Fp::new(i as u64, 101)).collect();
        let ratios: Vec<_> = (keep..3).map(|i| Fp::new(i as u64 + 2, 101)).collect();
        let mut points = scale.clone();
        let mut evaluator = GeometricEvaluator::new(&f, keep, &scale, &ratios);
        for _ in 0..6 {
            let mut expected = f.clone();
            for (i, (x, r)) in points.iter_mut().zip(&ratios).enumerate() {
                *x = *x * *r;
                expected = expected.eval_var(keep + i, x);
            }
            let terms = expected
                .terms
                .iter()
                .map(|(m, c)| (Monomial::new(m.exps()[..keep].to_vec()), *c))
                .collect();
            assert_eq!(evaluator.advance(), Poly::new(terms, keep, f.order.clone()));
        }
    }
}

#[test]
fn exact_division_preserves_order_and_integer_exactness() {
    for order in [
        Order::Lex,
        Order::GRevLex,
        Order::block(vec![(Order::GrLex, 1), (Order::Lex, 2)]),
    ] {
        let ring = Ring::new(["x", "y", "z"], order);
        let g = ring.parse("3*x^2+2*y*z+1").unwrap();
        let q = ring.parse("x*y-3*z^2+5").unwrap();
        let f = &g * &q;
        assert_eq!(f.exact(&g), Some(q.clone()));
        assert!(q.exact(&g).is_none());
        assert!(f.exact_with_budget(&g, 0).is_none());
        assert_eq!(f.exact_with_budget(&g, 1000), Some(q.clone()));
        let integer = |f: &Poly<num_rational::BigRational>| f.map(|c| c.numer().clone());
        assert_eq!(integer(&f).exact(&integer(&g)), Some(integer(&q)));
        let x = integer(&ring.parse("x").unwrap());
        let two_x = integer(&ring.parse("2*x").unwrap());
        assert!(x.exact(&two_x).is_none());
        let zero = Poly::<BigInt>::zero(3, f.order.clone());
        assert!(zero.exact(&zero).is_none());
        assert_eq!(zero.exact(&two_x), Some(zero));
    }
}

#[test]
fn sparse_division_handles_unpacked_exponents_and_canonicalizes_terms() {
    use polycore::division::exact_lex;
    let e = vec![1; 40];
    let z = vec![0; 40];
    let a = vec![
        (z.clone(), BigInt::from(1)),
        (e.clone(), BigInt::from(2)),
        (z.clone(), BigInt::from(1)),
    ];
    let b = vec![(e.clone(), BigInt::from(1)), (z.clone(), BigInt::from(1))];
    assert_eq!(
        exact_lex(40, &a, &b, None),
        Some(vec![(z, BigInt::from(2))])
    );
    assert!(exact_lex(39, &a, &b, None).is_none());
    assert!(exact_lex(40, &a, &b, Some(0)).is_none());
    let bad = vec![(e, BigInt::from(3))];
    assert!(exact_lex(40, &a, &bad, None).is_none());
}
