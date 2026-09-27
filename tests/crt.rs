use num_bigint::BigInt;
use num_rational::BigRational;
use polycore::crt::{self, CrtAccumulator, CrtError, WangContext};

#[test]
fn checked_updates_preserve_state_on_error() {
    let mut acc = CrtAccumulator::new(2);
    acc.add(5, &[1, 2]).unwrap();
    for (p, values, expected) in [
        (
            7,
            vec![1],
            CrtError::LengthMismatch {
                expected: 2,
                actual: 1,
            },
        ),
        (
            7,
            vec![1, 2, 3],
            CrtError::LengthMismatch {
                expected: 2,
                actual: 3,
            },
        ),
        (5, vec![1, 2], CrtError::NonCoprimeModuli),
        (1, vec![0, 0], CrtError::InvalidModulus),
        (7, vec![1, 7], CrtError::UnreducedResidue),
    ] {
        assert_eq!(acc.add(p, &values), Err(expected));
        assert_eq!(acc.residues(), &[BigInt::from(1), BigInt::from(2)]);
        assert_eq!(acc.modulus(), &BigInt::from(5));
        assert_eq!(acc.image_count(), 1);
    }
    let mut xs = vec![BigInt::from(1)];
    let mut m = BigInt::from(0);
    assert_eq!(
        crt::try_garner(&mut xs, &mut m, &[1], 7),
        Err(CrtError::InvalidModulus)
    );
    assert_eq!(xs, vec![BigInt::from(1)]);
    assert_eq!(m, BigInt::from(0));
}

#[test]
fn persistent_reconstruction_and_independent_verification() {
    let truth: Vec<_> = [(0, 1), (-7, 5), (123_456, 789), (1, 3)]
        .into_iter()
        .map(|(n, d)| BigRational::new(n.into(), d.into()))
        .collect();
    let mut acc = CrtAccumulator::new(truth.len());
    assert_eq!(acc.reconstruct(), None);
    for (i, p) in [101, 103, 107, 109, 113].into_iter().enumerate() {
        let values: Vec<_> = truth.iter().map(|q| crt::reduce(q, p).unwrap()).collect();
        acc.add(p, &values).unwrap();
        assert_eq!(acc.image_count(), i + 1);
        for (&v, x) in values.iter().zip(acc.residues()) {
            assert_eq!(crt::residue(x, p), v);
        }
        let _candidate = acc.reconstruct(); // An attempt must not consume accumulated state.
    }
    let got = acc.reconstruct().unwrap();
    assert_eq!(got, truth);
    for (a, b) in got.iter().zip(&truth) {
        assert_eq!(crt::reduce(a, 127), crt::reduce(b, 127));
    }
    acc.add(
        127,
        &truth
            .iter()
            .map(|q| crt::reduce(q, 127).unwrap())
            .collect::<Vec<_>>(),
    )
    .unwrap();
    assert_eq!(acc.reconstruct(), Some(truth));
}

#[test]
fn shared_wang_bounds_handle_signed_and_noncanonical_residues() {
    let m = BigInt::from(1009);
    let context = WangContext::new(&m).unwrap();
    for n in -12..=12 {
        for d in 1..=12 {
            let q = BigRational::new(n.into(), d.into());
            let x = BigInt::from(crt::reduce(&q, 1009).unwrap());
            for shifted in [&x - &m, x.clone(), &x + &m] {
                assert_eq!(context.reconstruct(&shifted), Some(q.clone()));
                assert_eq!(crt::wang(&shifted, &m), Some(q.clone()));
            }
        }
    }
    assert_eq!(context.reconstruct_many(&[]), Some(vec![]));
    assert_eq!(
        WangContext::new(&BigInt::from(5))
            .unwrap()
            .reconstruct_many(&[0.into(), 2.into()]),
        None
    );
    for m in [-1, 0, 1] {
        assert!(WangContext::new(&m.into()).is_none());
    }
}

#[test]
fn zero_width_accumulator_tracks_images() {
    let mut acc = CrtAccumulator::new(0);
    assert_eq!(acc.reconstruct(), None);
    acc.add(5, &[]).unwrap();
    assert_eq!(acc.reconstruct(), Some(vec![]));
    assert_eq!(acc.image_count(), 1);
}

#[test]
fn legacy_garner_rejects_lengths_before_mutation() {
    let mut xs = vec![BigInt::from(1), BigInt::from(2)];
    let mut m = BigInt::from(5);
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        crt::garner(&mut xs, &mut m, &[1], 7);
    }));
    assert!(result.is_err());
    assert_eq!(xs, vec![BigInt::from(1), BigInt::from(2)]);
    assert_eq!(m, BigInt::from(5));
}

#[test]
fn reconstruction_helpers_reject_inconsistent_lengths() {
    for voted in [false, true] {
        assert!(
            std::panic::catch_unwind(|| {
                let image = |p| Some(((), vec![0; if p == 5 { 2 } else { 1 }]));
                if voted {
                    crt::reconstruct_voted(|ps| ps.iter().map(|&p| image(p)).collect(), [5, 7]);
                } else {
                    crt::reconstruct_keyed(image, [5, 7]);
                }
            })
            .is_err()
        );
    }
    for len in [0, 2] {
        assert!(
            std::panic::catch_unwind(|| {
                crt::reconstruct_voted(|_| vec![Some(((), vec![0])); len], [5]);
            })
            .is_err()
        );
    }
}
