use num_rational::BigRational;
use polycore::{Echelon, Fp, Gf, Lead, Order, Primes, Ring, Uni, combination, crt};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ring = Ring::new(["x", "y"], Order::Lex);
    let f = ring.parse("x^2*y + x*y^2 + y^2")?;
    let gens = [ring.parse("x*y - 1")?, ring.parse("y^2 - 1")?];
    let (quots, rem) = f.divide(&gens);
    assert_eq!(&combination(&quots, &gens) + &rem, f);
    println!("f = {}", ring.show(&f));
    println!("  = ({}) * (x*y - 1)", ring.show(&quots[0]));
    println!("  + ({}) * (y^2 - 1)", ring.show(&quots[1]));
    println!("  + {}", ring.show(&rem));

    let prime = Primes::new().next().unwrap_or(2);
    let fp = f.map(|c| Fp::from_rational(c, prime).expect("integral"));
    let at = fp.eval(&[Fp::new(2, prime), Fp::new(3, prime)]);
    println!("f(2, 3) mod {prime} = {at}");

    let line = Ring::new(["x"], Order::Lex);
    let show = |u: &Uni<Gf<101>>| line.show(&u.to_poly(0, 1, Order::Lex));
    let lhs = Uni::new(vec![Gf::from_i64(-1), Gf::new(0), Gf::new(1)]);
    let rhs = Uni::new(vec![Gf::new(3), Gf::new(4), Gf::new(1)]);
    let (s, t, gcd) = lhs.bezout(&rhs);
    println!(
        "over GF(101): ({}) * (x^2 - 1) + ({}) * (x^2 + 4*x + 3) = {}",
        show(&s),
        show(&t),
        show(&gcd)
    );

    let secret = [
        BigRational::new(355.into(), 113.into()),
        BigRational::new((-7).into(), 3.into()),
    ];
    let image = |p| secret.iter().map(|q| crt::reduce(q, p)).collect();
    let got = crt::reconstruct(image, Primes::below(1 << 12)).unwrap_or_default();
    let got: Vec<String> = got.iter().map(ToString::to_string).collect();
    println!("reconstructed from primes below 4096: {}", got.join(", "));

    let mut ech = Echelon::new(Lead::Low);
    for row in [[1, 2, 3], [2, 4, 6], [1, 0, 1]] {
        let row: Vec<(usize, BigRational)> = row
            .iter()
            .map(|&v| BigRational::from_integer(v.into()))
            .enumerate()
            .collect();
        ech.insert(&row);
    }
    let null: Vec<String> = ech.nullspace(3)[0]
        .iter()
        .map(|(_, v)| v.to_string())
        .collect();
    println!("rank {}, null vector ({})", ech.rank(), null.join(", "));
    Ok(())
}
