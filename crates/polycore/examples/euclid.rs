use num_rational::BigRational;
use polycore::{Order, Ring, Uni};

fn main() -> Result<(), polycore::ParseError> {
    let ring = Ring::new(["x"], Order::Lex);
    let uni = |src| Ok::<_, polycore::ParseError>(Uni::from_poly(&ring.parse(src)?, 0).unwrap());
    let show = |u: &Uni<BigRational>| ring.show(&u.to_poly(0, 1, Order::Lex));
    let (f, g) = (uni("x^4 - 1")?, uni("x^6 - 1")?);
    let (s, t, gcd) = f.bezout(&g);
    println!("gcd(x^4 - 1, x^6 - 1) = {}", show(&gcd));
    println!(
        "({}) * (x^4 - 1) + ({}) * (x^6 - 1) = {}",
        show(&s),
        show(&t),
        show(&gcd)
    );
    Ok(())
}
