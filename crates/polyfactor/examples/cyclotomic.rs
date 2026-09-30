use num_rational::BigRational;
use polycore::{Order, Ring, Uni};
use polyfactor::factor;

fn main() {
    let q = |n: i64| BigRational::from_integer(n.into());
    let ring = Ring::new(["x"], Order::Lex);
    let mut c = vec![q(0); 13];
    (c[0], c[12]) = (q(-1), q(1));
    for (f, _) in factor(&Uni::new(c)) {
        println!("{}", ring.show(&f.to_poly(0, 1, Order::Lex)));
    }
}
