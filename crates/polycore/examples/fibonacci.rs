use num_rational::BigRational;
use polycore::interp::Massey;
use polycore::{Order, Ring};

fn main() {
    let mut massey = Massey::default();
    let (mut a, mut b) = (0i64, 1i64);
    for _ in 0..10 {
        massey.push(BigRational::from_integer(a.into()));
        (a, b) = (b, a + b);
    }
    let ring = Ring::new(["z"], Order::Lex);
    println!(
        "recurrence: {}",
        ring.show(&massey.generator().to_poly(0, 1, Order::Lex))
    );
}
