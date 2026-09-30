use num_rational::BigRational;
use polycore::interp::Newton;
use polycore::{Order, Ring};

fn main() {
    let q = |n: i64| BigRational::from_integer(n.into());
    let mut newton = Newton::default();
    let mut sum = 0;
    for n in 1.. {
        sum += n * n;
        if !newton.add(q(n), q(sum)) {
            break;
        }
    }
    let ring = Ring::new(["n"], Order::Lex);
    println!(
        "1^2 + ... + n^2 = {}",
        ring.show(&newton.poly().to_poly(0, 1, Order::Lex))
    );
}
