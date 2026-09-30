use num_rational::BigRational;
use polycore::interp::Thiele;
use polycore::{Order, Ring};

fn main() {
    let q = |n: i64| BigRational::from_integer(n.into());
    let f = |t: &BigRational| q(1) / (q(1) - t - t * t);
    let mut thiele = Thiele::default();
    for t in 0.. {
        if !thiele.add(q(t), f(&q(t))) {
            break;
        }
    }
    let ring = Ring::new(["t"], Order::Lex);
    let (num, den) = thiele.rational();
    let show = |u: &polycore::Uni<BigRational>| ring.show(&u.to_poly(0, 1, Order::Lex));
    println!("f(t) = ({}) / ({})", show(&num), show(&den));
}
