use num_rational::BigRational;
use num_traits::{One, Zero};
use polycore::{Order, RatFunc, Ring, Uni};

fn main() {
    let x = RatFunc::var();
    let k = |n: i64| RatFunc::from(Uni::constant(BigRational::from_integer(n.into())));
    let sum = (0..5).fold(RatFunc::zero(), |acc, n| {
        acc + RatFunc::one() / ((x.clone() + k(n)) * (x.clone() + k(n + 1)))
    });
    let ring = Ring::new(["x"], Order::Lex);
    let show = |u: &Uni<BigRational>| ring.show(&u.to_poly(0, 1, Order::Lex));
    println!(
        "sum 1/((x+k)(x+k+1)), k = 0..4 = ({}) / ({})",
        show(sum.num()),
        show(sum.den())
    );
}
