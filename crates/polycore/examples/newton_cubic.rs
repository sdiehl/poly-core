use num_rational::BigRational;
use num_traits::ToPrimitive;
use polycore::Uni;

fn main() {
    let q = |n: i64| BigRational::from_integer(n.into());
    let f = Uni::new(vec![q(-5), q(-2), q(0), q(1)]);
    let (mut lo, mut hi) = f.isolate()[0].clone();
    for _ in 0..40 {
        (lo, hi) = f.refine(&lo, &hi);
    }
    println!("real roots of x^3 - 2x - 5: {}", f.isolate().len());
    println!(
        "root in ({:.12}, {:.12}]",
        lo.to_f64().unwrap_or(0.0),
        hi.to_f64().unwrap_or(0.0)
    );
}
