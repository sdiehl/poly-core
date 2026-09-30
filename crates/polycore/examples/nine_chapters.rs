use num_rational::BigRational;
use polycore::dense;

fn main() {
    let q = |n: i64| BigRational::from_integer(n.into());
    let a = [
        vec![q(3), q(2), q(1)],
        vec![q(2), q(3), q(1)],
        vec![q(1), q(2), q(3)],
    ];
    let b = [q(39), q(34), q(26)];
    let x = dense::solve(&a, &b).expect("unique solution");
    println!("top = {}, middle = {}, bottom = {}", x[0], x[1], x[2]);
}
