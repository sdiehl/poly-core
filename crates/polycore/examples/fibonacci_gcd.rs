use num_bigint::BigUint;
use polycore::lehmer;

fn fib(n: usize) -> BigUint {
    let (mut a, mut b) = (BigUint::from(0u32), BigUint::from(1u32));
    for _ in 0..n {
        (a, b) = (b.clone(), a + b);
    }
    a
}

fn main() {
    let g = lehmer::gcd(&fib(300), &fib(200));
    println!("gcd(F_300, F_200) = {g}");
    println!("F_100             = {}", fib(100));
}
