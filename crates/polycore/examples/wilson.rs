use num_traits::One;
use polycore::Fp;

fn main() {
    for p in [5, 7, 11, 13, 101] {
        let fact = (1..p).fold(Fp::one(), |acc, k| acc * Fp::new(k, p));
        println!("({p} - 1)! = {fact} mod {p}");
    }
}
