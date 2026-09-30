use polycore::{Fp, Primes, Uni};

fn main() {
    let mut primes: Vec<u64> = Primes::below(60).collect();
    primes.reverse();
    for p in primes {
        let f = Uni::new(vec![Fp::new(1, p), Fp::new(0, p), Fp::new(1, p)]);
        let roots: Vec<String> = f.roots().iter().map(ToString::to_string).collect();
        println!("x^2 + 1 = 0 mod {p}: {{{}}}", roots.join(", "));
    }
}
