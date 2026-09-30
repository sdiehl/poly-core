use polycore::{Fp, Primes, crt};

fn main() {
    let image = |p| {
        let h = (1..=20).fold(Fp::new(0, p), |s, k| s + Fp::new(1, p) / Fp::new(k, p));
        Some(vec![h.value()])
    };
    let h = crt::reconstruct(image, Primes::below(100)).expect("enough primes");
    println!("H_20 = {}", h[0]);
}
