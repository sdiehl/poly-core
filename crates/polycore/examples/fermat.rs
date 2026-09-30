use polycore::{Gf, pow};

type F = Gf<101>;

fn main() {
    for a in [2, 3, 10, 57] {
        let x = F::new(a);
        println!("{a}^100 = {} mod 101", pow(&x, 100));
    }
}
