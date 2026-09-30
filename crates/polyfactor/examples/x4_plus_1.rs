use polycore::{Fp, Order, Ring, Uni};
use polyfactor::factor_mod;

fn main() {
    let ring = Ring::new(["x"], Order::Lex);
    for p in [2, 3, 5, 7, 11, 13, 17] {
        let f = Uni::new(vec![
            Fp::new(1, p),
            Fp::new(0, p),
            Fp::new(0, p),
            Fp::new(0, p),
            Fp::new(1, p),
        ]);
        let (_, parts) = factor_mod(&f);
        let parts: Vec<String> = parts
            .iter()
            .map(|(g, e)| {
                let g = ring.show(&g.to_poly(0, 1, Order::Lex));
                if *e == 1 {
                    format!("({g})")
                } else {
                    format!("({g})^{e}")
                }
            })
            .collect();
        println!("x^4 + 1 mod {p} = {}", parts.join(" * "));
    }
}
