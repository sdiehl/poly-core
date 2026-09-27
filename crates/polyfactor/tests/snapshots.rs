use std::fmt::Write;
use std::rc::Rc;

use num_rational::BigRational;
use num_traits::Zero;
use polycore::{Monomial, Order, Poly, Ring, Uni};
use polyfactor::{Alg, NumberField, factor, factor_mod, factor_over, is_irreducible};

type Q = BigRational;

struct Session {
    ring: Ring,
    p: u64,
    k: Option<Rc<NumberField>>,
    out: String,
}

impl Session {
    fn parse(&self, src: &str) -> Poly<Q> {
        self.ring
            .parse(src)
            .unwrap_or_else(|e| panic!("{src}: {e}"))
    }

    fn uni(&self, src: &str, k: usize) -> Uni<Q> {
        Uni::from_poly(&self.parse(src), k).expect("univariate")
    }

    fn alg(&self, src: &str, k: &Rc<NumberField>) -> Uni<Alg> {
        let mut cs = Vec::new();
        for (m, c) in self.parse(src).coeffs_in(1) {
            let i = m.exps()[0] as usize;
            cs.resize(cs.len().max(i + 1), Alg::zero());
            cs[i] = Alg::new(c, k);
        }
        Uni::new(cs)
    }

    fn show(&self, f: &Uni<Alg>) -> String {
        let groups = f.0.iter().enumerate().map(|(i, c)| {
            let i = u32::try_from(i).unwrap();
            (Monomial::new(vec![i, 0]), c.c.clone())
        });
        self.ring
            .show(&Poly::from_coeffs_in(1, 2, Order::Lex, groups))
    }

    fn line(&mut self, s: impl std::fmt::Display) {
        writeln!(self.out, "{s}").unwrap();
    }

    fn run(&mut self, cmd: &str) {
        let (head, rest) = cmd.split_once(' ').unwrap_or((cmd, ""));
        match head {
            "mod" => (self.p, self.k) = (rest.parse().unwrap(), None),
            "field" => {
                let (m, k) = rest
                    .split_once(';')
                    .map_or((rest, None), |(m, k)| (m, Some(k)));
                let mut field = NumberField::new(&self.uni(m, 1));
                if let Some(k) = k {
                    let k: usize = k.trim().parse().unwrap();
                    field.interval = Some(field.m.isolate().swap_remove(k - 1));
                }
                self.k = Some(Rc::new(field));
            }
            "irreducible" => {
                let f = self.uni(rest, 0);
                self.line(format!("{}", is_irreducible(&f)));
            }
            "alg" => {
                let k = self.k.clone().expect("a field");
                let a = self.alg(rest, &k).eval(&Alg::zero());
                let inv = Uni::constant(Alg::from(Q::from_integer(1.into())) / a.clone());
                self.line(format!(
                    "norm {}, trace {}, positive {}, inverse {}",
                    a.norm(&k),
                    a.trace(&k),
                    a.positive(),
                    self.show(&inv)
                ));
            }
            "factor" => {
                let parts: Vec<String> = if let Some(k) = self.k.clone() {
                    let parts = factor_over(&self.alg(rest, &k), &k);
                    parts
                        .iter()
                        .map(|(h, m)| format!("[{m}] {}", self.show(h)))
                        .collect()
                } else if self.p > 0 {
                    let f = self.uni(rest, 0).modp(self.p).expect("no inverse");
                    let (lc, parts) = factor_mod(&f);
                    let mut out = vec![format!("lc {}", lc.value())];
                    out.extend(
                        parts
                            .iter()
                            .map(|(h, m)| format!("[{m}] {}", self.show_q(&h.residues()))),
                    );
                    out
                } else {
                    let parts = factor(&self.uni(rest, 0));
                    parts
                        .iter()
                        .map(|(h, m)| format!("[{m}] {}", self.show_q(h)))
                        .collect()
                };
                for p in parts {
                    self.line(p);
                }
            }
            _ => panic!("unknown command {head}"),
        }
    }

    fn show_q(&self, f: &Uni<Q>) -> String {
        self.ring.show(&f.to_poly(0, 2, Order::Lex))
    }
}

#[test]
fn cases() {
    insta::glob!("cases/*.txt", |path| {
        let src = std::fs::read_to_string(path).unwrap();
        let mut s = Session {
            ring: Ring::new(["x", "a"], Order::Lex),
            p: 0,
            k: None,
            out: String::new(),
        };
        for cmd in src
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
        {
            s.line(format!("> {cmd}"));
            s.run(cmd);
        }
        insta::assert_snapshot!(s.out);
    });
}
