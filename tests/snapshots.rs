#![allow(clippy::many_single_char_names)]

use std::fmt::{Display, Write};

use num_rational::BigRational;
use poly_core::{combination, crt, Echelon, Field, Fp, Lead, Order, Poly, Primes, Ring, Uni};

type Q = BigRational;

struct Session {
    ring: Ring,
    p: u64,
    out: String,
}

impl Session {
    fn line(&mut self, s: impl Display) {
        writeln!(self.out, "{s}").unwrap();
    }

    fn parse(&self, src: &str) -> Poly<Q> {
        self.ring
            .parse(src)
            .unwrap_or_else(|e| panic!("{src}: {e}"))
    }

    fn rational(&self, src: &str) -> Q {
        let p = self.parse(src);
        assert!(p.is_constant(), "{src} is not a constant");
        p.lc().cloned().unwrap_or_default()
    }

    fn run(&mut self, cmd: &str) {
        let (head, rest) = cmd.split_once(' ').unwrap_or((cmd, ""));
        let args: Vec<&str> = rest.split(';').map(str::trim).collect();
        match head {
            "ring" => {
                let mut words = rest.split_whitespace();
                let order = words.next().unwrap().parse().unwrap();
                self.ring = Ring::new(words, order);
            }
            "order" => {
                let n = self.ring.nvars();
                let order = match rest.strip_prefix("elim ") {
                    Some(k) => {
                        Order::elimination(k.parse().unwrap(), n - k.parse::<usize>().unwrap())
                    }
                    None if rest == "weighted" => Order::weighted(vec![1; n], Order::GRevLex),
                    None => rest.parse().unwrap(),
                };
                self.ring = self.ring.with_order(order);
            }
            "mod" => self.p = rest.parse().unwrap(),
            "crt" => self.crt(&args),
            _ if self.p == 0 => self.field(head, &args, Q::clone),
            _ => {
                let p = self.p;
                self.field(head, &args, |q| {
                    Fp::from_rational(q, p).expect("p divides a denominator")
                });
            }
        }
    }

    fn field<F: Field + Display>(&mut self, head: &str, args: &[&str], conv: impl Fn(&Q) -> F) {
        let poly = |s: &Self, src: &str| s.parse(src).map(&conv);
        let show = |s: &Self, p: &Poly<F>| s.ring.show(p);
        match head {
            "show" => {
                let p = poly(self, args[0]);
                self.line(show(self, &p));
            }
            "eval" => {
                let f = poly(self, args[0]);
                let x: Vec<F> = args[1..].iter().map(|a| conv(&self.rational(a))).collect();
                self.line(f.eval(&x));
            }
            "subst" => {
                let f = poly(self, args[0]);
                let k = self.ring.names.iter().position(|v| v == args[1]).unwrap();
                let a = conv(&self.rational(args[2]));
                self.line(show(self, &f.eval_var(k, &a)));
            }
            "spoly" => {
                let (f, g) = (poly(self, args[0]), poly(self, args[1]));
                self.line(show(self, &f.spoly(&g)));
            }
            "divide" => {
                let f = poly(self, args[0]);
                let gs: Vec<Poly<F>> = args[1..].iter().map(|a| poly(self, a)).collect();
                let (q, r) = f.divide(&gs);
                assert_eq!(&combination(&q, &gs) + &r, f, "certificate");
                for (i, qi) in q.iter().enumerate() {
                    self.line(format!("q{} = {}", i + 1, show(self, qi)));
                }
                self.line(format!("r = {}", show(self, &r)));
            }
            "uni" => {
                let to_uni =
                    |s: &Self, src: &str| Uni::from_poly(&poly(s, src), 0).expect("univariate");
                let (f, g) = (to_uni(self, args[0]), to_uni(self, args[1]));
                let (s, t, h) = f.bezout(&g);
                assert_eq!(&(&s * &f) + &(&t * &g), h, "bezout");
                assert_eq!(h, f.gcd(&g));
                let n = self.ring.nvars();
                let back = |u: &Uni<F>| show(self, &u.to_poly(0, n, self.ring.order.clone()));
                let lines = [
                    format!("gcd = {}", back(&h)),
                    format!("s = {}", back(&s)),
                    format!("t = {}", back(&t)),
                    format!("res = {}", f.resultant(&g)),
                    format!("f' = {}", back(&f.derivative())),
                    format!("f / g = {}, rem {}", back(&(&f / &g)), back(&(&f % &g))),
                ];
                for l in lines {
                    self.line(l);
                }
            }
            "interp" => {
                let nums = |s: &Self, a: &str| -> Vec<F> {
                    a.split_whitespace().map(|x| conv(&s.rational(x))).collect()
                };
                let (xs, ys) = (nums(self, args[0]), nums(self, args[1]));
                let u = Uni::interpolate(&xs, &ys);
                assert!(xs.iter().zip(&ys).all(|(x, y)| u.eval(x) == *y));
                let n = self.ring.nvars();
                self.line(show(self, &u.to_poly(0, n, self.ring.order.clone())));
            }
            "echelon" => self.echelon(args, &conv),
            _ => panic!("unknown command {head}"),
        }
    }

    fn echelon<F: Field + Display>(&mut self, args: &[&str], conv: impl Fn(&Q) -> F) {
        let lead = if args[0] == "high" {
            Lead::High
        } else {
            Lead::Low
        };
        let rows: Vec<Vec<F>> = args[1..]
            .iter()
            .map(|r| {
                r.split_whitespace()
                    .map(|x| conv(&self.rational(x)))
                    .collect()
            })
            .collect();
        let ncols = rows[0].len();
        let sparse = |r: &[F]| {
            r.iter()
                .cloned()
                .enumerate()
                .filter(|t| !t.1.is_zero())
                .collect::<Vec<_>>()
        };
        let dense = |r: &[(usize, F)]| {
            let mut v = vec![F::zero(); ncols];
            for (j, a) in r {
                v[*j] = a.clone();
            }
            v.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(" ")
        };
        let mut ech = Echelon::recording(lead);
        let inserted: Vec<Option<usize>> = rows.iter().map(|r| ech.insert(&sparse(r))).collect();
        self.line(format!(
            "rank {} pivots {:?} inserted {inserted:?}",
            ech.rank(),
            ech.pivots()
        ));
        for c in ech
            .pivots()
            .into_iter()
            .filter(|&c| !ech.uses(c).is_empty())
        {
            self.line(format!("row {c} used {:?}", ech.uses(c)));
        }
        for r in ech.rref() {
            self.line(format!("rref [{}]", dense(&r)));
        }
        for v in ech.nullspace(ncols) {
            for r in &rows {
                let dot = v
                    .iter()
                    .fold(F::zero(), |acc, (j, a)| acc + r[*j].clone() * a.clone());
                assert!(dot.is_zero(), "null vector");
            }
            self.line(format!("null [{}]", dense(&v)));
        }
    }

    fn crt(&mut self, args: &[&str]) {
        let qs: Vec<Q> = args[1..].iter().map(|a| self.rational(a)).collect();
        let mut used = Vec::new();
        let image = |p: u64| {
            used.push(p);
            qs.iter().map(|q| crt::reduce(q, p)).collect()
        };
        let got = crt::reconstruct(image, Primes::below(args[0].parse().unwrap()));
        assert_eq!(got.as_ref(), Some(&qs));
        self.line(format!("{} primes {:?}", used.len(), used));
        let residues: Vec<(u64, u64)> = used
            .iter()
            .map(|&p| (crt::reduce(&qs[0], p).unwrap(), p))
            .collect();
        let (x, m) = crt::crt(&residues);
        self.line(format!("{} = {x} mod {m}", qs[0]));
        let back: Vec<String> = got.unwrap().iter().map(ToString::to_string).collect();
        self.line(back.join(", "));
    }
}

#[test]
fn cases() {
    insta::glob!("cases/*.txt", |path| {
        let src = std::fs::read_to_string(path).unwrap();
        let mut s = Session {
            ring: Ring::new(["x", "y", "z"], Order::GRevLex),
            p: 0,
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
