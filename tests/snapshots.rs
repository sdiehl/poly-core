#![allow(clippy::many_single_char_names)]

use std::fmt::{Display, Write};

use num_rational::BigRational;
use poly_core::interp::{self, Massey, Thiele};
use poly_core::{
    combination, crt, dense, modp, Echelon, Field, Fp, Lead, Order, Poly, Primes, RatFunc, Ring,
    Uni,
};

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

    fn nums<F>(&self, src: &str, conv: impl Fn(&Q) -> F) -> Vec<F> {
        src.split_whitespace()
            .map(|x| conv(&self.rational(x)))
            .collect()
    }

    fn univariate<F: Field + Display>(&self, u: &Uni<F>) -> String {
        let n = self.ring.nvars();
        self.ring.show(&u.to_poly(0, n, self.ring.order.clone()))
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
            "badring" => match Ring::try_new(rest.split(' '), Order::Lex) {
                Ok(r) => self.line(format!("ok {:?}", r.names)),
                Err(e) => self.line(e),
            },
            "mod" => self.p = rest.parse().unwrap(),
            "crt" => self.crt(&args),
            "prime" => {
                let ns: Vec<u64> = rest
                    .split_whitespace()
                    .map(|n| n.parse().unwrap())
                    .collect();
                let primes: Vec<_> = ns.iter().filter(|&&n| modp::is_prime(n)).collect();
                self.line(format!("primes {primes:?}"));
            }
            "many" => {
                for p in self.ring.parse_many(rest).unwrap() {
                    let s = self.ring.show(&p);
                    self.line(s);
                }
            }
            "latex" => {
                let s = self.ring.latex(&self.parse(rest));
                self.line(s);
            }
            "content" => {
                let f = self.parse(rest);
                let s = self.ring.show(&f.primitive());
                self.line(format!("{} * ({s})", f.content()));
            }
            "param" => self.param(&args),
            "roots" => {
                let p = self.p;
                let f = self.parse(rest).map(|q| Fp::from_rational(q, p).unwrap());
                let u = Uni::from_poly(&f, 0).expect("univariate");
                let roots = u.roots();
                assert!(roots.iter().all(|r| u.eval(r) == Fp::new(0, p)));
                let shown: Vec<String> = roots.iter().map(ToString::to_string).collect();
                self.line(format!("roots [{}]", shown.join(", ")));
            }
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
            _ => self.extra(head, args, &conv),
        }
    }

    fn extra<F: Field + Display>(&mut self, head: &str, args: &[&str], conv: impl Fn(&Q) -> F) {
        let poly = |s: &Self, src: &str| s.parse(src).map(&conv);
        let show = |s: &Self, p: &Poly<F>| s.ring.show(p);
        match head {
            "calc" => {
                let f = poly(self, args[0]);
                for (k, v) in self.ring.names.clone().iter().enumerate() {
                    self.line(format!("d/d{v} = {}", show(self, &f.derivative(k))));
                }
                self.line(format!("^2 = {}", show(self, &f.pow(2))));
                let perm: Vec<usize> = (0..self.ring.nvars()).rev().collect();
                self.line(format!("reversed = {}", show(self, &f.permute(&perm))));
                let (c, pp) = f.primitive_in(0);
                let c = self.univariate(&c);
                self.line(format!(
                    "content in x0 = {c}, primitive = {}",
                    show(self, &pp)
                ));
                for (m, u) in f.coeffs_in(0) {
                    self.line(format!("[{m}] {}", self.univariate(&u)));
                }
            }
            "unix" => {
                let to_uni =
                    |s: &Self, src: &str| Uni::from_poly(&poly(s, src), 0).expect("univariate");
                let (f, g) = (to_uni(self, args[0]), to_uni(self, args[1]));
                let fg = &f * &g;
                assert_eq!(fg.exact(&g), Some(f.clone()));
                let lines = [
                    format!("f^3 = {}", self.univariate(&f.pow(3))),
                    format!("f(g) = {}", self.univariate(&f.compose(&g))),
                    format!("f^5 mod g = {}", self.univariate(&f.powmod(5, &g))),
                    format!("int f = {}", self.univariate(&f.integral())),
                    format!(
                        "f / (f + 1) exact = {:?}",
                        f.exact(&(&f + &Uni::x())).is_some()
                    ),
                ];
                for l in lines {
                    self.line(l);
                }
            }
            "thiele" => {
                let (xs, ys) = (self.nums(args[0], &conv), self.nums(args[1], &conv));
                let mut th = Thiele::default();
                let used = xs
                    .into_iter()
                    .zip(ys)
                    .take_while(|(x, y)| th.add(x.clone(), y.clone()));
                let n = used.count();
                let (num, den) = th.rational();
                let (num, den) = (self.univariate(&num), self.univariate(&den));
                self.line(format!("{n} points: ({num}) / ({den})"));
            }
            "massey" => {
                let mut m = Massey::default();
                for a in self.nums(args[0], &conv) {
                    m.push(a);
                }
                let g = m.generator();
                self.line(format!(
                    "L = {} settled {}: {}",
                    m.complexity(),
                    m.settled(2),
                    self.univariate(&g)
                ));
            }
            "vander" => {
                let (vals, cs) = (self.nums(args[0], &conv), self.nums(args[1], &conv));
                let w: Vec<F> = (1..=vals.len() as u64)
                    .map(|j| {
                        vals.iter().zip(&cs).fold(F::zero(), |acc, (v, c)| {
                            acc + c.clone() * poly_core::pow(v, j)
                        })
                    })
                    .collect();
                let got = interp::solve(&vals, &interp::master(&vals), &w);
                assert_eq!(got, cs, "vandermonde");
                self.line(format!(
                    "master = {}",
                    self.univariate(&interp::master(&vals))
                ));
            }
            "dense" => self.dense(args, &conv),
            _ => panic!("unknown command {head}"),
        }
    }

    fn dense<F: Field + Display>(&mut self, args: &[&str], conv: impl Fn(&Q) -> F) {
        let rows: Vec<Vec<F>> = args.iter().map(|r| self.nums(r, &conv)).collect();
        let (a, b): (Vec<Vec<F>>, Vec<F>) = rows
            .iter()
            .map(|r| (r[..r.len() - 1].to_vec(), r[r.len() - 1].clone()))
            .unzip();
        let fmt = |r: &[F]| {
            r.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(" ")
        };
        let mut m = a.clone();
        let pivots = dense::rref(&mut m);
        self.line(format!("pivots {pivots:?}"));
        for v in dense::nullspace(&a) {
            self.line(format!("null [{}]", fmt(&v)));
        }
        match dense::solve(&a, &b) {
            Ok(x) => self.line(format!("x = [{}]", fmt(&x))),
            Err(e) => self.line(e),
        }
        if let Some(inv) = dense::invert(&a) {
            for r in &inv {
                self.line(format!("inv [{}]", fmt(r)));
            }
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

    fn param(&mut self, args: &[&str]) {
        let a = RatFunc::<Q>::var();
        let lift = |n| RatFunc::from(Uni::constant(Q::from_integer(n)));
        let f = self
            .ring
            .parse_with(args[1], &lift, &[(args[0], a)])
            .unwrap_or_else(|e| panic!("{e}"));
        let at = self.rational(args[2]);
        let special = f.map(|c| c.eval(&at).expect("pole"));
        let t = Ring::new([args[0]], Order::Lex);
        let side = |u: &Uni<Q>| t.show(&u.to_poly(0, 1, Order::Lex));
        for (m, c) in &f.terms {
            let m = self.ring.monomial(m);
            self.line(format!("[{m}] ({}) / ({})", side(c.num()), side(c.den())));
        }
        let s = self.ring.show(&special);
        self.line(format!("at {} = {at}: {s}", args[0]));
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
        let unlucky = used[1];
        let keyed = crt::reconstruct_keyed(
            |p| {
                let key = u8::from(p != unlucky);
                Some((
                    key,
                    qs.iter()
                        .map(|q| crt::reduce(q, p).map_or(0, |v| v ^ u64::from(key == 0)))
                        .collect(),
                ))
            },
            Primes::below(args[0].parse().unwrap()),
        );
        assert_eq!(keyed.map(|k| k.1).as_ref(), got.as_ref(), "keyed");
        self.line(format!("symmetric {}", crt::symmetric(&x, &m)));
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
