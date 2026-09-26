use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::{self, Display};
use std::ops::{Add, Mul, Neg, Sub};

use crate::field::{nat, pow, Field};
use crate::monomial::{Monomial, Order};
use crate::parse::Ring;
use crate::uni::Uni;

pub type Term<F> = (Monomial, F);

/// A sparse multivariate polynomial: nonzero terms in strictly descending `order`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly<F> {
    pub terms: Vec<Term<F>>,
    pub nvars: usize,
    pub order: Order,
}

impl<F: Field> Poly<F> {
    /// Sorts, combines like terms and drops zeros.
    pub fn new(mut terms: Vec<Term<F>>, nvars: usize, order: Order) -> Self {
        terms.sort_by(|a, b| order.compare(&b.0, &a.0));
        let mut out: Vec<Term<F>> = Vec::with_capacity(terms.len());
        for (m, c) in terms {
            match out.last_mut() {
                Some(last) if last.0 == m => last.1 = last.1.clone() + c,
                _ => out.push((m, c)),
            }
        }
        out.retain(|t| !t.1.is_zero());
        Self {
            terms: out,
            nvars,
            order,
        }
    }

    pub const fn zero(nvars: usize, order: Order) -> Self {
        Self {
            terms: Vec::new(),
            nvars,
            order,
        }
    }

    pub fn constant(c: F, nvars: usize, order: Order) -> Self {
        Self::new(vec![(Monomial::one(nvars), c)], nvars, order)
    }

    pub fn var(i: usize, nvars: usize, order: Order) -> Self {
        Self::monomial(Monomial::var(i, nvars), order)
    }

    /// The monomial `m` with coefficient one.
    pub fn monomial(m: Monomial, order: Order) -> Self {
        let n = m.nvars();
        Self::new(vec![(m, F::one())], n, order)
    }

    fn with(&self, terms: Vec<Term<F>>) -> Self {
        Self {
            terms,
            nvars: self.nvars,
            order: self.order.clone(),
        }
    }

    pub const fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn is_constant(&self) -> bool {
        self.terms.iter().all(|t| t.0.is_one())
    }

    pub fn lt(&self) -> Option<&Term<F>> {
        self.terms.first()
    }

    pub fn lm(&self) -> Option<&Monomial> {
        self.lt().map(|t| &t.0)
    }

    pub fn lc(&self) -> Option<&F> {
        self.lt().map(|t| &t.1)
    }

    pub fn total_degree(&self) -> u32 {
        self.terms.iter().map(|t| t.0.degree()).max().unwrap_or(0)
    }

    pub fn degree(&self, k: usize) -> u32 {
        self.terms.iter().map(|t| t.0.exps()[k]).max().unwrap_or(0)
    }

    #[must_use]
    pub fn reorder(&self, order: Order) -> Self {
        Self::new(self.terms.clone(), self.nvars, order)
    }

    pub fn map<G: Field>(&self, f: impl Fn(&F) -> G) -> Poly<G> {
        let terms = self.terms.iter().map(|(m, c)| (m.clone(), f(c))).collect();
        Poly::new(terms, self.nvars, self.order.clone())
    }

    /// [`Poly::map`] by a partial map, such as reduction modulo a prime that may divide a
    /// denominator.
    pub fn try_map<G: Field>(&self, f: impl Fn(&F) -> Option<G>) -> Option<Poly<G>> {
        let terms = self.terms.iter().map(|(m, c)| Some((m.clone(), f(c)?)));
        Some(Poly::new(
            terms.collect::<Option<_>>()?,
            self.nvars,
            self.order.clone(),
        ))
    }

    pub fn support(&self) -> Vec<Monomial> {
        self.terms.iter().map(|t| t.0.clone()).collect()
    }

    /// The largest monomial dividing every term, `1` for zero.
    pub fn min_exps(&self) -> Monomial {
        let mut it = self.terms.iter().map(|t| t.0.exps());
        let Some(first) = it.next() else {
            return Monomial::one(self.nvars);
        };
        let e = it.fold(first.to_vec(), |mut acc, e| {
            acc.iter_mut().zip(e).for_each(|(a, b)| *a = (*a).min(*b));
            acc
        });
        Monomial::new(e)
    }

    /// Variable `i` of the result is variable `perm[i]` of `self`.
    #[must_use]
    pub fn permute(&self, perm: &[usize]) -> Self {
        let terms = self.terms.iter().map(|(m, c)| {
            let e: Vec<u32> = perm.iter().map(|&i| m.exps()[i]).collect();
            (Monomial::new(e), c.clone())
        });
        Self::new(terms.collect(), perm.len(), self.order.clone())
    }
    #[must_use]
    pub fn scale(&self, c: &F) -> Self {
        if c.is_zero() {
            return self.with(Vec::new());
        }
        self.with(
            self.terms
                .iter()
                .map(|(m, a)| (m.clone(), a.clone() * c.clone()))
                .collect(),
        )
    }

    #[must_use]
    pub fn mul_term(&self, c: &F, m: &Monomial) -> Self {
        let shifted = self
            .terms
            .iter()
            .map(|(n, a)| (n * m, a.clone() * c.clone()));
        self.with(shifted.filter(|t| !t.1.is_zero()).collect())
    }

    #[must_use]
    pub fn monic(&self) -> Self {
        self.lc()
            .and_then(Field::inverse)
            .map_or_else(|| self.clone(), |l| self.scale(&l))
    }

    /// `self - c * m * g` in one merge, shifting `g` as it goes.
    #[must_use]
    pub fn sub_mul(&self, c: &F, m: &Monomial, g: &Self) -> Self {
        debug_assert_eq!(self.nvars, g.nvars);
        self.with(merge(
            &self.order,
            self.terms.iter().cloned(),
            shifted(g, &-c.clone(), m),
        ))
    }

    fn merge(&self, o: &Self) -> Self {
        debug_assert_eq!(self.nvars, o.nvars);
        self.with(merge(
            &self.order,
            self.terms.iter().cloned(),
            o.terms.iter().cloned(),
        ))
    }

    #[must_use]
    pub fn pow(&self, mut e: u32) -> Self {
        let (mut base, mut acc) = (
            self.clone(),
            Self::constant(F::one(), self.nvars, self.order.clone()),
        );
        while e > 0 {
            if e & 1 == 1 {
                acc = &acc * &base;
            }
            e >>= 1;
            if e > 0 {
                base = &base * &base;
            }
        }
        acc
    }

    /// The partial derivative in `x_k`.
    #[must_use]
    pub fn derivative(&self, k: usize) -> Self {
        let terms = self
            .terms
            .iter()
            .filter(|t| t.0.exps()[k] > 0)
            .map(|(m, c)| {
                let mut e = m.exps().to_vec();
                let d = e[k];
                e[k] -= 1;
                (Monomial::new(e), c.clone() * nat(d.into()))
            });
        Self::new(terms.collect(), self.nvars, self.order.clone())
    }

    pub fn eval(&self, x: &[F]) -> F {
        self.terms
            .iter()
            .fold(F::zero(), |acc, (m, c)| acc + c.clone() * m.eval(x))
    }

    /// Substitutes `x[i]` for every `x_i` with `i != k`, leaving a polynomial in `x_k`.
    pub fn eval_except(&self, k: usize, x: &[F]) -> Uni<F> {
        let mut d = vec![F::zero(); self.degree(k) as usize + 1];
        for (m, c) in &self.terms {
            let v = m
                .exps()
                .iter()
                .enumerate()
                .filter(|&(i, &e)| i != k && e > 0);
            let v = v.fold(c.clone(), |v, (i, &e)| v * pow(&x[i], e.into()));
            let j = m.exps()[k] as usize;
            d[j] = d[j].clone() + v;
        }
        Uni::new(d)
    }

    /// The restriction to the line `x_i = z_i * t + s_i`, as a polynomial in `t`.
    pub fn on_line(&self, z: &[F], s: &[F]) -> Uni<F> {
        let lines: Vec<Uni<F>> = z
            .iter()
            .zip(s)
            .map(|(a, b)| Uni::new(vec![b.clone(), a.clone()]))
            .collect();
        let mut powers: Vec<Vec<Uni<F>>> = vec![vec![Uni::constant(F::one())]; lines.len()];
        let mut out = Uni::zero();
        for (m, c) in &self.terms {
            let mut t = Uni::constant(c.clone());
            for (i, &e) in m.exps().iter().enumerate().filter(|t| *t.1 > 0) {
                while powers[i].len() <= e as usize {
                    let next = &powers[i][powers[i].len() - 1] * &lines[i];
                    powers[i].push(next);
                }
                t = &t * &powers[i][e as usize];
            }
            out = &out + &t;
        }
        out
    }

    /// The coefficients in `x_k`: for each monomial in the other variables, descending, the
    /// polynomial in `x_k` multiplying it.
    pub fn coeffs_in(&self, k: usize) -> Vec<(Monomial, Uni<F>)> {
        let mut groups: BTreeMap<Monomial, Vec<F>> = BTreeMap::new();
        for (m, c) in &self.terms {
            let mut e = m.exps().to_vec();
            let i = std::mem::take(&mut e[k]) as usize;
            let v = groups.entry(Monomial::new(e)).or_default();
            if v.len() <= i {
                v.resize(i + 1, F::zero());
            }
            v[i] = c.clone();
        }
        let mut out: Vec<(Monomial, Uni<F>)> =
            groups.into_iter().map(|(m, v)| (m, Uni::new(v))).collect();
        out.sort_by(|a, b| self.order.compare(&b.0, &a.0));
        out
    }

    /// The inverse of [`Poly::coeffs_in`].
    pub fn from_coeffs_in(
        k: usize,
        nvars: usize,
        order: Order,
        groups: impl IntoIterator<Item = (Monomial, Uni<F>)>,
    ) -> Self {
        let mut terms = Vec::new();
        for (m, u) in groups {
            for (i, c) in u.0.into_iter().enumerate() {
                let mut e = m.exps().to_vec();
                e[k] = i as u32;
                terms.push((Monomial::new(e), c));
            }
        }
        Self::new(terms, nvars, order)
    }

    /// The content in `F[x_k]`: the monic gcd of the coefficients in [`Poly::coeffs_in`].
    pub fn content_in(&self, k: usize) -> Uni<F> {
        self.coeffs_in(k)
            .iter()
            .fold(Uni::zero(), |g, (_, u)| g.gcd(u))
    }

    /// The content in `F[x_k]` and the primitive part.
    pub fn primitive_in(&self, k: usize) -> (Uni<F>, Self) {
        let c = self.content_in(k);
        if c.deg() == 0 {
            return (c, self.clone());
        }
        let groups = self.coeffs_in(k).into_iter().map(|(m, u)| (m, &u / &c));
        let p = Self::from_coeffs_in(k, self.nvars, self.order.clone(), groups);
        (c, p)
    }

    /// Substitutes `a` for `x_k`.
    #[must_use]
    pub fn eval_var(&self, k: usize, a: &F) -> Self {
        let terms = self.terms.iter().map(|(m, c)| {
            let mut e = m.exps().to_vec();
            let v = c.clone() * pow(a, std::mem::take(&mut e[k]).into());
            (Monomial::new(e), v)
        });
        Self::new(terms.collect(), self.nvars, self.order.clone())
    }

    /// The S-polynomial, zero if either is.
    #[must_use]
    pub fn spoly(&self, o: &Self) -> Self {
        let (Some((m1, c1)), Some((m2, c2))) = (self.lt(), o.lt()) else {
            return self.with(Vec::new());
        };
        let l = m1.lcm(m2);
        let a = self.mul_term(&(F::one() / c1.clone()), &l.quo(m1).expect("lcm"));
        a.sub_mul(&(F::one() / c2.clone()), &l.quo(m2).expect("lcm"), o)
    }

    /// Division by an ordered list: `self = sum q[i] * divs[i] + r`, no term of `r` divisible by a
    /// leading monomial, each step reducing by the first divisor that applies.
    pub fn divide(&self, divs: &[Self]) -> (Vec<Self>, Self) {
        let mut q: Vec<Vec<Term<F>>> = vec![Vec::new(); divs.len()];
        let (mut rem, mut work, mut i) = (Vec::new(), self.terms.clone(), 0);
        while let Some((m, c)) = work.get(i) {
            let hit = divs.iter().enumerate().find_map(|(j, g)| {
                let (gm, gc) = g.lt()?;
                Some((j, m.quo(gm)?, c.clone() / gc.clone()))
            });
            let Some((j, t, k)) = hit else {
                i += 1;
                continue;
            };
            rem.extend(work.drain(..i));
            work = merge(
                &self.order,
                work.into_iter(),
                shifted(&divs[j], &-k.clone(), &t),
            );
            i = 0;
            q[j].push((t, k));
        }
        rem.extend(work);
        let q = q.into_iter().map(|t| self.with(t)).collect();
        (q, self.with(rem))
    }

    /// The remainder of [`Poly::divide`].
    #[must_use]
    pub fn reduce(&self, divs: &[Self]) -> Self {
        self.divide(divs).1
    }
}

/// The terms of `c * m * g`, still descending.
fn shifted<'a, F: Field>(
    g: &'a Poly<F>,
    c: &'a F,
    m: &'a Monomial,
) -> impl Iterator<Item = Term<F>> + 'a {
    let live = !c.is_zero();
    g.terms
        .iter()
        .filter(move |_| live)
        .map(move |(n, a)| (n * m, a.clone() * c.clone()))
}

/// Merges two descending term streams, adding like terms and dropping cancellations.
fn merge<F: Field>(
    order: &Order,
    a: impl Iterator<Item = Term<F>>,
    b: impl Iterator<Item = Term<F>>,
) -> Vec<Term<F>> {
    let (mut a, mut b) = (a.peekable(), b.peekable());
    let mut out = Vec::with_capacity(a.size_hint().0 + b.size_hint().0);
    loop {
        let side = match (a.peek(), b.peek()) {
            (Some(x), Some(y)) => order.compare(&x.0, &y.0),
            (Some(_), None) => Ordering::Greater,
            (None, Some(_)) => Ordering::Less,
            (None, None) => return out,
        };
        match side {
            Ordering::Greater => out.extend(a.next()),
            Ordering::Less => out.extend(b.next()),
            Ordering::Equal => {
                if let (Some((m, x)), Some((_, y))) = (a.next(), b.next()) {
                    let c = x + y;
                    if !c.is_zero() {
                        out.push((m, c));
                    }
                }
            }
        }
    }
}

/// With variables named `x0, x1, ...`.
impl<F: Field + Display> Display for Poly<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let names = (0..self.nvars).map(|i| format!("x{i}"));
        f.write_str(&Ring::new(names, self.order.clone()).show(self))
    }
}

/// `sum h[i] * gens[i]`, the polynomial a certificate claims.
pub fn combination<F: Field>(h: &[Poly<F>], gens: &[Poly<F>]) -> Poly<F> {
    let zero = gens.first().map(|g| Poly::zero(g.nvars, g.order.clone()));
    h.iter().zip(gens).fold(
        zero.unwrap_or_else(|| Poly::zero(0, Order::Lex)),
        |acc, (a, g)| &acc + &(a * g),
    )
}

impl<F: Field> Add for &Poly<F> {
    type Output = Poly<F>;
    fn add(self, o: Self) -> Poly<F> {
        self.merge(o)
    }
}

impl<F: Field> Neg for &Poly<F> {
    type Output = Poly<F>;
    fn neg(self) -> Poly<F> {
        self.with(
            self.terms
                .iter()
                .map(|(m, c)| (m.clone(), -c.clone()))
                .collect(),
        )
    }
}

impl<F: Field> Sub for &Poly<F> {
    type Output = Poly<F>;
    fn sub(self, o: Self) -> Poly<F> {
        self.merge(&-o)
    }
}

impl<F: Field> Mul for &Poly<F> {
    type Output = Poly<F>;
    fn mul(self, o: Self) -> Poly<F> {
        o.terms.iter().fold(self.with(Vec::new()), |acc, (m, c)| {
            self.with(merge(
                &self.order,
                acc.terms.into_iter(),
                shifted(self, c, m),
            ))
        })
    }
}
