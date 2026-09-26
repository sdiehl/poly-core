use std::fmt::{self, Display};

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;

use crate::field::Field;
use crate::monomial::{Monomial, Order};
use crate::poly::Poly;

type Q = BigRational;

/// Variable names and an order: enough to read and write polynomials.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ring {
    pub names: Vec<String>,
    pub order: Order,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParseError {
    pub pos: usize,
    pub msg: String,
}

impl Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "at {}: {}", self.pos, self.msg)
    }
}

impl std::error::Error for ParseError {}

impl Ring {
    pub fn new<S: Into<String>>(names: impl IntoIterator<Item = S>, order: Order) -> Self {
        Self {
            names: names.into_iter().map(Into::into).collect(),
            order,
        }
    }

    pub const fn nvars(&self) -> usize {
        self.names.len()
    }

    #[must_use]
    pub fn with_order(&self, order: Order) -> Self {
        Self {
            names: self.names.clone(),
            order,
        }
    }

    /// Parses `+ - * / ^`, parentheses, integers and variable names over Q. Map the result into
    /// another field with [`Poly::map`].
    pub fn parse(&self, src: &str) -> Result<Poly<Q>, ParseError> {
        let mut p = Parser {
            ring: self,
            src: src.as_bytes(),
            pos: 0,
        };
        let e = p.expr()?;
        match p.peek() {
            None => Ok(e),
            Some(_) => p.fail("unexpected input"),
        }
    }

    pub fn monomial(&self, m: &Monomial) -> String {
        let vars = m.exps().iter().zip(&self.names).filter(|t| *t.0 > 0);
        let parts: Vec<String> = vars
            .map(|(&e, v)| {
                if e == 1 {
                    v.clone()
                } else {
                    format!("{v}^{e}")
                }
            })
            .collect();
        if parts.is_empty() {
            "1".into()
        } else {
            parts.join("*")
        }
    }

    pub fn show<F: Field + Display>(&self, p: &Poly<F>) -> String {
        let mut out = String::new();
        for (i, (m, c)) in p.terms.iter().enumerate() {
            let s = c.to_string();
            let (neg, mag) = s
                .strip_prefix('-')
                .map_or((false, s.as_str()), |r| (true, r));
            let body = match (m.is_one(), mag) {
                (true, _) => mag.to_string(),
                (false, "1") => self.monomial(m),
                (false, _) => format!("{mag}*{}", self.monomial(m)),
            };
            let sign = match (i, neg) {
                (0, false) => "",
                (0, true) => "-",
                (_, false) => " + ",
                (_, true) => " - ",
            };
            out.push_str(sign);
            out.push_str(&body);
        }
        if out.is_empty() {
            "0".into()
        } else {
            out
        }
    }
}

struct Parser<'a> {
    ring: &'a Ring,
    src: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn fail<T>(&self, msg: &str) -> Result<T, ParseError> {
        Err(ParseError {
            pos: self.pos,
            msg: msg.into(),
        })
    }

    fn peek(&mut self) -> Option<u8> {
        while self.src.get(self.pos).is_some_and(u8::is_ascii_whitespace) {
            self.pos += 1;
        }
        self.src.get(self.pos).copied()
    }

    fn eat(&mut self, b: u8) -> bool {
        let hit = self.peek() == Some(b);
        self.pos += usize::from(hit);
        hit
    }

    fn take(&mut self, f: impl Fn(u8) -> bool) -> &str {
        let start = self.pos;
        while self.src.get(self.pos).is_some_and(|&b| f(b)) {
            self.pos += 1;
        }
        std::str::from_utf8(&self.src[start..self.pos]).unwrap_or_default()
    }

    fn constant(&self, q: Q) -> Poly<Q> {
        Poly::constant(q, self.ring.nvars(), self.ring.order.clone())
    }

    fn expr(&mut self) -> Result<Poly<Q>, ParseError> {
        let mut acc = self.term()?;
        loop {
            if self.eat(b'+') {
                acc = &acc + &self.term()?;
            } else if self.eat(b'-') {
                acc = &acc - &self.term()?;
            } else {
                return Ok(acc);
            }
        }
    }

    fn term(&mut self) -> Result<Poly<Q>, ParseError> {
        let mut acc = self.factor()?;
        loop {
            if self.eat(b'*') {
                acc = &acc * &self.factor()?;
            } else if self.eat(b'/') {
                let d = self.factor()?;
                match (d.is_constant(), d.lc()) {
                    (true, Some(c)) => acc = acc.scale(&c.recip()),
                    _ => return self.fail("division by a nonconstant or zero"),
                }
            } else {
                return Ok(acc);
            }
        }
    }

    fn factor(&mut self) -> Result<Poly<Q>, ParseError> {
        if self.eat(b'-') {
            return Ok(-&self.factor()?);
        }
        let base = self.atom()?;
        if !self.eat(b'^') {
            return Ok(base);
        }
        self.peek();
        let Ok(e) = self.take(|b| b.is_ascii_digit()).parse::<u32>() else {
            return self.fail("expected an exponent");
        };
        Ok((0..e).fold(self.constant(Q::one()), |acc, _| &acc * &base))
    }

    fn atom(&mut self) -> Result<Poly<Q>, ParseError> {
        match self.peek() {
            Some(b'(') => {
                self.pos += 1;
                let e = self.expr()?;
                if self.eat(b')') {
                    Ok(e)
                } else {
                    self.fail("expected )")
                }
            }
            Some(b) if b.is_ascii_digit() => {
                let n: BigInt = self.take(|b| b.is_ascii_digit()).parse().expect("digits");
                Ok(self.constant(Q::from_integer(n)))
            }
            Some(b) if b.is_ascii_alphabetic() || b == b'_' => {
                let name = self
                    .take(|b| b.is_ascii_alphanumeric() || b == b'_')
                    .to_string();
                match self.ring.names.iter().position(|v| *v == name) {
                    Some(i) => Ok(Poly::var(i, self.ring.nvars(), self.ring.order.clone())),
                    None => self.fail(&format!("unknown variable {name}")),
                }
            }
            _ => self.fail("expected a number, variable or ("),
        }
    }
}
