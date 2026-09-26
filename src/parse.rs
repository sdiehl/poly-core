use std::fmt::{self, Display};

use num_bigint::BigInt;
use num_rational::BigRational;

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
    /// Panics on names [`Ring::try_new`] rejects.
    pub fn new<S: Into<String>>(names: impl IntoIterator<Item = S>, order: Order) -> Self {
        Self::try_new(names, order).unwrap_or_else(|e| panic!("{e}"))
    }

    /// Rejects empty, duplicate and non-identifier names.
    pub fn try_new<S: Into<String>>(
        names: impl IntoIterator<Item = S>,
        order: Order,
    ) -> Result<Self, ParseError> {
        let names: Vec<String> = names.into_iter().map(Into::into).collect();
        for (i, v) in names.iter().enumerate() {
            let fail = |msg: String| Err(ParseError { pos: i, msg });
            if !is_ident(v) {
                return fail(format!("invalid variable name {v:?}"));
            }
            if names[..i].contains(v) {
                return fail(format!("duplicate variable {v}"));
            }
        }
        Ok(Self { names, order })
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
        self.parse_with(src, &Q::from_integer, &[])
    }

    /// Several polynomials separated by `,` or `;`, blanks skipped.
    pub fn parse_many(&self, src: &str) -> Result<Vec<Poly<Q>>, ParseError> {
        let mut out = Vec::new();
        let mut start = 0;
        for piece in src.split([',', ';']) {
            if !piece.trim().is_empty() {
                let shift = |e: ParseError| ParseError {
                    pos: e.pos + start,
                    ..e
                };
                out.push(self.parse(piece).map_err(shift)?);
            }
            start += piece.len() + 1;
        }
        Ok(out)
    }

    /// Parses over any field, reading integers through `lift` and the names in `params` as the
    /// given constants, such as the generators of `Q(a)` as [`crate::RatFunc`]. Division is by
    /// nonzero constants only.
    pub fn parse_with<F: Field>(
        &self,
        src: &str,
        lift: &dyn Fn(BigInt) -> F,
        params: &[(&str, F)],
    ) -> Result<Poly<F>, ParseError> {
        let mut p = Parser {
            ring: self,
            lift,
            params,
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

    /// A monomial in LaTeX: `x_{1}^{2} y`.
    pub fn latex_monomial(&self, m: &Monomial) -> String {
        let vars = m.exps().iter().zip(&self.names).filter(|t| *t.0 > 0);
        let parts: Vec<String> = vars
            .map(|(&e, v)| match e {
                1 => latex_name(v),
                _ => format!("{}^{{{e}}}", latex_name(v)),
            })
            .collect();
        parts.join(" ")
    }

    /// `p` in LaTeX, fractions as `\frac`.
    pub fn latex<F: Field + Display>(&self, p: &Poly<F>) -> String {
        let mut out = String::new();
        for (i, (m, c)) in p.terms.iter().enumerate() {
            let s = c.to_string();
            let (neg, mag) = s
                .strip_prefix('-')
                .map_or((false, s.as_str()), |r| (true, r));
            let coeff = match mag.split_once('/') {
                Some((a, b)) => format!("\\frac{{{a}}}{{{b}}}"),
                None => mag.to_string(),
            };
            let body = match (m.is_one(), mag) {
                (true, _) => coeff,
                (false, "1") => self.latex_monomial(m),
                (false, _) => format!("{coeff} {}", self.latex_monomial(m)),
            };
            out.push_str(match (i, neg) {
                (0, false) => "",
                (0, true) => "-",
                (_, false) => " + ",
                (_, true) => " - ",
            });
            out.push_str(&body);
        }
        if out.is_empty() {
            "0".into()
        } else {
            out
        }
    }
}

fn is_ident(v: &str) -> bool {
    let mut b = v.bytes();
    b.next()
        .is_some_and(|c| c.is_ascii_alphabetic() || c == b'_')
        && b.all(|c| c.is_ascii_alphanumeric() || c == b'_')
}

/// Trailing digits, after an optional `_`, become a subscript; other underscores are escaped.
fn latex_name(v: &str) -> String {
    let stem = v.trim_end_matches(|c: char| c.is_ascii_digit());
    let digits = &v[stem.len()..];
    let stem = stem
        .strip_suffix('_')
        .filter(|_| !digits.is_empty())
        .unwrap_or(stem);
    let stem = stem.replace('_', "\\_");
    if digits.is_empty() || stem.is_empty() {
        format!("{stem}{digits}")
    } else {
        format!("{stem}_{{{digits}}}")
    }
}

struct Parser<'a, F> {
    ring: &'a Ring,
    lift: &'a dyn Fn(BigInt) -> F,
    params: &'a [(&'a str, F)],
    src: &'a [u8],
    pos: usize,
}

impl<F: Field> Parser<'_, F> {
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

    fn constant(&self, q: F) -> Poly<F> {
        Poly::constant(q, self.ring.nvars(), self.ring.order.clone())
    }

    fn expr(&mut self) -> Result<Poly<F>, ParseError> {
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

    fn term(&mut self) -> Result<Poly<F>, ParseError> {
        let mut acc = self.factor()?;
        loop {
            if self.eat(b'*') {
                acc = &acc * &self.factor()?;
            } else if self.eat(b'/') {
                let d = self.factor()?;
                match d.lc().filter(|_| d.is_constant()).and_then(F::inverse) {
                    Some(c) => acc = acc.scale(&c),
                    None => return self.fail("division by a nonconstant or zero"),
                }
            } else {
                return Ok(acc);
            }
        }
    }

    fn factor(&mut self) -> Result<Poly<F>, ParseError> {
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
        Ok(base.pow(e))
    }

    fn atom(&mut self) -> Result<Poly<F>, ParseError> {
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
                Ok(self.constant((self.lift)(n)))
            }
            Some(b) if b.is_ascii_alphabetic() || b == b'_' => {
                let name = self
                    .take(|b| b.is_ascii_alphanumeric() || b == b'_')
                    .to_string();
                if let Some(i) = self.ring.names.iter().position(|v| *v == name) {
                    return Ok(Poly::var(i, self.ring.nvars(), self.ring.order.clone()));
                }
                match self.params.iter().find(|(v, _)| *v == name) {
                    Some((_, c)) => Ok(self.constant(c.clone())),
                    None => self.fail(&format!("unknown variable {name}")),
                }
            }
            _ => self.fail("expected a number, variable or ("),
        }
    }
}
