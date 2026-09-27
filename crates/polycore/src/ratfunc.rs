use std::ops::{Add, Div, Mul, Neg, Sub};

use num_traits::{One, Zero};

use crate::field::Field;
use crate::uni::Uni;

/// A univariate rational function `num / den` in lowest terms with `den` monic: the field `F(a)`,
/// coefficients for polynomials with a parameter.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct RatFunc<F> {
    num: Uni<F>,
    den: Uni<F>,
}

impl<F: Field> RatFunc<F> {
    pub fn new(num: &Uni<F>, den: &Uni<F>) -> Self {
        assert!(!den.is_zero(), "zero denominator");
        let g = num.gcd(den);
        let (num, den) = (num / &g, den / &g);
        let l = den.lc().inverse().expect("nonzero denominator");
        Self {
            num: num.scale(&l),
            den: den.scale(&l),
        }
    }

    /// The parameter itself.
    pub fn var() -> Self {
        Uni::x().into()
    }

    pub const fn num(&self) -> &Uni<F> {
        &self.num
    }

    pub const fn den(&self) -> &Uni<F> {
        &self.den
    }

    /// The value at `a`, or `None` at a pole.
    pub fn eval(&self, a: &F) -> Option<F> {
        Some(self.num.eval(a) * self.den.eval(a).inverse()?)
    }

    /// Maps the coefficients, or `None` if `f` refuses one or the denominator vanishes.
    pub fn try_map<G: Field>(&self, f: impl Fn(&F) -> Option<G>) -> Option<RatFunc<G>> {
        let lift = |u: &Uni<F>| u.0.iter().map(&f).collect::<Option<Vec<G>>>().map(Uni::new);
        let (num, den) = (lift(&self.num)?, lift(&self.den)?);
        (!den.is_zero()).then(|| RatFunc::new(&num, &den))
    }
}

impl<F: Field> From<Uni<F>> for RatFunc<F> {
    fn from(num: Uni<F>) -> Self {
        Self {
            num,
            den: Uni::constant(F::one()),
        }
    }
}

impl<F: Field> Zero for RatFunc<F> {
    fn zero() -> Self {
        Uni::zero().into()
    }

    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

impl<F: Field> One for RatFunc<F> {
    fn one() -> Self {
        Uni::constant(F::one()).into()
    }
}

impl<F: Field> Add for RatFunc<F> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        if self.den == o.den {
            return Self::new(&(&self.num + &o.num), &self.den);
        }
        let num = &(&self.num * &o.den) + &(&o.num * &self.den);
        Self::new(&num, &(&self.den * &o.den))
    }
}

impl<F: Field> Neg for RatFunc<F> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            num: -&self.num,
            den: self.den,
        }
    }
}

impl<F: Field> Sub for RatFunc<F> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + -o
    }
}

impl<F: Field> Mul for RatFunc<F> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(&(&self.num * &o.num), &(&self.den * &o.den))
    }
}

impl<F: Field> Div for RatFunc<F> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        assert!(!o.is_zero(), "division by zero");
        Self::new(&(&self.num * &o.den), &(&self.den * &o.num))
    }
}
