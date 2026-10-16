//! Series coefficients: exact numbers (fast path) or general expressions.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Zero};

use crate::expr::{Expr, Kind};
use crate::number::Number;
use crate::ops::{add, add2, mul2, neg, pow};

/// A coefficient of a truncated series.
///
/// `N` holds plain numbers (the common case for elementary functions
/// expanded at rational points), `E` any other expression. An `E` never
/// holds a bare number: constructors normalise through [`C::from_expr`].
#[derive(Clone, Debug, PartialEq)]
pub enum C {
    N(Number),
    E(Expr),
}

impl C {
    #[inline]
    pub fn zero() -> C {
        C::N(Number::ZERO)
    }

    #[inline]
    pub fn one() -> C {
        C::N(Number::ONE)
    }

    #[inline]
    pub fn int(i: i64) -> C {
        C::N(Number::int(i))
    }

    pub fn from_expr(e: Expr) -> C {
        match e.kind() {
            Kind::Num(n) => C::N(n.clone()),
            _ => C::E(e),
        }
    }

    pub fn to_expr(&self) -> Expr {
        match self {
            C::N(n) => Expr::num(n.clone()),
            C::E(e) => e.clone(),
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        match self {
            C::N(n) => n.is_zero(),
            C::E(_) => false,
        }
    }

    #[inline]
    pub fn is_one(&self) -> bool {
        matches!(self, C::N(n) if n.is_one())
    }

    pub fn as_number(&self) -> Option<&Number> {
        match self {
            C::N(n) => Some(n),
            C::E(_) => None,
        }
    }

    pub fn add(&self, o: &C) -> C {
        match (self, o) {
            (C::N(a), C::N(b)) => C::N(nadd(a, b)),
            (C::N(a), _) if a.is_zero() => o.clone(),
            (_, C::N(b)) if b.is_zero() => self.clone(),
            _ => C::from_expr(add2(&self.to_expr(), &o.to_expr())),
        }
    }

    pub fn sub(&self, o: &C) -> C {
        self.add(&o.neg())
    }

    pub fn neg(&self) -> C {
        match self {
            C::N(a) => C::N(a.neg()),
            C::E(e) => C::from_expr(neg(e)),
        }
    }

    pub fn mul(&self, o: &C) -> C {
        match (self, o) {
            (C::N(a), C::N(b)) => C::N(nmul(a, b)),
            (C::N(a), C::E(e)) | (C::E(e), C::N(a)) => {
                if a.is_zero() {
                    C::zero()
                } else if a.is_one() {
                    C::E(e.clone())
                } else {
                    C::from_expr(mul2(&Expr::num(a.clone()), e))
                }
            }
            (C::E(a), C::E(b)) => C::from_expr(emul(a, b)),
        }
    }

    pub fn mul_num(&self, k: &Number) -> C {
        match self {
            C::N(a) => C::N(nmul(a, k)),
            C::E(e) => {
                if k.is_zero() {
                    C::zero()
                } else if k.is_one() {
                    self.clone()
                } else {
                    C::from_expr(mul2(&Expr::num(k.clone()), e))
                }
            }
        }
    }

    pub fn mul_int(&self, k: i64) -> C {
        match self {
            C::N(a) => match a {
                Number::Int(i) => match i.checked_mul(k) {
                    Some(v) => C::N(Number::Int(v)),
                    None => C::N(a.mul(&Number::int(k))),
                },
                _ => C::N(nmul(a, &Number::int(k))),
            },
            _ => self.mul_num(&Number::int(k)),
        }
    }

    /// `self / k` for a non-zero integer `k`.
    pub fn div_int(&self, k: i64) -> C {
        if k == 1 {
            return self.clone();
        }
        let q = Number::rational_i64(1, k).expect("division by zero");
        self.mul_num(&q)
    }

    /// Multiplicative inverse, `None` for zero.
    pub fn inv(&self) -> Option<C> {
        match self {
            C::N(n) => n.inv().map(C::N),
            C::E(e) => Some(C::from_expr(pow(e, &Expr::minus_one()))),
        }
    }

    /// `self ** r`.
    pub fn pow(&self, r: &C) -> C {
        match (self, r) {
            (C::N(a), C::N(k)) => {
                if let Some(i) = k.to_i64() {
                    if let Some(v) = a.pow_i64(i) {
                        return C::N(v);
                    }
                }
                C::from_expr(pow(&Expr::num(a.clone()), &Expr::num(k.clone())))
            }
            _ => C::from_expr(pow(&self.to_expr(), &r.to_expr())),
        }
    }
}

/// Product of two expressions distributed one level over sums, so that
/// coefficients stay in a canonical expanded form (needed to detect
/// cancellation).
pub fn emul(a: &Expr, b: &Expr) -> Expr {
    if !a.is_add() && !b.is_add() {
        return mul2(a, b);
    }
    let ta = a.args_add();
    let tb = b.args_add();
    let mut v = Vec::with_capacity(ta.len() * tb.len());
    for x in &ta {
        for y in &tb {
            v.push(mul2(x, y));
        }
    }
    add(&v)
}

trait AddArgs {
    fn args_add(&self) -> Vec<Expr>;
}

impl AddArgs for Expr {
    fn args_add(&self) -> Vec<Expr> {
        if self.is_add() {
            self.args()
        } else {
            vec![self.clone()]
        }
    }
}

/// `a + b` without the clones of `Number::add` for rationals.
pub fn nadd(a: &Number, b: &Number) -> Number {
    match (a, b) {
        (Number::Int(x), Number::Int(y)) => match x.checked_add(*y) {
            Some(v) => Number::Int(v),
            None => a.add(b),
        },
        (Number::Rational(x), Number::Rational(y)) => Number::from_ratio(x + y),
        (Number::Rational(x), Number::Int(i)) | (Number::Int(i), Number::Rational(x)) => {
            if *i == 0 {
                return Number::Rational(x.clone());
            }
            // (p + i q)/q is reduced when p/q is
            let num = x.numer() + x.denom() * BigInt::from(*i);
            Number::Rational(BigRational::new_raw(num, x.denom().clone()))
        }
        _ => a.add(b),
    }
}

/// `a * b` without the clones of `Number::mul` for rationals.
pub fn nmul(a: &Number, b: &Number) -> Number {
    match (a, b) {
        (Number::Int(x), Number::Int(y)) => match x.checked_mul(*y) {
            Some(v) => Number::Int(v),
            None => a.mul(b),
        },
        (Number::Rational(x), Number::Rational(y)) => Number::from_ratio(x * y),
        (Number::Rational(x), Number::Int(i)) | (Number::Int(i), Number::Rational(x)) => {
            if *i == 0 {
                return Number::ZERO;
            }
            if *i == 1 {
                return Number::Rational(x.clone());
            }
            let bi = BigInt::from(*i);
            let g = bi.gcd(x.denom());
            let num = x.numer() * (&bi / &g);
            let den = x.denom() / &g;
            if den.is_one() {
                Number::from_bigint(num)
            } else if num.is_zero() {
                Number::ZERO
            } else {
                Number::Rational(BigRational::new_raw(num, den))
            }
        }
        _ => a.mul(b),
    }
}
