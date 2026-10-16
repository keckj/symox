//! Rational functions `num/den` over `ZZ` (the fraction field `QQ(a, b, ...)`
//! used for coefficients with symbolic parameters, e.g. in `apart`).

use std::fmt;

use super::gcd::gcd_mz;
use super::int::{Coeff, Field, Int, Rat};
use super::mpoly::MPoly;

/// A reduced fraction of integer polynomials with a positive leading
/// coefficient in the denominator.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct RF {
    pub n: MPoly<Int>,
    pub d: MPoly<Int>,
}

impl fmt::Debug for RF {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:?})/({:?})", self.n, self.d)
    }
}

thread_local! {
    static NVARS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Number of variables of the rational functions created by `zero()`,
/// `one()` and `from_i64()` (set by the code using the field).
pub fn set_nvars(n: usize) -> usize {
    NVARS.with(|c| c.replace(n))
}

/// Restore the ambient coefficient dimension even when an operation exits
/// early. Nested polynomial operations may safely use a different field.
pub struct Scope(usize);
impl Scope {
    pub fn new(n: usize) -> Self {
        Self(set_nvars(n))
    }
}
impl Drop for Scope {
    fn drop(&mut self) {
        set_nvars(self.0);
    }
}

fn nv() -> usize {
    NVARS.with(|c| c.get())
}

impl RF {
    pub fn new(n: MPoly<Int>, d: MPoly<Int>) -> RF {
        assert!(!d.is_zero(), "zero denominator");
        if n.is_zero() {
            return RF {
                d: MPoly::one(n.nvars),
                n,
            };
        }
        let (g, a, b) = gcd_mz(&n, &d);
        let (a, b) = if g.is_one() { (n, d) } else { (a, b) };
        if b.lc().is_negative() {
            RF {
                n: a.neg(),
                d: b.neg(),
            }
        } else {
            RF { n: a, d: b }
        }
    }

    pub fn from_poly(p: MPoly<Int>) -> RF {
        let n = p.nvars;
        RF {
            n: p,
            d: MPoly::one(n),
        }
    }

    pub fn from_rat_poly(p: &MPoly<Rat>) -> RF {
        let (den, z) = p.clear_denoms();
        RF::new(z, MPoly::constant(den, p.nvars))
    }
}

impl Coeff for RF {
    fn zero() -> Self {
        RF {
            n: MPoly::zero(nv()),
            d: MPoly::one(nv()),
        }
    }
    fn one() -> Self {
        RF {
            n: MPoly::one(nv()),
            d: MPoly::one(nv()),
        }
    }
    fn from_i64(i: i64) -> Self {
        RF {
            n: MPoly::constant(Int::S(i), nv()),
            d: MPoly::one(nv()),
        }
    }
    fn is_zero(&self) -> bool {
        self.n.is_zero()
    }
    fn is_one(&self) -> bool {
        self.n.is_one() && self.d.is_one()
    }
    fn add(&self, o: &Self) -> Self {
        if self.n.is_zero() {
            return o.clone();
        }
        if o.n.is_zero() {
            return self.clone();
        }
        if self.d == o.d {
            return RF::new(self.n.add(&o.n), self.d.clone());
        }
        RF::new(self.n.mul(&o.d).add(&o.n.mul(&self.d)), self.d.mul(&o.d))
    }
    fn sub(&self, o: &Self) -> Self {
        self.add(&o.neg())
    }
    fn mul(&self, o: &Self) -> Self {
        if self.n.is_zero() || o.n.is_zero() {
            return RF {
                n: MPoly::zero(self.n.nvars),
                d: MPoly::one(self.n.nvars),
            };
        }
        RF::new(self.n.mul(&o.n), self.d.mul(&o.d))
    }
    fn neg(&self) -> Self {
        RF {
            n: self.n.neg(),
            d: self.d.clone(),
        }
    }
    fn is_negative(&self) -> bool {
        self.n.lc().is_negative()
    }
    fn div_exact(&self, o: &Self) -> Option<Self> {
        if o.n.is_zero() {
            None
        } else {
            Some(self.mul(&o.inv()))
        }
    }
    fn to_rat(&self) -> Rat {
        match (self.n.constant_value(), self.d.constant_value()) {
            (Some(a), Some(b)) => Rat::new(a, b),
            _ => panic!("not a constant rational function"),
        }
    }
}

impl Field for RF {
    fn inv(&self) -> Self {
        assert!(!self.n.is_zero(), "division by zero");
        if self.n.lc().is_negative() {
            RF {
                n: self.d.neg(),
                d: self.n.neg(),
            }
        } else {
            RF {
                n: self.d.clone(),
                d: self.n.clone(),
            }
        }
    }
}
