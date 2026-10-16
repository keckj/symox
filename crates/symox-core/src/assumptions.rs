//! Three-valued assumption queries (`is_positive`, `is_real`, ...).
//!
//! The inference propagates facts from symbol assumptions through sums,
//! products, powers and elementary functions, like SymPy's old assumption
//! system (`expr.is_positive`). `None` means "unknown".

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::symbol::Assumptions as A;

/// Possible signs of a real quantity, as a set: bit 0 = negative,
/// bit 1 = zero, bit 2 = positive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Signs(u8);

impl Signs {
    pub const NEG: Signs = Signs(1);
    pub const ZERO: Signs = Signs(2);
    pub const POS: Signs = Signs(4);
    pub const ANY: Signs = Signs(7);
    pub const NONNEG: Signs = Signs(6);
    pub const NONPOS: Signs = Signs(3);
    pub const NONZERO: Signs = Signs(5);

    fn union(self, o: Signs) -> Signs {
        Signs(self.0 | o.0)
    }
    fn subset(self, o: Signs) -> bool {
        self.0 & !o.0 == 0
    }
    fn has(self, o: Signs) -> bool {
        self.0 & o.0 != 0
    }

    fn mul(self, o: Signs) -> Signs {
        let mut r = 0u8;
        for a in [Signs::NEG, Signs::ZERO, Signs::POS] {
            if !self.has(a) {
                continue;
            }
            for b in [Signs::NEG, Signs::ZERO, Signs::POS] {
                if !o.has(b) {
                    continue;
                }
                r |= match (a.0, b.0) {
                    (2, _) | (_, 2) => 2,
                    (1, 1) | (4, 4) => 4,
                    _ => 1,
                };
            }
        }
        Signs(r)
    }

    fn add(self, o: Signs) -> Signs {
        // a + b
        if self == Signs::ZERO {
            return o;
        }
        if o == Signs::ZERO {
            return self;
        }
        let nonneg = self.subset(Signs::NONNEG) && o.subset(Signs::NONNEG);
        let nonpos = self.subset(Signs::NONPOS) && o.subset(Signs::NONPOS);
        if nonneg {
            if self == Signs::POS || o == Signs::POS || !(self.has(Signs::ZERO) && o.has(Signs::ZERO)) {
                return Signs::POS;
            }
            return Signs::NONNEG;
        }
        if nonpos {
            if self == Signs::NEG || o == Signs::NEG || !(self.has(Signs::ZERO) && o.has(Signs::ZERO)) {
                return Signs::NEG;
            }
            return Signs::NONPOS;
        }
        Signs::ANY
    }
}

fn number_sign(n: &Number) -> Signs {
    if n.is_positive() {
        Signs::POS
    } else if n.is_negative() {
        Signs::NEG
    } else {
        Signs::ZERO
    }
}

/// Sign set of `e` if `e` is known to be real (finite), else `None`.
pub fn real_signs(e: &Expr) -> Option<Signs> {
    match e.kind() {
        Kind::Num(n) => Some(number_sign(n)),
        Kind::Const(c) => match c {
            Constant::Pi | Constant::E | Constant::EulerGamma => Some(Signs::POS),
            _ => None,
        },
        Kind::Sym(s) => {
            let a = s.assumptions();
            if !a.has(A::REAL) {
                return None;
            }
            Some(if a.has(A::ZERO) {
                Signs::ZERO
            } else if a.has(A::POSITIVE) {
                Signs::POS
            } else if a.has(A::NEGATIVE) {
                Signs::NEG
            } else if a.has(A::NONNEGATIVE) {
                Signs::NONNEG
            } else if a.has(A::NONPOSITIVE) {
                Signs::NONPOS
            } else if a.has(A::NONZERO) {
                Signs::NONZERO
            } else {
                Signs::ANY
            })
        }
        Kind::Add(c, terms) => {
            let mut s = number_sign(c);
            for (t, k) in terms {
                let ts = real_signs(t)?.mul(number_sign(k));
                s = s.add(ts);
                if s == Signs::ANY {
                    // keep checking realness of the remaining terms
                    for (t2, _) in terms {
                        real_signs(t2)?;
                    }
                    return numeric_refine(e, Signs::ANY);
                }
            }
            Some(s)
        }
        Kind::Mul(c, facs) => {
            let mut s = number_sign(c);
            for (b, x) in facs {
                s = s.mul(pow_signs(b, x)?);
            }
            Some(s)
        }
        Kind::Pow(b, x) => pow_signs(b, x),
        Kind::Func(f, args) => func_signs(f, args),
        _ => None,
    }
}

/// Refine an unknown sign of a numeric expression by evaluating it.
fn numeric_refine(e: &Expr, s: Signs) -> Option<Signs> {
    if s != Signs::ANY || !crate::functions::is_numeric_value(e) {
        return Some(s);
    }
    match crate::evalf::eval_f64(e) {
        Some(v) if v > 1e-12 * (1.0 + v.abs()) => Some(Signs::POS),
        Some(v) if v < -1e-12 * (1.0 + v.abs()) => Some(Signs::NEG),
        _ => Some(s),
    }
}

fn pow_signs(b: &Expr, x: &Expr) -> Option<Signs> {
    let bs = real_signs(b);
    if let Some(k) = x.as_i64() {
        let bs = bs?;
        if k == 0 {
            return Some(Signs::POS);
        }
        if k < 0 && bs.has(Signs::ZERO) {
            // possibly division by zero
            return None;
        }
        return Some(if k % 2 == 0 {
            if bs.has(Signs::ZERO) {
                if bs == Signs::ZERO {
                    Signs::ZERO
                } else {
                    Signs::NONNEG
                }
            } else {
                Signs::POS
            }
        } else {
            bs
        });
    }
    // non-integer exponent: need a nonnegative base and a real exponent
    let bs = bs?;
    let xs = real_signs(x)?;
    if bs == Signs::POS {
        return Some(Signs::POS);
    }
    if bs.subset(Signs::NONNEG) {
        if xs == Signs::POS {
            return Some(bs);
        }
        if bs == Signs::ZERO {
            return None;
        }
    }
    None
}

fn func_signs(f: &Func, args: &[Expr]) -> Option<Signs> {
    let a = args.first()?;
    match f {
        Func::Exp => {
            real_signs(a)?;
            Some(Signs::POS)
        }
        Func::Abs => {
            let s = if is_nonzero(a) == Some(true) { Signs::POS } else { Signs::NONNEG };
            Some(s)
        }
        Func::Log => {
            let s = real_signs(a)?;
            if s != Signs::POS {
                return None;
            }
            if crate::functions::is_numeric_value(a) {
                return numeric_refine(&crate::functions::log(a), Signs::ANY);
            }
            Some(Signs::ANY)
        }
        Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc => {
            real_signs(a)?;
            numeric_refine(&Expr::from_kind(Kind::Func(f.clone(), args.to_vec())), Signs::ANY)
        }
        Func::Atan | Func::Sinh | Func::Tanh | Func::Asinh | Func::Erf | Func::Sign => real_signs(a),
        Func::Cosh => {
            real_signs(a)?;
            Some(Signs::POS)
        }
        Func::Floor | Func::Ceiling => real_signs(a).map(|_| Signs::ANY),
        Func::Re | Func::Im => Some(Signs::ANY),
        Func::Gamma => {
            if real_signs(a)? == Signs::POS {
                Some(Signs::POS)
            } else {
                None
            }
        }
        Func::Max | Func::Min => {
            let mut s: Option<Signs> = None;
            for x in args {
                let xs = real_signs(x)?;
                s = Some(match s {
                    None => xs,
                    Some(p) => {
                        if matches!(f, Func::Max) {
                            if xs == Signs::POS || p == Signs::POS {
                                Signs::POS
                            } else {
                                p.union(xs)
                            }
                        } else if xs == Signs::NEG || p == Signs::NEG {
                            Signs::NEG
                        } else {
                            p.union(xs)
                        }
                    }
                });
            }
            s
        }
        _ => None,
    }
}

// ----------------------------------------------------------------------
// public queries
// ----------------------------------------------------------------------

pub fn is_real(e: &Expr) -> Option<bool> {
    if real_signs(e).is_some() {
        return Some(true);
    }
    match e.kind() {
        Kind::Const(Constant::I) => Some(false),
        Kind::Const(c) if c.is_infinite() || *c == Constant::NaN => Some(false),
        Kind::Sym(s) if s.assumptions().has(A::NOT_REAL) => Some(false),
        Kind::Mul(c, facs) if facs.len() == 1 && facs[0].0.is_constant(Constant::I) && facs[0].1.is_one() && !c.is_zero() => {
            Some(false)
        }
        _ => {
            if crate::functions::is_numeric_value(e) {
                if let Some(z) = crate::evalf::eval_c64(e) {
                    if z.im.abs() > 1e-10 * (1.0 + z.re.abs()) {
                        return Some(false);
                    }
                }
            }
            None
        }
    }
}

pub fn is_positive(e: &Expr) -> Option<bool> {
    if e.is_constant(Constant::Infinity) {
        return Some(true);
    }
    let s = real_signs(e)?;
    if s == Signs::POS {
        Some(true)
    } else if !s.has(Signs::POS) {
        Some(false)
    } else {
        None
    }
}

pub fn is_negative(e: &Expr) -> Option<bool> {
    if e.is_constant(Constant::NegInfinity) {
        return Some(true);
    }
    let s = real_signs(e)?;
    if s == Signs::NEG {
        Some(true)
    } else if !s.has(Signs::NEG) {
        Some(false)
    } else {
        None
    }
}

pub fn is_nonnegative(e: &Expr) -> Option<bool> {
    let s = real_signs(e)?;
    if s.subset(Signs::NONNEG) {
        Some(true)
    } else if s == Signs::NEG {
        Some(false)
    } else {
        None
    }
}

pub fn is_nonpositive(e: &Expr) -> Option<bool> {
    let s = real_signs(e)?;
    if s.subset(Signs::NONPOS) {
        Some(true)
    } else if s == Signs::POS {
        Some(false)
    } else {
        None
    }
}

pub fn is_zero(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => return Some(n.is_zero()),
        Kind::Sym(s) if s.assumptions().has(A::NONZERO) => return Some(false),
        Kind::Sym(s) if s.assumptions().has(A::ZERO) => return Some(true),
        _ => {}
    }
    if let Some(s) = real_signs(e) {
        if s == Signs::ZERO {
            return Some(true);
        }
        if !s.has(Signs::ZERO) {
            return Some(false);
        }
    }
    if crate::functions::is_numeric_value(e) {
        if let Some(z) = crate::evalf::eval_c64(e) {
            if z.norm() > 1e-10 {
                return Some(false);
            }
        }
    }
    None
}

pub fn is_nonzero(e: &Expr) -> Option<bool> {
    is_zero(e).map(|z| !z)
}

pub fn is_integer(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(n.is_integer()),
        Kind::Sym(s) => {
            let a = s.assumptions();
            if a.has(A::INTEGER) {
                Some(true)
            } else if a.has(A::NOT_INTEGER) || a.has(A::IRRATIONAL) || a.has(A::NOT_REAL) {
                Some(false)
            } else {
                None
            }
        }
        Kind::Const(_) => Some(false),
        Kind::Add(c, terms) => {
            if !c.is_integer() {
                // integer + non-integer rational is not an integer
                if terms.iter().all(|(t, k)| k.is_integer() && is_integer(t) == Some(true)) {
                    return Some(false);
                }
                return None;
            }
            for (t, k) in terms {
                if !(k.is_integer() && is_integer(t) == Some(true)) {
                    return None;
                }
            }
            Some(true)
        }
        Kind::Mul(c, facs) => {
            if !c.is_integer() {
                return None;
            }
            for (b, x) in facs {
                let ok = is_integer(b) == Some(true) && x.as_i64().is_some_and(|k| k >= 0);
                if !ok {
                    return None;
                }
            }
            Some(true)
        }
        Kind::Pow(b, x) => {
            if is_integer(b) == Some(true) && x.as_i64().is_some_and(|k| k >= 0) {
                Some(true)
            } else {
                None
            }
        }
        Kind::Func(Func::Floor | Func::Ceiling, args) => {
            if is_real(&args[0]) == Some(true) {
                Some(true)
            } else {
                None
            }
        }
        Kind::Func(Func::Sign, args) if is_real(&args[0]) == Some(true) => Some(true),
        Kind::Func(Func::Abs, args) if is_integer(&args[0]) == Some(true) => Some(true),
        _ => None,
    }
}

pub fn is_rational(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(n.is_rational()),
        Kind::Sym(s) => {
            let a = s.assumptions();
            if a.has(A::RATIONAL) {
                Some(true)
            } else if a.has(A::IRRATIONAL) {
                Some(false)
            } else {
                None
            }
        }
        Kind::Const(Constant::Pi | Constant::E) => Some(false),
        _ => is_integer(e).filter(|b| *b),
    }
}

/// `Some(true)` if even, `Some(false)` if odd, `None` if unknown / not an
/// integer.
pub fn parity_even(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => {
            let b = n.to_bigint()?;
            Some(num_integer::Integer::is_even(&b))
        }
        Kind::Sym(s) => {
            let a = s.assumptions();
            if a.has(A::EVEN) {
                Some(true)
            } else if a.has(A::ODD) {
                Some(false)
            } else {
                None
            }
        }
        Kind::Mul(c, facs) => {
            if is_integer(e) != Some(true) {
                return None;
            }
            if c.to_bigint().is_some_and(|b| num_integer::Integer::is_even(&b)) {
                return Some(true);
            }
            let mut all_odd = c.to_bigint().is_some_and(|b| num_integer::Integer::is_odd(&b));
            for (b, x) in facs {
                match parity_even(b) {
                    Some(true) if x.as_i64().is_some_and(|k| k > 0) => return Some(true),
                    Some(false) => {}
                    _ => all_odd = false,
                }
            }
            if all_odd {
                Some(false)
            } else {
                None
            }
        }
        Kind::Add(c, terms) => {
            let mut even = parity_even(&Expr::num(c.clone()))?;
            for (t, k) in terms {
                let te = parity_even(&crate::ops::add_term(t, k))?;
                even = even == te;
            }
            Some(even)
        }
        Kind::Pow(b, x) => {
            let k = x.as_i64()?;
            if k <= 0 {
                return None;
            }
            parity_even(b)
        }
        _ => None,
    }
}

pub fn is_even(e: &Expr) -> Option<bool> {
    parity_even(e)
}

pub fn is_odd(e: &Expr) -> Option<bool> {
    parity_even(e).map(|b| !b)
}

pub fn is_finite(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(_) => Some(true),
        Kind::Const(c) => Some(!c.is_infinite() && *c != Constant::NaN),
        Kind::Sym(s) => {
            if s.assumptions().has(A::FINITE) {
                Some(true)
            } else {
                None
            }
        }
        _ => {
            if real_signs(e).is_some() {
                Some(true)
            } else {
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use crate::symbol::Assumptions;

    #[test]
    fn inference() {
        let p = Expr::symbol_with("p", Assumptions::POSITIVE);
        let n = Expr::symbol_with("n", Assumptions::NEGATIVE);
        let r = Expr::symbol_with("r", Assumptions::REAL);
        let k = Expr::symbol_with("k", Assumptions::INTEGER);
        let x = Expr::symbol("x");
        assert_eq!(is_positive(&(&p + 1)), Some(true));
        assert_eq!(is_positive(&(&p * &n)), Some(false));
        assert_eq!(is_negative(&(&p * &n)), Some(true));
        assert_eq!(is_nonnegative(&r.powi(2)), Some(true));
        assert_eq!(is_positive(&r.powi(2)), None);
        assert_eq!(is_positive(&(r.powi(2) + 1)), Some(true));
        assert_eq!(is_positive(&exp(&r)), Some(true));
        assert_eq!(is_real(&x), None);
        assert_eq!(is_real(&(&r + &p)), Some(true));
        assert_eq!(is_integer(&(2 * &k + 1)), Some(true));
        assert_eq!(is_odd(&(2 * &k + 1)), Some(true));
        assert_eq!(is_even(&(2 * &k)), Some(true));
        assert_eq!(is_positive(&sqrt(&p)), Some(true));
        assert_eq!(is_positive(&sub(&Expr::pi(), &Expr::int(3))), Some(true));
        assert_eq!(is_real(&Expr::i()), Some(false));
    }
}
