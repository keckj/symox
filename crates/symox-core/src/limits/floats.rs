//! `Limit.doit` replaces the floats of the expression by exact numbers
//! before computing (`nsimplify(e)`, which for an expression with free
//! symbols is `_real_to_rational(e, rational_conversion='base10')`): every
//! double `f` becomes `nsimplify(f, rational=False)` (`mpmath.identify` at
//! tolerance `1e-15`, whose first candidate is a rational `p/q` with
//! `|p|, |q| <= 1000`), or else the rational of its 15 significant digits.
//!
//! Only the decisions that are certain are reproduced: a double whose
//! residual `|q*f - p|` is near the tolerance, or with no small rational and
//! more than 12 significant digits (where `identify` may find an algebraic
//! or transcendental closed form), or any arbitrary precision float, leaves
//! the whole expression unchanged (the limit is then computed with the
//! floats, as before).

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, ToPrimitive, Zero};

use crate::expr::{Expr, Kind};
use crate::number::Number;
use crate::ops::{add, mul};
use crate::traverse::map_children;

/// `identify`'s `maxcoeff`.
const MAXCOEFF: i64 = 1000;
/// `nsimplify`'s tolerance for 53-bit floats.
const TOL: f64 = 1e-15;

/// The exact replacement of one double, `None` when undecided.
fn rationalize(x: f64) -> Option<Number> {
    if !x.is_finite() {
        return None;
    }
    if x == 0.0 {
        return Some(Number::Int(0));
    }
    let a = x.abs();
    if !(1e-12..1e15).contains(&a) {
        return None;
    }
    let exact = BigRational::from_float(x)?;
    // pslq([x, 1], tol, maxcoeff): a relation q*x - p with the residual
    // below tol*|(x, 1)|; the smallest residual for q <= maxcoeff is at a
    // convergent of the continued fraction of x.
    let scale = TOL * (x * x + 1.0).sqrt();
    for (p, q) in convergents(&exact) {
        if q > BigInt::from(MAXCOEFF) {
            break;
        }
        if p.abs() > BigInt::from(MAXCOEFF) {
            continue;
        }
        let residual = (BigRational::from_integer(q.clone()) * &exact - BigRational::from_integer(p.clone())).abs();
        let r = residual.to_f64()? / scale;
        if r < 0.5 {
            return Number::rational(p, q);
        }
        if r <= 2.0 {
            return None;
        }
    }
    // no small rational: the 15 significant digits, decided only for a
    // short decimal (its shortest representation has at most 12 digits)
    let shortest = format!("{:e}", x);
    let (mant, exp) = shortest.split_once('e')?;
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let digits = digits.trim_start_matches('0').trim_end_matches('0');
    if digits.len() > 12 {
        return None;
    }
    decimal_rational(mant, exp.parse().ok()?)
}

/// `m * 10**e` for a decimal mantissa `m`.
fn decimal_rational(mant: &str, e: i32) -> Option<Number> {
    let neg = mant.starts_with('-');
    let m = mant.trim_start_matches('-');
    let (ip, fp) = m.split_once('.').unwrap_or((m, ""));
    let n: BigInt = format!("{ip}{fp}").parse().ok()?;
    let e = e - fp.len() as i32;
    let ten = BigInt::from(10);
    let r = if e >= 0 {
        BigRational::from_integer(n * num_traits::pow(ten, e as usize))
    } else {
        BigRational::new(n, num_traits::pow(ten, (-e) as usize))
    };
    Some(Number::from_ratio(if neg { -r } else { r }))
}

/// The convergents `p/q` (q > 0) of a rational.
fn convergents(r: &BigRational) -> Vec<(BigInt, BigInt)> {
    let (mut num, mut den) = (r.numer().clone(), r.denom().clone());
    let (mut h0, mut h1) = (BigInt::from(0), BigInt::from(1));
    let (mut k0, mut k1) = (BigInt::from(1), BigInt::from(0));
    let mut out = Vec::new();
    while !den.is_zero() && out.len() < 64 {
        let a = num_integer::Integer::div_floor(&num, &den);
        let h2 = &a * &h1 + &h0;
        let k2 = &a * &k1 + &k0;
        out.push((h2.clone(), k2.clone()));
        if k2 > BigInt::from(MAXCOEFF) {
            break;
        }
        let rem = &num - &a * &den;
        num = den;
        den = rem;
        (h0, h1, k0, k1) = (h1, h2, k1, k2);
    }
    out
}

enum Scan {
    None,
    Floats,
    Undecided,
}

fn scan(e: &Expr) -> Scan {
    let mut found = false;
    let mut stack = vec![e.clone()];
    while let Some(n) = stack.pop() {
        if n.float_info().max_precision.is_none() {
            continue;
        }
        match n.kind() {
            Kind::Num(Number::Float(f)) => {
                if rationalize(*f).is_none() {
                    return Scan::Undecided;
                }
                found = true;
            }
            Kind::Num(Number::MpFloat(..)) => return Scan::Undecided,
            Kind::Num(_) => {}
            _ => stack.extend(n.args()),
        }
    }
    if found {
        Scan::Floats
    } else {
        Scan::None
    }
}

fn replace(e: &Expr) -> Expr {
    if e.float_info().max_precision.is_none() {
        return e.clone();
    }
    match e.kind() {
        Kind::Num(Number::Float(f)) => rationalize(*f).map(Expr::num).unwrap_or_else(|| e.clone()),
        Kind::Num(_) => e.clone(),
        // the numeric coefficients are arguments, not children
        Kind::Add(..) => add(&e.args().iter().map(replace).collect::<Vec<_>>()),
        Kind::Mul(..) => mul(&e.args().iter().map(replace).collect::<Vec<_>>()),
        _ => map_children(e, &mut |c| replace(c)),
    }
}

/// `e` with its floats replaced as `Limit.doit` does, or `None` when `e`
/// has no float or a replacement is undecided.
pub fn nsimplify_floats(e: &Expr) -> Option<Expr> {
    match scan(e) {
        Scan::Floats => Some(replace(e)),
        Scan::None | Scan::Undecided => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(x: f64) -> Option<String> {
        rationalize(x).map(|n| Expr::num(n).to_string())
    }

    #[test]
    fn identify_rationals_and_short_decimals() {
        assert_eq!(q(0.5).as_deref(), Some("1/2"));
        assert_eq!(q(-0.5).as_deref(), Some("-1/2"));
        assert_eq!(q(0.3).as_deref(), Some("3/10"));
        assert_eq!(q(0.03).as_deref(), Some("3/100"));
        assert_eq!(q(4.0).as_deref(), Some("4"));
        assert_eq!(q(1.0 / 3.0).as_deref(), Some("1/3"));
        assert_eq!(q(1234.5).as_deref(), Some("2469/2"));
        assert_eq!(q(0.0001234).as_deref(), Some("617/5000000"));
        // no small rational and 16 digits: identify may find a closed form
        assert_eq!(q(std::f64::consts::E), None);
        assert_eq!(q(std::f64::consts::SQRT_2), None);
        // residual at the tolerance: undecided
        assert_eq!(q(0.333333333333333), None);
    }
}
