//! Certified real poles of rational definite integrals.
//!
//! An ordinary improper integral sums independent one-sided integrals.
//! Opposing infinite contributions give NaN, never a Cauchy principal value.

use std::cmp::Ordering;

use num_bigint::BigInt;

use crate::expr::{Constant, Expr, Kind};
use crate::number::Number;
use crate::ops::{mul2, sub};
use crate::poly::{dense, Int, Rat};

pub(super) enum Analysis {
    NotRational,
    Regular(Expr),
    Unknown,
    Value(Expr),
}

#[derive(Clone)]
struct Interval {
    lo: Rat,
    hi: Rat,
}

impl Interval {
    fn exact(x: Rat) -> Self {
        Self {
            lo: x.clone(),
            hi: x,
        }
    }
    fn add(&self, other: &Self) -> Self {
        Self {
            lo: self.lo.add(&other.lo),
            hi: self.hi.add(&other.hi),
        }
    }
    fn mul(&self, other: &Self) -> Self {
        let values = [
            self.lo.mul(&other.lo),
            self.lo.mul(&other.hi),
            self.hi.mul(&other.lo),
            self.hi.mul(&other.hi),
        ];
        Self {
            lo: values.iter().min().unwrap().clone(),
            hi: values.iter().max().unwrap().clone(),
        }
    }
    fn reciprocal(&self) -> Option<Self> {
        if self.lo <= Rat::zero() && self.hi >= Rat::zero() {
            return None;
        }
        Some(Self {
            lo: Rat::one().div(&self.hi),
            hi: Rat::one().div(&self.lo),
        })
    }
    fn integer_power(&self, exponent: i64) -> Option<Self> {
        if exponent.unsigned_abs() > 64 {
            return None;
        }
        let mut result = Self::exact(Rat::one());
        let mut base = self.clone();
        let mut k = exponent.unsigned_abs();
        while k != 0 {
            if k & 1 != 0 {
                result = result.mul(&base);
            }
            k >>= 1;
            if k != 0 {
                base = base.mul(&base);
            }
        }
        if exponent < 0 {
            result.reciprocal()
        } else {
            Some(result)
        }
    }
    fn positive_root(&self, n: u32, bits: usize) -> Option<Self> {
        if self.lo.is_negative() || n == 0 || n > 16 {
            return None;
        }
        let scale = BigInt::from(1) << bits;
        let bound = |x: &Rat, upper: bool| {
            let scaled = (x.n.to_big() << (bits * n as usize)) / x.d.to_big();
            let floor = scaled.nth_root(n);
            Rat::new(
                Int::from_big(if upper { floor + 1 } else { floor }),
                Int::from_big(scale.clone()),
            )
        };
        Some(Self {
            lo: bound(&self.lo, false),
            hi: bound(&self.hi, true),
        })
    }
}

fn enclose(e: &Expr, bits: usize) -> Option<Interval> {
    match e.kind() {
        Kind::Num(n) => {
            if !n.to_f64().is_finite() {
                return None;
            }
            Some(Interval::exact(crate::poly::convert::number_to_rat(n)))
        }
        Kind::Const(Constant::Pi) => Some(Interval {
            lo: Rat::new(Int::from(3141592653589793), Int::from(1000000000000000)),
            hi: Rat::new(Int::from(3141592653589794), Int::from(1000000000000000)),
        }),
        Kind::Const(Constant::E) => Some(Interval {
            lo: Rat::new(Int::from(2718281828459045), Int::from(1000000000000000)),
            hi: Rat::new(Int::from(2718281828459046), Int::from(1000000000000000)),
        }),
        Kind::Add(c, terms) => {
            let mut out = Interval::exact(crate::poly::convert::number_to_rat(c));
            for (term, coefficient) in terms {
                out = out.add(&enclose(term, bits)?.mul(&Interval::exact(
                    crate::poly::convert::number_to_rat(coefficient),
                )));
            }
            Some(out)
        }
        Kind::Mul(c, factors) => {
            let mut out = Interval::exact(crate::poly::convert::number_to_rat(c));
            for (base, exponent) in factors {
                out = out.mul(&enclose_power(base, exponent, bits)?);
            }
            Some(out)
        }
        Kind::Pow(base, exponent) => enclose_power(base, exponent, bits),
        _ => None,
    }
}

fn enclose_power(base: &Expr, exponent: &Expr, bits: usize) -> Option<Interval> {
    let value = enclose(base, bits)?;
    if let Some(k) = exponent.as_i64() {
        return value.integer_power(k);
    }
    let Kind::Num(Number::Rational(r)) = exponent.kind() else {
        return None;
    };
    use num_traits::ToPrimitive;
    let p = r.numer().to_i64()?;
    let q = r.denom().to_u32()?;
    value.positive_root(q, bits)?.integer_power(p)
}

fn certified_sign(e: &Expr) -> Option<i8> {
    let e = crate::simplify::simplify(&crate::expand::expand(e));
    if e.is_zero() {
        return Some(0);
    }
    for bits in [32, 64, 128, 256] {
        let interval = enclose(&e, bits)?;
        if interval.lo > Rat::zero() {
            return Some(1);
        }
        if interval.hi < Rat::zero() {
            return Some(-1);
        }
        if interval.lo.is_zero() && interval.hi.is_zero() {
            return Some(0);
        }
    }
    None
}

pub(super) fn compare(a: &Expr, b: &Expr) -> Option<Ordering> {
    if a == b {
        return Some(Ordering::Equal);
    }
    if a.is_constant(Constant::NegInfinity) || b.is_constant(Constant::Infinity) {
        return Some(Ordering::Less);
    }
    if a.is_constant(Constant::Infinity) || b.is_constant(Constant::NegInfinity) {
        return Some(Ordering::Greater);
    }
    certified_sign(&sub(a, b)).map(|s| s.cmp(&0))
}

fn variations(signs: impl Iterator<Item = i8>) -> usize {
    let mut previous = 0;
    let mut changes = 0;
    for sign in signs.filter(|&s| s != 0) {
        if previous != 0 && sign != previous {
            changes += 1;
        }
        previous = sign;
    }
    changes
}

/// Exact number of distinct real roots, including multiple roots.
fn real_root_count(coefficients: &[Rat]) -> usize {
    let derivative = dense::derivative(coefficients);
    let gcd = dense::gcd_field(coefficients, &derivative);
    let squarefree = dense::divrem(coefficients, &gcd).0;
    if squarefree.len() < 2 {
        return 0;
    }
    let mut sturm = vec![squarefree.clone(), dense::derivative(&squarefree)];
    loop {
        let r = dense::neg(&dense::rem(&sturm[sturm.len() - 2], sturm.last().unwrap()));
        if r.is_empty() {
            break;
        }
        sturm.push(r);
    }
    let at_infinity = |negative: bool| {
        variations(sturm.iter().map(|p| {
            let sign = if p.last().unwrap().is_negative() {
                -1
            } else {
                1
            };
            if negative && (p.len() - 1) % 2 == 1 {
                -sign
            } else {
                sign
            }
        }))
    };
    at_infinity(true) - at_infinity(false)
}

fn evaluated(coefficients: &[Expr], root: &Expr) -> Expr {
    coefficients.iter().rev().fold(Expr::zero(), |sum, c| {
        crate::ops::add2(&mul2(&sum, root), c)
    })
}

/// Analyze all rational poles before evaluating an antiderivative. Exact
/// Sturm counts prevent a partial radical-root list from hiding real poles.
pub(super) fn analyze(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Analysis {
    let cancelled = crate::poly::api::cancel(f).unwrap_or_else(|_| f.clone());
    let (num, den) = crate::simplify::numer_denom(&crate::simplify::together(&cancelled));
    let Some(nc) = super::poly_coeffs(&num, x) else {
        return Analysis::NotRational;
    };
    let Some(dc) = super::poly_coeffs(&den, x) else {
        return Analysis::NotRational;
    };
    if cancelled.is_zero() {
        return Analysis::Value(Expr::zero());
    }
    let mut positive = false;
    let mut negative = false;
    let mut contribution = |sign: i8| {
        if sign > 0 {
            positive = true;
        } else {
            negative = true;
        }
    };
    if dc.len() > 1 {
        let rational: Option<Vec<Rat>> = dc
            .iter()
            .map(|e| e.as_number().and_then(Rat::from_number))
            .collect();
        let expected_real = rational.as_ref().map(|c| real_root_count(c));
        let roots = if expected_real == Some(0) {
            Vec::new()
        } else {
            match crate::poly::roots::roots(&den, std::slice::from_ref(x)) {
                Ok(r) => r,
                Err(_) => return Analysis::Unknown,
            }
        };
        let real: Vec<_> = roots
            .iter()
            .filter(|(r, _)| enclose(r, 64).is_some())
            .collect();
        if let Some(count) = expected_real {
            if real.len() != count {
                return Analysis::Unknown;
            }
        } else if real.iter().map(|(_, m)| *m as usize).sum::<usize>() != dc.len() - 1 {
            return Analysis::Unknown;
        }
        for (root, multiplicity) in real {
            let Some(left) = compare(root, a) else {
                return Analysis::Unknown;
            };
            let Some(right) = compare(root, b) else {
                return Analysis::Unknown;
            };
            if left == Ordering::Less || right == Ordering::Greater {
                continue;
            }
            let mut derivative = dc.clone();
            for _ in 0..*multiplicity {
                derivative = derivative
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(i, c)| mul2(&Expr::int(i as i64), c))
                    .collect();
            }
            let coefficient = mul2(&evaluated(&nc, root), &evaluated(&derivative, root));
            let Some(sign) = certified_sign(&coefficient) else {
                return Analysis::Unknown;
            };
            if sign == 0 {
                return Analysis::Unknown;
            }
            if left != Ordering::Equal {
                contribution(if multiplicity % 2 == 0 { sign } else { -sign });
            }
            if right != Ordering::Equal {
                contribution(sign);
            }
        }
    }
    // Every infinite tail of a rational function is determined by its
    // leading power; include it independently of finite-pole divergences.
    let degree = nc.len() as i64 - dc.len() as i64;
    if degree >= -1 && (a.is_constant(Constant::NegInfinity) || b.is_constant(Constant::Infinity)) {
        let Some(sign) = certified_sign(&mul2(nc.last().unwrap(), dc.last().unwrap())) else {
            return Analysis::Unknown;
        };
        if sign == 0 {
            return Analysis::Unknown;
        }
        if a.is_constant(Constant::NegInfinity) {
            contribution(if degree % 2 == 0 { sign } else { -sign });
        }
        if b.is_constant(Constant::Infinity) {
            contribution(sign);
        }
    }
    if positive && negative {
        Analysis::Value(Expr::nan())
    } else if positive {
        Analysis::Value(Expr::oo())
    } else if negative {
        Analysis::Value(Expr::neg_oo())
    } else {
        Analysis::Regular(cancelled)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    #[test]
    fn exact_sturm_counts_and_improper_integrals() {
        let x = Expr::symbol("x");
        for (input, expected) in [
            ("1/x", "nan"),
            ("1/x**2", "oo"),
            ("-1/(x**2-2)**2", "-oo"),
            ("1/(x**2-2)", "nan"),
        ] {
            let value = super::super::integrate_definite(
                &parse_str(input).unwrap(),
                &x,
                &Expr::int(-2),
                &Expr::int(2),
            );
            assert_eq!(value, parse_str(expected).unwrap(), "{input}");
        }
        assert_eq!(
            real_root_count(&[
                Rat::from_i64(1),
                Rat::zero(),
                Rat::from_i64(2),
                Rat::zero(),
                Rat::one()
            ]),
            0
        );
        assert_eq!(
            real_root_count(&[Rat::from_i64(-2), Rat::zero(), Rat::one()]),
            2
        );
    }

    #[test]
    fn endpoint_reversal_and_independent_tails() {
        let x = Expr::symbol("x");
        for (f, a, b, want) in [
            ("1/x", "-1", "0", "-oo"),
            ("1/x", "1", "0", "-oo"),
            ("1/x**2", "1", "0", "-oo"),
            ("1/(x-sqrt(2))", "0", "sqrt(2)", "-oo"),
            ("1/(x-sqrt(2))", "sqrt(2)", "2", "oo"),
            ("1/x**2-1", "0", "oo", "nan"),
            ("x", "-oo", "oo", "nan"),
            ("1/(x**2+1)", "-oo", "oo", "pi"),
        ] {
            let value = super::super::integrate_definite(
                &parse_str(f).unwrap(),
                &x,
                &parse_str(a).unwrap(),
                &parse_str(b).unwrap(),
            );
            assert_eq!(value, parse_str(want).unwrap(), "{f},({a},{b})");
        }
        let root = parse_str("sqrt(2)").unwrap();
        let upper = parse_str("sqrt(2)+1/1267650600228229401496703205376").unwrap();
        assert_eq!(compare(&root, &upper), Some(Ordering::Less));
        let lower =
            parse_str("1414213562373095048801688724209/1000000000000000000000000000000").unwrap();
        let upper =
            parse_str("1414213562373095048801688724210/1000000000000000000000000000000").unwrap();
        assert_eq!(compare(&lower, &root), Some(Ordering::Less));
        assert_eq!(compare(&root, &upper), Some(Ordering::Less));
    }
}
