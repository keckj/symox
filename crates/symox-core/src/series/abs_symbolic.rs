//! `Abs(g)` and `sign(g)` when the sign of the leading coefficient `c` of
//! `g` is not a decided number (`Abs(x + y)`, `Abs(x + p)` with `p > 0`).
//!
//! The reference `Abs._eval_nseries` is `(sign(c)*g).expand()` and
//! `sign._eval_nseries` is `sign(c)` (through `as_leading_term`), with
//! `c = g.leadterm(x)[0]`. `sign(c)` is `1` or `-1` when the assumptions
//! decide it and stays a symbolic `sign(c)` otherwise
//! (`Abs(x + y).series(x)` is `x*sign(y) + y*sign(y)`).
//!
//! A coefficient that is provably nonreal is declined: there the
//! reference evaluates `Abs` through its real and imaginary parts
//! (`Abs(x + I*p)` is a square root), which this rule does not give. A
//! coefficient holding `log(x)` is declined too (the reference substitutes
//! `logx`, the engine keeps a separate symbol).
//!
//! The reference expands on a positive dummy, where `Abs(x**k*r)` and
//! `sign(x**k*r)` with `r` free of `x` evaluate to `x**k*Abs(r)` and
//! `sign(r)` before the expansion (`Abs(x*y).series(x)` is `x*Abs(y)`);
//! `split_positive` gives that form first, from the sign of the leading
//! coefficient of `x` along the expansion parameter (`x = -t` from the
//! left gives `-x*Abs(y)`).

use super::coeff::C;
use super::engine::depends;
use crate::assumptions::{is_negative, is_positive, is_real, is_zero};
use crate::expr::{Expr, Func, Kind};
use crate::ops::{div, mul, mul2, mul_factor};

/// `sign(c)` for a leading coefficient whose numeric sign is undecided,
/// or `None` when the reference rule does not apply.
pub fn lead_sign(c: &C, x: &Expr, logt: &Expr) -> Option<Expr> {
    let e = c.to_expr();
    if depends(&e, x) || depends(&e, logt) || is_zero(&e) == Some(true) {
        return None;
    }
    if is_positive(&e) == Some(true) {
        return Some(Expr::one());
    }
    if is_negative(&e) == Some(true) {
        return Some(Expr::minus_one());
    }
    if is_real(&e) == Some(false) || crate::functions::is_numeric_value(&e) {
        return None;
    }
    Some(crate::functions::func(Func::Sign, vec![e]))
}

/// `x**k*Abs(r)` (`Abs`) or `sign(r)` (`sign`) when `arg = x**k*r` with
/// real numeric `k`, `r` free of `x` and `x > 0` along the expansion
/// (`xsign = 1`); for `x < 0` (`xsign = -1`, an expansion from the left)
/// `k` must be an integer and the result takes the factor `(-1)**k`.
/// `None` otherwise.
pub fn split_positive(f: Func, arg: &Expr, x: &Expr, xsign: Option<i32>) -> Option<Expr> {
    let Kind::Mul(_, facs) = arg.kind() else { return None };
    let xsign = xsign?;
    if xsign != 1 && xsign != -1 {
        return None;
    }
    let mut xs = Vec::new();
    let mut odd = false;
    for (b, p) in facs {
        let g = mul_factor(b, p);
        if !depends(&g, x) {
            continue;
        }
        if b != x || !crate::functions::is_numeric_value(p) || is_real(p) != Some(true) {
            return None;
        }
        if xsign == -1 {
            odd ^= p.as_i64()? % 2 != 0;
        }
        xs.push(g);
    }
    let unit = if odd { Expr::minus_one() } else { Expr::one() };
    if xs.is_empty() {
        return None;
    }
    let xpart = mul(&xs);
    let rest = div(arg, &xpart);
    if depends(&rest, x) {
        return None;
    }
    match f {
        Func::Abs => Some(mul(&[unit, xpart, crate::functions::func(Func::Abs, vec![rest])])),
        Func::Sign => Some(mul2(&unit, &crate::functions::func(Func::Sign, vec![rest]))),
        _ => None,
    }
}
