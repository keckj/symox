//! Limits with SymPy semantics.
//!
//! ```
//! use symox::prelude::*;
//! use symox::limits::{limit, Dir};
//!
//! let x = Expr::symbol("x");
//! let e = div(&sin(&x), &x);
//! assert_eq!(limit(&e, &x, &Expr::zero(), Dir::Plus).unwrap(), Expr::one());
//! ```
//!
//! Strategy: a leading term computed by the series engine answers most
//! limits at finite points (and rational functions at infinity) directly;
//! everything else is mapped to `x -> oo` and handed to the Gruntz
//! comparability class (mrv) algorithm in [`gruntz`]. Expressions with
//! bounded oscillating factors (`sin(x)/x` at `oo`) are handled by treating
//! the oscillating parts as bounded parameters.

pub mod gruntz;

use rustc_hash::FxHashMap;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::ops::{add, mul2, neg, pow};
use crate::series::engine::{degree_in, depends};
use crate::series::{leading_term_with, LeadTerm};
use crate::subs::xreplace;

/// Direction of a limit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dir {
    Plus,
    Minus,
    Both,
}

impl Dir {
    pub fn parse(s: &str) -> Option<Dir> {
        match s {
            "+" => Some(Dir::Plus),
            "-" => Some(Dir::Minus),
            "+-" | "-+" => Some(Dir::Both),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum LimitError {
    /// The algorithm could not decide (return an unevaluated `Limit`).
    Unknown(String),
    /// Two-sided limit with different one-sided values.
    Differ(Expr, Expr),
}

impl std::fmt::Display for LimitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LimitError::Unknown(s) => write!(f, "cannot compute the limit: {s}"),
            LimitError::Differ(l, r) => write!(
                f,
                "The limit does not exist since left hand limit = {l} and right hand limit = {r}"
            ),
        }
    }
}

pub type LRes<T> = Result<T, LimitError>;

pub(crate) fn unknown<T>(s: impl Into<String>) -> LRes<T> {
    Err(LimitError::Unknown(s.into()))
}

fn is_infinite(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(c) if c.is_infinite())
        || matches!(e.kind(), Kind::Mul(_, f) if f.iter().any(|(b, _)| b.is_constant(Constant::Infinity)))
}

/// `limit(e, x, x0, dir)`.
pub fn limit(e: &Expr, x: &Expr, x0: &Expr, dir: Dir) -> LRes<Expr> {
    if !depends(e, x) {
        return Ok(e.clone());
    }
    let infinite_point = x0.is_constant(Constant::Infinity) || x0.is_constant(Constant::NegInfinity);
    match dir {
        _ if infinite_point => limit_one_sided(e, x, x0, false),
        Dir::Plus => limit_one_sided(e, x, x0, false),
        Dir::Minus => limit_one_sided(e, x, x0, true),
        Dir::Both => {
            let r = limit_one_sided(e, x, x0, false)?;
            let l = limit_one_sided(e, x, x0, true)?;
            if r == l {
                return Ok(r);
            }
            if is_infinite(&r) && is_infinite(&l) {
                return Ok(Expr::zoo());
            }
            Err(LimitError::Differ(l, r))
        }
    }
}

fn limit_one_sided(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> LRes<Expr> {
    if x0.is_constant(Constant::ComplexInfinity) || x0.is_constant(Constant::NaN) {
        return unknown("bad limit point");
    }
    if let Some(r) = fast_path(e, x, x0, minus) {
        return Ok(r);
    }
    // map to x -> oo
    let t = Expr::symbol("_symox_lim_t");
    let inv = pow(&t, &Expr::minus_one());
    let v = if x0.is_constant(Constant::Infinity) {
        t.clone()
    } else if x0.is_constant(Constant::NegInfinity) {
        neg(&t)
    } else if minus {
        add(&[x0.clone(), neg(&inv)])
    } else {
        add(&[x0.clone(), inv])
    };
    let mut m = FxHashMap::default();
    m.insert(x.clone(), v);
    let et = xreplace(e, &m);
    let mut g = gruntz::Gruntz::new(&t);
    match g.limitinf(&et) {
        Ok(r) => Ok(r),
        Err(err) => match bounded_heuristic(&et, &t) {
            Some(r) => Ok(r),
            None => Err(err),
        },
    }
}

/// Value of the limit from a leading term `c * t**e` (with `t -> 0+`
/// and `log(t)` written as `lt.logt` in `lt.coeff_logt`).
pub(crate) fn from_leading(lt: &LeadTerm) -> Option<Expr> {
    if lt.exp.is_positive() {
        return Some(Expr::zero());
    }
    let c = &lt.coeff_logt;
    let k = degree_in(c, &lt.logt)?;
    if lt.exp.is_zero() && k == 0 {
        return Some(c.clone());
    }
    // infinite: sign of the dominant coefficient (log(t) -> -oo)
    let ck = if k == 0 { c.clone() } else { coeff_of_power(c, &lt.logt, k)? };
    let parity = if k % 2 == 0 { 1 } else { -1 };
    match crate::series::engine::sign_of(&crate::series::coeff::C::from_expr(ck.clone())) {
        Some(s) => Some(if s * parity > 0 { Expr::oo() } else { Expr::neg_oo() }),
        None => {
            if crate::functions::is_numeric_value(&ck) {
                return None; // complex
            }
            if ck.free_symbols_empty() {
                return None;
            }
            let s = crate::functions::sign(&ck);
            Some(mul2(&mul2(&Expr::int(parity as i64), &s), &Expr::oo()))
        }
    }
}

trait FreeEmpty {
    fn free_symbols_empty(&self) -> bool;
}

impl FreeEmpty for Expr {
    fn free_symbols_empty(&self) -> bool {
        crate::traverse::free_symbols(self).is_empty()
    }
}

/// Coefficient of `s**k` in a polynomial in `s`.
pub(crate) fn coeff_of_power(e: &Expr, s: &Expr, k: i64) -> Option<Expr> {
    let terms: Vec<Expr> = if e.is_add() { e.args() } else { vec![e.clone()] };
    let mut acc = Vec::new();
    for t in terms {
        if degree_in(&t, s)? == k {
            let sk = pow(s, &Expr::int(k));
            acc.push(crate::ops::div(&t, &sk));
        }
    }
    Some(add(&acc))
}

fn fast_path(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> Option<Expr> {
    if has_bad(e) {
        return None;
    }
    let logt = Expr::symbol("_symox_lim_logt");
    let lt = leading_term_with(e, x, x0, minus, &logt, &logt).ok()?;
    if depends(&lt.coeff_logt, x) {
        return None;
    }
    from_leading(&lt)
}

fn has_bad(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(c) => c.is_infinite() || *c == Constant::NaN,
        Kind::Func(Func::Order, _) => true,
        Kind::Num(_) | Kind::Sym(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(has_bad),
    }
}

/// Replace oscillating `sin(u)`/`cos(u)` (with `u -> ±oo`) by bounded
/// parameters; if the limit exists and does not depend on them, it is
/// the limit of the original expression.
fn bounded_heuristic(e: &Expr, t: &Expr) -> Option<Expr> {
    let mut map: FxHashMap<Expr, Expr> = FxHashMap::default();
    let mut g = gruntz::Gruntz::new(t);
    collect_oscillating(e, t, &mut g, &mut map);
    if map.is_empty() {
        return None;
    }
    let params: Vec<Expr> = map.values().cloned().collect();
    let e2 = xreplace(e, &map);
    let mut g2 = gruntz::Gruntz::new(t);
    let r = g2.limitinf(&e2).ok()?;
    if params.iter().any(|p| depends(&r, p)) {
        return None;
    }
    Some(r)
}

fn collect_oscillating(e: &Expr, t: &Expr, g: &mut gruntz::Gruntz, map: &mut FxHashMap<Expr, Expr>) {
    if !depends(e, t) || map.contains_key(e) {
        return;
    }
    if let Kind::Func(Func::Sin | Func::Cos, a) = e.kind() {
        if let Ok(l) = g.limitinf(&a[0]) {
            if is_infinite(&l) {
                let p = Expr::symbol(&format!("_symox_bounded{}", map.len()));
                map.insert(e.clone(), p);
                return;
            }
        }
    }
    for a in e.args() {
        collect_oscillating(&a, t, g, map);
    }
}

#[cfg(test)]
mod tests;
