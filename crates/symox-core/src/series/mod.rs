//! Truncated power series (Laurent/Puiseux, with `log` terms) and the
//! `Order` term.
//!
//! ```
//! use symox::prelude::*;
//! use symox::series::series;
//!
//! let x = Expr::symbol("x");
//! let s = series(&sin(&x), &x, &Expr::zero(), 6, false).unwrap();
//! assert_eq!(s.to_string(), "x - x**3/6 + x**5/120 + O(x**6)");
//! ```
//!
//! Expansions are computed by [`engine::Engine`] with `O(n²)` recurrences
//! (`exp`, `log`, powers, `sin`/`cos`, `tan`, inverse functions through
//! `F(a) = F(a0) + ∫ F'(a) a'`), never by repeated differentiation.

pub mod asym;
pub mod coeff;
pub mod engine;
pub mod order;
pub mod print;
pub mod ser;

use rustc_hash::FxHashMap;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, mul, mul2, mul_factor, neg, pow, sub};
use crate::subs::xreplace;
use coeff::C;
use engine::{degree_in, depends, var_series, Engine};
pub use order::{contains, get_o, getn, is_order, order, remove_o};
pub use ser::{SErr, SRes, Ser};

/// Maximal ramification index tried for Puiseux expansions.
const MAX_D: i64 = 1 << 12;

/// The symbol standing for `log(t^D)` inside the engine.
pub fn logt_symbol() -> Expr {
    Expr::symbol("_symox_logt")
}

/// Run `f` with an engine for the expansion point, increasing the
/// ramification when needed. Returns the result and `D`.
pub fn with_engine<T>(
    x: &Expr,
    x0: &Expr,
    minus: bool,
    logt: &Expr,
    mut f: impl FnMut(&mut Engine) -> SRes<T>,
) -> SRes<(T, i64)> {
    let mut d = 1i64;
    loop {
        let mut eng = Engine::new(x, var_series(x0, minus, d), d, logt);
        match f(&mut eng) {
            Err(SErr::NeedD(k)) if k > 1 && d * k <= MAX_D => d *= k,
            Err(SErr::NeedD(_)) => return Err(SErr::Unsupported("ramification too large".into())),
            Err(e) => return Err(e),
            Ok(v) => return Ok((v, d)),
        }
    }
}

/// `x - x0` (or the small variable at infinity) as an expression.
fn small_var(x: &Expr, x0: &Expr, minus: bool) -> Expr {
    if x0.is_constant(Constant::Infinity) {
        return pow(x, &Expr::minus_one());
    }
    if x0.is_constant(Constant::NegInfinity) {
        return neg(&pow(x, &Expr::minus_one()));
    }
    if minus {
        sub(x0, x)
    } else {
        sub(x, x0)
    }
}

/// Leading term of an expansion.
#[derive(Clone, Debug)]
pub struct LeadTerm {
    /// coefficient with `log` of the small variable written out
    pub coeff: Expr,
    /// coefficient in terms of [`LeadTerm::logt`]
    pub coeff_logt: Expr,
    pub logt: Expr,
    /// exponent of the small variable (`x - x0`, or `1/x` at infinity)
    pub exp: Number,
}

/// Leading term `c * (x - x0)**e` of `expr` at `x0` (from the right unless
/// `minus`).
pub fn leading_term(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> SRes<LeadTerm> {
    let logt = logt_symbol();
    leading_term_with(e, x, x0, minus, &logt, &crate::functions::log(&small_var(x, x0, minus)))
}

/// Whether `e` is exactly zero on the one-sided neighbourhood of `x0`
/// (`Heaviside(-x)` or `floor(x)` for `x -> 0+`), where the leading term
/// does not exist but the expansion and the limit are `0`.
pub fn vanishes_near(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> bool {
    let logt = logt_symbol();
    with_engine(x, x0, minus, &logt, |eng| Ok(eng.ser(e, 1)?.is_exact_zero())).is_ok_and(|(zero, _)| zero)
}

/// [`leading_term`] with a caller supplied `log` symbol and its value.
pub fn leading_term_with(e: &Expr, x: &Expr, x0: &Expr, minus: bool, logt: &Expr, logval: &Expr) -> SRes<LeadTerm> {
    let ((c, v), d) = with_engine(x, x0, minus, logt, |eng| {
        eng.logval = Some(logval.clone());
        eng.leading(e)
    })?;
    let ce = c.to_expr();
    let coeff = if depends(&ce, logt) { xreplace(&ce, &single(logt, logval)) } else { ce.clone() };
    Ok(LeadTerm { coeff, coeff_logt: ce, logt: logt.clone(), exp: Number::rational_i64(v, d).unwrap() })
}

fn single(a: &Expr, b: &Expr) -> FxHashMap<Expr, Expr> {
    let mut m = FxHashMap::default();
    m.insert(a.clone(), b.clone());
    m
}

/// Split `log(x)**k` factors off a product (SymPy keeps them outside the
/// expansion: `series(log(x)*sin(x))` has an `O(x**6*log(x))` term).
fn split_log_factor(e: &Expr, x: &Expr) -> Option<(Expr, Expr, i64)> {
    let Kind::Mul(c, facs) = e.kind() else { return None };
    let mut logs = Vec::new();
    let mut rest = vec![Expr::num(c.clone())];
    let mut k = 0i64;
    for (b, p) in facs {
        let is_logx = matches!(b.kind(), Kind::Func(Func::Log, a) if a[0] == *x);
        match (is_logx, p.as_i64()) {
            (true, Some(m)) if m > 0 => {
                k += m;
                logs.push(mul_factor(b, p));
            }
            _ => rest.push(mul_factor(b, p)),
        }
    }
    if logs.is_empty() {
        return None;
    }
    Some((mul(&logs), mul(&rest), k))
}

/// `series(expr, x, x0, n, dir)` with SymPy semantics: the expansion of
/// `expr` around `x = x0` up to (excluding) `(x - x0)**n`, followed by an
/// order term, or the exact result when the expansion terminates.
/// `minus` selects `dir='-'`.
pub fn series(e: &Expr, x: &Expr, x0: &Expr, n: i64, minus: bool) -> SRes<Expr> {
    if !depends(e, x) {
        return Ok(e.clone());
    }
    if e.may_have_order() {
        if let Some(o) = get_o(e) {
            let rest = remove_o(e);
            if !rest.may_have_order() {
                let s = series(&rest, x, x0, n, minus)?;
                return Ok(add(&[s, o]));
            }
        }
        if contains_order(e) {
            return Err(SErr::Unsupported("series of an expression containing O".into()));
        }
    }
    if x0.is_zero() && !minus {
        if let Some((l, rest, k)) = split_log_factor(e, x) {
            return series_impl(&rest, x, x0, n, minus, Some((l, k)));
        }
    }
    match series_impl(e, x, x0, n, minus, None) {
        // exp(1/x) and friends are returned unevaluated, like SymPy
        Err(SErr::Pole(_)) if e.is_func(&Func::Exp) => Ok(e.clone()),
        r => r,
    }
}

fn series_impl(e: &Expr, x: &Expr, x0: &Expr, n: i64, minus: bool, logfac: Option<(Expr, i64)>) -> SRes<Expr> {
    let logt = logt_symbol();
    let ((s, extra_logdeg), d) = with_engine(x, x0, minus, &logt, |eng| {
        let target = n.saturating_mul(eng.d);
        let mut m = target;
        let mut s = eng.ser(e, m)?;
        for _ in 0..6 {
            if s.exact || s.prec >= target {
                break;
            }
            m += (target - s.prec).max(1);
            s = eng.ser(e, m)?;
        }
        let s = s.trunc(target.min(s.prec.max(target)));
        // log degree of the order term
        let mut logdeg = 0i64;
        if !s.exact && s.c.iter().any(|c| matches!(c, C::E(e) if depends(e, &logt))) {
            let s2 = eng.ser(e, target + eng.d)?;
            let mut found = false;
            for k in target..target + eng.d {
                if let C::E(c) = s2.at(k) {
                    if let Some(dg) = degree_in(&c, &logt) {
                        logdeg = logdeg.max(dg);
                        found = true;
                    }
                } else if !s2.at(k).is_zero() {
                    found = true;
                }
            }
            if !found {
                for c in &s.c {
                    if let C::E(c) = c {
                        logdeg = logdeg.max(degree_in(c, &logt).unwrap_or(0));
                    }
                }
            }
        }
        Ok((s, logdeg))
    })?;
    let base = small_var(x, x0, minus);
    let logval = crate::functions::log(&base);
    let map = single(&logt, &logval);
    let (lfac, lk) = match &logfac {
        Some((l, k)) => (Some(l.clone()), *k),
        None => (None, 0),
    };
    let mut terms: Vec<Expr> = Vec::with_capacity(s.c.len() + 1);
    for (i, c) in s.c.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let k = s.val + i as i64;
        let q = Number::rational_i64(k, d).unwrap();
        let mut ce = c.to_expr();
        if depends(&ce, &logt) {
            ce = xreplace(&ce, &map);
        }
        let mut t = mul2(&ce, &pow(&base, &Expr::num(q)));
        if let Some(l) = &lfac {
            t = mul2(&t, l);
        }
        terms.push(t);
    }
    if !s.exact {
        let q = Number::rational_i64(s.prec, d).unwrap();
        // order term is always in terms of (x - x0), whatever the direction
        let obase = small_var(x, x0, false);
        let mut oe = pow(&obase, &Expr::num(q));
        let ld = extra_logdeg + lk;
        if ld > 0 {
            oe = mul2(&oe, &pow(&crate::functions::log(&obase), &Expr::int(ld)));
        }
        terms.push(order(&oe, &[(x.clone(), x0.clone())]));
    }
    Ok(add(&terms))
}

fn contains_order(e: &Expr) -> bool {
    e.may_have_order() && (is_order(e) || e.args().iter().any(contains_order))
}

/// SymPy's `taylor_term(n, x)`: the `x**n` term of the expansion at 0.
pub fn taylor_term(e: &Expr, x: &Expr, n: i64) -> SRes<Expr> {
    let logt = logt_symbol();
    let (c, _) = with_engine(x, &Expr::zero(), false, &logt, |eng| {
        let s = eng.ser(e, (n + 1) * eng.d)?;
        Ok(s.at(n * eng.d))
    })?;
    Ok(mul2(&c.to_expr(), &pow(x, &Expr::int(n))))
}

/// The expansion as a list of `(coefficient, exponent)` pairs of the small
/// variable plus the precision (`None` when exact). Coefficients may contain
/// `log(x - x0)`.
pub fn series_coeffs(e: &Expr, x: &Expr, x0: &Expr, n: i64, minus: bool) -> SRes<(Vec<(Expr, Number)>, Option<Number>)> {
    let logt = logt_symbol();
    let (s, d) = with_engine(x, x0, minus, &logt, |eng| eng.ser(e, n * eng.d))?;
    let logval = crate::functions::log(&small_var(x, x0, minus));
    let map = single(&logt, &logval);
    let mut v = Vec::new();
    for (i, c) in s.c.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let mut ce = c.to_expr();
        if depends(&ce, &logt) {
            ce = xreplace(&ce, &map);
        }
        v.push((ce, Number::rational_i64(s.val + i as i64, d).unwrap()));
    }
    let p = if s.exact { None } else { Some(Number::rational_i64(s.prec, d).unwrap()) };
    Ok((v, p))
}

#[cfg(test)]
mod tests;
