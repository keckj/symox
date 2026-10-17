//! Convergence of infinite series and products, after SymPy's
//! `Sum.is_convergent`, `Sum.is_absolutely_convergent` and
//! `Product.is_convergent` (sympy/concrete/summations.py and products.py).
//!
//! Exact decisions come first: finite ranges, constant and rational terms
//! (p-series on the degrees) and hypergeometric terms (ratio and Raabe tests
//! on the exact ratio). The remaining terms go through SymPy's test
//! sequence on the native limit engine: the divergence test, the
//! alternating series test, the bounded-times-absolutely-convergent test,
//! the comparison with `k**p*log(k)**q*log(log(k))**r` for terms of
//! eventually constant sign, and the ratio, Raabe and root tests. `Err` is
//! SymPy's NotImplementedError. Unlike SymPy, a p-series comparison never
//! claims divergence for an oscillating term (SymPy's `O()` drops bounded
//! factors, so it reports `Sum(sin(k)/k)` as divergent).

use std::cmp::Ordering;

use rustc_hash::FxHashMap;

use crate::binders::depends_on as has;
use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::functions::{abs, log};
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::subs::{subs1, xreplace};

const MANY_SYMBOLS: &str = "convergence checking for more than one symbol containing series is not handled";

fn oo() -> Expr {
    Expr::constant(Constant::Infinity)
}

fn neg_oo() -> Expr {
    Expr::constant(Constant::NegInfinity)
}

fn is_finite_limit(e: &Expr) -> bool {
    crate::assumptions::is_finite(e) == Some(true)
}

fn single_symbol(e: &Expr, k: &Expr) -> Result<(), String> {
    if crate::traverse::free_symbols(e).iter().any(|s| s != k) {
        return Err(MANY_SYMBOLS.into());
    }
    Ok(())
}

fn rational_in(e: &Expr, k: &Expr) -> bool {
    crate::inspection::classify(e, std::slice::from_ref(k), crate::inspection::PolynomialKind::Rational) == Some(true)
}

fn contains_limit(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Limit, _)) || e.args().iter().any(contains_limit)
}

/// `lim_{k -> oo} e` when it is a real number or `+-oo`.
fn lim(e: &Expr, k: &Expr) -> Option<Expr> {
    let r = crate::limits::limit(e, k, &oo(), crate::limits::Dir::Minus).ok()?;
    if r.is_constant(Constant::Infinity) || r.is_constant(Constant::NegInfinity) {
        return Some(r);
    }
    let real = crate::traverse::free_symbols(&r).is_empty() && !contains_limit(&r) && !r.is_constant(Constant::NaN)
        && crate::assumptions::is_real(&r) == Some(true) && is_finite_limit(&r);
    real.then_some(r)
}

/// How a real number (or `+-oo`) compares with the integer `c`.
fn compare(e: &Expr, c: i64) -> Option<Ordering> {
    if e.is_constant(Constant::Infinity) {
        return Some(Ordering::Greater);
    }
    if e.is_constant(Constant::NegInfinity) {
        return Some(Ordering::Less);
    }
    let d = sub(e, &Expr::int(c));
    if crate::assumptions::is_zero(&d) == Some(true) {
        return Some(Ordering::Equal);
    }
    if crate::assumptions::is_positive(&d) == Some(true) {
        return Some(Ordering::Greater);
    }
    if crate::assumptions::is_negative(&d) == Some(true) {
        return Some(Ordering::Less);
    }
    let v = crate::evalf::eval_f64(&d)?;
    (v.is_finite() && v.abs() > 1e-12).then(|| v.partial_cmp(&0.0)).flatten()
}

/// Degrees of numerator and denominator of a rational function of `k`,
/// the ratio of their leading coefficients and whether it is zero.
fn rational_degrees(f: &Expr, k: &Expr) -> Option<(u32, u32, Expr, bool)> {
    let cancelled = crate::poly::api::cancel(f).ok()?;
    let (p, q) = crate::simplify::numer_denom(&cancelled);
    let gens = std::slice::from_ref(k);
    let (p, q) = (crate::poly::polyobj::Poly::from_expr(&p, gens).ok()?, crate::poly::polyobj::Poly::from_expr(&q, gens).ok()?);
    let zero = p.as_expr().is_zero();
    Some((p.degree(0).unwrap_or(0), q.degree(0)?, div(&p.lc(), &q.lc()), zero))
}

/// The tail test for a term rational in `k`: a nonzero limit diverges,
/// otherwise `g ~ C/k**d` and the p-series test decides (`d > 1`).
fn rational_sum_tail(g: &Expr, k: &Expr) -> Option<bool> {
    let (p, q, _, zero) = rational_degrees(g, k)?;
    Some(zero || q >= p + 2)
}

/// SymPy's `k -> -k` transform of a range starting at `-oo`: `(-oo, oo)` is
/// split at zero and `(-oo, b)` becomes `(-b, oo)`.
fn negative_range(
    g: &Expr,
    k: &Expr,
    a: &Expr,
    b: &Expr,
    check: fn(&Expr, &Expr, &Expr, &Expr) -> Result<bool, String>,
) -> Option<Result<bool, String>> {
    if *a != neg_oo() {
        return None;
    }
    if *b == oo() {
        return Some(check(g, k, &Expr::zero(), &oo()).and_then(|first| Ok(first && check(g, k, &neg_oo(), &Expr::zero())?)));
    }
    let reflected = xreplace(g, &FxHashMap::from_iter([(k.clone(), neg(k))]));
    Some(check(&reflected, k, &neg(b), &oo()))
}

/// The term on a positive integer Dummy index, as SymPy's `is_convergent`
/// does before running its tests (so `Abs((-1)**k/k)` becomes `1/k`).
fn on_positive_integers(g: &Expr, k: &Expr) -> (Expr, Expr) {
    use crate::symbol::Assumptions as A;
    let Some(symbol) = k.as_sym() else { return (g.clone(), k.clone()) };
    let dummy = Expr::from_kind(Kind::Sym(crate::symbol::Symbol::dummy(symbol.name(), A::INTEGER.with(A::POSITIVE))));
    (xreplace(g, &FxHashMap::from_iter([(k.clone(), dummy.clone())])), dummy)
}

/// `Sum(g, (k, a, b)).is_convergent()`; `Err` is SymPy's NotImplementedError.
pub fn sum_is_convergent(g: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Result<bool, String> {
    single_symbol(g, k)?;
    if is_finite_limit(a) && is_finite_limit(b) {
        return Ok(true);
    }
    if let Some(halves) = negative_range(g, k, a, b, sum_is_convergent) {
        return halves;
    }
    let (term, index) = on_positive_integers(g, k);
    let term = drop_positive_abs(&term);
    series_tail(&term, &index, a, 0)
        .ok_or_else(|| format!("The algorithm to find the Sum convergence of {g} is not yet implemented"))
}

/// `Sum(g, (k, a, b)).is_absolutely_convergent()`: the convergence of
/// `Sum(Abs(g), (k, a, b))`.
pub fn sum_is_absolutely_convergent(g: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Result<bool, String> {
    sum_is_convergent(&abs(g), k, a, b)
}

/// Factorials, double factorials, binomials and rising factorials of
/// nonnegative integers are positive, as SymPy's `_eval_is_positive`
/// (the native assumptions do not see through these heads yet).
fn positive_special(e: &Expr) -> bool {
    use crate::special::SpecialFn as S;
    let nonneg_int = |x: &Expr| crate::assumptions::is_integer(x) == Some(true) && crate::assumptions::is_nonnegative(x) == Some(true);
    match e.kind() {
        Kind::Pow(b, x) => positive_special(b) && crate::assumptions::is_real(x) == Some(true),
        Kind::Mul(c, facs) => c.is_positive() && facs.iter().all(|(b, x)| positive_special(&crate::ops::mul_factor(b, x))),
        Kind::Func(Func::Special(S::Factorial | S::Factorial2), args) => nonneg_int(&args[0]),
        Kind::Func(Func::Special(S::Binomial), args) => nonneg_int(&args[1]) && nonneg_int(&sub(&args[0], &args[1])),
        Kind::Func(Func::Special(S::RisingFactorial), args) => crate::assumptions::is_positive(&args[0]) == Some(true) && nonneg_int(&args[1]),
        _ => crate::assumptions::is_positive(e) == Some(true),
    }
}

/// `Abs(f)` factors of a term whose argument is a positive special
/// function value lose their `Abs`, as SymPy's automatic `Abs` evaluation.
fn drop_positive_abs(t: &Expr) -> Expr {
    let strip = |f: &Expr| -> Expr {
        match f.kind() {
            Kind::Func(Func::Abs, args) if positive_special(&args[0]) => args[0].clone(),
            Kind::Pow(b, x) => match b.kind() {
                Kind::Func(Func::Abs, args) if positive_special(&args[0]) => pow(&args[0], x),
                _ => f.clone(),
            },
            _ => f.clone(),
        }
    };
    match t.kind() {
        Kind::Mul(c, facs) => {
            let mut parts = vec![Expr::num(c.clone())];
            parts.extend(facs.iter().map(|(b, x)| strip(&crate::ops::mul_factor(b, x))));
            mul(&parts)
        }
        _ => strip(t),
    }
}

/// Nested tests (absolute convergence of a factor) stop at this depth.
const DEPTH: usize = 3;

/// Convergence of `Sum(t, (k, a, oo))` for `k` a positive integer.
fn series_tail(t: &Expr, k: &Expr, a: &Expr, depth: usize) -> Option<bool> {
    if crate::assumptions::is_zero(t) == Some(true) {
        return Some(true);
    }
    if let Kind::Piecewise(pieces) = t.kind() {
        return piecewise_tail(pieces, k, a, depth);
    }
    if !has(t, k) {
        return crate::assumptions::is_zero(t);
    }
    if rational_in(t, k) {
        return rational_sum_tail(t, k);
    }
    if let Some(decided) = hypergeometric_tail(t, k) {
        return Some(decided);
    }
    general_tail(t, k, a, depth)
}

/// SymPy's Piecewise rule: the first piece whose condition reaches `oo`
/// decides; when none does the series is finite.
fn piecewise_tail(pieces: &[(Expr, Expr)], k: &Expr, a: &Expr, depth: usize) -> Option<bool> {
    for (piece, condition) in pieces {
        let reaches_oo = match condition.kind() {
            Kind::Bool(true) => true,
            Kind::Rel(RelOp::Gt | RelOp::Ge, l, r) => l == k && !has(r, k),
            Kind::Rel(RelOp::Lt | RelOp::Le, l, r) => r == k && !has(l, k),
            _ => return None,
        };
        if reaches_oo {
            return series_tail(piece, k, a, depth);
        }
    }
    Some(true)
}

/// Ratio and Raabe tests on the exact ratio of a hypergeometric term.
fn hypergeometric_tail(t: &Expr, k: &Expr) -> Option<bool> {
    let ratio = crate::recurr::hypersimp(t, k)?;
    if !rational_in(&ratio, k) {
        return None;
    }
    let (p, q, lc, zero) = rational_degrees(&ratio, k)?;
    if zero || p < q {
        return Some(true);
    }
    if p > q {
        return Some(false);
    }
    match compare(&abs(&lc), 1)? {
        Ordering::Less => Some(true),
        Ordering::Greater => Some(false),
        Ordering::Equal if lc.is_one() => raabe(&ratio, k),
        Ordering::Equal => None,
    }
}

/// Raabe's test for a ratio `t(k + 1)/t(k)` tending to one.
fn raabe(ratio: &Expr, k: &Expr) -> Option<bool> {
    let test = mul(&[k.clone(), sub(&div(&Expr::one(), ratio), &Expr::one())]);
    let test = crate::poly::api::cancel(&test).unwrap_or(test);
    match compare(&lim(&test, k)?, 1)? {
        Ordering::Greater => Some(true),
        Ordering::Less => Some(false),
        Ordering::Equal => None,
    }
}

fn factors(t: &Expr) -> Vec<Expr> {
    if t.is_mul() { t.args() } else { vec![t.clone()] }
}

/// The factor `(-1)**(k + c)` of an alternating term and the rest.
fn alternating_part(t: &Expr, k: &Expr) -> Option<Expr> {
    let all = factors(t);
    let position = all.iter().position(|f| match f.kind() {
        Kind::Pow(b, e) => *b == Expr::minus_one() && !has(&sub(e, k), k),
        _ => false,
    })?;
    let rest: Vec<Expr> = all.iter().enumerate().filter(|(i, _)| *i != position).map(|(_, f)| f.clone()).collect();
    Some(mul(&rest))
}

/// Whether `q` is nonincreasing on the real interval `[a, oo)`, as SymPy's
/// `is_decreasing(q, Interval(a, oo))`: the derivative is nonpositive.
fn decreasing(q: &Expr, k: &Expr, a: &Expr) -> bool {
    use crate::symbol::Assumptions as A;
    if !a.is_number() {
        return false;
    }
    let t = Expr::from_kind(Kind::Sym(crate::symbol::Symbol::dummy("t", A::NONNEGATIVE)));
    let d = subs1(&crate::diff::diff(q, k), k, &add(&[a.clone(), t]));
    if crate::assumptions::is_nonpositive(&d) == Some(true) {
        return true;
    }
    let cancelled = crate::poly::api::cancel(&d).unwrap_or(d);
    let (p, r) = crate::simplify::numer_denom(&cancelled);
    let (p, r) = (crate::expand::expand(&p), crate::expand::expand(&r));
    let sign = |e: &Expr| {
        if crate::assumptions::is_nonpositive(e) == Some(true) {
            Some(-1)
        } else if crate::assumptions::is_nonnegative(e) == Some(true) {
            Some(1)
        } else {
            None
        }
    };
    matches!((sign(&p), sign(&r)), (Some(-1), Some(1)) | (Some(1), Some(-1)))
        || (crate::assumptions::is_zero(&p) == Some(true))
}

/// A bounded oscillating factor: `sin`, `cos` of a real argument,
/// `(-1)**n` for an integer `n`, or the absolute value of one.
fn bounded(f: &Expr, k: &Expr) -> bool {
    if !has(f, k) {
        return false;
    }
    match f.kind() {
        Kind::Func(Func::Sin | Func::Cos, args) => crate::assumptions::is_real(&args[0]) == Some(true),
        Kind::Func(Func::Abs, args) => bounded(&args[0], k),
        Kind::Pow(b, e) if *b == Expr::minus_one() => crate::assumptions::is_integer(e) == Some(true),
        Kind::Pow(b, e) => bounded(b, k) && e.as_number().is_some_and(|n| !n.is_negative()),
        _ => false,
    }
}

/// SymPy's divergence, alternating, bounded-times-convergent, comparison,
/// ratio, Raabe and root tests on the native limit engine.
fn general_tail(t: &Expr, k: &Expr, a: &Expr, depth: usize) -> Option<bool> {
    // divergence test
    let nonzero = |l: &Expr| crate::assumptions::is_zero(l) == Some(false);
    if let Some(l) = lim(t, k) {
        if nonzero(&l) {
            return Some(false);
        }
    }
    let magnitude = abs(t);
    if let Some(l) = lim(&magnitude, k) {
        if nonzero(&l) {
            return Some(false);
        }
    }
    // alternating series test
    if let Some(q) = alternating_part(t, k) {
        if lim(&q, k).is_some_and(|l| l.is_zero()) && decreasing(&q, k, a) {
            return Some(true);
        }
    }
    // bounded factors times an absolutely convergent series
    if depth < DEPTH {
        let (osc, rest): (Vec<Expr>, Vec<Expr>) = factors(t).into_iter().partition(|f| bounded(f, k));
        if !osc.is_empty() && has(&mul(&rest), k) && series_tail(&abs(&mul(&rest)), k, a, depth + 1) == Some(true) {
            return Some(true);
        }
    }
    // comparison with k**p*log(k)**q*log(log(k))**r
    if let Some(decided) = comparison(t, k) {
        return Some(decided);
    }
    // ratio and Raabe tests
    let ratio = shift_ratio(t, k);
    let limit = lim(&ratio, k).map(|l| (ratio.clone(), l)).or_else(|| {
        let ratio = crate::simplify::simplify(&div(&subs1(t, k, &add(&[k.clone(), Expr::one()])), t));
        lim(&ratio, k).map(|l| (ratio, l))
    });
    if let Some((ratio, l)) = limit {
        match compare(&abs(&l), 1) {
            Some(Ordering::Less) => return Some(true),
            Some(Ordering::Greater) => return Some(false),
            Some(Ordering::Equal) if l.is_one() => {
                if let Some(decided) = raabe(&ratio, k) {
                    return Some(decided);
                }
            }
            _ => {}
        }
    }
    // root test
    if let Some(l) = lim(&pow(&magnitude, &div(&Expr::one(), k)), k) {
        match compare(&l, 1) {
            Some(Ordering::Less) => return Some(true),
            Some(Ordering::Greater) => return Some(false),
            _ => {}
        }
    }
    None
}

/// `t(k + 1)/t(k)` factor by factor: hypergeometric factors (factorials,
/// binomials, gamma functions) through their exact rational ratio, the
/// others shifted and recombined by `powsimp`.
fn shift_ratio(t: &Expr, k: &Expr) -> Expr {
    let factors: Vec<Expr> = match t.kind() {
        Kind::Mul(_, facs) => facs.iter().map(|(b, e)| crate::ops::mul_factor(b, e)).collect(),
        _ => vec![t.clone()],
    };
    let next = add(&[k.clone(), Expr::one()]);
    let parts: Vec<Expr> = factors.iter().filter(|f| has(f, k)).map(|f| {
        crate::recurr::hypersimp(f, k).filter(|r| rational_in(r, k))
            .unwrap_or_else(|| div(&subs1(f, k, &next), f))
    }).collect();
    crate::simp::powsimp(&mul(&parts))
}

/// The sign of `t` for large `k`, when it is eventually constant.
fn eventual_sign(t: &Expr, k: &Expr) -> Option<i32> {
    if positive_special(t) {
        return Some(1);
    }
    if let Some(l) = lim(&drop_positive_abs(&div(t, &abs(t))), k) {
        if l.is_one() {
            return Some(1);
        }
        if l == Expr::minus_one() {
            return Some(-1);
        }
    }
    if crate::traverse::free_symbols(t).len() != 1 || has_oscillation(t, k) {
        return None;
    }
    match crate::limits::gruntz::Gruntz::new(k).sign(t).ok()? {
        0 => None,
        s => Some(s),
    }
}

fn has_oscillation(e: &Expr, k: &Expr) -> bool {
    let own = match e.kind() {
        Kind::Func(Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc, args) => {
            // trigonometric functions of a vanishing argument do not oscillate
            has(&args[0], k) && lim(&args[0], k).map_or(true, |l| !is_finite_limit(&l))
        }
        Kind::Pow(b, x) => has(x, k) && compare(b, 0) != Some(Ordering::Greater),
        _ => false,
    };
    own || e.args().iter().any(|a| has_oscillation(a, k))
}

/// Comparison with `1/(k**p*log(k)**q*log(log(k))**r)` for a term of
/// eventually constant sign: SymPy's p-series and logarithmic tests.
fn comparison(t: &Expr, k: &Expr) -> Option<bool> {
    let s = eventual_sign(t, k)?;
    let u = if s < 0 { neg(t) } else { t.clone() };
    let log_u = log(&u);
    let log_k = log(k);
    let loglog_k = log(&log_k);
    let p = lim(&div(&log_u, &log_k), k)?;
    match compare(&p, -1)? {
        Ordering::Less => return Some(true),
        Ordering::Greater => return Some(false),
        Ordering::Equal => {}
    }
    // log(u*k) rather than log(u) + log(k): the product cancels exactly,
    // the sum of logarithms can leave an unrecognised zero numerator
    let u = mul(&[u, k.clone()]);
    let q = lim(&div(&log(&u), &loglog_k), k)?;
    match compare(&q, -1)? {
        Ordering::Less => return Some(true),
        Ordering::Greater => return Some(false),
        Ordering::Equal => {}
    }
    let r = lim(&div(&log(&mul(&[u, log_k])), &log(&loglog_k)), k)?;
    Some(compare(&r, -1)? == Ordering::Less)
}

/// `Product(f, (k, a, b)).is_convergent()`: the convergence of
/// `Sum(log(f))`, falling back to the absolute convergence of
/// `Sum(f - 1)`; decided exactly for factors rational in `k`.
pub fn is_convergent(f: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Result<bool, String> {
    single_symbol(f, k)?;
    if is_finite_limit(a) && is_finite_limit(b) {
        return Ok(true);
    }
    if let Some(halves) = negative_range(f, k, a, b, is_convergent) {
        return halves;
    }
    let unsupported = || format!("The algorithm to find the product convergence of {f} is not yet implemented");
    if let Some(decided) = rational_product(f, k) {
        return decided.ok_or_else(unsupported);
    }
    if let Kind::Func(Func::Exp, g) = f.kind() {
        if let Ok(decided) = sum_is_convergent(&g[0], k, a, b) {
            return Ok(decided);
        }
    }
    if let Ok(decided) = sum_is_convergent(&log(f), k, a, b) {
        return Ok(decided);
    }
    let excess = sub(f, &Expr::one());
    if sum_is_absolutely_convergent(&excess, k, a, b) == Ok(true) {
        return Ok(true);
    }
    // f - 1 of eventually constant sign: the product converges (to a
    // nonzero limit, as Sum(log(f)) does) exactly when Sum(f - 1) does
    let (tail, index) = on_positive_integers(&excess, k);
    if eventual_sign(&tail, &index).is_some() {
        if let Ok(decided) = sum_is_convergent(&excess, k, a, b) {
            return Ok(decided);
        }
    }
    Err(unsupported())
}

/// Constant and rational factors: `log(f)` tends to a nonzero limit unless
/// `f -> 1`, and then `log(f) ~ f - 1`. `None` when `f` is not rational.
fn rational_product(f: &Expr, k: &Expr) -> Option<Option<bool>> {
    if !has(f, k) {
        if crate::assumptions::is_zero(f) != Some(false) {
            return Some(None);
        }
        return Some(crate::assumptions::is_zero(&sub(f, &Expr::one())));
    }
    if !rational_in(f, k) {
        return None;
    }
    let Some((p, q, ratio, _)) = rational_degrees(f, k) else { return Some(None) };
    if p != q || crate::assumptions::is_zero(&sub(&ratio, &Expr::one())) == Some(false) {
        return Some(Some(false));
    }
    Some(rational_sum_tail(&sub(f, &Expr::one()), k))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn sum(f: &str) -> Result<bool, String> {
        sum_is_convergent(&parse_str(f).unwrap(), &Expr::symbol("k"), &Expr::one(), &oo())
    }

    #[test]
    fn series_tests_match_sympy() {
        for (f, expected) in [
            ("1/k**2", Some(true)), ("1/k", Some(false)), ("(-1)**k/k", Some(true)), ("(-1)**k", Some(false)),
            ("1/sqrt(k)", Some(false)), ("sin(k)/k**2", Some(true)), ("log(k)/k**2", Some(true)),
            ("exp(-k)", Some(true)), ("2**k/k", Some(false)), ("sin(1/k)", Some(false)),
        ] {
            assert_eq!(sum(f).ok(), expected, "{f}");
        }
    }
}
