//! The Landau order term `O(expr, (x, x0), ...)` and its absorbing rules.
//!
//! Representation: `Kind::Func(Func::Order, [expr, x1, p1, x2, p2, ...])`
//! with the variables sorted. `expr` is normalised to its leading monomial
//! (`O(2*x + x**2) = O(x)`, `O(x*y, x) = O(x)`), which makes the
//! containment test cheap: `t ∈ O(e)` iff `t/e` stays bounded at the point.
//!
//! The hooks [`add_hook`], [`mul_hook`] and [`pow_hook`] are called by the
//! canonical constructors in `ops` when an operand contains an order term:
//! `x + O(x**2)` keeps, `x**3 + O(x**2) = O(x**2)`, `x*O(x) = O(x**2)`,
//! `O(x)**2 = O(x**2)`.

use rustc_hash::FxHashMap;

use super::engine::{degree_in, depends};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, add_from_terms, as_coeff_mul, div, mul, mul2, mul_factor, neg, pow};
use crate::subs::xreplace;
use crate::traverse::free_symbols;

#[inline]
pub fn is_order(e: &Expr) -> bool {
    e.is_func(&Func::Order)
}

/// `(expr, [(var, point)])` of an order term.
pub fn order_parts(e: &Expr) -> Option<(Expr, Vec<(Expr, Expr)>)> {
    match e.kind() {
        Kind::Func(Func::Order, args) => {
            let vars = args[1..]
                .chunks(2)
                .map(|p| (p[0].clone(), p.get(1).cloned().unwrap_or_else(Expr::zero)))
                .collect();
            Some((args[0].clone(), vars))
        }
        _ => None,
    }
}

/// `Order.__new__` refuses variables at different points
/// (`O(x*y, (x, 0), (y, 1))`, also reached through `O(x*y).subs(x, x + 1)`).
pub const DIFFERENT_POINTS: &str = "Multivariable orders at different points are not supported.";

/// Whether the points of an order's variables are not all equal.
pub fn different_points(vars: &[(Expr, Expr)]) -> bool {
    vars.windows(2).any(|w| w[0].1 != w[1].1)
}

/// `Err(DIFFERENT_POINTS)` when `e` holds an order term whose variables
/// sit at different points (the reference raises `NotImplementedError`
/// when such a term would be built).
pub fn check_order_points(e: &Expr) -> Result<(), &'static str> {
    if !e.may_have_order() {
        return Ok(());
    }
    if let Some((_, vars)) = order_parts(e) {
        if different_points(&vars) {
            return Err(DIFFERENT_POINTS);
        }
    }
    e.args().iter().try_for_each(check_order_points)
}

fn raw(expr: Expr, vars: &[(Expr, Expr)]) -> Expr {
    let mut args = Vec::with_capacity(1 + 2 * vars.len());
    args.push(expr);
    for (v, p) in vars {
        args.push(v.clone());
        args.push(p.clone());
    }
    Expr::from_kind(Kind::Func(Func::Order, args))
}

/// Build from flat arguments `[expr, x1, p1, ...]`.
pub fn order_from_args(args: &[Expr]) -> Expr {
    if args.is_empty() {
        return raw(Expr::one(), &[]);
    }
    let vars: Vec<(Expr, Expr)> = args[1..]
        .chunks(2)
        .map(|p| (p[0].clone(), p.get(1).cloned().unwrap_or_else(Expr::zero)))
        .collect();
    order(&args[0], &vars)
}

fn one_map(x: &Expr, v: Expr) -> FxHashMap<Expr, Expr> {
    let mut m = FxHashMap::default();
    m.insert(x.clone(), v);
    m
}

/// Map the point `p` of `x` to `0`.
pub(crate) fn to_zero(e: &Expr, x: &Expr, p: &Expr) -> Expr {
    if p.is_zero() {
        return e.clone();
    }
    let v = if p.is_constant(Constant::Infinity) {
        pow(x, &Expr::minus_one())
    } else if p.is_constant(Constant::NegInfinity) {
        neg(&pow(x, &Expr::minus_one()))
    } else {
        add(&[x.clone(), p.clone()])
    };
    xreplace(e, &one_map(x, v))
}

/// Inverse of [`to_zero`].
pub(crate) fn from_zero(e: &Expr, x: &Expr, p: &Expr) -> Expr {
    if p.is_zero() {
        return e.clone();
    }
    let v = if p.is_constant(Constant::Infinity) {
        pow(x, &Expr::minus_one())
    } else if p.is_constant(Constant::NegInfinity) {
        neg(&pow(x, &Expr::minus_one()))
    } else {
        add(&[x.clone(), neg(p)])
    };
    xreplace(e, &one_map(x, v))
}

/// `e = c * x**a * log(x)**b` with `c` free of `x`: returns `(a, b)`.
pub(crate) fn monomial_x(e: &Expr, x: &Expr) -> Option<(Number, i64)> {
    if !depends(e, x) {
        return Some((Number::ZERO, 0));
    }
    if e == x {
        return Some((Number::ONE, 0));
    }
    match e.kind() {
        Kind::Pow(b, p) => {
            let k = p.as_number()?;
            if b == x && k.is_rational() {
                return Some((k.clone(), 0));
            }
            if is_log_of(b, x) {
                return Some((Number::ZERO, k.to_i64()?));
            }
            None
        }
        Kind::Func(Func::Log, a) if a[0] == *x => Some((Number::ZERO, 1)),
        Kind::Mul(_, facs) => {
            let mut a = Number::ZERO;
            let mut l = 0i64;
            for (b, p) in facs {
                let (da, dl) = monomial_x(&mul_factor(b, p), x)?;
                a = a.add(&da);
                l += dl;
            }
            Some((a, l))
        }
        _ => None,
    }
}

fn is_log_of(e: &Expr, x: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Log, a) if a[0] == *x)
}

fn build_monomial(x: &Expr, a: &Number, b: i64) -> Expr {
    let xa = pow(x, &Expr::num(a.clone()));
    if b == 0 {
        xa
    } else {
        mul2(&xa, &pow(&crate::functions::log(x), &Expr::int(b)))
    }
}

/// Leading behaviour `(a, b)` of `e ~ c x**a log(x)**b` as `x -> 0+`.
fn lead_exponents(e: &Expr, x: &Expr) -> Option<(Number, i64)> {
    if let Some(r) = monomial_x(e, x) {
        return Some(r);
    }
    if e.may_have_order() {
        return None;
    }
    let lt = super::leading_term(e, x, &Expr::zero(), false).ok()?;
    // a coefficient that is not a polynomial in `log(x)` (`log(log(x))`)
    // has no such leading behaviour
    let b = match degree_in(&lt.coeff_logt, &lt.logt) {
        Some(b) => b,
        None if depends(&lt.coeff_logt, &lt.logt) => return None,
        None => 0,
    };
    Some((lt.exp, b))
}

thread_local! {
    static TRIG_ACTIVE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static CONTAINS_ACTIVE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Runs `f` unless `flag` is already set on this thread (the limits taken
/// below build order terms themselves).
fn reentry_guarded<T>(flag: &'static std::thread::LocalKey<std::cell::Cell<bool>>, f: impl FnOnce() -> Option<T>) -> Option<T> {
    if flag.with(|a| a.replace(true)) {
        return None;
    }
    let r = f();
    flag.with(|a| a.set(false));
    r
}

/// `sin(a)`, `cos(a)` or `tan(a)` (to a nonzero integer power) whose
/// argument is infinite at `x = 0` (`a.subs(x, 0)` is `zoo` or `oo`) and
/// tends to `oo` or `-oo` as `x -> 0+`: the reference leading term is the
/// constant `AccumBounds(-1, 1)` (`AccumBounds(-oo, oo)` for `tan`), which
/// `Order` strips as a factor free of `x`.
fn bounded_trig_factor(b: &Expr, k: &Expr, x: &Expr) -> bool {
    let Kind::Func(Func::Sin | Func::Cos | Func::Tan, a) = b.kind() else { return false };
    if a.len() != 1 || !depends(&a[0], x) {
        return false;
    }
    if !k.as_number().is_some_and(|n| n.is_integer() && !n.is_zero()) {
        return false;
    }
    let a0 = crate::subs::subs1(&a[0], x, &Expr::zero());
    if !(a0.is_constant(Constant::ComplexInfinity) || a0.is_constant(Constant::Infinity) || a0.is_constant(Constant::NegInfinity)) {
        return false;
    }
    reentry_guarded(&TRIG_ACTIVE, || crate::limits::limit(&a[0], x, &Expr::zero(), crate::limits::Dir::Plus).ok())
        .is_some_and(|l| l.is_constant(Constant::Infinity) || l.is_constant(Constant::NegInfinity))
}

/// `e0` (a point moved to 0) without its bounded trigonometric factors
/// (`O(x*sin(x), (x, oo)) = O(x, (x, oo))`, `O(tan(x), (x, oo)) = O(1, (x, oo))`).
/// Admitted only when the remaining factors have a power or logarithmic
/// leading term, so `Mul._eval_as_leading_term` raises no `PoleError`;
/// otherwise (`sin(1/x)*exp(1/x)`) the expression is kept, as the
/// reference does when every factor is a function.
fn drop_bounded_trig(e0: &Expr, x: &Expr) -> Option<Expr> {
    if e0.may_have_order() {
        return None;
    }
    let facs: Vec<(Expr, Expr)> = match e0.kind() {
        Kind::Mul(c, facs) => {
            let mut v: Vec<(Expr, Expr)> = facs.to_vec();
            v.push((Expr::num(c.clone()), Expr::one()));
            v
        }
        Kind::Pow(b, k) => vec![(b.clone(), k.clone())],
        _ => vec![(e0.clone(), Expr::one())],
    };
    let mut dropped = false;
    let mut rest = Vec::with_capacity(facs.len());
    for (b, k) in &facs {
        if bounded_trig_factor(b, k, x) {
            dropped = true;
        } else {
            rest.push(mul_factor(b, k));
        }
    }
    if !dropped {
        return None;
    }
    let r = mul(&rest);
    let simple = !depends(&r, x) || lead_exponents(&r, x).is_some() || super::order_log::log_coefficient_order(&r, x).is_some();
    simple.then_some(r)
}

/// The `PoleError` branch of `Order.__new__` (no power leading term): a
/// function, or a product of functions only, is kept; otherwise every
/// factor is replaced by its own order, a power with an exponent depending
/// on `x` being written `exp(k*log(b))` and a bounded trigonometric factor
/// (`AccumBounds` leading term) dropping out. The factors free of `x` are
/// removed (`Expr.as_independent`).
fn pole_factors(e: &Expr, x: &Expr) -> Expr {
    let one = |f: Expr| super::order_dominant::pole_power(&f, x).unwrap_or(f);
    match e.kind() {
        Kind::Mul(c, facs) => {
            let all_functions = c.is_one() && facs.iter().all(|(b, k)| k.is_one() && matches!(b.kind(), Kind::Func(..)));
            let keep = |(b, k): &&(Expr, Expr)| (depends(b, x) || depends(k, x)) && (all_functions || !bounded_trig_factor(b, k, x));
            // a sum factor takes its own order (`Order(x + 1) = O(1)`)
            let own = |(b, k): &(Expr, Expr)| {
                if !all_functions && k.is_one() && b.is_add() && !b.may_have_order() {
                    order_parts(&order(b, &[(x.clone(), Expr::zero())])).map_or_else(|| b.clone(), |(oe, _)| oe)
                } else {
                    one(mul_factor(b, k))
                }
            };
            mul(&facs.iter().filter(keep).map(own).collect::<Vec<_>>())
        }
        _ => one(as_coeff_mul(e).1),
    }
}

/// Canonical constructor `O(expr, (x1, p1), ...)`.
pub fn order(expr: &Expr, vars: &[(Expr, Expr)]) -> Expr {
    let mut expr = expr.clone();
    let mut vars: Vec<(Expr, Expr)> = vars.to_vec();
    if let Some((inner, ivars)) = order_parts(&expr) {
        // Order.__new__ on an order: the variables of both, the inner ones
        // appended (`Order(Order(x), y) = O(x, x, y)`).
        if vars.iter().all(|v| ivars.iter().any(|w| w.0 == v.0)) {
            return expr;
        }
        expr = inner;
        for v in ivars {
            if !vars.iter().any(|w| w.0 == v.0) {
                vars.push(v);
            }
        }
    }
    if expr.may_have_order() {
        // An order term inside the expression (`(x + O(x**2))**2` built by
        // `O(x**2).subs(x, x + O(x**2))`): the expanded expression with
        // every order term replaced by its own expression keeps the same
        // leading terms, as `extract_leading_order` does on the reference.
        let ex = if expr.is_add() && expr.args().iter().all(|t| is_order(t) || !t.may_have_order()) {
            expr.clone()
        } else {
            crate::expand::expand(&expr)
        };
        // (`expand` builds the products of order terms again, so that an
        // order absorbs its factors: `x*O(x**2)` is `O(x**3)`)
        let terms: Vec<Expr> = add_args(&ex).into_iter().map(|t| order_parts(&t).map_or(t, |(oe, _)| oe)).collect();
        let r = add(&terms);
        expr = if r.may_have_order() { remove_o(&ex) } else { r };
    }
    if vars.is_empty() {
        vars = free_symbols(&expr).into_iter().map(|s| (s, Expr::zero())).collect();
    }
    vars.sort_by(|a, b| a.0.cmp_struct(&b.0));
    vars.dedup_by(|a, b| a.0 == b.0);
    if vars.iter().any(|(v, _)| !v.is_symbol()) {
        return raw(expr, &vars);
    }
    if expr.is_constant(Constant::NaN) {
        // Order.__new__: `O(nan)` is `nan`
        return Expr::nan();
    }
    if expr.is_zero() {
        return raw(expr, &vars);
    }
    if vars.is_empty() {
        return raw(Expr::one(), &[]);
    }
    if vars.len() == 1 {
        let (x, p) = &vars[0];
        let mut e0 = to_zero(&expr, x, p);
        if let Kind::Mul(_, facs) = e0.kind() {
            // factors whose `x` is bound (`Integral(x, (x, 1, 2))`) are
            // constants, dropped by `as_independent` as the free ones
            let facs: Vec<Expr> = facs.iter().map(|(b, q)| mul_factor(b, q)).collect();
            if facs.iter().any(|f| depends(f, x) && !free_symbols(f).contains(x)) {
                e0 = mul(&facs.into_iter().filter(|f| free_symbols(f).contains(x)).collect::<Vec<_>>());
            }
        }
        if let Some(d) = drop_bounded_trig(&e0, x) {
            e0 = d;
        }
        if let Some(lead) = super::order_log::log_coefficient_order(&e0, x) {
            return raw(from_zero(&lead, x, p), &vars);
        }
        let lead = match lead_exponents(&e0, x) {
            Some((a, b)) => build_monomial(x, &a, b),
            None => {
                let r = super::order_dominant::dominant_terms(&e0, x).unwrap_or(e0);
                pole_factors(&r, x)
            }
        };
        return raw(from_zero(&lead, x, p), &vars);
    }
    // several variables: the expanded expression (points moved to the
    // origin); in a product the factors free of all variables are dropped
    // and a factor of a single variable is replaced by its own leading
    // term (`O(x*y*sin(x)*Integral(x, (x, 1, 2))) = O(x**2*y, x, y)`).
    let syms: Vec<Expr> = vars.iter().map(|v| v.0.clone()).collect();
    let used = |f: &Expr| -> Vec<Expr> {
        if !syms.iter().any(|s| depends(f, s)) {
            return vec![];
        }
        let fs = free_symbols(f);
        syms.iter().filter(|s| fs.contains(s)).cloned().collect()
    };
    let factor_lead = |f: &Expr| -> Expr {
        let u = used(f);
        match u.len() {
            0 => Expr::one(),
            1 if monomial_x(f, &u[0]).is_none() && !f.may_have_order() => {
                order_parts(&order(f, &[(u[0].clone(), Expr::zero())])).map_or_else(|| f.clone(), |(oe, _)| oe)
            }
            _ => f.clone(),
        }
    };
    // the reference leading term of the whole product raises `PoleError`
    // when a factor has none (`exp(1/y)`), and a product of functions is
    // then kept: `O(exp(x)*exp(1/y)) = O(exp(x)*exp(1/y), x, y)`
    let pole = |f: &Expr| {
        let u = used(f);
        u.len() == 1
            && monomial_x(f, &u[0]).is_none()
            && !f.may_have_order()
            && matches!(super::leading_term(f, &u[0], &Expr::zero(), false), Err(super::ser::SErr::Pole(_)))
    };
    let strip = |t: &Expr| match t.kind() {
        Kind::Mul(_, facs) if facs.iter().all(|(b, p)| matches!(mul_factor(b, p).kind(), Kind::Func(..)) || used(&mul_factor(b, p)).is_empty())
            && facs.iter().any(|(b, p)| pole(&mul_factor(b, p))) =>
        {
            mul(&facs.iter().map(|(b, p)| mul_factor(b, p)).filter(|f| !used(f).is_empty()).collect::<Vec<_>>())
        }
        Kind::Mul(_, facs) => mul(&facs.iter().map(|(b, p)| factor_lead(&mul_factor(b, p))).collect::<Vec<_>>()),
        _ => factor_lead(&as_coeff_mul(t).1),
    };
    let mut e0 = expr.clone();
    for (x, p) in &vars {
        e0 = to_zero(&e0, x, p);
    }
    let e0 = crate::expand::expand(&e0);
    if !e0.is_add() {
        let mut e = strip(&e0);
        for (x, p) in &vars {
            e = from_zero(&e, x, p);
        }
        return raw(e, &vars);
    }
    // A sum keeps its leading terms only (reference `extract_leading_order`
    // on the expanded expression): `O(1 + x*y) = O(1)`,
    // `O(x + y + x*y) = O(x + y)`. A term whose ratio to another is not a
    // monomial is kept.
    let zvars: Vec<(Expr, Expr)> = vars.iter().map(|(x, _)| (x.clone(), Expr::zero())).collect();
    let mut kept: Vec<Expr> = Vec::new();
    for t in add_args(&e0) {
        let s = strip(&t);
        if s.is_zero() || kept.iter().any(|k| *k == s || contains_expr(k, &zvars, &s)) {
            continue;
        }
        kept.retain(|k| !contains_expr(&s, &zvars, k));
        kept.push(s);
    }
    let mut e = add(&kept);
    for (x, p) in &vars {
        e = from_zero(&e, x, p);
    }
    raw(e, &vars)
}

/// Is the term `t` contained in the order term `o` (i.e. `t = O(o)`)?
pub fn contains(o: &Expr, t: &Expr) -> bool {
    let Some((oe, vars)) = order_parts(o) else {
        return false;
    };
    if let Some((te, tvars)) = order_parts(t) {
        if tvars.len() != vars.len() || tvars.iter().zip(&vars).any(|(a, b)| a != b) {
            // The reference `Order.contains` takes the variables of the
            // order that has some when the other has none (`O(1)`): `O(1)`
            // holds `O(x)` (bounded at its point), `O(1/x)` holds `O(1)`.
            if vars.is_empty() && !tvars.is_empty() {
                return !oe.is_zero() && contains_expr(&Expr::one(), &tvars, &te);
            }
            if tvars.is_empty() && !vars.is_empty() {
                return contains_expr(&oe, &vars, &Expr::one());
            }
            if tvars.is_empty() && vars.is_empty() {
                return true;
            }
            return contains_multi(&oe, &vars, &te, &tvars, 0) == Some(true);
        }
        return contains_expr(&oe, &vars, &te) || (vars.len() > 1 && contains_multi(&oe, &vars, &te, &tvars, 0) == Some(true));
    }
    contains_expr(&oe, &vars, t) || (vars.len() > 1 && contains_plain_multi(&oe, &vars, t))
}

fn contains_plain_multi(oe: &Expr, vars: &[(Expr, Expr)], t: &Expr) -> bool {
    if t.may_have_order() {
        return false;
    }
    match order_parts(&order(t, vars)) {
        Some((te, tvars)) => contains_multi(oe, vars, &te, &tvars, 0) == Some(true),
        None => false,
    }
}

/// The reference `Order.contains` of an order term for several variables
/// (or different variable sets): equal expressions compare the variables,
/// a sum is held term by term, a sum order holds what one of its terms
/// holds, otherwise `limit(o/t, s, point) != 0` must hold alike for every
/// common variable `s`. `None` when undecided.
fn contains_multi(oe: &Expr, ovars: &[(Expr, Expr)], te: &Expr, tvars: &[(Expr, Expr)], depth: u32) -> Option<bool> {
    if te.is_zero() {
        return Some(true);
    }
    if depth > 4 || oe.may_have_order() || te.may_have_order() {
        return None;
    }
    let point = ovars.first().map_or_else(Expr::zero, |v| v.1.clone());
    if tvars.iter().chain(ovars).any(|v| v.1 != point) {
        return None;
    }
    if te == oe {
        return Some(tvars.iter().all(|v| ovars.contains(v)));
    }
    if te.is_add() {
        return Some(te.args().iter().all(|a| match order_parts(&order(a, ovars)) {
            Some((ae, av)) => contains_multi(oe, ovars, &ae, &av, depth + 1) == Some(true),
            None => false,
        }));
    }
    if oe.is_add() && point.is_zero() {
        return Some(oe.args().iter().any(|a| match order_parts(&order(a, ovars)) {
            Some((ae, av)) => contains_multi(&ae, &av, te, tvars, depth + 1) == Some(true),
            None => false,
        }));
    }
    let common: Vec<Expr> = if !ovars.is_empty() && !tvars.is_empty() {
        ovars.iter().filter(|v| tvars.iter().any(|w| w.0 == v.0)).map(|v| v.0.clone()).collect()
    } else if !ovars.is_empty() {
        ovars.iter().map(|v| v.0.clone()).collect()
    } else {
        tvars.iter().map(|v| v.0.clone()).collect()
    };
    if common.is_empty() {
        return None;
    }
    let ratio = div(oe, te);
    let mut r: Option<bool> = None;
    for s in &common {
        let l = reentry_guarded(&CONTAINS_ACTIVE, || {
            crate::limits::limit(&ratio, s, &point, crate::limits::Dir::Plus).ok()
        })?;
        if l.is_func(&Func::Limit) {
            return None;
        }
        let l = !l.is_zero();
        match r {
            None => r = Some(l),
            Some(prev) if prev != l => return None,
            _ => {}
        }
    }
    r
}

/// The reference `Order.contains` for terms without a power leading
/// ratio: `O(o)` contains `t` iff `limit(o/t, x, 0+)` is not `0` (`O(x)`
/// holds `exp(-1/x)`, `O(exp(1/x))` holds `x**-3`). `None` when the native
/// limit declines or leaves a symbol: the term is then kept, as before.
fn limit_contains(o0: &Expr, t0: &Expr, x: &Expr) -> Option<bool> {
    if o0.may_have_order() || t0.may_have_order() || t0.is_zero() {
        return None;
    }
    let l = reentry_guarded(&CONTAINS_ACTIVE, || crate::limits::limit(&div(o0, t0), x, &Expr::zero(), crate::limits::Dir::Plus).ok())?;
    if !free_symbols(&l).is_empty() || l.is_constant(Constant::NaN) {
        return None;
    }
    Some(!l.is_zero())
}

fn contains_expr(oe: &Expr, vars: &[(Expr, Expr)], t: &Expr) -> bool {
    if t.is_zero() {
        return true;
    }
    if vars.is_empty() {
        if oe.is_zero() {
            return false;
        }
        // O(1): bounded at the origin
        return free_symbols(t).iter().all(|s| matches!(monomial_x(t, s), Some((a, _)) if !a.is_negative()));
    }
    if vars.len() == 1 {
        let (x, p) = &vars[0];
        let t0 = to_zero(t, x, p);
        let o0 = to_zero(oe, x, p);
        if !depends(&t0, x) {
            return match monomial_x(&o0, x) {
                Some((a, b)) => a.is_negative() || (a.is_zero() && b >= 0),
                None => limit_contains(&o0, &t0, x).unwrap_or(false),
            };
        }
        let ratio = div(&t0, &o0);
        return match lead_exponents(&ratio, x) {
            Some((a, b)) => a.is_positive() || (a.is_zero() && b <= 0),
            None => limit_contains(&o0, &t0, x).unwrap_or(false),
        };
    }
    let ratio = div(t, oe);
    vars.iter().all(|(x, _)| matches!(monomial_x(&ratio, x), Some((a, b)) if a.is_positive() || (a.is_zero() && b <= 0)))
}

fn add_args(e: &Expr) -> Vec<Expr> {
    if e.is_add() {
        e.args()
    } else if e.is_zero() {
        vec![]
    } else {
        vec![e.clone()]
    }
}

/// `add` with order terms; `None` when no argument is an order term.
pub fn add_hook(args: &[Expr]) -> Option<Expr> {
    // cheap check first: the mask bit may also come from an ordinary symbol
    let has = args.iter().any(|a| match a.kind() {
        Kind::Func(Func::Order, _) => true,
        Kind::Add(_, terms) => a.may_have_order() && terms.iter().any(|(t, _)| is_order(t)),
        _ => false,
    });
    if !has {
        return None;
    }
    let mut orders: Vec<Expr> = Vec::new();
    let mut rest: Vec<Expr> = Vec::with_capacity(args.len());
    for a in args {
        if is_order(a) {
            orders.push(a.clone());
        } else if a.is_add() && a.may_have_order() {
            for t in a.args() {
                if is_order(&t) {
                    orders.push(t);
                } else {
                    rest.push(t);
                }
            }
        } else {
            rest.push(a.clone());
        }
    }
    if orders.is_empty() {
        return None;
    }
    let mut kept: Vec<Expr> = Vec::with_capacity(orders.len());
    for o in orders {
        // Add.flatten: an order of a zero expression (`O(w).subs(w, 0)`)
        // is dropped
        if order_parts(&o).is_some_and(|(oe, _)| oe.is_zero()) {
            continue;
        }
        if kept.iter().any(|k| *k == o || contains(k, &o)) {
            continue;
        }
        kept.retain(|k| !contains(&o, k));
        kept.push(o);
    }
    let s = add(&rest);
    if s.is_constant(Constant::NaN) {
        return Some(s);
    }
    let mut coeff = Number::ZERO;
    let mut pairs: Vec<(Expr, Number)> = Vec::new();
    for t in add_args(&s) {
        if kept.iter().any(|o| contains(o, &t)) {
            continue;
        }
        let (c, r) = as_coeff_mul(&t);
        if r.is_one() {
            coeff = coeff.add(&c);
        } else {
            pairs.push((r, c));
        }
    }
    pairs.extend(kept.into_iter().map(|o| (o, Number::ONE)));
    Some(add_from_terms(coeff, pairs))
}

/// `mul` with order factors; `None` when no argument is an order term.
pub fn mul_hook(args: &[Expr]) -> Option<Expr> {
    if !args.iter().any(is_order) {
        return None;
    }
    let mut rest: Vec<Expr> = Vec::new();
    let mut orders: Vec<(Expr, Vec<(Expr, Expr)>)> = Vec::new();
    for a in args {
        match order_parts(a) {
            Some(p) => orders.push(p),
            None => rest.push(a.clone()),
        }
    }
    let p = mul(&rest);
    if p.is_zero() || p.is_constant(Constant::NaN) {
        return Some(p);
    }
    let mut e = p;
    let mut vars: Vec<(Expr, Expr)> = Vec::new();
    for (oe, ov) in orders {
        e = mul2(&e, &oe);
        for v in ov {
            if !vars.iter().any(|w| w.0 == v.0) {
                vars.push(v);
            }
        }
    }
    Some(order(&e, &vars))
}

/// `O(e)**k = O(e**k)` for a positive number `k`.
pub fn pow_hook(b: &Expr, k: &Expr) -> Option<Expr> {
    // Order._eval_power: `O(e)**O(1) = O(e)`
    if order_parts(k).is_some_and(|(ke, kv)| ke.is_one() && kv.is_empty()) {
        return Some(b.clone());
    }
    let n = k.as_number()?;
    if !n.is_positive() {
        return None;
    }
    let (oe, vars) = order_parts(b)?;
    Some(order(&pow(&oe, k), &vars))
}

/// Drop the order terms of a sum (`removeO`).
pub fn remove_o(e: &Expr) -> Expr {
    if is_order(e) {
        return Expr::zero();
    }
    if e.is_add() && e.may_have_order() {
        let v: Vec<Expr> = e.args().into_iter().filter(|t| !is_order(t)).collect();
        return add(&v);
    }
    e.clone()
}

/// The order term of a sum (`getO`).
pub fn get_o(e: &Expr) -> Option<Expr> {
    if is_order(e) {
        return Some(e.clone());
    }
    if e.is_add() && e.may_have_order() {
        let v: Vec<Expr> = e.args().into_iter().filter(is_order).collect();
        if !v.is_empty() {
            return Some(add(&v));
        }
    }
    None
}

/// Exponent `n` of `O(x**n)` (SymPy `getn`), `None` if there is no
/// (univariate power) order term.
pub fn getn(e: &Expr) -> Option<Expr> {
    let o = get_o(e)?;
    let (oe, vars) = order_parts(&o)?;
    if vars.len() != 1 {
        return None;
    }
    let (x, p) = &vars[0];
    let e0 = to_zero(&oe, x, p);
    if e0.is_one() {
        return Some(Expr::zero());
    }
    if e0 == *x {
        return Some(Expr::one());
    }
    match e0.kind() {
        Kind::Pow(b, k) if b == x => Some(k.clone()),
        Kind::Mul(..) => monomial_x(&e0, x).map(|(a, _)| Expr::num(a)),
        _ => None,
    }
}
