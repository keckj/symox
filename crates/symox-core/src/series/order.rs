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
    let b = degree_in(&lt.coeff_logt, &lt.logt).unwrap_or(0);
    Some((lt.exp, b))
}

/// Canonical constructor `O(expr, (x1, p1), ...)`.
pub fn order(expr: &Expr, vars: &[(Expr, Expr)]) -> Expr {
    let mut expr = expr.clone();
    let mut vars: Vec<(Expr, Expr)> = vars.to_vec();
    if let Some((inner, ivars)) = order_parts(&expr) {
        expr = inner;
        if vars.is_empty() {
            vars = ivars;
        }
    }
    if expr.may_have_order() {
        expr = remove_o(&expr);
    }
    if vars.is_empty() {
        vars = free_symbols(&expr).into_iter().map(|s| (s, Expr::zero())).collect();
    }
    vars.sort_by(|a, b| a.0.cmp_struct(&b.0));
    vars.dedup_by(|a, b| a.0 == b.0);
    if vars.iter().any(|(v, _)| !v.is_symbol()) {
        return raw(expr, &vars);
    }
    if expr.is_zero() {
        return raw(expr, &vars);
    }
    if vars.is_empty() {
        return raw(Expr::one(), &[]);
    }
    if vars.len() == 1 {
        let (x, p) = &vars[0];
        let e0 = to_zero(&expr, x, p);
        let lead = match lead_exponents(&e0, x) {
            Some((a, b)) => build_monomial(x, &a, b),
            None => as_coeff_mul(&e0).1,
        };
        return raw(from_zero(&lead, x, p), &vars);
    }
    // several variables: drop the factors free of all of them
    let syms: Vec<Expr> = vars.iter().map(|v| v.0.clone()).collect();
    let dep = |f: &Expr| syms.iter().any(|s| depends(f, s));
    let strip = |t: &Expr| match t.kind() {
        Kind::Mul(_, facs) => mul(&facs.iter().map(|(b, p)| mul_factor(b, p)).filter(|f| dep(f)).collect::<Vec<_>>()),
        _ if !dep(t) => Expr::one(),
        _ => as_coeff_mul(t).1,
    };
    if !expr.is_add() {
        return raw(strip(&expr), &vars);
    }
    // A sum keeps its leading terms only (reference `extract_leading_order`
    // on the expanded expression, with every point moved to the origin):
    // `O(1 + x*y) = O(1)`, `O(x + y + x*y) = O(x + y)`. A term whose ratio
    // to another is not a monomial is kept.
    let mut e0 = expr.clone();
    for (x, p) in &vars {
        e0 = to_zero(&e0, x, p);
    }
    let e0 = crate::expand::expand(&e0);
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
            return tvars.is_empty() && vars.is_empty();
        }
        return contains_expr(&oe, &vars, &te);
    }
    contains_expr(&oe, &vars, t)
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
                None => false,
            };
        }
        let ratio = div(&t0, &o0);
        return match lead_exponents(&ratio, x) {
            Some((a, b)) => a.is_positive() || (a.is_zero() && b <= 0),
            None => false,
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
