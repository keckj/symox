//! `_laplace_transform_integration`: the Laplace transform as the definite
//! integral `int_0^oo f(t) exp(-s*t) dt`, used when no rule applies. The
//! convergence half plane and the auxiliary conditions are read from the
//! conditions of the integral (`process_conds`), as in the reference.

use super::util::*;
use crate::expr::{Expr, Func, Kind, RelOp};
use crate::functions::{abs, exp, max, min, re};
use crate::logic::LogicOp;
use crate::ops::{add, div, mul, neg, sub};
use crate::special::SpecialFn;

type R = (Expr, Expr, Expr);

fn has_func(e: &Expr, pred: &dyn Fn(&Func) -> bool) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(g, _) = s.kind() {
            found |= pred(g);
        }
    });
    found
}

pub(super) fn has_integral(e: &Expr) -> bool {
    has_func(e, &|g| matches!(g, Func::Integral))
}

/// `ft.atoms(AppliedUndef)` contains a function of `t`.
pub(crate) fn has_applied_undef_in(e: &Expr, t: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(Func::Undefined(_), _) = s.kind() {
            found |= has(s, t);
        }
    });
    found
}

pub(super) fn logic_args(e: &Expr, op: LogicOp) -> Vec<Expr> {
    match e.kind() {
        Kind::Func(Func::Logic(o), a) if *o == op => a.clone(),
        _ => vec![e.clone()],
    }
}

fn rel(op: RelOp, a: &Expr, b: &Expr) -> Expr {
    crate::relational::rel(op, a, b)
}

fn is_bool(e: &Expr, v: bool) -> bool {
    matches!(e.kind(), Kind::Bool(b) if *b == v)
}

/// `Relational.lts`: the lesser side.
pub(super) fn lts(op: RelOp, l: &Expr, r: &Expr) -> Expr {
    match op {
        RelOp::Gt | RelOp::Ge => r.clone(),
        _ => l.clone(),
    }
}

/// `cond.xreplace({r: r.canonical for r in cond.atoms(Relational)})`.
fn canonical_all(c: &Expr) -> Expr {
    match c.kind() {
        Kind::Rel(..) => crate::logic::canonical(c),
        Kind::Func(Func::Logic(op), args) => {
            let new: Vec<Expr> = args.iter().map(canonical_all).collect();
            if new == *args {
                c.clone()
            } else {
                crate::logic::eval(*op, new)
            }
        }
        _ => c.clone(),
    }
}

/// `ie.subs(s, i)` classified as true, false or undecided (`None`).
fn classify(ie: &Expr, t: &Expr, i: &Expr) -> Option<bool> {
    let v = crate::subs::subs(ie, &[(t.clone(), i.clone())]);
    match v.kind() {
        Kind::Bool(b) => Some(*b),
        _ => None,
    }
}

/// `_solve_inequality(ie, t)` for the relationals met here: `t` isolated
/// on the left when the relation is linear with a coefficient of known
/// sign, otherwise the univariate inequality solver.
pub(super) fn solve_inequality(ie: &Expr, t: &Expr) -> Expr {
    let Kind::Rel(op, l, r) = ie.kind() else { return ie.clone() };
    if !has(ie, t) {
        return ie.clone();
    }
    let (mut op, mut l, mut r) = (*op, l.clone(), r.clone());
    if r == *t {
        op = crate::logic::reversed_op(op);
        std::mem::swap(&mut l, &mut r);
    }
    if l == *t && !has(&r, t) {
        return rel(op, &l, &r);
    }
    let expr = sub(&l, &r);
    let degree = poly(&expr, t).map(|c| c.len() as i64 - 1);
    match degree {
        Some(0) => return rel(op, &crate::expand::expand(&expr), &Expr::zero()),
        Some(1) => {}
        _ => {
            return crate::sets::solve_univariate_inequality(ie, t, true, &crate::sets::reals())
                .unwrap_or_else(|_| ie.clone());
        }
    }
    let e = crate::expand::expand(&expr);
    let (b, e1) = split_add(&e, t);
    let rhs = neg(&b);
    let ef = crate::simptools::exprtools::factor_terms(&e1, false, false, false, true);
    let (mut a, mut lhs) = split_mul(&ef, t);
    let zero_unknown = crate::assumptions::is_zero(&a) != Some(false);
    let sign_unknown = crate::assumptions::is_negative(&a).is_none() && crate::assumptions::is_positive(&a).is_none();
    if zero_unknown || (sign_unknown && !matches!(op, RelOp::Eq | RelOp::Ne)) {
        lhs = ef;
        a = Expr::one();
    }
    let rhs = div(&rhs, &a);
    let rv = if crate::assumptions::is_positive(&a) == Some(true) {
        rel(op, &lhs, &rhs)
    } else {
        rel(crate::logic::reversed_op(op), &lhs, &rhs)
    };
    let mut conds = Vec::new();
    for (i, upper) in [(Expr::neg_oo(), false), (Expr::oo(), true)] {
        if classify(&rv, t, &i) == Some(true) && classify(ie, t, &i) != Some(true) {
            conds.push(if upper { rel(RelOp::Lt, t, &i) } else { rel(RelOp::Lt, &i, t) });
        }
    }
    conds.push(rv);
    crate::logic::and(conds)
}

/// `e.as_independent(t, as_Add=True)`.
fn split_add(e: &Expr, t: &Expr) -> (Expr, Expr) {
    let (mut c, mut v) = (Vec::new(), Vec::new());
    for term in terms(e) {
        if has(&term, t) {
            v.push(term)
        } else {
            c.push(term)
        }
    }
    (add(&c), add(&v))
}

/// The Wild match of `(s + w3)` (exponent and polar lift handled by the
/// caller): returns `w3`.
fn shift_of(b: &Expr, s: &Expr) -> Option<Expr> {
    let w3 = sub(b, s);
    (!has(&w3, s)).then_some(w3)
}

/// `X = B*q` with `q` free of `s`: `(q, B)`.
fn split_q(x: &Expr, s: &Expr) -> (Expr, Expr) {
    split_mul(x, s)
}

/// The six argument patterns of `process_conds`. Returns `(p, q, w2, w3)`.
fn match_arg_patterns(lhs: &Expr, rhs: &Expr, s: &Expr) -> Option<(Expr, Expr, Expr, Expr)> {
    if has(rhs, s) {
        return None;
    }
    let (p_out, core) = split_mul(lhs, s);
    let Kind::Func(Func::Abs, a) = core.kind() else { return None };
    match a[0].kind() {
        // p*Abs(arg((s + w3)*q))
        Kind::Func(Func::Arg, x) => {
            let (q, b) = split_q(&x[0], s);
            let w3 = shift_of(&b, s)?;
            Some((p_out, q, rhs.clone(), w3))
        }
        // Abs(periodic_argument((s + w3)**p*q, w1)), also with polar_lift
        Kind::Func(Func::Special(SpecialFn::PeriodicArgument), x) if p_out.is_one() => {
            if has(&x[1], s) {
                return None;
            }
            let (q, d) = split_q(&x[0], s);
            let (base, p) = match d.kind() {
                Kind::Pow(b, e) if !has(e, s) => (b.clone(), e.clone()),
                _ => (d.clone(), Expr::one()),
            };
            let base = match base.kind() {
                Kind::Func(Func::Special(SpecialFn::PolarLift), y) => y[0].clone(),
                _ => base,
            };
            let w3 = shift_of(&base, s)?;
            Some((p, q, rhs.clone(), w3))
        }
        _ => None,
    }
}

/// `x.expand().as_real_imag()[0]` for every `re(x)`.
pub(super) fn expand_re(e: &Expr) -> Expr {
    let e = crate::traverse::map_children(e, &mut |c| expand_re(c));
    if let Kind::Func(Func::Re, a) = e.kind() {
        let x = crate::expand::expand(&a[0]);
        if let Ok(Some((r, _))) = crate::complex_parts::as_real_imag(&x, true, None) {
            return r;
        }
    }
    e
}

/// `process_conds`: a half plane `Re(s) > a` and auxiliary conditions, or
/// `None` when a condition bounds `Re(s)` from above.
fn process_conds(conds: &Expr, s: &Expr, t: &Expr) -> Option<(Expr, Expr)> {
    let mut a = Expr::neg_oo();
    let mut aux = tru();
    let cnf = crate::logic::to_cnf(conds, false, false).unwrap_or_else(|_| conds.clone());
    for c in logic_args(&cnf, LogicOp::And) {
        let mut a_ = Expr::oo();
        let mut aux_: Vec<Expr> = Vec::new();
        for d0 in logic_args(&c, LogicOp::Or) {
            let mut d = d0.clone();
            if let Kind::Rel(op, l, r) = d.kind() {
                if has(r, s) {
                    d = Expr::from_kind(Kind::Rel(crate::logic::reversed_op(*op), r.clone(), l.clone()));
                }
            }
            if let Kind::Rel(op @ (RelOp::Gt | RelOp::Ge), l, r) = d.kind() {
                let op = if *op == RelOp::Gt { RelOp::Lt } else { RelOp::Le };
                d = rel(op, &neg(l), &neg(r));
            }
            if let Kind::Rel(RelOp::Lt | RelOp::Le, l, r) = d.kind() {
                if let Some((p, q, w2, w3)) = match_arg_patterns(l, r, s) {
                    let half = div(&w2, &p) == div(&Expr::pi(), &Expr::int(2));
                    if crate::assumptions::is_positive(&q) == Some(true) && half {
                        d = rel(RelOp::Lt, &neg(&re(&add(&[s.clone(), w3]))), &Expr::zero());
                    }
                }
            }
            let d_ = crate::subs::subs1(&expand_re(&d), &re(s), t);
            let relational = matches!(d.kind(), Kind::Rel(op, ..) if !matches!(op, RelOp::Eq | RelOp::Ne));
            if !relational || has(&d_, s) || !has(&d_, t) {
                aux_.push(d);
                continue;
            }
            let soln = solve_inequality(&d_, t);
            match soln.kind() {
                Kind::Rel(op, l, r) if !matches!(op, RelOp::Eq | RelOp::Ne) => {
                    let low = lts(*op, l, r);
                    if low == *t {
                        return None;
                    }
                    a_ = min(vec![low, a_]);
                }
                _ => {
                    aux_.push(d);
                    continue;
                }
            }
        }
        if a_ != Expr::oo() {
            a = max(vec![a_, a]);
        } else {
            aux = crate::logic::and(vec![aux, crate::logic::or(aux_)]);
        }
    }
    let aux = if matches!(aux.kind(), Kind::Rel(..)) { crate::logic::canonical(&aux) } else { aux };
    Some((a, aux))
}

/// `_simplifyconds(expr, s, a)`, given `Re(s) > a`.
fn simplify_conds(expr: &Expr, s: &Expr, a: &Expr) -> Expr {
    fn power(ex: &Expr, s: &Expr) -> Option<Expr> {
        if ex == s {
            return Some(Expr::one());
        }
        match ex.kind() {
            Kind::Pow(b, e) if b == s => Some(e.clone()),
            _ => None,
        }
    }
    fn strip_abs(e: &Expr) -> Expr {
        match e.kind() {
            Kind::Func(Func::Abs, a) => a[0].clone(),
            _ => e.clone(),
        }
    }
    fn bigger(ex1: &Expr, ex2: &Expr, s: &Expr, a: &Expr) -> Option<bool> {
        if has(ex1, s) && has(ex2, s) {
            return None;
        }
        let (ex1, ex2) = (strip_abs(ex1), strip_abs(ex2));
        if has(&ex1, s) {
            return bigger(&crate::ops::inv(&ex2), &crate::ops::inv(&ex1), s, a);
        }
        let n = power(&ex2, s)?;
        let bound = crate::ops::pow(&abs(a), &n);
        let is_pos = crate::assumptions::is_positive(&n);
        let is_neg = crate::assumptions::is_negative(&n);
        // `n > 0` raises TypeError (caught) for an undecided sign.
        if is_pos.is_none() {
            return None;
        }
        if is_pos == Some(true) && is_bool(&rel(RelOp::Le, &abs(&ex1), &bound), true) {
            return Some(false);
        }
        if is_neg.is_none() {
            return None;
        }
        if is_neg == Some(true) && is_bool(&rel(RelOp::Ge, &abs(&ex1), &bound), true) {
            return Some(true);
        }
        None
    }
    let pos_or_abs = |x: &Expr| crate::assumptions::is_positive(x) == Some(true) || matches!(x.kind(), Kind::Func(Func::Abs, _));
    let replie = |x: &Expr, y: &Expr| -> Expr {
        if !pos_or_abs(x) || !pos_or_abs(y) {
            return rel(RelOp::Lt, x, y);
        }
        match bigger(x, y, s, a) {
            Some(r) => Expr::boolean(!r),
            None => rel(RelOp::Lt, x, y),
        }
    };
    fn walk(e: &Expr, f: &dyn Fn(&Expr) -> Option<Expr>) -> Expr {
        let e = crate::traverse::map_children(e, &mut |c| walk(c, f));
        f(&e).unwrap_or(e)
    }
    if let Kind::Bool(_) = expr.kind() {
        return expr.clone();
    }
    let e = walk(expr, &|e| match e.kind() {
        Kind::Rel(RelOp::Lt, x, y) => Some(replie(x, y)),
        _ => None,
    });
    let e = walk(&e, &|e| match e.kind() {
        Kind::Rel(RelOp::Gt, x, y) => Some(replie(y, x)),
        _ => None,
    });
    walk(&e, &|e| match e.kind() {
        Kind::Rel(RelOp::Ne, x, y) => Some(match bigger(x, y, s, a) {
            Some(_) => tru(),
            None => rel(RelOp::Ne, x, y),
        }),
        _ => None,
    })
}

/// The transforms' `_simplify(expr, doit)`: `simplify(powdenest(expr,
/// polar=True))`, including the gamma simplification step of the reference
/// `simplify` (`combsimp` passes gamma functions of non-integer arguments to
/// `gammasimp`).
pub(super) fn simplify_result(e: &Expr, simplify: bool) -> Expr {
    if !simplify {
        return e.clone();
    }
    let p = crate::integrate::meijerint::powdenest_polar(e);
    let r = crate::simplify::simplify(&p);
    let mut gamma = false;
    crate::traverse::preorder(&r, &mut |n| {
        if let Kind::Func(Func::Gamma, a) = n.kind() {
            gamma |= crate::assumptions::is_integer(&a[0]) != Some(true);
        }
    });
    if !gamma {
        return r;
    }
    let g = crate::simplify::simplify(&crate::simptools::gammasimp::gammasimp(&r));
    let ops = |x: &Expr| crate::simptools::methods::count_ops(x);
    if ops(&g) <= ops(&r) { g } else { r }
}

/// `_laplace_transform_integration(f, t, s_)`.
pub(crate) fn laplace_integration(f: &Expr, t: &Expr, s_: &Expr, simplify: bool) -> Option<R> {
    if has_func(f, &|g| matches!(g, Func::DiracDelta)) {
        return None;
    }
    let s = Expr::dummy("s", Default::default());
    let integrand = mul(&[f.clone(), exp(&neg(&mul(&[s.clone(), t.clone()])))]);
    let big_f = crate::integrate::integrate_definite(&integrand, t, &Expr::zero(), &Expr::oo());
    let back = |e: &Expr| crate::subs::subs(e, &[(s.clone(), s_.clone())]);
    if !has_integral(&big_f) {
        return Some((simplify_result(&back(&big_f), simplify), Expr::neg_oo(), tru()));
    }
    let Kind::Piecewise(pieces) = big_f.kind() else { return None };
    let (value, cond) = pieces[0].clone();
    if has_integral(&value) {
        return None;
    }
    let mut conds = Vec::new();
    for c in logic_args(&cond, LogicOp::Or) {
        // The reference subscripts the `None` of an upper bound and raises
        // a TypeError; the transform is left unevaluated instead.
        conds.push(process_conds(&c, &s, t)?);
    }
    let mut conds2: Vec<(Expr, Expr)> = conds
        .iter()
        .filter(|(a, aux)| !is_bool(aux, false) && *a != Expr::neg_oo())
        .cloned()
        .collect();
    if conds2.is_empty() {
        conds2 = conds.iter().filter(|(_, aux)| !is_bool(aux, false)).cloned().collect();
    }
    if conds2.is_empty() {
        return None;
    }
    if conds2.len() > 1 {
        let mut keyed: Vec<(Expr, Expr)> = conds2.clone();
        keyed.sort_by(|x, y| crate::logic::ordered_cmp(&x.0, &y.0).then_with(|| crate::logic::ordered_cmp(&x.1, &y.1)));
        let plane = |c: &(Expr, Expr)| crate::evalf::eval_f64(&c.0).or(match c.0.kind() {
            Kind::Const(crate::expr::Constant::Infinity) => Some(f64::INFINITY),
            Kind::Const(crate::expr::Constant::NegInfinity) => Some(f64::NEG_INFINITY),
            _ => None,
        });
        if keyed.iter().all(|c| plane(c).is_some()) {
            let cnt = |e: &Expr| if matches!(e.kind(), Kind::Bool(_)) { 0 } else { crate::simptools::methods::count_ops(e) };
            keyed.sort_by(|x, y| {
                plane(y).unwrap().partial_cmp(&plane(x).unwrap()).unwrap_or(std::cmp::Ordering::Equal).then(cnt(&x.1).cmp(&cnt(&y.1)))
            });
        }
        conds2 = keyed;
    }
    let (a, aux) = conds2.swap_remove(0);
    let (value, aux) = if simplify { (simplify_conds(&value, &s, &a), simplify_conds(&aux, &s, &a)) } else { (value, aux) };
    Some((simplify_result(&back(&value), simplify), back(&a), canonical_all(&back(&aux))))
}
