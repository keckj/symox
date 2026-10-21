//! The reference's integration backends of the Mellin, Fourier type and
//! Hankel transforms (`_mellin_transform`, `_fourier_transform`,
//! `_sine_cosine_transform`, `_hankel_transform`): the transform is the
//! definite integral of the kernel times the function, and the convergence
//! conditions are read from the conditions of the integral. They are used
//! when the closed form tables do not apply.

use super::laplace_integration::{expand_re, has_integral, logic_args, lts, simplify_result, solve_inequality};
use super::util::*;
use crate::expr::{Expr, Kind, RelOp};
use crate::functions::{exp, max, min, re};
use crate::logic::LogicOp;
use crate::ops::{mul, pow, sub};

type V = (Expr, Expr);

/// The first piece of a piecewise result, with the reference's errors.
fn first_piece(big_f: &Expr) -> Result<V, &'static str> {
    let Kind::Piecewise(pieces) = big_f.kind() else { return Err("could not compute integral") };
    let (value, cond) = pieces[0].clone();
    if has_integral(&value) {
        return Err("integral in unexpected form");
    }
    Ok((value, cond))
}

/// `_fourier_transform(f, x, k, a, b)`: `a * int_{-oo}^{oo} f(x) exp(b*I*x*k) dx`.
pub(crate) fn fourier_integration(f: &Expr, x: &Expr, k: &Expr, a: &Expr, b: &Expr, simplify: bool) -> Result<V, &'static str> {
    let kernel = exp(&mul(&[b.clone(), Expr::i(), x.clone(), k.clone()]));
    let big_f = crate::integrate::integrate_definite(&mul(&[a.clone(), f.clone(), kernel]), x, &Expr::neg_oo(), &Expr::oo());
    if !has_integral(&big_f) {
        return Ok((simplify_result(&big_f, simplify), tru()));
    }
    let integral_f = crate::integrate::integrate_definite(f, x, &Expr::neg_oo(), &Expr::oo());
    if integral_f == Expr::oo() || integral_f == Expr::neg_oo() || integral_f == Expr::nan() || has_integral(&integral_f) {
        return Err("function not integrable on real axis");
    }
    let (value, cond) = first_piece(&big_f)?;
    Ok((simplify_result(&value, simplify), cond))
}

/// `_sine_cosine_transform(f, x, k, a, b, K)`: `a * int_0^oo f(x) K(b*x*k) dx`.
pub(crate) fn sine_cosine_integration(f: &Expr, x: &Expr, k: &Expr, a: &Expr, b: &Expr, sine: bool, simplify: bool) -> Result<V, &'static str> {
    let arg = mul(&[b.clone(), x.clone(), k.clone()]);
    let kernel = if sine { crate::functions::sin(&arg) } else { crate::functions::cos(&arg) };
    let big_f = crate::integrate::integrate_definite(&mul(&[a.clone(), f.clone(), kernel]), x, &Expr::zero(), &Expr::oo());
    if !has_integral(&big_f) {
        return Ok((simplify_result(&big_f, simplify), tru()));
    }
    let (value, cond) = first_piece(&big_f)?;
    Ok((simplify_result(&value, simplify), cond))
}

/// `_hankel_transform(f, r, k, nu)`: `int_0^oo f(r) J_nu(k*r) r dr`.
pub(crate) fn hankel_integration(f: &Expr, r: &Expr, k: &Expr, nu: &Expr, simplify: bool) -> Result<V, &'static str> {
    let integrand = mul(&[f.clone(), besselj(nu, &mul(&[k.clone(), r.clone()])), r.clone()]);
    let big_f = crate::integrate::integrate_definite(&integrand, r, &Expr::zero(), &Expr::oo());
    if !has_integral(&big_f) {
        return Ok((simplify_result(&big_f, simplify), tru()));
    }
    let (value, cond) = first_piece(&big_f)?;
    Ok((simplify_result(&value, simplify), cond))
}

/// `process_conds` of `_mellin_transform`: the strip `(a, b)` and the
/// auxiliary conditions of one disjunct.
fn mellin_process_conds(cond: &Expr, s: &Expr, t: &Expr) -> (Expr, Expr, Expr) {
    let mut a = Expr::neg_oo();
    let mut b = Expr::oo();
    let mut aux = tru();
    let cnf = crate::logic::to_cnf(cond, false, false).unwrap_or_else(|_| cond.clone());
    for c in logic_args(&cnf, LogicOp::And) {
        let mut a_ = Expr::oo();
        let mut b_ = Expr::neg_oo();
        let mut aux_: Vec<Expr> = Vec::new();
        for d in logic_args(&c, LogicOp::Or) {
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
                        let high = if matches!(op, RelOp::Gt | RelOp::Ge) { l.clone() } else { r.clone() };
                        b_ = max(vec![high, b_]);
                    } else {
                        a_ = min(vec![low, a_]);
                    }
                }
                _ => aux_.push(d),
            }
        }
        if a_ != Expr::oo() && a_ != b {
            a = max(vec![a_, a]);
        } else if b_ != Expr::neg_oo() && b_ != a {
            b = min(vec![b_, b]);
        } else {
            aux = crate::logic::and(vec![aux, crate::logic::or(aux_)]);
        }
    }
    (a, b, aux)
}

/// `_mellin_transform(f, x, s_)`: `int_0^oo x**(s-1) f(x) dx` with its
/// fundamental strip.
pub(crate) fn mellin_integration(f: &Expr, x: &Expr, s_: &Expr, simplify: bool) -> Result<(Expr, (Expr, Expr), Expr), &'static str> {
    mellin_integration_with(f, x, s_, simplify, &|g, x| crate::integrate::integrate_definite(g, x, &Expr::zero(), &Expr::oo()))
}

/// `_mellin_transform(f, x, s_, integrator)`: `integrator(g, x)` is the
/// integral of `g` over `(0, oo)`.
pub(crate) fn mellin_integration_with(
    f: &Expr,
    x: &Expr,
    s_: &Expr,
    simplify: bool,
    integrator: &dyn Fn(&Expr, &Expr) -> Expr,
) -> Result<(Expr, (Expr, Expr), Expr), &'static str> {
    let s = Expr::dummy("s", Default::default());
    let integrand = mul(&[pow(x, &sub(&s, &Expr::one())), f.clone()]);
    let big_f = integrator(&integrand, x);
    let back = |e: &Expr| crate::subs::subs(e, &[(s.clone(), s_.clone())]);
    if !has_integral(&big_f) {
        return Ok((simplify_result(&back(&big_f), simplify), (Expr::neg_oo(), Expr::oo()), tru()));
    }
    let (value, cond) = first_piece(&big_f)?;
    let t = Expr::dummy("t", crate::symbol::Assumptions::REAL);
    let mut conds: Vec<(Expr, Expr, Expr)> = logic_args(&cond, LogicOp::Or)
        .iter()
        .map(|c| mellin_process_conds(c, &s, &t))
        .filter(|c| truth(&c.2) != Some(false))
        .collect();
    if conds.is_empty() {
        return Err("no convergence found");
    }
    if conds.len() > 1 {
        // key (a - b, count_ops(aux)), compared numerically when possible
        let width = |c: &(Expr, Expr, Expr)| -> Option<f64> {
            let w = sub(&c.0, &c.1);
            crate::evalf::eval_f64(&w).or(match w.kind() {
                Kind::Const(crate::expr::Constant::Infinity) => Some(f64::INFINITY),
                Kind::Const(crate::expr::Constant::NegInfinity) => Some(f64::NEG_INFINITY),
                _ => None,
            })
        };
        if conds.iter().all(|c| width(c).is_some()) {
            let cnt = |e: &Expr| if matches!(e.kind(), Kind::Bool(_)) { 0 } else { crate::simptools::methods::count_ops(e) };
            conds.sort_by(|x, y| {
                width(x).unwrap().partial_cmp(&width(y).unwrap()).unwrap_or(std::cmp::Ordering::Equal).then(cnt(&x.2).cmp(&cnt(&y.2)))
            });
        }
    }
    let (a, b, aux) = conds.swap_remove(0);
    Ok((simplify_result(&back(&value), simplify), (back(&a), back(&b)), back(&aux)))
}
