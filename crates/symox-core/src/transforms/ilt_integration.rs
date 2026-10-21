//! The last stage of the inverse Laplace transform: the reference's
//! `_inverse_laplace_transform_integration`.
//!
//! `F(s)` is first read as an inverse Mellin transform of `f(exp(-t))` on
//! the strip `(None, oo)`; when that fails the Bromwich integral is computed
//! with Meijer G-functions (`meijerint_inversion`). Heaviside factors of
//! `exp(-t)` are brought back to `Heaviside(t - c)` and exponentials are
//! split into real and imaginary parts.
//!
//! The reference raises the errors of the inverse Mellin rewrite (a pole
//! order that cannot be decided, gammas across the strip); they leave the
//! transform unevaluated here.

use super::util::{has, split_mul, terms};
use super::SpecialFn;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions as fx;
use crate::ops::{add, div, mul, neg};

type R = (Expr, Expr);

fn exp_arg(e: &Expr) -> Option<Expr> {
    match e.kind() {
        Kind::Func(Func::Exp, a) => Some(a[0].clone()),
        Kind::Pow(base, a) if base.is_constant(Constant::E) => Some(a.clone()),
        _ => None,
    }
}

fn is_rational_function(e: &Expr, s: &Expr) -> bool {
    if !has(e, s) {
        return true;
    }
    match e.kind() {
        Kind::Sym(_) => true,
        Kind::Add(..) | Kind::Mul(..) => e.args().iter().all(|a| is_rational_function(a, s)),
        Kind::Pow(b, n) => n.is_integer_number() && is_rational_function(b, s),
        _ => false,
    }
}

/// `_simplify(expr, doit)` of the transforms module.
fn simplify_result(e: &Expr, doit: bool) -> Expr {
    if !doit {
        return e.clone();
    }
    crate::simplify::simplify(&crate::integrate::meijerint::powdenest_polar(e))
}

/// `pw_simp`: a three piece Piecewise from hyperexpand as Heaviside terms.
fn pw_simp(e: &Expr, t: &Expr) -> Expr {
    let Kind::Piecewise(pieces) = e.kind() else { return e.clone() };
    if pieces.len() != 3 {
        return e.clone();
    }
    let Kind::Func(Func::Special(SpecialFn::MeijerG), margs) = pieces[2].0.kind() else { return e.clone() };
    let Some((coeff, exponent)) = crate::integrate::meijerint::get_coeff_exp(&margs[2], t) else { return e.clone() };
    let te = crate::ops::pow(t, &exponent);
    let c = div(&Expr::one(), &fx::abs(&coeff));
    add(&[
        mul(&[fx::heaviside(&crate::ops::sub(&c, &te)), pieces[0].0.clone()]),
        mul(&[fx::heaviside(&crate::ops::sub(&te, &c)), pieces[1].0.clone()]),
    ])
}

/// `simp_heaviside`: `Heaviside(g(exp(-t)))` as `Heaviside(+-(t + k))`.
fn simp_heaviside(h: &Expr, t: &Expr, u: &Expr) -> Expr {
    let Kind::Func(Func::Heaviside, hargs) = h.kind() else { return h.clone() };
    let arg = &hargs[0];
    let a = crate::subs::subs1(arg, &fx::exp(&neg(t)), u);
    if has(&a, t) {
        return h.clone();
    }
    let rel = super::laplace_integration::solve_inequality(&crate::relational::gt(&a, &Expr::zero()), u);
    let Kind::Rel(op, l, r) = rel.kind() else { return h.clone() };
    // `rel.lts` and `rel.gts`
    let (lts, gts) = match op {
        crate::expr::RelOp::Lt | crate::expr::RelOp::Le => (l.clone(), r.clone()),
        crate::expr::RelOp::Gt | crate::expr::RelOp::Ge => (r.clone(), l.clone()),
        _ => return h.clone(),
    };
    let rebuild = |x: Expr| {
        let mut v = hargs.clone();
        v[0] = x;
        crate::functions::func(Func::Heaviside, v)
    };
    if &lts == u {
        rebuild(add(&[t.clone(), fx::log(&gts)]))
    } else {
        rebuild(neg(&add(&[t.clone(), fx::log(&lts)])))
    }
}

/// `expand_complex(exp(arg))`.
fn simp_exp(e: &Expr) -> Expr {
    let only = crate::expand::ExpandOptions {
        basic: false,
        mul: false,
        multinomial: false,
        power_base: false,
        power_exp: false,
        log: false,
        ..Default::default()
    };
    crate::complex_parts::expand_complex(e, &only).unwrap_or_else(|_| e.clone())
}

fn has_integral(e: &Expr) -> bool {
    super::util::has_head(e, &|k| matches!(k, Kind::Func(Func::Integral, _)))
}

/// `InverseMellinTransform(F, s, x, None, oo).doit(needeval=True)`.
fn inverse_mellin(f: &Expr, s: &Expr, x: &Expr) -> Option<R> {
    let args = [f.clone(), s.clone(), x.clone(), super::none_plane(), Expr::oo()];
    let v = super::transform_doit(SpecialFn::InverseMellinTransform, &args, Some(true), true).ok()?;
    match v.as_slice() {
        [f, c] => Some((f.clone(), c.clone())),
        _ => None,
    }
}

/// `_inverse_laplace_transform_integration(F, s, t, plane, simplify)`.
pub(super) fn integration(f: &Expr, s: &Expr, t_: &Expr, simplify: bool) -> Option<R> {
    let t = Expr::dummy("t", crate::symbol::Assumptions::REAL);
    let f = if is_rational_function(f, s) { crate::poly::api::apart(f, Some(s)).unwrap_or_else(|_| f.clone()) } else { f.clone() };
    if f.is_add() {
        // the reference passes `simplify` positionally here and raises
        // TypeError; this is the sum it means
        let mut out = Vec::new();
        let mut conds = Vec::new();
        for x in terms(&f) {
            let (k, g) = split_mul(&x, s);
            let (r, c) = integration(&g, s, &t, simplify)?;
            out.push(mul(&[k, r]));
            conds.push(c);
        }
        let res = crate::subs::subs1(&add(&out), &t, t_);
        return Some((simplify_result(&res, simplify), crate::logic::and(conds)));
    }
    let (mut res, cond) = match inverse_mellin(&f, s, &fx::exp(&neg(&t))) {
        Some(r) => r,
        None => {
            let r = crate::integrate::meijerint::meijerint_inversion(&f, s, &t)?;
            let (r, cond) = match r.kind() {
                Kind::Piecewise(pieces) => {
                    let (r, c) = pieces[0].clone();
                    if has_integral(&r) {
                        return None;
                    }
                    (r, c)
                }
                _ => (r.clone(), Expr::boolean(true)),
            };
            (crate::simptools::bottom_up(&r, &mut |e| pw_simp(e, &t)), cond)
        }
    };
    if matches!(res.kind(), Kind::Piecewise(_)) {
        return Some((crate::subs::subs1(&res, &t, t_), cond));
    }
    let u = Expr::dummy("u", crate::symbol::Assumptions::NONE);
    res = crate::simptools::bottom_up(&res, &mut |e| simp_heaviside(e, &t, &u));
    res = crate::simptools::bottom_up(&res, &mut |e| if exp_arg(e).is_some() { simp_exp(e) } else { e.clone() });
    Some((simplify_result(&crate::subs::subs1(&res, &t, t_), simplify), cond))
}
