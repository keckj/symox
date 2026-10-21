//! The recursive Mellin transform step of the reference `_rewrite_single`.
//!
//! A factor with no table entry is rewritten through its Mellin transform:
//! `F(s) = int_0^oo x**(s-1) f(x) dx` is computed with Meijer G integration
//! of products of two G-functions only (`_meijerint_definite_4` with
//! `only_double`), and the inverse transform of `F` is read as a sum of
//! G-functions of `x`. When the transform does not converge, `f(a*x)` is
//! tried with a dummy `a` set back to 1 (analytic continuation).
//!
//! The reference lets the errors other than `IntegralTransformError` of the
//! inverse transform retry propagate out of the integration: they are kept
//! in the pending error slot (`super::set_pending`), which
//! `meijerint_indefinite` raises; the definite paths still read them as
//! `None`.

use std::cell::Cell;

use super::conds::my_unpolarify;
use super::inversion::is_analytic;
use super::{definite_4_opt, get_coeff_exp, Terms, G};
use crate::expr::Constant;
use crate::ops::{mul, pow};
use crate::special::polar::{polarify, unpolarify};
use crate::special::SpecialFn;
use crate::{Expr, Func, Kind};

/// Nesting bound of the recursion (the reference is bounded by the Python
/// recursion limit only; the factors split at each level, so real inputs
/// stay far below it).
const MAX_DEPTH: u32 = 8;

thread_local! {
    static DEPTH: Cell<u32> = const { Cell::new(0) };
    /// `_dummy_('a', 'rewrite-single')`: one dummy for every call, so a
    /// continuation is never nested in another.
    static A: Expr = Expr::dummy("a", crate::symbol::Assumptions::NONE);
}

struct Guard;

impl Guard {
    fn enter() -> Option<Guard> {
        DEPTH.with(|d| {
            if d.get() >= MAX_DEPTH {
                return None;
            }
            d.set(d.get() + 1);
            Some(Guard)
        })
    }
}

impl Drop for Guard {
    fn drop(&mut self) {
        DEPTH.with(|d| d.set(d.get() - 1));
    }
}

/// `my_integrator`: the integral over `(0, oo)` by two G-functions, as a
/// Piecewise with the unevaluated integral as the otherwise branch.
fn integrator(f: &Expr, x: &Expr) -> Expr {
    let held = crate::integrate::integral_definite_unevaluated(f, x, &Expr::zero(), &Expr::oo());
    let Some((res, cond)) = definite_4_opt(f, x, true) else { return held };
    let res = crate::special::hyperexpand_options(&res, false, Some("nonrepsmall"), None).unwrap_or(res);
    crate::relational::piecewise(vec![(my_unpolarify(&res), cond), (held, Expr::boolean(true))])
}

/// `mellin_transform(...)` followed by `my_imt(...)`.
fn mellin_inverse(f: &Expr, x: &Expr, s: &Expr) -> Option<Expr> {
    let (big_f, (a, b), _) = crate::transforms::mellin_needeval(f, x, s, &integrator)?;
    match crate::transforms::inverse_mellin_meijerg(&big_f, s, x, &a, &b) {
        Ok(g) => g,
        Err((class, message)) => {
            super::set_pending(class, message);
            None
        }
    }
}

fn has_bad_value(e: &Expr) -> bool {
    let mut bad = false;
    crate::traverse::preorder(e, &mut |t| {
        if matches!(t.kind(), Kind::Const(Constant::Infinity | Constant::ComplexInfinity | Constant::NaN)) {
            bad = true;
        }
    });
    bad
}

/// The recursive part of `_rewrite_single(f, x, recursive=True)`.
pub(super) fn rewrite(f: &Expr, x: &Expr) -> Option<(Terms, Expr)> {
    let _guard = Guard::enter()?;
    let s = Expr::dummy("s", crate::symbol::Assumptions::NONE);
    let mut g = mellin_inverse(f, x, &s);
    if g.is_none() && !super::has_pending() {
        let a = A.with(|a| a.clone());
        if !crate::traverse::has(f, &a) && is_analytic(f, x) {
            let scaled = crate::subs::subs1(f, x, &mul(&[a.clone(), x.clone()]));
            g = mellin_inverse(&scaled, x, &s).map(|g| crate::subs::subs1(&g, &a, &Expr::one()));
        }
    }
    let g = g?;
    if has_bad_value(&g) {
        return None;
    }
    let mut res = Vec::new();
    for term in crate::simptools::add_args(&g) {
        let (c, m) = crate::decompose::as_coeff_mul(&term, std::slice::from_ref(x), false);
        let gm = match m.as_slice() {
            [gm] => gm,
            [] => {
                super::set_pending("IndexError", "list index out of range".into());
                return None;
            }
            _ => {
                super::set_pending("NotImplementedError", "Unexpected form...".into());
                return None;
            }
        };
        let Kind::Func(Func::Special(SpecialFn::MeijerG), _) = gm.kind() else { return None };
        let gf = G::from_expr(gm)?;
        let Some((a, b)) = get_coeff_exp(&gf.z, x) else {
            super::set_pending("ValueError", super::coeff_exp_error(&gf.z, x));
            return None;
        };
        let z = mul(&[unpolarify(&polarify(&a, true), true), pow(x, &b)]);
        match G::try_new(gf.an, gf.ao, gf.bm, gf.bo, z) {
            Ok(g) => res.push((c, Expr::zero(), g)),
            Err(message) => {
                super::set_pending("ValueError", message);
                return None;
            }
        }
    }
    Some((res, Expr::boolean(true)))
}
