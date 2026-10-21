//! Indefinite integration by Meijer G-functions: the reference
//! `meijerint_indefinite` and `_meijerint_indefinite_1`.
//!
//! The integrand is rewritten as `fac*x**c*sum(C*G(a*x**b))` (with the
//! recursive Mellin rewrite); each term `t**rho*G(t)` with `t = a*x**b` is
//! integrated by adding the parameters `1 - rho` and `-rho`, and the result
//! is expanded with `hyperexpand`. The condition of the rewrite guards the
//! value: `Piecewise((F, cond), (Integral(f, x), True))`.
//!
//! An error inside a term (`hyperexpand` failing, a non-admissible
//! G-function) aborts the method like the reference exception; the
//! reference lets those other than `NotImplementedError` propagate out of
//! `integrate`, the method gives `None` here.

use super::conds::my_unpolarify;
use super::{find_splitting_points, get_coeff_exp, has_head_in, hyperbolic_heads, rewrite1, G};
use crate::expr::Constant;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::special::SpecialFn;
use crate::integrate::piecewise::{piecewise_fold, raw_piecewise};
use crate::{Expr, Func, Kind};

/// `_has(res, hyper, meijerg)`: every piece of the folded result contains
/// a hypergeometric or Meijer G-function.
fn has_hyper(e: &Expr) -> bool {
    fn any(e: &Expr) -> bool {
        crate::special::has_head(e, SpecialFn::MeijerG) || crate::special::has_head(e, SpecialFn::Hyper)
    }
    match piecewise_fold(e, Some(true)).kind() {
        Kind::Piecewise(pieces) => pieces.iter().all(|(v, c)| any(v) || any(c)),
        _ => any(e),
    }
}

fn has_nan_or_zoo(e: &Expr) -> bool {
    let mut bad = false;
    crate::traverse::preorder(e, &mut |t| {
        if matches!(t.kind(), Kind::Const(Constant::NaN | Constant::ComplexInfinity)) {
            bad = true;
        }
    });
    bad
}

/// `_clean`: multiply out the powers of `x` and drop the constant terms.
fn clean(e: &Expr, x: &Expr) -> Expr {
    // the reference `cancel` expands the numerator and the denominator
    // (logarithms of polar factors included) while building its polynomials
    let (n, d) = crate::poly::api::as_numer_denom(e);
    let full = |t: &Expr| crate::simptools::expand_hints(t, true, true, true, true, true, true);
    let e = div(&full(&n), &full(&d));
    let c = crate::poly::api::cancel(&e).unwrap_or_else(|_| e.clone());
    let c = crate::simptools::expand_hints_shallow_mul(&c);
    let (_, dep) = crate::decompose::as_coeff_add(&c, std::slice::from_ref(x));
    add(&dep)
}

fn is_true(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Bool(true))
}

/// `_meijerint_indefinite_1`: no substitution of `x`. `Ok(None)` when the
/// integrand has no rewrite, `Err` for the errors that abort the method.
fn indefinite_1(f: &Expr, x: &Expr) -> Result<Option<Expr>, ()> {
    let Some((fac, po, gl, cond)) = rewrite1(f, x, true) else { return Ok(None) };
    let mut res = Vec::new();
    for (big_c, s, g) in gl {
        let (a, b) = get_coeff_exp(&g.z, x).ok_or(())?;
        let (_, c) = get_coeff_exp(&po, x).ok_or(())?;
        let c = add(&[c, s]);
        // the substitution t = a*x**b gives fac*t**rho*G(t)
        let fac_ = div(&mul(&[fac.clone(), big_c, pow(x, &add(&[Expr::one(), c.clone()]))]), &b);
        let rho = div(&add(&[c, Expr::one()]), &b);
        let z = mul(&[a, pow(x, &b)]);
        let one_rho = sub(&Expr::one(), &rho);
        let neg_rho = neg(&rho);
        let shifted_nonpositive = g.bm.iter().any(|bj| {
            let v = add(&[bj.clone(), rho.clone()]);
            crate::assumptions::fact(&v, "integer") == Some(true) && is_true(&crate::relational::le(&v, &Expr::zero()))
        });
        let r = if shifted_nonpositive {
            let mut ao = g.ao.clone();
            ao.push(one_rho);
            let mut bm = g.bm.clone();
            bm.push(neg_rho);
            neg(&G::new(g.an.clone(), ao, bm, g.bo.clone(), z).ok_or(())?.to_expr())
        } else {
            let mut an = g.an.clone();
            an.push(one_rho);
            let mut bo = g.bo.clone();
            bo.push(neg_rho);
            G::new(an, g.ao.clone(), g.bm.clone(), bo, z).ok_or(())?.to_expr()
        };
        // the antiderivative is most often expected near x = 0
        let place = (crate::assumptions::fact(&b, "extended_nonnegative") == Some(true)
            && !has_nan_or_zoo(&crate::subs::subs1(f, x, &Expr::zero())))
        .then(Expr::zero);
        let r = crate::special::hyperexpand_options(&r, false, Some("default"), place.as_ref()).map_err(|_| ())?;
        res.push(super::powdenest_polar(&mul(&[fac_, r])));
    }
    let res = piecewise_fold(&add(&res), None);
    let res = match res.kind() {
        Kind::Piecewise(pieces) => raw_piecewise(
            pieces.iter().map(|(e, c)| (my_unpolarify(&clean(e, x)), c.clone())).collect(),
        ),
        _ => my_unpolarify(&clean(&res, x)),
    };
    Ok(Some(crate::relational::piecewise(vec![
        (res, my_unpolarify(&cond)),
        (crate::integrate::integral_unevaluated(f, x), Expr::boolean(true)),
    ])))
}

/// `meijerint_indefinite(f, x)`.
pub(crate) fn meijerint_indefinite(f: &Expr, x: &Expr) -> Option<Expr> {
    let mut points = vec![Expr::zero()];
    find_splitting_points(f, x, &mut points);
    points.sort_by(crate::sortkey::cmp_default);
    points.dedup();
    let mut results = Vec::new();
    for a in points {
        let shifted = crate::subs::subs1(f, x, &add(&[x.clone(), a.clone()]));
        let Some(res) = indefinite_1(&shifted, x).ok()? else { continue };
        let res = crate::subs::subs1(&res, x, &sub(x, &a));
        if has_hyper(&res) {
            results.push(res);
        } else {
            return Some(res);
        }
    }
    if has_head_in(f, &hyperbolic_heads()) {
        let rewritten = crate::special::rewrite::rewrite(
            f,
            &["sinh", "cosh", "tanh", "coth", "sech", "csch"].map(String::from),
            "exp",
            true,
        );
        if let Some(rv) = meijerint_indefinite(&rewritten, x) {
            let mut exps = Vec::new();
            crate::traverse::preorder(&rv, &mut |s| {
                if s.is_func(&Func::Exp) && !exps.contains(s) {
                    exps.push(s.clone())
                }
            });
            return Some(crate::poly::api::collect(
                &crate::simptools::exprtools::factor_terms(&rv, false, false, false, true),
                &exps,
            ));
        }
    }
    results.into_iter().min_by(|r1, r2| {
        crate::simptools::nodes(r1)
            .cmp(&crate::simptools::nodes(r2))
            .then_with(|| crate::sortkey::cmp_default(r1, r2))
    })
}
