//! Method of variation of parameters, as the reference (SymPy 1.14,
//! `sympy/solvers/ode/nonhomogeneous.py`, `_solve_variation_of_parameters`)
//! builds it:
//!
//! ```text
//! wr = wronskian(roots, x)            # simplified when simplify_flag
//! psol = Σ_i (-1)**(order - i) * Integral(wronskian(roots \ {r_i}, x)*r[-1]/wr, x) * r_i / r[order]
//! sol  = homogen_sol.rhs + psol       # psol simplified when simplify_flag
//! ```
//!
//! The Wronskians, the integrands, the held integrals and the sums are built
//! here with the reference's binary operations in the reference's order
//! (`a*b`, `a/b`, `0 + t`, `psol + t`), so the expression is the reference's
//! one. The two simplifications (`simplify` then `trigsimp(deep=True,
//! recursive=True)` of the Wronskian, `simplify` then `trigsimp(deep=True)`
//! of `psol`) are supplied by the caller. The integrals stay unevaluated, as
//! in the reference.
//!
//! `Ok(None)` declines: the number of solutions differs from the order, the
//! Wronskian is zero (the reference raises `NotImplementedError`), the
//! Wronskian computation fails, or an integrand leaves the content on which
//! the reference `Integral` constructor is `binders::held` (a `Piecewise`
//! that it would fold, `nan`, a bare applied function whose class could
//! define `_eval_Integral`).

use crate::binders::held;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::matrix::wronskian;
use crate::ops::{add2, div, mul2, pow};

/// Content on which the reference `Integral(body, x)` constructor is
/// `binders::held`: no `Piecewise` folding, no `nan`, no user hook on a bare
/// applied function.
fn admissible_integrand(body: &Expr) -> bool {
    if body.is_constant(Constant::NaN) || matches!(body.kind(), Kind::Func(Func::Undefined(_), _)) {
        return false;
    }
    let mut piecewise = false;
    crate::traverse::preorder(body, &mut |s| {
        if matches!(s.kind(), Kind::Piecewise(_)) {
            piecewise = true;
        }
    });
    !piecewise
}

/// `_solve_variation_of_parameters`: the right-hand side
/// `homogen_sol.rhs + psol` of the returned equation. `forcing` is
/// `match_obj[-1]`, `lead` is `match_obj[order]`. `simplify_wronskian` and
/// `simplify_psol` are called only when `simplify_flag` is set.
#[allow(clippy::too_many_arguments)]
pub fn variation_of_parameters<E>(
    roots: &[Expr],
    x: &Expr,
    forcing: &Expr,
    lead: &Expr,
    order: usize,
    homogen_rhs: &Expr,
    simplify_flag: bool,
    simplify_wronskian: &mut dyn FnMut(&Expr) -> Result<Expr, E>,
    simplify_psol: &mut dyn FnMut(&Expr) -> Result<Expr, E>,
) -> Result<Option<Expr>, E> {
    if roots.len() != order || order == 0 || !x.is_symbol() {
        return Ok(None);
    }
    let Ok(mut wr) = wronskian(roots, x) else {
        return Ok(None);
    };
    if simplify_flag {
        wr = simplify_wronskian(&wr)?;
    }
    // `not wr`: only an exact zero is false; relational and boolean values
    // are left to the reference
    match wr.kind() {
        Kind::Num(n) if n.is_zero() => return Ok(None),
        Kind::Rel(..) | Kind::Bool(_) => return Ok(None),
        _ => {}
    }
    let mut negoneterm = pow(&Expr::int(-1), &Expr::int(order as i64));
    let mut psol = Expr::zero();
    for i in roots {
        let others: Vec<Expr> = roots.iter().filter(|s| *s != i).cloned().collect();
        let Ok(w) = wronskian(&others, x) else {
            return Ok(None);
        };
        let body = div(&mul2(&w, forcing), &wr);
        if !admissible_integrand(&body) {
            return Ok(None);
        }
        let Ok(integral) = held(Func::Integral, &body, vec![vec![x.clone()]]) else {
            return Ok(None);
        };
        let term = div(&mul2(&mul2(&negoneterm, &integral), i), lead);
        psol = add2(&psol, &term);
        negoneterm = mul2(&negoneterm, &Expr::int(-1));
    }
    if simplify_flag {
        psol = simplify_psol(&psol)?;
    }
    Ok(Some(add2(homogen_rhs, &psol)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(s: &str) -> Expr {
        parse_str(s).unwrap()
    }

    #[test]
    fn unsimplified_construction() {
        let x = p("x");
        let roots = vec![p("exp(x)"), p("x*exp(x)")];
        let mut id = |e: &Expr| -> Result<Expr, ()> { Ok(e.clone()) };
        let mut id2 = |e: &Expr| -> Result<Expr, ()> { Ok(e.clone()) };
        let rhs = variation_of_parameters(&roots, &x, &p("exp(x)/x"), &p("1"), 2, &p("C1*exp(x) + C2*x*exp(x)"), false, &mut id, &mut id2)
            .unwrap()
            .unwrap();
        let wr = wronskian(&roots, &x).unwrap();
        let i1 = held(Func::Integral, &div(&mul2(&p("x*exp(x)"), &p("exp(x)/x")), &wr), vec![vec![x.clone()]]).unwrap();
        let i2 = held(Func::Integral, &div(&mul2(&p("exp(x)"), &p("exp(x)/x")), &wr), vec![vec![x.clone()]]).unwrap();
        let want = add2(&p("C1*exp(x) + C2*x*exp(x)"), &add2(&mul2(&mul2(&Expr::int(1), &i1), &p("exp(x)")), &mul2(&mul2(&Expr::int(-1), &i2), &p("x*exp(x)"))));
        assert_eq!(rhs, want);
        // order mismatch and dependent solutions are declined
        assert!(variation_of_parameters(&roots, &x, &p("exp(x)"), &p("1"), 3, &p("0"), false, &mut id, &mut id2).unwrap().is_none());
        let dependent = vec![p("exp(x)"), p("2*exp(x)")];
        assert!(variation_of_parameters(&dependent, &x, &p("exp(x)"), &p("1"), 2, &p("0"), false, &mut id, &mut id2).unwrap().is_none());
    }
}
