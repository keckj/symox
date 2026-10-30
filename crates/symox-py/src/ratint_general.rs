//! Native `ratint_ratpart` over the algebraic (`QQ<a>`, `ZZ_I`, `QQ_I`),
//! rational function (`QQ(a, ...)`) and expression (`EX`) coefficient
//! fields. The caller has already checked the stock Poly protocols; this
//! module only selects stored representations and returns `None` for any
//! other input, so the reference algorithm then runs unchanged.
use pyo3::prelude::*;
use symox::expr::Expr;
use symox::poly::fraction::Polynomial;
use symox::poly::{dense, Coeff, MPoly};

use crate::poly::PyPoly;

/// The same representation family for numerator and denominator, each
/// univariate in `x`.
enum Family<'a> {
    Algebraic(&'a symox::poly::algebraic::AlgebraicPoly, &'a symox::poly::algebraic::AlgebraicPoly),
    Fraction(&'a symox::poly::fraction::FractionPoly, &'a symox::poly::fraction::FractionPoly),
    Expression(&'a symox::poly::expression_domain::ExpressionPoly, &'a symox::poly::expression_domain::ExpressionPoly),
}

fn family<'a>(f: &'a Polynomial, g: &'a Polynomial, x: &Expr) -> Option<Family<'a>> {
    if f.base.modulus.is_some() || g.base.modulus.is_some() { return None; }
    let univariate = |gens: &[Expr]| gens.len() == 1 && gens[0] == *x;
    match (&f.algebraic, &g.algebraic, &f.field, &g.field, &f.expression, &g.expression) {
        (Some(a), Some(b), None, None, None, None)
            if univariate(&a.gens) && univariate(&b.gens) && a.gaussian == b.gaussian
                && a.number_field.same_embedding(&b.number_field) => Some(Family::Algebraic(a, b)),
        (None, None, Some(a), Some(b), None, None)
            if univariate(&a.gens) && univariate(&b.gens) && a.parameters == b.parameters
                && a.display == b.display && !a.parameters.is_empty() => Some(Family::Fraction(a, b)),
        (None, None, None, None, Some(a), Some(b))
            if univariate(&a.gens) && univariate(&b.gens) && a.ground.is_none() && b.ground.is_none()
                => Some(Family::Expression(a, b)),
        _ => None,
    }
}

fn dense_of<C: Coeff>(p: &MPoly<C>) -> Option<Vec<C>> {
    (p.nvars == 1).then(|| dense::from_mpoly(p, 0))
}

/// `(A, u, B, v)` as expressions, in the reference order of construction.
type Parts = (Expr, Expr, Expr, Expr);

fn reduce(py: Python<'_>, f: &Polynomial, g: &Polynomial, x: &Expr) -> PyResult<Option<Parts>> {
    let Some(family) = family(f, g, x) else { return Ok(None) };
    // u, v, _ = g.cofactors(g.diff()), with the stock Poly methods.
    let derivative = match &g.expression {
        Some(e) => match e.diff(1, 0) { Ok(d) => Polynomial::from(d), Err(_) => return Ok(None) },
        None => g.diff(0),
    };
    let Ok((u, v, _)) = g.checked_cofactors(&derivative) else { return Ok(None) };
    let mut check = || py.check_signals();
    match family {
        Family::Algebraic(fa, ga) => {
            let (Some(ua), Some(va)) = (&u.algebraic, &v.algebraic) else { return Ok(None) };
            let _scope = symox::poly::algebraic::AlgebraicScope::new(ga.number_field.clone());
            let (Some(fd), Some(ud), Some(vd)) = (dense_of(&fa.rep), dense_of(&ua.rep), dense_of(&va.rep)) else { return Ok(None) };
            let Some((a, b)) = symox::integrate::ratint_general::solve(&fd, &ud, &vd, 96, &mut check)? else { return Ok(None) };
            let build = |p: &[symox::poly::algebraic::AN]| ga.with_rep(dense::to_mpoly(p, 0, 1)).as_expr();
            Ok(Some((build(&a), ua.as_expr(), build(&b), va.as_expr())))
        },
        Family::Fraction(fa, ga) => {
            let (Some(uf), Some(vf)) = (&u.field, &v.field) else { return Ok(None) };
            let _scope = symox::poly::ratfunc::Scope::new(ga.parameters.len());
            let (Some(fd), Some(ud), Some(vd)) = (dense_of(&fa.rep), dense_of(&uf.rep), dense_of(&vf.rep)) else { return Ok(None) };
            let Some((a, b)) = symox::integrate::ratint_general::solve(&fd, &ud, &vd, 48, &mut check)? else { return Ok(None) };
            let build = |p: &[symox::poly::ratfunc::RF]| {
                let mut q = ga.clone();
                q.rep = dense::to_mpoly(p, 0, 1);
                q.as_expr()
            };
            Ok(Some((build(&a), uf.as_expr(), build(&b), vf.as_expr())))
        },
        Family::Expression(fe, ge) => {
            let (Some(ue), Some(ve)) = (&u.expression, &v.expression) else { return Ok(None) };
            if let Some((a, b)) = parametric(&[fe.as_expr(), ue.as_expr(), ve.as_expr()], x, &mut check)? {
                return Ok(Some((a, ue.as_expr(), b, ve.as_expr())));
            }
            let solved = symox::poly::expression_domain::guarded(|| {
                let (Some(fd), Some(ud), Some(vd)) = (dense_of(&fe.to_sparse()), dense_of(&ue.to_sparse()), dense_of(&ve.to_sparse()))
                    else { return Ok(None) };
                symox::integrate::ratint_general::solve(&fd, &ud, &vd, 16, &mut check)
            });
            let Ok(solved) = solved else { return Ok(None) };
            let Some((a, b)) = solved? else { return Ok(None) };
            let build = |p: &[symox::poly::expression_domain::EX]| ge.with_sparse_rep(dense::to_mpoly(p, 0, 1)).as_expr();
            Ok(Some((build(&a), ue.as_expr(), build(&b), ve.as_expr())))
        },
    }
}

/// `EX` coefficients that are all rational functions of commuting symbols:
/// the unique solution is computed in the fraction field `QQ(a, ...)`
/// instead of with expression arithmetic. Only the values of `A`, `B` enter
/// the result, which `cancel(., x)` then puts in its canonical form.
fn parametric(polys: &[Expr; 3], x: &Expr, check: &mut impl FnMut() -> PyResult<()>)
    -> PyResult<Option<(Expr, Expr)>>
{
    use symox::poly::fraction::FractionPoly;
    let mut parameters: Vec<Expr> = Vec::new();
    for e in polys {
        for s in symox::traverse::free_symbols(e) {
            if s != *x && !parameters.contains(&s) { parameters.push(s); }
        }
    }
    if parameters.is_empty() || parameters.iter().any(|s| !matches!(s.kind(), symox::expr::Kind::Sym(_))
        || symox::decompose::commutativity(s) != Some(true)) {
        return Ok(None);
    }
    parameters.sort_by(|p, q| p.cmp_struct(q));
    let gens = std::slice::from_ref(x);
    let converted = polys.iter().map(|e| FractionPoly::from_expr(e, gens, &parameters, "QQ").ok())
        .collect::<Option<Vec<_>>>();
    let Some(converted) = converted else { return Ok(None) };
    let _scope = symox::poly::ratfunc::Scope::new(parameters.len());
    let (Some(fd), Some(ud), Some(vd)) = (dense_of(&converted[0].rep), dense_of(&converted[1].rep), dense_of(&converted[2].rep))
        else { return Ok(None) };
    let Some((a, b)) = symox::integrate::ratint_general::solve(&fd, &ud, &vd, 48, check)? else { return Ok(None) };
    let build = |p: &[symox::poly::ratfunc::RF]| converted[0].with_rep(dense::to_mpoly(p, 0, 1)).as_expr();
    Ok(Some((build(&a), build(&b))))
}

/// `cancel(A/u.as_expr(), x), cancel(B/v.as_expr(), x)` through the
/// `cancel` global of rationaltools, exactly as the reference body calls it.
pub(crate) fn ratpart(py: Python<'_>, f: &PyPoly, g: &PyPoly, x: &Expr, cancel: &Bound<'_, PyAny>)
    -> PyResult<Option<(PyObject, PyObject)>>
{
    if f.p.true_gens() != g.p.true_gens() { return Ok(None); }
    let Some((a, u, b, v)) = reduce(py, &f.p, &g.p, x)? else { return Ok(None) };
    py.check_signals()?;
    let variable = crate::wrap(x.clone()).into_pyobject(py)?.into_any().unbind();
    let cancelled = |p: &Expr, q: &Expr| -> PyResult<PyObject> {
        Ok(cancel.call1((crate::wrap(symox::ops::div(p, q)), variable.clone_ref(py)))?.unbind())
    };
    let rational = cancelled(&a, &u)?;
    let logarithmic = cancelled(&b, &v)?;
    Ok(Some((rational, logarithmic)))
}
