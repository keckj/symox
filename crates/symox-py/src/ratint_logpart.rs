//! Lazard-Rioboo-Trager orchestration of `ratint_logpart` for the coefficient
//! domains outside the plain QQ/ZZ kernel (`ZZ_I`, `QQ_I`, `ZZ[a]`, `QQ(a)`,
//! `EX`, ...).
//!
//! Every step is the reference step on the same objects, with the stock
//! native Poly methods that the caller has verified against their
//! descriptors. The public constructor calls of the reference body are
//! replaced by the backend calls the stock constructor makes for exactly
//! these arguments (`_from_normalized`, the term path `_from_terms`, the
//! `composite=False` EX selection and `to_field`); the float-domain pass of
//! the constructor is a no-op because inexact inputs are refused. The
//! domains of the results, the order of the pairs and the Dummy are
//! therefore those of the reference. Any input the reference would reject
//! (zero resultant, missing PRS degree, failed coefficient conversion)
//! returns `None` before anything is published, and the reference body then
//! runs and raises its own error.
use pyo3::prelude::*;
use pyo3::basic::CompareOp;
use pyo3::types::{PyDict, PyTuple};

use super::{w, PyPoly};
use symox::integrate::ratint_algebraic::{self, ExtensionRingPoly};
use symox::poly::algebraic::AlgebraicPoly;
use symox::poly::expression_domain::ExpressionPoly;

type Pairs = Vec<(PyObject, PyObject)>;

fn exact(p: &PyPoly) -> bool {
    if p.p.base.modulus.is_some() || p.p.base.domain.is_rr() { return false; }
    if p.p.algebraic.is_some() || p.p.field.is_some() { return true; }
    match &p.p.expression {
        Some(e) => e.ground.as_ref().is_none_or(|g| g.algebraic.is_some()) && e.as_expr().float_info().max_precision.is_none(),
        None => matches!(p.p.base.domain, symox::poly::polyobj::Domain::ZZ | symox::poly::polyobj::Domain::QQ),
    }
}

/// `Poly(rep, *gens)` for an expression or Poly `rep` without options.
fn construct(rep: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>) -> PyResult<Option<PyPoly>> {
    match PyPoly::_from_normalized(rep, gens, None, None, None) {
        Ok(p) => Ok(Some(p)),
        Err(e) if e.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(rep.py()) => Err(e),
        Err(_) => Ok(None),
    }
}

/// The square-free list as owned pairs.
fn listed(value: Bound<'_, PyTuple>) -> PyResult<Option<(Bound<'_, PyAny>, Vec<(PyPoly, Bound<'_, PyAny>)>)>> {
    if value.len() != 2 { return Ok(None); }
    let mut factors = Vec::new();
    for item in value.get_item(1)?.try_iter()? {
        let item = item?;
        let Ok(item) = item.downcast::<PyTuple>() else { return Ok(None) };
        if item.len() != 2 { return Ok(None); }
        let Ok(factor) = item.get_item(0)?.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
        factors.push((factor, item.get_item(1)?));
    }
    Ok(Some((value.get_item(0)?, factors)))
}

/// `_include_sign(c, sqf)` of the reference.
fn include_sign(c: &Bound<'_, PyAny>, sqf: &mut [(PyPoly, Bound<'_, PyAny>)]) -> PyResult<Option<()>> {
    if !c.getattr("is_extended_real")?.is_truthy()? { return Ok(Some(())); }
    if !c.rich_compare(0, CompareOp::Lt)?.eq(true)? { return Ok(Some(())); }
    let Some((h, _)) = sqf.first_mut() else { return Ok(None) };
    // c.as_poly(h.gens), then h*c_poly.
    let Some(c_poly) = construct(c, &h.gens(c.py())?)? else { return Ok(None) };
    let product = h.__mul__(&Bound::new(c.py(), c_poly)?.into_any())?;
    let Ok(product) = product.bind(c.py()).downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
    *h = product;
    Ok(Some(()))
}

fn int(value: &Bound<'_, PyAny>) -> Option<i64> { value.extract::<i64>().ok() }

/// The reference body after its admission guard, for the non-plain domains.
pub(crate) fn logpart(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>,
    t: &Bound<'_, PyAny>) -> PyResult<Option<Pairs>>
{
    let py = f.py();
    let (Ok(fp), Ok(gp)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) else { return Ok(None) };
    if !exact(fp.get()) || !exact(gp.get()) { return Ok(None); }
    let in_x = PyTuple::new(py, [x])?;
    let in_t = PyTuple::new(py, [t])?;
    let new = |p: PyPoly| Bound::new(py, p).map(|b| b.into_any());
    // f, g = Poly(f, x), Poly(g, x)
    let (Some(f), Some(g)) = (construct(f, &in_x)?, construct(g, &in_x)?) else { return Ok(None) };
    if let Some((fa, ga)) = extension_pair(&f, &g) {
        return extension_logpart(py, &g, &fa, &ga, x, t);
    }
    let g_object = new(g.clone())?;
    py.check_signals()?;
    // a, b = g, f - g.diff()*Poly(t, x)
    let Some(t_poly) = construct(t, &in_x)? else { return Ok(None) };
    let product = new(g.diff(&PyTuple::empty(py))?)?.mul(new(t_poly)?)?;
    let b = new(f)?.sub(product)?;
    if !b.is_exact_instance_of::<PyPoly>() { return Ok(None); }
    // res, R = resultant(a, b, includePRS=True)
    let value = g.resultant(&b, true)?;
    let Ok(value) = value.bind(py).downcast::<PyTuple>().cloned() else { return Ok(None) };
    if value.len() != 2 { return Ok(None); }
    // res = Poly(res, t, composite=False)
    let Some(res) = construct(&value.get_item(0)?, &in_t)? else { return Ok(None) };
    let ground = new(res.clone())?.getattr("domain")?;
    let composite = ground.getattr("is_PolynomialRing").map_or(Ok(false), |v| v.is_truthy())?
        || ground.getattr("is_FractionField").map_or(Ok(false), |v| v.is_truthy())?;
    let res = if composite {
        res.set_domain(&pyo3::types::PyString::new(py, "EX").into_any())?
    } else { res };
    // assert res: the reference raises its own error.
    if !res.__bool__() { return Ok(None); }
    py.check_signals()?;
    let mut map: Vec<(i64, PyPoly)> = Vec::new();
    for r in value.get_item(1)?.try_iter()? {
        let r = r?;
        let Ok(r) = r.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
        // A non-integer degree (-oo) is never looked up below.
        let Some(degree) = int(r.degree(py, None)?.bind(py)) else { continue };
        match map.iter_mut().find(|(d, _)| *d == degree) {
            Some(slot) => slot.1 = r,
            None => map.push((degree, r)),
        }
    }
    let Some((c, mut res_sqf)) = listed(res.sqf_list(py, false)?)? else { return Ok(None) };
    if include_sign(&c, &mut res_sqf)?.is_none() { return Ok(None); }
    let Some(g_degree) = int(g.degree(py, None)?.bind(py)) else { return Ok(None) };
    let one = py.import("symox.core.singleton")?.getattr("S")?.getattr("One")?;
    let from_terms = py.import("symox._poly_constructor")?.getattr("_from_terms")?;
    let mut pairs = Vec::new();
    for (q, i) in res_sqf {
        py.check_signals()?;
        let primitive = q.primitive(py)?;
        if primitive.len() != 2 { return Ok(None); }
        let Ok(q) = primitive.get_item(1)?.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
        let Some(i) = int(&i) else { return Ok(None) };
        let q_object = new(q.clone())?;
        if g_degree == i {
            pairs.push((g_object.clone().unbind(), q_object.unbind()));
            continue;
        }
        let Some((_, h)) = map.iter().find(|(d, _)| *d == i) else { return Ok(None) };
        let mut h = h.clone();
        // h_lc = Poly(h.LC(), t, field=True)
        let lc = h.LC().into_pyobject(py)?.into_any();
        let Some(h_lc) = construct(&lc, &in_t)? else { return Ok(None) };
        let h_lc = h_lc.to_field()?;
        let Some((c, mut h_lc_sqf)) = listed(h_lc.sqf_list(py, true)?)? else { return Ok(None) };
        if include_sign(&c, &mut h_lc_sqf)?.is_none() { return Ok(None); }
        for (a, j) in h_lc_sqf {
            py.check_signals()?;
            let Some(j) = int(&j).and_then(|j| u32::try_from(j).ok()) else { return Ok(None) };
            let power = new(a.gcd(&q_object)?.__pow__(j, None)?)?;
            let Some(divisor) = construct(&power, &in_x)? else { return Ok(None) };
            h = h.quo(&new(divisor)?, true)?;
        }
        let inv = h_lc.invert(&q_object, true)?;
        let inv_gens = inv.gens(py)?;
        let inv_object = new(inv)?;
        let terms = h.terms(py)?;
        let mut coeffs = vec![one.clone()];
        for term in terms.iter().skip(1) {
            py.check_signals()?;
            let Some(coeff) = construct(&term.get_item(1)?, &inv_gens)? else { return Ok(None) };
            let product = inv_object.mul(new(coeff)?)?;
            let Ok(product) = product.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
            coeffs.push(product.rem(&q_object, true)?.as_expr().into_pyobject(py)?.into_any());
        }
        let data = PyDict::new(py);
        for (term, coeff) in terms.iter().zip(coeffs) {
            data.set_item(term.get_item(0)?, coeff)?;
        }
        // Poly(dict(...), x): the stock constructor's term path (its
        // coefficient sympification is the identity on these expressions).
        let h = from_terms.call1((data, &in_x, py.None(), py.None(), py.None()))?;
        if !h.is_exact_instance_of::<PyPoly>() { return Ok(None); }
        pairs.push((h.unbind(), q_object.unbind()));
    }
    py.check_signals()?;
    Ok(Some(pairs))
}

/// A polynomial over `ZZ` or `QQ` without coefficient generators.
fn plain_rational(p: &PyPoly) -> bool {
    p.p.algebraic.is_none() && p.p.field.is_none() && p.p.expression.is_none()
        && p.p.base.modulus.is_none() && !p.p.base.has_domain_gens()
        && matches!(p.p.base.domain, symox::poly::polyobj::Domain::ZZ | symox::poly::polyobj::Domain::QQ)
}

/// `f, g` over one algebraic field `QQ<a>` (not the Gaussian domains), or
/// one of them over `ZZ`/`QQ`, which the reference unifies into `QQ<a>`.
fn extension_pair(f: &PyPoly, g: &PyPoly) -> Option<(AlgebraicPoly, AlgebraicPoly)> {
    let genuine = |p: &PyPoly| p.p.algebraic.as_ref().filter(|a| a.gaussian.is_none() && a.gens.len() == 1).cloned();
    let (fa, ga) = (genuine(f), genuine(g));
    let lift = |p: &PyPoly, other: &AlgebraicPoly| -> Option<AlgebraicPoly> {
        if !plain_rational(p) { return None; }
        other.convert(&p.p.as_expr()).ok().filter(|a| a.gens == other.gens)
    };
    let (fa, ga) = match (fa, ga) {
        (Some(fa), Some(ga)) => (fa, ga),
        (Some(fa), None) => { let ga = lift(g, &fa)?; (fa, ga) }
        (None, Some(ga)) => { let fa = lift(f, &ga)?; (fa, ga) }
        (None, None) => return None,
    };
    if fa.gens != ga.gens { return None; }
    if fa.number_field.same_embedding(&ga.number_field) { return Some((fa, ga)); }
    // Two different fields: the reference unifies `QQ<a>` with `QQ<b>`
    // into `QQ<b, a>` (`Domain.unify` swaps algebraic operands), a field of
    // one primitive element; the pairs come out over EX, independent of the
    // primitive element chosen.
    let field = symox::poly::algebraic::NumberField::from_generators(&[
        ga.number_field.generator.clone(), fa.number_field.generator.clone()]).ok()?;
    let convert = |p: &AlgebraicPoly| AlgebraicPoly::from_expr(&p.as_expr(), &p.gens, field.clone()).ok();
    Some((convert(&fa)?, convert(&ga)?))
}

/// A subresultant over `QQ<a>[t]`, until a division by an `EX` divisor
/// moves it to `EX` (the reference unification of `QQ<a>[t]` with `EX`).
enum Member {
    Ring(ExtensionRingPoly),
    Stock(PyPoly),
}

/// The reference body when the coefficients lie in `QQ<a>`: `b`, the
/// resultant and the subresultants over `QQ<a>[t]` are computed natively
/// (the stock Poly arithmetic has no `QQ<a>[t]` composite); every later
/// step is the generic one on stock polynomials.
fn extension_logpart(py: Python<'_>, g: &PyPoly, fa: &AlgebraicPoly, ga: &AlgebraicPoly,
    x: &Bound<'_, PyAny>, t: &Bound<'_, PyAny>) -> PyResult<Option<Pairs>>
{
    let in_x = PyTuple::new(py, [x])?;
    let in_t = PyTuple::new(py, [t])?;
    let new = |p: PyPoly| Bound::new(py, p).map(|b| b.into_any());
    let (xe, te) = (crate::to_expr(x)?, crate::to_expr(t)?);
    if ga.gens.first() != Some(&xe) { return Ok(None); }
    let Ok((res, sequence)) = ratint_algebraic::logpart_prs(fa, ga, &te) else { return Ok(None) };
    py.check_signals()?;
    // res = Poly(res, t, composite=False)
    let res_object = w(res).into_pyobject(py)?.into_any();
    let Some(res) = construct(&res_object, &in_t)? else { return Ok(None) };
    let ground = Bound::new(py, res.clone())?.getattr("domain")?;
    let composite = ground.getattr("is_PolynomialRing").map_or(Ok(false), |v| v.is_truthy())?
        || ground.getattr("is_FractionField").map_or(Ok(false), |v| v.is_truthy())?;
    let res = if composite {
        res.set_domain(&pyo3::types::PyString::new(py, "EX").into_any())?
    } else { res };
    if !res.__bool__() { return Ok(None); }
    let mut map: Vec<(u32, ExtensionRingPoly)> = Vec::new();
    for r in sequence {
        let Some(degree) = r.degree() else { continue };
        match map.iter_mut().find(|(d, _)| *d == degree) {
            Some(slot) => slot.1 = r,
            None => map.push((degree, r)),
        }
    }
    let Some((c, mut res_sqf)) = listed(res.sqf_list(py, false)?)? else { return Ok(None) };
    if include_sign(&c, &mut res_sqf)?.is_none() { return Ok(None); }
    let Some(g_degree) = int(g.degree(py, None)?.bind(py)) else { return Ok(None) };
    let g_object = new(g.clone())?;
    let one = py.import("symox.core.singleton")?.getattr("S")?.getattr("One")?;
    let from_terms = py.import("symox._poly_constructor")?.getattr("_from_terms")?;
    let mut pairs = Vec::new();
    for (q, i) in res_sqf {
        py.check_signals()?;
        let primitive = q.primitive(py)?;
        if primitive.len() != 2 { return Ok(None); }
        let Ok(q) = primitive.get_item(1)?.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
        let Some(i) = int(&i) else { return Ok(None) };
        let q_object = new(q.clone())?;
        if g_degree == i {
            pairs.push((g_object.clone().unbind(), q_object.unbind()));
            continue;
        }
        let Some((_, ring)) = map.iter().find(|(d, _)| i64::from(*d) == i) else { return Ok(None) };
        // h_lc = Poly(h.LC(), t, field=True)
        let lc = w(ring.lc_expr()).into_pyobject(py)?.into_any();
        let mut h = Member::Ring(ring.clone());
        let Some(h_lc) = construct(&lc, &in_t)? else { return Ok(None) };
        let h_lc = h_lc.to_field()?;
        let Some((c, mut h_lc_sqf)) = listed(h_lc.sqf_list(py, true)?)? else { return Ok(None) };
        if include_sign(&c, &mut h_lc_sqf)?.is_none() { return Ok(None); }
        for (a, j) in h_lc_sqf {
            py.check_signals()?;
            let Some(j) = int(&j).and_then(|j| u32::try_from(j).ok()) else { return Ok(None) };
            let power = new(a.gcd(&q_object)?.__pow__(j, None)?)?;
            let Some(divisor) = construct(&power, &in_x)? else { return Ok(None) };
            h = match h {
                Member::Stock(p) => Member::Stock(p.quo(&new(divisor)?, true)?),
                Member::Ring(r) if divisor.p.expression.is_some() => {
                    // unify(QQ<a>[t], EX) = EX, each coefficient being
                    // EX(QQ<a>[t].to_sympy(c)).
                    let terms = r.terms().into_iter().map(|(k, c)| (vec![k], c)).collect();
                    let Ok(p) = ExpressionPoly::from_terms(terms, std::slice::from_ref(&xe)) else { return Ok(None) };
                    Member::Stock(PyPoly { p: p.into() }.quo(&new(divisor)?, true)?)
                }
                Member::Ring(r) => {
                    // unify(QQ<a>[t], ZZ, QQ, ZZ[t] or QQ[t]) = QQ<a>[t]; the
                    // quotient over QQ<a>(t) must retract to the ring.
                    if divisor.p.algebraic.is_some() || divisor.p.field.is_some() || divisor.p.expression.is_some()
                        || divisor.p.base.modulus.is_some() || divisor.p.base.domain.is_rr() { return Ok(None); }
                    let Some(r) = r.quo_by_t_polynomial(&divisor.p.as_expr(), &xe) else { return Ok(None) };
                    Member::Ring(r)
                }
            };
        }
        let inv = h_lc.invert(&q_object, true)?;
        let inv_gens = inv.gens(py)?;
        let inv_object = new(inv)?;
        let terms: Vec<(Bound<'_, PyAny>, Bound<'_, PyAny>)> = match &h {
            Member::Ring(r) => {
                let terms = r.terms();
                let terms = if terms.is_empty() { vec![(0, symox::Expr::zero())] } else { terms };
                terms.into_iter()
                    .map(|(k, c)| Ok((PyTuple::new(py, [k])?.into_any(), w(c).into_pyobject(py)?.into_any())))
                    .collect::<PyResult<_>>()?
            }
            Member::Stock(p) => p.terms(py)?.into_iter()
                .map(|term| Ok((term.get_item(0)?, term.get_item(1)?)))
                .collect::<PyResult<_>>()?,
        };
        let mut coeffs = vec![one.clone()];
        for (_, coeff) in terms.iter().skip(1) {
            py.check_signals()?;
            let Some(coeff) = construct(coeff, &inv_gens)? else { return Ok(None) };
            let product = inv_object.mul(new(coeff)?)?;
            let Ok(product) = product.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
            coeffs.push(product.rem(&q_object, true)?.as_expr().into_pyobject(py)?.into_any());
        }
        let data = PyDict::new(py);
        for ((monom, _), coeff) in terms.iter().zip(coeffs) {
            data.set_item(monom, coeff)?;
        }
        // Poly(dict(...), x)
        let h = from_terms.call1((data, &in_x, py.None(), py.None(), py.None()))?;
        if !h.is_exact_instance_of::<PyPoly>() { return Ok(None); }
        pairs.push((h.unbind(), q_object.unbind()));
    }
    py.check_signals()?;
    Ok(Some(pairs))
}
