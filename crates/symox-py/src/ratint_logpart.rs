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

use super::PyPoly;

type Pairs = Vec<(PyObject, PyObject)>;

fn exact(p: &PyPoly) -> bool {
    if p.p.base.modulus.is_some() || p.p.base.domain.is_rr() { return false; }
    if p.p.algebraic.is_some() || p.p.field.is_some() { return true; }
    match &p.p.expression {
        Some(e) => e.ground.is_none() && e.as_expr().float_info().max_precision.is_none(),
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
