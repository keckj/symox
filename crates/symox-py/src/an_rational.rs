//! Coercion and return-type boundary for exact rational operations over AN.
#![allow(non_snake_case)]
use super::{
    err, expr_or_coeff_list, gens_of, w, AlgebraicPoly, NumberField,
    PolyValue, PyPoly,
};
use crate::{to_expr, PyExpr};
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use symox::expr::Expr;
use symox::poly::algebraic_rational as native;
use symox::poly::convert::PolyError;

pyo3::create_exception!(symox, NotInvertible, PyValueError);

pub(super) struct Selected {
    pub a: AlgebraicPoly,
    pub b: AlgebraicPoly,
    pub polys: bool,
}
fn selected_options<'py>(
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<(Option<Bound<'py, PyAny>>, Option<Bound<'py, PyAny>>)> {
    let extension = kwargs
        .map(|k| k.get_item("extension"))
        .transpose()?
        .flatten()
        .filter(|v| !v.is_none());
    let gaussian = match kwargs {
        Some(k) => super::division::gaussian_option(k)?,
        None => None,
    };
    let domain = match gaussian {
        Some(d) => Some(d),
        None => kwargs
            .map(|k| k.get_item("domain"))
            .transpose()?
            .flatten()
            .filter(|v| !v.is_none()),
    };
    Ok((extension, domain))
}
pub(super) fn select(
    f: &Bound<'_, PyAny>,
    g: &Bound<'_, PyAny>,
    gens: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Option<Selected>> {
    let (fe, fg) = expr_or_coeff_list(f)?;
    let (ge, gg) = expr_or_coeff_list(g)?;
    let (extension, domain) = selected_options(kwargs)?;
    let selected=super::domains::Selection::new(domain.as_ref(),None,extension.as_ref())?;
    if selected.native.as_ref().is_some_and(|d|!matches!(d.inner,symox::poly::domains::ExactDomain::ZZ|symox::poly::domains::ExactDomain::QQ|symox::poly::domains::ExactDomain::Algebraic(_)|symox::poly::domains::ExactDomain::Gaussian{..})){return Err(PyNotImplementedError::new_err("this rational operation currently supports QQ and certified algebraic domains"));}
    let mut field=selected.algebraic(extension.as_ref(),&[fe.clone(),ge.clone()])?;
    if field.is_none() && domain.is_none() {
        for obj in [f, g] {
            if let Ok(p) = obj.downcast::<PyPoly>() {
                if let Some(a) = &p.get().p.algebraic {
                    field = Some(match field {
                        None => a.number_field.clone(),
                        Some(previous) => {
                            if previous == a.number_field {
                                previous
                            } else {
                                NumberField::from_generators(&[
                                    previous.generator.clone(),
                                    a.number_field.generator.clone(),
                                ])
                                .map_err(err)?
                            }
                        }
                    });
                }
            }
        }
    }
    let mut default_gens = None;
    if field.is_none() && domain.is_none() && extension.is_none()
        && [f, g].iter().all(|x| x.downcast::<PyPoly>().map_or(true, |p| p.get().p.algebraic.as_ref().is_some_and(|a| a.gaussian.is_some())))
    {
        // The reference default domain of Gaussian rational coefficients.
        if let Some(selection) = symox::poly::gaussian::default_selection(&[fe.clone(), ge.clone()], &gens_of(gens)?) {
            field = Some(symox::poly::gaussian::field());
            default_gens = Some(selection.gens);
        }
    }
    let gaussian = selected.gaussian().is_some() || default_gens.is_some()
        || (domain.is_none() && extension.is_none()
            && [f, g].iter().any(|x| x.downcast::<PyPoly>().is_ok_and(|p| p.get().p.algebraic.as_ref().is_some_and(|a| a.gaussian.is_some())))
            && [f, g].iter().all(|x| x.downcast::<PyPoly>().map_or(true, |p| p.get().p.algebraic.as_ref().is_some_and(|a| a.gaussian.is_some()))));
    let Some(field) = field else { return Ok(None) };
    if kwargs
        .map(|k| k.get_item("modulus"))
        .transpose()?
        .flatten()
        .is_some_and(|v| !v.is_none())
    {
        return Err(PyValueError::new_err(
            "algebraic extensions cannot use a finite modulus",
        ));
    }
    for obj in [f, g] {
        if let Ok(p) = obj.downcast::<PyPoly>() {
            if p.get().p.modulus.is_some() || p.get().p.field.is_some() {
                return Err(PyValueError::new_err(
                    "incompatible algebraic, finite, or rational-function coefficient fields",
                ));
            }
        }
    }
    let mut generators = default_gens.unwrap_or(gens_of(gens)?);
    if generators.is_empty() {
        generators = symox::poly::polyobj::unify_gens(&fg, &gg);
    }
    if generators.is_empty() {
        generators = symox::poly::algebraic::infer_generators(&[fe.clone(), ge.clone()], &field);
    }
    if generators.is_empty() {
        generators.push(Expr::symbol("_x"));
    }
    let polys = kwargs
        .map(|k| k.get_item("polys"))
        .transpose()?
        .flatten()
        .map(|v| v.extract::<bool>())
        .transpose()?
        .unwrap_or_else(|| f.downcast::<PyPoly>().is_ok() || g.downcast::<PyPoly>().is_ok());
    let mut a = AlgebraicPoly::from_expr(&fe, &generators, field.clone()).map_err(err)?;
    let mut b = AlgebraicPoly::from_expr(&ge, &generators, field).map_err(err)?;
    // Rational operations run over the fraction field of a Gaussian ring.
    if gaussian {
        a.gaussian = Some(false);
        b.gaussian = Some(false);
    }
    Ok(Some(Selected { a, b, polys }))
}
pub(super) fn poly_output(py: Python<'_>, p: AlgebraicPoly, polys: bool) -> PyResult<PyObject> {
    if polys {
        Ok(PyPoly { p: p.into() }
            .into_pyobject(py)?
            .into_any()
            .unbind())
    } else {
        Ok(w(p.as_expr()).into_pyobject(py)?.into_any().unbind())
    }
}
pub(super) fn bezout(py: Python<'_>, selected: Selected, half: bool) -> PyResult<PyObject> {
    if !half && selected.b.rep.is_zero() {
        return Err(err(PolyError::Unsupported(
            "polynomial division by zero".into(),
        )));
    }
    let (s, t, h) = native::gcdex(&selected.a, &selected.b).map_err(err)?;
    let mut values = vec![poly_output(py, s, selected.polys)?];
    if !half {
        values.push(poly_output(py, t, selected.polys)?);
    }
    values.push(poly_output(py, h, selected.polys)?);
    Ok(PyTuple::new(py, values)?.into_any().unbind())
}
pub(super) fn inverse(py: Python<'_>, selected: Selected) -> PyResult<PyObject> {
    poly_output(
        py,
        native::invert(&selected.a, &selected.b).map_err(err)?,
        selected.polys,
    )
}
pub(super) fn checked_poly_pair(
    a: &PolyValue,
    b: &PolyValue,
) -> PyResult<Option<(AlgebraicPoly, AlgebraicPoly)>> {
    a.checked_algebraic_pair(b).map_err(err)
}
/// The expression-level cancel: over `ZZ_I`/`QQ_I` the reference
/// `PolyElement.cancel` with its canonical denominator unit, otherwise the
/// monic-gcd cancel.
fn expression_cancel(a: &AlgebraicPoly, b: &AlgebraicPoly) -> Result<(AlgebraicPoly, AlgebraicPoly), PolyError> {
    if a.gaussian.is_some() && b.gaussian.is_some() && a.gens == b.gens && a.number_field == b.number_field {
        let (p, q) = symox::poly::gaussian::expr_cancel(&a.rep, &b.rep)?;
        return Ok((a.with_rep(p), b.with_rep(q)));
    }
    native::cancel(a, b)
}
/// The reference `Poly.cancel` over `ZZ_I`/`QQ_I` (`dmp_cancel`):
/// `(cp, cq, p, q)` with `f/g = (cp*p)/(cq*q)`.
pub(super) fn gaussian_poly_cancel(
    a: &PolyValue,
    b: &PolyValue,
) -> PyResult<Option<(symox::poly::algebraic::AN, symox::poly::algebraic::AN, AlgebraicPoly, AlgebraicPoly)>> {
    let gaussian = |p: &PolyValue| p.algebraic.as_ref().is_some_and(|a| a.gaussian.is_some());
    if !gaussian(a) && !gaussian(b) {
        return Ok(None);
    }
    let Some((a, b)) = checked_poly_pair(a, b)? else { return Ok(None) };
    let Some(ring) = a.gaussian else { return Ok(None) };
    if b.rep.is_zero() {
        return Err(err(PolyError::Unsupported("polynomial division by zero".into())));
    }
    let (cp, cq, p, q) = symox::poly::gaussian::cancel(&a.rep, &b.rep, ring).map_err(err)?;
    Ok(Some((cp, cq, a.with_rep(p), b.with_rep(q))))
}
pub(super) fn cancel_input(
    py: Python<'_>,
    f: &Bound<'_, PyAny>,
    gens: &Bound<'_, PyTuple>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Option<PyObject>> {
    if let Ok(pair) = f.downcast::<PyTuple>() {
        if pair.len() != 2 {
            return Err(PyValueError::new_err(
                "cancel expects a numerator/denominator pair",
            ));
        }
        let a = pair.get_item(0)?;
        let b = pair.get_item(1)?;
        let Some(selected) = select(&a, &b, gens, kwargs)? else {
            return Ok(None);
        };
        let (p, q) = expression_cancel(&selected.a, &selected.b).map_err(err)?;
        if selected.polys && selected.a.gaussian.is_some() {
            // The reference rebuilds each side with `Poly(P, *R.symbols)`,
            // inferring its domain again.
            let poly = py.import("symox")?.getattr("Poly")?;
            let gens: Vec<PyObject> = selected.a.gens.iter().map(|x| w(x.clone()).into_pyobject(py).map(|v| v.into_any().unbind())).collect::<PyResult<_>>()?;
            let build = |e: Expr| -> PyResult<PyObject> {
                let mut args: Vec<PyObject> = vec![w(e).into_pyobject(py)?.into_any().unbind()];
                args.extend(gens.iter().map(|g| g.clone_ref(py)));
                Ok(poly.call1(PyTuple::new(py, args)?)?.unbind())
            };
            return Ok(Some(
                PyTuple::new(py, [w(Expr::one()).into_pyobject(py)?.into_any().unbind(), build(p.as_expr())?, build(q.as_expr())?])?
                    .into_any()
                    .unbind(),
            ));
        }
        return Ok(Some(
            PyTuple::new(
                py,
                [
                    w(Expr::one()).into_pyobject(py)?.into_any().unbind(),
                    poly_output(py, p, selected.polys)?,
                    poly_output(py, q, selected.polys)?,
                ],
            )?
            .into_any()
            .unbind(),
        ));
    }
    let value = to_expr(f)?;
    let (n, d) = symox::poly::api::as_numer_denom(&value);
    let n = w(n).into_pyobject(py)?;
    let d = w(d).into_pyobject(py)?;
    let selected = if let Ok(p) = f.downcast::<PyPoly>() {
        if p.get().p.algebraic.is_some() {
            select(f, d.as_any(), gens, kwargs)?
        } else {
            select(n.as_any(), d.as_any(), gens, kwargs)?
        }
    } else {
        select(n.as_any(), d.as_any(), gens, kwargs)?
    };
    let Some(selected) = selected else {
        return Ok(None);
    };
    let (p, q) = expression_cancel(&selected.a, &selected.b).map_err(err)?;
    Ok(Some(
        w(symox::ops::div(&p.as_expr(), &q.as_expr()))
            .into_pyobject(py)?
            .into_any()
            .unbind(),
    ))
}
pub(super) fn apart_input(
    f: &Bound<'_, PyAny>,
    x: Option<&Bound<'_, PyAny>>,
    full: bool,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Option<PyExpr>> {
    let py = f.py();
    let value = to_expr(f)?;
    let (n, d) = symox::poly::api::as_numer_denom(&value);
    let n = w(n).into_pyobject(py)?;
    let d = w(d).into_pyobject(py)?;
    let gens = PyTuple::new(py, x.into_iter())?;
    let selected = if let Ok(p) = f.downcast::<PyPoly>() {
        if p.get().p.algebraic.is_some() {
            select(f, d.as_any(), &gens, kwargs)?
        } else {
            select(n.as_any(), d.as_any(), &gens, kwargs)?
        }
    } else {
        select(n.as_any(), d.as_any(), &gens, kwargs)?
    };
    let Some(selected) = selected else {
        return Ok(None);
    };
    if full {
        return Err(PyNotImplementedError::new_err(
            "full=True algebraic partial fractions need root-isolation/RootSum representations",
        ));
    }
    Ok(Some(w(native::apart(&selected.a, &selected.b)
        .map_err(err)?
        .as_expr())))
}
pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("NotInvertible", m.py().get_type::<NotInvertible>())
}

/// The reference error of the extended Euclid over a ring with `auto=False`.
fn ring_euclid_error(half: bool, domain: &str) -> PyErr {
    super::domains::DomainError::new_err(if half {
        format!("Cannot compute half extended GCD over {domain}")
    } else {
        "ground domain must be a field".to_string()
    })
}
fn univariate_expected() -> PyErr {
    PyValueError::new_err("univariate polynomial expected")
}
fn ordinary_pair(
    a: &PolyValue,
    b: &PolyValue,
    auto: bool,
    half: bool,
) -> PyResult<Vec<symox::poly::polyobj::Poly>> {
    use symox::poly::polyobj::{Domain, Poly};
    if a.field.is_some() || b.field.is_some() || a.modulus.is_some() || b.modulus.is_some() {
        return Err(PyNotImplementedError::new_err(
            "these Poly methods currently support QQ and algebraic number fields",
        ));
    }
    let gens = symox::poly::polyobj::unify_gens(a.true_gens(), b.true_gens());
    if !auto && (a.domain == Domain::ZZ && b.domain == Domain::ZZ) {
        return Err(ring_euclid_error(half, "ZZ"));
    }
    if gens.len() != 1 {
        return Err(univariate_expected());
    }
    Poly::from_exprs(&[a.as_expr(), b.as_expr()], &gens).map_err(err)
}
/// `Poly.gcdex`/`invert` over ZZ_I: `auto` moves to QQ_I, otherwise the
/// reference DomainError.
fn gaussian_field_pair(
    mut a: AlgebraicPoly,
    mut b: AlgebraicPoly,
    auto: bool,
    message: &str,
) -> PyResult<(AlgebraicPoly, AlgebraicPoly)> {
    if a.gaussian == Some(true) {
        if !auto {
            return Err(super::domains::DomainError::new_err(message.to_string()));
        }
        a.gaussian = Some(false);
        b.gaussian = Some(false);
    }
    Ok((a, b))
}
/// The pair over a rational-function coefficient domain, promoted to its
/// fraction field when `auto` is set (the reference `Poly.gcdex`).
fn fraction_poly_pair(
    a: &PolyValue,
    b: &PolyValue,
    auto: bool,
    half: bool,
) -> PyResult<Option<(symox::poly::fraction::FractionPoly, symox::poly::fraction::FractionPoly)>> {
    if a.field.is_none() && b.field.is_none() {
        return Ok(None);
    }
    let Some((mut a, mut b)) = a.checked_fraction_pair(b).map_err(err)? else { return Ok(None) };
    if !a.coefficient_domain.is_field() && !auto {
        return Err(ring_euclid_error(half, &a.coefficient_domain.to_string()));
    }
    if a.gens.len() != 1 {
        return Err(univariate_expected());
    }
    if !a.coefficient_domain.is_field() {
        a.promote_field();
        b.promote_field();
    }
    Ok(Some((a, b)))
}
/// The reference `dmp_cancel` over a rational-function coefficient domain:
/// `(cp/cq, p, q)` with the ground fraction as an expression, or
/// `(cp*p, cq*q)` with `include`.
pub(super) fn fraction_poly_cancel(
    a: &PolyValue,
    b: &PolyValue,
    include: bool,
) -> PyResult<Option<(Option<Expr>, PolyValue, PolyValue)>> {
    if a.field.is_none() && b.field.is_none() {
        return Ok(None);
    }
    let Some((a, b)) = a.checked_fraction_pair(b).map_err(err)? else { return Ok(None) };
    if b.rep.is_zero() {
        return Err(err(PolyError::Unsupported("polynomial division by zero".into())));
    }
    let (cp, cq, p, q) = a.reference_cancel(&b);
    if include {
        let _s = symox::poly::ratfunc::Scope::new(a.parameters.len());
        let p = p.with_rep(p.rep.scale(&cp));
        let q = q.with_rep(q.rep.scale(&cq));
        return Ok(Some((None, p.into(), q.into())));
    }
    let ratio = symox::ops::div(
        &symox::poly::fraction::coefficient_expr(&cp, &a.parameters),
        &symox::poly::fraction::coefficient_expr(&cq, &a.parameters),
    );
    Ok(Some((Some(ratio), p.into(), q.into())))
}
pub(super) fn poly_bezout(
    a: &PolyValue,
    b: &PolyValue,
    half: bool,
    auto: bool,
) -> PyResult<Vec<PolyValue>> {
    if let Some((a, b)) = checked_poly_pair(a, b)? {
        let (a, b) = gaussian_field_pair(a, b, auto, if half { "Cannot compute half extended GCD over ZZ_I" } else { "ground domain must be a field" })?;
        if !half && b.rep.is_zero() {
            return Err(err(PolyError::Unsupported(
                "polynomial division by zero".into(),
            )));
        }
        let (s, t, h) = native::gcdex(&a, &b).map_err(err)?;
        return Ok(if half {
            vec![s.into(), h.into()]
        } else {
            vec![s.into(), t.into(), h.into()]
        });
    }
    if let Some((a, b)) = fraction_poly_pair(a, b, auto, half)? {
        if a.rep.is_zero() && b.rep.is_zero() || !half && b.rep.is_zero() {
            return Err(err(PolyError::Unsupported("polynomial division by zero".into())));
        }
        return Ok(if half {
            let (s, h) = a.half_gcdex(&b).map_err(err)?;
            vec![s.into(), h.into()]
        } else {
            let (s, t, h) = a.gcdex(&b).map_err(err)?;
            vec![s.into(), t.into(), h.into()]
        });
    }
    let p = ordinary_pair(a, b, auto, half)?;
    if p[0].rep.is_zero() && p[1].rep.is_zero() || !half && p[1].rep.is_zero() {
        return Err(err(PolyError::Unsupported(
            "polynomial division by zero".into(),
        )));
    }
    let (s, t, h) =
        symox::poly::api::gcdex(&p[0].as_expr(), &p[1].as_expr(), &p[0].gens).map_err(err)?;
    let values = if half { vec![s, h] } else { vec![s, t, h] };
    values
        .into_iter()
        .map(|e| {
            let mut q = symox::poly::polyobj::Poly::from_expr(&e, &p[0].gens).map_err(err)?;
            q.domain = symox::poly::polyobj::Domain::QQ;
            Ok(q.into())
        })
        .collect()
}
pub(super) fn poly_inverse(a: &PolyValue, b: &PolyValue, auto: bool) -> PyResult<PolyValue> {
    if let Some((a, b)) = checked_poly_pair(a, b)? {
        let (a, b) = gaussian_field_pair(a, b, auto, "Cannot compute half extended GCD over ZZ_I")?;
        return native::invert(&a, &b).map(Into::into).map_err(err);
    }
    if let Some((a, b)) = fraction_poly_pair(a, b, auto, true)? {
        return a.invert(&b).map(Into::into).map_err(err);
    }
    let p = ordinary_pair(a, b, auto, true)?;
    let e = symox::poly::api::invert(&p[0].as_expr(), &p[1].as_expr(), &p[0].gens).map_err(err)?;
    let mut p = symox::poly::polyobj::Poly::from_expr(&e, &p[0].gens).map_err(err)?;
    p.domain = symox::poly::polyobj::Domain::QQ;
    Ok(p.into())
}
pub(super) fn poly_cancel(a: &PolyValue, b: &PolyValue) -> PyResult<(PolyValue, PolyValue)> {
    if let Some((a, b)) = checked_poly_pair(a, b)? {
        let (p, q) = native::cancel(&a, &b).map_err(err)?;
        return Ok((p.into(), q.into()));
    }
    let p = ordinary_pair(a, b, true, true)?;
    if p[1].rep.is_zero() {
        return Err(err(PolyError::Unsupported(
            "polynomial division by zero".into(),
        )));
    }
    let (n, d) = symox::poly::api::cancel_polys(&p[0].rep, &p[1].rep);
    let mut a = p[0].clone();
    let mut b = p[1].clone();
    a.rep = n;
    b.rep = d;
    Ok((a.into(), b.into()))
}
