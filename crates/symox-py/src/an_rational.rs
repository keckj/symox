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
    let domain = kwargs
        .map(|k| k.get_item("domain"))
        .transpose()?
        .flatten()
        .filter(|v| !v.is_none());
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
    if selected.native.as_ref().is_some_and(|d|!matches!(d.inner,symox::poly::domains::ExactDomain::ZZ|symox::poly::domains::ExactDomain::QQ|symox::poly::domains::ExactDomain::Algebraic(_))){return Err(PyNotImplementedError::new_err("this rational operation currently supports QQ and certified algebraic domains"));}
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
    let mut generators = gens_of(gens)?;
    if generators.is_empty() {
        generators = symox::poly::polyobj::unify_gens(&fg, &gg);
    }
    if generators.is_empty() {
        generators = symox::poly::fraction::infer_parameters(&[fe.clone(), ge.clone()], &[]);
    }
    if generators.is_empty() {
        generators.push(Expr::symbol("_x"));
    }
    if generators.len() != 1 {
        return Err(PyNotImplementedError::new_err(
            "algebraic rational operations require one polynomial generator",
        ));
    }
    let polys = kwargs
        .map(|k| k.get_item("polys"))
        .transpose()?
        .flatten()
        .map(|v| v.extract::<bool>())
        .transpose()?
        .unwrap_or_else(|| f.downcast::<PyPoly>().is_ok() || g.downcast::<PyPoly>().is_ok());
    let a = AlgebraicPoly::from_expr(&fe, &generators, field.clone()).map_err(err)?;
    let b = AlgebraicPoly::from_expr(&ge, &generators, field).map_err(err)?;
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
        let (p, q) = native::cancel(&selected.a, &selected.b).map_err(err)?;
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
    let (p, q) = native::cancel(&selected.a, &selected.b).map_err(err)?;
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

fn ordinary_pair(
    a: &PolyValue,
    b: &PolyValue,
    auto: bool,
) -> PyResult<Vec<symox::poly::polyobj::Poly>> {
    use symox::poly::polyobj::{Domain, Poly};
    if a.field.is_some() || b.field.is_some() || a.modulus.is_some() || b.modulus.is_some() {
        return Err(PyNotImplementedError::new_err(
            "these Poly methods currently support QQ and algebraic number fields",
        ));
    }
    if !auto && (a.domain == Domain::ZZ || b.domain == Domain::ZZ) {
        return Err(PyValueError::new_err(
            "extended polynomial Euclid requires a coefficient field; use auto=True",
        ));
    }
    let gens = symox::poly::polyobj::unify_gens(a.true_gens(), b.true_gens());
    if gens.len() != 1 {
        return Err(PyNotImplementedError::new_err(
            "extended polynomial Euclid requires one generator",
        ));
    }
    Poly::from_exprs(&[a.as_expr(), b.as_expr()], &gens).map_err(err)
}
pub(super) fn poly_bezout(
    a: &PolyValue,
    b: &PolyValue,
    half: bool,
    auto: bool,
) -> PyResult<Vec<PolyValue>> {
    if let Some((a, b)) = checked_poly_pair(a, b)? {
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
    let p = ordinary_pair(a, b, auto)?;
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
        return native::invert(&a, &b).map(Into::into).map_err(err);
    }
    let p = ordinary_pair(a, b, auto)?;
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
    let p = ordinary_pair(a, b, true)?;
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
