//! Exact rational boundary for real-root analysis. No runtime SymPy calls.
use super::{division, domains, poly_obj, w, PolyValue, PyExpr, PyPoly};
use crate::to_expr;
use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};
use symox::expr::{Constant, Expr};
use symox::poly::real_roots::{
    self as native, Bound as RootBound, IsolationOptions, RefinementOptions, RootError,
};
use symox::poly::Rat;
pyo3::create_exception!(symox, RefinementFailed, PyValueError);

fn error(error: RootError) -> PyErr {
    match error {
        RootError::Domain(message) => domains::DomainError::new_err(message),
        RootError::Refinement(message) => RefinementFailed::new_err(message),
        RootError::Invalid(message) | RootError::Limit(message) => PyValueError::new_err(message),
    }
}
fn rational(value: &Bound<'_, PyAny>) -> PyResult<Rat> {
    let expression = to_expr(value)?;
    let number = expression.as_number().ok_or_else(|| {
        domains::CoercionFailed::new_err("expected a finite rational real-root bound or precision")
    })?;
    if matches!(number,symox::Number::Float(f) if !f.is_finite()) {
        return Err(domains::CoercionFailed::new_err(
            "real-root bounds must be finite rationals",
        ));
    }
    // Float arguments denote their exact stored binary rational; they never
    // supply a root estimate or a decision about a polynomial's sign.
    Ok(symox::poly::convert::number_to_rat(number))
}
fn optional_rational(value: Option<&Bound<'_, PyAny>>) -> PyResult<Option<Rat>> {
    value.filter(|v| !v.is_none()).map(rational).transpose()
}
fn bound(value: Option<&Bound<'_, PyAny>>, lower: bool) -> PyResult<RootBound> {
    let Some(value) = value.filter(|v| !v.is_none()) else {
        return Ok(if lower {
            RootBound::NegativeInfinity
        } else {
            RootBound::PositiveInfinity
        });
    };
    let expression = to_expr(value)?;
    if expression.is_constant(Constant::NegInfinity) {
        return Ok(RootBound::NegativeInfinity);
    }
    if expression.is_constant(Constant::Infinity) {
        return Ok(RootBound::PositiveInfinity);
    }
    if symox::assumptions::is_real(&expression) == Some(false) {
        return Err(PyNotImplementedError::new_err(
            "complex rectangle root counting is not implemented; use real rational bounds",
        ));
    }
    Ok(RootBound::Rational(rational(value)?))
}
fn isolation_bound(value: Option<&Bound<'_, PyAny>>, lower: bool) -> PyResult<Option<Rat>> {
    match bound(value, lower)? {
        RootBound::Rational(value) => Ok(Some(value)),
        RootBound::NegativeInfinity if lower => Ok(None),
        RootBound::PositiveInfinity if !lower => Ok(None),
        _ => Err(PyValueError::new_err(
            "isolation bounds must be in increasing order",
        )),
    }
}
fn unary<'py>(f: &Bound<'py, PyAny>) -> PyResult<PolyValue> {
    let gens = PyTuple::empty(f.py());
    Ok(division::prepare_unary(f, &gens, None)?.0)
}
pub(super) fn sturm_method(p: &PolyValue, auto: bool) -> PyResult<Vec<PyPoly>> {
    Ok(native::sturm(p, auto)
        .map_err(error)?
        .into_iter()
        .map(|p| PyPoly { p })
        .collect())
}
pub(super) fn count_method(
    p: &PolyValue,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExpr> {
    Ok(w(Expr::int(
        native::count_roots(p, bound(inf, true)?, bound(sup, false)?).map_err(error)? as i64,
    )))
}
pub(super) fn intervals_method(
    py: Python<'_>,
    p: &PolyValue,
    all: bool,
    eps: Option<&Bound<'_, PyAny>>,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
    fast: bool,
    sqf: bool,
) -> PyResult<PyObject> {
    if all {
        return Err(PyNotImplementedError::new_err(
            "complex root isolation (all=True) is not implemented",
        ));
    }
    let options = IsolationOptions {
        epsilon: optional_rational(eps)?,
        lower: isolation_bound(inf, true)?,
        upper: isolation_bound(sup, false)?,
        fast,
    };
    let result = PyList::empty(py);
    for interval in native::intervals(p, &options).map_err(error)? {
        let endpoints = PyTuple::new(
            py,
            [
                w(Expr::num(interval.lower.to_number())),
                w(Expr::num(interval.upper.to_number())),
            ],
        )?
        .into_any()
        .unbind();
        if sqf {
            result.append(endpoints)?;
        } else {
            result.append(PyTuple::new(
                py,
                [
                    endpoints,
                    interval.multiplicity.into_pyobject(py)?.into_any().unbind(),
                ],
            )?)?;
        }
    }
    Ok(result.into_any().unbind())
}
pub(super) fn refine_method(
    p: &PolyValue,
    s: &Bound<'_, PyAny>,
    t: &Bound<'_, PyAny>,
    eps: Option<&Bound<'_, PyAny>>,
    steps: Option<i64>,
    fast: bool,
    check_sqf: bool,
) -> PyResult<(PyExpr, PyExpr)> {
    let options = RefinementOptions {
        epsilon: optional_rational(eps)?,
        steps: steps.map(|n| n.max(0) as u64),
        fast,
        check_squarefree: check_sqf,
    };
    let (a, b) = native::refine_root(p, rational(s)?, rational(t)?, &options).map_err(error)?;
    Ok((w(Expr::num(a.to_number())), w(Expr::num(b.to_number()))))
}
#[pyfunction]
#[pyo3(signature=(f,*gens,**kwargs))]
fn sturm<'py>(
    f: &Bound<'py, PyAny>,
    gens: &Bound<'py, PyTuple>,
    kwargs: Option<&Bound<'py, PyDict>>,
) -> PyResult<PyObject> {
    let options = PyDict::new(f.py());
    let mut auto = true;
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs.iter() {
            if key.extract::<String>()? == "auto" {
                auto = division::boolean(&value, "auto")?.unwrap_or(true);
            } else {
                options.set_item(key, value)?;
            }
        }
    }
    let (p, polys) = division::prepare_unary(f, gens, Some(&options))?;
    let values = native::sturm(&p, auto).map_err(error)?;
    Ok(PyList::new(
        f.py(),
        values
            .into_iter()
            .map(|p| {
                if polys {
                    poly_obj(f.py(), p)
                } else {
                    Ok(w(p.as_expr()).into_pyobject(f.py())?.into_any().unbind())
                }
            })
            .collect::<PyResult<Vec<_>>>()?,
    )?
    .into_any()
    .unbind())
}
#[pyfunction]
#[pyo3(signature=(f,inf=None,sup=None))]
fn count_roots(
    f: &Bound<'_, PyAny>,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExpr> {
    count_method(&unary(f)?, inf, sup)
}
#[pyfunction]
#[pyo3(signature=(F,all=false,eps=None,inf=None,sup=None,strict=false,fast=false,sqf=false))]
#[allow(non_snake_case)]
fn intervals(
    F: &Bound<'_, PyAny>,
    all: bool,
    eps: Option<&Bound<'_, PyAny>>,
    inf: Option<&Bound<'_, PyAny>>,
    sup: Option<&Bound<'_, PyAny>>,
    strict: bool,
    fast: bool,
    sqf: bool,
) -> PyResult<PyObject> {
    if F.is_instance_of::<PyList>() || F.is_instance_of::<PyTuple>() || strict {
        return Err(PyNotImplementedError::new_err("simultaneous root isolation and strict=True are not implemented; supply one polynomial"));
    }
    let p = if F.downcast::<PyPoly>().is_err() {
        let expression = to_expr(F)?;
        if expression.as_number().is_some_and(|n| n.is_rational()) {
            symox::poly::polyobj::Poly::from_expr(
                &expression,
                &[Expr::symbol("_root_interval_generator")],
            )
            .map_err(super::err)?
            .into()
        } else {
            unary(F)?
        }
    } else {
        unary(F)?
    };
    intervals_method(F.py(), &p, all, eps, inf, sup, fast, sqf)
}
#[pyfunction]
#[pyo3(signature=(f,s,t,eps=None,steps=None,fast=false,check_sqf=false))]
fn refine_root(
    f: &Bound<'_, PyAny>,
    s: &Bound<'_, PyAny>,
    t: &Bound<'_, PyAny>,
    eps: Option<&Bound<'_, PyAny>>,
    steps: Option<i64>,
    fast: bool,
    check_sqf: bool,
) -> PyResult<(PyExpr, PyExpr)> {
    refine_method(&unary(f)?, s, t, eps, steps, fast, check_sqf)
}
pub(super) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("RefinementFailed", m.py().get_type::<RefinementFailed>())?;
    m.add_function(wrap_pyfunction!(sturm, m)?)?;
    m.add_function(wrap_pyfunction!(count_roots, m)?)?;
    m.add_function(wrap_pyfunction!(intervals, m)?)?;
    m.add_function(wrap_pyfunction!(refine_root, m)?)?;
    Ok(())
}
