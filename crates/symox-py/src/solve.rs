//! Python bindings for `solve` and `nsolve`.

use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::solve::{self as sv, SystemSolution};
use symox::Expr;

use crate::{to_expr, wrap};

fn expr_list(o: &Bound<'_, PyAny>) -> PyResult<(Vec<Expr>, bool)> {
    if let Ok(e) = to_expr(o) {
        return Ok((vec![e], false));
    }
    let mut v = Vec::new();
    for item in o.try_iter()? {
        v.push(to_expr(&item?)?);
    }
    Ok((v, true))
}

#[pyfunction]
#[pyo3(signature = (f, *symbols, **flags))]
fn solve(py: Python<'_>, f: &Bound<'_, PyAny>, symbols: &Bound<'_, PyTuple>, flags: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let as_dict = match flags {
        Some(d) => d.get_item("dict")?.map(|v| v.is_truthy()).transpose()?.unwrap_or(false),
        None => false,
    };
    let (eqs, is_system) = expr_list(f)?;
    let mut vars: Vec<Expr> = Vec::new();
    for s in symbols.iter() {
        let (v, _) = expr_list(&s)?;
        vars.extend(v);
    }
    if vars.is_empty() {
        vars = sv::unknowns(&eqs);
    }
    if vars.is_empty() {
        return Ok(PyList::empty(py).into_any().unbind());
    }
    if !is_system && eqs.len() == 1 && vars.len() == 1 {
        let x = &vars[0];
        let sols = sv::solve(&eqs[0], x).ok_or_else(|| {
            PyNotImplementedError::new_err(format!("could not solve {} for {}", eqs[0], x))
        })?;
        if as_dict {
            let l = PyList::empty(py);
            for s in sols {
                let d = PyDict::new(py);
                d.set_item(wrap(x.clone()), wrap(s))?;
                l.append(d)?;
            }
            return Ok(l.into_any().unbind());
        }
        return Ok(PyList::new(py, sols.into_iter().map(wrap))?.into_any().unbind());
    }
    let res = sv::solve_system(&eqs, &vars)
        .ok_or_else(|| PyNotImplementedError::new_err("could not solve the system"))?;
    match res {
        SystemSolution::NoSolution => Ok(PyList::empty(py).into_any().unbind()),
        SystemSolution::Map(m) => {
            let d = PyDict::new(py);
            for (k, v) in m {
                d.set_item(wrap(k), wrap(v))?;
            }
            if as_dict {
                return Ok(PyList::new(py, [d])?.into_any().unbind());
            }
            Ok(d.into_any().unbind())
        }
        SystemSolution::Tuples(ts) => {
            let l = PyList::empty(py);
            for t in ts {
                if as_dict {
                    let d = PyDict::new(py);
                    for (k, v) in vars.iter().zip(t) {
                        d.set_item(wrap(k.clone()), wrap(v))?;
                    }
                    l.append(d)?;
                } else {
                    l.append(PyTuple::new(py, t.into_iter().map(wrap))?)?;
                }
            }
            Ok(l.into_any().unbind())
        }
    }
}

#[pyfunction]
#[pyo3(signature = (f, x, x0, **_kwargs))]
fn nsolve(py: Python<'_>, f: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, x0: &Bound<'_, PyAny>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (eqs, sys) = expr_list(f)?;
    let (vars, _) = expr_list(x)?;
    let starts: Vec<f64> = if let Ok(v) = x0.extract::<f64>() {
        vec![v]
    } else {
        x0.extract::<Vec<f64>>()?
    };
    if vars.len() != starts.len() || vars.len() != eqs.len() {
        return Err(PyValueError::new_err("nsolve needs as many equations, unknowns and starting points"));
    }
    let r = sv::nsolve(&eqs, &vars, &starts).ok_or_else(|| PyValueError::new_err("Could not find root within given tolerance"))?;
    if !sys && r.len() == 1 {
        return Ok(wrap(Expr::float(r[0])).into_pyobject(py)?.into_any().unbind());
    }
    Ok(PyList::new(py, r.into_iter().map(|v| wrap(Expr::float(v))))?.into_any().unbind())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    m.add_function(wrap_pyfunction!(nsolve, m)?)?;
    Ok(())
}
