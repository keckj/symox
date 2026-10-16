//! Python bindings for `solve` and `nsolve`.

use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PySet, PyTuple};

use symox::solve::{self as sv, SystemSolution};
use symox::Expr;

use crate::{to_expr, wrap};

fn set_result(py: Python<'_>, variables: &[Expr], solutions: Vec<Vec<Expr>>) -> PyResult<PyObject> {
    let symbols = PyList::new(py, variables.iter().cloned().map(wrap))?;
    let rows = solutions.into_iter().map(|row| PyTuple::new(py, row.into_iter().map(wrap))).collect::<PyResult<Vec<_>>>()?;
    Ok((symbols, PySet::new(py, rows)?).into_pyobject(py)?.into_any().unbind())
}

fn solve_flags(flags: Option<&Bound<'_, PyDict>>) -> PyResult<(bool, bool)> {
    let (mut dict, mut set) = (false, false);
    if let Some(flags) = flags {
        for (key, value) in flags.iter() {
            let key: String = key.extract()?;
            match key.as_str() {
                "dict" => dict = value.is_truthy()?,
                "set" => set = value.is_truthy()?,
                "check" | "simplify" | "cubics" | "quartics" | "quintics" if value.is_truthy()? => {},
                "manual" | "warn" | "minimal" if !value.is_truthy()? => {},
                _ => return Err(PyNotImplementedError::new_err(format!("solve option '{key}' with this value is not implemented"))),
            }
        }
    }
    Ok((dict, set && !dict))
}

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
    let (as_dict, as_set) = solve_flags(flags)?;
    let (eqs, is_system) = expr_list(f)?;
    let mut vars: Vec<Expr> = Vec::new();
    for s in symbols.iter() {
        let (v, _) = expr_list(&s)?;
        vars.extend(v);
    }
    if vars.is_empty() {
        if as_set { return set_result(py, &vars, vec![]); }
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
        if as_set { return set_result(py, &vars, sols.into_iter().map(|s| vec![s]).collect()); }
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
        SystemSolution::NoSolution if as_set => set_result(py, &vars, vec![]),
        SystemSolution::NoSolution => Ok(PyList::empty(py).into_any().unbind()),
        SystemSolution::Map(m) => {
            if as_set {
                let row = vars.iter().map(|v| m.iter().find(|(k, _)| k == v).map(|(_, value)| value.clone()).unwrap_or_else(|| v.clone())).collect();
                return set_result(py, &vars, vec![row]);
            }
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
            if as_set { return set_result(py, &vars, ts); }
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
    let (mut tolerance, mut maxsteps, mut verify, mut as_dict) = (1e-14_f64, 100_usize, true, false);
    let mut solver = "newton".to_owned();
    if let Some(kwargs) = _kwargs {
        for (key, value) in kwargs.iter() {
            let key: String = key.extract()?;
            match key.as_str() {
                "tol" => tolerance = value.extract()?,
                "maxsteps" => maxsteps = value.extract()?,
                "verify" => verify = value.extract()?,
                "dict" => as_dict = value.is_truthy()?,
                "solver" => solver = value.extract()?,
                "prec" if (1..=15).contains(&value.extract::<u32>()?) => {},
                "verbose" if !value.is_truthy()? => {},
                _ => return Err(PyNotImplementedError::new_err(format!("nsolve option '{key}' with this value is not implemented"))),
            }
        }
    }
    if !tolerance.is_finite() || tolerance <= 0.0 || maxsteps == 0 {
        return Err(PyValueError::new_err("tol must be finite and positive, and maxsteps must be positive"));
    }
    let (eqs, sys) = expr_list(f)?;
    let (vars, _) = expr_list(x)?;
    let starts: Vec<f64> = if let Ok(v) = x0.extract::<f64>() {
        vec![v]
    } else {
        x0.extract::<Vec<f64>>()?
    };
    if vars.iter().any(|v| !v.is_symbol()) || vars.iter().enumerate().any(|(i, v)| vars[..i].contains(v)) {
        return Err(PyValueError::new_err("nsolve unknowns must be distinct symbols"));
    }
    let r = if solver == "bisect" {
        if vars.len() != 1 || eqs.len() != 1 || starts.len() != 2 {
            return Err(PyValueError::new_err("bisection requires one equation, one symbol and a two-endpoint bracket"));
        }
        sv::nsolve_bisect(&eqs[0], &vars[0], starts[0], starts[1], tolerance, maxsteps, verify).map(|r| vec![r])
    } else if solver == "newton" {
        if vars.len() != starts.len() || vars.len() != eqs.len() {
        return Err(PyValueError::new_err("nsolve needs as many equations, unknowns and starting points"));
        }
        sv::nsolve_with_options(&eqs, &vars, &starts, tolerance, maxsteps, verify)
    } else {
        return Err(PyNotImplementedError::new_err(format!("nsolve solver '{solver}' is not implemented")));
    }.ok_or_else(|| PyValueError::new_err("Could not find root within given tolerance"))?;
    if as_dict {
        let mapping = PyDict::new(py);
        for (symbol, value) in vars.iter().zip(r) { mapping.set_item(wrap(symbol.clone()), wrap(Expr::float(value)))?; }
        return Ok(PyList::new(py, [mapping])?.into_any().unbind());
    }
    if !sys && r.len() == 1 {
        return Ok(wrap(Expr::float(r[0])).into_pyobject(py)?.into_any().unbind());
    }
    let values = PyList::new(py, r.into_iter().map(|v| wrap(Expr::float(v))))?;
    Ok(py.import("symox")?.getattr("Matrix")?.call1((vars.len(), 1, values))?.unbind())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    m.add_function(wrap_pyfunction!(nsolve, m)?)?;
    Ok(())
}
