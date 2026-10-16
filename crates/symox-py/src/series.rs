//! Python bindings for series expansions, `Order` and limits.

use pyo3::exceptions::{PyNotImplementedError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use symox::expr::Expr;
use symox::series::{self as sr, SErr};
use symox::traverse::free_symbols;

use crate::{to_expr, PyExpr};

fn w(e: Expr) -> PyExpr {
    PyExpr { e }
}

fn serr(e: SErr) -> PyErr {
    match e {
        SErr::Pole(s) => PyValueError::new_err(format!("PoleError: {s}")),
        e => PyNotImplementedError::new_err(e.to_string()),
    }
}

fn the_symbol(e: &Expr, x: Option<&Bound<'_, PyAny>>) -> PyResult<Expr> {
    match x {
        Some(x) if !x.is_none() => to_expr(x),
        _ => {
            let fs = free_symbols(e);
            match fs.len() {
                1 => Ok(fs[0].clone()),
                0 => Ok(Expr::symbol("x")),
                _ => Err(PyValueError::new_err("x must be given for multivariate functions.")),
            }
        }
    }
}

/// `series(expr, x=None, x0=0, n=6, dir='+')`
#[pyfunction]
#[pyo3(signature = (expr, x=None, x0=None, n=None, dir="+"))]
fn _series(
    expr: &Bound<'_, PyAny>,
    x: Option<&Bound<'_, PyAny>>,
    x0: Option<&Bound<'_, PyAny>>,
    n: Option<i64>,
    dir: &str,
) -> PyResult<PyExpr> {
    let e = to_expr(expr)?;
    let x = the_symbol(&e, x)?;
    let x0 = match x0 {
        Some(v) if !v.is_none() => to_expr(v)?,
        _ => Expr::zero(),
    };
    let minus = match dir {
        "+" => false,
        "-" => true,
        _ => return Err(PyValueError::new_err("Dir must be '+' or '-'")),
    };
    sr::series(&e, &x, &x0, n.unwrap_or(6), minus).map(w).map_err(serr)
}

/// `Order(expr, *args)`; args are symbols or `(symbol, point)` tuples.
#[pyfunction]
#[pyo3(signature = (expr, *args))]
fn _order(expr: &Bound<'_, PyAny>, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let e = to_expr(expr)?;
    let mut vars = Vec::new();
    let mut pending_syms: Vec<Expr> = Vec::new();
    for a in args.iter() {
        if let Ok(t) = a.downcast::<PyTuple>() {
            if t.len() != 2 {
                return Err(PyValueError::new_err("expected (symbol, point) tuples"));
            }
            vars.push((to_expr(&t.get_item(0)?)?, to_expr(&t.get_item(1)?)?));
        } else {
            pending_syms.push(to_expr(&a)?);
        }
    }
    // O(x**2, x, y): all at the point given last (default 0)
    for s in pending_syms {
        vars.push((s, Expr::zero()));
    }
    Ok(w(sr::order(&e, &vars)))
}

#[pyfunction]
fn _remove_o(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(w(sr::remove_o(&to_expr(e)?)))
}

#[pyfunction]
fn _get_o(e: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    Ok(sr::get_o(&to_expr(e)?).map(w))
}

#[pyfunction]
fn _getn(e: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    Ok(sr::getn(&to_expr(e)?).map(w))
}

#[pyfunction]
fn _is_order(e: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(sr::is_order(&to_expr(e)?))
}

#[pyfunction]
fn _order_contains(o: &Bound<'_, PyAny>, t: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(sr::contains(&to_expr(o)?, &to_expr(t)?))
}

#[pyfunction]
fn _taylor_term(e: &Bound<'_, PyAny>, n: i64, x: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    sr::taylor_term(&to_expr(e)?, &to_expr(x)?, n).map(w).map_err(serr)
}

#[pyfunction]
#[pyo3(signature = (e, x, x0=None, dir="+"))]
fn _leadterm<'py>(
    py: Python<'py>,
    e: &Bound<'py, PyAny>,
    x: &Bound<'py, PyAny>,
    x0: Option<&Bound<'py, PyAny>>,
    dir: &str,
) -> PyResult<Bound<'py, PyTuple>> {
    let x0 = match x0 {
        Some(v) if !v.is_none() => to_expr(v)?,
        _ => Expr::zero(),
    };
    let lt = sr::leading_term(&to_expr(e)?, &to_expr(x)?, &x0, dir == "-").map_err(serr)?;
    PyTuple::new(py, [w(lt.coeff), w(Expr::num(lt.exp))])
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_series, m)?)?;
    m.add_function(wrap_pyfunction!(_order, m)?)?;
    m.add_function(wrap_pyfunction!(_remove_o, m)?)?;
    m.add_function(wrap_pyfunction!(_get_o, m)?)?;
    m.add_function(wrap_pyfunction!(_getn, m)?)?;
    m.add_function(wrap_pyfunction!(_is_order, m)?)?;
    m.add_function(wrap_pyfunction!(_order_contains, m)?)?;
    m.add_function(wrap_pyfunction!(_taylor_term, m)?)?;
    m.add_function(wrap_pyfunction!(_leadterm, m)?)?;
    Ok(())
}
