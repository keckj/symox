//! `_derivative_orders` and `_derivative_order_plan`: the native variable
//! parsing and order plan of the symbolic-order derivatives of matrices,
//! matrix expressions and arrays (`tensor/array/array_derivatives.py`).
//! `NotImplemented` means the kernel declined and the Python code runs.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyInt, PyList, PyTuple};
use symox::diff::order_plan::{self as op, OrderItem, OrderPlan};
use symox::{Expr, Kind};

fn native(obj: &Bound<'_, PyAny>) -> Option<Expr> {
    if obj.is_instance_of::<crate::PyExprNode>() {
        crate::to_expr(obj).ok()
    } else {
        None
    }
}

fn integer(e: &Expr) -> Option<num_bigint::BigInt> {
    match e.kind() {
        Kind::Num(n) if n.is_integer() => n.to_bigint(),
        _ => None,
    }
}

/// One argument as `_parse` reads it; `Err(true)` when `_parse` gives
/// None, `Err(false)` when the Python code must decide.
fn item(v: &Bound<'_, PyAny>) -> Result<OrderItem, bool> {
    if v.is_instance_of::<PyInt>() || v.is_instance_of::<PyBool>() {
        return v.extract::<num_bigint::BigInt>().map(OrderItem::Order).map_err(|_| false);
    }
    let one = || Expr::one();
    if let Some(e) = native(v) {
        if let Some(n) = integer(&e) {
            return Ok(OrderItem::Order(n));
        }
        if !v.getattr("is_Tuple").and_then(|t| t.is_truthy()).map_err(|_| false)? {
            return Ok(OrderItem::Variable(e, one()));
        }
    } else if !v.is_instance_of::<PyTuple>() {
        return Err(false);
    }
    // a (variable, count) pair
    let n = v.len().map_err(|_| false)?;
    if n != 2 {
        return Err(true);
    }
    let (variable, count) = (v.get_item(0).map_err(|_| false)?, v.get_item(1).map_err(|_| false)?);
    let count = if count.is_instance_of::<PyInt>() && !count.is_instance_of::<PyBool>() {
        Expr::num(symox::Number::from_bigint(count.extract().map_err(|_| false)?))
    } else {
        native(&count).ok_or(false)?
    };
    Ok(OrderItem::Variable(native(&variable).ok_or(false)?, count))
}

fn pair_list<'py>(py: Python<'py>, pairs: Vec<(Expr, Expr)>) -> PyResult<Bound<'py, PyList>> {
    let items = pairs
        .into_iter()
        .map(|(v, c)| PyTuple::new(py, [crate::expr_object(py, v)?, crate::expr_object(py, c)?]))
        .collect::<PyResult<Vec<_>>>()?;
    PyList::new(py, items)
}

/// `_parse(variables)`: the merged `(Symbol, count)` pairs, None, or
/// NotImplemented; raises the ValueError of a negative order.
#[pyfunction]
fn _derivative_orders(py: Python<'_>, variables: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
    let mut items = Vec::with_capacity(variables.len());
    for v in variables.iter() {
        match item(&v) {
            Ok(i) => items.push(i),
            Err(true) => return Ok(py.None()),
            Err(false) => return Ok(py.NotImplemented()),
        }
    }
    match op::parse_orders(items) {
        Ok(Some(pairs)) => Ok(pair_list(py, pairs)?.into_any().unbind()),
        Ok(None) => Ok(py.None()),
        Err(message) => Err(PyValueError::new_err(message)),
    }
}

/// The plan of `_evaluate`: True for the zero exit, `(ordered, split)`
/// otherwise (the Integer orders of `ordered[:split]` are derived, the rest
/// held), None without pairs, NotImplemented when declined.
#[pyfunction]
fn _derivative_order_plan(py: Python<'_>, pairs: &Bound<'_, PyAny>, free: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let mut native_pairs = Vec::new();
    for pair in pairs.try_iter()? {
        let pair = pair?;
        let (Some(v), Some(c)) = (native(&pair.get_item(0)?), native(&pair.get_item(1)?)) else {
            return Ok(py.NotImplemented());
        };
        native_pairs.push((v, c));
    }
    let mut symbols = Vec::new();
    for s in free.try_iter()? {
        match native(&s?) {
            Some(s) => symbols.push(s),
            None => return Ok(py.NotImplemented()),
        }
    }
    Ok(match op::plan_orders(native_pairs, &symbols) {
        OrderPlan::Unchanged => py.None(),
        OrderPlan::Zero => PyBool::new(py, true).to_owned().into_any().unbind(),
        OrderPlan::Derive(ordered, split) => (pair_list(py, ordered)?, split).into_pyobject(py)?.into_any().unbind(),
    })
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_derivative_orders, m)?)?;
    m.add_function(wrap_pyfunction!(_derivative_order_plan, m)?)
}
