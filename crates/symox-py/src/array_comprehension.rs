//! `ArrayComprehension` and `ArrayComprehensionMap`
//! (sympy/tensor/array/array_comprehension.py): limit validation, shape,
//! loop size, free symbols and the expansion of the entries. The Python
//! classes only hold the arguments and pass them here.

use num_bigint::BigInt;
use num_traits::{One, ToPrimitive};
use pyo3::exceptions::{PyAttributeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PySet, PyTuple};
use symox::expr::{Expr, Kind, RelOp};
use symox::special::{tuple, tuple_items};

use crate::{checked_subs, expr_object, matrix_or_expr, non_expr, to_expr};

/// The positions an expansion may create, as for dense arrays.
const MAX_ENTRIES: usize = 100_000;

/// `Basic.atoms()` of an expression: its leaves.
fn visit_atoms(e: &Expr, f: &mut dyn FnMut(&Expr) -> bool) -> bool {
    let args = e.args();
    if args.is_empty() { return f(e); }
    args.iter().all(|a| visit_atoms(a, f))
}

/// `i.atoms(Symbol, Integer) != i.atoms()`: a leaf other than a symbol or
/// an integer.
fn bound_atoms_valid(e: &Expr) -> bool {
    visit_atoms(e, &mut |leaf| leaf.is_symbol() || leaf.is_integer_number())
}

/// `Basic(inf, sup).atoms(Symbol)` is not empty.
fn has_symbol_atom(e: &Expr) -> bool {
    !visit_atoms(e, &mut |leaf| !leaf.is_symbol())
}

fn limit_parts(limit: &Expr) -> PyResult<(Expr, Expr, Expr)> {
    match tuple_items(limit) {
        Some([v, a, b]) => Ok((v.clone(), a.clone(), b.clone())),
        _ => Err(PyValueError::new_err("ArrayComprehension requires values lower and upper bound for the expression")),
    }
}

fn limits_of(limits: &Bound<'_, PyTuple>) -> PyResult<Vec<(Expr, Expr, Expr)>> {
    limits.iter().map(|l| limit_parts(&to_expr(&l)?)).collect()
}

/// `ArrayComprehension._check_limits_validity`: `limits` holds the
/// sympified `(var, inf, sup)` triples (a list `sup` already a Tuple);
/// returns them as `Tuple`s.
#[pyfunction]
fn _array_comprehension_limits(py: Python<'_>, limits: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
    let mut out = Vec::with_capacity(limits.len());
    for limit in limits.iter() {
        let items = limit.downcast::<PyTuple>()?;
        let (var, inf, sup) = (to_expr(&items.get_item(0)?)?, to_expr(&items.get_item(1)?)?, to_expr(&items.get_item(2)?)?);
        out.push(tuple(vec![var.clone(), inf.clone(), sup.clone()]));
        if [&inf, &sup].iter().any(|b| non_expr(b).is_some() || !bound_atoms_valid(b)) {
            return Err(PyTypeError::new_err("Bounds should be an Expression(combination of Integer and Symbol)"));
        }
        if matches!(symox::relational::rel(RelOp::Gt, &inf, &sup).kind(), Kind::Bool(true)) {
            return Err(PyValueError::new_err("Lower bound should be inferior to upper bound"));
        }
        let free = |e: &Expr| symox::traverse::free_symbols(e).contains(&var);
        if free(&inf) || free(&sup) {
            return Err(PyValueError::new_err("Variable should not be part of its bounds"));
        }
    }
    Ok(PyTuple::new(py, out.into_iter().map(|e| expr_object(py, e)).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind())
}

/// `_calculate_shape_from_limits`: `sup - inf + 1` per limit.
fn shape(limits: &[(Expr, Expr, Expr)]) -> Vec<Expr> {
    limits.iter().map(|(_, a, b)| symox::ops::add2(&symox::ops::sub(b, a), &Expr::one())).collect()
}

#[pyfunction]
fn _array_comprehension_shape(py: Python<'_>, limits: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
    let shape = shape(&limits_of(limits)?);
    Ok(PyTuple::new(py, shape.into_iter().map(|e| expr_object(py, e)).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind())
}

/// `__len__`: the loop size, the product of the shape (the int 0, which
/// has no `free_symbols`, without limits).
#[pyfunction]
fn _array_comprehension_len(limits: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
    let shape = shape(&limits_of(limits)?);
    if shape.is_empty() {
        return Err(PyAttributeError::new_err("'int' object has no attribute 'free_symbols'"));
    }
    let size = shape.iter().fold(Expr::one(), |p, l| symox::ops::mul2(&p, l));
    if !symox::traverse::free_symbols(&size).is_empty() {
        return Err(PyValueError::new_err("Symbolic length is not supported"));
    }
    let py = limits.py();
    let n = size.as_number().and_then(|n| n.to_bigint())
        .ok_or_else(|| PyTypeError::new_err(format!("'{}' object cannot be interpreted as an integer", size.class_name())))?;
    Ok(n.into_pyobject(py)?.into_any().unbind())
}

/// `is_shape_numeric`: no symbol among the atoms of the bounds.
#[pyfunction]
fn _array_comprehension_is_shape_numeric(limits: &Bound<'_, PyTuple>) -> PyResult<bool> {
    Ok(!limits_of(limits)?.iter().any(|(_, a, b)| has_symbol_atom(a) || has_symbol_atom(b)))
}

/// `free_symbols`: those of `function` (for `ArrayComprehensionMap`, the
/// first limit), each variable removed in turn and the bounds' added.
#[pyfunction]
fn _array_comprehension_free_symbols(py: Python<'_>, function: &Bound<'_, PyAny>, limits: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
    let symbols = PySet::empty(py)?;
    match to_expr(function) {
        Ok(f) => for s in symox::traverse::free_symbols(&f) { symbols.add(expr_object(py, s)?)?; },
        Err(_) => for s in function.getattr("free_symbols")?.try_iter()? { symbols.add(s?)?; },
    }
    for (var, inf, sup) in limits_of(limits)? {
        symbols.discard(expr_object(py, var)?)?;
        for s in symox::traverse::free_symbols(&inf).into_iter().chain(symox::traverse::free_symbols(&sup)) {
            symbols.add(expr_object(py, s)?)?;
        }
    }
    Ok(symbols.into_any().unbind())
}

/// The integer ranges `range(inf, sup + 1)` of numeric limits.
fn ranges(limits: &[(Expr, Expr, Expr)]) -> PyResult<Vec<(BigInt, usize)>> {
    let integer = |e: &Expr| e.as_number().filter(|n| n.is_integer()).and_then(|n| n.to_bigint())
        .ok_or_else(|| PyTypeError::new_err(format!("'{}' object cannot be interpreted as an integer", e.class_name())));
    let mut total: usize = 1;
    let mut out = Vec::with_capacity(limits.len());
    for (_, a, b) in limits {
        let (a, b) = (integer(a)?, integer(b)?);
        let width = if b < a { 0 } else { (&b - &a + 1u32).to_usize().unwrap_or(usize::MAX) };
        total = total.saturating_mul(width);
        out.push((a, width));
    }
    if total > MAX_ENTRIES {
        return Err(PyValueError::new_err("Comprehension expansion exceeds 100000 positions"));
    }
    Ok(out)
}

/// `itertools.product` of the ranges, the last index fastest (one empty
/// point without ranges).
fn for_each_point(ranges: &[(BigInt, usize)], f: &mut dyn FnMut(&[BigInt]) -> PyResult<()>) -> PyResult<()> {
    if ranges.iter().any(|(_, w)| *w == 0) { return Ok(()); }
    let mut offsets = vec![0usize; ranges.len()];
    loop {
        let point = ranges.iter().zip(&offsets).map(|((a, _), o)| a + BigInt::from(*o)).collect::<Vec<_>>();
        f(&point)?;
        let mut axis = ranges.len();
        loop {
            if axis == 0 { return Ok(()); }
            axis -= 1;
            offsets[axis] += 1;
            if offsets[axis] < ranges[axis].1 { break; }
            offsets[axis] = 0;
        }
    }
}

/// `_expand_array`'s entries, in row-major order. `ArrayComprehension`
/// substitutes the variables one at a time (`temp.subs(var, val)`);
/// `ArrayComprehensionMap` (`mapped`) calls the lambda with no argument or
/// with the product of the indices, and keeps the lambda itself when it takes
/// more arguments (the array constructor then rejects it).
#[pyfunction]
#[pyo3(signature = (function, limits, mapped=false))]
fn _array_comprehension_entries(py: Python<'_>, function: &Bound<'_, PyAny>, limits: &Bound<'_, PyTuple>, mapped: bool) -> PyResult<PyObject> {
    let limits = limits_of(limits)?;
    let ranges = ranges(&limits)?;
    let entries = PyList::empty(py);
    if mapped {
        let argcount: usize = function.getattr("__code__")?.getattr("co_argcount")?.extract()?;
        for_each_point(&ranges, &mut |point| {
            let value = match argcount {
                0 => function.call0()?,
                1 => function.call1((point.iter().fold(BigInt::one(), |p, v| p * v),))?,
                _ => function.clone(),
            };
            entries.append(value)
        })?;
        return Ok(entries.into_any().unbind());
    }
    let native = to_expr(function).ok();
    for_each_point(&ranges, &mut |point| {
        match &native {
            Some(f) => {
                let mut temp = f.clone();
                for ((var, _, _), value) in limits.iter().zip(point) {
                    temp = checked_subs(py, &temp, &[(var.clone(), Expr::num(symox::number::Number::from_bigint(value.clone())))], false)?;
                }
                entries.append(matrix_or_expr(py, temp)?)
            }
            None => {
                let mut temp = function.clone();
                for ((var, _, _), value) in limits.iter().zip(point) {
                    temp = temp.call_method1("subs", (expr_object(py, var.clone())?, expr_object(py, Expr::num(symox::number::Number::from_bigint(value.clone())))?))?;
                }
                entries.append(temp)
            }
        }
    })?;
    Ok(entries.into_any().unbind())
}

/// `subs` with the indices bound: a pair whose old value is an index is not
/// applied to the function, and an index occurring free in a replacement is
/// first renamed to an integer `Dummy` (no capture). Returns the new
/// function (unchanged for a map) and limits.
#[pyfunction]
#[pyo3(signature = (function, limits, args, kwargs=None, mapped=false))]
fn _array_comprehension_subs(py: Python<'_>, function: &Bound<'_, PyAny>, limits: &Bound<'_, PyTuple>, args: &Bound<'_, PyTuple>,
                             kwargs: Option<&Bound<'_, pyo3::types::PyDict>>, mapped: bool) -> PyResult<PyObject> {
    let (args, options) = crate::substitution_arguments(args, kwargs)?;
    let pairs = crate::subs_pairs(&args)?;
    let simultaneous = crate::substitution_options(options.as_ref())?;
    let limits = limits_of(limits)?;
    let replacements = pairs.iter().flat_map(|(_, new)| symox::traverse::free_symbols(new)).collect::<Vec<_>>();
    let mut rename = rustc_hash::FxHashMap::default();
    for (var, _, _) in &limits {
        if replacements.contains(var) {
            let name = var.as_symbol().unwrap_or("i").to_string();
            rename.insert(var.clone(), Expr::dummy(&name, symox::Assumptions::INTEGER));
        }
    }
    let new_function = if mapped { function.clone().unbind() } else {
        let f = to_expr(function)?;
        let protected = pairs.iter().filter(|(old, _)| !limits.iter().any(|(v, _, _)| v == old)).cloned().collect::<Vec<_>>();
        let renamed = symox::subs::xreplace(&f, &rename);
        matrix_or_expr(py, if protected.is_empty() { renamed } else { checked_subs(py, &renamed, &protected, simultaneous)? })?
    };
    let mut new_limits = Vec::with_capacity(limits.len());
    for (var, inf, sup) in &limits {
        let bound = |e: &Expr| if pairs.is_empty() { Ok(e.clone()) } else { checked_subs(py, e, &pairs, simultaneous) };
        let var = rename.get(var).cloned().unwrap_or_else(|| var.clone());
        new_limits.push(expr_object(py, tuple(vec![var, bound(inf)?, bound(sup)?]))?);
    }
    Ok(PyTuple::new(py, [new_function, PyTuple::new(py, new_limits)?.into_any().unbind()])?.into_any().unbind())
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_array_comprehension_limits, m)?)?;
    m.add_function(wrap_pyfunction!(_array_comprehension_shape, m)?)?;
    m.add_function(wrap_pyfunction!(_array_comprehension_len, m)?)?;
    m.add_function(wrap_pyfunction!(_array_comprehension_is_shape_numeric, m)?)?;
    m.add_function(wrap_pyfunction!(_array_comprehension_free_symbols, m)?)?;
    m.add_function(wrap_pyfunction!(_array_comprehension_entries, m)?)?;
    m.add_function(wrap_pyfunction!(_array_comprehension_subs, m)?)?;
    Ok(())
}
