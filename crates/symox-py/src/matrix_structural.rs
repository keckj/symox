//! Shared native kernels and thin dense bindings for structural transforms.
use super::{err, PyMatrix};
use num_bigint::BigInt;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use symox::{
    matrix::{
        structural::{self, IndexPermutation},
        Matrix,
    },
    Expr,
};

pub fn plan(
    py: Python<'_>,
    specification: Vec<Vec<BigInt>>,
    dimension: BigInt,
    cycles: bool,
    backward: bool,
) -> PyResult<Bound<'_, PyDict>> {
    let plan = IndexPermutation::new(&dimension, &specification, cycles, backward)
        .map_err(|e| err(py, e))?;
    let output = PyDict::new(py);
    for (a, b) in plan.destinations() {
        output.set_item(a, b)?;
    }
    Ok(output)
}
pub fn limit(py: Python<'_>, count: BigInt) -> PyResult<()> {
    structural::check_count(&count).map_err(|e| err(py, e))
}
pub fn values<'py>(
    py: Python<'py>,
    add: bool,
    coefficient: &Bound<'py, PyAny>,
    targets: &Bound<'py, PyAny>,
    sources: &Bound<'py, PyAny>,
) -> PyResult<Bound<'py, PyList>> {
    let k = crate::scalar_expr(coefficient)?;
    let extract = |v: &Bound<'_, PyAny>| -> PyResult<Vec<Expr>> {
        // Input is a concrete compact vector assembled by the sparse adapter;
        // consume boundedly even if this private entry point is called directly.
        let mut output = Vec::new();
        for item in v.try_iter()? {
            if output.len() >= structural::MAX_ITEMS {
                return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                    "Matrix structural transformation resource limit exceeded",
                ));
            }
            output.push(crate::scalar_expr(&item?)?);
        }
        Ok(output)
    };
    let a = extract(targets)?;
    let b = extract(sources)?;
    let output = structural::elementary_values(add, &k, &a, &b).map_err(|e| err(py, e))?;
    PyList::new(py, output.into_iter().map(crate::wrap))
}

pub fn permute(
    py: Python<'_>,
    matrix: &Matrix,
    perm: &Bound<'_, PyAny>,
    orientation: &str,
    direction: &str,
) -> PyResult<Matrix> {
    let rows = match orientation {
        "rows" => true,
        "cols" | "columns" => false,
        _ => {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "orientation must be 'rows' or 'cols'",
            ))
        }
    };
    let dimension = if rows { matrix.rows() } else { matrix.cols() };
    let normalized = py
        .import("symox._matrix_structural")?
        .getattr("normalize_permutation")?
        .call1((perm, dimension, orientation, direction))?;
    let (specification, cycles, backward, _): (Vec<Vec<BigInt>>, bool, bool, bool) =
        normalized.extract()?;
    IndexPermutation::new(&BigInt::from(dimension), &specification, cycles, backward)
        .and_then(|p| p.apply(matrix, rows))
        .map_err(|e| err(py, e))
}
pub fn elementary(
    py: Python<'_>,
    matrix: &Matrix,
    rows: bool,
    op: &str,
    target: Option<&Bound<'_, PyAny>>,
    k: Option<&Bound<'_, PyAny>>,
    first: Option<&Bound<'_, PyAny>>,
    second: Option<&Bound<'_, PyAny>>,
) -> PyResult<Matrix> {
    let size = if rows { matrix.rows() } else { matrix.cols() };
    let normalized = py
        .import("symox._matrix_structural")?
        .getattr("normalize_elementary")?
        .call1((op, size, target, k, first, second))?;
    let (kind, target, source, k): (u8, usize, usize, Bound<'_, PyAny>) = normalized.extract()?;
    let coefficient = crate::scalar_expr(&k)?;
    matrix
        .elementary_transform(rows, kind, target, source, &coefficient)
        .map_err(|e| err(py, e))
}

/// Release the Matrix borrow around callbacks: they may read or mutate that
/// very matrix, and SymPy applies updates sequentially rather than atomically.
pub fn callback(
    py: Python<'_>,
    slf: Py<PyMatrix>,
    index: isize,
    source: Option<isize>,
    function: &Bound<'_, PyAny>,
    rows: bool,
) -> PyResult<()> {
    let (target, source, width) = {
        let matrix = slf.borrow(py);
        matrix.require_mutable()?;
        let dimension = if rows {
            matrix.m.rows()
        } else {
            matrix.m.cols()
        };
        let width = if rows {
            matrix.m.cols()
        } else {
            matrix.m.rows()
        };
        // An empty line invokes no callback and accesses no index, as in SymPy.
        if width == 0 {
            return Ok(());
        }
        let target = super::norm_index(index, dimension)?;
        let source = source
            .map(|i| super::norm_index(i, dimension))
            .transpose()?;
        structural::check_count(&BigInt::from(width)).map_err(|e| err(py, e))?;
        (target, source, width)
    };
    for p in 0..width {
        let (a, b) = {
            let matrix = slf.borrow(py);
            let dimension = if rows {
                matrix.m.rows()
            } else {
                matrix.m.cols()
            };
            let current_width = if rows {
                matrix.m.cols()
            } else {
                matrix.m.rows()
            };
            if target >= dimension || source.is_some_and(|i| i >= dimension) || p >= current_width {
                return Err(pyo3::exceptions::PyIndexError::new_err(
                    "Matrix shape changed during callback",
                ));
            }
            let get = |i| {
                if rows {
                    matrix.m.get(i, p)
                } else {
                    matrix.m.get(p, i)
                }
                .clone()
            };
            (get(target), source.map(get))
        };
        let result = if let Some(b) = b {
            function.call1((crate::wrap(a), crate::wrap(b)))?
        } else {
            function.call1((crate::wrap(a), p))?
        };
        let value = crate::scalar_expr(&result)?;
        let mut matrix = slf.borrow_mut(py);
        // A callback may change shape. Validate again instead of indexing with
        // stale dimensions; this also preserves earlier successful updates.
        if (rows && (target >= matrix.m.rows() || p >= matrix.m.cols()))
            || (!rows && (p >= matrix.m.rows() || target >= matrix.m.cols()))
        {
            return Err(pyo3::exceptions::PyIndexError::new_err(
                "Matrix shape changed during callback",
            ));
        }
        if rows {
            matrix.m.set(target, p, value);
        } else {
            matrix.m.set(p, target, value);
        }
    }
    Ok(())
}
