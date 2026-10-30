//! Rank-four matrix derivatives (`symox._matrix_derivatives`): the native
//! matrix lines of `crates/symox-core/src/matrixexpr/lines.rs` and SymPy's
//! `_rank_four` array form, computed in Rust.
//!
//! The binding returns the addends as `(a, b, permutation)` triples
//! (`permutation` is None for a bare `ArrayTensorProduct(a, b)`); the Python
//! wrapper only builds the array-expression objects. `None` when the kernel
//! declines, and the Python algorithm then runs unchanged.
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use crate::matrixexpr::{matrix_object, PyMatrixExpr};

#[pyfunction]
fn _matrix_lines_rank_four(py: Python<'_>, expr: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    use symox::matrixexpr::lines::{matrix_lines_rank_four, RankFourTerm, Value};
    let Ok(x) = x.downcast::<PyMatrixExpr>() else { return Ok(None) };
    let value = if let Ok(m) = expr.downcast::<PyMatrixExpr>() {
        Value::Matrix(m.get().m.clone())
    } else if let Ok(e) = expr.downcast::<crate::PyExprNode>() {
        Value::Scalar(e.get().e.clone())
    } else {
        return Ok(None);
    };
    let Some(terms) = matrix_lines_rank_four(&value, &x.get().m) else { return Ok(None) };
    let out = PyList::empty(py);
    for term in terms {
        let item = match term {
            RankFourTerm::Permuted(a, b, permutation) => {
                let permutation = match permutation { Some(p) => PyList::new(py, p)?.into_any(), None => py.None().into_bound(py) };
                PyTuple::new(py, [matrix_object(py, &a)?, matrix_object(py, &b)?, permutation])?
            },
            RankFourTerm::Product(a, b) => PyTuple::new(py, [matrix_object(py, &a)?, matrix_object(py, &b)?, py.None().into_bound(py)])?,
        };
        out.append(item)?;
    }
    Ok(Some(out.into_any().unbind()))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_matrix_lines_rank_four, m)?)?;
    Ok(())
}
