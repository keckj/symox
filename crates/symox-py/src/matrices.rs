//! Python bindings for symbolic matrices (`symox.Matrix`).
#[path = "matrix_structural.rs"]
mod structural;
#[path = "matrix_auxiliary.rs"]
mod auxiliary;

/// Install the host `simplify` used by the native pivot search.
pub fn register_pivot_host() {
    auxiliary::register_pivot_host();
}

use pyo3::exceptions::{PyIndexError, PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyInt, PyList, PySlice, PyTuple};
use rustc_hash::FxHashMap;
use num_traits::ToPrimitive;

use symox::expr::Expr;
use symox::matrix::{self, DetMethod, Matrix, MatrixError};

use crate::{diff_vars, subs_pairs, to_expr, wrap, PyExpr};

pub(crate) fn err(py: Python<'_>, e: MatrixError) -> PyErr {
    let (cls, msg) = match &e {
        MatrixError::Shape(m) => ("ShapeError", m.clone()),
        MatrixError::NonSquare(m) => ("NonSquareMatrixError", m.clone()),
        MatrixError::NonInvertible(m) => ("NonInvertibleMatrixError", m.clone()),
        MatrixError::Other(m) => ("MatrixError", m.clone()),
        MatrixError::NonPositiveDefinite(m) => ("NonPositiveDefiniteMatrixError", m.clone()),
        MatrixError::NotImplemented(m) => return PyNotImplementedError::new_err(m.clone()),
        MatrixError::Value(m) => return PyValueError::new_err(m.clone()),
        MatrixError::Index(m) => return PyIndexError::new_err(m.clone()),
        MatrixError::Type(m) => return PyTypeError::new_err(m.clone()),
    };
    match py.import("symox").and_then(|m| m.getattr(cls)) {
        Ok(c) => match c.call1((msg.clone(),)) {
            Ok(inst) => PyErr::from_value(inst),
            Err(_) => PyValueError::new_err(msg),
        },
        Err(_) => PyValueError::new_err(msg),
    }
}


fn modular_modulus(py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<num_bigint::BigInt> {
    if value.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err("inv_mod: modulus m must be an integer"));
    }
    py.import("operator")?.getattr("index")?.call1((value,))
        .and_then(|integer| integer.extract::<num_bigint::BigInt>())
        .map_err(|_| PyTypeError::new_err("inv_mod: modulus m must be an integer"))
}

struct LllDelta(Expr);
impl Default for LllDelta {
    fn default() -> Self { Self(Expr::rational(3,4)) }
}
impl<'py> FromPyObject<'py> for LllDelta {
    fn extract_bound(value: &Bound<'py, PyAny>) -> PyResult<Self> {
        crate::scalar_expr(value).map(Self)
    }
}
fn lll_error(py: Python<'_>, error: matrix::lll::LllError) -> PyErr {
    use matrix::lll::LllError as E;
    let message=error.to_string();
    let name=match error {
        E::Delta=>"DMValueError", E::Shape=>"DMShapeError", E::Rank=>"DMRankError",
        E::Parameter|E::Domain=>return crate::poly::domains::CoercionFailed::new_err(message),
        E::Resource|E::Certificate=>return PyNotImplementedError::new_err(message),
    };
    match py.import("symox.polys.matrices.exceptions").and_then(|m|m.getattr(name)).and_then(|c|c.call1((&message,))) {
        Ok(value)=>PyErr::from_value(value), Err(error)=>error,
    }
}

// Keep omitted/default cancel distinct from explicitly disabled None/False.
// A callable is an explicit matrix postprocessor, never a symbolic fallback.
enum LogSimplifier {
    Native,
    Disabled,
    Callback(Py<PyAny>),
}

enum PivotSimplifier {
    Native,
    Callback(Py<PyAny>),
}

impl<'py> FromPyObject<'py> for PivotSimplifier {
    fn extract_bound(value:&Bound<'py,PyAny>)->PyResult<Self> {
        if value.is_none() || auxiliary::boolean_scalar(value)?.is_some() {Ok(Self::Native)}
        else {
            auxiliary::checked_callback(Some(value),"matrix simplify")?;
            Ok(Self::Callback(value.clone().unbind()))
        }
    }
}

impl PivotSimplifier {
    fn callback<'a,'py>(&'a self,py:Python<'py>)->Option<&'a Bound<'py,PyAny>> {
        match self {Self::Native=>None,Self::Callback(callback)=>Some(callback.bind(py))}
    }
}

impl<'py> FromPyObject<'py> for LogSimplifier {
    fn extract_bound(value: &Bound<'py, PyAny>) -> PyResult<Self> {
        if !value.is_truthy()? {
            Ok(Self::Disabled)
        } else if value.is_callable() {
            if value.is(&value.py().import("symox._core")?.getattr("cancel")?) {
                Ok(Self::Native)
            } else {
                Ok(Self::Callback(value.clone().unbind()))
            }
        } else {
            Err(PyTypeError::new_err("Matrix.log: simplify must be a callable, False or None"))
        }
    }
}

fn reject_callback(name: &str, callback: Option<&Bound<'_, PyAny>>) -> PyResult<()> {
    if callback.is_some_and(|c| !c.is_none()) {
        Err(PyNotImplementedError::new_err(format!("custom {name} callbacks are not implemented for matrices")))
    } else {
        Ok(())
    }
}

fn reject_kwargs(kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    if let Some(kwargs) = kwargs {
        if let Some((key, _)) = kwargs.iter().next() {
            return Err(PyNotImplementedError::new_err(format!("matrix option {key} is not implemented")));
        }
    }
    Ok(())
}

/// Applies SymPy's `cls=` option of the special constructors: the native
/// classes keep the native result (immutable when asked), any other class
/// receives the explicit result through its own constructor.
fn finish_special(py: Python<'_>, mut value: PyMatrix, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let mut cls = None;
    if let Some(kwargs) = kwargs {
        for (key, item) in kwargs.iter() {
            if key.extract::<String>().ok().as_deref() == Some("cls") {
                if !item.is_none() { cls = Some(item); }
            } else {
                return Err(PyNotImplementedError::new_err(format!("matrix option {key} is not implemented")));
            }
        }
    }
    let Some(cls) = cls else { return Ok(value.into_pyobject(py)?.into_any().unbind()) };
    let base = py.get_type::<PyMatrix>();
    let native = cls.is(&base)
        || MUTABLE_TYPE.get(py).is_some_and(|t| cls.is(t.bind(py)))
        || IMMUTABLE_TYPE.get(py).is_some_and(|t| cls.is(t.bind(py)));
    if native {
        if IMMUTABLE_TYPE.get(py).is_some_and(|t| cls.is(t.bind(py))) { value.immutable = true; }
        return Ok(value.into_pyobject(py)?.into_any().unbind());
    }
    Ok(cls.call1((value,))?.unbind())
}

trait IntoPy<T> {
    fn py(self, py: Python<'_>) -> PyResult<T>;
}

impl<T> IntoPy<T> for Result<T, MatrixError> {
    fn py(self, py: Python<'_>) -> PyResult<T> {
        self.map_err(|e| err(py, e))
    }
}

/// A mutable dense symbolic matrix.
/// Instances carry a `__dict__` like SymPy matrices (downstream caches such as
/// SageMath's `_sage_object` are stored on the matrix).
#[pyclass(name = "DenseMatrix", module = "symox.matrices.dense", weakref, subclass, dict, extends = pyo3::PyAny)]
#[derive(Clone)]
pub struct PyMatrix {
    pub m: Matrix,
    immutable: bool,
}

/// The Python classes `MutableDenseMatrix` and `ImmutableDenseMatrix`
/// (subclasses of the native `DenseMatrix` that add no storage). Every
/// native matrix result is allocated as one of them, according to its
/// mutability, as SymPy's `cls._new` does.
static MUTABLE_TYPE: pyo3::sync::GILOnceCell<Py<pyo3::types::PyType>> = pyo3::sync::GILOnceCell::new();
static IMMUTABLE_TYPE: pyo3::sync::GILOnceCell<Py<pyo3::types::PyType>> = pyo3::sync::GILOnceCell::new();

#[pyfunction]
fn _register_matrix_types(py: Python<'_>, mutable: Bound<'_, pyo3::types::PyType>, immutable: Bound<'_, pyo3::types::PyType>) -> PyResult<()> {
    let base = py.get_type::<PyMatrix>();
    if !mutable.is_subclass(&base)? || !immutable.is_subclass(&base)? {
        return Err(PyTypeError::new_err("matrix classes must derive from DenseMatrix"));
    }
    let _ = MUTABLE_TYPE.set(py, mutable.unbind());
    let _ = IMMUTABLE_TYPE.set(py, immutable.unbind());
    Ok(())
}

impl<'py> IntoPyObject<'py> for PyMatrix {
    type Target = PyMatrix;
    type Output = Bound<'py, PyMatrix>;
    type Error = PyErr;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Bound<'py, PyMatrix>> {
        let cell = if self.immutable { &IMMUTABLE_TYPE } else { &MUTABLE_TYPE };
        match cell.get(py) {
            Some(class) => {
                use pyo3::impl_::pyclass_init::PyObjectInit;
                let initializer = pyo3::PyClassInitializer::from(self);
                // SAFETY: the registered class is a Python subclass of
                // DenseMatrix with the same layout (checked at registration).
                unsafe {
                    let obj = initializer.into_new_object(py, class.bind(py).as_type_ptr())?;
                    Ok(Bound::from_owned_ptr(py, obj).downcast_into_unchecked())
                }
            }
            None => Bound::new(py, self),
        }
    }
}

pub(crate) fn wrapm(m: Matrix) -> PyMatrix {
    PyMatrix { m, immutable: false }
}

pub(crate) fn wrapm_immutable(m: Matrix) -> PyMatrix {
    PyMatrix { m, immutable: true }
}

pub(crate) fn matrix_value_hash(rows: usize, cols: usize, mut entries: Vec<((usize, usize), Expr)>) -> isize {
    use std::hash::{Hash, Hasher};
    entries.sort_unstable_by_key(|(index, _)| *index);
    let mut hash = std::collections::hash_map::DefaultHasher::new();
    (rows, cols).hash(&mut hash);
    entries.hash(&mut hash);
    let value = hash.finish() as isize;
    if value == -1 { -2 } else { value }
}

/// Shared value hash for dense and sparse immutable matrices, without
/// constructing implicit sparse zeros or a dense temporary.
#[pyfunction]
fn _sparse_matrix_hash(rows: usize, cols: usize, entries: &Bound<'_, PyDict>) -> PyResult<isize> {
    let mut values = Vec::with_capacity(entries.len());
    for (key, value) in entries.iter() {
        let (i, j): (usize, usize) = key.extract()?;
        if i >= rows || j >= cols {
            return Err(PyValueError::new_err("Sparse coordinate is outside matrix shape"));
        }
        let value = to_expr(&value)?;
        if !value.is_zero() { values.push(((i, j), value)); }
    }
    Ok(matrix_value_hash(rows, cols, values))
}

impl PyMatrix {
    fn wrapped(&self, m: Matrix) -> PyMatrix {
        PyMatrix { m, immutable: self.immutable }
    }

    fn with_peer(&self, peer: &Bound<'_, PyAny>, m: Matrix) -> PyMatrix {
        let immutable = self.immutable || peer.downcast::<PyMatrix>()
            .is_ok_and(|p| p.borrow().immutable)
            || (peer.getattr("_symox_sparse").is_ok_and(|value| value.extract::<bool>().unwrap_or(false))
                && peer.getattr("_immutable").is_ok_and(|value| value.extract::<bool>().unwrap_or(false)));
        PyMatrix { m, immutable }
    }

    fn require_mutable(&self) -> PyResult<()> {
        if self.immutable {
            // SymPy's MatrixBase wording names the immutable class.
            let class = Python::with_gil(|py| py.import("symox")
                .and_then(|module| module.getattr("ImmutableDenseMatrix"))
                .and_then(|class| class.repr().map(|text| text.to_string())))
                .unwrap_or_else(|_| "an immutable matrix".to_owned());
            Err(PyTypeError::new_err(format!("Cannot set values of {class}")))
        } else {
            Ok(())
        }
    }
}

fn exprs_list<'py>(py: Python<'py>, v: &[Expr]) -> PyResult<Bound<'py, PyList>> {
    PyList::new(py, v.iter().map(|e| wrap(e.clone())))
}

fn is_sympy_obj(obj: &Bound<'_, PyAny>) -> bool {
    crate::defining_module(obj.get_type().as_any()).starts_with("sympy")
}

/// Convert a Python object to a matrix (symox/SymPy matrices, nested
/// sequences, numpy arrays).
pub fn to_matrix(obj: &Bound<'_, PyAny>) -> PyResult<Matrix> {
    if obj.is_instance_of::<crate::matrixexpr::PyBlockGrid>() {
        return Err(PyTypeError::new_err("A block grid contains matrix-valued cells and cannot be converted to a scalar matrix"));
    }
    if let Ok(matrix)=obj.downcast::<crate::matrixexpr::PyMatrixExpr>() {
        return matrix.get().m.as_explicit().map_err(|e|err(obj.py(),e));
    }
    if let Ok(m) = obj.downcast::<PyMatrix>() {
        return Ok(m.borrow().m.clone());
    }
    if is_sympy_obj(obj) && obj.getattr("is_MatrixExpr").is_ok_and(|v|v.is_truthy().unwrap_or(false)) {
        return crate::matrixexpr::to_matrixexpr(obj)?.as_explicit().map_err(|e|err(obj.py(),e));
    }
    if obj.hasattr("tolist")? && (is_sympy_obj(obj) || obj.hasattr("shape")?) {
        let shape: Option<(usize, usize)> = obj.getattr("shape").ok().and_then(|s| s.extract().ok());
        let l = obj.call_method0("tolist")?;
        let m = from_nested(&l)?;
        if let Some((r, c)) = shape {
            if r * c == 0 {
                return Ok(Matrix::zeros(r, c));
            }
        }
        return Ok(m);
    }
    from_nested(obj)
}

/// A scalar matrix index (an integer or a native scalar expression), as
/// SymPy's `isinstance(i, Expr)` test sees it.
fn scalar_index(k: &Bound<'_, PyAny>) -> Option<Expr> {
    if k.downcast::<PySlice>().is_ok() || is_seq(k) || k.is_instance_of::<PyBool>() {
        return None;
    }
    to_expr(k).ok().filter(|e| !matches!(e.kind(), symox::expr::Kind::Bool(_) | symox::expr::Kind::Rel(..) | symox::expr::Kind::Func(symox::Func::MatrixProxy(_), _)))
}

/// A Python `int` (not a `bool`): the only index type SymPy's `(j < 0) is
/// True` boundary test can see.
fn host_int(k: &Bound<'_, PyAny>) -> bool {
    k.is_instance_of::<PyInt>() && !k.is_instance_of::<PyBool>()
}

fn is_seq(obj: &Bound<'_, PyAny>) -> bool {
    obj.downcast::<PyList>().is_ok() || obj.downcast::<PyTuple>().is_ok() || is_sequence_row(obj)
}

/// Host sequences SymPy's `is_sequence` accepts as matrix rows, such as a
/// geometry `Point` (its class sets `_symox_matrix_row`).
fn is_sequence_row(obj: &Bound<'_, PyAny>) -> bool {
    !obj.is_instance_of::<crate::PyExprNode>()
        && obj.downcast::<PyMatrix>().is_err()
        && obj.get_type().getattr("_symox_matrix_row").is_ok_and(|v| v.is_truthy().unwrap_or(false))
}

fn from_nested(obj: &Bound<'_, PyAny>) -> PyResult<Matrix> {
    let items: Vec<Bound<'_, PyAny>> = obj.try_iter()?.collect::<PyResult<_>>()?;
    if items.is_empty() {
        return Ok(Matrix::zeros(0, 0));
    }
    let is_matrix = |item: &Bound<'_, PyAny>| -> bool {
        item.downcast::<PyMatrix>().is_ok()
            || item.is_instance_of::<crate::matrixexpr::PyMatrixExpr>()
            || (is_sympy_obj(item) && (item.getattr("is_MatrixExpr").is_ok_and(|v|v.is_truthy().unwrap_or(false))
                || (item.hasattr("shape").unwrap_or(false) && item.hasattr("tolist").unwrap_or(false))))
    };
    // A flat list of matrices is stacked vertically. Nested rows of
    // matrices are block rows and are joined horizontally first.
    if items.iter().any(&is_matrix) {
        let blocks: Vec<Matrix> = items.iter().map(|item| {
            if is_matrix(item) { to_matrix(item) }
            else { Ok(Matrix::column(vec![scalar_matrix_entry(item)?])) }
        }).collect::<PyResult<_>>()?;
        // As in the reference, a 0x0 block adds no rows and sets no width:
        // Matrix([q, Matrix([])]) is q.
        let refs: Vec<&Matrix> = blocks.iter().filter(|m| m.rows() != 0 || m.cols() != 0).collect();
        return Matrix::vstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()));
    }
    if items.iter().all(|i| is_seq(i) || (i.hasattr("tolist").unwrap_or(false) && !is_sympy_obj(i) && i.hasattr("shape").unwrap_or(false))) {
        let mut rows = Vec::with_capacity(items.len());
        let mut plain = true;
        for it in &items {
            let it = if is_seq(it) { it.clone() } else { it.call_method0("tolist")? };
            let entries: Vec<Bound<'_, PyAny>> = it.try_iter()?.collect::<PyResult<_>>()?;
            if entries.iter().any(&is_matrix) {
                plain = false;
                let blocks: Vec<Matrix> = entries.iter().map(|entry| {
                    if is_matrix(entry) { to_matrix(entry) }
                    else { Ok(Matrix::column(vec![scalar_matrix_entry(entry)?])) }
                }).collect::<PyResult<_>>()?;
                let refs: Vec<&Matrix> = blocks.iter().collect();
                rows.push(Matrix::hstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()))?);
            } else {
                let row: Vec<Expr> = entries.iter().map(scalar_matrix_entry).collect::<PyResult<_>>()?;
                rows.push(Matrix::new(1, row.len(), row).map_err(|e| PyValueError::new_err(e.to_string()))?);
            }
        }
        // The reference rejects plain rows of unequal lengths.
        if plain && rows.windows(2).any(|w| w[0].cols() != w[1].cols()) {
            return Err(PyValueError::new_err("mismatched dimensions"));
        }
        let refs: Vec<&Matrix> = rows.iter().collect();
        return Matrix::vstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()));
    }
    // A sequence mixed with scalars is not a list of lists.
    if items.iter().any(is_seq) {
        return Err(PyValueError::new_err("expecting list of lists"));
    }
    let col: Vec<Expr> = items.iter().map(scalar_matrix_entry).collect::<PyResult<_>>()?;
    Ok(Matrix::column(col))
}

impl PyMatrix {
    /// SymPy's `minor`/`cofactor` checks: squareness first, then the
    /// `minor_submatrix` bounds (negative indices wrap once).
    fn minor_index(&self, py: Python<'_>, i: isize, j: isize, nonempty: bool) -> PyResult<(usize, usize)> {
        let (rows, cols) = (self.m.rows() as isize, self.m.cols() as isize);
        if rows != cols || (nonempty && rows < 1) {
            return Err(err(py, MatrixError::NonSquare(String::new())));
        }
        let (i, j) = (if i < 0 { i + rows } else { i }, if j < 0 { j + cols } else { j });
        if !(0..rows).contains(&i) || !(0..cols).contains(&j) {
            return Err(PyValueError::new_err(format!(
                "`i` and `j` must satisfy 0 <= i < ``M.rows`` ({rows})and 0 <= j < ``M.cols`` ({cols})."
            )));
        }
        Ok((i as usize, j as usize))
    }
}

fn norm_index(i: isize, n: usize) -> PyResult<usize> {
    let j = if i < 0 { i + n as isize } else { i };
    if j < 0 || j >= n as isize {
        return Err(PyIndexError::new_err(format!("Index out of range: a[{i}]")));
    }
    Ok(j as usize)
}

/// Indices selected by an int, a slice or a sequence of ints.
fn key_indices(k: &Bound<'_, PyAny>, n: usize) -> PyResult<(Vec<usize>, bool)> {
    if let Ok(s) = k.downcast::<PySlice>() {
        let ind = s.indices(n as isize)?;
        let mut v = Vec::new();
        let mut i = ind.start;
        while (ind.step > 0 && i < ind.stop) || (ind.step < 0 && i > ind.stop) {
            v.push(i as usize);
            i += ind.step;
        }
        return Ok((v, true));
    }
    if k.is_instance_of::<PyInt>() || k.hasattr("__index__")? && !is_seq(k) {
        let i: isize = k.extract()?;
        return Ok((vec![norm_index(i, n)?], false));
    }
    if is_seq(k) {
        let items: Vec<Bound<'_, PyAny>> = k.try_iter()?.collect::<PyResult<_>>()?;
        // SymPy's ``extract``: a list made only of bools is a mask.
        if !items.is_empty() && items.iter().all(|x| x.is_instance_of::<pyo3::types::PyBool>()) {
            let mut v = Vec::new();
            for (i, x) in items.iter().enumerate() {
                if x.is_truthy()? {
                    v.push(norm_index(i as isize, n)?);
                }
            }
            return Ok((v, true));
        }
        let mut v = Vec::new();
        for x in items {
            let i: isize = x.extract()?;
            v.push(norm_index(i, n)?);
        }
        return Ok((v, true));
    }
    Err(PyTypeError::new_err("invalid index"))
}

/// An N-dim array operand: SymPy's NDimArray has the higher `_op_priority`,
/// so `call_highest_priority` hands the operation to its reflected method,
/// which refuses a matrix (TypeError for `+`/`-`).
fn ndim_array(o: &Bound<'_, PyAny>) -> bool {
    !o.is_instance_of::<PyMatrix>() && o.getattr(pyo3::intern!(o.py(), "is_NDimArray")).is_ok_and(|v| v.is_truthy().unwrap_or(false))
}

fn scalar_or_matrix(o: &Bound<'_, PyAny>) -> PyResult<Result<Matrix, Expr>> {
    if let Ok(m) = o.downcast::<PyMatrix>() {
        return Ok(Ok(m.borrow().m.clone()));
    }
    if is_seq(o) {
        return Ok(Ok(from_nested(o)?));
    }
    if o.hasattr("shape")? && o.hasattr("tolist")? {
        // Two-dimensional array operands denote matrices. Other ranks retain
        // NumPy's scalar/broadcast dispatch instead of becoming column vectors.
        let shape = o.getattr("shape")?;
        if is_sympy_obj(o) || shape.extract::<(usize, usize)>().is_ok() {
            return Ok(Ok(to_matrix(o)?));
        }
    }
    Ok(Err(scalar_matrix_entry(o)?))
}

fn scalar_matrix_entry(value: &Bound<'_, PyAny>) -> PyResult<Expr> {
    if value.is_instance_of::<crate::matrixexpr::PyMatrixExpr>() || value.is_instance_of::<crate::matrixexpr::PyBlockGrid>()
        || value.getattr("is_Matrix").is_ok_and(|v|v.is_truthy().unwrap_or(false)) {
        return Err(PyTypeError::new_err("Dense matrices require scalar entries; materialize the matrix explicitly or use a typed block grid"));
    }
    // RepMatrix._dod_to_DomainMatrix: a Poly is not an Expr. SymPy keeps it
    // and warns; the native matrix stores its expression.
    if value.is_instance_of::<crate::poly::PyPoly>() { warn_non_expr_entry(value.py())?; }
    to_expr(value)
}

fn warn_non_expr_entry(py: Python<'_>) -> PyResult<()> {
    let kwargs = PyDict::new(py);
    kwargs.set_item("deprecated_since_version", "1.9")?;
    kwargs.set_item("active_deprecations_target", "deprecated-non-expr-in-matrix")?;
    kwargs.set_item("stacklevel", 2)?;
    py.import("symox.utilities.exceptions")?.getattr("sympy_deprecation_warning")?.call((
        "\n                non-Expr objects in a Matrix is deprecated. Matrix represents\n                a mathematical matrix. To represent a container of non-numeric\n                entities, Use a list of lists, TableForm, NumPy array, or some\n                other data structure instead.\n                ",), Some(&kwargs))?;
    Ok(())
}

fn has_matrix_entry(obj: &Bound<'_, PyAny>, depth: usize) -> bool {
    let is_matrix = |v: &Bound<'_, PyAny>| v.getattr("is_Matrix").is_ok_and(|m| m.is_truthy().unwrap_or(false))
        || v.getattr("is_MatrixExpr").is_ok_and(|m| m.is_truthy().unwrap_or(false));
    if !is_seq(obj) { return false; }
    obj.try_iter().map(|items| items.flatten().any(|item| is_matrix(&item) || (depth > 0 && has_matrix_entry(&item, depth - 1))))
        .unwrap_or(false)
}

/// SymPy forwards constructor keywords to ``_handle_creation_inputs``: only
/// ``evaluate`` and ``copy`` have meaning there and other keywords are ignored.
pub(crate) fn matrix_creation_options(args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
    let Some(kwargs) = kwargs else { return Ok(()) };
    if let Some(copy) = kwargs.get_item("copy")? {
        if !copy.is_truthy()? && args.len() != 3 {
            return Err(PyTypeError::new_err("'copy=False' requires a matrix be initialized as rows,cols,[list]"));
        }
    }
    if let Some(evaluate) = kwargs.get_item("evaluate")? {
        if !evaluate.is_truthy()? && args.len() == 1 && has_matrix_entry(&args.get_item(0)?, 1) {
            return Err(PyNotImplementedError::new_err("Matrix(..., evaluate=False) with matrix entries is not implemented"));
        }
    }
    Ok(())
}

impl PyMatrix {
    fn py_new_value(py: Python<'_>, args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        matrix_creation_options(args, kwargs)?;
        match args.len() {
            0 => Ok(wrapm(Matrix::zeros(0, 0))),
            1 => Ok(wrapm(to_matrix(&args.get_item(0)?)?)),
            3 => {
                let (r, c) = creation_dims(&args.get_item(0)?, Some(&args.get_item(1)?))?;
                let src = args.get_item(2)?;
                if src.is_callable() && !is_seq(&src) {
                    let mut data = Vec::with_capacity(r * c);
                    for i in 0..r {
                        for j in 0..c {
                            data.push(scalar_matrix_entry(&src.call1((i, j))?)?);
                        }
                    }
                    return Ok(wrapm(Matrix::new(r, c, data).py(py)?));
                }
                let mut flat: Vec<Expr> = Vec::with_capacity(r * c);
                for it in src.try_iter()? {
                    let it = it?;
                    if is_seq(&it) {
                        for x in it.try_iter()? {
                            flat.push(scalar_matrix_entry(&x?)?);
                        }
                    } else {
                        flat.push(scalar_matrix_entry(&it)?);
                    }
                }
                Ok(wrapm(Matrix::new(r, c, flat).py(py)?))
            }
            // SymPy has no two-argument form: Matrix(2, 2) is a TypeError.
            2 => Err(PyTypeError::new_err(
                "\nData type not understood; expecting list of lists or lists of values.")),
            _ => Err(PyTypeError::new_err("Matrix() takes at most 3 arguments")),
        }
    }
}

#[pymethods]
impl PyMatrix {
    #[new]
    #[classmethod]
    #[pyo3(signature = (*args, **kwargs))]
    fn py_new(cls: &Bound<'_, pyo3::types::PyType>, py: Python<'_>, args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<pyo3::PyClassInitializer<Self>> {
        let mut value = Self::py_new_value(py, args, kwargs)?;
        if let Some(immutable) = IMMUTABLE_TYPE.get(py) {
            if cls.is_subclass(immutable.bind(py))? { value.immutable = true; }
        }
        // The native base class itself builds the public mutable class.
        if cls.is(&py.get_type::<PyMatrix>()) {
            return Ok(pyo3::PyClassInitializer::from(value.into_pyobject(py)?.unbind()));
        }
        Ok(pyo3::PyClassInitializer::from(value))
    }

    // ------------------------------------------------------------------
    // printing / conversion
    // ------------------------------------------------------------------

    fn __str__(&self) -> String {
        self.m.to_str()
    }

    fn __repr__(&self) -> String {
        self.m.to_repr()
    }

    fn _repr_latex_(&self) -> String {
        format!("$\\displaystyle {}$", self.m.to_str())
    }

    fn tolist<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        let rows: Vec<Bound<'py, PyList>> =
            self.m.to_rows().iter().map(|r| exprs_list(py, r)).collect::<PyResult<_>>()?;
        PyList::new(py, rows)
    }

    /// Convert to a SymPy ``Matrix``.
    fn to_sympy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let (r, c) = self.m.shape();
        let name = if self.immutable { "ImmutableMatrix" } else { "Matrix" };
        let entries = self.m.data().iter().map(|e| crate::expr_to_sympy(py, e)).collect::<PyResult<Vec<_>>>()?;
        py.import("sympy")?.getattr(name)?.call1((r, c, PyList::new(py, entries)?))
    }

    fn _sympy_<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.to_sympy(py)
    }

    // ------------------------------------------------------------------
    // shape
    // ------------------------------------------------------------------

    #[getter]
    fn shape(&self) -> (usize, usize) {
        self.m.shape()
    }

    #[getter]
    fn rows(&self) -> usize {
        self.m.rows()
    }

    #[getter]
    fn cols(&self) -> usize {
        self.m.cols()
    }

    #[getter]
    fn is_square(&self) -> bool {
        self.m.is_square()
    }

    #[getter]
    fn is_zero_matrix(&self) -> Option<bool> {
        let mut unknown = false;
        for entry in self.m.data() {
            match symox::assumptions::is_zero(entry)
                .or_else(|| matrix::ratfunc::is_zero_assume(entry)) {
                Some(false) => return Some(false),
                None => unknown = true,
                Some(true) => {},
            }
        }
        if unknown { None } else { Some(true) }
    }

    #[getter]
    fn free_symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PySet>> {
        pyo3::types::PySet::new(py, self.m.free_symbols().into_iter().map(wrap))
    }

    fn is_symmetric(&self) -> bool {
        self.m.is_symmetric()
    }

    fn is_diagonal(&self) -> bool {
        self.m.is_diagonal()
    }

    #[getter]
    fn is_upper(&self) -> bool {
        self.m.is_upper()
    }

    #[getter]
    fn is_lower(&self) -> bool {
        self.m.is_lower()
    }

    fn __len__(&self) -> usize {
        self.m.len()
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        exprs_list(py, self.m.data())?.try_iter().map(|i| i.into_any())
    }

    fn copy(&self) -> PyMatrix {
        self.clone()
    }

    #[pyo3(signature = (*patterns))]
    fn has(&self, patterns: &Bound<'_, PyTuple>) -> PyResult<bool> {
        for p in patterns.iter() {
            let p = to_expr(&p)?;
            if self.m.data().iter().any(|e| symox::traverse::has(e, &p)) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    // ------------------------------------------------------------------
    // indexing
    // ------------------------------------------------------------------

    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        let (r, c) = self.m.shape();
        if let Ok(t) = key.downcast::<PyTuple>() {
            if t.len() != 2 {
                return Err(PyIndexError::new_err("expected (row, col)"));
            }
            let (a, b) = (t.get_item(0)?, t.get_item(1)?);
            let (ri, rs, ci, cs) = match (key_indices(&a, r), key_indices(&b, c)) {
                (Ok((ri, rs)), Ok((ci, cs))) => (ri, rs, ci, cs),
                (first, second) => {
                    let error = first.err().or(second.err()).expect("one index failed");
                    // A symbolic index gives `MatrixElement(M, i, j)` (SymPy's
                    // RepMatrix.__getitem__); any other failure is re-raised.
                    return match (scalar_index(&a), scalar_index(&b)) {
                        (Some(i), Some(j)) if !i.is_number() || !j.is_number() => {
                            let element = symox::matrixexpr::MatrixExpr::explicit_symbolic_element(
                                self.m.clone(), &i, &j, (host_int(&a), host_int(&b))).map_err(|e| err(py, e))?;
                            Ok(wrap(element).into_pyobject(py)?.into_any().unbind())
                        }
                        _ => Err(error),
                    };
                }
            };
            if !rs && !cs {
                return Ok(wrap(self.m.get(ri[0], ci[0]).clone()).into_pyobject(py)?.into_any().unbind());
            }
            let sub = self.m.extract(&ri, &ci).py(py)?;
            return Ok(self.wrapped(sub).into_pyobject(py)?.into_any().unbind());
        }
        let (ix, is_slice) = key_indices(key, self.m.len())?;
        if is_slice {
            let v: Vec<Expr> = ix.iter().map(|&i| self.m.data()[i].clone()).collect();
            return Ok(exprs_list(py, &v)?.into_any().unbind());
        }
        Ok(wrap(self.m.data()[ix[0]].clone()).into_pyobject(py)?.into_any().unbind())
    }

    fn __setitem__(&mut self, key: &Bound<'_, PyAny>, value: &Bound<'_, PyAny>) -> PyResult<()> {
        self.require_mutable()?;
        let (r, c) = self.m.shape();
        let (ri, ci, flat) = if let Ok(t) = key.downcast::<PyTuple>() {
            if t.len() != 2 {
                return Err(PyIndexError::new_err("expected (row, col)"));
            }
            let (ri, _) = key_indices(&t.get_item(0)?, r)?;
            let (ci, _) = key_indices(&t.get_item(1)?, c)?;
            (ri, ci, false)
        } else {
            let (ix, _) = key_indices(key, self.m.len())?;
            (ix, vec![], true)
        };
        let val = scalar_or_matrix(value)?;
        if flat {
            match val {
                Err(e) => {
                    for i in ri {
                        self.m.set(i / c, i % c, e.clone());
                    }
                }
                Ok(mv) => {
                    if mv.len() != ri.len() {
                        return Err(PyValueError::new_err("size mismatch in assignment"));
                    }
                    for (k, i) in ri.into_iter().enumerate() {
                        self.m.set(i / c, i % c, mv.data()[k].clone());
                    }
                }
            }
            return Ok(());
        }
        match val {
            Err(e) => {
                for &i in &ri {
                    for &j in &ci {
                        self.m.set(i, j, e.clone());
                    }
                }
            }
            Ok(mv) => {
                let mv = if mv.shape() == (ri.len(), ci.len()) {
                    mv
                } else if mv.len() == ri.len() * ci.len() {
                    mv.reshape(ri.len(), ci.len()).map_err(|e| PyValueError::new_err(e.to_string()))?
                } else {
                    return Err(PyValueError::new_err("size mismatch in assignment"));
                };
                for (a, &i) in ri.iter().enumerate() {
                    for (b, &j) in ci.iter().enumerate() {
                        self.m.set(i, j, mv.get(a, b).clone());
                    }
                }
            }
        }
        Ok(())
    }

    fn row(&self, py: Python<'_>, i: isize) -> PyResult<PyMatrix> {
        let i = norm_index(i, self.m.rows()).map_err(|_| err(py, MatrixError::Index("Index out of range".into())))?;
        Ok(self.wrapped(self.m.row(i)))
    }

    fn col(&self, j: isize) -> PyResult<PyMatrix> {
        let j = norm_index(j, self.m.cols())?;
        Ok(self.wrapped(self.m.col(j)))
    }

    fn extract(&self, py: Python<'_>, rowsList: &Bound<'_, PyAny>, colsList: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        let (ri, _) = key_indices(rowsList, self.m.rows())?;
        let (ci, _) = key_indices(colsList, self.m.cols())?;
        Ok(self.wrapped(self.m.extract(&ri, &ci).py(py)?))
    }

    fn reshape(&self, py: Python<'_>, rows: usize, cols: usize) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.reshape(rows, cols).py(py)?))
    }

    fn row_del(&mut self, py: Python<'_>, row: isize) -> PyResult<()> {
        self.require_mutable()?;
        let i = norm_index(row, self.m.rows())?;
        self.m.row_del(i).py(py)
    }

    fn col_del(&mut self, py: Python<'_>, col: isize) -> PyResult<()> {
        self.require_mutable()?;
        let j = norm_index(col, self.m.cols())?;
        self.m.col_del(j).py(py)
    }

    fn row_insert(&self, py: Python<'_>, pos: isize, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.row_insert(pos, &to_matrix(other)?).py(py)?))
    }

    fn col_insert(&self, py: Python<'_>, pos: isize, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.col_insert(pos, &to_matrix(other)?).py(py)?))
    }

    fn row_join(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.row_join(&to_matrix(other)?).py(py)?))
    }

    fn col_join(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.col_join(&to_matrix(other)?).py(py)?))
    }

    fn row_swap(&mut self, i: isize, j: isize) -> PyResult<()> {
        self.require_mutable()?;
        let n = self.m.rows();
        self.m.row_swap(norm_index(i, n)?, norm_index(j, n)?);
        Ok(())
    }

    fn col_swap(&mut self, i: isize, j: isize) -> PyResult<()> {
        self.require_mutable()?;
        let n = self.m.cols();
        self.m.col_swap(norm_index(i, n)?, norm_index(j, n)?);
        Ok(())
    }

    #[pyo3(signature = (perm, orientation="rows", direction="forward"))]
    fn permute(&self, py: Python<'_>, perm: &Bound<'_, PyAny>, orientation: &str, direction: &str) -> PyResult<PyMatrix> {
        Ok(self.wrapped(structural::permute(py, &self.m, perm, orientation, direction)?))
    }
    #[pyo3(signature = (swaps, direction="forward"))]
    fn permute_rows(&self, py: Python<'_>, swaps: &Bound<'_, PyAny>, direction: &str) -> PyResult<PyMatrix> {
        self.permute(py, swaps, "rows", direction)
    }
    #[pyo3(signature = (swaps, direction="forward"))]
    fn permute_cols(&self, py: Python<'_>, swaps: &Bound<'_, PyAny>, direction: &str) -> PyResult<PyMatrix> {
        self.permute(py, swaps, "cols", direction)
    }
    #[pyo3(name = "permuteFwd")]
    fn permute_fwd(&self, py: Python<'_>, perm: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        self.permute(py, perm, "rows", "forward")
    }
    #[pyo3(name = "permuteBkwd")]
    fn permute_bkwd(&self, py: Python<'_>, perm: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        self.permute(py, perm, "rows", "backward")
    }
    #[pyo3(signature = (op="n->kn", row=None, k=None, row1=None, row2=None))]
    fn elementary_row_op(&self, py: Python<'_>, op: &str, row: Option<&Bound<'_, PyAny>>, k: Option<&Bound<'_, PyAny>>, row1: Option<&Bound<'_, PyAny>>, row2: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(structural::elementary(py, &self.m, true, op, row, k, row1, row2)?))
    }
    #[pyo3(signature = (op="n->kn", col=None, k=None, col1=None, col2=None))]
    fn elementary_col_op(&self, py: Python<'_>, op: &str, col: Option<&Bound<'_, PyAny>>, k: Option<&Bound<'_, PyAny>>, col1: Option<&Bound<'_, PyAny>>, col2: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(structural::elementary(py, &self.m, false, op, col, k, col1, col2)?))
    }
    fn row_op(slf: Py<Self>, py: Python<'_>, i: isize, f: &Bound<'_, PyAny>) -> PyResult<()> {
        structural::callback(py,slf,i,None,f,true)
    }
    fn col_op(slf: Py<Self>, py: Python<'_>, j: isize, f: &Bound<'_, PyAny>) -> PyResult<()> {
        structural::callback(py,slf,j,None,f,false)
    }
    fn zip_row_op(slf: Py<Self>, py: Python<'_>, i: isize, k: isize, f: &Bound<'_, PyAny>) -> PyResult<()> {
        structural::callback(py,slf,i,Some(k),f,true)
    }
    #[staticmethod]
    fn _structural_plan(py: Python<'_>, specification: Vec<Vec<num_bigint::BigInt>>, dimension: num_bigint::BigInt, cycles: bool, backward: bool) -> PyResult<Bound<'_, PyDict>> {
        structural::plan(py,specification,dimension,cycles,backward)
    }
    #[staticmethod]
    fn _structural_limit(py: Python<'_>, count: num_bigint::BigInt) -> PyResult<()> { structural::limit(py,count) }
    #[staticmethod]
    fn _structural_values<'py>(py: Python<'py>, add: bool, coefficient: &Bound<'py, PyAny>, targets: &Bound<'py, PyAny>, sources: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyList>> {
        structural::values(py,add,coefficient,targets,sources)
    }

    // ------------------------------------------------------------------
    // auxiliary MatrixBase methods (kernels in symox::matrix::auxiliary)
    #[pyo3(signature = (k=None))]
    fn lower_triangular(&self, py: Python<'_>, k: Option<&Bound<'_,PyAny>>) -> PyResult<PyMatrix> { Ok(self.wrapped(auxiliary::triangular(py,&self.m,k,false)?)) }
    #[pyo3(signature = (k=None))]
    fn upper_triangular(&self, py: Python<'_>, k: Option<&Bound<'_,PyAny>>) -> PyResult<PyMatrix> { Ok(self.wrapped(auxiliary::triangular(py,&self.m,k,true)?)) }
    #[pyo3(signature = (k=1))]
    fn rot90(&self, k: i64) -> PyMatrix { self.wrapped(self.m.rot90(k)) }
    fn _vech(&self, py: Python<'_>, diagonal: bool) -> PyResult<PyMatrix> { Ok(self.wrapped(self.m.vech(diagonal).py(py)?)) }
    fn hat(&self, py: Python<'_>) -> PyResult<PyMatrix> { Ok(self.wrapped(self.m.hat().py(py)?)) }
    fn vee(&self, py: Python<'_>) -> PyResult<PyMatrix> { Ok(self.wrapped(self.m.vee().py(py)?)) }
    fn dual(&self, py: Python<'_>) -> PyResult<PyMatrix> { Ok(wrapm(self.m.dual().py(py)?)) }
    fn per(&self, py: Python<'_>) -> PyResult<PyExpr> { Ok(wrap(self.m.permanent().py(py)?)) }
    fn berkowitz<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> { auxiliary::berkowitz(py, &self.m) }
    fn berkowitz_minors<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.m.berkowitz_minors().py(py)?.into_iter().map(wrap))
    }
    fn _lu_ff(&self, py: Python<'_>) -> PyResult<(PyMatrix, PyMatrix, PyMatrix, PyMatrix)> {
        let (p, l, d, u) = self.m.lu_decomposition_ff().py(py)?;
        Ok((wrapm(p), wrapm(l), wrapm(d), wrapm(u)))
    }
    #[pyo3(signature = (iszerofunc=None, simpfunc=None, rankcheck=false))]
    fn LUdecomposition_Simple<'py>(&self, py: Python<'py>, iszerofunc: Option<&Bound<'_, PyAny>>, simpfunc: Option<&Bound<'_, PyAny>>, rankcheck: bool) -> PyResult<(PyMatrix, Bound<'py, PyList>)> {
        let (lu, swaps) = auxiliary::lu_simple(py,&self.m,iszerofunc,simpfunc,rankcheck)?;
        let lu = if self.m.data().is_empty() { self.wrapped(lu) } else { wrapm(lu) };
        Ok((lu, PyList::new(py, swaps.iter().map(|(a, b)| vec![*a, *b]))?))
    }
    fn is_nilpotent(&self, py: Python<'_>) -> PyResult<bool> { self.m.is_nilpotent().py(py) }
    fn is_symbolic(&self) -> bool { self.m.is_symbolic() }
    fn get_diag_blocks(&self) -> PyResult<Vec<PyMatrix>> {
        if self.m.rows() != 0 && self.m.cols() == 0 {
            return Err(PyIndexError::new_err("Index out of range"));
        }
        Ok(self.m.diagonal_blocks().into_iter().map(|b| self.wrapped(b)).collect())
    }
    fn connected_components(&self, py: Python<'_>) -> PyResult<Vec<Vec<usize>>> { self.m.connected_components().py(py) }
    #[pyo3(signature = (iszerofunc=None))]
    fn inverse_BLOCK(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        let _ = iszerofunc; // SymPy ignores it and always uses _iszero.
        Ok(self.wrapped(self.m.inverse_block().py(py)?))
    }
    fn upper_hessenberg_decomposition(&self, py: Python<'_>) -> PyResult<(PyMatrix, PyMatrix)> {
        let (h, p) = self.m.upper_hessenberg_decomposition().py(py)?;
        Ok((wrapm(h), wrapm(p)))
    }
    #[pyo3(signature = (upper=None))]
    fn bidiagonal_decomposition(&self, py: Python<'_>, upper: Option<&Bound<'_, PyAny>>) -> PyResult<(PyMatrix, PyMatrix, PyMatrix)> {
        let upper = match upper { Some(u) => auxiliary::upper_flag(u)?, None => true };
        let (u, b, v) = auxiliary::bidiagonal(py, &self.m, upper)?;
        Ok((wrapm(u), wrapm(b), wrapm(v)))
    }
    #[pyo3(signature = (upper=None))]
    fn bidiagonalize(&self, py: Python<'_>, upper: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        let upper = match upper { Some(u) => auxiliary::upper_flag(u)?, None => true };
        Ok(wrapm(auxiliary::bidiagonalize(py, &self.m, upper)?))
    }
    fn rref_rhs(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>) -> PyResult<(PyMatrix, PyMatrix)> {
        let right = to_matrix(rhs)?;
        let rational = self.m.data().iter().chain(right.data()).all(|e| e.as_number().is_some_and(|n| !n.is_float()));
        let (a, b) = self.m.rref_rhs(&right).py(py)?;
        let peer = !rational && self.with_peer(rhs, Matrix::zeros(0, 0)).immutable;
        Ok((PyMatrix { m: a, immutable: peer }, PyMatrix { m: b, immutable: peer }))
    }
    #[pyo3(signature = (iszerofunc=None))]
    fn normalized(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(auxiliary::normalized(py, &self.m, iszerofunc)?))
    }
    #[pyo3(signature = (simplify=None))]
    fn is_anti_symmetric(&self, py: Python<'_>, simplify: Option<&Bound<'_, PyAny>>) -> PyResult<Option<bool>> {
        match simplify {
            Some(s) => auxiliary::is_anti_symmetric(py, &self.m, s),
            None => auxiliary::is_anti_symmetric(py, &self.m, PyBool::new(py, true).as_any()),
        }
    }
    #[pyo3(signature = (rhs, det_method=None))]
    fn cramer_solve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>, det_method: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        let laplace = pyo3::types::PyString::new(py, "laplace");
        let method = det_method.filter(|m| !m.is_none()).cloned().unwrap_or_else(|| laplace.into_any());
        Ok(self.wrapped(auxiliary::cramer_solve(py, self, rhs, &method)?))
    }
    #[pyo3(signature = (other, failing_expression=false))]
    fn equals(&self, py: Python<'_>, other: &Bound<'_, PyAny>, failing_expression: bool) -> PyResult<PyObject> {
        auxiliary::equals(py, &self.m, other, failing_expression)
    }
    #[staticmethod]
    fn _companion(py: Python<'_>, coeffs: &Bound<'_, PyAny>) -> PyResult<PyMatrix> { Ok(wrapm(auxiliary::companion(py, coeffs)?)) }
    #[staticmethod]
    fn _jordan_block(py: Python<'_>, size: usize, eigenvalue: &Bound<'_, PyAny>, lower: bool) -> PyResult<PyMatrix> {
        Ok(wrapm(Matrix::jordan_block(size, &crate::scalar_expr(eigenvalue)?, lower).py(py)?))
    }
    #[staticmethod]
    fn _wilkinson(py: Python<'_>, n: usize) -> PyResult<(PyMatrix, PyMatrix)> { let (a, b) = Matrix::wilkinson(n).py(py)?; Ok((wrapm(a), wrapm(b))) }
    #[staticmethod]
    fn _orthogonalize(py: Python<'_>, vecs: &Bound<'_, PyAny>, normalize: bool, rankcheck: bool) -> PyResult<Vec<PyMatrix>> {
        Ok(auxiliary::orthogonalize(py, vecs, normalize, rankcheck)?.into_iter().map(wrapm).collect())
    }
    #[staticmethod]
    fn _irregular(py: Python<'_>, ntop: usize, matrices: &Bound<'_, PyList>) -> PyResult<PyMatrix> { Ok(wrapm(auxiliary::irregular(py, ntop, matrices)?)) }
    #[staticmethod]
    #[pyo3(signature = (seqs, n, zero=true))]
    fn _casoratian(py: Python<'_>, seqs: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>, zero: bool) -> PyResult<PyExpr> { Ok(wrap(auxiliary::casoratian(py, seqs, n, zero)?)) }

    #[staticmethod]
    #[pyo3(signature = (*args))]
    fn hstack(py: Python<'_>, args: &Bound<'_, PyTuple>) -> PyResult<PyMatrix> {
        let ms: Vec<Matrix> = args.iter().map(|a| to_matrix(&a)).collect::<PyResult<_>>()?;
        let refs: Vec<&Matrix> = ms.iter().collect();
        let immutable = args.iter().any(|arg| arg.downcast::<PyMatrix>().is_ok_and(|m| m.borrow().immutable));
        Ok(PyMatrix { m: Matrix::hstack(&refs).py(py)?, immutable })
    }

    #[staticmethod]
    #[pyo3(signature = (*args))]
    fn vstack(py: Python<'_>, args: &Bound<'_, PyTuple>) -> PyResult<PyMatrix> {
        let ms: Vec<Matrix> = args.iter().map(|a| to_matrix(&a)).collect::<PyResult<_>>()?;
        let refs: Vec<&Matrix> = ms.iter().collect();
        let immutable = args.iter().any(|arg| arg.downcast::<PyMatrix>().is_ok_and(|m| m.borrow().immutable));
        Ok(PyMatrix { m: Matrix::vstack(&refs).py(py)?, immutable })
    }

    // ------------------------------------------------------------------
    // structure
    // ------------------------------------------------------------------

    #[getter]
    #[pyo3(name = "T")]
    fn t(&self) -> PyMatrix {
        self.wrapped(self.m.transpose())
    }

    #[getter]
    #[pyo3(name = "H")]
    fn h(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(noncommutative_adjoint(py, &self.m)?))
    }

    #[getter]
    #[pyo3(name = "C")]
    fn c(&self) -> PyMatrix {
        self.wrapped(self.m.conjugate())
    }

    fn transpose(&self) -> PyMatrix {
        self.wrapped(self.m.transpose())
    }

    fn conjugate(&self) -> PyMatrix {
        self.wrapped(self.m.conjugate())
    }

    #[pyo3(signature = (deep=true, **hints))]
    fn as_real_imag(&self,deep:bool,hints:Option<&Bound<'_,PyDict>>)->PyResult<(PyMatrix,PyMatrix)>{
        crate::complex_parts_options(hints,false)?;
        let parts=symox::complex_parts::as_real_imag_many(self.m.data(),deep)
            .map_err(|error|PyNotImplementedError::new_err(error.to_string()))?;
        let(real,imag):(Vec<_>,Vec<_>)=parts.into_iter().unzip();
        let real=Matrix::new(self.m.rows(),self.m.cols(),real).map_err(|error|PyValueError::new_err(error.to_string()))?;
        let imag=Matrix::new(self.m.rows(),self.m.cols(),imag).map_err(|error|PyValueError::new_err(error.to_string()))?;
        Ok((self.wrapped(real),self.wrapped(imag)))
    }

    fn adjoint(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(noncommutative_adjoint(py, &self.m)?))
    }

    fn trace(&self, py: Python<'_>) -> PyResult<PyExpr> {
        Ok(wrap(self.m.trace().py(py)?))
    }

    #[pyo3(signature = (k=None))]
    fn diagonal(&self, k: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        let k = match k {
            None => 0i64,
            Some(k) => {
                let index = k.py().import("operator")?.call_method1("index", (k,));
                match index {
                    Ok(i) => i.extract::<i64>()?,
                    Err(e) if e.is_instance_of::<PyTypeError>(k.py()) => {
                        return Err(PyValueError::new_err(format!("{} is not an integer", k.str()?)))
                    }
                    Err(e) => return Err(e),
                }
            }
        };
        let (rows, cols) = (self.m.rows() as i64, self.m.cols() as i64);
        if k < 1 - rows || k > cols - 1 {
            return Err(PyValueError::new_err(format!("\nThe {k} diagonal is out of range [{}, {}]",1 - rows,cols - 1)));
        }
        let (r0, c0) = if k > 0 { (0, k) } else { (-k, 0) };
        let n = (rows - r0).min(cols - c0);
        if n <= 0 {
            return Err(PyValueError::new_err(format!(
                "\nThe {k} diagonal is out of range [{}, {}]",
                1 - rows,
                cols - 1
            )));
        }
        let items = (0..n).map(|i| self.m.get((r0 + i) as usize, (c0 + i) as usize).clone()).collect();
        Ok(self.wrapped(Matrix::row_vector(items)))
    }

    // ------------------------------------------------------------------
    // arithmetic
    // ------------------------------------------------------------------

    fn __neg__(&self) -> PyMatrix {
        self.wrapped(self.m.neg())
    }

    fn __pos__(&self) -> PyMatrix {
        self.clone()
    }

    fn __add__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if ndim_array(o) { return Ok(py.NotImplemented()); }
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, self.m.add(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __radd__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if ndim_array(o) { return Ok(py.NotImplemented()); }
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, m.add(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __sub__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if ndim_array(o) { return Ok(py.NotImplemented()); }
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, self.m.sub(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __rsub__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if ndim_array(o) { return Ok(py.NotImplemented()); }
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, m.sub(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __mul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if o.is_instance_of::<crate::matrixexpr::PyMatrixExpr>() { return Ok(py.NotImplemented()); }
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, self.m.matmul(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            Ok(Err(e)) => Ok(self.with_peer(o, self.m.scale(&e)).into_pyobject(py)?.into_any().unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }

    fn __rmul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        if o.is_instance_of::<crate::matrixexpr::PyMatrixExpr>() { return Ok(py.NotImplemented()); }
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, m.matmul(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            Ok(Err(e)) => Ok(self.with_peer(o, self.m.scale_left(&e)).into_pyobject(py)?.into_any().unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }

    fn __matmul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, self.m.matmul(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __rmatmul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, m.matmul(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __truediv__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Err(e)) => {
                let inv = symox::ops::pow(&e, &Expr::minus_one());
                Ok(self.wrapped(self.m.scale(&inv)).into_pyobject(py)?.into_any().unbind())
            }
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __pow__(&self, py: Python<'_>, o: &Bound<'_, PyAny>, _m: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        if _m.is_some_and(|m| !m.is_none()) {
            return Err(PyTypeError::new_err("matrix powers do not support a modulus"));
        }
        let e = to_expr(o)?;
        let n = match e.as_number() {
            Some(n) if n.is_integer() => n.to_i64().ok_or_else(|| PyValueError::new_err("exponent too large"))?,
            Some(n) if n.is_float() && n.to_f64().fract() == 0.0 => {
                let value = n.to_f64();
                if value < i64::MIN as f64 || value >= -(i64::MIN as f64) {
                    return Err(PyValueError::new_err("exponent too large"));
                }
                value as i64
            }
            _ => {
                return Err(PyNotImplementedError::new_err(
                    "only integer powers of matrices are implemented",
                ))
            }
        };
        Ok(self.wrapped(self.m.powi(n).py(py)?))
    }

    fn multiply_elementwise(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.with_peer(o, self.m.multiply_elementwise(&to_matrix(o)?).py(py)?))
    }

    fn multiply(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.with_peer(o, self.m.matmul(&to_matrix(o)?).py(py)?))
    }

    fn __richcmp__(&self, py: Python<'_>, other: &Bound<'_, PyAny>, op: pyo3::basic::CompareOp) -> PyResult<PyObject> {
        use pyo3::basic::CompareOp::*;
        if other.is_instance_of::<crate::matrixexpr::PyMatrixExpr>() { return Ok(py.NotImplemented()); }
        if other.getattr("_symox_sparse").is_ok_and(|value| value.extract::<bool>().unwrap_or(false)) {
            return Ok(py.NotImplemented());
        }
        let eq = match other.downcast::<PyMatrix>() {
            Ok(o) => o.borrow().m == self.m,
            Err(_) => {
                if is_sympy_obj(other) && other.hasattr("shape")? {
                    to_matrix(other).map(|m| m == self.m).unwrap_or(false)
                } else {
                    // RepMatrix.__eq__ returns NotImplemented for anything
                    // that is not a matrix, so the other operand's __eq__
                    // runs (a RawMatrix, a NumPy array, ...).
                    return Ok(py.NotImplemented());
                }
            }
        };
        Ok(match op {
            Eq => eq.into_pyobject(py)?.to_owned().into_any().unbind(),
            Ne => (!eq).into_pyobject(py)?.to_owned().into_any().unbind(),
            _ => py.NotImplemented(),
        })
    }

    fn __hash__(&self) -> PyResult<isize> {
        if !self.immutable {
            return Err(PyTypeError::new_err("unhashable type: 'Matrix'"));
        }
        let entries = self.m.data().iter().enumerate().filter(|(_, value)| !value.is_zero())
            .map(|(index, value)| ((index / self.m.cols(), index % self.m.cols()), value.clone())).collect();
        Ok(matrix_value_hash(self.m.rows(), self.m.cols(), entries))
    }

    #[getter]
    fn _immutable(&self) -> bool {
        self.immutable
    }

    // ------------------------------------------------------------------
    // element-wise
    // ------------------------------------------------------------------

    fn applyfunc(&self, f: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.try_applyfunc(|e| scalar_matrix_entry(&f.call1((wrap(e.clone()),))?))?))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn subs(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        if self.m.is_empty() {return Ok(self.wrapped(self.m.clone()));}
        let (args,options)=crate::substitution_arguments(args,_kwargs)?;
        let pairs = subs_pairs(&args)?;
        let simultaneous=crate::substitution_options(options.as_ref())?;
        Ok(self.wrapped(self.m.try_applyfunc(|entry|crate::checked_subs(args.py(),entry,&pairs,simultaneous))?))
    }

    fn xreplace(&self, rule: &Bound<'_, PyDict>) -> PyResult<PyMatrix> {
        let mut m = FxHashMap::default();
        for (k, v) in rule.iter() {
            m.insert(to_expr(&k)?, to_expr(&v)?);
        }
        crate::with_derivative_errors(||Ok(self.wrapped(self.m.xreplace(&m))))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn diff(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        let vars = diff_vars(args)?;
        crate::with_derivative_errors(||Ok(self.wrapped(self.m.diff_many(&vars))))
    }

    #[pyo3(signature = (deep=true, **hints))]
    fn expand(&self, deep: bool, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        let (options,mut budget,complex)=crate::parse_scalar_expand_options(deep,hints)?;
        budget.claim_terms(self.m.data().len()).map_err(|e|PyValueError::new_err(e.to_string()))?;
        Ok(self.wrapped(self.m.try_applyfunc(|e| crate::expand_scalar(e,&options,&mut budget,complex))?))
    }

    #[pyo3(signature = (n=15, subs=None, **_kwargs))]
    fn evalf(&self, n: i64, subs: Option<&Bound<'_, PyDict>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        let precision=crate::decimal_precision(n)?;
        let settings=crate::evalf_options(_kwargs,precision)?;
        Ok(self.wrapped(self.m.try_applyfunc(|e| crate::evalf_prepared(e,precision,subs,settings))?))
    }

    #[pyo3(signature = (n=15, subs=None, **_kwargs))]
    fn n(&self, n: i64, subs: Option<&Bound<'_, PyDict>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        self.evalf(n, subs, _kwargs)
    }

    /// Simplify the entries in place (like SymPy's mutable matrices).
    #[pyo3(signature = (**_kwargs))]
    fn simplify(&mut self, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Option<PyMatrix>> {
        reject_kwargs(_kwargs)?;
        Ok(if self.immutable {
            Some(self.wrapped(self.m.simplify()))
        } else {
            self.m = self.m.simplify();
            None
        })
    }

    fn cancel(&self) -> PyMatrix {
        self.wrapped(self.m.cancel())
    }

    // ------------------------------------------------------------------
    // linear algebra
    // ------------------------------------------------------------------

    fn _has_rational_entries(&self) -> bool {
        matrix::domain::has_rational_entries(&self.m)
    }

    fn _has_rational_polynomial_entries(&self) -> bool {
        matrix::domain::has_rational_polynomial_entries(&self.m)
    }

    #[pyo3(signature = (method="bareiss", iszerofunc=None))]
    fn det(&self, py: Python<'_>, method: &str, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        let m = DetMethod::from_name(method).py(py)?;
        Ok(wrap(auxiliary::determinant(py,&self.m,m,iszerofunc)?))
    }

    #[pyo3(signature = (method=None, iszerofunc=None, try_block_diag=false))]
    fn inv(
        &self,
        py: Python<'_>,
        method: Option<&str>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        try_block_diag: bool,
    ) -> PyResult<PyMatrix> {
        Ok(self.wrapped(auxiliary::inverse(py,&self.m,method,iszerofunc,try_block_diag)?))
    }

    fn _inexact_domain(&self) -> Option<(bool,u32)> {
        matrix::inexact::domain(&self.m).map(|domain|(domain.complex,domain.precision))
    }

    fn _inverse_inexact(&self, py: Python<'_>, complex: bool, precision: u32) -> PyResult<PyMatrix> {
        let domain=matrix::inexact::InexactDomain{complex,precision};
        let result=matrix::inexact::inverse(&self.m,domain,||py.check_signals());
        match result {
            Ok(matrix)=>Ok(self.wrapped(matrix)),
            Err(matrix::inexact::InexactError::Matrix(error))=>Err(err(py,error)),
            Err(matrix::inexact::InexactError::Interrupted(error))=>Err(error),
        }
    }

    #[pyo3(signature=(delta=LllDelta::default()), text_signature="($self, delta=0.75)")]
    fn lll(&self, py: Python<'_>, delta: LllDelta) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.lll(&delta.0).map_err(|e|lll_error(py,e))?))
    }

    #[pyo3(signature=(delta=LllDelta::default()), text_signature="($self, delta=0.75)")]
    fn lll_transform<'py>(&self, py: Python<'py>, delta: LllDelta) -> PyResult<Bound<'py, PyTuple>> {
        let (basis,transform)=self.m.lll_transform(&delta.0).map_err(|e|lll_error(py,e))?;
        PyTuple::new(py,[self.wrapped(basis),self.wrapped(transform)])
    }

    #[staticmethod]
    #[pyo3(signature=(rows, cols, delta=LllDelta::default()))]
    fn _lll_preflight(py: Python<'_>, rows: num_bigint::BigInt, cols: num_bigint::BigInt, delta: LllDelta) -> PyResult<()> {
        matrix::lll::preflight(0,0,&delta.0).map_err(|e|lll_error(py,e))?;
        if rows>cols { return Err(lll_error(py,matrix::lll::LllError::Shape)); }
        let limited=||lll_error(py,matrix::lll::LllError::Resource);
        let rows=rows.to_usize().ok_or_else(limited)?;
        let cols=cols.to_usize().ok_or_else(limited)?;
        matrix::lll::preflight(rows,cols,&delta.0).map_err(|e|lll_error(py,e))
    }

    fn inv_mod(&self, py: Python<'_>, m: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        if !self.m.is_square() { return Err(err(py, MatrixError::NonSquare(String::new()))); }
        let modulus = modular_modulus(py,m)?;
        Ok(self.wrapped(self.m.inv_mod(&modulus).py(py)?))
    }

    #[staticmethod]
    fn _inv_mod_preflight(py: Python<'_>, rows: num_bigint::BigInt, cols: num_bigint::BigInt, m: &Bound<'_, PyAny>) -> PyResult<()> {
        if rows != cols { return Err(err(py, MatrixError::NonSquare(String::new()))); }
        let modulus = modular_modulus(py,m)?;
        matrix::inv_mod::preflight(0,0,&modulus).py(py)?;
        let n = rows.to_usize().ok_or_else(|| PyNotImplementedError::new_err("Exact modular matrix inverse resource limit exceeded"))?;
        matrix::inv_mod::preflight(n,n,&modulus).py(py)
    }

    #[pyo3(signature = (iszerofunc=None))]
    fn inverse_GE(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        self.inv(py,Some("GE"),iszerofunc,false)
    }

    #[pyo3(signature = (iszerofunc=None))]
    fn inverse_LU(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        self.inv(py,Some("LU"),iszerofunc,false)
    }

    #[pyo3(signature = (iszerofunc=None))]
    fn inverse_ADJ(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        self.inv(py,Some("ADJ"),iszerofunc,false)
    }

    #[pyo3(signature = (iszerofunc=None))]
    fn inverse_CH(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        self.inv(py, Some("CH"), iszerofunc, false)
    }

    #[pyo3(signature = (iszerofunc=None))]
    fn inverse_LDL(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        self.inv(py, Some("LDL"), iszerofunc, false)
    }

    #[pyo3(signature = (iszerofunc=None))]
    fn inverse_QR(&self, py: Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        self.inv(py, Some("QR"), iszerofunc, false)
    }

    #[pyo3(signature = (method="berkowitz"))]
    fn adjugate(&self, py: Python<'_>, method: &str) -> PyResult<PyMatrix> {
        let _ = method;
        Ok(self.wrapped(self.m.adjugate().py(py)?))
    }

    #[pyo3(signature = (i, j, method="berkowitz"))]
    fn minor(&self, py: Python<'_>, i: isize, j: isize, method: &str) -> PyResult<PyExpr> {
        let _ = method;
        let (i, j) = self.minor_index(py, i, j, false)?;
        Ok(wrap(self.m.minor(i, j).py(py)?))
    }

    #[pyo3(signature = (i, j, method="berkowitz"))]
    fn cofactor(&self, py: Python<'_>, i: isize, j: isize, method: &str) -> PyResult<PyExpr> {
        let _ = method;
        let (i, j) = self.minor_index(py, i, j, true)?;
        Ok(wrap(self.m.cofactor(i, j).py(py)?))
    }

    #[pyo3(signature = (method="berkowitz"))]
    fn cofactor_matrix(&self, py: Python<'_>, method: &str) -> PyResult<PyMatrix> {
        let _ = method;
        Ok(self.wrapped(self.m.cofactor_matrix().py(py)?))
    }

    fn minor_submatrix(&self, i: isize, j: isize) -> PyResult<PyMatrix> {
        let (i, j) = (norm_index(i, self.m.rows())?, norm_index(j, self.m.cols())?);
        Ok(self.wrapped(self.m.minor_submatrix(i, j)))
    }

    #[pyo3(signature = (iszerofunc=None, simplify=PivotSimplifier::Native, pivots=true, normalize_last=true),
        text_signature="($self, iszerofunc=None, simplify=False, pivots=True, normalize_last=True)")]
    fn rref(
        &self,
        py: Python<'_>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        simplify: PivotSimplifier,
        pivots: bool,
        normalize_last: bool,
    ) -> PyResult<PyObject> {
        let custom=simplify.callback(py);
        let (r,p)=if iszerofunc.is_none_or(|callback|callback.is_none()) && custom.is_none() && normalize_last {
            self.m.rref()
        } else {
            let (matrix,pivots,_)=auxiliary::row_reduce(py,&self.m,iszerofunc,custom,normalize_last,true,true)?;
            (matrix,pivots)
        };
        if pivots {
            let t = PyTuple::new(py, p)?;
            Ok(PyTuple::new(py, [self.wrapped(r).into_pyobject(py)?.into_any(), t.into_any()])?.into_any().unbind())
        } else {
            Ok(self.wrapped(r).into_pyobject(py)?.into_any().unbind())
        }
    }

    #[pyo3(signature = (iszerofunc=None, simplify=PivotSimplifier::Native, with_pivots=false),
        text_signature="($self, iszerofunc=None, simplify=False, with_pivots=False)")]
    fn echelon_form(
        &self,
        py: Python<'_>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        simplify: PivotSimplifier,
        with_pivots: bool,
    ) -> PyResult<PyObject> {
        let custom=simplify.callback(py);
        let (r,p)=if iszerofunc.is_none_or(|callback|callback.is_none()) && custom.is_none() {
            self.m.echelon_form()
        } else {
            let (matrix,pivots,_)=auxiliary::row_reduce(py,&self.m,iszerofunc,custom,true,false,false)?;
            (matrix,pivots)
        };
        if with_pivots {
            let t = PyTuple::new(py, p)?;
            Ok(PyTuple::new(py, [self.wrapped(r).into_pyobject(py)?.into_any(), t.into_any()])?.into_any().unbind())
        } else {
            Ok(self.wrapped(r).into_pyobject(py)?.into_any().unbind())
        }
    }

    #[pyo3(signature = (iszerofunc=None, simplify=PivotSimplifier::Native),
        text_signature="($self, iszerofunc=None, simplify=False)")]
    fn rank(&self, py:Python<'_>, iszerofunc: Option<&Bound<'_, PyAny>>, simplify: PivotSimplifier) -> PyResult<usize> {
        let custom=simplify.callback(py);
        if iszerofunc.is_none_or(|callback|callback.is_none()) && custom.is_none() {return Ok(self.m.rank());}
        Ok(auxiliary::row_reduce(py,&self.m,iszerofunc,custom,true,false,false)?.1.len())
    }

    #[pyo3(signature = (simplify=PivotSimplifier::Native, iszerofunc=None),
        text_signature="($self, simplify=False, iszerofunc=None)")]
    fn nullspace(&self, py:Python<'_>, simplify: PivotSimplifier, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<Vec<PyMatrix>> {
        let custom=simplify.callback(py);
        if iszerofunc.is_none_or(|callback|callback.is_none()) && custom.is_none() {
            return Ok(self.m.nullspace().into_iter().map(|m| self.wrapped(m)).collect());
        }
        let (reduced,pivots,_)=auxiliary::row_reduce(py,&self.m,iszerofunc,custom,true,true,true)?;
        let mut basis=Vec::new();
        for free in (0..self.m.cols()).filter(|column|!pivots.contains(column)) {
            let mut entries=vec![Expr::zero();self.m.cols()];entries[free]=Expr::one();
            for (row,column) in pivots.iter().enumerate() {entries[*column]=symox::ops::neg(reduced.get(row,free));}
            basis.push(self.wrapped(Matrix::column(entries)));
        }
        Ok(basis)
    }

    #[pyo3(signature = (simplify=PivotSimplifier::Native),text_signature="($self, simplify=False)")]
    fn columnspace(&self, py:Python<'_>, simplify: PivotSimplifier) -> PyResult<Vec<PyMatrix>> {
        match simplify.callback(py) {
            None=>Ok(self.m.columnspace().into_iter().map(|m|self.wrapped(m)).collect()),
            custom=>{
                let pivots=auxiliary::row_reduce(py,&self.m,None,custom,true,false,false)?.1;
                Ok(pivots.iter().map(|column|self.wrapped(self.m.col(*column))).collect())
            }
        }
    }

    #[pyo3(signature = (simplify=PivotSimplifier::Native),text_signature="($self, simplify=False)")]
    fn rowspace(&self, py:Python<'_>, simplify: PivotSimplifier) -> PyResult<Vec<PyMatrix>> {
        match simplify.callback(py) {
            None=>Ok(self.m.rowspace().into_iter().map(|m|self.wrapped(m)).collect()),
            custom=>{
                let (reduced,pivots,_)=auxiliary::row_reduce(py,&self.m,None,custom,true,false,false)?;
                Ok((0..pivots.len()).map(|row|self.wrapped(reduced.row(row))).collect())
            }
        }
    }

    #[pyo3(signature = (iszerofunc=None, simpfunc=None, rankcheck=false))]
    fn LUdecomposition<'py>(
        &self,
        py: Python<'py>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        simpfunc: Option<&Bound<'_, PyAny>>,
        rankcheck: bool,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let (l, u, p) = auxiliary::lu_decomposition(py,&self.m,iszerofunc,simpfunc,rankcheck)?;
        let perm = PyList::new(py, p.iter().map(|(a, b)| vec![*a, *b]))?;
        PyTuple::new(
            py,
            [
                self.wrapped(l).into_pyobject(py)?.into_any(),
                self.wrapped(u).into_pyobject(py)?.into_any(),
                perm.into_any(),
            ],
        )
    }

    #[pyo3(signature = (rhs, iszerofunc=None))]
    fn LUsolve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(auxiliary::lu_solve(py,&self.m,&to_matrix(rhs)?,iszerofunc)?))
    }

    #[pyo3(signature = (iszerofunc=None, simplify=PivotSimplifier::Native),
        text_signature="($self, iszerofunc=None, simplify=False)")]
    fn rank_decomposition<'py>(&self, py: Python<'py>, iszerofunc: Option<&Bound<'_, PyAny>>, simplify: PivotSimplifier) -> PyResult<Bound<'py, PyTuple>> {
        let custom=simplify.callback(py);
        let (c,f)=if iszerofunc.is_none_or(|callback|callback.is_none()) && custom.is_none() {
            self.m.rank_decomposition().py(py)?
        } else {
            let (reduced,pivots,_)=auxiliary::row_reduce(py,&self.m,iszerofunc,custom,true,true,true)?;
            let rows=(0..self.m.rows()).collect::<Vec<_>>();
            let columns=(0..self.m.cols()).collect::<Vec<_>>();
            (self.m.extract(&rows,&pivots).py(py)?,reduced.extract(&(0..pivots.len()).collect::<Vec<_>>(),&columns).py(py)?)
        };
        PyTuple::new(py, [self.wrapped(c), self.wrapped(f)])
    }

    #[pyo3(signature = (method="RD"))]
    fn pinv(&self, py: Python<'_>, method: &str) -> PyResult<PyMatrix> {
        let result = match method {
            "RD" => self.m.pinv(),
            "ED" => self.m.pinv_ed(),
            _ => return Err(PyValueError::new_err(format!("invalid pinv method {method}"))),
        };
        Ok(self.wrapped(result.py(py)?))
    }

    #[pyo3(signature = (B, arbitrary_matrix=None))]
    fn pinv_solve(&self, py: Python<'_>, B: &Bound<'_, PyAny>, arbitrary_matrix: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        let b = to_matrix(B)?;
        let arbitrary = arbitrary_matrix.filter(|w| !w.is_none()).map(to_matrix).transpose()?;
        Ok(self.wrapped(self.m.pinv_solve(&b, arbitrary.as_ref()).py(py)?))
    }

    fn lower_triangular_solve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.lower_triangular_solve(&to_matrix(rhs)?).py(py)?))
    }

    fn upper_triangular_solve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.upper_triangular_solve(&to_matrix(rhs)?).py(py)?))
    }

    fn diagonal_solve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.diagonal_solve(&to_matrix(rhs)?).py(py)?))
    }

    fn QRsolve(&self, py: Python<'_>, b: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.qr_solve(&to_matrix(b)?).py(py)?))
    }

    fn cholesky_solve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.cholesky_solve(&to_matrix(rhs)?).py(py)?))
    }

    fn LDLsolve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.ldl_solve(&to_matrix(rhs)?).py(py)?))
    }

    #[pyo3(signature = (rhs, method="CH"))]
    fn solve_least_squares(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>, method: &str) -> PyResult<PyMatrix> {
        let b = to_matrix(rhs)?;
        let result = match method.to_uppercase().as_str() {
            "CH" => self.m.cholesky_solve(&b),
            "LDL" => self.m.ldl_solve(&b),
            "QR" => self.m.qr_solve(&b),
            "PINV" => self.m.pinv_solve(&b, None),
            "GJ" | "GE" | "LU" | "ADJ" | "DM" | "DMNC" => {
                let h = self.m.adjoint();
                let normal = h.matmul(&self.m).py(py)?;
                let rhs = h.matmul(&b).py(py)?;
                match method.to_uppercase().as_str() {
                    "LU" => normal.lu_solve(&rhs),
                    "ADJ" => normal.inverse_adj().and_then(|i| i.matmul(&rhs)),
                    "DM" | "DMNC" => normal.inv().and_then(|i| i.matmul(&rhs)),
                    _ => normal.solve(&rhs),
                }
            }
            _ => return Err(PyNotImplementedError::new_err(format!("matrix least-squares method {method} is not implemented"))),
        };
        Ok(self.wrapped(result.py(py)?))
    }

    #[pyo3(signature = (rhs, method="GJ"))]
    fn solve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>, method: &str) -> PyResult<PyMatrix> {
        let b = to_matrix(rhs)?;
        let r = match method.to_uppercase().as_str() {
            "GJ" | "GE" => self.m.solve(&b),
            "LU" => self.m.lu_solve(&b),
            "CH" => self.m.cholesky_solve(&b),
            "LDL" => self.m.ldl_solve(&b),
            "QR" => self.m.qr_solve(&b),
            "PINV" => self.m.pinv_solve(&b, None),
            "CRAMER" => return self.cramer_solve(py, rhs, None),
            "ADJ" => self.m.inverse_adj().and_then(|i| i.matmul(&b)),
            "DM" | "DMNC" => self.m.inv().and_then(|i| i.matmul(&b)),
            "BLOCK" => self.inverse_BLOCK(py, None)?.m.matmul(&b),
            _ => return Err(PyValueError::new_err("Inversion method unrecognized")),
        };
        Ok(self.wrapped(r.py(py)?))
    }

    #[pyo3(signature = (B, freevar=false))]
    fn gauss_jordan_solve<'py>(&self, py: Python<'py>, B: &Bound<'_, PyAny>, freevar: bool) -> PyResult<Bound<'py, PyTuple>> {
        let (s, t) = self.m.gauss_jordan_solve(&to_matrix(B)?).py(py)?;
        let mut values = vec![self.wrapped(s).into_pyobject(py)?.into_any(), self.wrapped(t).into_pyobject(py)?.into_any()];
        if freevar {
            let pivots = self.m.rref().1;
            let free: Vec<usize> = (0..self.m.cols()).filter(|c| !pivots.contains(c)).collect();
            values.push(PyList::new(py, free)?.into_any());
        }
        PyTuple::new(py, values)
    }

    /// Characteristic polynomial ``det(lambda*I - M)`` as a ``PurePoly``.
    #[pyo3(signature = (x=None, simplify=None))]
    fn charpoly<'py>(
        &self,
        py: Python<'py>,
        x: Option<&Bound<'_, PyAny>>,
        simplify: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        reject_callback("simplify", simplify)?;
        // Preserve a supplied symbol's identity unless its name clashes
        // with a matrix coefficient. The native constructor chooses a fresh
        // name in that case so coefficients cannot capture the generator.
        let mut gen: Option<Expr> = None;
        let name = match x {
            None => "lambda".to_string(),
            Some(o) => {
                if let Ok(s) = o.extract::<String>() {
                    s
                } else {
                    let e = to_expr(o)?;
                    if e.as_symbol().is_none() {
                        return Err(PyTypeError::new_err("expected a symbol"));
                    }
                    let name = e.as_sym().unwrap().name().to_string();
                    gen = Some(e);
                    name
                }
            }
        };
        let cp = self.m.charpoly(&name).py(py)?;
        let coeffs = exprs_list(py, &cp.coeffs)?;
        let var = match gen {
            Some(g) if cp.var.as_sym().unwrap().name() == g.as_sym().unwrap().name() =>
                g,
            _ => cp.var.clone(),
        };
        let cls = py.import("symox")?.getattr("PurePoly")?;
        cls.call_method1("_from_coeffs", (coeffs, wrap(var), cp.domain.clone()))
    }

    #[pyo3(signature = (error_when_incomplete=true, multiple=false, **_flags))]
    fn eigenvals<'py>(&self, py: Python<'py>, error_when_incomplete: bool, multiple: bool, _flags: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
        reject_kwargs(_flags)?;
        let _ = error_when_incomplete;
        let ev = self.m.eigenvals().py(py)?;
        if multiple {
            let values: Vec<PyExpr> = ev.into_iter().flat_map(|(v, k)| std::iter::repeat_n(v, k)).map(wrap).collect();
            return Ok(PyList::new(py, values)?.into_any());
        }
        let d = PyDict::new(py);
        for (v, k) in ev {
            d.set_item(wrap(v), k)?;
        }
        Ok(d.into_any())
    }

    #[pyo3(signature = (error_when_incomplete=true, **_flags))]
    fn eigenvects<'py>(&self, py: Python<'py>, error_when_incomplete: bool, _flags: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyList>> {
        reject_kwargs(_flags)?;
        let _ = error_when_incomplete;
        let ev = self.m.eigenvects().py(py)?;
        let items: Vec<Bound<'py, PyTuple>> = ev
            .into_iter()
            .map(|(v, k, vs)| {
                let l = PyList::new(py, vs.into_iter().map(|m| self.wrapped(m)))?;
                PyTuple::new(
                    py,
                    [
                        wrap(v).into_pyobject(py)?.into_any(),
                        k.into_pyobject(py)?.into_any(),
                        l.into_any(),
                    ],
                )
            })
            .collect::<PyResult<_>>()?;
        PyList::new(py, items)
    }

    fn jacobian(&self, py: Python<'_>, X: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        let vars = to_matrix(X)?;
        Ok(self.wrapped(self.m.jacobian(vars.data()).py(py)?))
    }

    #[pyo3(signature = (ord=None))]
    fn norm(&self, py: Python<'_>, ord: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        let o: Option<String> = match ord {
            None => None,
            Some(o) if o.is_none() => None,
            Some(o) => {
                if let Ok(value) = o.extract::<String>() {
                    if self.m.is_vector() {
                        return Err(PyValueError::new_err("Expected order to be Number, Symbol, oo"));
                    }
                    let value = value.to_lowercase();
                    if !matches!(value.as_str(), "f" | "fro" | "frobenius" | "vector") {
                        return Err(PyNotImplementedError::new_err("Matrix Norms under development"));
                    }
                    Some(value)
                } else {
                    let order = to_expr(o)?;
                    if self.m.is_vector() && !order.is_constant(symox::Constant::Infinity)
                        && !order.is_constant(symox::Constant::NegInfinity) {
                        let values = self.m.values();
                        let values = if values.is_empty() { vec![Expr::zero()] } else { values };
                        let terms: Vec<Expr> = values.iter().map(|e| {
                            let abs = symox::functions::func(symox::Func::Abs, vec![e.clone()]);
                            symox::ops::pow(&abs, &order)
                        }).collect();
                        let inverse = symox::ops::pow(&order, &Expr::minus_one());
                        return Ok(wrap(symox::ops::pow(&symox::ops::add(&terms), &inverse)));
                    }
                    Some(order.to_string())
                }
            }
        };
        Ok(wrap(self.m.norm_ord(o.as_deref()).py(py)?))
    }

    #[pyo3(signature = (hermitian=true))]
    fn cholesky(&self, py: Python<'_>, hermitian: bool) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.cholesky(hermitian).py(py)?))
    }

    #[pyo3(signature = (hermitian=true))]
    fn LDLdecomposition<'py>(&self, py: Python<'py>, hermitian: bool) -> PyResult<Bound<'py, PyTuple>> {
        let (l, d) = self.m.ldl_decomposition(hermitian).py(py)?;
        PyTuple::new(py, [self.wrapped(l), self.wrapped(d)])
    }

    fn QRdecomposition<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (q, r) = self.m.qr_decomposition().py(py)?;
        PyTuple::new(py, [self.wrapped(q), self.wrapped(r)])
    }

    fn singular_value_decomposition<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (u, s, v) = self.m.singular_value_decomposition().py(py)?;
        PyTuple::new(py, [self.wrapped(u), self.wrapped(s), self.wrapped(v)])
    }

    fn singular_values(&self, py: Python<'_>) -> PyResult<Vec<PyExpr>> {
        Ok(self.m.singular_values().py(py)?.into_iter().map(wrap).collect())
    }

    fn condition_number(&self, py: Python<'_>) -> PyResult<PyExpr> {
        Ok(wrap(self.m.condition_number().py(py)?))
    }

    fn exp(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.exp().py(py)?))
    }

    /// Matrix logarithm using scalar log branches. Explicit simplify callables
    /// postprocess the matrix once; False/None omit optional final cancellation.
    #[pyo3(signature = (simplify=LogSimplifier::Native), text_signature = "($self, simplify=cancel)")]
    fn log(&self, py: Python<'_>, simplify: LogSimplifier) -> PyResult<PyMatrix> {
        let result = self.m.log_with_simplification(matches!(simplify, LogSimplifier::Native)).py(py)?;
        match simplify {
            LogSimplifier::Callback(callback) => {
                let value = callback.call1(py, (self.wrapped(result),))?;
                Ok(self.wrapped(to_matrix(value.bind(py))?))
            }
            _ => Ok(self.wrapped(result)),
        }
    }

    fn analytic_func(&self, py: Python<'_>, f: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        let function=crate::scalar_expr(f)?;
        let parameter=crate::scalar_expr(x)?;
        Ok(self.wrapped(self.m.analytic_func(&function,&parameter).py(py)?))
    }

    #[pyo3(signature = (calc_transform=true, **kwargs))]
    fn jordan_form<'py>(&self, py: Python<'py>, calc_transform: bool, kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyAny>> {
        reject_kwargs(kwargs)?;
        let (p, j, _) = self.m.jordan_decomposition().py(py)?;
        if calc_transform {
            Ok(PyTuple::new(py, [self.wrapped(p), self.wrapped(j)])?.into_any())
        } else {
            Ok(self.wrapped(j).into_pyobject(py)?.into_any())
        }
    }

    #[pyo3(signature = (calc_transformation=true))]
    fn jordan_cells<'py>(&self, py: Python<'py>, calc_transformation: bool) -> PyResult<Bound<'py, PyTuple>> {
        let _ = calc_transformation;
        let (p, _, blocks) = self.m.jordan_decomposition().py(py)?;
        let cells = PyList::new(py, blocks.into_iter().map(|b|self.wrapped(b)))?;
        PyTuple::new(py, [self.wrapped(p).into_pyobject(py)?.into_any(), cells.into_any()])
    }

    #[pyo3(signature = (reals_only=false, sort=false, normalize=false))]
    fn diagonalize<'py>(&self, py: Python<'py>, reals_only: bool, sort: bool, normalize: bool) -> PyResult<Bound<'py, PyTuple>> {
        let (p, d) = self.m.diagonalize_with_options(reals_only, sort, normalize).py(py)?;
        PyTuple::new(py, [self.wrapped(p), self.wrapped(d)])
    }

    #[pyo3(signature = (reals_only=false, **_kwargs))]
    fn is_diagonalizable(&self, py: Python<'_>, reals_only: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<bool> {
        reject_kwargs(_kwargs)?;
        self.m.is_diagonalizable(reals_only).py(py)
    }

    #[getter]
    fn is_hermitian(&self) -> Option<bool> {
        self.m.is_hermitian()
    }

    fn vec(&self) -> PyMatrix {
        self.wrapped(self.m.vec())
    }

    fn values(&self) -> Vec<PyExpr> {
        self.m.values().into_iter().map(wrap).collect()
    }

    fn flat(&self) -> Vec<PyExpr> {
        self.m.data().iter().cloned().map(wrap).collect()
    }

    fn as_mutable(&self) -> PyMatrix {
        wrapm(self.m.clone())
    }

    fn as_immutable(&self) -> PyMatrix {
        PyMatrix { m: self.m.clone(), immutable: true }
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, Bound<'py, PyTuple>)> {
        if self.m.data().iter().any(Expr::has_matrix) {
            return crate::matrix_dag::reduce_dense(py,&self.m,self.immutable);
        }
        let cls = py.import("symox")?.getattr(if self.immutable { "ImmutableMatrix" } else { "Matrix" })?;
        let data: Vec<PyExpr> = self.m.data().iter().cloned().map(wrap).collect();
        Ok((cls, (self.m.rows(), self.m.cols(), data).into_pyobject(py)?))
    }

    #[staticmethod]
    #[pyo3(name = "eye", signature = (n, m=None, **_kwargs))]
    fn eye_s(n: &Bound<'_, PyAny>, m: Option<&Bound<'_, PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
        eye(n, m, _kwargs)
    }

    #[staticmethod]
    #[pyo3(name = "zeros", signature = (r, c=None, **_kwargs))]
    fn zeros_s(r: &Bound<'_, PyAny>, c: Option<&Bound<'_, PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
        zeros(r, c, _kwargs)
    }

    #[staticmethod]
    #[pyo3(name = "ones", signature = (r, c=None, **_kwargs))]
    fn ones_s(r: &Bound<'_, PyAny>, c: Option<&Bound<'_, PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
        ones(r, c, _kwargs)
    }

    #[staticmethod]
    #[pyo3(name = "diag", signature = (*values, **_kwargs))]
    fn diag_s(values: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        diag(values, _kwargs)
    }

    #[pyo3(signature = (b, hermitian=None, conjugate_convention=None))]
    fn dot(
        &self,
        py: Python<'_>,
        b: &Bound<'_, PyAny>,
        hermitian: Option<bool>,
        conjugate_convention: Option<&str>,
    ) -> PyResult<PyExpr> {
        let other = to_matrix(b)?;
        let hermitian = hermitian.unwrap_or(conjugate_convention.is_some());
        let result = if hermitian {
            match conjugate_convention.unwrap_or("maths") {
                "math" | "maths" | "left" => self.m.conjugate().dot(&other),
                "physics" | "right" => self.m.dot(&other.conjugate()),
                _ => return Err(PyValueError::new_err("Unknown conjugate_convention; expected math, maths, left, physics or right")),
            }
        } else {
            self.m.dot(&other)
        };
        Ok(wrap(result.py(py)?))
    }

    fn cross(&self, py: Python<'_>, b: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.cross(&to_matrix(b)?).py(py)?))
    }

    fn strongly_connected_components(&self, py: Python<'_>) -> PyResult<Vec<Vec<usize>>> {
        if !self.m.is_square() {
            return Err(err(py, MatrixError::NonSquare(String::new())));
        }
        Ok(matrix::linalg::scc(&self.m))
    }
}

// ----------------------------------------------------------------------
// module functions
// ----------------------------------------------------------------------

/// SymPy's ``as_int`` on a matrix dimension.
fn as_int_dim(o: &Bound<'_, PyAny>) -> PyResult<isize> {
    let py = o.py();
    match py.import("operator")?.getattr("index")?.call1((o,)) {
        Ok(v) => v.extract::<isize>(),
        Err(e) if e.is_instance_of::<PyTypeError>(py) => {
            Err(PyValueError::new_err(format!("{} is not an integer", o.str()?)))
        }
        Err(e) => Err(e),
    }
}

/// Dimensions of ``zeros``/``ones``/``eye``/``Matrix(r, c, data)``: integers
/// (``as_int``) and non-negative, with SymPy's error messages.
fn creation_dims(r: &Bound<'_, PyAny>, c: Option<&Bound<'_, PyAny>>) -> PyResult<(usize, usize)> {
    let r = as_int_dim(r)?;
    let c = match c {
        Some(c) if !c.is_none() => as_int_dim(c)?,
        _ => r,
    };
    if r < 0 || c < 0 {
        return Err(PyValueError::new_err(format!(
            "Cannot create a {} x {} matrix. Both dimensions must be positive", r, c)));
    }
    Ok((r as usize, c as usize))
}

#[pyfunction]
#[pyo3(signature = (n, m=None, **_kwargs))]
fn eye(n: &Bound<'_, PyAny>, m: Option<&Bound<'_, PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (r, c) = creation_dims(n, m)?;
    finish_special(n.py(), wrapm(Matrix::from_fn(r, c, |i, j| if i == j { Expr::one() } else { Expr::zero() })), _kwargs)
}

#[pyfunction]
#[pyo3(signature = (r, c=None, **_kwargs))]
fn zeros(r: &Bound<'_, PyAny>, c: Option<&Bound<'_, PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let py = r.py();
    let (r, c) = creation_dims(r, c)?;
    finish_special(py, wrapm(Matrix::zeros(r, c)), _kwargs)
}

#[pyfunction]
#[pyo3(signature = (r, c=None, **_kwargs))]
fn ones(r: &Bound<'_, PyAny>, c: Option<&Bound<'_, PyAny>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let py = r.py();
    let (r, c) = creation_dims(r, c)?;
    finish_special(py, wrapm(Matrix::ones(r, c)), _kwargs)
}

#[pyfunction]
#[pyo3(signature = (*values, **_kwargs))]
fn diag(values: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
    reject_kwargs(_kwargs)?;
    let mut blocks: Vec<Matrix> = Vec::new();
    for v in values.iter() {
        if v.downcast::<PyMatrix>().is_ok() || is_seq(&v) || v.getattr("is_Matrix").is_ok_and(|v|v.is_truthy().unwrap_or(false)) || (is_sympy_obj(&v) && v.hasattr("shape")?) {
            blocks.push(to_matrix(&v)?);
        } else {
            blocks.push(Matrix::new(1, 1, vec![scalar_matrix_entry(&v)?]).unwrap());
        }
    }
    Ok(wrapm(Matrix::diag_blocks(&blocks)))
}

fn expr_seq(obj: &Bound<'_, PyAny>) -> PyResult<Vec<Expr>> {
    if let Ok(m) = obj.downcast::<PyMatrix>() {
        return Ok(m.borrow().m.data().to_vec());
    }
    if is_sympy_obj(obj) && obj.hasattr("shape")? {
        return Ok(to_matrix(obj)?.data().to_vec());
    }
    if let Ok(it) = obj.try_iter() {
        return it.map(|x| to_expr(&x?)).collect();
    }
    Ok(vec![to_expr(obj)?])
}

#[pyfunction]
#[pyo3(name = "jacobian")]
fn py_jacobian(exprs: &Bound<'_, PyAny>, vars: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
    let fs = expr_seq(exprs)?;
    let vs = expr_seq(vars)?;
    Ok(wrapm(matrix::jacobian(&fs, &vs)))
}

#[pyfunction]
#[pyo3(name = "hessian", signature = (f, varlist, constraints=None))]
fn py_hessian(f: &Bound<'_, PyAny>, varlist: &Bound<'_, PyAny>, constraints: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
    if constraints.is_some_and(|c| !c.is_none()) {
        return Err(PyNotImplementedError::new_err("bordered Hessians are not implemented"));
    }
    let vs = expr_seq(varlist)?;
    Ok(wrapm(matrix::hessian(&to_expr(f)?, &vs)))
}

#[pyfunction]
#[pyo3(name = "wronskian", signature = (functions, var, method="bareiss"))]
fn py_wronskian(py: Python<'_>, functions: &Bound<'_, PyAny>, var: &Bound<'_, PyAny>, method: &str) -> PyResult<PyExpr> {
    let _ = method;
    let fs = expr_seq(functions)?;
    Ok(wrap(matrix::wronskian(&fs, &to_expr(var)?).py(py)?))
}

/// `cancel` of a single expression (used by the Python layer).
#[pyfunction]
#[pyo3(name = "_mat_cancel")]
fn py_cancel(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(matrix::ratfunc::cancel(&to_expr(e)?)))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMatrix>()?;
    // Internal modules reach the native base as `_core.Matrix`.
    m.add("Matrix", m.py().get_type::<PyMatrix>())?;
    m.add_function(wrap_pyfunction!(_register_matrix_types, m)?)?;
    m.add_function(wrap_pyfunction!(eye, m)?)?;
    m.add_function(wrap_pyfunction!(zeros, m)?)?;
    m.add_function(wrap_pyfunction!(ones, m)?)?;
    m.add_function(wrap_pyfunction!(diag, m)?)?;
    m.add_function(wrap_pyfunction!(py_jacobian, m)?)?;
    m.add_function(wrap_pyfunction!(py_hessian, m)?)?;
    m.add_function(wrap_pyfunction!(py_wronskian, m)?)?;
    m.add_function(wrap_pyfunction!(py_cancel, m)?)?;
    m.add_function(wrap_pyfunction!(_sparse_matrix_hash, m)?)?;
    Ok(())
}


/// Conjugate transpose; noncommutative entries take the reference
/// `adjoint(entry)` (reversed products) instead of an entrywise conjugate.
fn noncommutative_adjoint(py: Python<'_>, m: &Matrix) -> PyResult<Matrix> {
    if !m.data().iter().any(symox::decompose::has_noncommutative_atoms) {
        return Ok(m.adjoint());
    }
    let adjoint = py.import("symox._nc_adjoint")?.getattr("adjoint_of")?;
    m.transpose().try_applyfunc(|e| crate::to_expr(&adjoint.call1((crate::wrap(e.clone()),))?))
}
