//! Python bindings for symbolic matrices (`symox.Matrix`).

use pyo3::exceptions::{PyIndexError, PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyInt, PyList, PySlice, PyTuple};
use rustc_hash::FxHashMap;

use symox::expr::Expr;
use symox::matrix::{self, DetMethod, Matrix, MatrixError};
use symox::printer;

use crate::{diff_vars, subs_pairs, to_expr, wrap, PyExpr};

fn err(py: Python<'_>, e: MatrixError) -> PyErr {
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

trait IntoPy<T> {
    fn py(self, py: Python<'_>) -> PyResult<T>;
}

impl<T> IntoPy<T> for Result<T, MatrixError> {
    fn py(self, py: Python<'_>) -> PyResult<T> {
        self.map_err(|e| err(py, e))
    }
}

/// A mutable dense symbolic matrix.
#[pyclass(name = "Matrix", module = "symox._core", weakref)]
#[derive(Clone)]
pub struct PyMatrix {
    pub m: Matrix,
    immutable: bool,
}

fn wrapm(m: Matrix) -> PyMatrix {
    PyMatrix { m, immutable: false }
}

fn matrix_value_hash(rows: usize, cols: usize, mut entries: Vec<((usize, usize), Expr)>) -> isize {
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
            Err(PyTypeError::new_err("Cannot set values of an immutable matrix"))
        } else {
            Ok(())
        }
    }
}

fn exprs_list<'py>(py: Python<'py>, v: &[Expr]) -> PyResult<Bound<'py, PyList>> {
    PyList::new(py, v.iter().map(|e| wrap(e.clone())))
}

fn is_sympy_obj(obj: &Bound<'_, PyAny>) -> bool {
    obj.get_type()
        .module()
        .map(|m| m.to_string().starts_with("sympy"))
        .unwrap_or(false)
}

/// Convert a Python object to a matrix (symox/SymPy matrices, nested
/// sequences, numpy arrays).
pub fn to_matrix(obj: &Bound<'_, PyAny>) -> PyResult<Matrix> {
    if let Ok(m) = obj.downcast::<PyMatrix>() {
        return Ok(m.borrow().m.clone());
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

fn is_seq(obj: &Bound<'_, PyAny>) -> bool {
    obj.downcast::<PyList>().is_ok() || obj.downcast::<PyTuple>().is_ok()
}

fn from_nested(obj: &Bound<'_, PyAny>) -> PyResult<Matrix> {
    let items: Vec<Bound<'_, PyAny>> = obj.try_iter()?.collect::<PyResult<_>>()?;
    if items.is_empty() {
        return Ok(Matrix::zeros(0, 0));
    }
    let is_matrix = |item: &Bound<'_, PyAny>| -> bool {
        item.downcast::<PyMatrix>().is_ok()
            || (is_sympy_obj(item) && item.hasattr("shape").unwrap_or(false)
                && item.hasattr("tolist").unwrap_or(false))
    };
    // A flat list of matrices is stacked vertically. Nested rows of
    // matrices are block rows and are joined horizontally first.
    if items.iter().any(&is_matrix) {
        let blocks: Vec<Matrix> = items.iter().map(|item| {
            if is_matrix(item) { to_matrix(item) }
            else { Ok(Matrix::column(vec![to_expr(item)?])) }
        }).collect::<PyResult<_>>()?;
        let refs: Vec<&Matrix> = blocks.iter().collect();
        return Matrix::vstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()));
    }
    if items.iter().all(|i| is_seq(i) || (i.hasattr("tolist").unwrap_or(false) && !is_sympy_obj(i) && i.hasattr("shape").unwrap_or(false))) {
        let mut rows = Vec::with_capacity(items.len());
        for it in &items {
            let it = if is_seq(it) { it.clone() } else { it.call_method0("tolist")? };
            let entries: Vec<Bound<'_, PyAny>> = it.try_iter()?.collect::<PyResult<_>>()?;
            if entries.iter().any(&is_matrix) {
                let blocks: Vec<Matrix> = entries.iter().map(|entry| {
                    if is_matrix(entry) { to_matrix(entry) }
                    else { Ok(Matrix::column(vec![to_expr(entry)?])) }
                }).collect::<PyResult<_>>()?;
                let refs: Vec<&Matrix> = blocks.iter().collect();
                rows.push(Matrix::hstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()))?);
            } else {
                let row: Vec<Expr> = entries.iter().map(to_expr).collect::<PyResult<_>>()?;
                rows.push(Matrix::new(1, row.len(), row).map_err(|e| PyValueError::new_err(e.to_string()))?);
            }
        }
        let refs: Vec<&Matrix> = rows.iter().collect();
        return Matrix::vstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()));
    }
    let col: Vec<Expr> = items.iter().map(to_expr).collect::<PyResult<_>>()?;
    Ok(Matrix::column(col))
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
        let mut v = Vec::new();
        for x in k.try_iter()? {
            let i: isize = x?.extract()?;
            v.push(norm_index(i, n)?);
        }
        return Ok((v, true));
    }
    Err(PyTypeError::new_err("invalid index"))
}

fn scalar_or_matrix(o: &Bound<'_, PyAny>) -> PyResult<Result<Matrix, Expr>> {
    if let Ok(m) = o.downcast::<PyMatrix>() {
        return Ok(Ok(m.borrow().m.clone()));
    }
    if is_seq(o) {
        return Ok(Ok(from_nested(o)?));
    }
    if is_sympy_obj(o) && o.hasattr("shape")? && o.hasattr("tolist")? {
        return Ok(Ok(to_matrix(o)?));
    }
    Ok(Err(to_expr(o)?))
}

#[pymethods]
impl PyMatrix {
    #[new]
    #[pyo3(signature = (*args))]
    fn py_new(py: Python<'_>, args: &Bound<'_, PyTuple>) -> PyResult<Self> {
        match args.len() {
            0 => Ok(wrapm(Matrix::zeros(0, 0))),
            1 => Ok(wrapm(to_matrix(&args.get_item(0)?)?)),
            3 => {
                let r: usize = args.get_item(0)?.extract()?;
                let c: usize = args.get_item(1)?.extract()?;
                let src = args.get_item(2)?;
                if src.is_callable() && !is_seq(&src) {
                    let mut data = Vec::with_capacity(r * c);
                    for i in 0..r {
                        for j in 0..c {
                            data.push(to_expr(&src.call1((i, j))?)?);
                        }
                    }
                    return Ok(wrapm(Matrix::new(r, c, data).py(py)?));
                }
                let mut flat: Vec<Expr> = Vec::with_capacity(r * c);
                for it in src.try_iter()? {
                    let it = it?;
                    if is_seq(&it) {
                        for x in it.try_iter()? {
                            flat.push(to_expr(&x?)?);
                        }
                    } else {
                        flat.push(to_expr(&it)?);
                    }
                }
                Ok(wrapm(Matrix::new(r, c, flat).py(py)?))
            }
            2 => {
                let r: usize = args.get_item(0)?.extract()?;
                let c: usize = args.get_item(1)?.extract()?;
                Ok(wrapm(Matrix::zeros(r, c)))
            }
            _ => Err(PyTypeError::new_err("Matrix() takes at most 3 arguments")),
        }
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
        let ns = PyDict::new(py);
        ns.set_item("__builtins__", py.import("builtins")?)?;
        py.run(&std::ffi::CString::new("from sympy import *")?, Some(&ns), None)?;
        let (r, c) = self.m.shape();
        let name = if self.immutable { "ImmutableMatrix" } else { "Matrix" };
        let mut s = format!("{name}({r}, {c}, [");
        for (i, e) in self.m.data().iter().enumerate() {
            if i > 0 {
                s.push_str(", ");
            }
            s.push_str(&printer::srepr(e));
        }
        s.push_str("])");
        py.eval(&std::ffi::CString::new(s)?, Some(&ns), None)
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
            let (ri, rs) = key_indices(&t.get_item(0)?, r)?;
            let (ci, cs) = key_indices(&t.get_item(1)?, c)?;
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

    fn row_del(&mut self, py: Python<'_>, i: isize) -> PyResult<()> {
        self.require_mutable()?;
        let i = norm_index(i, self.m.rows())?;
        self.m.row_del(i).py(py)
    }

    fn col_del(&mut self, py: Python<'_>, j: isize) -> PyResult<()> {
        self.require_mutable()?;
        let j = norm_index(j, self.m.cols())?;
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
    fn h(&self) -> PyMatrix {
        self.wrapped(self.m.adjoint())
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

    fn adjoint(&self) -> PyMatrix {
        self.wrapped(self.m.adjoint())
    }

    fn trace(&self, py: Python<'_>) -> PyResult<PyExpr> {
        Ok(wrap(self.m.trace().py(py)?))
    }

    fn diagonal(&self) -> PyMatrix {
        self.wrapped(self.m.diagonal())
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
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, self.m.add(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __radd__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, m.add(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __sub__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, self.m.sub(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __rsub__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, m.sub(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __mul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(self.with_peer(o, self.m.matmul(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            Ok(Err(e)) => Ok(self.with_peer(o, self.m.scale(&e)).into_pyobject(py)?.into_any().unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }

    fn __rmul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
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
        if other.getattr("_symox_sparse").is_ok_and(|value| value.extract::<bool>().unwrap_or(false)) {
            return Ok(py.NotImplemented());
        }
        let eq = match other.downcast::<PyMatrix>() {
            Ok(o) => o.borrow().m == self.m,
            Err(_) => {
                if is_sympy_obj(other) && other.hasattr("shape")? {
                    to_matrix(other).map(|m| m == self.m).unwrap_or(false)
                } else {
                    false
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
        Ok(self.wrapped(self.m.try_applyfunc(|e| to_expr(&f.call1((wrap(e.clone()),))?))?))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn subs(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        reject_kwargs(_kwargs)?;
        let pairs = subs_pairs(args)?;
        Ok(self.wrapped(self.m.subs(&pairs)))
    }

    fn xreplace(&self, rule: &Bound<'_, PyDict>) -> PyResult<PyMatrix> {
        let mut m = FxHashMap::default();
        for (k, v) in rule.iter() {
            m.insert(to_expr(&k)?, to_expr(&v)?);
        }
        Ok(self.wrapped(self.m.xreplace(&m)))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn diff(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        let vars = diff_vars(args)?;
        Ok(self.wrapped(self.m.diff_many(&vars)))
    }

    #[pyo3(signature = (deep=true, **hints))]
    fn expand(&self, deep: bool, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        if !deep {
            return Err(PyNotImplementedError::new_err("matrix expand(deep=False) is not implemented"));
        }
        Ok(self.wrapped(self.m.try_applyfunc(|e| crate::expand_hints(e, hints))?))
    }

    #[pyo3(signature = (n=15, subs=None, **_kwargs))]
    fn evalf(&self, n: u32, subs: Option<&Bound<'_, PyDict>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        crate::decimal_precision(n)?;
        crate::evalf_options(_kwargs)?;
        Ok(self.wrapped(self.m.try_applyfunc(|e| crate::evalf_expr(e, n, subs, _kwargs))?))
    }

    #[pyo3(signature = (n=15, subs=None, **_kwargs))]
    fn n(&self, n: u32, subs: Option<&Bound<'_, PyDict>>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
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

    #[pyo3(signature = (method="bareiss", iszerofunc=None))]
    fn det(&self, py: Python<'_>, method: &str, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        reject_callback("iszerofunc", iszerofunc)?;
        let m = DetMethod::from_name(method).py(py)?;
        Ok(wrap(self.m.det_method(m).py(py)?))
    }

    #[pyo3(signature = (method=None, iszerofunc=None, try_block_diag=false))]
    fn inv(
        &self,
        py: Python<'_>,
        method: Option<&str>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        try_block_diag: bool,
    ) -> PyResult<PyMatrix> {
        reject_callback("iszerofunc", iszerofunc)?;
        let _ = try_block_diag;
        if !self.m.is_square() {
            return Err(err(py, MatrixError::NonSquare("A Matrix must be square to invert.".into())));
        }
        let r = match method.map(|s| s.to_uppercase()) {
            None => self.m.inv(),
            Some(s) if s == "GE" => self.m.inverse_ge(),
            Some(s) if s == "LU" => self.m.inverse_lu(),
            Some(s) if s == "ADJ" => self.m.inverse_adj(),
            Some(s) if s == "CH" => self.m.cholesky_solve(&Matrix::eye(self.m.rows())),
            Some(s) if s == "LDL" => self.m.ldl_solve(&Matrix::eye(self.m.rows())),
            Some(s) if s == "QR" => self.m.qr_solve(&Matrix::eye(self.m.rows())),
            Some(s) if s == "DM" || s == "DMNC" => self.m.inv(),
            Some(s) => return Err(PyValueError::new_err(format!("Inversion method unrecognized: {s}"))),
        };
        Ok(self.wrapped(r.py(py)?))
    }

    fn inverse_GE(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.inverse_ge().py(py)?))
    }

    fn inverse_LU(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.inverse_lu().py(py)?))
    }

    fn inverse_ADJ(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(self.wrapped(self.m.inverse_adj().py(py)?))
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
        let (i, j) = (norm_index(i, self.m.rows())?, norm_index(j, self.m.cols())?);
        Ok(wrap(self.m.minor(i, j).py(py)?))
    }

    #[pyo3(signature = (i, j, method="berkowitz"))]
    fn cofactor(&self, py: Python<'_>, i: isize, j: isize, method: &str) -> PyResult<PyExpr> {
        let _ = method;
        let (i, j) = (norm_index(i, self.m.rows())?, norm_index(j, self.m.cols())?);
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

    #[pyo3(signature = (iszerofunc=None, simplify=false, pivots=true, normalize_last=true))]
    fn rref(
        &self,
        py: Python<'_>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        simplify: bool,
        pivots: bool,
        normalize_last: bool,
    ) -> PyResult<PyObject> {
        reject_callback("iszerofunc", iszerofunc)?;
        let _ = (simplify, normalize_last);
        let (r, p) = self.m.rref();
        if pivots {
            let t = PyTuple::new(py, p)?;
            Ok(PyTuple::new(py, [self.wrapped(r).into_pyobject(py)?.into_any(), t.into_any()])?.into_any().unbind())
        } else {
            Ok(self.wrapped(r).into_pyobject(py)?.into_any().unbind())
        }
    }

    #[pyo3(signature = (iszerofunc=None, simplify=false, with_pivots=false))]
    fn echelon_form(
        &self,
        py: Python<'_>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        simplify: bool,
        with_pivots: bool,
    ) -> PyResult<PyObject> {
        reject_callback("iszerofunc", iszerofunc)?;
        let _ = simplify;
        let (r, p) = self.m.echelon_form();
        if with_pivots {
            let t = PyTuple::new(py, p)?;
            Ok(PyTuple::new(py, [self.wrapped(r).into_pyobject(py)?.into_any(), t.into_any()])?.into_any().unbind())
        } else {
            Ok(self.wrapped(r).into_pyobject(py)?.into_any().unbind())
        }
    }

    #[pyo3(signature = (iszerofunc=None, simplify=false))]
    fn rank(&self, iszerofunc: Option<&Bound<'_, PyAny>>, simplify: bool) -> PyResult<usize> {
        reject_callback("iszerofunc", iszerofunc)?;
        let _ = simplify;
        Ok(self.m.rank())
    }

    #[pyo3(signature = (simplify=false, iszerofunc=None))]
    fn nullspace(&self, simplify: bool, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<Vec<PyMatrix>> {
        reject_callback("iszerofunc", iszerofunc)?;
        let _ = simplify;
        Ok(self.m.nullspace().into_iter().map(|m| self.wrapped(m)).collect())
    }

    #[pyo3(signature = (simplify=false))]
    fn columnspace(&self, simplify: bool) -> Vec<PyMatrix> {
        let _ = simplify;
        self.m.columnspace().into_iter().map(|m| self.wrapped(m)).collect()
    }

    #[pyo3(signature = (simplify=false))]
    fn rowspace(&self, simplify: bool) -> Vec<PyMatrix> {
        let _ = simplify;
        self.m.rowspace().into_iter().map(|m| self.wrapped(m)).collect()
    }

    #[pyo3(signature = (iszerofunc=None, simpfunc=None, rankcheck=false))]
    fn LUdecomposition<'py>(
        &self,
        py: Python<'py>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        simpfunc: Option<&Bound<'_, PyAny>>,
        rankcheck: bool,
    ) -> PyResult<Bound<'py, PyTuple>> {
        reject_callback("iszerofunc", iszerofunc)?;
        reject_callback("simpfunc", simpfunc)?;
        let (l, u, p) = self.m.lu_decomposition_with_rankcheck(rankcheck).py(py)?;
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
        reject_callback("iszerofunc", iszerofunc)?;
        Ok(self.wrapped(self.m.lu_solve(&to_matrix(rhs)?).py(py)?))
    }

    #[pyo3(signature = (iszerofunc=None, simplify=false))]
    fn rank_decomposition<'py>(&self, py: Python<'py>, iszerofunc: Option<&Bound<'_, PyAny>>, simplify: bool) -> PyResult<Bound<'py, PyTuple>> {
        reject_callback("iszerofunc", iszerofunc)?;
        let _ = simplify;
        let (c, f) = self.m.rank_decomposition().py(py)?;
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
            "ADJ" => self.m.inverse_adj().and_then(|i| i.matmul(&b)),
            "DM" | "DMNC" => self.m.inv().and_then(|i| i.matmul(&b)),
            _ => return Err(PyNotImplementedError::new_err(format!("matrix solve method {method} is not implemented"))),
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
        let name = match x {
            None => "lambda".to_string(),
            Some(o) => {
                if let Ok(s) = o.extract::<String>() {
                    s
                } else {
                    let e = to_expr(o)?;
                    e.as_symbol().map(|s| s.to_string()).ok_or_else(|| PyTypeError::new_err("expected a symbol"))?
                }
            }
        };
        let cp = self.m.charpoly(&name).py(py)?;
        let coeffs = exprs_list(py, &cp.coeffs)?;
        let cls = py.import("symox")?.getattr("PurePoly")?;
        cls.call_method1("_from_coeffs", (coeffs, wrap(cp.var.clone()), cp.domain.clone(), wrap(cp.as_expr())))
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

    fn singular_values(&self, py: Python<'_>) -> PyResult<Vec<PyExpr>> {
        Ok(self.m.singular_values().py(py)?.into_iter().map(wrap).collect())
    }

    fn condition_number(&self, py: Python<'_>) -> PyResult<PyExpr> {
        Ok(wrap(self.m.condition_number().py(py)?))
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

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, (usize, usize, Vec<PyExpr>))> {
        let cls = py.import("symox")?.getattr(if self.immutable { "ImmutableMatrix" } else { "Matrix" })?;
        // Serialize entries through Expr.__reduce__ to preserve assumptions,
        // held nodes and exact floating point payloads.
        let data: Vec<PyExpr> = self.m.data().iter().cloned().map(wrap).collect();
        Ok((cls, (self.m.rows(), self.m.cols(), data)))
    }

    #[staticmethod]
    #[pyo3(name = "eye", signature = (n, m=None, **_kwargs))]
    fn eye_s(n: usize, m: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        eye(n, m, _kwargs)
    }

    #[staticmethod]
    #[pyo3(name = "zeros", signature = (r, c=None, **_kwargs))]
    fn zeros_s(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        zeros(r, c, _kwargs)
    }

    #[staticmethod]
    #[pyo3(name = "ones", signature = (r, c=None, **_kwargs))]
    fn ones_s(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
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

    fn strongly_connected_components(&self) -> Vec<Vec<usize>> {
        matrix::linalg::scc(&self.m)
    }
}

// ----------------------------------------------------------------------
// module functions
// ----------------------------------------------------------------------

fn dims(r: usize, c: Option<usize>) -> (usize, usize) {
    (r, c.unwrap_or(r))
}

#[pyfunction]
#[pyo3(signature = (n, m=None, **_kwargs))]
fn eye(n: usize, m: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
    reject_kwargs(_kwargs)?;
    let (r, c) = dims(n, m);
    Ok(wrapm(Matrix::from_fn(r, c, |i, j| if i == j { Expr::one() } else { Expr::zero() })))
}

#[pyfunction]
#[pyo3(signature = (r, c=None, **_kwargs))]
fn zeros(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
    reject_kwargs(_kwargs)?;
    let (r, c) = dims(r, c);
    Ok(wrapm(Matrix::zeros(r, c)))
}

#[pyfunction]
#[pyo3(signature = (r, c=None, **_kwargs))]
fn ones(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
    reject_kwargs(_kwargs)?;
    let (r, c) = dims(r, c);
    Ok(wrapm(Matrix::ones(r, c)))
}

#[pyfunction]
#[pyo3(signature = (*values, **_kwargs))]
fn diag(values: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
    reject_kwargs(_kwargs)?;
    let mut blocks: Vec<Matrix> = Vec::new();
    for v in values.iter() {
        if v.downcast::<PyMatrix>().is_ok() || is_seq(&v) || (is_sympy_obj(&v) && v.hasattr("shape")?) {
            blocks.push(to_matrix(&v)?);
        } else {
            blocks.push(Matrix::new(1, 1, vec![to_expr(&v)?]).unwrap());
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
