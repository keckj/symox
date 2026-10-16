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
}

fn wrapm(m: Matrix) -> PyMatrix {
    PyMatrix { m }
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
    // list of matrices: stack them vertically (rows)
    if items.iter().all(|i| i.downcast::<PyMatrix>().is_ok()) {
        let ms: Vec<Matrix> = items.iter().map(|i| i.downcast::<PyMatrix>().unwrap().borrow().m.clone()).collect();
        if ms.iter().all(|m| m.rows() == 1) || ms.iter().all(|m| m.cols() != 1) {
            let refs: Vec<&Matrix> = ms.iter().collect();
            return Matrix::vstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()));
        }
        let refs: Vec<&Matrix> = ms.iter().collect();
        return Matrix::hstack(&refs).map_err(|e| PyValueError::new_err(e.to_string()));
    }
    if items.iter().all(|i| is_seq(i) || (i.hasattr("tolist").unwrap_or(false) && !is_sympy_obj(i) && i.hasattr("shape").unwrap_or(false))) {
        let mut rows = Vec::with_capacity(items.len());
        for it in &items {
            let it = if is_seq(it) { it.clone() } else { it.call_method0("tolist")? };
            let row: Vec<Expr> = it.try_iter()?.map(|x| to_expr(&x?)).collect::<PyResult<_>>()?;
            rows.push(row);
        }
        let c = rows[0].len();
        if rows.iter().any(|r| r.len() != c) {
            return Err(PyValueError::new_err("mismatched dimensions"));
        }
        if c == 0 {
            return Ok(Matrix::zeros(rows.len(), 0));
        }
        return Matrix::from_rows(rows).map_err(|e| PyValueError::new_err(e.to_string()));
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
        let mut s = format!("Matrix({r}, {c}, [");
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
    fn is_zero_matrix(&self) -> bool {
        self.m.is_zero_matrix()
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
            return Ok(wrapm(sub).into_pyobject(py)?.into_any().unbind());
        }
        let (ix, is_slice) = key_indices(key, self.m.len())?;
        if is_slice {
            let v: Vec<Expr> = ix.iter().map(|&i| self.m.data()[i].clone()).collect();
            return Ok(exprs_list(py, &v)?.into_any().unbind());
        }
        Ok(wrap(self.m.data()[ix[0]].clone()).into_pyobject(py)?.into_any().unbind())
    }

    fn __setitem__(&mut self, key: &Bound<'_, PyAny>, value: &Bound<'_, PyAny>) -> PyResult<()> {
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
        Ok(wrapm(self.m.row(i)))
    }

    fn col(&self, j: isize) -> PyResult<PyMatrix> {
        let j = norm_index(j, self.m.cols())?;
        Ok(wrapm(self.m.col(j)))
    }

    fn extract(&self, py: Python<'_>, rowsList: &Bound<'_, PyAny>, colsList: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        let (ri, _) = key_indices(rowsList, self.m.rows())?;
        let (ci, _) = key_indices(colsList, self.m.cols())?;
        Ok(wrapm(self.m.extract(&ri, &ci).py(py)?))
    }

    fn reshape(&self, py: Python<'_>, rows: usize, cols: usize) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.reshape(rows, cols).py(py)?))
    }

    fn row_del(&mut self, py: Python<'_>, i: isize) -> PyResult<()> {
        let i = norm_index(i, self.m.rows())?;
        self.m.row_del(i).py(py)
    }

    fn col_del(&mut self, py: Python<'_>, j: isize) -> PyResult<()> {
        let j = norm_index(j, self.m.cols())?;
        self.m.col_del(j).py(py)
    }

    fn row_insert(&self, py: Python<'_>, pos: isize, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.row_insert(pos, &to_matrix(other)?).py(py)?))
    }

    fn col_insert(&self, py: Python<'_>, pos: isize, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.col_insert(pos, &to_matrix(other)?).py(py)?))
    }

    fn row_join(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.row_join(&to_matrix(other)?).py(py)?))
    }

    fn col_join(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.col_join(&to_matrix(other)?).py(py)?))
    }

    fn row_swap(&mut self, i: isize, j: isize) -> PyResult<()> {
        let n = self.m.rows();
        self.m.row_swap(norm_index(i, n)?, norm_index(j, n)?);
        Ok(())
    }

    fn col_swap(&mut self, i: isize, j: isize) -> PyResult<()> {
        let n = self.m.cols();
        self.m.col_swap(norm_index(i, n)?, norm_index(j, n)?);
        Ok(())
    }

    #[staticmethod]
    #[pyo3(signature = (*args))]
    fn hstack(py: Python<'_>, args: &Bound<'_, PyTuple>) -> PyResult<PyMatrix> {
        let ms: Vec<Matrix> = args.iter().map(|a| to_matrix(&a)).collect::<PyResult<_>>()?;
        let refs: Vec<&Matrix> = ms.iter().collect();
        Ok(wrapm(Matrix::hstack(&refs).py(py)?))
    }

    #[staticmethod]
    #[pyo3(signature = (*args))]
    fn vstack(py: Python<'_>, args: &Bound<'_, PyTuple>) -> PyResult<PyMatrix> {
        let ms: Vec<Matrix> = args.iter().map(|a| to_matrix(&a)).collect::<PyResult<_>>()?;
        let refs: Vec<&Matrix> = ms.iter().collect();
        Ok(wrapm(Matrix::vstack(&refs).py(py)?))
    }

    // ------------------------------------------------------------------
    // structure
    // ------------------------------------------------------------------

    #[getter]
    #[pyo3(name = "T")]
    fn t(&self) -> PyMatrix {
        wrapm(self.m.transpose())
    }

    #[getter]
    #[pyo3(name = "H")]
    fn h(&self) -> PyMatrix {
        wrapm(self.m.adjoint())
    }

    #[getter]
    #[pyo3(name = "C")]
    fn c(&self) -> PyMatrix {
        wrapm(self.m.conjugate())
    }

    fn transpose(&self) -> PyMatrix {
        wrapm(self.m.transpose())
    }

    fn conjugate(&self) -> PyMatrix {
        wrapm(self.m.conjugate())
    }

    fn adjoint(&self) -> PyMatrix {
        wrapm(self.m.adjoint())
    }

    fn trace(&self, py: Python<'_>) -> PyResult<PyExpr> {
        Ok(wrap(self.m.trace().py(py)?))
    }

    fn diagonal(&self) -> PyMatrix {
        wrapm(self.m.diagonal())
    }

    // ------------------------------------------------------------------
    // arithmetic
    // ------------------------------------------------------------------

    fn __neg__(&self) -> PyMatrix {
        wrapm(self.m.neg())
    }

    fn __pos__(&self) -> PyMatrix {
        self.clone()
    }

    fn __add__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(self.m.add(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __radd__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(m.add(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __sub__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(self.m.sub(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __rsub__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(m.sub(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __mul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(self.m.matmul(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            Ok(Err(e)) => Ok(wrapm(self.m.scale(&e)).into_pyobject(py)?.into_any().unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }

    fn __rmul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(m.matmul(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            Ok(Err(e)) => Ok(wrapm(self.m.scale_left(&e)).into_pyobject(py)?.into_any().unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }

    fn __matmul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(self.m.matmul(&m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __rmatmul__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Ok(m)) => Ok(wrapm(m.matmul(&self.m).py(py)?).into_pyobject(py)?.into_any().unbind()),
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __truediv__(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        match scalar_or_matrix(o) {
            Ok(Err(e)) => {
                let inv = symox::ops::pow(&e, &Expr::minus_one());
                Ok(wrapm(self.m.scale(&inv)).into_pyobject(py)?.into_any().unbind())
            }
            _ => Ok(py.NotImplemented()),
        }
    }

    fn __pow__(&self, py: Python<'_>, o: &Bound<'_, PyAny>, _m: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        let e = to_expr(o)?;
        let n = match e.as_number() {
            Some(n) if n.is_integer() => n.to_i64().ok_or_else(|| PyValueError::new_err("exponent too large"))?,
            Some(n) if n.is_float() && n.to_f64().fract() == 0.0 => n.to_f64() as i64,
            _ => {
                return Err(PyNotImplementedError::new_err(
                    "only integer powers of matrices are implemented",
                ))
            }
        };
        Ok(wrapm(self.m.powi(n).py(py)?))
    }

    fn multiply_elementwise(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.multiply_elementwise(&to_matrix(o)?).py(py)?))
    }

    fn multiply(&self, py: Python<'_>, o: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.matmul(&to_matrix(o)?).py(py)?))
    }

    fn __richcmp__(&self, py: Python<'_>, other: &Bound<'_, PyAny>, op: pyo3::basic::CompareOp) -> PyResult<PyObject> {
        use pyo3::basic::CompareOp::*;
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

    #[classattr]
    const __hash__: Option<PyObject> = None;

    // ------------------------------------------------------------------
    // element-wise
    // ------------------------------------------------------------------

    fn applyfunc(&self, f: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.try_applyfunc(|e| to_expr(&f.call1((wrap(e.clone()),))?))?))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn subs(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        let pairs = subs_pairs(args)?;
        Ok(wrapm(self.m.subs(&pairs)))
    }

    fn xreplace(&self, rule: &Bound<'_, PyDict>) -> PyResult<PyMatrix> {
        let mut m = FxHashMap::default();
        for (k, v) in rule.iter() {
            m.insert(to_expr(&k)?, to_expr(&v)?);
        }
        Ok(wrapm(self.m.xreplace(&m)))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn diff(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        let vars = diff_vars(args)?;
        Ok(wrapm(self.m.diff_many(&vars)))
    }

    #[pyo3(signature = (*_args, **_kwargs))]
    fn expand(&self, _args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
        wrapm(self.m.expand())
    }

    #[pyo3(signature = (n=15, **_kwargs))]
    fn evalf(&self, n: u32, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
        let _ = n;
        wrapm(self.m.evalf())
    }

    #[pyo3(signature = (n=15, **_kwargs))]
    fn n(&self, n: u32, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
        let _ = n;
        wrapm(self.m.evalf())
    }

    /// Simplify the entries in place (like SymPy's mutable matrices).
    #[pyo3(signature = (**_kwargs))]
    fn simplify(&mut self, _kwargs: Option<&Bound<'_, PyDict>>) {
        self.m = self.m.simplify();
    }

    fn cancel(&self) -> PyMatrix {
        wrapm(self.m.cancel())
    }

    // ------------------------------------------------------------------
    // linear algebra
    // ------------------------------------------------------------------

    #[pyo3(signature = (method="bareiss", iszerofunc=None))]
    fn det(&self, py: Python<'_>, method: &str, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        let _ = iszerofunc;
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
        let _ = (iszerofunc, try_block_diag);
        let r = match method.map(|s| s.to_uppercase()) {
            None => self.m.inv(),
            Some(s) if s == "GE" => self.m.inverse_ge(),
            Some(s) if s == "LU" => self.m.inverse_lu(),
            Some(s) if s == "ADJ" => self.m.inverse_adj(),
            Some(s) if s == "DM" || s == "DMNC" => self.m.inv(),
            Some(s) => return Err(PyValueError::new_err(format!("Inversion method unrecognized: {s}"))),
        };
        Ok(wrapm(r.py(py)?))
    }

    fn inverse_GE(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.inverse_ge().py(py)?))
    }

    fn inverse_LU(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.inverse_lu().py(py)?))
    }

    fn inverse_ADJ(&self, py: Python<'_>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.inverse_adj().py(py)?))
    }

    #[pyo3(signature = (method="berkowitz"))]
    fn adjugate(&self, py: Python<'_>, method: &str) -> PyResult<PyMatrix> {
        let _ = method;
        Ok(wrapm(self.m.adjugate().py(py)?))
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
        Ok(wrapm(self.m.cofactor_matrix().py(py)?))
    }

    fn minor_submatrix(&self, i: isize, j: isize) -> PyResult<PyMatrix> {
        let (i, j) = (norm_index(i, self.m.rows())?, norm_index(j, self.m.cols())?);
        Ok(wrapm(self.m.minor_submatrix(i, j)))
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
        let _ = (iszerofunc, simplify, normalize_last);
        let (r, p) = self.m.rref();
        if pivots {
            let t = PyTuple::new(py, p)?;
            Ok(PyTuple::new(py, [wrapm(r).into_pyobject(py)?.into_any(), t.into_any()])?.into_any().unbind())
        } else {
            Ok(wrapm(r).into_pyobject(py)?.into_any().unbind())
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
        let _ = (iszerofunc, simplify);
        let (r, p) = self.m.echelon_form();
        if with_pivots {
            let t = PyTuple::new(py, p)?;
            Ok(PyTuple::new(py, [wrapm(r).into_pyobject(py)?.into_any(), t.into_any()])?.into_any().unbind())
        } else {
            Ok(wrapm(r).into_pyobject(py)?.into_any().unbind())
        }
    }

    #[pyo3(signature = (iszerofunc=None, simplify=false))]
    fn rank(&self, iszerofunc: Option<&Bound<'_, PyAny>>, simplify: bool) -> usize {
        let _ = (iszerofunc, simplify);
        self.m.rank()
    }

    #[pyo3(signature = (simplify=false, iszerofunc=None))]
    fn nullspace(&self, simplify: bool, iszerofunc: Option<&Bound<'_, PyAny>>) -> Vec<PyMatrix> {
        let _ = (iszerofunc, simplify);
        self.m.nullspace().into_iter().map(wrapm).collect()
    }

    #[pyo3(signature = (simplify=false))]
    fn columnspace(&self, simplify: bool) -> Vec<PyMatrix> {
        let _ = simplify;
        self.m.columnspace().into_iter().map(wrapm).collect()
    }

    #[pyo3(signature = (simplify=false))]
    fn rowspace(&self, simplify: bool) -> Vec<PyMatrix> {
        let _ = simplify;
        self.m.rowspace().into_iter().map(wrapm).collect()
    }

    #[pyo3(signature = (iszerofunc=None, simpfunc=None, rankcheck=false))]
    fn LUdecomposition<'py>(
        &self,
        py: Python<'py>,
        iszerofunc: Option<&Bound<'_, PyAny>>,
        simpfunc: Option<&Bound<'_, PyAny>>,
        rankcheck: bool,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let _ = (iszerofunc, simpfunc, rankcheck);
        let (l, u, p) = self.m.lu_decomposition().py(py)?;
        let perm = PyList::new(py, p.iter().map(|(a, b)| vec![*a, *b]))?;
        PyTuple::new(
            py,
            [
                wrapm(l).into_pyobject(py)?.into_any(),
                wrapm(u).into_pyobject(py)?.into_any(),
                perm.into_any(),
            ],
        )
    }

    #[pyo3(signature = (rhs, iszerofunc=None))]
    fn LUsolve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>, iszerofunc: Option<&Bound<'_, PyAny>>) -> PyResult<PyMatrix> {
        let _ = iszerofunc;
        Ok(wrapm(self.m.lu_solve(&to_matrix(rhs)?).py(py)?))
    }

    #[pyo3(signature = (rhs, method="GJ"))]
    fn solve(&self, py: Python<'_>, rhs: &Bound<'_, PyAny>, method: &str) -> PyResult<PyMatrix> {
        let b = to_matrix(rhs)?;
        let r = match method.to_uppercase().as_str() {
            "GJ" | "GE" => self.m.solve(&b),
            "LU" => self.m.lu_solve(&b),
            _ => self.m.inv().and_then(|i| i.matmul(&b)),
        };
        Ok(wrapm(r.py(py)?))
    }

    #[pyo3(signature = (B, freevar=false))]
    fn gauss_jordan_solve<'py>(&self, py: Python<'py>, B: &Bound<'_, PyAny>, freevar: bool) -> PyResult<Bound<'py, PyTuple>> {
        let _ = freevar;
        let (s, t) = self.m.gauss_jordan_solve(&to_matrix(B)?).py(py)?;
        PyTuple::new(py, [wrapm(s), wrapm(t)])
    }

    /// Characteristic polynomial ``det(lambda*I - M)`` as a ``PurePoly``.
    #[pyo3(signature = (x=None, simplify=None))]
    fn charpoly<'py>(
        &self,
        py: Python<'py>,
        x: Option<&Bound<'_, PyAny>>,
        simplify: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let _ = simplify;
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

    #[pyo3(signature = (error_when_incomplete=true, **_flags))]
    fn eigenvals<'py>(&self, py: Python<'py>, error_when_incomplete: bool, _flags: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyDict>> {
        let _ = error_when_incomplete;
        let ev = self.m.eigenvals().py(py)?;
        let d = PyDict::new(py);
        for (v, k) in ev {
            d.set_item(wrap(v), k)?;
        }
        Ok(d)
    }

    #[pyo3(signature = (error_when_incomplete=true, **_flags))]
    fn eigenvects<'py>(&self, py: Python<'py>, error_when_incomplete: bool, _flags: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyList>> {
        let _ = error_when_incomplete;
        let ev = self.m.eigenvects().py(py)?;
        let items: Vec<Bound<'py, PyTuple>> = ev
            .into_iter()
            .map(|(v, k, vs)| {
                let l = PyList::new(py, vs.into_iter().map(wrapm))?;
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
        Ok(wrapm(self.m.jacobian(vars.data()).py(py)?))
    }

    #[pyo3(signature = (ord=None))]
    fn norm(&self, py: Python<'_>, ord: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        let o: Option<String> = match ord {
            None => None,
            Some(o) if o.is_none() => None,
            Some(o) => Some(if let Ok(s) = o.extract::<String>() { s } else { to_expr(o)?.to_string() }),
        };
        Ok(wrap(self.m.norm_ord(o.as_deref()).py(py)?))
    }

    #[pyo3(signature = (hermitian=true))]
    fn cholesky(&self, py: Python<'_>, hermitian: bool) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.cholesky(hermitian).py(py)?))
    }

    #[pyo3(signature = (hermitian=true))]
    fn LDLdecomposition<'py>(&self, py: Python<'py>, hermitian: bool) -> PyResult<Bound<'py, PyTuple>> {
        let (l, d) = self.m.ldl_decomposition(hermitian).py(py)?;
        PyTuple::new(py, [wrapm(l), wrapm(d)])
    }

    fn QRdecomposition<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (q, r) = self.m.qr_decomposition().py(py)?;
        PyTuple::new(py, [wrapm(q), wrapm(r)])
    }

    #[pyo3(signature = (reals_only=false, sort=false, normalize=false))]
    fn diagonalize<'py>(&self, py: Python<'py>, reals_only: bool, sort: bool, normalize: bool) -> PyResult<Bound<'py, PyTuple>> {
        let _ = (reals_only, sort, normalize);
        let (p, d) = self.m.diagonalize().py(py)?;
        PyTuple::new(py, [wrapm(p), wrapm(d)])
    }

    #[pyo3(signature = (reals_only=false, **_kwargs))]
    fn is_diagonalizable(&self, reals_only: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> bool {
        let _ = reals_only;
        self.m.is_diagonalizable()
    }

    #[getter]
    fn is_hermitian(&self) -> Option<bool> {
        self.m.is_hermitian()
    }

    fn vec(&self) -> PyMatrix {
        wrapm(self.m.vec())
    }

    fn values(&self) -> Vec<PyExpr> {
        self.m.values().into_iter().map(wrap).collect()
    }

    fn flat(&self) -> Vec<PyExpr> {
        self.m.data().iter().cloned().map(wrap).collect()
    }

    fn as_mutable(&self) -> PyMatrix {
        self.clone()
    }

    fn as_immutable(&self) -> PyMatrix {
        self.clone()
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, (usize, usize, Vec<String>))> {
        let cls = py.import("symox")?.getattr("Matrix")?;
        let data: Vec<String> = self.m.data().iter().map(printer::srepr).collect();
        Ok((cls, (self.m.rows(), self.m.cols(), data)))
    }

    #[staticmethod]
    #[pyo3(name = "eye", signature = (n, m=None, **_kwargs))]
    fn eye_s(n: usize, m: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
        eye(n, m, None)
    }

    #[staticmethod]
    #[pyo3(name = "zeros", signature = (r, c=None, **_kwargs))]
    fn zeros_s(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
        zeros(r, c, None)
    }

    #[staticmethod]
    #[pyo3(name = "ones", signature = (r, c=None, **_kwargs))]
    fn ones_s(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
        ones(r, c, None)
    }

    #[staticmethod]
    #[pyo3(name = "diag", signature = (*values, **_kwargs))]
    fn diag_s(values: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
        diag(values, None)
    }

    #[pyo3(signature = (b, hermitian=None, conjugate_convention=None))]
    fn dot(
        &self,
        py: Python<'_>,
        b: &Bound<'_, PyAny>,
        hermitian: Option<bool>,
        conjugate_convention: Option<&str>,
    ) -> PyResult<PyExpr> {
        let _ = (hermitian, conjugate_convention);
        Ok(wrap(self.m.dot(&to_matrix(b)?).py(py)?))
    }

    fn cross(&self, py: Python<'_>, b: &Bound<'_, PyAny>) -> PyResult<PyMatrix> {
        Ok(wrapm(self.m.cross(&to_matrix(b)?).py(py)?))
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
fn eye(n: usize, m: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
    let (r, c) = dims(n, m);
    wrapm(Matrix::from_fn(r, c, |i, j| if i == j { Expr::one() } else { Expr::zero() }))
}

#[pyfunction]
#[pyo3(signature = (r, c=None, **_kwargs))]
fn zeros(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
    let (r, c) = dims(r, c);
    wrapm(Matrix::zeros(r, c))
}

#[pyfunction]
#[pyo3(signature = (r, c=None, **_kwargs))]
fn ones(r: usize, c: Option<usize>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyMatrix {
    let (r, c) = dims(r, c);
    wrapm(Matrix::ones(r, c))
}

#[pyfunction]
#[pyo3(signature = (*values, **_kwargs))]
fn diag(values: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyMatrix> {
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
    Ok(())
}
