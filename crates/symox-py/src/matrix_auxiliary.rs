//! Bindings for the auxiliary `MatrixBase` kernels in
//! `symox::matrix::auxiliary` (Berkowitz vectors, fraction-free LU,
//! Householder reductions, Cramer's rule, `equals` and constructors).
use super::{err, to_matrix, PyMatrix};
use crate::{to_expr, wrap};
use pyo3::exceptions::{PyValueError,PyTypeError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyList, PyTuple};
use std::cell::RefCell;
use symox::expr::Expr;
use symox::matrix::auxiliary::{self, Equality};
use symox::matrix::{DetMethod, MResult, Matrix, MatrixError};

/// Boolean scalars, including NumPy/Tensor scalar containers via item(),
/// without importing their libraries or treating arbitrary integers as bools.
pub fn boolean_scalar(value:&Bound<'_,PyAny>)->PyResult<Option<bool>> {
    if value.is_instance_of::<PyBool>() {return Ok(Some(value.extract()?));}
    if let Ok(expression)=value.extract::<PyRef<'_,crate::PyExpr>>() {
        if let symox::expr::Kind::Bool(value)=expression.e.kind() {return Ok(Some(*value));}
    }
    if value.hasattr("item")? {
        let scalar=value.call_method0("item")?;
        if scalar.is_instance_of::<PyBool>() {return Ok(Some(scalar.extract()?));}
    }
    Ok(None)
}

/// Runs a native kernel whose callbacks may raise Python exceptions; the
/// first Python error wins over the placeholder matrix error it caused.
struct Callbacks {
    error: RefCell<Option<PyErr>>,
}

impl Callbacks {
    fn new() -> Self {
        Callbacks {
            error: RefCell::new(None),
        }
    }
    fn trap<T>(&self, value: PyResult<T>) -> MResult<T> {
        value.map_err(|e| {
            self.error.borrow_mut().get_or_insert(e);
            MatrixError::Other("python callback failed".into())
        })
    }
    fn zero(&self,value:&Expr,callback:Option<&Bound<'_,PyAny>>)->MResult<Option<bool>> {
        match callback {
            None=>Ok(symox::matrix::ratfunc::is_zero_assume(value)),
            Some(callback)=>self.trap(callback.call1((wrap(value.clone()),)).and_then(|result|{
                if result.is_none() {return Ok(None);}
                if let Some(value)=boolean_scalar(&result)? {return Ok(Some(value));}
                match to_expr(&result)?.kind() {
                    symox::expr::Kind::Bool(value)=>Ok(Some(*value)),
                    _=>Err(PyTypeError::new_err("iszerofunc must return True, False or None")),
                }
            })),
        }
    }
    fn transform(&self,value:&Expr,callback:Option<&Bound<'_,PyAny>>)->MResult<Expr> {
        match callback {
            None=>Ok(value.clone()),
            Some(callback)=>self.trap(callback.call1((wrap(value.clone()),)).and_then(|result|to_expr(&result))),
        }
    }
    fn finish<T>(self, py: Python<'_>, value: MResult<T>) -> PyResult<T> {
        match (value, self.error.into_inner()) {
            (Err(_), Some(e)) => Err(e),
            (value, _) => value.map_err(|e| err(py, e)),
        }
    }
}

pub(super) fn checked_callback<'a,'py>(callback:Option<&'a Bound<'py,PyAny>>,name:&str)->PyResult<Option<&'a Bound<'py,PyAny>>> {
    let callback=callback.filter(|value|!value.is_none());
    let valid=|value:&Bound<'py,PyAny>| {
        if let Ok(expression)=value.extract::<PyRef<'_,crate::PyExpr>>() {
            matches!(expression.e.kind(),symox::expr::Kind::Func(symox::expr::Func::Set(symox::sets::SetOp::Lambda|symox::sets::SetOp::IdentityFunction),_))
        } else {value.is_callable()}
    };
    if callback.is_some_and(|value|!valid(value)) {
        return Err(PyTypeError::new_err(format!("{name} must be callable or None")));
    }
    Ok(callback)
}

pub fn row_reduce(
    py:Python<'_>,m:&Matrix,iszerofunc:Option<&Bound<'_,PyAny>>,
    simplify:Option<&Bound<'_,PyAny>>,normalize_last:bool,normalize:bool,zero_above:bool,
)->PyResult<(Matrix,Vec<usize>,Vec<(usize,usize)>)> {
    let zero=checked_callback(iszerofunc,"iszerofunc")?;
    let callbacks=Callbacks::new();
    let iszero=|value:&Expr|->MResult<Option<bool>> {
        callbacks.zero(value,zero)
    };
    let simp=|value:&Expr|->MResult<Expr> {
        match simplify {
            None=>Ok(symox::matrix::linalg::pivot_simp(value)),
            Some(callback)=>callbacks.transform(value,Some(callback)),
        }
    };
    let result=symox::matrix::linalg::row_reduce_with(m,normalize_last,normalize,zero_above,&iszero,Some(&simp));
    callbacks.finish(py,result)
}

pub fn lu_simple(
    py:Python<'_>,m:&Matrix,iszerofunc:Option<&Bound<'_,PyAny>>,
    simplify:Option<&Bound<'_,PyAny>>,rankcheck:bool,
)->PyResult<(Matrix,Vec<(usize,usize)>)> {
    let zero=checked_callback(iszerofunc,"iszerofunc")?;
    let simplify=checked_callback(simplify,"simpfunc")?;
    let callbacks=Callbacks::new();
    let iszero=|value:&Expr|callbacks.zero(value,zero);
    let simp=|value:&Expr|callbacks.transform(value,simplify);
    let result=symox::matrix::linalg::lu_simple_with(m,rankcheck,&iszero,if simplify.is_some() {Some(&simp)} else {None});
    callbacks.finish(py,result)
}

pub fn lu_decomposition(
    py:Python<'_>,m:&Matrix,iszerofunc:Option<&Bound<'_,PyAny>>,
    simplify:Option<&Bound<'_,PyAny>>,rankcheck:bool,
)->PyResult<(Matrix,Matrix,Vec<(usize,usize)>)> {
    let zero=checked_callback(iszerofunc,"iszerofunc")?;
    let simplify=checked_callback(simplify,"simpfunc")?;
    let callbacks=Callbacks::new();
    let iszero=|value:&Expr|callbacks.zero(value,zero);
    let simp=|value:&Expr|callbacks.transform(value,simplify);
    let result=m.lu_decomposition_with(rankcheck,&iszero,if simplify.is_some() {Some(&simp)} else {None});
    callbacks.finish(py,result)
}

pub fn lu_solve(py:Python<'_>,m:&Matrix,rhs:&Matrix,iszerofunc:Option<&Bound<'_,PyAny>>)->PyResult<Matrix> {
    let zero=checked_callback(iszerofunc,"iszerofunc")?;
    let callbacks=Callbacks::new();
    let result=m.lu_solve_with(rhs,&|value|callbacks.zero(value,zero));
    callbacks.finish(py,result)
}

pub fn determinant(py:Python<'_>,m:&Matrix,method:DetMethod,iszerofunc:Option<&Bound<'_,PyAny>>)->PyResult<Expr> {
    let zero=checked_callback(iszerofunc,"iszerofunc")?;
    let callbacks=Callbacks::new();
    let result=if zero.is_some() {m.det_method_with(method,&|value|callbacks.zero(value,zero))} else {m.det_method(method)};
    callbacks.finish(py,result)
}

pub fn inverse(
    py:Python<'_>,m:&Matrix,method:Option<&str>,iszerofunc:Option<&Bound<'_,PyAny>>,try_block_diag:bool,
)->PyResult<Matrix> {
    let zero=checked_callback(iszerofunc,"iszerofunc")?;
    let method=method.map(str::to_uppercase);
    if let Some(method)=method.as_deref() {
        if !["GE","LU","ADJ","CH","LDL","QR","DM","DMNC"].contains(&method) {
            return Err(PyValueError::new_err(format!("Inversion method unrecognized: {method}")));
        }
    }
    let callbacks=Callbacks::new();
    fn run(m:&Matrix,method:Option<&str>,zero:Option<&dyn Fn(&Expr)->MResult<Option<bool>>>,blocks:bool)->MResult<Matrix> {
        if !m.is_square() {return Err(MatrixError::NonSquare("A Matrix must be square to invert.".into()));}
        if blocks {
            let parts=m.diagonal_blocks();
            if parts.len()>1 {
                let inverses=parts.iter().map(|part|run(part,method,zero,false)).collect::<MResult<Vec<_>>>()?;
                return Ok(Matrix::diag_blocks(&inverses));
            }
        }
        match (method,zero) {
            (None,Some(zero))|(Some("GE"),Some(zero))=>m.inverse_ge_with(zero),
            (Some("LU"),Some(zero))=>m.inverse_lu_with(zero),
            (Some("ADJ"),Some(zero))=>m.inverse_adj_with(zero),
            (None,None)|(Some("DM"|"DMNC"),_)=>m.inv(),
            (Some("GE"),None)=>m.inverse_ge(),
            (Some("LU"),None)=>m.inverse_lu(),
            (Some("ADJ"),None)=>m.inverse_adj(),
            (Some(method),zero)=>{
                if let Some(zero)=zero {m.verify_invertible_with(zero)?;}
                let identity=Matrix::eye(m.rows());
                match method {"CH"=>m.cholesky_solve(&identity),"LDL"=>m.ldl_solve(&identity),"QR"=>m.qr_solve(&identity),_=>unreachable!()}
            },
        }
    }
    let predicate=|value:&Expr|callbacks.zero(value,zero);
    let result=run(m,method.as_deref(),if zero.is_some() {Some(&predicate)} else {None},try_block_diag);
    callbacks.finish(py,result)
}

pub fn upper_flag(upper: &Bound<'_, PyAny>) -> PyResult<bool> {
    if !upper.is_instance_of::<PyBool>() {
        return Err(PyValueError::new_err("upper must be a boolean"));
    }
    upper.extract()
}

pub fn bidiagonal(py: Python<'_>, m: &Matrix, upper: bool) -> PyResult<(Matrix, Matrix, Matrix)> {
    if upper {
        m.householder_bidiagonal(true).map_err(|e| err(py, e))
    } else {
        let (u, b, v) = m
            .adjoint()
            .householder_bidiagonal(true)
            .map_err(|e| err(py, e))?;
        Ok((v.adjoint(), b.adjoint(), u.adjoint()))
    }
}

pub fn bidiagonalize(py: Python<'_>, m: &Matrix, upper: bool) -> PyResult<Matrix> {
    if upper {
        Ok(m.householder_bidiagonal(false).map_err(|e| err(py, e))?.1)
    } else {
        Ok(m.adjoint()
            .householder_bidiagonal(false)
            .map_err(|e| err(py, e))?
            .1
            .adjoint())
    }
}

/// `berkowitz()` as nested tuples; the leading `(1,)` is a Python int as in
/// SymPy.
pub fn berkowitz<'py>(py: Python<'py>, m: &Matrix) -> PyResult<Bound<'py, PyTuple>> {
    let vectors = m.berkowitz_vectors().map_err(|e| err(py, e))?;
    let mut items = Vec::with_capacity(vectors.len());
    for (k, v) in vectors.into_iter().enumerate() {
        if k == 0 {
            items.push(PyTuple::new(py, [1i64])?);
        } else {
            items.push(PyTuple::new(py, v.into_iter().map(wrap))?);
        }
    }
    PyTuple::new(py, items)
}

pub fn is_anti_symmetric(
    py: Python<'_>,
    m: &Matrix,
    simplify: &Bound<'_, PyAny>,
) -> PyResult<Option<bool>> {
    let callbacks = Callbacks::new();
    let custom = !simplify.is_instance_of::<PyBool>() && simplify.is_callable();
    let enabled = custom || simplify.is_truthy()?;
    let result = m.is_anti_symmetric_with(&mut |e| {
        if custom {
            callbacks.trap(simplify.call1((wrap(e.clone()),)).and_then(|r| to_expr(&r)))
        } else if enabled {
            Ok(symox::simplify::simplify(e))
        } else {
            Ok(e.clone())
        }
    });
    callbacks.finish(py, result)
}

pub fn normalized(
    py: Python<'_>,
    m: &Matrix,
    iszerofunc: Option<&Bound<'_, PyAny>>,
) -> PyResult<Matrix> {
    let callbacks = Callbacks::new();
    let result = m.normalized_with(&mut |norm| match iszerofunc {
        Some(f) if !f.is_none() => {
            callbacks.trap(f.call1((wrap(norm.clone()),)).and_then(|r| r.is_truthy()))
        }
        _ => Ok(symox::assumptions::is_zero(norm) == Some(true)),
    });
    callbacks.finish(py, result)
}

pub fn cramer_solve(
    py: Python<'_>,
    this: &PyMatrix,
    rhs: &Bound<'_, PyAny>,
    det_method: &Bound<'_, PyAny>,
) -> PyResult<Matrix> {
    let rhs = to_matrix(rhs)?;
    let callbacks = Callbacks::new();
    let named = det_method.extract::<String>().ok();
    let method = match &named {
        Some(name) => Some(DetMethod::from_name(name).map_err(|e| err(py, e))?),
        None => None,
    };
    let result = this.m.cramer_solve(&rhs, &mut |matrix| match method {
        Some(method) => matrix.det_method(method),
        None => callbacks.trap(
            det_method
                .call1((this.wrapped(matrix.clone()),))
                .and_then(|r| to_expr(&r)),
        ),
    });
    callbacks.finish(py, result)
}

pub fn equals(
    py: Python<'_>,
    m: &Matrix,
    other: &Bound<'_, PyAny>,
    failing_expression: bool,
) -> PyResult<PyObject> {
    let shape: Option<(usize, usize)> = other.getattr("shape").ok().and_then(|s| s.extract().ok());
    if shape != Some(m.shape()) {
        return Ok(false.into_pyobject(py)?.to_owned().into_any().unbind());
    }
    let other = match to_matrix(other) {
        Ok(o) => o,
        Err(_) => return Ok(false.into_pyobject(py)?.to_owned().into_any().unbind()),
    };
    Ok(match m.equals(&other) {
        Equality::True => true.into_pyobject(py)?.to_owned().into_any().unbind(),
        Equality::False => false.into_pyobject(py)?.to_owned().into_any().unbind(),
        Equality::Unknown(diff) if failing_expression => {
            wrap(diff).into_pyobject(py)?.into_any().unbind()
        }
        Equality::Unknown(_) => py.None(),
    })
}

pub fn companion(py: Python<'_>, coeffs: &Bound<'_, PyAny>) -> PyResult<Matrix> {
    let coeffs: Vec<Expr> = coeffs
        .try_iter()?
        .map(|c| to_expr(&c?))
        .collect::<PyResult<_>>()?;
    if coeffs.len() < 2 {
        return Err(err(
            py,
            MatrixError::Value("companion matrices need degree >= 1".into()),
        ));
    }
    Matrix::companion(&coeffs).map_err(|e| err(py, e))
}

pub fn orthogonalize<'py>(
    py: Python<'py>,
    vecs: &Bound<'py, PyAny>,
    normalize: bool,
    rankcheck: bool,
) -> PyResult<Vec<Matrix>> {
    let vecs: Vec<Matrix> = vecs
        .try_iter()?
        .map(|v| to_matrix(&v?))
        .collect::<PyResult<_>>()?;
    Matrix::orthogonalize(&vecs, normalize, rankcheck).map_err(|e| err(py, e))
}

pub fn irregular(py: Python<'_>, ntop: usize, matrices: &Bound<'_, PyList>) -> PyResult<Matrix> {
    let matrices: Vec<Matrix> = matrices
        .iter()
        .map(|m| to_matrix(&m))
        .collect::<PyResult<_>>()?;
    Matrix::irregular(ntop, &matrices).map_err(|e| err(py, e))
}

pub fn casoratian(
    py: Python<'_>,
    seqs: &Bound<'_, PyAny>,
    n: &Bound<'_, PyAny>,
    zero: bool,
) -> PyResult<Expr> {
    let seqs: Vec<Expr> = seqs
        .try_iter()?
        .map(|s| to_expr(&s?))
        .collect::<PyResult<_>>()?;
    let n = to_expr(n)?;
    auxiliary::casoratian(&seqs, &n, zero).map_err(|e| err(py, e))
}

/// Integer diagonal offsets are compared with arbitrary scalar boundaries.
/// Rounding in the appropriate direction also accepts SymPy's fractional k.
pub fn triangular(
    py: Python<'_>,
    m: &Matrix,
    k: Option<&Bound<'_, PyAny>>,
    upper: bool,
) -> PyResult<Matrix> {
    if m.data().is_empty() {
        return Ok(m.clone());
    }
    let offset = match k {
        Some(k) => to_expr(k)?,
        None => Expr::zero(),
    };
    let rounded = symox::functions::func(
        if upper {
            symox::expr::Func::Ceiling
        } else {
            symox::expr::Func::Floor
        },
        vec![offset],
    );
    let number = rounded.as_number().ok_or_else(|| {
        pyo3::exceptions::PyTypeError::new_err("cannot determine a diagonal offset")
    })?;
    let offset = number.to_i64().unwrap_or_else(|| {
        if number.is_negative() {
            i64::MIN
        } else {
            i64::MAX
        }
    });
    let _ = py;
    Ok(m.triangular_part(offset, upper))
}
