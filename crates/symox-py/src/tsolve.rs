//! Python bindings for the native `_tsolve` steps
//! (`symox::solve::tsolve`) and the `sympy.solvers.bivariate` helpers.
//! The recursive `solve` calls, `nsimplify`, `equals` and the inverse of
//! host function classes call back into Python; everything else is native.

use pyo3::exceptions::{PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use symox::solve::tsolve::{self as ts, Fail, Flags, FuncClass, Host, Options};
use symox::Expr;

use crate::{to_expr, wrap};

struct PyHost<'py> {
    py: Python<'py>,
    solve: Bound<'py, PyAny>,
    flags: Bound<'py, PyDict>,
}

type F = Fail<PyErr>;

fn host_err(py: Python<'_>, e: PyErr) -> F {
    if e.is_instance_of::<PyNotImplementedError>(py) {
        Fail::NotImplemented(Some(e), String::new())
    } else if e.is_instance_of::<PyValueError>(py) {
        Fail::Value(Some(e), String::new())
    } else {
        Fail::Host(e)
    }
}

fn to_py_err(e: F) -> PyErr {
    match e {
        Fail::NotImplemented(Some(e), _) | Fail::Value(Some(e), _) | Fail::Host(e) => e,
        Fail::NotImplemented(None, m) => PyNotImplementedError::new_err(m),
        Fail::Value(None, m) => PyValueError::new_err(m),
        Fail::Type(m) => PyTypeError::new_err(m),
        Fail::Poly(e) => crate::poly::err(e),
    }
}

impl<'py> PyHost<'py> {
    fn new(py: Python<'py>, solve: Option<Bound<'py, PyAny>>, flags: Option<Bound<'py, PyDict>>) -> PyResult<Self> {
        let solve = match solve {
            Some(s) => s,
            None => py.import("symox")?.getattr("solve")?,
        };
        let flags = flags.unwrap_or_else(|| PyDict::new(py));
        Ok(PyHost { py, solve, flags })
    }
    fn exprs(&self, v: &Bound<'py, PyAny>) -> Result<Vec<Expr>, F> {
        let mut out = Vec::new();
        for item in v.try_iter().map_err(Fail::Host)? {
            out.push(to_expr(&item.map_err(Fail::Host)?).map_err(Fail::Host)?);
        }
        Ok(out)
    }
}

impl Host for PyHost<'_> {
    type Error = PyErr;
    fn solve(&mut self, eq: &Expr, x: &Expr, flags: Flags) -> Result<Vec<Expr>, F> {
        let py = self.py;
        let kw = match flags {
            Flags::Plain => None,
            Flags::Same => Some(self.flags.copy().map_err(Fail::Host)?),
            Flags::NoBivariate => {
                let d = self.flags.copy().map_err(Fail::Host)?;
                d.set_item("bivariate", false).map_err(Fail::Host)?;
                Some(d)
            }
        };
        let r = self.solve.call((wrap(eq.clone()), wrap(x.clone())), kw.as_ref()).map_err(|e| host_err(py, e))?;
        self.exprs(&r)
    }
    fn nsimplify(&mut self, e: &Expr) -> Result<Expr, F> {
        let py = self.py;
        let f = py.import("symox").and_then(|m| m.getattr("nsimplify")).map_err(Fail::Host)?;
        let r = f.call1((wrap(e.clone()),)).map_err(|e| host_err(py, e))?;
        to_expr(&r).map_err(Fail::Host)
    }
    fn equals_zero(&mut self, e: &Expr) -> Result<Option<bool>, F> {
        let py = self.py;
        let obj = wrap(e.clone()).into_pyobject(py).map_err(Fail::Host)?;
        let r = obj.call_method1("equals", (0,)).map_err(|e| host_err(py, e))?;
        if r.is_none() {
            return Ok(None);
        }
        Ok(Some(r.is_truthy().map_err(Fail::Host)?))
    }
    fn count_ops(&mut self, e: &Expr) -> usize {
        use num_traits::ToPrimitive;
        crate::operation_count(e).to_usize().unwrap_or(usize::MAX)
    }
    fn inverse(&mut self, f: &Expr, rhs: &Expr) -> Result<Option<Expr>, F> {
        let py = self.py;
        let helper = py
            .import("symox.solvers._tsolve_steps")
            .and_then(|m| m.getattr("_inverse"))
            .map_err(Fail::Host)?;
        let inv = helper.call1((wrap(f.clone()),)).map_err(|e| host_err(py, e))?;
        if inv.is_none() {
            return Ok(None);
        }
        let r = inv.call1((wrap(rhs.clone()),)).map_err(|e| host_err(py, e))?;
        to_expr(&r).map(Some).map_err(Fail::Host)
    }
}

fn flag(flags: &Bound<'_, PyDict>, key: &str, default: bool) -> PyResult<bool> {
    match flags.get_item(key)? {
        Some(v) => v.is_truthy(),
        None => Ok(default),
    }
}

fn exprs_out(py: Python<'_>, v: Vec<Expr>) -> PyResult<PyObject> {
    Ok(PyList::new(py, v.into_iter().map(wrap))?.into_any().unbind())
}

/// The `_tsolve` steps for `f = 0` in the Symbol `symbol`: a list of
/// solutions, or None when no step applies.
#[pyfunction]
fn _tsolve_step<'py>(
    py: Python<'py>, solve: Bound<'py, PyAny>, f: &Bound<'py, PyAny>, symbol: &Bound<'py, PyAny>, flags: Bound<'py, PyDict>,
) -> PyResult<Option<PyObject>> {
    let (f, sym) = (to_expr(f)?, to_expr(symbol)?);
    let opts = Options {
        force: flag(&flags, "force", true)?,
        bivariate: flag(&flags, "bivariate", true)?,
        check: flag(&flags, "check", true)?,
        numerical: flag(&flags, "numerical", true)?,
    };
    let mut host = PyHost::new(py, Some(solve), Some(flags))?;
    match ts::tsolve_step(&mut host, &f, &sym, opts) {
        Ok(Some(v)) => exprs_out(py, v).map(Some),
        Ok(None) => Ok(None),
        // the reference `_solve` catches PolynomialError around `_tsolve`
        Err(Fail::Poly(_)) => Ok(None),
        Err(e) => Err(to_py_err(e)),
    }
}

#[pyfunction]
fn _tsolve_invert(py: Python<'_>, eq: &Bound<'_, PyAny>, sym: &Bound<'_, PyAny>) -> PyResult<(crate::PyExpr, crate::PyExpr)> {
    let mut host = PyHost::new(py, None, None)?;
    let (r, l) = ts::invert(&mut host, &to_expr(eq)?, &to_expr(sym)?).map_err(to_py_err)?;
    Ok((wrap(r), wrap(l)))
}

#[pyfunction]
#[pyo3(signature = (expr, recursive=false))]
fn _bivariate_mexpand(expr: &Bound<'_, PyAny>, recursive: bool) -> PyResult<crate::PyExpr> {
    Ok(wrap(ts::mexpand(&to_expr(expr)?, recursive)))
}

#[pyfunction]
fn _bivariate_filtered_gens(gens: Vec<Bound<'_, PyAny>>, symbol: &Bound<'_, PyAny>) -> PyResult<Vec<crate::PyExpr>> {
    let gens = gens.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    Ok(ts::filtered_gens(&gens, &to_expr(symbol)?).into_iter().map(wrap).collect())
}

#[pyfunction]
#[pyo3(signature = (lhs, func, x=None))]
fn _bivariate_mostfunc(lhs: &Bound<'_, PyAny>, func: &str, x: Option<&Bound<'_, PyAny>>) -> PyResult<Option<crate::PyExpr>> {
    let class = match func {
        "log" => FuncClass::Log,
        "exp" => FuncClass::Exp,
        "Pow" => FuncClass::Pow,
        _ => return Err(PyValueError::new_err(format!("unsupported class {func}"))),
    };
    let x = x.map(to_expr).transpose()?;
    Ok(ts::mostfunc(&to_expr(lhs)?, class, x.as_ref()).map(wrap))
}

#[pyfunction]
fn _bivariate_linab(arg: &Bound<'_, PyAny>, symbol: &Bound<'_, PyAny>) -> PyResult<(crate::PyExpr, crate::PyExpr, crate::PyExpr)> {
    let (a, b, x) = ts::linab(&to_expr(arg)?, &to_expr(symbol)?);
    Ok((wrap(a), wrap(b), wrap(x)))
}

#[pyfunction]
fn _bivariate_lambert(py: Python<'_>, eq: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let mut host = PyHost::new(py, None, None)?;
    let v = ts::lambert(&mut host, &to_expr(eq)?, &to_expr(x)?).map_err(to_py_err)?;
    exprs_out(py, v)
}

#[pyfunction]
fn _bivariate_solve_lambert(py: Python<'_>, f: &Bound<'_, PyAny>, symbol: &Bound<'_, PyAny>, gens: Vec<Bound<'_, PyAny>>) -> PyResult<PyObject> {
    let gens = gens.iter().map(to_expr).collect::<PyResult<Vec<_>>>()?;
    let mut host = PyHost::new(py, None, None)?;
    let v = ts::solve_lambert(&mut host, &to_expr(f)?, &to_expr(symbol)?, &gens).map_err(to_py_err)?;
    exprs_out(py, v)
}

#[pyfunction]
fn _bivariate_type(f: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, y: &Bound<'_, PyAny>) -> PyResult<Option<(crate::PyExpr, crate::PyExpr, crate::PyExpr)>> {
    let r = ts::bivariate_type(&to_expr(f)?, &to_expr(x)?, &to_expr(y)?).map_err(crate::poly::err)?;
    Ok(r.map(|(g, p, u)| (wrap(g), wrap(p), wrap(u))))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_tsolve_step, m)?)?;
    m.add_function(wrap_pyfunction!(_tsolve_invert, m)?)?;
    m.add_function(wrap_pyfunction!(_bivariate_mexpand, m)?)?;
    m.add_function(wrap_pyfunction!(_bivariate_filtered_gens, m)?)?;
    m.add_function(wrap_pyfunction!(_bivariate_mostfunc, m)?)?;
    m.add_function(wrap_pyfunction!(_bivariate_linab, m)?)?;
    m.add_function(wrap_pyfunction!(_bivariate_lambert, m)?)?;
    m.add_function(wrap_pyfunction!(_bivariate_solve_lambert, m)?)?;
    m.add_function(wrap_pyfunction!(_bivariate_type, m)?)?;
    Ok(())
}
