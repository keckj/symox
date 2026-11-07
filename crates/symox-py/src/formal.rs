//! Python bindings for the native formal power series
//! (`symox::series::formal`, a port of `sympy/series/formal.py`).
//!
//! The public operations the reference calls on intermediate expressions
//! (`simplify`, `limit`, `integrate`, `apart`, `factor`, `cancel`, `doit`,
//! `leadterm`, `roots`, `linsolve`, `Coeff`) are the host's public entry
//! points; their exceptions propagate unchanged. The reference exceptions
//! raised by the algorithms map to the same Python classes.

use pyo3::exceptions::{PyIndexError, PyKeyError, PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::expr::{Expr, Func, Kind};
use symox::series::formal::{self as fm, Error, Exc, Fps, Hooks, Options, Seq, SimpleDe};

use crate::{to_expr, wrap as w};

struct Host<'py> {
    py: Python<'py>,
    symox: Bound<'py, PyModule>,
}

impl<'py> Host<'py> {
    fn new(py: Python<'py>) -> PyResult<Self> {
        Ok(Host { py, symox: py.import("symox")? })
    }

    fn obj(&self, e: &Expr) -> PyResult<Bound<'py, PyAny>> {
        Ok(w(e.clone()).into_pyobject(self.py)?.into_any())
    }

    fn method0(&self, e: &Expr, name: &str) -> PyResult<Expr> {
        to_expr(&self.obj(e)?.call_method0(name)?)
    }

    fn fun(&self, name: &str, args: &[&Expr]) -> PyResult<Expr> {
        let args: Vec<Bound<'py, PyAny>> = args.iter().map(|a| self.obj(a)).collect::<PyResult<_>>()?;
        to_expr(&self.symox.getattr(name)?.call1(PyTuple::new(self.py, args)?)?)
    }
}

impl Hooks for Host<'_> {
    type Error = PyErr;

    fn simplify(&mut self, e: &Expr) -> PyResult<Expr> {
        self.method0(e, "simplify")
    }

    fn limit0(&mut self, e: &Expr, x: &Expr) -> PyResult<Expr> {
        to_expr(&self.obj(e)?.call_method1("limit", (w(x.clone()), 0))?)
    }

    fn integrate(&mut self, e: &Expr, x: &Expr) -> PyResult<Expr> {
        self.fun("integrate", &[e, x])
    }

    fn apart(&mut self, e: &Expr, x: &Expr, full: bool) -> PyResult<Expr> {
        let kw = PyDict::new(self.py);
        kw.set_item("full", full)?;
        to_expr(&self.symox.getattr("apart")?.call((w(e.clone()), w(x.clone())), Some(&kw))?)
    }

    fn doit(&mut self, e: &Expr) -> PyResult<Expr> {
        self.method0(e, "doit")
    }

    fn factor(&mut self, e: &Expr) -> PyResult<Expr> {
        self.method0(e, "factor")
    }

    fn cancel(&mut self, e: &Expr) -> PyResult<Expr> {
        self.method0(e, "cancel")
    }

    fn leadterm_coeff(&mut self, e: &Expr, x: &Expr) -> PyResult<Expr> {
        to_expr(&self.obj(e)?.call_method1("leadterm", (w(x.clone()),))?.get_item(0)?)
    }

    fn roots(&mut self, e: &Expr, x: &Expr) -> PyResult<Vec<(Expr, Expr)>> {
        let found = self.symox.getattr("roots")?.call1((w(e.clone()), w(x.clone())))?;
        let found = found.downcast::<PyDict>()?;
        found.iter().map(|(k, v)| Ok((to_expr(&k)?, to_expr(&v)?))).collect()
    }

    fn linsolve(&mut self, eqs: &[Expr], syms: &[Expr]) -> PyResult<Vec<Expr>> {
        let eqs = PyList::new(self.py, eqs.iter().cloned().map(w))?;
        let syms = PyList::new(self.py, syms.iter().cloned().map(w))?;
        let sol = self.symox.getattr("linsolve")?.call1((eqs, syms))?;
        let mut out = Vec::new();
        for s in sol.try_iter()? {
            for i in s?.try_iter()? {
                out.push(to_expr(&i?)?);
            }
        }
        Ok(out)
    }

    fn coeff_node(&mut self, p: &Expr, x: &Expr, n: &Expr) -> PyResult<Expr> {
        let coeff = self.py.import("symox.series.formal")?.getattr("Coeff")?;
        to_expr(&coeff.call1((w(p.clone()), w(x.clone()), w(n.clone())))?)
    }

    fn tick(&mut self) -> PyResult<()> {
        self.py.check_signals()
    }
}

fn exc(e: Exc) -> PyErr {
    match e {
        Exc::Value(m) => PyValueError::new_err(m),
        Exc::NotImplemented(m) => PyNotImplementedError::new_err(m),
        Exc::Type(m) => PyTypeError::new_err(m),
        Exc::Index(m) => PyIndexError::new_err(m),
        Exc::Key(m) => PyKeyError::new_err(m),
    }
}

fn err(e: Error<PyErr>) -> PyErr {
    match e {
        Error::Exc(e) => exc(e),
        Error::Hook(e) => e,
    }
}

fn tuple<'py>(py: Python<'py>, items: Vec<Expr>) -> PyResult<PyObject> {
    Ok(PyTuple::new(py, items.into_iter().map(w))?.into_any().unbind())
}

fn opt3(py: Python<'_>, r: Option<(Expr, Expr, Expr)>) -> PyResult<PyObject> {
    match r {
        Some((a, b, c)) => tuple(py, vec![a, b, c]),
        None => Ok(py.None()),
    }
}

fn head(e: &Expr) -> PyResult<Func> {
    match e.kind() {
        Kind::Func(h, _) => Ok(h.clone()),
        _ => Err(PyTypeError::new_err(format!("{e} is not a function application"))),
    }
}

fn node(obj: &Bound<'_, PyAny>) -> PyResult<Fps> {
    Fps::from_node(&to_expr(obj)?).map_err(exc)
}

fn int(obj: &Bound<'_, PyAny>) -> PyResult<i64> {
    obj.extract::<i64>()
}

fn parts(py: Python<'_>, f: &Fps) -> PyResult<PyObject> {
    tuple(py, vec![f.function.clone(), f.parts_tuple().map_err(exc)?])
}

/// `rational_algorithm(f, x, k, order, full)`
#[pyfunction]
fn fps_rational_algorithm(py: Python<'_>, f: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>, order: i64, full: bool) -> PyResult<PyObject> {
    let mut h = Host::new(py)?;
    let r = fm::rational_algorithm(&mut h, &to_expr(f)?, &to_expr(x)?, &to_expr(k)?, order, full).map_err(err)?;
    opt3(py, r)
}

/// `rational_independent(terms, x)`
#[pyfunction]
fn fps_rational_independent(py: Python<'_>, terms: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let terms: Vec<Expr> = terms.try_iter()?.map(|t| to_expr(&t?)).collect::<PyResult<_>>()?;
    let mut h = Host::new(py)?;
    let r = fm::rational_independent(&mut h, &terms, &to_expr(x)?).map_err(err)?;
    Ok(PyList::new(py, r.into_iter().map(w))?.into_any().unbind())
}

/// The lazy `simpleDE(f, x, g, order)` generator; `g` is passed applied.
#[pyclass(unsendable)]
struct FpsSimpleDe {
    inner: SimpleDe,
}

#[pymethods]
impl FpsSimpleDe {
    #[new]
    fn new(f: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, gx: &Bound<'_, PyAny>, order: i64) -> PyResult<Self> {
        let g = head(&to_expr(gx)?)?;
        Ok(FpsSimpleDe { inner: SimpleDe::new(&to_expr(f)?, &to_expr(x)?, &g, order) })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PyObject>> {
        let mut h = Host::new(py)?;
        match self.inner.next(&mut h).map_err(err)? {
            Some((de, k)) => Ok(Some(PyTuple::new(py, [w(de).into_pyobject(py)?.into_any(), k.into_pyobject(py)?.into_any()])?.into_any().unbind())),
            None => Ok(None),
        }
    }
}

/// `exp_re(DE, r, k)`; `rk` is `r(k)`.
#[pyfunction]
fn fps_exp_re(de: &Bound<'_, PyAny>, rk: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>) -> PyResult<PyExprOut> {
    let r = head(&to_expr(rk)?)?;
    fm::exp_re(&to_expr(de)?, &r, &to_expr(k)?).map(w).map_err(exc)
}

/// `hyper_re(DE, r, k)`; `rk` is `r(k)`.
#[pyfunction]
fn fps_hyper_re(de: &Bound<'_, PyAny>, rk: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>) -> PyResult<PyExprOut> {
    let r = head(&to_expr(rk)?)?;
    fm::hyper_re(&to_expr(de)?, &r, &to_expr(k)?).map(w).map_err(exc)
}

type PyExprOut = crate::PyExpr;

/// `rsolve_hypergeometric(f, x, P, Q, k, m)`
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn fps_rsolve_hypergeometric(
    py: Python<'_>,
    f: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    p: &Bound<'_, PyAny>,
    q: &Bound<'_, PyAny>,
    k: &Bound<'_, PyAny>,
    m: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let mut h = Host::new(py)?;
    let r = fm::rsolve_hypergeometric(&mut h, &to_expr(f)?, &to_expr(x)?, &to_expr(p)?, &to_expr(q)?, &to_expr(k)?, &to_expr(m)?).map_err(err)?;
    opt3(py, r)
}

/// `solve_de(f, x, DE, order, g, k)`; `gx` is `g(x)`.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn fps_solve_de(
    py: Python<'_>,
    f: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    de: &Bound<'_, PyAny>,
    order: i64,
    gx: &Bound<'_, PyAny>,
    k: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let g = head(&to_expr(gx)?)?;
    let mut h = Host::new(py)?;
    let r = fm::solve_de(&mut h, &to_expr(f)?, &to_expr(x)?, &to_expr(de)?, order, &g, &to_expr(k)?).map_err(err)?;
    opt3(py, r)
}

/// `hyper_algorithm(f, x, k, order)`
#[pyfunction]
fn fps_hyper_algorithm(py: Python<'_>, f: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>, order: i64) -> PyResult<PyObject> {
    let mut h = Host::new(py)?;
    let r = fm::hyper_algorithm(&mut h, &to_expr(f)?, &to_expr(x)?, &to_expr(k)?, order).map_err(err)?;
    opt3(py, r)
}

fn parts_of(py: Python<'_>, r: Option<fm::Parts>) -> PyResult<PyObject> {
    match r {
        Some((ak, xk, ind)) => {
            let t = |s: &Seq| s.to_sum().map_err(exc);
            tuple(py, vec![t(&ak)?, t(&xk)?, ind])
        }
        None => Ok(py.None()),
    }
}

/// `compute_fps(f, x, x0, dir, hyper, order, rational, full)` with `dir`
/// already mapped to `1`/`-1` (or left as given, then rejected): `None` or
/// `(Sum(ak), Sum(xk), ind)`.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn fps_compute(
    py: Python<'_>,
    f: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    x0: &Bound<'_, PyAny>,
    dir: &Bound<'_, PyAny>,
    hyper: bool,
    order: i64,
    rational: bool,
    full: bool,
) -> PyResult<PyObject> {
    let mut h = Host::new(py)?;
    let o = Options { hyper, order, rational, full };
    let r = fm::compute_fps(&mut h, &to_expr(f)?, &to_expr(x)?, &to_expr(x0)?, &to_expr(dir)?, o).map_err(err)?;
    parts_of(py, r)
}

/// `fps._eval_term(pt)`
#[pyfunction]
fn fps_term(py: Python<'_>, s: &Bound<'_, PyAny>, pt: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    let f = node(s)?;
    let mut h = Host::new(py)?;
    f.term(&mut h, &to_expr(pt)?).map(w).map_err(err)
}

/// `fps.polynomial(n)`
#[pyfunction]
fn fps_polynomial(py: Python<'_>, s: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    let f = node(s)?;
    let mut h = Host::new(py)?;
    f.polynomial(&mut h, &to_expr(n)?).map(w).map_err(err)
}

/// `fps.truncate(n)` for a non-`None` `n`.
#[pyfunction]
fn fps_truncate(py: Python<'_>, s: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    let f = node(s)?;
    let mut h = Host::new(py)?;
    f.truncate(&mut h, &to_expr(n)?).map(w).map_err(err)
}

/// `fps.infinite`
#[pyfunction]
fn fps_infinite(s: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    node(s)?.infinite().map(w).map_err(exc)
}

/// `fps._eval_as_leading_term`
#[pyfunction]
fn fps_leading_term(py: Python<'_>, s: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    let f = node(s)?;
    let mut h = Host::new(py)?;
    f.leading_term(&mut h).map(w).map_err(err)
}

/// `fps.free_symbols`
#[pyfunction]
fn fps_free_symbols(py: Python<'_>, s: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    Ok(PyList::new(py, node(s)?.free_symbols().into_iter().map(w))?.into_any().unbind())
}

/// `self + other` for two series: `("fps", function, parts)` or
/// `("expr", f)`.
#[pyfunction]
fn fps_add_fps(py: Python<'_>, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    match node(a)?.add_fps(&node(b)?).map_err(exc)? {
        Ok(f) => parts(py, &f),
        Err(e) => tuple(py, vec![e]),
    }
}

/// `self + other` for `other` free of `x`: `(function, parts)`, or `None`
/// for `Add(self, other)`.
#[pyfunction]
fn fps_add_expr(py: Python<'_>, a: &Bound<'_, PyAny>, other: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    match node(a)?.add_expr(&to_expr(other)?) {
        Some(f) => parts(py, &f),
        None => Ok(py.None()),
    }
}

/// `-self`: `(function, parts)`.
#[pyfunction]
fn fps_neg(py: Python<'_>, a: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    parts(py, &node(a)?.neg())
}

/// `self*other`: `(function, parts)`, or `None` for `Mul(self, other)`.
#[pyfunction]
fn fps_mul_expr(py: Python<'_>, a: &Bound<'_, PyAny>, other: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    match node(a)?.mul_expr(&to_expr(other)?) {
        Some(f) => parts(py, &f),
        None => Ok(py.None()),
    }
}

/// `fps._eval_derivative(v)`: `(function, parts)`.
#[pyfunction]
fn fps_derivative(py: Python<'_>, a: &Bound<'_, PyAny>, v: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    parts(py, &node(a)?.derivative(&to_expr(v)?).map_err(exc)?)
}

/// `fps.integrate(v)`: `(function, parts)`.
#[pyfunction]
fn fps_integrate(py: Python<'_>, a: &Bound<'_, PyAny>, v: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let f = node(a)?;
    let mut h = Host::new(py)?;
    let r = f.integrate(&mut h, &to_expr(v)?).map_err(err)?;
    parts(py, &r)
}

/// The checks of `product` (`compose = False`) or `compose`.
#[pyfunction]
fn fps_check_pair(py: Python<'_>, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>, compose: bool) -> PyResult<()> {
    let (a, b) = (node(a)?, node(b)?);
    if compose {
        let mut h = Host::new(py)?;
        fm::check_compose(&mut h, &a, &b).map_err(err)
    } else {
        fm::check_pair(&a, &b).map_err(exc)
    }
}

/// The check of `inverse`.
#[pyfunction]
fn fps_check_inverse(py: Python<'_>, a: &Bound<'_, PyAny>) -> PyResult<()> {
    let a = node(a)?;
    let mut h = Host::new(py)?;
    fm::check_inverse(&mut h, &a).map_err(err)
}

/// `_eval_terms(n)` of the finite series: `kind` is `"product"`,
/// `"compose"` or `"inverse"` (then `b` is ignored).
#[pyfunction]
fn fps_finite_terms(py: Python<'_>, kind: &str, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    let a = node(a)?;
    let n = int(n)?;
    let mut h = Host::new(py)?;
    let r = match kind {
        "product" => fm::product_terms(&mut h, &a, &node(b)?, n),
        "compose" => fm::compose_terms(&mut h, &a, &node(b)?, n),
        "inverse" => fm::inverse_terms(&mut h, &a, n),
        _ => return Err(PyValueError::new_err(format!("unknown finite series {kind}"))),
    };
    r.map(w).map_err(err)
}

/// `FiniteFormalPowerSeries.truncate(n)` given `polynomial(n)`.
#[pyfunction]
fn fps_finite_truncate(a: &Bound<'_, PyAny>, terms: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    fm::finite_truncate(&node(a)?, &to_expr(terms)?, &to_expr(n)?).map(w).map_err(exc)
}

/// `Coeff.eval(p, x, n)`
#[pyfunction]
fn fps_coeff_eval(py: Python<'_>, p: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    match fm::coeff_eval(&to_expr(p)?, &to_expr(x)?, &to_expr(n)?).map_err(exc)? {
        Some(e) => Ok(w(e).into_pyobject(py)?.into_any().unbind()),
        None => Ok(py.None()),
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<FpsSimpleDe>()?;
    m.add_function(wrap_pyfunction!(fps_rational_algorithm, m)?)?;
    m.add_function(wrap_pyfunction!(fps_rational_independent, m)?)?;
    m.add_function(wrap_pyfunction!(fps_exp_re, m)?)?;
    m.add_function(wrap_pyfunction!(fps_hyper_re, m)?)?;
    m.add_function(wrap_pyfunction!(fps_rsolve_hypergeometric, m)?)?;
    m.add_function(wrap_pyfunction!(fps_solve_de, m)?)?;
    m.add_function(wrap_pyfunction!(fps_hyper_algorithm, m)?)?;
    m.add_function(wrap_pyfunction!(fps_compute, m)?)?;
    m.add_function(wrap_pyfunction!(fps_term, m)?)?;
    m.add_function(wrap_pyfunction!(fps_polynomial, m)?)?;
    m.add_function(wrap_pyfunction!(fps_truncate, m)?)?;
    m.add_function(wrap_pyfunction!(fps_infinite, m)?)?;
    m.add_function(wrap_pyfunction!(fps_leading_term, m)?)?;
    m.add_function(wrap_pyfunction!(fps_free_symbols, m)?)?;
    m.add_function(wrap_pyfunction!(fps_add_fps, m)?)?;
    m.add_function(wrap_pyfunction!(fps_add_expr, m)?)?;
    m.add_function(wrap_pyfunction!(fps_neg, m)?)?;
    m.add_function(wrap_pyfunction!(fps_mul_expr, m)?)?;
    m.add_function(wrap_pyfunction!(fps_derivative, m)?)?;
    m.add_function(wrap_pyfunction!(fps_integrate, m)?)?;
    m.add_function(wrap_pyfunction!(fps_check_pair, m)?)?;
    m.add_function(wrap_pyfunction!(fps_check_inverse, m)?)?;
    m.add_function(wrap_pyfunction!(fps_finite_terms, m)?)?;
    m.add_function(wrap_pyfunction!(fps_finite_truncate, m)?)?;
    m.add_function(wrap_pyfunction!(fps_coeff_eval, m)?)?;
    Ok(())
}
