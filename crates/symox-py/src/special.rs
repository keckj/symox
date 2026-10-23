//! Python bindings for the special functions (`symox::special`).

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict, PyTuple};

use symox::special::{self as sp, SpecialFn};
use symox::Expr;

use crate::{to_expr, wrap, PyExpr};

fn nonneg(o: &Bound<'_, PyAny>, what: &str) -> PyResult<u64> {
    let e = to_expr(o)?;
    e.as_number()
        .and_then(|n| n.to_bigint())
        .and_then(|b| u64::try_from(b).ok())
        .ok_or_else(|| PyValueError::new_err(format!("{what} must be a non-negative integer, got {e}")))
}

/// `stirling(n, k, d=None, kind=2, signed=False)`.
#[pyfunction]
#[pyo3(signature = (n, k, d=None, kind=2, signed=false))]
fn stirling(n: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>, d: Option<&Bound<'_, PyAny>>, kind: u32, signed: bool) -> PyResult<PyExpr> {
    let n = nonneg(n, "n")?;
    let k = nonneg(k, "k")?;
    let d = match d {
        Some(o) if !o.is_none() => Some(nonneg(o, "d")?),
        _ => None,
    };
    if kind != 1 && kind != 2 {
        return Err(PyValueError::new_err("kind must be 1 or 2"));
    }
    Ok(wrap(Expr::num(symox::Number::from_bigint(sp::stirling(n, k, d, kind, signed)))))
}

/// `multinomial_coefficients(m, n)` as a dict `{(k1, ..., km): coeff}`.
#[pyfunction]
fn multinomial_coefficients<'py>(py: Python<'py>, m: usize, n: u64) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    for (ks, c) in sp::multinomial_coefficients(m, n) {
        let t = PyTuple::new(py, ks)?;
        d.set_item(t, c)?;
    }
    Ok(d)
}

#[pyfunction]
#[pyo3(signature = (expr, deep=true))]
fn expand_func(expr: &Bound<'_, PyAny>, deep: bool) -> PyResult<PyExpr> {
    Ok(wrap(sp::expand_func_deep(&to_expr(expr)?, deep)))
}

/// Native engine behind `Expr.rewrite(*patterns, rule, deep=True)`.
#[pyfunction]
#[pyo3(signature = (e, patterns, rule, deep=true))]
fn _rewrite(e: &Bound<'_, PyAny>, patterns: Vec<String>, rule: &str, deep: bool) -> PyResult<PyExpr> {
    Ok(wrap(sp::rewrite::rewrite(&to_expr(e)?, &patterns, rule, deep)))
}

#[pyfunction]
fn besselsimp(expr: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(sp::besselsimp(&to_expr(expr)?)))
}

#[pyfunction]
#[pyo3(signature = (f, allow_hyper=false, rewrite=Some("default".to_string()), place=None))]
fn hyperexpand(
    f: &Bound<'_, PyAny>,
    allow_hyper: bool,
    rewrite: Option<String>,
    place: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyExpr> {
    let place = match place {
        Some(p) if !p.is_none() => Some(to_expr(p)?),
        _ => None,
    };
    match sp::hyperexpand_options(&to_expr(f)?, allow_hyper, rewrite.as_deref(), place.as_ref()) {
        Ok(r) => Ok(wrap(r)),
        Err(("ValueError", m)) => Err(pyo3::exceptions::PyValueError::new_err(m)),
        Err(("TypeError", m)) => Err(pyo3::exceptions::PyTypeError::new_err(m)),
        Err(("PolynomialError", m)) => {
            let class = f.py().import("symox")?.getattr("PolynomialError")?;
            Err(PyErr::from_value(class.call1((m,))?))
        }
        Err((_, m)) => Err(pyo3::exceptions::PyNotImplementedError::new_err(m)),
    }
}

/// Arbitrary precision value of a numeric expression as a decimal string
/// with `dps` significant digits (`None` if not supported).
#[pyfunction]
#[pyo3(name = "_nstr")]
fn nstr(e: &Bound<'_, PyAny>, dps: u32) -> PyResult<Option<String>> {
    let e = to_expr(e)?;
    let prec = symox::bigfloat::dps_to_prec(dps);
    Ok(symox::evalf_mp::eval_mp(&e, prec).map(|v| {
        let re = v.re.to_string_dps(dps, false);
        if v.im.is_zero() {
            re
        } else {
            format!("{re} + {}*I", v.im.to_string_dps(dps, false))
        }
    }))
}

/// Direct double precision evaluation `_numeric(name, *args)` (complex in,
/// complex out) for testing and benchmarking.
#[pyfunction]
#[pyo3(name = "_numeric", signature = (name, *args))]
fn numeric<'py>(py: Python<'py>, name: &str, args: &Bound<'py, PyTuple>) -> PyResult<Option<Bound<'py, PyComplex>>> {
    let f = SpecialFn::from_name(name).ok_or_else(|| PyValueError::new_err(format!("unknown special function {name}")))?;
    if !f.nargs().contains(&args.len()) {
        return Err(PyTypeError::new_err(format!(
            "{name} accepts argument counts {:?}, got {}", f.nargs(), args.len()
        )));
    }
    let mut v = Vec::new();
    for a in args.iter() {
        let c: num_complex::Complex64 = if let Ok(x) = a.extract::<f64>() {
            num_complex::Complex64::new(x, 0.0)
        } else if let Ok(z) = a.downcast::<PyComplex>() {
            num_complex::Complex64::new(z.real(), z.imag())
        } else {
            return Err(PyTypeError::new_err("numeric arguments expected"));
        };
        v.push(c);
    }
    Ok(symox::special::numeric::apply(f, &v).map(|z| PyComplex::from_doubles(py, z.re, z.im)))
}

/// Vectorised real evaluation `_numeric_many(name, xs)` of a unary
/// function (benchmark helper).
#[pyfunction]
#[pyo3(name = "_numeric_many")]
fn numeric_many(name: &str, xs: Vec<f64>) -> PyResult<Vec<f64>> {
    let f = SpecialFn::from_name(name).ok_or_else(|| PyValueError::new_err(format!("unknown special function {name}")))?;
    if !f.nargs().contains(&1) {
        return Err(PyTypeError::new_err(format!("{name} does not accept one argument")));
    }
    Ok(xs.into_iter().map(|x| symox::special::numeric::real1(f, x)).collect())
}

/// Names of all special function heads.
#[pyfunction]
#[pyo3(name = "_special_names")]
fn special_names() -> Vec<&'static str> {
    SpecialFn::ALL.iter().map(|s| s.name()).collect()
}

/// Partial derivative with respect to a special-function argument.
#[pyfunction]
fn _special_partial(e: &Bound<'_, PyAny>, index: usize) -> PyResult<Option<PyExpr>> {
    let e=to_expr(e)?;
    if let symox::expr::Kind::Func(symox::expr::Func::Special(s),args)=e.kind(){
        Ok(sp::diff::partial(*s,args,index).map(wrap))
    }else{Ok(None)}
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    symox::evalf_mp::register_special_host(mpmath_special);
    m.add_function(wrap_pyfunction!(_special_partial, m)?)?;
    m.add_function(wrap_pyfunction!(stirling, m)?)?;
    m.add_function(wrap_pyfunction!(multinomial_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(expand_func, m)?)?;
    m.add_function(wrap_pyfunction!(_rewrite, m)?)?;
    m.add_function(wrap_pyfunction!(besselsimp, m)?)?;
    m.add_function(wrap_pyfunction!(hyperexpand, m)?)?;
    m.add_function(wrap_pyfunction!(nstr, m)?)?;
    m.add_function(wrap_pyfunction!(numeric, m)?)?;
    m.add_function(wrap_pyfunction!(numeric_many, m)?)?;
    m.add_function(wrap_pyfunction!(special_names, m)?)?;
    Ok(())
}

static MPMATH_EVALUATE: pyo3::sync::GILOnceCell<PyObject> = pyo3::sync::GILOnceCell::new();

fn host_arg_object(py: Python<'_>, arg: &symox::evalf_mp::HostArg) -> PyResult<PyObject> {
    use symox::evalf_mp::HostArg;
    Ok(match arg {
        HostArg::Number(value) => (
            value.re.man.clone(), value.re.exp, value.im.man.clone(), value.im.exp,
        ).into_pyobject(py)?.into_any().unbind(),
        HostArg::Tuple(items) => pyo3::types::PyList::new(
            py, items.iter().map(|item| host_arg_object(py, item)).collect::<PyResult<Vec<_>>>()?,
        )?.into_any().unbind(),
    })
}

fn mpmath_evaluate(py: Python<'_>, name: &str, args: &[symox::evalf_mp::HostArg], prec: u32)
    -> PyResult<Option<symox::evalf_mp::MpComplex>> {
    let evaluate = MPMATH_EVALUATE.get_or_try_init(py, || -> PyResult<PyObject> {
        Ok(py.import("symox._mpmath_special")?.getattr("evaluate")?.unbind())
    })?;
    let values = args.iter().map(|arg| host_arg_object(py, arg)).collect::<PyResult<Vec<_>>>()?;
    let result = evaluate.call1(py, (name, pyo3::types::PyList::new(py, values)?, prec))?;
    if result.is_none(py) {
        return Ok(None);
    }
    let (re_man, re_exp, im_man, im_exp): (num_bigint::BigInt, num_bigint::BigInt, num_bigint::BigInt, num_bigint::BigInt) =
        result.extract(py)?;
    // mpmath exponents are unbounded; a value whose exponent does not fit
    // is beyond MAX_EXPONENT and is unavailable, as for native kernels.
    let (Ok(re_exp), Ok(im_exp)) = (i128::try_from(re_exp), i128::try_from(im_exp)) else {
        return Ok(None);
    };
    let part = |man: num_bigint::BigInt, exp: i128| symox::bigfloat::BigFloat { man, exp };
    Ok(Some(symox::evalf_mp::MpComplex { re: part(re_man, re_exp), im: part(im_man, im_exp) }))
}

/// Multiprecision special functions through mpmath, as SymPy evaluates them.
fn mpmath_special(name: &str, args: &[symox::evalf_mp::HostArg], prec: u32)
    -> Result<Option<symox::evalf_mp::MpComplex>, symox::evalf_mp::HostError> {
    Python::with_gil(|py| mpmath_evaluate(py, name, args, prec))
        .map_err(|error| symox::evalf_mp::HostError(Box::new(error)))
}

/// Run `f` with host errors recorded, and raise the first one as SymPy's
/// `evalf` propagates mpmath's exceptions.
pub(crate) fn with_host_errors<T>(f: impl FnOnce() -> PyResult<T>) -> PyResult<T> {
    let value = {
        let _scope = symox::evalf_mp::HostErrorScope::enter(true);
        let value = f();
        (value, symox::evalf_mp::take_host_error())
    };
    match value {
        (_, Some(error)) => match error.0.downcast::<PyErr>() {
            Ok(error) => Err(*error),
            Err(_) => Err(PyValueError::new_err("special function evaluation failed")),
        },
        (value, None) => value,
    }
}
