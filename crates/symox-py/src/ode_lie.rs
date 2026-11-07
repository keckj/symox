//! Bindings of the native `lie_group` solver (`symox::ode::lie_group`).
//!
//! The public callables `simplify`, `expand`, `cancel`, `integrate`,
//! `solve`, `div`, `doit` and `dsolve` arrive in a dict; an error they raise
//! propagates unchanged. A result the port does not model (or a native
//! decline) makes the binding return `None`, and the caller keeps the
//! reference body.

use pyo3::exceptions::{PyIndexError, PyKeyError, PyNotImplementedError, PyTypeError, PyUnboundLocalError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PySet, PyTuple};
use symox::functions::undefined;
use symox::ode::lie_group::{self, Dsolved, Exc, Hooks, Inf, LieError, LieMatch, Params, Solved};
use symox::Expr;

use crate::{to_expr, wrap as w};

pub(crate) struct LieHooks<'a, 'py> {
    pub(crate) d: &'a Bound<'py, PyDict>,
    pub(crate) foreign: bool,
}

impl<'py> LieHooks<'_, 'py> {
    fn get(&mut self, name: &str) -> PyResult<Bound<'py, PyAny>> {
        match self.d.get_item(name)? {
            Some(f) => Ok(f),
            None => {
                self.foreign = true;
                Err(PyTypeError::new_err("native lie_group hook missing"))
            }
        }
    }

    fn foreign<T>(&mut self, r: PyResult<T>) -> PyResult<T> {
        r.inspect_err(|_| self.foreign = true)
    }

    fn expr_of(&mut self, v: &Bound<'_, PyAny>) -> PyResult<Expr> {
        let r = to_expr(v);
        self.foreign(r)
    }

    fn one(&mut self, name: &str, e: &Expr) -> PyResult<Expr> {
        let f = self.get(name)?;
        let v = f.call1((w(e.clone()),))?;
        self.expr_of(&v)
    }

    fn exprs(&mut self, v: &Bound<'_, PyAny>) -> PyResult<Vec<Expr>> {
        let Ok(list) = v.downcast::<PyList>() else {
            self.foreign = true;
            return Err(PyTypeError::new_err("not a list"));
        };
        let r = list.iter().map(|item| to_expr(&item)).collect::<PyResult<Vec<_>>>();
        self.foreign(r)
    }

    fn dict(&mut self, v: &Bound<'_, PyAny>) -> PyResult<Vec<(Expr, Expr)>> {
        let Ok(d) = v.downcast::<PyDict>() else {
            self.foreign = true;
            return Err(PyTypeError::new_err("not a dict"));
        };
        let r = d.iter().map(|(k, v)| Ok((to_expr(&k)?, to_expr(&v)?))).collect::<PyResult<Vec<_>>>();
        self.foreign(r)
    }

    fn list_of_exprs<'b>(&self, py: Python<'b>, v: &[Expr]) -> PyResult<Bound<'b, PyList>> {
        PyList::new(py, v.iter().cloned().map(w))
    }
}

impl Hooks for LieHooks<'_, '_> {
    type Error = PyErr;

    fn exc(&self, e: &PyErr) -> Exc {
        if self.foreign {
            return Exc::Other;
        }
        Python::with_gil(|py| {
            if e.is_instance_of::<PyNotImplementedError>(py) {
                Exc::NotImplemented
            } else if e.is_instance_of::<PyValueError>(py) {
                Exc::Value
            } else if e.is_instance_of::<PyTypeError>(py) {
                Exc::Type
            } else if e.is_instance_of::<PyIndexError>(py) {
                Exc::Index
            } else if e.is_instance_of::<PyUnboundLocalError>(py) {
                Exc::UnboundLocal
            } else if e.is_instance_of::<PyKeyError>(py) {
                Exc::Key
            } else {
                Exc::Other
            }
        })
    }
    fn simplify(&mut self, e: &Expr) -> PyResult<Expr> {
        self.one("simplify", e)
    }
    fn expand(&mut self, e: &Expr) -> PyResult<Expr> {
        self.one("expand", e)
    }
    fn cancel(&mut self, e: &Expr) -> PyResult<Expr> {
        self.one("cancel", e)
    }
    fn integrate(&mut self, e: &Expr, x: Option<&Expr>) -> PyResult<Expr> {
        let f = self.get("integrate")?;
        let v = match x {
            Some(x) => f.call1((w(e.clone()), w(x.clone())))?,
            None => f.call1((w(e.clone()),))?,
        };
        self.expr_of(&v)
    }
    fn solve(&mut self, eq: &Expr, sym: &Expr) -> PyResult<Vec<Expr>> {
        let f = self.get("solve")?;
        let v = f.call1((w(eq.clone()), w(sym.clone())))?;
        self.exprs(&v)
    }
    fn solve_dicts(&mut self, eqs: &[Expr], syms: &[Expr]) -> PyResult<Vec<Vec<(Expr, Expr)>>> {
        let f = self.get("solve")?;
        let py = f.py();
        let mut args = vec![self.list_of_exprs(py, eqs)?.into_any()];
        args.extend(syms.iter().map(|s| w(s.clone()).into_pyobject(py)).collect::<PyResult<Vec<_>>>()?);
        let kw = PyDict::new(py);
        kw.set_item("dict", true)?;
        let v = f.call(PyTuple::new(py, args)?, Some(&kw))?;
        let Ok(list) = v.downcast::<PyList>() else {
            self.foreign = true;
            return Err(PyTypeError::new_err("not a list"));
        };
        let items: Vec<_> = list.iter().collect();
        items.iter().map(|d| self.dict(d)).collect()
    }
    fn solve_system(&mut self, eqs: &[Expr], syms: &[Expr], unpacked_set: bool) -> PyResult<Solved> {
        let f = self.get("solve")?;
        let py = f.py();
        let eqs = self.list_of_exprs(py, eqs)?.into_any();
        let v = if unpacked_set {
            // solve(eqs, *set(syms)): the set iteration order is the
            // reference one.
            let set = PySet::new(py, syms.iter().cloned().map(w).map(|e| e.into_pyobject(py)).collect::<PyResult<Vec<_>>>()?)?;
            let mut args = vec![eqs];
            args.extend(set.iter());
            f.call1(PyTuple::new(py, args)?)?
        } else {
            f.call1((eqs, self.list_of_exprs(py, syms)?))?
        };
        if v.downcast::<PyDict>().is_ok() {
            return Ok(Solved::Dict(self.dict(&v)?));
        }
        let Ok(list) = v.downcast::<PyList>() else {
            self.foreign = true;
            return Err(PyTypeError::new_err("not a dict or a list"));
        };
        let items: Vec<_> = list.iter().collect();
        Ok(Solved::Dicts(items.iter().map(|d| self.dict(d)).collect::<PyResult<Vec<_>>>()?))
    }
    fn div(&mut self, f: &Expr, g: &Expr) -> PyResult<(Expr, Expr)> {
        let h = self.get("div")?;
        let v = h.call1((w(f.clone()), w(g.clone())))?;
        let r = v.extract::<(Bound<'_, PyAny>, Bound<'_, PyAny>)>().and_then(|(q, r)| Ok((to_expr(&q)?, to_expr(&r)?)));
        self.foreign(r)
    }
    fn doit(&mut self, e: &Expr) -> PyResult<Expr> {
        self.one("doit", e)
    }
    fn dsolve(&mut self, eq: &Expr, func: &Expr, hint: Option<&str>) -> PyResult<Dsolved> {
        let f = self.get("dsolve")?;
        let py = f.py();
        let v = match hint {
            Some(h) => {
                let kw = PyDict::new(py);
                kw.set_item("hint", h)?;
                f.call((w(eq.clone()), w(func.clone())), Some(&kw))?
            }
            None => f.call1((w(eq.clone()), w(func.clone())))?,
        };
        if v.downcast::<PyList>().is_ok() {
            return Ok(Dsolved::Many(self.exprs(&v)?));
        }
        Ok(Dsolved::One(self.expr_of(&v)?))
    }
}

/// The value of a native result: `Ok(None)` for a decline or a hook result
/// the port does not model, else the reference exception.
pub(crate) fn settle<T>(hooks: &LieHooks<'_, '_>, r: Result<T, LieError<PyErr>>) -> PyResult<Option<T>> {
    match r {
        Ok(v) => Ok(Some(v)),
        Err(_) if hooks.foreign => Ok(None),
        Err(LieError::Decline(_)) => Ok(None),
        Err(LieError::Hook(e)) => Err(e),
        Err(LieError::KeyMissing(k)) => Err(PyKeyError::new_err(w(k))),
        Err(LieError::Raise(kind, msg)) => Err(match kind {
            Exc::NotImplemented => PyNotImplementedError::new_err(msg),
            Exc::Value => PyValueError::new_err(msg),
            Exc::Type => PyTypeError::new_err(msg),
            Exc::Index => PyIndexError::new_err(msg),
            Exc::UnboundLocal => PyUnboundLocalError::new_err(msg),
            Exc::Key => PyKeyError::new_err(msg),
            Exc::Other => return Ok(None),
        }),
    }
}

/// The reference dicts `{xi(x, f(x)): xi, eta(x, f(x)): eta}`.
fn inf_list(py: Python<'_>, infs: Vec<Inf>, fx: &Expr, x: &Expr) -> PyResult<PyObject> {
    let xi = undefined("xi", vec![x.clone(), fx.clone()]);
    let eta = undefined("eta", vec![x.clone(), fx.clone()]);
    let out = PyList::empty(py);
    for inf in infs {
        let d = PyDict::new(py);
        if inf.eta_first {
            d.set_item(w(eta.clone()), w(inf.eta))?;
            d.set_item(w(xi.clone()), w(inf.xi))?;
        } else {
            d.set_item(w(xi.clone()), w(inf.xi))?;
            d.set_item(w(eta.clone()), w(inf.eta))?;
        }
        out.append(d)?;
    }
    Ok(out.into_any().unbind())
}

/// `infinitesimals(eq, func, order=1, hint, match)` after the argument
/// checks; `h`/`y` are those of a supplied `match`. `None` declines.
#[pyfunction]
#[pyo3(signature = (eq, fx, x, hint, h, y, hooks))]
#[allow(clippy::too_many_arguments)]
fn ode_lie_infinitesimals(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    fx: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    hint: &str,
    h: Option<&Bound<'_, PyAny>>,
    y: Option<&Bound<'_, PyAny>>,
    hooks: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    let (eq, fx, x) = (to_expr(eq)?, to_expr(fx)?, to_expr(x)?);
    let given = match (h, y) {
        (Some(h), Some(y)) => match (to_expr(h), to_expr(y)) {
            (Ok(h), Ok(y)) => Some((h, y)),
            _ => return Ok(py.None()),
        },
        _ => None,
    };
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = lie_group::infinitesimals(&eq, &fx, &x, hint, given, &mut hk);
    match settle(&hk, r)? {
        Some(v) => inf_list(py, v, &fx, &x),
        None => Ok(py.None()),
    }
}

/// `lie_heuristic_<name>(match, comp)`: `(result,)` with `result` the list
/// of infinitesimals or `None`; `None` declines.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn ode_lie_heuristic(
    py: Python<'_>,
    name: &str,
    h: &Bound<'_, PyAny>,
    hx: &Bound<'_, PyAny>,
    hy: &Bound<'_, PyAny>,
    y: &Bound<'_, PyAny>,
    hinv: &Bound<'_, PyAny>,
    fx: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    comp: bool,
    hooks: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    let m = LieMatch { h: to_expr(h)?, hx: to_expr(hx)?, hy: to_expr(hy)?, y: to_expr(y)?, hinv: to_expr(hinv)?, fx: to_expr(fx)?, x: to_expr(x)? };
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = lie_group::heuristic(name, &m, comp, &mut hk);
    match settle(&hk, r)? {
        Some(Some(v)) => Ok(PyTuple::new(py, [inf_list(py, v, &m.fx, &m.x)?])?.into_any().unbind()),
        Some(None) => Ok(PyTuple::new(py, [py.None()])?.into_any().unbind()),
        None => Ok(py.None()),
    }
}

/// `checkinfsol(eq, infinitesimals, func, order=1)` for the `(xi, eta)`
/// pairs: the list of `(bool, expr)`, `None` declines.
#[pyfunction]
fn ode_lie_checkinfsol(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    fx: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    infs: Vec<(Bound<'_, PyAny>, Bound<'_, PyAny>)>,
    hooks: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    let (eq, fx, x) = (to_expr(eq)?, to_expr(fx)?, to_expr(x)?);
    let mut pairs = Vec::new();
    for (a, b) in infs {
        match (to_expr(&a), to_expr(&b)) {
            (Ok(a), Ok(b)) => pairs.push((a, b)),
            _ => return Ok(py.None()),
        }
    }
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = lie_group::checkinfsol(&eq, &fx, &x, &pairs, &mut hk);
    match settle(&hk, r)? {
        Some(v) => {
            let out = PyList::empty(py);
            // SymPy appends `(True, 0)` with the Python int `0`.
            for (ok, e) in v {
                let second = if ok { 0i64.into_pyobject(py)?.into_any().unbind() } else { w(e).into_pyobject(py)?.into_any().unbind() };
                out.append(PyTuple::new(py, [ok.into_pyobject(py)?.to_owned().into_any().unbind(), second])?)?;
            }
            Ok(out.into_any().unbind())
        }
        None => Ok(py.None()),
    }
}

/// `_lie_group_remove(coords)`.
#[pyfunction]
fn ode_lie_group_remove(coords: &Bound<'_, PyAny>) -> PyResult<crate::PyExpr> {
    Ok(w(lie_group::lie_group_remove(&to_expr(coords)?)))
}

/// `LieGroup._matches` split: `(d, e)` with `f(x)` replaced by `y`, `False`
/// without a match, `None` declines.
#[pyfunction]
fn ode_lie_group_split(py: Python<'_>, eq: &Bound<'_, PyAny>, fx: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, y: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let (Ok(eq), fx, x, y) = (to_expr(eq), to_expr(fx)?, to_expr(x)?, to_expr(y)?) else { return Ok(py.None()) };
    Ok(match lie_group::lie_group_split(&eq, &fx, &x, &y) {
        Ok(Some((d, e))) => PyTuple::new(py, [w(d), w(e)])?.into_any().unbind(),
        Ok(None) => false.into_pyobject(py)?.to_owned().into_any().unbind(),
        Err(_) => py.None(),
    })
}

/// `LieGroup._get_general_solution`: the list of solutions, `None` declines.
/// `user` is the supplied `(xi, eta)`, `split` the `(y, d, e)` of the match.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn ode_lie_group_solve(
    py: Python<'_>,
    eq: &Bound<'_, PyAny>,
    fx: &Bound<'_, PyAny>,
    x: &Bound<'_, PyAny>,
    user: Option<(Bound<'_, PyAny>, Bound<'_, PyAny>)>,
    split: Option<(Bound<'_, PyAny>, Bound<'_, PyAny>, Bound<'_, PyAny>)>,
    hooks: &Bound<'_, PyDict>,
) -> PyResult<PyObject> {
    let (eq, fx, x) = (to_expr(eq)?, to_expr(fx)?, to_expr(x)?);
    let user = match user {
        Some((a, b)) => match (to_expr(&a), to_expr(&b)) {
            (Ok(a), Ok(b)) => Some((a, b)),
            _ => return Ok(py.None()),
        },
        None => None,
    };
    let split = match split {
        Some((y, d, e)) => match (to_expr(&y), to_expr(&d), to_expr(&e)) {
            (Ok(y), Ok(d), Ok(e)) => Some((y, d, e)),
            _ => return Ok(py.None()),
        },
        None => None,
    };
    let mut hk = LieHooks { d: hooks, foreign: false };
    let r = lie_group::lie_group_solve(&eq, &fx, &x, Params { user, split }, &mut hk);
    match settle(&hk, r)? {
        Some(v) => Ok(PyList::new(py, v.into_iter().map(w))?.into_any().unbind()),
        None => Ok(py.None()),
    }
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ode_lie_infinitesimals, m)?)?;
    m.add_function(wrap_pyfunction!(ode_lie_heuristic, m)?)?;
    m.add_function(wrap_pyfunction!(ode_lie_checkinfsol, m)?)?;
    m.add_function(wrap_pyfunction!(ode_lie_group_remove, m)?)?;
    m.add_function(wrap_pyfunction!(ode_lie_group_split, m)?)?;
    m.add_function(wrap_pyfunction!(ode_lie_group_solve, m)?)?;
    Ok(())
}
