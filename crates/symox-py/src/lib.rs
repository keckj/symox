//! Python bindings: a SymPy compatible API on top of the Rust core.

// SymPy's API uses `is_Add`, `as_coeff_Mul`, ...
#![allow(non_snake_case)]

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use num_bigint::BigInt;
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1};
use pyo3::exceptions::{PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyComplex, PyDict, PyFloat, PyInt, PyList, PySet, PyString, PyTuple};
use rustc_hash::FxHashMap;

mod series;

use symox::expr::{Constant, Expr, Func, Kind, RelOp};
use symox::lambdify::Compiled;
use symox::number::Number;
mod printing;

use symox::{diff, evalf, expand, functions as fx, ops, parser, printer, relational, simplify, subs, traverse};

mod poly;
mod ntheory;
mod simp;
mod ode;
mod solve;

/// A symbolic expression.
#[pyclass(name = "Expr", module = "symox._core", frozen, weakref)]
#[derive(Clone)]
pub struct PyExpr {
    pub e: Expr,
}

pub(crate) fn wrap(e: Expr) -> PyExpr {
    PyExpr { e }
}

/// Convert any supported Python object into an expression.
pub(crate) fn to_expr(obj: &Bound<'_, PyAny>) -> PyResult<Expr> {
    if let Ok(p) = obj.downcast::<PyExpr>() {
        return Ok(p.get().e.clone());
    }
    if obj.is_instance_of::<pyo3::types::PyBool>() {
        return Ok(Expr::boolean(obj.extract::<bool>()?));
    }
    if obj.is_instance_of::<PyInt>() {
        if let Ok(i) = obj.extract::<i64>() {
            return Ok(Expr::int(i));
        }
        let b: BigInt = obj.extract()?;
        return Ok(Expr::num(Number::from_bigint(b)));
    }
    if obj.is_instance_of::<PyFloat>() {
        return Ok(Expr::float(obj.extract::<f64>()?));
    }
    if let Ok(c) = obj.downcast::<PyComplex>() {
        let (re, im) = (c.real(), c.imag());
        return Ok(evalf::c64_to_expr(num_complex_c64(re, im)));
    }
    if obj.is_instance_of::<PyString>() {
        let s: String = obj.extract()?;
        return parser::parse_str(&s).map_err(|e| PyValueError::new_err(e.to_string()));
    }
    // fractions.Fraction
    let ty = obj.get_type();
    let tname = ty.name()?.to_string();
    let module = ty.module().map(|m| m.to_string()).unwrap_or_default();
    if tname == "Fraction" {
        let n: BigInt = obj.getattr("numerator")?.extract()?;
        let d: BigInt = obj.getattr("denominator")?.extract()?;
        return Number::rational(n, d)
            .map(Expr::num)
            .ok_or_else(|| PyZeroDivisionError::new_err("zero denominator"));
    }
    // numpy scalars
    if module == "numpy" {
        if let Ok(i) = obj.extract::<i64>() {
            if tname.contains("int") {
                return Ok(Expr::int(i));
            }
        }
        if let Ok(f) = obj.extract::<f64>() {
            return Ok(Expr::float(f));
        }
    }
    // SymPy objects (and anything with a _symox_ hook)
    if let Ok(hook) = obj.getattr("_symox_") {
        return to_expr(&hook.call0()?);
    }
    if module.starts_with("sympy") {
        let py = obj.py();
        let srepr = py.import("sympy")?.getattr("srepr")?.call1((obj,))?;
        let s: String = srepr.extract()?;
        return parser::parse_str(&s).map_err(|e| {
            PyValueError::new_err(format!("cannot convert SymPy object {s}: {e}"))
        });
    }
    Err(PyTypeError::new_err(format!(
        "cannot convert object of type {tname} to a symox expression"
    )))
}

fn num_complex_c64(re: f64, im: f64) -> num_complex::Complex64 {
    num_complex::Complex64::new(re, im)
}

fn opt_expr(obj: &Bound<'_, PyAny>) -> Option<Expr> {
    to_expr(obj).ok()
}

/// Parse the variable list of `diff(expr, x, y, (z, 2), x, 3)`.
fn diff_vars(args: &Bound<'_, PyTuple>) -> PyResult<Vec<(Expr, u32)>> {
    let mut vars: Vec<(Expr, u32)> = Vec::new();
    let items: Vec<Bound<'_, PyAny>> = args.iter().collect();
    let mut i = 0;
    while i < items.len() {
        let it = &items[i];
        if let Ok(t) = it.downcast::<PyTuple>() {
            if t.len() != 2 {
                return Err(PyValueError::new_err("expected (symbol, order) tuples"));
            }
            let v = to_expr(&t.get_item(0)?)?;
            let n: u32 = t.get_item(1)?.extract()?;
            vars.push((v, n));
            i += 1;
            continue;
        }
        if it.is_instance_of::<PyInt>() {
            let n: u32 = it.extract()?;
            match vars.last_mut() {
                Some(last) => last.1 = last.1 + n - 1,
                None => return Err(PyValueError::new_err("derivative order without a variable")),
            }
            i += 1;
            continue;
        }
        vars.push((to_expr(it)?, 1));
        i += 1;
    }
    Ok(vars)
}

fn do_diff(e: &Expr, args: &Bound<'_, PyTuple>) -> PyResult<Expr> {
    let vars = diff_vars(args)?;
    if vars.is_empty() {
        let fs = traverse::free_symbols(e);
        return match fs.len() {
            0 => Ok(Expr::zero()),
            1 => Ok(diff::diff(e, &fs[0])),
            _ => Err(PyValueError::new_err(
                "Since there is more than one variable in the expression, the variable(s) of differentiation must be supplied",
            )),
        };
    }
    Ok(diff::diff_many(e, &vars))
}

fn subs_pairs(args: &Bound<'_, PyTuple>) -> PyResult<Vec<(Expr, Expr)>> {
    let mut pairs = Vec::new();
    match args.len() {
        1 => {
            let a = args.get_item(0)?;
            if let Ok(d) = a.downcast::<PyDict>() {
                for (k, v) in d.iter() {
                    pairs.push((to_expr(&k)?, to_expr(&v)?));
                }
            } else {
                for item in a.try_iter()? {
                    let item = item?;
                    let t = item.downcast::<PyTuple>().map_err(|_| {
                        PyTypeError::new_err("subs expects a dict, an iterable of pairs or (old, new)")
                    })?;
                    pairs.push((to_expr(&t.get_item(0)?)?, to_expr(&t.get_item(1)?)?));
                }
            }
        }
        2 => pairs.push((to_expr(&args.get_item(0)?)?, to_expr(&args.get_item(1)?)?)),
        _ => return Err(PyTypeError::new_err("subs expects 1 or 2 arguments")),
    }
    Ok(pairs)
}

#[pymethods]
impl PyExpr {
    #[new]
    fn py_new(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(wrap(to_expr(obj)?))
    }

    fn __str__(&self) -> String {
        printer::to_str(&self.e)
    }

    fn __repr__(&self) -> String {
        printer::to_str(&self.e)
    }

    fn _repr_latex_(&self) -> String {
        printing::repr_latex(self)
    }

    fn srepr(&self) -> String {
        printer::srepr(&self.e)
    }

    fn __hash__(&self) -> u64 {
        let mut h = DefaultHasher::new();
        self.e.hash(&mut h);
        h.finish()
    }

    fn __richcmp__(&self, other: &Bound<'_, PyAny>, op: pyo3::basic::CompareOp) -> PyResult<PyObject> {
        let py = other.py();
        use pyo3::basic::CompareOp::*;
        let Some(o) = opt_expr(other) else {
            return Ok(match op {
                Eq => false.into_pyobject(py)?.to_owned().into_any().unbind(),
                Ne => true.into_pyobject(py)?.to_owned().into_any().unbind(),
                _ => py.NotImplemented(),
            });
        };
        Ok(match op {
            Eq => (self.e == o).into_pyobject(py)?.to_owned().into_any().unbind(),
            Ne => (self.e != o).into_pyobject(py)?.to_owned().into_any().unbind(),
            Lt => wrap(relational::lt(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
            Le => wrap(relational::le(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
            Gt => wrap(relational::gt(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
            Ge => wrap(relational::ge(&self.e, &o)).into_pyobject(py)?.into_any().unbind(),
        })
    }

    fn __bool__(&self) -> PyResult<bool> {
        match self.e.kind() {
            Kind::Bool(b) => Ok(*b),
            Kind::Num(n) => Ok(!n.is_zero()),
            Kind::Rel(..) => Err(PyTypeError::new_err("cannot determine truth value of Relational")),
            _ => Ok(true),
        }
    }

    // arithmetic -------------------------------------------------------

    fn __add__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::add2(&self.e, &to_expr(o)?)))
    }
    fn __radd__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::add2(&to_expr(o)?, &self.e)))
    }
    fn __sub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::sub(&self.e, &to_expr(o)?)))
    }
    fn __rsub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::sub(&to_expr(o)?, &self.e)))
    }
    fn __mul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::mul2(&self.e, &to_expr(o)?)))
    }
    fn __rmul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::mul2(&to_expr(o)?, &self.e)))
    }
    fn __truediv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::div(&self.e, &to_expr(o)?)))
    }
    fn __rtruediv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(wrap(ops::div(&to_expr(o)?, &self.e)))
    }
    fn __pow__(&self, o: &Bound<'_, PyAny>, m: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        if m.is_some_and(|m| !m.is_none()) {
            return Err(PyTypeError::new_err("modular pow is not supported"));
        }
        Ok(wrap(ops::pow(&self.e, &to_expr(o)?)))
    }
    fn __rpow__(&self, o: &Bound<'_, PyAny>, m: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        if m.is_some_and(|m| !m.is_none()) {
            return Err(PyTypeError::new_err("modular pow is not supported"));
        }
        Ok(wrap(ops::pow(&to_expr(o)?, &self.e)))
    }
    fn __neg__(&self) -> PyExpr {
        wrap(ops::neg(&self.e))
    }
    fn __pos__(&self) -> PyExpr {
        self.clone()
    }
    fn __abs__(&self) -> PyExpr {
        wrap(fx::abs(&self.e))
    }

    fn __float__(&self) -> PyResult<f64> {
        evalf::eval_f64(&self.e).ok_or_else(|| PyTypeError::new_err(format!("cannot convert {} to float", self.e)))
    }

    fn __complex__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyComplex>> {
        let z = evalf::eval_c64(&self.e)
            .ok_or_else(|| PyTypeError::new_err(format!("cannot convert {} to complex", self.e)))?;
        Ok(PyComplex::from_doubles(py, z.re, z.im))
    }

    fn __int__(&self) -> PyResult<BigInt> {
        match self.e.kind() {
            Kind::Num(n) => {
                if let Some(b) = n.to_bigint() {
                    return Ok(b);
                }
                let f = n.to_f64().trunc();
                if let Some((num, den)) = n.numer_denom() {
                    return Ok(&num / &den);
                }
                Ok(BigInt::from(f as i64))
            }
            _ => {
                let f = self.__float__()?;
                Ok(BigInt::from(f.trunc() as i64))
            }
        }
    }

    fn __index__(&self) -> PyResult<BigInt> {
        self.e
            .as_number()
            .and_then(|n| n.to_bigint())
            .ok_or_else(|| PyTypeError::new_err("only integers can be used as an index"))
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<(Bound<'py, PyAny>, (String,))> {
        let f = py.import("symox")?.getattr("sympify")?;
        Ok((f, (printer::srepr(&self.e),)))
    }

    // structure --------------------------------------------------------

    #[getter]
    fn args<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.e.args().into_iter().map(wrap))
    }

    /// SymPy class name of the node (`Add`, `Mul`, `Symbol`, `sin`, ...).
    #[getter]
    fn class_name(&self) -> String {
        self.e.class_name().to_string()
    }

    /// The head of the expression as a callable rebuilding it from args.
    #[getter]
    fn func(&self, py: Python<'_>) -> PyResult<PyObject> {
        let head = Head::of(&self.e);
        Ok(Py::new(py, head)?.into_any())
    }

    #[getter]
    fn name(&self) -> PyResult<String> {
        match self.e.kind() {
            Kind::Sym(s) => Ok(s.to_string()),
            Kind::Func(f, _) => Ok(f.name().to_string()),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("expression has no name")),
        }
    }

    #[getter]
    fn p(&self) -> PyResult<BigInt> {
        self.e
            .as_number()
            .and_then(|n| n.numer_denom())
            .map(|x| x.0)
            .ok_or_else(|| pyo3::exceptions::PyAttributeError::new_err("not a rational number"))
    }

    #[getter]
    fn q(&self) -> PyResult<BigInt> {
        self.e
            .as_number()
            .and_then(|n| n.numer_denom())
            .map(|x| x.1)
            .ok_or_else(|| pyo3::exceptions::PyAttributeError::new_err("not a rational number"))
    }

    #[getter]
    fn free_symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PySet>> {
        PySet::new(py, traverse::free_symbols(&self.e).into_iter().map(wrap))
    }

    #[getter]
    fn is_Atom(&self) -> bool {
        self.e.is_atom()
    }
    #[getter]
    fn is_Symbol(&self) -> bool {
        self.e.is_symbol()
    }
    #[getter]
    fn is_Number(&self) -> bool {
        self.e.is_number()
    }
    #[getter]
    fn is_Integer(&self) -> bool {
        self.e.is_integer_number()
    }
    #[getter]
    fn is_Rational(&self) -> bool {
        self.e.is_rational_number()
    }
    #[getter]
    fn is_Float(&self) -> bool {
        matches!(self.e.as_number(), Some(Number::Float(_)))
    }
    #[getter]
    fn is_Add(&self) -> bool {
        self.e.is_add()
    }
    #[getter]
    fn is_Mul(&self) -> bool {
        self.e.is_mul()
    }
    #[getter]
    fn is_Pow(&self) -> bool {
        self.e.is_pow()
    }
    #[getter]
    fn is_Function(&self) -> bool {
        matches!(self.e.kind(), Kind::Func(..))
    }
    #[getter]
    fn is_Derivative(&self) -> bool {
        matches!(self.e.kind(), Kind::Derivative(..))
    }
    #[getter]
    fn is_Relational(&self) -> bool {
        matches!(self.e.kind(), Kind::Rel(..))
    }
    #[getter]
    fn is_Piecewise(&self) -> bool {
        matches!(self.e.kind(), Kind::Piecewise(..))
    }
    /// True if the expression has no free symbols (it is a number).
    #[getter]
    fn is_number(&self) -> bool {
        fx::is_numeric_value(&self.e)
    }
    #[getter]
    fn is_zero(&self) -> Option<bool> {
        symox::assumptions::is_zero(&self.e)
    }
    #[getter]
    fn is_nonzero(&self) -> Option<bool> {
        symox::assumptions::is_nonzero(&self.e)
    }
    #[getter]
    fn is_integer(&self) -> Option<bool> {
        symox::assumptions::is_integer(&self.e)
    }
    #[getter]
    fn is_rational(&self) -> Option<bool> {
        symox::assumptions::is_rational(&self.e)
    }
    #[getter]
    fn is_real(&self) -> Option<bool> {
        symox::assumptions::is_real(&self.e)
    }
    #[getter]
    fn is_positive(&self) -> Option<bool> {
        symox::assumptions::is_positive(&self.e)
    }
    #[getter]
    fn is_negative(&self) -> Option<bool> {
        symox::assumptions::is_negative(&self.e)
    }
    #[getter]
    fn is_nonnegative(&self) -> Option<bool> {
        symox::assumptions::is_nonnegative(&self.e)
    }
    #[getter]
    fn is_nonpositive(&self) -> Option<bool> {
        symox::assumptions::is_nonpositive(&self.e)
    }
    #[getter]
    fn is_even(&self) -> Option<bool> {
        symox::assumptions::is_even(&self.e)
    }
    #[getter]
    fn is_odd(&self) -> Option<bool> {
        symox::assumptions::is_odd(&self.e)
    }
    #[getter]
    fn is_finite(&self) -> Option<bool> {
        symox::assumptions::is_finite(&self.e)
    }
    #[getter]
    fn is_complex(&self) -> Option<bool> {
        symox::assumptions::is_finite(&self.e)
    }
    /// The assumptions of a symbol as a dict (SymPy's `assumptions0`).
    #[getter]
    fn assumptions0<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        if let Some(s) = self.e.as_sym() {
            for (n, f) in symox::Assumptions::NAMES {
                if s.assumptions().has(*f) {
                    match n.strip_prefix('!') {
                        Some(m) => d.set_item(m, false)?,
                        None => d.set_item(*n, true)?,
                    }
                }
            }
        }
        Ok(d)
    }
    /// Left/right hand sides of a relational.
    #[getter]
    fn lhs(&self) -> PyResult<PyExpr> {
        match self.e.kind() {
            Kind::Rel(_, a, _) => Ok(wrap(a.clone())),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("not a relational")),
        }
    }
    #[getter]
    fn rhs(&self) -> PyResult<PyExpr> {
        match self.e.kind() {
            Kind::Rel(_, _, b) => Ok(wrap(b.clone())),
            _ => Err(pyo3::exceptions::PyAttributeError::new_err("not a relational")),
        }
    }

    // algorithms -------------------------------------------------------

    #[pyo3(signature = (*args))]
    fn diff(&self, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
        Ok(wrap(do_diff(&self.e, args)?))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn subs(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
        let pairs = subs_pairs(args)?;
        Ok(wrap(subs::subs(&self.e, &pairs)))
    }

    /// `expr.match(pattern)`: dict of wildcard bindings or None.
    #[pyo3(name = "match")]
    fn py_match<'py>(&self, py: Python<'py>, pattern: &Bound<'py, PyAny>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let p = to_expr(pattern)?;
        match symox::pattern::pattern_match(&self.e, &p) {
            None => Ok(None),
            Some(m) => {
                let d = PyDict::new(py);
                for (k, v) in m {
                    d.set_item(wrap(k), wrap(v))?;
                }
                Ok(Some(d))
            }
        }
    }

    /// `expr.replace(pattern, value)`: pattern may contain Wilds, value may
    /// be an expression or a callable receiving the bindings by name; a
    /// function head (e.g. `sin`) may be replaced by another callable.
    fn replace(&self, py: Python<'_>, pattern: &Bound<'_, PyAny>, value: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        let pat = to_expr(pattern);
        match pat {
            Ok(p) => {
                if let Ok(v) = to_expr(value) {
                    return Ok(wrap(symox::pattern::replace(&self.e, &p, &v)));
                }
                if !value.is_callable() {
                    return Err(PyTypeError::new_err("replacement must be an expression or a callable"));
                }
                let mut err: Option<PyErr> = None;
                let r = symox::pattern::replace_with(&self.e, &p, &mut |bind, _| {
                    let kw = PyDict::new(py);
                    for (k, v) in bind {
                        let name = k.as_sym().map(|s| s.name().to_string()).unwrap_or_default();
                        if kw.set_item(name, wrap(v.clone())).is_err() {
                            return None;
                        }
                    }
                    match value.call((), Some(&kw)).and_then(|r| to_expr(&r)) {
                        Ok(e) => Some(e),
                        Err(e) => {
                            err = Some(e);
                            None
                        }
                    }
                });
                if let Some(e) = err {
                    return Err(e);
                }
                Ok(wrap(r))
            }
            Err(_) => {
                // function head replacement: replace(sin, cos)
                let name: String = pattern.getattr("__name__")?.extract()?;
                let mut err: Option<PyErr> = None;
                let r = replace_heads(&self.e, &name, &mut |args: Vec<Expr>| {
                    let t = match PyTuple::new(py, args.into_iter().map(wrap)) {
                        Ok(t) => t,
                        Err(e) => {
                            err = Some(e);
                            return None;
                        }
                    };
                    match value.call1(t).and_then(|r| to_expr(&r)) {
                        Ok(e) => Some(e),
                        Err(e) => {
                            err = Some(e);
                            None
                        }
                    }
                });
                if let Some(e) = err {
                    return Err(e);
                }
                Ok(wrap(r))
            }
        }
    }

    /// Sub-expressions matching a pattern.
    fn find<'py>(&self, py: Python<'py>, pattern: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PySet>> {
        let p = to_expr(pattern)?;
        PySet::new(py, symox::pattern::find(&self.e, &p).into_iter().map(wrap))
    }

    fn xreplace(&self, rule: &Bound<'_, PyDict>) -> PyResult<PyExpr> {
        let mut m = FxHashMap::default();
        for (k, v) in rule.iter() {
            m.insert(to_expr(&k)?, to_expr(&v)?);
        }
        Ok(wrap(subs::xreplace(&self.e, &m)))
    }

    #[pyo3(signature = (deep=true, **hints))]
    fn expand(&self, deep: bool, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
        let _ = deep;
        Ok(wrap(expand_hints(&self.e, hints)?))
    }

    fn trigsimp(&self) -> PyExpr {
        wrap(symox::simp::trigsimp(&self.e))
    }

    fn powsimp(&self) -> PyExpr {
        wrap(symox::simp::powsimp(&self.e))
    }

    fn nsimplify(&self) -> PyExpr {
        wrap(symox::simp::nsimplify(&self.e, &[]))
    }

    #[pyo3(signature = (n=15, subs=None, **_kwargs))]
    fn evalf(
        &self,
        n: u32,
        subs: Option<&Bound<'_, PyDict>>,
        _kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyExpr> {
        let _ = n;
        let e = match subs {
            Some(d) => {
                let mut pairs = Vec::new();
                for (k, v) in d.iter() {
                    pairs.push((to_expr(&k)?, to_expr(&v)?));
                }
                symox::subs::subs(&self.e, &pairs)
            }
            None => self.e.clone(),
        };
        Ok(wrap(evalf::evalf(&e)))
    }

    #[pyo3(signature = (n=15, subs=None, **kwargs))]
    fn n(
        &self,
        n: u32,
        subs: Option<&Bound<'_, PyDict>>,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<PyExpr> {
        self.evalf(n, subs, kwargs)
    }

    fn simplify(&self) -> PyExpr {
        wrap(simplify::simplify(&self.e))
    }

    #[pyo3(signature = (*args, **_kwargs))]
    fn integrate(&self, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
        Ok(wrap(do_integrate(&self.e, args)?))
    }

    fn together(&self) -> PyExpr {
        wrap(simplify::together(&self.e))
    }

    fn radsimp(&self) -> PyExpr {
        wrap(simplify::radsimp(&self.e))
    }

    #[pyo3(signature = (*patterns))]
    fn has(&self, patterns: &Bound<'_, PyTuple>) -> PyResult<bool> {
        for p in patterns.iter() {
            if traverse::has(&self.e, &to_expr(&p)?) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn as_numer_denom<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (n, d) = simplify::numer_denom(&self.e);
        PyTuple::new(py, [wrap(n), wrap(d)])
    }

    fn as_base_exp<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (b, x) = match self.e.kind() {
            Kind::Func(Func::Exp, a) => (Expr::e(), a[0].clone()),
            _ => ops::as_base_exp(&self.e),
        };
        PyTuple::new(py, [wrap(b), wrap(x)])
    }

    fn as_coeff_Mul<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (c, r) = ops::as_coeff_mul(&self.e);
        PyTuple::new(py, [wrap(Expr::num(c)), wrap(r)])
    }

    fn as_coeff_Add<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (c, r) = ops::as_coeff_add(&self.e);
        PyTuple::new(py, [wrap(Expr::num(c)), wrap(r)])
    }

    fn as_ordered_terms(&self) -> Vec<PyExpr> {
        printer::ordered_terms(&self.e).into_iter().map(wrap).collect()
    }

    fn as_coefficients_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        let push = |t: Expr, k: Number| -> PyResult<()> {
            let key = wrap(t);
            let prev = d.get_item(key.clone())?;
            let val = match prev {
                Some(p) => ops::add2(&p.downcast::<PyExpr>()?.get().e, &Expr::num(k)),
                None => Expr::num(k),
            };
            d.set_item(key, wrap(val))
        };
        match self.e.kind() {
            Kind::Add(c, terms) => {
                for (t, k) in terms {
                    push(t.clone(), k.clone())?;
                }
                if !c.is_zero() {
                    push(Expr::one(), c.clone())?;
                }
            }
            _ => {
                let (c, r) = ops::as_coeff_mul(&self.e);
                push(r, c)?;
            }
        }
        Ok(d)
    }

    /// Coefficient of `x**n` in the expanded expression.
    #[pyo3(signature = (x, n=None))]
    fn coeff(&self, x: &Bound<'_, PyAny>, n: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
        let x = to_expr(x)?;
        let n = match n {
            Some(n) => to_expr(n)?,
            None => Expr::one(),
        };
        Ok(wrap(coeff_of(&self.e, &x, &n)))
    }

    /// Convert to a SymPy expression.
    fn to_sympy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let sympy = py.import("sympy")?;
        let s = printer::srepr(&self.e);
        let ns = PyDict::new(py);
        ns.set_item("__builtins__", py.import("builtins")?)?;
        py.run(&std::ffi::CString::new("from sympy import *")?, Some(&ns), None)?;
        let _ = sympy;
        py.eval(&std::ffi::CString::new(s)?, Some(&ns), None)
    }

    fn _sympy_<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.to_sympy(py)
    }

    fn __getnewargs__(&self) -> (String,) {
        (printer::srepr(&self.e),)
    }
}

/// Replace applications of the function named `name` bottom-up.
fn replace_heads(e: &Expr, name: &str, f: &mut dyn FnMut(Vec<Expr>) -> Option<Expr>) -> Expr {
    let rebuilt = if e.is_atom() {
        e.clone()
    } else {
        traverse::map_children(e, &mut |c| replace_heads(c, name, f))
    };
    if let Kind::Func(g, args) = rebuilt.kind() {
        if g.name() == name {
            if let Some(r) = f(args.clone()) {
                return r;
            }
        }
    }
    rebuilt
}

/// Coefficient of `x**n` in `e` (SymPy `Expr.coeff` semantics: the sum
/// of the cofactors of terms containing exactly `x**n` and nothing else
/// depending on `x`; `n = 0` gives the terms independent of `x`).
fn coeff_of(e: &Expr, x: &Expr, n: &Expr) -> Expr {
    let terms: Vec<Expr> = match e.kind() {
        Kind::Add(..) => e.args(),
        _ => vec![e.clone()],
    };
    let mut out = Vec::new();
    if n.is_zero() {
        for t in &terms {
            if !traverse::has(t, x) {
                out.push(t.clone());
            }
        }
        return ops::add(&out);
    }
    let target = ops::pow(x, n);
    for t in &terms {
        if *t == target {
            out.push(Expr::one());
            continue;
        }
        let factors = match t.kind() {
            Kind::Mul(..) => t.args(),
            _ => continue,
        };
        if let Some(i) = factors.iter().position(|f| *f == target) {
            let rest: Vec<Expr> = factors
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, f)| f.clone())
                .collect();
            let rest = ops::mul(&rest);
            if !traverse::has(&rest, x) {
                out.push(rest);
            }
        }
    }
    ops::add(&out)
}

/// The head (`func`) of an expression: calling it with arguments rebuilds a
/// node of the same kind.
#[pyclass(name = "Head", module = "symox._core", frozen)]
pub struct Head {
    kind: HeadKind,
}

enum HeadKind {
    Add,
    Mul,
    Pow,
    Func(Func),
    Rel(RelOp),
    Atom(Expr),
    Derivative,
    Piecewise,
}

impl Head {
    fn of(e: &Expr) -> Head {
        let kind = match e.kind() {
            Kind::Add(..) => HeadKind::Add,
            Kind::Mul(..) => HeadKind::Mul,
            Kind::Pow(..) => HeadKind::Pow,
            Kind::Func(f, _) => HeadKind::Func(f.clone()),
            Kind::Rel(op, ..) => HeadKind::Rel(*op),
            Kind::Derivative(..) => HeadKind::Derivative,
            Kind::Piecewise(..) => HeadKind::Piecewise,
            _ => HeadKind::Atom(e.clone()),
        };
        Head { kind }
    }
}

#[pymethods]
impl Head {
    #[pyo3(signature = (*args))]
    fn __call__(&self, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
        let v: Vec<Expr> = args.iter().map(|a| to_expr(&a)).collect::<PyResult<_>>()?;
        Ok(wrap(match &self.kind {
            HeadKind::Add => ops::add(&v),
            HeadKind::Mul => ops::mul(&v),
            HeadKind::Pow => ops::pow(&v[0], &v[1]),
            HeadKind::Func(f) => fx::func(f.clone(), v),
            HeadKind::Rel(op) => relational::rel(*op, &v[0], &v[1]),
            HeadKind::Atom(e) => e.clone(),
            HeadKind::Derivative => {
                let vars = v[1..].iter().map(|x| (x.clone(), 1)).collect();
                diff::derivative_unevaluated(&v[0], vars)
            }
            HeadKind::Piecewise => {
                let pairs = v.chunks(2).map(|c| (c[0].clone(), c[1].clone())).collect();
                relational::piecewise(pairs)
            }
        }))
    }

    fn __repr__(&self) -> String {
        match &self.kind {
            HeadKind::Add => "Add".into(),
            HeadKind::Mul => "Mul".into(),
            HeadKind::Pow => "Pow".into(),
            HeadKind::Func(f) => f.name().into(),
            HeadKind::Rel(op) => format!("Rel({})", op.symbol()),
            HeadKind::Atom(e) => e.class_name().to_string(),
            HeadKind::Derivative => "Derivative".into(),
            HeadKind::Piecewise => "Piecewise".into(),
        }
    }

    #[getter]
    fn __name__(&self) -> String {
        self.__repr__()
    }
}

// ----------------------------------------------------------------------
// module level functions
// ----------------------------------------------------------------------

#[pyfunction]
#[pyo3(name = "Symbol", signature = (name, **assumptions))]
fn py_symbol(name: &str, assumptions: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let Some(d) = assumptions else {
        return Ok(wrap(Expr::symbol(name)));
    };
    let mut kw: Vec<(String, bool)> = Vec::new();
    for (k, v) in d.iter() {
        let k: String = k.extract()?;
        if v.is_none() {
            continue;
        }
        kw.push((k, v.is_truthy()?));
    }
    let a = symox::Assumptions::from_kwargs(kw.iter().map(|(k, v)| (k.as_str(), *v)))
        .map_err(PyValueError::new_err)?;
    Ok(wrap(Expr::symbol_with(name, a)))
}

#[pyfunction]
#[pyo3(name = "Wild", signature = (name, exclude=None, properties=None))]
fn py_wild(name: &str, exclude: Option<&Bound<'_, PyAny>>, properties: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let _ = properties;
    let mut ex = Vec::new();
    if let Some(it) = exclude {
        for x in it.try_iter()? {
            ex.push(to_expr(&x?)?);
        }
    }
    Ok(wrap(Expr::wild(name, ex)))
}

#[pyfunction]
fn sympify(obj: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(to_expr(obj)?))
}

#[pyfunction]
#[pyo3(name = "Integer")]
fn py_integer(obj: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let e = to_expr(obj)?;
    if !e.is_integer_number() {
        return Err(PyValueError::new_err(format!("{e} is not an integer")));
    }
    Ok(wrap(e))
}

#[pyfunction]
#[pyo3(name = "Rational", signature = (p, q=None))]
fn py_rational(p: &Bound<'_, PyAny>, q: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let p = to_expr(p)?;
    let r = match q {
        Some(q) => {
            let q = to_expr(q)?;
            if q.is_zero() {
                return Err(PyZeroDivisionError::new_err("Rational with zero denominator"));
            }
            ops::div(&p, &q)
        }
        None => p,
    };
    Ok(wrap(r))
}

#[pyfunction]
#[pyo3(name = "Float")]
fn py_float(obj: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(evalf::evalf(&to_expr(obj)?)))
}

#[pyfunction]
#[pyo3(name = "srepr")]
fn py_srepr(obj: &Bound<'_, PyAny>) -> PyResult<String> {
    Ok(printer::srepr(&to_expr(obj)?))
}

#[pyfunction]
#[pyo3(name = "diff", signature = (f, *args))]
fn py_diff(f: &Bound<'_, PyAny>, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    Ok(wrap(do_diff(&to_expr(f)?, args)?))
}

#[pyfunction]
#[pyo3(name = "Derivative", signature = (f, *args))]
fn py_derivative(f: &Bound<'_, PyAny>, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let vars = diff_vars(args)?;
    Ok(wrap(diff::derivative_unevaluated(&to_expr(f)?, vars)))
}

/// `integrate(f, x)`, `integrate(f, (x, a, b))`, several variables in turn.
fn do_integrate(f: &Expr, args: &Bound<'_, PyTuple>) -> PyResult<Expr> {
    let mut r = f.clone();
    if args.is_empty() {
        let fs = traverse::free_symbols(f);
        if fs.len() != 1 {
            return Err(PyValueError::new_err("specify the integration variable"));
        }
        return Ok(symox::integrate::integrate(f, &fs[0]));
    }
    for a in args.iter() {
        if let Ok(t) = a.downcast::<PyTuple>() {
            let x = to_expr(&t.get_item(0)?)?;
            match t.len() {
                1 => r = symox::integrate::integrate(&r, &x),
                3 => {
                    let lo = to_expr(&t.get_item(1)?)?;
                    let hi = to_expr(&t.get_item(2)?)?;
                    r = symox::integrate::integrate_definite(&r, &x, &lo, &hi);
                }
                _ => return Err(PyValueError::new_err("expected (x,) or (x, a, b)")),
            }
        } else {
            r = symox::integrate::integrate(&r, &to_expr(&a)?);
        }
    }
    Ok(r)
}

#[pyfunction]
#[pyo3(name = "summation", signature = (f, *limits))]
fn py_summation(f: &Bound<'_, PyAny>, limits: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let mut r = to_expr(f)?;
    for l in limits.iter() {
        let t = l
            .downcast::<PyTuple>()
            .map_err(|_| PyTypeError::new_err("summation expects (k, a, b) tuples"))?;
        if t.len() != 3 {
            return Err(PyValueError::new_err("summation expects (k, a, b) tuples"));
        }
        let (k, a, b) = (to_expr(&t.get_item(0)?)?, to_expr(&t.get_item(1)?)?, to_expr(&t.get_item(2)?)?);
        r = symox::summation::summation(&r, &k, &a, &b);
    }
    Ok(wrap(r))
}

#[pyfunction]
#[pyo3(name = "bernoulli")]
fn py_bernoulli(n: usize) -> PyExpr {
    wrap(symox::summation::bernoulli(n))
}

#[pyfunction]
#[pyo3(name = "integrate", signature = (f, *args, **_kwargs))]
fn py_integrate(f: &Bound<'_, PyAny>, args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(wrap(do_integrate(&to_expr(f)?, args)?))
}

/// `expand` with SymPy hints (`trig`, `log`, `force`, `power_base`, `mul`,
/// `multinomial`, `power_exp`).
fn expand_hints(e: &Expr, hints: Option<&Bound<'_, PyDict>>) -> PyResult<Expr> {
    let get = |k: &str, default: bool| -> PyResult<bool> {
        match hints {
            Some(h) => match h.get_item(k)? {
                Some(v) => v.is_truthy(),
                None => Ok(default),
            },
            None => Ok(default),
        }
    };
    let force = get("force", false)?;
    let mut r = e.clone();
    if get("trig", false)? {
        r = symox::simp::expand_trig(&r);
    }
    if get("log", true)? && (force || get("log", false)?) {
        r = symox::simp::expand_log(&r, force);
    } else if get("log", true)? {
        r = symox::simp::expand_log(&r, false);
    }
    if get("power_base", true)? {
        r = symox::simp::expand_power_base(&r, force);
    }
    let full = get("mul", true)? && get("multinomial", true)?;
    if full {
        r = expand::expand(&r);
    } else if get("mul", true)? {
        r = symox::simp::expand_mul(&r);
    }
    Ok(r)
}

#[pyfunction]
#[pyo3(name = "expand", signature = (e, deep=true, **hints))]
fn py_expand(e: &Bound<'_, PyAny>, deep: bool, hints: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(wrap(expand_hints(&to_expr(e)?, hints)?))
}

#[pyfunction]
#[pyo3(name = "simplify", signature = (e, **_kwargs))]
fn py_simplify(e: &Bound<'_, PyAny>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(wrap(simplify::simplify(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "together")]
fn py_together(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(simplify::together(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "radsimp")]
fn py_radsimp(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(simplify::radsimp(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "N", signature = (e, n=15))]
fn py_n(e: &Bound<'_, PyAny>, n: u32) -> PyResult<PyExpr> {
    let _ = n;
    Ok(wrap(evalf::evalf(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "sqrt")]
fn py_sqrt(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(ops::sqrt(&to_expr(e)?)))
}

#[pyfunction]
#[pyo3(name = "cbrt")]
fn py_cbrt(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(ops::pow(&to_expr(e)?, &Expr::rational(1, 3))))
}

#[pyfunction]
#[pyo3(name = "root")]
fn py_root(e: &Bound<'_, PyAny>, n: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(ops::pow(&to_expr(e)?, &ops::inv(&to_expr(n)?))))
}

#[pyfunction]
#[pyo3(name = "log", signature = (x, base=None))]
fn py_log(x: &Bound<'_, PyAny>, base: Option<&Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let x = to_expr(x)?;
    Ok(wrap(match base {
        Some(b) => fx::log_base(&x, &to_expr(b)?),
        None => fx::log(&x),
    }))
}

/// Apply a known function by name (used by the Python layer).
#[pyfunction]
#[pyo3(name = "_apply", signature = (name, *args))]
fn py_apply(name: &str, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let v: Vec<Expr> = args.iter().map(|a| to_expr(&a)).collect::<PyResult<_>>()?;
    match Func::from_name(name) {
        Some(f) => {
            if let Some(n) = f.arity() {
                if v.len() != n {
                    return Err(PyTypeError::new_err(format!("{name} takes exactly {n} argument(s)")));
                }
            }
            Ok(wrap(fx::func(f, v)))
        }
        None => Ok(wrap(fx::undefined(name, v))),
    }
}

#[pyfunction]
#[pyo3(name = "Rel")]
fn py_rel(op: &str, a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    let op = match op {
        "==" => RelOp::Eq,
        "!=" => RelOp::Ne,
        "<" => RelOp::Lt,
        "<=" => RelOp::Le,
        ">" => RelOp::Gt,
        ">=" => RelOp::Ge,
        _ => return Err(PyValueError::new_err(format!("unknown relational {op}"))),
    };
    Ok(wrap(relational::rel(op, &to_expr(a)?, &to_expr(b)?)))
}

#[pyfunction]
#[pyo3(name = "Piecewise", signature = (*pieces))]
fn py_piecewise(pieces: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
    let mut v = Vec::new();
    for p in pieces.iter() {
        let t = p
            .downcast::<PyTuple>()
            .map_err(|_| PyTypeError::new_err("Piecewise expects (expr, cond) tuples"))?;
        v.push((to_expr(&t.get_item(0)?)?, to_expr(&t.get_item(1)?)?));
    }
    Ok(wrap(relational::piecewise(v)))
}

#[pyfunction]
#[pyo3(name = "free_symbols")]
fn py_free_symbols<'py>(py: Python<'py>, e: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PySet>> {
    PySet::new(py, traverse::free_symbols(&to_expr(e)?).into_iter().map(wrap))
}

#[pyfunction]
#[pyo3(name = "constant")]
fn py_constant(name: &str) -> PyResult<PyExpr> {
    Ok(wrap(match name {
        "pi" => Expr::pi(),
        "E" => Expr::e(),
        "I" => Expr::i(),
        "oo" => Expr::oo(),
        "zoo" => Expr::zoo(),
        "nan" => Expr::nan(),
        "EulerGamma" => Expr::constant(Constant::EulerGamma),
        "true" => Expr::boolean(true),
        "false" => Expr::boolean(false),
        _ => return Err(PyValueError::new_err(format!("unknown constant {name}"))),
    }))
}

// ----------------------------------------------------------------------
// lambdify
// ----------------------------------------------------------------------

/// A compiled numerical function.
#[pyclass(name = "CompiledFunction", module = "symox._core", frozen)]
pub struct PyCompiled {
    c: Compiled,
    scalar_out: bool,
    text: String,
}

#[pymethods]
impl PyCompiled {
    #[pyo3(signature = (*args))]
    fn __call__(&self, py: Python<'_>, args: &Bound<'_, PyTuple>) -> PyResult<PyObject> {
        if args.len() != self.c.n_args() {
            return Err(PyTypeError::new_err(format!(
                "expected {} argument(s), got {}",
                self.c.n_args(),
                args.len()
            )));
        }
        // scalar fast path
        let scalars: Option<Vec<f64>> = args
            .iter()
            .map(|a| if a.is_instance_of::<PyFloat>() || a.is_instance_of::<PyInt>() { a.extract::<f64>().ok() } else { None })
            .collect();
        if let Some(s) = scalars {
            let mut out = vec![0.0; self.c.n_outputs()];
            self.c.call(&s, &mut out);
            return if self.scalar_out {
                Ok(out[0].into_pyobject(py)?.into_any().unbind())
            } else {
                Ok(PyList::new(py, out)?.into_any().unbind())
            };
        }
        // vectorised path: broadcast scalars against 1-d arrays
        let np = py.import("numpy")?;
        let mut arrays: Vec<Vec<f64>> = Vec::with_capacity(args.len());
        let mut n: Option<usize> = None;
        let mut shape: Option<PyObject> = None;
        let mut scalar_flags = Vec::new();
        for a in args.iter() {
            let arr = np.call_method1("asarray", (a, "float64"))?;
            let ndim: usize = arr.getattr("ndim")?.extract()?;
            if ndim == 0 {
                scalar_flags.push(true);
                arrays.push(vec![arr.extract::<f64>()?]);
                continue;
            }
            scalar_flags.push(false);
            let size: usize = arr.getattr("size")?.extract()?;
            if let Some(m) = n {
                if m != size {
                    return Err(PyValueError::new_err("array arguments must have the same size"));
                }
            }
            n = Some(size);
            if shape.is_none() {
                shape = Some(arr.getattr("shape")?.unbind());
            }
            let flat = arr.call_method1("reshape", (-1,))?;
            let flat = np.call_method1("ascontiguousarray", (flat,))?;
            let ro: PyReadonlyArray1<f64> = flat.extract()?;
            arrays.push(ro.as_slice()?.to_vec());
        }
        let n = n.unwrap_or(1);
        for (i, s) in scalar_flags.iter().enumerate() {
            if *s {
                let v = arrays[i][0];
                arrays[i] = vec![v; n];
            }
        }
        let mut outs: Vec<Vec<f64>> = vec![vec![0.0; n]; self.c.n_outputs()];
        {
            let arg_slices: Vec<&[f64]> = arrays.iter().map(|v| v.as_slice()).collect();
            let mut out_slices: Vec<&mut [f64]> = outs.iter_mut().map(|v| v.as_mut_slice()).collect();
            py.allow_threads(|| self.c.call_many(&arg_slices, &mut out_slices));
        }
        let shape = shape.unwrap();
        let mut results = Vec::new();
        for o in outs {
            let a = PyArray1::from_vec(py, o);
            results.push(a.reshape_with_order(
                shape.extract::<Vec<usize>>(py)?,
                numpy::npyffi::NPY_ORDER::NPY_CORDER,
            )?.into_any());
        }
        if self.scalar_out {
            Ok(results.pop().unwrap().unbind())
        } else {
            Ok(PyList::new(py, results)?.into_any().unbind())
        }
    }

    fn __repr__(&self) -> String {
        format!("<compiled {}>", self.text)
    }

    #[getter]
    fn nargs(&self) -> usize {
        self.c.n_args()
    }
}

#[pyfunction]
#[pyo3(name = "lambdify", signature = (args, expr, modules=None, **_kwargs))]
fn py_lambdify(
    args: &Bound<'_, PyAny>,
    expr: &Bound<'_, PyAny>,
    modules: Option<&Bound<'_, PyAny>>,
    _kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<PyCompiled> {
    let _ = modules;
    let arg_list: Vec<Expr> = if let Ok(e) = to_expr(args) {
        vec![e]
    } else {
        args.try_iter()?.map(|a| to_expr(&a?)).collect::<PyResult<_>>()?
    };
    let (exprs, scalar_out): (Vec<Expr>, bool) = if let Ok(e) = to_expr(expr) {
        (vec![e], true)
    } else {
        (expr.try_iter()?.map(|a| to_expr(&a?)).collect::<PyResult<_>>()?, false)
    };
    let c = Compiled::new(&arg_list, &exprs).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let text = exprs.iter().map(|e| e.to_string()).collect::<Vec<_>>().join(", ");
    Ok(PyCompiled { c, scalar_out, text })
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyExpr>()?;
    m.add_class::<Head>()?;
    m.add_class::<PyCompiled>()?;
    m.add_function(wrap_pyfunction!(py_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(sympify, m)?)?;
    m.add_function(wrap_pyfunction!(py_wild, m)?)?;
    m.add_function(wrap_pyfunction!(py_integer, m)?)?;
    m.add_function(wrap_pyfunction!(py_rational, m)?)?;
    m.add_function(wrap_pyfunction!(py_float, m)?)?;
    m.add_function(wrap_pyfunction!(py_srepr, m)?)?;
    m.add_function(wrap_pyfunction!(py_diff, m)?)?;
    m.add_function(wrap_pyfunction!(py_derivative, m)?)?;
    m.add_function(wrap_pyfunction!(py_expand, m)?)?;
    m.add_function(wrap_pyfunction!(py_integrate, m)?)?;
    m.add_function(wrap_pyfunction!(py_summation, m)?)?;
    m.add_function(wrap_pyfunction!(py_bernoulli, m)?)?;
    m.add_function(wrap_pyfunction!(py_simplify, m)?)?;
    m.add_function(wrap_pyfunction!(py_together, m)?)?;
    m.add_function(wrap_pyfunction!(py_radsimp, m)?)?;
    m.add_function(wrap_pyfunction!(py_n, m)?)?;
    m.add_function(wrap_pyfunction!(py_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(py_cbrt, m)?)?;
    m.add_function(wrap_pyfunction!(py_root, m)?)?;
    m.add_function(wrap_pyfunction!(py_log, m)?)?;
    m.add_function(wrap_pyfunction!(py_apply, m)?)?;
    m.add_function(wrap_pyfunction!(py_rel, m)?)?;
    m.add_function(wrap_pyfunction!(py_piecewise, m)?)?;
    m.add_function(wrap_pyfunction!(py_free_symbols, m)?)?;
    m.add_function(wrap_pyfunction!(py_constant, m)?)?;
    m.add_function(wrap_pyfunction!(py_lambdify, m)?)?;
    poly::register(m)?;
    ntheory::register(m)?;
    simp::register(m)?;
    solve::register(m)?;
    ode::register(m)?;
    series::register(m)?;
    printing::register(m)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
