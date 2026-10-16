//! Python bindings of the polynomial module: the `Poly` class and the
//! functions `factor`, `factor_list`, `sqf_list`, `gcd`, `lcm`, `cancel`,
//! `apart`, `collect`, `div`, `rem`, `quo`, `resultant`, `discriminant`,
//! `roots`, `nroots`, `degree`, ...

#![allow(non_snake_case)]

use pyo3::exceptions::{PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyTuple};

use symox::expr::Expr;
use symox::poly::api;
use symox::poly::convert::PolyError;
use symox::poly::polyobj::{Domain, Poly};
use symox::poly::roots as proots;
use symox::poly::Rat;

use crate::{to_expr, PyExpr};

fn w(e: Expr) -> PyExpr {
    PyExpr { e }
}

fn err(e: PolyError) -> PyErr {
    match &e {
        PolyError::Unsupported(s) if s.contains("division by zero") => PyZeroDivisionError::new_err(s.clone()),
        _ => PyValueError::new_err(e.to_string()),
    }
}

fn gens_of(args: &Bound<'_, PyTuple>) -> PyResult<Vec<Expr>> {
    let mut v = Vec::new();
    for a in args.iter() {
        if let Ok(t) = a.downcast::<PyTuple>() {
            for b in t.iter() {
                v.push(to_expr(&b)?);
            }
        } else if let Ok(l) = a.downcast::<PyList>() {
            for b in l.iter() {
                v.push(to_expr(&b)?);
            }
        } else {
            v.push(to_expr(&a)?);
        }
    }
    Ok(v)
}

fn rat_expr(r: &Rat) -> Expr {
    Expr::num(r.to_number())
}

/// A univariate polynomial given as a coefficient list (SymPy accepts
/// `roots([1, 0, -1])`).
fn expr_or_coeff_list(obj: &Bound<'_, PyAny>) -> PyResult<(Expr, Vec<Expr>)> {
    if let Ok(p) = obj.downcast::<PyPoly>() {
        let p = &p.get().p;
        return Ok((p.as_expr(), p.true_gens().to_vec()));
    }
    if let Ok(l) = obj.downcast::<PyList>() {
        let x = Expr::symbol("x");
        let n = l.len();
        let mut terms = Vec::new();
        for (i, c) in l.iter().enumerate() {
            let c = to_expr(&c)?;
            terms.push(symox::ops::mul(&[c, symox::ops::pow(&x, &Expr::int((n - 1 - i) as i64))]));
        }
        return Ok((symox::ops::add(&terms), vec![x]));
    }
    Ok((to_expr(obj)?, vec![]))
}

// ----------------------------------------------------------------------
// Poly
// ----------------------------------------------------------------------

/// A polynomial with explicit generators (SymPy's `Poly`).
#[pyclass(name = "Poly", module = "symox._core", frozen)]
#[derive(Clone)]
pub struct PyPoly {
    pub p: Poly,
}

fn set_domain(mut p: Poly, domain: Option<&str>) -> PyResult<Poly> {
    match domain {
        None => Ok(p),
        Some("ZZ") => {
            if !p.rep.is_integral() {
                return Err(PyValueError::new_err("coefficients are not integers"));
            }
            p.domain = Domain::ZZ;
            Ok(p)
        }
        Some("QQ") => {
            p.domain = Domain::QQ;
            Ok(p)
        }
        Some("RR") => {
            p.domain = Domain::RR;
            Ok(p)
        }
        Some(d) => Err(PyValueError::new_err(format!("unsupported domain {d}"))),
    }
}

impl PyPoly {
    fn other(&self, o: &Bound<'_, PyAny>) -> Option<Poly> {
        if let Ok(q) = o.downcast::<PyPoly>() {
            return Some(q.get().p.clone());
        }
        let e = to_expr(o).ok()?;
        Poly::from_expr(&e, self.p.true_gens()).ok().map(|mut q| {
            if q.domain < self.p.domain {
                q.domain = self.p.domain;
            }
            q
        })
    }

    fn binop(&self, o: &Bound<'_, PyAny>, f: impl Fn(&Poly, &Poly) -> Poly, g: impl Fn(&Expr, &Expr) -> Expr) -> PyResult<PyObject> {
        let py = o.py();
        match self.other(o) {
            Some(q) => Ok(PyPoly { p: f(&self.p, &q) }.into_pyobject(py)?.into_any().unbind()),
            None => {
                let e = to_expr(o)?;
                Ok(w(g(&self.p.as_expr(), &e)).into_pyobject(py)?.into_any().unbind())
            }
        }
    }

    fn gen_index(&self, g: Option<&Bound<'_, PyAny>>) -> PyResult<usize> {
        match g {
            None => Ok(0),
            Some(g) => {
                if let Ok(i) = g.extract::<usize>() {
                    return Ok(i);
                }
                let e = to_expr(g)?;
                self.p
                    .gen_index(&e)
                    .ok_or_else(|| PyValueError::new_err(format!("{} is not a generator", symox::printer::to_str(&e))))
            }
        }
    }
}

fn poly_obj(py: Python<'_>, p: Poly) -> PyResult<PyObject> {
    Ok(PyPoly { p }.into_pyobject(py)?.into_any().unbind())
}

fn degree_obj(py: Python<'_>, d: Option<u32>) -> PyResult<PyObject> {
    match d {
        Some(k) => Ok(k.into_pyobject(py)?.into_any().unbind()),
        None => Ok(w(Expr::neg_oo()).into_pyobject(py)?.into_any().unbind()),
    }
}

#[pymethods]
impl PyPoly {
    #[new]
    #[pyo3(signature = (expr, *gens, domain=None, **_kwargs))]
    fn py_new(expr: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, domain: Option<&str>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let gens = gens_of(gens)?;
        if let Ok(q) = expr.downcast::<PyPoly>() {
            let q = &q.get().p;
            if gens.is_empty() {
                return Ok(PyPoly { p: set_domain(q.clone(), domain)? });
            }
            let p = Poly::from_expr(&q.as_expr(), &gens).map_err(err)?;
            return Ok(PyPoly { p: set_domain(p, domain)? });
        }
        let (e, g2) = expr_or_coeff_list(expr)?;
        let gens = if gens.is_empty() { g2 } else { gens };
        let p = Poly::from_expr(&e, &gens).map_err(err)?;
        Ok(PyPoly { p: set_domain(p, domain)? })
    }

    fn __repr__(&self) -> String {
        self.p.to_sympy_string()
    }

    fn __str__(&self) -> String {
        self.p.to_sympy_string()
    }

    fn _symox_(&self) -> PyExpr {
        w(self.p.as_expr())
    }

    fn as_expr(&self) -> PyExpr {
        w(self.p.as_expr())
    }

    #[getter]
    fn expr(&self) -> PyExpr {
        w(self.p.as_expr())
    }

    #[getter]
    fn gens<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(py, self.p.true_gens().iter().map(|g| w(g.clone())))
    }

    #[getter]
    fn gen(&self) -> PyExpr {
        w(self.p.true_gens()[0].clone())
    }

    #[getter]
    fn domain(&self) -> String {
        self.p.domain_str()
    }

    #[getter]
    fn free_symbols<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PySet>> {
        let fs = symox::traverse::free_symbols(&self.p.as_expr());
        let mut v: Vec<PyExpr> = Vec::new();
        for g in self.p.true_gens() {
            for s in symox::traverse::free_symbols(g) {
                if !v.iter().any(|x| x.e == s) {
                    v.push(w(s));
                }
            }
        }
        for s in fs {
            if !v.iter().any(|x| x.e == s) {
                v.push(w(s));
            }
        }
        pyo3::types::PySet::new(py, v)
    }

    #[pyo3(signature = (gen=None))]
    fn degree(&self, py: Python<'_>, gen: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
        let i = self.gen_index(gen)?;
        degree_obj(py, self.p.degree(i))
    }

    fn degree_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let v: Vec<i64> = (0..self.p.ngens).map(|i| self.p.rep.degree(i)).collect();
        PyTuple::new(py, v)
    }

    fn total_degree(&self, py: Python<'_>) -> PyResult<PyObject> {
        degree_obj(py, self.p.total_degree())
    }

    fn coeffs(&self) -> Vec<PyExpr> {
        self.p.coeffs().into_iter().map(w).collect()
    }

    fn all_coeffs(&self) -> PyResult<Vec<PyExpr>> {
        Ok(self.p.all_coeffs().map_err(err)?.into_iter().map(w).collect())
    }

    fn monoms<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyTuple>>> {
        self.p.monoms().into_iter().map(|m| PyTuple::new(py, m)).collect()
    }

    fn terms<'py>(&self, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyTuple>>> {
        let t = self.p.terms();
        let t = if t.is_empty() { vec![(vec![0; self.p.ngens], Expr::zero())] } else { t };
        t.into_iter()
            .map(|(m, c)| {
                let mt = PyTuple::new(py, m)?;
                PyTuple::new(py, [mt.into_any(), w(c).into_pyobject(py)?.into_any()])
            })
            .collect()
    }

    fn as_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        for (m, c) in self.p.terms() {
            d.set_item(PyTuple::new(py, m)?, w(c))?;
        }
        Ok(d)
    }

    fn LC(&self) -> PyExpr {
        w(self.p.lc())
    }

    fn TC(&self) -> PyExpr {
        w(self.p.tc())
    }

    fn LM<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let m = self.p.monoms().into_iter().next().unwrap_or_default();
        PyTuple::new(py, m)
    }

    #[pyo3(signature = (*n))]
    fn nth(&self, n: Vec<u32>) -> PyExpr {
        w(self.p.coeff_monomial(&n))
    }

    fn length(&self) -> usize {
        self.p.terms().len()
    }

    #[getter]
    fn is_zero(&self) -> bool {
        self.p.is_zero()
    }

    #[getter]
    fn is_one(&self) -> bool {
        self.p.rep.is_one()
    }

    #[getter]
    fn is_ground(&self) -> bool {
        self.p.rep.iter().all(|(e, _)| e[..self.p.ngens].iter().all(|&k| k == 0))
    }

    #[getter]
    fn is_monomial(&self) -> bool {
        self.p.terms().len() <= 1
    }

    #[getter]
    fn is_univariate(&self) -> bool {
        self.p.ngens == 1
    }

    #[getter]
    fn is_multivariate(&self) -> bool {
        self.p.ngens > 1
    }

    #[getter]
    fn is_monic(&self) -> bool {
        self.p.lc().is_one()
    }

    #[getter]
    fn is_irreducible(&self) -> bool {
        let (_, fs) = self.p.factor_list();
        fs.len() == 1 && fs[0].1 == 1
    }

    #[getter]
    fn is_sqf(&self) -> bool {
        let (_, fs) = self.p.sqf_list();
        fs.iter().all(|(_, k)| *k == 1)
    }

    fn __add__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.add(b), |a, b| symox::ops::add2(a, b))
    }
    fn __radd__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.add(a), |a, b| symox::ops::add2(b, a))
    }
    fn __sub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.sub(b), |a, b| symox::ops::sub(a, b))
    }
    fn __rsub__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.sub(a), |a, b| symox::ops::sub(b, a))
    }
    fn __mul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| a.mul(b), |a, b| symox::ops::mul2(a, b))
    }
    fn __rmul__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyObject> {
        self.binop(o, |a, b| b.mul(a), |a, b| symox::ops::mul2(b, a))
    }
    fn __neg__(&self) -> PyPoly {
        PyPoly { p: self.p.neg() }
    }
    fn __pos__(&self) -> PyPoly {
        self.clone()
    }
    fn __pow__(&self, k: u32, _m: Option<&Bound<'_, PyAny>>) -> PyPoly {
        PyPoly { p: self.p.pow(k) }
    }
    fn __truediv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        Ok(w(symox::ops::div(&self.p.as_expr(), &to_expr(o)?)))
    }
    fn __floordiv__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: self.p.quo(&q).map_err(err)? })
    }
    fn __mod__(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: self.p.rem(&q).map_err(err)? })
    }
    fn __divmod__<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        self.div(py, o)
    }

    fn __eq__(&self, o: &Bound<'_, PyAny>) -> bool {
        if let Ok(q) = o.downcast::<PyPoly>() {
            let q = &q.get().p;
            return self.p.true_gens() == q.true_gens() && self.p.as_expr() == q.as_expr();
        }
        match self.other(o) {
            Some(q) => self.p.as_expr() == q.as_expr(),
            None => false,
        }
    }

    fn __ne__(&self, o: &Bound<'_, PyAny>) -> bool {
        !self.__eq__(o)
    }

    fn __hash__(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        self.p.as_expr().hash(&mut h);
        h.finish()
    }

    fn __bool__(&self) -> bool {
        !self.p.is_zero()
    }

    fn div<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        let (a, b) = self.p.div(&q).map_err(err)?;
        PyTuple::new(py, [PyPoly { p: a }, PyPoly { p: b }])
    }

    fn rem(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        self.__mod__(o)
    }

    fn quo(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        self.__floordiv__(o)
    }

    fn exquo(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: self.p.exquo(&q).map_err(err)? })
    }

    fn gcd(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: self.p.gcd(&q) })
    }

    fn lcm(&self, o: &Bound<'_, PyAny>) -> PyResult<PyPoly> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(PyPoly { p: self.p.lcm(&q) })
    }

    fn cofactors<'py>(&self, py: Python<'py>, o: &Bound<'_, PyAny>) -> PyResult<Bound<'py, PyTuple>> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        let (h, a, b) = self.p.cofactors(&q);
        PyTuple::new(py, [PyPoly { p: h }, PyPoly { p: a }, PyPoly { p: b }])
    }

    fn factor_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (c, fs) = self.p.factor_list();
        let l = PyList::empty(py);
        for (f, k) in fs {
            l.append(PyTuple::new(py, [poly_obj(py, f)?, k.into_pyobject(py)?.into_any().unbind()])?)?;
        }
        PyTuple::new(py, [w(rat_expr(&c)).into_pyobject(py)?.into_any(), l.into_any()])
    }

    fn sqf_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (c, fs) = self.p.sqf_list();
        let l = PyList::empty(py);
        for (f, k) in fs {
            l.append(PyTuple::new(py, [poly_obj(py, f)?, k.into_pyobject(py)?.into_any().unbind()])?)?;
        }
        PyTuple::new(py, [w(rat_expr(&c)).into_pyobject(py)?.into_any(), l.into_any()])
    }

    fn primitive<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let (c, p) = self.p.primitive();
        PyTuple::new(py, [w(rat_expr(&c)).into_pyobject(py)?.into_any(), PyPoly { p }.into_pyobject(py)?.into_any()])
    }

    fn content(&self) -> PyExpr {
        w(rat_expr(&self.p.primitive().0))
    }

    fn monic(&self) -> PyResult<PyPoly> {
        if self.p.is_zero() {
            return Ok(self.clone());
        }
        let lc = self.p.rep.lc();
        let mut p = self.p.mul_ground(&lc.inv());
        if p.domain == Domain::ZZ && !lc.is_one() {
            p.domain = Domain::QQ;
        }
        Ok(PyPoly { p })
    }

    #[pyo3(signature = (*gens))]
    fn diff(&self, gens: &Bound<'_, PyTuple>) -> PyResult<PyPoly> {
        let mut p = self.p.clone();
        if gens.is_empty() {
            return Ok(PyPoly { p: p.diff(0) });
        }
        for g in gens.iter() {
            let i = self.gen_index(Some(&g))?;
            p = p.diff(i);
        }
        Ok(PyPoly { p })
    }

    #[pyo3(signature = (x, a=None))]
    fn eval(&self, py: Python<'_>, x: &Bound<'_, PyAny>, a: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
        let (i, val) = match a {
            Some(a) => (self.gen_index(Some(x))?, to_expr(a)?),
            None => (0, to_expr(x)?),
        };
        let g = self.p.true_gens()[i].clone();
        let r = symox::subs::subs(&self.p.as_expr(), &[(g, val)]);
        let r = symox::expand::expand(&r);
        if self.p.ngens > 1 {
            let mut gens = self.p.true_gens().to_vec();
            gens.remove(i);
            if let Ok(p) = Poly::from_expr(&r, &gens) {
                return poly_obj(py, p);
            }
        }
        Ok(w(r).into_pyobject(py)?.into_any().unbind())
    }

    #[pyo3(signature = (*args))]
    fn __call__(&self, args: &Bound<'_, PyTuple>) -> PyResult<PyExpr> {
        let mut pairs = Vec::new();
        for (i, a) in args.iter().enumerate() {
            pairs.push((self.p.true_gens()[i].clone(), to_expr(&a)?));
        }
        Ok(w(symox::expand::expand(&symox::subs::subs(&self.p.as_expr(), &pairs))))
    }

    #[pyo3(signature = (n=15, **_kwargs))]
    fn nroots(&self, n: u32, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
        Ok(proots::nroots_poly(&self.p, n).map_err(err)?.into_iter().map(w).collect())
    }

    fn all_roots(&self) -> PyResult<Vec<PyExpr>> {
        let rs = proots::roots_poly(&self.p).map_err(err)?;
        let mut v = Vec::new();
        for (r, k) in rs {
            for _ in 0..k {
                v.push(w(r.clone()));
            }
        }
        Ok(v)
    }

    fn as_poly(&self) -> PyPoly {
        self.clone()
    }

    fn resultant(&self, o: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
        let q = self.other(o).ok_or_else(|| PyTypeError::new_err("expected a polynomial"))?;
        Ok(w(api::resultant(&self.p.as_expr(), &q.as_expr(), self.p.true_gens()).map_err(err)?))
    }

    fn discriminant(&self) -> PyResult<PyExpr> {
        Ok(w(api::discriminant(&self.p.as_expr(), self.p.true_gens()).map_err(err)?))
    }
}

// ----------------------------------------------------------------------
// functions
// ----------------------------------------------------------------------

fn pair<'py>(py: Python<'py>, c: Expr, fs: Vec<(Expr, Expr)>) -> PyResult<Bound<'py, PyTuple>> {
    let l = PyList::empty(py);
    for (f, k) in fs {
        let kk: PyObject = match k.as_i64() {
            Some(i) => i.into_pyobject(py)?.into_any().unbind(),
            None => w(k).into_pyobject(py)?.into_any().unbind(),
        };
        l.append(PyTuple::new(py, [w(f).into_pyobject(py)?.into_any().unbind(), kk])?)?;
    }
    PyTuple::new(py, [w(c).into_pyobject(py)?.into_any(), l.into_any()])
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn factor(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = gens;
    let e = to_expr(f)?;
    Ok(w(api::factor(&e).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn sqf(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = gens;
    Ok(w(api::sqf(&to_expr(f)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn factor_list<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    if let Ok(p) = f.downcast::<PyPoly>() {
        return p.get().factor_list(py);
    }
    let gens = gens_of(gens)?;
    let e = to_expr(f)?;
    if !gens.is_empty() {
        let (c, fs) = api::poly_factor_list(&e, &gens).map_err(err)?;
        return pair(py, c, fs.into_iter().map(|(f, k)| (f, Expr::int(k as i64))).collect());
    }
    let (c, fs) = api::factor_list(&e).map_err(err)?;
    pair(py, c, fs)
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn sqf_list<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    if let Ok(p) = f.downcast::<PyPoly>() {
        return p.get().sqf_list(py);
    }
    let _ = gens;
    let (c, fs) = api::sqf_list(&to_expr(f)?).map_err(err)?;
    pair(py, c, fs)
}

fn list_items(f: &Bound<'_, PyAny>) -> Option<Vec<Expr>> {
    if f.downcast::<PyList>().is_ok() || f.downcast::<PyTuple>().is_ok() {
        let mut v = Vec::new();
        for x in f.try_iter().ok()? {
            v.push(to_expr(&x.ok()?).ok()?);
        }
        return Some(v);
    }
    None
}

fn poly_result(py: Python<'_>, e: Expr, as_poly: bool, gens: &[Expr]) -> PyResult<PyObject> {
    if as_poly {
        let p = Poly::from_expr(&e, gens).map_err(err)?;
        return poly_obj(py, p);
    }
    Ok(w(e).into_pyobject(py)?.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, g=None, *gens, **_kwargs))]
fn gcd(py: Python<'_>, f: &Bound<'_, PyAny>, g: Option<&Bound<'_, PyAny>>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let as_poly = f.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if let Some(items) = list_items(f) {
        if let Some(g) = g {
            gv.insert(0, to_expr(g)?);
        }
        let mut acc = Expr::zero();
        for x in items {
            acc = api::gcd(&acc, &x, &gv).map_err(err)?;
        }
        return Ok(w(acc).into_pyobject(py)?.into_any().unbind());
    }
    let g = g.ok_or_else(|| PyTypeError::new_err("gcd() takes 2 arguments or a sequence"))?;
    let a = to_expr(f)?;
    let b = to_expr(g)?;
    let r = api::gcd(&a, &b, &gv).map_err(err)?;
    let gens = if as_poly { f.downcast::<PyPoly>().unwrap().get().p.true_gens().to_vec() } else { gv };
    poly_result(py, r, as_poly, &gens)
}

#[pyfunction]
#[pyo3(signature = (f, g=None, *gens, **_kwargs))]
fn lcm(py: Python<'_>, f: &Bound<'_, PyAny>, g: Option<&Bound<'_, PyAny>>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let as_poly = f.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if let Some(items) = list_items(f) {
        if let Some(g) = g {
            gv.insert(0, to_expr(g)?);
        }
        let mut acc = Expr::one();
        for x in items {
            acc = api::lcm(&acc, &x, &gv).map_err(err)?;
        }
        return Ok(w(acc).into_pyobject(py)?.into_any().unbind());
    }
    let g = g.ok_or_else(|| PyTypeError::new_err("lcm() takes 2 arguments or a sequence"))?;
    let r = api::lcm(&to_expr(f)?, &to_expr(g)?, &gv).map_err(err)?;
    let gens = if as_poly { f.downcast::<PyPoly>().unwrap().get().p.true_gens().to_vec() } else { gv };
    poly_result(py, r, as_poly, &gens)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn cofactors<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    let (h, a, b) = api::cofactors(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?;
    PyTuple::new(py, [w(h), w(a), w(b)])
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn cancel(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = gens;
    Ok(w(api::cancel(&to_expr(f)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, x=None, full=false, **_kwargs))]
fn apart(f: &Bound<'_, PyAny>, x: Option<&Bound<'_, PyAny>>, full: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = full;
    let x = match x {
        Some(x) => Some(to_expr(x)?),
        None => None,
    };
    Ok(w(api::apart(&to_expr(f)?, x.as_ref()).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (expr, syms, func=None, evaluate=true, **_kwargs))]
fn collect(expr: &Bound<'_, PyAny>, syms: &Bound<'_, PyAny>, func: Option<&Bound<'_, PyAny>>, evaluate: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let _ = evaluate;
    let py = expr.py();
    let e = to_expr(expr)?;
    let syms = match list_items(syms) {
        Some(v) => v,
        None => vec![to_expr(syms)?],
    };
    let r = api::collect(&e, &syms);
    match func {
        Some(f) => {
            // apply func to each collected coefficient
            let parts: Vec<Expr> = match r.kind() {
                symox::Kind::Add(..) => r.args(),
                _ => vec![r.clone()],
            };
            let mut out = Vec::new();
            for t in parts {
                out.push(to_expr(&f.call1((w(t),))?)?);
            }
            Ok(w(symox::ops::add(&out)).into_pyobject(py)?.into_any().unbind())
        }
        None => Ok(w(r).into_pyobject(py)?.into_any().unbind()),
    }
}

fn div_like(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>) -> PyResult<(Expr, Expr, Vec<Expr>, bool)> {
    let as_poly = f.downcast::<PyPoly>().is_ok() || g.downcast::<PyPoly>().is_ok();
    let mut gv = gens_of(gens)?;
    if gv.is_empty() {
        if let Ok(p) = f.downcast::<PyPoly>() {
            gv = p.get().p.true_gens().to_vec();
        }
    }
    let (q, r) = api::pdiv(&to_expr(f)?, &to_expr(g)?, &gv).map_err(err)?;
    Ok((q, r, gv, as_poly))
}

#[pyfunction]
#[pyo3(name = "div", signature = (f, g, *gens, **_kwargs))]
fn py_div<'py>(py: Python<'py>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Bound<'py, PyTuple>> {
    if let (Ok(a), Ok(_)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
        return a.get().div(py, g);
    }
    let (q, r, gv, as_poly) = div_like(f, g, gens)?;
    PyTuple::new(py, [poly_result(py, q, as_poly, &gv)?, poly_result(py, r, as_poly, &gv)?])
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn rem(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if let (Ok(a), Ok(_)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
        return Ok(a.get().rem(g)?.into_pyobject(py)?.into_any().unbind());
    }
    let (_, r, gv, as_poly) = div_like(f, g, gens)?;
    poly_result(py, r, as_poly, &gv)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn quo(py: Python<'_>, f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    if let (Ok(a), Ok(_)) = (f.downcast::<PyPoly>(), g.downcast::<PyPoly>()) {
        return Ok(a.get().quo(g)?.into_pyobject(py)?.into_any().unbind());
    }
    let (q, _, gv, as_poly) = div_like(f, g, gens)?;
    poly_result(py, q, as_poly, &gv)
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn exquo(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::exquo(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, g, *gens, **_kwargs))]
fn resultant(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::resultant(&to_expr(f)?, &to_expr(g)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn discriminant(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    Ok(w(api::discriminant(&to_expr(f)?, &gens_of(gens)?).map_err(err)?))
}

#[pyfunction]
#[pyo3(signature = (f, gen=None))]
fn degree(py: Python<'_>, f: &Bound<'_, PyAny>, gen: Option<&Bound<'_, PyAny>>) -> PyResult<PyObject> {
    if let Ok(p) = f.downcast::<PyPoly>() {
        let p = p.get();
        let i = p.gen_index(gen)?;
        return match p.p.degree(i) {
            Some(k) => Ok(w(Expr::int(k as i64)).into_pyobject(py)?.into_any().unbind()),
            None => Ok(w(Expr::neg_oo()).into_pyobject(py)?.into_any().unbind()),
        };
    }
    let e = to_expr(f)?;
    let g = match gen {
        None => None,
        Some(g) => {
            if g.extract::<i64>().is_ok() {
                None
            } else {
                Some(to_expr(g)?)
            }
        }
    };
    let d = if let Some(g) = &g {
        if !symox::traverse::has(&e, g) {
            Some(if e.is_zero() { None } else { Some(0) })
        } else {
            None
        }
    } else {
        None
    };
    let d = match d {
        Some(d) => d,
        None => api::degree(&e, g.as_ref()).map_err(err)?,
    };
    Ok(w(match d {
        Some(k) => Expr::int(k as i64),
        None => Expr::neg_oo(),
    })
    .into_pyobject(py)?
    .into_any()
    .unbind())
}

#[pyfunction]
#[pyo3(signature = (f, *gens, multiple=false, filter=None, **_kwargs))]
fn roots(py: Python<'_>, f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, multiple: bool, filter: Option<&str>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyObject> {
    let (e, g2) = expr_or_coeff_list(f)?;
    let gv = gens_of(gens)?;
    let gv = if gv.is_empty() { g2 } else { gv };
    let mut rs = proots::roots(&e, &gv).map_err(err)?;
    if let Some(flt) = filter {
        rs.retain(|(r, _)| match flt {
            "Z" => r.is_integer_number(),
            "Q" => r.is_rational_number(),
            "R" => r.mask() == 0 && symox::evalf::eval_c64(r).is_some_and(|z| z.im == 0.0) && !symox::traverse::has(r, &Expr::i()),
            "I" => r.mask() == 0 && symox::evalf::eval_c64(r).is_some_and(|z| z.re.abs() < 1e-300 && z.im != 0.0),
            _ => true,
        });
    }
    if multiple {
        let l = PyList::empty(py);
        for (r, k) in rs {
            for _ in 0..k {
                l.append(w(r.clone()))?;
            }
        }
        return Ok(l.into_any().unbind());
    }
    let d = PyDict::new(py);
    for (r, k) in rs {
        d.set_item(w(r), k)?;
    }
    Ok(d.into_any().unbind())
}

#[pyfunction]
#[pyo3(signature = (f, n=15, **_kwargs))]
fn nroots(f: &Bound<'_, PyAny>, n: u32, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
    let (e, g2) = expr_or_coeff_list(f)?;
    Ok(proots::nroots(&e, &g2, n).map_err(err)?.into_iter().map(w).collect())
}

#[pyfunction]
#[pyo3(signature = (f, *gens, **_kwargs))]
fn real_roots(f: &Bound<'_, PyAny>, gens: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<PyExpr>> {
    let (e, g2) = expr_or_coeff_list(f)?;
    let gv = gens_of(gens)?;
    let gv = if gv.is_empty() { g2 } else { gv };
    let rs = proots::roots(&e, &gv).map_err(err)?;
    let mut out: Vec<(f64, Expr, u32)> = Vec::new();
    for (r, k) in rs {
        if r.mask() == 0 && !symox::traverse::has(&r, &Expr::i()) {
            if let Some(v) = symox::evalf::eval_f64(&r) {
                out.push((v, r, k));
            }
        }
    }
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut v = Vec::new();
    for (_, r, k) in out {
        for _ in 0..k {
            v.push(w(r.clone()));
        }
    }
    Ok(v)
}

#[pyfunction]
#[pyo3(name = "_poly_together", signature = (f, deep=false, **_kwargs))]
fn py_together(f: &Bound<'_, PyAny>, deep: bool, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<PyExpr> {
    let _ = deep;
    Ok(w(api::together(&to_expr(f)?)))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPoly>()?;
    m.add_function(wrap_pyfunction!(factor, m)?)?;
    m.add_function(wrap_pyfunction!(sqf, m)?)?;
    m.add_function(wrap_pyfunction!(factor_list, m)?)?;
    m.add_function(wrap_pyfunction!(sqf_list, m)?)?;
    m.add_function(wrap_pyfunction!(gcd, m)?)?;
    m.add_function(wrap_pyfunction!(lcm, m)?)?;
    m.add_function(wrap_pyfunction!(cofactors, m)?)?;
    m.add_function(wrap_pyfunction!(cancel, m)?)?;
    m.add_function(wrap_pyfunction!(apart, m)?)?;
    m.add_function(wrap_pyfunction!(collect, m)?)?;
    m.add_function(wrap_pyfunction!(py_div, m)?)?;
    m.add_function(wrap_pyfunction!(rem, m)?)?;
    m.add_function(wrap_pyfunction!(quo, m)?)?;
    m.add_function(wrap_pyfunction!(exquo, m)?)?;
    m.add_function(wrap_pyfunction!(resultant, m)?)?;
    m.add_function(wrap_pyfunction!(discriminant, m)?)?;
    m.add_function(wrap_pyfunction!(degree, m)?)?;
    m.add_function(wrap_pyfunction!(roots, m)?)?;
    m.add_function(wrap_pyfunction!(nroots, m)?)?;
    m.add_function(wrap_pyfunction!(real_roots, m)?)?;
    m.add_function(wrap_pyfunction!(py_together, m)?)?;
    Ok(())
}
