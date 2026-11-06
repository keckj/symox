//! The reference `RootSum` constructor (`sympy.polys.rootoftools`):
//! `_transform` (`preprocess_roots`), `__new__` after the argument checks,
//! `new`, `_rational_case`, `doit` and `_eval_evalf`.
//!
//! The Python class keeps the argument conversion (`PurePoly(expr, x,
//! greedy=False)`, the univariate and function checks) and the side table
//! of the polynomials of the raw nodes. Every computing step is here, on
//! the native Poly methods of the reference calls (`clear_denoms`,
//! `primitive`, the `retract`, `inject`, `eval` and `eject` constructions,
//! `factor_list`, `all_coeffs`, `nroots`, the native `roots`) and on the
//! core closed forms, symmetric reduction and `factor`.
//!
//! A function returns `None` (the Python body then runs) only for inputs
//! whose reference steps leave the native representation: the decline is
//! never a different answer.
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};

use super::{err, w, PyPoly};
use crate::{to_expr, PyExpr};
use symox::expr::Expr;
use symox::ops::{add, mul};
use symox::poly::roots::{self as proots, RootDomain};

fn is_error(py: Python<'_>, e: &PyErr, kind: &str) -> PyResult<bool> {
    let class = py.import("symox.polys.polyerrors")?.getattr(kind)?;
    Ok(e.value(py).is_instance(&class)?)
}

fn domain_flag(py: Python<'_>, p: &PyPoly, flag: &str) -> PyResult<bool> {
    p.domain(py)?.bind(py).getattr(flag)?.is_truthy()
}

fn root_domain(py: Python<'_>, p: &PyPoly) -> PyResult<RootDomain> {
    Ok(if domain_flag(py, p, "is_Numerical")? {
        RootDomain::Numerical
    } else if domain_flag(py, p, "is_Composite")? {
        RootDomain::Composite
    } else {
        RootDomain::Expression
    })
}

fn degree(p: &PyPoly) -> i64 { p.p.degree(0).map_or(-1, i64::from) }

fn poly_expr(p: &PyPoly) -> Expr { symox::poly::structure::sympy_expr(&p.p) }

fn construct(py: Python<'_>, e: &Expr, gens: &[Expr], domain: Option<&Bound<'_, PyAny>>) -> PyResult<PyPoly> {
    let rep = w(e.clone()).into_pyobject(py)?.into_any();
    let gens = PyTuple::new(py, gens.iter().map(|g| w(g.clone())))?;
    PyPoly::_from_normalized(&rep, &gens, domain, None, None)
}

/// `Poly.retract()`: the polynomial rebuilt from its expression.
fn retract(py: Python<'_>, p: &PyPoly) -> PyResult<PyPoly> {
    if p.p.base.modulus.is_some() { return Ok(p.clone()); }
    construct(py, &poly_expr(p), p.p.true_gens(), None)
}

fn integer(e: &Expr) -> Option<num_bigint::BigInt> {
    let n = e.as_number()?;
    if !n.is_integer() { return None; }
    n.to_bigint()
}

fn big(n: num_bigint::BigInt) -> Expr { Expr::num(symox::number::Number::from_bigint(n)) }

/// The reference `preprocess_roots(poly)`: `(coeff, poly)`, `None` when a
/// step leaves the native representation.
fn preprocess(py: Python<'_>, poly: &PyPoly) -> PyResult<Option<(Expr, PyPoly)>> {
    let mut coeff = Expr::one();
    let cleared = match poly.clear_denoms(py, true) {
        Ok(pair) => pair,
        Err(e) if is_error(py, &e, "DomainError")? => return Ok(Some((coeff, poly.clone()))),
        Err(e) => return Err(e),
    };
    let Ok(p) = cleared.get_item(1)?.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
    let Ok(p) = p.primitive(py)?.get_item(1)?.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
    let mut p = retract(py, &p)?;
    let domain = p.domain(py)?.into_bound(py);
    if domain.getattr("is_Poly")?.is_truthy()? {
        if p.p.true_gens().len() != 1 { return Ok(None); }
        let x = p.p.true_gens()[0].clone();
        let symbols: Vec<Expr> = domain.getattr("symbols")?.try_iter()?.map(|s| to_expr(&s?)).collect::<PyResult<_>>()?;
        let ground = domain.getattr("dom")?;
        let coeffs = p.p.all_coeffs().map_err(err)?;
        let n = coeffs.len() - 1;
        // `all(c.is_term for c in poly.rep.coeffs())`: every nonzero
        // coefficient is a single term of the ring.
        let mut rows: Vec<(u32, Vec<u32>)> = Vec::new();
        let mut terms = true;
        for (i, c) in coeffs.iter().enumerate() {
            if c.is_zero() { continue; }
            if c.is_add() { terms = false; break; }
            let Ok(converted) = symox::poly::convert::from_exprs_gens(std::slice::from_ref(c), &symbols) else { return Ok(None) };
            if converted.polys[0].len() != 1 || converted.float { return Ok(None); }
            let (e, _) = converted.polys[0].iter().next().unwrap();
            rows.push(((n - i) as u32, e[..symbols.len()].to_vec()));
        }
        if terms {
            // `poly.inject()`: the monomials over `(x, *symbols)` in lex
            // order, one per degree of `x`, highest first.
            let base: Vec<u32> = rows.iter().map(|r| r.0).collect();
            let mut remaining: Vec<usize> = (0..symbols.len()).collect();
            let mut evaluated: Vec<(Expr, Expr)> = Vec::new();
            for j in 0..symbols.len() {
                let mut strip: Vec<u32> = rows.iter().map(|r| r.1[j]).collect();
                let reverse = strip[0] < strip[strip.len() - 1];
                if reverse { strip.reverse(); }
                let mut ratio: Option<i64> = None;
                let mut complete = true;
                for (&a, &b) in base.iter().zip(&strip) {
                    if a == 0 && b == 0 { continue; }
                    if a == 0 || b == 0 || b % a != 0 { complete = false; break; }
                    let r = (b / a) as i64;
                    if ratio.is_some_and(|q| q != r) { complete = false; break; }
                    ratio = Some(r);
                }
                if !complete { continue; }
                // `ratio` stays `None` only for an all-zero strip, where the
                // reference fails on `-None`.
                let Some(mut ratio) = ratio else { return Ok(None) };
                if reverse { ratio = -ratio; }
                coeff = mul(&[coeff, symox::ops::pow(&symbols[j], &Expr::int(-ratio))]);
                evaluated.push((symbols[j].clone(), Expr::one()));
                remaining.retain(|&k| k != j);
            }
            let e = symox::subs::subs(&poly_expr(&p), &evaluated);
            let kept: Vec<Expr> = remaining.iter().map(|&k| symbols[k].clone()).collect();
            let target = if kept.is_empty() {
                ground.clone()
            } else {
                ground.call_method1("poly_ring", PyTuple::new(py, kept.iter().map(|g| w(g.clone())))?)?
            };
            p = construct(py, &e, std::slice::from_ref(&x), Some(&target))?;
        }
    }
    let domain = p.domain(py)?.into_bound(py);
    if p.p.true_gens().len() == 1 && domain.getattr("is_ZZ")?.is_truthy()? {
        let coeffs = p.p.all_coeffs().map_err(err)?;
        let n = coeffs.len() - 1;
        let mut terms = Vec::new();
        for (i, c) in coeffs.iter().enumerate() {
            if c.is_zero() { continue; }
            let Some(v) = integer(c) else { return Ok(None) };
            terms.push(((n - i) as u32, v));
        }
        if let Some(basis) = proots::integer_basis_of(&terms) {
            let parts: Vec<Expr> = terms.iter().map(|(k, c)| {
                let scaled = num_integer_floor(c, &num_traits::pow(basis.clone(), n - *k as usize));
                mul(&[big(scaled), symox::ops::pow(&p.p.true_gens()[0], &Expr::int(*k as i64))])
            }).collect();
            p = construct(py, &add(&parts), p.p.true_gens(), Some(&domain))?;
            coeff = mul(&[coeff, big(basis)]);
        }
    }
    Ok(Some((coeff, p)))
}

/// Python's `c // d` on integers.
fn num_integer_floor(c: &num_bigint::BigInt, d: &num_bigint::BigInt) -> num_bigint::BigInt {
    use num_traits::{Signed, Zero};
    let (q, r) = (c / d, c % d);
    if !r.is_zero() && (r.is_negative() != d.is_negative()) { q - 1 } else { q }
}

/// `RootSum._transform(expr, x)` after `PurePoly(expr, x, greedy=False)`.
#[pyfunction]
pub(crate) fn _root_sum_transform(py: Python<'_>, poly: &Bound<'_, PyAny>) -> PyResult<Option<(PyExpr, PyPoly)>> {
    let Ok(p) = poly.downcast::<PyPoly>() else { return Ok(None) };
    Ok(preprocess(py, p.get())?.map(|(c, p)| (w(c), p)))
}

/// `Lambda(var, expr)(value)`.
fn apply(var: &Expr, body: &Expr, value: &Expr) -> Expr {
    let mut map = rustc_hash::FxHashMap::default();
    map.insert(var.clone(), value.clone());
    symox::subs::xreplace(body, &map)
}

fn is_rational(expr: &Expr, var: &Expr) -> bool {
    symox::inspection::classify(expr, std::slice::from_ref(var), symox::inspection::PolynomialKind::Rational) == Some(true)
}

fn public_expand(e: &Expr) -> Option<Expr> { crate::expand_hints(e, true, None).ok() }

/// `RootSum._rational_case(poly, Lambda(var, expr))`.
fn rational_case(poly: &PyPoly, var: &Expr, expr: &Expr) -> PyResult<Option<Expr>> {
    let coeffs = poly.p.all_coeffs().map_err(err)?;
    if coeffs.len() < 2 { return Ok(None); }
    let lc = coeffs[0].clone();
    let mut sign = -1i64;
    let mut viete = Vec::new();
    for c in &coeffs[1..] {
        viete.push(mul(&[Expr::int(sign), symox::ops::div(c, &lc)]));
        sign = -sign;
    }
    Ok(symox::poly::rootsum::rational_case(expr, var, &viete, &public_expand))
}

type Node = (PyExpr, PyPoly, bool);

/// `RootSum._new(poly, func, auto)`: the raw node, recorded for the side
/// table of the Python class.
fn raw(poly: &PyPoly, var: &Expr, expr: &Expr, auto: bool, nodes: &mut Vec<Node>) -> PyResult<Option<Expr>> {
    let gen = poly.p.true_gens()[0].clone();
    let Ok(signature) = symox::sets::lambda_checked(&gen, &poly_expr(poly)) else { return Ok(None) };
    let Ok(fun) = symox::sets::lambda_checked(var, expr) else { return Ok(None) };
    let node = symox::functions::undefined("RootSum", vec![signature, fun]);
    nodes.push((w(node.clone()), poly.clone(), auto));
    Ok(Some(node))
}

/// `RootSum.new(poly, func, auto)`.
fn new_sum(poly: &PyPoly, var: &Expr, expr: &Expr, auto: bool, nodes: &mut Vec<Node>) -> PyResult<Option<Expr>> {
    if !symox::traverse::has(expr, var) { return Ok(Some(expr.clone())); }
    if !is_rational(expr, var) || !auto { raw(poly, var, expr, auto, nodes) } else { rational_case(poly, var, expr) }
}

fn nodes_list<'py>(py: Python<'py>, value: Expr, nodes: Vec<Node>) -> PyResult<Bound<'py, PyTuple>> {
    let items = nodes.into_iter().map(|(n, p, a)| PyTuple::new(py, [n.into_pyobject(py)?.into_any(), p.into_pyobject(py)?.into_any(), a.into_pyobject(py)?.to_owned().into_any()]))
        .collect::<PyResult<Vec<_>>>()?;
    PyTuple::new(py, [w(value).into_pyobject(py)?.into_any(), PyList::new(py, items)?.into_any()])
}

/// `RootSum.new(poly, Lambda(var, expr), auto)`: `(value, raw nodes)`.
#[pyfunction]
pub(crate) fn _root_sum_new<'py>(py: Python<'py>, poly: &Bound<'py, PyAny>, var: &Bound<'py, PyAny>, expr: &Bound<'py, PyAny>, auto: bool) -> PyResult<Option<Bound<'py, PyTuple>>> {
    let Ok(p) = poly.downcast::<PyPoly>() else { return Ok(None) };
    let (var, expr) = (to_expr(var)?, to_expr(expr)?);
    let mut nodes = Vec::new();
    let Some(value) = new_sum(p.get(), &var, &expr, auto, &mut nodes)? else { return Ok(None) };
    Ok(Some(nodes_list(py, value, nodes)?))
}

/// `RootSum.__new__` after `_transform` and the argument checks, for
/// `func = Lambda(var, expr)`: `(value, raw nodes)`.
#[pyfunction]
#[pyo3(signature = (coeff, poly, var, expr, auto=true, quadratic=false))]
pub(crate) fn _root_sum_build<'py>(py: Python<'py>, coeff: &Bound<'py, PyAny>, poly: &Bound<'py, PyAny>, var: &Bound<'py, PyAny>, expr: &Bound<'py, PyAny>, auto: bool, quadratic: bool) -> PyResult<Option<Bound<'py, PyTuple>>> {
    let Ok(p) = poly.downcast::<PyPoly>() else { return Ok(None) };
    let poly = p.get();
    let (coeff, var, mut expr) = (to_expr(coeff)?, to_expr(var)?, to_expr(expr)?);
    if poly.p.true_gens().len() != 1 { return Ok(None); }
    if !(coeff.as_number().is_some_and(|n| n.is_one())) {
        expr = symox::subs::subs(&expr, &[(var.clone(), mul(&[coeff.clone(), var.clone()]))]);
    }
    let deg = degree(poly);
    if !symox::traverse::has(&expr, &var) {
        return Ok(Some(nodes_list(py, mul(&[Expr::int(deg), expr]), Vec::new())?));
    }
    let deps = std::slice::from_ref(&var);
    let add_const = if expr.is_add() {
        let (independent, dependent) = symox::decompose::as_independent(&expr, deps, None);
        expr = dependent;
        independent
    } else { Expr::zero() };
    let mul_const = if expr.is_mul() {
        let (independent, dependent) = symox::decompose::as_independent(&expr, deps, None);
        expr = dependent;
        independent
    } else { Expr::one() };
    let rational = is_rational(&expr, &var);
    let factors = poly.factor_list(py)?;
    let mut nodes = Vec::new();
    let mut terms = Vec::new();
    for item in factors.get_item(1)?.try_iter()? {
        let item = item?;
        let Ok(factor) = item.get_item(0)?.downcast::<PyPoly>().map(|p| p.get().clone()) else { return Ok(None) };
        let k: i64 = item.get_item(1)?.extract()?;
        let d = degree(&factor);
        let term = if d == 1 {
            let c = factor.p.all_coeffs().map_err(err)?;
            let r = proots::roots_linear_in(&c, root_domain(py, &factor)?);
            apply(&var, &expr, &r[0])
        } else if quadratic && d == 2 {
            let c = factor.p.all_coeffs().map_err(err)?;
            let r = proots::roots_quadratic_in(&c, root_domain(py, &factor)?);
            add(&r.iter().map(|r| apply(&var, &expr, r)).collect::<Vec<_>>())
        } else if !rational || !auto {
            let Some(node) = raw(&factor, &var, &expr, auto, &mut nodes)? else { return Ok(None) };
            node
        } else {
            let Some(value) = rational_case(&factor, &var, &expr)? else { return Ok(None) };
            value
        };
        terms.push(mul(&[Expr::int(k), term]));
    }
    let value = add(&[mul(&[mul_const, add(&terms)]), mul(&[Expr::int(deg), add_const])]);
    Ok(Some(nodes_list(py, value, nodes)?))
}

/// `RootSum.doit()`: the node itself when the roots are not all found.
#[pyfunction]
pub(crate) fn _root_sum_doit(py: Python<'_>, node: &Bound<'_, PyAny>, poly: &Bound<'_, PyAny>, fun: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    let Ok(p) = poly.downcast::<PyPoly>() else { return Ok(None) };
    let fun = to_expr(fun)?;
    let Some((vars, body)) = symox::sets::lambda_parts(&fun) else { return Ok(None) };
    if vars.len() != 1 { return Ok(None); }
    let found = super::roots(py, poly, &PyTuple::empty(py), true, None, None)?;
    let found: Vec<Expr> = found.bind(py).try_iter()?.map(|r| to_expr(&r?)).collect::<PyResult<_>>()?;
    if (found.len() as i64) < degree(p.get()) { return Ok(Some(w(to_expr(node)?))); }
    Ok(Some(w(add(&found.iter().map(|r| apply(&vars[0], &body, r)).collect::<Vec<_>>()))))
}

/// mpmath's `prec_to_dps`.
fn prec_to_dps(prec: i64) -> i64 { ((prec as f64 / 3.3219280948873626).round() as i64 - 1).max(1) }

/// `RootSum._eval_evalf(prec)` with the `as_real_imag` split of the
/// evalf fallback: `Ok(None)` keeps the node, `Err` of the decline kind
/// is never produced (the caller's body is the same algorithm).
#[pyfunction]
pub(crate) fn _root_sum_evalf(py: Python<'_>, node: &Bound<'_, PyAny>, poly: &Bound<'_, PyAny>, fun: &Bound<'_, PyAny>, prec: i64) -> PyResult<Option<PyExpr>> {
    let Ok(p) = poly.downcast::<PyPoly>() else { return Ok(None) };
    let fun = to_expr(fun)?;
    let Some((vars, body)) = symox::sets::lambda_parts(&fun) else { return Ok(None) };
    if vars.len() != 1 { return Ok(None); }
    let dps = prec_to_dps(prec);
    let found = match p.get().nroots(py, dps as u32, 50, true) {
        Ok(r) => r,
        Err(e) if is_error(py, &e, "DomainError")? || is_error(py, &e, "PolynomialError")? => return Ok(Some(w(to_expr(node)?))),
        Err(e) => return Err(e),
    };
    let terms: Vec<Expr> = found.iter().map(|r| apply(&vars[0], &body, &r.e)).collect();
    let value = add(&terms);
    for second in [false, true] {
        let candidate = if second {
            let mut evaluated = Vec::new();
            for t in &terms {
                match crate::evalf_expr(t, dps, None, None) {
                    Ok(v) => evaluated.push(v),
                    Err(_) => return Ok(Some(w(value))),
                }
            }
            add(&evaluated)
        } else { value.clone() };
        let (re, im) = match symox::complex_parts::as_real_imag(&candidate, true, None) {
            Ok(Some(pair)) => pair,
            Ok(None) => return Ok(Some(w(value))),
            // The public `as_real_imag` raises NotImplementedError, which
            // the reference does not catch: the Python body raises it.
            Err(_) => return Ok(None),
        };
        if w(im).into_pyobject(py)?.eq(0)? { return Ok(Some(w(re))); }
    }
    Ok(Some(w(value)))
}
