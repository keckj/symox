//! Python bindings for the ported `sympy.simplify` helpers (`radsimp`,
//! `exprtools`, `powsimp`, `sqrtdenest`, `gammasimp`, `ratsimp`, `fu`).
//! The public signatures live in `python/symox/_simptools.py`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use symox::expr::Expr;
use symox::simptools as st;

use crate::{to_expr, wrap, PyExpr};

type Pairs = Vec<(PyExpr, PyExpr)>;

fn exprs(v: &[Bound<'_, PyAny>]) -> PyResult<Vec<Expr>> {
    v.iter().map(to_expr).collect()
}

fn pairs(v: Vec<(Expr, Expr)>) -> Pairs {
    v.into_iter().map(|(a, b)| (wrap(a), wrap(b))).collect()
}

#[pyfunction]
fn _st_fraction(e: &Bound<'_, PyAny>, exact: bool) -> PyResult<(PyExpr, PyExpr)> {
    let (n, d) = st::radsimp::fraction(&to_expr(e)?, exact);
    Ok((wrap(n), wrap(d)))
}

#[pyfunction]
fn _st_collect_const(e: &Bound<'_, PyAny>, vars: Vec<Bound<'_, PyAny>>, numbers: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::radsimp::collect_const(&to_expr(e)?, &exprs(&vars)?, numbers)))
}

#[pyfunction]
fn _st_collect_sqrt(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::radsimp::collect_sqrt(&to_expr(e)?)))
}

#[pyfunction]
fn _st_rcollect(e: &Bound<'_, PyAny>, vars: Vec<Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    Ok(wrap(st::radsimp::rcollect(&to_expr(e)?, &exprs(&vars)?)))
}

#[pyfunction]
fn _st_radsimp(e: &Bound<'_, PyAny>, symbolic: bool, max_terms: usize) -> PyResult<PyExpr> {
    Ok(wrap(st::radsimp::radsimp(&to_expr(e)?, symbolic, max_terms)))
}

#[pyfunction]
fn _st_rad_rationalize(num: &Bound<'_, PyAny>, den: &Bound<'_, PyAny>) -> PyResult<(PyExpr, PyExpr)> {
    let (n, d) = st::radsimp::rad_rationalize(&to_expr(num)?, &to_expr(den)?);
    Ok((wrap(n), wrap(d)))
}

#[pyfunction]
fn _st_split_surds(e: &Bound<'_, PyAny>) -> PyResult<(PyExpr, PyExpr, PyExpr)> {
    let (g, a, b) = st::radsimp::split_surds(&to_expr(e)?);
    Ok((wrap(g), wrap(a), wrap(b)))
}

#[pyfunction]
fn _st_signsimp(e: &Bound<'_, PyAny>, evaluate: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::misc::signsimp(&to_expr(e)?, evaluate)))
}

#[pyfunction]
fn _st_posify(e: &Bound<'_, PyAny>) -> PyResult<(PyExpr, Pairs)> {
    let (out, reps) = st::misc::posify(&to_expr(e)?);
    Ok((wrap(out), pairs(reps)))
}

#[pyfunction]
fn _st_posify_seq(v: Vec<Bound<'_, PyAny>>) -> PyResult<(Vec<PyExpr>, Pairs)> {
    let (out, reps) = st::misc::posify_seq(&exprs(&v)?);
    Ok((out.into_iter().map(wrap).collect(), pairs(reps)))
}

#[pyfunction]
fn _st_separatevars(e: &Bound<'_, PyAny>, force: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::misc::separatevars(&to_expr(e)?, force)))
}

#[pyfunction]
#[pyo3(signature = (e, symbols, force))]
fn _st_separatevars_dict(e: &Bound<'_, PyAny>, symbols: Option<Vec<Bound<'_, PyAny>>>, force: bool) -> PyResult<Option<Vec<(Option<PyExpr>, PyExpr)>>> {
    let e = to_expr(e)?;
    let out = match symbols {
        None => Some(st::misc::separatevars_dict_none(&e, force)),
        Some(s) => {
            let s = exprs(&s)?;
            if !s.iter().all(Expr::is_atom) {
                return Err(PyValueError::new_err("symbols must be Atoms."));
            }
            st::misc::separatevars_dict(&e, Some(&s), force)
        }
    };
    Ok(out.map(|v| v.into_iter().map(|(k, v)| (k.map(wrap), wrap(v))).collect()))
}

#[pyfunction]
fn _st_hypersimp(f: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    Ok(st::misc::hypersimp(&to_expr(f)?, &to_expr(k)?).map(wrap))
}

#[pyfunction]
fn _st_hypersimilar(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>, k: &Bound<'_, PyAny>) -> PyResult<Option<bool>> {
    Ok(st::misc::hypersimilar(&to_expr(f)?, &to_expr(g)?, &to_expr(k)?))
}

#[pyfunction]
fn _st_powsimp(e: &Bound<'_, PyAny>, deep: bool, combine: &str, force: bool) -> PyResult<PyExpr> {
    let combine = match combine {
        "all" => st::powsimp::Combine::All,
        "exp" => st::powsimp::Combine::Exp,
        "base" => st::powsimp::Combine::Base,
        _ => return Err(PyValueError::new_err("combine must be one of ('all', 'exp', 'base').")),
    };
    Ok(wrap(st::powsimp::powsimp(&to_expr(e)?, &st::powsimp::PowsimpOptions { deep, combine, force })))
}

#[pyfunction]
fn _st_powdenest(e: &Bound<'_, PyAny>, force: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::powsimp::powdenest(&to_expr(e)?, force, false)))
}

#[pyfunction]
fn _st_factor_terms(e: &Bound<'_, PyAny>, radical: bool, clear: bool, fraction: bool, sign: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::exprtools::factor_terms(&to_expr(e)?, radical, clear, fraction, sign)))
}

#[pyfunction]
fn _st_gcd_terms(e: &Bound<'_, PyAny>, isprimitive: bool, clear: bool, fraction: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::exprtools::gcd_terms(&to_expr(e)?, isprimitive, clear, fraction)))
}

#[pyfunction]
fn _st_gcd_terms_list(v: Vec<Bound<'_, PyAny>>, isprimitive: bool, clear: bool, fraction: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::exprtools::gcd_terms_list(&exprs(&v)?, isprimitive, clear, fraction)))
}

#[pyfunction]
fn _st_factor_nc(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::exprtools::factor_nc(&to_expr(e)?)))
}

#[pyfunction]
fn _has_noncommutative_atoms(e: &Bound<'_, PyAny>) -> bool {
    to_expr(e).map(|e| symox::decompose::has_noncommutative_atoms(&e)).unwrap_or(false)
}

#[pyfunction]
fn _st_primitive(e: &Bound<'_, PyAny>) -> PyResult<(PyExpr, PyExpr)> {
    let (c, p) = st::primitive(&to_expr(e)?);
    Ok((wrap(Expr::num(c)), wrap(p)))
}

#[pyfunction]
fn _st_extract_additively(e: &Bound<'_, PyAny>, c: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    Ok(st::extract_additively(&to_expr(e)?, &to_expr(c)?).map(wrap))
}

#[pyfunction]
fn _st_extract_multiplicatively(e: &Bound<'_, PyAny>, c: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    Ok(st::extract_multiplicatively(&to_expr(e)?, &to_expr(c)?).map(wrap))
}

#[pyfunction]
fn _st_as_coefficient(e: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>) -> PyResult<Option<PyExpr>> {
    Ok(st::as_coefficient(&to_expr(e)?, &to_expr(x)?).map(wrap))
}

#[pyfunction]
fn _st_coeff(e: &Bound<'_, PyAny>, x: &Bound<'_, PyAny>, n: i64) -> PyResult<PyExpr> {
    Ok(wrap(st::coeff(&to_expr(e)?, &to_expr(x)?, n)))
}

#[pyfunction]
fn _st_normal(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::normal(&to_expr(e)?)))
}

#[pyfunction]
fn _st_equals(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<Option<bool>> {
    Ok(st::methods::equals(&to_expr(a)?, &to_expr(b)?))
}

#[pyfunction]
fn _st_is_constant(e: &Bound<'_, PyAny>, wrt: Vec<Bound<'_, PyAny>>, simplify: bool) -> PyResult<Option<bool>> {
    Ok(st::methods::is_constant_with(&to_expr(e)?, &exprs(&wrt)?, simplify))
}

#[pyfunction]
fn _st_count_ops(e: &Bound<'_, PyAny>) -> PyResult<usize> {
    Ok(st::methods::count_ops(&to_expr(e)?))
}

#[pyfunction]
fn _st_sqrtdenest(e: &Bound<'_, PyAny>, max_iter: usize) -> PyResult<PyExpr> {
    Ok(wrap(st::sqrtdenest::sqrtdenest(&to_expr(e)?, max_iter)))
}

#[pyfunction]
fn _st_gammasimp(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::gammasimp::gammasimp(&to_expr(e)?)))
}

#[pyfunction]
fn _st_combsimp(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::gammasimp::combsimp(&to_expr(e)?)))
}

#[pyfunction]
fn _st_compare(a: &Bound<'_, PyAny>, b: &Bound<'_, PyAny>) -> PyResult<i32> {
    Ok(st::methods::compare(&to_expr(a)?, &to_expr(b)?))
}

#[pyfunction]
fn _st_ratsimp(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::ratsimp::ratsimp(&to_expr(e)?)))
}

#[pyfunction]
fn _st_ratsimpmodprime(
    e: &Bound<'_, PyAny>,
    g: Vec<Bound<'_, PyAny>>,
    gens: Vec<Bound<'_, PyAny>>,
    quick: bool,
    polynomial: bool,
    order: &str,
    modulus: u64,
    force_field: bool,
) -> PyResult<PyExpr> {
    let order = symox::poly::groebner::MonomialOrder::parse(order)
        .ok_or_else(|| PyValueError::new_err(format!("unknown monomial order {:?}", order)))?;
    st::ratsimp::ratsimpmodprime(&to_expr(e)?, &exprs(&g)?, &exprs(&gens)?, quick, polynomial, order,
        (modulus > 1).then_some(modulus), force_field)
        .map(wrap)
        .map_err(|err| {
            let msg = err.to_string();
            if msg.contains("not implemented") {
                pyo3::exceptions::PyNotImplementedError::new_err(msg)
            } else {
                PyValueError::new_err(msg)
            }
        })
}

// --- fu --------------------------------------------------------------------

use std::cell::RefCell;
use symox::expr::Func;
use symox::simptools::fu;

fn trig_func(name: &str) -> PyResult<Func> {
    Ok(match name {
        "sin" => Func::Sin,
        "cos" => Func::Cos,
        "tan" => Func::Tan,
        "cot" => Func::Cot,
        "sec" => Func::Sec,
        "csc" => Func::Csc,
        _ => return Err(PyValueError::new_err(format!("unsupported function {name}"))),
    })
}

#[pyfunction]
#[pyo3(signature = (name, e, max=None, pow=false, half=false, first=true, base=None))]
fn _st_tr(
    name: &str,
    e: &Bound<'_, PyAny>,
    max: Option<Bound<'_, PyAny>>,
    pow: bool,
    half: bool,
    first: bool,
    base: Option<Bound<'_, PyAny>>,
) -> PyResult<PyExpr> {
    let e = to_expr(e)?;
    let max = match max {
        Some(m) => to_expr(&m)?,
        None => Expr::int(4),
    };
    let base = base.map(|b| to_expr(&b)).transpose()?;
    let r = match name {
        "TR0" => fu::tr0(&e),
        "TR1" => fu::tr1(&e),
        "TR2" => fu::tr2(&e),
        "TR2i" => fu::tr2i(&e, half),
        "TR3" => fu::tr3(&e),
        "TR4" => fu::tr4(&e),
        "TR5" => fu::tr5(&e, &max, pow),
        "TR6" => fu::tr6(&e, &max, pow),
        "TR7" => fu::tr7(&e),
        "TR8" => fu::tr8(&e, first),
        "TR9" => fu::tr9(&e),
        "TR10" => fu::tr10(&e, first),
        "TR10i" => fu::tr10i(&e),
        "TR11" => fu::tr11(&e, base.as_ref()),
        "_TR11" => fu::tr11_ratio(&e),
        "TR12" => fu::tr12(&e, first),
        "TR12i" => fu::tr12i(&e),
        "TR13" => fu::tr13(&e),
        "TR14" => fu::tr14(&e, first),
        "TR15" => fu::tr15(&e, &max, pow),
        "TR16" => fu::tr16(&e, &max, pow),
        "TR111" => fu::tr111(&e),
        "TR22" => fu::tr22(&e, &max, pow),
        "TRmorrie" => fu::trmorrie(&e),
        "TRpower" => fu::trpower(&e),
        "sincos_to_sum" => fu::sincos_to_sum(&e),
        _ => return Err(PyValueError::new_err(format!("unknown rule {name}"))),
    };
    Ok(wrap(r))
}

#[pyfunction]
fn _st_fu_l(e: &Bound<'_, PyAny>) -> PyResult<usize> {
    Ok(fu::l_count(&to_expr(e)?))
}

/// Run `f` with a comparator built from an optional Python `measure`.
fn with_less<R>(measure: Option<&Bound<'_, PyAny>>, f: impl FnOnce(fu::Less) -> R) -> PyResult<R> {
    let Some(m) = measure else { return Ok(f(&fu::default_less)) };
    let err: RefCell<Option<PyErr>> = RefCell::new(None);
    let less = |a: &Expr, b: &Expr| -> bool {
        if err.borrow().is_some() {
            return false;
        }
        let r = (|| -> PyResult<bool> {
            let ma = m.call1((wrap(a.clone()),))?;
            let mb = m.call1((wrap(b.clone()),))?;
            ma.lt(mb)
        })();
        r.unwrap_or_else(|e| {
            *err.borrow_mut() = Some(e);
            false
        })
    };
    let out = f(&less);
    match err.into_inner() {
        Some(e) => Err(e),
        None => Ok(out),
    }
}

#[pyfunction]
#[pyo3(signature = (e, measure=None))]
fn _st_fu(e: &Bound<'_, PyAny>, measure: Option<Bound<'_, PyAny>>) -> PyResult<PyExpr> {
    let e = to_expr(e)?;
    with_less(measure.as_ref(), |less| wrap(fu::fu(&e, less)))
}

#[pyfunction]
fn _st_tr56(
    e: &Bound<'_, PyAny>,
    f: &str,
    g: &str,
    h: &Bound<'_, PyAny>,
    max: &Bound<'_, PyAny>,
    pow: bool,
) -> PyResult<PyExpr> {
    let (f, g) = (trig_func(f)?, trig_func(g)?);
    let err: RefCell<Option<PyErr>> = RefCell::new(None);
    let hh = |x: &Expr| -> Expr {
        let r = h.call1((wrap(x.clone()),)).and_then(|r| to_expr(&r));
        r.unwrap_or_else(|e| {
            *err.borrow_mut() = Some(e);
            x.clone()
        })
    };
    let out = fu::tr56(&to_expr(e)?, f, g, &hh, &to_expr(max)?, pow);
    match err.into_inner() {
        Some(e) => Err(e),
        None => Ok(wrap(out)),
    }
}

#[pyfunction]
#[pyo3(signature = (a, b, two=false))]
fn _st_trig_split(
    a: &Bound<'_, PyAny>,
    b: &Bound<'_, PyAny>,
    two: bool,
) -> PyResult<Option<(PyExpr, i64, i64, PyExpr, PyExpr, bool)>> {
    Ok(fu::trig_split(&to_expr(a)?, &to_expr(b)?, two).map(|(g, n1, n2, a, b, t)| (wrap(g), n1, n2, wrap(a), wrap(b), t)))
}

#[pyfunction]
fn _st_as_f_sign_1(e: &Bound<'_, PyAny>) -> PyResult<Option<(PyExpr, PyExpr, PyExpr)>> {
    Ok(fu::as_f_sign_1(&to_expr(e)?).map(|(g, f, s)| (wrap(g), wrap(f), wrap(s))))
}

#[pyfunction]
fn _st_osborne(e: &Bound<'_, PyAny>, d: &Bound<'_, PyAny>, inverse: bool) -> PyResult<PyExpr> {
    let (e, d) = (to_expr(e)?, to_expr(d)?);
    Ok(wrap(if inverse { fu::osbornei(&e, &d) } else { fu::osborne(&e, &d) }))
}

#[pyfunction]
fn _st_hyper_as_trig(e: &Bound<'_, PyAny>) -> PyResult<(PyExpr, PyExpr, Pairs)> {
    let (t, d, back) = fu::hyper_as_trig(&to_expr(e)?);
    Ok((wrap(t), wrap(d), pairs(back)))
}

#[pyfunction]
fn _st_hyper_as_trig_back(x: &Bound<'_, PyAny>, d: &Bound<'_, PyAny>, back: Vec<(Bound<'_, PyAny>, Bound<'_, PyAny>)>) -> PyResult<PyExpr> {
    let back = back.iter().map(|(a, b)| Ok((to_expr(a)?, to_expr(b)?))).collect::<PyResult<Vec<_>>>()?;
    Ok(wrap(fu::hyper_as_trig_back(&to_expr(x)?, &to_expr(d)?, &back)))
}

#[pyfunction]
#[pyo3(signature = (e, hyper=true))]
fn _st_futrig(e: &Bound<'_, PyAny>, hyper: bool) -> PyResult<PyExpr> {
    Ok(wrap(st::trigsimp::futrig(&to_expr(e)?, hyper)))
}

#[pyfunction]
fn _st_exptrigsimp(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::trigsimp::exptrigsimp(&to_expr(e)?)))
}

#[pyfunction]
fn _st_trigsimp_inverse(e: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(wrap(st::trigsimp::trigsimp_inverse(&to_expr(e)?)))
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_st_futrig, m)?)?;
    m.add_function(wrap_pyfunction!(_st_exptrigsimp, m)?)?;
    m.add_function(wrap_pyfunction!(_st_trigsimp_inverse, m)?)?;
    m.add_function(wrap_pyfunction!(_st_tr, m)?)?;
    m.add_function(wrap_pyfunction!(_st_fu_l, m)?)?;
    m.add_function(wrap_pyfunction!(_st_fu, m)?)?;
    m.add_function(wrap_pyfunction!(_st_tr56, m)?)?;
    m.add_function(wrap_pyfunction!(_st_trig_split, m)?)?;
    m.add_function(wrap_pyfunction!(_st_as_f_sign_1, m)?)?;
    m.add_function(wrap_pyfunction!(_st_osborne, m)?)?;
    m.add_function(wrap_pyfunction!(_st_hyper_as_trig, m)?)?;
    m.add_function(wrap_pyfunction!(_st_hyper_as_trig_back, m)?)?;
    m.add_function(wrap_pyfunction!(_st_ratsimp, m)?)?;
    m.add_function(wrap_pyfunction!(_st_ratsimpmodprime, m)?)?;
    m.add_function(wrap_pyfunction!(_st_fraction, m)?)?;
    m.add_function(wrap_pyfunction!(_st_collect_const, m)?)?;
    m.add_function(wrap_pyfunction!(_st_collect_sqrt, m)?)?;
    m.add_function(wrap_pyfunction!(_st_rcollect, m)?)?;
    m.add_function(wrap_pyfunction!(_st_radsimp, m)?)?;
    m.add_function(wrap_pyfunction!(_st_rad_rationalize, m)?)?;
    m.add_function(wrap_pyfunction!(_st_split_surds, m)?)?;
    m.add_function(wrap_pyfunction!(_st_signsimp, m)?)?;
    m.add_function(wrap_pyfunction!(_st_posify, m)?)?;
    m.add_function(wrap_pyfunction!(_st_posify_seq, m)?)?;
    m.add_function(wrap_pyfunction!(_st_separatevars, m)?)?;
    m.add_function(wrap_pyfunction!(_st_separatevars_dict, m)?)?;
    m.add_function(wrap_pyfunction!(_st_hypersimp, m)?)?;
    m.add_function(wrap_pyfunction!(_st_hypersimilar, m)?)?;
    m.add_function(wrap_pyfunction!(_st_powsimp, m)?)?;
    m.add_function(wrap_pyfunction!(_st_powdenest, m)?)?;
    m.add_function(wrap_pyfunction!(_st_factor_terms, m)?)?;
    m.add_function(wrap_pyfunction!(_st_gcd_terms, m)?)?;
    m.add_function(wrap_pyfunction!(_st_gcd_terms_list, m)?)?;
    m.add_function(wrap_pyfunction!(_st_factor_nc, m)?)?;
    m.add_function(wrap_pyfunction!(_has_noncommutative_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(_st_primitive, m)?)?;
    m.add_function(wrap_pyfunction!(_st_extract_additively, m)?)?;
    m.add_function(wrap_pyfunction!(_st_extract_multiplicatively, m)?)?;
    m.add_function(wrap_pyfunction!(_st_as_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(_st_coeff, m)?)?;
    m.add_function(wrap_pyfunction!(_st_normal, m)?)?;
    m.add_function(wrap_pyfunction!(_st_equals, m)?)?;
    m.add_function(wrap_pyfunction!(_st_is_constant, m)?)?;
    m.add_function(wrap_pyfunction!(_st_count_ops, m)?)?;
    m.add_function(wrap_pyfunction!(_st_sqrtdenest, m)?)?;
    m.add_function(wrap_pyfunction!(_st_compare, m)?)?;
    m.add_function(wrap_pyfunction!(_st_gammasimp, m)?)?;
    m.add_function(wrap_pyfunction!(_st_combsimp, m)?)?;
    Ok(())
}
