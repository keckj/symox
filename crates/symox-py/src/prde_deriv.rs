//! `is_deriv_k(fa, fd, DE)` of SymPy `prde.py` (structure theorem for the
//! logarithmic derivative `Df/f`) over rational function fields.
//!
//! The stock body is reproduced step by step: `Df/f` is computed with the
//! native `derivation` loop and the Rust Poly `*`, `-` and `cancel`; the
//! `E_part + L_part` row and the `Df/f` entry get the PolyMatrix domains of
//! `PolyMatrix([...], dum)`; `constant_system` is the native fraction-field
//! kernel; the constant and rational checks of `u` are exact rational
//! function computations (`D(u_i) == 0` over the active levels, `u_i` a
//! rational number); the answer `(ans, u, const)` is built with the native
//! `Add`/`Mul`/`Pow` constructors, `as_numer_denom`, `sqf_list` and
//! `cancel`. Any step the kernel does not reproduce returns `Ok(None)`.
use pyo3::exceptions::PyNotImplementedError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyInt, PyList, PyString, PyStringMethods, PyTuple};
use symox::expr::Expr;
use symox::integrate::prde_const as k;
use symox::ops;
use crate::poly::PyPoly;
use crate::prde_const::{derivation_poly, exact_poly, field_system, infer_domain, items};

fn expr_of(v: &Bound<'_, PyAny>) -> PyResult<Option<Expr>> {
    if !v.is_instance_of::<crate::PyExprNode>() {
        return Ok(None);
    }
    let e = crate::to_expr(v)?;
    Ok((!e.has_host() && !e.has_wild_properties()).then_some(e))
}

enum Outcome {
    /// `not all(derivation(i, DE, basic=True).is_zero for i in u)`.
    NotConstant,
    /// `not all(i.is_Rational for i in u)`.
    NotRational,
    Rationals(Vec<Expr>),
}

/// `args = (fa, fd, dum, exts, extargs)`, `t = DE.t`, `tower = DE.T`,
/// `derivations = DE.D`, `level = DE.level`, with `len(DE.exts) ==
/// len(DE.D)` checked by the admission. Returns `(result,)` where `result`
/// is the stock return value (`None` or `(ans, u, const)`), raises the
/// stock `NotImplementedError` for non-rational coefficients, and returns
/// `None` to request the Python body.
pub(crate) fn is_deriv_k(py: Python<'_>, args: &Bound<'_, PyTuple>, t: &Bound<'_, PyAny>,
    tower: &Bound<'_, PyAny>, derivations: &Bound<'_, PyAny>, level: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    if args.len() != 5 {
        return Ok(None);
    }
    let (fa, fd, dum) = (args.get_item(0)?, args.get_item(1)?, args.get_item(2)?);
    let (exts, extargs) = (args.get_item(3)?, args.get_item(4)?);
    let (Ok(exts), Ok(extargs)) = (exts.downcast_exact::<PyList>(), extargs.downcast_exact::<PyList>()) else { return Ok(None) };
    let (Ok(tl), Ok(dl)) = (tower.downcast_exact::<PyList>(), derivations.downcast_exact::<PyList>()) else { return Ok(None) };
    if exts.len() != dl.len() || !exact_poly(&fa)? || !exact_poly(&fd)? || !level.is_exact_instance_of::<PyInt>() {
        return Ok(None);
    }
    let Ok(lv) = level.extract::<isize>() else { return Ok(None) };
    // dfa, dfd = fd*derivation(fa, DE) - fa*derivation(fd, DE), fd*fa
    let Some(da) = derivation_poly(py, &(&fa, false).into_pyobject(py)?, t, tower, derivations, level)? else { return Ok(None) };
    let Some(dd) = derivation_poly(py, &(&fd, false).into_pyobject(py)?, t, tower, derivations, level)? else { return Ok(None) };
    let num = fd.mul(da.bind(py))?.sub(fa.mul(dd.bind(py))?)?;
    let den = fd.mul(&fa)?;
    if !num.is_exact_instance_of::<PyPoly>() || !den.is_exact_instance_of::<PyPoly>() {
        return Ok(None);
    }
    let kw = PyDict::new(py);
    kw.set_item("include", true)?;
    let pair = num.call_method("cancel", (&den,), Some(&kw))?;
    let Ok(pair) = pair.downcast_exact::<PyTuple>() else { return Ok(None) };
    if pair.len() != 2 {
        return Ok(None);
    }
    // E_part and L_part (DE.indices over every extension).
    let (mut exp_idx, mut log_idx) = (Vec::new(), Vec::new());
    for (i, e) in exts.iter().enumerate() {
        if let Ok(s) = e.downcast_exact::<PyString>() {
            match s.to_cow()?.as_ref() {
                "exp" => exp_idx.push(i),
                "log" => log_idx.push(i),
                _ => {}
            }
        } else if !e.is_none() {
            return Ok(None);
        }
    }
    let mut row = Vec::with_capacity(exp_idx.len() + log_idx.len());
    for &i in &exp_idx {
        if i >= tl.len() {
            return Ok(None);
        }
        let (d, v) = (dl.get_item(i)?, tl.get_item(i)?);
        if !d.is_exact_instance_of::<PyPoly>() || expr_of(&v)?.is_none() {
            return Ok(None);
        }
        let pv = PyPoly::_from_normalized(&v, &PyTuple::new(py, [&v])?, None, None, None)?;
        row.push(d.call_method1("quo", (pv,))?.call_method0("as_expr")?);
    }
    for &i in &log_idx {
        let d = dl.get_item(i)?;
        if !d.is_exact_instance_of::<PyPoly>() {
            return Ok(None);
        }
        row.push(d.call_method0("as_expr")?);
    }
    let rhs = pair.get_item(0)?.call_method0("as_expr")?.div(pair.get_item(1)?.call_method0("as_expr")?)?;
    let n = row.len();
    let (Some(ai), Some(ui)) = (items(PyList::new(py, &row)?.as_any(), n)?, items(PyList::new(py, [&rhs])?.as_any(), 1)?) else { return Ok(None) };
    let Some(dlu) = infer_domain(py, &dum, &ui)? else { return Ok(None) };
    if n == 0 {
        // Empty E_part + L_part: the stock constant_system returns (A, u)
        // unchanged and `not A` gives None.
        return Ok(Some((py.None(),).into_pyobject(py)?.into_any().unbind()));
    }
    let Some(dla) = infer_domain(py, &dum, &ai)? else { return Ok(None) };
    let outcome = field_system(py, &dum, (1, n, ai, dla), (1, 1, ui, dlu), tl, dl, lv, |expr, _a, u, tw| {
        if u.values.iter().any(|v| !k::derivation(v, &tw.active).n.is_zero()) {
            return Outcome::NotConstant;
        }
        if !u.values.iter().all(|v| v.n.is_constant() && v.d.is_constant()) {
            return Outcome::NotRational;
        }
        Outcome::Rationals(u.values.iter().map(expr).collect())
    })?;
    let us = match outcome {
        None => return Ok(None),
        Some(Outcome::NotConstant) => return Ok(Some((py.None(),).into_pyobject(py)?.into_any().unbind())),
        Some(Outcome::NotRational) => return Err(PyNotImplementedError::new_err("Cannot work with non-rational coefficients in this case.")),
        Some(Outcome::Rationals(us)) => us,
    };
    if us.iter().any(|u| !u.as_number().is_some_and(|q| q.is_rational())) {
        return Ok(None);
    }
    // terms = [DE.extargs[i] for exp] + [DE.T[i] for log]; ans = list(zip(terms, u))
    let mut terms = Vec::new();
    let mut argterms = Vec::new();
    for &i in &exp_idx {
        terms.push(extargs.get_item(i)?);
        argterms.push(tl.get_item(i)?);
    }
    for &i in &log_idx {
        terms.push(tl.get_item(i)?);
        argterms.push(extargs.get_item(i)?);
    }
    let ans = PyList::empty(py);
    let mut summands = Vec::new();
    for (term, j) in terms.iter().zip(&us) {
        let Some(te) = expr_of(term)? else { return Ok(None) };
        summands.push(ops::mul(&[te, j.clone()]));
        ans.append((term, crate::wrap(j.clone())))?;
    }
    let result = ops::add(&summands);
    let no_gens = PyTuple::empty(py);
    let mut l = Vec::new();
    let mut ld = Vec::new();
    for (arg, j) in argterms.iter().zip(&us) {
        if expr_of(arg)?.is_none() {
            return Ok(None);
        }
        let nd = arg.call_method0("as_numer_denom")?;
        let Ok(nd) = nd.downcast_exact::<PyTuple>() else { return Ok(None) };
        if nd.len() != 2 {
            return Ok(None);
        }
        for (part, out) in [(nd.get_item(0)?, &mut l), (nd.get_item(1)?, &mut ld)] {
            let sqf = crate::poly::sqf_list(&part, &no_gens, None, None)?.into_bound(py);
            let Ok((coeff, factors)) = sqf.extract::<(Bound<'_, PyAny>, Bound<'_, PyList>)>() else { return Ok(None) };
            let Some(c) = expr_of(&coeff)? else { return Ok(None) };
            let mut factor = vec![ops::pow(&c, j)];
            for f in factors.iter() {
                let Ok((b, e)) = f.extract::<(Bound<'_, PyAny>, Bound<'_, PyAny>)>() else { return Ok(None) };
                let Some(b) = expr_of(&b)? else { return Ok(None) };
                if !e.is_exact_instance_of::<PyInt>() {
                    return Ok(None);
                }
                let e = crate::to_expr(&e)?;
                factor.push(ops::pow(&b, &ops::mul(&[e, j.clone()])));
            }
            out.push(ops::mul(&factor));
        }
    }
    // const = cancel(fa.as_expr()/fd.as_expr()/Mul(*l)*Mul(*ld))
    let (Some(fae), Some(fde)) = (expr_of(&fa.call_method0("as_expr")?)?, expr_of(&fd.call_method0("as_expr")?)?) else { return Ok(None) };
    let q = ops::mul(&[ops::div(&ops::div(&fae, &fde), &ops::mul(&l)), ops::mul(&ld)]);
    let q = crate::wrap(q).into_pyobject(py)?.into_any();
    let constant = crate::poly::cancel(&q, &no_gens, None)?;
    let out = (ans, crate::wrap(result), constant);
    Ok(Some((out,).into_pyobject(py)?.into_any().unbind()))
}
