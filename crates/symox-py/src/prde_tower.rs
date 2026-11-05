//! `is_log_deriv_k_t_radical_in_field` and `parametric_log_deriv_heu` on
//! exp/primitive towers (see `symox::integrate::prde_tower`).
//!
//! Admission: stock Poly arguments in `DE.t` over ZZ, QQ, ZZ[...], QQ[...],
//! ZZ(...) or QQ(...) of plain symbols, at an active level `>= 1`, with
//! `DE.T` plain symbols, the active `DE.D` stock Poly objects in their
//! `DE.T` entry and `DE.cases` a list of strings whose base level is
//! `'base'` with `D = 1`. The answer `u` (or `v`) is the stock
//! `cancel(...)` of the computed rational function.
use pyo3::prelude::*;
use pyo3::types::{PyList, PyString, PyTuple};
use symox::expr::Expr;
use symox::integrate::prde_tower::{self as k, Level, Stop, Tower, UP};
use symox::poly::domains::{DomainElement, ElementValue, ExactDomain};
use symox::poly::ratfunc::{Scope, RF};
use crate::poly::PyPoly;

fn plain_symbol(e: &Expr) -> bool {
    e.is_symbol() && !e.has_host() && !e.has_wild_properties()
        && symox::assumptions::fact(e, "commutative") == Some(true)
}

fn function(space: &ExactDomain, e: &Expr) -> Option<RF> {
    match space.convert(e).ok()?.value {
        ElementValue::Function(r) => Some(r),
        _ => None,
    }
}

/// The expression of a stock Poly in `gen` with an exact ground domain.
fn poly_expr(p: &Bound<'_, PyAny>, gen: &Expr) -> PyResult<Option<Expr>> {
    if !crate::prde_const::exact_poly(p)? {
        return Ok(None);
    }
    let p = p.downcast::<PyPoly>()?.get();
    if p.p.true_gens() != [gen.clone()] {
        return Ok(None);
    }
    Ok(Some(p.p.as_expr()))
}

fn strings(v: &Bound<'_, PyAny>) -> Option<Vec<String>> {
    let l = v.downcast_exact::<PyList>().ok()?;
    l.iter().map(|s| if s.is_exact_instance_of::<PyString>() { s.extract::<String>().ok() } else { None }).collect()
}

enum Answer {
    Log(Option<(symox::poly::Int, RF)>),
    Heu(Option<(symox::poly::Int, symox::poly::Int, RF)>),
}

/// `args = (fa, fd, case, cases)` (`is_log_deriv_tower`) or `(fa, fd, wa,
/// wd, cases)` (`heu_tower`). Returns `None` (stock body) or `(status,
/// payload, trace)`: status 0 with the stock answer, 1 with the message of
/// a `NotImplementedError`, 2 with the message of a `ValueError`; `trace`
/// lists the names of the Dummies the stock body creates, in order.
#[allow(clippy::too_many_arguments)]
pub(crate) fn kernel(py: Python<'_>, kind: &str, args: &Bound<'_, PyTuple>, t: &Bound<'_, PyAny>, d: &Bound<'_, PyAny>,
    tower: &Bound<'_, PyAny>, derivations: &Bound<'_, PyAny>, level: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    let heu = kind == "heu_tower";
    if args.len() != if heu { 5 } else { 4 } {
        return Ok(None);
    }
    let (Ok(tower), Ok(derivations)) = (tower.downcast_exact::<PyList>(), derivations.downcast_exact::<PyList>()) else { return Ok(None) };
    if !level.is_exact_instance_of::<pyo3::types::PyInt>() {
        return Ok(None);
    }
    let Ok(level) = level.extract::<isize>() else { return Ok(None) };
    let nt = tower.len() as isize;
    let active = nt + level + 1;
    if active < 2 || active > nt || tower.len() != derivations.len() || tower.len() > 16 {
        return Ok(None);
    }
    let active = active as usize;
    let top = active - 1;
    if !t.is(&tower.get_item(top)?) || !d.is(&derivations.get_item(top)?) {
        return Ok(None);
    }
    let Some(cases) = strings(&args.get_item(args.len() - 1)?) else { return Ok(None) };
    if cases.len() < active || cases[0] != "base" {
        return Ok(None);
    }
    let case = if heu {
        None
    } else {
        let c = args.get_item(2)?;
        if !c.is_exact_instance_of::<PyString>() {
            return Ok(None);
        }
        let Ok(c) = c.extract::<String>() else { return Ok(None) };
        Some(c)
    };
    let mut t_vars = Vec::with_capacity(active);
    for v in tower.iter().take(active) {
        if !v.is_instance_of::<crate::PyExprNode>() {
            return Ok(None);
        }
        let e = crate::to_expr(&v)?;
        if !plain_symbol(&e) || t_vars.contains(&e) {
            return Ok(None);
        }
        t_vars.push(e);
    }
    let mut d_exprs = Vec::with_capacity(active);
    for (i, dv) in derivations.iter().take(active).enumerate() {
        let Some(e) = poly_expr(&dv, &t_vars[i])? else { return Ok(None) };
        d_exprs.push(e);
    }
    if d_exprs[0] != Expr::one() {
        return Ok(None);
    }
    let nargs = if heu { 4 } else { 2 };
    let mut inputs = Vec::with_capacity(nargs);
    for i in 0..nargs {
        let Some(e) = poly_expr(&args.get_item(i)?, &t_vars[top])? else { return Ok(None) };
        inputs.push(e);
    }
    let mut space: Vec<Expr> = t_vars.clone();
    for e in d_exprs.iter().chain(&inputs) {
        for s in symox::traverse::free_symbols(e) {
            if !plain_symbol(&s) {
                return Ok(None);
            }
            space.push(s);
        }
    }
    space.sort_by_cached_key(symox::sortkey::sort_key);
    let mut unique: Vec<Expr> = Vec::with_capacity(space.len());
    for s in space {
        if !unique.contains(&s) {
            unique.push(s);
        }
    }
    let space = unique;
    if space.len() > 32 {
        return Ok(None);
    }
    let index = |e: &Expr| space.iter().position(|s| s == e).unwrap();
    let Ok(field) = ExactDomain::functions(false, false, space.clone(), "lex") else { return Ok(None) };
    let _scope = Scope::new(space.len());
    let mut interrupted: Option<PyErr> = None;
    let outcome = symox::poly::mpoly::checked_operation(|| -> symox::poly::convert::PResult<Option<(Result<Answer, Stop>, Vec<&'static str>)>> {
        let convert = |es: &[Expr]| es.iter().map(|e| function(&field, e)).collect::<Option<Vec<_>>>();
        let (Some(dv), Some(iv)) = (convert(&d_exprs), convert(&inputs)) else { return Ok(None) };
        let levels = t_vars.iter().zip(dv).zip(&cases)
            .map(|((v, d), c)| Level { var: index(v), d, case: c.clone() }).collect();
        let var = index(&t_vars[top]);
        let Some(polys) = iv.iter().map(|f| k::up_of(f, var)).collect::<Option<Vec<UP>>>() else { return Ok(None) };
        if polys.iter().any(|p| p.len() > 65) {
            return Ok(None);
        }
        let mut check = || match py.check_signals() {
            Ok(()) => true,
            Err(e) => { interrupted = Some(e); false }
        };
        let Some(mut tw) = Tower::new(levels, space.len(), &mut check) else { return Ok(None) };
        let r = if heu {
            tw.heu_top(&polys[0], &polys[1], &polys[2], &polys[3]).map(Answer::Heu)
        } else {
            tw.log_deriv_top(&polys[0], &polys[1], case.as_deref().unwrap()).map(Answer::Log)
        };
        let trace = std::mem::take(&mut tw.trace);
        Ok(Some((r, trace)))
    });
    if let Some(e) = interrupted {
        return Err(e);
    }
    py.check_signals()?;
    let Ok(Some((r, trace))) = outcome else { return Ok(None) };
    let expr = |v: &RF| DomainElement { domain: field.clone(), value: ElementValue::Function(v.clone()) }.as_expr();
    let no_gens = PyTuple::empty(py);
    let cancelled = |v: &RF| -> PyResult<PyObject> {
        let obj = crate::wrap(expr(v)).into_pyobject(py)?.into_any();
        crate::poly::cancel(&obj, &no_gens, None)
    };
    let int = |n: symox::poly::Int| crate::wrap(symox::poly::convert::rat_to_expr(&symox::poly::Rat::from(n)));
    let (status, payload): (u8, PyObject) = match r {
        Err(Stop::Decline) => return Ok(None),
        Err(Stop::NotImplemented(m)) => (1, m.into_pyobject(py)?.into_any().unbind()),
        Err(Stop::Value(m)) => (2, m.into_pyobject(py)?.into_any().unbind()),
        Ok(Answer::Log(None)) | Ok(Answer::Heu(None)) => (0, py.None()),
        Ok(Answer::Log(Some((n, u)))) => (0, (int(n), cancelled(&u)?).into_pyobject(py)?.into_any().unbind()),
        Ok(Answer::Heu(Some((n, m, v)))) => (0, (int(n), int(m), cancelled(&v)?).into_pyobject(py)?.into_any().unbind()),
    };
    let trace = PyList::new(py, trace)?;
    py.check_signals()?;
    Ok(Some((status, payload, trace).into_pyobject(py)?.into_any().unbind()))
}
