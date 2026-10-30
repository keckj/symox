//! `constant_system(A, u, DE)` over rational function fields (see
//! `symox::integrate::prde_const`).
//!
//! Admission: stock PolyMatrix items (expressions without floats) over QQ or
//! a fraction field `ZZ(...)`/`QQ(...)` of plain commutative symbols, a plain
//! symbol generator, `DE.T` plain symbols and the active `DE.D` stock Poly
//! objects whose expressions are rational functions of plain symbols. The
//! PolyMatrix domain inference calls the Poly constructor itself.
use pyo3::prelude::*;
use pyo3::types::{PyInt, PyList, PyTuple};
use symox::expr::Expr;
use symox::integrate::prde_const::{self as k, Host, Matrix, Tower};
use symox::poly::domains::{DomainElement, ElementValue, ExactDomain};
use symox::poly::ratfunc::{Scope, RF};
use crate::poly::domains::PyDomain;
use crate::poly::PyPoly;

fn plain_symbol(e: &Expr) -> bool {
    e.is_symbol() && !e.has_host() && !e.has_wild_properties()
        && symox::assumptions::fact(e, "commutative") == Some(true)
}

/// QQ or a fraction field over plain symbols.
fn admitted(d: &ExactDomain) -> bool {
    match d {
        ExactDomain::QQ => true,
        ExactDomain::Functions { ring: false, parameters, .. } => parameters.iter().all(plain_symbol),
        _ => false,
    }
}

struct PyHost<'py> {
    py: Python<'py>,
    gens: Bound<'py, PyTuple>,
    space: ExactDomain,
}

impl PyHost<'_> {
    fn domain_of(&self, e: &Expr, ground: Option<&Bound<'_, PyAny>>) -> Option<PyPoly> {
        let obj = crate::wrap(e.clone()).into_pyobject(self.py).ok()?.into_any();
        PyPoly::_from_normalized(&obj, &self.gens, ground, None, None).ok()
    }
}

impl Host for PyHost<'_> {
    type D = ExactDomain;
    fn infer(&mut self, items: &[Expr]) -> Option<ExactDomain> {
        // parallel_poly_from_expr(items, gens, field=True)
        let mut domains = Vec::with_capacity(items.len());
        for e in items {
            let p = self.domain_of(e, None)?.to_field().ok()?;
            domains.push(p.domain_descriptor().ok()?.inner);
        }
        let first = domains.first()?.clone();
        let label = first.to_string();
        if domains.iter().all(|d| d.to_string() == label) {
            return admitted(&first).then_some(first);
        }
        // _field_joint: construct_domain(EX coefficients, field=True).
        let ex = self.py.import("symox._core").ok()?.getattr("EX").ok()?;
        let mut coefficients = Vec::new();
        for e in items {
            coefficients.extend(self.domain_of(e, Some(&ex))?.p.coeffs());
        }
        if coefficients.iter().all(|c| c.as_number().is_some_and(|n| n.is_rational())) {
            return Some(ExactDomain::QQ);
        }
        let mut symbols: Vec<Expr> = Vec::new();
        for c in &coefficients {
            if c.as_number().is_some_and(|n| n.is_float()) {
                return None;
            }
            for s in symox::traverse::free_symbols(c) {
                if !plain_symbol(&s) {
                    return None;
                }
                if !symbols.contains(&s) {
                    symbols.push(s);
                }
            }
        }
        if symbols.is_empty() {
            return None;
        }
        symbols.sort_by_cached_key(symox::sortkey::sort_key);
        let joint = ExactDomain::functions(true, false, symbols, "lex").ok()?;
        for c in &coefficients {
            joint.convert(c).ok()?;
        }
        Some(joint)
    }
    fn unify(&mut self, a: &ExactDomain, b: &ExactDomain) -> ExactDomain {
        if a == b {
            return a.clone();
        }
        crate::poly::domains::unify_native(a, b).unwrap_or_else(|_| a.clone())
    }
    fn element(&mut self, d: &ExactDomain, v: &RF) -> Option<Expr> {
        if !admitted(d) {
            return None;
        }
        let e = self.expr(v);
        d.convert(&e).ok().map(|x| x.as_expr())
    }
    fn is_field(&mut self, d: &ExactDomain) -> bool {
        admitted(d) && d.is_field()
    }
    fn expr(&mut self, v: &RF) -> Expr {
        DomainElement { domain: self.space.clone(), value: ElementValue::Function(v.clone()) }.as_expr()
    }
}

fn items(v: &Bound<'_, PyAny>, len: usize) -> PyResult<Option<Vec<Expr>>> {
    let Ok(l) = v.downcast_exact::<PyList>() else { return Ok(None) };
    if l.len() != len {
        return Ok(None);
    }
    let mut out = Vec::with_capacity(len);
    for item in l.iter() {
        if !item.is_instance_of::<crate::PyExprNode>() {
            return Ok(None);
        }
        let e = crate::to_expr(&item)?;
        let mut bad = e.has_wild_properties() || e.has_host();
        symox::traverse::preorder(&e, &mut |x| {
            if x.as_number().is_some_and(|n| n.is_float()) {
                bad = true;
            }
        });
        if bad {
            return Ok(None);
        }
        out.push(e);
    }
    Ok(Some(out))
}

fn size(v: &Bound<'_, PyAny>) -> Option<usize> {
    if !v.is_exact_instance_of::<PyInt>() {
        return None;
    }
    v.extract::<usize>().ok().filter(|n| *n <= 256)
}

fn domain(v: &Bound<'_, PyAny>) -> PyResult<Option<ExactDomain>> {
    let Ok(d) = v.downcast_exact::<PyDomain>() else { return Ok(None) };
    let d = d.get().inner.clone();
    Ok(admitted(&d).then_some(d))
}

/// Python slice `seq[:k]` length.
fn prefix(len: usize, k: isize) -> usize {
    if k >= 0 { (k as usize).min(len) } else { len.saturating_sub(k.unsigned_abs()) }
}

fn function(space: &ExactDomain, e: &Expr) -> Option<RF> {
    match space.convert(e).ok()?.value {
        ElementValue::Function(r) => Some(r),
        _ => None,
    }
}

/// `args = (A_rows, A_cols, A_items, A_domain, u_rows, u_cols, u_items,
/// u_domain)`; returns `((rows, cols, items, domain), (rows, 1, items,
/// domain))` or `None`.
pub(crate) fn constant_system_field(py: Python<'_>, args: &Bound<'_, PyTuple>, gen: &Bound<'_, PyAny>,
    tower: &Bound<'_, PyAny>, derivations: &Bound<'_, PyAny>, level: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    if args.len() != 8 || !gen.is_instance_of::<crate::PyExprNode>() {
        return Ok(None);
    }
    let g = crate::to_expr(gen)?;
    if !plain_symbol(&g) {
        return Ok(None);
    }
    let (Ok(tower), Ok(derivations)) = (tower.downcast_exact::<PyList>(), derivations.downcast_exact::<PyList>()) else { return Ok(None) };
    if !level.is_exact_instance_of::<PyInt>() {
        return Ok(None);
    }
    let Ok(level) = level.extract::<isize>() else { return Ok(None) };
    if level.unsigned_abs() > 1 << 16 {
        return Ok(None);
    }
    let arg = |i: usize| args.get_item(i);
    let (Some(ar), Some(ac), Some(ur), Some(uc)) = (size(&arg(0)?), size(&arg(1)?), size(&arg(4)?), size(&arg(5)?)) else { return Ok(None) };
    let (Some(da), Some(du)) = (domain(&arg(3)?)?, domain(&arg(7)?)?) else { return Ok(None) };
    let (Some(ai), Some(ui)) = (items(&arg(2)?, ar * ac)?, items(&arg(6)?, ur * uc)?) else { return Ok(None) };
    let mut space: Vec<Expr> = vec![g.clone()];
    let mut t_vars = Vec::with_capacity(tower.len());
    for v in tower.iter() {
        if !v.is_instance_of::<crate::PyExprNode>() {
            return Ok(None);
        }
        let e = crate::to_expr(&v)?;
        if !plain_symbol(&e) {
            return Ok(None);
        }
        t_vars.push(e);
    }
    let nd = prefix(derivations.len(), derivations.len() as isize + level + 1);
    let nt = prefix(t_vars.len(), t_vars.len() as isize + level + 1);
    let mut d_exprs = Vec::new();
    for d in derivations.iter().take(nd.min(nt)) {
        if !d.is_exact_instance_of::<PyPoly>() {
            return Ok(None);
        }
        let e = d.downcast::<PyPoly>()?.get().p.as_expr();
        for s in symox::traverse::free_symbols(&e) {
            if !plain_symbol(&s) {
                return Ok(None);
            }
            space.push(s);
        }
        d_exprs.push(e);
    }
    space.extend(t_vars.iter().cloned());
    space.extend(da.flat_parameters());
    space.extend(du.flat_parameters());
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
    let index = |e: &Expr| space.iter().position(|s| s == e);
    let Ok(field) = ExactDomain::functions(false, false, space.clone(), "lex") else { return Ok(None) };
    let _scope = Scope::new(space.len());
    let mut interrupted: Option<PyErr> = None;
    let result = symox::poly::mpoly::checked_operation(|| -> symox::poly::convert::PResult<Option<_>> {
        let convert = |es: &[Expr]| es.iter().map(|e| function(&field, e)).collect::<Option<Vec<_>>>();
        let (Some(av), Some(uv), Some(dv)) = (convert(&ai), convert(&ui), convert(&d_exprs)) else { return Ok(None) };
        let active = t_vars.iter().zip(dv).map(|(v, d)| (index(v).unwrap(), d)).collect();
        let tw = Tower { active, tower: t_vars.iter().map(|v| index(v).unwrap()).collect(), gen: index(&g) };
        let gens = PyTuple::new(py, [gen.clone()]).map_err(|_| symox::poly::convert::PolyError::Unsupported("gens".into()))?;
        let mut host = PyHost { py, gens, space: field.clone() };
        let a = Matrix { rows: ar, cols: ac, values: av, items: ai.clone(), domain: da };
        let u = Matrix { rows: ur, cols: uc, values: uv, items: ui.clone(), domain: du };
        Ok(k::constant_system(&mut host, a, u, &tw, &mut || match py.check_signals() {
            Ok(()) => true,
            Err(e) => { interrupted = Some(e); false }
        }))
    });
    if let Some(e) = interrupted {
        return Err(e);
    }
    py.check_signals()?;
    let Ok(Some((a, u))) = result else { return Ok(None) };
    let out = |m: Matrix<ExactDomain>| -> PyResult<Bound<'_, PyAny>> {
        let items = PyList::new(py, m.items.into_iter().map(crate::wrap))?;
        let d = PyDomain { inner: m.domain }.object(py)?;
        Ok((m.rows, m.cols, items, d).into_pyobject(py)?.into_any())
    };
    Ok(Some((out(a)?, out(u)?).into_pyobject(py)?.into_any().unbind()))
}
