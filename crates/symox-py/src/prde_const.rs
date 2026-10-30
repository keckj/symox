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

pub(crate) fn items(v: &Bound<'_, PyAny>, len: usize) -> PyResult<Option<Vec<Expr>>> {
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
    if args.len() != 8 {
        return Ok(None);
    }
    let (Ok(tower), Ok(derivations)) = (tower.downcast_exact::<PyList>(), derivations.downcast_exact::<PyList>()) else { return Ok(None) };
    if !level.is_exact_instance_of::<PyInt>() {
        return Ok(None);
    }
    let Ok(level) = level.extract::<isize>() else { return Ok(None) };
    let arg = |i: usize| args.get_item(i);
    let (Some(ar), Some(ac), Some(ur), Some(uc)) = (size(&arg(0)?), size(&arg(1)?), size(&arg(4)?), size(&arg(5)?)) else { return Ok(None) };
    let (Some(da), Some(du)) = (domain(&arg(3)?)?, domain(&arg(7)?)?) else { return Ok(None) };
    let (Some(ai), Some(ui)) = (items(&arg(2)?, ar * ac)?, items(&arg(6)?, ur * uc)?) else { return Ok(None) };
    let result = field_system(py, gen, (ar, ac, ai, da), (ur, uc, ui, du), tower, derivations, level, |_, a, u, _| (a, u))?;
    let Some((a, u)) = result else { return Ok(None) };
    let out = |m: Matrix<ExactDomain>| -> PyResult<Bound<'_, PyAny>> {
        let items = PyList::new(py, m.items.into_iter().map(crate::wrap))?;
        let d = PyDomain { inner: m.domain }.object(py)?;
        Ok((m.rows, m.cols, items, d).into_pyobject(py)?.into_any())
    };
    Ok(Some((out(a)?, out(u)?).into_pyobject(py)?.into_any().unbind()))
}

/// The domain of `PolyMatrix(items, gen)` (`parallel_poly_from_expr(items,
/// gen, field=True)`) when it is QQ or a fraction field over plain symbols.
pub(crate) fn infer_domain(py: Python<'_>, gen: &Bound<'_, PyAny>, items: &[Expr]) -> PyResult<Option<ExactDomain>> {
    let gens = PyTuple::new(py, [gen.clone()])?;
    let mut host = PyHost { py, gens, space: ExactDomain::QQ };
    Ok(host.infer(items).filter(admitted))
}

/// `constant_system(A, u, DE)` over the fraction field of every symbol
/// involved, `A = (rows, cols, items, domain)` and the same for `u`. The
/// `finish` callback runs in the rational function scope with the result
/// matrices (values and items) and the derivation of the active levels.
/// `Ok(None)` when the kernel does not reproduce the stock body.
#[allow(clippy::too_many_arguments)]
pub(crate) fn field_system<R>(py: Python<'_>, gen: &Bound<'_, PyAny>, a: (usize, usize, Vec<Expr>, ExactDomain),
    u: (usize, usize, Vec<Expr>, ExactDomain), tower: &Bound<'_, PyList>, derivations: &Bound<'_, PyList>, level: isize,
    finish: impl FnOnce(&mut dyn FnMut(&RF) -> Expr, Matrix<ExactDomain>, Matrix<ExactDomain>, &Tower) -> R) -> PyResult<Option<R>> {
    if !gen.is_instance_of::<crate::PyExprNode>() || level.unsigned_abs() > 1 << 16 {
        return Ok(None);
    }
    let g = crate::to_expr(gen)?;
    if !plain_symbol(&g) {
        return Ok(None);
    }
    let (ar, ac, ai, da) = a;
    let (ur, uc, ui, du) = u;
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
        let solved = k::constant_system(&mut host, a, u, &tw, &mut || match py.check_signals() {
            Ok(()) => true,
            Err(e) => { interrupted = Some(e); false }
        });
        Ok(solved.map(|(a, u)| finish(&mut |v: &RF| host.expr(v), a, u, &tw)))
    });
    if let Some(e) = interrupted {
        return Err(e);
    }
    py.check_signals()?;
    Ok(result.ok().flatten())
}

/// `ZZ`, `QQ`, `ZZ[...]`, `QQ[...]`, `ZZ(...)` or `QQ(...)` over plain
/// symbols: the Poly expressions have no floats, so the constructor's
/// `_zero_float_domain` step is the identity on every term below.
fn exact_ground(d: &ExactDomain) -> bool {
    match d {
        ExactDomain::ZZ | ExactDomain::QQ => true,
        ExactDomain::Functions { parameters, .. } => parameters.iter().all(plain_symbol),
        _ => false,
    }
}

pub(crate) fn exact_poly(p: &Bound<'_, PyAny>) -> PyResult<bool> {
    let Ok(p) = p.downcast_exact::<PyPoly>() else { return Ok(false) };
    Ok(exact_ground(&p.borrow().domain_descriptor()?.inner))
}

/// `derivation(p, DE, coefficientD)` with `basic=False` (risch.py): the
/// stock loop over the active `(D, T)` pairs, with the two Poly
/// constructions (`p.as_poly(v)` and `term.as_poly(t)`) run by the native
/// constructor `Poly._from_normalized` that `symox._poly_constructor` calls
/// for these arguments, and the Poly/expression `diff`, `as_expr`, `*` and
/// `+` of the stock body. The term expressions keep the stock unexpanded
/// form, so each term Poly has the stock domain and the sum unifies them
/// in the stock order, starting from `Poly(0, t)` over ZZ.
/// `Ok(None)` requests the Python body (any failure it would report).
pub(crate) fn derivation_poly(py: Python<'_>, args: &Bound<'_, PyTuple>, t: &Bound<'_, PyAny>,
    tower: &Bound<'_, PyAny>, derivations: &Bound<'_, PyAny>, level: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    if args.len() != 2 || !t.is_instance_of::<crate::PyExprNode>() {
        return Ok(None);
    }
    let p = args.get_item(0)?;
    let coefficient_d = args.get_item(1)?.is_truthy()?;
    let (Ok(tower), Ok(derivations)) = (tower.downcast_exact::<PyList>(), derivations.downcast_exact::<PyList>()) else { return Ok(None) };
    if !level.is_exact_instance_of::<PyInt>() || !exact_poly(&p)? {
        return Ok(None);
    }
    let Ok(mut level) = level.extract::<isize>() else { return Ok(None) };
    let gen_t = PyTuple::new(py, [t])?;
    let zero = crate::wrap(Expr::zero()).into_pyobject(py)?.into_any();
    let mut r = PyPoly::_from_normalized(&zero, &gen_t, None, None, None)?.into_pyobject(py)?.into_any();
    let (nd, nt) = (derivations.len() as isize, tower.len() as isize);
    if coefficient_d {
        if level <= -nt {
            return Ok(Some(r.unbind()));
        }
        level -= 1;
    }
    let (sd, st) = (nd + level + 1, nt + level + 1);
    if sd < 0 || st < 0 || sd > nd || st > nt {
        return Ok(None);
    }
    let errors = py.import("symox.polys.polyerrors")?;
    let (poly_error, gens_needed) = (errors.getattr("PolynomialError")?, errors.getattr("GeneratorsNeeded")?);
    for i in 0..(sd.min(st) as usize) {
        let (d, v) = (derivations.get_item(i)?, tower.get_item(i)?);
        if !v.is_instance_of::<crate::PyExprNode>() || !exact_poly(&d)? {
            return Ok(None);
        }
        let gen_v = PyTuple::new(py, [&v])?;
        let derivative = match PyPoly::_from_normalized(&p, &gen_v, None, None, None) {
            Ok(pv) => pv.into_pyobject(py)?.into_any().call_method1("diff", (&v,))?.call_method0("as_expr")?,
            Err(e) if e.is_instance(py, &poly_error) || e.is_instance(py, &gens_needed) =>
                p.call_method0("as_expr")?.call_method1("diff", (&v,))?,
            Err(_) => return Ok(None),
        };
        let term = d.call_method0("as_expr")?.mul(&derivative)?;
        if !term.is_instance_of::<crate::PyExprNode>() {
            return Ok(None);
        }
        let Ok(term) = PyPoly::_from_normalized(&term, &gen_t, None, None, None) else { return Ok(None) };
        r = r.add(term.into_pyobject(py)?)?;
    }
    Ok(Some(r.unbind()))
}
