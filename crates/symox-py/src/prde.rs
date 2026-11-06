//! Base-case parametric Risch DE kernels (`symox.integrals.prde`).
//!
//! Admission only: exact stock `Poly` instances in the single active
//! variable `x` with `D = [Poly(1, x)]`, plain ZZ/QQ coefficients and bounded
//! sizes. Anything else returns `None` and the Python body runs unchanged.
use pyo3::prelude::*;
use pyo3::types::{PyInt, PyList, PyTuple};
use symox::expr::{Expr, Kind};
use symox::integrate::prde::{self as k, Mat, P, ParamRisch};
use symox::integrate::prde_log as log;
use symox::poly::{dense, polyobj::Domain, Rat};
use crate::poly::PyPoly;

fn plain(p: &symox::poly::fraction::Polynomial) -> bool {
    p.field.is_none() && p.algebraic.is_none() && p.expression.is_none()
        && p.base.modulus.is_none() && !p.base.has_domain_gens() && p.base.ngens == 1
        && p.base.gens.len() == 1 && matches!(p.base.domain, Domain::ZZ | Domain::QQ)
        && p.base.rep.degree(0) <= 64
        && p.base.rep.iter().all(|(_, v)| v.n.bits() <= 4096 && v.d.bits() <= 4096)
}

fn symbol(value: &Bound<'_, PyAny>) -> PyResult<Option<Expr>> {
    if !value.is_instance_of::<crate::PyExprNode>() {
        return Ok(None);
    }
    let e = crate::to_expr(value)?;
    let Kind::Sym(s) = e.kind() else { return Ok(None) };
    if s.host().is_some() || s.scalar_identity().is_some() || s.is_wild()
        || symox::decompose::commutativity(&e) != Some(true)
        || !crate::node_type(value.py(), &e).is_some_and(|c| value.get_type().is(&c)) {
        return Ok(None);
    }
    Ok(Some(e))
}

pub(crate) fn stock_poly_type(py: Python<'_>) -> bool {
    use pyo3::impl_::pyclass::PyClassImpl;
    let expected = <PyPoly as PyClassImpl>::items_iter().flat_map(|i| i.slots)
        .find(|slot| slot.slot == pyo3::ffi::Py_tp_new).map(|slot| slot.pfunc);
    let current = unsafe { pyo3::ffi::PyType_GetSlot(py.get_type::<PyPoly>().as_type_ptr(), pyo3::ffi::Py_tp_new) };
    expected.is_some_and(|p| !p.is_null() && p == current)
}

struct Ctx {
    t: Expr,
    template: symox::poly::polyobj::Poly,
}

impl Ctx {
    fn poly(&self, v: &Bound<'_, PyAny>) -> PyResult<Option<P>> {
        if !v.is_exact_instance_of::<PyPoly>() {
            return Ok(None);
        }
        let p = v.downcast::<PyPoly>()?.get();
        if !plain(&p.p) || p.p.base.gens != [self.t.clone()] {
            return Ok(None);
        }
        Ok(Some(P { c: dense::from_mpoly(&p.p.base.rep, 0), qq: p.p.base.domain == Domain::QQ }))
    }
    fn polys(&self, v: &Bound<'_, PyAny>) -> PyResult<Option<Vec<P>>> {
        let items: Vec<Bound<'_, PyAny>> = if let Ok(l) = v.downcast_exact::<PyList>() {
            l.iter().collect()
        } else if let Ok(t) = v.downcast_exact::<PyTuple>() {
            t.iter().collect()
        } else {
            return Ok(None);
        };
        if items.is_empty() || items.len() > 64 {
            return Ok(None);
        }
        let mut out = Vec::with_capacity(items.len());
        for item in &items {
            let Some(p) = self.poly(item)? else { return Ok(None) };
            out.push(p);
        }
        Ok(Some(out))
    }
    fn pairs(&self, v: &Bound<'_, PyAny>) -> PyResult<Option<Vec<(P, P)>>> {
        let Ok(l) = v.downcast_exact::<PyList>() else { return Ok(None) };
        if l.is_empty() {
            return Ok(None);
        }
        self.pairs_or_empty(v)
    }
    /// `pairs` that also admits the empty list (`limited_integrate`'s `G`).
    fn pairs_or_empty(&self, v: &Bound<'_, PyAny>) -> PyResult<Option<Vec<(P, P)>>> {
        let Ok(l) = v.downcast_exact::<PyList>() else { return Ok(None) };
        if l.len() > 64 {
            return Ok(None);
        }
        let mut out = Vec::with_capacity(l.len());
        for item in l.iter() {
            let Ok(t) = item.downcast_exact::<PyTuple>() else { return Ok(None) };
            if t.len() != 2 {
                return Ok(None);
            }
            let (Some(a), Some(d)) = (self.poly(&t.get_item(0)?)?, self.poly(&t.get_item(1)?)?) else { return Ok(None) };
            out.push((a, d));
        }
        Ok(Some(out))
    }
    fn out(&self, py: Python<'_>, p: &P) -> PyResult<PyObject> {
        let mut base = self.template.clone();
        base.rep = dense::to_mpoly(&p.c, 0, 1);
        base.domain = if p.qq { Domain::QQ } else { Domain::ZZ };
        Ok(PyPoly { p: base.into() }.into_pyobject(py)?.into_any().unbind())
    }
    fn outs(&self, py: Python<'_>, ps: &[P]) -> PyResult<Vec<PyObject>> {
        ps.iter().map(|p| self.out(py, p)).collect()
    }
    fn pair_list(&self, py: Python<'_>, ps: &[(P, P)]) -> PyResult<PyObject> {
        let items = ps.iter().map(|(a, d)| Ok(PyTuple::new(py, [self.out(py, a)?, self.out(py, d)?])?.into_any()))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyList::new(py, items)?.into_any().unbind())
    }
}

fn int_expr(n: symox::poly::Int) -> Expr {
    symox::poly::convert::rat_to_expr(&Rat::from(n))
}

/// `(n, u)` of `is_log_deriv_k_t_radical_in_field`, or `None`.
fn log_answer(py: Python<'_>, answer: Option<(symox::poly::Int, Expr)>) -> PyResult<PyObject> {
    Ok(match answer {
        None => py.None(),
        Some((n, u)) => (crate::wrap(int_expr(n)), crate::wrap(u)).into_pyobject(py)?.into_any().unbind(),
    })
}

fn mat_out(py: Python<'_>, m: &Mat) -> PyResult<PyObject> {
    let items = m.e.iter().map(|r| Ok(crate::wrap(symox::poly::convert::rat_to_expr(r)).into_pyobject(py)?.into_any()))
        .collect::<PyResult<Vec<_>>>()?;
    Ok((m.rows, m.cols, PyList::new(py, items)?).into_pyobject(py)?.into_any().unbind())
}

fn mat_in(rows: &Bound<'_, PyAny>, cols: &Bound<'_, PyAny>, items: &Bound<'_, PyAny>) -> PyResult<Option<Mat>> {
    if !rows.is_exact_instance_of::<PyInt>() || !cols.is_exact_instance_of::<PyInt>() {
        return Ok(None);
    }
    let (Ok(r), Ok(c)) = (rows.extract::<usize>(), cols.extract::<usize>()) else { return Ok(None) };
    let Ok(items) = items.downcast_exact::<PyList>() else { return Ok(None) };
    if r > 256 || c > 256 || items.len() != r * c {
        return Ok(None);
    }
    let mut e = Vec::with_capacity(r * c);
    for item in items.iter() {
        if !item.is_instance_of::<crate::PyExprNode>() {
            return Ok(None);
        }
        let x = crate::to_expr(&item)?;
        let Some(v) = x.as_number().and_then(Rat::from_number) else { return Ok(None) };
        e.push(v);
    }
    Ok(Some(Mat { rows: r, cols: c, e }))
}

fn int_arg(v: &Bound<'_, PyAny>) -> Option<i64> {
    if !v.is_exact_instance_of::<PyInt>() {
        return None;
    }
    v.extract::<i64>().ok().filter(|n| (-(1 << 20)..=(1 << 20)).contains(n))
}

/// Validates the base-case derivation `D = [Poly(1, x)]` at the active level.
fn context(t: &Bound<'_, PyAny>, d: &Bound<'_, PyAny>, variables: &Bound<'_, PyAny>,
    derivations: &Bound<'_, PyAny>, level: &Bound<'_, PyAny>) -> PyResult<Option<Ctx>> {
    let py = t.py();
    if !stock_poly_type(py) {
        return Ok(None);
    }
    if variables.is_none() {
        // poly_linear_constraints: no derivation, only the generator.
        let Some(t) = symbol(t)? else { return Ok(None) };
        if !d.is_exact_instance_of::<PyPoly>() {
            return Ok(None);
        }
        let dp = d.downcast::<PyPoly>()?.get();
        if !plain(&dp.p) || dp.p.base.gens != [t.clone()] {
            return Ok(None);
        }
        let template = dp.p.base.clone();
        return Ok(Some(Ctx { t, template }));
    }
    if !d.is_exact_instance_of::<PyPoly>() || !level.is_exact_instance_of::<PyInt>() {
        return Ok(None);
    }
    let (Ok(variables), Ok(derivations)) = (variables.downcast_exact::<PyList>(), derivations.downcast_exact::<PyList>()) else { return Ok(None) };
    let Ok(level) = level.extract::<isize>() else { return Ok(None) };
    if variables.is_empty() || variables.len() != derivations.len() || level >= 0 || level < -(variables.len() as isize) {
        return Ok(None);
    }
    let active = (variables.len() as isize + level + 1) as usize;
    if active != 1 || !t.is(&variables.get_item(0)?) || !d.is(&derivations.get_item(0)?) {
        return Ok(None);
    }
    let Some(t) = symbol(t)? else { return Ok(None) };
    let dp = d.downcast::<PyPoly>()?.get();
    if !plain(&dp.p) || dp.p.base.gens != [t.clone()] || !dp.p.base.rep.is_one() {
        return Ok(None);
    }
    let template = dp.p.base.clone();
    Ok(Some(Ctx { t, template }))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub(crate) fn _prde_kernel(kind: &str, args: &Bound<'_, PyTuple>, t: &Bound<'_, PyAny>, d: &Bound<'_, PyAny>,
    variables: &Bound<'_, PyAny>, derivations: &Bound<'_, PyAny>, level: &Bound<'_, PyAny>,
    py: Python<'_>) -> PyResult<Option<PyObject>> {
    py.check_signals()?;
    if kind == "constant_system" {
        // Any DE case: rational entries never involve `DE.T` (see admission).
        return constant_system_any(py, args, t, variables);
    }
    if kind == "derivation" {
        return crate::prde_const::derivation_poly(py, args, t, variables, derivations, level);
    }
    if kind == "is_deriv_k" {
        return crate::prde_deriv::is_deriv_k(py, args, t, variables, derivations, level);
    }
    if kind == "constant_system_field" {
        return crate::prde_const::constant_system_field(py, args, t, variables, derivations, level);
    }
    if kind == "is_log_deriv_tower" || kind == "heu_tower" {
        if !stock_poly_type(py) {
            return Ok(None);
        }
        return crate::prde_tower::kernel(py, kind, args, t, d, variables, derivations, level);
    }
    if matches!(kind, "param_rischDE_tower" | "param_poly_rischDE_tower" | "prde_cancel_liouvillian_tower" | "limited_integrate_tower") {
        if !stock_poly_type(py) {
            return Ok(None);
        }
        return crate::prde_param::kernel(py, kind, args, t, d, variables, derivations, level);
    }
    let Some(cx) = context(t, d, variables, derivations, level)? else { return Ok(None) };
    let arg = |i: usize| args.get_item(i);
    let result: Option<PyObject> = match (kind, args.len()) {
        ("param_rischDE", 3) => {
            let (Some(fa), Some(fd), Some(g)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?, cx.pairs(&arg(2)?)?) else { return Ok(None) };
            match k::param_rischde(&fa, &fd, &g) {
                None => None,
                Some(ParamRisch::NoSolution(m)) => Some((PyList::empty(py), mat_out(py, &Mat::eye(m))?).into_pyobject(py)?.into_any().unbind()),
                Some(ParamRisch::Solution(h, c)) => Some((cx.pair_list(py, &h)?, mat_out(py, &c)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("param_poly_rischDE", 4) => {
            let (Some(a), Some(b), Some(q), Some(n)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?, cx.polys(&arg(2)?)?, int_arg(&arg(3)?)) else { return Ok(None) };
            match k::param_poly_rischde(&a, &b, &q, n, 0) {
                None => None,
                Some((h, m)) => Some((PyList::new(py, cx.outs(py, &h)?)?, mat_out(py, &m)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("prde_no_cancel_b_large" | "prde_no_cancel_b_small", 3) => {
            let (Some(b), Some(mut q), Some(n)) = (cx.poly(&arg(0)?)?, cx.polys(&arg(1)?)?, int_arg(&arg(2)?)) else { return Ok(None) };
            let r = if kind == "prde_no_cancel_b_large" { k::prde_no_cancel_b_large(&b, &mut q, n) } else { k::prde_no_cancel_b_small(&b, &mut q, n) };
            match r {
                None => None,
                Some((h, m)) => Some((PyList::new(py, cx.outs(py, &h)?)?, mat_out(py, &m)?, PyList::new(py, cx.outs(py, &q)?)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("prde_spde", 4) => {
            let (Some(a), Some(b), Some(q), Some(n)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?, cx.polys(&arg(2)?)?, int_arg(&arg(3)?)) else { return Ok(None) };
            if a.is_zero() {
                return Ok(None);
            }
            match k::prde_spde(&a, &b, &q, n) {
                None => None,
                Some(s) => Some((cx.out(py, &s.b)?, PyList::new(py, cx.outs(py, &s.q)?)?, PyList::new(py, cx.outs(py, &s.r)?)?, s.n).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("prde_normal_denom", 3) => {
            let (Some(fa), Some(fd), Some(g)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?, cx.pairs(&arg(2)?)?) else { return Ok(None) };
            match k::prde_normal_denom(&fa, &fd, &g) {
                None => None,
                Some(r) => Some((cx.out(py, &r.a)?, (cx.out(py, &r.ba)?, cx.out(py, &r.bd)?), cx.pair_list(py, &r.g)?, cx.out(py, &r.h)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("prde_special_denom", 2) => {
            let (Some(ba), Some(bd)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?) else { return Ok(None) };
            match k::prde_special_denom(&ba, &bd) {
                None => None,
                Some((b, h)) => Some((cx.out(py, &b)?, cx.out(py, &h)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("prde_linear_constraints", 1) => {
            let Some(g) = cx.pairs(&arg(0)?)? else { return Ok(None) };
            match k::prde_linear_constraints(&g) {
                None => None,
                Some((q, m)) => Some((PyTuple::new(py, cx.outs(py, &q)?)?, mat_out(py, &m)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("poly_linear_constraints", 1) => {
            let (Some(p), Some(dd)) = (cx.polys(&arg(0)?)?, cx.poly(d)?) else { return Ok(None) };
            match k::poly_linear_constraints(&p, &dd) {
                None => None,
                Some((q, m)) => Some((PyTuple::new(py, cx.outs(py, &q)?)?, mat_out(py, &m)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("limited_integrate_reduce", 3) => {
            let (Some(fa), Some(fd), Some(g)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?, cx.pairs(&arg(2)?)?) else { return Ok(None) };
            match symox::integrate::prde_limited::limited_integrate_reduce(&fa, &fd, &g) {
                None => None,
                Some(r) => Some((cx.out(py, &r.a)?, cx.out(py, &r.b)?, r.n, (cx.out(py, &r.g.0)?, cx.out(py, &r.g.1)?), cx.pair_list(py, &r.v)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("is_log_deriv_in_field", 2) => {
            let (Some(fa), Some(fd)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?) else { return Ok(None) };
            let mut trace = Vec::new();
            match log::is_log_deriv_traced(&fa, &fd, &cx.t, &mut trace) {
                None => None,
                Some(answer) => Some((log_answer(py, answer)?, PyList::new(py, trace)?).into_pyobject(py)?.into_any().unbind()),
            }
        }
        ("parametric_log_deriv_heu", 4) => {
            let (Some(fa), Some(fd), Some(wa), Some(wd)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?, cx.poly(&arg(2)?)?, cx.poly(&arg(3)?)?) else { return Ok(None) };
            if fd.is_zero() || wd.is_zero() {
                return Ok(None);
            }
            match log::parametric_log_deriv_heu(&fa, &fd, &wa, &wd, &cx.t) {
                None => None,
                Some(heu) => {
                    let answer: PyObject = match heu.answer {
                        None => py.None(),
                        Some((n, m, v)) => (crate::wrap(int_expr(n)), crate::wrap(int_expr(m)), crate::wrap(v)).into_pyobject(py)?.into_any().unbind(),
                    };
                    Some((answer, PyList::new(py, heu.trace)?).into_pyobject(py)?.into_any().unbind())
                }
            }
        }
        ("limited_integrate", 3) => {
            let (Some(fa), Some(fd), Some(g)) = (cx.poly(&arg(0)?)?, cx.poly(&arg(1)?)?, cx.pairs_or_empty(&arg(2)?)?) else { return Ok(None) };
            match log::limited_integrate(&fa, &fd, &g, &cx.t) {
                None => None,
                Some(None) => Some((py.None(),).into_pyobject(py)?.into_any().unbind()),
                Some(Some(((ya, yd), cs))) => {
                    let c = cs.into_iter().map(|v| cx.out(py, &P::new(vec![v], true))).collect::<PyResult<Vec<_>>>()?;
                    let y = (cx.out(py, &ya)?, cx.out(py, &yd)?);
                    Some(((y, PyList::new(py, c)?),).into_pyobject(py)?.into_any().unbind())
                }
            }
        }
        _ => None,
    };
    py.check_signals()?;
    Ok(result)
}

/// `constant_system(A, u, DE)` for QQ matrices with rational entries in any
/// DE state: `t` is the matrix generator and `variables` is `DE.T`; both must
/// be plain stock symbols so `rref`'s ground check and `has(*DE.T)` behave
/// as proved (ground entries, no derivation step).
pub(crate) fn constant_system_any(py: Python<'_>, args: &Bound<'_, PyTuple>, gen: &Bound<'_, PyAny>,
    tower: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    if args.len() != 6 || symbol(gen)?.is_none() {
        return Ok(None);
    }
    let Ok(tower) = tower.downcast_exact::<PyList>() else { return Ok(None) };
    for v in tower.iter() {
        if symbol(&v)?.is_none() {
            return Ok(None);
        }
    }
    let arg = |i: usize| args.get_item(i);
    let (Some(a), Some(u)) = (mat_in(&arg(0)?, &arg(1)?, &arg(2)?)?, mat_in(&arg(3)?, &arg(4)?, &arg(5)?)?) else { return Ok(None) };
    if a.rows * a.cols == 0 {
        return Ok(None);
    }
    let Some((a, u)) = k::constant_system(&a, &u) else { return Ok(None) };
    let out = (mat_out(py, &a)?, mat_out(py, &u)?).into_pyobject(py)?.into_any().unbind();
    py.check_signals()?;
    Ok(Some(out))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_prde_kernel, m)?)
}
