//! The parametric Risch differential equation above the base level:
//! `param_rischDE`, `param_poly_rischDE`, `prde_cancel_liouvillian` and
//! `limited_integrate` (prde.py) on exp/primitive towers, with their
//! substeps (`weak_normalizer`, `splitfactor`, `gcdex_diophantine`,
//! `prde_normal_denom`, `prde_special_denom`, `prde_linear_constraints`,
//! `poly_linear_constraints`, `bound_degree(parametric=True)`, `prde_spde`,
//! `prde_no_cancel_b_large`, `constant_system`, the PolyMatrix operations).
//!
//! The stock bodies are transliterated step by step on native Poly objects
//! (`symox._core.Poly` methods and operators), with the native derivation,
//! `constant_system` (QQ and fraction-field kernels), PolyMatrix nullspace
//! and the native lower-level kernels: the base-case `param_rischDE`,
//! `parametric_log_deriv_heu` and `is_log_deriv_k_t_radical_in_field`
//! (`prde.rs`) and their exp/primitive tower forms (`prde_tower.rs`).
//!
//! Every step that the stock body would leave through an exception, every
//! case it does not implement (tan, nonlinear, `b_small`, `b_equal`), every
//! lower kernel that declines and every unproved PolyMatrix domain returns
//! `None`, so the stock Python body runs and reports its own answer or
//! error. The stock `NotImplementedError`/`ValueError` of the lower tower
//! kernels are reported with status 1/2 (`bound_degree` catches the first,
//! as the stock body does). `trace` lists the Dummies the stock body
//! creates, in its order.
use std::cell::{Cell, RefCell};
use pyo3::exceptions::PyException;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyInt, PyList, PyString, PyTuple};
use symox::expr::Expr;
use crate::poly::domains::PyDomain;
use crate::poly::PyPoly;

type O<'py> = Bound<'py, PyAny>;

enum Stop {
    /// The stock body must run (failure, unsupported or unproved step).
    Decline,
    /// A stock exception of a lower kernel: 1 NotImplementedError, 2 ValueError.
    Raise(u8, String),
    /// KeyboardInterrupt and other non-Exception errors.
    Interrupt(PyErr),
}

type R<T> = Result<T, Stop>;

thread_local! {
    /// Where the last decline came from (`SYMOX_PRDE_DEBUG` diagnostics).
    static WHERE: RefCell<String> = const { RefCell::new(String::new()) };
}

#[track_caller]
fn note(what: &str) {
    let at = std::panic::Location::caller();
    WHERE.with(|w| *w.borrow_mut() = format!("{}:{} {}", at.file(), at.line(), what));
}

impl From<PyErr> for Stop {
    #[track_caller]
    fn from(e: PyErr) -> Self {
        let plain = Python::with_gil(|py| {
            let plain = e.is_instance_of::<PyException>(py);
            if plain {
                note(&e.to_string());
            }
            plain
        });
        if plain { Stop::Decline } else { Stop::Interrupt(e) }
    }
}

impl From<std::convert::Infallible> for Stop {
    fn from(e: std::convert::Infallible) -> Self {
        match e {}
    }
}

impl<'a, 'py> From<pyo3::DowncastError<'a, 'py>> for Stop {
    #[track_caller]
    fn from(e: pyo3::DowncastError<'a, 'py>) -> Self {
        note(&e.to_string());
        Stop::Decline
    }
}

impl<'py> From<pyo3::DowncastIntoError<'py>> for Stop {
    #[track_caller]
    fn from(e: pyo3::DowncastIntoError<'py>) -> Self {
        note(&e.to_string());
        Stop::Decline
    }
}

#[track_caller]
fn decline<T>() -> R<T> {
    note("");
    Err(Stop::Decline)
}

/// `PolyMatrix` state: `rows x cols` expressions, generators and domain.
#[derive(Clone)]
struct PM<'py> {
    r: usize,
    c: usize,
    e: Vec<O<'py>>,
    gens: Bound<'py, PyTuple>,
    dom: O<'py>,
}

const STEPS: u64 = 1_000_000;

struct Cx<'py> {
    py: Python<'py>,
    qq: O<'py>,
    tl: Bound<'py, PyList>,
    dl: Bound<'py, PyList>,
    cases: Vec<String>,
    cases_obj: Bound<'py, PyList>,
    nt: usize,
    act: Cell<usize>,
    trace: RefCell<Vec<String>>,
    steps: Cell<u64>,
    z: O<'py>,
    inc: Bound<'py, PyDict>,
    fld: Bound<'py, PyDict>,
}

impl<'py> Cx<'py> {
    fn tick(&self) -> R<()> {
        let s = self.steps.get() + 1;
        self.steps.set(s);
        if s > STEPS {
            return decline();
        }
        if s % 32 == 0 {
            self.py.check_signals().map_err(Stop::Interrupt)?;
        }
        Ok(())
    }
    fn push(&self, name: &str) {
        self.trace.borrow_mut().push(name.to_string());
    }
    fn t(&self) -> R<O<'py>> {
        Ok(self.tl.get_item(self.act.get() - 1)?)
    }
    fn d(&self) -> R<O<'py>> {
        Ok(self.dl.get_item(self.act.get() - 1)?)
    }
    fn level(&self) -> R<O<'py>> {
        let l = self.act.get() as isize - self.nt as isize - 1;
        Ok(l.into_pyobject(self.py)?.into_any())
    }
    fn case(&self) -> String {
        self.cases[self.act.get() - 1].clone()
    }
    /// The body of `with DecrementLevel(DE)`.
    fn lower<T>(&self, f: impl FnOnce(&Self) -> R<T>) -> R<T> {
        let a = self.act.get();
        if a < 2 {
            return decline();
        }
        self.act.set(a - 1);
        let r = f(self);
        self.act.set(a);
        r
    }
    fn gt(&self) -> R<Bound<'py, PyTuple>> {
        Ok(PyTuple::new(self.py, [self.t()?])?)
    }
    fn int(&self, i: i64) -> R<O<'py>> {
        Ok(crate::wrap(Expr::int(i)).into_pyobject(self.py)?.into_any())
    }
    fn pyint(&self, i: i64) -> R<O<'py>> {
        Ok(i.into_pyobject(self.py)?.into_any())
    }
    /// `Poly(rep, *gens)`.
    fn mk(&self, rep: &O<'py>, gens: &Bound<'py, PyTuple>) -> R<O<'py>> {
        Ok(PyPoly::_from_normalized(rep, gens, None, None, None)?.into_pyobject(self.py)?.into_any())
    }
    /// `Poly(rep, *gens, **kwargs)` through the stock constructor.
    fn cls_call(&self, args: Vec<O<'py>>, kwargs: Option<&Bound<'py, PyDict>>) -> R<O<'py>> {
        let args = PyTuple::new(self.py, args)?;
        Ok(self.py.get_type::<PyPoly>().call(args, kwargs)?)
    }
    /// `Poly(rep, DE.t)`.
    fn pt(&self, rep: &O<'py>) -> R<O<'py>> {
        self.mk(rep, &self.gt()?)
    }
    fn pti(&self, i: i64) -> R<O<'py>> {
        self.pt(&self.int(i)?)
    }
    fn lc(&self, p: &O<'py>) -> R<O<'py>> {
        Ok(p.call_method0("LC")?)
    }
    fn nth(&self, p: &O<'py>, i: i64) -> R<O<'py>> {
        Ok(p.call_method1("nth", (i,))?)
    }
    fn ex(&self, p: &O<'py>) -> R<O<'py>> {
        Ok(p.call_method0("as_expr")?)
    }
    fn is_zero(&self, p: &O<'py>) -> R<bool> {
        Ok(p.getattr("is_zero")?.is_truthy()?)
    }
    fn ground(&self, p: &O<'py>) -> R<bool> {
        Ok(p.getattr("is_ground")?.is_truthy()?)
    }
    /// `p.degree()`, `None` for `-oo`.
    fn deg(&self, p: &O<'py>) -> R<Option<i64>> {
        if self.is_zero(p)? {
            return Ok(None);
        }
        Ok(Some(p.call_method0("degree")?.extract::<i64>()?))
    }
    fn m1(&self, p: &O<'py>, name: &str, o: &O<'py>) -> R<O<'py>> {
        self.tick()?;
        Ok(p.call_method1(name, (o,))?)
    }
    fn gcd(&self, a: &O<'py>, b: &O<'py>) -> R<O<'py>> {
        self.m1(a, "gcd", b)
    }
    fn lcm(&self, a: &O<'py>, b: &O<'py>) -> R<O<'py>> {
        self.m1(a, "lcm", b)
    }
    fn quo(&self, a: &O<'py>, b: &O<'py>) -> R<O<'py>> {
        self.m1(a, "quo", b)
    }
    fn exquo(&self, a: &O<'py>, b: &O<'py>) -> R<O<'py>> {
        self.m1(a, "exquo", b)
    }
    fn div(&self, a: &O<'py>, b: &O<'py>) -> R<(O<'py>, O<'py>)> {
        let r = self.m1(a, "div", b)?;
        Ok((r.get_item(0)?, r.get_item(1)?))
    }
    fn diff(&self, p: &O<'py>) -> R<O<'py>> {
        self.m1(p, "diff", &self.t()?)
    }
    fn to_field(&self, p: &O<'py>) -> R<O<'py>> {
        Ok(p.call_method0("to_field")?)
    }
    fn monic(&self, p: &O<'py>) -> R<O<'py>> {
        Ok(p.call_method0("monic")?)
    }
    /// `a.cancel(b, include=True)`.
    fn cancel(&self, a: &O<'py>, b: &O<'py>) -> R<(O<'py>, O<'py>)> {
        self.tick()?;
        let r = a.call_method("cancel", (b,), Some(&self.inc))?;
        Ok((r.get_item(0)?, r.get_item(1)?))
    }
    /// `derivation(p, DE)` (native stock loop).
    fn deriv(&self, p: &O<'py>) -> R<O<'py>> {
        self.tick()?;
        let args = PyTuple::new(self.py, [p.clone(), false.into_pyobject(self.py)?.to_owned().into_any()])?;
        let r = crate::prde_const::derivation_poly(self.py, &args, &self.t()?, self.tl.as_any(), self.dl.as_any(), &self.level()?)?;
        match r {
            Some(r) => Ok(r.into_bound(self.py)),
            None => decline(),
        }
    }
    /// The stock `cancel(expr)`.
    fn ecancel(&self, e: &O<'py>) -> R<O<'py>> {
        self.tick()?;
        Ok(crate::poly::cancel(e, &PyTuple::empty(self.py), None)?.into_bound(self.py))
    }
    /// `frac_in(f, gen, **kwargs)` (`field=True` when `field`).
    fn frac_in(&self, f: &O<'py>, gen: &O<'py>, field: bool) -> R<(O<'py>, O<'py>)> {
        self.tick()?;
        let e = if f.is_instance_of::<PyPoly>() { self.ex(f)? } else { f.clone() };
        let nd = e.call_method0("as_expr")?.call_method0("as_numer_denom")?;
        let (n, d) = (nd.get_item(0)?, nd.get_item(1)?);
        let kw = if field { Some(&self.fld) } else { None };
        let a = n.call_method("as_poly", (gen,), kw)?;
        let b = d.call_method("as_poly", (gen,), kw)?;
        if a.is_none() || b.is_none() {
            return decline();
        }
        Ok((a, b))
    }
    /// `int(e)` of an Integer answer, otherwise the stock body.
    fn integer(&self, e: &O<'py>) -> R<i64> {
        if !e.getattr("is_Integer")?.is_truthy()? {
            return decline();
        }
        Ok(e.call_method0("__int__")?.extract::<i64>()?)
    }

    // ---- PolyMatrix ----------------------------------------------------

    fn pm(&self, r: usize, c: usize, e: Vec<O<'py>>, gens: Bound<'py, PyTuple>, dom: O<'py>) -> PM<'py> {
        PM { r, c, e, gens, dom }
    }
    fn norm_dom(&self, d: O<'py>) -> R<O<'py>> {
        if !d.is(&self.qq) && d.eq(&self.qq)? { Ok(self.qq.clone()) } else { Ok(d) }
    }
    fn eye(&self, n: usize) -> R<PM<'py>> {
        let (one, zero) = (self.int(1)?, self.int(0)?);
        let mut e = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                e.push(if i == j { one.clone() } else { zero.clone() });
            }
        }
        Ok(self.pm(n, n, e, self.gt()?, self.qq.clone()))
    }
    fn zeros(&self, r: usize, c: usize) -> R<PM<'py>> {
        let zero = self.int(0)?;
        Ok(self.pm(r, c, vec![zero; r * c], self.gt()?, self.qq.clone()))
    }
    /// `Matrix(rows, cols, items)` of Polys: the domain and generators of
    /// the unified items (`p.unify(p2)` fold); no items: `QQ`, `(DE.t,)`.
    fn from_polys(&self, r: usize, c: usize, items: Vec<O<'py>>) -> R<PM<'py>> {
        if items.is_empty() {
            return Ok(self.pm(r, c, vec![], self.gt()?, self.qq.clone()));
        }
        let mut p = items[0].clone();
        let (mut pg, mut pd) = (p.getattr("gens")?, p.getattr("domain")?);
        for q in &items[1..] {
            self.tick()?;
            // Same generators and domain: `unify` leaves both unchanged
            // (only they are read from the fold).
            if q.getattr("gens")?.eq(&pg)? && q.getattr("domain")?.eq(&pd)? {
                continue;
            }
            p = p.call_method1("unify", (q,))?.get_item(0)?;
            (pg, pd) = (p.getattr("gens")?, p.getattr("domain")?);
        }
        let gens = pg.downcast_into::<PyTuple>()?;
        if !gens.as_any().eq(self.gt()?)? {
            return decline();
        }
        let dom = self.norm_dom(p.getattr("domain")?)?;
        let e = items.iter().map(|x| self.ex(x)).collect::<R<Vec<_>>>()?;
        Ok(self.pm(r, c, e, gens, dom))
    }
    /// `Matrix(rows, cols, exprs, DE.t)`: `parallel_poly_from_expr(exprs,
    /// t, field=True)` domain; no items: `QQ`.
    fn from_exprs(&self, r: usize, c: usize, items: Vec<O<'py>>) -> R<PM<'py>> {
        let gens = self.gt()?;
        if items.is_empty() {
            return Ok(self.pm(r, c, vec![], gens, self.qq.clone()));
        }
        let mut exprs = Vec::with_capacity(items.len());
        for i in &items {
            if !i.is_instance_of::<crate::PyExprNode>() {
                return decline();
            }
            exprs.push(crate::to_expr(i)?);
        }
        let Some(d) = crate::prde_const::infer_domain(self.py, &self.t()?, &exprs)? else { return decline() };
        let dom = self.norm_dom(PyDomain { inner: d }.object(self.py)?.into_bound(self.py))?;
        Ok(self.pm(r, c, items, gens, dom))
    }
    fn getitem(&self, m: &PM<'py>, i: usize) -> R<O<'py>> {
        let Some(e) = m.e.get(i) else { return decline() };
        // An exact rational over QQ or a fraction field of plain symbols:
        // the stock constructor only parses the options before the native
        // `_from_normalized` (as `_polymatrix_rational.ground_poly`).
        if e.is_instance_of::<crate::PyExprNode>() {
            if let Ok(d) = m.dom.downcast::<PyDomain>() {
                let inner = &d.get().inner;
                let shared = match inner {
                    symox::poly::domains::ExactDomain::Functions { parameters, .. } => {
                        let mut s = false;
                        for g in m.gens.iter() {
                            s |= !g.is_instance_of::<crate::PyExprNode>() || parameters.contains(&crate::to_expr(&g)?);
                        }
                        s
                    }
                    _ => false,
                };
                if !shared && crate::prde_const::admitted(inner)
                    && crate::to_expr(e)?.as_number().is_some_and(|n| n.is_rational())
                {
                    return Ok(PyPoly::_from_normalized(e, &m.gens, Some(&m.dom), None, None)?.into_pyobject(self.py)?.into_any());
                }
            }
        }
        let kw = PyDict::new(self.py);
        kw.set_item("domain", &m.dom)?;
        let mut args = vec![e.clone()];
        args.extend(m.gens.iter());
        self.cls_call(args, Some(&kw))
    }
    fn slice(&self, m: &PM<'py>, lo: usize, hi: usize) -> R<Vec<O<'py>>> {
        (lo..hi.min(m.e.len())).map(|i| self.getitem(m, i)).collect()
    }
    fn unified(&self, a: &PM<'py>, b: &PM<'py>) -> R<O<'py>> {
        if a.dom.eq(&b.dom)? {
            return Ok(a.dom.clone());
        }
        match a.dom.call_method1("unify", (&b.dom,)) {
            Ok(d) => self.norm_dom(d),
            Err(e) => {
                let s: Stop = e.into();
                match s {
                    Stop::Decline => Ok(a.dom.clone()),
                    other => Err(other),
                }
            }
        }
    }
    fn mul(&self, a: &PM<'py>, b: &PM<'py>) -> R<PM<'py>> {
        if a.c != b.r {
            return decline();
        }
        let mut items = Vec::with_capacity(a.r * b.c);
        for i in 0..a.r {
            for j in 0..b.c {
                let mut acc = self.int(0)?;
                for k in 0..a.c {
                    self.tick()?;
                    acc = acc.add(a.e[i * a.c + k].mul(&b.e[k * b.c + j])?)?;
                }
                items.push(self.ecancel(&acc)?);
            }
        }
        let dom = self.unified(a, b)?;
        Ok(self.pm(a.r, b.c, items, a.gens.clone(), dom))
    }
    /// `M / poly`.
    fn pm_div(&self, m: &PM<'py>, o: &O<'py>) -> R<PM<'py>> {
        let o = if o.is_instance_of::<PyPoly>() { self.ex(o)? } else { o.clone() };
        if o.getattr("is_zero")?.is_truthy()? {
            return decline();
        }
        let e = m.e.iter().map(|x| self.ecancel(&x.div(&o)?)).collect::<R<Vec<_>>>()?;
        Ok(self.pm(m.r, m.c, e, m.gens.clone(), m.dom.clone()))
    }
    fn neg(&self, m: &PM<'py>) -> R<PM<'py>> {
        let e = m.e.iter().map(|x| Ok(x.neg()?)).collect::<R<Vec<_>>>()?;
        Ok(self.pm(m.r, m.c, e, m.gens.clone(), m.dom.clone()))
    }
    fn row_join(&self, a: &PM<'py>, b: &PM<'py>) -> R<PM<'py>> {
        if a.r != b.r {
            return decline();
        }
        let mut e = Vec::with_capacity(a.r * (a.c + b.c));
        for i in 0..a.r {
            e.extend_from_slice(&a.e[i * a.c..(i + 1) * a.c]);
            e.extend_from_slice(&b.e[i * b.c..(i + 1) * b.c]);
        }
        let dom = self.unified(a, b)?;
        Ok(self.pm(a.r, a.c + b.c, e, a.gens.clone(), dom))
    }
    fn col_join(&self, a: &PM<'py>, b: &PM<'py>) -> R<PM<'py>> {
        if a.c != b.c {
            return decline();
        }
        let mut e = a.e.clone();
        e.extend_from_slice(&b.e);
        let dom = self.unified(a, b)?;
        Ok(self.pm(a.r + b.r, a.c, e, a.gens.clone(), dom))
    }
    /// `M.set_gens(DE.t)`: `from_Matrix(to_Matrix(), t)`.
    fn set_gens(&self, m: &PM<'py>) -> R<PM<'py>> {
        if m.r * m.c == 0 {
            return Ok(self.pm(m.r, m.c, vec![], self.gt()?, self.qq.clone()));
        }
        let e = m.e.iter().map(|x| self.ecancel(x)).collect::<R<Vec<_>>>()?;
        self.from_exprs(m.r, m.c, e)
    }
    /// `M.nullspace()` for ground entries over QQ or a fraction field (the
    /// native `_polymatrix_ground` route); anything else is the stock body.
    fn nullspace(&self, m: &PM<'py>) -> R<Vec<PM<'py>>> {
        self.tick()?;
        let Ok(dom) = m.dom.downcast::<PyDomain>() else { return decline() };
        let inner = dom.get().inner.clone();
        let mut gens = Vec::with_capacity(m.gens.len());
        for g in m.gens.iter() {
            if !g.is_instance_of::<crate::PyExprNode>() {
                return decline();
            }
            gens.push(crate::to_expr(&g)?);
        }
        if inner.flat_parameters().iter().any(|p| gens.contains(p)) {
            return decline();
        }
        let mut entries = Vec::with_capacity(m.e.len());
        for x in &m.e {
            if !x.is_instance_of::<crate::PyExprNode>() {
                return decline();
            }
            let e = crate::to_expr(x)?;
            if symox::traverse::free_symbols(&e).iter().any(|s| gens.contains(s)) {
                return decline();
            }
            entries.push(e);
        }
        let Some((vectors, _)) = crate::poly::ddm::ground_entries(entries, m.r, m.c, &inner, true, self.py)? else { return decline() };
        let vectors = vectors.into_bound(self.py).downcast_into::<PyList>()?;
        let mut out = Vec::with_capacity(vectors.len());
        for v in vectors.iter() {
            let v = v.downcast_into::<PyList>()?;
            let e: Vec<O<'py>> = v.iter().collect();
            if e.len() != m.c {
                return decline();
            }
            out.push(self.pm(m.c, 1, e, m.gens.clone(), m.dom.clone()));
        }
        Ok(out)
    }
    /// `constant_system(A, u, DE)`: the native QQ or fraction-field kernel.
    fn constant_system(&self, a: &PM<'py>, u: &PM<'py>) -> R<(PM<'py>, PM<'py>)> {
        self.tick()?;
        if a.r * a.c == 0 {
            return Ok((a.clone(), u.clone()));
        }
        if a.gens.len() != 1 {
            return decline();
        }
        let gen = a.gens.get_item(0)?;
        let py = self.py;
        let list = |m: &PM<'py>| PyList::new(py, m.e.iter());
        if a.dom.is(&self.qq) && u.dom.is(&self.qq) {
            let args = (a.r, a.c, list(a)?, u.r, u.c, list(u)?).into_pyobject(py)?;
            let Some(r) = crate::prde::constant_system_any(py, &args, &gen, self.tl.as_any())? else { return decline() };
            let r = r.into_bound(py);
            let spec = |s: O<'py>| -> R<PM<'py>> {
                let (rr, cc, ii) = (s.get_item(0)?.extract::<usize>()?, s.get_item(1)?.extract::<usize>()?, s.get_item(2)?);
                let e: Vec<O<'py>> = ii.downcast_into::<PyList>()?.iter().collect();
                Ok(self.pm(rr, cc, e, a.gens.clone(), self.qq.clone()))
            };
            return Ok((spec(r.get_item(0)?)?, spec(r.get_item(1)?)?));
        }
        if !a.dom.is_exact_instance_of::<PyDomain>() || !u.dom.is_exact_instance_of::<PyDomain>() {
            return decline();
        }
        let args = (a.r, a.c, list(a)?, &a.dom, u.r, u.c, list(u)?, &u.dom).into_pyobject(py)?;
        let Some(r) = crate::prde_const::constant_system_field(py, &args, &gen, self.tl.as_any(), self.dl.as_any(), &self.level()?)? else { return decline() };
        let r = r.into_bound(py);
        let spec = |s: O<'py>| -> R<PM<'py>> {
            let (rr, cc, ii) = (s.get_item(0)?.extract::<usize>()?, s.get_item(1)?.extract::<usize>()?, s.get_item(2)?);
            let e: Vec<O<'py>> = ii.downcast_into::<PyList>()?.iter().collect();
            let dom = self.norm_dom(s.get_item(3)?)?;
            Ok(self.pm(rr, cc, e, a.gens.clone(), dom))
        };
        Ok((spec(r.get_item(0)?)?, spec(r.get_item(1)?)?))
    }

    // ---- lower-level kernels ---------------------------------------------

    fn tower_trace(&self, trace: &O<'py>) -> R<()> {
        for s in trace.downcast::<PyList>()?.iter() {
            self.push(&s.extract::<String>()?);
        }
        Ok(())
    }
    /// The tower kernels of `prde_tower.rs` at the current level (>= 2).
    fn tower_kernel(&self, kind: &str, args: Bound<'py, PyTuple>) -> R<O<'py>> {
        let r = crate::prde_tower::kernel(self.py, kind, &args, &self.t()?, &self.d()?, self.tl.as_any(), self.dl.as_any(), &self.level()?)?;
        let Some(r) = r else { return decline() };
        let r = r.into_bound(self.py);
        let (status, payload, trace) = (r.get_item(0)?.extract::<u8>()?, r.get_item(1)?, r.get_item(2)?);
        self.tower_trace(&trace)?;
        match status {
            0 => Ok(payload),
            1 | 2 => Err(Stop::Raise(status, payload.extract::<String>()?)),
            _ => decline(),
        }
    }
    /// The base-case kernels of `prde.rs` (current level 1).
    fn base_kernel(&self, kind: &str, args: Bound<'py, PyTuple>) -> R<O<'py>> {
        let level = (-(self.nt as isize)).into_pyobject(self.py)?.into_any();
        let r = crate::prde::_prde_kernel(kind, &args, &self.t()?, &self.d()?, self.tl.as_any(), self.dl.as_any(), &level, self.py)?;
        match r {
            Some(r) => Ok(r.into_bound(self.py)),
            None => {
                note(&format!("base kernel {kind}"));
                Err(Stop::Decline)
            }
        }
    }
    /// `parametric_log_deriv(fa, fd, wa, wd, DE)` at the current level:
    /// `None` or the first two answer entries.
    fn heu(&self, fa: &O<'py>, fd: &O<'py>, wa: &O<'py>, wd: &O<'py>) -> R<Option<(O<'py>, O<'py>)>> {
        self.tick()?;
        let answer = if self.act.get() == 1 {
            let r = self.base_kernel("parametric_log_deriv_heu", PyTuple::new(self.py, [fa, fd, wa, wd])?)?;
            self.tower_trace(&r.get_item(1)?)?;
            r.get_item(0)?
        } else {
            let args = (fa, fd, wa, wd, &self.cases_obj).into_pyobject(self.py)?;
            self.tower_kernel("heu_tower", args)?
        };
        if answer.is_none() {
            return Ok(None);
        }
        Ok(Some((answer.get_item(0)?, answer.get_item(1)?)))
    }
    /// `is_log_deriv_k_t_radical_in_field(fa, fd, DE)` at the current level:
    /// `None` or the answer `(n, u, ...)`.
    fn log_deriv(&self, fa: &O<'py>, fd: &O<'py>) -> R<Option<O<'py>>> {
        self.tick()?;
        let answer = if self.act.get() == 1 {
            let r = self.base_kernel("is_log_deriv_in_field", PyTuple::new(self.py, [fa, fd])?)?;
            self.tower_trace(&r.get_item(1)?)?;
            r.get_item(0)?
        } else {
            let auto = PyString::new(self.py, "auto");
            let args = (fa, fd, auto, &self.cases_obj).into_pyobject(self.py)?;
            self.tower_kernel("is_log_deriv_tower", args)?
        };
        if answer.is_none() {
            return Ok(None);
        }
        Ok(Some(answer))
    }
    /// `derivation(z, DE)` of an expression `z` at the current level: the
    /// stock loop over `DE.D`/`DE.T` (`z.as_poly(v)` or `z`, `diff(v)`,
    /// `(d*dz).as_poly(t)`); a term that is not a polynomial in `t` makes
    /// the stock body add `None` (its exception), left to it.
    fn deriv_expr(&self, z: &O<'py>) -> R<O<'py>> {
        self.tick()?;
        let t = self.t()?;
        let mut r = self.pti(0)?;
        for k in 0..self.act.get() {
            let (d, v) = (self.dl.get_item(k)?, self.tl.get_item(k)?);
            let pv = z.call_method1("as_poly", (&v,))?;
            let pv = if pv.is_none() { z.call_method0("as_expr")? } else { pv };
            let dv = pv.call_method1("diff", (&v,))?.call_method0("as_expr")?;
            let term = d.call_method0("as_expr")?.mul(dv)?.call_method1("as_poly", (&t,))?;
            if term.is_none() {
                return decline();
            }
            r = r.add(term)?;
        }
        Ok(r)
    }
    /// `param_rischDE(fa, fd, G, DE)` at the current level.
    fn param_risch_at(&self, fa: &O<'py>, fd: &O<'py>, g: &[(O<'py>, O<'py>)]) -> R<(Vec<(O<'py>, O<'py>)>, PM<'py>)> {
        if self.act.get() >= 2 {
            return self.param_risch(fa, fd, g);
        }
        if g.is_empty() {
            return decline();
        }
        let pairs = g.iter().map(|(a, b)| PyTuple::new(self.py, [a, b])).collect::<PyResult<Vec<_>>>()?;
        let args = (fa, fd, PyList::new(self.py, pairs)?).into_pyobject(self.py)?;
        let r = match self.base_kernel("param_rischDE", args) {
            Ok(r) => r,
            // The base kernel covers the constant-`d1` case: the other base
            // inputs run the generic transliteration below at level 1.
            Err(Stop::Decline) => return self.param_risch(fa, fd, g),
            Err(e) => return Err(e),
        };
        // Dummy('z') of the base-case weak_normalizer.
        self.push("z");
        let mut h = Vec::new();
        for p in r.get_item(0)?.downcast_into::<PyList>()?.iter() {
            h.push((p.get_item(0)?, p.get_item(1)?));
        }
        let spec = r.get_item(1)?;
        let (rr, cc) = (spec.get_item(0)?.extract::<usize>()?, spec.get_item(1)?.extract::<usize>()?);
        let e: Vec<O<'py>> = spec.get_item(2)?.downcast_into::<PyList>()?.iter().collect();
        if e.len() != rr * cc {
            return decline();
        }
        Ok((h, self.pm(rr, cc, e, self.gt()?, self.qq.clone())))
    }

    // ---- risch.py / rde.py substeps --------------------------------------

    fn splitfactor(&self, p: &O<'py>) -> R<(O<'py>, O<'py>)> {
        self.tick()?;
        let kw = PyDict::new(self.py);
        kw.set_item("domain", p.call_method0("get_domain")?)?;
        let one = self.cls_call(vec![self.int(1)?, self.t()?], Some(&kw))?;
        let dp = self.deriv(p)?;
        if self.is_zero(p)? {
            return Ok((p.clone(), one));
        }
        if self.ground(p)? {
            let mut kinv = Vec::new();
            for k in 0..self.act.get() - 1 {
                kinv.push(self.pyint(1)?.div(self.tl.get_item(k)?)?);
            }
            let kinv = PyTuple::new(self.py, kinv)?;
            let a = p.call_method1("as_poly", &kinv)?;
            let b = dp.call_method1("as_poly", &kinv)?;
            if a.is_none() || b.is_none() {
                return decline();
            }
            let s = self.gcd(&a, &b)?.call_method1("as_poly", (self.t()?,))?;
            if s.is_none() {
                return decline();
            }
            return Ok((self.exquo(p, &s)?, s));
        }
        if !self.is_zero(&dp)? {
            let h = self.to_field(&self.gcd(p, &dp)?)?;
            let g = self.to_field(&self.gcd(p, &self.diff(p)?)?)?;
            let s = self.exquo(&h, &g)?;
            if self.deg(&s)? == Some(0) {
                return Ok((p.clone(), one));
            }
            let (q0, q1) = self.splitfactor(&self.exquo(p, &s)?)?;
            return Ok((q0, q1.mul(&s)?));
        }
        Ok((p.clone(), one))
    }
    fn gcdex_diophantine(&self, a: &O<'py>, b: &O<'py>, c: &O<'py>) -> R<(O<'py>, O<'py>)> {
        let sg = self.m1(a, "half_gcdex", b)?;
        let (s, g) = (sg.get_item(0)?, sg.get_item(1)?);
        let mut s = s.mul(self.exquo(c, &g)?)?;
        if !self.is_zero(&s)? && self.deg(&s)? >= self.deg(b)? {
            s = self.div(&s, b)?.1;
        }
        let t = self.exquo(&c.sub(s.mul(a)?)?, b)?;
        Ok((s, t))
    }
    fn weak_normalizer(&self, a: &O<'py>, d: &O<'py>) -> R<(O<'py>, (O<'py>, O<'py>))> {
        self.push("z");
        let (dn, _) = self.splitfactor(d)?;
        let g = self.gcd(&dn, &self.diff(&dn)?)?;
        let dsq = self.quo(&dn, &g)?;
        let d1 = self.quo(&dsq, &self.gcd(&dsq, &g)?)?;
        self.gcdex_diophantine(&self.quo(d, &d1)?, &d1, a)?;
        let dd1 = self.deriv(&d1)?;
        let zt = self.pt(&self.z)?;
        let r = self.m1(&a.sub(zt.mul(&dd1)?)?, "resultant", &d1)?;
        let r = self.cls_call(vec![r, self.z.clone()], None)?;
        if self.ground(&r)? {
            return Ok((self.pti(1)?, (a.clone(), d.clone())));
        }
        // Dummy() of Poly.real_roots (CRootOf).
        self.push("");
        let inner = r.getattr("domain")?.downcast_into::<PyDomain>()?.get().inner.clone();
        let r = if matches!(inner, symox::poly::domains::ExactDomain::ZZ | symox::poly::domains::ExactDomain::QQ) {
            r
        } else {
            // CRootOf._preprocess_roots: `preprocess_roots` starts with
            // clear_denoms(convert=True), primitive() and retract(); when that
            // already gives ZZ the remaining steps keep the roots and CRootOf
            // proceeds. Any other outcome (symbolic stripping, the
            // "sorted roots not supported" NotImplementedError) is the stock
            // body's.
            let kw = PyDict::new(self.py);
            kw.set_item("convert", true)?;
            let p = r.call_method("clear_denoms", (), Some(&kw))?.get_item(1)?;
            let p = p.call_method0("primitive")?.get_item(1)?;
            let p = self.cls_call(vec![self.ex(&p)?, self.z.clone()], None)?;
            let pd = p.getattr("domain")?.downcast_into::<PyDomain>()?.get().inner.clone();
            if !matches!(pd, symox::poly::domains::ExactDomain::ZZ) {
                return decline();
            }
            p
        };
        let fl = r.call_method0("factor_list")?;
        let mut roots: Vec<i64> = Vec::new();
        for fk in fl.get_item(1)?.try_iter()? {
            let fk = fk?;
            let (f, k) = (fk.get_item(0)?, fk.get_item(1)?.extract::<usize>()?);
            if self.deg(&f)? != Some(1) {
                continue;
            }
            let root = self.nth(&f, 0)?.neg()?.div(self.nth(&f, 1)?)?;
            if root.getattr("is_Integer")?.is_truthy()? {
                let n = self.integer(&root)?;
                if n > 0 {
                    roots.extend(std::iter::repeat_n(n, k));
                }
            }
        }
        roots.sort_unstable();
        let mut q = self.pti(1)?;
        for n in roots {
            let f = self.gcd(&a.sub(self.pti(n)?.mul(&dd1)?)?, &d1)?;
            q = q.mul(f)?;
        }
        let dq = self.deriv(&q)?;
        let sn = q.mul(a)?.sub(d.mul(&dq)?)?;
        let sd = q.mul(d)?;
        let (sn, sd) = self.cancel(&sn, &sd)?;
        Ok((q, (sn, sd)))
    }
    #[allow(clippy::type_complexity)]
    fn normal_denom(&self, fa: &O<'py>, fd: &O<'py>, g: &[(O<'py>, O<'py>)]) -> R<(O<'py>, (O<'py>, O<'py>), Vec<(O<'py>, O<'py>)>, O<'py>)> {
        let (dn, _) = self.splitfactor(fd)?;
        let mut gd = self.pti(1)?;
        for (_, d) in g {
            gd = self.lcm(&gd, d)?;
        }
        let (en, _) = self.splitfactor(&gd)?;
        let p = self.gcd(&dn, &en)?;
        let h = self.quo(&self.gcd(&en, &self.diff(&en)?)?, &self.gcd(&p, &self.diff(&p)?)?)?;
        let a = dn.mul(&h)?;
        let c = a.mul(&h)?;
        let ba = a.mul(fa)?.sub(dn.mul(self.deriv(&h)?)?.mul(fd)?)?;
        let (ba, bd) = self.cancel(&ba, fd)?;
        let g = g.iter().map(|(x, y)| self.cancel(&c.mul(x)?, y)).collect::<R<Vec<_>>>()?;
        Ok((a, (ba, bd), g, h))
    }
    /// `order_at(f, Poly(t, t), t)`, `None` for `oo`.
    fn order_t(&self, f: &O<'py>) -> R<Option<i64>> {
        let Some(n) = self.deg(f)? else { return Ok(None) };
        for i in 0..=n {
            self.tick()?;
            if !self.nth(f, i)?.getattr("is_zero")?.is_truthy()? {
                return Ok(Some(i));
            }
        }
        decline()
    }
    #[allow(clippy::type_complexity)]
    fn special_denom(&self, a: &O<'py>, ba: &O<'py>, bd: &O<'py>, g: Vec<(O<'py>, O<'py>)>) -> R<(O<'py>, O<'py>, Vec<(O<'py>, O<'py>)>, O<'py>)> {
        let case = self.case();
        if case == "primitive" || case == "base" {
            return Ok((a.clone(), self.quo(ba, bd)?, g, self.pti(1)?));
        }
        if case != "exp" {
            return decline();
        }
        let t = self.t()?;
        let p = self.pt(&t)?;
        if self.is_zero(bd)? || g.is_empty() {
            return decline();
        }
        let nb = match self.order_t(ba)? {
            None => None,
            Some(x) => Some(x - self.order_t(bd)?.ok_or(Stop::Decline)?),
        };
        let mut nc: Option<i64> = None;
        for (ga, gd) in &g {
            let term = match self.order_t(ga)? {
                None => None,
                Some(x) => match self.order_t(gd)? {
                    None => return decline(),
                    Some(y) => Some(x - y),
                },
            };
            if let Some(v) = term {
                nc = Some(nc.map_or(v, |c| c.min(v)));
            }
        }
        let mut n = match nc {
            None => 0,
            Some(c) => 0.min(c - 0.min(nb.unwrap_or(0))),
        };
        if nb == Some(0) {
            let dcoeff = self.quo(&self.d()?, &p)?;
            let e = ba.call_method1("eval", (0,))?.neg()?.div(bd.call_method1("eval", (0,))?)?.div(a.call_method1("eval", (0,))?)?;
            let m = self.lower(|cx| {
                let tl = cx.t()?;
                let (aa, ad) = cx.frac_in(&e, &tl, false)?;
                let (ea, ed) = cx.frac_in(&dcoeff, &tl, false)?;
                match cx.heu(&aa, &ad, &ea, &ed)? {
                    Some((q, m)) if q.eq(1i64)? => Ok(Some(cx.integer(&m)?)),
                    _ => Ok(None),
                }
            })?;
            if let Some(m) = m {
                n = n.min(m);
            }
        }
        let big_n = nb.map_or(0, |b| 0.max(-b));
        let pn_big = p.pow(big_n, self.py.None())?;
        let pn = p.pow(-n, self.py.None())?;
        let aa = a.mul(&pn_big)?;
        let bb = ba.mul(self.quo(&pn_big, bd)?)?
            .add(self.pti(n)?.mul(a)?.mul(self.quo(&self.deriv(&p)?, &p)?)?.mul(&pn_big)?)?;
        let g = g.iter().map(|(x, y)| self.cancel(&x.mul(&pn_big)?.mul(&pn)?, y)).collect::<R<Vec<_>>>()?;
        Ok((aa, bb, g, pn))
    }
    fn linear_constraints(&self, g: &[(O<'py>, O<'py>)]) -> R<(Vec<O<'py>>, PM<'py>)> {
        let m = g.len();
        let Some((_, d0)) = g.first() else { return decline() };
        let mut d = d0.clone();
        for (_, y) in &g[1..] {
            d = self.lcm(&d, y)?;
        }
        let d = self.cls_call(vec![d], Some(&self.fld))?;
        let q = g.iter().map(|(x, y)| self.div(&x.mul(self.quo(&d, y)?)?, &d)).collect::<R<Vec<_>>>()?;
        let qs: Vec<O<'py>> = q.iter().map(|x| x.0.clone()).collect();
        let rs: Vec<O<'py>> = q.into_iter().map(|x| x.1).collect();
        Ok((qs, self.remainder_matrix(&rs, m)?))
    }
    /// `Matrix(N + 1, m, lambda i, j: r[j].nth(i), t)` with `N` the largest
    /// remainder degree, or `Matrix(0, m, [], t)` when every one is zero.
    fn remainder_matrix(&self, rs: &[O<'py>], m: usize) -> R<PM<'py>> {
        let mut n: Option<i64> = None;
        for r in rs {
            n = n.max(self.deg(r)?);
        }
        let Some(n) = n else { return Ok(self.pm(0, m, vec![], self.gt()?, self.qq.clone())) };
        let rows = (n + 1) as usize;
        let mut items = Vec::with_capacity(rows * m);
        for i in 0..rows {
            for r in rs {
                items.push(self.nth(r, i as i64)?);
            }
        }
        self.from_exprs(rows, m, items)
    }
    fn poly_linear_constraints(&self, p: &[O<'py>], d: &O<'py>) -> R<(Vec<O<'py>>, PM<'py>)> {
        let gens = d.getattr("gens")?.downcast_into::<PyTuple>()?;
        if !gens.as_any().eq(self.gt()?)? {
            return decline();
        }
        let q = p.iter().map(|x| self.div(x, d)).collect::<R<Vec<_>>>()?;
        let qs: Vec<O<'py>> = q.iter().map(|x| x.0.clone()).collect();
        let rs: Vec<O<'py>> = q.into_iter().map(|x| x.1).collect();
        Ok((qs, self.remainder_matrix(&rs, p.len())?))
    }
    /// `bound_degree(a, b, cQ, DE, parametric=True)`.
    fn bound_degree(&self, a: &O<'py>, b: &O<'py>, cq: &[O<'py>]) -> R<i64> {
        let (da, db) = (self.deg(a)?, self.deg(b)?);
        let mut dc: Option<i64> = None;
        let mut first = true;
        for q in cq {
            let x = self.deg(q)?;
            dc = if first { x } else { dc.max(x) };
            first = false;
        }
        if first {
            return decline();
        }
        let alpha = self.ecancel(&self.ex(&self.lc(b)?)?.neg()?.div(self.ex(&self.lc(a)?)?)?)?;
        match self.case().as_str() {
            "base" => {
                // n = max(0, dc - max(db, da - 1)); a != 0.
                let Some(da) = da else { return decline() };
                let mx = db.unwrap_or(i64::MIN).max(da - 1);
                let mut n = dc.map_or(0, |c| 0.max(c - mx));
                if db == Some(da - 1) && alpha.getattr("is_Integer")?.is_truthy()? {
                    let al = self.integer(&alpha)?;
                    n = 0.max(al).max(dc.map_or(i64::MIN, |c| c - (da - 1)));
                }
                if n > 1 << 20 {
                    return decline();
                }
                Ok(n)
            }
            "exp" => {
                let Some(mx) = da.max(db) else { return decline() };
                let mut n = dc.map_or(0, |c| 0.max(c - mx));
                if da == db {
                    let tq = self.quo(&self.d()?, &self.pt(&self.t()?)?)?;
                    let below = self.tl.get_item(self.act.get() - 2)?;
                    let (ea, ed) = self.frac_in(&tq, &below, false)?;
                    let m = self.lower(|cx| {
                        let (aa, ad) = cx.frac_in(&alpha, &cx.t()?, false)?;
                        match cx.heu(&aa, &ad, &ea, &ed)? {
                            Some((q, m)) if q.eq(1i64)? => Ok(Some(cx.integer(&m)?)),
                            _ => Ok(None),
                        }
                    })?;
                    if let Some(m) = m {
                        n = n.max(m);
                    }
                }
                Ok(n)
            }
            "primitive" => {
                let n = if db > da {
                    dc.map_or(0, |c| 0.max(c - db.unwrap()))
                } else {
                    match (dc, da) {
                        (None, Some(_)) => 0,
                        (Some(c), Some(x)) => 0.max(c - x + 1),
                        _ => return decline(),
                    }
                };
                let below = self.tl.get_item(self.act.get() - 2)?;
                let eta = self.frac_in(&self.d()?, &below, false)?;
                let t1 = self.t()?;
                let one_less = match (db, da) {
                    (Some(x), Some(y)) => x == y - 1,
                    (None, None) => true,
                    _ => false,
                };
                self.lower(|cx| {
                    let (aa, ad) = cx.frac_in(&alpha, &cx.t()?, false)?;
                    if one_less {
                        // The stock body unpacks `limited_integrate(...)` one
                        // level down and takes `max(n, m[0])` of a Poly: a
                        // TypeError either way (None or Poly vs int), unless
                        // NonElementaryIntegralException. Left to it.
                        return decline();
                    }
                    if db == da {
                        if let Some(ans) = cx.log_deriv(&aa, &ad)? {
                            let z = ans.get_item(1)?;
                            if ans.get_item(0)?.eq(1i64)? {
                                // beta = -(a*Dz + b*z).LC()/(z*a.LC())
                                let dz = cx.deriv_expr(&z)?.call_method1("as_poly", (&t1,))?;
                                let z1 = z.call_method1("as_poly", (&t1,))?;
                                if dz.is_none() || z1.is_none() {
                                    return decline();
                                }
                                let lcs = a.mul(dz)?.add(b.mul(z1)?)?.call_method0("LC")?;
                                let beta = lcs.neg()?.div(z.call_method0("as_expr")?.mul(a.call_method0("LC")?)?)?;
                                let (ba, bd) = cx.frac_in(&beta, &cx.t()?, false)?;
                                // `None` makes the stock body raise its
                                // TypeError (unpacking), left to it.
                                let Some((_, m)) = cx.limited_integrate(&ba, &bd, std::slice::from_ref(&eta))? else { return decline() };
                                if m.len() != 1 {
                                    return Err(Stop::Raise(2, "Length of m should be 1".into()));
                                }
                                let m0 = cx.ex(&m[0])?;
                                if m0.getattr("is_Integer")?.is_truthy()? {
                                    let m0 = cx.integer(&m0)?;
                                    if m0 > 1 << 20 {
                                        return decline();
                                    }
                                    return Ok(n.max(m0));
                                }
                                // A larger non-integer would leave the stock
                                // body with a rational degree bound.
                                if m0.gt(n)? {
                                    return decline();
                                }
                            }
                        }
                    }
                    Ok(n)
                })
            }
            _ => decline(),
        }
    }
    #[allow(clippy::type_complexity)]
    fn spde(&self, a: &O<'py>, b: &O<'py>, q: &[O<'py>], n: i64) -> R<(O<'py>, O<'py>, Vec<O<'py>>, Vec<O<'py>>, i64)> {
        if self.is_zero(a)? || q.is_empty() {
            return decline();
        }
        let rz = q.iter().map(|qi| self.gcdex_diophantine(b, a, qi)).collect::<R<Vec<_>>>()?;
        let bb = b.add(self.deriv(a)?)?;
        let qq = rz.iter().map(|(r, z)| Ok(z.sub(self.deriv(r)?)?)).collect::<R<Vec<_>>>()?;
        let r: Vec<O<'py>> = rz.into_iter().map(|x| x.0).collect();
        let Some(da) = self.deg(a)? else { return decline() };
        Ok((a.clone(), bb, qq, r, n - da))
    }
    fn b_large(&self, b: &O<'py>, mut q: Vec<O<'py>>, n: i64) -> R<(Vec<O<'py>>, PM<'py>)> {
        let Some(db) = self.deg(b)? else { return decline() };
        let m = q.len();
        let mut h = vec![self.pti(0)?; m];
        let t = self.t()?;
        let lcb = self.lc(b)?;
        let mut nn = n;
        while nn >= 0 {
            for i in 0..m {
                self.tick()?;
                let si = self.nth(&q[i], nn + db)?.div(&lcb)?;
                let sitn = self.pt(&si.mul(t.pow(nn, self.py.None())?)?)?;
                h[i] = h[i].add(&sitn)?;
                q[i] = q[i].sub(self.deriv(&sitn)?)?.sub(b.mul(&sitn)?)?;
            }
            nn -= 1;
        }
        let mut dc: i64 = -1;
        for qi in &q {
            if let Some(x) = self.deg(qi)? {
                dc = dc.max(x);
            }
        }
        let rows = (dc + 1) as usize;
        let mut items = Vec::with_capacity(rows * m);
        for i in 0..rows {
            for qj in &q {
                items.push(self.nth(qj, i as i64)?);
            }
        }
        let mm = self.from_exprs(rows, m, items)?;
        let (a, _) = self.constant_system(&mm, &self.zeros(rows, 1)?)?;
        let c = self.eye(m)?;
        let a = self.col_join(&self.row_join(&a, &self.zeros(a.r, m)?)?, &self.row_join(&c, &self.neg(&c)?)?)?;
        Ok((h, a))
    }
    /// `prde_no_cancel_b_small(b, Q, n, DE)` in the base case (`b == 0`,
    /// `D == d/dt`).
    fn b_small_base(&self, b: &O<'py>, mut q: Vec<O<'py>>, n: i64) -> R<(Vec<O<'py>>, PM<'py>)> {
        if !self.is_zero(b)? || self.case() != "base" {
            return decline();
        }
        let m = q.len();
        let mut h = vec![self.pti(0)?; m];
        let t = self.t()?;
        let d = self.d()?;
        let Some(dd) = self.deg(&d)? else { return decline() };
        let lcd = self.lc(&d)?;
        for nn in (1..=n).rev() {
            for i in 0..m {
                self.tick()?;
                let si = self.nth(&q[i], nn + dd - 1)?.div(self.pyint(nn)?.mul(&lcd)?)?;
                let sitn = self.pt(&si.mul(t.pow(nn, self.py.None())?)?)?;
                h[i] = h[i].add(&sitn)?;
                q[i] = q[i].sub(self.deriv(&sitn)?)?.sub(b.mul(&sitn)?)?;
            }
        }
        // b == 0: f = [Poly(1, t, field=True)], B = [[qi.TC() ...] + [0]].
        let f = self.cls_call(vec![self.int(1)?, t.clone()], Some(&self.fld))?;
        let mut tcs = Vec::with_capacity(m + 1);
        for qi in &q {
            tcs.push(qi.call_method0("TC")?);
        }
        tcs.push(self.int(0)?);
        let bm = self.from_exprs(1, m + 1, tcs)?;
        let mut dmax: Option<i64> = None;
        for qi in &q {
            dmax = dmax.max(self.deg(qi)?);
        }
        let a = match dmax {
            Some(dg) if dg > 0 => {
                let rows = dg as usize;
                let mut items = Vec::with_capacity(rows * m);
                for i in 0..rows {
                    for qj in &q {
                        items.push(self.nth(qj, i as i64 + 1)?);
                    }
                }
                let mm = self.from_exprs(rows, m, items)?;
                self.constant_system(&mm, &self.zeros(rows, 1)?)?.0
            }
            _ => self.pm(0, m, vec![], self.gt()?, self.qq.clone()),
        };
        let r = 1;
        let id = self.eye(m)?;
        let a = self.row_join(&a, &self.zeros(a.r, r + m)?)?;
        let bm = self.row_join(&bm, &self.zeros(bm.r, m)?)?;
        let c = self.row_join(&self.row_join(&id, &self.zeros(m, r)?)?, &self.neg(&id)?)?;
        let mut out = vec![f];
        out.extend(h);
        Ok((out, self.col_join(&self.col_join(&a, &bm)?, &c)?))
    }
    fn cancel_liouvillian(&self, b: &O<'py>, q: Vec<O<'py>>, n: i64) -> R<(Vec<O<'py>>, PM<'py>)> {
        if n < 0 {
            return decline();
        }
        let case = self.case();
        if case != "exp" && case != "primitive" {
            return decline();
        }
        let mut q = q;
        let mut h: Vec<O<'py>> = Vec::new();
        let mut bab: Option<(O<'py>, O<'py>)> = None;
        if case == "primitive" {
            bab = Some(self.lower(|cx| cx.frac_in(b, &cx.t()?, true))?);
        }
        let gens = b.getattr("gens")?.downcast_into::<PyTuple>()?;
        let t = self.t()?;
        let mut mat: Option<PM<'py>> = None;
        for i in (0..=n).rev() {
            self.tick()?;
            if case == "exp" {
                bab = Some(self.lower(|cx| {
                    let tl = cx.t()?;
                    let dt = cx.deriv(&cx.pt(&tl)?)?;
                    let e = cx.pyint(i)?.mul(dt.div(&tl)?)?;
                    let p = e.call_method1("as_poly", (&gens,))?;
                    if p.is_none() {
                        return decline();
                    }
                    cx.frac_in(&b.add(p)?, &tl, true)
                })?);
            }
            let (ba, bd) = bab.clone().ok_or(Stop::Decline)?;
            let (fi, ai) = self.lower(|cx| {
                let tl = cx.t()?;
                let qy = q.iter().map(|qq| cx.frac_in(&cx.nth(qq, i)?, &tl, true)).collect::<R<Vec<_>>>()?;
                cx.param_risch_at(&ba, &bd, &qy)
            })?;
            let fi = fi.iter().map(|(fa, fd)| self.cls_call(vec![self.ex(fa)?.div(self.ex(fd)?)?, t.clone()], Some(&self.fld))).collect::<R<Vec<_>>>()?;
            let ai = self.set_gens(&ai)?;
            let ri = fi.len();
            mat = Some(match mat {
                None => ai,
                Some(m) => self.col_join(&ai, &self.row_join(&m, &self.zeros(m.r, ri)?)?)?,
            });
            let ti = t.pow(i, self.py.None())?;
            let mut fs = Vec::with_capacity(ri);
            for f in &fi {
                let fg = f.getattr("gens")?;
                let tp = ti.call_method1("as_poly", (fg,))?;
                if tp.is_none() {
                    return decline();
                }
                let hji = f.mul(tp)?;
                fs.push(self.deriv(&hji)?.sub(b.mul(&hji)?)?.neg()?);
                h.push(hji);
            }
            q.extend(fs);
        }
        Ok((h, mat.ok_or(Stop::Decline)?))
    }

    // ---- prde.py entry points --------------------------------------------

    fn param_poly(&self, a: &O<'py>, b: &O<'py>, q: Vec<O<'py>>, n: i64) -> R<(Vec<O<'py>>, PM<'py>)> {
        self.tick()?;
        let m = q.len();
        if n < 0 {
            let mut all_zero = true;
            for qi in &q {
                all_zero &= self.is_zero(qi)?;
            }
            if all_zero {
                return Ok((vec![], self.zeros(1, m)?));
            }
            let mm = self.remainder_matrix(&q, m)?;
            let (a, _) = self.constant_system(&mm, &self.zeros(mm.r, 1)?)?;
            return Ok((vec![], a));
        }
        if self.ground(a)? {
            let c = self.lc(a)?;
            if c.getattr("is_zero")?.is_truthy()? {
                return decline();
            }
            let b = self.to_field(b)?.call_method1("exquo_ground", (&c,))?;
            let q = q.iter().map(|qi| Ok(self.to_field(qi)?.call_method1("exquo_ground", (&c,))?)).collect::<R<Vec<_>>>()?;
            let Some(dd) = self.deg(&self.d()?)? else { return decline() };
            let degb = self.deg(&b)?;
            let base = self.case() == "base";
            if degb.is_some() && (base || degb.is_some_and(|x| x > 0.max(dd - 1))) {
                return self.b_large(&b, q, n);
            }
            if degb.is_none_or(|x| x < dd - 1) && base {
                return self.b_small_base(&b, q, n);
            }
            if dd >= 2 {
                // prde_no_cancel_b_small / b_equal / hypertangent: stock body.
                return decline();
            }
            return match self.case().as_str() {
                "exp" | "primitive" => self.cancel_liouvillian(&b, q, n),
                _ => decline(),
            };
        }
        let mut a = a.clone();
        let mut b = b.clone();
        let mut q = q;
        let mut n = n;
        let kw = PyDict::new(self.py);
        kw.set_item("domain", a.getattr("domain")?)?;
        let mut cargs = vec![self.int(1)?];
        cargs.extend(a.getattr("gens")?.downcast_into::<PyTuple>()?.iter());
        let mut alpha = self.cls_call(cargs.clone(), Some(&kw))?;
        cargs[0] = self.int(0)?;
        let zero = self.cls_call(cargs, Some(&kw))?;
        let mut beta = vec![zero; m];
        let mut d;
        loop {
            self.tick()?;
            let (a1, b1, q1, r, n1) = self.spde(&a, &b, &q, n)?;
            a = a1;
            b = b1;
            q = q1;
            n = n1;
            for (bi, ri) in beta.iter_mut().zip(&r) {
                *bi = bi.add(alpha.mul(ri)?)?;
            }
            alpha = alpha.mul(&a)?;
            d = self.gcd(&a, &b)?;
            if !self.ground(&d)? || n < 0 {
                break;
            }
        }
        let (qq, mm) = self.poly_linear_constraints(&q, &d)?;
        let (aa, _) = self.constant_system(&mm, &self.zeros(mm.r, 1)?)?;
        let v = self.nullspace(&aa)?;
        if v.is_empty() {
            return Ok((vec![], self.eye(m)?));
        }
        let mqq = self.from_polys(1, qq.len(), qq)?;
        let r = v.iter().map(|vj| self.getitem(&self.mul(&mqq, vj)?, 0)).collect::<R<Vec<_>>>()?;
        let mbeta = self.from_polys(1, beta.len(), beta)?;
        let f = v.iter().map(|vj| self.getitem(&self.mul(&mbeta, vj)?, 0)).collect::<R<Vec<_>>>()?;
        let (g, bm) = self.param_poly(&self.quo(&a, &d)?, &self.quo(&b, &d)?, r, n)?;
        let mut h = f;
        for gk in &g {
            h.push(alpha.mul(gk)?);
        }
        let am = self.assemble(m, &v, g.len(), &bm)?;
        Ok((h, am))
    }
    /// `A = -eye(m)`, `row_join` each `vj`, `zeros(m, len)`, then
    /// `col_join(zeros(B.rows, m).row_join(B))`.
    fn assemble(&self, m: usize, v: &[PM<'py>], len: usize, bm: &PM<'py>) -> R<PM<'py>> {
        let mut a = self.neg(&self.eye(m)?)?;
        for vj in v {
            a = self.row_join(&a, vj)?;
        }
        a = self.row_join(&a, &self.zeros(m, len)?)?;
        self.col_join(&a, &self.row_join(&self.zeros(bm.r, m)?, bm)?)
    }
    fn param_risch(&self, fa: &O<'py>, fd: &O<'py>, g: &[(O<'py>, O<'py>)]) -> R<(Vec<(O<'py>, O<'py>)>, PM<'py>)> {
        self.tick()?;
        let m = g.len();
        if m == 0 {
            return decline();
        }
        let (q, (fa, fd)) = self.weak_normalizer(fa, fd)?;
        let mut gamma = q.clone();
        let g = g.iter().map(|(x, y)| self.cancel(&q.mul(x)?, y)).collect::<R<Vec<_>>>()?;
        let (a, (ba, bd), g, hn) = self.normal_denom(&fa, &fd, &g)?;
        gamma = gamma.mul(&hn)?;
        let (aa, bb, g, hs) = self.special_denom(&a, &ba, &bd, g)?;
        gamma = gamma.mul(&hs)?;
        let gg = self.gcd(&aa, &bb)?;
        let a = self.quo(&aa, &gg)?;
        let b = self.quo(&bb, &gg)?;
        let g = g.iter().map(|(x, y)| self.cancel(x, &y.mul(&gg)?)).collect::<R<Vec<_>>>()?;
        let (qs, mm) = self.linear_constraints(&g)?;
        let (mm, _) = self.constant_system(&mm, &self.zeros(mm.r, 1)?)?;
        let v = self.nullspace(&mm)?;
        if v.is_empty() {
            return Ok((vec![], self.eye(m)?));
        }
        let mq = self.from_polys(1, qs.len(), qs)?;
        let r = v.iter().map(|vj| self.getitem(&self.mul(&mq, vj)?, 0)).collect::<R<Vec<_>>>()?;
        let n = match self.bound_degree(&a, &b, &r) {
            Ok(n) => n,
            Err(Stop::Raise(1, _)) => 5,
            Err(e) => return Err(e),
        };
        let (h, bm) = self.param_poly(&a, &b, r, n)?;
        let am = self.assemble(m, &v, h.len(), &bm)?;
        let w = self.nullspace(&am)?;
        let nv = h.len();
        if nv == 0 {
            return decline();
        }
        let mut items = Vec::new();
        for wl in &w {
            let len = wl.e.len();
            items.extend(self.slice(wl, 0, m)?);
            items.extend(self.slice(wl, len.saturating_sub(nv), len)?);
        }
        let mw = self.from_polys(w.len(), m + nv, items)?;
        let nn = self.nullspace(&mw)?;
        let c = if nn.is_empty() {
            self.pm(0, 0, vec![], self.gt()?, self.qq.clone())
        } else {
            let cols = nn[0].e.len();
            let mut items = Vec::new();
            for ni in &nn {
                if ni.e.len() != cols {
                    return decline();
                }
                items.extend(self.slice(ni, 0, ni.e.len())?);
            }
            self.from_polys(nn.len(), cols, items)?
        };
        let h = h.iter().map(|hk| self.cancel(hk, &gamma)).collect::<R<Vec<_>>>()?;
        Ok((h, c))
    }
    #[allow(clippy::type_complexity)]
    fn limited_integrate(&self, fa: &O<'py>, fd: &O<'py>, g: &[(O<'py>, O<'py>)]) -> R<Option<((O<'py>, O<'py>), Vec<O<'py>>)>> {
        let fa = fa.mul(self.pt(&self.pyint(1)?.div(self.lc(fd)?)?)?)?;
        let fd = self.monic(fd)?;
        let mut gg = vec![(fa, fd)];
        gg.extend(g.iter().cloned());
        let (h, a) = self.param_risch_at(&self.pti(0)?, &self.pti(1)?, &gg)?;
        let mut v = Vec::new();
        for x in self.nullspace(&a)? {
            if self.getitem(&x, 0)?.ne(0i64)? {
                v.push(x);
            }
        }
        let Some(v0) = v.into_iter().next() else { return Ok(None) };
        let c0 = self.getitem(&v0, 0)?;
        let v = self.pm_div(&v0, &c0.neg()?)?;
        let r = h.len();
        if r == 0 {
            return decline();
        }
        let len = v.r * v.c;
        if len < r + 1 {
            return decline();
        }
        let m = len - r - 1;
        let c = self.slice(&v, 1, m + 1)?;
        let mut acc = self.pyint(0)?;
        for (i, (ha, hd)) in h.iter().enumerate() {
            acc = acc.add(self.getitem(&v, m + 1 + i)?.mul(self.ex(ha)?)?.div(self.ex(hd)?)?)?;
        }
        let y = acc.neg()?;
        let nd = y.call_method0("as_numer_denom")?;
        let ya = self.pt(&nd.get_item(0)?)?;
        let yd = self.pt(&nd.get_item(1)?)?;
        let ya = ya.mul(self.pt(&self.pyint(1)?.div(self.lc(&yd)?)?)?)?;
        Ok(Some(((ya, self.monic(&yd)?), c)))
    }
}

fn plain_symbol(e: &Expr) -> bool {
    e.is_symbol() && !e.has_host() && !e.has_wild_properties()
        && symox::assumptions::fact(e, "commutative") == Some(true)
}

fn strings(v: &Bound<'_, PyAny>) -> Option<Vec<String>> {
    let l = v.downcast_exact::<PyList>().ok()?;
    l.iter().map(|s| if s.is_exact_instance_of::<PyString>() { s.extract::<String>().ok() } else { None }).collect()
}

/// A stock Poly over an exact admitted ground domain with generators `(t,)`.
fn poly_in(p: &Bound<'_, PyAny>, t: &Expr) -> PyResult<bool> {
    if !crate::prde_const::exact_poly(p)? {
        return Ok(false);
    }
    let p = p.downcast::<PyPoly>()?.get();
    Ok(p.p.true_gens() == [t.clone()])
}

fn pairs<'py>(v: &O<'py>, t: &Expr) -> PyResult<Option<Vec<(O<'py>, O<'py>)>>> {
    let Ok(l) = v.downcast_exact::<PyList>() else { return Ok(None) };
    if l.len() > 64 {
        return Ok(None);
    }
    let mut out = Vec::with_capacity(l.len());
    for item in l.iter() {
        let Ok(tp) = item.downcast_exact::<PyTuple>() else { return Ok(None) };
        if tp.len() != 2 {
            return Ok(None);
        }
        let (a, b) = (tp.get_item(0)?, tp.get_item(1)?);
        if !poly_in(&a, t)? || !poly_in(&b, t)? {
            return Ok(None);
        }
        out.push((a, b));
    }
    Ok(Some(out))
}

fn polys<'py>(v: &O<'py>, t: &Expr) -> PyResult<Option<Vec<O<'py>>>> {
    let Ok(l) = v.downcast_exact::<PyList>() else { return Ok(None) };
    if l.len() > 64 {
        return Ok(None);
    }
    let mut out = Vec::with_capacity(l.len());
    for item in l.iter() {
        if !poly_in(&item, t)? {
            return Ok(None);
        }
        out.push(item);
    }
    Ok(Some(out))
}

fn small_int(v: &O<'_>) -> Option<i64> {
    if !v.is_exact_instance_of::<PyInt>() {
        return None;
    }
    v.extract::<i64>().ok().filter(|n| n.abs() <= 1 << 20)
}

/// `kind` with `args`:
/// - `param_rischDE_tower`: `(fa, fd, G, cases)` -> `(h, spec)`;
/// - `param_poly_rischDE_tower`: `(a, b, q, n, cases)` -> `(h, spec)`;
/// - `prde_cancel_liouvillian_tower`: `(b, Q, n, cases)` -> `(H, spec)`;
/// - `limited_integrate_tower`: `(fa, fd, G, cases)` -> `None` or
///   `((Ya, Yd), C)`;
///
/// `spec = (rows, cols, items, gens, domain)` of the PolyMatrix. Returns
/// `None` (stock body) or `(status, payload, trace)` as `prde_tower::kernel`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn kernel(py: Python<'_>, kind: &str, args: &Bound<'_, PyTuple>, t: &Bound<'_, PyAny>, d: &Bound<'_, PyAny>,
    tower: &Bound<'_, PyAny>, derivations: &Bound<'_, PyAny>, level: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    let nargs = match kind {
        "param_poly_rischDE_tower" => 5,
        _ => 4,
    };
    if args.len() != nargs {
        return Ok(None);
    }
    let (Ok(tl), Ok(dl)) = (tower.downcast_exact::<PyList>(), derivations.downcast_exact::<PyList>()) else { return Ok(None) };
    if !level.is_exact_instance_of::<PyInt>() {
        return Ok(None);
    }
    let Ok(level) = level.extract::<isize>() else { return Ok(None) };
    let nt = tl.len();
    let act = nt as isize + level + 1;
    if act < 2 || act > nt as isize || tl.len() != dl.len() || nt > 16 {
        return Ok(None);
    }
    let act = act as usize;
    if !t.is(&tl.get_item(act - 1)?) || !d.is(&dl.get_item(act - 1)?) {
        return Ok(None);
    }
    let cases_obj = args.get_item(nargs - 1)?;
    let Some(cases) = strings(&cases_obj) else { return Ok(None) };
    if cases.len() != nt || cases[0] != "base" || !cases[1..act].iter().all(|c| c == "exp" || c == "primitive") {
        return Ok(None);
    }
    let mut t_vars: Vec<Expr> = Vec::with_capacity(act);
    for v in tl.iter().take(act) {
        if !v.is_instance_of::<crate::PyExprNode>() {
            return Ok(None);
        }
        let e = crate::to_expr(&v)?;
        if !plain_symbol(&e) || t_vars.contains(&e) {
            return Ok(None);
        }
        t_vars.push(e);
    }
    for (i, dv) in dl.iter().take(act).enumerate() {
        if !poly_in(&dv, &t_vars[i])? {
            return Ok(None);
        }
        if i == 0 && dv.downcast::<PyPoly>()?.get().p.as_expr() != Expr::one() {
            return Ok(None);
        }
    }
    let top = &t_vars[act - 1];
    let arg = |i: usize| args.get_item(i);
    enum Input<'py> {
        Risch(O<'py>, O<'py>, Vec<(O<'py>, O<'py>)>),
        Poly(O<'py>, O<'py>, Vec<O<'py>>, i64),
        Cancel(O<'py>, Vec<O<'py>>, i64),
        Limited(O<'py>, O<'py>, Vec<(O<'py>, O<'py>)>),
    }
    let input = match kind {
        "param_rischDE_tower" | "limited_integrate_tower" => {
            let (fa, fd) = (arg(0)?, arg(1)?);
            if !poly_in(&fa, top)? || !poly_in(&fd, top)? {
                return Ok(None);
            }
            let Some(g) = pairs(&arg(2)?, top)? else { return Ok(None) };
            if kind == "param_rischDE_tower" { Input::Risch(fa, fd, g) } else { Input::Limited(fa, fd, g) }
        }
        "param_poly_rischDE_tower" => {
            let (a, b) = (arg(0)?, arg(1)?);
            if !poly_in(&a, top)? || !poly_in(&b, top)? {
                return Ok(None);
            }
            let (Some(q), Some(n)) = (polys(&arg(2)?, top)?, small_int(&arg(3)?)) else { return Ok(None) };
            Input::Poly(a, b, q, n)
        }
        "prde_cancel_liouvillian_tower" => {
            let b = arg(0)?;
            if !poly_in(&b, top)? {
                return Ok(None);
            }
            let (Some(q), Some(n)) = (polys(&arg(1)?, top)?, small_int(&arg(2)?)) else { return Ok(None) };
            Input::Cancel(b, q, n)
        }
        _ => return Ok(None),
    };
    let polymatrix = py.import("symox.polys.polymatrix")?;
    let inc = PyDict::new(py);
    inc.set_item("include", true)?;
    let fld = PyDict::new(py);
    fld.set_item("field", true)?;
    let z = Expr::dummy_with_index("z", Default::default(), num_bigint::BigInt::from(-(1i64 << 62)));
    let cx = Cx {
        py,
        qq: polymatrix.getattr("QQ")?,
        tl: tl.clone(),
        dl: dl.clone(),
        cases,
        cases_obj: cases_obj.downcast_into::<PyList>()?,
        nt,
        act: Cell::new(act),
        trace: RefCell::new(Vec::new()),
        steps: Cell::new(0),
        z: crate::wrap(z).into_pyobject(py)?.into_any(),
        inc,
        fld,
    };
    let spec = |m: &PM<'_>| -> PyResult<PyObject> {
        let items = PyList::new(py, m.e.iter())?;
        Ok((m.r, m.c, items, &m.gens, &m.dom).into_pyobject(py)?.into_any().unbind())
    };
    let result: R<PyObject> = (|| -> R<PyObject> {
        Ok(match input {
            Input::Risch(fa, fd, g) => {
                let (h, c) = cx.param_risch(&fa, &fd, &g)?;
                let h = h.iter().map(|(x, y)| PyTuple::new(py, [x, y])).collect::<PyResult<Vec<_>>>()?;
                (PyList::new(py, h)?, spec(&c)?).into_pyobject(py)?.into_any().unbind()
            }
            Input::Poly(a, b, q, n) => {
                let (h, c) = cx.param_poly(&a, &b, q, n)?;
                (PyList::new(py, h)?, spec(&c)?).into_pyobject(py)?.into_any().unbind()
            }
            Input::Cancel(b, q, n) => {
                let (h, c) = cx.cancel_liouvillian(&b, q, n)?;
                (PyList::new(py, h)?, spec(&c)?).into_pyobject(py)?.into_any().unbind()
            }
            Input::Limited(fa, fd, g) => match cx.limited_integrate(&fa, &fd, &g)? {
                None => py.None(),
                Some(((ya, yd), c)) => (PyTuple::new(py, [ya, yd])?, PyList::new(py, c)?).into_pyobject(py)?.into_any().unbind(),
            },
        })
    })();
    let (status, payload): (u8, PyObject) = match result {
        Ok(p) => (0, p),
        Err(Stop::Decline) => {
            if std::env::var_os("SYMOX_PRDE_DEBUG").is_some() {
                WHERE.with(|w| eprintln!("prde_param {kind} declined at {}", w.borrow()));
            }
            return Ok(None);
        }
        Err(Stop::Interrupt(e)) => return Err(e),
        Err(Stop::Raise(k, m)) => (k, m.into_pyobject(py)?.into_any().unbind()),
    };
    py.check_signals()?;
    let trace = PyList::new(py, cx.trace.borrow().iter())?;
    Ok(Some((status, payload, trace).into_pyobject(py)?.into_any().unbind()))
}
