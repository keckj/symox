//! Native `lie_group` solver.
//!
//! Port of SymPy 1.14 `sympy/solvers/ode/lie_group.py` (`infinitesimals`,
//! the nine `lie_heuristic_*` functions, `_lie_group_remove`,
//! `_ode_lie_group`, `_ode_lie_group_try_heuristic`), of `checkinfsol`
//! (`sympy/solvers/ode/ode.py`) and of `LieGroup._matches` /
//! `LieGroup._get_general_solution` (`sympy/solvers/ode/single.py`).
//!
//! The algorithm bodies run here: the determining equation and its
//! derivatives, the heuristic assumptions and their tests, the bivariate and
//! `chi` polynomial ansatzes and their coefficient systems, the canonical
//! coordinates (`_lie_group_remove`), the separation in canonical
//! coordinates, the quadrature equation and the heuristic order. `collect`,
//! `Expr.match`, `separatevars`, `powsimp`, `factor_terms`, `Poly`,
//! `as_numer_denom`, `diff`, `subs` and `doit` are the native ones.
//! `simplify`, `expand`, `cancel`, `integrate`, `solve`, `div` and `pdsolve`
//! are the public callables supplied by the binding through [`Hooks`]; an
//! error they raise propagates unchanged, and the `except` clauses of the
//! reference catch it by class ([`Hooks::exc`]). [`LieError::Decline`] leaves
//! the computation to the reference body.

use rustc_hash::FxHashSet;

use crate::diff::{diff, diff_n};
use crate::expr::{Expr, Func, Kind, RelOp};
use crate::functions::{func, undefined};
use crate::matches::{reference_has, reference_match};
use crate::ops::{add, add2, div, inv, mul, mul2, neg, pow, sub};
use crate::poly::collect::{self, Collection, Options};
use crate::relational::rel;
use crate::subs::{subs, subs1};
use crate::symbol::Assumptions;
use crate::traverse::{free_symbols, preorder};

/// `lie_heuristics`, in the reference order.
pub const HEURISTICS: [&str; 9] = [
    "abaco1_simple",
    "abaco1_product",
    "abaco2_similar",
    "abaco2_unique_unknown",
    "abaco2_unique_general",
    "linear",
    "function_sum",
    "bivariate",
    "chi",
];

/// The Python class of an exception, as far as the reference catches it or
/// raises it itself.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Exc {
    NotImplemented,
    Value,
    Type,
    Index,
    UnboundLocal,
    Other,
}

/// The result of `solve(eqs, *syms)` / `solve(eqs, syms)`: a dict or a list
/// of dicts (other shapes are not modelled and decline in the binding).
pub enum Solved {
    Dict(Vec<(Expr, Expr)>),
    Dicts(Vec<Vec<(Expr, Expr)>>),
}

/// The public callables of the reference bodies.
pub trait Hooks {
    type Error;
    /// The class of a raised error, for the reference `except` clauses.
    fn exc(&self, e: &Self::Error) -> Exc;
    fn simplify(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    fn expand(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    fn cancel(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    /// `integrate(e, x)`, or `integrate(e)` without a variable.
    fn integrate(&mut self, e: &Expr, x: Option<&Expr>) -> Result<Expr, Self::Error>;
    /// `solve(eq, sym)`: the list of solutions.
    fn solve(&mut self, eq: &Expr, sym: &Expr) -> Result<Vec<Expr>, Self::Error>;
    /// `solve(eqs, *syms, dict=True)`.
    fn solve_dicts(&mut self, eqs: &[Expr], syms: &[Expr]) -> Result<Vec<Vec<(Expr, Expr)>>, Self::Error>;
    /// `solve(eqs, *set(syms))` when `unpacked_set`, else `solve(eqs, syms)`.
    fn solve_system(&mut self, eqs: &[Expr], syms: &[Expr], unpacked_set: bool) -> Result<Solved, Self::Error>;
    /// `div(f, g)`: `(q, r)`.
    fn div(&mut self, f: &Expr, g: &Expr) -> Result<(Expr, Expr), Self::Error>;
    /// `pdsolve(eq, func=fxy).rhs`.
    fn pdsolve(&mut self, eq: &Expr, fxy: &Expr) -> Result<Expr, Self::Error>;
}

/// An error of the reference body: a hook error, an exception the reference
/// raises itself (its class and message), the `KeyError` of a missing key,
/// or a decline (the port does not reproduce the reference here).
#[derive(Debug)]
pub enum LieError<E> {
    Hook(E),
    Raise(Exc, String),
    KeyMissing(Expr),
    Decline(&'static str),
}

type R<T, H> = Result<T, LieError<<H as Hooks>::Error>>;

fn class<H: Hooks>(hooks: &H, e: &LieError<H::Error>) -> Exc {
    match e {
        LieError::Hook(e) => hooks.exc(e),
        LieError::Raise(k, _) => *k,
        LieError::KeyMissing(_) | LieError::Decline(_) => Exc::Other,
    }
}

/// `try: v = call() except NotImplementedError: None`.
fn unless_nie<H: Hooks, T>(hooks: &H, r: Result<T, H::Error>) -> R<Option<T>, H> {
    match r {
        Ok(v) => Ok(Some(v)),
        Err(e) if hooks.exc(&e) == Exc::NotImplemented => Ok(None),
        Err(e) => Err(LieError::Hook(e)),
    }
}

fn hk<H: Hooks, T>(r: Result<T, H::Error>) -> R<T, H> {
    r.map_err(LieError::Hook)
}

/// The reference truth value of an expression: false only for a numeric zero.
fn truthy(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => !n.is_zero(),
        _ => true,
    }
}

/// `e.diff(v, 2)`: the reference `Derivative` ends a derivative of order
/// above one with `factor_terms(signsimp(expr))`.
fn diff2(e: &Expr, v: &Expr) -> Expr {
    let r = diff_n(e, v, 2);
    if r.is_zero() || r.has_deferred_numeric_power() {
        return r;
    }
    crate::simptools::exprtools::factor_terms(&crate::simptools::misc::signsimp(&r, true), false, false, false, true)
}

fn has(e: &Expr, s: &Expr) -> bool {
    reference_has(e, s)
}

fn exp(e: &Expr) -> Expr {
    func(Func::Exp, vec![e.clone()])
}

fn index_error<E>() -> LieError<E> {
    LieError::Raise(Exc::Index, "list index out of range".into())
}

/// Content the port does not model: matrices, host classes and symbolic
/// derivatives.
fn unmodelled(e: &Expr) -> bool {
    e.has_matrix() || e.has_host() || e.has_symbolic_derivative()
}

/// Held operations other than derivatives (`Integral`, `Sum`, `Product`,
/// `Limit`, `Subs`, integral transforms): their `doit` stays with the
/// reference.
fn held(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Func(Func::Integral | Func::Sum | Func::Product | Func::Limit | Func::Subs, _)) || crate::transforms::as_held(s).is_some() {
            found = true;
        }
    });
    found
}

/// One infinitesimal `{xi(x, f(x)): xi, eta(x, f(x)): eta}`; `eta_first`
/// records the reference dict order. Equality is dict equality.
#[derive(Clone, Debug)]
pub struct Inf {
    pub xi: Expr,
    pub eta: Expr,
    pub eta_first: bool,
}

impl PartialEq for Inf {
    fn eq(&self, o: &Self) -> bool {
        self.xi == o.xi && self.eta == o.eta
    }
}

impl Inf {
    fn xe(xi: Expr, eta: Expr) -> Inf {
        Inf { xi, eta, eta_first: false }
    }
    fn ex(eta: Expr, xi: Expr) -> Inf {
        Inf { xi, eta, eta_first: true }
    }
}

/// The `match` dict of the heuristics.
#[derive(Clone, Debug)]
pub struct LieMatch {
    pub h: Expr,
    pub hx: Expr,
    pub hy: Expr,
    pub y: Expr,
    pub hinv: Expr,
    pub fx: Expr,
    pub x: Expr,
}

impl LieMatch {
    /// The dict `infinitesimals` builds from `h` (in `x` and the dummy `y`).
    pub fn new(h: &Expr, y: &Expr, fx: &Expr, x: &Expr) -> LieMatch {
        let u = Expr::dummy("u", Assumptions::default());
        let hinv = subs1(&subs(&inv(h), &[(x.clone(), u.clone()), (y.clone(), x.clone())]), &u, y);
        LieMatch { h: h.clone(), hx: diff(h, x), hy: diff(h, y), y: y.clone(), hinv, fx: fx.clone(), x: x.clone() }
    }
}

fn dummy(name: &str) -> Expr {
    Expr::dummy(name, Assumptions::default())
}

/// `-simplify(b/a)` of `collect(expand(eq), df).match(a*df + b)`, else the
/// first `solve(eq, df)`: the right-hand side `h` of `df = h` (still in
/// `f(x)`), as `infinitesimals` and `checkinfsol` compute it.
fn first_order_h<H: Hooks>(eq: &Expr, fx: &Expr, x: &Expr, hooks: &mut H) -> R<Expr, H> {
    let df = diff(fx, x);
    let expanded = hk::<H, _>(hooks.expand(eq))?;
    if unmodelled(&expanded) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let collected = match collect::run::<()>(&expanded, std::slice::from_ref(&df), Options { exact: false, discover: false, distribute_order_term: true }, true, None) {
        Ok(Collection::Expression(v)) => v,
        _ => return Err(LieError::Decline("collect")),
    };
    let a = Expr::wild("a", vec![df.clone()]);
    let b = Expr::wild("b", vec![df.clone()]);
    match reference_match(&collected, &add2(&mul2(&a, &df), &b), false) {
        Ok(Some(m)) if !m.is_empty() => {
            let (Some(va), Some(vb)) = (m.get(&a), m.get(&b)) else { return Err(LieError::Decline("match")) };
            let q = div(vb, va);
            Ok(neg(&hk::<H, _>(hooks.simplify(&q))?))
        }
        Ok(_) => match hooks.solve(eq, &df) {
            Ok(v) => v.into_iter().next().ok_or_else(index_error),
            Err(e) if hooks.exc(&e) == Exc::NotImplemented => {
                Err(LieError::Raise(Exc::NotImplemented, "Infinitesimals for the first order ODE could not be found".into()))
            }
            Err(e) => Err(LieError::Hook(e)),
        },
        Err(_) => Err(LieError::Decline("match")),
    }
}

/// `infinitesimals(eq, func, order=1, hint, match)` after the argument
/// checks: `given` is the `(h, y)` of a supplied `match`.
pub fn infinitesimals<H: Hooks>(eq: &Expr, fx: &Expr, x: &Expr, hint: &str, given: Option<(Expr, Expr)>, hooks: &mut H) -> R<Vec<Inf>, H> {
    let (h, y) = match given {
        Some(v) => v,
        None => {
            let h = first_order_h(eq, fx, x, hooks)?;
            let y = dummy("y");
            (subs1(&h, fx, &y), y)
        }
    };
    if unmodelled(&h) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let m = LieMatch::new(&h, &y, fx, x);
    match hint {
        "all" => {
            let mut xieta: Vec<Inf> = Vec::new();
            for name in HEURISTICS {
                if let Some(list) = heuristic(name, &m, true, hooks)? {
                    for inf in list {
                        if !xieta.contains(&inf) {
                            xieta.push(inf);
                        }
                    }
                }
            }
            if xieta.is_empty() {
                return Err(LieError::Raise(Exc::NotImplemented, "Infinitesimals could not be found for the given ODE".into()));
            }
            Ok(xieta)
        }
        "default" => {
            for name in HEURISTICS {
                if let Some(list) = heuristic(name, &m, false, hooks)? {
                    if !list.is_empty() {
                        return Ok(list);
                    }
                }
            }
            Err(LieError::Raise(Exc::NotImplemented, "Infinitesimals could not be found for the given ODE".into()))
        }
        name if !HEURISTICS.contains(&name) => Err(LieError::Raise(Exc::Value, format!("Heuristic not recognized: {name}"))),
        name => match heuristic(name, &m, true, hooks)? {
            Some(list) if !list.is_empty() => Ok(list),
            _ => Err(LieError::Raise(Exc::Value, "Infinitesimals could not be found using the given heuristic".into())),
        },
    }
}

/// `lie_heuristic_<name>(match, comp)`.
pub fn heuristic<H: Hooks>(name: &str, m: &LieMatch, comp: bool, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    match name {
        "abaco1_simple" => abaco1_simple(m, comp, hooks),
        "abaco1_product" => abaco1_product(m, comp, hooks),
        "bivariate" => bivariate(m, hooks),
        "chi" => chi(m, hooks),
        "function_sum" => function_sum(m, comp, hooks),
        "abaco2_similar" => abaco2_similar(m, hooks),
        "abaco2_unique_unknown" => abaco2_unique_unknown(m, hooks),
        "abaco2_unique_general" => abaco2_unique_general(m, hooks),
        "linear" => linear(m, hooks),
        _ => Err(LieError::Decline("unknown heuristic")),
    }
}

fn record(xieta: &mut Vec<Inf>, inf: Inf) {
    if !xieta.contains(&inf) {
        xieta.push(inf);
    }
}

fn abaco1_simple<H: Hooks>(m: &LieMatch, comp: bool, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hx, hy, fx) = (&m.x, &m.y, &m.h, &m.hx, &m.hy, &m.fx);
    let mut xieta = Vec::new();
    if !free_symbols(hy).contains(y) {
        let r = hooks.integrate(hy, Some(x));
        if let Some(v) = unless_nie(hooks, r)? {
            let inf = Inf::xe(Expr::zero(), exp(&v));
            if !comp {
                return Ok(Some(vec![inf]));
            }
            record(&mut xieta, inf);
        }
    }
    let factor = div(hy, h);
    if !free_symbols(&factor).contains(x) {
        let r = hooks.integrate(&factor, Some(y));
        if let Some(v) = unless_nie(hooks, r)? {
            let inf = Inf::xe(Expr::zero(), subs1(&exp(&v), y, fx));
            if !comp {
                return Ok(Some(vec![inf]));
            }
            record(&mut xieta, inf);
        }
    }
    let factor = div(&neg(hx), h);
    if !free_symbols(&factor).contains(y) {
        let r = hooks.integrate(&factor, Some(x));
        if let Some(v) = unless_nie(hooks, r)? {
            let inf = Inf::xe(exp(&v), Expr::zero());
            if !comp {
                return Ok(Some(vec![inf]));
            }
            record(&mut xieta, inf);
        }
    }
    let factor = div(&neg(hx), &pow(h, &Expr::int(2)));
    if !free_symbols(&factor).contains(x) {
        let r = hooks.integrate(&factor, Some(y));
        if let Some(v) = unless_nie(hooks, r)? {
            let inf = Inf::xe(subs1(&exp(&v), y, fx), Expr::zero());
            if !comp {
                return Ok(Some(vec![inf]));
            }
            record(&mut xieta, inf);
        }
    }
    Ok(if xieta.is_empty() { None } else { Some(xieta) })
}

/// `separatevars(e, dict=True, symbols=syms)`: `None` for the reference
/// `None`, else the `(key, value)` entries (`None` key for `'coeff'`).
fn separated(e: &Expr, syms: &[Expr]) -> Option<Vec<(Option<Expr>, Expr)>> {
    crate::simptools::misc::separatevars_dict(e, Some(syms), false)
}

fn entry<E>(sep: &[(Option<Expr>, Expr)], key: Option<&Expr>) -> Result<Expr, LieError<E>> {
    sep.iter().find(|(k, _)| k.as_ref() == key).map(|(_, v)| v.clone()).ok_or(LieError::Decline("separatevars key"))
}

fn abaco1_product<H: Hooks>(m: &LieMatch, comp: bool, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hinv, fx) = (&m.x, &m.y, &m.h, &m.hinv, &m.fx);
    let mut xieta = Vec::new();
    let two = Expr::int(2);
    let xy = [x.clone(), y.clone()];
    let lh = diff(&diff(&func(Func::Log, vec![h.clone()]), y), x);
    if let Some(inf) = separated(&div(&lh, &pow(h, &two)), &xy) {
        if !inf.is_empty() && truthy(&entry(&inf, None)?) {
            let fxv = entry(&inf, Some(x))?;
            let gy = hk::<H, _>(hooks.simplify(&mul2(&fxv, &diff(&inv(&mul2(&fxv, h)), x))))?;
            if !free_symbols(&gy).contains(x) {
                let gy = exp(&hk::<H, _>(hooks.integrate(&gy, Some(y)))?);
                let inf = Inf::ex(Expr::zero(), subs1(&mul2(&fxv, &gy), y, fx));
                if !comp {
                    return Ok(Some(vec![inf]));
                }
                record(&mut xieta, inf);
            }
        }
    }
    let u1 = dummy("u1");
    let lh = diff(&diff(&func(Func::Log, vec![hinv.clone()]), y), x);
    if let Some(inf) = separated(&div(&lh, &pow(hinv, &two)), &xy) {
        if !inf.is_empty() && truthy(&entry(&inf, None)?) {
            let fxv = entry(&inf, Some(x))?;
            let gy = hk::<H, _>(hooks.simplify(&mul2(&fxv, &diff(&inv(&mul2(&fxv, hinv)), x))))?;
            if !free_symbols(&gy).contains(x) {
                let gy = exp(&hk::<H, _>(hooks.integrate(&gy, Some(y)))?);
                let etaval = mul2(&fxv, &gy);
                let etaval = subs1(&subs(&etaval, &[(x.clone(), u1.clone()), (y.clone(), x.clone())]), &u1, y);
                let inf = Inf::ex(subs1(&etaval, y, fx), Expr::zero());
                if !comp {
                    return Ok(Some(vec![inf]));
                }
                record(&mut xieta, inf);
            }
        }
    }
    Ok(if xieta.is_empty() { None } else { Some(xieta) })
}

fn is_rational_function(e: &Expr) -> bool {
    crate::inspection::classify(e, &free_symbols(e), crate::inspection::PolynomialKind::Rational) == Some(true)
}

fn is_polynomial_in(e: &Expr, gens: &[Expr]) -> bool {
    crate::inspection::classify(e, gens, crate::inspection::PolynomialKind::Polynomial) == Some(true)
}

/// `Poly(e, *gens)`; a construction the native ring does not cover declines.
fn poly<E>(e: &Expr, gens: &[Expr]) -> Result<crate::poly::polyobj::Poly, LieError<E>> {
    crate::poly::polyobj::Poly::from_expr(e, gens).map_err(|_| LieError::Decline("Poly"))
}

/// `Poly(e, *gens).total_degree()`.
fn total_degree<E>(e: &Expr, gens: &[Expr]) -> Result<u64, LieError<E>> {
    poly(e, gens)?.total_degree().ok_or(LieError::Decline("total degree of zero"))
}

/// `Add(*[Symbol(prefix + '_' + p + '_' + (i-p)) * x**p * y**(i-p) for p in range(i+1)])`.
fn homogeneous_ansatz(prefix: &str, i: u64, x: &Expr, y: &Expr) -> Expr {
    let terms: Vec<Expr> = (0..=i)
        .map(|p| {
            let c = Expr::symbol(&format!("{prefix}_{p}_{}", i - p));
            mul(&[c, pow(x, &Expr::int(p as i64)), pow(y, &Expr::int((i - p) as i64))])
        })
        .collect();
    add(&terms)
}

/// The dict of a `solve(...)` result after `if isinstance(soldict, list):
/// soldict = soldict[0]`.
fn first_dict<E>(s: Solved) -> Result<Vec<(Expr, Expr)>, LieError<E>> {
    match s {
        Solved::Dict(d) => Ok(d),
        Solved::Dicts(l) => l.into_iter().next().ok_or_else(index_error),
    }
}

fn subs_pairs(e: &Expr, pairs: &[(Expr, Expr)]) -> Expr {
    subs(e, pairs)
}

fn ones(syms: &[Expr]) -> Vec<(Expr, Expr)> {
    syms.iter().map(|s| (s.clone(), Expr::one())).collect()
}

/// `xieq.free_symbols.union(etaeq.free_symbols) - {x, y}`, in first-seen
/// order (the binding iterates it as a Python set).
fn unknowns(parts: &[&Expr], x: &Expr, y: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    for p in parts {
        for s in free_symbols(p) {
            if s != *x && s != *y && !out.contains(&s) {
                out.push(s);
            }
        }
    }
    out
}

fn bivariate<H: Hooks>(m: &LieMatch, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hx, hy, fx) = (&m.x, &m.y, &m.h, &m.hx, &m.hy, &m.fx);
    if !is_rational_function(h) {
        return Ok(None);
    }
    let two = Expr::int(2);
    let s = |n: &str| Expr::symbol(n);
    let (etax, etay, etad, xix, xiy, xid) = (s("etax"), s("etay"), s("etad"), s("xix"), s("xiy"), s("xid"));
    let ipde = add(&[
        etax,
        mul2(&sub(&etay, &xix), h),
        neg(&mul2(&xiy, &pow(h, &two))),
        neg(&mul2(&xid, hx)),
        neg(&mul2(&etad, hy)),
    ]);
    let canceled = hk::<H, _>(hooks.cancel(&ipde))?;
    let (num, _) = crate::simptools::as_numer_denom(&canceled);
    let xy = [x.clone(), y.clone()];
    let deg = total_degree(&num, &xy)?;
    let mut xieq = s("xi0");
    let mut etaeq = s("eta0");
    let mut polyy: Option<Vec<Expr>> = None;
    for i in 0..=deg {
        if i > 0 {
            xieq = add2(&xieq, &homogeneous_ansatz("xi", i, x, y));
            etaeq = add2(&etaeq, &homogeneous_ansatz("eta", i, x, y));
        }
        // ipde.subs({dxi: xieq, deta: etaeq}).doit(): the derivatives of the
        // polynomial ansatz (h is a rational function, nothing else is held)
        let val = add(&[
            diff(&etaeq, x),
            mul2(&sub(&diff(&etaeq, y), &diff(&xieq, x)), h),
            neg(&mul2(&diff(&xieq, y), &pow(h, &two))),
            neg(&mul2(&xieq, hx)),
            neg(&mul2(&etaeq, hy)),
        ]);
        let (pden, _) = crate::simptools::as_numer_denom(&val);
        let pden = hk::<H, _>(hooks.expand(&pden))?;
        if is_polynomial_in(&pden, &xy) && pden.is_add() {
            polyy = Some(poly(&pden, &xy)?.terms().into_iter().map(|(_, c)| c).collect());
        }
        let Some(values) = &polyy else {
            return Err(LieError::Raise(Exc::UnboundLocal, "local variable 'polyy' referenced before assignment".into()));
        };
        if !values.is_empty() {
            let symset = unknowns(&[&xieq, &etaeq], x, y);
            let solved = hk::<H, _>(hooks.solve_system(values, &symset, true))?;
            let soldict = first_dict(solved)?;
            if soldict.iter().any(|(_, v)| truthy(v)) {
                let xired = subs_pairs(&xieq, &soldict);
                let etared = subs_pairs(&etaeq, &soldict);
                let one = ones(&symset);
                let inf = Inf::ex(subs1(&subs(&etared, &one), y, fx), subs1(&subs(&xired, &one), y, fx));
                return Ok(Some(vec![inf]));
            }
        }
    }
    Ok(None)
}

fn chi<H: Hooks>(m: &LieMatch, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hy, fx) = (&m.x, &m.y, &m.h, &m.hy, &m.fx);
    if !is_rational_function(h) {
        return Ok(None);
    }
    let s = |n: &str| Expr::symbol(n);
    let (schi, schix, schiy) = (s("schi"), s("schix"), s("schiy"));
    let cpde = add(&[schix, mul2(h, &schiy), neg(&mul2(hy, &schi))]);
    let canceled = hk::<H, _>(hooks.cancel(&cpde))?;
    let (num, _) = crate::simptools::as_numer_denom(&canceled);
    let xy = [x.clone(), y.clone()];
    let deg = total_degree(&num, &xy)?;
    let mut chieq = s("chi");
    for i in 1..=deg {
        chieq = add2(&chieq, &homogeneous_ansatz("chi", i, x, y));
        let val = add(&[diff(&chieq, x), mul2(h, &diff(&chieq, y)), neg(&mul2(hy, &chieq))]);
        let canceled = hk::<H, _>(hooks.cancel(&val))?;
        let (cnum, _) = crate::simptools::as_numer_denom(&canceled);
        let cnum = hk::<H, _>(hooks.expand(&cnum))?;
        if is_polynomial_in(&cnum, &xy) && cnum.is_add() {
            let cpoly: Vec<Expr> = poly(&cnum, &xy)?.terms().into_iter().map(|(_, c)| c).collect();
            if !cpoly.is_empty() {
                let solsyms = unknowns(&[&chieq], x, y);
                let solved = hk::<H, _>(hooks.solve_system(&cpoly, &solsyms, true))?;
                let soldict = first_dict(solved)?;
                if soldict.iter().any(|(_, v)| truthy(v)) {
                    let chieq = subs(&subs_pairs(&chieq, &soldict), &ones(&solsyms));
                    let (xic, etac) = hk::<H, _>(hooks.div(&chieq, h))?;
                    let inf = Inf::ex(subs1(&etac, y, fx), neg(&subs1(&xic, y, fx)));
                    return Ok(Some(vec![inf]));
                }
            }
        }
    }
    Ok(None)
}

fn function_sum<H: Hooks>(m: &LieMatch, comp: bool, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hinv, fx) = (&m.x, &m.y, &m.h, &m.hinv, &m.fx);
    let mut xieta = Vec::new();
    let xy = [x.clone(), y.clone()];
    for odefac in [h, hinv] {
        let factor = mul2(odefac, &diff2(&div(&Expr::one(), odefac), x));
        let sep = separated(&diff(&div(&Expr::one(), &factor), y), &xy);
        if let Some(sep) = sep.filter(|s| !s.is_empty()) {
            let coeff = entry(&sep, None)?;
            let (sx, sy) = (entry(&sep, Some(x))?, entry(&sep, Some(y))?);
            if truthy(&coeff) && has(&sx, x) && has(&sy, y) {
                let k = dummy("k");
                let r = hooks.integrate(&sy, Some(y));
                if let Some(int) = unless_nie(hooks, r)? {
                    let mut gy = mul2(&k, &int);
                    let fdd = inv(&mul(&[k.clone(), sx.clone(), coeff.clone()]));
                    let mut fxv = hk::<H, _>(hooks.simplify(&sub(&div(&fdd, &factor), &gy)))?;
                    let check = hk::<H, _>(hooks.simplify(&sub(&diff2(&fxv, x), &fdd)))?;
                    if truthy(&fxv) {
                        if !truthy(&check) {
                            fxv = subs1(&fxv, &k, &Expr::one());
                            gy = div(&gy, &k);
                        } else {
                            let sol = hk::<H, _>(hooks.solve(&check, &k))?;
                            let Some(sol) = sol.into_iter().next() else { continue };
                            fxv = subs1(&fxv, &k, &sol);
                            gy = mul2(&div(&gy, &k), &sol);
                        }
                        let swapped = odefac == hinv;
                        if swapped {
                            fxv = subs1(&fxv, x, y);
                            gy = subs1(&gy, y, x);
                        }
                        let mut etaval = crate::simptools::exprtools::factor_terms(&add2(&fxv, &gy), false, false, false, true);
                        if etaval.is_mul() {
                            let kept: Vec<Expr> = etaval.args().into_iter().filter(|a| has(a, x) || has(a, y)).collect();
                            etaval = mul(&kept);
                        }
                        let inf = if swapped { Inf::ex(subs1(&etaval, y, fx), Expr::zero()) } else { Inf::xe(subs1(&etaval, y, fx), Expr::zero()) };
                        if !comp {
                            return Ok(Some(vec![inf]));
                        }
                        xieta.push(inf);
                    }
                }
            }
        }
        if !xieta.is_empty() {
            return Ok(Some(xieta));
        }
    }
    Ok(None)
}

fn abaco2_similar<H: Hooks>(m: &LieMatch, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hx, hy, hinv, fx) = (&m.x, &m.y, &m.h, &m.hx, &m.hy, &m.hinv, &m.fx);
    for swapped in [false, true] {
        let g = if swapped { hinv } else { h };
        let (gx, gy) = if swapped { (diff(hinv, x), diff(hinv, y)) } else { (hx.clone(), hy.clone()) };
        let factor = hk::<H, _>(hooks.cancel(&div(&diff(g, y), &diff2(g, y))))?;
        let factorx = diff(&factor, x);
        let factory = diff(&factor, y);
        if !has(&factor, x) && !has(&factor, y) {
            // The reference matches h (not hinv) in both branches and then
            // calls exp with two arguments: only a NotImplementedError of
            // the integration lets the heuristic continue.
            let a = Expr::wild("A", vec![y.clone()]);
            let b = Expr::wild("B", vec![y.clone()]);
            let c = Expr::wild("C", vec![x.clone(), y.clone()]);
            let pattern = add2(&a, &mul2(&b, &exp(&div(y, &c))));
            let found = match reference_match(h, &pattern, false) {
                Ok(v) => v,
                Err(_) => return Err(LieError::Decline("match")),
            };
            let Some(found) = found else {
                return Err(LieError::Raise(Exc::Type, "'NoneType' object is not subscriptable".into()));
            };
            // `match[A]/match[C]`: a missing key raises KeyError in that order.
            let va = found.get(&a).ok_or_else(|| LieError::KeyMissing(a.clone()))?;
            let vc = found.get(&c).ok_or_else(|| LieError::KeyMissing(c.clone()))?;
            let r = hooks.integrate(&div(va, vc), None);
            if unless_nie(hooks, r)?.is_some() {
                return Err(LieError::Raise(Exc::Type, "exp takes exactly 1 argument (2 given)".into()));
            }
        } else {
            let gamma = hk::<H, _>(hooks.cancel(&div(&factorx, &factory)))?;
            if !has(&gamma, y) {
                let num = sub(&sub(&mul2(&gamma, &gy), &diff(&gamma, x)), &gx);
                let tauint = hk::<H, _>(hooks.cancel(&div(&num, &add2(g, &gamma))))?;
                if !has(&tauint, y) {
                    let r = hooks.integrate(&tauint, Some(x));
                    if let Some(int) = unless_nie(hooks, r)? {
                        let tau = exp(&int);
                        let gxv = neg(&mul2(&tau, &gamma));
                        let inf = if swapped { Inf::ex(subs1(&tau, x, fx), subs1(&gxv, x, fx)) } else { Inf::xe(tau, gxv) };
                        return Ok(Some(vec![inf]));
                    }
                }
            }
        }
    }
    Ok(None)
}

/// `h.atoms(Pow)` / `h.atoms(AppliedUndef)` in first preorder occurrence
/// (the reference iterates a set, whose order depends on the hash seed).
fn atoms_where(e: &Expr, pred: impl Fn(&Expr) -> bool) -> Vec<Expr> {
    let mut seen = FxHashSet::default();
    let mut out = Vec::new();
    preorder(e, &mut |s| {
        if pred(s) && seen.insert(s.clone()) {
            out.push(s.clone());
        }
    });
    out
}

fn abaco2_unique_unknown<H: Hooks>(m: &LieMatch, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hx, hy, fx) = (&m.x, &m.y, &m.h, &m.hx, &m.hy, &m.fx);
    let two = Expr::int(2);
    let mut funclist = Vec::new();
    for atom in atoms_where(h, |s| s.is_pow()) {
        let (base, e) = crate::decompose::as_base_exp(&atom);
        if has(&base, x) && has(&base, y) && !e.is_integer_number() {
            funclist.push(atom);
        }
    }
    for function in atoms_where(h, |s| matches!(s.kind(), Kind::Func(Func::Undefined(_), _))) {
        let syms = free_symbols(&function);
        if syms.contains(x) && syms.contains(y) {
            funclist.push(function);
        }
    }
    let xy = [x.clone(), y.clone()];
    for f in funclist {
        let frac = hk::<H, _>(hooks.cancel(&div(&diff(&f, y), &diff(&f, x))))?;
        let sep = separated(&frac, &xy).filter(|s| !s.is_empty());
        let coeff = match &sep {
            Some(s) => Some(entry(s, None)?),
            None => None,
        };
        if let (Some(sep), Some(coeff)) = (&sep, coeff.filter(truthy)) {
            let xitry1 = entry(sep, Some(x))?;
            let etatry1 = neg(&inv(&mul2(&entry(sep, Some(y))?, &coeff)));
            let pde1 = add(&[
                mul2(&diff(&etatry1, y), h),
                neg(&mul2(&diff(&xitry1, x), h)),
                neg(&mul2(&xitry1, hx)),
                neg(&mul2(&etatry1, hy)),
            ]);
            if !truthy(&hk::<H, _>(hooks.simplify(&pde1))?) {
                return Ok(Some(vec![Inf::xe(xitry1, subs1(&etatry1, y, fx))]));
            }
            let xitry2 = inv(&etatry1);
            let etatry2 = inv(&xitry1);
            let pde2 = add(&[
                diff(&etatry2, x),
                neg(&mul2(&diff(&xitry2, y), &pow(h, &two))),
                neg(&mul2(&xitry2, hx)),
                neg(&mul2(&etatry2, hy)),
            ]);
            let expanded = hk::<H, _>(hooks.expand(&pde2))?;
            if !truthy(&hk::<H, _>(hooks.simplify(&expanded))?) {
                return Ok(Some(vec![Inf::xe(subs1(&xitry2, y, fx), etatry2)]));
            }
        } else {
            let etatry = neg(&inv(&frac));
            let pde = add(&[diff(&etatry, x), mul2(&diff(&etatry, y), h), neg(hx), neg(&mul2(&etatry, hy))]);
            if !truthy(&hk::<H, _>(hooks.simplify(&pde))?) {
                return Ok(Some(vec![Inf::xe(Expr::one(), subs1(&etatry, y, fx))]));
            }
            let xitry = neg(&frac);
            let pde = add(&[
                neg(&mul2(&diff(&xitry, x), h)),
                neg(&mul2(&diff(&xitry, y), &pow(h, &two))),
                neg(&mul2(&xitry, hx)),
                neg(hy),
            ]);
            let expanded = hk::<H, _>(hooks.expand(&pde))?;
            if !truthy(&hk::<H, _>(hooks.simplify(&expanded))?) {
                return Ok(Some(vec![Inf::xe(subs1(&xitry, y, fx), Expr::one())]));
            }
        }
    }
    Ok(None)
}

/// The `TypeError` of the reference `x not in etaval`.
fn not_iterable<E>(e: &Expr) -> LieError<E> {
    LieError::Raise(Exc::Type, format!("argument of type '{}' is not iterable", e.class_name()))
}

fn abaco2_unique_general<H: Hooks>(m: &LieMatch, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, hx, hy) = (&m.x, &m.y, &m.hx, &m.hy);
    let n = Expr::int;
    let two = n(2);
    let a = diff(hx, y);
    let b = add2(&diff(hy, y), &pow(hy, &two));
    let c = sub(&diff(hx, x), &pow(hx, &two));
    if !(truthy(&a) && truthy(&b) && truthy(&c)) {
        return Ok(None);
    }
    let ax = diff(&a, x);
    let ay = diff(&a, y);
    let axy = diff(&ax, y);
    let axx = diff(&ax, x);
    let ayy = diff(&ay, y);
    let inner = add(&[mul2(&two, &axy), mul2(hx, &ay), neg(&mul2(&ax, hy)), mul2(&add2(&mul2(hx, hy), &mul2(&two, &a)), &a)]);
    let d = sub(&mul2(&hk::<H, _>(hooks.simplify(&inner))?, &a), &mul(&[n(3), ax.clone(), ay.clone()]));
    let a3 = pow(&a, &n(3));
    if !truthy(&d) {
        let e1 = hk::<H, _>(hooks.simplify(&add2(
            &mul2(&n(3), &pow(&ax, &two)),
            &mul2(&sub(&mul2(&add2(&pow(hx, &two), &mul2(&two, &c)), &a), &mul2(&two, &axx)), &a),
        )))?;
        if truthy(&e1) {
            let e2 = hk::<H, _>(hooks.simplify(&sub(
                &mul2(&add2(&mul2(&two, &ayy), &mul2(&sub(&mul2(&two, &b), &pow(hy, &two)), &a)), &a),
                &mul2(&n(3), &pow(&ay, &two)),
            )))?;
            if !truthy(&e2) {
                let left = mul2(&e1, &sub(&mul2(&add2(&mul2(&n(28), &ax), &mul(&[n(4), hx.clone(), a.clone()])), &a3), &mul2(&e1, &add2(&mul2(hy, &a), &ay))));
                let e3 = hk::<H, _>(hooks.simplify(&sub(&left, &mul(&[diff(&e1, x), n(8), pow(&a, &n(4))]))))?;
                if !truthy(&e3) {
                    let num = add2(&mul(&[n(4), a3.clone(), sub(&ax, &mul2(hx, &a))]), &mul2(&e1, &sub(&mul2(hy, &a), &ay)));
                    let etaval = hk::<H, _>(hooks.cancel(&div(&num, &mul(&[two.clone(), a.clone(), e1.clone()]))))?;
                    return Err(not_iterable(&etaval));
                }
            }
        }
    } else {
        let e1 = hk::<H, _>(hooks.simplify(&sub(
            &mul2(&add2(&mul2(&two, &ayy), &mul2(&sub(&mul2(&two, &b), &pow(hy, &two)), &a)), &a),
            &mul2(&n(3), &pow(&ay, &two)),
        )))?;
        if truthy(&e1) {
            let tail = sub(&mul2(&sub(&mul2(&two, &axx), &mul2(&add2(&pow(hx, &two), &mul2(&two, &c)), &a)), &a), &mul2(&n(3), &pow(&ax, &two)));
            let e2 = hk::<H, _>(hooks.simplify(&add(&[mul(&[n(4), a3.clone(), d.clone()]), neg(&pow(&d, &two)), mul2(&e1, &tail)])))?;
            if !truthy(&e2) {
                let first = mul2(&neg(&mul2(&a, &d)), &diff(&e1, y));
                let bracket = add(&[
                    mul2(&sub(&diff(&e1, x), &mul2(hy, &d)), &a),
                    mul(&[n(3), ay.clone(), d.clone()]),
                    mul2(&sub(&mul2(&a, hx), &mul2(&n(3), &ax)), &e1),
                ]);
                let e3 = hk::<H, _>(hooks.simplify(&add2(&first, &mul2(&bracket, &e1))))?;
                if !truthy(&e3) {
                    let num = sub(&mul2(&sub(&mul2(&a, hx), &ax), &e1), &mul2(&add2(&ay, &mul2(&a, hy)), &d));
                    let etaval = hk::<H, _>(hooks.cancel(&div(&num, &mul(&[two.clone(), a.clone(), d.clone()]))))?;
                    return Err(not_iterable(&etaval));
                }
            }
        }
    }
    Ok(None)
}

fn linear<H: Hooks>(m: &LieMatch, hooks: &mut H) -> R<Option<Vec<Inf>>, H> {
    let (x, y, h, hx, hy, fx) = (&m.x, &m.y, &m.h, &m.hx, &m.hy, &m.fx);
    // `[next(symbols) for _ in islice(symbols, 6)]` takes every other
    // numbered dummy: c1, c3, ..., c11.
    let syms: Vec<Expr> = (0..6).map(|i| dummy(&format!("c{}", 2 * i + 1))).collect();
    let (c0, c1, c2, c3, c4, c5) = (&syms[0], &syms[1], &syms[2], &syms[3], &syms[4], &syms[5]);
    let pde = add(&[
        c3.clone(),
        mul2(&sub(c4, c0), h),
        neg(&mul2(&add(&[mul2(c0, x), mul2(c1, y), c2.clone()]), hx)),
        neg(&mul2(&add(&[mul2(c3, x), mul2(c4, y), c5.clone()]), hy)),
        neg(&mul2(c1, &pow(h, &Expr::int(2)))),
    ]);
    let (pde, _) = crate::simptools::as_numer_denom(&pde);
    let pde = crate::simptools::powsimp::powsimp(&hk::<H, _>(hooks.expand(&pde))?, &Default::default());
    let mut coeffdict: Vec<(Expr, Expr)> = Vec::new();
    let mut accumulate = |key: Expr, value: Expr| match coeffdict.iter_mut().find(|(k, _)| *k == key) {
        Some((_, v)) => *v = add2(v, &value),
        None => coeffdict.push((key, value)),
    };
    if pde.is_add() {
        for term in pde.args() {
            if term.is_mul() {
                let rem = mul(&term.args().into_iter().filter(|f| !(has(f, x) || has(f, y))).collect::<Vec<_>>());
                accumulate(div(&term, &rem), rem);
            } else {
                accumulate(term, Expr::one());
            }
        }
    }
    let sollist: Vec<Expr> = coeffdict.into_iter().map(|(_, v)| v).collect();
    let soldict = match hk::<H, _>(hooks.solve_system(&sollist, &syms, false))? {
        Solved::Dict(d) if d.is_empty() => return Ok(None),
        Solved::Dicts(l) if l.is_empty() => return Ok(None),
        s => first_dict(s)?,
    };
    if !soldict.iter().any(|(_, v)| truthy(v)) {
        return Ok(None);
    }
    let one = ones(&syms);
    let xival = subs(&subs_pairs(&add(&[mul2(c0, x), mul2(c1, fx), c2.clone()]), &soldict), &one);
    let etaval = subs(&subs_pairs(&add(&[mul2(c3, x), mul2(c4, fx), c5.clone()]), &soldict), &one);
    Ok(Some(vec![Inf::xe(xival, etaval)]))
}

/// `_lie_group_remove(coords)`.
pub fn lie_group_remove(coords: &Expr) -> Expr {
    match coords.kind() {
        Kind::Func(Func::Undefined(_), args) => args[0].clone(),
        Kind::Add(..) => {
            let funcs = atoms_where(coords, |s| matches!(s.kind(), Kind::Func(Func::Undefined(_), _)));
            let mut out = coords.clone();
            for f in funcs {
                out = subs1(&out, &f, &Expr::zero());
            }
            out
        }
        Kind::Pow(..) => {
            let (base, e) = crate::decompose::as_base_exp(coords);
            pow(&lie_group_remove(&base), &lie_group_remove(&e))
        }
        Kind::Mul(..) => mul(&coords.args().iter().map(lie_group_remove).collect::<Vec<_>>()),
        _ => coords.clone(),
    }
}

/// `checkinfsol(eq, infinitesimals, func, order=1)` after the argument
/// checks, for the `(xi, eta)` values of each infinitesimal.
pub fn checkinfsol<H: Hooks>(eq: &Expr, fx: &Expr, x: &Expr, infs: &[(Expr, Expr)], hooks: &mut H) -> R<Vec<(bool, Expr)>, H> {
    let h = first_order_h(eq, fx, x, hooks)?;
    let y = dummy("y");
    let h = subs1(&h, fx, &y);
    if unmodelled(&h) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let xi = undefined("xi", vec![x.clone(), y.clone()]);
    let eta = undefined("eta", vec![x.clone(), y.clone()]);
    let pde = add(&[
        diff(&eta, x),
        mul2(&sub(&diff(&eta, &y), &diff(&xi, x)), &h),
        neg(&mul2(&diff(&xi, &y), &pow(&h, &Expr::int(2)))),
        neg(&mul2(&xi, &diff(&h, x))),
        neg(&mul2(&eta, &diff(&h, &y))),
    ]);
    let mut out = Vec::new();
    for (vxi, veta) in infs {
        let (vxi, veta) = (subs1(vxi, fx, &y), subs1(veta, fx, &y));
        if unmodelled(&vxi) || unmodelled(&veta) || has(&vxi, &xi) || has(&vxi, &eta) || has(&veta, &xi) || has(&veta, &eta) {
            return Err(LieError::Decline("unmodelled infinitesimal"));
        }
        // pde.subs({xi: ..., eta: ...}).doit()
        let substituted = subs(&pde, &[(eta.clone(), veta), (xi.clone(), vxi)]);
        if unmodelled(&substituted) || held(&substituted) {
            return Err(LieError::Decline("unmodelled content"));
        }
        let done = crate::functions::doit(&substituted, true);
        let sol = hk::<H, _>(hooks.simplify(&done))?;
        out.push(if truthy(&sol) { (false, subs1(&sol, &y, fx)) } else { (true, Expr::zero()) });
    }
    Ok(out)
}

/// `LieGroup._matches` after the order and parameter checks: the values of
/// `d` and `e` in `collect(eq, df, exact=True).match(d + e*df)` with `f(x)`
/// replaced by `y`, `Ok(None)` without a match.
pub fn lie_group_split(eq: &Expr, fx: &Expr, x: &Expr, y: &Expr) -> Result<Option<(Expr, Expr)>, &'static str> {
    if unmodelled(eq) {
        return Err("unmodelled content");
    }
    let df = diff(fx, x);
    let d2 = diff_n(fx, x, 2);
    let collected = match collect::run::<()>(eq, std::slice::from_ref(&df), Options { exact: true, discover: false, distribute_order_term: true }, true, None) {
        Ok(Collection::Expression(v)) => v,
        _ => return Err("collect"),
    };
    let d = Expr::wild("d", vec![df.clone(), d2]);
    let e = Expr::wild("e", vec![df.clone()]);
    match reference_match(&collected, &add2(&d, &mul2(&e, &df)), false) {
        Ok(Some(m)) if !m.is_empty() => match (m.get(&d), m.get(&e)) {
            (Some(vd), Some(ve)) => Ok(Some((subs1(vd, fx, y), subs1(ve, fx, y)))),
            _ => Err("match"),
        },
        Ok(_) => Ok(None),
        Err(_) => Err("match"),
    }
}

/// The function `f` of `f(x)` applied to `args`.
fn apply_same(fx: &Expr, args: Vec<Expr>) -> Option<Expr> {
    match fx.kind() {
        Kind::Func(head @ Func::Undefined(_), _) => Some(func(head.clone(), args)),
        _ => None,
    }
}

/// `_ode_lie_group_try_heuristic(Eq(df, s), heuristic, func, {h, y}, inf)`.
fn try_heuristic<H: Hooks>(name: &str, fx: &Expr, x: &Expr, h: &Expr, y: &Expr, inf: &[Inf], hooks: &mut H) -> R<Option<Vec<Expr>>, H> {
    let infs = if inf.is_empty() {
        match infinitesimals(&Expr::zero(), fx, x, name, Some((h.clone(), y.clone())), hooks) {
            Ok(v) => v,
            Err(e) if class(hooks, &e) == Exc::Value => return Ok(None),
            Err(e) => return Err(e),
        }
    } else {
        inf.to_vec()
    };
    let fxy = apply_same(fx, vec![x.clone(), y.clone()]).ok_or(LieError::Decline("function head"))?;
    let c1 = Expr::symbol("C1");
    let mut tempsol: Vec<Expr> = Vec::new();
    for infsim in infs {
        let xiinf = subs1(&infsim.xi, fx, y);
        let etainf = subs1(&infsim.eta, fx, y);
        if hk::<H, _>(hooks.simplify(&div(&etainf, &xiinf)))? == *h {
            continue;
        }
        let rpde = add2(&mul2(&diff(&fxy, x), &xiinf), &mul2(&diff(&fxy, y), &etainf));
        let r = hk::<H, _>(hooks.pdsolve(&rpde, &fxy))?;
        let s = hk::<H, _>(hooks.pdsolve(&sub(&rpde, &Expr::one()), &fxy))?;
        let rcoord = lie_group_remove(&r);
        let scoord = lie_group_remove(&s);
        let (rd, sd) = (dummy("r"), dummy("s"));
        let system = [sub(&rd, &rcoord), sub(&sd, &scoord)];
        let r = hooks.solve_dicts(&system, &[x.clone(), y.clone()]);
        let Some(sol) = unless_nie(hooks, r)? else { continue };
        let Some(sol) = sol.into_iter().next() else { continue };
        let lookup = |k: &Expr| sol.iter().find(|(key, _)| key == k).map(|(_, v)| v.clone()).ok_or_else(|| LieError::KeyMissing(k.clone()));
        let xsub = lookup(x)?;
        let ysub = lookup(y)?;
        let num = hk::<H, _>(hooks.simplify(&add2(&diff(&scoord, x), &mul2(&diff(&scoord, y), h))))?;
        let denom = hk::<H, _>(hooks.simplify(&add2(&diff(&rcoord, x), &mul2(&diff(&rcoord, y), h))))?;
        if truthy(&num) && truthy(&denom) {
            let diffeq = hk::<H, _>(hooks.simplify(&subs(&div(&num, &denom), &[(x.clone(), xsub), (y.clone(), ysub)])))?;
            if let Some(sep) = separated(&diffeq, &[rd.clone(), sd.clone()]).filter(|s| !s.is_empty()) {
                let (ss, sr, coeff) = (entry(&sep, Some(&sd))?, entry(&sep, Some(&rd))?, entry(&sep, None)?);
                let is = hk::<H, _>(hooks.integrate(&inv(&ss), Some(&sd)))?;
                let ir = hk::<H, _>(hooks.integrate(&mul2(&coeff, &sr), Some(&rd)))?;
                let deq = sub(&add2(&is, &c1), &ir);
                let deq = subs(&deq, &[(rd.clone(), rcoord.clone()), (sd.clone(), scoord.clone())]);
                let r = hooks.solve(&deq, y);
                match unless_nie(hooks, r)? {
                    None => tempsol.push(deq),
                    Some(sdeq) => return Ok(Some(sdeq.iter().map(|s| rel(RelOp::Eq, fx, s)).collect())),
                }
            }
        } else if truthy(&denom) {
            let sol = hk::<H, _>(hooks.solve(&sub(&scoord, &c1), y))?;
            return Ok(Some(vec![rel(RelOp::Eq, fx, &sol.into_iter().next().ok_or_else(index_error)?)]));
        } else if truthy(&num) {
            let sol = hk::<H, _>(hooks.solve(&sub(&rcoord, &c1), y))?;
            return Ok(Some(vec![rel(RelOp::Eq, fx, &sol.into_iter().next().ok_or_else(index_error)?)]));
        }
    }
    if tempsol.is_empty() {
        return Ok(None);
    }
    Ok(Some(tempsol.iter().map(|s| rel(RelOp::Eq, &subs1(s, y, fx), &Expr::zero())).collect()))
}

/// The user infinitesimal `(xi, eta)` of the solver parameters and the
/// `(y, d, e)` entries of `LieGroup._matches`.
pub struct Params {
    pub user: Option<(Expr, Expr)>,
    pub split: Option<(Expr, Expr, Expr)>,
}

/// `_ode_lie_group(s, func, order=1, match)`: the solutions, `Ok(None)` for
/// the reference `None` (an empty list stands for the reference `[]`).
fn ode_lie_group<H: Hooks>(s: &Expr, fx: &Expr, x: &Expr, user: Option<&(Expr, Expr)>, split: Option<(Expr, Expr, Expr)>, hooks: &mut H) -> R<Option<Vec<Expr>>, H> {
    let df = diff(fx, x);
    let (y, h) = match split {
        Some((y, d, e)) => {
            let h = neg(&hk::<H, _>(hooks.simplify(&div(&d, &e)))?);
            (y, h)
        }
        None => {
            let y = dummy("y");
            let h = subs1(s, fx, &y);
            (y, h)
        }
    };
    if unmodelled(&h) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let mut inf = Vec::new();
    if let Some((xi, eta)) = user {
        inf.push(Inf::xe(xi.clone(), eta.clone()));
        // checkinfsol(Eq(df, s), inf, func=f(x), order=1) only prepends
        // 'user_defined' to the heuristic names; with a supplied
        // infinitesimal every name tries that same infinitesimal, so one
        // attempt decides them all. Its errors still propagate.
        checkinfsol(&sub(&df, s), fx, x, &[(xi.clone(), eta.clone())], hooks)?;
        return try_heuristic("user_defined", fx, x, &h, &y, &inf, hooks);
    }
    let mut sol = None;
    for name in HEURISTICS {
        sol = try_heuristic(name, fx, x, &h, &y, &inf, hooks)?;
        if sol.as_ref().is_some_and(|v| !v.is_empty()) {
            return Ok(sol);
        }
    }
    Ok(sol)
}

/// `LieGroup._get_general_solution`: the solution equations, or the
/// reference `NotImplementedError` when none was found.
pub fn lie_group_solve<H: Hooks>(eq: &Expr, fx: &Expr, x: &Expr, params: Params, hooks: &mut H) -> R<Vec<Expr>, H> {
    if unmodelled(eq) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let df = diff(fx, x);
    let r = hooks.solve(eq, &df);
    let eqsol = unless_nie(hooks, r)?.unwrap_or_default();
    let mut split = params.split;
    let mut desols = Vec::new();
    for s in &eqsol {
        if let Some(sol) = ode_lie_group(s, fx, x, params.user.as_ref(), split.take(), hooks)? {
            desols.extend(sol);
        }
    }
    if desols.is_empty() {
        return Err(LieError::Raise(
            Exc::NotImplemented,
            format!("The given ODE {} cannot be solved by the lie group method", crate::printer::to_str(eq)),
        ));
    }
    Ok(desols)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn remove_drops_arbitrary_functions() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let f = undefined("F", vec![sub(&x, &y)]);
        assert_eq!(lie_group_remove(&f), sub(&x, &y));
        assert_eq!(lie_group_remove(&add2(&x, &f)), x.clone());
        assert_eq!(lie_group_remove(&mul2(&y, &f)), mul2(&y, &sub(&x, &y)));
    }

    #[test]
    fn linear_ansatz_uses_alternate_dummies() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let a = homogeneous_ansatz("xi", 1, &x, &y);
        assert_eq!(free_symbols(&a).len(), 4);
    }

    #[test]
    fn function_sum_separates_like_the_reference_derivative() {
        // `x**(n*(m + 1) - m)*f' - a*f**n - b*x**(n*(m + 1))`: the second
        // derivative of `1/hinv` is factored as in the reference, so the
        // separation of `(1/factor).diff(y)` succeeds on the inverse ODE.
        let s = |n: &str| Expr::symbol(n);
        let (x, y, a, b, m, n, u) = (s("x"), s("y"), s("a"), s("b"), s("m"), s("n"), s("u"));
        let e = mul2(&n, &add2(&m, &Expr::one()));
        let h = div(&add2(&mul2(&a, &pow(&y, &n)), &mul2(&b, &pow(&x, &e))), &pow(&x, &sub(&e, &m)));
        let hinv = subs1(&subs(&inv(&h), &[(x.clone(), u.clone()), (y.clone(), x.clone())]), &u, &y);
        let xy = [x.clone(), y.clone()];
        let separate = |odefac: &Expr| {
            let factor = mul2(odefac, &diff2(&div(&Expr::one(), odefac), &x));
            separated(&diff(&div(&Expr::one(), &factor), &y), &xy)
        };
        assert!(separate(&h).is_none());
        let sep = separate(&hinv).expect("separable");
        assert_eq!(entry::<()>(&sep, Some(&x)).ok(), Some(pow(&x, &sub(&Expr::int(2), &n))));
    }
}
