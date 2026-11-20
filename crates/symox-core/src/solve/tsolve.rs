//! SymPy 1.14's `_tsolve` steps for one equation in one symbol, taken by
//! `solve` after the native solver declined (raised NotImplementedError),
//! and the helpers of `sympy.solvers.bivariate` that they use.
//!
//! * `_invert` (`integer_power=False`);
//! * the factor, `logcombine` and `f(x)**e = rhs` steps;
//! * the `bivariate` block: the `x**spow` change of variable,
//!   `_solve_lambert` (`_lambert`, `_linab`, `_mostfunc`,
//!   `_filtered_gens`) and `bivariate_type`;
//! * the closing `posify` retry;
//! * SymPy's `_solve` check of the `_tsolve` candidates
//!   (`checksol(f_num, sol) is not False`).
//!
//! The recursive `solve` calls, `nsimplify`, `equals` and the `inverse()`
//! of host functions go through a [`Host`], so that the public `solve`
//! entry (and its own fallbacks) is used exactly as by the reference.

use crate::expr::{Constant, Expr, Func, Kind};
use crate::ops::{add, div, mul, neg, pow, sqrt, sub};
use crate::poly::convert::PolyError;
use crate::subs::{subs1, xreplace};
use crate::traverse::{free_symbols, has};
use rustc_hash::FxHashMap;
use std::cell::Cell;

/// Which flags a recursive `solve` call receives.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flags {
    /// The flags of the current call.
    Same,
    /// No flags (`solve(eq, x)`).
    Plain,
    /// The flags of the current call with `bivariate=False`.
    NoBivariate,
    /// The flags of the current call once the bivariate block has popped
    /// `bivariate` (the posify retry: `bivariate` is back to its default).
    Popped,
}

/// Why a step stopped. `NotImplemented` and `Value` are the exceptions the
/// reference catches at some places; `None` means the error was raised
/// here (the binding builds the exception), `Some` carries a host error.
pub enum Fail<E> {
    NotImplemented(Option<E>, String),
    Value(Option<E>, String),
    Type(String),
    Poly(PolyError),
    Host(E),
}

impl<E> Fail<E> {
    fn ni(msg: &str) -> Self {
        Fail::NotImplemented(None, msg.to_string())
    }
    fn value(msg: &str) -> Self {
        Fail::Value(None, msg.to_string())
    }
    fn is_ni(&self) -> bool {
        matches!(self, Fail::NotImplemented(..))
    }
    /// The reference `except (NotImplementedError, ValueError)`; the
    /// polynomial errors are ValueError subclasses there.
    fn is_ni_or_value(&self) -> bool {
        matches!(self, Fail::NotImplemented(..) | Fail::Value(..) | Fail::Poly(_))
    }
}

impl<E> From<PolyError> for Fail<E> {
    fn from(e: PolyError) -> Self {
        Fail::Poly(e)
    }
}

type R<T, E> = Result<T, Fail<E>>;

/// The calls back into the public API.
pub trait Host {
    type Error;
    /// `solve(eq, x, **flags)` (a list of solutions).
    fn solve(&mut self, eq: &Expr, x: &Expr, flags: Flags) -> R<Vec<Expr>, Self::Error>;
    /// `nsimplify(e)`.
    fn nsimplify(&mut self, e: &Expr) -> R<Expr, Self::Error>;
    /// `e.equals(0)`.
    fn equals_zero(&mut self, e: &Expr) -> R<Option<bool>, Self::Error>;
    /// `e.count_ops()`.
    fn count_ops(&mut self, e: &Expr) -> usize;
    /// `f.inverse()(rhs)` for a one-argument function that the native
    /// table does not know (a host function class); None without one.
    fn inverse(&mut self, f: &Expr, rhs: &Expr) -> R<Option<Expr>, Self::Error>;
}

/// The flags of the call that this module reads.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub force: bool,
    pub bivariate: bool,
    pub check: bool,
    pub numerical: bool,
}

thread_local! {
    /// Depth of the inner `solve` calls of the posify retry, which the
    /// reference makes with `force=False`.
    static NO_FORCE: Cell<usize> = const { Cell::new(0) };
}

struct NoForce;
impl NoForce {
    fn enter() -> Self {
        NO_FORCE.with(|c| c.set(c.get() + 1));
        NoForce
    }
}
impl Drop for NoForce {
    fn drop(&mut self) {
        NO_FORCE.with(|c| c.set(c.get() - 1));
    }
}

fn no_force() -> bool {
    NO_FORCE.with(|c| c.get() > 0)
}

// ---------------------------------------------------------------- helpers

fn is_exp(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Exp, _) => true,
        Kind::Pow(b, _) => b.is_constant(Constant::E),
        _ => false,
    }
}

fn is_log(e: &Expr) -> bool {
    e.is_func(&Func::Log)
}

/// `e.is_Pow or isinstance(e, exp)`.
fn ispow(e: &Expr) -> bool {
    e.is_pow() || is_exp(e)
}

/// The reference `is_Function`.
fn is_function(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(f, _) if !matches!(f, Func::GenericExpr | Func::SymbolicDerivative | Func::Subs | Func::Limit | Func::Integral | Func::Sum | Func::Product | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _)) && !matches!(f, Func::Undefined(head) if head.is_internal("Tuple")))
}

fn func_args(e: &Expr) -> &[Expr] {
    match e.kind() {
        Kind::Func(_, a) => a,
        _ => &[],
    }
}

fn same_head(a: &Expr, b: &Expr) -> bool {
    matches!((a.kind(), b.kind()), (Kind::Func(f, _), Kind::Func(g, _)) if f == g)
}

/// Pre-order walk over the public arguments.
fn walk(e: &Expr, f: &mut dyn FnMut(&Expr)) {
    crate::traverse::preorder(e, f)
}

/// `e.atoms(cls)`, in first-seen order.
fn atoms(e: &Expr, pred: &dyn Fn(&Expr) -> bool) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    let mut seen = rustc_hash::FxHashSet::default();
    walk(e, &mut |n| {
        if pred(n) && seen.insert(n.clone()) {
            out.push(n.clone());
        }
    });
    out
}

/// `e.count(cls)`.
fn count(e: &Expr, pred: &dyn Fn(&Expr) -> bool) -> usize {
    let mut n = 0;
    walk(e, &mut |x| {
        if pred(x) {
            n += 1;
        }
    });
    n
}

fn has_pred(e: &Expr, pred: &dyn Fn(&Expr) -> bool) -> bool {
    count(e, pred) > 0
}

fn in_free(s: &Expr, e: &Expr) -> bool {
    free_symbols(e).contains(s)
}

/// `e.subs(dict)`: larger patterns first, like the binding.
fn subs_dict(e: &Expr, pairs: &[(Expr, Expr)]) -> Expr {
    let mut pairs = pairs.to_vec();
    pairs.sort_by_cached_key(|(old, _)| (std::cmp::Reverse(crate::traverse::node_count(old)), crate::sortkey::sort_key(old)));
    let mut r = e.clone();
    for (o, n) in &pairs {
        r = subs1(&r, o, n);
    }
    r
}

fn xrep(e: &Expr, old: &Expr, new: &Expr) -> Expr {
    let mut m = FxHashMap::default();
    m.insert(old.clone(), new.clone());
    xreplace(e, &m)
}

fn log(e: &Expr) -> Expr {
    crate::functions::func(Func::Log, vec![e.clone()])
}

fn exp(e: &Expr) -> Expr {
    crate::functions::func(Func::Exp, vec![e.clone()])
}

/// `expand_log(e, force=force)`.
fn expand_log(e: &Expr, force: bool) -> Expr {
    let mut v = e.clone();
    for _ in 0..16 {
        let n = crate::simp::expand_log_deep(&v, force, true);
        if n == v {
            break;
        }
        v = n;
    }
    v
}

fn expand_power_exp(e: &Expr) -> Expr {
    crate::simp::expand_power_exp_deep(e, true)
}

fn powsimp(e: &Expr) -> Expr {
    crate::simptools::powsimp::powsimp(e, &crate::simptools::powsimp::PowsimpOptions::default())
}

fn powdenest(e: &Expr) -> Expr {
    crate::simptools::powsimp::powdenest(e, false, false)
}

fn as_independent(e: &Expr, x: &Expr, as_add: Option<bool>) -> (Expr, Expr) {
    crate::decompose::as_independent(e, std::slice::from_ref(x), as_add)
}

fn is_zero(e: &Expr) -> Option<bool> {
    crate::assumptions::settle_nonreal(e, crate::assumptions::is_zero(e))
}

fn is_real(e: &Expr) -> Option<bool> {
    crate::assumptions::settle_nonreal(e, crate::assumptions::is_real(e))
}

fn make_args_add(e: &Expr) -> Vec<Expr> {
    if e.is_add() {
        e.args()
    } else {
        vec![e.clone()]
    }
}

/// `e.has(nan, zoo, oo, -oo)`.
pub fn has_illegal(e: &Expr) -> bool {
    has(e, &Expr::nan()) || has(e, &Expr::zoo()) || has(e, &Expr::oo()) || has(e, &Expr::neg_oo())
}

/// `_mexpand(e, recursive)`.
pub fn mexpand(e: &Expr, recursive: bool) -> Expr {
    let mut e = e.clone();
    loop {
        let n = crate::simptools::mexpand(&e);
        if n == e || !recursive {
            return n;
        }
        e = n;
    }
}

/// `_uniq`.
fn uniq(v: Vec<Expr>) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    for s in v {
        if !out.contains(&s) {
            out.push(s);
        }
    }
    out
}

/// The generators of `Poly(e)`, None when that fails.
fn poly_gens(e: &Expr) -> Option<Vec<Expr>> {
    let p = crate::poly::polyobj::Poly::from_expr(e, &[]).ok()?;
    Some(p.gens.clone())
}

// ---------------------------------------------------------------- _invert

const INVERSES: &[(&str, &str)] = &[
    ("acos", "cos"), ("acosh", "cosh"), ("acot", "cot"), ("acoth", "coth"),
    ("acsc", "csc"), ("acsch", "csch"), ("asec", "sec"), ("asech", "sech"),
    ("asin", "sin"), ("asinh", "sinh"), ("atan", "tan"), ("atanh", "tanh"),
    ("conjugate", "conjugate"), ("cot", "acot"), ("coth", "acoth"),
    ("erf", "erfinv"), ("erfc", "erfcinv"), ("erfcinv", "erfc"),
    ("erfinv", "erf"), ("exp", "log"), ("exp_polar", "log"), ("log", "exp"),
    ("sinh", "asinh"), ("tan", "atan"), ("tanh", "atanh"),
];

/// `lhs.inverse()(rhs)` for a one-argument function, None when it has none.
fn apply_inverse<H: Host>(h: &mut H, lhs: &Expr, rhs: &Expr) -> R<Option<Expr>, H::Error> {
    let Kind::Func(f, _) = lhs.kind() else { return Ok(None) };
    if let Func::Undefined(_) = f {
        return h.inverse(lhs, rhs);
    }
    let name = f.name();
    let Some((_, inv)) = INVERSES.iter().find(|(n, _)| *n == name) else { return Ok(None) };
    match Func::from_name(inv) {
        Some(g) => Ok(Some(crate::functions::func(g, vec![rhs.clone()]))),
        None => h.inverse(lhs, rhs),
    }
}

/// SymPy's `_invert(eq, sym)` for one symbol: `(rhs, lhs)`.
pub fn invert<H: Host>(h: &mut H, eq: &Expr, sym: &Expr) -> R<(Expr, Expr), H::Error> {
    // the reference returns `(eq, 0)` when `eq` is free of `sym`
    if !in_free(sym, eq) {
        return Ok((eq.clone(), Expr::zero()));
    }
    let mut lhs = eq.clone();
    let mut rhs = Expr::zero();
    loop {
        let was = lhs.clone();
        loop {
            let (indep, dep) = as_independent(&lhs, sym, None);
            if lhs.is_add() {
                if is_zero(&indep) == Some(true) {
                    break;
                }
                lhs = dep;
                rhs = sub(&rhs, &indep);
            } else {
                if indep.is_one() {
                    break;
                }
                lhs = dep;
                rhs = div(&rhs, &indep);
            }
        }
        if lhs.is_add() {
            let mut terms: Vec<(Expr, Vec<Expr>)> = Vec::new();
            for a in lhs.args() {
                let (i, d) = as_independent(&a, sym, None);
                match terms.iter_mut().find(|(k, _)| *k == d) {
                    Some((_, v)) => v.push(i),
                    None => terms.push((d, vec![i])),
                }
            }
            if terms.iter().any(|(_, v)| v.len() > 1) {
                let args: Vec<Expr> = terms
                    .iter()
                    .map(|(d, i)| if i.len() > 1 { mul(&[add(i), d.clone()]) } else { mul(&[i[0].clone(), d.clone()]) })
                    .collect();
                lhs = add(&args);
            }
        }
        if lhs.is_add()
            && rhs.is_zero()
            && lhs.args().len() == 2
            && crate::inspection::classify(&lhs, std::slice::from_ref(sym), crate::inspection::PolynomialKind::Polynomial) != Some(true)
        {
            let ab = crate::simptools::ordered(lhs.args());
            let (a, b) = (&ab[0], &ab[1]);
            let (ai, ad) = as_independent(a, sym, None);
            let (bi, bd) = as_independent(b, sym, None);
            if ispow(&ad) || ispow(&bd) {
                let (a_base, a_exp) = crate::decompose::as_base_exp(&ad);
                let (b_base, b_exp) = crate::decompose::as_base_exp(&bd);
                if a_base == b_base && crate::simptools::extract_additively(&a_exp, &b_exp).is_none() {
                    lhs = powsimp(&powdenest(&div(&ad, &bd)));
                    rhs = neg(&div(&bi, &ai));
                } else {
                    let rat = div(&ad, &bd);
                    let l = powsimp(&rat);
                    if l != rat {
                        lhs = l;
                        rhs = neg(&div(&bi, &ai));
                    }
                }
            } else if ai == neg(&bi) && is_function(&ad) && same_head(&ad, &bd) {
                let (na, nb) = (func_args(&ad).len(), func_args(&bd).len());
                if na == 1 && nb == 1 {
                    lhs = sub(&func_args(&ad)[0], &func_args(&bd)[0]);
                } else if na == nb {
                    return Err(Fail::ni("equal function with more than 1 argument"));
                } else {
                    return Err(Fail::value("function with different numbers of args"));
                }
            }
        } else if lhs.is_mul() && lhs.args().iter().any(ispow) {
            lhs = powsimp(&powdenest(&lhs));
        }
        if is_function(&lhs) {
            let inverse = if func_args(&lhs).len() == 1 { apply_inverse(h, &lhs, &rhs)? } else { None };
            if let Some(r) = inverse {
                rhs = r;
                lhs = func_args(&lhs)[0].clone();
            } else if lhs.is_func(&Func::Atan2) {
                let (y, x) = (func_args(&lhs)[0].clone(), func_args(&lhs)[1].clone());
                let s = sqrt(&add(&[pow(&x, &Expr::int(2)), pow(&y, &Expr::int(2))]));
                lhs = mul(&[Expr::int(2), crate::functions::func(Func::Atan, vec![div(&y, &add(&[s, x]))])]);
            } else if same_head(&lhs, &rhs) {
                let (nl, nr) = (func_args(&lhs).len(), func_args(&rhs).len());
                if nl == 1 && nr == 1 {
                    let (l, r) = (func_args(&lhs)[0].clone(), func_args(&rhs)[0].clone());
                    lhs = l;
                    rhs = r;
                } else if nl == nr {
                    return Err(Fail::ni("equal function with more than 1 argument"));
                } else {
                    return Err(Fail::value("function with different numbers of args"));
                }
            }
        }
        if !rhs.is_zero() {
            if let Kind::Pow(_, x) = lhs.kind() {
                if x.is_integer_number() && x.as_number().is_some_and(|n| n.is_negative()) {
                    lhs = div(&Expr::one(), &lhs);
                    rhs = div(&Expr::one(), &rhs);
                }
            }
        }
        if lhs == was {
            break;
        }
    }
    Ok((rhs, lhs))
}

// ---------------------------------------------------------------- bivariate.py

/// `_filtered_gens`: the generators with `symbol`; of `g` and `1/g` only
/// one is kept.
pub fn filtered_gens(gens: &[Expr], symbol: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    for g in gens {
        if in_free(symbol, g) && !out.contains(g) {
            out.push(g.clone());
        }
    }
    for g in out.clone() {
        let ag = div(&Expr::one(), &g);
        if out.contains(&g) && out.contains(&ag) {
            let g = if !crate::simptools::as_numer_denom(&ag).1.is_one() { ag } else { g };
            out.retain(|x| *x != g);
        }
    }
    out
}

/// The kind of terms `_mostfunc` looks for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FuncClass {
    Log,
    Exp,
    Pow,
}

fn class_pred(c: FuncClass) -> impl Fn(&Expr) -> bool {
    move |e: &Expr| match c {
        FuncClass::Log => is_log(e),
        FuncClass::Exp => e.is_func(&Func::Exp),
        FuncClass::Pow => e.is_pow(),
    }
}

/// `_mostfunc(lhs, func, X)`.
pub fn mostfunc(lhs: &Expr, func: FuncClass, x: Option<&Expr>) -> Option<Expr> {
    let pred = class_pred(func);
    let fterms: Vec<Expr> = atoms(lhs, &pred)
        .into_iter()
        .filter(|t| match x {
            None => true,
            Some(x) if x.is_symbol() => in_free(x, t),
            Some(x) => has(t, x),
        })
        .collect();
    match fterms.len() {
        0 => None,
        1 => Some(fterms[0].clone()),
        _ => {
            let mut best: Option<(usize, Expr)> = None;
            for t in crate::simptools::ordered(fterms) {
                let n = count(&t, &pred);
                if best.as_ref().is_none_or(|(m, _)| n > *m) {
                    best = Some((n, t));
                }
            }
            best.map(|(_, t)| t)
        }
    }
}

/// `_linab(arg, symbol)`: `(a, b, x)` with `arg = a*x + b`.
pub fn linab(arg: &Expr, symbol: &Expr) -> (Expr, Expr, Expr) {
    let arg = crate::simptools::exprtools::factor_terms(&crate::expand::expand(arg), false, false, false, true);
    let (ind, dep) = as_independent(&arg, symbol, None);
    if arg.is_mul() && dep.is_add() {
        let (a, b, x) = linab(&dep, symbol);
        return (mul(&[ind.clone(), a]), mul(&[ind, b]), x);
    }
    let (mut a, b, mut x) = if !arg.is_add() {
        (ind, Expr::zero(), dep)
    } else {
        let (a, x) = as_independent(&crate::simptools::misc::separatevars(&dep, false), symbol, Some(false));
        (a, ind, x)
    };
    if crate::decompose::could_extract_minus_sign(&x) {
        a = neg(&a);
        x = neg(&x);
    }
    (a, b, x)
}

/// `_lambert(eq, x)`: the solutions of `a*log(b*X + c) + d*X + f = 0`.
pub fn lambert<H: Host>(h: &mut H, eq: &Expr, x: &Expr) -> R<Vec<Expr>, H::Error> {
    let mut eq = mexpand(&expand_log(eq, false), false);
    let Some(mut mainlog) = mostfunc(&eq, FuncClass::Log, Some(x)) else { return Ok(vec![]) };
    let mut other = subs1(&eq, &mainlog, &Expr::zero());
    let mother = neg(&other);
    if is_log(&mother) {
        let inner = func_args(&mainlog)[0].clone();
        eq = subs1(&sub(&eq, &other), &mainlog, &inner);
        mainlog = inner;
        if !is_log(&mainlog) {
            return Ok(vec![]);
        }
        other = neg(&func_args(&mother)[0]);
        eq = add(&[eq, other.clone()]);
    }
    if !in_free(x, &other) {
        return Ok(vec![]);
    }
    let (d, f, x2) = linab(&other, x);
    let logterm = crate::poly::api::collect(&sub(&eq, &other), std::slice::from_ref(&mainlog));
    let a = match crate::simptools::as_coefficient(&logterm, &mainlog) {
        Some(a) if !in_free(x, &a) => a,
        _ => return Ok(vec![]),
    };
    let logarg = func_args(&mainlog)[0].clone();
    let (b, c, x1) = linab(&logarg, x);
    if x1 != x2 {
        return Ok(vec![]);
    }
    let u = Expr::dummy("rhs", crate::Assumptions::NONE);
    let xusolns = h.solve(&sub(&x1, &u), x, Flags::Plain)?;
    let ratio = div(&div(&sub(&mul(&[c.clone(), d.clone()]), &mul(&[b.clone(), f])), &a), &b);
    let (num, den) = crate::simptools::as_numer_denom(&ratio);
    let (p, den) = crate::simptools::as_coeff_mul(&den, false);
    let e = exp(&div(&num, &den));
    let t = Expr::dummy("t", crate::Assumptions::NONE);
    let rts = crate::poly::roots::roots(&sub(&pow(&t, &p), &e), std::slice::from_ref(&t))?;
    let mut sol = Vec::new();
    let scale = div(&d, &mul(&[a.clone(), b.clone()]));
    for (r, _) in rts {
        let arg = mul(&[scale.clone(), r]);
        for k in [-1i64, 0] {
            let w = if k == 0 { crate::special::lambertw(&arg) } else { crate::special::lambertw_k(&arg, &Expr::int(k)) };
            if k != 0 && is_real(&w) != Some(true) {
                continue;
            }
            let rhs = add(&[neg(&div(&c, &b)), mul(&[div(&a, &d), w])]);
            sol.extend(xusolns.iter().map(|xu| subs1(xu, &u, &rhs)));
        }
    }
    Ok(sol)
}

/// `_solve_lambert(f, symbol, gens)`: NotImplemented without a Lambert form.
pub fn solve_lambert<H: Host>(h: &mut H, f: &Expr, symbol: &Expr, gens: &[Expr]) -> R<Vec<Expr>, H::Error> {
    let (nrhs, lhs0) = as_independent(f, symbol, Some(true));
    let mut rhs = neg(&nrhs);
    let mut lhs = lhs0;
    let lamcheck = gens.iter().any(|g| is_exp(g) || is_log(g) || matches!(g.kind(), Kind::Pow(_, x) if in_free(symbol, x)));
    if !lamcheck {
        return Err(Fail::ni(""));
    }
    if lhs.is_add() || lhs.is_mul() {
        let t = match symbol.as_sym() {
            Some(s) => Expr::from_kind(Kind::Sym(s.fresh_dummy_like("t"))),
            None => Expr::dummy("t", crate::Assumptions::NONE),
        };
        let even = |i: &Expr| -> Option<Expr> {
            match i.kind() {
                Kind::Pow(b, x) if b == symbol && crate::assumptions::fact(x, "even") == Some(true) => Some(pow(&t, x)),
                _ => None,
            }
        };
        lhs = crate::simptools::bottom_up(&lhs, &mut |i| even(i).unwrap_or_else(|| i.clone()));
        let even_degree = |h: &mut H, expr: &Expr| -> R<Vec<Expr>, H::Error> {
            let nl = xrep(expr, &t, &neg(symbol));
            let pl = xrep(expr, &t, symbol);
            let mut sols = solve_lambert(h, &nl, symbol, gens)?;
            if pl != nl {
                sols.extend(solve_lambert(h, &pl, symbol, gens)?);
            }
            Ok(uniq(sols))
        };
        if lhs.is_add() && has(&lhs, &t) {
            let t_indep = subs1(&lhs, &t, &Expr::zero());
            let t_term = sub(&lhs, &t_indep);
            let rhs2 = sub(&rhs, &t_indep);
            if !t_term.is_add() && !rhs2.is_zero() && !(has(&t_term, &Expr::zoo()) || has(&t_term, &Expr::nan())) {
                let eq = expand_log(&sub(&log(&t_term), &log(&rhs2)), false);
                return even_degree(h, &eq);
            }
        } else if lhs.is_mul() && !rhs.is_zero() {
            lhs = expand_log(&log(&lhs), true);
            rhs = log(&rhs);
            if has(&lhs, &t) && lhs.is_add() {
                let eq = sub(&lhs, &rhs);
                return even_degree(h, &eq);
            }
        }
        lhs = xrep(&lhs, &t, symbol);
    }
    lhs = powsimp(&crate::poly::api::factor_deep(&lhs)?);
    let r = Expr::dummy("r", crate::Assumptions::NONE);
    let (i, l) = invert(h, &sub(&lhs, &r), symbol)?;
    lhs = l;
    rhs = xrep(&i, &r, &rhs);

    let mut soln: Vec<Expr> = Vec::new();
    if let Some(mainlog) = mostfunc(&lhs, FuncClass::Log, Some(symbol)) {
        if lhs.is_mul() && !rhs.is_zero() {
            soln = lambert(h, &sub(&log(&lhs), &log(&rhs)), symbol)?;
        } else if lhs.is_add() {
            let other = subs1(&lhs, &mainlog, &Expr::zero());
            if !other.is_zero() && !other.is_add() && atoms(&other, &|e: &Expr| e.is_pow()).iter().any(|p| in_free(symbol, p)) {
                let diff = if rhs.is_zero() {
                    sub(&log(&other), &log(&sub(&other, &lhs)))
                } else {
                    sub(&log(&sub(&lhs, &other)), &log(&sub(&rhs, &other)))
                };
                soln = lambert(h, &expand_log(&diff, false), symbol)?;
            } else {
                soln = lambert(h, &sub(&lhs, &rhs), symbol)?;
            }
        }
    }
    for class in [FuncClass::Exp, FuncClass::Pow] {
        if !soln.is_empty() {
            break;
        }
        let Some(main) = mostfunc(&lhs, class, Some(symbol)) else { continue };
        if class == FuncClass::Pow {
            let Kind::Pow(_, x) = main.kind() else { continue };
            if !in_free(symbol, x) {
                continue;
            }
        }
        lhs = crate::poly::api::collect(&lhs, std::slice::from_ref(&main));
        if lhs.is_mul() && !rhs.is_zero() {
            soln = lambert(h, &expand_log(&sub(&log(&lhs), &log(&rhs)), false), symbol)?;
        } else if lhs.is_add() {
            let other = subs1(&lhs, &main, &Expr::zero());
            let mut mainterm = sub(&lhs, &other);
            rhs = sub(&rhs, &other);
            if class == FuncClass::Exp
                && crate::decompose::could_extract_minus_sign(&mainterm)
                && crate::decompose::could_extract_minus_sign(&rhs)
            {
                mainterm = neg(&mainterm);
                rhs = neg(&rhs);
            }
            let diff = sub(&log(&mainterm), &log(&rhs));
            soln = lambert(h, &expand_log(&diff, false), symbol)?;
        }
    }
    if soln.is_empty() {
        return Err(Fail::NotImplemented(
            None,
            format!("{} does not appear to have a solution in terms of LambertW", crate::printer::to_str(f)),
        ));
    }
    Ok(crate::simptools::ordered(soln))
}

/// `bivariate_type(f, x, y)`: `(g, p, u)` with `f(x, y) = p(u)`, `u = g(x, y)`.
pub fn bivariate_type(f: &Expr, x: &Expr, y: &Expr) -> Result<Option<(Expr, Expr, Expr)>, PolyError> {
    use crate::poly::polyobj::Poly;
    let f = Poly::from_expr(f, &[x.clone(), y.clone()])?.as_expr();
    let dx = Expr::dummy("_x", crate::Assumptions::NONE);
    let dy = Expr::dummy("_y", crate::Assumptions::NONE);
    let fd = subs_dict(&f, &[(x.clone(), dx.clone()), (y.clone(), dy.clone())]);
    let p = Poly::from_expr(&fd, &[dx.clone(), dy.clone()])?;
    let Some((g, q, u)) = bivariate_inner(&p, &dx, &dy) else { return Ok(None) };
    let mut m = FxHashMap::default();
    m.insert(dx, x.clone());
    m.insert(dy, y.clone());
    Ok(Some((xreplace(&g, &m), xreplace(&q, &m), u)))
}

fn bivariate_inner(p: &crate::poly::polyobj::Poly, x0: &Expr, y0: &Expr) -> Option<(Expr, Expr, Expr)> {
    let u = Expr::dummy("u", crate::symbol::Assumptions::POSITIVE);
    let f = p.as_expr();
    let free_xy = |e: &Expr, x: &Expr, y: &Expr| {
        let fs = free_symbols(e);
        fs.contains(x) || fs.contains(y)
    };
    // f(x*y)
    let mut new = Vec::new();
    let mut all = true;
    for a in make_args_add(&f) {
        let a = mexpand(&subs1(&a, x0, &div(&u, y0)), false);
        if free_xy(&a, x0, y0) {
            all = false;
            break;
        }
        new.push(a);
    }
    if all {
        return Some((mul(&[x0.clone(), y0.clone()]), add(&new), u));
    }
    let ok = |v: &Expr, c: &Expr| -> Option<Expr> {
        let n = mexpand(&subs1(&f, v, c), false);
        (!free_xy(&n, x0, y0)).then_some(n)
    };
    let (ix, iy) = (p.gen_index(x0)?, p.gen_index(y0)?);
    let ngens = p.true_gens().len();
    let mono = |pairs: &[(usize, u32)]| {
        let mut m = vec![0u32; ngens];
        for (i, d) in pairs {
            m[*i] = *d;
        }
        m
    };
    let root = |c: Expr, d: u32| pow(&c, &div(&Expr::one(), &Expr::int(d as i64)));
    // f(a*x + b*y)
    let d = p.degree(ix)?;
    if p.degree(iy) == Some(d) {
        let a = root(p.coeff_monomial(&mono(&[(ix, d)])), d);
        let b = root(p.coeff_monomial(&mono(&[(iy, d)])), d);
        if let Some(n) = ok(x0, &div(&sub(&u, &mul(&[b.clone(), y0.clone()])), &a)) {
            return Some((add(&[mul(&[a, x0.clone()]), mul(&[b, y0.clone()])]), n, u));
        }
    }
    // f(a*x*y + b*y)
    if p.degree(iy) == Some(d) {
        let (mut x, mut y, mut jx, mut jy) = (x0.clone(), y0.clone(), ix, iy);
        for _ in 0..2 {
            let a = root(p.coeff_monomial(&mono(&[(jx, d), (jy, d)])), d);
            let b = root(p.coeff_monomial(&mono(&[(jy, d)])), d);
            if let Some(n) = ok(&x, &div(&div(&sub(&u, &mul(&[b.clone(), y.clone()])), &a), &y)) {
                return Some((add(&[mul(&[a, x.clone(), y.clone()]), mul(&[b, y])]), n, u));
            }
            std::mem::swap(&mut x, &mut y);
            std::mem::swap(&mut jx, &mut jy);
        }
    }
    None
}

// ---------------------------------------------------------------- checksol

/// SymPy's `checksol(f, {sym: sol})` for one expression (`minimal=False`,
/// `simplify=False` as after `_solve` simplified its results); see
/// `solve::checksol`.
pub fn checksol(f: &Expr, sym: &Expr, sol: &Expr, numerical: bool, force: bool) -> Option<bool> {
    crate::solve::checksol::checksol(f, sym, sol, crate::solve::checksol::CheckFlags { numerical, force, simplify: false, minimal: false })
}

// ---------------------------------------------------------------- _tsolve

/// `_bivariate_eq`: the reference rewrite of `eq` before `_solve_lambert`.
fn bivariate_eq<H>(eq: &Expr, sym: &Expr) -> R<Expr, H> {
    let Some(gens) = poly_gens(eq) else { return Err(Fail::ni("no polynomial form")) };
    let g = filtered_gens(&gens, sym);
    let mut up_or_log: Vec<Expr> = Vec::new();
    for gi in g {
        if is_exp(&gi) || is_log(&gi) {
            up_or_log.push(gi);
        } else if gi.is_pow() {
            let s = powdenest(&expand_power_exp(&gi));
            if matches!(s.kind(), Kind::Pow(_, x) if in_free(sym, x)) {
                up_or_log.push(gi);
            }
        }
    }
    let pairs: Vec<(Expr, Expr)> = up_or_log.into_iter().map(|g| (g, Expr::zero())).collect();
    let eq_down = subs_dict(&expand_log(&expand_power_exp(eq), false), &pairs);
    let factored = crate::poly::api::factor_deep(&eq_down)?;
    Ok(expand_power_exp(&add(&[factored, sub(eq, &eq_down)])))
}

fn numeric_min(v: &[Expr]) -> R<Expr, ()> {
    let mut best = v[0].clone();
    for x in &v[1..] {
        match crate::simptools::num_sign(&sub(x, &best)) {
            Some(std::cmp::Ordering::Less) => best = x.clone(),
            Some(_) => {}
            None => return Err(Fail::Type("cannot determine truth value of Relational".into())),
        }
    }
    Ok(best)
}

/// The `bivariate` block of `_tsolve`. Ok(None) when it gives nothing.
fn bivariate_step<H: Host>(h: &mut H, eq: &Expr, sym: &Expr) -> R<Option<Vec<Expr>>, H::Error> {
    let mut exps: Vec<Expr> = Vec::new();
    for j in atoms(eq, &is_log) {
        for i in atoms(&j, &|e: &Expr| e.is_pow()) {
            if let Kind::Pow(b, x) = i.kind() {
                if b == sym && !exps.contains(x) {
                    exps.push(x.clone());
                }
            }
        }
    }
    let spow = match exps.len() {
        0 => Expr::one(),
        1 => exps[0].clone(),
        _ => numeric_min(&exps).map_err(|e| match e {
            Fail::Type(m) => Fail::Type(m),
            _ => Fail::ni(""),
        })?,
    };
    if !spow.is_one() {
        let p = pow(sym, &spow);
        let u = Expr::dummy("bivariate-cov", crate::Assumptions::NONE);
        let ueq = subs1(eq, &p, &u);
        if !in_free(sym, &ueq) {
            let sol = match h.solve(&ueq, &u, Flags::Same) {
                Ok(s) => s,
                Err(e) if e.is_ni() => return Ok(None),
                Err(e) => return Err(e),
            };
            let inv = match h.solve(&sub(&p, &u), sym, Flags::Plain) {
                Ok(s) => s,
                Err(e) if e.is_ni() => return Ok(None),
                Err(e) => return Err(e),
            };
            let mut out = Vec::new();
            for i in &inv {
                for s in &sol {
                    out.push(subs1(i, &u, s));
                }
            }
            return Ok(Some(out));
        }
    }
    let (rhs, lhs) = match bivariate_eq(eq, sym).and_then(|e| invert(h, &e, sym)) {
        Ok(v) => v,
        Err(e) if e.is_ni_or_value() => return Ok(None),
        Err(e) => return Err(e),
    };
    if !has(&lhs, sym) {
        return Ok(None);
    }
    let Some(gens) = poly_gens(&lhs) else { return Ok(None) };
    let g = filtered_gens(&gens, sym);
    let eq2 = sub(&lhs, &rhs);
    let lambert = (|| -> R<Vec<Expr>, H::Error> {
        let mut sols = solve_lambert(h, &eq2, sym, &g)?;
        for s in sols.iter_mut() {
            let ns = h.nsimplify(s)?;
            if ns != *s && h.count_ops(&ns) <= h.count_ops(s) {
                let ok = match crate::solve::checksol::checksol(&eq2, sym, &ns, Default::default()) {
                    Some(v) => Some(v),
                    None => h.equals_zero(&subs1(&eq2, sym, &ns))?,
                };
                if ok == Some(true) {
                    *s = ns;
                }
            }
        }
        Ok(sols)
    })();
    match lambert {
        Ok(s) => return Ok(Some(s)),
        Err(e) if e.is_ni() => {}
        Err(e) => return Err(e),
    }
    if g.len() != 2 {
        return Ok(None);
    }
    let (gg, p, u) = match bivariate_type(&eq2, &g[0], &g[1]) {
        Ok(Some(v)) => v,
        _ => return Ok(None),
    };
    let inversion = match h.solve(&sub(&gg, &u), sym, Flags::NoBivariate) {
        Ok(s) => s,
        Err(e) if e.is_ni() => return Ok(None),
        Err(e) => return Err(e),
    };
    if inversion.is_empty() {
        return Ok(None);
    }
    let sol = match h.solve(&p, &u, Flags::Same) {
        Ok(s) => s,
        Err(e) if e.is_ni() => return Ok(None),
        Err(e) => return Err(e),
    };
    let mut out = Vec::new();
    for i in &inversion {
        for s in &sol {
            out.push(subs1(i, &u, s));
        }
    }
    Ok(Some(uniq(out)))
}

/// The `_tsolve` steps for `f = 0` in `sym` (a Symbol). Ok(None) when no
/// step applies (the caller then re-raises NotImplementedError).
pub fn tsolve_step<H: Host>(h: &mut H, f: &Expr, sym: &Expr, opts: Options) -> R<Option<Vec<Expr>>, H::Error> {
    let (num, den) = crate::simptools::as_numer_denom(f);
    let f = if has(&den, sym) { num.clone() } else { f.clone() };
    let result = steps(h, &f, sym, opts)?;
    let Some(result) = result else { return Ok(None) };
    let result = uniq(result);
    if !opts.check {
        return Ok(Some(result));
    }
    // the reference `_solve` keeps the candidates whose check is not False
    let force = opts.force && !no_force();
    Ok(Some(result.into_iter().filter(|r| checksol(&num, sym, r, opts.numerical, force) != Some(false)).collect()))
}

/// The steps of `_tsolve` between `_invert` and the bivariate block: the
/// `Add` (factor, logcombine), `Pow`, `sin`/`cos`/`LambertW` and
/// `rewrite(exp)` branches. Ok(None) when none of them returned.
fn inverted_steps<H: Host>(h: &mut H, f: &Expr, sym: &Expr, lhs: &Expr, rhs: &Expr, opts: Options) -> R<Option<Vec<Expr>>, H::Error> {
    if lhs.is_add() {
        // factoring (after powdenest) may split the equation
        let g = crate::poly::api::factor(&powdenest(&sub(lhs, rhs)))?;
        if g.is_mul() {
            let (indep, dep) = as_independent(&g, sym, Some(false));
            if !indep.is_one() && !indep.is_zero() && has(&dep, sym) {
                return h.solve(&dep, sym, Flags::Same).map(Some);
            }
            if dep.is_mul() {
                return h.solve(&g, sym, Flags::Same).map(Some);
            }
        }
        if !rhs.is_zero() && has_pred(lhs, &is_log) {
            let g = crate::simp::logcombine(lhs, opts.force && !no_force());
            if count(&g, &is_log) != count(lhs, &is_log) {
                if is_log(&g) {
                    return h.solve(&sub(&func_args(&g)[0], &exp(rhs)), sym, Flags::Same).map(Some);
                }
                return h.solve(&sub(&g, rhs), sym, Flags::Same).map(Some);
            }
        }
    } else if let (Kind::Pow(b, e), false) = (lhs.kind(), is_exp(lhs)) {
        return pow_step(h, f, sym, lhs, rhs, b, e).map(Some);
    } else if let Kind::Func(func, args) = lhs.kind() {
        if args.len() == 1 {
            // `multi_inverses`
            let inverses = match func {
                Func::Sin => {
                    let a = crate::functions::func(Func::Asin, vec![rhs.clone()]);
                    Some(vec![a.clone(), sub(&Expr::pi(), &a)])
                }
                Func::Cos => {
                    let a = crate::functions::func(Func::Acos, vec![rhs.clone()]);
                    Some(vec![a.clone(), sub(&mul(&[Expr::int(2), Expr::pi()]), &a)])
                }
                _ => None,
            };
            if let Some(inverses) = inverses {
                let mut out = Vec::new();
                for i in inverses {
                    out.extend(h.solve(&sub(&args[0], &i), sym, Flags::Same)?);
                }
                return Ok(Some(uniq(out)));
            }
            if func.name() == "LambertW" {
                return h.solve(&sub(&args[0], &mul(&[rhs.clone(), exp(rhs)])), sym, Flags::Same).map(Some);
            }
        }
    }
    // `rebuild(lhs.rewrite(exp))`: the rewrite keeps unevaluated nodes
    // (`exp(log(x)/2)`), which the reference re-evaluates bottom-up.
    let rewritten = rl_rebuild(&crate::special::rewrite::rewrite(lhs, &[], "exp", true));
    if rewritten != *lhs {
        return h.solve(&sub(&rewritten, rhs), sym, Flags::Same).map(Some);
    }
    Ok(None)
}

/// `sympy.strategies.rl.rebuild`: every node rebuilt with the evaluating
/// constructors, bottom-up.
fn rl_rebuild(e: &Expr) -> Expr {
    use crate::expr::{Func, Kind};
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Rel(..) => {
            let new: Vec<Expr> = crate::simptools::args(e).iter().map(rl_rebuild).collect();
            crate::simptools::rebuild(e, &new)
        }
        Kind::Func(f, _) if !matches!(f, Func::SymbolicDerivative | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(..) | Func::Subs | Func::Integral | Func::Sum | Func::Limit) => {
            let new: Vec<Expr> = crate::simptools::args(e).iter().map(rl_rebuild).collect();
            crate::simptools::rebuild(e, &new)
        }
        _ => e.clone(),
    }
}

/// `checksol(f, {sym: s})`, then `f.subs(sym, s).equals(0)` when undecided.
fn checked<H: Host>(h: &mut H, f: &Expr, sym: &Expr, s: &Expr) -> R<bool, H::Error> {
    match crate::solve::checksol::checksol(f, sym, s, Default::default()) {
        Some(v) => Ok(v),
        None => Ok(h.equals_zero(&subs1(f, sym, s))? == Some(true)),
    }
}

/// The `lhs.is_Pow` branch of `_tsolve` for `lhs = b**e` and `rhs`.
fn pow_step<H: Host>(h: &mut H, f: &Expr, sym: &Expr, lhs: &Expr, rhs: &Expr, b: &Expr, e: &Expr) -> R<Vec<Expr>, H::Error> {
    if e.is_integer_number() && sub(lhs, rhs) != *f {
        return h.solve(&sub(lhs, rhs), sym, Flags::Same);
    }
    if !in_free(sym, e) {
        return h.solve(&sub(b, &pow(rhs, &div(&Expr::one(), e))), sym, Flags::Same);
    }
    // `_tsolve` calls this with a Dummy before passing the actual number
    if free_symbols(rhs).iter().any(|s| s.as_sym().is_some_and(|s| s.is_dummy())) {
        return Err(Fail::ni(""));
    }
    if rhs.is_zero() {
        // f(x)**g(x) = 0 where f(x) = 0 and g(x) != 0
        let sol_base = h.solve(b, sym, Flags::Same)?;
        return Ok(sol_base.into_iter().filter(|s| !subs1(e, sym, s).is_zero()).collect());
    }
    if !has(b, sym) {
        if b.is_zero() {
            return h.solve(e, sym, Flags::Same);
        }
        let (rb, re) = crate::ops::as_base_exp(rhs);
        let mut sol = if *b == rb {
            h.solve(&sub(e, &re), sym, Flags::Same)?
        } else {
            let g = sub(&exp(&mul(&[log(b), e.clone()])), &exp(&log(rhs)));
            h.solve(&g, sym, Flags::Same)?
        };
        // guess a rational exponent
        let ratio = div(&log(&crate::functions::func(Func::Abs, vec![rhs.clone()])), &log(&crate::functions::func(Func::Abs, vec![b.clone()])));
        let e_rat = h.nsimplify(&ratio)?;
        let e_rat = crate::simplify::simplify(&crate::simptools::misc::posify(&e_rat).0);
        let (n, d) = crate::simptools::radsimp::fraction(&e_rat, false);
        if crate::expand::expand(&sub(&pow(b, &n), &pow(rhs, &d))).is_zero() {
            let mut kept = Vec::new();
            for s in sol {
                if !equal(h, &subs1(e, sym, &s), &e_rat)? {
                    kept.push(s);
                }
            }
            sol = kept;
            sol.extend(h.solve(&sub(e, &e_rat), sym, Flags::Same)?);
        }
        return Ok(uniq(sol));
    }
    // f(x)**g(x) = c
    let mut sol = Vec::new();
    let logform = sub(&mul(&[e.clone(), log(b)]), &log(rhs));
    // Whether the general (logarithmic) form went unsolved: the guesses
    // below are then not a complete answer.
    let mut logform_unsolved = false;
    if logform != sub(lhs, rhs) {
        match h.solve(&logform, sym, Flags::Same) {
            Ok(v) => sol.extend(v),
            Err(err) if err.is_ni() => logform_unsolved = true,
            Err(err) => return Err(err),
        }
    }
    let mut check: Vec<Expr> = Vec::new();
    if rhs.is_one() {
        check.extend(h.solve(e, sym, Flags::Same)?);
        check.extend(h.solve(&sub(b, &Expr::one()), sym, Flags::Same)?);
        check.extend(h.solve(&add(&[b.clone(), Expr::one()]), sym, Flags::Same)?);
    } else if let Some(r) = rhs.as_number().filter(|n| n.is_rational()) {
        use num_traits::Signed;
        let Some((p, q)) = r.numer_denom() else { return Ok(uniq(sol)) };
        let one = num_bigint::BigInt::from(1);
        for dv in crate::ntheory::divisors(&p.abs()) {
            if dv == one {
                continue;
            }
            let (ex, exact) = crate::ntheory::misc::integer_log(&p, &dv).map_err(|_| Fail::value("integer_log"))?;
            if !exact {
                continue;
            }
            for sv in crate::ntheory::divisors(&q.abs()) {
                if num_traits::pow::pow(sv.clone(), ex as usize) == q {
                    let Some(r) = crate::number::Number::rational(dv.clone(), sv.clone()) else { continue };
                    let r = Expr::num(r);
                    check.extend(h.solve(&sub(b, &r), sym, Flags::Same)?);
                    check.extend(h.solve(&add(&[b.clone(), r]), sym, Flags::Same)?);
                    check.extend(h.solve(&sub(e, &Expr::int(ex as i64)), sym, Flags::Same)?);
                }
            }
        }
    } else if crate::assumptions::is_real(rhs) == Some(true) && crate::assumptions::is_rational(rhs) == Some(false) {
        let (b_l, e_l) = crate::ops::as_base_exp(b);
        let (_, d) = crate::simplify::numer_denom(&mul(&[e_l.clone(), e.clone()]));
        let n = crate::simplify::numer_denom(&mul(&[e_l.clone(), e.clone()])).0;
        let (rb, re) = crate::ops::as_base_exp(&crate::simptools::sqrtdenest::sqrtdenest(rhs, 3));
        for i in h.solve(&sub(b, &rb), sym, Flags::Same)? {
            check.push(crate::simptools::sqrtdenest::sqrtdenest(&i, 3));
        }
        for i in h.solve(&sub(e, &re), sym, Flags::Same)? {
            check.push(crate::simptools::sqrtdenest::sqrtdenest(&i, 3));
        }
        if !mul(&[e_l.clone(), d.clone()]).is_one() {
            check.extend(h.solve(&sub(&pow(&b_l, &n), &pow(rhs, &mul(&[e_l, d]))), sym, Flags::Same)?);
        }
    }
    for s in check {
        if checked(h, f, sym, &s)? {
            sol.push(s);
        }
    }
    if logform_unsolved && sol.is_empty() {
        // The reference returns `[]` once its integer guesses fail
        // (`solve(x**x - 2, x, bivariate=False, force=False)`), although
        // `exp(LambertW(log(2)))` is a root: an empty list here is not a
        // proof that there is no solution.
        return Err(Fail::ni("f(x)**g(x) = c: the logarithmic form is not solved"));
    }
    Ok(uniq(sol))
}

/// The `equal` helper of the reference `_tsolve`.
fn equal<H: Host>(h: &mut H, a: &Expr, b: &Expr) -> R<bool, H::Error> {
    let d = Expr::dummy("equal", crate::Assumptions::NONE);
    match crate::solve::checksol::checksol(&sub(a, &d), &d, b, Default::default()) {
        Some(v) => Ok(v),
        None => {
            if h.nsimplify(a)? != h.nsimplify(b)? {
                return Ok(false);
            }
            Ok(h.equals_zero(&sub(a, b))? == Some(true))
        }
    }
}

fn steps<H: Host>(h: &mut H, f: &Expr, sym: &Expr, opts: Options) -> R<Option<Vec<Expr>>, H::Error> {
    let (rhs, lhs) = match invert(h, f, sym) {
        Ok(v) => v,
        Err(e) if e.is_ni_or_value() => return Ok(None),
        Err(e) => return Err(e),
    };
    if lhs == *sym {
        return Ok(Some(vec![rhs]));
    }
    // the reference wraps these steps in `except NotImplementedError: pass`
    match inverted_steps(h, f, sym, &lhs, &rhs, opts) {
        Ok(Some(r)) => return Ok(Some(r)),
        Ok(None) => {}
        Err(e) if e.is_ni() => {}
        Err(e) => return Err(e),
    }
    let mut candidates = vec![(rhs, lhs)];
    if opts.bivariate {
        if let Some(r) = bivariate_step(h, f, sym)? {
            return Ok(Some(r));
        }
        // the reference posify retry works on the rewritten equation; the
        // equation before the rewrite is tried next
        match bivariate_eq(f, sym).and_then(|e| invert(h, &e, sym)) {
            Ok(c) => candidates.insert(0, c),
            Err(e) if e.is_ni_or_value() => {}
            Err(e) => return Err(e),
        }
    }
    // the reference ends with a posify retry (`force` is then False)
    if opts.force && !no_force() {
        let mut seen: Vec<Expr> = Vec::new();
        for (rhs, lhs) in candidates {
            let eq = sub(&lhs, &rhs);
            if seen.contains(&eq) {
                continue;
            }
            seen.push(eq.clone());
            let (pos, reps) = crate::simptools::misc::posify(&eq);
            let u = reps.iter().find(|(_, s)| s == sym).map(|(d, _)| d.clone()).unwrap_or_else(|| sym.clone());
            if !has(&pos, &u) {
                continue;
            }
            let soln = {
                let _guard = NoForce::enter();
                h.solve(&pos, &u, Flags::Popped)
            };
            match soln {
                Ok(s) => return Ok(Some(s.iter().map(|s| subs_dict(s, &reps)).collect())),
                Err(e) if e.is_ni() => continue,
                Err(e) => return Err(e),
            }
        }
    }
    Ok(None)
}

