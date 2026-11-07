//! Native port of SymPy 1.14 `sympy/solvers/pde.py`.
//!
//! `pdsolve` (with the `_desolve` hint logic of `sympy/solvers/deutils.py`
//! and `_handle_Integral`), `classify_pde`, `checkpdesol`, the solvers
//! `1st_linear_constant_coeff_homogeneous`, `1st_linear_constant_coeff`
//! (and its `_Integral` form), `1st_linear_variable_coeff` with
//! `_simplify_variable_coeff`, and `pde_separate` / `_separate`, together
//! with `ode_order` and `_preprocess` of `deutils`.
//!
//! The algorithm bodies run here: preprocessing and function detection, the
//! order, the reduction by the smallest power of `f(x, y)`, the matches of
//! the linear forms (native `collect` and `Expr.match`), the hint selection
//! and its errors, the characteristic formulas, the reductions of the
//! variable coefficient method to ODEs, the substitution check and the
//! separation of variables. `simplify`, `expand`, `integrate`, `solve`,
//! `doit` and `dsolve` are the public callables supplied through
//! [`Hooks`] (shared with the `lie_group` port, which calls [`pdsolve`]
//! directly); an error they raise propagates unchanged. A
//! [`LieError::Decline`] leaves the computation to the reference body.
//!
//! The symox deviations of the Python port are kept: the default arbitrary
//! function takes a free name when the equation applies `F`
//! ([`default_solvefun`]), the characteristic coordinates `xi`/`eta` become
//! dummies when the equation contains such symbols, integration variables
//! captured by `eta = funcarg` are renamed, and `checkpdesol` combines held
//! integrals with the same limits.

use crate::diff::diff;
use crate::expr::{Expr, Func, Kind, RelOp};
use crate::functions::undefined;
use crate::matches::reference_match_ordered;
use crate::number::Number;
use crate::ops::{add, add2, div, inv, mul, mul2, neg, pow, sub};
use crate::poly::collect::{self, Collection, Options};
use crate::printer::to_str;
use crate::relational::rel;
use crate::subs::{subs, subs1};
use crate::traverse::{free_symbols, has_structural, preorder};

use super::lie_group::{apply_same, dummy, exp, hk, index_error, truthy, unmodelled, Dsolved, Exc, Hooks, LieError, Solved, R};

/// `allhints` of `pde.py`, in the reference order.
pub const ALLHINTS: [&str; 4] = [
    "1st_linear_constant_coeff_homogeneous",
    "1st_linear_constant_coeff",
    "1st_linear_constant_coeff_Integral",
    "1st_linear_variable_coeff",
];

fn raise<E>(k: Exc, m: impl Into<String>) -> LieError<E> {
    LieError::Raise(k, m.into())
}

/// `e is S.Zero`.
fn is_int_zero(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(Number::Int(0)))
}

fn has(e: &Expr, s: &Expr) -> bool {
    has_structural(e, s)
}

fn has_any(e: &Expr, ss: &[Expr]) -> bool {
    ss.iter().any(|s| has(e, s))
}

fn push_unique(v: &mut Vec<Expr>, e: Expr) {
    if !v.contains(&e) {
        v.push(e);
    }
}

/// An application of a (non-internal, non-host) undefined function: an
/// `AppliedUndef` of the reference.
fn applied_undef(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Undefined(h), _) if !h.internal && h.host().is_none())
}

fn derivatives(e: &Expr) -> Vec<Expr> {
    let mut out = Vec::new();
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Derivative(..)) {
            push_unique(&mut out, s.clone());
        }
    });
    out
}

fn applied_undefs(e: &Expr) -> Vec<Expr> {
    let mut out = Vec::new();
    preorder(e, &mut |s| {
        if applied_undef(s) {
            push_unique(&mut out, s.clone());
        }
    });
    out
}

fn has_class(e: &Expr, pred: impl Fn(&Expr) -> bool) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if !found && pred(s) {
            found = true;
        }
    });
    found
}

fn has_integral(e: &Expr) -> bool {
    has_class(e, |s| matches!(s.kind(), Kind::Func(Func::Integral, _)))
}

fn has_piecewise(e: &Expr) -> bool {
    has_class(e, |s| matches!(s.kind(), Kind::Piecewise(_)))
}

fn eq_rel(a: &Expr, b: &Expr) -> Expr {
    rel(RelOp::Eq, a, b)
}

/// `Eq(a, b)` that must stay an equation (a boolean result declines).
fn equation<E>(a: &Expr, b: &Expr) -> Result<Expr, LieError<E>> {
    let r = eq_rel(a, b);
    if matches!(r.kind(), Kind::Rel(RelOp::Eq, ..)) {
        Ok(r)
    } else {
        Err(LieError::Decline("equation evaluated"))
    }
}

/// `ode_order(expr, func)` of `deutils`.
pub fn ode_order(e: &Expr, f: &Expr) -> u32 {
    if !has(e, f) {
        return 0;
    }
    match e.kind() {
        Kind::Derivative(body, vars) => {
            let n: u32 = vars.iter().map(|(_, c)| *c).sum();
            if body == f {
                n
            } else {
                n + body.args().iter().map(|a| ode_order(a, f)).max().unwrap_or(0)
            }
        }
        _ => e.args().iter().map(|a| ode_order(a, f)).max().unwrap_or(0),
    }
}

/// `_preprocess(expr, func)` of `deutils` (hint `'_Integral'`): the
/// expression with the derivatives of `func` evaluated, and the function.
/// Without `reps`, only the function is computed (the callers that discard
/// the expression).
pub fn preprocess<H: Hooks>(expr: &Expr, func: Option<&Expr>, reps: bool, hooks: &mut H) -> R<(Expr, Expr), H> {
    let mut expr = expr.clone();
    if let Kind::Pow(b, e) = expr.kind() {
        if crate::assumptions::is_positive(e) == Some(true) {
            expr = b.clone();
        }
    }
    let derivs = derivatives(&expr);
    let func = match func {
        Some(f) => f.clone(),
        None => {
            let mut funcs = Vec::new();
            for d in &derivs {
                for a in applied_undefs(d) {
                    push_unique(&mut funcs, a);
                }
            }
            if funcs.len() != 1 {
                return Err(raise(Exc::Value, format!("The function cannot be automatically detected for {}.", to_str(&expr))));
            }
            funcs.pop().expect("one function")
        }
    };
    if !reps {
        return Ok((expr, func));
    }
    let fvars = func.args();
    let mut pairs: Vec<(Expr, Expr)> = Vec::new();
    for d in &derivs {
        let Kind::Derivative(body, vars) = d.kind() else { continue };
        if !(has(d, &func) || vars.iter().any(|(v, _)| fvars.contains(v))) {
            continue;
        }
        // `_doit_derivative`: a derivative of an applied function in its own
        // symbol arguments is held.
        let args = body.args();
        if applied_undef(body) && vars.iter().all(|(v, _)| v.is_symbol() && args.contains(v)) {
            continue;
        }
        let v = hk::<H, _>(hooks.doit(d))?;
        if v != *d {
            pairs.push((d.clone(), v));
        }
    }
    // The reference substitutes in the iteration order of a set; nested
    // replaced derivatives would depend on it.
    for (i, (a, _)) in pairs.iter().enumerate() {
        if pairs.iter().enumerate().any(|(j, (b, _))| i != j && has(a, b)) {
            return Err(LieError::Decline("nested derivative replacements"));
        }
    }
    Ok((subs(&expr, &pairs), func))
}

/// The match dict of a hint: the Wild bindings (in binding order) and the
/// Wilds stored under `'b'`, `'c'`, `'d'` (and `'e'`).
#[derive(Clone, Debug)]
pub struct Match {
    pub entries: Vec<(Expr, Expr)>,
    pub b: Expr,
    pub c: Expr,
    pub d: Expr,
    pub e: Option<Expr>,
}

impl Match {
    fn get<E>(&self, w: &Expr) -> Result<Expr, LieError<E>> {
        self.entries.iter().find(|(k, _)| k == w).map(|(_, v)| v.clone()).ok_or_else(|| LieError::KeyMissing(w.clone()))
    }
}

/// The result of `classify_pde`: `tuple_only` when the reference returns a
/// tuple whatever `dict` (the recursion on `Eq(lhs, rhs)`). All matching
/// hints share `matched`.
#[derive(Clone, Debug)]
pub struct Classified {
    pub order: u32,
    pub tuple_only: bool,
    pub hints: Vec<&'static str>,
    pub matched: Option<Match>,
}

fn lookup<E>(m: &[(Expr, Expr)], w: &Expr) -> Result<Expr, LieError<E>> {
    m.iter().find(|(k, _)| k == w).map(|(_, v)| v.clone()).ok_or_else(|| LieError::KeyMissing(w.clone()))
}

fn matched(e: &Expr, pattern: &Expr) -> Result<Option<Vec<(Expr, Expr)>>, &'static str> {
    match reference_match_ordered(e, pattern, false) {
        Ok(Some(b)) => Ok(Some(b.into_vec())),
        Ok(None) => Ok(None),
        Err(_) => Err("match"),
    }
}

/// The function `f` of `f(...)` and its arguments `x`, `y`.
fn split_func<E>(func: &Expr) -> Result<(Expr, Expr), LieError<E>> {
    if !applied_undef(func) {
        return Err(LieError::Decline("function head"));
    }
    let args = func.args();
    if args.len() < 2 {
        return Err(raise(Exc::Index, "tuple index out of range"));
    }
    if !args[0].is_symbol() || !args[1].is_symbol() {
        return Err(LieError::Decline("function arguments"));
    }
    Ok((args[0].clone(), args[1].clone()))
}

/// `classify_pde(eq, func, prep=prep)`.
pub fn classify_pde<H: Hooks>(eq: &Expr, func: Option<&Expr>, prep: bool, hooks: &mut H) -> R<Classified, H> {
    if unmodelled(eq) || func.is_some_and(unmodelled) {
        return Err(LieError::Decline("unmodelled content"));
    }
    if let Some(f) = func {
        if f.args().len() != 2 {
            return Err(raise(Exc::NotImplemented, "Right now only partial differential equations of two variables are supported"));
        }
    }
    let mut func = func.cloned();
    if prep || func.is_none() {
        let (_, f) = preprocess(eq, func.as_ref(), false, hooks)?;
        if func.is_none() {
            func = Some(f);
        }
    }
    let func = func.expect("function");
    let eq = match eq.kind() {
        Kind::Rel(RelOp::Eq, l, r) => {
            if !is_int_zero(r) {
                let mut c = classify_pde(&sub(l, r), Some(&func), true, hooks)?;
                c.tuple_only = true;
                return Ok(c);
            }
            l.clone()
        }
        Kind::Rel(..) | Kind::Bool(_) => return Err(LieError::Decline("relational")),
        _ => eq.clone(),
    };
    let (x, y) = split_func(&func)?;
    let fxy = apply_same(&func, vec![x.clone(), y.clone()]).ok_or(LieError::Decline("function head"))?;
    let fx = diff(&fxy, &x);
    let fy = diff(&fxy, &y);
    let order = ode_order(&eq, &fxy);
    let mut out = Classified { order, tuple_only: false, hints: Vec::new(), matched: None };
    if order == 0 {
        return Ok(out);
    }
    let eq = hk::<H, _>(hooks.expand(&eq))?;
    if unmodelled(&eq) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let a = Expr::wild("a", vec![fxy.clone()]);
    let strict = vec![fxy.clone(), fx.clone(), fy.clone(), x.clone(), y.clone()];
    let b = Expr::wild("b", strict.clone());
    let c = Expr::wild("c", strict.clone());
    let d = Expr::wild("d", strict);
    let e = Expr::wild("e", vec![fxy.clone(), fx.clone(), fy.clone()]);
    let n = Expr::wild("n", vec![x.clone(), y.clone()]);
    let mut reduced = eq.clone();
    if matches!(eq.kind(), Kind::Add(..)) {
        // The minimum of the powers `n` of `a*f**n` over the coefficients of
        // the highest derivatives (a set iteration in the reference: ties
        // and undecidable comparisons decline).
        let pattern = mul2(&a, &pow(&fxy, &n));
        let mut candidates: Vec<Expr> = Vec::new();
        for i in (0..=order).rev() {
            let mut target = fxy.clone();
            for _ in 0..i {
                target = diff(&target, &x);
            }
            for _ in i..order {
                target = diff(&target, &y);
            }
            let co = super::classify::coeff(&eq, &target);
            if co == Expr::one() {
                continue;
            }
            let Some(m) = matched(&co, &pattern).map_err(LieError::Decline)? else { continue };
            if m.is_empty() {
                continue;
            }
            let va = lookup::<()>(&m, &a).map_err(|_| LieError::Decline("power match"))?;
            if truthy(&va) {
                candidates.push(lookup::<()>(&m, &n).map_err(|_| LieError::Decline("power match"))?);
            }
        }
        let mut power: Option<Expr> = None;
        for (i, p) in candidates.iter().enumerate() {
            let minimal = candidates.iter().enumerate().all(|(j, q)| {
                i == j || p == q || matches!(rel(RelOp::Lt, p, q).kind(), Kind::Bool(true))
            });
            let comparable = candidates.iter().all(|q| p == q || matches!(rel(RelOp::Lt, p, q).kind(), Kind::Bool(_)));
            if !comparable {
                return Err(LieError::Decline("power comparison"));
            }
            if minimal {
                power = Some(p.clone());
                break;
            }
        }
        if !candidates.is_empty() && power.is_none() {
            return Err(LieError::Decline("power comparison"));
        }
        if let Some(p) = power.filter(truthy) {
            let den = pow(&fxy, &p);
            reduced = add(&eq.args().iter().map(|t| div(t, &den)).collect::<Vec<_>>());
        }
    }
    if order == 1 {
        let collected = match collect::run::<()>(&reduced, std::slice::from_ref(&fxy), Options { exact: false, discover: false, distribute_order_term: true }, true, None) {
            Ok(Collection::Expression(v)) => v,
            _ => return Err(LieError::Decline("collect")),
        };
        let pattern = |b: &Expr, c: &Expr, d: &Expr| add(&[mul2(b, &fx), mul2(c, &fy), mul2(d, &fxy), e.clone()]);
        match matched(&collected, &pattern(&b, &c, &d)).map_err(LieError::Decline)? {
            Some(r) if !r.is_empty() => {
                let re = lookup(&r, &e)?;
                if !truthy(&re) {
                    out.hints.push(ALLHINTS[0]);
                    out.matched = Some(Match { entries: r, b, c, d, e: None });
                } else {
                    let (rb, rc) = (lookup(&r, &b)?, lookup(&r, &c)?);
                    if !is_int_zero(&add2(&pow(&rb, &Expr::int(2)), &pow(&rc, &Expr::int(2)))) {
                        out.hints.push(ALLHINTS[1]);
                        out.hints.push(ALLHINTS[2]);
                        out.matched = Some(Match { entries: r, b, c, d, e: Some(e) });
                    }
                }
            }
            _ => {
                let loose = vec![fxy.clone(), fx.clone(), fy.clone()];
                let b = Expr::wild("b", loose.clone());
                let c = Expr::wild("c", loose.clone());
                let d = Expr::wild("d", loose);
                if let Some(r) = matched(&collected, &pattern(&b, &c, &d)).map_err(LieError::Decline)?.filter(|r| !r.is_empty()) {
                    out.hints.push(ALLHINTS[3]);
                    out.matched = Some(Match { entries: r, b, c, d, e: Some(e) });
                }
            }
        }
    }
    Ok(out)
}

/// The arbitrary function `solvefun` applied to `arg`.
fn apply_fun<E>(solvefun: &Expr, arg: Expr) -> Result<Expr, LieError<E>> {
    apply_same(solvefun, vec![arg]).ok_or(LieError::Decline("solvefun"))
}

/// The default arbitrary function of the symox `pdsolve`: `F`, or the first
/// of `F_`, `F__`, ... that the equation does not apply, as a sample
/// application `F(x)` (only its head is used).
pub fn default_solvefun(eq: &Expr, func: &Expr) -> Expr {
    let mut names: Vec<String> = Vec::new();
    preorder(eq, &mut |s| {
        if let Kind::Func(Func::Undefined(h), _) = s.kind() {
            if !h.internal {
                names.push(h.name.to_string());
            }
        }
    });
    let mut name = String::from("F");
    while names.contains(&name) {
        name.push('_');
    }
    let arg = func.args().into_iter().next().unwrap_or_else(|| Expr::symbol("x"));
    undefined(&name, vec![arg])
}

/// `pde_<hint>(eq, func, order, match, solvefun)` for the base hint (the
/// `_Integral` suffix removed), with `b`, `c`, `d`, `e` the matched values.
#[allow(clippy::too_many_arguments)]
pub fn solve_hint<H: Hooks>(hint: &str, eq: &Expr, func: &Expr, b: &Expr, c: &Expr, d: &Expr, e: Option<&Expr>, solvefun: &Expr, hooks: &mut H) -> R<Expr, H> {
    if [eq, func, b, c, d].iter().any(|v| unmodelled(v)) || e.is_some_and(unmodelled) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let (x, y) = split_func(func)?;
    let fxy = apply_same(func, vec![x.clone(), y.clone()]).ok_or(LieError::Decline("function head"))?;
    match hint {
        "1st_linear_constant_coeff_homogeneous" => {
            let den = add2(&pow(b, &Expr::int(2)), &pow(c, &Expr::int(2)));
            let arg = mul2(&div(&neg(d), &den), &add2(&mul2(b, &x), &mul2(c, &y)));
            let rhs = mul2(&exp(&arg), &apply_fun(solvefun, sub(&mul2(c, &x), &mul2(b, &y)))?);
            Ok(eq_rel(&fxy, &rhs))
        }
        "1st_linear_constant_coeff" => {
            let e = neg(e.ok_or_else(|| LieError::KeyMissing(Expr::symbol("e")))?);
            constant_coeff(eq, &fxy, &x, &y, b, c, d, &e, solvefun, hooks)
        }
        "1st_linear_variable_coeff" => {
            let e = neg(e.ok_or_else(|| LieError::KeyMissing(Expr::symbol("e")))?);
            variable_coeff(&fxy, &x, &y, b, c, d, &e, solvefun, hooks)
        }
        _ => Err(LieError::Decline("hint")),
    }
}

#[allow(clippy::too_many_arguments)]
fn constant_coeff<H: Hooks>(eq: &Expr, fxy: &Expr, x: &Expr, y: &Expr, b: &Expr, c: &Expr, d: &Expr, e: &Expr, solvefun: &Expr, hooks: &mut H) -> R<Expr, H> {
    let (mut xi, mut eta) = (Expr::symbol("xi"), Expr::symbol("eta"));
    // symox: the characteristic coordinates of the reference capture symbols
    // named xi or eta of the equation; fresh dummies are used in that case.
    if has(eq, &xi) || has(eq, &eta) {
        xi = dummy("xi");
        eta = dummy("eta");
    }
    let den = add2(&pow(b, &Expr::int(2)), &pow(c, &Expr::int(2)));
    let expterm = exp(&mul2(&div(&neg(d), &den), &xi));
    let functerm = apply_fun(solvefun, eta.clone())?;
    let s1 = add2(&mul2(b, x), &mul2(c, y));
    let s2 = sub(&mul2(c, x), &mul2(b, y));
    let system = [sub(&s1, &xi), sub(&s2, &eta)];
    let pairs = match hk::<H, _>(hooks.solve_system(&system, &[x.clone(), y.clone()], false))? {
        Solved::Dict(p) => p,
        Solved::Dicts(v) if v.is_empty() => Vec::new(),
        Solved::Dicts(_) => return Err(LieError::Decline("solve result")),
    };
    // The values do not contain x or y: the dict order of the reference
    // substitution is immaterial.
    if pairs.iter().any(|(_, v)| has(v, x) || has(v, y)) {
        return Err(LieError::Decline("solve result"));
    }
    let body = subs(&mul2(&inv(&expterm), e), &pairs);
    if has_piecewise(&body) {
        return Err(LieError::Decline("piecewise integrand"));
    }
    let integral = crate::binders::held(Func::Integral, &body, vec![vec![xi.clone(), s1.clone()]]).map_err(|_| LieError::Decline("integral"))?;
    let genterm = mul2(&inv(&den), &integral);
    let held = crate::subs::held_subs(&mul2(&expterm, &add2(&functerm, &genterm)), vec![eta, xi], vec![s2, s1]);
    Ok(eq_rel(fxy, &held))
}

fn free_minus(e: &Expr, others: &[&Expr], x: &Expr, y: &Expr) -> Vec<Expr> {
    let mut drop: Vec<Expr> = vec![x.clone(), y.clone()];
    for o in others {
        drop.extend(free_symbols(o));
    }
    free_symbols(e).into_iter().filter(|s| !drop.contains(s)).collect()
}

fn one_solution<E>(s: Dsolved) -> Result<Expr, LieError<E>> {
    match s {
        Dsolved::One(e) => Ok(e),
        Dsolved::Many(_) => Err(LieError::Decline("dsolve list")),
    }
}

fn rhs_of<E>(e: &Expr) -> Result<Expr, LieError<E>> {
    match e.kind() {
        Kind::Rel(RelOp::Eq, _, r) => Ok(r.clone()),
        _ => Err(LieError::Decline("dsolve result")),
    }
}

#[allow(clippy::too_many_arguments)]
fn variable_coeff<H: Hooks>(fxy: &Expr, x: &Expr, y: &Expr, b: &Expr, c: &Expr, d: &Expr, e: &Expr, solvefun: &Expr, hooks: &mut H) -> R<Expr, H> {
    let eta = Expr::symbol("eta");
    let integrate = |hooks: &mut H, v: &Expr, s: &Expr| -> R<Expr, H> {
        match hooks.integrate(v, Some(s)) {
            Ok(t) => Ok(t),
            Err(err) if hooks.exc(&err) == Exc::NotImplemented => {
                Err(raise(Exc::NotImplemented, "Unable to find a solution due to inability of integrate"))
            }
            Err(err) => Err(LieError::Hook(err)),
        }
    };
    if !truthy(d) && !(truthy(b) && truthy(c)) {
        if truthy(c) {
            let tsol = integrate(hooks, &div(e, c), y)?;
            return Ok(eq_rel(fxy, &add2(&apply_fun(solvefun, x.clone())?, &tsol)));
        }
        if truthy(b) {
            let tsol = integrate(hooks, &div(e, b), x)?;
            return Ok(eq_rel(fxy, &add2(&apply_fun(solvefun, y.clone())?, &tsol)));
        }
    }
    for (t, coef, funcarg) in [(x, c, y), (y, b, x)] {
        // `not c`: b*u_x + d*u = e is a linear ODE in x (`not b`: in y).
        if truthy(coef) {
            continue;
        }
        let other = if t == x { b } else { c };
        let ft = apply_same(fxy, vec![t.clone()]).ok_or(LieError::Decline("function head"))?;
        let plode = sub(&add2(&mul2(&diff(&ft, t), other), &mul2(d, &ft)), e);
        let sol = one_solution(hk::<H, _>(hooks.dsolve(&plode, &ft, None))?)?;
        let syms = free_minus(&sol, &[&plode], x, y);
        let rhs = simplify_variable_coeff(&rhs_of(&sol)?, &syms, solvefun, funcarg, hooks)?;
        return Ok(eq_rel(fxy, &rhs));
    }
    let dfun = undefined("d", vec![x.clone()]);
    let h = subs1(&div(c, b), y, &dfun);
    let sol = match hk::<H, _>(hooks.dsolve(&sub(&diff(&dfun, x), &h), &dfun, None))? {
        Dsolved::One(s) => s,
        Dsolved::Many(v) => v.into_iter().next().ok_or_else(index_error)?,
    };
    let solsym = free_minus(&sol, &[&h], x, y);
    if solsym.len() != 1 {
        return Err(raise(Exc::NotImplemented, "Cannot solve the partial differential equation due to inability of constantsimp"));
    }
    let etat = subs1(&hk::<H, _>(hooks.solve(&sol, &solsym[0]))?.into_iter().next().ok_or_else(index_error)?, &dfun, y);
    let ysub = hk::<H, _>(hooks.solve(&sub(&eta, &etat), y))?.into_iter().next().ok_or_else(index_error)?;
    let fx1 = apply_same(fxy, vec![x.clone()]).ok_or(LieError::Decline("function head"))?;
    let deq = subs1(&sub(&add2(&mul2(b, &diff(&fx1, x)), &mul2(d, &fx1)), e), y, &ysub);
    let fin = rhs_of(&one_solution(hk::<H, _>(hooks.dsolve(&deq, &fx1, Some("1st_linear")))?)?)?;
    let finsyms = free_minus(&fin, &[&deq], x, y);
    let rhs = simplify_variable_coeff(&fin, &finsyms, solvefun, &etat, hooks)?;
    Ok(eq_rel(fxy, &rhs))
}

/// `_simplify_variable_coeff(sol, syms, func, funcarg)`.
fn simplify_variable_coeff<H: Hooks>(sol: &Expr, syms: &[Expr], solvefun: &Expr, funcarg: &Expr, hooks: &mut H) -> R<Expr, H> {
    let eta = Expr::symbol("eta");
    let fin = match syms {
        [] => return Err(raise(Exc::UnboundLocal, "local variable 'final' referenced before assignment")),
        [sym] => subs1(sol, sym, &apply_fun(solvefun, funcarg.clone())?),
        // The reference keeps the substitution of the last symbol of a set.
        _ => return Err(LieError::Decline("several constants")),
    };
    let fin = rename_captured(&fin, &eta, funcarg);
    hk::<H, _>(hooks.simplify(&subs1(&fin, &eta, funcarg)))
}

/// symox `_rename_captured_integration_variables`: in the held integrals
/// whose free symbols contain `eta`, an integration variable that occurs in
/// `funcarg` becomes a dummy of the same name (an indefinite limit `(v,)`
/// becomes `(v_dummy, v)`), so that `eta = funcarg` does not capture it.
fn rename_captured(e: &Expr, eta: &Expr, funcarg: &Expr) -> Expr {
    let outer = free_symbols(funcarg);
    fn rebuild(node: &Expr, eta: &Expr, outer: &[Expr]) -> Expr {
        let node = crate::traverse::map_children(node, &mut |a: &Expr| rebuild(a, eta, outer));
        let Kind::Func(Func::Integral, args) = node.kind() else { return node };
        if !free_symbols(&node).contains(eta) {
            return node;
        }
        let mut body = args[0].clone();
        let mut limits = Vec::new();
        for limit in crate::binders::limits(args) {
            let variable = limit[0].clone();
            if outer.contains(&variable) {
                let fresh = dummy(variable.as_symbol().unwrap_or("x"));
                let mut map = rustc_hash::FxHashMap::default();
                map.insert(variable.clone(), fresh.clone());
                body = crate::subs::xreplace(&body, &map);
                if limit.len() == 1 {
                    limits.push(vec![fresh, variable]);
                } else {
                    let mut l = vec![fresh];
                    l.extend(limit[1..].iter().cloned());
                    limits.push(l);
                }
            } else {
                limits.push(limit);
            }
        }
        crate::binders::held(Func::Integral, &body, limits).unwrap_or(node)
    }
    rebuild(e, eta, &outer)
}

/// The result of `pdsolve`: one solution, the dict of `hint='all'` (each
/// hint with its solution or its `NotImplementedError`), or the
/// `NotImplementedError` of the classification step (which the binding
/// raises, the symox `_native_fallback` handling it for `hint='default'`).
pub enum Pdsolved<E> {
    One(Expr),
    All { order: u32, default: Option<&'static str>, entries: Vec<(&'static str, Result<Expr, LieError<E>>)> },
    Unclassified(String),
}

/// `_helper_simplify` and `_handle_Integral` for one hint.
fn helper_simplify<H: Hooks>(hint: &str, eq: &Expr, func: &Expr, m: &Match, solvefun: &Expr, hooks: &mut H) -> R<Expr, H> {
    let base = hint.strip_suffix("_Integral").unwrap_or(hint);
    let (b, c, d) = (m.get(&m.b)?, m.get(&m.c)?, m.get(&m.d)?);
    let e = match &m.e {
        Some(w) => Some(m.get(w)?),
        None => None,
    };
    let expr = solve_hint(base, eq, func, &b, &c, &d, e.as_ref(), solvefun, hooks)?;
    if hint == "1st_linear_constant_coeff" {
        let v = hk::<H, _>(hooks.doit(&expr))?;
        return hk::<H, _>(hooks.simplify(&v));
    }
    Ok(expr)
}

/// `pdsolve(eq, func, hint, solvefun=solvefun, prep=prep)`: `solvefun` is a
/// sample application of the arbitrary function.
pub fn pdsolve<H: Hooks>(eq: &Expr, func: Option<&Expr>, hint: &str, prep: bool, solvefun: &Expr, hooks: &mut H) -> R<Pdsolved<H::Error>, H> {
    if unmodelled(eq) || func.is_some_and(unmodelled) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let eq = match eq.kind() {
        Kind::Rel(RelOp::Eq, l, r) => sub(l, r),
        Kind::Rel(..) | Kind::Bool(_) => return Err(LieError::Decline("relational")),
        _ => eq.clone(),
    };
    let (eq, func) = match func {
        Some(f) if !prep => (eq, f.clone()),
        _ => preprocess(&eq, func, true, hooks)?,
    };
    let cls = match classify_pde(&eq, Some(&func), false, hooks) {
        Err(LieError::Raise(Exc::NotImplemented, m)) => return Ok(Pdsolved::Unclassified(m)),
        r => r?,
    };
    let Some(&default) = cls.hints.first() else {
        if !ALLHINTS.contains(&hint) && hint != "default" {
            return Err(raise(Exc::Value, format!("Hint not recognized: {hint}")));
        }
        if cls.order == 0 {
            // `hints['ordered_hints']` is evaluated first and is missing.
            return Err(raise(Exc::Key, "ordered_hints"));
        }
        if hint != "default" {
            return Err(raise(Exc::Value, format!("PDE {} does not match hint {}", to_str(&eq), hint)));
        }
        return Ok(Pdsolved::Unclassified(format!("psolve: Cannot solve {}", to_str(&eq))));
    };
    let m = cls.matched.clone().ok_or(LieError::Decline("match"))?;
    match hint {
        "default" => Ok(Pdsolved::One(helper_simplify(default, &eq, &func, &m, solvefun, hooks)?)),
        "all" | "all_Integral" | "best" => {
            let mut gethints: Vec<&'static str> = cls.hints.clone();
            if hint == "all_Integral" {
                for h in &cls.hints {
                    if let Some(base) = h.strip_suffix("_Integral") {
                        gethints.retain(|g| *g != base);
                    }
                }
            }
            let g = classify_pde(&eq, None, true, hooks)?;
            let mut entries = Vec::new();
            for h in gethints {
                match helper_simplify(h, &eq, &func, &m, solvefun, hooks) {
                    Ok(v) => entries.push((h, Ok(v))),
                    Err(err) if super::lie_group::class(hooks, &err) == Exc::NotImplemented => entries.push((h, Err(err))),
                    Err(err) => return Err(err),
                }
            }
            Ok(Pdsolved::All { order: g.order, default: g.hints.first().copied(), entries })
        }
        _ if !ALLHINTS.contains(&hint) => Err(raise(Exc::Value, format!("Hint not recognized: {hint}"))),
        _ if !cls.hints.contains(&hint) => Err(raise(Exc::Value, format!("PDE {} does not match hint {}", to_str(&eq), hint))),
        _ => Ok(Pdsolved::One(helper_simplify(hint, &eq, &func, &m, solvefun, hooks)?)),
    }
}

/// A `checkpdesol` solution argument: one expression or a sequence.
pub enum Sol {
    One(Expr),
    Seq(Vec<Expr>),
}

/// The result of `checkpdesol`: `(bool, expr)` or a list of results.
#[derive(Debug)]
pub enum Checked {
    One(bool, Expr),
    List(Vec<Checked>),
}

/// The raw message of the final `NotImplementedError` of `checkpdesol`,
/// before `filldedent`.
pub const UNABLE_PREFIX: &str = "\n        Unable to test if ";

/// `checkpdesol(pde, sol, func, solve_for_func)`.
pub fn checkpdesol<H: Hooks>(pde: &Expr, sol: &Sol, func: Option<&Expr>, solve_for_func: bool, hooks: &mut H) -> R<Checked, H> {
    let sols: Vec<&Expr> = match sol {
        Sol::One(s) => vec![s],
        Sol::Seq(v) => v.iter().collect(),
    };
    if unmodelled(pde) || func.is_some_and(unmodelled) || sols.iter().any(|s| unmodelled(s)) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let pde = match pde.kind() {
        Kind::Rel(RelOp::Eq, ..) => pde.clone(),
        Kind::Rel(..) | Kind::Bool(_) => return Err(LieError::Decline("relational")),
        _ => equation(pde, &Expr::zero())?,
    };
    let Kind::Rel(_, lhs, _) = pde.kind() else { unreachable!() };
    let func = match func {
        Some(f) => f.clone(),
        None => match preprocess(lhs, None, false, hooks) {
            Ok((_, f)) => f,
            Err(LieError::Raise(Exc::Value, _)) => {
                return Err(if matches!(sol, Sol::Seq(v) if v.is_empty()) {
                    raise(Exc::Value, "must pass func arg to checkpdesol for this case.")
                } else {
                    raise(Exc::Type, "unhashable type: 'set'")
                });
            }
            Err(e) => return Err(e),
        },
    };
    match sol {
        Sol::One(s) => check_one(&pde, s, &func, solve_for_func, hooks),
        Sol::Seq(v) => Ok(Checked::List(v.iter().map(|s| check_one(&pde, s, &func, solve_for_func, hooks)).collect::<Result<_, _>>()?)),
    }
}

fn check_one<H: Hooks>(pde: &Expr, sol: &Expr, func: &Expr, solve_for_func: bool, hooks: &mut H) -> R<Checked, H> {
    let sol = match sol.kind() {
        Kind::Rel(RelOp::Eq, l, r) => {
            if r == func {
                Expr::from_kind(Kind::Rel(RelOp::Eq, r.clone(), l.clone()))
            } else {
                sol.clone()
            }
        }
        Kind::Rel(..) | Kind::Bool(_) => return Err(LieError::Decline("relational solution")),
        _ => equation(func, sol)?,
    };
    let Kind::Rel(_, lhs, rhs) = sol.kind() else { unreachable!() };
    let solved = lhs == func && !has(rhs, func);
    if solve_for_func && !solved {
        let roots = hk::<H, _>(hooks.solve(&sol, func))?;
        match roots.len() {
            0 => {}
            1 => return check_one(pde, &equation(func, &roots[0])?, func, false, hooks),
            _ => {
                let mut out = Vec::new();
                for t in &roots {
                    out.push(check_one(pde, &equation(func, t)?, func, false, hooks)?);
                }
                return Ok(Checked::List(out));
            }
        }
    }
    if lhs == func {
        let Kind::Rel(_, pl, pr) = pde.kind() else { unreachable!() };
        let p = sub(pl, pr);
        let v = hk::<H, _>(hooks.doit(&subs1(&p, func, rhs)))?;
        let s = hk::<H, _>(hooks.simplify(&v))?;
        if !is_int_zero(&s) && has_integral(&s) {
            // symox: held integrals with the same limits are combined by
            // linearity, so their integrands can cancel.
            if is_int_zero(&combine_integrals(&s, hooks)?) {
                return Ok(Checked::One(true, Expr::zero()));
            }
        }
        return Ok(Checked::One(is_int_zero(&s), s));
    }
    Err(raise(Exc::NotImplemented, format!("{UNABLE_PREFIX}{} is a solution to {}.", to_str(&sol), to_str(pde))))
}

/// symox `_combine_integrals`.
fn combine_integrals<H: Hooks>(e: &Expr, hooks: &mut H) -> R<Expr, H> {
    let expanded = hk::<H, _>(hooks.expand(e))?;
    let terms = if matches!(expanded.kind(), Kind::Add(..)) { expanded.args() } else { vec![expanded] };
    let mut groups: Vec<(Vec<Expr>, Vec<Expr>)> = Vec::new();
    let mut rest = Vec::new();
    for term in terms {
        let factors = if matches!(term.kind(), Kind::Mul(..)) { term.args() } else { vec![term.clone()] };
        let integrals: Vec<&Expr> = factors.iter().filter(|f| matches!(f.kind(), Kind::Func(Func::Integral, _))).collect();
        if integrals.len() != 1 {
            rest.push(term);
            continue;
        }
        let integral = integrals[0];
        let Kind::Func(_, args) = integral.kind() else { unreachable!() };
        let coefficient = div(&term, integral);
        let limits = args[1..].to_vec();
        let bound: Vec<Expr> = limits.iter().filter_map(|l| l.args().into_iter().next()).collect();
        if free_symbols(&coefficient).iter().any(|s| bound.contains(s)) {
            rest.push(term);
            continue;
        }
        let part = mul2(&coefficient, &args[0]);
        match groups.iter_mut().find(|(k, _)| *k == limits) {
            Some((_, parts)) => parts.push(part),
            None => groups.push((limits, vec![part])),
        }
    }
    let mut combined = Vec::new();
    for (limits, parts) in groups {
        let body = hk::<H, _>(hooks.simplify(&add(&parts)))?;
        if has_piecewise(&body) {
            return Err(LieError::Decline("piecewise integrand"));
        }
        if is_int_zero(&body) {
            continue;
        }
        let lims = limits.iter().map(|l| l.args()).collect();
        combined.push(crate::binders::held(Func::Integral, &body, lims).map_err(|_| LieError::Decline("integral"))?);
    }
    hk::<H, _>(hooks.simplify(&add2(&add(&rest), &add(&combined))))
}

/// The arguments of an expression iterated by `_separate`, when they are
/// expressions (held and relational nodes decline).
fn plain_args<E>(e: &Expr) -> Result<Vec<Expr>, LieError<E>> {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => Ok(e.args()),
        Kind::Func(f, _) if !matches!(f, Func::Integral | Func::Sum | Func::Product | Func::Subs | Func::Limit) && !crate::binders::is_calculus(f) => {
            if matches!(f, Func::Undefined(h) if h.internal) {
                return Err(LieError::Decline("internal function"));
            }
            Ok(e.args())
        }
        _ => Err(LieError::Decline("separate arguments")),
    }
}

const REDUCE_EMPTY: &str = "reduce() of empty iterable with no initial value";

/// `pde_separate(eq, fun, sep, strategy)`: `Ok(None)` for the reference
/// `None` (no separation).
pub fn pde_separate<H: Hooks>(eq: &Expr, fun: &Expr, sep: &[Expr], strategy: &str, hooks: &mut H) -> R<Option<[Expr; 2]>, H> {
    let do_add = match strategy {
        "add" => true,
        "mul" => false,
        _ => return Err(raise(Exc::Value, format!("Unknown strategy: {strategy}"))),
    };
    if unmodelled(eq) || unmodelled(fun) || sep.iter().any(unmodelled) {
        return Err(LieError::Decline("unmodelled content"));
    }
    let lhs = match eq.kind() {
        Kind::Rel(RelOp::Eq, l, r) => {
            if !is_int_zero(r) {
                return pde_separate(&equation(&sub(l, r), &Expr::zero())?, fun, sep, strategy, hooks);
            }
            l.clone()
        }
        Kind::Rel(..) | Kind::Bool(_) => return Err(LieError::Decline("relational")),
        _ => return pde_separate(&equation(eq, &Expr::zero())?, fun, sep, strategy, hooks),
    };
    let orig_args = fun.args();
    let subs_args: Vec<Expr> = sep.iter().flat_map(|s| s.args()).collect();
    if sep.is_empty() {
        return Err(raise(Exc::Type, REDUCE_EMPTY));
    }
    let functions = if do_add { add(sep) } else { mul(sep) };
    if subs_args.len() != orig_args.len() {
        return Err(raise(Exc::Value, "Variable counts do not match"));
    }
    let mut seen = Vec::new();
    for a in &subs_args {
        if seen.contains(a) {
            return Err(raise(Exc::Value, "Duplicate substitution arguments detected"));
        }
        seen.push(a.clone());
    }
    if orig_args.iter().any(|a| !subs_args.contains(a)) || subs_args.iter().any(|a| !orig_args.contains(a)) {
        return Err(raise(Exc::Value, "Arguments do not match"));
    }
    let mut result = hk::<H, _>(hooks.doit(&subs1(&lhs, fun, &functions)))?;
    if unmodelled(&result) {
        return Err(LieError::Decline("unmodelled content"));
    }
    if !do_add {
        let args = plain_args(&result)?;
        if args.is_empty() {
            return Err(LieError::Decline("separate arguments"));
        }
        result = add(&args.iter().map(|a| div(a, &functions)).collect::<Vec<_>>());
    }
    separate(&result, &subs_args[0], &subs_args[1..], hooks)
}

/// `_separate(eq, dep, others)`.
fn separate<H: Hooks>(eq: &Expr, dep: &Expr, others: &[Expr], hooks: &mut H) -> R<Option<[Expr; 2]>, H> {
    let args = plain_args(eq)?;
    let mut terms: Vec<Expr> = Vec::new();
    for term in &args {
        match term.kind() {
            Kind::Mul(..) => {
                for i in term.args() {
                    if matches!(i.kind(), Kind::Derivative(..)) && !has_any(&i, others) {
                        push_unique(&mut terms, term.clone());
                    }
                }
            }
            Kind::Derivative(..) if !has_any(term, others) => push_unique(&mut terms, term.clone()),
            _ => {}
        }
    }
    let mut divs: Vec<Expr> = Vec::new();
    for term in &terms {
        let expanded = hk::<H, _>(hooks.expand(term))?;
        let (ext, sep) = crate::decompose::as_independent(&expanded, std::slice::from_ref(dep), None);
        if has_any(&sep, others) {
            return Ok(None);
        }
        push_unique(&mut divs, ext);
    }
    let mut eq = eq.clone();
    if !divs.is_empty() {
        let mut parts = Vec::new();
        for term in &args {
            let s = add(&divs.iter().map(|i| div(term, i)).collect::<Vec<_>>());
            parts.push(hk::<H, _>(hooks.simplify(&s))?);
        }
        eq = add(&parts);
    }
    let mut divs: Vec<Expr> = Vec::new();
    let (mut lhs, mut rhs) = (Expr::zero(), Expr::zero());
    for term in plain_args(&eq)? {
        if !has_any(&term, others) {
            lhs = add2(&lhs, &term);
            continue;
        }
        let expanded = hk::<H, _>(hooks.expand(&term))?;
        let (_, sep) = crate::decompose::as_independent(&expanded, std::slice::from_ref(dep), None);
        if has_any(&sep, others) {
            return Ok(None);
        }
        push_unique(&mut divs, sep);
        let expanded = hk::<H, _>(hooks.expand(&term))?;
        rhs = sub(&rhs, &expanded);
    }
    if divs.is_empty() {
        return Err(raise(Exc::Type, REDUCE_EMPTY));
    }
    let fulldiv = add(&divs);
    let l = hk::<H, _>(hooks.simplify(&div(&lhs, &fulldiv)))?;
    let lhs = hk::<H, _>(hooks.expand(&l))?;
    let r = hk::<H, _>(hooks.simplify(&div(&rhs, &fulldiv)))?;
    let rhs = hk::<H, _>(hooks.expand(&r))?;
    if has_any(&lhs, others) || has(&rhs, dep) {
        return Ok(None);
    }
    Ok(Some([lhs, rhs]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_counts_held_derivatives() {
        let (x, y) = (Expr::symbol("x"), Expr::symbol("y"));
        let f = undefined("f", vec![x.clone(), y.clone()]);
        let fxy = diff(&diff(&f, &x), &y);
        assert_eq!(ode_order(&add2(&fxy, &x), &f), 2);
        assert_eq!(ode_order(&x, &f), 0);
        assert_eq!(ode_order(&mul2(&f, &diff(&f, &x)), &f), 1);
    }

    #[test]
    fn default_solvefun_avoids_the_equation_functions() {
        let (x, y) = (Expr::symbol("x"), Expr::symbol("y"));
        let f = undefined("f", vec![x.clone(), y.clone()]);
        let big_f = undefined("F", vec![sub(&x, &y)]);
        let s = default_solvefun(&add2(&diff(&f, &x), &big_f), &f);
        assert!(matches!(s.kind(), Kind::Func(Func::Undefined(h), _) if &*h.name == "F_"));
        let s = default_solvefun(&diff(&f, &x), &f);
        assert!(matches!(s.kind(), Kind::Func(Func::Undefined(h), _) if &*h.name == "F"));
    }
}
