//! `piecewise_fold` and the integration of Piecewise functions
//! (`Piecewise._eval_integral`, `Piecewise._eval_interval`).

use crate::{Expr, Kind};

/// `Piecewise(*args, evaluate=False)`: false conditions are dropped, the
/// arguments stop at the first true condition, and a single true piece is
/// its value.
pub(crate) fn raw_piecewise(pieces: Vec<(Expr, Expr)>) -> Expr {
    let mut out = Vec::with_capacity(pieces.len());
    for (e, c) in pieces {
        match c.kind() {
            Kind::Bool(false) => continue,
            Kind::Bool(true) => {
                out.push((e, c));
                break;
            }
            _ => out.push((e, c)),
        }
    }
    if out.len() == 1 && matches!(out[0].1.kind(), Kind::Bool(true)) {
        return out.pop().unwrap().0;
    }
    Expr::from_kind(Kind::Piecewise(out))
}

fn pieces(e: &Expr) -> Option<&[(Expr, Expr)]> {
    match e.kind() {
        Kind::Piecewise(p) => Some(p),
        _ => None,
    }
}

pub(crate) fn has_piecewise(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| found |= matches!(s.kind(), Kind::Piecewise(_)));
    found
}

/// `piecewise_fold(expr, evaluate)`: `Some(true)` evaluates the result,
/// `Some(false)` keeps the raw form, `None` drops duplicate conditions only.
pub(crate) fn piecewise_fold(expr: &Expr, evaluate: Option<bool>) -> Expr {
    if !has_piecewise(expr) {
        return expr.clone();
    }
    let mut new_args: Vec<(Expr, Expr)> = Vec::new();
    if let Some(ps) = pieces(expr) {
        for (e, c) in ps {
            let e = if pieces(e).is_some() { e.clone() } else { piecewise_fold(e, Some(true)) };
            match pieces(&e) {
                Some(inner) => new_args.extend(
                    inner.iter().map(|(ei, ci)| (piecewise_fold(ei, Some(true)), crate::logic::and(vec![ci.clone(), c.clone()]))),
                ),
                None => new_args.push((e, c.clone())),
            }
        }
    } else {
        let args: Vec<Expr> = if matches!(expr.kind(), Kind::Add(..) | Kind::Mul(..)) && crate::assumptions::fact(expr, "commutative") != Some(false) {
            group_common_conditions(expr)
        } else {
            crate::traverse::ordered_args(expr)
        };
        let folded: Vec<Expr> = args.iter().map(|a| piecewise_fold(a, Some(true))).collect();
        let mut combos: Vec<(Vec<Expr>, Vec<Expr>)> = vec![(Vec::new(), Vec::new())];
        for a in &folded {
            let choices: Vec<(Expr, Expr)> = match pieces(a) {
                Some(p) => p.to_vec(),
                None => vec![(a.clone(), Expr::boolean(true))],
            };
            let mut next = Vec::with_capacity(combos.len() * choices.len());
            for (vals, conds) in &combos {
                for (v, c) in &choices {
                    let mut vals = vals.clone();
                    vals.push(v.clone());
                    let mut conds = conds.clone();
                    conds.push(c.clone());
                    next.push((vals, conds));
                }
            }
            combos = next;
        }
        for (vals, conds) in combos {
            new_args.push((crate::simptools::rebuild(expr, &vals), crate::logic::and(conds)));
        }
    }
    if evaluate.is_none() {
        // keep the first value of each condition, in the order of the last
        // occurrences (the reference dictionary round trip)
        let mut keyed: Vec<(Expr, Expr)> = Vec::new();
        for (e, c) in new_args.iter().rev() {
            match keyed.iter_mut().find(|(k, _)| k == c) {
                Some(slot) => slot.1 = e.clone(),
                None => keyed.push((c.clone(), e.clone())),
            }
        }
        new_args = keyed.into_iter().rev().map(|(c, e)| (e, c)).collect();
    }
    let rv = if evaluate == Some(true) { crate::relational::piecewise(new_args) } else { raw_piecewise(new_args) };
    if evaluate.is_none() {
        if let Some([(e, c)]) = pieces(&rv) {
            if matches!(c.kind(), Kind::Bool(true)) {
                return e.clone();
            }
        }
    }
    let mut nested = false;
    crate::traverse::preorder(&rv, &mut |s| {
        if let Some(p) = pieces(s) {
            if p.iter().any(|(e, _)| has_piecewise(e)) {
                nested = true;
            }
        }
    });
    if nested {
        return piecewise_fold(&rv, Some(true));
    }
    rv
}

/// The grouping of `piecewise_fold` for a commutative Add or Mul:
/// Piecewise arguments with the same conditions are combined piece by
/// piece over their common condition prefix.
fn group_common_conditions(expr: &Expr) -> Vec<Expr> {
    let all = crate::traverse::ordered_args(expr);
    let (p, mut args): (Vec<Expr>, Vec<Expr>) = all.into_iter().partition(|a| pieces(a).is_some());
    // sift by the tuple of conditions, keys in `ordered` order
    let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
    for pw in p {
        let key = crate::special::tuple(pieces(&pw).unwrap().iter().map(|(_, c)| c.clone()).collect());
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some(g) => g.1.push(pw),
            None => groups.push((key, vec![pw])),
        }
    }
    let keys = crate::simptools::ordered(groups.iter().map(|(k, _)| k.clone()).collect());
    for key in keys {
        let members = groups.iter().find(|(k, _)| *k == key).map(|(_, m)| m.clone()).unwrap_or_default();
        if members.len() <= 1 {
            args.extend(members);
            continue;
        }
        let pargs: Vec<Vec<(Expr, Expr)>> = members.iter().map(|m| pieces(m).unwrap().to_vec()).collect();
        // the common prefix of the condition lists
        let shortest = pargs.iter().map(Vec::len).min().unwrap_or(0);
        let n = (0..shortest).take_while(|&i| pargs.iter().all(|a| a[i].1 == pargs[0][i].1)).count();
        let mut collected: Vec<(Expr, Expr)> = (0..n)
            .map(|i| {
                let vals: Vec<Expr> = pargs.iter().map(|a| a[i].0.clone()).collect();
                (crate::simptools::rebuild(expr, &vals), pargs[0][i].1.clone())
            })
            .collect();
        let mut remains = Vec::new();
        for a in &pargs {
            if n == a.len() {
                continue;
            }
            if matches!(a[n].1.kind(), Kind::Bool(true)) {
                remains.push(a[n].0.clone());
            } else {
                remains.push(raw_piecewise(a[n..].to_vec()));
            }
        }
        if !remains.is_empty() {
            collected.push((crate::simptools::rebuild(expr, &remains), Expr::boolean(true)));
        }
        args.push(raw_piecewise(collected));
    }
    args
}

// --------------------------------------------------------------------
// Integration of Piecewise functions
// --------------------------------------------------------------------

use crate::expr::{Constant, Func, RelOp};
use crate::logic::LogicOp;
use crate::ops::{add, neg, sub};
use crate::relational::{lt, rel};
use rustc_hash::FxHashMap;

/// The exceptions of the reference methods: `NotImplementedError` (an
/// undecided limit) and `ValueError` (an integral across an undefined
/// region). The callers present both as an unevaluated integral.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PwError {
    NotImplemented,
    Value,
}

fn is_true(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Bool(true))
}

fn is_false(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Bool(false))
}

fn has_x(e: &Expr, x: &Expr) -> bool {
    crate::traverse::free_symbols(e).contains(x)
}

/// `e.is_comparable`: a real number (the signed infinities included).
pub(crate) fn is_comparable(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity))
        || (crate::functions::is_numeric_value(e)
            && crate::assumptions::fact(e, "extended_real") != Some(false)
            && crate::evalf::eval_f64(e).is_some())
}

/// `e.has(*_illegal)`: nan or an infinity anywhere.
fn has_illegal(e: &Expr) -> bool {
    let mut bad = false;
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Const(c) = s.kind() {
            bad |= c.is_infinite() || matches!(c, Constant::NaN);
        }
    });
    bad
}

fn xrep(e: &Expr, pairs: &[(Expr, Expr)]) -> Expr {
    let map: FxHashMap<Expr, Expr> = pairs.iter().cloned().collect();
    crate::subs::xreplace(e, &map)
}

/// The relationals of `e`, each once, in traversal order.
fn relationals(e: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    crate::traverse::preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Rel(..)) && !out.contains(s) {
            out.push(s.clone());
        }
    });
    out
}

/// `classify(ie, s, i)` of `_solve_inequality`: the truth value of `ie` at
/// `s = i`, `None` when it does not evaluate (nan included).
fn classify(ie: &Expr, s: &Expr, i: &Expr) -> Option<bool> {
    match crate::subs::subs1(ie, s, i).kind() {
        Kind::Bool(b) => Some(*b),
        _ => None,
    }
}

/// `solvers.denoms`: the denominator factors of every subexpression.
fn denoms(e: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    crate::traverse::preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Rel(..) | Kind::Bool(_) | Kind::Func(Func::Logic(_), _)) {
            return;
        }
        let (_, d) = crate::simptools::radsimp::fraction(s, false);
        if d.is_one() {
            return;
        }
        for f in crate::simptools::mul_args(&d) {
            if !out.contains(&f) {
                out.push(f);
            }
        }
    });
    out
}

/// `_solve_inequality(ie, s, linear)`: the relational `ie` solved for `s`.
/// `Err` is the reference `NotImplementedError`.
pub(crate) fn solve_inequality(ie: &Expr, s: &Expr, linear: bool) -> Result<Expr, String> {
    let Kind::Rel(op, l, r) = ie.kind() else { return Ok(ie.clone()) };
    if !has_x(ie, s) {
        return Ok(ie.clone());
    }
    let (mut op, mut l, mut r) = (*op, l.clone(), r.clone());
    let mut ie = ie.clone();
    if r == *s {
        op = crate::logic::reversed_op(op);
        std::mem::swap(&mut l, &mut r);
        ie = rel(op, &l, &r);
    }
    if l == *s && !has_x(&r, s) {
        return Ok(ie);
    }
    let oo = Expr::oo();
    let noo = Expr::neg_oo();
    let expr = sub(&l, &r);
    let degree = crate::integrate::poly_coeffs(&expr, s).map(|c| {
        let n = c.iter().rposition(|t| !t.is_zero()).map(|i| i as i64);
        n.unwrap_or(-1)
    });
    let rv = (degree == Some(0)).then(|| rel(op, &crate::expand::expand(&expr), &Expr::zero()));
    let nonlinear = !linear && match degree {
        None => true,
        Some(d) => d > 1,
    };
    if nonlinear {
        let mut v = crate::sets::solve_univariate_inequality(&ie, s, true, &crate::sets::reals())?;
        // remove the restrictions at +-oo that the sets may have added
        let okoo = classify(&ie, s, &oo);
        if okoo == Some(true) && classify(&v, s, &oo) == Some(false) {
            v = crate::subs::subs1(&v, &lt(s, &oo), &Expr::boolean(true));
        }
        let oknoo = classify(&ie, s, &noo);
        if oknoo == Some(true) && classify(&v, s, &noo) == Some(false) {
            v = crate::subs::subs1(&v, &lt(&noo, s), &Expr::boolean(true));
            v = crate::subs::subs1(&v, &rel(RelOp::Gt, s, &noo), &Expr::boolean(true));
        }
        if is_true(&v) {
            v = if okoo == Some(true) { rel(RelOp::Le, s, &oo) } else { lt(s, &oo) };
            if oknoo != Some(true) {
                v = crate::logic::and(vec![lt(&noo, s), v]);
            }
        }
        return Ok(v);
    }
    let mut conds = Vec::new();
    let rv = match rv {
        Some(v) => v,
        None => {
            // a safe inversion: the terms free of `s` go to the right side,
            // and the division is by a factor of known sign (or nonzero
            // for Eq and Ne)
            let e = crate::expand::expand(&expr);
            let (b, _) = crate::decompose::as_independent(&e, std::slice::from_ref(s), Some(true));
            let e = sub(&e, &b);
            let rhs = neg(&b);
            let ef = crate::simptools::exprtools::factor_terms(&e, false, false, false, true);
            let (mut a, mut e) = crate::decompose::as_independent(&ef, std::slice::from_ref(s), Some(false));
            let zero_unknown = crate::assumptions::is_zero(&a) != Some(false);
            let sign_unknown = crate::assumptions::is_negative(&a).is_none() && crate::assumptions::is_positive(&a).is_none();
            if zero_unknown || (sign_unknown && !matches!(op, RelOp::Eq | RelOp::Ne)) {
                e = ef;
                a = Expr::one();
            }
            let rhs = crate::ops::div(&rhs, &a);
            let rv = if crate::assumptions::is_positive(&a) == Some(true) {
                rel(op, &e, &rhs)
            } else {
                rel(crate::logic::reversed_op(op), &e, &rhs)
            };
            // the conditions under which the value is valid
            let mut beginning = denoms(&l);
            for d in denoms(&r) {
                if !beginning.contains(&d) {
                    beginning.push(d);
                }
            }
            let current = denoms(&rv);
            for d in beginning.iter().filter(|d| !current.contains(d)) {
                let c = solve_inequality(&crate::relational::eq(d, &Expr::zero()), s, linear)?;
                if let Kind::Rel(RelOp::Eq, cl, cr) = c.kind() {
                    if cl == s && classify(&rv, s, cr) == Some(true) {
                        conds.push(crate::logic::not(&c));
                    }
                }
            }
            for (i, upper) in [(noo.clone(), false), (oo.clone(), true)] {
                if classify(&rv, s, &i) == Some(true) && classify(&ie, s, &i) != Some(true) {
                    conds.push(if upper { lt(s, &i) } else { lt(&i, s) });
                }
            }
            rv
        }
    };
    conds.push(rv);
    Ok(crate::logic::and(conds))
}

/// `(cond.lts, cond.gts)`: the lesser and the greater side.
fn lts_gts(op: RelOp, l: &Expr, r: &Expr) -> (Expr, Expr) {
    match op {
        RelOp::Gt | RelOp::Ge => (r.clone(), l.clone()),
        _ => (l.clone(), r.clone()),
    }
}

/// `_solve_relational` of `Piecewise._intervals`.
fn solve_relational(r: &Expr, x: &Expr) -> Result<Expr, String> {
    if !has_x(r, x) {
        return Err(format!("A condition not involving {x} appeared: {r}"));
    }
    let rv = solve_inequality(r, x, false).map_err(|_| format!("Unable to solve relational {r} for {x}."))?;
    if let Kind::Rel(op, l, rhs) = rv.kind() {
        if l != x || has_x(rhs, x) {
            return Err(format!("Unable to solve relational {r} for {x}."));
        }
        return Ok(match op {
            RelOp::Eq => Expr::boolean(false),
            RelOp::Ne => {
                if crate::assumptions::fact(rhs, "extended_real") == Some(false) || rhs.is_constant(Constant::NaN) {
                    Expr::boolean(true)
                } else {
                    crate::logic::or(vec![lt(x, rhs), rel(RelOp::Gt, x, rhs)])
                }
            }
            _ => rv.clone(),
        });
    }
    if rv == crate::logic::and(vec![lt(&Expr::neg_oo(), x), lt(x, &Expr::oo())]) {
        return Ok(Expr::boolean(true));
    }
    Ok(rv)
}

/// One `(a, b, e, i)` of `Piecewise._intervals`: the piece `e`, argument
/// `i`, is the value over `(a, b)`.
type Abei = (Expr, Expr, Expr, isize);

/// `Piecewise._intervals(x)`; `Err` carries the reference message.
fn intervals(pw: &Expr, x: &Expr) -> Result<Vec<Abei>, String> {
    let Some(args) = pieces(pw) else { return Err("not a Piecewise".into()) };
    let mut reps = Vec::new();
    for r in relationals(pw) {
        let s = solve_relational(&r, x)?;
        reps.push((r, s));
    }
    let args: Vec<(Expr, Expr)> = args.iter().map(|(e, c)| (xrep(e, &reps), xrep(c, &reps))).collect();
    let mut expr_cond: Vec<(isize, Expr, Expr)> = Vec::new();
    let mut default: Option<(Expr, isize)> = None;
    for (i, (expr, cond)) in args.iter().enumerate() {
        let i = i as isize;
        if is_false(cond) {
            continue;
        }
        if is_true(cond) {
            default = Some((expr.clone(), i));
            break;
        }
        if matches!(cond.kind(), Kind::Rel(RelOp::Eq, ..)) {
            continue; // zero width interval
        }
        let mut cond = crate::logic::to_cnf(cond, false, false)?;
        if matches!(cond.kind(), Kind::Func(Func::Logic(LogicOp::And), _)) {
            cond = crate::logic::distribute_or_over_and(&cond);
        }
        match cond.kind() {
            Kind::Func(Func::Logic(LogicOp::Or), os) => expr_cond.extend(
                os.iter().filter(|o| !matches!(o.kind(), Kind::Rel(RelOp::Eq, ..))).map(|o| (i, expr.clone(), o.clone())),
            ),
            Kind::Bool(false) => {}
            _ => expr_cond.push((i, expr.clone(), cond.clone())),
        }
    }
    let nonsym = |c: &Expr| format!("A condition not involving {x} appeared: {c}");
    let mut int_expr: Vec<Abei> = Vec::new();
    let mut k = 0;
    while k < expr_cond.len() {
        let (iarg, expr, cond) = expr_cond[k].clone();
        k += 1;
        let (lower, upper);
        match cond.kind() {
            Kind::Func(Func::Logic(LogicOp::And), cargs) => {
                let mut lo = Expr::neg_oo();
                let mut up = Expr::oo();
                let mut exclude = Vec::new();
                for cond2 in cargs {
                    let Kind::Rel(op, l, r) = cond2.kind() else { return Err("expecting only Relationals".into()) };
                    match op {
                        RelOp::Eq => {
                            lo = up.clone(); // ignore
                            break;
                        }
                        RelOp::Ne => {
                            if l == x {
                                exclude.push(r.clone());
                            } else if r == x {
                                exclude.push(l.clone());
                            } else {
                                return Err(nonsym(cond2));
                            }
                            continue;
                        }
                        _ => {}
                    }
                    let (lts, gts) = lts_gts(*op, l, r);
                    if lts == *x {
                        up = crate::functions::min(vec![gts, up]);
                    } else if gts == *x {
                        lo = crate::functions::max(vec![lts, lo]);
                    } else {
                        return Err(nonsym(cond2));
                    }
                }
                if !exclude.is_empty() {
                    let exclude = crate::simptools::ordered(exclude);
                    let mut newcond: Vec<(Option<Expr>, Expr)> = Vec::new();
                    for e in &exclude {
                        // `e < lower == True or e > upper == True` is a
                        // chained comparison that never holds, but the truth
                        // value of an undecided relational is an error
                        for t in [lt(e, &lo), rel(RelOp::Gt, e, &up)] {
                            if !matches!(t.kind(), Kind::Bool(_)) {
                                return Err("cannot determine truth value of Relational".into());
                            }
                        }
                        if newcond.is_empty() {
                            newcond.push((None, lo.clone()));
                        }
                        let last = newcond.last().unwrap().1.clone();
                        newcond.push((Some(last), e.clone()));
                    }
                    let last = newcond.last().unwrap().1.clone();
                    newcond.push((Some(last), up.clone()));
                    newcond.remove(0);
                    for (a, b) in newcond {
                        let a = a.expect("primer removed");
                        expr_cond.push((iarg, expr.clone(), crate::logic::and(vec![lt(&a, x), lt(x, &b)])));
                    }
                    continue;
                }
                lower = lo;
                upper = up;
            }
            Kind::Rel(op, l, r) if *op != RelOp::Ne => {
                let (lts, gts) = lts_gts(*op, l, r);
                if lts == *x {
                    lower = Expr::neg_oo();
                    upper = gts;
                } else if gts == *x {
                    lower = lts;
                    upper = Expr::oo();
                } else {
                    return Err(nonsym(&cond));
                }
            }
            _ => return Err(format!("unrecognized condition: {cond}")),
        }
        let upper = crate::functions::max(vec![lower.clone(), upper]);
        if !is_true(&crate::relational::ge(&lower, &upper)) {
            int_expr.push((lower, upper, expr, iarg));
        }
    }
    if let Some((e, i)) = default {
        int_expr.push((Expr::neg_oo(), Expr::oo(), e, i));
    }
    let mut out: Vec<Abei> = Vec::new();
    for t in int_expr {
        if !out.contains(&t) {
            out.push(t);
        }
    }
    Ok(out)
}

/// `_clip(A, B, k)`: the parts of `B` outside `A` keep the index -1, the
/// overlap receives `k`.
fn clip(a_: (&Expr, &Expr), b_: (&Expr, &Expr), k: isize) -> Vec<(Expr, Expr, isize)> {
    use crate::functions::{max, min};
    let (a, b) = (b_.0.clone(), b_.1.clone());
    let (c, d) = (a_.0, a_.1);
    let c = min(vec![max(vec![c.clone(), a.clone()]), b.clone()]);
    let d = min(vec![max(vec![d.clone(), a.clone()]), b.clone()]);
    let a = min(vec![a, b.clone()]);
    let mut p: Vec<(Expr, Expr, isize)> = Vec::new();
    if a != c {
        p.push((a.clone(), c.clone(), -1));
    }
    if c != d {
        p.push((c.clone(), d.clone(), k));
    }
    if b != d {
        if d == c && p.last().map(|t| t.2) == Some(-1) {
            let first = p.last().unwrap().0.clone();
            *p.last_mut().unwrap() = (first, b, -1);
        } else {
            p.push((d, b, -1));
        }
    }
    p
}

/// The covering loop shared by `_eval_integral` and `_eval_interval`:
/// every interval starts unassigned (-1) and the pieces claim their parts
/// in order.
fn cover(abei: &[Abei], start: (Expr, Expr)) -> Vec<(Expr, Expr, isize)> {
    let mut done: Vec<(Expr, Expr, isize)> = vec![(start.0, start.1, -1)];
    let (noo, oo) = (Expr::neg_oo(), Expr::oo());
    for (k, (pa, pb, _, _)) in abei.iter().enumerate() {
        let k = k as isize;
        if *pa == noo && *pb == oo {
            for t in done.iter_mut() {
                if t.2 == -1 {
                    t.2 = k;
                }
            }
            break; // nothing else to consider
        }
        let mut j = done.len();
        while j > 0 {
            j -= 1;
            if done[j].2 == -1 {
                let (a, b) = (done[j].0.clone(), done[j].1.clone());
                let parts = clip((pa, pb), (&a, &b), k);
                done.splice(j..j + 1, parts);
            }
        }
    }
    done.retain(|(a, b, _)| a != b);
    done
}

/// `Piecewise._handle_irel(x, handler)`: `None` when the conditions depend
/// on `x` only; otherwise the handled Piecewise for each truth assignment
/// of the relationals free of `x`, under the matching condition.
fn handle_irel(
    pw: &Expr,
    x: &Expr,
    handler: &mut dyn FnMut(&Expr) -> Result<Expr, PwError>,
) -> Result<Option<Expr>, PwError> {
    let irel = crate::simptools::ordered(relationals(pw).into_iter().filter(|r| !has_x(r, x)).collect());
    if irel.is_empty() {
        return Ok(None);
    }
    let n = irel.len();
    if n > 16 {
        // 2**n assignments, each integrated
        return Err(PwError::NotImplemented);
    }
    let mut args: Vec<(Expr, Vec<Expr>)> = Vec::new();
    let mut exprinorder: Vec<Expr> = Vec::new();
    let mut last = Expr::nan();
    for m in 0..(1usize << n) {
        // `product((1, 0), repeat=n)`: all ones first
        let truth: Vec<bool> = (0..n).map(|j| (m >> (n - 1 - j)) & 1 == 0).collect();
        let reps: Vec<(Expr, Expr)> = irel.iter().zip(&truth).map(|(r, t)| (r.clone(), Expr::boolean(*t))).collect();
        let mut cond = if !truth.contains(&true) {
            None // flag this one so it doesn't get combined
        } else {
            let andargs: Vec<Expr> = irel.iter().zip(&truth).filter(|(_, t)| **t).map(|(r, _)| r.clone()).collect();
            let mut free: Vec<Expr> = Vec::new();
            for a in &andargs {
                for s in crate::traverse::free_symbols(a) {
                    if !free.contains(&s) {
                        free.push(s);
                    }
                }
            }
            let t = if free.len() == 1 {
                match crate::sets::reduce_inequalities(&andargs, Some(&free[0])) {
                    Ok(t) => t,
                    Err(_) => crate::logic::and(
                        andargs
                            .iter()
                            .map(|a| solve_inequality(a, &free[0], true))
                            .collect::<Result<Vec<_>, _>>()
                            .map_err(|_| PwError::NotImplemented)?,
                    ),
                }
            } else {
                crate::logic::and(andargs)
            };
            if is_false(&t) {
                continue; // an impossible combination
            }
            Some(t)
        };
        let mut expr = handler(&xrep(pw, &reps))?;
        if let Some([(e, c)]) = pieces(&expr) {
            let (e, c) = (e.clone(), c.clone());
            cond = Some(crate::logic::and(vec![c, cond.unwrap_or_else(|| Expr::boolean(true))]));
            expr = e;
        }
        if let Some(c) = cond {
            match args.iter_mut().find(|(k, _)| *k == expr) {
                Some(slot) => slot.1.push(c),
                None => args.push((expr.clone(), vec![c])),
            }
            exprinorder.push(expr.clone());
        }
        last = expr;
    }
    let mut out: Vec<(Expr, Expr)> = Vec::new();
    let mut seen: Vec<Expr> = Vec::new();
    for e in exprinorder {
        if seen.contains(&e) {
            continue;
        }
        let conds = args.iter().find(|(k, _)| *k == e).map(|(_, c)| c.clone()).unwrap_or_default();
        out.push((e.clone(), crate::logic::or(conds)));
        seen.push(e);
    }
    out.push((last, Expr::boolean(true)));
    Ok(Some(crate::relational::piecewise(out)))
}

/// `integrate(e, x, meijerg=...)` for a piece.
fn integrate_piece(e: &Expr, x: &Expr, meijerg: Option<bool>) -> Expr {
    super::integrate_hint(e, x, meijerg)
}

/// `Piecewise.piecewise_integrate(x)`: each piece integrated on its own.
pub(crate) fn piecewise_integrate(pw: &Expr, x: &Expr, meijerg: Option<bool>) -> Expr {
    let Some(ps) = pieces(pw) else { return integrate_piece(pw, x, meijerg) };
    crate::relational::piecewise(ps.iter().map(|(e, c)| (integrate_piece(e, x, meijerg), c.clone())).collect())
}

/// `Piecewise._eval_integral(x)`: an antiderivative whose value at `x` is
/// the integral up to that point, continuous across the pieces.
pub(crate) fn eval_integral(pw: &Expr, x: &Expr, first: bool, meijerg: Option<bool>) -> Result<Expr, PwError> {
    if first {
        let mut handler = |ipw: &Expr| -> Result<Expr, PwError> {
            if pieces(ipw).is_some() {
                eval_integral(ipw, x, false, meijerg)
            } else {
                Ok(integrate_piece(ipw, x, meijerg))
            }
        };
        if let Some(irv) = handle_irel(pw, x, &mut handler)? {
            return Ok(irv);
        }
    }
    let Ok(mut abei) = intervals(pw, x) else {
        return Ok(super::integral_unevaluated(pw, x));
    };
    let done = cover(&abei, (Expr::neg_oo(), Expr::oo()));
    if done.iter().any(|t| t.2 == -1) {
        abei.push((Expr::neg_oo(), Expr::oo(), Expr::nan(), -1));
    }
    let orig = pieces(pw).expect("a Piecewise");
    let at = |i: isize| -> usize { if i < 0 { abei.len() - 1 } else { i as usize } };
    let mut out = Vec::new();
    let mut sum: Option<Expr> = None;
    for (a, b, i) in done {
        let piece = &abei[at(i)];
        let anti = integrate_piece(&piece.2, x, meijerg);
        sum = Some(match sum {
            None => anti.clone(),
            Some(s) => {
                let s = crate::subs::subs1(&s, x, &a);
                let e = eval_interval_any(&anti, x, &a, x)?;
                if has_illegal(&s) || has_illegal(&e) {
                    anti.clone()
                } else {
                    add(&[s, e])
                }
            }
        });
        let cond = if b == Expr::oo() {
            Expr::boolean(true)
        } else {
            let idx = if piece.3 < 0 { orig.len() - 1 } else { piece.3 as usize };
            if is_false(&crate::subs::subs1(&orig[idx].1, x, &b)) {
                lt(x, &b)
            } else {
                crate::relational::le(x, &b)
            }
        };
        out.push((sum.clone().unwrap(), cond));
    }
    Ok(crate::relational::piecewise(out))
}

/// `_eval_interval` of any expression: the Piecewise method or the
/// generic one.
pub(crate) fn eval_interval_any(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Result<Expr, PwError> {
    if pieces(f).is_some() {
        eval_interval(f, x, a, b, true)
    } else {
        expr_eval_interval(f, x, a, b)
    }
}

/// `Expr._eval_interval(x, a, b)`: the difference of the endpoint values,
/// through one-sided limits at a singular endpoint, corrected for the jumps
/// at the poles and logarithmic branch points inside comparable bounds.
pub(crate) fn expr_eval_interval(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Result<Expr, PwError> {
    use crate::limits::Dir;
    if a == b {
        return Ok(Expr::zero());
    }
    let increasing = !is_false(&lt(a, b));
    let endpoint = |left: bool| -> Result<Expr, PwError> {
        let c = if left { a } else { b };
        let value = crate::subs::subs1(f, x, c);
        if !has_illegal(&value) {
            return Ok(value);
        }
        let dir = if left == increasing { Dir::Plus } else { Dir::Minus };
        let value = crate::limits::limit(f, x, c, dir).map_err(|_| PwError::NotImplemented)?;
        // a limit left with an infinity inside (exp(-oo*s), oo**(n + 1))
        // is not decided for every parameter value: never a value
        if value.is_constant(Constant::NaN) || super::settled_bound_value(&value) && !super::has_bad(&value) {
            Ok(value)
        } else {
            Err(PwError::NotImplemented)
        }
    };
    let big_a = endpoint(true)?;
    if big_a.is_constant(Constant::NaN) {
        return Ok(big_a);
    }
    let big_b = endpoint(false)?;
    let value = sub(&big_b, &big_a);
    if is_comparable(a) && is_comparable(b) {
        return super::interval_singularities(f, x, a, b, value).ok_or(PwError::NotImplemented);
    }
    Ok(value)
}

/// `Piecewise._eval_interval(x, lo, hi, _first)`.
pub(crate) fn eval_interval(pw: &Expr, x: &Expr, lo: &Expr, hi: &Expr, first: bool) -> Result<Expr, PwError> {
    let (oo, noo) = (Expr::oo(), Expr::neg_oo());
    if first {
        let mut handler = |ipw: &Expr| -> Result<Expr, PwError> {
            if pieces(ipw).is_some() {
                eval_interval(ipw, x, lo, hi, false)
            } else {
                expr_eval_interval(ipw, x, lo, hi)
            }
        };
        if let Some(irv) = handle_irel(pw, x, &mut handler)? {
            return Ok(irv);
        }
        let lo_hi = lt(lo, hi);
        if is_false(&lo_hi) || *lo == oo || *hi == noo {
            let rv = eval_interval(pw, x, hi, lo, false)?;
            return Ok(match pieces(&rv) {
                Some(ps) => crate::relational::piecewise(ps.iter().map(|(e, c)| (neg(e), c.clone())).collect()),
                None => neg(&rv),
            });
        }
        if !(is_true(&lo_hi) || *hi == oo || *lo == noo) {
            let da = Expr::dummy("lo", crate::symbol::Assumptions::NONE);
            let db = Expr::dummy("hi", crate::symbol::Assumptions::NONE);
            let a = if is_comparable(lo) { lo.clone() } else { da.clone() };
            let b = if is_comparable(hi) { hi.clone() } else { db.clone() };
            let pos = eval_interval(pw, x, &a, &b, false)?;
            let (neg_, pos) = if a == da && b == db {
                (
                    neg(&xrep(&pos, &[(da.clone(), hi.clone()), (db.clone(), lo.clone())])),
                    xrep(&pos, &[(da.clone(), lo.clone()), (db.clone(), hi.clone())]),
                )
            } else {
                (
                    neg(&eval_interval(pw, x, hi, lo, false)?),
                    xrep(&pos, &[(da.clone(), lo.clone()), (db.clone(), hi.clone())]),
                )
            };
            let p = Expr::dummy("", crate::symbol::Assumptions::POSITIVE);
            let (mut pos, mut neg_) = (pos, neg_);
            if lo.is_symbol() {
                pos = xrep(&xrep(&pos, &[(lo.clone(), sub(hi, &p))]), &[(p.clone(), sub(hi, lo))]);
                neg_ = xrep(&xrep(&neg_, &[(lo.clone(), add(&[hi.clone(), p.clone()]))]), &[(p.clone(), sub(lo, hi))]);
            } else if hi.is_symbol() {
                pos = xrep(&xrep(&pos, &[(hi.clone(), add(&[lo.clone(), p.clone()]))]), &[(p.clone(), sub(hi, lo))]);
                neg_ = xrep(&xrep(&neg_, &[(hi.clone(), sub(lo, &p))]), &[(p.clone(), sub(lo, hi))]);
            }
            let neg_ = touch_min_max(&neg_);
            let pos = touch_min_max(&pos);
            let rv = if a == da {
                crate::relational::piecewise(vec![(pos.clone(), lt(lo, hi)), (neg_.clone(), Expr::boolean(true))])
            } else {
                crate::relational::piecewise(vec![(neg_.clone(), lt(hi, lo)), (pos.clone(), Expr::boolean(true))])
            };
            if rv.is_constant(Constant::NaN) {
                return Err(PwError::Value); // can't integrate across an undefined region
            }
            if pieces(&pos).is_some() || pieces(&neg_).is_some() {
                return Ok(piecewise_fold(&rv, Some(true)));
            }
            return Ok(rv);
        }
    }
    let Ok(abei) = intervals(pw, x) else {
        return Ok(super::integral_definite_unevaluated(&crate::diff::diff(pw, x), x, lo, hi));
    };
    let done = cover(&abei, (lo.clone(), hi.clone()));
    let mut sum = Expr::zero();
    let mut upto: Option<Expr> = None;
    for (a, b, i) in done {
        if i == -1 {
            return Ok(match upto {
                None => Expr::nan(),
                Some(u) => crate::relational::piecewise(vec![
                    (sum, crate::relational::le(hi, &u)),
                    (Expr::nan(), Expr::boolean(true)),
                ]),
            });
        }
        let e = &abei[i as usize].2;
        sum = add(&[sum, eval_interval_any(e, x, &a, &b)?]);
        upto = Some(b);
    }
    Ok(sum)
}

/// `touch`: every Min and Max rebuilt so that it evaluates again.
fn touch_min_max(e: &Expr) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    let e = crate::traverse::map_children(e, &mut |c| touch_min_max(c));
    match e.kind() {
        Kind::Func(Func::Min, args) => crate::functions::min(args.clone()),
        Kind::Func(Func::Max, args) => crate::functions::max(args.clone()),
        _ => e,
    }
}

fn has_func(e: &Expr, heads: &[Func]) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(f, _) = s.kind() {
            found |= heads.contains(f);
        }
    });
    found
}

/// The rewrites of `Integral.doit` before an antiderivative is sought:
/// `Abs` and `sign` become Piecewise when the variable (indefinite) or the
/// bounds (definite, finite) are extended real, `Min` and `Max` otherwise,
/// and a Piecewise inside the integrand is folded to the top.
pub(crate) fn doit_rewrite(f: &Expr, x: &Expr, bounds: Option<(&Expr, &Expr)>) -> Expr {
    use crate::assumptions::fact;
    let real = |e: &Expr| fact(e, "extended_real") == Some(true);
    let gate = match bounds {
        None => real(x),
        Some((a, b)) => [a, b].iter().all(|e| real(e) && fact(e, "infinite") != Some(true)),
    };
    let mut g = f.clone();
    if has_func(&g, &[Func::Abs, Func::Sign]) && gate {
        let xr = Expr::dummy("xr", crate::symbol::Assumptions::REAL);
        let r = crate::special::rewrite::rewrite(&xrep(&g, &[(x.clone(), xr.clone())]), &[], "Piecewise", true);
        g = xrep(&r, &[(xr, x.clone())]);
    } else if has_func(&g, &[Func::Min, Func::Max]) {
        g = crate::special::rewrite::rewrite(&g, &[], "Piecewise", true);
    }
    if has_piecewise(&g) && pieces(&g).is_none() {
        g = piecewise_fold(&g, Some(true));
    }
    g
}

/// Whether `e` is a Piecewise.
pub(crate) fn is_piecewise(e: &Expr) -> bool {
    pieces(e).is_some()
}

/// The definite evaluation of `Integral.doit` for an antiderivative with
/// Piecewise terms: those terms are folded together and evaluated by
/// `Piecewise._eval_interval`, the others by `others` (the generic
/// evaluation). `None` when an evaluation is undecided, or when a piece is
/// left as an unevaluated integral.
pub(crate) fn eval_split(
    big_f: &Expr,
    x: &Expr,
    a: &Expr,
    b: &Expr,
    others: &dyn Fn(&Expr) -> Option<Expr>,
) -> Option<Expr> {
    if has_func(big_f, &[Func::Integral]) {
        return None;
    }
    let (mut rest, mut pws) = (Vec::new(), Vec::new());
    for t in crate::simptools::add_args(big_f) {
        if crate::simptools::mul_args(&t).iter().any(|g| pieces(g).is_some()) {
            pws.push(piecewise_fold(&t, Some(true)));
        } else {
            rest.push(t);
        }
    }
    let evalued = if rest.is_empty() { Expr::zero() } else { others(&add(&rest))? };
    let evalued_pw = if pws.is_empty() {
        Expr::zero()
    } else {
        eval_interval_any(&piecewise_fold(&add(&pws), Some(true)), x, a, b).ok()?
    };
    Some(add(&[evalued, evalued_pw]))
}
