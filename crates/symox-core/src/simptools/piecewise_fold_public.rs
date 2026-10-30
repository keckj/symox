//! The public `piecewise_fold(expr, evaluate)` and the `piecewise_fold` of
//! the `Sum`/`Product`/`Integral` constructors (`ExprWithLimits._common_new`),
//! step for step as the Python API code of symox they replace
//! (`_public_calculus.piecewise_fold`, `_fold_piecewise_function`).
//!
//! The integration code keeps its own `integrate::piecewise::piecewise_fold`
//! (it differs on bound conditions, ITE conditions and non-commutative
//! products). Heads other than Add, Mul and Pow are rebuilt by the caller
//! (`expr.func(*args)` of the Python API), which also converts ITE
//! conditions; every other step is computed here.

use std::cmp::Ordering;

use crate::expr::Func;
use crate::logic::LogicOp;
use crate::{Expr, Kind};
use rustc_hash::FxHashMap;

/// The operations of the public API that stay with the caller.
pub trait FoldHost {
    type Error;
    /// `expr.func(*args)` for a head other than Add, Mul and Pow.
    fn rebuild(&mut self, expr: &Expr, args: Vec<Expr>) -> Result<Expr, Self::Error>;
    /// `simplify_logic(c.to_nnf(), form='cnf')` of an ITE condition.
    fn ite_condition(&mut self, condition: &Expr) -> Result<Expr, Self::Error>;
    /// An input this kernel does not reproduce exactly: the caller runs the
    /// Python algorithm.
    fn decline(&mut self) -> Self::Error;
    /// A fresh `Dummy()`.
    fn dummy(&mut self) -> Expr;
}

/// Whether the kernel reproduces the Python code on `expr`: no host
/// object, no matrix, and no Piecewise whose condition holds a Piecewise
/// (`ExprCondPair.__new__` folds those conditions in Python).
pub fn admitted(expr: &Expr) -> bool {
    if expr.has_host() || expr.has_matrix() {
        return false;
    }
    piecewise_atoms(expr)
        .iter()
        .all(|p| pieces(p).unwrap_or(&[]).iter().all(|(_, c)| !contains_piecewise(c)))
}

/// `expr.has(Piecewise)` / `_has_piecewise(expr)`.
pub fn contains_piecewise(expr: &Expr) -> bool {
    // Visited nodes stay alive (see `Expr::dag_id`).
    let mut seen = rustc_hash::FxHashMap::default();
    let mut pending = vec![expr.clone()];
    while let Some(node) = pending.pop() {
        if matches!(node.kind(), Kind::Piecewise(_)) {
            return true;
        }
        if seen.insert(node.dag_id(), node.clone()).is_some() {
            continue;
        }
        pending.extend(crate::inspection::scalar_children(&node));
    }
    false
}

/// `expr.atoms(Piecewise)`, without repetition.
fn piecewise_atoms(expr: &Expr) -> Vec<Expr> {
    crate::inspection::scalar_nodes(expr)
        .into_iter()
        .filter_map(|(node, _)| matches!(node.kind(), Kind::Piecewise(_)).then_some(node))
        .collect()
}

/// `expr.args` as the Python API presents them.
fn public_args(e: &Expr) -> Vec<Expr> {
    if let Kind::Func(Func::Heaviside, args) = e.kind() {
        if args.len() == 1 {
            return vec![args[0].clone(), Expr::rational(1, 2)];
        }
    }
    crate::traverse::ordered_args(e)
}

fn pieces(e: &Expr) -> Option<&[(Expr, Expr)]> {
    match e.kind() {
        Kind::Piecewise(p) => Some(p),
        _ => None,
    }
}

fn is_true(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Bool(true))
}

/// `expr.func(*args)`: Add, Mul and Pow are the native evaluated
/// constructors of `Add(*args)`, `Mul(*args)`, `Pow(b, e)`.
fn rebuild<H: FoldHost>(host: &mut H, expr: &Expr, args: Vec<Expr>) -> Result<Expr, H::Error> {
    match expr.kind() {
        Kind::Add(..) => Ok(crate::ops::add(&args)),
        Kind::Mul(..) => Ok(crate::ops::mul(&args)),
        Kind::Pow(..) if args.len() == 2 => Ok(crate::ops::pow(&args[0], &args[1])),
        _ => {
            let rebuilt = host.rebuild(expr, args)?;
            if !admitted(&rebuilt) {
                return Err(host.decline());
            }
            Ok(rebuilt)
        }
    }
}

/// `And(*conditions)`.
fn and<H: FoldHost>(host: &mut H, conditions: Vec<Expr>) -> Result<Expr, H::Error> {
    if conditions.iter().any(|c| !crate::logic::is_boolean(c) || matches!(c.kind(), Kind::Num(_))) {
        return Err(host.decline());
    }
    Ok(crate::logic::and(conditions))
}

/// `Piecewise(*pieces, evaluate=evaluate)` of the public constructor: the
/// native `Piecewise` when evaluated; held, the native held Piecewise with
/// the false branches dropped, stopped at the first true branch, and a
/// single true branch replaced by its value.
fn piecewise<H: FoldHost>(host: &mut H, pieces: Vec<(Expr, Expr)>, evaluate: bool) -> Result<Expr, H::Error> {
    let mut converted = Vec::with_capacity(pieces.len());
    for (e, c) in pieces {
        let c = crate::logic::boolean_piecewise_to_ite(&c).unwrap_or(c);
        if !crate::logic::is_settled_condition(&c) || e.has_matrix() {
            return Err(host.decline());
        }
        converted.push((e, c));
    }
    if evaluate {
        return Ok(crate::relational::piecewise(converted));
    }
    let mut branches = Vec::with_capacity(converted.len());
    for (e, c) in converted {
        if matches!(c.kind(), Kind::Bool(false)) {
            continue;
        }
        let last = is_true(&c);
        branches.push((e, c));
        if last {
            break;
        }
    }
    if branches.len() == 1 && is_true(&branches[0].1) {
        return Ok(branches.pop().unwrap().0);
    }
    Ok(Expr::from_unevaluated_kind(Kind::Piecewise(branches)))
}

/// `_node_count` of the Python API (a Float counts one half), doubled.
fn node_count2(e: &Expr) -> u64 {
    if let Kind::Num(n) = e.kind() {
        if n.is_float() {
            return 1;
        }
    }
    if let Some(p) = pieces(e) {
        // ExprCondPair arguments
        return 2 + p.iter().map(|(v, c)| 2 + node_count2(v) + node_count2(c)).sum::<u64>();
    }
    2 + public_args(e).iter().map(node_count2).sum::<u64>()
}

/// `ordered(keys)` for Python tuples of conditions: by `_nodes`, then by
/// `default_sort_key` (length, then the items), stable on ties.
fn ordered_condition_tuples(mut keys: Vec<Vec<Expr>>) -> Vec<Vec<Expr>> {
    let nodes = |k: &Vec<Expr>| 2 + k.iter().map(node_count2).sum::<u64>();
    keys.sort_by(|a, b| {
        nodes(a).cmp(&nodes(b)).then_with(|| {
            a.len().cmp(&b.len()).then_with(|| {
                a.iter()
                    .zip(b.iter())
                    .map(|(x, y)| crate::sortkey::cmp_default(x, y))
                    .find(|o| *o != Ordering::Equal)
                    .unwrap_or(Ordering::Equal)
            })
        })
    });
    keys
}

/// `_bound_piecewise_variables(expr)`: `None` unless `expr` is a Sum,
/// Product or Integral, else whether a Piecewise condition of its function
/// depends on one of its bound variables.
fn bound_piecewise_variables(expr: &Expr) -> Option<bool> {
    let Kind::Func(Func::Integral | Func::Sum | Func::Product, _) = expr.kind() else {
        return None;
    };
    let args = public_args(expr);
    let Some(function) = args.first() else { return Some(false) };
    let mut variables = Vec::new();
    for limit in &args[1..] {
        let items = public_args(limit);
        variables.push(items.first().cloned().unwrap_or_else(|| limit.clone()));
    }
    for piece in piecewise_atoms(function) {
        for (_, condition) in pieces(&piece).unwrap_or(&[]) {
            let free = crate::traverse::free_symbols(condition);
            if variables.iter().any(|v| free.contains(v)) {
                return Some(true);
            }
        }
    }
    Some(false)
}

/// The public `piecewise_fold(expr, evaluate)`; `evaluate` is `None`,
/// `Some(true)` or `Some(false)` (the truth value of the argument).
pub fn piecewise_fold<H: FoldHost>(host: &mut H, expr: &Expr, evaluate: Option<bool>) -> Result<Expr, H::Error> {
    if !contains_piecewise(expr) {
        return Ok(expr.clone());
    }
    let held = evaluate.unwrap_or(false);
    let mut new_args: Vec<(Expr, Expr)> = Vec::new();
    if let Some(ps) = pieces(expr) {
        for (e, c) in ps {
            let e = if pieces(e).is_some() { e.clone() } else { piecewise_fold(host, e, Some(true))? };
            let c = if matches!(c.kind(), Kind::Func(Func::Logic(LogicOp::ITE), _)) {
                host.ite_condition(c)?
            } else {
                c.clone()
            };
            match pieces(&e) {
                Some(inner) => {
                    for (ei, ci) in inner.to_vec() {
                        let value = piecewise_fold(host, &ei, Some(true))?;
                        let condition = and(host, vec![ci, c.clone()])?;
                        new_args.push((value, condition));
                    }
                }
                None => new_args.push((e, c)),
            }
        }
    } else {
        if bound_piecewise_variables(expr) == Some(true) {
            // A condition on a summation or integration variable cannot
            // leave its binder: fold the function only.
            let mut args = public_args(expr);
            args[0] = piecewise_fold(host, &args[0], evaluate)?;
            return rebuild(host, expr, args);
        }
        let commutative_product = matches!(expr.kind(), Kind::Mul(..)) && crate::decompose::commutativity(expr) == Some(true);
        let args = if matches!(expr.kind(), Kind::Add(..)) || commutative_product {
            group_common_conditions(host, expr)?
        } else {
            public_args(expr)
        };
        let mut folded = Vec::with_capacity(args.len());
        for a in &args {
            folded.push(piecewise_fold(host, a, Some(true))?);
        }
        let choices: Vec<Vec<(Expr, Expr)>> = folded
            .iter()
            .map(|a| match pieces(a) {
                Some(p) => p.to_vec(),
                None => vec![(a.clone(), Expr::boolean(true))],
            })
            .collect();
        // itertools.product order: the last argument varies fastest
        let total: usize = choices.iter().map(Vec::len).product();
        let mut index = vec![0usize; choices.len()];
        for _ in 0..total {
            let values: Vec<Expr> = index.iter().zip(&choices).map(|(&i, c)| c[i].0.clone()).collect();
            let conditions: Vec<Expr> = index.iter().zip(&choices).map(|(&i, c)| c[i].1.clone()).collect();
            let value = rebuild(host, expr, values)?;
            let condition = and(host, conditions)?;
            new_args.push((value, condition));
            for k in (0..index.len()).rev() {
                index[k] += 1;
                if index[k] < choices[k].len() {
                    break;
                }
                index[k] = 0;
            }
        }
    }
    if evaluate.is_none() {
        // {c: e for e, c in reversed(new_args)}, items reversed
        let mut keyed: Vec<(Expr, Expr)> = Vec::new();
        for (e, c) in new_args.iter().rev() {
            match keyed.iter_mut().find(|(k, _)| k == c) {
                Some(slot) => slot.1 = e.clone(),
                None => keyed.push((c.clone(), e.clone())),
            }
        }
        new_args = keyed.into_iter().rev().map(|(c, e)| (e, c)).collect();
    }
    let rv = piecewise(host, new_args, held)?; // held is the evaluate flag
    if evaluate.is_none() {
        match pieces(&rv) {
            Some([(e, c)]) if is_true(c) => return Ok(e.clone()),
            Some(_) => {}
            // `_pairs(rv)[0]` of a one-argument result: the Python code decides
            None if public_args(&rv).len() == 1 => return Err(host.decline()),
            None => {}
        }
    }
    if pieces(&rv).is_none() {
        return Ok(rv);
    }
    let nested = piecewise_atoms(&rv)
        .iter()
        .any(|q| pieces(q).unwrap_or(&[]).iter().any(|(e, _)| contains_piecewise(e)));
    if nested {
        return piecewise_fold(host, &rv, Some(true));
    }
    Ok(rv)
}

/// The grouping of a commutative Add or Mul: Piecewise arguments with the
/// same tuple of conditions are combined piece by piece.
fn group_common_conditions<H: FoldHost>(host: &mut H, expr: &Expr) -> Result<Vec<Expr>, H::Error> {
    let all = public_args(expr);
    let (p, mut args): (Vec<Expr>, Vec<Expr>) = all.into_iter().partition(|a| pieces(a).is_some());
    let mut groups: Vec<(Vec<Expr>, Vec<Expr>)> = Vec::new();
    for pw in p {
        let key: Vec<Expr> = pieces(&pw).unwrap().iter().map(|(_, c)| c.clone()).collect();
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some(g) => g.1.push(pw),
            None => groups.push((key, vec![pw])),
        }
    }
    let keys = ordered_condition_tuples(groups.iter().map(|(k, _)| k.clone()).collect());
    for key in keys {
        let members = groups.iter().find(|(k, _)| *k == key).map(|(_, m)| m.clone()).unwrap_or_default();
        if members.len() <= 1 {
            args.extend(members);
            continue;
        }
        let pargs: Vec<Vec<(Expr, Expr)>> = members.iter().map(|m| pieces(m).unwrap().to_vec()).collect();
        let shortest = pargs.iter().map(Vec::len).min().unwrap_or(0);
        let n = (0..shortest).take_while(|&i| pargs.iter().all(|a| a[i].1 == pargs[0][i].1)).count();
        let mut collected = Vec::with_capacity(n + 1);
        for i in 0..n {
            let values: Vec<Expr> = pargs.iter().map(|a| a[i].0.clone()).collect();
            collected.push((rebuild(host, expr, values)?, pargs[0][i].1.clone()));
        }
        let mut remains = Vec::new();
        for a in &pargs {
            if n == a.len() {
                continue;
            }
            if is_true(&a[n].1) {
                remains.push(a[n].0.clone());
            } else {
                remains.push(piecewise(host, a[n..].to_vec(), false)?);
            }
        }
        if !remains.is_empty() {
            collected.push((rebuild(host, expr, remains)?, Expr::boolean(true)));
        }
        args.push(piecewise(host, collected, false)?);
    }
    Ok(args)
}

/// `_fold_piecewise_function(function, limits)` of the Sum, Product and
/// Integral constructors with the limit variables `variables` (the free
/// symbols of `function` when there are no limits): Piecewise atoms free of
/// the variables are protected by Dummies while `function` is folded.
pub fn fold_limits_function<H: FoldHost>(host: &mut H, function: &Expr, variables: &[Expr]) -> Result<Expr, H::Error> {
    if !contains_piecewise(function) {
        return Ok(function.clone());
    }
    let mut reps: FxHashMap<Expr, Expr> = FxHashMap::default();
    for p in piecewise_atoms(function) {
        if !variables.iter().any(|v| crate::traverse::has_structural(&p, v)) {
            let d = host.dummy();
            reps.insert(p, d);
        }
    }
    let protected = if reps.is_empty() { function.clone() } else { crate::subs::xreplace(function, &reps) };
    let folded = piecewise_fold(host, &protected, Some(true))?;
    if reps.is_empty() {
        return Ok(folded);
    }
    let back: FxHashMap<Expr, Expr> = reps.into_iter().map(|(k, v)| (v, k)).collect();
    Ok(crate::subs::xreplace(&folded, &back))
}
