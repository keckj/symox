//! Normal forms: negation normal form, CNF and DNF.

use rustc_hash::FxHashSet;

use super::{and, as_logic, not, or, LogicOp};
use crate::expr::{Expr, Func, Kind};

/// SymPy's `is_literal`: an atom, a relational between atoms, or the
/// negation of a literal.
pub fn is_literal(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Logic(LogicOp::Not), a) => is_literal(&a[0]),
        Kind::Func(Func::Logic(_), _) => false,
        _ if e.is_atom() => true,
        _ => e.args().iter().all(|a| a.is_atom()),
    }
}

/// Is `e` in negation normal form? With `simplified`, also require that no
/// `And`/`Or` contains a literal together with its negation.
pub fn is_nnf(e: &Expr, simplified: bool) -> bool {
    if is_literal(e) {
        return true;
    }
    let mut stack = vec![e.clone()];
    while let Some(x) = stack.pop() {
        if is_literal(&x) {
            continue;
        }
        match as_logic(&x) {
            Some((LogicOp::And | LogicOp::Or, args)) => {
                if simplified {
                    let set: FxHashSet<&Expr> = args.iter().collect();
                    if args.iter().any(|a| set.contains(&not(a))) {
                        return false;
                    }
                }
                stack.extend(args.iter().cloned());
            }
            _ => return false,
        }
    }
    true
}

fn is_form(e: &Expr, outer: LogicOp, inner: LogicOp) -> bool {
    let vals: Vec<Expr> = match as_logic(e) {
        Some((op, a)) if op == outer => a.to_vec(),
        _ => vec![e.clone()],
    };
    for lit in &vals {
        match as_logic(lit) {
            Some((op, a)) if op == inner => {
                if !a.iter().all(is_literal) {
                    return false;
                }
            }
            _ => {
                if !is_literal(lit) {
                    return false;
                }
            }
        }
    }
    true
}

/// Is `e` in conjunctive normal form?
pub fn is_cnf(e: &Expr) -> bool {
    is_form(e, LogicOp::And, LogicOp::Or)
}

/// Is `e` in disjunctive normal form?
pub fn is_dnf(e: &Expr) -> bool {
    is_form(e, LogicOp::Or, LogicOp::And)
}

/// `And._to_nnf` / `Or._to_nnf`: convert operands, optionally removing
/// complementary pairs.
fn lattice_nnf(op: LogicOp, args: Vec<Expr>, simplify: bool) -> Expr {
    let zero = op == LogicOp::Or; // the absorbing element
    let mut set: Vec<Expr> = Vec::with_capacity(args.len());
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    for a in args {
        let a = if is_literal(&a) { a } else { to_nnf(&a, simplify) };
        if simplify {
            let parts: Vec<Expr> = match as_logic(&a) {
                Some((o, sub)) if o == op => sub.to_vec(),
                _ => vec![a],
            };
            for p in parts {
                if seen.contains(&not(&p)) {
                    return Expr::boolean(zero);
                }
                if seen.insert(p.clone()) {
                    set.push(p);
                }
            }
        } else if seen.insert(a.clone()) {
            set.push(a);
        }
    }
    if op == LogicOp::And {
        and(set)
    } else {
        or(set)
    }
}

/// Call `f` on every `k`-subset of `0..n` (lexicographic order).
pub(crate) fn combinations(n: usize, k: usize, f: &mut dyn FnMut(&[usize])) {
    if k > n {
        return;
    }
    let mut idx: Vec<usize> = (0..k).collect();
    loop {
        f(&idx);
        let mut i = k as isize - 1;
        while i >= 0 && idx[i as usize] == n - k + i as usize {
            i -= 1;
        }
        if i < 0 {
            return;
        }
        let i = i as usize;
        idx[i] += 1;
        for t in i + 1..k {
            idx[t] = idx[t - 1] + 1;
        }
    }
}

/// Clauses `Or(...)` of `args` where the members of each subset of size
/// `k` (`k` running over `start, start+2, ...`) are negated.
fn xor_clauses(args: &[Expr], start: usize) -> Vec<Expr> {
    let n = args.len();
    let mut out = Vec::new();
    let mut k = start;
    while k <= n {
        combinations(n, k, &mut |idx| {
            let clause: Vec<Expr> = args
                .iter()
                .enumerate()
                .map(|(i, a)| if idx.contains(&i) { not(a) } else { a.clone() })
                .collect();
            out.push(or(clause));
        });
        k += 2;
    }
    out
}

/// Negation normal form (SymPy's `to_nnf`, `simplify=True` by default
/// there).
pub fn to_nnf(e: &Expr, simplify: bool) -> Expr {
    let Some((op, args)) = as_logic(e) else { return e.clone() };
    match op {
        LogicOp::And | LogicOp::Or => lattice_nnf(op, args.to_vec(), simplify),
        LogicOp::Not => nnf_not(&args[0], simplify),
        LogicOp::Xor => lattice_nnf(LogicOp::And, xor_clauses(args, 0), simplify),
        LogicOp::Implies => lattice_nnf(LogicOp::Or, vec![not(&args[0]), args[1].clone()], simplify),
        LogicOp::Equivalent => {
            let n = args.len();
            let mut v = Vec::with_capacity(n);
            for i in 0..n - 1 {
                v.push(or(vec![not(&args[i]), args[i + 1].clone()]));
            }
            v.push(or(vec![not(&args[n - 1]), args[0].clone()]));
            lattice_nnf(LogicOp::And, v, simplify)
        }
        LogicOp::ITE => {
            let (a, b, c) = (&args[0], &args[1], &args[2]);
            lattice_nnf(
                LogicOp::And,
                vec![or(vec![not(a), b.clone()]), or(vec![a.clone(), c.clone()])],
                simplify,
            )
        }
    }
}

fn nnf_not(x: &Expr, simplify: bool) -> Expr {
    let Some((op, args)) = as_logic(x) else { return not(x) };
    match op {
        LogicOp::Not => to_nnf(&args[0], simplify),
        LogicOp::And => lattice_nnf(LogicOp::Or, args.iter().map(not).collect(), simplify),
        LogicOp::Or => lattice_nnf(LogicOp::And, args.iter().map(not).collect(), simplify),
        LogicOp::Implies => lattice_nnf(LogicOp::And, vec![args[0].clone(), not(&args[1])], simplify),
        LogicOp::Equivalent => lattice_nnf(
            LogicOp::And,
            vec![or(args.to_vec()), or(args.iter().map(not).collect())],
            simplify,
        ),
        LogicOp::Xor => lattice_nnf(LogicOp::And, xor_clauses(args, 1), simplify),
        LogicOp::ITE => {
            let (a, b, c) = (&args[0], &args[1], &args[2]);
            lattice_nnf(
                LogicOp::And,
                vec![or(vec![a.clone(), not(c)]), or(vec![not(a), not(b)])],
                simplify,
            )
        }
    }
}

/// Replace `Implies`/`Equivalent`/`Xor`/`ITE` by `And`/`Or`/`Not`.
pub fn eliminate_implications(e: &Expr) -> Expr {
    to_nnf(e, false)
}

fn distribute(e: &Expr, outer: LogicOp, inner: LogicOp) -> Expr {
    let build = |op: LogicOp, v: Vec<Expr>| if op == LogicOp::And { and(v) } else { or(v) };
    match as_logic(e) {
        Some((op, args)) if op == inner => {
            let Some(pos) = args.iter().position(|a| matches!(as_logic(a), Some((o, _)) if o == outer)) else {
                return e.clone();
            };
            let conj = &args[pos];
            let rest: Vec<Expr> = args.iter().enumerate().filter(|(i, _)| *i != pos).map(|x| x.1.clone()).collect();
            let rest = build(inner, rest);
            let (_, cargs) = as_logic(conj).unwrap();
            let parts: Vec<Expr> = cargs
                .iter()
                .map(|c| distribute(&build(inner, vec![c.clone(), rest.clone()]), outer, inner))
                .collect();
            build(outer, parts)
        }
        Some((op, args)) if op == outer => {
            let parts: Vec<Expr> = args.iter().map(|a| distribute(a, outer, inner)).collect();
            build(outer, parts)
        }
        _ => e.clone(),
    }
}

/// Distribute `And` over `Or` (the CNF step).
pub fn distribute_and_over_or(e: &Expr) -> Expr {
    distribute(e, LogicOp::And, LogicOp::Or)
}

/// Distribute `Or` over `And` (the DNF step).
pub fn distribute_or_over_and(e: &Expr) -> Expr {
    distribute(e, LogicOp::Or, LogicOp::And)
}

fn is_boolean_function(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Logic(_), _))
}

/// Error for `simplify=True` on more than 8 variables without `force`.
pub const TOO_MANY_VARIABLES: &str =
    "To simplify expressions with more than 8 variables force=True must be used";

/// Conjunctive normal form. With `simplify`, a minimal CNF is computed
/// ([`super::simplify_logic`]); that is refused (`Err`) above 8 variables
/// unless `force`.
pub fn to_cnf(e: &Expr, simplify: bool, force: bool) -> Result<Expr, String> {
    if !is_boolean_function(e) {
        return Ok(e.clone());
    }
    if simplify {
        if super::atoms(e).len() > 20 {
            return Err("truth-table operations support at most 20 variables".into());
        }
        if !force && super::atoms(e).len() > 8 {
            return Err(TOO_MANY_VARIABLES.into());
        }
        return Ok(super::simplify_logic(e, Some(super::LogicForm::Cnf), true, force));
    }
    if is_cnf(e) {
        return Ok(e.clone());
    }
    Ok(distribute_and_over_or(&eliminate_implications(e)))
}

/// Disjunctive normal form (see [`to_cnf`]).
pub fn to_dnf(e: &Expr, simplify: bool, force: bool) -> Result<Expr, String> {
    if !is_boolean_function(e) {
        return Ok(e.clone());
    }
    if simplify {
        if super::atoms(e).len() > 20 {
            return Err("truth-table operations support at most 20 variables".into());
        }
        if !force && super::atoms(e).len() > 8 {
            return Err(TOO_MANY_VARIABLES.into());
        }
        return Ok(super::simplify_logic(e, Some(super::LogicForm::Dnf), true, force));
    }
    if is_dnf(e) {
        return Ok(e.clone());
    }
    Ok(distribute_or_over_and(&eliminate_implications(e)))
}

