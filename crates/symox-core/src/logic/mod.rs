//! Boolean algebra: `And`, `Or`, `Not`, `Xor`, `Implies`, `Equivalent`,
//! `ITE` with SymPy's automatic evaluation, normal forms ([`to_cnf`],
//! [`to_dnf`], [`to_nnf`]), two-level minimisation ([`simplify_logic`],
//! Quine–McCluskey with an exact cover), a CDCL SAT solver
//! ([`satisfiable`]) and truth tables.
//!
//! Boolean expressions are ordinary [`Expr`] nodes
//! `Kind::Func(Func::Logic(op), args)`. Commutative connectives keep their
//! arguments in SymPy's `ordered()` order (node count, then
//! `default_sort_key`), which is also the printing order.

mod minimize;
mod normal;
mod anf;
pub mod sat;

pub use minimize::{simplify_logic, sop_form, pos_form, truth_table, LogicForm};
pub use minimize::simplify_logic_dontcare;
pub use anf::{to_anf, is_anf, anf_form, anf_coeffs, anf_from_coeffs};
pub use normal::{distribute_and_over_or, distribute_or_over_and, eliminate_implications, is_cnf, is_dnf, is_literal, is_nnf, to_cnf, to_dnf, to_nnf};
pub use sat::{satisfiable, satisfiable_all, Model};

use std::cmp::Ordering;

use rustc_hash::FxHashSet;

use crate::expr::{Constant, Expr, Func, Kind, RelOp};
// Relational.canonical uses Expr.could_extract_minus_sign, whose ties are
// broken by sort_key.
use crate::decompose::could_extract_minus_sign;
use crate::ops::neg;
use crate::relational::rel;
use crate::sortkey::{sort_key, Key};

/// Boolean connectives.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum LogicOp {
    And,
    Or,
    Not,
    Xor,
    Implies,
    Equivalent,
    ITE,
}

impl LogicOp {
    pub fn name(self) -> &'static str {
        match self {
            LogicOp::And => "And",
            LogicOp::Or => "Or",
            LogicOp::Not => "Not",
            LogicOp::Xor => "Xor",
            LogicOp::Implies => "Implies",
            LogicOp::Equivalent => "Equivalent",
            LogicOp::ITE => "ITE",
        }
    }

    pub fn from_name(name: &str) -> Option<LogicOp> {
        Some(match name {
            "And" => LogicOp::And,
            "Or" => LogicOp::Or,
            "Not" => LogicOp::Not,
            "Xor" => LogicOp::Xor,
            "Implies" => LogicOp::Implies,
            "Equivalent" => LogicOp::Equivalent,
            "ITE" => LogicOp::ITE,
            _ => return None,
        })
    }
}

// ----------------------------------------------------------------------
// small helpers
// ----------------------------------------------------------------------

#[inline]
pub(crate) fn mk(op: LogicOp, args: Vec<Expr>) -> Expr {
    Expr::from_kind(Kind::Func(Func::Logic(op), args))
}

/// The connective and arguments of a boolean function node.
#[inline]
pub fn as_logic(e: &Expr) -> Option<(LogicOp, &[Expr])> {
    match e.kind() {
        Kind::Func(Func::Logic(op), a) => Some((*op, a.as_slice())),
        _ => None,
    }
}

/// Is `e` a boolean valued expression (SymPy's `Boolean`)? Symbols count
/// as boolean variables.
pub fn is_boolean(e: &Expr) -> bool {
    match e.kind() {
        Kind::Bool(_) | Kind::Sym(_) | Kind::Rel(..) | Kind::Func(Func::Logic(_) | Func::AppliedPredicate, _) => true,
        Kind::Func(Func::Set(crate::sets::SetOp::Contains), _) => true,
        Kind::Piecewise(p) => p.iter().all(|(x, _)| is_boolean(x)),
        _ => false,
    }
}

/// Twice SymPy's `_nodes` (floats count 1/2).
pub fn nodes2(e: &Expr) -> u64 {
    match e.kind() {
        Kind::Num(crate::number::Number::Float(_)) => 1,
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => 2,
        // `Lambda` stores its signature as a `Tuple`
        Kind::Func(Func::Set(crate::sets::SetOp::Lambda), a) => {
            2 + a.iter().map(nodes2).sum::<u64>()
        }
        _ => 2 + e.args().iter().map(nodes2).sum::<u64>(),
    }
}

/// SymPy's `ordered()` key: node count then `default_sort_key`.
pub fn ordered_key(e: &Expr) -> (u64, Key) {
    (nodes2(e), sort_key(e))
}

/// Cheap ordering keys for literals (`x`, `~x`) and clauses of literals,
/// equivalent to [`ordered_key`] on them.
#[derive(PartialEq, Eq)]
enum Light<'a> {
    /// (negated, name)
    Lit(bool, &'a str),
    /// (connective class rank, literals)
    Clause(u8, Vec<(bool, &'a str)>),
}

fn light_lit(e: &Expr) -> Option<(bool, &str)> {
    match e.kind() {
        Kind::Sym(s) if !s.is_wild() && !s.name().starts_with('_') => Some((false, s.name())),
        Kind::Func(Func::Logic(LogicOp::Not), a) => match a[0].kind() {
            Kind::Sym(s) if !s.is_wild() && !s.name().starts_with('_') => Some((true, s.name())),
            _ => None,
        },
        _ => None,
    }
}

fn light(e: &Expr) -> Option<Light<'_>> {
    if let Some((n, s)) = light_lit(e) {
        return Some(Light::Lit(n, s));
    }
    match e.kind() {
        Kind::Func(Func::Logic(op @ (LogicOp::And | LogicOp::Or)), a) => {
            let rank = if *op == LogicOp::And { 0 } else { 2 };
            let lits: Option<Vec<(bool, &str)>> = a.iter().map(light_lit).collect();
            Some(Light::Clause(rank, lits?))
        }
        _ => None,
    }
}

impl Light<'_> {
    fn nodes(&self) -> u64 {
        match self {
            Light::Lit(n, _) => 2 + 2 * (*n as u64),
            Light::Clause(_, l) => 2 + l.iter().map(|(n, _)| 2 + 2 * (*n as u64)).sum::<u64>(),
        }
    }
    /// class key rank: Symbol < And < Not < Or
    fn class(&self) -> u8 {
        match self {
            Light::Lit(false, _) => 0,
            Light::Clause(0, _) => 1,
            Light::Lit(true, _) => 2,
            _ => 3,
        }
    }
}

fn lit_cmp(a: &(bool, &str), b: &(bool, &str)) -> Ordering {
    a.0.cmp(&b.0).then_with(|| a.1.cmp(b.1))
}

fn light_cmp(a: &Light, b: &Light) -> Ordering {
    a.nodes().cmp(&b.nodes()).then_with(|| a.class().cmp(&b.class())).then_with(|| match (a, b) {
        (Light::Lit(_, x), Light::Lit(_, y)) => x.cmp(y),
        (Light::Clause(_, x), Light::Clause(_, y)) => x.len().cmp(&y.len()).then_with(|| {
            for (p, q) in x.iter().zip(y) {
                let r = lit_cmp(p, q);
                if r != Ordering::Equal {
                    return r;
                }
            }
            Ordering::Equal
        }),
        _ => Ordering::Equal,
    })
}

/// Sort like SymPy's `ordered()` (ties broken structurally).
pub fn ordered(v: &mut [Expr]) {
    if v.len() < 2 {
        return;
    }
    // fast path: literals and clauses of literals
    if let Some(mut keyed) = v.iter().map(|e| light(e).map(|l| (l, e.clone()))).collect::<Option<Vec<_>>>() {
        keyed.sort_by(|a, b| light_cmp(&a.0, &b.0).then_with(|| a.1.cmp_struct(&b.1)));
        let sorted: Vec<Expr> = keyed.into_iter().map(|x| x.1).collect();
        for (slot, e) in v.iter_mut().zip(sorted) {
            *slot = e;
        }
        return;
    }
    let mut keyed: Vec<((u64, Key), Expr)> = v.iter().map(|e| (ordered_key(e), e.clone())).collect();
    keyed.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp_struct(&b.1)));
    for (slot, (_, e)) in v.iter_mut().zip(keyed) {
        *slot = e;
    }
}

/// Compare two expressions like SymPy's `ordered()`.
pub fn ordered_cmp(a: &Expr, b: &Expr) -> Ordering {
    ordered_key(a).cmp(&ordered_key(b)).then_with(|| a.cmp_struct(b))
}

// ----------------------------------------------------------------------
// relationals
// ----------------------------------------------------------------------

/// `a op b` -> `b op' a`.
pub fn reversed_op(op: RelOp) -> RelOp {
    match op {
        RelOp::Lt => RelOp::Gt,
        RelOp::Gt => RelOp::Lt,
        RelOp::Le => RelOp::Ge,
        RelOp::Ge => RelOp::Le,
        o => o,
    }
}

/// Logical negation of a relational operator.
pub fn negated_op(op: RelOp) -> RelOp {
    match op {
        RelOp::Lt => RelOp::Ge,
        RelOp::Ge => RelOp::Lt,
        RelOp::Le => RelOp::Gt,
        RelOp::Gt => RelOp::Le,
        RelOp::Eq => RelOp::Ne,
        RelOp::Ne => RelOp::Eq,
    }
}

/// Free of symbols (SymPy's `is_number`, which includes `oo` and `I`).
pub(crate) fn is_number_like(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) | Kind::Const(_) => true,
        Kind::Sym(_) | Kind::Bool(_) | Kind::Rel(..) | Kind::Piecewise(_) | Kind::Derivative(..) => false,
        Kind::Func(Func::Predicate(_) | Func::AppliedPredicate | Func::Undefined(_) | Func::Logic(_) | Func::Set(_), _) => false,
        _ => e.mask() == 0 && e.args().iter().all(is_number_like),
    }
}

fn is_number_atom(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(_))
        || matches!(e.kind(), Kind::Const(c) if matches!(c, Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN))
}

fn num_gt(a: &Expr, b: &Expr) -> bool {
    match (a.kind(), b.kind()) {
        (Kind::Num(x), Kind::Num(y)) => x.cmp_value(y) == Some(Ordering::Greater),
        _ => {
            let f = |e: &Expr| match e.kind() {
                Kind::Num(n) => Some(n.to_f64()),
                Kind::Const(Constant::Infinity) => Some(f64::INFINITY),
                Kind::Const(Constant::NegInfinity) => Some(f64::NEG_INFINITY),
                _ => None,
            };
            matches!((f(a), f(b)), (Some(x), Some(y)) if x > y)
        }
    }
}

/// SymPy's `Relational.canonical`.
pub fn canonical(r: &Expr) -> Expr {
    let Kind::Rel(op, l, rr) = r.kind() else { return r.clone() };
    let (mut op, mut l, mut rr) = (*op, l.clone(), rr.clone());
    let rev = |op: &mut RelOp, l: &mut Expr, r: &mut Expr| {
        std::mem::swap(l, r);
        *op = reversed_op(*op);
    };
    if is_number_like(&rr) {
        if is_number_atom(&rr) && is_number_atom(&l) && num_gt(&l, &rr) {
            rev(&mut op, &mut l, &mut rr);
        }
    } else if is_number_like(&l) {
        rev(&mut op, &mut l, &mut rr);
    } else if ordered_cmp(&l, &rr) == Ordering::Greater {
        rev(&mut op, &mut l, &mut rr);
    }
    if matches!(l.kind(), Kind::Bool(_)) || matches!(rr.kind(), Kind::Bool(_)) {
        return rel(op, &l, &rr);
    }
    if could_extract_minus_sign(&l) {
        // reversedsign
        return rel(reversed_op(op), &neg(&l), &neg(&rr));
    } else if !is_number_like(&rr) && could_extract_minus_sign(&rr) {
        let nr = neg(&rr);
        if ordered_cmp(&nr, &l) == Ordering::Less {
            // reversed, then reversedsign: `a op b` -> `-b op -a`
            return rel(op, &nr, &neg(&l));
        }
    }
    let out = rel(op, &l, &rr);
    if out == *r {
        r.clone()
    } else {
        out
    }
}

/// Negation of a relational (`x < 1` -> `x >= 1`).
pub fn negated_rel(r: &Expr) -> Expr {
    match r.kind() {
        Kind::Rel(op, a, b) => rel(negated_op(*op), a, b),
        _ => not(r),
    }
}

// ----------------------------------------------------------------------
// evaluation
// ----------------------------------------------------------------------

/// Evaluate `op(args)` (entry point used by [`crate::functions::func`]).
pub fn eval(op: LogicOp, args: Vec<Expr>) -> Expr {
    // Substitution can introduce integer 0/1 after Python's boolean argument
    // conversion. Rebuilding the connective must apply the same identities.
    let args = args.into_iter().map(|a| {
        if a.is_integer_number() && (a.is_zero() || a.is_one()) {
            Expr::boolean(a.is_one())
        } else { a }
    }).collect();
    match op {
        LogicOp::And => lattice(true, args),
        LogicOp::Or => lattice(false, args),
        LogicOp::Not => {
            if args.len() == 1 {
                not(&args[0])
            } else {
                mk(LogicOp::Not, args)
            }
        }
        LogicOp::Xor => eval_xor(args),
        LogicOp::Implies => {
            if args.len() == 2 {
                implies(&args[0], &args[1])
            } else {
                mk(LogicOp::Implies, args)
            }
        }
        LogicOp::Equivalent => eval_equivalent(args),
        LogicOp::ITE => {
            if args.len() == 3 {
                ite(&args[0], &args[1], &args[2])
            } else {
                mk(LogicOp::ITE, args)
            }
        }
    }
}

/// `And` (`is_and`) / `Or` evaluation: flattening, constants, idempotence,
/// complementary relationals.
fn lattice(is_and: bool, args: Vec<Expr>) -> Expr {
    let (op, zero) = if is_and { (LogicOp::And, false) } else { (LogicOp::Or, true) };
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let mut flat: Vec<Expr> = Vec::with_capacity(args.len());
    let mut has_rel = false;
    let mut stack: Vec<Expr> = args;
    stack.reverse();
    while let Some(a) = stack.pop() {
        match a.kind() {
            Kind::Bool(b) => {
                if *b == zero {
                    return Expr::boolean(zero);
                }
            }
            Kind::Func(Func::Logic(o), sub) if *o == op => {
                for s in sub.iter().rev() {
                    stack.push(s.clone());
                }
            }
            _ => {
                if seen.insert(a.clone()) {
                    has_rel |= matches!(a.kind(), Kind::Rel(..));
                    flat.push(a);
                }
            }
        }
    }
    ordered(&mut flat);
    if has_rel {
        let mut rels: FxHashSet<Expr> = FxHashSet::default();
        let mut out = Vec::with_capacity(flat.len());
        for x in flat {
            if matches!(x.kind(), Kind::Rel(..)) {
                let c = canonical(&x);
                if rels.contains(&c) {
                    continue;
                }
                let nc = canonical(&negated_rel(&c));
                if rels.contains(&nc) {
                    return Expr::boolean(zero);
                }
                rels.insert(c);
            }
            out.push(x);
        }
        flat = out;
    }
    match flat.len() {
        0 => Expr::boolean(!zero),
        1 => flat.pop().unwrap(),
        _ => mk(op, flat),
    }
}

fn eval_xor(args: Vec<Expr>) -> Expr {
    // toggle set (keeps insertion order for determinism)
    let mut set: Vec<Expr> = Vec::new();
    let mut idx: FxHashSet<Expr> = FxHashSet::default();
    let mut truth = false;
    let toggle = |a: Expr, set: &mut Vec<Expr>, idx: &mut FxHashSet<Expr>| {
        if idx.remove(&a) {
            set.retain(|x| *x != a);
        } else {
            idx.insert(a.clone());
            set.push(a);
        }
    };
    for a in args {
        match a.kind() {
            Kind::Bool(b) => truth ^= *b,
            Kind::Func(Func::Logic(LogicOp::Xor), sub) => {
                for s in sub.clone() {
                    toggle(s, &mut set, &mut idx);
                }
            }
            _ => toggle(a, &mut set, &mut idx),
        }
    }
    // complementary / equal relationals
    let rels: Vec<(usize, Expr, Expr)> = set
        .iter()
        .enumerate()
        .filter(|(_, r)| matches!(r.kind(), Kind::Rel(..)))
        .map(|(i, r)| (i, canonical(r), canonical(&negated_rel(r))))
        .collect();
    let mut removed: FxHashSet<usize> = FxHashSet::default();
    for i in 0..rels.len() {
        if removed.contains(&rels[i].0) {
            continue;
        }
        for j in i + 1..rels.len() {
            if removed.contains(&rels[j].0) {
                continue;
            }
            if rels[j].1 == rels[i].2 {
                truth = !truth;
            } else if rels[j].1 != rels[i].1 {
                continue;
            }
            removed.insert(rels[i].0);
            removed.insert(rels[j].0);
            break;
        }
    }
    let mut rest: Vec<Expr> = set
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !removed.contains(i))
        .map(|x| x.1)
        .collect();
    if rest.is_empty() {
        return Expr::boolean(truth);
    }
    if truth {
        return not(&eval_xor(rest));
    }
    if rest.len() == 1 {
        return rest.pop().unwrap();
    }
    ordered(&mut rest);
    mk(LogicOp::Xor, rest)
}

fn eval_equivalent(args: Vec<Expr>) -> Expr {
    let mut set: Vec<Expr> = Vec::new();
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let (mut has_true, mut has_false) = (false, false);
    for a in args {
        match a.kind() {
            Kind::Bool(true) => has_true = true,
            Kind::Bool(false) => has_false = true,
            _ => {
                if seen.insert(a.clone()) {
                    set.push(a);
                }
            }
        }
    }
    let rels: Vec<(usize, Expr, Expr)> = set
        .iter()
        .enumerate()
        .filter(|(_, r)| matches!(r.kind(), Kind::Rel(..)))
        .map(|(i, r)| (i, canonical(r), canonical(&negated_rel(r))))
        .collect();
    let mut removed: FxHashSet<usize> = FxHashSet::default();
    for i in 0..rels.len() {
        for j in i + 1..rels.len() {
            if rels[j].1 == rels[i].2 {
                return Expr::boolean(false);
            } else if rels[j].1 == rels[i].1 {
                removed.insert(rels[i].0);
                removed.insert(rels[j].0);
                has_true = true;
                break;
            }
        }
    }
    let rest: Vec<Expr> = set
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !removed.contains(i))
        .map(|x| x.1)
        .collect();
    let n = rest.len() + has_true as usize + has_false as usize;
    if n <= 1 {
        return Expr::boolean(true);
    }
    if has_true {
        let mut v = rest;
        if has_false {
            v.push(Expr::boolean(false));
        }
        return and(v);
    }
    if has_false {
        return and(rest.iter().map(not).collect());
    }
    let mut rest = rest;
    ordered(&mut rest);
    mk(LogicOp::Equivalent, rest)
}

// ----------------------------------------------------------------------
// public constructors
// ----------------------------------------------------------------------

/// `And(*args)`.
pub fn and(args: Vec<Expr>) -> Expr {
    lattice(true, args)
}

/// `Or(*args)`.
pub fn or(args: Vec<Expr>) -> Expr {
    lattice(false, args)
}

/// `a & b`.
pub fn and2(a: &Expr, b: &Expr) -> Expr {
    lattice(true, vec![a.clone(), b.clone()])
}

/// `a | b`.
pub fn or2(a: &Expr, b: &Expr) -> Expr {
    lattice(false, vec![a.clone(), b.clone()])
}

/// `Not(a)` (`~a`).
pub fn not(a: &Expr) -> Expr {
    match a.kind() {
        Kind::Bool(b) => Expr::boolean(!b),
        Kind::Func(Func::Logic(LogicOp::Not), x) => x[0].clone(),
        Kind::Rel(op, l, r) => rel(negated_op(*op), l, r),
        _ => mk(LogicOp::Not, vec![a.clone()]),
    }
}

/// `Xor(*args)`.
pub fn xor(args: Vec<Expr>) -> Expr {
    eval_xor(args)
}

/// `Nand(*args) = Not(And(*args))`.
pub fn nand(args: Vec<Expr>) -> Expr {
    not(&and(args))
}

/// `Nor(*args) = Not(Or(*args))`.
pub fn nor(args: Vec<Expr>) -> Expr {
    not(&or(args))
}

/// `Xnor(*args) = Not(Xor(*args))`.
pub fn xnor(args: Vec<Expr>) -> Expr {
    not(&xor(args))
}

/// `Implies(a, b)`.
pub fn implies(a: &Expr, b: &Expr) -> Expr {
    if matches!(a.kind(), Kind::Bool(_)) || matches!(b.kind(), Kind::Bool(_)) {
        return or2(&not(a), b);
    }
    if a == b {
        return Expr::boolean(true);
    }
    if matches!(a.kind(), Kind::Rel(..)) && matches!(b.kind(), Kind::Rel(..)) {
        let (ca, cb) = (canonical(a), canonical(b));
        if ca == cb {
            return Expr::boolean(true);
        }
        if canonical(&negated_rel(a)) == cb {
            return b.clone();
        }
    }
    mk(LogicOp::Implies, vec![a.clone(), b.clone()])
}

/// `Equivalent(*args)`.
pub fn equivalent(args: Vec<Expr>) -> Expr {
    eval_equivalent(args)
}

/// Native nodes belonging to the Boolean category (distinct from is_Boolean).
pub fn is_boolean_category(expression: &Expr) -> bool {
    matches!(expression.kind(), Kind::Sym(_) | Kind::Rel(..) | Kind::Bool(_)
        | Kind::Func(Func::Logic(_) | Func::Predicate(_) | Func::AppliedPredicate
            | Func::Set(crate::sets::SetOp::Contains), _))
}

/// Whether a Piecewise node occurs anywhere in `expression` (shared DAG
/// nodes are visited once).
pub fn contains_piecewise(expression: &Expr) -> bool {
    // Visited nodes stay alive (see `Expr::dag_id`).
    let mut seen = rustc_hash::FxHashMap::default();
    let mut pending = vec![expression.clone()];
    while let Some(node) = pending.pop() {
        if seen.insert(node.dag_id(), node.clone()).is_some() { continue; }
        if matches!(node.kind(), Kind::Piecewise(_)) { return true; }
        pending.extend(node.args());
    }
    false
}

/// A Piecewise condition the ExprCondPair constructor keeps unchanged: a
/// Boolean-category node free of Piecewise (no fold or ITE rewrite applies).
pub fn is_settled_condition(condition: &Expr) -> bool {
    is_boolean_category(condition) && !condition.has_matrix() && !contains_piecewise(condition)
}

/// Convert Boolean-valued Piecewise expressions with an explicit default.
/// Missing coverage or unsupported branches retain the caller's fallback.
pub fn boolean_piecewise_to_ite(expression: &Expr) -> Option<Expr> {
    fn convert(expression: &Expr, depth: usize,
        memo: &mut rustc_hash::FxHashMap<usize, Option<Expr>>) -> Option<Expr> {
        if let Some(value) = memo.get(&expression.dag_id()) { return value.clone(); }
        if depth > 128 || memo.len() >= 100_000 { return None; }
        let result = if let Kind::Piecewise(pieces) = expression.kind() {
            let last = pieces.iter().position(|(_, condition)| matches!(condition.kind(), Kind::Bool(true)))?;
            let mut result = convert(&pieces[last].0, depth+1, memo)?;
            for (value, condition) in pieces[..last].iter().rev() {
                let value = convert(value, depth+1, memo)?;
                let condition = convert(condition, depth+1, memo)?;
                result = ite(&condition, &value, &result);
            }
            Some(result)
        } else if is_boolean_category(expression) { Some(expression.clone()) } else { None };
        memo.insert(expression.dag_id(), result.clone());
        result
    }
    if !matches!(expression.kind(), Kind::Piecewise(_)) {
        return is_boolean_category(expression).then(|| expression.clone());
    }
    convert(expression, 0, &mut rustc_hash::FxHashMap::default())
}

/// `ITE(c, t, e)`: if-then-else.
pub fn ite(c: &Expr, t: &Expr, e: &Expr) -> Expr {
    match c.kind() {
        Kind::Bool(true) => return t.clone(),
        Kind::Bool(false) => return e.clone(),
        _ => {}
    }
    if t == e {
        return t.clone();
    }
    match (t.kind(), e.kind()) {
        (Kind::Bool(true), Kind::Bool(false)) => return c.clone(),
        (Kind::Bool(false), Kind::Bool(true)) => return not(c),
        _ => {}
    }
    mk(LogicOp::ITE, vec![c.clone(), t.clone(), e.clone()])
}

// ----------------------------------------------------------------------
// atoms / evaluation under an assignment
// ----------------------------------------------------------------------

/// The propositional atoms of a boolean expression (symbols, relationals,
/// other non-boolean-function leaves), in SymPy's `ordered()` order.
pub fn atoms(e: &Expr) -> Vec<Expr> {
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let mut out = Vec::new();
    fn rec(e: &Expr, seen: &mut FxHashSet<Expr>, out: &mut Vec<Expr>) {
        match e.kind() {
            Kind::Func(Func::Logic(_), a) => {
                for x in a {
                    rec(x, seen, out);
                }
            }
            Kind::Bool(_) => {}
            _ => {
                if seen.insert(e.clone()) {
                    out.push(e.clone());
                }
            }
        }
    }
    rec(e, &mut seen, &mut out);
    ordered(&mut out);
    out
}

/// Evaluate a boolean expression for an assignment of its atoms. Returns
/// `None` when an atom is not assigned.
pub fn eval_bool(e: &Expr, val: &dyn Fn(&Expr) -> Option<bool>) -> Option<bool> {
    Some(match e.kind() {
        Kind::Bool(b) => *b,
        Kind::Func(Func::Logic(op), a) => match op {
            LogicOp::And => {
                for x in a {
                    if !eval_bool(x, val)? {
                        return Some(false);
                    }
                }
                true
            }
            LogicOp::Or => {
                for x in a {
                    if eval_bool(x, val)? {
                        return Some(true);
                    }
                }
                false
            }
            LogicOp::Not => !eval_bool(&a[0], val)?,
            LogicOp::Xor => {
                let mut r = false;
                for x in a {
                    r ^= eval_bool(x, val)?;
                }
                r
            }
            LogicOp::Implies => !eval_bool(&a[0], val)? || eval_bool(&a[1], val)?,
            LogicOp::Equivalent => {
                let first = eval_bool(&a[0], val)?;
                for x in &a[1..] {
                    if eval_bool(x, val)? != first {
                        return Some(false);
                    }
                }
                true
            }
            LogicOp::ITE => {
                if eval_bool(&a[0], val)? {
                    eval_bool(&a[1], val)?
                } else {
                    eval_bool(&a[2], val)?
                }
            }
        },
        _ => val(e)?,
    })
}

// ----------------------------------------------------------------------
// printing
// ----------------------------------------------------------------------

/// SymPy precedence of a boolean operand.
fn bprec(e: &Expr) -> u32 {
    match e.kind() {
        Kind::Func(Func::Logic(op), _) => match op {
            LogicOp::And => 30,
            LogicOp::Or => 20,
            LogicOp::Xor | LogicOp::Implies | LogicOp::Equivalent => 10,
            LogicOp::Not => 100,
            LogicOp::ITE => 1000,
        },
        Kind::Rel(RelOp::Eq | RelOp::Ne, ..) => 1000,
        Kind::Rel(..) => 35,
        Kind::Add(..) => 40,
        Kind::Mul(c, _) if c.is_negative() => 40,
        Kind::Mul(..) => 50,
        Kind::Pow(..) => 60,
        _ => 1000,
    }
}

fn write_join(args: &[Expr], sep: &str, level: u32, out: &mut String) {
    for (i, a) in args.iter().enumerate() {
        if i > 0 {
            out.push_str(sep);
        }
        if bprec(a) <= level {
            out.push('(');
            crate::printer::write_expr(a, out);
            out.push(')');
        } else {
            crate::printer::write_expr(a, out);
        }
    }
}

fn canonical_rhs_is_neg_oo(e: &Expr) -> bool {
    if !matches!(e.kind(), Kind::Rel(..)) {
        return false;
    }
    matches!(canonical(e).kind(), Kind::Rel(_, _, r) if r.is_constant(Constant::NegInfinity))
}

/// `str()` of a boolean function node.
pub(crate) fn write_logic(op: LogicOp, args: &[Expr], out: &mut String) {
    match op {
        LogicOp::And => {
            let mut v: Vec<Expr> = args.to_vec();
            for j in 0..v.len() {
                if canonical_rhs_is_neg_oo(&v[j]) {
                    let a = v.remove(j);
                    v.insert(0, a);
                }
            }
            write_join(&v, " & ", 38, out)
        }
        LogicOp::Or => write_join(args, " | ", 36, out),
        LogicOp::Xor => write_join(args, " ^ ", 37, out),
        LogicOp::Not if args.len() == 1 => {
            out.push('~');
            write_join(args, "", 100, out);
        }
        _ => {
            out.push_str(op.name());
            out.push('(');
            for (i, a) in args.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                crate::printer::write_expr(a, out);
            }
            out.push(')');
        }
    }
}

#[cfg(test)]
mod tests;
