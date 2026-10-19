//! Common subexpression elimination with SymPy's `cse` semantics.
//!
//! The algorithm follows SymPy's two stages:
//!
//! 1. *opt_cse*: find optimisation opportunities — negative coefficients
//!    (`-x*y` shares `x*y`), negative powers (`x**-2` shares `x**2`) and
//!    common argument subsets of `Add`s and `Mul`s (`x + y + z` and
//!    `x + y + w` share `x + y`).
//! 2. *tree_cse*: count occurrences of every sub-expression of the
//!    (optimised) DAG and pull out the ones seen at least twice, numbering
//!    them `x0, x1, ...` in post-order.
//!
//! Every distinct node is visited once: nodes are keyed by their cached
//! structural hash, so the cost is linear in the size of the expression DAG
//! (plus the argument-subset matching, which uses an inverted index).

use std::collections::BTreeSet;
use std::rc::Rc;

use num_bigint::BigUint;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::expr::{Expr, Kind};
use crate::ops::{add, mul, neg, pow};

/// Options of [`cse_with`].
#[derive(Clone, Debug)]
pub struct CseOptions {
    /// Prefix of the generated symbols (`x` gives `x0, x1, ...`).
    pub symbol_prefix: String,
    /// Run the `opt_cse` pre-pass (minus signs, negative powers, common
    /// `Add`/`Mul` argument subsets). Without it only identical
    /// sub-trees are shared.
    pub optimizations: bool,
    /// Symbols (names) whose sub-expressions must not be extracted.
    pub ignore: Vec<String>,
    /// Extra names the generated symbols must not take (e.g. function
    /// arguments that do not appear in the expressions).
    pub excluded_names: Vec<String>,
    /// Sort rebuilt sums/products with the reference node-count/key order.
    pub canonical_order: bool,
}

impl Default for CseOptions {
    fn default() -> Self {
        CseOptions {
            symbol_prefix: "x".into(),
            optimizations: true,
            ignore: vec![],
            excluded_names: vec![],
            canonical_order: true,
        }
    }
}

/// SymPy's `cse(exprs)`: `(replacements, reduced_exprs)`.
pub fn cse(exprs: &[Expr]) -> (Vec<(Expr, Expr)>, Vec<Expr>) {
    cse_with(exprs, &CseOptions::default())
}

/// `cse` with options.
pub fn cse_with(exprs: &[Expr], opts: &CseOptions) -> (Vec<(Expr, Expr)>, Vec<Expr>) {
    let mut st = Cse {
        canonical_order: opts.canonical_order,
        ..Default::default()
    };
    if opts.optimizations {
        st.opt_cse(exprs);
    }
    st.tree_cse(exprs, opts)
}

// ----------------------------------------------------------------------
// internal representation
// ----------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Op {
    Add,
    Mul,
    Pow,
}

/// A value: an expression or an unevaluated `Add`/`Mul`/`Pow` (by index).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Val {
    E(Expr),
    U(usize),
}

#[derive(Debug)]
struct Uneval {
    op: Op,
    args: Vec<Val>,
}

#[derive(Default)]
struct Cse {
    keys: super::order::KeyCache,
    args_cache: FxHashMap<Expr, Rc<[Expr]>>,
    unevals: Vec<Uneval>,
    opt_subs: FxHashMap<Expr, usize>,
    nodes: FxHashMap<Expr, BigUint>,
    canonical_order: bool,
}

/// SymPy's `could_extract_minus_sign` (ties of a sum are broken with
/// `default_sort_key`, like SymPy).
fn could_extract_minus_sign(e: &Expr, keys: &mut super::order::KeyCache) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_negative(),
        Kind::Const(crate::expr::Constant::NegInfinity) => true,
        Kind::Mul(c, _) => c.is_negative(),
        Kind::Add(c, terms) => {
            let mut negs = terms.iter().filter(|(_, k)| k.is_negative()).count();
            let mut n = terms.len();
            if !c.is_zero() {
                n += 1;
                if c.is_negative() {
                    negs += 1;
                }
            }
            let pos = n - negs;
            if pos != negs {
                return negs > pos;
            }
            keys.key(e) < keys.key(&neg(e))
        }
        _ => false,
    }
}

/// Leaves of the CSE tree. Derivatives are opaque (their inner expression
/// must not be replaced by a symbol).
fn is_leaf(e: &Expr) -> bool {
    e.is_atom() || matches!(e.kind(), Kind::Derivative(..))
}

fn rebuild_kind(e: &Expr, args: &[Expr]) -> Expr {
    match e.kind() {
        Kind::Add(..) => add(args),
        Kind::Mul(..) => mul(args),
        Kind::Pow(..) => pow(&args[0], &args[1]),
        Kind::Func(f, _) => crate::functions::func(f.clone(), args.to_vec()),
        Kind::Rel(op, ..) => crate::relational::rel(*op, &args[0], &args[1]),
        Kind::Piecewise(_) => crate::relational::piecewise(
            args.chunks(2).map(|c| (c[0].clone(), c[1].clone())).collect(),
        ),
        _ => e.clone(),
    }
}

fn build(op: Op, args: &[Expr]) -> Expr {
    match op {
        Op::Add => add(args),
        Op::Mul => mul(args),
        Op::Pow => pow(&args[0], &args[1]),
    }
}

impl Cse {
    fn args(&mut self, e: &Expr) -> Rc<[Expr]> {
        if let Some(a) = self.args_cache.get(e) {
            return a.clone();
        }
        let a: Rc<[Expr]> = match e.kind() {
            Kind::Derivative(..) => Rc::from(Vec::new()),
            _ => Rc::from(crate::traverse::ordered_args(e)),
        };
        self.args_cache.insert(e.clone(), a.clone());
        a
    }

    fn new_uneval(&mut self, op: Op, args: Vec<Val>) -> usize {
        self.unevals.push(Uneval { op, args });
        self.unevals.len() - 1
    }

    // ------------------------------------------------------------------
    // opt_cse
    // ------------------------------------------------------------------

    fn opt_cse(&mut self, exprs: &[Expr]) {
        let mut seen: FxHashSet<Expr> = FxHashSet::default();
        let mut adds: Vec<Expr> = Vec::new();
        let mut muls: Vec<Expr> = Vec::new();
        let mut in_list: FxHashSet<Expr> = FxHashSet::default();
        // iterative post-order traversal (deep expressions must not
        // overflow the stack)
        for root in exprs {
            let mut stack: Vec<(Expr, bool)> = vec![(root.clone(), false)];
            while let Some((e, expanded)) = stack.pop() {
                if is_leaf(&e) {
                    continue;
                }
                if !expanded {
                    if !seen.insert(e.clone()) {
                        continue;
                    }
                    stack.push((e.clone(), true));
                    let a = self.args(&e);
                    for c in a.iter().rev() {
                        if !is_leaf(c) && !seen.contains(c) {
                            stack.push((c.clone(), false));
                        }
                    }
                    continue;
                }
                let mut expr = e.clone();
                if matches!(e.kind(), Kind::Add(..) | Kind::Mul(..) | Kind::Num(_)) && could_extract_minus_sign(&e, &mut self.keys) {
                    let ne = neg(&e);
                    if !ne.is_atom() {
                        let u = self.new_uneval(Op::Mul, vec![Val::E(Expr::minus_one()), Val::E(ne.clone())]);
                        self.opt_subs.insert(e.clone(), u);
                        seen.insert(ne.clone());
                        expr = ne;
                    }
                }
                match expr.kind() {
                    Kind::Mul(..) => {
                        if in_list.insert(expr.clone()) {
                            muls.push(expr);
                        }
                    }
                    Kind::Add(..) => {
                        if in_list.insert(expr.clone()) {
                            adds.push(expr);
                        }
                    }
                    Kind::Pow(b, x) => {
                        if could_extract_minus_sign(x, &mut self.keys) {
                            let p = pow(b, &neg(x));
                            let u = self.new_uneval(Op::Pow, vec![Val::E(p), Val::E(Expr::minus_one())]);
                            self.opt_subs.insert(expr.clone(), u);
                        }
                    }
                    _ => {}
                }
            }
        }
        self.match_common_args(Op::Add, adds);
        self.match_common_args(Op::Mul, muls);
    }

    fn match_common_args(&mut self, op: Op, funcs: Vec<Expr>) {
        if funcs.len() < 2 {
            return;
        }
        let mut funcs: Vec<(usize, Expr)> = funcs
            .into_iter()
            .map(|f| {
                let n = self.args(&f).len();
                (n, f)
            })
            .collect();
        funcs.sort_by_key(|(n, _)| *n);
        let funcs: Vec<Expr> = funcs.into_iter().map(|(_, f)| f).collect();

        let mut t = ArgTracker::default();
        for (fi, f) in funcs.iter().enumerate() {
            let mut set: Vec<usize> = Vec::new();
            for a in self.args(f).iter() {
                let n = t.number(Val::E(a.clone()));
                set.push(n);
                t.arg_to_funcs[n].insert(fi);
            }
            set.sort_unstable();
            set.dedup();
            t.func_to_args.push(set);
        }
        let mut changed = vec![false; funcs.len()];
        let mut remaining = vec![false; funcs.len()];
        for i in 0..funcs.len() {
            let counts = t.common_candidates(i, i + 1);
            let mut cands: Vec<(usize, usize)> = counts.into_iter().map(|(k, c)| (c, k)).collect();
            cands.sort_unstable();
            for &(_, j) in &cands {
                remaining[j] = true;
            }
            for &(_, j) in &cands {
                if !remaining[j] {
                    continue;
                }
                remaining[j] = false;
                let com = sorted_inter(&t.func_to_args[i], &t.func_to_args[j]);
                if com.len() <= 1 {
                    continue;
                }
                let diff_i = sorted_diff(&t.func_to_args[i], &com);
                let com_num = if !diff_i.is_empty() {
                    let args: Vec<Val> = com.iter().map(|n| t.values[*n].clone()).collect();
                    let u = self.new_uneval(op, args);
                    let num = t.number(Val::U(u));
                    t.update(i, sorted_insert(diff_i, num));
                    changed[i] = true;
                    num
                } else {
                    t.number(Val::E(funcs[i].clone()))
                };
                let diff_j = sorted_diff(&t.func_to_args[j], &com);
                t.update(j, sorted_insert(diff_j, com_num));
                changed[j] = true;
                for k in t.subset_candidates(&com, &remaining) {
                    let diff_k = sorted_diff(&t.func_to_args[k], &com);
                    t.update(k, sorted_insert(diff_k, com_num));
                    changed[k] = true;
                }
            }
            for &(_, j) in &cands {
                remaining[j] = false;
            }
            if changed[i] {
                let args: Vec<Val> = t.func_to_args[i].iter().map(|n| t.values[*n].clone()).collect();
                let u = self.new_uneval(op, args);
                self.opt_subs.insert(funcs[i].clone(), u);
            }
            t.stop_tracking(i);
        }
    }

    // ------------------------------------------------------------------
    // tree_cse
    // ------------------------------------------------------------------

    fn children(&mut self, v: &Val) -> Vec<Val> {
        match v {
            Val::E(e) => {
                if let Some(&u) = self.opt_subs.get(e) {
                    self.unevals[u].args.clone()
                } else {
                    self.args(e).iter().map(|a| Val::E(a.clone())).collect()
                }
            }
            Val::U(u) => self.unevals[*u].args.clone(),
        }
    }

    fn tree_cse(&mut self, exprs: &[Expr], opts: &CseOptions) -> (Vec<(Expr, Expr)>, Vec<Expr>) {
        let mut seen: FxHashSet<Val> = FxHashSet::default();
        let mut to_eliminate: FxHashSet<Val> = FxHashSet::default();
        let mut excluded: FxHashSet<String> = opts.excluded_names.iter().cloned().collect();
        let ignore: Vec<Expr> = opts.ignore.iter().map(|s| Expr::symbol(s)).collect();

        // find repeated (iterative pre-order)
        for root in exprs {
            let mut stack = vec![Val::E(root.clone())];
            while let Some(v) = stack.pop() {
                if let Val::E(e) = &v {
                    if is_leaf(e) {
                        if let Some(n) = e.as_symbol() {
                            excluded.insert(n.to_string());
                        }
                        if !matches!(e.kind(), Kind::Derivative(..)) {
                            continue;
                        }
                    }
                }
                if seen.contains(&v) {
                    let ignored = !ignore.is_empty()
                        && match &v {
                            Val::E(e) => ignore.iter().any(|s| crate::traverse::has(e, s)),
                            Val::U(_) => {
                                let mut it = self.children(&v);
                                let mut found = false;
                                while let Some(c) = it.pop() {
                                    match c {
                                        Val::E(e) => {
                                            if ignore.iter().any(|s| crate::traverse::has(&e, s)) {
                                                found = true;
                                                break;
                                            }
                                        }
                                        u => it.extend(self.children(&u)),
                                    }
                                }
                                found
                            }
                        };
                    if !ignored {
                        to_eliminate.insert(v);
                        continue;
                    }
                }
                seen.insert(v.clone());
                let ch = self.children(&v);
                for c in ch.into_iter().rev() {
                    stack.push(c);
                }
            }
        }

        // rebuild (post-order)
        let mut counter = 0usize;
        let mut next_symbol = || loop {
            let name = format!("{}{}", opts.symbol_prefix, counter);
            counter += 1;
            if !excluded.contains(&name) {
                return Expr::symbol(&name);
            }
        };
        let mut subs: FxHashMap<Val, Expr> = FxHashMap::default();
        let mut replacements: Vec<(Expr, Expr)> = Vec::new();
        let mut reduced = Vec::with_capacity(exprs.len());
        for root in exprs {
            // explicit stack: Some(children) once the children are pushed
            let mut stack: Vec<(Val, Option<Vec<Val>>)> = vec![(Val::E(root.clone()), None)];
            let mut results: Vec<Expr> = Vec::new();
            while let Some((v, done)) = stack.pop() {
                let Some(old) = done else {
                    if let Val::E(e) = &v {
                        if e.is_atom() {
                            results.push(e.clone());
                            continue;
                        }
                    }
                    if let Some(s) = subs.get(&v) {
                        results.push(s.clone());
                        continue;
                    }
                    let ch = self.ordered_children(&v);
                    let pushed: Vec<Val> = ch.iter().rev().cloned().collect();
                    stack.push((v, Some(ch)));
                    for c in pushed {
                        stack.push((c, None));
                    }
                    continue;
                };
                let new_args: Vec<Expr> = results.split_off(results.len() - old.len());
                let new_expr = match &v {
                    Val::U(u) => build(self.unevals[*u].op, &new_args),
                    Val::E(e) => {
                        if let Some(&u) = self.opt_subs.get(e) {
                            build(self.unevals[u].op, &new_args)
                        } else {
                            let same = old.iter().zip(&new_args).all(|(o, n)| matches!(o, Val::E(x) if x.ptr_eq(n)));
                            if same {
                                e.clone()
                            } else {
                                rebuild_kind(e, &new_args)
                            }
                        }
                    }
                };
                if to_eliminate.contains(&v) {
                    let sym = next_symbol();
                    subs.insert(v, sym.clone());
                    replacements.push((sym.clone(), new_expr));
                    results.push(sym);
                } else {
                    results.push(new_expr);
                }
            }
            reduced.push(results.pop().unwrap());
        }
        (replacements, reduced)
    }

    /// Children in SymPy's canonical processing order (`Add`/`Mul`
    /// arguments sorted, everything else as is).
    fn ordered_children(&mut self, v: &Val) -> Vec<Val> {
        let mut ch = self.children(v);
        if let Val::E(e) = v {
            if self.canonical_order
                && !self.opt_subs.contains_key(e)
                && matches!(e.kind(), Kind::Add(..) | Kind::Mul(..))
            {
                // ordered() first compares expanded node counts, then keys.
                let mut keyed: Vec<(BigUint, Rc<super::order::SortKey>, Val)> = ch
                    .into_iter()
                    .map(|v| match &v {
                        Val::E(x) => (self.node_size(x), self.keys.key(x), v),
                        Val::U(_) => unreachable!(),
                    })
                    .collect();
                keyed.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
                ch = keyed.into_iter().map(|(_, _, v)| v).collect();
            }
        }
        ch
    }

    // Count the mathematical tree using its DAG, with no depth/size cutoff.
    // Big counts distinguish shared trees larger than a machine counter.
    fn node_size(&mut self, expression: &Expr) -> BigUint {
        let mut pending = vec![(expression.clone(), false)];
        while let Some((node, expanded)) = pending.pop() {
            if self.nodes.contains_key(&node) {
                continue;
            }
            let children = match node.kind() {
                Kind::Derivative(body, _) => vec![body.clone()],
                _ => self.args(&node).to_vec(),
            };
            if !expanded {
                pending.push((node, true));
                pending.extend(children.into_iter().rev().map(|child| (child, false)));
                continue;
            }
            let size = if let Kind::Derivative(body, variables) = node.kind() {
                variables
                    .iter()
                    .fold(self.nodes[body].clone(), |size, (_, order)| {
                        size + BigUint::from(*order)
                    })
            } else {
                children
                    .iter()
                    .fold(BigUint::from(1u32), |size, child| size + &self.nodes[child])
            };
            self.nodes.insert(node, size);
        }
        self.nodes[expression].clone()
    }

}

fn sorted_inter(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

fn sorted_diff(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len());
    let mut j = 0;
    for &x in a {
        while j < b.len() && b[j] < x {
            j += 1;
        }
        if j < b.len() && b[j] == x {
            continue;
        }
        out.push(x);
    }
    out
}

fn sorted_insert(mut v: Vec<usize>, x: usize) -> Vec<usize> {
    if let Err(p) = v.binary_search(&x) {
        v.insert(p, x);
    }
    v
}

fn is_subset(sub: &[usize], sup: &[usize]) -> bool {
    let mut j = 0;
    for &x in sub {
        while j < sup.len() && sup[j] < x {
            j += 1;
        }
        if j >= sup.len() || sup[j] != x {
            return false;
        }
        j += 1;
    }
    true
}

/// `FuncArgTracker`: functions ↔ argument value numbers.
#[derive(Default)]
struct ArgTracker {
    numbers: FxHashMap<Val, usize>,
    values: Vec<Val>,
    arg_to_funcs: Vec<BTreeSet<usize>>,
    /// sorted argument numbers of every function
    func_to_args: Vec<Vec<usize>>,
    /// scratch counters of `common_candidates`
    counts: Vec<u32>,
}

impl ArgTracker {
    fn number(&mut self, v: Val) -> usize {
        if let Some(&n) = self.numbers.get(&v) {
            return n;
        }
        let n = self.values.len();
        self.numbers.insert(v.clone(), n);
        self.values.push(v);
        self.arg_to_funcs.push(BTreeSet::new());
        n
    }

    fn stop_tracking(&mut self, fi: usize) {
        for k in 0..self.func_to_args[fi].len() {
            let a = self.func_to_args[fi][k];
            self.arg_to_funcs[a].remove(&fi);
        }
    }

    /// Functions `>= min` sharing at least two arguments with function
    /// `fi`, with the number of shared arguments.
    fn common_candidates(&mut self, fi: usize, min: usize) -> Vec<(usize, usize)> {
        let argset = &self.func_to_args[fi];
        if argset.is_empty() {
            return Vec::new();
        }
        if self.counts.len() < self.func_to_args.len() {
            self.counts.resize(self.func_to_args.len(), 0);
        }
        // the largest function set is only used to bump existing counts
        let largest = argset.iter().max_by_key(|a| self.arg_to_funcs[**a].len()).copied().unwrap();
        let mut touched: Vec<usize> = Vec::new();
        for a in argset {
            if *a == largest {
                continue;
            }
            for &f in self.arg_to_funcs[*a].range(min..) {
                if self.counts[f] == 0 {
                    touched.push(f);
                }
                self.counts[f] += 1;
            }
        }
        let big = &self.arg_to_funcs[largest];
        for f in &touched {
            if big.contains(f) {
                self.counts[*f] += 1;
            }
        }
        let mut out = Vec::new();
        for f in touched {
            let c = self.counts[f] as usize;
            self.counts[f] = 0;
            if c >= 2 {
                out.push((f, c));
            }
        }
        out
    }

    /// Remaining candidate functions whose arguments contain `argset`.
    fn subset_candidates(&self, argset: &[usize], remaining: &[bool]) -> Vec<usize> {
        let smallest = argset.iter().min_by_key(|a| self.arg_to_funcs[**a].len()).copied().unwrap();
        self.arg_to_funcs[smallest]
            .iter()
            .copied()
            .filter(|f| remaining[*f] && is_subset(argset, &self.func_to_args[*f]))
            .collect()
    }

    fn update(&mut self, fi: usize, new: Vec<usize>) {
        let old = std::mem::take(&mut self.func_to_args[fi]);
        for d in sorted_diff(&old, &new) {
            self.arg_to_funcs[d].remove(&fi);
        }
        for a in sorted_diff(&new, &old) {
            self.arg_to_funcs[a].insert(fi);
        }
        self.func_to_args[fi] = new;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    use crate::subs::xreplace;

    #[test]
    fn ordering_counts_shared_trees_beyond_machine_size() {
        let head = crate::expr::Func::Undefined(crate::expr::UndefinedFunction::new(
            "f",
            crate::Assumptions::NONE,
            None,
        ));
        let mut expression = Expr::symbol("x");
        for _ in 0..80 {
            expression = crate::functions::func(head.clone(), vec![expression.clone(), expression]);
        }
        let mut state = Cse::default();
        assert_eq!(
            state.node_size(&expression),
            (BigUint::from(1u32) << 81usize) - BigUint::from(1u32)
        );
    }

    fn back(reps: &[(Expr, Expr)], e: &Expr) -> Expr {
        let mut r = e.clone();
        for (s, v) in reps.iter().rev() {
            let mut m = FxHashMap::default();
            m.insert(s.clone(), v.clone());
            r = xreplace(&r, &m);
        }
        r
    }

    #[test]
    fn simple() {
        let e = parse_str("(x + y)**2 + sin(x + y)").unwrap();
        let (r, red) = cse(&[e.clone()]);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].0.to_string(), "x0");
        assert_eq!(r[0].1.to_string(), "x + y");
        assert_eq!(red[0].to_string(), "x0**2 + sin(x0)");
        assert_eq!(back(&r, &red[0]), e);
    }

    #[test]
    fn common_add_args() {
        let a = parse_str("x + y + z").unwrap();
        let b = parse_str("sin(x + y + w)").unwrap();
        let (r, red) = cse(&[a.clone(), b.clone()]);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].1.to_string(), "x + y");
        assert_eq!(back(&r, &red[0]), a);
        assert_eq!(back(&r, &red[1]), b);
    }

    #[test]
    fn negatives() {
        let a = parse_str("x**-2 + x**2*sin(y)").unwrap();
        let (r, red) = cse(&[a.clone()]);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].1.to_string(), "x**2");
        assert_eq!(back(&r, &red[0]), a);
        let a = parse_str("exp(-x*y) + cos(x*y)").unwrap();
        let (r, red) = cse(&[a.clone()]);
        assert_eq!(r.len(), 1);
        assert_eq!(back(&r, &red[0]), a);
    }

    #[test]
    fn excluded_names() {
        let a = parse_str("(x0 + x1)**2 + sin(x0 + x1)").unwrap();
        let (r, _) = cse(&[a]);
        assert_eq!(r[0].0.to_string(), "x2");
    }
}
