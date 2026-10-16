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

use rustc_hash::{FxHashMap, FxHashSet};

use crate::expr::{Expr, Kind};
use crate::ops::{add, could_extract_minus_sign, mul, neg, pow};

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
}

impl Default for CseOptions {
    fn default() -> Self {
        CseOptions { symbol_prefix: "x".into(), optimizations: true, ignore: vec![], excluded_names: vec![] }
    }
}

/// SymPy's `cse(exprs)`: `(replacements, reduced_exprs)`.
pub fn cse(exprs: &[Expr]) -> (Vec<(Expr, Expr)>, Vec<Expr>) {
    cse_with(exprs, &CseOptions::default())
}

/// `cse` with options.
pub fn cse_with(exprs: &[Expr], opts: &CseOptions) -> (Vec<(Expr, Expr)>, Vec<Expr>) {
    let mut st = Cse::default();
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
    args_cache: FxHashMap<Expr, Rc<[Expr]>>,
    unevals: Vec<Uneval>,
    opt_subs: FxHashMap<Expr, usize>,
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
            _ => Rc::from(e.args()),
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
                if matches!(e.kind(), Kind::Add(..) | Kind::Mul(..) | Kind::Num(_)) && could_extract_minus_sign(&e) {
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
                        if could_extract_minus_sign(x) {
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
            let mut set = BTreeSet::new();
            for a in self.args(f).iter() {
                let n = t.number(Val::E(a.clone()));
                set.insert(n);
                t.arg_to_funcs[n].insert(fi);
            }
            t.func_to_args.push(set);
        }
        let mut changed: FxHashSet<usize> = FxHashSet::default();
        for i in 0..funcs.len() {
            let counts = t.common_candidates(i, i + 1);
            let mut cands: Vec<(usize, usize)> = counts.into_iter().map(|(k, c)| (c, k)).collect();
            cands.sort_unstable();
            // ordered set of remaining candidates
            let mut remaining: BTreeSet<usize> = cands.iter().map(|(_, k)| *k).collect();
            for &(_, j) in &cands {
                if !remaining.remove(&j) {
                    continue;
                }
                let com: BTreeSet<usize> = t.func_to_args[i].intersection(&t.func_to_args[j]).copied().collect();
                if com.len() <= 1 {
                    continue;
                }
                let diff_i: BTreeSet<usize> = t.func_to_args[i].difference(&com).copied().collect();
                let com_num = if !diff_i.is_empty() {
                    let args: Vec<Val> = com.iter().map(|n| t.values[*n].clone()).collect();
                    let u = self.new_uneval(op, args);
                    let num = t.number(Val::U(u));
                    let mut ns = diff_i;
                    ns.insert(num);
                    t.update(i, ns);
                    changed.insert(i);
                    num
                } else {
                    t.number(Val::E(funcs[i].clone()))
                };
                let mut diff_j: BTreeSet<usize> = t.func_to_args[j].difference(&com).copied().collect();
                diff_j.insert(com_num);
                t.update(j, diff_j);
                changed.insert(j);
                for k in t.subset_candidates(&com, &remaining) {
                    let mut diff_k: BTreeSet<usize> = t.func_to_args[k].difference(&com).copied().collect();
                    diff_k.insert(com_num);
                    t.update(k, diff_k);
                    changed.insert(k);
                }
            }
            if changed.contains(&i) {
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
            // explicit stack: (value, children pushed?)
            let mut stack: Vec<(Val, bool)> = vec![(Val::E(root.clone()), false)];
            let mut results: Vec<Expr> = Vec::new();
            while let Some((v, expanded)) = stack.pop() {
                if !expanded {
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
                    stack.push((v, true));
                    for c in ch.into_iter().rev() {
                        stack.push((c, false));
                    }
                    continue;
                }
                let ch_len = self.ordered_children_len(&v);
                let new_args: Vec<Expr> = results.split_off(results.len() - ch_len);
                let new_expr = match &v {
                    Val::U(u) => build(self.unevals[*u].op, &new_args),
                    Val::E(e) => {
                        if let Some(&u) = self.opt_subs.get(e) {
                            build(self.unevals[u].op, &new_args)
                        } else {
                            let old = self.ordered_children(&v);
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
            if !self.opt_subs.contains_key(e) && matches!(e.kind(), Kind::Add(..) | Kind::Mul(..)) {
                ch.sort_by(|a, b| match (a, b) {
                    (Val::E(x), Val::E(y)) => x.cmp_struct(y),
                    _ => std::cmp::Ordering::Equal,
                });
            }
        }
        ch
    }

    fn ordered_children_len(&mut self, v: &Val) -> usize {
        match v {
            Val::E(e) => match self.opt_subs.get(e) {
                Some(&u) => self.unevals[u].args.len(),
                None => self.args(e).len(),
            },
            Val::U(u) => self.unevals[*u].args.len(),
        }
    }
}

/// `FuncArgTracker`: functions ↔ argument value numbers.
#[derive(Default)]
struct ArgTracker {
    numbers: FxHashMap<Val, usize>,
    values: Vec<Val>,
    arg_to_funcs: Vec<BTreeSet<usize>>,
    func_to_args: Vec<BTreeSet<usize>>,
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
        let args: Vec<usize> = self.func_to_args[fi].iter().copied().collect();
        for a in args {
            self.arg_to_funcs[a].remove(&fi);
        }
    }

    /// Functions `>= min` sharing at least two arguments with function
    /// `fi`, with the number of shared arguments.
    fn common_candidates(&self, fi: usize, min: usize) -> FxHashMap<usize, usize> {
        let mut count: FxHashMap<usize, usize> = FxHashMap::default();
        let argset = &self.func_to_args[fi];
        if argset.is_empty() {
            return count;
        }
        let largest = argset.iter().max_by_key(|a| self.arg_to_funcs[**a].len()).copied().unwrap();
        for a in argset {
            if *a == largest {
                continue;
            }
            for &f in self.arg_to_funcs[*a].range(min..) {
                *count.entry(f).or_insert(0) += 1;
            }
        }
        let big = &self.arg_to_funcs[largest];
        if big.len() < count.len() {
            for &f in big.range(min..) {
                if let Some(c) = count.get_mut(&f) {
                    *c += 1;
                }
            }
        } else {
            for (f, c) in count.iter_mut() {
                if big.contains(f) {
                    *c += 1;
                }
            }
        }
        count.retain(|_, c| *c >= 2);
        count
    }

    fn subset_candidates(&self, argset: &BTreeSet<usize>, restrict: &BTreeSet<usize>) -> Vec<usize> {
        let mut it = argset.iter();
        let first = *it.next().unwrap();
        let mut idx: BTreeSet<usize> = self.arg_to_funcs[first].intersection(restrict).copied().collect();
        for a in it {
            idx = idx.intersection(&self.arg_to_funcs[*a]).copied().collect();
            if idx.is_empty() {
                break;
            }
        }
        idx.into_iter().collect()
    }

    fn update(&mut self, fi: usize, new: BTreeSet<usize>) {
        let old = std::mem::take(&mut self.func_to_args[fi]);
        for d in old.difference(&new) {
            self.arg_to_funcs[*d].remove(&fi);
        }
        for a in new.difference(&old) {
            self.arg_to_funcs[*a].insert(fi);
        }
        self.func_to_args[fi] = new;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    use crate::subs::xreplace;

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
