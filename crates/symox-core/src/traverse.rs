//! Generic tree traversal helpers.

use rustc_hash::FxHashSet;

use crate::expr::{Expr, Kind};
use crate::functions::func;
use crate::ops::{add, mul, mul2, pow};
use crate::relational::{piecewise, rel};

/// Rebuild `e` with every direct child mapped through `f`, using the
/// canonical constructors. Returns `e` itself when nothing changed.
pub fn map_children(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => e.clone(),
        Kind::Add(c, terms) => {
            let mut changed = false;
            let mut v: Vec<Expr> = Vec::with_capacity(terms.len() + 1);
            let mut keep: Vec<(Expr, crate::number::Number)> = Vec::new();
            for (t, k) in terms {
                let nt = f(t);
                if nt.ptr_eq(t) || nt == *t {
                    keep.push((t.clone(), k.clone()));
                } else {
                    changed = true;
                    v.push(mul2(&Expr::num(k.clone()), &nt));
                }
            }
            if !changed {
                return e.clone();
            }
            v.push(crate::ops::add_from_terms(c.clone(), keep));
            add(&v)
        }
        Kind::Mul(c, facs) => {
            let mut changed = false;
            let mut v: Vec<Expr> = Vec::with_capacity(facs.len() + 1);
            v.push(Expr::num(c.clone()));
            for (b, x) in facs {
                let nb = f(b);
                let nx = if x.is_number() { x.clone() } else { f(x) };
                if !(nb == *b && nx == *x) {
                    changed = true;
                }
                v.push(if x.is_one() && nx.is_one() { nb } else { pow(&nb, &nx) });
            }
            if !changed {
                return e.clone();
            }
            mul(&v)
        }
        Kind::Pow(b, x) => {
            let (nb, nx) = (f(b), f(x));
            if nb == *b && nx == *x {
                return e.clone();
            }
            pow(&nb, &nx)
        }
        Kind::Func(fun, args) => {
            let na: Vec<Expr> = args.iter().map(|a| f(a)).collect();
            if na.iter().zip(args).all(|(a, b)| a == b) {
                return e.clone();
            }
            func(fun.clone(), na)
        }
        Kind::Derivative(inner, vars) => {
            let ni = f(inner);
            if ni == *inner {
                return e.clone();
            }
            let mut r = ni;
            for (v, n) in vars {
                r = crate::diff::diff_n(&r, v, *n);
            }
            r
        }
        Kind::Rel(op, a, b) => {
            let (na, nb) = (f(a), f(b));
            if na == *a && nb == *b {
                return e.clone();
            }
            rel(*op, &na, &nb)
        }
        Kind::Piecewise(p) => {
            let np: Vec<(Expr, Expr)> = p.iter().map(|(x, c)| (f(x), f(c))).collect();
            if np == *p {
                return e.clone();
            }
            piecewise(np)
        }
    }
}

/// Pre-order iterator over all sub-expressions (with repetitions).
pub fn preorder(e: &Expr, f: &mut dyn FnMut(&Expr)) {
    f(e);
    for a in e.args() {
        preorder(&a, f);
    }
}

/// Set of free symbols.
pub fn free_symbols(e: &Expr) -> Vec<Expr> {
    e.cached_free_symbols(|| collect_free_symbols(e))
}

fn collect_free_symbols(e: &Expr) -> Vec<Expr> {
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let mut out = Vec::new();
    fn rec(e: &Expr, seen: &mut FxHashSet<Expr>, out: &mut Vec<Expr>, visited: &mut FxHashSet<Expr>) {
        if e.mask() == 0 {
            return;
        }
        if e.is_symbol() {
            if seen.insert(e.clone()) {
                out.push(e.clone());
            }
            return;
        }
        if !visited.insert(e.clone()) {
            return;
        }
        match e.kind() {
            Kind::Func(crate::expr::Func::Set(op), args)
                if matches!(op, crate::sets::SetOp::Lambda | crate::sets::SetOp::ConditionSet)
                    && args.len() >= 2 => {
                // The binder applies to the body/condition, but not the base set.
                for symbol in free_symbols(&args[1]).into_iter().filter(|s| *s != args[0]) {
                    if seen.insert(symbol.clone()) {
                        out.push(symbol);
                    }
                }
                for base in &args[2..] {
                    rec(base, seen, out, visited);
                }
            }
            Kind::Add(_, terms) => {
                for (t, _) in terms {
                    rec(t, seen, out, visited)
                }
            }
            Kind::Mul(_, facs) => {
                for (b, x) in facs {
                    rec(b, seen, out, visited);
                    rec(x, seen, out, visited);
                }
            }
            _ => {
                for a in e.args() {
                    rec(&a, seen, out, visited)
                }
            }
        }
    }
    let mut visited = FxHashSet::default();
    rec(e, &mut seen, &mut out, &mut visited);
    out.sort_by(|a, b| a.cmp_struct(b));
    out
}

/// Does `e` contain the sub-expression `s` (structurally)?
pub fn has(e: &Expr, s: &Expr) -> bool {
    if s.mask() & !e.mask() != 0 {
        return false;
    }
    if e == s {
        return true;
    }
    match e.kind() {
        Kind::Add(c, terms) => {
            if let Kind::Add(sc, subterms) = s.kind() {
                // `x + y` is contained in `x + y + z`
                // Structural containment does not allow the sign changes or
                // constant offsets supported by algebraic substitution.
                if (sc.is_zero() || c == sc)
                    && subterms.iter().all(|term| terms.contains(term))
                {
                    return true;
                }
            }
            terms.iter().any(|(t, _)| has(t, s))
        }
        Kind::Mul(_, facs) => facs.iter().any(|(b, x)| has(b, s) || has(x, s)),
        _ => e.args().iter().any(|a| has(a, s)),
    }
}

/// Number of nodes (SymPy's `count_ops`-like size measure).
pub fn node_count(e: &Expr) -> usize {
    1 + e.args().iter().map(node_count).sum::<usize>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbol_cache_preserves_assumptions_and_replacements() {
        let x = Expr::symbol("x");
        let y = Expr::symbol_with("x", crate::symbol::Assumptions::from_kwargs([("positive", true)]).unwrap());
        let e = add(&[x.clone(), y.clone()]);
        assert_eq!(free_symbols(&e).len(), 2);
        let mut returned = free_symbols(&e);
        returned.clear();
        assert_eq!(free_symbols(&e).len(), 2);
        let map = [(x, Expr::one())].into_iter().collect();
        let replaced = crate::subs::xreplace(&e, &map);
        assert_eq!(free_symbols(&replaced), vec![y]);
    }

    #[test]
    fn structural_sum_containment_differs_from_substitution() {
        let parse = |s| crate::parser::parse_str(s).unwrap();
        assert!(has(&parse("x+y+z"), &parse("x+y")));
        assert!(!has(&parse("x+y"), &parse("x+y+1")));
        assert!(!has(&parse("x-y"), &parse("y-x")));
        assert!(!has(&parse("2*x+2*y"), &parse("x+y")));
    }
}
