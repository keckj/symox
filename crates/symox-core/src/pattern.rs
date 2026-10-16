//! Pattern matching with wildcards: `match`, `replace`, `find`.
//!
//! Sums and products are matched modulo commutativity with a small
//! backtracking search; a lone wildcard term (or factor) absorbs whatever
//! is left (possibly 0 or 1), which gives SymPy's behaviour on the usual
//! cases (`(2*x + 3).match(a*x + b)` gives `{a: 2, b: 3}`).

use rustc_hash::FxHashMap;

use crate::expr::{Expr, Kind};
use crate::ops::{add, div, mul, sub};
use crate::subs::xreplace;
use crate::traverse::{has, map_children};

pub type Bindings = Vec<(Expr, Expr)>;

fn is_wild(e: &Expr) -> bool {
    e.as_sym().is_some_and(|s| s.is_wild())
}

/// Does the pattern contain wildcards?
pub fn has_wild(p: &Expr) -> bool {
    if p.mask() == 0 {
        return false;
    }
    match p.kind() {
        Kind::Sym(s) => s.is_wild(),
        Kind::Add(_, terms) => terms.iter().any(|(t, _)| has_wild(t)),
        Kind::Mul(_, facs) => facs.iter().any(|(b, e)| has_wild(b) || has_wild(e)),
        _ => p.args().iter().any(has_wild),
    }
}

fn lookup<'a>(b: &'a Bindings, w: &Expr) -> Option<&'a Expr> {
    b.iter().find(|(k, _)| k == w).map(|(_, v)| v)
}

fn bind(w: &Expr, e: &Expr, b: &mut Bindings) -> bool {
    if let Some(v) = lookup(b, w) {
        return v == e;
    }
    let s = w.as_sym().unwrap();
    if s.exclude().iter().any(|x| has(e, x)) {
        return false;
    }
    b.push((w.clone(), e.clone()));
    true
}

/// Match `e` against pattern `p`; returns the wildcard bindings.
pub fn pattern_match(e: &Expr, p: &Expr) -> Option<FxHashMap<Expr, Expr>> {
    let mut b = Bindings::new();
    if m(e, p, &mut b) {
        Some(b.into_iter().collect())
    } else {
        None
    }
}

fn m(e: &Expr, p: &Expr, b: &mut Bindings) -> bool {
    if !has_wild(p) {
        return e == p;
    }
    if is_wild(p) {
        return bind(p, e, b);
    }
    match p.kind() {
        Kind::Add(..) => match_add(e, p, b),
        Kind::Mul(..) => match_mul(e, p, b),
        Kind::Pow(pb, pe) => {
            if let Kind::Pow(eb, ee) = e.kind() {
                let save = b.len();
                if m(eb, pb, b) && m(ee, pe, b) {
                    return true;
                }
                b.truncate(save);
            }
            // x matches x**a with a = 1
            let save = b.len();
            if m(e, pb, b) && m(&Expr::one(), pe, b) {
                return true;
            }
            b.truncate(save);
            false
        }
        Kind::Func(f, pargs) => match e.kind() {
            Kind::Func(g, eargs) if f == g && pargs.len() == eargs.len() => {
                let save = b.len();
                for (x, y) in eargs.iter().zip(pargs) {
                    if !m(x, y, b) {
                        b.truncate(save);
                        return false;
                    }
                }
                true
            }
            _ => false,
        },
        _ => {
            let (ea, pa) = (e.args(), p.args());
            if e.class_name() != p.class_name() || ea.len() != pa.len() {
                return false;
            }
            let save = b.len();
            for (x, y) in ea.iter().zip(&pa) {
                if !m(x, y, b) {
                    b.truncate(save);
                    return false;
                }
            }
            true
        }
    }
}

/// Assign pattern items to expression items with backtracking. `rest`
/// combines the leftover expression items for the final lone wildcard.
fn assign(
    items: &[Expr],
    used: &mut Vec<bool>,
    pats: &[Expr],
    lone: &[Expr],
    neutral: &Expr,
    combine: &dyn Fn(&[Expr]) -> Expr,
    b: &mut Bindings,
) -> bool {
    if let Some((p0, prest)) = pats.split_first() {
        for i in 0..items.len() {
            if used[i] {
                continue;
            }
            let save = b.len();
            used[i] = true;
            if m(&items[i], p0, b) && assign(items, used, prest, lone, neutral, combine, b) {
                return true;
            }
            used[i] = false;
            b.truncate(save);
        }
        // the pattern item may match the neutral element (0 or 1)
        let save = b.len();
        if m(neutral, p0, b) && assign(items, used, prest, lone, neutral, combine, b) {
            return true;
        }
        b.truncate(save);
        return false;
    }
    let mut left: Vec<Expr> = items.iter().zip(used.iter()).filter(|(_, u)| !**u).map(|(x, _)| x.clone()).collect();
    // lone wildcards already bound by the other pattern items consume an
    // equal item (or stand for the neutral element)
    let mut unbound: Vec<Expr> = Vec::new();
    for w in lone {
        match lookup(b, w).cloned() {
            Some(v) => {
                if let Some(i) = left.iter().position(|x| *x == v) {
                    left.remove(i);
                } else if v != *neutral {
                    return false;
                }
            }
            None => unbound.push(w.clone()),
        }
    }
    match unbound.split_last() {
        None => left.is_empty(),
        Some((last, firsts)) => {
            // all but the last free wildcard take one item each
            if firsts.is_empty() {
                return m(&combine(&left), last, b);
            }
            let save = b.len();
            for i in 0..left.len() {
                if m(&left[i], &firsts[0], b) {
                    let mut rest = left.clone();
                    rest.remove(i);
                    let mut used2 = vec![false; rest.len()];
                    if assign(&rest, &mut used2, &[], &unbound[1..], neutral, combine, b) {
                        return true;
                    }
                }
                b.truncate(save);
            }
            false
        }
    }
}

fn match_add(e: &Expr, p: &Expr, b: &mut Bindings) -> bool {
    let items: Vec<Expr> = match e.kind() {
        Kind::Add(..) => e.args(),
        _ => vec![e.clone()],
    };
    let pargs = p.args();
    let (fixed, wilds): (Vec<Expr>, Vec<Expr>) = pargs.into_iter().partition(|t| !has_wild(t));
    // fixed pattern terms must be present (numbers are subtracted)
    let mut remaining = add(&items);
    for f in &fixed {
        remaining = sub(&remaining, f);
    }
    let items: Vec<Expr> = match remaining.kind() {
        Kind::Add(..) => remaining.args(),
        _ if remaining.is_zero() => vec![],
        _ => vec![remaining.clone()],
    };
    let (lone, pats): (Vec<Expr>, Vec<Expr>) = wilds.into_iter().partition(is_wild);
    let mut used = vec![false; items.len()];
    let save = b.len();
    let ok = assign(&items, &mut used, &pats, &lone, &Expr::zero(), &|v: &[Expr]| add(v), b);
    if !ok {
        b.truncate(save);
    }
    ok
}

fn match_mul(e: &Expr, p: &Expr, b: &mut Bindings) -> bool {
    let pargs = p.args();
    let (fixed, wilds): (Vec<Expr>, Vec<Expr>) = pargs.into_iter().partition(|t| !has_wild(t));
    // divide out the fixed factors; they must actually be present
    let fixed_prod = mul(&fixed);
    let q = div(e, &fixed_prod);
    if !fixed.iter().all(|f| f.is_number()) {
        // a fixed symbolic factor must not reappear inverted in q
        for f in &fixed {
            if f.is_number() {
                continue;
            }
            let (fb, _) = crate::ops::as_base_exp(f);
            if has(&q, &fb) && !has(e, &fb) {
                return false;
            }
            if has(&q, &crate::ops::inv(f)) {
                return false;
            }
        }
        // q must not have negative powers of fixed bases that e lacks
        if let Kind::Mul(_, facs) = q.kind() {
            for (qb, qe) in facs {
                if fixed.iter().any(|f| crate::ops::as_base_exp(f).0 == *qb)
                    && qe.as_number().is_some_and(|n| n.is_negative())
                {
                    return false;
                }
            }
        } else if let Kind::Pow(qb, qe) = q.kind() {
            if fixed.iter().any(|f| crate::ops::as_base_exp(f).0 == *qb)
                && qe.as_number().is_some_and(|n| n.is_negative())
            {
                return false;
            }
        }
    }
    let items: Vec<Expr> = match q.kind() {
        Kind::Mul(..) => q.args(),
        _ if q.is_one() => vec![],
        _ => vec![q.clone()],
    };
    let (lone, pats): (Vec<Expr>, Vec<Expr>) = wilds.into_iter().partition(is_wild);
    let mut used = vec![false; items.len()];
    let save = b.len();
    let ok = assign(&items, &mut used, &pats, &lone, &Expr::one(), &|v: &[Expr]| mul(v), b);
    if !ok {
        b.truncate(save);
    }
    ok
}

/// Replace every sub-expression matching `pattern` (bottom-up) by
/// `f(bindings)`.
pub fn replace_with(e: &Expr, pattern: &Expr, f: &mut dyn FnMut(&FxHashMap<Expr, Expr>, &Expr) -> Option<Expr>) -> Expr {
    let rebuilt = if e.is_atom() {
        e.clone()
    } else {
        map_children(e, &mut |c| replace_with(c, pattern, f))
    };
    if let Some(bind) = pattern_match(&rebuilt, pattern) {
        if let Some(r) = f(&bind, &rebuilt) {
            return r;
        }
    }
    rebuilt
}

/// `expr.replace(pattern, template)`: substitute the bindings into the
/// template.
pub fn replace(e: &Expr, pattern: &Expr, template: &Expr) -> Expr {
    if !has_wild(pattern) {
        return crate::subs::subs1(e, pattern, template);
    }
    replace_with(e, pattern, &mut |bind, _| Some(xreplace(template, bind)))
}

/// All sub-expressions matching `pattern`.
pub fn find(e: &Expr, pattern: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    crate::traverse::preorder(e, &mut |s| {
        if pattern_match(s, pattern).is_some() && !out.contains(s) {
            out.push(s.clone());
        }
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    fn mm(e: &str, p: &str) -> String {
        let a = Expr::wild("a", vec![]);
        let b = Expr::wild("b", vec![]);
        let x = Expr::symbol("x");
        let ax = Expr::wild("ax", vec![x.clone()]);
        let bx = Expr::wild("bx", vec![x.clone()]);
        let env = [("a", a), ("b", b), ("ax", ax), ("bx", bx)];
        let mut pe = parse_str(p).unwrap();
        for (n, w) in &env {
            pe = subs(&pe, &[(Expr::symbol(n), w.clone())]);
        }
        let r = pattern_match(&parse_str(e).unwrap(), &pe);
        match r {
            None => "None".into(),
            Some(m) => {
                let mut v: Vec<String> = m.iter().map(|(k, v)| format!("{k}: {v}")).collect();
                v.sort();
                v.join(", ")
            }
        }
    }

    #[test]
    fn matching() {
        assert_eq!(mm("2*x+3", "a*x+b"), "a_: 2, b_: 3");
        assert_eq!(mm("x", "a*x+b"), "a_: 1, b_: 0");
        assert_eq!(mm("2*x+3", "ax*x+bx"), "ax_: 2, bx_: 3");
        assert_eq!(mm("x**2", "x**a"), "a_: 2");
        assert_eq!(mm("x", "x**a"), "a_: 1");
        assert_eq!(mm("sin(x)**2", "sin(a)**b"), "a_: x, b_: 2");
        assert_eq!(mm("3*x*y", "a*x"), "a_: 3*y");
        assert_eq!(mm("3*x*y", "ax*x"), "ax_: 3*y");
        assert_eq!(mm("exp(2*x)", "exp(ax*x)"), "ax_: 2");
        assert_eq!(mm("2*x**3+x", "ax*x**3+bx*x"), "ax_: 2, bx_: 1");
        assert_eq!(mm("y*x+z", "ax*x+bx"), "ax_: y, bx_: z");
        assert_eq!(mm("x+y", "x+a"), "a_: y");
        assert_eq!(mm("x", "y"), "None");
        assert_eq!(mm("2*x", "a"), "a_: 2*x");
        assert_eq!(mm("y", "ax*x"), "None");
    }
}
