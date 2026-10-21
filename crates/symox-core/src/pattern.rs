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
pub type MatchError = crate::symbol::WildcardError;
#[derive(Default)]
struct Context { error: Option<MatchError> }

fn is_wild(e: &Expr) -> bool {
    e.as_sym().is_some_and(|s| s.is_wild())
}

/// Does the pattern contain wildcards?
pub fn has_wild(p: &Expr) -> bool {
    fn visit(p:&Expr,seen:&mut rustc_hash::FxHashSet<usize>)->bool{
        if p.mask()==0||!seen.insert(p.dag_id()){return false;}
        match p.kind(){
            Kind::Sym(s)=>s.is_wild(),
            Kind::Add(_,terms)=>terms.iter().any(|(t,_)|visit(t,seen)),
            Kind::Mul(_,factors)=>factors.iter().any(|(b,e)|visit(b,seen)||visit(e,seen)),
            _=>p.args().iter().any(|arg|visit(arg,seen)),
        }
    }
    visit(p,&mut rustc_hash::FxHashSet::default())
}

fn lookup<'a>(b: &'a Bindings, w: &Expr) -> Option<&'a Expr> {
    b.iter().find(|(k, _)| k == w).map(|(_, v)| v)
}

fn bind(w: &Expr, e: &Expr, b: &mut Bindings, context: &mut Context) -> bool {
    if let Some(v) = lookup(b, w) {
        return v == e;
    }
    let s = w.as_sym().unwrap();
    if s.exclude().iter().any(|x| has(e, x)) {
        return false;
    }
    if let Some(properties)=s.wild_properties() {
        match properties.evaluate(e) {
            Ok(true)=>{}, Ok(false)=>return false,
            Err(error)=>{context.error=Some(error);return false;},
        }
    }
    b.push((w.clone(), e.clone()));
    true
}

/// Match `e` against pattern `p`; returns the wildcard bindings.
pub fn pattern_match(e: &Expr, p: &Expr) -> Result<Option<FxHashMap<Expr, Expr>>,MatchError> {
    let mut b = Bindings::new();
    let mut context=Context::default();
    let matched=m(e,p,&mut b,&mut context);
    if let Some(error)=context.error { return Err(error); }
    Ok(matched.then(||b.into_iter().collect()))
}

fn m(e: &Expr, p: &Expr, b: &mut Bindings, context: &mut Context) -> bool {
    if context.error.is_some(){return false;}
    if !has_wild(p) {
        return e == p;
    }
    if is_wild(p) {
        return bind(p, e, b, context);
    }
    match p.kind() {
        Kind::Add(..) => match_add(e, p, b, context),
        Kind::Mul(..) => match_mul(e, p, b, context),
        Kind::Pow(pb, pe) => {
            if let Kind::Pow(eb, ee) = e.kind() {
                let save = b.len();
                if m(eb, pb, b, context) && m(ee, pe, b, context) {
                    return true;
                }
                b.truncate(save);
            }
            // x matches x**a with a = 1
            let save = b.len();
            if (has_wild(pe)||pe.is_one()) && m(e, pb, b, context) && m(&Expr::one(), pe, b, context) {
                return true;
            }
            b.truncate(save);
            // Canonical numeric reciprocals are numbers rather than Pow nodes.
            // Zero has no reciprocal match in the reference convention.
            if pe.is_minus_one()&&!e.is_zero(){
                let inverse=crate::ops::inv(e);
                if m(&inverse,pb,b,context){return true;}
                b.truncate(save);
            }
            false
        }
        Kind::Func(f, pargs) => match e.kind() {
            Kind::Func(g, eargs) if f == g && pargs.len() == eargs.len() => {
                let save = b.len();
                for (x, y) in eargs.iter().zip(pargs) {
                    if !m(x, y, b, context) {
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
                if !m(x, y, b, context) {
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
    context: &mut Context,
    scan: &mut ScanCache,
) -> bool {
    if context.error.is_some(){return false;}
    if let Some((p0, prest)) = pats.split_first() {
        // Terms that cannot match `p0` under the current bindings are not
        // visited; the others are visited in their original order.
        let candidates = if !scan.unfiltered && neutral.is_zero() && items.len() >= DETERMINED_SCAN_MIN {
            match determined_target(p0, b) {
                Some(target) => Some(determined_candidates(items, &target, scan)),
                None => factor_candidates(items, p0, pats.len(), scan, b),
            }
        } else { None };
        if let Some(candidates) = &candidates {
            if check_determined_filter() {
                for i in 0..items.len() {
                    if used[i] || candidates.binary_search(&i).is_ok() { continue; }
                    let save = b.len();
                    assert!(!m(&items[i], p0, b, context), "scan filter skipped a match: {} vs {}", items[i], p0);
                    b.truncate(save);
                }
            }
        }
        let mut cursor = 0;
        loop {
            let i = match &candidates {
                Some(candidates) => match candidates.get(cursor) { Some(&i) => i, None => break },
                None if cursor < items.len() => cursor,
                None => break,
            };
            cursor += 1;
            if context.error.is_some(){return false;}
            if used[i] {
                continue;
            }
            let save = b.len();
            used[i] = true;
            if m(&items[i], p0, b, context) && assign(items, used, prest, lone, neutral, combine, b, context, scan) {
                return true;
            }
            used[i] = false;
            b.truncate(save);
        }
        // the pattern item may match the neutral element (0 or 1)
        let save = b.len();
        if m(neutral, p0, b, context) && assign(items, used, prest, lone, neutral, combine, b, context, scan) {
            return true;
        }
        b.truncate(save);
        return false;
    }
    // A bound lone wildcard without an equal unused item (and other than
    // the neutral element) fails below whatever the other items are.
    if !scan.unfiltered && items.len() >= DETERMINED_SCAN_MIN && lone.iter().any(|w| lookup(b, w).is_some()) {
        let values = scan.values(items);
        let missing = lone.iter().any(|w| lookup(b, w).is_some_and(|v| {
            v != neutral && !values.get(v).is_some_and(|positions| positions.iter().any(|&j| !used[j]))
        }));
        if missing {
            if check_determined_filter() {
                let save = b.len();
                assert!(!assign(items, used, pats, lone, neutral, combine, b, context, &mut ScanCache { unfiltered: true, ..ScanCache::default() }),
                    "lone wildcard filter skipped a match");
                b.truncate(save);
            }
            return false;
        }
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
                return m(&combine(&left), last, b, context);
            }
            let save = b.len();
            for i in 0..left.len() {
                if context.error.is_some(){return false;}
                if m(&left[i], &firsts[0], b, context) {
                    let mut rest = left.clone();
                    rest.remove(i);
                    let mut used2 = vec![false; rest.len()];
                    if assign(&rest, &mut used2, &[], &unbound[1..], neutral, combine, b, context, &mut ScanCache::default()) {
                        return true;
                    }
                }
                b.truncate(save);
            }
            if m(neutral,&firsts[0],b,context) {
                let mut used=vec![false;left.len()];
                if assign(&left,&mut used,&[],&unbound[1..],neutral,combine,b,context,&mut ScanCache::default()){return true;}
            }
            b.truncate(save);
            // A constrained wildcard can need a group of terms, rather than
            // a singleton. Enumerate subsets lazily without a size cutoff.
            for size in 2..=left.len(){
                let mut chosen=Vec::with_capacity(size);
                if assign_group(&left,&unbound,size,0,&mut chosen,neutral,combine,b,context){return true;}
                if context.error.is_some(){return false;}
            }
            false
        }
    }
}

/// Sums (or products, for bound lone wildcards) at least this long use
/// the indexes of `ScanCache` instead of visiting every term.
const DETERMINED_SCAN_MIN: usize = 8;

fn check_determined_filter() -> bool {
    static CHECK: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *CHECK.get_or_init(|| std::env::var_os("SYMOX_CHECK_DETERMINED_FILTER").is_some())
}

/// Exact arithmetic without held, bound, host or floating-point nodes: the
/// rebuilding arithmetic of the matcher (quotients by fixed factors,
/// differences with fixed terms, reciprocals) round-trips on such terms.
fn plain_term(e: &Expr) -> bool {
    !e.has_held_arithmetic() && !e.has_binding_identity() && !e.has_host() && e.float_info().max_precision.is_none()
}

/// The term a sum pattern term rebuilds to when the current bindings bind
/// every wildcard in it. Matching such a term binds nothing new and
/// succeeds only on an equal term: each matcher step accepts a sub-term
/// only when the corresponding pattern part rebuilds to it, so a plain term
/// that differs from the target cannot match and is not visited. The
/// equal terms are found by value, so the cost does not grow with the
/// length of the sum for every binding of the earlier pattern terms.
fn determined_target(p: &Expr, b: &Bindings) -> Option<Expr> {
    if !has_wild(p) {
        return None;
    }
    let mut map: FxHashMap<Expr, Expr> = FxHashMap::default();
    let mut complete = true;
    crate::traverse::preorder(p, &mut |node| {
        if complete && is_wild(node) {
            match lookup(b, node) {
                Some(v) => {
                    map.insert(node.clone(), v.clone());
                }
                None => complete = false,
            }
        }
    });
    if !complete {
        return None;
    }
    let target = xreplace(p, &map);
    plain_term(&target).then_some(target)
}

/// Per-search indexes of the sum terms, one for each position of the
/// pattern terms (keyed by the number of pattern terms still to place, which
/// identifies the position within one search over the same terms).
#[derive(Default)]
struct ScanCache {
    factors: FxHashMap<usize, Option<FactorIndex>>,
    /// Positions of the terms by value.
    values: Option<FxHashMap<Expr, Vec<usize>>>,
    /// Positions of the terms that are not `plain_term`.
    irregular: Option<Vec<usize>>,
    /// Visit every term (used to check the filters).
    unfiltered: bool,
}

impl ScanCache {
    fn values(&mut self, items: &[Expr]) -> &FxHashMap<Expr, Vec<usize>> {
        self.values.get_or_insert_with(|| {
            let mut values: FxHashMap<Expr, Vec<usize>> = FxHashMap::default();
            for (j, item) in items.iter().enumerate() {
                values.entry(item.clone()).or_default().push(j);
            }
            values
        })
    }
}

/// The terms (ascending) that can match a pattern term whose wildcards are
/// all bound: those equal to its rebuilt `target` and those that are not
/// plain (see `determined_target`).
fn determined_candidates(items: &[Expr], target: &Expr, scan: &mut ScanCache) -> Vec<usize> {
    let mut candidates = scan.values(items).get(target).cloned().unwrap_or_default();
    let irregular = scan.irregular.get_or_insert_with(|| (0..items.len()).filter(|&j| !plain_term(&items[j])).collect());
    if !irregular.is_empty() {
        candidates.extend(irregular.iter().copied());
        candidates.sort_unstable();
        candidates.dedup();
    }
    candidates
}

/// For a product pattern term with bare wildcard factors: the factors each
/// sum term leaves once divided by the fixed factors, exactly as
/// `match_mul` computes them.
struct FactorIndex {
    lone: Vec<Expr>,
    by_factor: FxHashMap<Expr, Vec<usize>>,
}

fn factor_index(items: &[Expr], p: &Expr) -> Option<FactorIndex> {
    // Wildcard properties are evaluated (and can fail with an error) while
    // matching: such patterns are scanned term by term.
    if !matches!(p.kind(), Kind::Mul(..)) || p.has_wild_properties() {
        return None;
    }
    let (fixed, wilds): (Vec<Expr>, Vec<Expr>) = p.args().into_iter().partition(|t| !has_wild(t));
    let lone: Vec<Expr> = wilds.into_iter().filter(is_wild).collect();
    if lone.is_empty() {
        return None;
    }
    let fixed_prod = mul(&fixed);
    let mut by_factor: FxHashMap<Expr, Vec<usize>> = FxHashMap::default();
    for (j, item) in items.iter().enumerate() {
        let q = div(item, &fixed_prod);
        let factors = match q.kind() {
            Kind::Mul(..) => q.args(),
            _ if q.is_one() => vec![],
            _ => vec![q],
        };
        for f in factors {
            let list = by_factor.entry(f).or_default();
            if list.last() != Some(&j) {
                list.push(j);
            }
        }
    }
    Some(FactorIndex { lone, by_factor })
}

/// The sum terms (ascending) that can match the product pattern term `p`
/// under the current bindings, when a bare wildcard factor of `p` is bound
/// to a value other than 1. `match_mul` matches the other pattern factors
/// against the factors of the term divided by the fixed factors and then
/// requires a bound bare wildcard to consume an equal leftover factor (or
/// to be 1), so a term lacking that factor cannot match. Returns `None`
/// when no such filter applies.
fn factor_candidates(items: &[Expr], p: &Expr, key: usize, scan: &mut ScanCache, b: &Bindings) -> Option<Vec<usize>> {
    let index = scan.factors.entry(key).or_insert_with(|| factor_index(items, p)).as_ref()?;
    let mut best: Option<&[usize]> = None;
    for w in &index.lone {
        if let Some(v) = lookup(b, w) {
            if v.is_one() {
                continue;
            }
            let list = index.by_factor.get(v).map(Vec::as_slice).unwrap_or(&[]);
            if best.is_none_or(|current| list.len() < current.len()) {
                best = Some(list);
            }
        }
    }
    best.map(<[usize]>::to_vec)
}

fn assign_group(items:&[Expr],wilds:&[Expr],remaining:usize,start:usize,chosen:&mut Vec<usize>,neutral:&Expr,
    combine:&dyn Fn(&[Expr])->Expr,b:&mut Bindings,context:&mut Context)->bool {
    if context.error.is_some(){return false;}
    if remaining==0 {
        let selected=chosen.iter().map(|&i|items[i].clone()).collect::<Vec<_>>();
        let save=b.len();
        if m(&combine(&selected),&wilds[0],b,context){
            let rest=items.iter().enumerate().filter(|(i,_)|chosen.binary_search(i).is_err()).map(|(_,e)|e.clone()).collect::<Vec<_>>();
            let mut used=vec![false;rest.len()];
            if assign(&rest,&mut used,&[],&wilds[1..],neutral,combine,b,context,&mut ScanCache::default()){return true;}
        }
        b.truncate(save);return false;
    }
    for index in start..=items.len()-remaining {
        chosen.push(index);
        if assign_group(items,wilds,remaining-1,index+1,chosen,neutral,combine,b,context){return true;}
        chosen.pop();
        if context.error.is_some(){return false;}
    }
    false
}

fn match_add(e: &Expr, p: &Expr, b: &mut Bindings, context: &mut Context) -> bool {
    if context.error.is_some(){return false;}
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
    let ok = assign(&items, &mut used, &pats, &lone, &Expr::zero(), &|v: &[Expr]| add(v), b, context, &mut ScanCache::default());
    if !ok {
        b.truncate(save);
    }
    ok
}

fn match_mul(e: &Expr, p: &Expr, b: &mut Bindings, context: &mut Context) -> bool {
    if context.error.is_some(){return false;}
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
    let ok = assign(&items, &mut used, &pats, &lone, &Expr::one(), &|v: &[Expr]| mul(v), b, context, &mut ScanCache::default());
    if !ok {
        b.truncate(save);
    }
    ok
}

/// Replace every sub-expression matching `pattern` (bottom-up) by
/// `f(bindings)`.
pub fn replace_with(e: &Expr, pattern: &Expr, f: &mut dyn FnMut(&FxHashMap<Expr, Expr>, &Expr) -> Result<Option<Expr>,MatchError>) -> Result<Expr,MatchError> {
    let mut error=None;
    let rebuilt = if e.is_atom() {
        e.clone()
    } else {
        map_children(e, &mut |c| {
            if error.is_some(){return c.clone();}
            match replace_with(c,pattern,f){Ok(value)=>value,Err(value)=>{error=Some(value);c.clone()}}
        })
    };
    if let Some(error)=error{return Err(error);}
    if let Some(bind) = pattern_match(&rebuilt, pattern)? {
        if let Some(r) = f(&bind, &rebuilt)? {
            return Ok(r);
        }
    }
    Ok(rebuilt)
}

/// `expr.replace(pattern, template)`: substitute the bindings into the
/// template.
pub fn replace(e: &Expr, pattern: &Expr, template: &Expr) -> Result<Expr,MatchError> {
    if !has_wild(pattern) {
        return Ok(crate::subs::subs1(e, pattern, template));
    }
    replace_with(e, pattern, &mut |bind, _| Ok(Some(xreplace(template, bind))))
}

/// All sub-expressions matching `pattern`.
pub fn find(e: &Expr, pattern: &Expr) -> Result<Vec<Expr>,MatchError> {
    let mut out: Vec<Expr> = Vec::new();
    let mut pending=vec![e.clone()];
    let mut seen=rustc_hash::FxHashSet::default();
    while let Some(node)=pending.pop(){
        if !seen.insert(node.dag_id()){continue;}
        if pattern_match(&node,pattern)?.is_some()&&!out.contains(&node){out.push(node.clone());}
        pending.extend(node.args().into_iter().rev());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[derive(Debug)]
    struct Predicate {
        id: u64,
        check: fn(&Expr) -> Result<bool, MatchError>,
    }

    impl crate::symbol::WildcardProperties for Predicate {
        fn identity(&self) -> u64 { self.id }
        fn as_any(&self) -> &dyn std::any::Any { self }
        fn evaluate(&self, value: &Expr) -> Result<bool, MatchError> { (self.check)(value) }
    }

    fn constrained(name: &str, id: u64, check: fn(&Expr) -> Result<bool, MatchError>) -> Expr {
        Expr::from_kind(Kind::Sym(crate::symbol::Symbol::wild_with_properties(
            name, vec![], crate::assumption_facts::FactSet::EMPTY,
            std::sync::Arc::new(Predicate { id, check }),
        ).unwrap()))
    }

    #[test]
    fn constrained_commutative_partitions_and_neutral_elements() {
        let a = constrained("a", 1, |v| Ok(has(v, &Expr::symbol("x")) && has(v, &Expr::symbol("y"))));
        let b = constrained("b", 2, |v| Ok(v.is_number()));
        let grouped = parse_str("x+y+1").unwrap();
        let bindings = pattern_match(&grouped, &add(&[a.clone(), b.clone()])).unwrap().unwrap();
        assert_eq!(bindings[&a], parse_str("x+y").unwrap());
        assert_eq!(bindings[&b], Expr::one());
        let grouped = parse_str("3*x*y").unwrap();
        let bindings = pattern_match(&grouped, &mul(&[a.clone(), b.clone()])).unwrap().unwrap();
        assert_eq!(bindings[&a], parse_str("x*y").unwrap());
        assert_eq!(bindings[&b], Expr::int(3));

        let zero = constrained("a", 3, |v| Ok(v.is_zero()));
        let rest = Expr::wild("z", vec![]);
        let value = parse_str("x+y").unwrap();
        let bindings = pattern_match(&value, &add(&[zero.clone(), rest.clone()])).unwrap().unwrap();
        assert_eq!(bindings[&zero], Expr::zero());
        assert_eq!(bindings[&rest], value);
    }

    #[test]
    fn predicate_errors_survive_matching_and_traversal() {
        let a = constrained("a", 4, |_| Err(std::io::Error::other("predicate sentinel").into()));
        let value = parse_str("sin(x)+2").unwrap();
        assert_eq!(pattern_match(&value, &a).unwrap_err().to_string(), "predicate sentinel");
        assert_eq!(find(&value, &a).unwrap_err().to_string(), "predicate sentinel");
        assert_eq!(replace(&value, &a, &Expr::zero()).unwrap_err().to_string(), "predicate sentinel");
    }

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
        let r = pattern_match(&parse_str(e).unwrap(), &pe).unwrap();
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
