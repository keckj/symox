//! Substitution: `subs` (algebraic, sequential) and `xreplace` (exact).

use rustc_hash::FxHashMap;

use crate::expr::{Expr, Kind};
use crate::number::Number;
use crate::ops::{add, add_from_terms, div, mul, mul2, mul_factor, pow};
use crate::traverse::map_children;

/// Simultaneous exact structural replacement.
pub fn xreplace(e: &Expr, map: &FxHashMap<Expr, Expr>) -> Expr {
    if map.is_empty() {
        return e.clone();
    }
    let mask: u64 = map.keys().fold(0, |m, k| m | k.mask());
    let all_masked = map.keys().all(|k| k.mask() != 0);
    let mut cache: FxHashMap<Expr, Expr> = FxHashMap::default();
    xr(e, map, mask, all_masked, &mut cache)
}

fn xr(
    e: &Expr,
    map: &FxHashMap<Expr, Expr>,
    mask: u64,
    all_masked: bool,
    cache: &mut FxHashMap<Expr, Expr>,
) -> Expr {
    if let Some(v) = map.get(e) {
        return v.clone();
    }
    if e.is_atom() || (all_masked && e.mask() & mask == 0) {
        return e.clone();
    }
    if let Some(v) = cache.get(e) {
        return v.clone();
    }
    let r = map_children(e, &mut |c| xr(c, map, mask, all_masked, cache));
    cache.insert(e.clone(), r.clone());
    r
}

/// SymPy's `subs`: apply each `(old, new)` pair in turn.
pub fn subs(e: &Expr, pairs: &[(Expr, Expr)]) -> Expr {
    let mut r = e.clone();
    for (old, new) in pairs {
        r = subs1(&r, old, new);
    }
    r
}

/// Substitute a single `old -> new` (with algebraic matching of sums,
/// products and powers).
pub fn subs1(e: &Expr, old: &Expr, new: &Expr) -> Expr {
    let mut cache = FxHashMap::default();
    s1(e, old, new, &mut cache)
}

fn s1(e: &Expr, old: &Expr, new: &Expr, cache: &mut FxHashMap<Expr, Expr>) -> Expr {
    if e == old {
        return new.clone();
    }
    if e.is_atom() || (old.mask() & !e.mask() != 0) {
        return e.clone();
    }
    if let Some(v) = cache.get(e) {
        return v.clone();
    }
    let r = match (e.kind(), old.kind()) {
        (Kind::Add(..), Kind::Add(..)) => match add_contains(e, old) {
            Some((k, rest)) => {
                let rest = s1(&rest, old, new, cache);
                add(&[mul2(&Expr::num(k), new), rest])
            }
            None => map_children(e, &mut |c| s1(c, old, new, cache)),
        },
        (Kind::Mul(..), Kind::Mul(..)) => match mul_contains(e, old) {
            Some(rest) => {
                let rest = s1(&rest, old, new, cache);
                mul2(new, &rest)
            }
            None => map_children(e, &mut |c| s1(c, old, new, cache)),
        },
        (Kind::Pow(b, x), _) => {
            // x**4 subs x**2 -> y  gives y**2;  x**3 subs x -> y gives y**3
            let (ob, ox) = crate::ops::as_base_exp(old);
            if *b == ob && !old.is_symbol() {
                let ratio = div(x, &ox);
                if ratio.is_integer_number() {
                    return pow(new, &ratio);
                }
            }
            map_children(e, &mut |c| s1(c, old, new, cache))
        }
        _ => map_children(e, &mut |c| s1(c, old, new, cache)),
    };
    cache.insert(e.clone(), r.clone());
    r
}

/// If `e = k*s + rest` term-wise (with every term of `s` present in `e`
/// with the same ratio `k`), return `(k, rest)`.
pub fn add_contains(e: &Expr, s: &Expr) -> Option<(Number, Expr)> {
    let (Kind::Add(ce, te), Kind::Add(cs, ts)) = (e.kind(), s.kind()) else {
        return None;
    };
    if ts.len() > te.len() {
        return None;
    }
    let lookup: FxHashMap<&Expr, &Number> = te.iter().map(|(t, k)| (t, k)).collect();
    let mut ratio: Option<Number> = None;
    for (t, k) in ts {
        let ke = lookup.get(t)?;
        let r = ke.div(k)?;
        match &ratio {
            None => ratio = Some(r),
            Some(q) if *q == r => {}
            _ => return None,
        }
    }
    let k = ratio?;
    let c_rest = ce.sub(&cs.mul(&k));
    if !cs.is_zero() && ce.is_zero() {
        return None;
    }
    let sub_terms: FxHashMap<&Expr, ()> = ts.iter().map(|(t, _)| (t, ())).collect();
    let rest: Vec<(Expr, Number)> = te
        .iter()
        .filter(|(t, _)| !sub_terms.contains_key(t))
        .cloned()
        .collect();
    Some((k, add_from_terms(c_rest, rest)))
}

/// If all factors of `s` appear in `e` (same exponents), return `e / s`.
fn mul_contains(e: &Expr, s: &Expr) -> Option<Expr> {
    let (Kind::Mul(ce, fe), Kind::Mul(cs, fs)) = (e.kind(), s.kind()) else {
        return None;
    };
    if fs.len() > fe.len() {
        return None;
    }
    for f in fs {
        if !fe.contains(f) {
            return None;
        }
    }
    let k = ce.div(cs)?;
    let mut v: Vec<Expr> = vec![Expr::num(k)];
    for f in fe {
        if !fs.contains(f) {
            v.push(mul_factor(&f.0, &f.1));
        }
    }
    Some(mul(&v))
}
