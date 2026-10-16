//! Substitution: `subs` (algebraic, sequential) and `xreplace` (exact).

use num_integer::Integer;
use num_traits::ToPrimitive;
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
    // substituting a differentiation variable: keep an unevaluated Subs
    if let Kind::Derivative(_, vars) = e.kind() {
        if vars.iter().any(|(v, _)| v == old) {
            return crate::functions::undefined("Subs", vec![e.clone(), old.clone(), new.clone()]);
        }
    }
    let r = match (e.kind(), old.kind()) {
        (Kind::Add(..), Kind::Add(..)) => match add_contains(e, old) {
            Some((k, rest)) => {
                let rest = s1(&rest, old, new, cache);
                add(&[mul2(&Expr::num(k), new), rest])
            }
            None => subs_children(e, old, new, cache),
        },
        (Kind::Mul(..), Kind::Mul(..)) => match mul_contains(e, old) {
            Some((count, rest)) => {
                let rest = s1(&rest, old, new, cache);
                mul2(&pow(new, &Expr::int(count)), &rest)
            }
            None => {
                // Canonical products store denominator factors with negative
                // powers. Give a positive product in that denominator the
                // same opportunity to match as a product in the numerator.
                let (numerator, denominator) = crate::simplify::numer_denom(e);
                if !denominator.is_one() {
                    let n = s1(&numerator, old, new, cache);
                    let d = s1(&denominator, old, new, cache);
                    if n != numerator || d != denominator {
                        div(&n, &d)
                    } else {
                        subs_children(e, old, new, cache)
                    }
                } else {
                    subs_children(e, old, new, cache)
                }
            }
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
            subs_children(e, old, new, cache)
        }
        _ => subs_children(e, old, new, cache),
    };
    cache.insert(e.clone(), r.clone());
    r
}

/// Substitution sees the public expression tree, including coefficients
/// in Add terms and whole powers in Mul factors. map_children deliberately
/// visits the compact internal representation instead, which would miss
/// patterns such as 2*x in 2*x+y, or x**2 in x**2*y.
fn subs_children(e: &Expr, old: &Expr, new: &Expr, cache: &mut FxHashMap<Expr, Expr>) -> Expr {
    if e.is_add() || e.is_mul() {
        let args = e.args();
        let replaced: Vec<Expr> = args.iter().map(|c| s1(c, old, new, cache)).collect();
        if replaced == args {
            return e.clone();
        }
        return if e.is_add() { add(&replaced) } else { mul(&replaced) };
    }
    map_children(e, &mut |c| s1(c, old, new, cache))
}

/// Match a subset of additive terms, including the opposite sign.
/// SymPy matches the coefficients of each term literally: `2*x + 2*y`
/// does not contain `x + y`. Rational constant offsets are independent.
pub fn add_contains(e: &Expr, s: &Expr) -> Option<(Number, Expr)> {
    let (Kind::Add(ce, te), Kind::Add(cs, ts)) = (e.kind(), s.kind()) else {
        return None;
    };
    if ts.len() > te.len() || (!(ce.is_rational() && cs.is_rational()) && ce != cs) {
        return None;
    }
    let lookup: FxHashMap<&Expr, &Number> = te.iter().map(|(t, k)| (t, k)).collect();
    let mut sign: Option<Number> = None;
    for (t, k) in ts {
        let ke = lookup.get(t)?;
        let r = ke.div(k)?;
        if !r.is_one() && !r.is_minus_one() {
            return None;
        }
        match &sign {
            None => sign = Some(r),
            Some(q) if *q == r => {}
            _ => return None,
        }
    }
    let k = sign?;
    let c_rest = ce.sub(&cs.mul(&k));
    let sub_terms: FxHashMap<&Expr, ()> = ts.iter().map(|(t, _)| (t, ())).collect();
    let rest = te.iter().filter(|(t, _)| !sub_terms.contains_key(t)).cloned().collect();
    Some((k, add_from_terms(c_rest, rest)))
}

/// Extract an integer number of copies of a rational exponent, preserving
/// SymPy's extractive divisibility rule for fractional powers.
fn exponent_copies(available: &Expr, wanted: &Expr) -> Option<i64> {
    let a = available.as_number()?.to_ratio()?;
    let b = wanted.as_number()?.to_ratio()?;
    if a.is_integer() && b.is_integer() {
        return (a / b).to_integer().to_i64().filter(|&n| n > 0);
    }
    let aq = a.denom();
    let bq = b.denom();
    if !aq.is_multiple_of(bq) && !bq.is_multiple_of(aq) {
        return None;
    }
    (a / b).to_integer().to_i64().filter(|&n| n > 0)
}

/// Multiplicity of an exact rational factor in a rational coefficient.
/// Numerator and denominator multiplicities are counted independently, so
/// e.g. `3/2` contains one `1/2` and contains no `1/4`.
fn coefficient_copies(available: &Number, wanted: &Number, limit: i64) -> Option<i64> {
    let a = available.to_ratio()?;
    let b = wanted.to_ratio()?;
    let mut numerator = a.numer().clone();
    let mut denominator = a.denom().clone();
    let mut count = 0;
    while count < limit
        && numerator.is_multiple_of(b.numer())
        && denominator.is_multiple_of(b.denom())
    {
        numerator /= b.numer();
        denominator /= b.denom();
        count += 1;
    }
    Some(count)
}

/// Extract as many copies of `s` as occur in `e`, allowing larger powers.
fn mul_contains(e: &Expr, s: &Expr) -> Option<(i64, Expr)> {
    let (Kind::Mul(ce, fe), Kind::Mul(cs, fs)) = (e.kind(), s.kind()) else {
        return None;
    };
    if fs.len() > fe.len() || (cs.is_negative() && !ce.is_negative()) {
        return None;
    }
    let lookup: FxHashMap<&Expr, &Expr> = fe.iter().map(|(b, x)| (b, x)).collect();
    let mut count = i64::MAX;
    for (b, x) in fs {
        count = count.min(exponent_copies(lookup.get(b)?, x)?);
    }
    if count == i64::MAX {
        return None;
    }
    let abs_old = cs.abs();
    if !abs_old.is_one() {
        count = coefficient_copies(ce, &abs_old, count)?;
        if count == 0 {
            return None;
        }
    }
    // Respect the sign when a negative pattern is extracted repeatedly.
    let residual = ce.div(&cs.pow_i64(count)?)?;
    let wanted: FxHashMap<&Expr, &Expr> = fs.iter().map(|(b, x)| (b, x)).collect();
    let mut rest = vec![Expr::num(residual)];
    for (b, x) in fe {
        let exponent = match wanted.get(b) {
            Some(old_x) => crate::ops::sub(x, &mul2(&Expr::int(count), old_x)),
            None => x.clone(),
        };
        rest.push(mul_factor(b, &exponent));
    }
    Some((count, mul(&rest)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn s(text: &str) -> Expr {
        parse_str(text).unwrap()
    }

    #[test]
    fn additive_subsets_and_constant_offsets() {
        for (expr, old, expected) in [
            ("x + y + z", "x + y", "w + z"),
            ("x - y + z", "y - x", "z - w"),
            ("x + y", "x + y + 1", "w - 1"),
            ("x + y + 2", "x + y - 1", "w + 3"),
            ("2*x + 2*y + z", "x + y", "2*x + 2*y + z"),
        ] {
            assert_eq!(subs1(&s(expr), &s(old), &s("w")), s(expected), "{expr} subs {old}");
        }
    }

    #[test]
    fn multiplicative_subsets_and_powers() {
        for (expr, old, expected) in [
            ("x**2*y", "x**2", "w*y"),
            ("x**4*y", "x**2", "w**2*y"),
            ("2*x + y", "2*x", "w + y"),
            ("2*x*y", "2*x", "w*y"),
            ("x**2*y*z", "x*y", "w*x*z"),
            ("x**3*y**2*z", "x*y", "w**2*x*z"),
            ("6*x**3*y**2", "2*x*y", "3*w*x**2*y"),
            ("8*x**3*y**3", "2*x*y", "w**3"),
            ("x**2*y**3", "2*x*y", "x**2*y**3"),
            ("3*x**3*y**3/2", "x*y/2", "3*w*x**2*y**2"),
            ("3*x**3*y**3/2", "x*y/4", "3*x**3*y**3/2"),
            ("2*x**3*y**3", "x*y/2", "2*x**3*y**3"),
            ("-x*y*z", "-x*y", "w*z"),
            ("x*y*z", "-x*y", "x*y*z"),
            ("z/(x*y)", "x*y", "z/w"),
            ("z/(x**3*y**2)", "x*y", "z/(w**2*x)"),
            ("x*z/y", "x/y", "w*z"),
        ] {
            assert_eq!(subs1(&s(expr), &s(old), &s("w")), s(expected), "{expr} subs {old}");
        }
    }
}
