//! Set properties: membership, bounds, measure, boundary, relational form
//! and subset tests.

use std::cmp::Ordering;

use super::*;
use crate::expr::RelOp;
use crate::logic::{and, not, or};
use crate::relational::rel;

fn complex_inf(e: &Expr) -> bool {
    e.is_constant(Constant::ComplexInfinity) || e.is_constant(Constant::NaN)
}

/// Membership of `x` in `s`: a boolean expression (`true`/`false` when
/// decided).
pub fn contains(s: &Expr, x: &Expr) -> Expr {
    let Some((op, a)) = as_set(s) else {
        return mk(SetOp::Contains, vec![x.clone(), s.clone()]);
    };
    let t = || Expr::boolean(true);
    let f = || Expr::boolean(false);
    let unevaluated = || mk(SetOp::Contains, vec![x.clone(), s.clone()]);
    let is_real = crate::assumptions::is_real(x);
    match op {
        SetOp::EmptySet => f(),
        SetOp::UniversalSet => t(),
        SetOp::Interval | SetOp::Reals => {
            if is_set(x) || matches!(x.kind(), Kind::Bool(_)) || complex_inf(x) || is_real == Some(false) {
                return f();
            }
            if op == SetOp::Reals {
                if let Some(r) = is_real {
                    return Expr::boolean(r);
                }
            }
            as_relational(s, x)
        }
        SetOp::FiniteSet => {
            if a.iter().any(|e| e == x) {
                return t();
            }
            or(a.iter().map(|e| rel(RelOp::Eq, e, x)).collect())
        }
        SetOp::Union => or(a.iter().map(|p| contains(p, x)).collect()),
        SetOp::Intersection => and(a.iter().map(|p| contains(p, x)).collect()),
        SetOp::Complement => and(vec![contains(&a[0], x), not(&contains(&a[1], x))]),
        SetOp::Integers | SetOp::Naturals | SetOp::Naturals0 | SetOp::Rationals => {
            use crate::assumptions as asm;
            if is_set(x) || complex_inf(x) || is_infinite(x) {
                return f();
            }
            let base = if op == SetOp::Rationals { asm::is_rational(x) } else { asm::is_integer(x) };
            let r = match op {
                SetOp::Naturals => match base {
                    Some(true) => asm::is_positive(x),
                    b => b.filter(|b| !b),
                },
                SetOp::Naturals0 => match base {
                    Some(true) => asm::is_nonnegative(x),
                    b => b.filter(|b| !b),
                },
                _ => base,
            };
            match r {
                Some(b) => Expr::boolean(b),
                None => {
                    if op == SetOp::Rationals && (is_real == Some(false) || is_irrational_surd(x)) {
                        return f();
                    }
                    unevaluated()
                }
            }
        }
        SetOp::Complexes => {
            if crate::logic::is_number_like(x) && !is_infinite(x) && !complex_inf(x) {
                t()
            } else if is_infinite(x) || complex_inf(x) || is_set(x) {
                f()
            } else {
                unevaluated()
            }
        }
        SetOp::Range => {
            let integer = |e: &Expr| e.as_number().and_then(|n| n.to_bigint());
            if let (Some(st), Some(sp), Some(stp)) = (integer(&a[0]), integer(&a[1]), integer(&a[2])) {
                if let Some(v) = x.as_number() {
                    let Some(v) = v.to_bigint() else { return f() };
                    let zero = num_bigint::BigInt::from(0);
                    let inside = if stp > zero { v >= st && v < sp } else { v <= st && v > sp };
                    return Expr::boolean(inside && stp != zero && (v - st) % stp == zero);
                }
            }
            if let Some(v) = x.as_number() {
                if !v.is_integer() {
                    return f();
                }
                // half infinite ranges
                let lo_ok = if a[0].is_constant(Constant::NegInfinity) { Some(true) } else { cmp3(&a[0], x).map(|o| o != Ordering::Greater) };
                let hi_ok = if a[1].is_constant(Constant::Infinity) { Some(true) } else { lt3(x, &a[1]) };
                if a[2].is_one() {
                    if let (Some(l), Some(h)) = (lo_ok, hi_ok) {
                        return Expr::boolean(l && h);
                    }
                }
            }
            // reference `Range._contains` of an integer expression in a
            // numeric range: the bounds as relations when `(ref - x) % step`
            // is zero, false when it is a nonzero integer
            if !crate::logic::is_number_like(x) && crate::assumptions::is_integer(x) == Some(true) {
                if let Some(p) = range_progression(s) {
                    let Some(p) = p else { return f() };
                    let int = |v: &num_bigint::BigInt| Expr::num(crate::number::Number::from_bigint(v.clone()));
                    let residue = match (&p.lo, &p.hi) {
                        _ if num_traits::One::is_one(&p.step) => Some(num_bigint::BigInt::from(0)),
                        (Some(r), _) | (None, Some(r)) => integer_residue(&crate::ops::sub(&int(r), x), &p.step),
                        (None, None) => Some(num_bigint::BigInt::from(0)),
                    };
                    match (&p.lo, &p.hi) {
                        (None, None) => return t(),
                        (Some(l), Some(h)) if l == h => return rel(RelOp::Eq, x, &int(l)),
                        (lo, hi) => match residue {
                            Some(r) if num_traits::Zero::is_zero(&r) => {
                                return and(vec![
                                    lo.as_ref().map_or_else(t, |l| crate::relational::ge(x, &int(l))),
                                    hi.as_ref().map_or_else(t, |h| crate::relational::le(x, &int(h))),
                                ])
                            }
                            Some(_) => return f(),
                            None => {}
                        },
                    }
                }
            }
            unevaluated()
        }
        SetOp::ImageSet => {
            if let Some((variables, body)) = lambda_parts(&a[0]) {
                if variables.len()==1 && lambda_arity(&a[0])==Some(1) {
                    let variable = &variables[0];
                    let body = &body;
                    let coefficient = crate::diff::diff(body, variable);
                    if !crate::traverse::has(&coefficient, variable)
                        && crate::assumptions::is_zero(&coefficient) == Some(false) {
                        let offset = crate::subs::subs1(body, variable, &Expr::zero());
                        let residual = crate::ops::sub(body, &crate::ops::add2(&crate::ops::mul2(&coefficient, variable), &offset));
                        if crate::expand::expand(&residual).is_zero() {
                            let inverse = crate::simplify::simplify(&crate::ops::div(&crate::ops::sub(x, &offset), &coefficient));
                            // reference `ImageSet._contains`: only the
                            // emptiness of the solution set is decided, an
                            // undecided membership stays `Contains(x, self)`
                            let r = contains(&a[1], &inverse);
                            if matches!(r.kind(), Kind::Bool(_)) {
                                return r;
                            }
                            return unevaluated();
                        }
                    }
                }
            }
            unevaluated()
        }
        SetOp::ProductSet => {
            // reference `ProductSet._contains`: the conjunction of the
            // component memberships of a tuple of the right length
            if matches!(x.kind(), Kind::Sym(_)) {
                return unevaluated();
            }
            match crate::special::tuple_items(x) {
                Some(items) if items.len() == a.len() => and(a.iter().zip(items).map(|(p, e)| contains(p, e)).collect()),
                _ => f(),
            }
        }
        SetOp::ConditionSet => {
            let c = if crate::special::tuple_items(&a[0]).is_some() {
                let Some(values) = crate::special::tuple_items(x) else {
                    return if x.is_symbol() { unevaluated() } else { f() };
                };
                let Ok(predicate) = lambda_checked(&a[0], &a[1]) else { return unevaluated(); };
                let Ok(condition) = lambda_apply(&predicate, values) else { return f(); };
                condition
            } else { crate::subs::subs1(&a[1], &a[0], x) };
            and(vec![c, contains(&a[2], x)])
        }
        _ => unevaluated(),
    }
}

/// `e mod step` in `0..step` when every non-constant term of the expanded
/// `e` is an integer expression with a coefficient divisible by `step` (the
/// cases where the reference `Mod` evaluates to an integer); `None` otherwise.
fn integer_residue(e: &Expr, step: &num_bigint::BigInt) -> Option<num_bigint::BigInt> {
    use num_integer::Integer;
    let e = crate::expand::expand(e);
    let (konst, terms): (num_bigint::BigInt, Vec<(Expr, crate::number::Number)>) = match e.kind() {
        Kind::Num(n) => (n.to_bigint()?, vec![]),
        Kind::Add(c, t) => (c.to_bigint()?, t.clone()),
        _ => {
            let (c, t) = crate::ops::as_coeff_mul(&e);
            (num_bigint::BigInt::from(0), vec![(t, c)])
        }
    };
    for (term, c) in &terms {
        let c = c.to_bigint()?;
        if crate::assumptions::is_integer(term) != Some(true) || !num_traits::Zero::is_zero(&c.mod_floor(step)) {
            return None;
        }
    }
    Some(konst.mod_floor(step))
}

/// `q * n**(p/r)` with rational `q != 0` and an integer `n` whose
/// automatically simplified root is irrational.
fn is_irrational_surd(x: &Expr) -> bool {
    let surd = |b: &Expr, e: &Expr| {
        b.is_integer_number() && e.as_number().is_some_and(|n| n.is_rational() && !n.is_integer())
    };
    match x.kind() {
        Kind::Pow(b, e) => surd(b, e),
        Kind::Mul(c, facs) => !c.is_zero() && c.is_rational() && facs.len() == 1 && surd(&facs[0].0, &facs[0].1),
        _ => false,
    }
}

/// `S.as_relational(x)`.
pub fn as_relational(s: &Expr, x: &Expr) -> Expr {
    let Some((op, a)) = as_set(s) else {
        return mk(SetOp::Contains, vec![x.clone(), s.clone()]);
    };
    match op {
        SetOp::EmptySet => Expr::boolean(false),
        SetOp::UniversalSet => Expr::boolean(true),
        SetOp::Interval | SetOp::Reals => {
            let (st, en, l, r) = interval_parts(s).unwrap();
            let real_x = crate::assumptions::is_real(x) == Some(true);
            let left = if real_x && st.is_constant(Constant::NegInfinity) {
                Expr::boolean(true)
            } else {
                rel(if l { RelOp::Lt } else { RelOp::Le }, &st, x)
            };
            let right = if real_x && en.is_constant(Constant::Infinity) {
                Expr::boolean(true)
            } else {
                rel(if r { RelOp::Lt } else { RelOp::Le }, x, &en)
            };
            and(vec![left, right])
        }
        SetOp::FiniteSet => or(a.iter().map(|e| rel(RelOp::Eq, x, e)).collect()),
        SetOp::Union => {
            if a.len() == 2 {
                if let (Some((s1, e1, l1, r1)), Some((s2, e2, l2, r2))) = (interval_parts(&a[0]), interval_parts(&a[1])) {
                    if e1 == s2 && r1 && l2 {
                        let real_x = crate::assumptions::is_real(x) == Some(true);
                        let mincond = if real_x && s1.is_constant(Constant::NegInfinity) {
                            Expr::boolean(true)
                        } else {
                            rel(if l1 { RelOp::Gt } else { RelOp::Ge }, x, &s1)
                        };
                        let maxcond = if real_x && e2.is_constant(Constant::Infinity) {
                            Expr::boolean(true)
                        } else {
                            rel(if r2 { RelOp::Lt } else { RelOp::Le }, x, &e2)
                        };
                        let necond = rel(RelOp::Ne, x, &e1);
                        return and(vec![necond, mincond, maxcond]);
                    }
                }
            }
            or(a.iter().map(|p| as_relational(p, x)).collect())
        }
        SetOp::Intersection => and(a.iter().map(|p| as_relational(p, x)).collect()),
        SetOp::Complement => and(vec![as_relational(&a[0], x), not(&as_relational(&a[1], x))]),
        _ => contains(s, x),
    }
}

/// Infimum (`None` when not known).
pub fn inf(s: &Expr) -> Option<Expr> {
    let (op, a) = as_set(s)?;
    Some(match op {
        SetOp::Interval => a[0].clone(),
        SetOp::FiniteSet => crate::functions::min(a.to_vec()),
        SetOp::Union => crate::functions::min(a.iter().map(inf).collect::<Option<Vec<_>>>()?),
        SetOp::Reals | SetOp::Integers | SetOp::Rationals => Expr::neg_oo(),
        SetOp::Naturals => Expr::one(),
        SetOp::Naturals0 => Expr::zero(),
        // `Range._inf`: an empty or symbolic Range is decided by the
        // Python layer (the reference raises there).
        SetOp::Range if a[0] == a[1] || a.iter().any(|v| !crate::traverse::free_symbols(v).is_empty()) => return None,
        SetOp::Range => {
            if a[2].as_number().is_some_and(|n| n.is_positive()) {
                a[0].clone()
            } else {
                crate::ops::sub(&a[1], &a[2])
            }
        }
        SetOp::Intersection | SetOp::Complement => return None,
        _ => return None,
    })
}

/// Supremum (`None` when not known).
pub fn sup(s: &Expr) -> Option<Expr> {
    let (op, a) = as_set(s)?;
    Some(match op {
        SetOp::Interval => a[1].clone(),
        SetOp::FiniteSet => crate::functions::max(a.to_vec()),
        SetOp::Union => crate::functions::max(a.iter().map(sup).collect::<Option<Vec<_>>>()?),
        SetOp::Reals | SetOp::Integers | SetOp::Rationals | SetOp::Naturals | SetOp::Naturals0 => Expr::oo(),
        SetOp::Range if a[0] == a[1] || a.iter().any(|v| !crate::traverse::free_symbols(v).is_empty()) => return None,
        SetOp::Range => {
            if a[2].as_number().is_some_and(|n| n.is_positive()) {
                if a[1].is_constant(Constant::Infinity) {
                    a[1].clone()
                } else {
                    crate::ops::sub(&a[1], &a[2])
                }
            } else {
                a[0].clone()
            }
        }
        _ => return None,
    })
}

/// Lebesgue measure (`None` when not known).
pub fn measure(s: &Expr) -> Option<Expr> {
    let (op, a) = as_set(s)?;
    Some(match op {
        SetOp::Interval => crate::ops::sub(&a[1], &a[0]),
        SetOp::FiniteSet | SetOp::EmptySet => Expr::zero(),
        SetOp::Reals => Expr::oo(),
        SetOp::Union => {
            // the evaluated union of numeric intervals is disjoint
            let ms: Option<Vec<Expr>> = a.iter().map(measure).collect();
            crate::ops::add(&ms?)
        }
        _ => return None,
    })
}

/// Topological boundary.
pub fn boundary(s: &Expr) -> Option<Expr> {
    let (op, a) = as_set(s)?;
    Some(match op {
        SetOp::Interval => finite_set(a[..2].iter().filter(|p| !is_infinite(p)).cloned().collect()),
        SetOp::FiniteSet => s.clone(),
        SetOp::EmptySet | SetOp::Reals | SetOp::UniversalSet | SetOp::Complexes => empty_set(),
        SetOp::Integers | SetOp::Naturals | SetOp::Naturals0 => s.clone(),
        SetOp::Rationals => reals(),
        SetOp::Union => {
            let mut parts = Vec::new();
            for (i, p) in a.iter().enumerate() {
                let mut b = boundary(p)?;
                for (j, q) in a.iter().enumerate() {
                    if i != j {
                        b = complement(&b, &interior(q)?);
                    }
                }
                parts.push(b);
            }
            union(parts)
        }
        _ => return None,
    })
}

/// Interior.
pub fn interior(s: &Expr) -> Option<Expr> {
    let (op, a) = as_set(s)?;
    Some(match op {
        SetOp::Interval => interval(&a[0], &a[1], true, true),
        SetOp::FiniteSet | SetOp::EmptySet | SetOp::Integers | SetOp::Naturals | SetOp::Naturals0 | SetOp::Rationals => {
            empty_set()
        }
        SetOp::Reals | SetOp::Complexes | SetOp::UniversalSet => s.clone(),
        SetOp::Union => union(a.iter().map(interior).collect::<Option<Vec<_>>>()?),
        _ => return None,
    })
}

/// Closure.
pub fn closure(s: &Expr) -> Option<Expr> {
    let (op, a) = as_set(s)?;
    Some(match op {
        SetOp::Interval => interval(&a[0], &a[1], false, false),
        SetOp::FiniteSet | SetOp::EmptySet | SetOp::Integers | SetOp::Naturals | SetOp::Naturals0 => s.clone(),
        SetOp::Reals | SetOp::Complexes | SetOp::UniversalSet => s.clone(),
        SetOp::Rationals => reals(),
        SetOp::Union => union(a.iter().map(closure).collect::<Option<Vec<_>>>()?),
        _ => return None,
    })
}

fn std_rank(op: SetOp) -> Option<u8> {
    Some(match op {
        SetOp::Naturals => 0,
        SetOp::Naturals0 => 1,
        SetOp::Integers => 2,
        SetOp::Rationals => 3,
        SetOp::Reals => 4,
        SetOp::Complexes => 5,
        _ => return None,
    })
}

/// Is `a ⊆ b`? (`None` when undecided.)
pub fn is_subset(a: &Expr, b: &Expr) -> Option<bool> {
    if a == b || is_empty_set(a) || is_op(b, SetOp::UniversalSet) {
        return Some(true);
    }
    if is_empty_set(b) {
        return Some(false);
    }
    let (oa, aa) = as_set(a)?;
    let (ob, ba) = as_set(b)?;
    if let (Some(ra), Some(rb)) = (std_rank(oa), std_rank(ob)) {
        return Some(ra <= rb);
    }
    match oa {
        SetOp::FiniteSet => {
            let mut all = Some(true);
            for x in aa {
                match contains3(b, x) {
                    Some(true) => {}
                    Some(false) => return Some(false),
                    None => all = None,
                }
            }
            return all;
        }
        SetOp::Union => {
            let mut all = Some(true);
            for p in aa {
                match is_subset(p, b) {
                    Some(true) => {}
                    Some(false) => return Some(false),
                    None => all = None,
                }
            }
            return all;
        }
        SetOp::Range if ob == SetOp::Integers || ob == SetOp::Rationals || ob == SetOp::Reals || ob == SetOp::Complexes => {
            return Some(true);
        }
        _ => {}
    }
    if oa == SetOp::Interval {
        match ob {
            SetOp::Reals | SetOp::Complexes => return Some(true),
            SetOp::FiniteSet | SetOp::Integers | SetOp::Naturals | SetOp::Naturals0 | SetOp::Rationals | SetOp::Range => {
                return Some(false)
            }
            _ => {}
        }
    }
    if let (Some((s1, e1, l1, r1)), Some((s2, e2, l2, r2))) = (interval_parts(a), interval_parts(b)) {
        let left = match cmp3(&s2, &s1)? {
            Ordering::Less => true,
            Ordering::Equal => !l2 || l1,
            Ordering::Greater => false,
        };
        let right = match cmp3(&e1, &e2)? {
            Ordering::Less => true,
            Ordering::Equal => !r2 || r1,
            Ordering::Greater => false,
        };
        return Some(left && right);
    }
    if ob == SetOp::Union && ba.iter().any(|p| is_subset(a, p) == Some(true)) {
        return Some(true);
    }
    if oa == SetOp::Reals && matches!(ob, SetOp::Interval | SetOp::FiniteSet) {
        return Some(false);
    }
    // the reference fallback `self.intersect(other) == self` (True or None)
    if oa == SetOp::ProductSet && ob == SetOp::ProductSet {
        return (super::intersection(vec![a.clone(), b.clone()]) == *a).then_some(true);
    }
    if oa == SetOp::ImageSet && ob == SetOp::Integers {
        return (super::lattice_intersection(a, b).as_ref() == Some(a)).then_some(true);
    }
    None
}
