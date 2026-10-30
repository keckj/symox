//! Sets: `Interval`, `FiniteSet`, `Union`, `Intersection`, `Complement`,
//! `ProductSet`, `ImageSet`, `ConditionSet`, `Range` and the named sets
//! (`EmptySet`, `UniversalSet`, `Reals`, `Integers`, `Naturals`,
//! `Naturals0`, `Complexes`, `Rationals`), with SymPy's automatic
//! evaluation and printing; membership ([`contains`]), `inf`/`sup`,
//! `measure`, `boundary`, `as_relational`, and [`solveset`] /
//! inequality solving.
//!
//! Sets are `Kind::Func(Func::Set(op), args)` nodes; `Interval(a, b, l, r)`
//! stores its openness flags as boolean atoms like SymPy.

mod props;
mod solveset;
mod solveset_steps;
mod solveset_modular;
pub use solveset_modular::{solveset_modular, ModularError};
mod solveset_radical;
pub use solveset_radical::{negative_base_guard, negative_base_step, radical_step, StepHost};
mod lambda;
pub use lambda::{identity_function,lambda_checked,lambda_signature,lambda_arity,lambda_parts,lambda_apply,validate_lambda,validate_lambda_tree};

pub use props::{as_relational, boundary, closure, contains, inf, interior, is_subset, measure, sup};
pub use solveset::{reduce_inequalities, register_periodic_host, solve_univariate_inequality, solveset, solveset_real, with_periodic_window, PeriodicHost, PeriodicWindow};

use std::cmp::Ordering;

use rustc_hash::FxHashSet;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::logic::ordered;

/// Set heads (and the helpers `Lambda`, `Contains`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SetOp {
    Interval,
    FiniteSet,
    Union,
    Intersection,
    Complement,
    ProductSet,
    ImageSet,
    ConditionSet,
    Range,
    Lambda,
    IdentityFunction,
    Contains,
    EmptySet,
    UniversalSet,
    Reals,
    Integers,
    Naturals,
    Naturals0,
    Complexes,
    Rationals,
}

impl SetOp {
    pub fn name(self) -> &'static str {
        match self {
            SetOp::Interval => "Interval",
            SetOp::FiniteSet => "FiniteSet",
            SetOp::Union => "Union",
            SetOp::Intersection => "Intersection",
            SetOp::Complement => "Complement",
            SetOp::ProductSet => "ProductSet",
            SetOp::ImageSet => "ImageSet",
            SetOp::ConditionSet => "ConditionSet",
            SetOp::Range => "Range",
            SetOp::Lambda => "Lambda",
            SetOp::IdentityFunction => "IdentityFunction",
            SetOp::Contains => "Contains",
            SetOp::EmptySet => "EmptySet",
            SetOp::UniversalSet => "UniversalSet",
            SetOp::Reals => "Reals",
            SetOp::Integers => "Integers",
            SetOp::Naturals => "Naturals",
            SetOp::Naturals0 => "Naturals0",
            SetOp::Complexes => "Complexes",
            SetOp::Rationals => "Rationals",
        }
    }

    pub fn from_name(name: &str) -> Option<SetOp> {
        Some(match name {
            "Interval" => SetOp::Interval,
            "FiniteSet" => SetOp::FiniteSet,
            "Union" => SetOp::Union,
            "Intersection" => SetOp::Intersection,
            "Complement" => SetOp::Complement,
            "ProductSet" => SetOp::ProductSet,
            "ImageSet" => SetOp::ImageSet,
            "ConditionSet" => SetOp::ConditionSet,
            "Range" => SetOp::Range,
            "Lambda" => SetOp::Lambda,
            "IdentityFunction" => SetOp::IdentityFunction,
            "Contains" => SetOp::Contains,
            "EmptySet" => SetOp::EmptySet,
            "UniversalSet" => SetOp::UniversalSet,
            "Reals" => SetOp::Reals,
            "Integers" => SetOp::Integers,
            "Naturals" => SetOp::Naturals,
            "Naturals0" => SetOp::Naturals0,
            "Complexes" => SetOp::Complexes,
            "Rationals" => SetOp::Rationals,
            _ => return None,
        })
    }

    /// Named (argument free) sets.
    pub fn is_named(self) -> bool {
        matches!(
            self,
            SetOp::EmptySet
                | SetOp::UniversalSet
                | SetOp::Reals
                | SetOp::Integers
                | SetOp::Naturals
                | SetOp::Naturals0
                | SetOp::Complexes
                | SetOp::Rationals
        )
    }
}

// ----------------------------------------------------------------------
// basic helpers
// ----------------------------------------------------------------------

#[inline]
pub(crate) fn mk(op: SetOp, args: Vec<Expr>) -> Expr {
    Expr::from_kind(Kind::Func(Func::Set(op), args))
}

/// Head and arguments of a set (or `Lambda`/`Contains`) node.
#[inline]
pub fn as_set(e: &Expr) -> Option<(SetOp, &[Expr])> {
    match e.kind() {
        Kind::Func(Func::Set(op), a) => Some((*op, a.as_slice())),
        _ => None,
    }
}

/// Is `e` a set?
pub fn is_set(e: &Expr) -> bool {
    matches!(as_set(e), Some((op, _)) if !matches!(op, SetOp::Lambda | SetOp::IdentityFunction | SetOp::Contains))
}

#[inline]
pub(crate) fn is_op(e: &Expr, op: SetOp) -> bool {
    matches!(e.kind(), Kind::Func(Func::Set(o), _) if *o == op)
}

pub fn empty_set() -> Expr {
    mk(SetOp::EmptySet, vec![])
}
pub fn universal_set() -> Expr {
    mk(SetOp::UniversalSet, vec![])
}
pub fn reals() -> Expr {
    mk(SetOp::Reals, vec![])
}
pub fn integers() -> Expr {
    mk(SetOp::Integers, vec![])
}
pub fn naturals() -> Expr {
    mk(SetOp::Naturals, vec![])
}
pub fn naturals0() -> Expr {
    mk(SetOp::Naturals0, vec![])
}
pub fn complexes() -> Expr {
    mk(SetOp::Complexes, vec![])
}
pub fn rationals() -> Expr {
    mk(SetOp::Rationals, vec![])
}

/// A named set by SymPy name (`"Reals"`, `"EmptySet"`, ...).
pub fn named_set(name: &str) -> Option<Expr> {
    let op = SetOp::from_name(name)?;
    op.is_named().then(|| mk(op, vec![]))
}

pub(crate) fn is_empty_set(e: &Expr) -> bool {
    is_op(e, SetOp::EmptySet)
}

/// `(start, end, left_open, right_open)` of an interval; `Reals` counts
/// as `(-oo, oo)`.
pub fn interval_parts(e: &Expr) -> Option<(Expr, Expr, bool, bool)> {
    match as_set(e)? {
        (SetOp::Interval, a) => Some((a[0].clone(), a[1].clone(), a[2] == Expr::boolean(true), a[3] == Expr::boolean(true))),
        (SetOp::Reals, _) => Some((Expr::neg_oo(), Expr::oo(), true, true)),
        _ => None,
    }
}

/// Elements of a finite set.
pub fn finite_elements(e: &Expr) -> Option<&[Expr]> {
    match as_set(e)? {
        (SetOp::FiniteSet, a) => Some(a),
        _ => None,
    }
}

// ----------------------------------------------------------------------
// comparisons of extended reals
// ----------------------------------------------------------------------

/// Three-valued comparison of (extended) real numbers or expressions.
pub(crate) fn cmp3(a: &Expr, b: &Expr) -> Option<Ordering> {
    if a == b {
        return Some(Ordering::Equal);
    }
    let inf = |e: &Expr| match e.kind() {
        Kind::Const(Constant::Infinity) => Some(1i8),
        Kind::Const(Constant::NegInfinity) => Some(-1),
        _ => None,
    };
    match (inf(a), inf(b)) {
        (Some(x), Some(y)) => return Some(x.cmp(&y)),
        (Some(x), None) => {
            return finite_real(b).then_some(if x > 0 { Ordering::Greater } else { Ordering::Less });
        }
        (None, Some(y)) => {
            return finite_real(a).then_some(if y > 0 { Ordering::Less } else { Ordering::Greater });
        }
        _ => {}
    }
    if let (Kind::Num(x), Kind::Num(y)) = (a.kind(), b.kind()) {
        return x.cmp_value(y);
    }
    let d = crate::ops::sub(a, b);
    if let Kind::Num(n) = d.kind() {
        return n.cmp_value(&crate::number::Number::ZERO);
    }
    if crate::functions::is_numeric_value(&d) {
        if let Some(z) = crate::evalf::eval_c64(&d) {
            if z.im.abs() > 1e-12 * (1.0 + z.re.abs()) {
                return None;
            }
            if z.re.abs() > 1e-12 {
                return z.re.partial_cmp(&0.0);
            }
            if crate::expand::expand(&d).is_zero() {
                return Some(Ordering::Equal);
            }
            return None;
        }
    }
    use crate::assumptions as asm;
    if asm::is_zero(&d) == Some(true) {
        return Some(Ordering::Equal);
    }
    if asm::is_positive(&d) == Some(true) {
        return Some(Ordering::Greater);
    }
    if asm::is_negative(&d) == Some(true) {
        return Some(Ordering::Less);
    }
    None
}

fn finite_real(e: &Expr) -> bool {
    crate::assumptions::is_real(e) == Some(true)
}

pub(crate) fn lt3(a: &Expr, b: &Expr) -> Option<bool> {
    cmp3(a, b).map(|o| o == Ordering::Less)
}

pub(crate) fn is_infinite(e: &Expr) -> bool {
    e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity)
}

// ----------------------------------------------------------------------
// constructors with evaluation
// ----------------------------------------------------------------------

/// `Interval(a, b, left_open, right_open)`.
pub fn interval(a: &Expr, b: &Expr, lopen: bool, ropen: bool) -> Expr {
    if a.is_constant(Constant::Infinity) || b.is_constant(Constant::NegInfinity) {
        return empty_set();
    }
    let lopen = lopen || a.is_constant(Constant::NegInfinity);
    let ropen = ropen || b.is_constant(Constant::Infinity);
    match cmp3(b, a) {
        Some(Ordering::Less) => return empty_set(),
        Some(Ordering::Equal) => {
            return if lopen || ropen { empty_set() } else { finite_set(vec![a.clone()]) };
        }
        _ => {}
    }
    mk(SetOp::Interval, vec![a.clone(), b.clone(), Expr::boolean(lopen), Expr::boolean(ropen)])
}

/// `FiniteSet(*elements)`.
pub fn finite_set(elems: Vec<Expr>) -> Expr {
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let mut v: Vec<Expr> = Vec::with_capacity(elems.len());
    for e in elems {
        if seen.insert(e.clone()) {
            v.push(e);
        }
    }
    if v.is_empty() {
        return empty_set();
    }
    ordered(&mut v);
    mk(SetOp::FiniteSet, v)
}

fn finite_contains_structural(fs: &[Expr], x: &Expr) -> bool {
    fs.iter().any(|e| e == x)
}

/// Membership as a three valued boolean.
pub(crate) fn contains3(s: &Expr, x: &Expr) -> Option<bool> {
    match contains(s, x).kind() {
        Kind::Bool(b) => Some(*b),
        _ => None,
    }
}

/// The reference `ordered(args, Set._infimum_key)` of the arguments of an
/// unevaluated `Union` or `Intersection`: by the numeric infimum (`oo` when
/// it is unknown or not comparable), ties in the default `ordered` order.
pub(crate) fn ordered_by_infimum(v: &mut [Expr]) {
    ordered(v);
    let key = |s: &Expr| -> f64 {
        match inf(s) {
            Some(i) if i.is_constant(Constant::NegInfinity) => f64::NEG_INFINITY,
            Some(i) if crate::integrate::piecewise::is_comparable(&i) => {
                crate::evalf::eval_f64(&i).unwrap_or(f64::INFINITY)
            }
            _ => f64::INFINITY,
        }
    };
    let mut keyed: Vec<(f64, Expr)> = v.iter().map(|e| (key(e), e.clone())).collect();
    keyed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    for (slot, (_, e)) in v.iter_mut().zip(keyed) {
        *slot = e;
    }
}

/// `Union(*sets)`.
pub fn union(sets: Vec<Expr>) -> Expr {
    let mut args: Vec<Expr> = Vec::new();
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let mut stack: Vec<Expr> = sets;
    stack.reverse();
    let mut finite: Vec<Expr> = Vec::new();
    while let Some(s) = stack.pop() {
        match as_set(&s) {
            Some((SetOp::Union, sub)) => {
                for t in sub.iter().rev() {
                    stack.push(t.clone());
                }
            }
            Some((SetOp::EmptySet, _)) => {}
            Some((SetOp::UniversalSet, _)) => return universal_set(),
            Some((SetOp::FiniteSet, el)) => finite.extend(el.iter().cloned()),
            _ => {
                if seen.insert(s.clone()) {
                    args.push(s);
                }
            }
        }
    }
    if !finite.is_empty() {
        args.push(finite_set(finite));
    }
    // pairwise simplification until nothing changes
    'outer: loop {
        for i in 0..args.len() {
            for j in 0..args.len() {
                if i == j {
                    continue;
                }
                if let Some(rep) = union_pair(&args[i], &args[j]) {
                    let (hi, lo) = if i > j { (i, j) } else { (j, i) };
                    args.remove(hi);
                    args.remove(lo);
                    for r in rep {
                        match as_set(&r) {
                            Some((SetOp::Union, sub)) => args.extend(sub.iter().cloned()),
                            Some((SetOp::EmptySet, _)) => {}
                            _ => args.push(r),
                        }
                    }
                    continue 'outer;
                }
            }
        }
        break;
    }
    match args.len() {
        0 => empty_set(),
        1 => args.pop().unwrap(),
        _ => {
            ordered_by_infimum(&mut args);
            mk(SetOp::Union, args)
        }
    }
}

/// Simplify `a ∪ b` if a rule applies: the reference `union_sets` rules,
/// dispatched on the ordered pair (the caller tries both orders).
fn union_pair(a: &Expr, b: &Expr) -> Option<Vec<Expr>> {
    let (oa, _) = as_set(a)?;
    let ob = as_set(b).map(|(op, _)| op)?;
    let keep_a = || Some(vec![a.clone()]);
    match (oa, ob) {
        (SetOp::Naturals0, SetOp::Naturals)
        | (SetOp::Rationals, SetOp::Naturals | SetOp::Naturals0)
        | (SetOp::Reals, SetOp::Naturals | SetOp::Naturals0 | SetOp::Rationals) => return keep_a(),
        _ => {}
    }
    match oa {
        SetOp::Integers => {
            let common = intersection(vec![a.clone(), b.clone()]);
            if common == *a {
                return Some(vec![b.clone()]);
            }
            if common == *b {
                return keep_a();
            }
            return None;
        }
        SetOp::Complexes => {
            // `b` real: `ComplexRegion.from_real(b)` is absorbed by the plane.
            return (is_subset(b, &reals()) == Some(true)).then(|| vec![a.clone()]);
        }
        SetOp::ProductSet => {
            if is_subset(b, a) == Some(true) {
                return keep_a();
            }
            if let (Some((_, pa)), Some((SetOp::ProductSet, pb))) = (as_set(a), as_set(b)) {
                if pa.len() == 2 && pb.len() == 2 {
                    if pa[0] == pb[0] {
                        return Some(vec![product_set(vec![pa[0].clone(), union(vec![pa[1].clone(), pb[1].clone()])])]);
                    }
                    if pa[1] == pb[1] {
                        return Some(vec![product_set(vec![union(vec![pa[0].clone(), pb[0].clone()]), pa[1].clone()])]);
                    }
                }
            }
            return None;
        }
        _ => {}
    }
    if let Some((s1, e1, l1, r1)) = interval_parts(a) {
        if let Some((s2, e2, l2, r2)) = interval_parts(b) {
            // merged only when every endpoint is comparable
            let comparable = [&s1, &e1, &s2, &e2].iter().all(|e| crate::integrate::piecewise::is_comparable(e));
            if !comparable {
                return None;
            }
            let c = cmp3(&s1, &s2)?;
            let ((s1, e1, l1, r1), (s2, e2, l2, r2)) = if c == Ordering::Greater {
                ((s2, e2, l2, r2), (s1, e1, l1, r1))
            } else {
                ((s1, e1, l1, r1), (s2, e2, l2, r2))
            };
            // disjoint with a gap?
            match cmp3(&e1, &s2)? {
                Ordering::Less => return None,
                Ordering::Equal if r1 && l2 => return None,
                _ => {}
            }
            let lopen = if c == Ordering::Equal { l1 && l2 } else { l1 };
            let (end, ropen) = match cmp3(&e1, &e2)? {
                Ordering::Less => (e2, r2),
                Ordering::Greater => (e1, r1),
                Ordering::Equal => (e1, r1 && r2),
            };
            return Some(vec![interval(&s1, &end, lopen, ropen)]);
        }
        // open finite endpoints contained in `b` are closed
        let finite = |e: &Expr| crate::assumptions::fact(e, "finite") == Some(true);
        let close_l = l1 && contains3(b, &s1) == Some(true) && finite(&s1);
        let close_r = r1 && contains3(b, &e1) == Some(true) && finite(&e1);
        if close_l || close_r {
            let lopen = l1 && contains3(b, &s1) != Some(true);
            let ropen = r1 && contains3(b, &e1) != Some(true);
            return Some(vec![interval(&s1, &e1, lopen, ropen), b.clone()]);
        }
        return None;
    }
    if let Some(el) = finite_elements(a) {
        // elements of `a` already in `b` are dropped
        if el.iter().any(|x| contains3(b, x) == Some(true)) {
            let rest: Vec<Expr> = el.iter().filter(|x| contains3(b, x) != Some(true)).cloned().collect();
            return Some(vec![finite_set(rest), b.clone()]);
        }
    }
    None
}

/// `Intersection(*sets)`.
pub fn intersection(sets: Vec<Expr>) -> Expr {
    let mut args: Vec<Expr> = Vec::new();
    let mut seen: FxHashSet<Expr> = FxHashSet::default();
    let mut stack: Vec<Expr> = sets;
    stack.reverse();
    while let Some(s) = stack.pop() {
        match as_set(&s) {
            Some((SetOp::Intersection, sub)) => {
                for t in sub.iter().rev() {
                    stack.push(t.clone());
                }
            }
            Some((SetOp::EmptySet, _)) => return empty_set(),
            Some((SetOp::UniversalSet, _)) => {}
            _ => {
                if seen.insert(s.clone()) {
                    args.push(s);
                }
            }
        }
    }
    match args.len() {
        0 => return universal_set(),
        1 => return args.pop().unwrap(),
        _ => {}
    }
    // finite sets
    if args.iter().any(|a| is_op(a, SetOp::FiniteSet)) {
        return intersect_finite(args);
    }
    // distribute over unions
    if let Some(pos) = args.iter().position(|a| is_op(a, SetOp::Union)) {
        let u = args.remove(pos);
        let other = intersection(args);
        let (_, parts) = as_set(&u).unwrap();
        return union(parts.iter().map(|p| intersection(vec![p.clone(), other.clone()])).collect());
    }
    // unevaluated complements
    if let Some(pos) = args.iter().position(|a| is_op(a, SetOp::Complement)) {
        let c = args.remove(pos);
        let (_, ca) = as_set(&c).unwrap();
        args.push(ca[0].clone());
        return complement(&intersection(args), &ca[1]);
    }
    'outer: loop {
        for i in 0..args.len() {
            for j in 0..args.len() {
                if i == j {
                    continue;
                }
                if let Some(r) = intersect_pair(&args[i], &args[j]) {
                    let (hi, lo) = if i > j { (i, j) } else { (j, i) };
                    args.remove(hi);
                    args.remove(lo);
                    if is_empty_set(&r) {
                        return r;
                    }
                    args.push(r);
                    if args.len() == 1 {
                        break 'outer;
                    }
                    continue 'outer;
                }
            }
        }
        break;
    }
    if args.len() == 1 {
        return args.pop().unwrap();
    }
    ordered_by_infimum(&mut args);
    mk(SetOp::Intersection, args)
}

fn intersect_finite(args: Vec<Expr>) -> Expr {
    let (fs, others): (Vec<Expr>, Vec<Expr>) = args.into_iter().partition(|a| is_op(a, SetOp::FiniteSet));
    let all: Vec<Expr> = fs.iter().chain(others.iter()).cloned().collect();
    // prune elements definitely outside some other set
    let mut pruned: Vec<Vec<Expr>> = Vec::new();
    let mut fully_definite: Option<usize> = None;
    for (k, f) in fs.iter().enumerate() {
        let el = finite_elements(f).unwrap();
        let mut keep = Vec::new();
        let mut all_def = true;
        for x in el {
            let mut status = Some(true);
            for (m, s) in all.iter().enumerate() {
                if m == k {
                    continue;
                }
                let c = match finite_elements(s) {
                    Some(sel) if finite_contains_structural(sel, x) => Some(true),
                    _ => contains3(s, x),
                };
                match c {
                    Some(false) => {
                        status = Some(false);
                        break;
                    }
                    None => status = None,
                    _ => {}
                }
            }
            match status {
                Some(false) => {}
                Some(true) => keep.push(x.clone()),
                None => {
                    all_def = false;
                    keep.push(x.clone());
                }
            }
        }
        if keep.is_empty() {
            return empty_set();
        }
        if all_def && fully_definite.is_none() {
            fully_definite = Some(k);
        }
        pruned.push(keep);
    }
    if let Some(k) = fully_definite {
        return finite_set(pruned.swap_remove(k));
    }
    let mut v: Vec<Expr> = pruned.into_iter().map(finite_set).collect();
    v.extend(others);
    ordered_by_infimum(&mut v);
    mk(SetOp::Intersection, v)
}

fn subset_chain_rank(op: SetOp) -> Option<u8> {
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

fn ceil_int(e: &Expr) -> Option<Expr> {
    let f = crate::functions::ceiling(e);
    f.is_integer_number().then_some(f)
}

fn floor_int(e: &Expr) -> Option<Expr> {
    let f = crate::functions::floor(e);
    f.is_integer_number().then_some(f)
}

/// Integers in an interval as a `Range` (or `EmptySet`).
fn integers_in(s: &Expr, e: &Expr, l: bool, r: bool, min: Option<i64>) -> Option<Expr> {
    use crate::ops::add2;
    let one = Expr::one();
    let mut lo = if s.is_constant(Constant::NegInfinity) {
        s.clone()
    } else {
        let c = ceil_int(s)?;
        if l && c == *s {
            add2(&c, &one)
        } else {
            c
        }
    };
    if let Some(m) = min {
        let m = Expr::int(m);
        if lo.is_constant(Constant::NegInfinity) || lt3(&lo, &m) == Some(true) {
            lo = m;
        }
    }
    let hi = if e.is_constant(Constant::Infinity) {
        e.clone()
    } else {
        let f = floor_int(e)?;
        let f = if r && f == *e { crate::ops::sub(&f, &one) } else { f };
        add2(&f, &one)
    };
    if !is_infinite(&lo) && !is_infinite(&hi) && lt3(&lo, &hi) != Some(true) {
        return Some(empty_set());
    }
    Some(mk(SetOp::Range, vec![lo, hi, one]))
}

/// Points of `ImageSet(Lambda(n, p*n + c), Integers)` in a bounded
/// numeric interval.
fn image_in_interval(img: &Expr, iv: &Expr) -> Option<Expr> {
    let (SetOp::ImageSet, a) = as_set(img)? else { return None };
    if !is_op(&a[1], SetOp::Integers) {
        return None;
    }
    if lambda_arity(&a[0])!=Some(1){return None;}
    let (variables, body)=lambda_parts(&a[0])?;
    if variables.len()!=1{return None;}
    let n=&variables[0];
    let (s, e, _, _) = interval_parts(iv)?;
    if is_infinite(&s) || is_infinite(&e) {
        return None;
    }
    let body = &body;
    let p = crate::diff::diff(body, n);
    let c = crate::expand::expand(&crate::subs::subs1(body, n, &Expr::zero()));
    if crate::traverse::has(&p, n) || p.is_zero() {
        return None;
    }
    let (pv, cv) = (crate::evalf::eval_f64(&p)?, crate::evalf::eval_f64(&c)?);
    let (sv, ev) = (crate::evalf::eval_f64(&s)?, crate::evalf::eval_f64(&e)?);
    let (k1, k2) = {
        let a = (sv - cv) / pv;
        let b = (ev - cv) / pv;
        (a.min(b).floor() as i64 - 1, a.max(b).ceil() as i64 + 1)
    };
    if k2 - k1 > 10_000 {
        return None;
    }
    let mut pts = Vec::new();
    for k in k1..=k2 {
        let v = crate::expand::expand(&crate::subs::subs1(body, n, &Expr::int(k)));
        match contains3(iv, &v) {
            Some(true) => pts.push(v),
            Some(false) => {}
            None => return None,
        }
    }
    Some(finite_set(pts))
}


/// `ImageSet(Lambda(n, p*n + q), Integers)` with exact rational `p != 0`
/// and `q`: the variable and the coefficients.
fn integer_lattice(img: &Expr) -> Option<(Expr, num_rational::BigRational, num_rational::BigRational)> {
    let (SetOp::ImageSet, ia) = as_set(img)? else { return None };
    if !is_op(&ia[1], SetOp::Integers) || lambda_arity(&ia[0]) != Some(1) {
        return None;
    }
    let (variables, body) = lambda_parts(&ia[0])?;
    if variables.len() != 1 {
        return None;
    }
    let n = variables[0].clone();
    let p = crate::diff::diff(&body, &n);
    if crate::traverse::has(&p, &n) {
        return None;
    }
    let q = crate::expand::expand(&crate::subs::subs1(&body, &n, &Expr::zero()));
    let (p, q) = (p.as_number()?.to_ratio()?, q.as_number()?.to_ratio()?);
    if num_traits::Zero::is_zero(&p) {
        return None;
    }
    Some((n, p, q))
}

/// Linear case of the reference `intersection_sets(ImageSet, Set)` with
/// `Integers` or another integer image set: the solutions of the linear
/// diophantine equation `p*n + q = r*m + s`, as the canonical image set the
/// reference `imageset` builds (`L*n + (b mod L)`, `Integers` for `n`).
pub(crate) fn lattice_intersection(a: &Expr, b: &Expr) -> Option<Expr> {
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_rational::BigRational;
    use num_traits::{One, Signed, Zero};
    let (n, p, q) = integer_lattice(a)?;
    let (r, s) = if is_op(b, SetOp::Integers) {
        (BigRational::one(), BigRational::zero())
    } else {
        let (_, r, s) = integer_lattice(b)?;
        (r, s)
    };
    // scale by the common denominator: P*n - R*m = S - Q over the integers
    let d = [&p, &q, &r, &s].iter().fold(BigInt::one(), |acc, v| acc.lcm(v.denom()));
    let scale = |v: &BigRational| (v * BigRational::from_integer(d.clone())).to_integer();
    let (pp, qq, rr, ss) = (scale(&p), scale(&q), scale(&r), scale(&s));
    let eg = pp.extended_gcd(&rr);
    let g = eg.gcd.abs();
    let rhs = &ss - &qq;
    if !(&rhs % &g).is_zero() {
        return Some(empty_set());
    }
    // pp*x + rr*y = gcd(pp, rr) (up to the sign of the gcd)
    let n0 = &eg.x * (&rhs / &eg.gcd);
    let den = BigRational::from_integer(d);
    let offset = (BigRational::from_integer(&pp * &n0 + &qq)) / &den;
    let step = BigRational::from_integer((&pp * &rr / &g).abs()) / &den;
    let modulo = &offset - (&offset / &step).floor() * &step;
    let num = |v: BigRational| Expr::num(crate::number::Number::from_ratio(v));
    let body = if step.is_one() {
        if modulo.is_zero() || modulo.is_integer() {
            return Some(integers());
        }
        crate::ops::add2(&n, &num(modulo))
    } else {
        crate::ops::add2(&crate::ops::mul2(&num(step), &n), &num(modulo))
    };
    Some(image_set(&lambda(&n, &body), &integers()))
}

/// An integer arithmetic progression `{c + k*step}` bounded by `lo` and
/// `hi` (`None` for an infinite end).
#[derive(Clone)]
struct Progression {
    lo: Option<num_bigint::BigInt>,
    hi: Option<num_bigint::BigInt>,
    step: num_bigint::BigInt,
}

impl Progression {
    fn point(&self) -> num_bigint::BigInt {
        self.lo.clone().or_else(|| self.hi.clone()).unwrap_or_default()
    }
}

/// The elements of a `Range` whose arguments are all numbers: `None` for
/// symbolic arguments, `Some(None)` for an empty range.
fn range_progression(r: &Expr) -> Option<Option<Progression>> {
    use num_integer::Integer;
    use num_traits::{One, Signed, Zero};
    let (SetOp::Range, a) = as_set(r)? else { return None };
    if a.len() != 3 {
        return None;
    }
    let val = |e: &Expr| -> Option<Result<num_bigint::BigInt, i8>> {
        if e.is_constant(Constant::Infinity) {
            return Some(Err(1));
        }
        if e.is_constant(Constant::NegInfinity) {
            return Some(Err(-1));
        }
        e.as_number()?.to_bigint().map(Ok)
    };
    let (start, stop, step) = (val(&a[0])?, val(&a[1])?, val(&a[2])?.ok()?);
    if step.is_zero() {
        return None;
    }
    let pos = step.is_positive();
    let st = step.abs();
    let p = |lo, hi| Some(Some(Progression { lo, hi, step: st.clone() }));
    match (start, stop) {
        (Ok(s), Ok(e)) => {
            let n = (&e - &s).div_ceil(&step);
            if !n.is_positive() {
                return Some(None);
            }
            let last = &s + (n - 1) * &step;
            if pos { p(Some(s), Some(last)) } else { p(Some(last), Some(s)) }
        }
        (Ok(s), Err(d)) => match ((d > 0) == pos, pos) {
            (false, _) => Some(None),
            (true, true) => p(Some(s), None),
            (true, false) => p(None, Some(s)),
        },
        (Err(d), Ok(e)) => {
            // the elements are `stop - step`, `stop - 2*step`, ...
            let first = &e - &step;
            match ((d < 0) == pos, pos) {
                (false, _) => Some(None),
                (true, true) => p(None, Some(first)),
                (true, false) => p(Some(first), None),
            }
        }
        (Err(x), Err(y)) => (x != y && st.is_one()).then(|| Some(Progression { lo: None, hi: None, step: st })),
    }
}

/// The reference `intersection_sets(Range, Range)` on progressions: the
/// common elements as a `Range` with a positive step.
fn progression_intersection(p1: &Progression, p2: &Progression) -> Expr {
    use num_bigint::BigInt;
    use num_integer::Integer;
    let (c1, c2) = (p1.point(), p2.point());
    let eg = p1.step.extended_gcd(&p2.step);
    let g = eg.gcd.clone();
    let d = &c2 - &c1;
    if !d.is_multiple_of(&g) {
        return empty_set();
    }
    let l = &p1.step / &g * &p2.step;
    let c = (&c1 + &p1.step * (&eg.x * (&d / &g))).mod_floor(&l);
    let lo = match (&p1.lo, &p2.lo) {
        (Some(x), Some(y)) => Some(x.max(y).clone()),
        (Some(x), None) | (None, Some(x)) => Some(x.clone()),
        _ => None,
    };
    let hi = match (&p1.hi, &p2.hi) {
        (Some(x), Some(y)) => Some(x.min(y).clone()),
        (Some(x), None) | (None, Some(x)) => Some(x.clone()),
        _ => None,
    };
    let first = lo.map(|lo| &lo + (&c - &lo).mod_floor(&l));
    let last = hi.map(|hi| &hi - (&hi - &c).mod_floor(&l));
    let int = |v: BigInt| Expr::num(crate::number::Number::from_bigint(v));
    let step = int(l.clone());
    match (first, last) {
        (Some(f), Some(t)) => {
            if f > t {
                return empty_set();
            }
            mk(SetOp::Range, vec![int(f), int(t + &l), step])
        }
        (Some(f), None) => mk(SetOp::Range, vec![int(f), Expr::oo(), step]),
        (None, Some(t)) => mk(SetOp::Range, vec![Expr::neg_oo(), int(t + &l), step]),
        (None, None) => mk(SetOp::Range, vec![Expr::neg_oo(), Expr::oo(), step]),
    }
}

/// The reference `intersection_sets(Range, Range)`.
fn range_range(a: &Expr, b: &Expr) -> Option<Expr> {
    let (pa, pb) = (range_progression(a)?, range_progression(b)?);
    let (Some(pa), Some(pb)) = (pa, pb) else { return Some(empty_set()) };
    if pa.lo.is_none() && pa.hi.is_none() {
        return Some(b.clone());
    }
    if pb.lo.is_none() && pb.hi.is_none() {
        return Some(a.clone());
    }
    Some(progression_intersection(&pa, &pb))
}

/// The reference `intersection_sets(Range, Interval)`: the interval is
/// trimmed to the integers it contains and intersected as a `Range`.
fn range_interval(r: &Expr, s: &Expr, e: &Expr, l: bool, ro: bool) -> Option<Expr> {
    use num_traits::One;
    if !(crate::logic::is_number_like(s) || is_infinite(s)) || !(crate::logic::is_number_like(e) || is_infinite(e)) {
        return None;
    }
    let pr = range_progression(r)?;
    let Some(pr) = pr else { return Some(empty_set()) };
    let to_int = |v: &Expr| v.as_number().and_then(|n| n.to_bigint());
    let lo = if s.is_constant(Constant::NegInfinity) {
        None
    } else if s.is_constant(Constant::Infinity) {
        return Some(empty_set());
    } else {
        let c = ceil_int(s)?;
        let ci = to_int(&c)?;
        Some(if l && cmp3(&c, s)? == Ordering::Equal { ci + 1 } else { ci })
    };
    let hi = if e.is_constant(Constant::Infinity) {
        None
    } else if e.is_constant(Constant::NegInfinity) {
        return Some(empty_set());
    } else {
        let f = floor_int(e)?;
        let fi = to_int(&f)?;
        Some(if ro && cmp3(&f, e)? == Ordering::Equal { fi - 1 } else { fi })
    };
    if let (Some(x), Some(y)) = (&lo, &hi) {
        if x > y {
            return Some(empty_set());
        }
    }
    if pr.lo.is_none() && pr.hi.is_none() && lo.is_none() && hi.is_none() {
        return Some(r.clone());
    }
    Some(progression_intersection(&pr, &Progression { lo, hi, step: num_bigint::BigInt::one() }))
}

/// Simplify `a ∩ b` if a rule applies.
fn intersect_pair(a: &Expr, b: &Expr) -> Option<Expr> {
    if let (Some((s1, e1, l1, r1)), Some((s2, e2, l2, r2))) = (interval_parts(a), interval_parts(b)) {
        // (-oo, oo) with any interval gives that interval, even when its
        // bounds are not comparable (reference Interval/Interval rule).
        let whole = |s: &Expr, e: &Expr| *s == Expr::neg_oo() && *e == Expr::oo();
        if whole(&s1, &e1) {
            return Some(b.clone());
        }
        if whole(&s2, &e2) {
            return Some(a.clone());
        }
        let (start, lopen) = match cmp3(&s1, &s2)? {
            Ordering::Less => (s2, l2),
            Ordering::Greater => (s1, l1),
            Ordering::Equal => (s1, l1 || l2),
        };
        let (end, ropen) = match cmp3(&e1, &e2)? {
            Ordering::Less => (e1, r1),
            Ordering::Greater => (e2, r2),
            Ordering::Equal => (e1, r1 || r2),
        };
        return Some(interval(&start, &end, lopen, ropen));
    }
    let (oa, ob) = (as_set(a)?.0, as_set(b)?.0);
    // reference `intersection_sets(ProductSet, ProductSet)`
    if let (Some((SetOp::ProductSet, pa)), Some((SetOp::ProductSet, pb))) = (as_set(a), as_set(b)) {
        if pa.len() != pb.len() {
            return Some(empty_set());
        }
        return Some(product_set(pa.iter().zip(pb).map(|(i, j)| intersection(vec![i.clone(), j.clone()])).collect()));
    }
    // reference `Range` rules (`Range`, `Interval`, `Naturals`, `Integers`,
    // `Rationals`)
    if oa == SetOp::Range && ob == SetOp::Range {
        if let Some(r) = range_range(a, b) {
            return Some(r);
        }
    }
    for (x, y) in [(a, b), (b, a)] {
        if !is_op(x, SetOp::Range) {
            continue;
        }
        if let Some((s, e, l, r)) = interval_parts(y) {
            if let Some(v) = range_interval(x, &s, &e, l, r) {
                return Some(v);
            }
        }
        match as_set(y).map(|t| t.0) {
            Some(SetOp::Naturals) => return range_interval(x, &Expr::one(), &Expr::oo(), false, true),
            Some(SetOp::Naturals0) => return range_interval(x, &Expr::zero(), &Expr::oo(), false, true),
            Some(SetOp::Integers | SetOp::Rationals) => return Some(x.clone()),
            _ => {}
        }
    }
    for (x, y) in [(a, b), (b, a)] {
        if let Some(r) = lattice_intersection(x, y) {
            return Some(r);
        }
        // reference `intersection_sets(ConditionSet, Set)`
        if let (Some((SetOp::ConditionSet, ca)), false) = (as_set(x), is_op(y, SetOp::ConditionSet)) {
            if ca.len() == 3 {
                return Some(condition_set(&ca[0], &ca[1], &intersection(vec![ca[2].clone(), y.clone()])));
            }
        }
    }
    // nested standard sets
    if let (Some(ra), Some(rb)) = (subset_chain_rank(oa), subset_chain_rank(ob)) {
        return Some(if ra <= rb { a.clone() } else { b.clone() });
    }
    // integers / naturals with an interval
    for (x, y) in [(a, b), (b, a)] {
        if let (Some(op), Some((s, e, l, r))) = (as_set(x).map(|t| t.0), interval_parts(y)) {
            let min = match op {
                SetOp::Integers => Some(None),
                SetOp::Naturals => Some(Some(1)),
                SetOp::Naturals0 => Some(Some(0)),
                _ => None,
            };
            if let Some(min) = min {
                if let Some(r) = integers_in(&s, &e, l, r, min) {
                    return Some(r);
                }
            }
        }
        // everything here is a set of complex numbers
        if is_op(x, SetOp::Complexes) && interval_parts(y).is_some() {
            return Some(y.clone());
        }
        if let Some(r) = image_in_interval(x, y) {
            return Some(r);
        }
    }
    if is_subset(a, b) == Some(true) {
        return Some(a.clone());
    }
    if is_subset(b, a) == Some(true) {
        return Some(b.clone());
    }
    None
}

/// `Complement(a, b)` (`a - b`).
pub fn complement(a: &Expr, b: &Expr) -> Expr {
    if is_empty_set(a) {
        return empty_set();
    }
    if is_empty_set(b) {
        return a.clone();
    }
    if is_op(b, SetOp::UniversalSet) || is_subset(a, b) == Some(true) {
        return empty_set();
    }
    if let Some((SetOp::Union, parts)) = as_set(b) {
        return intersection(parts.iter().map(|p| complement(a, p)).collect());
    }
    // reference `Set._complement` of two products:
    // (A x B) - (C x D) = ((A - C) x B) U (A x (B - D))
    if let (Some((SetOp::ProductSet, pa)), Some((SetOp::ProductSet, pb))) = (as_set(a), as_set(b)) {
        if pa.len() != pb.len() {
            return a.clone();
        }
        let overlaps = (0..pa.len())
            .map(|k| product_set(pa.iter().zip(pb).enumerate().map(|(i, (o, s))| if i != k { o.clone() } else { complement(o, s) }).collect()))
            .collect();
        return union(overlaps);
    }
    // B-specific rules
    if let Some((s, e, l, r)) = interval_parts(b) {
        if is_op(a, SetOp::Reals) {
            return union(vec![
                interval(&Expr::neg_oo(), &s, true, !l),
                interval(&e, &Expr::oo(), !r, true),
            ]);
        }
    }
    if let Some(bel) = finite_elements(b) {
        if interval_parts(a).is_some() {
            let mut nums = Vec::new();
            let mut syms = Vec::new();
            for m in bel {
                let real = crate::assumptions::is_real(m);
                if crate::logic::is_number_like(m) && real == Some(true) {
                    nums.push(m.clone());
                } else if real == Some(false) {
                } else {
                    syms.push(m.clone());
                }
            }
            if is_op(a, SetOp::Reals) && !nums.is_empty() {
                nums.sort_by(|x, y| cmp3(x, y).unwrap_or(Ordering::Equal));
                let mut ivs = vec![interval(&Expr::neg_oo(), &nums[0], true, true)];
                for w in nums.windows(2) {
                    ivs.push(interval(&w[0], &w[1], true, true));
                }
                ivs.push(interval(nums.last().unwrap(), &Expr::oo(), true, true));
                let u = union(ivs);
                return if syms.is_empty() { u } else { mk(SetOp::Complement, vec![u, finite_set(syms)]) };
            } else if nums.is_empty() {
                return if syms.is_empty() { a.clone() } else { mk(SetOp::Complement, vec![a.clone(), finite_set(syms)]) };
            }
        }
        if let Some(ael) = finite_elements(a) {
            let unk: Vec<Expr> = bel
                .iter()
                .filter(|x| !finite_contains_structural(ael, x) && contains3(a, x).is_none())
                .cloned()
                .collect();
            if unk.len() != bel.len() {
                let not_true: Vec<Expr> = ael
                    .iter()
                    .filter(|x| !finite_contains_structural(bel, x) && contains3(b, x) != Some(true))
                    .cloned()
                    .collect();
                let nt = finite_set(not_true);
                return if unk.is_empty() { nt } else { complement_raw(nt, finite_set(unk)) };
            }
        }
    }
    // generic rules on A
    if interval_parts(a).is_some() && (interval_parts(b).is_some() || is_op(b, SetOp::FiniteSet)) {
        return intersection(vec![a.clone(), complement(&reals(), b)]);
    }
    if let Some((SetOp::Union, parts)) = as_set(a) {
        return union(parts.iter().map(|p| complement(p, b)).collect());
    }
    if let Some((SetOp::Complement, ca)) = as_set(a) {
        return mk(SetOp::Complement, vec![ca[0].clone(), union(vec![ca[1].clone(), b.clone()])]);
    }
    if let Some(ael) = finite_elements(a) {
        let mut keep = Vec::new();
        let mut unk = Vec::new();
        for x in ael {
            match contains3(b, x) {
                Some(true) => {}
                Some(false) => keep.push(x.clone()),
                None => unk.push(x.clone()),
            }
        }
        let mut parts = vec![finite_set(keep)];
        if !unk.is_empty() {
            parts.push(mk(SetOp::Complement, vec![finite_set(unk), b.clone()]));
        }
        return union(parts);
    }
    mk(SetOp::Complement, vec![a.clone(), b.clone()])
}

fn complement_raw(a: Expr, b: Expr) -> Expr {
    if is_empty_set(&a) {
        return a;
    }
    mk(SetOp::Complement, vec![a, b])
}

/// `ProductSet(*sets)`.
pub fn product_set(sets: Vec<Expr>) -> Expr {
    if sets.iter().any(is_empty_set) {
        return empty_set();
    }
    if sets.is_empty() {
        return finite_set(vec![Expr::from_kind(Kind::Func(Func::Undefined("Tuple".into()), vec![]))]);
    }
    mk(SetOp::ProductSet, sets)
}

/// Single-variable compatibility constructor; new fallible APIs use lambda_checked.
pub fn lambda(var: &Expr, body: &Expr) -> Expr {
    lambda_checked(var,body).expect("valid internal Lambda signature")
}

/// Single-variable compatibility application.
pub fn apply_lambda(lam:&Expr,x:&Expr)->Option<Expr>{lambda_apply(lam,std::slice::from_ref(x)).ok()}

/// `ImageSet(Lambda(var, body), base)`.
pub fn image_set(lam: &Expr, base: &Expr) -> Expr {
    if is_empty_set(base) {
        return empty_set();
    }
    if is_op(lam,SetOp::IdentityFunction) { return base.clone(); }
    mk(SetOp::ImageSet, vec![lam.clone(), base.clone()])
}

/// `ConditionSet(sym, cond, base)`.
pub fn condition_set(sym: &Expr, cond: &Expr, base: &Expr) -> Expr {
    match cond.kind() {
        Kind::Bool(true) => return base.clone(),
        Kind::Bool(false) => return empty_set(),
        _ => {}
    }
    if is_empty_set(base) {
        return empty_set();
    }
    if let Some(el) = finite_elements(base) {
        let mut keep = Vec::new();
        let mut unk = Vec::new();
        for x in el {
            match crate::subs::subs1(cond, sym, x).kind() {
                Kind::Bool(true) => keep.push(x.clone()),
                Kind::Bool(false) => {}
                _ => unk.push(x.clone()),
            }
        }
        if unk.len() == el.len() {
            return mk(SetOp::ConditionSet, vec![sym.clone(), cond.clone(), base.clone()]);
        }
        let mut parts = vec![finite_set(keep)];
        if !unk.is_empty() {
            parts.push(mk(SetOp::ConditionSet, vec![sym.clone(), cond.clone(), finite_set(unk)]));
        }
        return union(parts);
    }
    mk(SetOp::ConditionSet, vec![sym.clone(), cond.clone(), base.clone()])
}

/// `Range(start, stop, step)` of integers.
pub fn range(start: &Expr, stop: &Expr, step: &Expr) -> Expr {
    if let (Some(a), Some(b), Some(s)) = (start.as_i64(), stop.as_i64(), step.as_i64()) {
        if s == 0 {
            return mk(SetOp::Range, vec![start.clone(), stop.clone(), step.clone()]);
        }
        // Widen before subtracting: valid i64 endpoints can span more than i64.
        let (a, b, s) = (a as i128, b as i128, s as i128);
        let n = (b - a + s - s.signum()) / s; // ceil((b - a) / s)
        if n <= 0 {
            return mk(SetOp::Range, vec![Expr::zero(), Expr::zero(), Expr::one()]);
        }
        let integer = |v: i128| Expr::num(crate::number::Number::from_bigint(v.into()));
        return mk(SetOp::Range, vec![integer(a), integer(a + n * s), integer(s)]);
    }
    mk(SetOp::Range, vec![start.clone(), stop.clone(), step.clone()])
}

/// `Contains(x, S)`: evaluates to a boolean when membership is decided.
pub fn contains_expr(x: &Expr, s: &Expr) -> Expr {
    let r = contains(s, x);
    if matches!(r.kind(), Kind::Bool(_)) {
        return r;
    }
    mk(SetOp::Contains, vec![x.clone(), s.clone()])
}

/// Entry point used by [`crate::functions::func`].
pub fn eval(op: SetOp, args: Vec<Expr>) -> Expr {
    let flag = |i: usize| args.get(i).is_some_and(|a| *a == Expr::boolean(true));
    match op {
        SetOp::Interval if (2..=4).contains(&args.len()) => interval(&args[0], &args[1], flag(2), flag(3)),
        SetOp::FiniteSet => finite_set(args),
        SetOp::Union => union(args),
        SetOp::Intersection => intersection(args),
        SetOp::Complement if args.len() == 2 => complement(&args[0], &args[1]),
        SetOp::ProductSet => product_set(args),
        SetOp::ImageSet if args.len() == 2 => image_set(&args[0], &args[1]),
        SetOp::ConditionSet if args.len() == 3 => condition_set(&args[0], &args[1], &args[2]),
        SetOp::ConditionSet if args.len() == 2 => condition_set(&args[0], &args[1], &universal_set()),
        SetOp::Range => match args.len() {
            1 => range(&Expr::zero(), &args[0], &Expr::one()),
            2 => range(&args[0], &args[1], &Expr::one()),
            3 => range(&args[0], &args[1], &args[2]),
            _ => mk(op, args),
        },
        SetOp::Contains if args.len() == 2 => contains_expr(&args[0], &args[1]),
        SetOp::Lambda if args.len()==2=>lambda_checked(&args[0],&args[1]).unwrap_or_else(|_|mk(op,args)),
        SetOp::IdentityFunction if args.is_empty()=>identity_function(),
        _ => mk(op, args),
    }
}

// ----------------------------------------------------------------------
// printing
// ----------------------------------------------------------------------

fn has_finite_set(e: &Expr) -> bool {
    is_op(e, SetOp::FiniteSet) || e.args().iter().any(has_finite_set)
}

/// `str()` of a set node.
pub(crate) fn write_set(op: SetOp, args: &[Expr], out: &mut String) {
    use crate::printer::write_expr;
    let list = |args: &[Expr], out: &mut String| {
        for (i, a) in args.iter().enumerate() {
            if i > 0 {
                out.push_str(", ");
            }
            write_expr(a, out);
        }
    };
    match op {
        SetOp::IdentityFunction => out.push_str("Lambda(_x, _x)"),
        SetOp::Lambda if args.len()==2=>{
            out.push_str("Lambda(");
            let signature=lambda_signature(&mk(op,args.to_vec())).unwrap();
            let vars=signature.args();
            if vars.len()==1 && vars[0].is_symbol(){write_expr(&vars[0],out);}else{write_expr(&signature,out);}
            out.push_str(", ");write_expr(&args[1],out);out.push(')');
        },
        SetOp::Interval if args.len() == 4 => {
            let (a, b) = (&args[0], &args[1]);
            let (l, r) = (args[2] == Expr::boolean(true), args[3] == Expr::boolean(true));
            let m = if (is_infinite(a) && is_infinite(b)) || (is_infinite(a) && !r) || (is_infinite(b) && !l) || (!l && !r) {
                ""
            } else if l && r {
                ".open"
            } else if l {
                ".Lopen"
            } else {
                ".Ropen"
            };
            out.push_str("Interval");
            out.push_str(m);
            out.push('(');
            // SymPy formats the end points with `str()` (full precision floats)
            for (i, p) in args[..2].iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                out.push_str(&crate::printer::to_str(p));
            }
            out.push(')');
        }
        SetOp::Union | SetOp::Intersection => {
            // SymPy's `_print_LatticeOp`: sorted by `default_sort_key`
            let mut v: Vec<Expr> = args.to_vec();
            v.sort_by_cached_key(crate::sortkey::sort_key);
            out.push_str(op.name());
            out.push('(');
            list(&v, out);
            out.push(')');
        }
        SetOp::FiniteSet => {
            let mut v: Vec<Expr> = args.to_vec();
            v.sort_by_cached_key(crate::sortkey::sort_key);
            if v.iter().any(has_finite_set) {
                out.push_str("FiniteSet(");
                list(&v, out);
                out.push(')');
            } else {
                out.push('{');
                list(&v, out);
                out.push('}');
            }
        }
        SetOp::ConditionSet if args.len() == 3 && is_op(&args[2], SetOp::UniversalSet) => {
            out.push_str("ConditionSet(");
            list(&args[..2], out);
            out.push(')');
        }
        _ if op.is_named() && args.is_empty() => out.push_str(op.name()),
        _ => {
            out.push_str(op.name());
            out.push('(');
            list(args, out);
            out.push(')');
        }
    }
}

#[cfg(test)]
mod tests;
