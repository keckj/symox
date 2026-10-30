//! The `solveset` steps for radicals and negative bases that the native
//! solver does not take: SymPy 1.14's `_solve_radical` (after its `unrad`
//! test in `_solveset`), and the real inversion of `b**g(x)` for a
//! negative constant base `b` (`_invert_real` with `integer_log`), with the
//! guard that keeps a ConditionSet where the principal logarithm of the
//! native solver would miss the integer exponents.
//!
//! The recursive `solveset` calls, `imageset` and the intersection with a
//! domain call back into the host (they are the public set algebra); the
//! candidates of `_solve_radical` are checked with the native `checksol`
//! of `solve::tsolve`, and an undecided candidate stays in a ConditionSet.

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::sets::{as_set, complement, condition_set, empty_set, finite_elements, finite_set, integers, is_subset, lambda, reals, union, SetOp};
use crate::solve::unrad::{UFail, UnradHost};
use crate::traverse::{free_symbols, has, preorder};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive};

type R<T, E> = Result<T, UFail<E>>;

/// The public calls of the steps.
pub trait StepHost: UnradHost {
    /// `solveset(f, x, domain)`, a set.
    fn solveset(&mut self, f: &Expr, x: &Expr, domain: &Expr) -> R<Expr, Self::Error>;
    /// `imageset(lam, base)`.
    fn imageset(&mut self, lam: &Expr, base: &Expr) -> R<Expr, Self::Error>;
    /// `a.intersect(b)`.
    fn intersect(&mut self, a: &Expr, b: &Expr) -> R<Expr, Self::Error>;
    /// Whether a host error is a ValueError.
    fn is_value_error(&self, e: &Self::Error) -> bool;
    /// The established elimination for an equation where the reference
    /// `unrad` raises NotImplementedError: `None` or `(eq, cov)` whose roots
    /// are a superset of the equation's (every candidate is checked).
    fn unrad_fallback(&mut self, _f: &Expr, _x: &Expr) -> R<Option<(Expr, Vec<Expr>)>, Self::Error> {
        Ok(None)
    }
}

fn is_zero(e: &Expr) -> Option<bool> {
    crate::assumptions::is_zero(e)
}

fn is_finite(e: &Expr) -> Option<bool> {
    crate::assumptions::is_finite(e)
}

fn is_exp(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(Func::Exp, _) => true,
        Kind::Pow(b, _) => b.is_constant(Constant::E),
        _ => false,
    }
}

fn exp_arg(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Func(Func::Exp, a) => a[0].clone(),
        Kind::Pow(_, x) => x.clone(),
        _ => unreachable!(),
    }
}

fn is_rational(e: &Expr) -> bool {
    e.as_number().is_some_and(|n| n.is_rational())
}

fn pows(f: &Expr) -> Vec<(Expr, Expr)> {
    let mut out: Vec<(Expr, Expr)> = Vec::new();
    preorder(f, &mut |n| {
        if let Kind::Pow(b, x) = n.kind() {
            let p = (b.clone(), x.clone());
            if !out.contains(&p) {
                out.push(p);
            }
        }
    });
    out
}

fn is_set_op(e: &Expr, op: SetOp) -> bool {
    as_set(e).is_some_and(|(o, _)| o == op)
}

fn is_empty(e: &Expr) -> bool {
    is_set_op(e, SetOp::EmptySet) || finite_elements(e).is_some_and(|v| v.is_empty())
}

/// `A - B` as the public `Complement`: an unevaluated result is empty when
/// `A` is a subset of `B`.
fn set_minus(a: &Expr, b: &Expr) -> Expr {
    let r = complement(a, b);
    if let Some((SetOp::Complement, args)) = as_set(&r) {
        if is_subset(&args[0], &args[1]) == Some(true) {
            return empty_set();
        }
    }
    r
}

/// The coefficient normalization at the start of `_solveset`.
fn normalize(f: &Expr, x: &Expr) -> Expr {
    if f.is_add() {
        let (a, h) = crate::decompose::as_independent(f, std::slice::from_ref(x), Some(true));
        let (m, h) = crate::decompose::as_independent(&h, std::slice::from_ref(x), Some(false));
        let bad = m.is_zero()
            || m.is_constant(Constant::ComplexInfinity)
            || m.is_constant(Constant::Infinity)
            || m.is_constant(Constant::NegInfinity);
        if !bad {
            return add(&[div(&a, &m), h]);
        }
    }
    f.clone()
}

fn has_radical(f: &Expr, x: &Expr) -> bool {
    pows(f).iter().any(|(b, e)| has(b, x) && is_rational(e) && !e.is_integer_number())
}

/// SymPy's `denoms(f, x)`: the factors of the denominators of every node
/// that contain `x`.
pub fn denoms(f: &Expr, x: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    preorder(f, &mut |n| {
        if crate::sets::is_set(n) || matches!(n.kind(), Kind::Bool(_) | Kind::Rel(..)) {
            return;
        }
        let d = crate::simptools::radsimp::fraction(n, false).1;
        if !d.is_one() {
            for m in if d.is_mul() { crate::traverse::ordered_args(&d) } else { vec![d] } {
                if !out.contains(&m) {
                    out.push(m);
                }
            }
        }
    });
    out.retain(|d| free_symbols(d).contains(x));
    out
}

/// `_invert_real` of `dep = -indep` for `dep = g**(num/den)` with an even
/// `den` and an odd `num`; None keeps the native result.
fn invert_real_radical<H: StepHost>(h: &mut H, indep: &Expr, dep: &Expr, x: &Expr, domain: &Expr) -> R<Option<Expr>, H::Error> {
    if is_subset(domain, &reals()) != Some(true) {
        return Ok(None);
    }
    let Kind::Pow(base, expo) = dep.kind() else { return Ok(None) };
    if has(expo, x) || !is_rational(expo) || expo.is_integer_number() {
        return Ok(None);
    }
    let Some((num, den)) = expo.as_number().and_then(|n| n.numer_denom()) else { return Ok(None) };
    if !(den.is_even() && num.is_odd()) {
        return Ok(None);
    }
    let rhs = neg(indep);
    if crate::assumptions::is_negative(&rhs) == Some(true) {
        return Ok(Some(empty_set()));
    }
    if crate::assumptions::is_nonnegative(&rhs) != Some(true) {
        return Ok(None);
    }
    let inv = Expr::num(Number::from_ratio(BigRational::new(den, num)));
    h.solveset(&sub(base, &pow(&rhs, &inv)), x, domain).map(Some)
}

fn check_finiteset(f: &Expr, x: &Expr, sols: &[Expr]) -> Expr {
    let (mut yes, mut undecided) = (Vec::new(), Vec::new());
    for s in sols {
        // the reference keeps a CRootOf unchecked. Only a candidate that
        // `checksol` proves is kept outright; the others go to the
        // ConditionSet, whose constructor removes those where the equation
        // is provably nonzero. A `False` of the native check is not trusted
        // to drop a candidate: with a parameter it can be valid for some
        // values (test_issue_11534's `-y/sqrt(y**2 + 1)` at y = 0, where the
        // reference keeps a ConditionSet).
        let rootof = matches!(s.kind(), Kind::Func(Func::Undefined(n), a) if crate::poly::rootof::is_rootof_head(n, a.len()));
        if rootof || crate::solve::tsolve::checksol(f, x, s, true, false) == Some(true) {
            yes.push(s.clone());
        } else {
            undecided.push(s.clone());
        }
    }
    let cond = condition_set(x, &crate::relational::eq(f, &Expr::zero()), &finite_set(undecided));
    union(vec![finite_set(yes), cond])
}

fn check_set(f: &Expr, x: &Expr, s: &Expr) -> Expr {
    let Some((op, args)) = as_set(s) else { return s.clone() };
    match op {
        SetOp::EmptySet | SetOp::ConditionSet => s.clone(),
        SetOp::FiniteSet => check_finiteset(f, x, args),
        SetOp::Complement => set_minus(&check_set(f, x, &args[0]), &args[1]),
        SetOp::Union => union(args.iter().map(|a| check_set(f, x, a)).collect()),
        _ => s.clone(),
    }
}

/// The reference branch after an `unrad` without result (and no `Abs`):
/// `_solve_as_rational`, the numerator's solutions minus the
/// denominator's when the denominator has the symbol (a symbol-free
/// denominator is the native polynomial solve, which already ran). A
/// ConditionSet keeps the native result (the reference then tries
/// `_transolve`, which the native solver did).
fn solve_as_rational<H: StepHost>(h: &mut H, f: &Expr, x: &Expr, domain: &Expr) -> R<Option<Expr>, H::Error> {
    let mut has_abs = false;
    preorder(f, &mut |n| has_abs |= matches!(n.kind(), Kind::Func(Func::Abs, _)));
    if has_abs {
        return Ok(None);
    }
    let t = crate::poly::api::together(&crate::solve::tsolve::mexpand(f, true));
    let (g, d) = crate::simptools::radsimp::fraction(&t, false);
    if !has(&d, x) {
        return Ok(None);
    }
    let valid = h.solveset(&g, x, domain)?;
    let invalid = h.solveset(&d, x, domain)?;
    let r = set_minus(&valid, &invalid);
    Ok(if is_set_op(&r, SetOp::ConditionSet) { None } else { Some(r) })
}

/// `radical_step`: SymPy's `_solve_radical` for an equation `f = 0` with a
/// radical of `x` that the native solver leaves as a ConditionSet, after
/// the reference inversion of a single radical or a product. None keeps
/// the native result.
pub fn radical_step<H: StepHost>(h: &mut H, f: &Expr, x: &Expr, domain: &Expr) -> R<Option<Expr>, H::Error> {
    if !x.is_symbol() {
        return Ok(None);
    }
    let f = normalize(f, x);
    if !has_radical(&f, x) {
        return Ok(None);
    }
    // `_invert` only moves the symbol-free terms of a sum, so the
    // reference reaches `unrad` with `equation == f` exactly when the
    // symbol-dependent part is itself a sum.
    let (indep, dep) = crate::decompose::as_independent(&f, std::slice::from_ref(x), Some(true));
    if dep.is_mul() {
        let (g, rest) = crate::decompose::as_independent(&dep, std::slice::from_ref(x), Some(false));
        if !g.is_one() {
            if is_zero(&g) != Some(false) || is_finite(&g) != Some(true) {
                return Ok(None);
            }
            return h.solveset(&add(&[rest, div(&indep, &g)]), x, domain).map(Some);
        }
    } else if !dep.is_add() {
        return invert_real_radical(h, &indep, &dep, x, domain);
    }
    let u = match crate::solve::unrad::unrad(h, &f, std::slice::from_ref(x), false, false) {
        Ok(u) => u,
        // SymPy raises here; the established resultant elimination keeps
        // the roots it finds (they are checked below).
        Err(UFail::NotImplemented(_)) => match h.unrad_fallback(&f, x) {
            Ok(Some(u)) => Some(u),
            Ok(None) | Err(UFail::NotImplemented(_) | UFail::Poly(_)) => return Ok(None),
            Err(UFail::Host(e)) if h.is_value_error(&e) => return Ok(None),
            Err(e) => return Err(e),
        },
        Err(UFail::Poly(_)) => return Ok(None),
        Err(UFail::Host(e)) if h.is_value_error(&e) => return Ok(None),
        Err(e) => return Err(e),
    };
    let Some((mut eq, cov)) = u else { return solve_as_rational(h, &f, x, domain) };
    let result = if cov.is_empty() {
        let main = h.solveset(&eq, x, domain)?;
        let mut sing = Vec::new();
        for g in denoms(&f, x) {
            sing.push(h.solveset(&g, x, domain)?);
        }
        set_minus(&main, &union(sing))
    } else {
        let (mut y, mut yeq) = (cov[0].clone(), cov[1].clone());
        let probe = h.solveset(&sub(&y, &Expr::i()), &y, domain)?;
        if is_empty(&probe) {
            let yreal = Expr::dummy("yreal", crate::symbol::Assumptions::REAL);
            let mut rep = rustc_hash::FxHashMap::default();
            rep.insert(y.clone(), yreal.clone());
            yeq = crate::subs::xreplace(&yeq, &rep);
            eq = crate::subs::xreplace(&eq, &rep);
            y = yreal;
        }
        let g_y_s = h.solveset(&yeq, x, domain)?;
        let f_y_sols = h.solveset(&eq, &y, domain)?;
        // iterating a set that is not finite is not supported here
        let Some(gs) = (if is_empty(&g_y_s) { Some(Vec::new()) } else { finite_elements(&g_y_s).map(|v| v.to_vec()) }) else {
            return Ok(None);
        };
        let mut parts = Vec::new();
        for g in gs {
            parts.push(h.imageset(&lambda(&y, &g), &f_y_sols)?);
        }
        union(parts)
    };
    Ok(Some(check_set(&f, x, &result)))
}

// ---------------------------------------------------------------- negative bases

/// The powers `b**g(x)` of `f` with a constant base that is not known to be
/// positive (or zero).
fn negative_base_pows(f: &Expr, x: &Expr) -> bool {
    pows(f).iter().any(|(b, e)| {
        has(e, x) && !has(b, x) && crate::assumptions::is_positive(b) != Some(true) && is_zero(b) != Some(true)
    })
}

/// SymPy's `_is_finite_with_finite_vars`.
fn finite_with_finite_vars(f: &Expr, domain_real: bool) -> Option<bool> {
    use crate::symbol::Assumptions as A;
    let mut rep = rustc_hash::FxHashMap::default();
    for s in free_symbols(f) {
        let Kind::Sym(sym) = s.kind() else { continue };
        let known = sym.assumptions();
        let mut a = known.with(A::FINITE);
        if domain_real {
            if !known.has(A::REAL) && !known.has(A::NOT_REAL) {
                a = a.with(A::REAL);
            }
        } else if !known.has(A::COMPLEX) {
            a = a.with(A::COMPLEX);
        }
        rep.insert(s.clone(), Expr::dummy("_", a));
    }
    is_finite(&crate::subs::xreplace(f, &rep))
}

fn ln_big(n: &BigInt) -> f64 {
    let bits = n.bits();
    if bits <= 1000 {
        n.to_f64().unwrap_or(f64::INFINITY).ln()
    } else {
        let shift = bits - 64;
        (n >> shift).to_f64().unwrap_or(f64::INFINITY).ln() + shift as f64 * std::f64::consts::LN_2
    }
}

/// The integer `k` with `a**k == m` for positive Rationals `a != 1` and `m`.
fn rational_log(a: &BigRational, m: &BigRational) -> Option<i64> {
    if m.is_one() {
        return Some(0);
    }
    let la = ln_big(a.numer()) - ln_big(a.denom());
    let lm = ln_big(m.numer()) - ln_big(m.denom());
    let k0 = (lm / la).round();
    if !k0.is_finite() || k0.abs() > 1e7 {
        return None;
    }
    let k0 = k0 as i64;
    for k in [k0, k0 - 1, k0 + 1] {
        let p = if k >= 0 { num_traits::pow(a.clone(), k as usize) } else { num_traits::pow(a.recip(), (-k) as usize) };
        if &p == m {
            return Some(k);
        }
    }
    None
}

enum Exps {
    Ints(Vec<BigInt>),
    Parity(i64),
}

fn rat_of(e: &Expr) -> Option<BigRational> {
    let (p, q) = e.as_number().filter(|n| n.is_rational())?.numer_denom()?;
    Some(BigRational::new(p, q))
}

fn has_float(e: &Expr) -> bool {
    let mut f = false;
    preorder(e, &mut |n| {
        if matches!(n.kind(), Kind::Num(m) if !m.is_rational()) {
            f = true;
        }
    });
    f
}

/// The real exponents `e` with `base**e == rhs` for a negative constant
/// `base`, or None when it is not decided.
fn negative_base_exponents(base: &Expr, rhs: &Expr) -> Option<Exps> {
    if !free_symbols(rhs).is_empty() || !free_symbols(base).is_empty() || has(rhs, &Expr::nan()) || has(rhs, &Expr::zoo()) {
        return None;
    }
    if is_zero(rhs) == Some(true) {
        return Some(Exps::Ints(vec![]));
    }
    if is_zero(rhs) != Some(false) || is_finite(rhs) != Some(true) {
        return None;
    }
    if crate::assumptions::is_real(rhs) != Some(true) {
        return None;
    }
    if has_float(base) || has_float(rhs) {
        return None;
    }
    if base == &Expr::int(-1) {
        if rhs.is_one() {
            return Some(Exps::Parity(0));
        }
        if rhs == &Expr::int(-1) {
            return Some(Exps::Parity(1));
        }
        if is_rational(rhs) {
            return Some(Exps::Ints(vec![]));
        }
        return None;
    }
    let a = neg(base);
    let abs_rhs = crate::functions::func(Func::Abs, vec![rhs.clone()]);
    let k: Option<BigInt> = match (rat_of(&a), rat_of(rhs)) {
        (Some(ar), Some(mr)) => rational_log(&ar, &mr.abs()).map(BigInt::from),
        _ => {
            let la = crate::simp::expand_log(&crate::functions::func(Func::Log, vec![abs_rhs.clone()]), true);
            let lb = crate::simp::expand_log(&crate::functions::func(Func::Log, vec![a.clone()]), true);
            let k = crate::poly::api::cancel(&div(&la, &lb)).ok()?;
            let kn = k.as_number().filter(|n| n.is_rational())?.clone();
            if !kn.is_integer() {
                return Some(Exps::Ints(vec![]));
            }
            let k = kn.to_bigint()?;
            let d = sub(&pow(&a, &Expr::num(Number::from_bigint(k.clone()))), &abs_rhs);
            if crate::simptools::methods::equals(&d, &Expr::zero()) != Some(true) {
                return None;
            }
            Some(k)
        }
    };
    let Some(k) = k else { return Some(Exps::Ints(vec![])) };
    let sign = crate::assumptions::is_positive(rhs)?;
    Some(Exps::Ints(if k.is_even() == sign { vec![k] } else { vec![] }))
}

fn domain_real(domain: &Expr) -> bool {
    is_subset(domain, &reals()) == Some(true)
}

/// `negative_base_step`: `_invert_real` of an equation whose symbol sits in
/// the exponent of a power with a negative constant base on a real domain.
/// None keeps the native result.
pub fn negative_base_step<H: StepHost>(h: &mut H, f: &Expr, x: &Expr, domain: &Expr) -> R<Option<Expr>, H::Error> {
    if !crate::sets::is_set(domain) || !x.is_symbol() || !domain_real(domain) {
        return Ok(None);
    }
    if !negative_base_pows(f, x) {
        return Ok(None);
    }
    let unsolved = condition_set(x, &crate::relational::eq(f, &Expr::zero()), domain);
    if f.is_mul() {
        let args = crate::traverse::ordered_args(f);
        if args.iter().all(|m| finite_with_finite_vars(m, true) == Some(true)) {
            let mut parts = Vec::new();
            for m in &args {
                parts.push(h.solveset(m, x, domain)?);
            }
            return Ok(Some(union(parts)));
        }
    }
    let mut lhs = f.clone();
    let mut ys = vec![Expr::zero()];
    loop {
        if is_exp(&lhs) {
            lhs = exp_arg(&lhs);
            ys = ys.iter().map(|y| crate::functions::func(Func::Log, vec![y.clone()])).collect();
        } else if matches!(lhs.kind(), Kind::Func(Func::Log, a) if a.len() == 1) {
            lhs = lhs.args()[0].clone();
            ys = ys.iter().map(|y| crate::functions::func(Func::Exp, vec![y.clone()])).collect();
        } else if lhs.is_add() && !crate::decompose::as_independent(&lhs, std::slice::from_ref(x), None).0.is_zero() {
            let (g, rest) = crate::decompose::as_independent(&lhs, std::slice::from_ref(x), None);
            lhs = rest;
            ys = ys.iter().map(|y| sub(y, &g)).collect();
        } else if lhs.is_mul() && !crate::decompose::as_independent(&lhs, std::slice::from_ref(x), None).0.is_one() {
            let (g, rest) = crate::decompose::as_independent(&lhs, std::slice::from_ref(x), None);
            if is_zero(&g) != Some(false) || is_finite(&g) != Some(true) {
                return Ok(None);
            }
            lhs = rest;
            ys = ys.iter().map(|y| div(y, &g)).collect();
        } else {
            break;
        }
    }
    let (base, expo) = match lhs.kind() {
        Kind::Pow(b, e) if !b.is_constant(Constant::E) && !has(b, x) && has(e, x) => (b.clone(), e.clone()),
        _ => return Ok(None),
    };
    if crate::assumptions::is_negative(&base) != Some(true) {
        let keep = crate::assumptions::is_positive(&base) != Some(true) && is_zero(&base) != Some(true);
        return Ok(if keep { Some(unsolved) } else { None });
    }
    let mut parts = Vec::new();
    for y in &ys {
        let Some(found) = negative_base_exponents(&base, y) else { return Ok(Some(unsolved)) };
        match found {
            Exps::Ints(ks) => {
                for k in ks {
                    parts.push(h.solveset(&sub(&expo, &Expr::num(Number::from_bigint(k))), x, domain)?);
                }
            }
            Exps::Parity(j) => {
                let c = crate::decompose::as_independent(&expo, std::slice::from_ref(x), Some(true)).0;
                let coeff = div(&sub(&expo, &c), x);
                if has(&coeff, x)
                    || crate::assumptions::is_real(&coeff) != Some(true)
                    || is_zero(&coeff) != Some(false)
                    || crate::assumptions::is_real(&c) != Some(true)
                {
                    return Ok(Some(unsolved));
                }
                let n = Expr::dummy("n", crate::symbol::Assumptions::INTEGER);
                let body = div(&sub(&add(&[mul(&[Expr::int(2), n.clone()]), Expr::int(j)]), &c), &coeff);
                let lattice = h.imageset(&lambda(&n, &body), &integers())?;
                parts.push(if domain == &reals() { lattice } else { h.intersect(domain, &lattice)? });
            }
        }
    }
    Ok(Some(if parts.is_empty() { empty_set() } else { union(parts) }))
}

/// `negative_base_guard`: the native result becomes a ConditionSet when
/// `f` has a power of a constant base not known to be positive, on a real
/// domain (the principal logarithm misses the integer exponents of a
/// negative base).
pub fn negative_base_guard(f: &Expr, x: &Expr, domain: &Expr, result: &Expr) -> Expr {
    if !crate::sets::is_set(domain) || !x.is_symbol() || is_set_op(result, SetOp::ConditionSet) {
        return result.clone();
    }
    if !domain_real(domain) || !negative_base_pows(f, x) {
        return result.clone();
    }
    condition_set(x, &crate::relational::eq(f, &Expr::zero()), domain)
}
