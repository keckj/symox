//! `A - Mod(B, C) = 0` over the integers: `_solve_modular` and
//! `_invert_modular` of the reference `solveset` (SymPy 1.14,
//! `sympy/solvers/solveset.py`).
//!
//! The reference raises `TypeError` when `Abs(rhs) >= Abs(C)` cannot be
//! decided and `NotInvertible` when the coefficient of the symbol has no
//! inverse modulo `C`, and fails with an `AttributeError` for a modular
//! term inside a product or a non-numeric base of a power of the symbol;
//! all are reported as [`ModularError`].
//!
//! Intentional corrections (the reference result is wrong there, the
//! equation is left unsolved as a ConditionSet): a symbolic coefficient of
//! the symbol (`Mod(k*x, 5) = 1` is not `5*Integers`: `k = x = 1` solves
//! it) and a non-integer rational base (`Mod((1/2)**x, 5) = 1` is not
//! `4*Naturals0`: the reference uses the numerator of the base).

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed};

use super::*;
use crate::expr::RelOp;
use crate::special::SpecialFn;
use crate::traverse::has;

/// The exceptions the reference raises out of `_invert_modular`.
#[derive(Debug, Clone, PartialEq)]
pub enum ModularError {
    Type(String),
    NotInvertible(String),
    Attribute(String),
}

fn is_mod(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Special(SpecialFn::Mod), a) if a.len() == 2)
}

fn mod_args(e: &Expr) -> Option<(Expr, Expr)> {
    match e.kind() {
        Kind::Func(Func::Special(SpecialFn::Mod), a) if a.len() == 2 => Some((a[0].clone(), a[1].clone())),
        _ => None,
    }
}

fn add_args(e: &Expr) -> Vec<Expr> {
    if e.is_add() { e.args() } else { vec![e.clone()] }
}

fn mul_args(e: &Expr) -> Vec<Expr> {
    if e.is_mul() { e.args() } else { vec![e.clone()] }
}

fn term_factors(f: &Expr) -> Vec<Expr> {
    add_args(f).iter().flat_map(mul_args).collect()
}

fn mod_atoms(f: &Expr) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    crate::traverse::preorder(f, &mut |s| {
        if is_mod(s) && !out.contains(s) {
            out.push(s.clone());
        }
    });
    out
}

/// The reference `_is_modular`.
fn is_modular(f: &Expr, x: &Expr) -> bool {
    let terms = mod_atoms(f);
    if terms.len() != 1 {
        return false;
    }
    let Some((a, m)) = mod_args(&terms[0]) else { return false };
    has(&a, x) && crate::assumptions::is_integer(&m) == Some(true) && term_factors(f).iter().any(is_mod)
}

/// The coefficient normalization at the start of the reference `_solveset`.
fn normalize(f: &Expr, x: &Expr) -> Expr {
    if !f.is_add() {
        return f.clone();
    }
    let deps = [x.clone()];
    let (a, h) = crate::decompose::as_independent(f, &deps, Some(true));
    let (m, h) = crate::decompose::as_independent(&h, &deps, Some(false));
    let special = m.is_zero() || m == Expr::oo() || m == Expr::neg_oo() || m == Expr::zoo();
    if special {
        return f.clone();
    }
    crate::ops::add2(&crate::ops::div(&a, &m), &h)
}

fn int_of(e: &Expr) -> Option<BigInt> {
    e.as_number().filter(|n| n.is_integer()).and_then(|n| n.to_bigint())
}

fn int(v: BigInt) -> Expr {
    Expr::num(crate::number::Number::from_bigint(v))
}

/// `ImageSet(Lambda(n, body), base)`; a body free of `n` gives
/// `FiniteSet(body)` as the reference constructor does.
fn image(n: &Expr, body: &Expr, base: &Expr) -> Expr {
    if !has(body, n) && !is_empty_set(base) {
        return finite_set(vec![body.clone()]);
    }
    image_set(&lambda(n, body), base)
}

enum Inverted {
    Unsolved,
    Solved(Expr, Expr),
}

/// Python's `pow(b, e, m)` for `m > 0`.
fn pow_mod(b: &BigInt, e: &BigInt, m: &BigInt) -> BigInt {
    b.mod_floor(m).modpow(e, m)
}

/// The reference `_invert_modular(Mod(a, m), rhs, n, symbol)`: the pair
/// `(f_x, g_n)` with `f_x` in `g_n`, or Unsolved for `(modterm, rhs)`.
fn invert_modular(modterm: &Expr, rhs: &Expr, n: &Expr, x: &Expr) -> Result<Inverted, ModularError> {
    let Some((a, m)) = mod_args(modterm) else { return Ok(Inverted::Unsolved) };
    if crate::assumptions::is_integer(rhs) == Some(false) {
        return Ok(Inverted::Solved(x.clone(), empty_set()));
    }
    if crate::assumptions::is_real(rhs) == Some(false)
        || term_factors(&a).iter().any(|t| crate::assumptions::is_real(t) == Some(false))
    {
        return Ok(Inverted::Unsolved);
    }
    let bound = crate::relational::ge(&crate::functions::abs(rhs), &crate::functions::abs(&m));
    match bound.kind() {
        Kind::Bool(true) => return Ok(Inverted::Solved(x.clone(), empty_set())),
        Kind::Bool(false) => {}
        _ => return Err(ModularError::Type(format!("cannot determine truth value of Relational: {bound}"))),
    }
    if a == *x {
        let body = crate::ops::add2(&crate::ops::mul2(&m, n), rhs);
        return Ok(Inverted::Solved(x.clone(), image(n, &body, &integers())));
    }
    let deps = [x.clone()];
    if a.is_add() {
        let (g, h) = crate::decompose::as_independent(&a, &deps, None);
        if !g.is_zero() {
            let indep = crate::ops::sub(rhs, &crate::special::modulo(&g, &m));
            return invert_modular(&crate::special::modulo(&h, &m), &crate::special::modulo(&indep, &m), n, x);
        }
    }
    if a.is_mul() {
        let (g, h) = crate::decompose::as_independent(&a, &deps, None);
        if !g.is_one() {
            let inverse = match (int_of(&g), int_of(&m)) {
                (Some(g), Some(m)) => crate::ntheory::mod_inverse(&g, &m)
                    .ok_or_else(|| ModularError::NotInvertible("zero divisor".into()))?,
                // `invert` over QQ: a fraction has no inverse modulo m
                _ if g.is_rational_number() && m.is_integer_number() => {
                    return Err(ModularError::NotInvertible("zero divisor".into()))
                }
                _ => return Ok(Inverted::Unsolved),
            };
            let indep = crate::ops::mul2(rhs, &int(inverse));
            return invert_modular(&crate::special::modulo(&h, &m), &crate::special::modulo(&indep, &m), n, x);
        }
    }
    if let Kind::Pow(base, expo) = a.kind() {
        if has(expo, x) && !has(base, x) {
            if !base.is_rational_number() {
                let class = if base.is_symbol() { "Symbol" } else if base.is_pow() { "Pow" } else if base.is_mul() { "Mul" } else if base.is_add() { "Add" } else { "Expr" };
                return Err(ModularError::Attribute(format!("'{class}' object has no attribute 'p'")));
            }
            let (Some(mi), Some(ri), Some(bi)) = (int_of(&m), int_of(rhs), int_of(base)) else {
                return Ok(Inverted::Unsolved);
            };
            let mdiv = mi.div_floor(&mi.gcd(&ri));
            if !mdiv.is_positive() {
                return Ok(Inverted::Unsolved);
            }
            let remainder = match crate::ntheory::residue::discrete_log(&mdiv, &ri.mod_floor(&mdiv), &bi.mod_floor(&mdiv), None, None) {
                Ok(r) => r,
                Err(_) => return Ok(Inverted::Unsolved),
            };
            if !mi.is_positive() {
                return Ok(Inverted::Unsolved);
            }
            let mut period = crate::ntheory::totient(&mi);
            let reduced = &mi / mi.gcd(&bi);
            if reduced > BigInt::one() {
                for p in crate::ntheory::divisors(&period) {
                    if pow_mod(&bi, &p, &reduced).is_one() {
                        period = p;
                        break;
                    }
                }
            }
            let body = crate::ops::add2(&crate::ops::mul2(&int(period), n), &int(remainder));
            return Ok(Inverted::Solved(expo.clone(), image(n, &body, &naturals0())));
        }
        if has(base, x) && !has(expo, x) {
            let (Some(mi), Some(ri), Some(ei)) = (int_of(&m), int_of(rhs), int_of(expo)) else {
                return Ok(Inverted::Unsolved);
            };
            if !ei.is_positive() || !mi.is_positive() {
                return Ok(Inverted::Unsolved);
            }
            let a = ri.mod_floor(&mi);
            let roots = if ei == BigInt::from(2) {
                let mut r = crate::ntheory::residue::sqrt_mod_iter(&a, &mi);
                r.sort();
                r.dedup();
                r
            } else {
                match crate::ntheory::residue::nthroot_mod_all(&a, &ei, &mi) {
                    Ok(r) => r,
                    Err(_) => return Ok(Inverted::Unsolved),
                }
            };
            if roots.is_empty() {
                return Ok(Inverted::Solved(x.clone(), empty_set()));
            }
            let parts = roots
                .into_iter()
                .map(|r| image(n, &crate::ops::add2(&crate::ops::mul2(&m, n), &int(r)), &integers()))
                .collect();
            return Ok(Inverted::Solved(base.clone(), union(parts)));
        }
    }
    Ok(Inverted::Unsolved)
}

/// Whether the coefficient of `term` in `f.as_coefficients_dict()` is
/// negative, None when the term is not a bare term of `f`.
fn coefficient_is_negative(f: &Expr, term: &Expr) -> Option<bool> {
    if f == term {
        return Some(false);
    }
    match f.kind() {
        Kind::Add(_, terms) => terms.iter().find(|(t, _)| t == term).map(|(_, c)| c.is_negative()),
        _ => None,
    }
}

/// The reference `_solve_modular` for a normalized modular `f`.
fn solve_modular(f: &Expr, x: &Expr, domain: &Expr) -> Result<Expr, ModularError> {
    let unsolved = || condition_set(x, &crate::relational::eq(f, &Expr::zero()), domain);
    let modterm = mod_atoms(f).remove(0);
    let mut map = rustc_hash::FxHashMap::default();
    map.insert(modterm.clone(), Expr::zero());
    let mut rhs = crate::ops::neg(&crate::subs::xreplace(f, &map));
    match coefficient_is_negative(f, &modterm) {
        Some(true) => rhs = crate::ops::neg(&rhs),
        Some(false) => {}
        // `as_coefficients_dict()[modterm]` is the int 0 there
        None => return Err(ModularError::Attribute("'int' object has no attribute 'is_negative'".into())),
    }
    if is_subset(domain, &integers()) != Some(true) || has(&rhs, x) {
        return Ok(unsolved());
    }
    let n = super::solveset::dummy_n();
    let (f_x, g_n) = match invert_modular(&modterm, &rhs, &n, x)? {
        Inverted::Unsolved => return Ok(unsolved()),
        Inverted::Solved(f_x, g_n) => (f_x, g_n),
    };
    if f_x == *x {
        if is_op(domain, SetOp::Integers) {
            return Ok(g_n);
        }
        return Ok(intersection(vec![domain.clone(), g_n]));
    }
    let Some((SetOp::ImageSet, parts)) = as_set(&g_n) else { return Ok(unsolved()) };
    let (lam, base) = (parts[0].clone(), parts[1].clone());
    let Some((vars, body)) = lambda_parts(&lam) else { return Ok(unsolved()) };
    let [var] = vars.as_slice() else { return Ok(unsolved()) };
    let var = var.clone();
    let sol_set = match super::solveset::solveset(&crate::ops::sub(&f_x, &body), x, &integers()) {
        Ok(s) => s,
        Err(_) => return Ok(unsolved()),
    };
    let sol_set = match finite_elements(&sol_set) {
        Some(elements) => union(elements.iter().map(|s| image(&var, s, &base)).collect()),
        // `ImageSet(Lambda(n, S), base)` with a set `S` free of `n` is
        // `{S}`, which has no integer element
        None if !has(&sol_set, &var) => return Ok(empty_set()),
        None => image(&var, &sol_set, &base),
    };
    Ok(intersection(vec![domain.clone(), sol_set]))
}

/// The reference result of a modular equation `f` (an expression or an
/// `Eq`), or None when `f` is not modular in `x`.
pub fn solveset_modular(f: &Expr, x: &Expr, domain: &Expr) -> Result<Option<Expr>, ModularError> {
    let f = match f.kind() {
        Kind::Rel(RelOp::Eq, a, b) => crate::ops::sub(a, b),
        Kind::Rel(..) | Kind::Bool(_) => return Ok(None),
        _ => f.clone(),
    };
    if !x.is_symbol() || f.is_mul() || !has(&f, x) || matches!(f.kind(), Kind::Func(Func::Set(_) | Func::Logic(_), _)) {
        return Ok(None);
    }
    let f = normalize(&f, x);
    if !is_modular(&f, x) {
        return Ok(None);
    }
    solve_modular(&f, x, domain).map(Some)
}
