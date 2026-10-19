//! Elementary heads completing SymPy's core set: `Mod`, `Rem`, `sech`,
//! `csch`, `asec`, `acsc`, `asech`, `acsch`, `acoth`, `sinc`, `frac`,
//! `KroneckerDelta` and `LeviCivita`.
//!
//! `Mod.eval`, the inverse function constant tables and the `sinc`, `frac`,
//! `KroneckerDelta` and `LeviCivita` rules are ported from SymPy
//! (BSD-3-Clause, see LICENSE-SYMPY).

use crate::decompose::could_extract_minus_sign;
use super::{
    extract_i, is_integer, is_negative, is_nonnegative, is_nonpositive, is_positive,
    is_zero, raw, sf, tuple_items, SpecialFn,
};
use crate::assumptions::{self, Signs};
use crate::bigfloat::BigFloat;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{acos, acosh, acot, asin, asinh, atanh, cosh, floor, im, log, sinh};
use crate::number::Number;
use crate::ops::{add, add2, as_coeff_mul, div, inv, mul, mul2, neg, sqrt, sub};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use std::cell::Cell;

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        Mod => guarded(|| eval_mod(&a[0], &a[1])),
        Rem => eval_rem(&a[0], &a[1]),
        Sech => reciprocal(s, Func::Cosh, true, SpecialFn::Asech, &a[0]),
        Csch => reciprocal(s, Func::Sinh, false, SpecialFn::Acsch, &a[0]),
        Asec => eval_asec(&a[0]),
        Acsc => eval_acsc(&a[0]),
        Asech => eval_asech(&a[0]),
        Acsch => eval_acsch(&a[0]),
        Acoth => eval_acoth(&a[0]),
        Sinc => eval_sinc(&a[0]),
        Frac => eval_frac(&a[0]),
        KroneckerDelta => eval_kronecker(a),
        LeviCivita => eval_levi_civita(a),
        _ => None,
    }
}

thread_local! {static MOD_DEPTH: Cell<u32> = const { Cell::new(0) };}

/// Bound the mutual recursion of the `Mod` canonicalisation.
fn guarded(f: impl FnOnce() -> Option<Expr>) -> Option<Expr> {
    let depth = MOD_DEPTH.with(|d| d.get());
    if depth > 48 {
        return None;
    }
    MOD_DEPTH.with(|d| d.set(depth + 1));
    struct Restore(u32);
    impl Drop for Restore {
        fn drop(&mut self) {
            MOD_DEPTH.with(|d| d.set(self.0));
        }
    }
    let _restore = Restore(depth);
    f()
}

// ----------------------------------------------------------------------
// numbers
// ----------------------------------------------------------------------

fn exact_ratio(n: &Number) -> Option<BigRational> {
    match n {
        Number::Float(f) => f.is_finite().then(|| BigFloat::from_f64(*f).to_ratio_exact()),
        Number::MpFloat(v, _) => Some(v.to_ratio_exact()),
        _ => n.to_ratio(),
    }
}

fn float_prec(p: &Number, q: &Number) -> u32 {
    p.float_precision().unwrap_or(0).max(q.float_precision().unwrap_or(0))
}

fn round_ratio(m: &BigRational, prec: u32) -> Number {
    if m.is_zero() {
        return if prec == 53 { Number::float(0.0) } else { Number::mp_float(BigFloat::from_i64(0), prec) };
    }
    Number::mp_float(BigFloat::from_ratio(m.numer(), m.denom(), prec), prec)
}

/// `p % q` with the sign of `q` (SymPy's `Number.__mod__`). Floats are
/// handled exactly on their binary values and rounded once.
pub(crate) fn num_mod(p: &Number, q: &Number) -> Option<Number> {
    if q.is_zero() {
        return None;
    }
    let (pr, qr) = (exact_ratio(p)?, exact_ratio(q)?);
    let m = &pr - &qr * (&pr / &qr).floor();
    let prec = float_prec(p, q);
    Some(if prec == 0 { Number::from_ratio(m) } else { round_ratio(&m, prec) })
}

/// Truncated quotient `int(p/q)`.
fn trunc_div(p: &Number, q: &Number) -> Option<BigInt> {
    let r = p.div(q)?;
    Some(exact_ratio(&r)?.trunc().to_integer())
}

/// `p - int(p/q)*q` (SymPy's `Rem` on numbers).
pub(crate) fn num_rem(p: &Number, q: &Number) -> Option<Number> {
    if q.is_zero() {
        return None;
    }
    if !p.is_float() && !q.is_float() {
        let (pr, qr) = (p.to_ratio()?, q.to_ratio()?);
        return Some(Number::from_ratio(&pr - &qr * (&pr / &qr).trunc()));
    }
    let t = trunc_div(p, q)?;
    Some(p.sub(&Number::from_bigint(t).mul(q)))
}

fn equal_valued(n: &Number, v: i64) -> bool {
    n.cmp_value(&Number::int(v)) == Some(std::cmp::Ordering::Equal)
}

fn has_float(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_float(),
        _ if e.is_atom() => false,
        _ => e.args().iter().any(has_float),
    }
}

/// Numeric `int(e)` (truncation toward zero) of a real numeric expression.
fn int_value(e: &Expr) -> Option<BigInt> {
    if let Some(n) = e.as_number() {
        return exact_ratio(n).map(|r| r.trunc().to_integer());
    }
    if !crate::functions::is_numeric_value(e) {
        return None;
    }
    let v = crate::evalf_mp::eval_mp(e, 160)?;
    if !v.im.is_zero() {
        return None;
    }
    Some(v.re.to_int_trunc())
}

fn real_f64(e: &Expr) -> Option<f64> {
    if !crate::functions::is_numeric_value(e) {
        return None;
    }
    let z = crate::evalf::eval_c64(e)?;
    (z.im == 0.0 && z.re.is_finite()).then_some(z.re)
}

// ----------------------------------------------------------------------
// Mod
// ----------------------------------------------------------------------

fn md(p: &Expr, q: &Expr) -> Expr {
    sf(SpecialFn::Mod, vec![p.clone(), q.clone()])
}

fn as_mod(e: &Expr) -> Option<(&Expr, &Expr)> {
    match e.kind() {
        Kind::Func(Func::Special(SpecialFn::Mod), a) if a.len() == 2 => Some((&a[0], &a[1])),
        _ => None,
    }
}

fn count_mod(e: &Expr) -> usize {
    let own = usize::from(as_mod(e).is_some());
    if e.is_atom() {
        return own;
    }
    own + e.args().iter().map(count_mod).sum::<usize>()
}

fn is_signed_inf(e: &Expr) -> bool {
    e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity)
}

fn is_nonfinite(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) => false,
        Kind::Const(c) => c.is_infinite(),
        // `oo*I` and similar directed infinities
        Kind::Mul(c, facs)
            if !c.is_zero()
                && facs.iter().any(|(b, x)| is_signed_inf(b) && x.is_one())
                && facs.iter().all(|(b, x)| (is_signed_inf(b) && x.is_one()) || b.is_constant(Constant::I)) =>
        {
            true
        }
        _ => assumptions::is_finite(e) == Some(false),
    }
}

fn same_value(a: &Expr, b: &Expr) -> bool {
    match (a.as_number(), b.as_number()) {
        (Some(x), Some(y)) => x.cmp_value(y) == Some(std::cmp::Ordering::Equal),
        _ => a == b,
    }
}

fn pow_mod(p: &Expr, q: &Expr) -> Option<Expr> {
    let Kind::Pow(b, e) = p.kind() else { return None };
    if !(is_integer(e) == Some(true) && is_positive(e) == Some(true)) {
        return None;
    }
    if is_integer(q) == Some(true) && md(b, q).is_zero() {
        return Some(Expr::zero());
    }
    let (bi, ei, qi) = (super::as_int(b)?, super::as_int(e)?, super::as_int(q)?);
    if qi.is_zero() || ei.is_negative() {
        return None;
    }
    let mut r = bi.modpow(&ei, &qi);
    if !r.is_zero() && (r.is_negative() != qi.is_negative()) {
        r += &qi;
    }
    Some(Expr::num(Number::from_bigint(r)))
}

fn factorial_mod(p: &Expr, q: &Expr) -> Option<Expr> {
    let Kind::Func(Func::Special(SpecialFn::Factorial), a) = p.kind() else { return None };
    let n = a.first()?;
    if !(is_integer(n) == Some(true) && is_nonnegative(n) == Some(true) && is_integer(q) == Some(true)) {
        return None;
    }
    let aq = crate::functions::abs(q);
    let d = sub(&aq, n);
    if is_nonpositive(&d) == Some(true) {
        return Some(Expr::zero());
    }
    if d.is_one() {
        let prime = match super::as_int(&aq) {
            Some(v) => Some(crate::ntheory::isprime(&v)),
            None => assumptions::fact(&aq, "prime"),
        };
        if prime == Some(true) {
            return Some(md(&Expr::minus_one(), q));
        }
        if prime == Some(false) && is_nonnegative(&sub(&aq, &Expr::int(6))) == Some(true) {
            return Some(Expr::zero());
        }
    }
    if let (Some(n), Some(modulus)) = (super::as_int(n), super::as_int(q)) {
        let positive_modulus = modulus.abs();
        let gap = &positive_modulus - &n;
        let residue = if &gap - 1 < n && crate::ntheory::isprime(&positive_modulus) {
            let complement = crate::ntheory::congruence::factorial_mod(&(&gap - 1), &positive_modulus)?;
            let inverse = crate::ntheory::mod_inverse(&complement, &positive_modulus)?;
            if gap.is_odd() { -inverse } else { inverse }
        } else {
            crate::ntheory::congruence::factorial_mod(&n, &positive_modulus)?
        };
        return Some(Expr::num(Number::from_bigint(residue.mod_floor(&modulus))));
    }
    None
}

fn binomial_mod(p: &Expr, q: &Expr) -> Option<Expr> {
    let Kind::Func(Func::Special(SpecialFn::Binomial), args) = p.kind() else { return None };
    let (n, k, modulus) = (super::as_int(&args[0])?, super::as_int(&args[1])?, super::as_int(q)?);
    if modulus.is_zero() { return None; }
    if k.is_negative() { return Some(Expr::zero()); }
    let negative = n.is_negative() && k.is_odd();
    let n = if n.is_negative() { &k - n - 1 } else { n };
    let mut residue = crate::ntheory::congruence::binomial_mod(&n, &k, &modulus.abs()).ok()?;
    if negative { residue = -residue; }
    Some(Expr::num(Number::from_bigint(residue.mod_floor(&modulus))))
}

#[cfg(test)]
mod modular_combinatorial_tests {
    use super::*;

    fn held(head: SpecialFn, args: Vec<Expr>) -> Expr {
        Expr::from_unevaluated_kind(Kind::Func(Func::Special(head), args))
    }

    #[test]
    fn held_binomials_use_modular_arithmetic_with_signed_arguments() {
        for n in -8..20 {
            for k in -2..25 {
                let value = held(SpecialFn::Binomial, vec![Expr::int(n), Expr::int(k)]);
                let exact = crate::ntheory::binomial(&n.into(), &k.into());
                for modulus in [-25, -12, -7, -1, 1, 7, 12, 25] {
                    let expected = exact.mod_floor(&modulus.into());
                    assert_eq!(eval_mod(&value, &Expr::int(modulus)),
                        Some(Expr::num(Number::from_bigint(expected))), "{n} {k} {modulus}");
                }
            }
        }
    }

    #[test]
    fn huge_held_combinatorial_examples_do_not_expand_integer_results() {
        let factorial = held(SpecialFn::Factorial, vec![Expr::int(1_000_000_000)]);
        assert_eq!(eval_mod(&factorial, &Expr::int(1_000_000_009)), Some(Expr::int(712_524_808)));
        let binomial = held(SpecialFn::Binomial, vec![Expr::int(1_000_000_000_000_000_000), Expr::int(1_000_000_000_000)]);
        assert_eq!(eval_mod(&binomial, &Expr::int(10_000_600_009)), Some(Expr::int(3_744_312_326)));
    }
}

/// SymPy's `Mod.eval.number_eval`. `None` also stands for a zero divisor
/// (the Python layer raises `ZeroDivisionError` before reaching here).
fn mod_number_eval(p: &Expr, q: &Expr) -> Option<Expr> {
    if q.is_zero() || is_zero(q) == Some(true) {
        return None;
    }
    if p.is_constant(Constant::NaN) || q.is_constant(Constant::NaN) || is_nonfinite(p) || is_nonfinite(q) {
        return Some(Expr::nan());
    }
    if p.is_zero()
        || same_value(p, q)
        || same_value(p, &neg(q))
        || (is_integer(p) == Some(true) && q.as_number().is_some_and(|n| equal_valued(n, 1)))
    {
        return Some(Expr::zero());
    }
    if let Some(qn) = q.as_number() {
        if let Some(pn) = p.as_number() {
            return num_mod(pn, qn).map(Expr::num);
        }
        if equal_valued(qn, 2) {
            if assumptions::is_even(p) == Some(true) {
                return Some(Expr::zero());
            }
            if assumptions::is_odd(p) == Some(true) {
                return Some(Expr::one());
            }
        }
    }
    if let Some(v) = pow_mod(p, q).or_else(|| factorial_mod(p, q)).or_else(|| binomial_mod(p, q)) {
        return Some(v);
    }
    let r = div(p, q);
    if is_integer(&r) == Some(true) {
        return Some(Expr::zero());
    }
    if let Some(d) = int_value(&r) {
        let rv = sub(p, &mul2(&Expr::num(Number::from_bigint(d)), q));
        if is_negative(&mul2(&rv, q)) == Some(true) {
            return Some(add2(&rv, q));
        }
        return Some(rv);
    }
    // by difference: p - k*q for a small k
    let qpos = if is_positive(q) == Some(true) {
        true
    } else if is_negative(q) == Some(true) {
        false
    } else {
        return None;
    };
    // comp1(a, b): a <= b (q > 0) or a >= b (q < 0); comp2 strict
    let comp1 = |a: &Expr, b: &Expr| {
        let d = if qpos { sub(b, a) } else { sub(a, b) };
        is_nonnegative(&d).or_else(|| integer_bound(&d, false))
    };
    let comp2 = |a: &Expr, b: &Expr| {
        let d = if qpos { sub(b, a) } else { sub(a, b) };
        is_positive(&d).or_else(|| integer_bound(&d, true))
    };
    let mut ls = mul2(&Expr::int(-2), q);
    let r = sub(p, q);
    for _ in 0..4 {
        if comp1(&ls, p) != Some(true) {
            return None;
        }
        if comp2(&r, &ls) == Some(true) {
            return Some(sub(p, &ls));
        }
        ls = add2(&ls, q);
    }
    None
}

/// Lower bound test for `c + sum(k_i*s_i)` with `k_i > 0` and integer
/// symbols `s_i` known positive (>= 1) or nonnegative (>= 0), in the spirit
/// of SymPy's `_monotonic_sign`.  Returns `Some(true)` when the bound proves
/// `e > 0` (`strict`) or `e >= 0`.
fn integer_bound(e: &Expr, strict: bool) -> Option<bool> {
    let Kind::Add(c, terms) = e.kind() else { return None };
    let mut low = exact_ratio(c)?;
    for (t, k) in terms {
        let k = exact_ratio(k)?;
        if !k.is_positive() || !matches!(t.kind(), Kind::Sym(_)) || is_integer(t) != Some(true) {
            return None;
        }
        if is_positive(t) == Some(true) {
            low += k;
        } else if is_nonnegative(t) != Some(true) {
            return None;
        }
    }
    let ok = if strict { low.is_positive() } else { !low.is_negative() };
    ok.then_some(true)
}

fn contains_add(e: &Expr) -> bool {
    match e.kind() {
        Kind::Add(..) => true,
        Kind::Mul(_, f) => f.iter().any(|(b, _)| b.is_add()),
        Kind::Pow(b, _) => b.is_add(),
        _ => false,
    }
}

/// `gcd_terms(e/g, clear=False, fraction=False)` for the shapes `Mod`
/// produces: term-wise division unless `g` involves a sum.
fn gcd_div(e: &Expr, g: &Expr) -> Expr {
    if !contains_add(g) && e.is_add() {
        return add(&e.args().iter().map(|t| div(t, g)).collect::<Vec<_>>());
    }
    div(e, g)
}

fn eval_mod(p0: &Expr, q0: &Expr) -> Option<Expr> {
    if q0.is_zero() {
        return None;
    }
    if let Some(v) = mod_number_eval(p0, q0) {
        return Some(v);
    }
    let mut p = p0.clone();
    let mut q = q0.clone();
    let minus_p = neg(&p);
    if let Some((pin, qin)) = as_mod(&p) {
        if md(qin, &q).is_zero() {
            return Some(md(pin, &q));
        }
        if is_nonnegative(&mul2(qin, &sub(&q, qin))) == Some(true) {
            return Some(p.clone());
        }
    } else if let Some((pin, qin)) = as_mod(&minus_p) {
        if md(qin, &q).is_zero() {
            return Some(md(&neg(pin), &q));
        }
        if is_nonpositive(&mul2(qin, &add2(&q, qin))) == Some(true) {
            return Some(p.clone());
        }
    } else if p.is_add() {
        let (mods, non): (Vec<Expr>, Vec<Expr>) = p.args().into_iter().partition(|t| as_mod(t).is_some());
        if !mods.is_empty() && mods.iter().all(|m| as_mod(m).is_some_and(|(_, qq)| *qq == q)) {
            let inner: Vec<Expr> = mods.iter().map(|m| as_mod(m).unwrap().0.clone()).collect();
            return Some(md(&add2(&add(&non), &add(&inner)), &q));
        }
    } else if p.is_mul() {
        let args = p.args();
        let (mod_l, mut non_mod_l): (Vec<Expr>, Vec<Expr>) = args.iter().cloned().partition(|t| as_mod(t).is_some());
        let all_int = || args.iter().all(|t| is_integer(t) == Some(true));
        if !mod_l.is_empty()
            && mod_l.iter().all(|m| as_mod(m).is_some_and(|(_, qq)| *qq == q))
            && all_int()
            && is_integer(&q) == Some(true)
        {
            let mut mods = Vec::new();
            let mut non = Vec::new();
            for x in &non_mod_l {
                let j = md(x, &q);
                match as_mod(&j) {
                    Some((inner, _)) => mods.push(inner.clone()),
                    None => non.push(j),
                }
            }
            let prod_mod1 = mul(&mod_l.iter().map(|m| as_mod(m).unwrap().0.clone()).collect::<Vec<_>>());
            let net = mul2(&prod_mod1, &mul(&mods));
            return Some(mul2(&mul(&non), &md(&net, &q)));
        }
        if q.is_integer_number() && !q.is_one() && all_int() {
            let qn = q.as_number().unwrap().clone();
            non_mod_l = args
                .iter()
                .map(|i| match i.as_number() {
                    Some(n) if i.is_integer_number() => Expr::num(num_mod(n, &qn).unwrap()),
                    _ => i.clone(),
                })
                .collect();
            if non_mod_l.iter().any(|t| t.is_zero()) {
                return Some(Expr::zero());
            }
        }
        let mut all = non_mod_l;
        all.extend(mod_l);
        p = mul(&all);
    }

    let mut g = if has_float(&p) || has_float(&q) {
        Expr::float(1.0)
    } else {
        crate::poly::api::gcd(&p, &q, &[]).unwrap_or_else(|_| Expr::one())
    };
    if !g.as_number().is_some_and(|n| equal_valued(n, 1)) {
        p = gcd_div(&p, &g);
        q = gcd_div(&q, &g);
    }
    let (pwas, qwas) = (p.clone(), q.clone());
    if p.is_add() {
        let args = p.args();
        let new: Vec<Expr> = args
            .iter()
            .map(|i| {
                let a = md(i, &q);
                if count_mod(&a) > count_mod(i) { i.clone() } else { a }
            })
            .collect();
        if new != args {
            p = add(&new);
        }
    } else {
        let (cp, pr) = as_coeff_mul(&p);
        let (cq, qr) = as_coeff_mul(&q);
        if (cp.is_float() || cq.is_float()) && num_mod(&cp, &cq).is_some_and(|r| r.is_zero()) {
            if let Some(t) = trunc_div(&cp, &cq) {
                g = mul2(&g, &Expr::num(cq.clone()));
                p = mul2(&pr, &Expr::num(Number::from_bigint(t)));
                q = qr;
            }
        }
    }
    if could_extract_minus_sign(&p) && could_extract_minus_sign(&q) {
        g = neg(&g);
        p = neg(&p);
        q = neg(&q);
    }
    if let Some(rv) = mod_number_eval(&p, &q) {
        return Some(mul2(&rv, &g));
    }
    if g.as_number().is_some_and(|n| n.is_float() && equal_valued(n, 1)) {
        return Some(raw(SpecialFn::Mod, vec![mul2(&p, &g), q]));
    }
    if let Kind::Mul(c, facs) = g.kind() {
        if c.is_float() && equal_valued(c, 1) {
            p = mul2(&Expr::num(c.clone()), &p);
            g = mul(&facs.iter().map(|(b, x)| crate::ops::mul_factor(b, x)).collect::<Vec<_>>());
        }
    }
    let m = if p != pwas || q != qwas { md(&p, &q) } else { raw(SpecialFn::Mod, vec![p, q]) };
    Some(mul2(&g, &m))
}

// ----------------------------------------------------------------------
// Rem
// ----------------------------------------------------------------------

fn eval_rem(p: &Expr, q: &Expr) -> Option<Expr> {
    if q.is_zero() || is_zero(q) == Some(true) {
        return None;
    }
    if p.is_constant(Constant::NaN) || q.is_constant(Constant::NaN) || is_nonfinite(p) || is_nonfinite(q) {
        return Some(Expr::nan());
    }
    if p.is_zero()
        || same_value(p, q)
        || same_value(p, &neg(q))
        || (is_integer(p) == Some(true) && q.as_number().is_some_and(|n| equal_valued(n, 1)))
    {
        return Some(Expr::zero());
    }
    match (p.as_number(), q.as_number()) {
        (Some(a), Some(b)) => num_rem(a, b).map(Expr::num),
        _ => None,
    }
}

// ----------------------------------------------------------------------
// sech, csch
// ----------------------------------------------------------------------

fn reciprocal(s: SpecialFn, base: Func, even: bool, inverse: SpecialFn, arg: &Expr) -> Option<Expr> {
    if could_extract_minus_sign(arg) {
        let m = sf(s, vec![neg(arg)]);
        return Some(if even { m } else { neg(&m) });
    }
    match arg.kind() {
        // cosh/sinh(oo) = oo, cosh/sinh(zoo) = nan
        Kind::Const(Constant::Infinity) => return Some(Expr::zero()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        _ => {}
    }
    if let Kind::Func(Func::Special(g), a) = arg.kind() {
        if *g == inverse && a.len() == 1 {
            return Some(a[0].clone());
        }
    }
    let t = if base == Func::Cosh { cosh(arg) } else { sinh(arg) };
    if let Kind::Func(f, a) = t.kind() {
        if *f == base && a.len() == 1 && a[0] == *arg {
            return None;
        }
    }
    Some(inv(&t))
}

// ----------------------------------------------------------------------
// inverse functions: constant tables
// ----------------------------------------------------------------------

type Table = Vec<(Expr, Expr)>;

fn build(rows: &[(&str, &str)]) -> Table {
    rows.iter()
        .filter_map(|(k, v)| Some((crate::parser::parse_str(k).ok()?, crate::parser::parse_str(v).ok()?)))
        .collect()
}

fn lookup(t: &Table, e: &Expr) -> Option<Expr> {
    t.iter().find(|(k, _)| k == e).map(|(_, v)| v.clone())
}

thread_local! {
    static ACSC_TABLE: Table = build(&[
        ("2*sqrt(3)/3", "pi/3"), ("sqrt(2)", "pi/4"), ("sqrt(2*sqrt(5)/5 + 2)", "pi/5"),
        ("1/sqrt(5/8 - sqrt(5)/8)", "pi/5"), ("sqrt(2 - 2*sqrt(5)/5)", "2*pi/5"),
        ("1/sqrt(sqrt(5)/8 + 5/8)", "2*pi/5"), ("2", "pi/6"), ("sqrt(2*sqrt(2) + 4)", "pi/8"),
        ("2/sqrt(2 - sqrt(2))", "pi/8"), ("sqrt(4 - 2*sqrt(2))", "3*pi/8"),
        ("2/sqrt(sqrt(2) + 2)", "3*pi/8"), ("1 + sqrt(5)", "pi/10"), ("-1 + sqrt(5)", "3*pi/10"),
        ("1 - sqrt(5)", "-3*pi/10"), ("sqrt(2) + sqrt(6)", "pi/12"), ("-sqrt(2) + sqrt(6)", "5*pi/12"),
        ("-sqrt(6) + sqrt(2)", "-5*pi/12"),
    ]);
    static ASECH_TABLE: Table = build(&[
        ("I", "log(1 + sqrt(2)) - I*pi/2"), ("-I", "log(1 + sqrt(2)) + I*pi/2"),
        ("-sqrt(2) + sqrt(6)", "pi/12"), ("-sqrt(6) + sqrt(2)", "11*pi/12"),
        ("sqrt(2 - 2*sqrt(5)/5)", "pi/10"), ("-sqrt(2 - 2*sqrt(5)/5)", "9*pi/10"),
        ("2/sqrt(sqrt(2) + 2)", "pi/8"), ("-2/sqrt(sqrt(2) + 2)", "7*pi/8"),
        ("2*sqrt(3)/3", "pi/6"), ("-2*sqrt(3)/3", "5*pi/6"), ("-1 + sqrt(5)", "pi/5"),
        ("1 - sqrt(5)", "4*pi/5"), ("sqrt(2)", "pi/4"), ("-sqrt(2)", "3*pi/4"),
        ("sqrt(2*sqrt(5)/5 + 2)", "3*pi/10"), ("-sqrt(2*sqrt(5)/5 + 2)", "7*pi/10"),
        ("2", "pi/3"), ("-2", "2*pi/3"), ("sqrt(2*sqrt(2) + 4)", "3*pi/8"),
        ("-sqrt(2*sqrt(2) + 4)", "5*pi/8"), ("1 + sqrt(5)", "2*pi/5"), ("-sqrt(5) - 1", "3*pi/5"),
        ("sqrt(2) + sqrt(6)", "5*pi/12"), ("-sqrt(6) - sqrt(2)", "7*pi/12"),
        ("oo*I", "-I*pi/2"), ("-oo*I", "I*pi/2"),
    ]);
    static ACSCH_TABLE: Table = build(&[
        ("I", "-pi/2"), ("I*(sqrt(2) + sqrt(6))", "-pi/12"), ("I*(1 + sqrt(5))", "-pi/10"),
        ("2*I/sqrt(2 - sqrt(2))", "-pi/8"), ("2*I", "-pi/6"), ("I*sqrt(2*sqrt(5)/5 + 2)", "-pi/5"),
        ("sqrt(2)*I", "-pi/4"), ("I*(-1 + sqrt(5))", "-3*pi/10"), ("2*sqrt(3)*I/3", "-pi/3"),
        ("2*I/sqrt(sqrt(2) + 2)", "-3*pi/8"), ("I*sqrt(2 - 2*sqrt(5)/5)", "-2*pi/5"),
        ("I*(-sqrt(2) + sqrt(6))", "-5*pi/12"), ("2", "-I*log(1/2 + sqrt(5)/2)"),
    ]);
}

fn half_pi() -> Expr {
    div(&Expr::pi(), &Expr::int(2))
}

fn is_unsigned_inf(e: &Expr) -> bool {
    e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity) || e.is_constant(Constant::ComplexInfinity)
}

fn func_arg(e: &Expr, f: Func) -> Option<&Expr> {
    match e.kind() {
        Kind::Func(g, a) if *g == f && a.len() == 1 => Some(&a[0]),
        _ => None,
    }
}

fn gt(a: &Expr, b: &Expr) -> bool {
    matches!((real_f64(a), real_f64(b)), (Some(x), Some(y)) if x > y)
}

fn comparable(e: &Expr) -> bool {
    real_f64(e).is_some()
}

fn eval_asec(arg: &Expr) -> Option<Expr> {
    if is_zero(arg) == Some(true) {
        return Some(Expr::zoo());
    }
    if arg.is_one() {
        return Some(Expr::zero());
    }
    if arg.is_minus_one() {
        return Some(Expr::pi());
    }
    if is_unsigned_inf(arg) {
        return Some(half_pi());
    }
    if crate::functions::is_numeric_value(arg) {
        if let Some(v) = ACSC_TABLE.with(|t| lookup(t, arg)) {
            return Some(sub(&half_pi(), &v));
        }
        if let Some(v) = ACSC_TABLE.with(|t| lookup(t, &neg(arg))) {
            return Some(add2(&half_pi(), &v));
        }
    }
    if is_nonfinite(arg) {
        return Some(half_pi());
    }
    let (narg, minus) = match arg.kind() {
        Kind::Mul(c, f) if c.is_minus_one() && f.len() == 1 && f[0].1.is_one() => (f[0].0.clone(), true),
        _ => (arg.clone(), false),
    };
    if let Some(ang) = func_arg(&narg, Func::Sec) {
        if comparable(ang) {
            let mut ang = if minus { sub(&Expr::pi(), ang) } else { ang.clone() };
            ang = md(&ang, &mul2(&Expr::int(2), &Expr::pi()));
            if gt(&ang, &Expr::pi()) {
                ang = sub(&mul2(&Expr::int(2), &Expr::pi()), &ang);
            }
            return Some(ang);
        }
    }
    if let Some(ang) = func_arg(&narg, Func::Csc) {
        if comparable(ang) {
            return Some(sub(&half_pi(), &sf(SpecialFn::Acsc, vec![narg.clone()])));
        }
    }
    None
}

fn eval_acsc(arg: &Expr) -> Option<Expr> {
    if is_zero(arg) == Some(true) {
        return Some(Expr::zoo());
    }
    if arg.is_one() {
        return Some(half_pi());
    }
    if arg.is_minus_one() {
        return Some(neg(&half_pi()));
    }
    if is_unsigned_inf(arg) {
        return Some(Expr::zero());
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&sf(SpecialFn::Acsc, vec![neg(arg)])));
    }
    if is_nonfinite(arg) {
        return Some(Expr::zero());
    }
    if crate::functions::is_numeric_value(arg) {
        if let Some(v) = ACSC_TABLE.with(|t| lookup(t, arg)) {
            return Some(v);
        }
    }
    if let Some(ang) = func_arg(arg, Func::Csc) {
        if comparable(ang) {
            let pi = Expr::pi();
            let mut ang = md(ang, &mul2(&Expr::int(2), &pi));
            if gt(&ang, &pi) {
                ang = sub(&pi, &ang);
            }
            if gt(&ang, &half_pi()) {
                ang = sub(&pi, &ang);
            }
            if gt(&neg(&half_pi()), &ang) {
                ang = sub(&neg(&pi), &ang);
            }
            return Some(ang);
        }
    }
    if let Some(ang) = func_arg(arg, Func::Sec) {
        if comparable(ang) {
            return Some(sub(&half_pi(), &sf(SpecialFn::Asec, vec![arg.clone()])));
        }
    }
    None
}

fn eval_asech(arg: &Expr) -> Option<Expr> {
    let i_half_pi = || mul2(&Expr::i(), &half_pi());
    if arg.is_number() {
        if arg.is_constant(Constant::Infinity) || arg.is_constant(Constant::NegInfinity) {
            return Some(i_half_pi());
        }
        if arg.is_zero() {
            return Some(Expr::oo());
        }
        if arg.is_one() {
            return Some(Expr::zero());
        }
        if arg.is_minus_one() {
            return Some(mul2(&Expr::pi(), &Expr::i()));
        }
        if arg.as_number().is_some_and(|n| n.is_zero()) {
            return Some(Expr::oo());
        }
    }
    if arg.is_constant(Constant::Infinity) || arg.is_constant(Constant::NegInfinity) {
        return Some(i_half_pi());
    }
    if crate::functions::is_numeric_value(arg) || is_imag_inf(arg) {
        if let Some(v) = ASECH_TABLE.with(|t| lookup(t, arg)) {
            return Some(if assumptions::is_real(arg) == Some(true) { mul2(&v, &Expr::i()) } else { v });
        }
    }
    if is_zero(arg) == Some(true) {
        return Some(Expr::oo());
    }
    None
}

fn is_imag_inf(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Mul(_, f) if f.iter().any(|(b, _)| is_signed_inf(b)))
}

fn eval_acsch(arg: &Expr) -> Option<Expr> {
    let l = || log(&add2(&Expr::one(), &sqrt(&Expr::int(2))));
    if arg.is_constant(Constant::Infinity) || arg.is_constant(Constant::NegInfinity) {
        return Some(Expr::zero());
    }
    if arg.as_number().is_some_and(|n| n.is_zero()) {
        return Some(Expr::zoo());
    }
    if arg.is_one() {
        return Some(l());
    }
    if arg.is_minus_one() {
        return Some(neg(&l()));
    }
    if crate::functions::is_numeric_value(arg) {
        if let Some(v) = ACSCH_TABLE.with(|t| lookup(t, arg)) {
            return Some(mul2(&v, &Expr::i()));
        }
    }
    if arg.is_constant(Constant::ComplexInfinity) || is_nonfinite(arg) {
        return Some(Expr::zero());
    }
    if is_zero(arg) == Some(true) {
        return Some(Expr::zoo());
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&sf(SpecialFn::Acsch, vec![neg(arg)])));
    }
    None
}

fn eval_acoth(arg: &Expr) -> Option<Expr> {
    let i_half_pi = || mul2(&Expr::i(), &half_pi());
    if let Some(n) = arg.as_number() {
        if n.is_zero() {
            return Some(i_half_pi());
        }
        if arg.is_one() {
            return Some(Expr::oo());
        }
        if arg.is_minus_one() {
            return Some(Expr::neg_oo());
        }
        if n.is_negative() {
            return Some(neg(&sf(SpecialFn::Acoth, vec![Expr::num(n.neg())])));
        }
        return None;
    }
    if arg.is_constant(Constant::Infinity) || arg.is_constant(Constant::NegInfinity) {
        return Some(Expr::zero());
    }
    if arg.is_constant(Constant::ComplexInfinity) {
        return Some(Expr::zero());
    }
    if assumptions::is_real(arg) != Some(true) {
        if let Some(c) = extract_i(arg) {
            if !has_i(&c) {
                return Some(neg(&mul2(&Expr::i(), &acot(&c))));
            }
        }
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&sf(SpecialFn::Acoth, vec![neg(arg)])));
    }
    if is_zero(arg) == Some(true) {
        return Some(i_half_pi());
    }
    None
}

fn has_i(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::I) => true,
        _ if e.is_atom() => false,
        _ => e.args().iter().any(has_i),
    }
}

// ----------------------------------------------------------------------
// sinc
// ----------------------------------------------------------------------

/// `(c, x)` with `e == c*pi*x`, `c` rational (floats with half-integer
/// multiples converted) and `x` free of `pi`.
fn pi_split(e: &Expr) -> Option<(Number, Expr)> {
    match e.kind() {
        Kind::Const(Constant::Pi) => Some((Number::ONE, Expr::one())),
        Kind::Mul(c, facs) => {
            let idx = facs.iter().position(|(b, x)| b.is_constant(Constant::Pi) && x.is_one())?;
            let rest: Vec<Expr> =
                facs.iter().enumerate().filter(|(j, _)| *j != idx).map(|(_, (b, x))| crate::ops::mul_factor(b, x)).collect();
            let c = if c.is_float() {
                let twice = exact_ratio(&c.mul(&Number::int(2)))?;
                if !twice.is_integer() {
                    return None;
                }
                Number::from_ratio(twice / BigInt::from(2))
            } else {
                c.clone()
            };
            Some((c, mul(&rest)))
        }
        _ => None,
    }
}

fn eval_sinc(arg: &Expr) -> Option<Expr> {
    if is_zero(arg) == Some(true) {
        return Some(Expr::one());
    }
    if arg.is_constant(Constant::Infinity) || arg.is_constant(Constant::NegInfinity) {
        return Some(Expr::zero());
    }
    if arg.is_constant(Constant::ComplexInfinity) {
        return Some(Expr::nan());
    }
    if could_extract_minus_sign(arg) {
        return Some(sf(SpecialFn::Sinc, vec![neg(arg)]));
    }
    let (c, x) = pi_split(arg)?;
    if x.is_one() {
        if c.is_integer() {
            return Some(Expr::zero());
        }
        let twice = c.mul(&Number::int(2));
        if twice.is_integer() {
            // (-1)**(c - 1/2) / arg
            let k = exact_ratio(&c)? - BigRational::new(BigInt::one(), BigInt::from(2));
            let odd = (k.to_integer() % BigInt::from(2)) != BigInt::zero();
            let s = if odd { Expr::minus_one() } else { Expr::one() };
            return Some(div(&s, arg));
        }
        return None;
    }
    if c.is_integer() && is_integer(&x) == Some(true) && is_zero(arg) == Some(false) {
        return Some(Expr::zero());
    }
    None
}

// ----------------------------------------------------------------------
// frac
// ----------------------------------------------------------------------

fn eval_frac(arg: &Expr) -> Option<Expr> {
    let terms: Vec<Expr> = if arg.is_add() { arg.args() } else { vec![arg.clone()] };
    let mut real = Vec::new();
    let mut imag = Vec::new();
    for t in terms {
        let imaginary = assumptions::fact(&t, "imaginary") == Some(true)
            || assumptions::is_real(&mul2(&Expr::i(), &t)) == Some(true);
        if imaginary {
            let i = im(&t);
            if !has_i(&i) {
                imag.push(i);
            } else {
                real.push(t);
            }
        } else {
            real.push(t);
        }
    }
    fn part(a: &Expr) -> Option<Expr> {
        if a.is_constant(Constant::Infinity) || a.is_constant(Constant::NegInfinity) {
            return None;
        }
        if is_integer(a) == Some(true) {
            return Some(Expr::zero());
        }
        if a.is_constant(Constant::NaN) || a.is_constant(Constant::ComplexInfinity) {
            return Some(Expr::nan());
        }
        if crate::functions::is_numeric_value(a) {
            return Some(sub(a, &floor(a)));
        }
        Some(raw(SpecialFn::Frac, vec![a.clone()]))
    }
    let r = part(&add(&real))?;
    let i = part(&add(&imag))?;
    let out = add2(&r, &mul2(&Expr::i(), &i));
    if let Kind::Func(Func::Special(SpecialFn::Frac), a) = out.kind() {
        if a[0] == *arg {
            return None;
        }
    }
    Some(out)
}

// ----------------------------------------------------------------------
// KroneckerDelta, LeviCivita
// ----------------------------------------------------------------------

fn eval_kronecker(a: &[Expr]) -> Option<Expr> {
    let (i, j) = (&a[0], &a[1]);
    if let Some(r) = a.get(2).and_then(|r| tuple_items(r)) {
        if r.len() == 2 {
            let (lo, hi) = (&r[0], &r[1]);
            if is_positive(&sub(lo, i)) == Some(true)
                || is_positive(&sub(lo, j)) == Some(true)
                || is_negative(&sub(hi, i)) == Some(true)
                || is_negative(&sub(hi, j)) == Some(true)
            {
                return Some(Expr::zero());
            }
        }
    }
    let d = sub(i, j);
    match is_zero(&d) {
        Some(true) => return Some(Expr::one()),
        Some(false) => return Some(Expr::zero()),
        None => {}
    }
    if crate::printing::order::sort_key(j) < crate::printing::order::sort_key(i) {
        let mut v = vec![j.clone(), i.clone()];
        v.extend(a.get(2).cloned());
        return Some(sf(SpecialFn::KroneckerDelta, v));
    }
    None
}

fn eval_levi_civita(a: &[Expr]) -> Option<Expr> {
    if a.iter().all(|x| x.is_integer_number()) {
        let n = a.len();
        let mut terms = Vec::new();
        for i in 0..n {
            let mut f: Vec<Expr> = (i + 1..n).map(|j| sub(&a[j], &a[i])).collect();
            f.push(inv(&super::factorial(&Expr::int(i as i64))));
            terms.push(mul(&f));
        }
        return Some(mul(&terms));
    }
    for i in 0..a.len() {
        if a[i + 1..].contains(&a[i]) {
            return Some(Expr::zero());
        }
    }
    None
}

// ----------------------------------------------------------------------
// assumptions hooks
// ----------------------------------------------------------------------

/// Sign set of a real, finite application (see `assumptions::real_signs`).
pub(crate) fn signs(s: SpecialFn, a: &[Expr]) -> Option<Signs> {
    use SpecialFn::*;
    let x = a.first()?;
    match s {
        Mod if a.len() == 2 => {
            if is_positive(&a[1]) == Some(true) {
                Some(Signs::NONNEG)
            } else if is_negative(&a[1]) == Some(true) {
                Some(Signs::NONPOS)
            } else {
                None
            }
        }
        Sech => assumptions::real_signs(x).map(|_| Signs::POS),
        Csch | Acsch => match assumptions::real_signs(x)? {
            sg if sg == Signs::POS => Some(if s == Csch { Signs::POS } else { Signs::ANY }),
            sg if sg == Signs::NEG => Some(if s == Csch { Signs::NEG } else { Signs::ANY }),
            _ => None,
        },
        Sinc => assumptions::real_signs(x).map(|_| Signs::ANY),
        Frac => assumptions::real_signs(x).map(|_| Signs::NONNEG),
        KroneckerDelta | LeviCivita => Some(Signs::ANY),
        _ => None,
    }
}

pub(crate) fn integer(s: SpecialFn, a: &[Expr]) -> Option<bool> {
    use SpecialFn::*;
    match s {
        Mod if a.len() == 2 => {
            (is_integer(&a[0]) == Some(true) && is_integer(&a[1]) == Some(true) && is_zero(&a[1]) == Some(false))
                .then_some(true)
        }
        Frac => is_integer(a.first()?),
        KroneckerDelta | LeviCivita => Some(true),
        _ => None,
    }
}

pub(crate) fn finite(s: SpecialFn, a: &[Expr]) -> Option<bool> {
    use SpecialFn::*;
    match s {
        Frac | KroneckerDelta | LeviCivita => Some(true),
        Asech => {
            let x = a.first()?;
            (is_positive(x) == Some(true) || is_negative(x) == Some(true)).then_some(true)
        }
        _ => None,
    }
}

// ----------------------------------------------------------------------
// numerics
// ----------------------------------------------------------------------

/// Arbitrary precision evaluation: rewrite to heads `evalf_mp` evaluates.
pub(crate) fn mp_rewrite(s: SpecialFn, args: &[Expr], prec: u32) -> Option<Expr> {
    use SpecialFn::*;
    let vals: Vec<Expr> = args
        .iter()
        .map(|a| crate::evalf_mp::ev(a, prec + 20).map(|v| crate::evalf_mp::to_expr(v, prec + 20)))
        .collect::<Option<_>>()?;
    let x = &vals[0];
    let one = Expr::one();
    let real = x.as_number().cloned();
    let big = real.as_ref().is_some_and(|n| n.abs().cmp_value(&Number::ONE) != Some(std::cmp::Ordering::Less));
    let w = || inv(x);
    let s1 = |w: &Expr| sqrt(&sub(&one, &mul2(w, w)));
    Some(match s {
        Mod => {
            let (p, q) = (&vals[0], &vals[1]);
            sub(p, &mul2(q, &floor(&div(p, q))))
        }
        Rem => {
            let (p, q) = (&vals[0], &vals[1]);
            let r = div(p, q);
            let t = if is_negative(&r) == Some(true) { neg(&floor(&neg(&r))) } else { floor(&r) };
            sub(p, &mul2(q, &t))
        }
        Sech => inv(&cosh(x)),
        Csch => inv(&sinh(x)),
        Asec if big => acos(&w()),
        Asec => {
            let w = w();
            mul2(&neg(&Expr::i()), &log(&add2(&w, &mul2(&Expr::i(), &s1(&w)))))
        }
        Acsc if big => asin(&w()),
        Acsc => {
            let w = w();
            mul2(&neg(&Expr::i()), &log(&add2(&mul2(&Expr::i(), &w), &s1(&w))))
        }
        Asech if real.as_ref().is_some_and(|n| n.is_positive() && !big || n.is_one()) => acosh(&w()),
        Asech => {
            let w = w();
            log(&add2(&w, &mul2(&sqrt(&add2(&w, &one)), &sqrt(&sub(&w, &one)))))
        }
        Acsch if real.is_some() => asinh(&w()),
        Acsch => {
            let w = w();
            log(&add2(&w, &sqrt(&add2(&mul2(&w, &w), &one))))
        }
        Acoth if real.is_some() && big => atanh(&w()),
        Acoth => {
            let w = w();
            div(&sub(&log(&add2(&one, &w)), &log(&sub(&one, &w))), &Expr::int(2))
        }
        Sinc if x.is_zero() => one,
        Sinc => div(&crate::functions::sin(x), x),
        Frac => {
            let (re, imv) = (crate::functions::re(x), im(x));
            add2(&sub(&re, &floor(&re)), &mul2(&Expr::i(), &sub(&imv, &floor(&imv))))
        }
        KroneckerDelta => {
            if vals.len() > 2 {
                return None;
            }
            if sub(&vals[0], &vals[1]).is_zero() { one } else { Expr::zero() }
        }
        _ => return None,
    })
}
