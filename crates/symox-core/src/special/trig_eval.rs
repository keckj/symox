//! Automatic evaluation of the trigonometric, hyperbolic, inverse,
//! exponential and logarithmic heads.
//!
//! The rules are a port of the `eval` classmethods of SymPy's
//! `sympy.functions.elementary` (`sin`, `cos`, `tan`, `cot`, the
//! reciprocal functions, `asin` ... `acot`, `atan2`, `sinh` ... `coth`,
//! `asinh`, `acosh`, `atanh`, `exp`, `log`) together with the helpers
//! `_pi_coeff`, `_peeloff_pi`, `_peeloff_ipi` and the constant tables
//! (BSD-3-Clause, see LICENSE-SYMPY).

use super::SpecialFn;
use crate::assumptions;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{
    acos, acot, asin, asinh, atan, atanh, cos, cosh, cot, exp, func, heaviside, is_numeric_value, log, sin,
    sinh, tan, tanh, coth, abs, floor, re, im,
};
use crate::number::Number;
use crate::decompose::could_extract_minus_sign;
use crate::ops::{add, add2, div, inv, mul, mul2, neg, pow, sqrt, sub};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};
use std::cmp::Ordering;

// ----------------------------------------------------------------------
// small helpers
// ----------------------------------------------------------------------

fn pi() -> Expr {
    Expr::pi()
}

fn ii() -> Expr {
    Expr::i()
}

fn num(k: &Number) -> Expr {
    Expr::num(k.clone())
}

fn q(p: i64, d: i64) -> Number {
    Number::rational_i64(p, d).unwrap()
}

fn kpi(k: &Number) -> Expr {
    mul2(&num(k), &pi())
}

fn ipi() -> Expr {
    mul2(&ii(), &pi())
}

fn half_pi() -> Expr {
    div(&pi(), &Expr::int(2))
}

/// `k mod m` for rationals, in `[0, m)`.
fn nmod(k: &Number, m: &Number) -> Number {
    k.sub(&k.div(m).unwrap().floor().mul(m))
}

fn ncmp(a: &Number, b: &Number) -> Option<Ordering> {
    a.cmp_value(b)
}

fn yes(v: Option<bool>) -> bool {
    v == Some(true)
}

fn is_int(e: &Expr) -> bool {
    int3(e) == Some(true)
}

/// Integrality with the parity reasoning SymPy's assumptions perform on
/// sums and products of integers (`n/2` with `n` odd is not an integer).
fn int3(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) if n.is_float() => Some(false),
        Kind::Num(n) => Some(n.is_integer()),
        Kind::Add(c, terms) if !c.is_float() => {
            let mut all = true;
            for (t, k) in terms {
                if !(k.is_integer() && int3(t) == Some(true)) {
                    all = false;
                    break;
                }
            }
            if all {
                return Some(c.is_integer());
            }
            assumptions::is_integer(e)
        }
        Kind::Mul(c, facs) if c.is_rational() && !c.is_integer() => {
            let fs: Vec<Expr> = facs.iter().map(|(b, x)| crate::ops::mul_factor(b, x)).collect();
            if fs.iter().all(|f| int3(f) == Some(true)) {
                let (_, d) = c.numer_denom()?;
                if fs.iter().all(|f| even3(f) == Some(false)) && d.is_even() {
                    return Some(false);
                }
                if d == BigInt::from(2) && fs.iter().any(|f| even3(f) == Some(true)) {
                    return Some(true);
                }
            }
            assumptions::is_integer(e)
        }
        _ => assumptions::is_integer(e),
    }
}

/// SymPy's `is_even` with the same structural reasoning as [`int3`].
fn even3(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) if n.is_float() => None,
        Kind::Num(n) if !n.is_integer() => Some(false),
        Kind::Num(n) => n.to_bigint().map(|b| b.is_even()),
        _ => {
            if int3(e) == Some(false) {
                return Some(false);
            }
            match e.kind() {
                Kind::Add(c, terms) if c.is_integer() => {
                    let mut odd_count = if c.to_bigint().is_some_and(|b| b.is_odd()) { 1 } else { 0 };
                    for (t, k) in terms {
                        if !k.is_integer() {
                            return assumptions::is_even(e);
                        }
                        let kb = k.to_bigint()?;
                        if kb.is_even() {
                            if int3(t) != Some(true) {
                                return assumptions::is_even(e);
                            }
                            continue;
                        }
                        match even3(t) {
                            Some(true) => {}
                            Some(false) if int3(t) == Some(true) => odd_count += 1,
                            _ => return assumptions::is_even(e),
                        }
                    }
                    Some(odd_count % 2 == 0)
                }
                Kind::Mul(c, facs) if c.is_integer() => {
                    let fs: Vec<Expr> = facs.iter().map(|(b, x)| crate::ops::mul_factor(b, x)).collect();
                    if !fs.iter().all(|f| int3(f) == Some(true)) {
                        return assumptions::is_even(e);
                    }
                    if c.to_bigint().is_some_and(|b| b.is_even()) || fs.iter().any(|f| even3(f) == Some(true)) {
                        return Some(true);
                    }
                    if fs.iter().all(|f| even3(f) == Some(false)) {
                        return Some(false);
                    }
                    assumptions::is_even(e)
                }
                _ => assumptions::is_even(e),
            }
        }
    }
}

/// SymPy's `is_even` (three valued).
fn even(e: &Expr) -> Option<bool> {
    even3(e)
}

/// SymPy's `is_odd` (three valued).
fn odd(e: &Expr) -> Option<bool> {
    match int3(e) {
        Some(true) => even3(e).map(|b| !b),
        Some(false) => Some(false),
        None => assumptions::is_odd(e),
    }
}

/// The only argument of `e` when `e` is an `f` node.
fn arg_of<'a>(e: &'a Expr, f: &Func) -> Option<&'a Expr> {
    match e.kind() {
        Kind::Func(g, a) if g == f && a.len() == 1 => Some(&a[0]),
        _ => None,
    }
}

fn special_arg(e: &Expr, s: SpecialFn) -> Option<&Expr> {
    arg_of(e, &Func::Special(s))
}

fn is_head(e: &Expr, f: &Func) -> bool {
    matches!(e.kind(), Kind::Func(g, _) if g == f)
}

fn call(f: &Func, a: &Expr) -> Expr {
    func(f.clone(), vec![a.clone()])
}

/// `cls.eval(arg)`: the constructor result, `None` when it stays `cls(arg)`.
fn eval_of(f: &Func, a: &Expr) -> Option<Expr> {
    let r = call(f, a);
    match r.kind() {
        Kind::Func(g, args) if g == f && args.len() == 1 && args[0] == *a => None,
        _ => Some(r),
    }
}

/// SymPy's `is_number` restricted to finite values.
fn is_num(e: &Expr) -> bool {
    is_numeric_value(e)
}

/// Real numeric value of a comparable expression.
fn real_f64(e: &Expr) -> Option<f64> {
    if !is_num(e) {
        return None;
    }
    crate::evalf::eval_f64(e).filter(|v| v.is_finite())
}

fn comparable(e: &Expr) -> bool {
    real_f64(e).is_some()
}

fn gt(a: &Expr, b: &Expr) -> bool {
    if a == b {
        return false;
    }
    matches!((real_f64(a), real_f64(b)), (Some(x), Some(y)) if x > y)
}

fn lt(a: &Expr, b: &Expr) -> bool {
    gt(b, a)
}

fn md(a: &Expr, m: &Expr) -> Expr {
    super::modulo(a, m)
}

fn two_pi() -> Expr {
    mul2(&Expr::int(2), &pi())
}

/// `Expr.as_coefficient(I)` (`_imaginary_unit_as_coefficient`).
fn i_coeff(arg: &Expr) -> Option<Expr> {
    fn term(t: &Expr) -> Option<Expr> {
        match t.kind() {
            Kind::Const(Constant::I) => Some(Expr::one()),
            Kind::Mul(c, facs) => {
                let idx = facs.iter().position(|(b, x)| b.is_constant(Constant::I) && x.is_one())?;
                let mut v: Vec<Expr> = Vec::with_capacity(facs.len());
                v.push(Expr::num(c.clone()));
                for (j, (b, x)) in facs.iter().enumerate() {
                    if j != idx {
                        v.push(crate::ops::mul_factor(b, x));
                    }
                }
                Some(mul(&v))
            }
            _ => None,
        }
    }
    let r = match arg.kind() {
        Kind::Add(c, terms) => {
            if !c.is_zero() {
                return None;
            }
            let mut v = Vec::with_capacity(terms.len());
            for (t, k) in terms {
                v.push(mul2(&term(t)?, &num(k)));
            }
            add(&v)
        }
        _ => term(arg)?,
    };
    if crate::functions::contains_i(&r) {
        return None;
    }
    Some(r)
}

/// `as_coefficient(pi*I)` of a product.
fn ipi_coeff(arg: &Expr) -> Option<Expr> {
    let Kind::Mul(c, facs) = arg.kind() else { return None };
    let i = facs.iter().position(|(b, x)| b.is_constant(Constant::I) && x.is_one())?;
    let p = facs.iter().position(|(b, x)| b.is_constant(Constant::Pi) && x.is_one())?;
    let mut v = vec![Expr::num(c.clone())];
    for (j, (b, x)) in facs.iter().enumerate() {
        if j != i && j != p {
            v.push(crate::ops::mul_factor(b, x));
        }
    }
    let r = mul(&v);
    if crate::functions::contains_i(&r) || crate::traverse::has(&r, &pi()) {
        return None;
    }
    Some(r)
}

/// Recast an exact binary fraction float to a rational (`_pi_coeff`).
fn recast_float(c: &Number) -> Number {
    let f = c.to_f64();
    let frac = f.abs() % 1.0;
    if frac != 0.0 {
        let p = -(frac.log2().round()) as i32;
        if (0..=60).contains(&p) {
            let m = 2f64.powi(p);
            let cm = f * m;
            let i = cm.trunc();
            if i == cm {
                return Number::rational_i64(i as i64, m as i64).unwrap_or_else(|| c.clone());
            }
        }
        c.clone()
    } else {
        Number::int(f.trunc() as i64)
    }
}

/// SymPy's `_pi_coeff`: `arg/pi` normalized when `arg` is a multiple of pi.
pub(crate) fn pi_coeff(arg: &Expr) -> Option<Expr> {
    if arg.is_constant(Constant::Pi) {
        return Some(Expr::one());
    }
    if arg.is_zero() {
        return Some(Expr::zero());
    }
    let Kind::Mul(c, facs) = arg.kind() else { return None };
    let idx = facs.iter().position(|(b, x)| b.is_constant(Constant::Pi) && x.is_one())?;
    let rest: Vec<Expr> =
        facs.iter().enumerate().filter(|(j, _)| *j != idx).map(|(_, (b, x))| crate::ops::mul_factor(b, x)).collect();
    let x = mul(&rest);
    let c = if c.is_float() { recast_float(c) } else { c.clone() };
    if c.is_float() {
        return Some(mul2(&num(&c), &x));
    }
    if is_int(&x) {
        let c2 = nmod(&c, &Number::int(2));
        if c2.is_one() {
            return Some(x);
        }
        if c2.is_zero() {
            return Some(if even(&x).is_some() { Expr::zero() } else { Expr::int(2) });
        }
        return Some(mul2(&num(&c2), &x));
    }
    Some(mul2(&num(&c), &x))
}

/// SymPy's `_peeloff_pi`: `arg = x + m*pi` with `m` a multiple of 1/2.
fn peeloff_pi(arg: &Expr) -> Option<(Expr, Expr)> {
    let Kind::Add(c, terms) = arg.kind() else { return None };
    let mut rest: Vec<Expr> = vec![num(c)];
    let mut coeffs: Vec<Expr> = Vec::new();
    for (t, k) in terms {
        let kk = if t.is_constant(Constant::Pi) {
            Some(num(k))
        } else {
            pi_factor_quotient(t).map(|q| mul2(&q, &num(k)))
        };
        match kk {
            Some(kk) if !kk.is_zero() && rational3(&kk) => coeffs.push(kk),
            _ => rest.push(mul2(t, &num(k))),
        }
    }
    if coeffs.is_empty() {
        return None;
    }
    let total = add(&coeffs);
    let m1 = match total.as_number() {
        Some(k) => num(&nmod(k, &q(1, 2))),
        None => md(&total, &Expr::rational(1, 2)),
    };
    let m2 = sub(&total, &m1);
    if m2.is_zero() {
        return None;
    }
    if is_int(&m2) || (is_int(&mul2(&Expr::int(2), &m2)) && even(&m2) == Some(false)) {
        rest.push(mul2(&m1, &pi()));
        return Some((add(&rest), m2));
    }
    None
}

fn rational3(e: &Expr) -> bool {
    match e.as_number() {
        Some(n) => n.is_rational(),
        None => int3(e) == Some(true) || yes(assumptions::is_rational(e)),
    }
}

/// `t.coeff(pi)` of a product with a plain `pi` factor.
fn pi_factor_quotient(t: &Expr) -> Option<Expr> {
    let Kind::Mul(c, facs) = t.kind() else { return None };
    let idx = facs.iter().position(|(b, x)| b.is_constant(Constant::Pi) && x.is_one())?;
    let mut v = vec![Expr::num(c.clone())];
    for (j, (b, x)) in facs.iter().enumerate() {
        if j != idx {
            v.push(crate::ops::mul_factor(b, x));
        }
    }
    Some(mul(&v))
}

/// SymPy's `_peeloff_ipi`: `arg = x + m*I*pi` with `m` a multiple of 1/2.
fn peeloff_ipi(arg: &Expr) -> Option<(Expr, Number)> {
    let Kind::Add(_, terms) = arg.kind() else { return None };
    let target = ipi();
    let k = terms.iter().find(|(t, k)| *t == target && k.is_rational()).map(|(_, k)| k.clone())?;
    let m1 = nmod(&k, &q(1, 2));
    let m2 = k.sub(&m1);
    if m2.is_zero() {
        return None;
    }
    Some((sub(arg, &mul2(&num(&m2), &target)), m2))
}

/// `match_real_imag`: `(re, im)` of `r + i*I` with both parts real.
fn match_real_imag(z: &Expr) -> Option<(Expr, Expr)> {
    let (r, i) = match z.kind() {
        Kind::Add(c, terms) => {
            let mut rs = vec![num(c)];
            let mut is = Vec::new();
            for (t, k) in terms {
                let term = mul2(t, &num(k));
                if crate::functions::contains_i(t) {
                    is.push(term);
                } else {
                    rs.push(term);
                }
            }
            (add(&rs), add(&is))
        }
        _ if crate::functions::contains_i(z) => (Expr::zero(), z.clone()),
        _ => (z.clone(), Expr::zero()),
    };
    if i.is_zero() {
        return yes(assumptions::is_real(&r)).then_some((r, Expr::zero()));
    }
    let ic = i_coeff(&i)?;
    if yes(assumptions::is_real(&ic)) && yes(assumptions::is_real(&r)) {
        return Some((r, ic));
    }
    None
}

/// The constant tables, parsed once per thread.
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
    static ASIN_TABLE: Table = build(&[
        ("sqrt(3)/2", "pi/3"), ("sqrt(2)/2", "pi/4"), ("sqrt(5/8 - sqrt(5)/8)", "pi/5"),
        ("sqrt(2)*sqrt(5 - sqrt(5))/4", "pi/5"), ("sqrt(sqrt(5)/8 + 5/8)", "2*pi/5"),
        ("sqrt(2)*sqrt(sqrt(5) + 5)/4", "2*pi/5"), ("1/2", "pi/6"), ("sqrt(2 - sqrt(2))/2", "pi/8"),
        ("sqrt(1/2 - sqrt(2)/4)", "pi/8"), ("sqrt(sqrt(2) + 2)/2", "3*pi/8"),
        ("sqrt(sqrt(2)/4 + 1/2)", "3*pi/8"), ("-1/4 + sqrt(5)/4", "pi/10"), ("1/4 - sqrt(5)/4", "-pi/10"),
        ("1/4 + sqrt(5)/4", "3*pi/10"), ("-sqrt(2)/4 + sqrt(6)/4", "pi/12"), ("-sqrt(6)/4 + sqrt(2)/4", "-pi/12"),
        ("sqrt(2)*(-1 + sqrt(3))/4", "pi/12"), ("sqrt(2)*(1 - sqrt(3))/4", "-pi/12"),
        ("sqrt(2)/4 + sqrt(6)/4", "5*pi/12"), ("sqrt(2)*(1 + sqrt(3))/4", "5*pi/12"),
    ]);
    static ATAN_TABLE: Table = build(&[
        ("sqrt(3)/3", "pi/6"), ("sqrt(3)", "pi/3"), ("-1 + sqrt(2)", "pi/8"), ("1 - sqrt(2)", "-pi/8"),
        ("1 + sqrt(2)", "3*pi/8"), ("sqrt(5 - 2*sqrt(5))", "pi/5"), ("sqrt(2*sqrt(5) + 5)", "2*pi/5"),
        ("sqrt(1 - 2*sqrt(5)/5)", "pi/10"), ("sqrt(2*sqrt(5)/5 + 1)", "3*pi/10"), ("2 - sqrt(3)", "pi/12"),
        ("-2 + sqrt(3)", "-pi/12"), ("sqrt(3) + 2", "5*pi/12"),
    ]);
    static ACOSH_TABLE: Table = build(&[
        ("I", "log(1 + sqrt(2)) + I*pi/2"), ("-I", "log(1 + sqrt(2)) - I*pi/2"), ("1/2", "pi/3"),
        ("-1/2", "2*pi/3"), ("sqrt(2)/2", "pi/4"), ("-sqrt(2)/2", "3*pi/4"), ("sqrt(3)/2", "pi/6"),
        ("-sqrt(3)/2", "5*pi/6"), ("sqrt(2)*(-1 + sqrt(3))/4", "5*pi/12"),
        ("sqrt(2)*(1 - sqrt(3))/4", "7*pi/12"), ("sqrt(sqrt(2) + 2)/2", "pi/8"),
        ("-sqrt(sqrt(2) + 2)/2", "7*pi/8"), ("sqrt(2 - sqrt(2))/2", "3*pi/8"),
        ("-sqrt(2 - sqrt(2))/2", "5*pi/8"), ("sqrt(2)*(1 + sqrt(3))/4", "pi/12"),
        ("sqrt(2)*(-sqrt(3) - 1)/4", "11*pi/12"), ("1/4 + sqrt(5)/4", "pi/5"), ("-sqrt(5)/4 - 1/4", "4*pi/5"),
    ]);
    static LOG_ATAN_TABLE: Table = build(&[
        ("sqrt(3)", "pi/3"), ("1", "pi/4"), ("sqrt(5 - 2*sqrt(5))", "pi/5"),
        ("sqrt(2)*sqrt(5 - sqrt(5))/(1 + sqrt(5))", "pi/5"), ("sqrt(2*sqrt(5) + 5)", "2*pi/5"),
        ("sqrt(2)*sqrt(sqrt(5) + 5)/(-1 + sqrt(5))", "2*pi/5"), ("sqrt(3)/3", "pi/6"), ("-1 + sqrt(2)", "pi/8"),
        ("sqrt(2 - sqrt(2))/sqrt(sqrt(2) + 2)", "pi/8"), ("1 + sqrt(2)", "3*pi/8"),
        ("sqrt(sqrt(2) + 2)/sqrt(2 - sqrt(2))", "3*pi/8"), ("sqrt(1 - 2*sqrt(5)/5)", "pi/10"),
        ("(-sqrt(2) + sqrt(10))/(2*sqrt(sqrt(5) + 5))", "pi/10"), ("sqrt(2*sqrt(5)/5 + 1)", "3*pi/10"),
        ("(sqrt(2) + sqrt(10))/(2*sqrt(5 - sqrt(5)))", "3*pi/10"), ("2 - sqrt(3)", "pi/12"),
        ("(-1 + sqrt(3))/(1 + sqrt(3))", "pi/12"), ("sqrt(3) + 2", "5*pi/12"),
        ("(1 + sqrt(3))/(-1 + sqrt(3))", "5*pi/12"),
    ]);
}

fn table(t: &'static std::thread::LocalKey<Table>, e: &Expr) -> Option<Expr> {
    t.with(|t| lookup(t, e))
}

// ----------------------------------------------------------------------
// sin, cos, tan, cot
// ----------------------------------------------------------------------

/// `_table2`: `q -> (a, b)` with `cos(p*pi/q)` from `p*pi/a` and `p*pi/b`.
fn table2(q: i64) -> Option<(i64, i64)> {
    Some(match q {
        12 => (3, 4),
        20 => (4, 5),
        30 => (5, 6),
        15 => (6, 10),
        24 => (6, 8),
        40 => (8, 10),
        60 => (20, 30),
        120 => (40, 60),
        _ => return None,
    })
}

fn pq(k: &Number) -> Option<(BigInt, i64)> {
    let (p, d) = k.numer_denom()?;
    Some((p, d.to_i64()?))
}

pub(crate) fn sin_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return None,
        _ => {}
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&sin(&neg(arg))));
    }
    if let Some(c) = i_coeff(arg) {
        return Some(mul2(&ii(), &sinh(&c)));
    }
    if let Some(pc) = pi_coeff(arg) {
        if is_int(&pc) {
            return Some(Expr::zero());
        }
        if is_int(&mul2(&Expr::int(2), &pc)) && even(&pc) == Some(false) {
            return Some(pow(&Expr::minus_one(), &sub(&pc, &Expr::rational(1, 2))));
        }
        let Some(k) = pc.as_number().filter(|k| k.is_rational()).cloned() else {
            let narg = mul2(&pc, &pi());
            return (narg != *arg).then(|| sin(&narg));
        };
        let x = nmod(&k, &Number::int(2));
        if ncmp(&x, &Number::ONE) == Some(Ordering::Greater) {
            return Some(neg(&sin(&kpi(&nmod(&x, &Number::ONE)))));
        }
        if ncmp(&x.mul(&Number::int(2)), &Number::ONE) == Some(Ordering::Greater) {
            return Some(sin(&kpi(&Number::ONE.sub(&x))));
        }
        let narg = kpi(&nmod(&k.add(&q(3, 2)), &Number::int(2)));
        let result = cos(&narg);
        if !is_head(&result, &Func::Cos) {
            return Some(result);
        }
        let own = kpi(&k);
        return (own != *arg).then(|| sin(&own));
    }
    if let Some((x, m)) = peeloff_pi(arg) {
        let m = mul2(&m, &pi());
        return Some(add2(&mul2(&sin(&m), &cos(&x)), &mul2(&cos(&m), &sin(&x))));
    }
    trig_of_inverse(Func::Sin, arg)
}

pub(crate) fn cos_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::one()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return None,
        _ => {}
    }
    if could_extract_minus_sign(arg) {
        return Some(cos(&neg(arg)));
    }
    if let Some(c) = i_coeff(arg) {
        return Some(cosh(&c));
    }
    if let Some(pc) = pi_coeff(arg) {
        if is_int(&pc) {
            return Some(pow(&Expr::minus_one(), &pc));
        }
        if is_int(&mul2(&Expr::int(2), &pc)) && even(&pc) == Some(false) {
            return Some(Expr::zero());
        }
        let Some(k) = pc.as_number().filter(|k| k.is_rational()).cloned() else {
            let narg = mul2(&pc, &pi());
            return (narg != *arg).then(|| cos(&narg));
        };
        return cos_rational(&k);
    }
    if let Some((x, m)) = peeloff_pi(arg) {
        let m = mul2(&m, &pi());
        return Some(sub(&mul2(&cos(&m), &cos(&x)), &mul2(&sin(&m), &sin(&x))));
    }
    trig_of_inverse(Func::Cos, arg)
}

fn cos_rational(k: &Number) -> Option<Expr> {
    let (p0, qd) = pq(k)?;
    let qb = BigInt::from(qd);
    let p = p0.mod_floor(&(&qb * 2));
    if p > qb {
        return Some(neg(&cos(&kpi(&k.sub(&Number::ONE)))));
    }
    if &p * 2 > qb {
        return Some(neg(&cos(&kpi(&Number::ONE.sub(k)))));
    }
    let pn = Number::from_bigint(p.clone());
    if let Some((a, b)) = table2(qd) {
        let a = kpi(&pn.div(&Number::int(a)).unwrap());
        let b = kpi(&pn.div(&Number::int(b)).unwrap());
        return Some(add2(
            &mul2(&cos(&a), &cos(&b)),
            &mul2(&cos(&sub(&half_pi(), &a)), &cos(&sub(&half_pi(), &b))),
        ));
    }
    if qd > 12 {
        return None;
    }
    let s5 = sqrt(&Expr::int(5));
    match (qd, p.to_i64()?) {
        (3, _) => return Some(chebyshev_t(p.to_i64()?, &Expr::rational(1, 2))),
        (5, _) => {
            let c = div(&add2(&s5, &Expr::one()), &Expr::int(4));
            return Some(chebyshev_t(p.to_i64()?, &c));
        }
        _ => {}
    }
    if qd % 2 == 0 {
        let nval = cos(&kpi(&k.mul(&Number::int(2))));
        let x = k.mul(&Number::int(2)).add(&Number::ONE).div(&Number::int(2)).unwrap();
        let t = x.abs().floor().to_i64()?;
        let sgn = if x.is_negative() { -t } else { t };
        let v = sqrt(&div(&add2(&Expr::one(), &nval), &Expr::int(2)));
        return Some(if sgn.rem_euclid(2) == 1 { neg(&v) } else { v });
    }
    None
}

/// `chebyshevt(n, x).expand()` for a small non-negative `n`.
pub(crate) fn chebyshev_t(n: i64, x: &Expr) -> Expr {
    let (mut a, mut b) = (Expr::one(), x.clone());
    if n == 0 {
        return a;
    }
    for _ in 1..n {
        let c = crate::expand::expand(&sub(&mul(&[Expr::int(2), x.clone(), b.clone()]), &a));
        a = b;
        b = c;
    }
    crate::expand::expand(&b)
}

const TABLE10: [&str; 4] = ["sqrt(1 - 2*sqrt(5)/5)", "sqrt(5 - 2*sqrt(5))", "sqrt(1 + 2*sqrt(5)/5)", "sqrt(5 + 2*sqrt(5))"];

fn table10(n: i64) -> Expr {
    crate::parser::parse_str(TABLE10[(n - 1) as usize]).unwrap()
}

pub(crate) fn tan_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return None,
        _ => {}
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&tan(&neg(arg))));
    }
    if let Some(c) = i_coeff(arg) {
        return Some(mul2(&ii(), &tanh(&c)));
    }
    if let Some(pc) = pi_coeff(arg) {
        if is_int(&pc) {
            return Some(Expr::zero());
        }
        let Some(k) = pc.as_number().filter(|k| k.is_rational()).cloned() else {
            let narg = mul2(&pc, &pi());
            return (narg != *arg).then(|| tan(&narg));
        };
        let (p0, qd) = pq(&k)?;
        let p = p0.mod_floor(&BigInt::from(qd));
        if qd == 5 || qd == 10 {
            let n = ((&p * BigInt::from(10)) / BigInt::from(qd)).to_i64()?;
            return Some(if n > 5 { neg(&table10(10 - n)) } else { table10(n) });
        }
        if qd % 2 == 0 {
            let narg = kpi(&k.mul(&Number::int(2)));
            let (c, s) = (cos(&narg), cos(&sub(&narg, &half_pi())));
            if !is_head(&c, &Func::Cos) && !is_head(&s, &Func::Cos) {
                if s.is_zero() {
                    return Some(Expr::zoo());
                }
                return Some(sub(&inv(&s), &div(&c, &s)));
            }
        }
        let pn = Number::from_bigint(p);
        if let Some((a, b)) = table2(qd) {
            let ta = tan(&kpi(&pn.div(&Number::int(a)).unwrap()));
            let tb = tan(&kpi(&pn.div(&Number::int(b)).unwrap()));
            return Some(div(&sub(&ta, &tb), &add2(&Expr::one(), &mul2(&ta, &tb))));
        }
        let narg = kpi(&nmod(&k.add(&q(1, 2)), &Number::ONE).sub(&q(1, 2)));
        let (c, s) = (cos(&narg), cos(&sub(&narg, &half_pi())));
        if !is_head(&c, &Func::Cos) && !is_head(&s, &Func::Cos) {
            if c.is_zero() {
                return Some(Expr::zoo());
            }
            return Some(div(&s, &c));
        }
        return (narg != *arg).then(|| tan(&narg));
    }
    if let Some((x, m)) = peeloff_pi(arg) {
        let tanm = tan(&mul2(&m, &pi()));
        return Some(if tanm.is_constant(Constant::ComplexInfinity) { neg(&cot(&x)) } else { tan(&x) });
    }
    trig_of_inverse(Func::Tan, arg)
}

pub(crate) fn cot_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::zoo()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return None,
        _ => {}
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&cot(&neg(arg))));
    }
    if let Some(c) = i_coeff(arg) {
        return Some(neg(&mul2(&ii(), &coth(&c))));
    }
    if let Some(pc) = pi_coeff(arg) {
        if is_int(&pc) {
            return Some(Expr::zoo());
        }
        let Some(k) = pc.as_number().filter(|k| k.is_rational()).cloned() else {
            let narg = mul2(&pc, &pi());
            return (narg != *arg).then(|| cot(&narg));
        };
        let (p0, qd) = pq(&k)?;
        if qd == 5 || qd == 10 {
            return Some(tan(&sub(&half_pi(), arg)));
        }
        if qd > 2 && qd % 2 == 0 {
            let narg = kpi(&k.mul(&Number::int(2)));
            let (c, s) = (cos(&narg), cos(&sub(&narg, &half_pi())));
            if !is_head(&c, &Func::Cos) && !is_head(&s, &Func::Cos) {
                return Some(add2(&inv(&s), &div(&c, &s)));
            }
        }
        let p = Number::from_bigint(p0.mod_floor(&BigInt::from(qd)));
        if let Some((a, b)) = table2(qd) {
            let ca = cot(&kpi(&p.div(&Number::int(a)).unwrap()));
            let cb = cot(&kpi(&p.div(&Number::int(b)).unwrap()));
            return Some(div(&add2(&Expr::one(), &mul2(&ca, &cb)), &sub(&cb, &ca)));
        }
        let narg = kpi(&nmod(&k.add(&q(1, 2)), &Number::ONE).sub(&q(1, 2)));
        let (c, s) = (cos(&narg), cos(&sub(&narg, &half_pi())));
        if !is_head(&c, &Func::Cos) && !is_head(&s, &Func::Cos) {
            if s.is_zero() {
                return Some(Expr::zoo());
            }
            return Some(div(&c, &s));
        }
        return (narg != *arg).then(|| cot(&narg));
    }
    if let Some((x, m)) = peeloff_pi(arg) {
        let cotm = cot(&mul2(&m, &pi()));
        return Some(if cotm.is_constant(Constant::ComplexInfinity) { cot(&x) } else { neg(&tan(&x)) });
    }
    trig_of_inverse(Func::Cot, arg)
}

/// `f(g(x))` for a trigonometric `f` and an inverse trigonometric `g`.
fn trig_of_inverse(f: Func, arg: &Expr) -> Option<Expr> {
    let one = Expr::one();
    let sq = |x: &Expr| x.powi(2);
    let s1mx2 = |x: &Expr| sqrt(&sub(&one, &sq(x)));
    let s1mix2 = |x: &Expr| sqrt(&sub(&one, &inv(&sq(x))));
    if let Kind::Func(Func::Atan2, a) = arg.kind() {
        let (y, x) = (&a[0], &a[1]);
        let r = sqrt(&add2(&sq(x), &sq(y)));
        return Some(match f {
            Func::Sin => div(y, &r),
            Func::Cos => div(x, &r),
            Func::Tan => div(y, x),
            _ => div(x, y),
        });
    }
    let (g, x) = match arg.kind() {
        Kind::Func(g, a) if a.len() == 1 => (g.clone(), &a[0]),
        _ => return None,
    };
    let v = match (&f, &g) {
        (Func::Sin, Func::Asin) | (Func::Cos, Func::Acos) | (Func::Tan, Func::Atan) | (Func::Cot, Func::Acot) => x.clone(),
        (Func::Sin, Func::Atan) => div(x, &sqrt(&add2(&one, &sq(x)))),
        (Func::Sin, Func::Acos) | (Func::Cos, Func::Asin) => s1mx2(x),
        (Func::Sin, Func::Acot) => inv(&mul2(&sqrt(&add2(&one, &inv(&sq(x)))), x)),
        (Func::Sin, Func::Special(SpecialFn::Acsc)) | (Func::Cos, Func::Special(SpecialFn::Asec)) => inv(x),
        (Func::Sin, Func::Special(SpecialFn::Asec)) | (Func::Cos, Func::Special(SpecialFn::Acsc)) => s1mix2(x),
        (Func::Cos, Func::Atan) => inv(&sqrt(&add2(&one, &sq(x)))),
        (Func::Cos, Func::Acot) => inv(&sqrt(&add2(&one, &inv(&sq(x))))),
        (Func::Tan, Func::Asin) | (Func::Cot, Func::Acos) => div(x, &s1mx2(x)),
        (Func::Tan, Func::Acos) | (Func::Cot, Func::Asin) => div(&s1mx2(x), x),
        (Func::Tan, Func::Acot) | (Func::Cot, Func::Atan) => inv(x),
        (Func::Tan, Func::Special(SpecialFn::Acsc)) | (Func::Cot, Func::Special(SpecialFn::Asec)) => {
            inv(&mul2(&s1mix2(x), x))
        }
        (Func::Tan, Func::Special(SpecialFn::Asec)) | (Func::Cot, Func::Special(SpecialFn::Acsc)) => mul2(&s1mix2(x), x),
        _ => return None,
    };
    Some(v)
}

// ----------------------------------------------------------------------
// sec, csc (ReciprocalTrigonometricFunction.eval)
// ----------------------------------------------------------------------

pub(crate) fn sec_eval(arg: &Expr) -> Option<Expr> {
    reciprocal_trig(Func::Sec, Func::Cos, true, SpecialFn::Asec, arg)
}

pub(crate) fn csc_eval(arg: &Expr) -> Option<Expr> {
    reciprocal_trig(Func::Csc, Func::Sin, false, SpecialFn::Acsc, arg)
}

fn reciprocal_trig(f: Func, base: Func, is_even: bool, inverse: SpecialFn, arg: &Expr) -> Option<Expr> {
    if could_extract_minus_sign(arg) {
        let m = call(&f, &neg(arg));
        return Some(if is_even { m } else { neg(&m) });
    }
    if let Some(pc) = pi_coeff(arg) {
        if let Some(k) = pc.as_number().filter(|k| k.is_rational() && !k.mul(&Number::int(2)).is_integer()) {
            let (p0, qd) = pq(k)?;
            let qb = BigInt::from(qd);
            let p = p0.mod_floor(&(&qb * 2));
            if p > qb {
                return Some(neg(&call(&f, &kpi(&k.sub(&Number::ONE)))));
            }
            if &p * 2 > qb {
                let r = call(&f, &kpi(&Number::ONE.sub(k)));
                return Some(if is_even { neg(&r) } else { r });
            }
        }
    }
    if let Some(x) = special_arg(arg, inverse) {
        return Some(x.clone());
    }
    let t = eval_of(&base, arg)?;
    let nt = neg(&t);
    for (cand, sign) in [(&t, false), (&nt, true)] {
        let r = if let Some(y) = arg_of(cand, &Func::Cos) {
            Some(call(&Func::Sec, y))
        } else {
            arg_of(cand, &Func::Sin).map(|y| call(&Func::Csc, y))
        };
        if let Some(r) = r {
            return Some(if sign { neg(&r) } else { r });
        }
    }
    Some(inv(&t))
}

// ----------------------------------------------------------------------
// asin, acos, atan, acot, atan2
// ----------------------------------------------------------------------

pub(crate) fn asin_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity) => return Some(mul2(&Expr::neg_oo(), &ii())),
        Kind::Const(Constant::NegInfinity) => return Some(mul2(&Expr::oo(), &ii())),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::zoo()),
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        Kind::Num(n) if n.is_one() => return Some(half_pi()),
        Kind::Num(n) if n.is_minus_one() => return Some(neg(&half_pi())),
        _ => {}
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&asin(&neg(arg))));
    }
    if is_num(arg) {
        if let Some(v) = table(&ASIN_TABLE, arg) {
            return Some(v);
        }
    }
    if let Some(c) = i_coeff(arg) {
        return Some(mul2(&ii(), &asinh(&c)));
    }
    if let Some(ang) = arg_of(arg, &Func::Sin) {
        if comparable(ang) {
            let mut ang = md(ang, &two_pi());
            if gt(&ang, &pi()) {
                ang = sub(&pi(), &ang);
            }
            if gt(&ang, &half_pi()) {
                ang = sub(&pi(), &ang);
            }
            if lt(&ang, &neg(&half_pi())) {
                ang = sub(&neg(&pi()), &ang);
            }
            return Some(ang);
        }
    }
    if let Some(ang) = arg_of(arg, &Func::Cos) {
        if comparable(ang) {
            return Some(sub(&half_pi(), &acos(arg)));
        }
    }
    None
}

pub(crate) fn acos_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity) => return Some(mul2(&Expr::oo(), &ii())),
        Kind::Const(Constant::NegInfinity) => return Some(mul2(&Expr::neg_oo(), &ii())),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::zoo()),
        Kind::Num(n) if n.is_zero() => return Some(half_pi()),
        Kind::Num(n) if n.is_one() => return Some(Expr::zero()),
        Kind::Num(n) if n.is_minus_one() => return Some(pi()),
        _ => {}
    }
    if is_num(arg) {
        if let Some(v) = table(&ASIN_TABLE, arg) {
            return Some(sub(&half_pi(), &v));
        }
        if let Some(v) = table(&ASIN_TABLE, &neg(arg)) {
            return Some(add2(&half_pi(), &v));
        }
    }
    if i_coeff(arg).is_some() {
        return Some(sub(&half_pi(), &asin(arg)));
    }
    let (narg, minus) = match arg.kind() {
        Kind::Mul(c, facs) if c.is_minus_one() && facs.len() == 1 => {
            (crate::ops::mul_factor(&facs[0].0, &facs[0].1), true)
        }
        _ => (arg.clone(), false),
    };
    if let Some(ang) = arg_of(&narg, &Func::Cos) {
        if comparable(ang) {
            let mut ang = if minus { sub(&pi(), ang) } else { ang.clone() };
            ang = md(&ang, &two_pi());
            if gt(&ang, &pi()) {
                ang = sub(&two_pi(), &ang);
            }
            return Some(ang);
        }
    }
    if let Some(ang) = arg_of(&narg, &Func::Sin) {
        if comparable(ang) {
            let a = asin(&narg);
            return Some(if minus { add2(&half_pi(), &a) } else { sub(&half_pi(), &a) });
        }
    }
    None
}

pub(crate) fn atan_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity) => return Some(half_pi()),
        Kind::Const(Constant::NegInfinity) => return Some(neg(&half_pi())),
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        Kind::Num(n) if n.is_one() => return Some(div(&pi(), &Expr::int(4))),
        Kind::Num(n) if n.is_minus_one() => return Some(div(&pi(), &Expr::int(-4))),
        _ => {}
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&atan(&neg(arg))));
    }
    if is_num(arg) {
        if let Some(v) = table(&ATAN_TABLE, arg) {
            return Some(v);
        }
    }
    if let Some(c) = i_coeff(arg) {
        return Some(mul2(&ii(), &atanh(&c)));
    }
    if let Some(ang) = arg_of(arg, &Func::Tan) {
        if comparable(ang) {
            let mut ang = md(ang, &pi());
            if gt(&ang, &half_pi()) {
                ang = sub(&ang, &pi());
            }
            return Some(ang);
        }
    }
    if let Some(ang) = arg_of(arg, &Func::Cot) {
        if comparable(ang) {
            let mut ang = sub(&half_pi(), &acot(arg));
            if gt(&ang, &half_pi()) {
                ang = sub(&ang, &pi());
            }
            return Some(ang);
        }
    }
    None
}

pub(crate) fn acot_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity) => {
            return Some(Expr::zero())
        }
        Kind::Num(n) if n.is_zero() => return Some(half_pi()),
        Kind::Num(n) if n.is_one() => return Some(div(&pi(), &Expr::int(4))),
        Kind::Num(n) if n.is_minus_one() => return Some(div(&pi(), &Expr::int(-4))),
        _ => {}
    }
    if could_extract_minus_sign(arg) {
        return Some(neg(&acot(&neg(arg))));
    }
    if is_num(arg) {
        if let Some(v) = table(&ATAN_TABLE, arg) {
            let mut ang = sub(&half_pi(), &v);
            if gt(&ang, &half_pi()) {
                ang = sub(&ang, &pi());
            }
            return Some(ang);
        }
    }
    if let Some(c) = i_coeff(arg) {
        return Some(neg(&mul2(&ii(), &super::acoth(&c))));
    }
    if let Some(ang) = arg_of(arg, &Func::Cot) {
        if comparable(ang) {
            let mut ang = md(ang, &pi());
            if gt(&ang, &half_pi()) {
                ang = sub(&ang, &pi());
            }
            return Some(ang);
        }
    }
    if let Some(ang) = arg_of(arg, &Func::Tan) {
        if comparable(ang) {
            let mut ang = sub(&half_pi(), &atan(arg));
            if gt(&ang, &half_pi()) {
                ang = sub(&ang, &pi());
            }
            return Some(ang);
        }
    }
    None
}

fn ext_real(e: &Expr) -> bool {
    e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity) || yes(assumptions::is_real(e))
}

/// A purely imaginary number (`is_imaginary` for numbers).
fn imaginary_number(e: &Expr) -> bool {
    is_num(e) && re(e).is_zero() && {
        let i = im(e);
        i.is_number() && !i.is_zero() || yes(assumptions::is_nonzero(&i))
    }
}

pub(crate) fn atan2_eval(y0: &Expr, x0: &Expr) -> Option<Expr> {
    if x0.is_constant(Constant::NegInfinity) {
        if yes(assumptions::is_zero(y0)) {
            return Some(pi());
        }
        return Some(sub(&mul(&[Expr::int(2), pi(), heaviside(&re(y0))]), &pi()));
    }
    if x0.is_constant(Constant::Infinity) {
        return Some(Expr::zero());
    }
    let (mut y, mut x) = (y0.clone(), x0.clone());
    if imaginary_number(&x) && imaginary_number(&y) {
        x = im(&x);
        y = im(&y);
    }
    if ext_real(&x) && ext_real(&y) {
        if yes(assumptions::is_positive(&x)) {
            return Some(atan(&div(&y, &x)));
        } else if yes(assumptions::is_negative(&x)) {
            if yes(assumptions::is_negative(&y)) {
                return Some(sub(&atan(&div(&y, &x)), &pi()));
            } else if yes(assumptions::is_nonnegative(&y)) {
                return Some(add2(&atan(&div(&y, &x)), &pi()));
            }
        } else if yes(assumptions::is_zero(&x)) {
            if yes(assumptions::is_positive(&y)) {
                return Some(half_pi());
            } else if yes(assumptions::is_negative(&y)) {
                return Some(neg(&half_pi()));
            } else if yes(assumptions::is_zero(&y)) {
                return Some(Expr::nan());
            }
        }
    }
    if yes(assumptions::is_zero(&y)) {
        if ext_real(&x) && yes(assumptions::is_nonzero(&x)) {
            return Some(mul2(&pi(), &sub(&Expr::one(), &heaviside(&x))));
        }
        if is_num(&x) {
            use crate::relational as rel;
            return Some(rel::piecewise(vec![
                (pi(), rel::lt(&re(&x), &Expr::zero())),
                (Expr::zero(), rel::ne(&x, &Expr::zero())),
                (Expr::nan(), Expr::boolean(true)),
            ]));
        }
    }
    if is_num(&x) && is_num(&y) {
        let z = div(&add2(&x, &mul2(&ii(), &y)), &sqrt(&add2(&x.powi(2), &y.powi(2))));
        return Some(neg(&mul2(&ii(), &log(&z))));
    }
    None
}

// ----------------------------------------------------------------------
// sinh, cosh, tanh, coth
// ----------------------------------------------------------------------

fn neg_number(arg: &Expr) -> bool {
    arg.as_number().is_some_and(|n| n.is_negative())
}

pub(crate) fn sinh_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity) => return Some(Expr::oo()),
        Kind::Const(Constant::NegInfinity) => return Some(Expr::neg_oo()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        _ => {}
    }
    if neg_number(arg) {
        return Some(neg(&sinh(&neg(arg))));
    }
    if !arg.is_number() {
        if let Some(c) = i_coeff(arg) {
            return Some(mul2(&ii(), &sin(&c)));
        }
        if could_extract_minus_sign(arg) {
            return Some(neg(&sinh(&neg(arg))));
        }
        if let Some((x, m)) = peeloff_ipi(arg) {
            let m = mul(&[num(&m), pi(), ii()]);
            return Some(add2(&mul2(&sinh(&m), &cosh(&x)), &mul2(&cosh(&m), &sinh(&x))));
        }
    }
    hyp_of_inverse(Func::Sinh, arg)
}

pub(crate) fn cosh_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return Some(Expr::oo()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Num(n) if n.is_zero() => return Some(Expr::one()),
        _ => {}
    }
    if neg_number(arg) {
        return Some(cosh(&neg(arg)));
    }
    if !arg.is_number() {
        if let Some(c) = i_coeff(arg) {
            return Some(cos(&c));
        }
        if could_extract_minus_sign(arg) {
            return Some(cosh(&neg(arg)));
        }
        if let Some((x, m)) = peeloff_ipi(arg) {
            let m = mul(&[num(&m), pi(), ii()]);
            return Some(add2(&mul2(&cosh(&m), &cosh(&x)), &mul2(&sinh(&m), &sinh(&x))));
        }
    }
    hyp_of_inverse(Func::Cosh, arg)
}

pub(crate) fn tanh_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity) => return Some(Expr::one()),
        Kind::Const(Constant::NegInfinity) => return Some(Expr::minus_one()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        _ => {}
    }
    if neg_number(arg) {
        return Some(neg(&tanh(&neg(arg))));
    }
    if !arg.is_number() {
        if let Some(c) = i_coeff(arg) {
            if could_extract_minus_sign(&c) {
                return Some(neg(&mul2(&ii(), &tan(&neg(&c)))));
            }
            return Some(mul2(&ii(), &tan(&c)));
        }
        if could_extract_minus_sign(arg) {
            return Some(neg(&tanh(&neg(arg))));
        }
        if let Some((x, m)) = peeloff_ipi(arg) {
            // tanh(m*pi*I) is zoo for half-integers and zero for integers
            return Some(if m.is_integer() { tanh(&x) } else { coth(&x) });
        }
    }
    hyp_of_inverse(Func::Tanh, arg)
}

pub(crate) fn coth_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity) => return Some(Expr::one()),
        Kind::Const(Constant::NegInfinity) => return Some(Expr::minus_one()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Num(n) if n.is_zero() => return Some(Expr::zoo()),
        _ => {}
    }
    if neg_number(arg) {
        return Some(neg(&coth(&neg(arg))));
    }
    if !arg.is_number() {
        if let Some(c) = i_coeff(arg) {
            if could_extract_minus_sign(&c) {
                return Some(mul2(&ii(), &cot(&neg(&c))));
            }
            return Some(neg(&mul2(&ii(), &cot(&c))));
        }
        if could_extract_minus_sign(arg) {
            return Some(neg(&coth(&neg(arg))));
        }
        if let Some((x, m)) = peeloff_ipi(arg) {
            // coth(m*pi*I) is zoo for integers and zero for half-integers
            return Some(if m.is_integer() { coth(&x) } else { tanh(&x) });
        }
    }
    hyp_of_inverse(Func::Coth, arg)
}

fn hyp_of_inverse(f: Func, arg: &Expr) -> Option<Expr> {
    let (g, x) = match arg.kind() {
        Kind::Func(g, a) if a.len() == 1 => (g.clone(), &a[0]),
        _ => return None,
    };
    let one = Expr::one();
    let sq = |x: &Expr| x.powi(2);
    let pm = |x: &Expr| mul2(&sqrt(&sub(x, &one)), &sqrt(&add2(x, &one)));
    let acoth = Func::Special(SpecialFn::Acoth);
    let v = match (&f, &g) {
        (Func::Sinh, Func::Asinh) | (Func::Cosh, Func::Acosh) | (Func::Tanh, Func::Atanh) => x.clone(),
        (Func::Sinh, Func::Acosh) => pm(x),
        (Func::Sinh, Func::Atanh) => div(x, &sqrt(&sub(&one, &sq(x)))),
        (Func::Cosh, Func::Asinh) => sqrt(&add2(&one, &sq(x))),
        (Func::Cosh, Func::Atanh) => inv(&sqrt(&sub(&one, &sq(x)))),
        (Func::Tanh, Func::Asinh) => div(x, &sqrt(&add2(&one, &sq(x)))),
        (Func::Tanh, Func::Acosh) => div(&pm(x), x),
        (Func::Coth, Func::Asinh) => div(&sqrt(&add2(&one, &sq(x))), x),
        (Func::Coth, Func::Acosh) => div(x, &pm(x)),
        (Func::Cosh, _) if g == acoth => div(x, &pm(x)),
        (Func::Coth, Func::Atanh) => inv(x),
        (Func::Sinh, _) if g == acoth => inv(&pm(x)),
        (Func::Tanh, _) if g == acoth => inv(x),
        (Func::Coth, _) if g == acoth => x.clone(),
        _ => return None,
    };
    Some(v)
}

// ----------------------------------------------------------------------
// asinh, acosh, atanh
// ----------------------------------------------------------------------

pub(crate) fn asinh_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity) => return Some(Expr::oo()),
        Kind::Const(Constant::NegInfinity) => return Some(Expr::neg_oo()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::zoo()),
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        Kind::Num(n) if n.is_one() => return Some(log(&add2(&sqrt(&Expr::int(2)), &Expr::one()))),
        Kind::Num(n) if n.is_minus_one() => return Some(log(&sub(&sqrt(&Expr::int(2)), &Expr::one()))),
        _ => {}
    }
    if neg_number(arg) {
        return Some(neg(&asinh(&neg(arg))));
    }
    if !arg.is_number() {
        if let Some(c) = i_coeff(arg) {
            return Some(mul2(&ii(), &asin(&c)));
        }
        if could_extract_minus_sign(arg) {
            return Some(neg(&asinh(&neg(arg))));
        }
    }
    if let Some(z) = arg_of(arg, &Func::Sinh).filter(|z| is_num(z)) {
        if yes(assumptions::is_real(z)) {
            return Some(z.clone());
        }
        if let Some((_, i)) = match_real_imag(z) {
            let f = floor(&div(&add2(&i, &half_pi()), &pi()));
            let m = sub(z, &mul(&[ii(), pi(), f.clone()]));
            match even(&f) {
                Some(true) => return Some(m),
                Some(false) => return Some(neg(&m)),
                None => {}
            }
        }
    }
    None
}

pub(crate) fn acosh_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return Some(Expr::oo()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::zoo()),
        Kind::Num(n) if n.is_zero() => return Some(mul2(&ii(), &half_pi())),
        Kind::Num(n) if n.is_one() => return Some(Expr::zero()),
        Kind::Num(n) if n.is_minus_one() => return Some(mul2(&ii(), &pi())),
        _ => {}
    }
    if is_num(arg) {
        if let Some(v) = table(&ACOSH_TABLE, arg) {
            return Some(if yes(assumptions::is_real(arg)) { mul2(&v, &ii()) } else { v });
        }
    }
    let ioo = mul2(&ii(), &Expr::oo());
    if *arg == ioo {
        return Some(add2(&Expr::oo(), &mul2(&ii(), &half_pi())));
    }
    if *arg == neg(&ioo) {
        return Some(sub(&Expr::oo(), &mul2(&ii(), &half_pi())));
    }
    if let Some(z) = arg_of(arg, &Func::Cosh).filter(|z| is_num(z)) {
        if yes(assumptions::is_real(z)) {
            return Some(abs(z));
        }
        if let Some((r, i)) = match_real_imag(z) {
            let f = floor(&div(&i, &pi()));
            let mut m = sub(z, &mul(&[ii(), pi(), f.clone()]));
            match even(&f) {
                Some(true) => {
                    if yes(assumptions::is_nonnegative(&r)) {
                        return Some(m);
                    } else if yes(assumptions::is_negative(&r)) {
                        return Some(neg(&m));
                    }
                }
                Some(false) => {
                    m = sub(&m, &ipi());
                    if yes(assumptions::is_nonpositive(&r)) {
                        return Some(neg(&m));
                    } else if yes(assumptions::is_positive(&r)) {
                        return Some(m);
                    }
                }
                None => {}
            }
        }
    }
    None
}

pub(crate) fn atanh_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::zero()),
        Kind::Num(n) if n.is_one() => return Some(Expr::oo()),
        Kind::Num(n) if n.is_minus_one() => return Some(Expr::neg_oo()),
        Kind::Const(Constant::Infinity) => return Some(neg(&mul2(&ii(), &atan(arg)))),
        Kind::Const(Constant::NegInfinity) => return Some(mul2(&ii(), &atan(&neg(arg)))),
        _ => {}
    }
    if neg_number(arg) {
        return Some(neg(&atanh(&neg(arg))));
    }
    if !arg.is_number() {
        if let Some(c) = i_coeff(arg) {
            return Some(mul2(&ii(), &atan(&c)));
        }
        if could_extract_minus_sign(arg) {
            return Some(neg(&atanh(&neg(arg))));
        }
    }
    if let Some(z) = arg_of(arg, &Func::Tanh).filter(|z| is_num(z)) {
        if yes(assumptions::is_real(z)) {
            return Some(z.clone());
        }
        if let Some((_, i)) = match_real_imag(z) {
            let f = floor(&div(&mul2(&Expr::int(2), &i), &pi()));
            let m = sub(z, &div(&mul(&[ii(), f.clone(), pi()]), &Expr::int(2)));
            match even(&f) {
                Some(true) => return Some(m),
                Some(false) => return Some(sub(&m, &mul2(&ii(), &half_pi()))),
                None => {}
            }
        }
    }
    None
}

// ----------------------------------------------------------------------
// exp, log
// ----------------------------------------------------------------------

pub(crate) fn exp_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::one()),
        Kind::Num(n) if n.is_one() => return Some(Expr::e()),
        Kind::Const(Constant::Infinity) => return Some(Expr::oo()),
        Kind::Const(Constant::NegInfinity) => return Some(Expr::zero()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Func(Func::Log, x) if x.len() == 1 => return Some(x[0].clone()),
        Kind::Mul(c, facs) => {
            if let Some(coeff) = ipi_coeff(arg) {
                if is_int(&mul2(&Expr::int(2), &coeff)) {
                    if even(&coeff) == Some(true) {
                        return Some(Expr::one());
                    } else if odd(&coeff) == Some(true) {
                        return Some(Expr::minus_one());
                    }
                    let h = add2(&coeff, &Expr::rational(1, 2));
                    if even(&h) == Some(true) {
                        return Some(neg(&ii()));
                    } else if odd(&h) == Some(true) {
                        return Some(ii());
                    }
                } else if let Some(k) = coeff.as_number().filter(|k| k.is_rational()) {
                    let mut nc = nmod(k, &Number::int(2));
                    if ncmp(&nc, &Number::ONE) == Some(Ordering::Greater) {
                        nc = nc.sub(&Number::int(2));
                    }
                    if nc != *k {
                        return Some(exp(&mul(&[num(&nc), pi(), ii()])));
                    }
                }
            }
            // a single log factor: exp(c*log(x)) = x**c
            if c.is_float() && facs.iter().any(|(b, _)| b.is_constant(Constant::Infinity)) {
                return None;
            }
            let mut coeffs = vec![num(c)];
            let mut log_term: Option<Expr> = None;
            for (b, x) in facs {
                if b.is_constant(Constant::Infinity) || b.is_constant(Constant::NegInfinity) {
                    return None;
                }
                let term = crate::ops::mul_factor(b, x);
                if let Some(l) = arg_of(&term, &Func::Log) {
                    if log_term.is_some() {
                        return None;
                    }
                    log_term = Some(l.clone());
                } else if comparable(&term) {
                    coeffs.push(term);
                } else {
                    return None;
                }
            }
            return log_term.map(|l| pow(&l, &mul(&coeffs)));
        }
        Kind::Add(..) => {
            let mut out: Vec<Expr> = Vec::new();
            let mut acc: Vec<Expr> = Vec::new();
            let mut changed = false;
            for a in arg.args() {
                if a.is_one() {
                    acc.push(a);
                    continue;
                }
                let newa = exp(&a);
                match arg_of(&newa, &Func::Exp) {
                    Some(b) => {
                        if *b != a {
                            acc.push(b.clone());
                            changed = true;
                        } else {
                            acc.push(a);
                        }
                    }
                    None => out.push(newa),
                }
            }
            if !out.is_empty() || changed {
                let rest = add(&acc);
                let e = if rest.is_add() { Expr::from_kind(Kind::Func(Func::Exp, vec![rest])) } else { exp(&rest) };
                out.push(e);
                return Some(mul(&out));
            }
        }
        _ => {}
    }
    None
}

pub(crate) fn log_eval(arg: &Expr) -> Option<Expr> {
    match arg.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::zoo()),
        Kind::Num(n) if n.is_one() => return Some(Expr::zero()),
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return Some(Expr::oo()),
        Kind::Num(n) if n.is_rational() && !n.is_integer() && n.numer_denom().is_some_and(|(p, _)| p.is_one()) => {
            let (_, d) = n.numer_denom()?;
            return Some(neg(&log(&Expr::num(Number::from_bigint(d)))));
        }
        _ => {}
    }
    let exp_arg = match arg.kind() {
        Kind::Func(Func::Exp, x) if x.len() == 1 => Some(&x[0]),
        Kind::Pow(b, x) if b.is_constant(Constant::E) => Some(x),
        _ => None,
    };
    if let Some(x) = exp_arg {
        if ext_real(x) {
            return Some(x.clone());
        }
        if is_num(x) && is_head(arg, &Func::Exp) {
            if let Some((r, i)) = match_real_imag(x) {
                if !i.is_zero() && comparable(&i) {
                    let mut i = md(&i, &two_pi());
                    if gt(&i, &pi()) {
                        i = sub(&i, &two_pi());
                    }
                    return Some(add2(&r, &crate::simp::expand_mul_deep(&mul2(&i, &ii()), false)));
                }
            }
            return None;
        }
    }
    if is_num(arg) {
        if yes(assumptions::is_negative(arg)) {
            return Some(add2(&ipi(), &log(&neg(arg))));
        }
    }
    match arg.kind() {
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::zoo()),
        Kind::Const(Constant::E) => return Some(Expr::one()),
        _ => {}
    }
    if !arg.is_add() {
        if let Some(c) = i_coeff(arg) {
            if c.is_constant(Constant::Infinity) || c.is_constant(Constant::NegInfinity) {
                return Some(Expr::oo());
            }
            if let Some(k) = c.as_number().filter(|k| k.is_rational()) {
                let h = mul2(&ii(), &half_pi());
                return Some(if !k.is_negative() {
                    add2(&h, &log(&c))
                } else {
                    sub(&log(&neg(&c)), &h)
                });
            }
        }
    }
    if is_num(arg) && crate::functions::contains_i(arg) && yes(assumptions::fact(arg, "algebraic")) {
        return log_algebraic(arg);
    }
    None
}

/// `log(coeff*(r + i*I))` with a tabulated `atan(i/r)`.
fn log_algebraic(arg: &Expr) -> Option<Expr> {
    let (mut coeff, mut a) = match arg.kind() {
        Kind::Mul(c, facs) => {
            let (dep, ind): (Vec<_>, Vec<_>) =
                facs.iter().map(|(b, x)| crate::ops::mul_factor(b, x)).partition(crate::functions::contains_i);
            let mut ind = ind;
            ind.push(num(c));
            (mul(&ind), mul(&dep))
        }
        _ => (Expr::one(), arg.clone()),
    };
    if yes(assumptions::is_negative(&coeff)) {
        coeff = neg(&coeff);
        a = neg(&a);
    }
    let a = crate::simp::expand_mul_deep(&a, false);
    let (r, i) = match a.kind() {
        Kind::Add(c, terms) => {
            let mut rs = vec![num(c)];
            let mut is = Vec::new();
            for (t, k) in terms {
                let term = mul2(t, &num(k));
                if crate::functions::contains_i(t) {
                    is.push(term);
                } else {
                    rs.push(term);
                }
            }
            (add(&rs), add(&is))
        }
        _ => (Expr::zero(), a.clone()),
    };
    let i = i_coeff(&i)?;
    if !(yes(assumptions::is_real(&coeff)) && !i.is_zero() && yes(assumptions::is_real(&i)) && yes(assumptions::is_real(&r))) {
        return None;
    }
    if r.is_zero() {
        let h = mul2(&ii(), &half_pi());
        if yes(assumptions::is_positive(&i)) {
            return Some(add2(&h, &log(&mul2(&coeff, &i))));
        } else if yes(assumptions::is_negative(&i)) {
            return Some(sub(&log(&mul2(&coeff, &neg(&i))), &h));
        }
        return None;
    }
    let t = div(&i, &r);
    let t1 = neg(&t);
    let candidates = |t: &Expr| -> Option<Expr> {
        table(&LOG_ATAN_TABLE, t).or_else(|| {
            let c = crate::matrix::ratfunc::cancel(t);
            if c != *t { table(&LOG_ATAN_TABLE, &c) } else { None }
        })
    };
    let modulus = || {
        let m = abs(&a);
        if is_head(&m, &Func::Abs) {
            let sq = crate::simp::expand_mul_deep(&add2(&r.powi(2), &i.powi(2)), true);
            mul2(&coeff, &sqrt(&crate::expand::expand(&sq)))
        } else {
            mul2(&coeff, &m)
        }
    };
    if let Some(v) = candidates(&t) {
        let im_part = if yes(assumptions::is_positive(&r)) { v } else { sub(&v, &pi()) };
        return Some(add2(&log(&modulus()), &mul2(&ii(), &im_part)));
    }
    if let Some(v) = candidates(&t1) {
        let im_part = if yes(assumptions::is_positive(&r)) { neg(&v) } else { sub(&pi(), &v) };
        return Some(add2(&log(&modulus()), &mul2(&ii(), &im_part)));
    }
    None
}

// ----------------------------------------------------------------------
// log(x, b)
// ----------------------------------------------------------------------

fn trailing(n: &BigInt) -> i64 {
    n.trailing_zeros().map(|t| t as i64).unwrap_or(0)
}

/// SymPy's `multiplicity(p, n)` for integers and rationals.
fn multiplicity_int(p: &BigInt, n: &BigInt) -> Option<i64> {
    if n.is_zero() {
        return None;
    }
    let two = BigInt::from(2);
    if *p == two {
        return Some(trailing(n));
    }
    if *p < two {
        return None;
    }
    if p == n {
        return Some(1);
    }
    let mut m = 0i64;
    let mut n = n.clone();
    loop {
        let (d, r) = n.div_rem(p);
        if !r.is_zero() {
            return Some(m);
        }
        m += 1;
        n = d;
    }
}

fn multiplicity(p: &Number, n: &Number) -> Option<i64> {
    if p.is_float() || n.is_float() {
        return None;
    }
    let (pp, pq) = p.numer_denom()?;
    let (np, nq) = n.numer_denom()?;
    if pq.is_one() && nq.is_one() {
        return multiplicity_int(&pp, &np);
    }
    if pq.is_one() {
        if np.is_one() {
            return Some(-multiplicity_int(&pp, &nq)?);
        }
        return Some(multiplicity_int(&pp, &np)? - multiplicity_int(&pp, &nq)?);
    }
    if pp.is_one() {
        return multiplicity_int(&pq, &nq);
    }
    let like = multiplicity_int(&pp, &np)?.min(multiplicity_int(&pq, &nq)?);
    let cross = multiplicity_int(&pq, &np)?.min(multiplicity_int(&pp, &nq)?);
    Some(like - cross)
}

/// `log(arg, base)` (`log.eval` with a base).
pub(crate) fn log_base(arg: &Expr, base: &Expr) -> Expr {
    if base.is_one() {
        return if arg.is_one() { Expr::nan() } else { Expr::zoo() };
    }
    if let (Some(b), Some(a)) = (base.as_number(), arg.as_number()) {
        if let Some(n) = multiplicity(b, a) {
            if n != 0 {
                let n = Expr::int(n);
                return add2(&n, &div(&log(&div(arg, &pow(base, &n))), &log(base)));
            }
            return div(&log(arg), &log(base));
        }
    }
    if base.is_constant(Constant::E) {
        return log(arg);
    }
    let (la, lb) = (log(arg), log(base));
    // zoo absorbs a finite nonzero real factor (Mul.flatten)
    if la.is_constant(Constant::ComplexInfinity) && is_num(&lb) && yes(assumptions::is_real(&lb)) && yes(assumptions::is_nonzero(&lb)) {
        return Expr::zoo();
    }
    div(&la, &lb)
}

/// Exact evaluation SymPy performs before the float evaluation of a call:
/// `sin(1.5*pi) = -1` (binary fraction multiples of pi) and
/// `log(-0.5) = log(0.5) + I*pi`.
pub(crate) fn float_exact(f: &Func, arg: &Expr) -> Option<Expr> {
    match f {
        Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc => {
            let Kind::Mul(c, facs) = arg.kind() else { return None };
            if !c.is_float() || facs.len() != 1 || !facs[0].0.is_constant(Constant::Pi) || !facs[0].1.is_one() {
                return None;
            }
            let k = recast_float(c);
            if k.is_float() {
                return None;
            }
            Some(call(f, &kpi(&k)))
        }
        Func::Log => {
            let n = arg.as_number().filter(|n| n.is_float() && n.is_negative())?;
            Some(add2(&log(&num(&n.neg())), &ipi()))
        }
        _ => None,
    }
}
