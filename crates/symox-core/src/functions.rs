//! Elementary and special functions with SymPy-like automatic evaluation.

use std::sync::Arc;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{
    add, add2, could_extract_minus_sign, div, mul, mul2, neg, pow, sqrt, sub,
};

/// Apply function `f` to `args`, performing automatic simplifications.
pub fn func(f: Func, args: Vec<Expr>) -> Expr {
    match f {
        Func::Logic(op) => return crate::logic::eval(op, args),
        Func::Set(op) => return crate::sets::eval(op, args),
        Func::Order => return crate::series::order::order_from_args(&args),
        _ => {}
    }
    if args.iter().any(|a| a.is_constant(Constant::NaN)) && !matches!(f, Func::Undefined(_)) {
        return Expr::nan();
    }
    if let Func::Special(s) = f {
        return crate::special::apply(s, args);
    }
    // numeric evaluation when an argument is a float
    if !matches!(f, Func::Undefined(_) | Func::Max | Func::Min)
        && !args.is_empty()
        && args.iter().all(is_numeric_value)
        && args.iter().any(has_float)
    {
        if let Some(v) = crate::evalf::eval_func_numeric(&f, &args) {
            return v;
        }
    }
    let r = match &f {
        Func::Exp => eval_exp(&args[0]),
        Func::Log => eval_log(&args[0]),
        Func::Sin => eval_sin(&args[0]),
        Func::Cos => eval_cos(&args[0]),
        Func::Tan => eval_tan(&args[0]),
        Func::Cot => eval_cot(&args[0]),
        Func::Sec => eval_sec(&args[0]),
        Func::Csc => eval_csc(&args[0]),
        Func::Asin => eval_asin(&args[0]),
        Func::Acos => eval_acos(&args[0]),
        Func::Atan => eval_atan(&args[0]),
        Func::Acot => eval_acot(&args[0]),
        Func::Atan2 => eval_atan2(&args[0], &args[1]),
        Func::Sinh => eval_sinh(&args[0]),
        Func::Cosh => eval_cosh(&args[0]),
        Func::Tanh => eval_tanh(&args[0]),
        Func::Coth => eval_odd(&Func::Coth, &args[0], |_| None),
        Func::Asinh => eval_odd(&Func::Asinh, &args[0], |a| {
            if a.is_zero() {
                Some(Expr::zero())
            } else {
                None
            }
        }),
        Func::Acosh => {
            if args[0].is_one() {
                Some(Expr::zero())
            } else {
                None
            }
        }
        Func::Atanh => eval_odd(&Func::Atanh, &args[0], |a| {
            if a.is_zero() {
                Some(Expr::zero())
            } else if a.is_one() {
                Some(Expr::oo())
            } else {
                None
            }
        }),
        Func::Abs => eval_abs(&args[0]),
        Func::Sign => eval_sign(&args[0]),
        Func::Floor | Func::Ceiling => eval_floor_ceil(&f, &args[0]),
        Func::Gamma => crate::special::eval_gamma(&args[0]),
        Func::LogGamma => crate::special::eval_loggamma(&args[0]),
        Func::Erf => crate::special::eval_erf(&args[0]),
        Func::Erfc => crate::special::eval_erfc(&args[0]),
        Func::Re => eval_re_im(&args[0], true),
        Func::Im => eval_re_im(&args[0], false),
        Func::Conjugate => eval_conjugate(&args[0]),
        Func::Arg => match args[0].kind() {
            Kind::Num(n) if n.is_positive() => Some(Expr::zero()),
            Kind::Num(n) if n.is_negative() => Some(Expr::pi()),
            Kind::Const(Constant::I) => Some(div(&Expr::pi(), &Expr::int(2))),
            _ => None,
        },
        Func::Heaviside => match args[0].kind() {
            Kind::Num(n) if n.is_positive() => Some(Expr::one()),
            Kind::Num(n) if n.is_negative() => Some(Expr::zero()),
            Kind::Num(_) => Some(Expr::rational(1, 2)),
            _ => None,
        },
        Func::DiracDelta => match args[0].kind() {
            Kind::Num(n) if !n.is_zero() => Some(Expr::zero()),
            _ => None,
        },
        Func::Max | Func::Min => return eval_minmax(&f, args),
        Func::Special(_) | Func::Order | Func::Undefined(_) | Func::Logic(_) | Func::Set(_) => None,
    };
    match r {
        Some(v) => v,
        None => Expr::from_kind(Kind::Func(f, args)),
    }
}

fn has_float(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_float(),
        _ => e.args().iter().any(has_float),
    }
}

/// Free of symbols and undefined functions (so it can be evaluated).
pub fn is_numeric_value(e: &Expr) -> bool {
    // any symbol (or Order term) inside sets a mask bit: O(1) rejection
    if e.mask() != 0 {
        return false;
    }
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(c) => !c.is_infinite() && *c != Constant::NaN,
        Kind::Sym(_) | Kind::Derivative(..) | Kind::Rel(..) | Kind::Piecewise(_) | Kind::Bool(_) => {
            false
        }
        Kind::Func(Func::Logic(_) | Func::Set(_), _) => false,
        Kind::Func(Func::Undefined(n), a) => &**n == "Tuple" && a.iter().all(is_numeric_value),
        _ => e.args().iter().all(is_numeric_value),
    }
}

// --------------------------------------------------------------------
// constructors
// --------------------------------------------------------------------

macro_rules! unary {
    ($($name:ident => $f:expr),* $(,)?) => {
        $(
            pub fn $name(x: &Expr) -> Expr { func($f, vec![x.clone()]) }
        )*
    };
}

unary! {
    exp => Func::Exp, log => Func::Log,
    sin => Func::Sin, cos => Func::Cos, tan => Func::Tan,
    cot => Func::Cot, sec => Func::Sec, csc => Func::Csc,
    asin => Func::Asin, acos => Func::Acos, atan => Func::Atan, acot => Func::Acot,
    sinh => Func::Sinh, cosh => Func::Cosh, tanh => Func::Tanh, coth => Func::Coth,
    asinh => Func::Asinh, acosh => Func::Acosh, atanh => Func::Atanh,
    abs => Func::Abs, sign => Func::Sign, floor => Func::Floor, ceiling => Func::Ceiling,
    gamma => Func::Gamma, loggamma => Func::LogGamma, erf => Func::Erf, erfc => Func::Erfc,
    re => Func::Re, im => Func::Im, conjugate => Func::Conjugate, arg => Func::Arg,
    heaviside => Func::Heaviside, dirac_delta => Func::DiracDelta,
}

pub fn atan2(y: &Expr, x: &Expr) -> Expr {
    func(Func::Atan2, vec![y.clone(), x.clone()])
}

pub fn max(args: Vec<Expr>) -> Expr {
    func(Func::Max, args)
}

pub fn min(args: Vec<Expr>) -> Expr {
    func(Func::Min, args)
}

/// Logarithm in base `b`: `log(x)/log(b)`.
pub fn log_base(x: &Expr, b: &Expr) -> Expr {
    div(&log(x), &log(b))
}

/// Applied undefined function `name(args...)`.
pub fn undefined(name: &str, args: Vec<Expr>) -> Expr {
    func(Func::Undefined(Arc::from(name)), args)
}

// --------------------------------------------------------------------
// helpers
// --------------------------------------------------------------------

/// If `e == k*pi` with rational `k`, return `k`.
fn pi_coeff(e: &Expr) -> Option<Number> {
    match e.kind() {
        Kind::Num(n) if n.is_zero() => Some(Number::ZERO),
        Kind::Const(Constant::Pi) => Some(Number::ONE),
        Kind::Mul(c, facs)
            if c.is_rational() && facs.len() == 1 && facs[0].0.is_constant(Constant::Pi) && facs[0].1.is_one() =>
        {
            Some(c.clone())
        }
        _ => None,
    }
}

/// Split `rest + k*pi` (k rational), returning `(rest, k)`.
fn split_pi_shift(e: &Expr) -> Option<(Expr, Number)> {
    if let Kind::Add(c, terms) = e.kind() {
        for (i, (t, k)) in terms.iter().enumerate() {
            if t.is_constant(Constant::Pi) && k.is_rational() {
                let mut rest: Vec<(Expr, Number)> = terms.clone();
                rest.remove(i);
                let rest = crate::ops::add_from_terms(c.clone(), rest);
                return Some((rest, k.clone()));
            }
        }
    }
    None
}

/// `e / I` if `e` is syntactically `I * something`.
fn extract_i(e: &Expr) -> Option<Expr> {
    match e.kind() {
        Kind::Const(Constant::I) => Some(Expr::one()),
        Kind::Mul(c, facs) => {
            let idx = facs
                .iter()
                .position(|(b, x)| b.is_constant(Constant::I) && x.is_one())?;
            let mut fs = facs.clone();
            fs.remove(idx);
            let mut v: Vec<Expr> = fs.iter().map(|(b, x)| crate::ops::mul_factor(b, x)).collect();
            v.push(Expr::num(c.clone()));
            Some(mul(&v))
        }
        _ => None,
    }
}

/// Reduce a rational `k` modulo `m` into `[0, m)`.
fn rmod(k: &Number, m: i64) -> Number {
    let m = Number::int(m);
    let q = k.div(&m).unwrap().floor();
    k.sub(&q.mul(&m))
}

fn q(n: i64, d: i64) -> Number {
    Number::rational_i64(n, d).unwrap()
}

fn half() -> Expr {
    Expr::rational(1, 2)
}

fn sqrt_i(n: i64) -> Expr {
    sqrt(&Expr::int(n))
}

fn eval_odd(f: &Func, a: &Expr, special: impl Fn(&Expr) -> Option<Expr>) -> Option<Expr> {
    if let Some(v) = special(a) {
        return Some(v);
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&func(f.clone(), vec![neg(a)])));
    }
    None
}

// --------------------------------------------------------------------
// exp / log
// --------------------------------------------------------------------

fn eval_exp(a: &Expr) -> Option<Expr> {
    match a.kind() {
        Kind::Num(n) if n.is_zero() => return Some(Expr::one()),
        Kind::Num(n) if n.is_one() => return Some(Expr::e()),
        Kind::Const(Constant::Infinity) => return Some(Expr::oo()),
        Kind::Const(Constant::NegInfinity) => return Some(Expr::zero()),
        Kind::Const(Constant::ComplexInfinity) => return Some(Expr::nan()),
        Kind::Func(Func::Log, x) => return Some(x[0].clone()),
        Kind::Mul(c, facs) if c.is_rational() && facs.len() == 1 => {
            // exp(k*log(x)) = x**k
            if let (Kind::Func(Func::Log, x), true) = (facs[0].0.kind(), facs[0].1.is_one()) {
                return Some(pow(&x[0], &Expr::num(c.clone())));
            }
        }
        _ => {}
    }
    // exp(k*I*pi) for k multiple of 1/2
    if let Some(r) = extract_i(a) {
        if let Some(k) = pi_coeff(&r) {
            let k2 = k.mul(&Number::int(2));
            if k2.is_integer() {
                let m = rmod(&k2, 4).to_i64().unwrap();
                return Some(match m {
                    0 => Expr::one(),
                    1 => Expr::i(),
                    2 => Expr::minus_one(),
                    _ => neg(&Expr::i()),
                });
            }
        }
    }
    None
}

fn eval_log(a: &Expr) -> Option<Expr> {
    match a.kind() {
        Kind::Num(n) => {
            if n.is_one() {
                return Some(Expr::zero());
            }
            if n.is_zero() {
                return Some(Expr::zoo());
            }
            if n.is_negative() {
                // log(-n) = log(n) + I*pi
                return Some(add2(&log(&Expr::num(n.neg())), &mul2(&Expr::i(), &Expr::pi())));
            }
            if let Number::Rational(r) = n {
                if r.numer() == &num_bigint::BigInt::from(1) {
                    return Some(neg(&log(&Expr::num(Number::from_bigint(r.denom().clone())))));
                }
            }
            None
        }
        Kind::Const(Constant::E) => Some(Expr::one()),
        Kind::Func(Func::Exp, x) if crate::assumptions::is_real(&x[0]) == Some(true) => Some(x[0].clone()),
        Kind::Const(Constant::I) => Some(mul(&[Expr::i(), Expr::pi(), half()])),
        Kind::Const(Constant::Infinity) | Kind::Const(Constant::NegInfinity) => Some(Expr::oo()),
        Kind::Const(Constant::ComplexInfinity) => Some(Expr::zoo()),
        _ => None,
    }
}

// --------------------------------------------------------------------
// trigonometric
// --------------------------------------------------------------------

/// sin(k*pi) for rational k already reduced into [0, 1/2].
fn sin_table(k: &Number) -> Option<Expr> {
    let v = if k.is_zero() {
        Expr::zero()
    } else if *k == q(1, 6) {
        half()
    } else if *k == q(1, 4) {
        div(&sqrt_i(2), &Expr::int(2))
    } else if *k == q(1, 3) {
        div(&sqrt_i(3), &Expr::int(2))
    } else if *k == q(1, 2) {
        Expr::one()
    } else if *k == q(1, 12) {
        div(&sub(&sqrt_i(6), &sqrt_i(2)), &Expr::int(4))
    } else if *k == q(5, 12) {
        div(&add2(&sqrt_i(6), &sqrt_i(2)), &Expr::int(4))
    } else {
        return None;
    };
    Some(v)
}

/// sin(k*pi) for any rational k: exact value or canonical reduced form.
fn sin_pi(k: &Number) -> Expr {
    let mut k = rmod(k, 2);
    let mut sign = false;
    if k.cmp_value(&Number::ONE) != Some(std::cmp::Ordering::Less) {
        k = k.sub(&Number::ONE);
        sign = true;
    }
    if k.cmp_value(&q(1, 2)) == Some(std::cmp::Ordering::Greater) {
        k = Number::ONE.sub(&k);
    }
    let v = match sin_table(&k) {
        Some(v) => v,
        None => Expr::from_kind(Kind::Func(Func::Sin, vec![mul2(&Expr::num(k), &Expr::pi())])),
    };
    if sign {
        neg(&v)
    } else {
        v
    }
}

fn cos_pi(k: &Number) -> Expr {
    // cos(k pi) = sin((1/2 - k) pi) but keep a `cos` node when not exact
    let mut k = rmod(k, 2);
    if k.cmp_value(&Number::ONE) == Some(std::cmp::Ordering::Greater) {
        k = Number::int(2).sub(&k);
    }
    let mut sign = false;
    if k.cmp_value(&q(1, 2)) == Some(std::cmp::Ordering::Greater) {
        k = Number::ONE.sub(&k);
        sign = true;
    }
    let s = q(1, 2).sub(&k);
    let v = match sin_table(&s) {
        Some(v) => v,
        None => Expr::from_kind(Kind::Func(Func::Cos, vec![mul2(&Expr::num(k), &Expr::pi())])),
    };
    if sign {
        neg(&v)
    } else {
        v
    }
}

/// `k` with `a = k*pi` and `k` a known (symbolic) integer.
fn integer_pi_multiple(a: &Expr) -> Option<Expr> {
    if let Kind::Mul(_, facs) = a.kind() {
        if facs.iter().any(|(b, e)| b.is_constant(Constant::Pi) && e.is_one()) {
            let k = div(a, &Expr::pi());
            if !k.is_number() && crate::assumptions::is_integer(&k) == Some(true) {
                return Some(k);
            }
        }
    }
    None
}

fn eval_sin(a: &Expr) -> Option<Expr> {
    if integer_pi_multiple(a).is_some() {
        return Some(Expr::zero());
    }
    if let Some(k) = pi_coeff(a) {
        let v = sin_pi(&k);
        return Some(v);
    }
    match a.kind() {
        Kind::Func(Func::Asin, x) => return Some(x[0].clone()),
        Kind::Func(Func::Acos, x) => return Some(sqrt(&sub(&Expr::one(), &x[0].powi(2)))),
        Kind::Func(Func::Atan, x) => {
            return Some(div(&x[0], &sqrt(&add2(&x[0].powi(2), &Expr::one()))))
        }
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => return None,
        _ => {}
    }
    if let Some(r) = extract_i(a) {
        return Some(mul2(&Expr::i(), &sinh(&r)));
    }
    if let Some((rest, k)) = split_pi_shift(a) {
        let k = rmod(&k, 2);
        if k.is_zero() {
            return Some(sin(&rest));
        } else if k.is_one() {
            return Some(neg(&sin(&rest)));
        } else if k == q(1, 2) {
            return Some(cos(&rest));
        } else if k == q(3, 2) {
            return Some(neg(&cos(&rest)));
        }
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&sin(&neg(a))));
    }
    None
}

fn eval_cos(a: &Expr) -> Option<Expr> {
    if let Some(k) = integer_pi_multiple(a) {
        return Some(pow(&Expr::minus_one(), &k));
    }
    if let Some(k) = pi_coeff(a) {
        return Some(cos_pi(&k));
    }
    match a.kind() {
        Kind::Func(Func::Acos, x) => return Some(x[0].clone()),
        Kind::Func(Func::Asin, x) => return Some(sqrt(&sub(&Expr::one(), &x[0].powi(2)))),
        Kind::Func(Func::Atan, x) => {
            return Some(pow(&add2(&x[0].powi(2), &Expr::one()), &Expr::rational(-1, 2)))
        }
        _ => {}
    }
    if let Some(r) = extract_i(a) {
        return Some(cosh(&r));
    }
    if let Some((rest, k)) = split_pi_shift(a) {
        let k = rmod(&k, 2);
        if k.is_zero() {
            return Some(cos(&rest));
        } else if k.is_one() {
            return Some(neg(&cos(&rest)));
        } else if k == q(1, 2) {
            return Some(neg(&sin(&rest)));
        } else if k == q(3, 2) {
            return Some(sin(&rest));
        }
    }
    if could_extract_minus_sign(a) {
        return Some(cos(&neg(a)));
    }
    None
}

fn eval_tan(a: &Expr) -> Option<Expr> {
    if integer_pi_multiple(a).is_some() {
        return Some(Expr::zero());
    }
    if let Some(k) = pi_coeff(a) {
        let k = rmod(&k, 1);
        let s = sin_pi(&k);
        let c = cos_pi(&k);
        if c.is_zero() {
            return Some(Expr::zoo());
        }
        if let (Kind::Func(..), _) | (_, Kind::Func(..)) = (s.kind(), c.kind()) {
            let k = if k.cmp_value(&q(1, 2)) == Some(std::cmp::Ordering::Greater) {
                // tan(k pi) = -tan((1-k) pi)
                return Some(neg(&Expr::from_kind(Kind::Func(
                    Func::Tan,
                    vec![mul2(&Expr::num(Number::ONE.sub(&k)), &Expr::pi())],
                ))));
            } else {
                k
            };
            return Some(Expr::from_kind(Kind::Func(Func::Tan, vec![mul2(&Expr::num(k), &Expr::pi())])));
        }
        return Some(crate::simplify::radsimp(&div(&s, &c)));
    }
    if let Kind::Func(Func::Atan, x) = a.kind() {
        return Some(x[0].clone());
    }
    if let Some(r) = extract_i(a) {
        return Some(mul2(&Expr::i(), &tanh(&r)));
    }
    if let Some((rest, k)) = split_pi_shift(a) {
        let k = rmod(&k, 1);
        if k.is_zero() {
            return Some(tan(&rest));
        } else if k == q(1, 2) {
            return Some(neg(&cot(&rest)));
        }
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&tan(&neg(a))));
    }
    None
}

fn eval_cot(a: &Expr) -> Option<Expr> {
    if let Some(k) = pi_coeff(a) {
        let k = rmod(&k, 1);
        let s = sin_pi(&k);
        let c = cos_pi(&k);
        if s.is_zero() {
            return Some(Expr::zoo());
        }
        if !matches!(s.kind(), Kind::Func(..)) && !matches!(c.kind(), Kind::Func(..)) {
            return Some(crate::simplify::radsimp(&div(&c, &s)));
        }
    }
    if let Kind::Func(Func::Acot, x) = a.kind() {
        return Some(x[0].clone());
    }
    if let Some((rest, k)) = split_pi_shift(a) {
        let k = rmod(&k, 1);
        if k.is_zero() {
            return Some(cot(&rest));
        } else if k == q(1, 2) {
            return Some(neg(&tan(&rest)));
        }
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&cot(&neg(a))));
    }
    None
}

fn eval_sec(a: &Expr) -> Option<Expr> {
    if let Some(k) = pi_coeff(a) {
        let c = cos_pi(&k);
        if c.is_zero() {
            return Some(Expr::zoo());
        }
        if !matches!(c.kind(), Kind::Func(..)) {
            return Some(crate::simplify::radsimp(&div(&Expr::one(), &c)));
        }
    }
    if could_extract_minus_sign(a) {
        return Some(sec(&neg(a)));
    }
    None
}

fn eval_csc(a: &Expr) -> Option<Expr> {
    if let Some(k) = pi_coeff(a) {
        let s = sin_pi(&k);
        if s.is_zero() {
            return Some(Expr::zoo());
        }
        if !matches!(s.kind(), Kind::Func(..)) {
            return Some(crate::simplify::radsimp(&div(&Expr::one(), &s)));
        }
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&csc(&neg(a))));
    }
    None
}

/// asin of special algebraic values, as a multiple of pi.
fn asin_table(a: &Expr) -> Option<Number> {
    let candidates = [
        (Expr::zero(), q(0, 1)),
        (half(), q(1, 6)),
        (div(&sqrt_i(2), &Expr::int(2)), q(1, 4)),
        (div(&sqrt_i(3), &Expr::int(2)), q(1, 3)),
        (Expr::one(), q(1, 2)),
        (div(&sub(&sqrt_i(6), &sqrt_i(2)), &Expr::int(4)), q(1, 12)),
        (div(&add2(&sqrt_i(6), &sqrt_i(2)), &Expr::int(4)), q(5, 12)),
    ];
    for (v, k) in candidates {
        if *a == v {
            return Some(k);
        }
    }
    None
}

fn eval_asin(a: &Expr) -> Option<Expr> {
    if let Some(k) = asin_table(a) {
        return Some(mul2(&Expr::num(k), &Expr::pi()));
    }
    if let Some(r) = extract_i(a) {
        return Some(mul2(&Expr::i(), &asinh(&r)));
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&asin(&neg(a))));
    }
    None
}

fn eval_acos(a: &Expr) -> Option<Expr> {
    if let Some(k) = asin_table(a) {
        return Some(mul2(&Expr::num(q(1, 2).sub(&k)), &Expr::pi()));
    }
    if could_extract_minus_sign(a) {
        let n = neg(a);
        if let Some(k) = asin_table(&n) {
            return Some(mul2(&Expr::num(q(1, 2).add(&k)), &Expr::pi()));
        }
    }
    None
}

fn eval_atan(a: &Expr) -> Option<Expr> {
    let candidates = [
        (Expr::zero(), q(0, 1)),
        (div(&sqrt_i(3), &Expr::int(3)), q(1, 6)),
        (Expr::one(), q(1, 4)),
        (sqrt_i(3), q(1, 3)),
        (sub(&Expr::int(2), &sqrt_i(3)), q(1, 12)),
        (add2(&Expr::int(2), &sqrt_i(3)), q(5, 12)),
    ];
    for (v, k) in candidates {
        if *a == v {
            return Some(mul2(&Expr::num(k), &Expr::pi()));
        }
    }
    match a.kind() {
        Kind::Const(Constant::Infinity) => return Some(div(&Expr::pi(), &Expr::int(2))),
        Kind::Const(Constant::NegInfinity) => return Some(div(&Expr::pi(), &Expr::int(-2))),
        _ => {}
    }
    if let Some(r) = extract_i(a) {
        return Some(mul2(&Expr::i(), &atanh(&r)));
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&atan(&neg(a))));
    }
    None
}

fn eval_acot(a: &Expr) -> Option<Expr> {
    if a.is_zero() {
        return Some(div(&Expr::pi(), &Expr::int(2)));
    }
    if a.is_one() {
        return Some(div(&Expr::pi(), &Expr::int(4)));
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&acot(&neg(a))));
    }
    None
}

fn eval_atan2(y: &Expr, x: &Expr) -> Option<Expr> {
    let (yn, xn) = (y.as_number()?, x.as_number()?);
    use std::cmp::Ordering::*;
    let pi = Expr::pi();
    match (yn.cmp_value(&Number::ZERO)?, xn.cmp_value(&Number::ZERO)?) {
        (Equal, Equal) => Some(Expr::nan()),
        (Equal, Greater) => Some(Expr::zero()),
        (Equal, Less) => Some(pi),
        (Greater, Equal) => Some(div(&pi, &Expr::int(2))),
        (Less, Equal) => Some(div(&pi, &Expr::int(-2))),
        (_, Greater) => Some(atan(&div(y, x))),
        (Greater, Less) => Some(add2(&atan(&div(y, x)), &pi)),
        (Less, Less) => Some(sub(&atan(&div(y, x)), &pi)),
    }
}

// --------------------------------------------------------------------
// hyperbolic
// --------------------------------------------------------------------

fn eval_sinh(a: &Expr) -> Option<Expr> {
    if a.is_zero() {
        return Some(Expr::zero());
    }
    if let Kind::Func(Func::Asinh, x) = a.kind() {
        return Some(x[0].clone());
    }
    if let Some(r) = extract_i(a) {
        return Some(mul2(&Expr::i(), &sin(&r)));
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&sinh(&neg(a))));
    }
    None
}

fn eval_cosh(a: &Expr) -> Option<Expr> {
    if a.is_zero() {
        return Some(Expr::one());
    }
    if let Kind::Func(Func::Acosh, x) = a.kind() {
        return Some(x[0].clone());
    }
    if let Some(r) = extract_i(a) {
        return Some(cos(&r));
    }
    if could_extract_minus_sign(a) {
        return Some(cosh(&neg(a)));
    }
    None
}

fn eval_tanh(a: &Expr) -> Option<Expr> {
    if a.is_zero() {
        return Some(Expr::zero());
    }
    if let Kind::Func(Func::Atanh, x) = a.kind() {
        return Some(x[0].clone());
    }
    if let Some(r) = extract_i(a) {
        return Some(mul2(&Expr::i(), &tan(&r)));
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&tanh(&neg(a))));
    }
    None
}

// --------------------------------------------------------------------
// misc
// --------------------------------------------------------------------

fn eval_abs(a: &Expr) -> Option<Expr> {
    if !a.is_number() {
        if crate::assumptions::is_nonnegative(a) == Some(true) {
            return Some(a.clone());
        }
        if crate::assumptions::is_nonpositive(a) == Some(true) {
            return Some(neg(a));
        }
    }
    match a.kind() {
        Kind::Num(n) => Some(Expr::num(n.abs())),
        Kind::Const(Constant::I) => Some(Expr::one()),
        Kind::Const(Constant::Pi | Constant::E | Constant::EulerGamma) => Some(a.clone()),
        Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity) => {
            Some(Expr::oo())
        }
        Kind::Func(Func::Abs, _) => Some(a.clone()),
        Kind::Mul(c, facs) => {
            // Abs(c * x) = |c| * Abs(x),  Abs(I*x) = Abs(x)
            let mut rest: Vec<Expr> = Vec::new();
            let mut changed = !c.is_one();
            for (b, e) in facs {
                if b.is_constant(Constant::I) && e.is_one() {
                    changed = true;
                    continue;
                }
                if matches!(b.kind(), Kind::Const(Constant::Pi | Constant::E))
                    || (b.is_number() && b.as_number().unwrap().is_positive() && e.is_rational_number())
                    || crate::assumptions::is_positive(&crate::ops::mul_factor(b, e)) == Some(true)
                {
                    // positive factor
                    return Some(mul2(
                        &crate::ops::mul_factor(b, e),
                        &abs(&div(a, &crate::ops::mul_factor(b, e))),
                    ));
                }
                rest.push(crate::ops::mul_factor(b, e));
            }
            if changed {
                Some(mul2(&Expr::num(c.abs()), &abs(&mul(&rest))))
            } else {
                None
            }
        }
        _ => {
            if could_extract_minus_sign(a) {
                Some(abs(&neg(a)))
            } else {
                None
            }
        }
    }
}

fn eval_sign(a: &Expr) -> Option<Expr> {
    if !a.is_number() {
        if crate::assumptions::is_positive(a) == Some(true) {
            return Some(Expr::one());
        }
        if crate::assumptions::is_negative(a) == Some(true) {
            return Some(Expr::minus_one());
        }
        if crate::assumptions::is_zero(a) == Some(true) {
            return Some(Expr::zero());
        }
    }
    match a.kind() {
        Kind::Num(n) => Some(Expr::int(if n.is_positive() {
            1
        } else if n.is_negative() {
            -1
        } else {
            0
        })),
        Kind::Const(Constant::Pi | Constant::E | Constant::Infinity) => Some(Expr::one()),
        Kind::Const(Constant::NegInfinity) => Some(Expr::minus_one()),
        Kind::Const(Constant::I) => Some(Expr::i()),
        _ => {
            if could_extract_minus_sign(a) {
                Some(neg(&sign(&neg(a))))
            } else {
                None
            }
        }
    }
}

fn eval_floor_ceil(f: &Func, a: &Expr) -> Option<Expr> {
    let ceil = matches!(f, Func::Ceiling);
    if !a.is_number() && crate::assumptions::is_integer(a) == Some(true) {
        return Some(a.clone());
    }
    match a.kind() {
        Kind::Num(n) => {
            let fl = match n {
                Number::Float(x) => {
                    let v = if ceil { x.ceil() } else { x.floor() };
                    if v.abs() < 9.0e15 {
                        return Some(Expr::int(v as i64));
                    }
                    return Some(Expr::float(v));
                }
                _ => n.floor(),
            };
            if ceil && fl != *n {
                Some(Expr::num(fl.add(&Number::ONE)))
            } else {
                Some(Expr::num(fl))
            }
        }
        _ => {
            if crate::functions::is_numeric_value(a) {
                // evaluate numerically with a safety margin
                let v = crate::evalf::eval_f64(a)?;
                let r = if ceil { v.ceil() } else { v.floor() };
                if (v - r).abs() > 1e-9 && r.abs() < 9.0e15 {
                    return Some(Expr::int(r as i64));
                }
            }
            None
        }
    }
}

/// Split a numeric value into exact real and imaginary parts when trivially
/// possible.
fn eval_re_im(a: &Expr, want_re: bool) -> Option<Expr> {
    if !a.is_number() && crate::assumptions::is_real(a) == Some(true) {
        return Some(if want_re { a.clone() } else { Expr::zero() });
    }
    match a.kind() {
        Kind::Num(_) | Kind::Const(Constant::Pi | Constant::E | Constant::EulerGamma) => {
            Some(if want_re { a.clone() } else { Expr::zero() })
        }
        Kind::Const(Constant::I) => Some(if want_re { Expr::zero() } else { Expr::one() }),
        Kind::Add(c, terms) => {
            let mut v = vec![if want_re { Expr::num(c.clone()) } else { Expr::zero() }];
            for (t, k) in terms {
                let tt = crate::ops::add_term(t, k);
                v.push(func(if want_re { Func::Re } else { Func::Im }, vec![tt]));
            }
            Some(add(&v))
        }
        Kind::Mul(c, _) if c.is_float() || !c.is_one() => {
            let rest = div(a, &Expr::num(c.clone()));
            let inner = func(if want_re { Func::Re } else { Func::Im }, vec![rest]);
            Some(mul2(&Expr::num(c.clone()), &inner))
        }
        Kind::Mul(_, _) => {
            if let Some(r) = extract_i(a) {
                // re(I*x) = -im(x), im(I*x) = re(x)
                return Some(if want_re { neg(&im(&r)) } else { re(&r) });
            }
            None
        }
        Kind::Func(Func::Re | Func::Im | Func::Abs, _) => {
            Some(if want_re { a.clone() } else { Expr::zero() })
        }
        _ => None,
    }
}

fn eval_conjugate(a: &Expr) -> Option<Expr> {
    if !a.is_number() && crate::assumptions::is_real(a) == Some(true) {
        return Some(a.clone());
    }
    match a.kind() {
        Kind::Num(_) | Kind::Const(Constant::Pi | Constant::E | Constant::EulerGamma) => {
            Some(a.clone())
        }
        Kind::Const(Constant::I) => Some(neg(a)),
        Kind::Func(Func::Conjugate, x) => Some(x[0].clone()),
        Kind::Pow(base, exponent) if crate::assumptions::is_integer(exponent) == Some(true) => {
            Some(pow(&conjugate(base), exponent))
        }
        Kind::Add(..) | Kind::Mul(..) => {
            let args: Vec<Expr> = a.args().iter().map(conjugate).collect();
            Some(if a.is_add() { add(&args) } else { mul(&args) })
        }
        _ => None,
    }
}

pub(crate) fn contains_i(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::I) => true,
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(contains_i),
    }
}

fn eval_minmax(f: &Func, args: Vec<Expr>) -> Expr {
    let is_max = matches!(f, Func::Max);
    let mut flat: Vec<Expr> = Vec::new();
    for a in args {
        match a.kind() {
            Kind::Func(g, inner) if g == f => flat.extend(inner.iter().cloned()),
            _ => flat.push(a),
        }
    }
    // fold numeric values
    let mut best: Option<(f64, Expr)> = None;
    let mut rest: Vec<Expr> = Vec::new();
    for a in flat {
        let v = if is_numeric_value(&a) && !crate::functions::contains_i(&a) {
            crate::evalf::eval_f64(&a)
        } else {
            None
        };
        match v {
            Some(v) => {
                let better = match &best {
                    None => true,
                    Some((b, _)) => {
                        if is_max {
                            v > *b
                        } else {
                            v < *b
                        }
                    }
                };
                if better {
                    best = Some((v, a));
                }
            }
            None => rest.push(a),
        }
    }
    if let Some((_, b)) = best {
        rest.push(b);
    }
    rest.sort();
    rest.dedup();
    // drop arguments dominated by another one (known sign of the difference)
    let mut i = 0;
    while i < rest.len() {
        let dominated = (0..rest.len()).any(|j| {
            j != i && {
                let d = sub(&rest[j], &rest[i]);
                if is_max {
                    crate::assumptions::is_positive(&d) == Some(true)
                } else {
                    crate::assumptions::is_negative(&d) == Some(true)
                }
            }
        });
        if dominated {
            rest.remove(i);
        } else {
            i += 1;
        }
    }
    if rest.len() == 1 {
        return rest.pop().unwrap();
    }
    if rest.is_empty() {
        return if is_max { Expr::neg_oo() } else { Expr::oo() };
    }
    Expr::from_kind(Kind::Func(f.clone(), rest))
}
