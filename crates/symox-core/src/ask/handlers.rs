//! Predicate handlers, ported from SymPy 1.14 `sympy/assumptions/handlers`
//! and `sympy/assumptions/relation`.
//!
//! SymPy dispatches on the class of the argument: the most specific
//! registered handler wins and only `MDNotImplementedError` falls through to
//! a more generic one (`Expr`, then `Basic`, then `object`). Each function
//! below tries the class specific handlers first and returns; the code that
//! follows the class `match` is the generic fallback. The `Expr` fallback is
//! the old assumption (`expr.is_<name>`). Numerical tests evaluate with
//! native arbitrary precision at two working precisions and only conclude
//! when both agree, so cancellations give `None` like SymPy's `_prec == 1`.

use std::cmp::Ordering;
use std::sync::Arc;

use num_integer::Integer as _;

use super::engine::{is_boolean_kind, is_matrix, Cx, R};
use super::predicates::Pred;
use super::{applied_unchecked, as_applied, q, AskError};
use crate::bigfloat::BigFloat;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::matrixexpr::{MatrixExpr, MatrixKind};
use crate::special::SpecialFn;

type O = Option<bool>;

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

fn tru() -> Expr {
    Expr::boolean(true)
}

fn not(e: Expr) -> Expr {
    crate::logic::mk(crate::logic::LogicOp::Not, vec![e])
}

fn and(v: Vec<Expr>) -> Expr {
    crate::logic::mk(crate::logic::LogicOp::And, v)
}

fn or(v: Vec<Expr>) -> Expr {
    crate::logic::mk(crate::logic::LogicOp::Or, v)
}

fn fnot(x: O) -> O {
    x.map(|v| !v)
}

/// `_fuzzy_group(..., quick_exit=True)` over `ask(key(arg))`.
fn closed_group(cx: &Cx, p: Pred, args: &[Expr], a: &Expr) -> R {
    let mut saw_other = false;
    for x in args {
        match cx.askp(p, x, a)? {
            Some(true) => {}
            None => return Ok(None),
            Some(false) => {
                if saw_other {
                    return Ok(None);
                }
                saw_other = true;
            }
        }
    }
    Ok(Some(!saw_other))
}

/// `fuzzy_and(ask(q) for q in queries)`, lazily.
fn ask_all(cx: &Cx, queries: &[Expr], a: &Expr) -> R {
    let mut rv = Some(true);
    for qq in queries {
        match cx.ask(qq, a)? {
            Some(false) => return Ok(Some(false)),
            None => rv = None,
            _ => {}
        }
    }
    Ok(rv)
}

/// `fuzzy_or(ask(q) for q in queries)`, lazily.
fn ask_any(cx: &Cx, queries: &[Expr], a: &Expr) -> R {
    let mut rv = Some(false);
    for qq in queries {
        match cx.ask(qq, a)? {
            Some(true) => return Ok(Some(true)),
            None => rv = None,
            _ => {}
        }
    }
    Ok(rv)
}

/// The old assumption `expr.is_<name>` (the `Expr` handlers).
fn old(e: &Expr, p: Pred) -> O {
    let name = p.old_assumption()?;
    if is_matrix(e) {
        // MatrixExpr is noncommutative: every fact implying commutativity is False.
        return match name {
            "finite" | "infinite" | "hermitian" | "antihermitian" => None,
            _ => Some(false),
        };
    }
    if !is_expr(e) {
        return None;
    }
    crate::assumptions::fact(e, name)
}

/// SymPy `isinstance(e, Expr)`.
fn is_expr(e: &Expr) -> bool {
    !is_boolean_kind(e) && !matches!(e.kind(), Kind::Func(Func::Predicate(_) | Func::Logic(_), _))
}

fn is_number(e: &Expr) -> bool {
    !is_matrix(e) && is_expr(e) && crate::symbolic_number::is_symbolic_number(e)
}

fn args_of(e: &Expr) -> Vec<Expr> {
    crate::traverse::ordered_args(e)
}

fn is_nan(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::NaN) => true,
        Kind::Num(n) => n.is_nan(),
        _ => false,
    }
}

/// SymPy `Integer`.
fn is_int(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(n) if !n.is_float() && !n.is_nan() && n.is_integer())
}

/// SymPy `Rational` (including `Integer`).
fn is_rat(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(n) if !n.is_float() && !n.is_nan())
}

fn is_float(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(n) if n.is_float() && !n.is_nan())
}

fn cst(e: &Expr) -> Option<Constant> {
    match e.kind() {
        Kind::Const(c) => Some(*c),
        _ => None,
    }
}

/// SymPy `Number` (numbers, `oo`, `-oo` and `nan`).
fn is_number_class(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(_))
        || matches!(
            cst(e),
            Some(Constant::Infinity | Constant::NegInfinity | Constant::NaN)
        )
}

/// SymPy `NumberSymbol` (`pi`, `E`, `EulerGamma`).
fn is_number_symbol(e: &Expr) -> bool {
    matches!(
        cst(e),
        Some(Constant::Pi | Constant::E | Constant::EulerGamma)
    )
}

fn func_args<'e>(e: &'e Expr, f: &Func) -> Option<&'e [Expr]> {
    match e.kind() {
        Kind::Func(g, args) if g == f => Some(args),
        _ => None,
    }
}

fn arg1<'e>(e: &'e Expr, f: Func) -> Option<&'e Expr> {
    func_args(e, &f).and_then(|a| a.first())
}

/// The base of `exp(x)` style powers: `x` for `exp(x)`.
fn exp_arg(e: &Expr) -> Option<&Expr> {
    arg1(e, Func::Exp)
}

fn sub1(x: &Expr) -> Expr {
    crate::ops::sub(x, &Expr::one())
}

fn add1(x: &Expr) -> Expr {
    crate::ops::add2(x, &Expr::one())
}

fn over_i_pi(x: &Expr) -> Expr {
    crate::ops::div(x, &crate::ops::mul2(&Expr::i(), &Expr::pi()))
}

fn conjuncts(a: &Expr) -> Vec<Expr> {
    match crate::logic::as_logic(a) {
        Some((crate::logic::LogicOp::And, args)) => args.iter().flat_map(conjuncts).collect(),
        _ => vec![a.clone()],
    }
}

fn in_conjuncts(item: &Expr, a: &Expr) -> bool {
    conjuncts(a).iter().any(|c| c == item)
}

// ---------------------------------------------------------------------------
// certified numerical tests
// ---------------------------------------------------------------------------

const LO: u32 = 96;
const HI: u32 = 256;

fn mp(e: &Expr, prec: u32) -> Option<crate::evalf_mp::MpComplex> {
    if !is_number(e) || is_nan(e) {
        return None;
    }
    crate::evalf_mp::eval_mp_default(e, prec)
}

/// The sign of a value computed at two precisions; `None` when they do not
/// agree to a few bits (cancellation or unsupported magnitude).
fn sure_sign(lo: &BigFloat, hi: &BigFloat) -> Option<i32> {
    if lo.is_zero() && hi.is_zero() {
        return Some(0);
    }
    if lo.is_zero() || hi.is_zero() || lo.sign() != hi.sign() {
        return None;
    }
    let d = hi.sub(lo, HI + 32).abs();
    if d.mul_i64(64, HI + 32).cmp(&hi.abs()) != Ordering::Less {
        return None;
    }
    Some(hi.sign())
}

/// Certified signs of the real and imaginary parts of a number.
fn num_signs(e: &Expr) -> Option<(Option<i32>, Option<i32>)> {
    let lo = mp(e, LO)?;
    let hi = mp(e, HI)?;
    Some((sure_sign(&lo.re, &hi.re), sure_sign(&lo.im, &hi.im)))
}

/// `as_real_imag()[1]` is structurally zero.
fn im_structurally_zero(e: &Expr, im: Option<i32>) -> bool {
    im == Some(0) || crate::assumptions::fact(e, "real") == Some(true)
}

/// `_PositivePredicate_number` (`sense = 1`) / `_NegativePredicate_number` (`-1`).
fn order_number(e: &Expr, sense: i32) -> O {
    let (re, im) = num_signs(e)?;
    if im_structurally_zero(e, im) {
        return re.map(|s| s == sense);
    }
    match im {
        Some(0) => re.map(|s| s == sense),
        Some(_) => Some(false),
        None => None,
    }
}

/// The nearest integer of a real number and whether the number provably
/// differs from it (`Some(false)`), provably equals it (`Some(true)`) or
/// neither (`None`). `None` overall for nonreal or non evaluable numbers.
fn round_number(e: &Expr) -> Option<(num_bigint::BigInt, O, bool)> {
    let lo = mp(e, LO)?;
    let hi = mp(e, HI)?;
    let im = sure_sign(&lo.im, &hi.im);
    if matches!(im, Some(s) if s != 0) {
        return Some((num_bigint::BigInt::from(0), Some(false), true));
    }
    let n = hi.re.to_int_round();
    let nb = BigFloat::from_int(&n);
    let dl = lo.re.sub(&nb, HI + 32);
    let dh = hi.re.sub(&nb, HI + 32);
    if matches!(sure_sign(&dl, &dh), Some(s) if s != 0) {
        return Some((n, Some(false), false));
    }
    let diff = crate::ops::sub(e, &Expr::num(crate::number::Number::from_bigint(n.clone())));
    let proven = diff.is_zero()
        || crate::assumptions::is_zero(&diff) == Some(true)
        || crate::simplify::simplify(&diff).is_zero();
    Some((n, if proven { Some(true) } else { None }, false))
}

fn has_float(e: &Expr) -> bool {
    if is_float(e) {
        return true;
    }
    match e.kind() {
        Kind::Add(c, _) | Kind::Mul(c, _) if c.is_float() => true,
        _ => e.args().iter().any(has_float),
    }
}

fn integer_number(e: &Expr) -> O {
    match round_number(e) {
        Some((_, Some(false), _)) => Some(false),
        Some((_, Some(true), _)) => Some(true),
        _ => None,
    }
}

fn even_number(e: &Expr) -> O {
    if let Kind::Num(n) = e.kind() {
        if n.is_float() {
            return if n.is_integer() { None } else { Some(false) };
        }
    }
    match round_number(e) {
        Some((_, Some(false), _)) => Some(false),
        Some((n, Some(true), _)) => Some(n.is_even()),
        _ => None,
    }
}

fn prime_number(e: &Expr) -> O {
    let exact = !has_float(e);
    match round_number(e) {
        Some((_, Some(false), _)) => Some(false),
        Some((n, _, _)) if exact => Some(crate::ntheory::isprime(&n)),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// entry point
// ---------------------------------------------------------------------------

/// `key(*args)._eval_ask(assumptions)`.
pub(super) fn eval(cx: &Cx, p: Pred, args: &[Expr], a: &Expr) -> R {
    if p.is_binrel() {
        return binrel(cx, p, &args[0], &args[1], a);
    }
    if p.is_undefined() {
        return match args {
            [x] => cx.external(p, x, a),
            _ => Ok(None),
        };
    }
    let [e] = args else { return Ok(None) };
    use Pred::*;
    match p {
        Positive => positive(cx, e, a),
        Negative => negative(cx, e, a),
        Nonnegative => nonneg(cx, e, a, Negative),
        Nonpositive => nonneg(cx, e, a, Positive),
        Nonzero => nonzero(cx, e, a),
        Zero => zero(cx, e, a),
        ExtendedPositive => cx.ask(&or(vec![q(Positive, e), q(PositiveInfinite, e)]), a),
        ExtendedNegative => cx.ask(&or(vec![q(Negative, e), q(NegativeInfinite, e)]), a),
        ExtendedNonzero => cx.ask(
            &or(vec![
                q(NegativeInfinite, e),
                q(Negative, e),
                q(Positive, e),
                q(PositiveInfinite, e),
            ]),
            a,
        ),
        ExtendedNonpositive => cx.ask(
            &or(vec![q(NegativeInfinite, e), q(Negative, e), q(Zero, e)]),
            a,
        ),
        ExtendedNonnegative => cx.ask(
            &or(vec![q(Zero, e), q(Positive, e), q(PositiveInfinite, e)]),
            a,
        ),
        Integer => integer(cx, e, a),
        Rational => rational(cx, e, a),
        Irrational => irrational(cx, e, a),
        Real => real(cx, e, a),
        ExtendedReal => extended_real(cx, e, a),
        Hermitian => hermitian(cx, e, a, false),
        Antihermitian => hermitian(cx, e, a, true),
        Complex => complex(cx, e, a),
        Imaginary => imaginary(cx, e, a),
        Algebraic => algebraic(cx, e, a),
        Prime => prime(cx, e, a),
        Composite => composite(cx, e, a),
        Even => even(cx, e, a),
        Odd => odd(cx, e, a),
        Finite => finite(cx, e, a),
        Infinite => infinite(cx, e, a),
        PositiveInfinite => Ok(match cst(e) {
            Some(Constant::Infinity) => Some(true),
            Some(Constant::NegInfinity | Constant::ComplexInfinity) => Some(false),
            _ => None,
        }),
        NegativeInfinite => Ok(match cst(e) {
            Some(Constant::NegInfinity) => Some(true),
            Some(Constant::Infinity | Constant::ComplexInfinity) => Some(false),
            _ => None,
        }),
        Commutative => commutative(cx, e, a),
        IsTrue => is_true(cx, e, a),
        _ if p.is_matrix() => matrix_pred(cx, p, e, a),
        _ => Ok(None),
    }
}

// ---------------------------------------------------------------------------
// order
// ---------------------------------------------------------------------------

fn positive(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    match e.kind() {
        Kind::Add(..) => {
            if is_number(e) {
                return Ok(order_number(e, 1));
            }
            let r = cx.askp(Real, e, a)?;
            if r != Some(true) {
                return Ok(r);
            }
            let args = args_of(e);
            let mut nonneg = 0;
            for x in &args {
                if cx.askp(Positive, x, a)? != Some(true) {
                    if cx.askp(Negative, x, a)? == Some(false) {
                        nonneg += 1;
                    } else {
                        return Ok(None);
                    }
                }
            }
            return Ok(if nonneg < args.len() {
                Some(true)
            } else {
                None
            });
        }
        Kind::Mul(..) => {
            if is_number(e) {
                return Ok(order_number(e, 1));
            }
            let mut result = true;
            for x in &args_of(e) {
                if cx.is(Positive, x, a)? {
                    continue;
                } else if cx.is(Negative, x, a)? {
                    result = !result;
                } else {
                    return Ok(None);
                }
            }
            return Ok(Some(result));
        }
        Kind::Pow(b, x) => {
            if is_number(e) {
                return Ok(order_number(e, 1));
            }
            if cx.is(Positive, b, a)? && cx.is(Real, x, a)? {
                return Ok(Some(true));
            }
            if cx.is(Negative, b, a)? {
                if cx.is(Even, x, a)? {
                    return Ok(Some(true));
                }
                if cx.is(Odd, x, a)? {
                    return Ok(Some(false));
                }
            }
            return Ok(None);
        }
        Kind::Const(Constant::I) => return Ok(Some(false)),
        Kind::Const(Constant::NaN) => return Ok(None),
        Kind::Func(f, args) => match f {
            Func::Exp => {
                let x = &args[0];
                if cx.is(Real, x, a)? {
                    return Ok(Some(true));
                }
                if cx.is(Imaginary, x, a)? {
                    return cx.askp(Even, &over_i_pi(x), a);
                }
                return Ok(None);
            }
            Func::Log => {
                let x = &args[0];
                let r = cx.askp(Real, x, a)?;
                if r != Some(true) {
                    return Ok(r);
                }
                if cx.is(Positive, &sub1(x), a)? {
                    return Ok(Some(true));
                }
                if cx.is(Negative, &sub1(x), a)? {
                    return Ok(Some(false));
                }
                return Ok(None);
            }
            Func::Special(SpecialFn::Factorial) => {
                let x = &args[0];
                if cx.ask(&and(vec![q(Integer, x), q(Positive, x)]), a)? == Some(true) {
                    return Ok(Some(true));
                }
                return Ok(None);
            }
            Func::Abs => return cx.askp(Nonzero, e, a),
            Func::MatrixTrace(m) | Func::MatrixDeterminant(m) => {
                return Ok(if cx.is(PositiveDefinite, &mx(m), a)? {
                    Some(true)
                } else {
                    None
                });
            }
            Func::MatrixElement(m) if args.len() == 2 => {
                return Ok(
                    if args[0] == args[1] && cx.is(PositiveDefinite, &mx(m), a)? {
                        Some(true)
                    } else {
                        None
                    },
                );
            }
            Func::Atan => return cx.askp(Positive, &args[0], a),
            Func::Asin => {
                let x = &args[0];
                if cx.ask(&and(vec![q(Positive, x), q(Nonpositive, &sub1(x))]), a)? == Some(true) {
                    return Ok(Some(true));
                }
                if cx.ask(&and(vec![q(Negative, x), q(Nonnegative, &add1(x))]), a)? == Some(true) {
                    return Ok(Some(false));
                }
                return Ok(None);
            }
            Func::Acos => {
                let x = &args[0];
                if cx.ask(
                    &and(vec![q(Nonpositive, &sub1(x)), q(Nonnegative, &add1(x))]),
                    a,
                )? == Some(true)
                {
                    return Ok(Some(true));
                }
                return Ok(None);
            }
            Func::Acot => return cx.askp(Real, &args[0], a),
            _ => {}
        },
        _ => {}
    }
    if is_nan(e) {
        return Ok(None);
    }
    if let Some(v) = old(e, Positive) {
        return Ok(Some(v));
    }
    Ok(if is_number(e) {
        order_number(e, 1)
    } else {
        None
    })
}

fn negative(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    match e.kind() {
        Kind::Add(..) => {
            if is_number(e) {
                return Ok(order_number(e, -1));
            }
            let r = cx.askp(Real, e, a)?;
            if r != Some(true) {
                return Ok(r);
            }
            let args = args_of(e);
            let mut nonpos = 0;
            for x in &args {
                if cx.askp(Negative, x, a)? != Some(true) {
                    if cx.askp(Positive, x, a)? == Some(false) {
                        nonpos += 1;
                    } else {
                        return Ok(None);
                    }
                }
            }
            return Ok(if nonpos < args.len() {
                Some(true)
            } else {
                None
            });
        }
        Kind::Mul(..) => {
            if is_number(e) {
                return Ok(order_number(e, -1));
            }
            let mut result = false;
            for x in &args_of(e) {
                if cx.is(Negative, x, a)? {
                    result = !result;
                } else if cx.is(Positive, x, a)? {
                } else {
                    return Ok(None);
                }
            }
            return Ok(Some(result));
        }
        Kind::Pow(b, x) => {
            if is_number(e) {
                return Ok(order_number(e, -1));
            }
            if cx.is(Real, b, a)? {
                if cx.is(Positive, b, a)? && cx.is(Real, x, a)? {
                    return Ok(Some(false));
                }
                if cx.is(Even, x, a)? {
                    return Ok(Some(false));
                }
                if cx.is(Odd, x, a)? {
                    return cx.askp(Negative, b, a);
                }
            }
            return Ok(None);
        }
        Kind::Const(Constant::I) => return Ok(Some(false)),
        Kind::Func(Func::Abs, _) => return Ok(Some(false)),
        Kind::Func(Func::Exp, args) => {
            if cx.is(Real, &args[0], a)? {
                return Ok(Some(false));
            }
        }
        _ => {}
    }
    if let Some(v) = old(e, Negative) {
        return Ok(Some(v));
    }
    Ok(if is_number(e) {
        order_number(e, -1)
    } else {
        None
    })
}

/// Nonnegative (`opposite = Negative`) / Nonpositive (`opposite = Positive`).
fn nonneg(cx: &Cx, e: &Expr, a: &Expr, opposite: Pred) -> R {
    let me = if opposite == Pred::Negative {
        Pred::Nonnegative
    } else {
        Pred::Nonpositive
    };
    if let Some(v) = old(e, me) {
        return Ok(Some(v));
    }
    if is_number(e) {
        let sense = if opposite == Pred::Negative { -1 } else { 1 };
        let notopp = fnot(order_number(e, sense));
        if notopp == Some(true) {
            return cx.askp(Pred::Real, e, a);
        }
        return Ok(notopp);
    }
    Ok(None)
}

fn nonzero(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    match e.kind() {
        Kind::Add(..) => {
            let args = args_of(e);
            let mut allpos = true;
            for x in &args {
                if !cx.is(Positive, x, a)? {
                    allpos = false;
                    break;
                }
            }
            if allpos {
                return Ok(Some(true));
            }
            for x in &args {
                if !cx.is(Negative, x, a)? {
                    return Ok(None);
                }
            }
            return Ok(Some(true));
        }
        Kind::Mul(..) => {
            for x in &args_of(e) {
                let r = cx.askp(Nonzero, x, a)?;
                if r == Some(true) {
                    continue;
                }
                return Ok(r);
            }
            return Ok(Some(true));
        }
        Kind::Pow(b, _) => return cx.askp(Nonzero, b, a),
        Kind::Func(Func::Abs, args) => return cx.askp(Nonzero, &args[0], a),
        _ if is_nan(e) => return Ok(None),
        _ => {}
    }
    if let Some(v) = old(e, Nonzero) {
        return Ok(Some(v));
    }
    if cx.askp(Real, e, &tru())? == Some(false) {
        return Ok(Some(false));
    }
    if is_number(e) {
        let Some((re, im)) = num_signs(e) else {
            return Ok(None);
        };
        let nz = |s: Option<i32>| s.map(|s| s != 0);
        return Ok(match (nz(re), nz(im)) {
            (Some(true), _) | (_, Some(true)) => Some(true),
            (Some(false), Some(false)) => Some(false),
            _ => None,
        });
    }
    Ok(None)
}

fn zero(cx: &Cx, e: &Expr, a: &Expr) -> R {
    if let Kind::Mul(..) = e.kind() {
        return ask_any(
            cx,
            &args_of(e)
                .iter()
                .map(|x| q(Pred::Zero, x))
                .collect::<Vec<_>>(),
            a,
        );
    }
    if let Some(v) = old(e, Pred::Zero) {
        return Ok(Some(v));
    }
    let nz = fnot(cx.askp(Pred::Nonzero, e, a)?);
    if nz == Some(false) {
        return Ok(Some(false));
    }
    let r = cx.askp(Pred::Real, e, a)?;
    Ok(match (nz, r) {
        (_, Some(false)) => Some(false),
        (Some(true), Some(true)) => Some(true),
        _ => None,
    })
}

// ---------------------------------------------------------------------------
// sets
// ---------------------------------------------------------------------------

fn trig_rational_like(cx: &Cx, p: Pred, e: &Expr, a: &Expr) -> Option<R> {
    let Kind::Func(f, args) = e.kind() else {
        return None;
    };
    let x = args.first()?;
    Some((|| match f {
        Func::Asin | Func::Atan | Func::Cos | Func::Sin | Func::Tan | Func::Exp => {
            if cx.is(p, x, a)? {
                return cx.ask(&not(q(Pred::Nonzero, x)), a);
            }
            Ok(None)
        }
        Func::Acot | Func::Cot => Ok(if cx.is(p, x, a)? { Some(false) } else { None }),
        Func::Acos | Func::Log => {
            if cx.is(p, x, a)? {
                return cx.ask(&not(q(Pred::Nonzero, &sub1(x))), a);
            }
            Ok(None)
        }
        _ => unreachable!(),
    })())
    .filter(|_| {
        matches!(
            f,
            Func::Asin
                | Func::Atan
                | Func::Cos
                | Func::Sin
                | Func::Tan
                | Func::Exp
                | Func::Acot
                | Func::Cot
                | Func::Acos
                | Func::Log
        )
    })
}

fn matrix_scalar_parent(e: &Expr) -> Option<MatrixExpr> {
    match e.kind() {
        Kind::Func(
            Func::MatrixTrace(m) | Func::MatrixDeterminant(m) | Func::MatrixElement(m),
            _,
        ) => Some((**m).clone()),
        _ => None,
    }
}

fn integer(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_int(e) {
        return Ok(Some(true));
    }
    if is_rat(e)
        || matches!(
            cst(e),
            Some(
                Constant::E
                    | Constant::I
                    | Constant::Infinity
                    | Constant::NegInfinity
                    | Constant::Pi
            )
        )
    {
        return Ok(Some(false));
    }
    match e.kind() {
        Kind::Add(..) => {
            if is_number(e) {
                return Ok(integer_number(e));
            }
            return closed_group(cx, Integer, &args_of(e), a);
        }
        Kind::Pow(b, x) => {
            if is_number(e) {
                return Ok(integer_number(e));
            }
            if ask_all(cx, &[not(q(Zero, b)), q(Finite, b), q(Zero, x)], a)? == Some(true) {
                return Ok(Some(true));
            }
            if ask_all(cx, &[q(Integer, b), q(Integer, x)], a)? == Some(true)
                && ask_any(
                    cx,
                    &[
                        q(Positive, x),
                        and(vec![q(Nonnegative, x), not(q(Zero, b))]),
                        q(Zero, &sub1(b)),
                        q(Zero, &add1(b)),
                    ],
                    a,
                )? == Some(true)
            {
                return Ok(Some(true));
            }
            return Ok(None);
        }
        Kind::Mul(..) => {
            if is_number(e) {
                return Ok(integer_number(e));
            }
            let mut output = true;
            for x in &args_of(e) {
                if !cx.is(Integer, x, a)? {
                    if is_rat(x) {
                        let den = x.as_number().map(|n| n.denom());
                        if den.as_ref().is_some_and(|d| d.to_i64() == Some(2)) {
                            return cx.askp(Even, &crate::ops::mul2(&Expr::int(2), e), a);
                        }
                        return Ok(None);
                    } else if cx.is(Irrational, x, a)? {
                        if output {
                            output = false;
                        } else {
                            return Ok(None);
                        }
                    } else {
                        return Ok(None);
                    }
                }
            }
            return Ok(Some(output));
        }
        Kind::Func(Func::Abs, args) => {
            return Ok(if cx.is(Integer, &args[0], a)? {
                Some(true)
            } else {
                None
            })
        }
        _ => {}
    }
    if let Some(m) = matrix_scalar_parent(e) {
        return cx.askp(IntegerElements, &mx(&m), a);
    }
    Ok(old(e, Integer))
}

fn rational(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_rat(e) {
        return Ok(Some(true));
    }
    if is_float(e) {
        return Ok(None);
    }
    if matches!(
        cst(e),
        Some(Constant::E | Constant::I | Constant::Infinity | Constant::NegInfinity | Constant::Pi)
    ) {
        return Ok(Some(false));
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) => {
            if is_number(e) {
                if let Some((_, im)) = num_signs(e) {
                    if matches!(im, Some(s) if s != 0) && !im_structurally_zero(e, im) {
                        return Ok(Some(false));
                    }
                }
            }
            return closed_group(cx, Rational, &args_of(e), a);
        }
        Kind::Pow(b, x) => {
            let is_exp_integer = cx.askp(Integer, x, a)?;
            if is_exp_integer == Some(true) {
                if cx.is(Rational, b, a)? {
                    let is_base_zero = cx.askp(Zero, b, a)?;
                    if is_base_zero == Some(false) {
                        return Ok(Some(true));
                    }
                    if is_base_zero == Some(true) && cx.is(Positive, x, &tru())? {
                        return Ok(Some(true));
                    }
                }
                if cx.askp(Algebraic, b, a)? == Some(false) {
                    return cx.askp(Zero, x, a);
                }
                if cx.is(Irrational, b, a)?
                    && cx.ask(
                        &applied_unchecked(Eq, vec![x.clone(), Expr::minus_one()]),
                        &tru(),
                    )? == Some(true)
                {
                    return Ok(Some(false));
                }
                return Ok(None);
            } else if cx.is(Rational, x, a)? {
                if cx.is(Prime, b, a)? && is_exp_integer == Some(false) {
                    return Ok(Some(false));
                }
                if cx.is(Zero, b, &tru())? && cx.is(Positive, x, &tru())? {
                    return Ok(Some(true));
                }
                if cx.ask(&applied_unchecked(Eq, vec![b.clone(), Expr::one()]), &tru())?
                    == Some(true)
                {
                    return Ok(Some(true));
                }
            }
            return Ok(None);
        }
        _ => {}
    }
    if let Some(r) = trig_rational_like(cx, Rational, e, a) {
        return r;
    }
    Ok(old(e, Rational))
}

fn irrational(cx: &Cx, e: &Expr, a: &Expr) -> R {
    if let Some(v) = old(e, Pred::Irrational) {
        return Ok(Some(v));
    }
    let r = cx.askp(Pred::Real, e, a)?;
    if r == Some(true) {
        let rat = cx.askp(Pred::Rational, e, a)?;
        return Ok(fnot(rat));
    }
    Ok(r)
}

/// `_RealPredicate_number`.
fn real_number(e: &Expr) -> O {
    let (_, im) = num_signs(e)?;
    im.map(|s| s == 0)
}

fn real(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_rat(e) || is_float(e) || matches!(cst(e), Some(Constant::E | Constant::Pi)) {
        return Ok(Some(true));
    }
    if matches!(
        cst(e),
        Some(Constant::I | Constant::Infinity | Constant::NegInfinity)
    ) {
        return Ok(Some(false));
    }
    match e.kind() {
        Kind::Func(Func::Abs | Func::Re | Func::Im, _) => return Ok(Some(true)),
        Kind::Add(..) => {
            if is_number(e) {
                return Ok(real_number(e));
            }
            return closed_group(cx, Real, &args_of(e), a);
        }
        Kind::Mul(..) => {
            if is_number(e) {
                return Ok(real_number(e));
            }
            let mut result = true;
            for x in &args_of(e) {
                if cx.is(Real, x, a)? {
                } else if cx.is(Imaginary, x, a)? {
                    result = !result;
                } else {
                    return Ok(None);
                }
            }
            return Ok(Some(result));
        }
        Kind::Pow(b, x) => {
            if is_number(e) {
                return Ok(real_number(e));
            }
            if let Some(bx) = exp_arg(b) {
                if cx.is(Imaginary, bx, a)? && cx.is(Imaginary, x, a)? {
                    return Ok(Some(true));
                }
                let i = over_i_pi(bx);
                if cx.is(Integer, &crate::ops::mul2(&Expr::int(2), &i), a)? {
                    let p = crate::ops::pow(&crate::ops::pow(&Expr::minus_one(), &i), x);
                    return cx.askp(Real, &p, a);
                }
                return Ok(None);
            }
            if cx.is(Imaginary, b, a)? && cx.is(Integer, x, a)? {
                let odd = cx.askp(Odd, x, a)?;
                return Ok(fnot(odd));
            }
            if cx.is(Imaginary, x, a)? {
                let imlog = cx.askp(
                    Imaginary,
                    &crate::functions::func(Func::Log, vec![b.clone()]),
                    a,
                )?;
                if imlog.is_some() {
                    return Ok(imlog);
                }
            }
            if cx.is(Real, b, a)? && cx.is(Real, x, a)? {
                if cx.askp(Zero, b, a)? != Some(false) {
                    if cx.is(Positive, x, a)? {
                        return Ok(Some(true));
                    }
                    return Ok(None);
                }
                if let Some(n) = x.as_number().filter(|_| is_rat(x)) {
                    if cx.is(Even, &Expr::num(n.denom()), a)? {
                        return cx.askp(Positive, b, a);
                    }
                }
                if cx.is(Integer, x, a)? {
                    return Ok(Some(true));
                }
                if cx.is(Positive, b, a)? {
                    return Ok(Some(true));
                }
            }
            return Ok(None);
        }
        Kind::Func(Func::Cos | Func::Sin, args) => {
            return Ok(if cx.is(Real, &args[0], a)? {
                Some(true)
            } else {
                None
            })
        }
        Kind::Func(Func::Exp, args) => {
            let x = &args[0];
            return cx.ask(&or(vec![q(Integer, &over_i_pi(x)), q(Real, x)]), a);
        }
        Kind::Func(Func::Log, args) => return cx.askp(Positive, &args[0], a),
        _ => {}
    }
    if let Some(m) = matrix_scalar_parent(e) {
        return cx.askp(RealElements, &mx(&m), a);
    }
    Ok(old(e, Real))
}

fn extended_real(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if matches!(cst(e), Some(Constant::Infinity | Constant::NegInfinity)) {
        return Ok(Some(true));
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) => closed_group(cx, ExtendedReal, &args_of(e), a),
        Kind::Pow(b, x) => closed_group(cx, ExtendedReal, &[b.clone(), x.clone()], a),
        _ => cx.ask(
            &or(vec![
                q(NegativeInfinite, e),
                q(Negative, e),
                q(Zero, e),
                q(Positive, e),
                q(PositiveInfinite, e),
            ]),
            a,
        ),
    }
}

/// Elementwise test of an explicit matrix against `conjugate` (or `-conjugate`).
fn explicit_hermitian(m: &crate::matrix::Matrix, anti: bool) -> O {
    let mut ret = Some(true);
    for i in 0..m.rows() {
        for j in i..m.cols() {
            if j >= m.rows() || i >= m.cols() {
                return Some(false);
            }
            let c = crate::functions::func(Func::Conjugate, vec![m.get(j, i).clone()]);
            let c = if anti { crate::ops::neg(&c) } else { c };
            match crate::relational::eq(m.get(i, j), &c).kind() {
                Kind::Bool(true) => {}
                Kind::Bool(false) => return Some(false),
                _ => ret = None,
            }
        }
    }
    ret
}

/// Hermitian (`anti = false`) and antihermitian handlers.
fn hermitian(cx: &Cx, e: &Expr, a: &Expr, anti: bool) -> R {
    use Pred::*;
    let me = if anti { Antihermitian } else { Hermitian };
    let num = is_number(e);
    match e.kind() {
        Kind::Add(..) if !num => return closed_group(cx, me, &args_of(e), a),
        Kind::Mul(..) if !num => {
            let mut nccount = 0;
            let mut result = !anti;
            for x in &args_of(e) {
                if cx.is(Antihermitian, x, a)? {
                    result = !result;
                } else if !cx.is(Hermitian, x, a)? {
                    return Ok(None);
                }
                if cx.ask(&not(q(Commutative, x)), a)? == Some(true) {
                    nccount += 1;
                    if nccount > 1 {
                        return Ok(None);
                    }
                }
            }
            return Ok(Some(result));
        }
        Kind::Pow(b, x) if !num => {
            if !anti {
                if cx.is(Hermitian, b, a)? && cx.is(Integer, x, a)? {
                    return Ok(Some(true));
                }
            } else if cx.is(Hermitian, b, a)? {
                if cx.is(Integer, x, a)? {
                    return Ok(Some(false));
                }
            } else if cx.is(Antihermitian, b, a)? {
                if cx.is(Even, x, a)? {
                    return Ok(Some(false));
                } else if cx.is(Odd, x, a)? {
                    return Ok(Some(true));
                }
            }
        }
        Kind::Func(Func::Cos | Func::Sin | Func::Exp, args) if !anti => {
            if cx.is(Hermitian, &args[0], a)? {
                return Ok(Some(true));
            }
        }
        Kind::Func(Func::MatrixProxy(m), _) => {
            if let MatrixKind::Explicit(mat) = m.kind() {
                // MatrixBase handler, then the object handler returns None.
                return Ok(explicit_hermitian(mat, anti).filter(|_| true));
            }
        }
        _ => {}
    }
    // object handler
    if anti {
        if cx.is(Zero, e, a)? {
            return Ok(Some(true));
        }
        return cx.askp(Imaginary, e, a);
    }
    cx.askp(Real, e, a)
}

fn complex(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_nan(e) {
        return Ok(None);
    }
    if matches!(cst(e), Some(Constant::Infinity | Constant::NegInfinity)) {
        return Ok(Some(false));
    }
    if is_number_class(e) || is_number_symbol(e) || matches!(cst(e), Some(Constant::I)) {
        return Ok(Some(true));
    }
    match e.kind() {
        Kind::Func(
            Func::Abs | Func::Cos | Func::Exp | Func::Im | Func::Log | Func::Re | Func::Sin,
            _,
        ) => return Ok(Some(true)),
        Kind::Add(..) | Kind::Mul(..) => return closed_group(cx, Complex, &args_of(e), a),
        Kind::Pow(b, x) => return closed_group(cx, Complex, &[b.clone(), x.clone()], a),
        _ => {}
    }
    if let Some(m) = matrix_scalar_parent(e) {
        return cx.askp(ComplexElements, &mx(&m), a);
    }
    Ok(old(e, Complex))
}

fn imaginary_number(e: &Expr) -> O {
    let (re, _) = num_signs(e)?;
    re.map(|s| s == 0)
}

fn imaginary(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if matches!(cst(e), Some(Constant::I)) {
        return Ok(Some(true));
    }
    if is_nan(e) {
        return Ok(None);
    }
    if is_number_class(e) || is_number_symbol(e) {
        // `not expr.as_real_imag()[1] == 0`
        return Ok(Some(false));
    }
    match e.kind() {
        Kind::Add(..) => {
            if is_number(e) {
                return Ok(imaginary_number(e));
            }
            let args = args_of(e);
            let mut reals = 0;
            for x in &args {
                if cx.is(Imaginary, x, a)? {
                } else if cx.is(Real, x, a)? {
                    reals += 1;
                } else {
                    return Ok(None);
                }
            }
            if reals == 0 {
                return Ok(Some(true));
            }
            if reals == 1 || reals == args.len() {
                return Ok(Some(false));
            }
            return Ok(None);
        }
        Kind::Mul(..) => {
            if is_number(e) {
                return Ok(imaginary_number(e));
            }
            let mut result = false;
            for x in &args_of(e) {
                if cx.is(Imaginary, x, a)? {
                    result = !result;
                } else if !cx.is(Real, x, a)? {
                    return Ok(None);
                }
            }
            return Ok(Some(result));
        }
        Kind::Pow(b, x) => {
            if is_number(e) {
                return Ok(imaginary_number(e));
            }
            if let Some(bx) = exp_arg(b) {
                if cx.is(Imaginary, bx, a)? {
                    if cx.is(Imaginary, x, a)? {
                        return Ok(Some(false));
                    }
                    let i = over_i_pi(bx);
                    if cx.is(Integer, &crate::ops::mul2(&Expr::int(2), &i), a)? {
                        let p = crate::ops::pow(&crate::ops::pow(&Expr::minus_one(), &i), x);
                        return cx.askp(Imaginary, &p, a);
                    }
                }
            }
            if cx.is(Imaginary, b, a)? && cx.is(Integer, x, a)? {
                return cx.askp(Odd, x, a);
            }
            if cx.is(Imaginary, x, a)? {
                let imlog = cx.askp(
                    Imaginary,
                    &crate::functions::func(Func::Log, vec![b.clone()]),
                    a,
                )?;
                if imlog.is_some() {
                    return Ok(Some(false));
                }
            }
            if cx.ask(&and(vec![q(Real, b), q(Real, x)]), a)? == Some(true) {
                if cx.is(Positive, b, a)? {
                    return Ok(Some(false));
                }
                let rat = cx.askp(Rational, x, a)?;
                if rat != Some(true) {
                    return Ok(rat);
                }
                if cx.is(Integer, x, a)? {
                    return Ok(Some(false));
                }
                let half = cx.askp(Integer, &crate::ops::mul2(&Expr::int(2), x), a)?;
                if half == Some(true) {
                    return cx.askp(Negative, b, a);
                }
                return Ok(half);
            }
            return Ok(None);
        }
        Kind::Func(Func::Log, args) => {
            let x = &args[0];
            if cx.is(Real, x, a)? {
                return Ok(if cx.is(Positive, x, a)? {
                    Some(false)
                } else {
                    None
                });
            }
            if let Some(ex) = exp_arg(x) {
                if *ex == Expr::i() || *ex == crate::ops::neg(&Expr::i()) {
                    return Ok(Some(true));
                }
            }
            return Ok(if cx.askp(Imaginary, x, a)? == Some(false) {
                Some(false)
            } else {
                None
            });
        }
        Kind::Func(Func::Exp, args) => {
            let i = over_i_pi(&args[0]);
            return cx.ask(
                &and(vec![
                    q(Integer, &crate::ops::mul2(&Expr::int(2), &i)),
                    not(q(Integer, &i)),
                ]),
                a,
            );
        }
        _ => {}
    }
    Ok(old(e, Imaginary))
}

fn algebraic(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_float(e) || matches!(cst(e), Some(Constant::I)) || is_rat(e) {
        return Ok(Some(true));
    }
    if matches!(
        cst(e),
        Some(
            Constant::ComplexInfinity
                | Constant::E
                | Constant::Infinity
                | Constant::NegInfinity
                | Constant::Pi
        )
    ) {
        return Ok(Some(false));
    }
    match e.kind() {
        Kind::Add(..) | Kind::Mul(..) => return closed_group(cx, Algebraic, &args_of(e), a),
        Kind::Pow(b, x) => {
            if matches!(cst(b), Some(Constant::Pi)) {
                if cx.is(Integer, x, a)? && cx.is(Positive, x, a)? {
                    return Ok(Some(false));
                }
                return Ok(None);
            }
            let exp_rational = cx.askp(Rational, x, a)?;
            let base_algebraic = cx.askp(Algebraic, b, a)?;
            let exp_algebraic = cx.askp(Algebraic, x, a)?;
            if base_algebraic == Some(true) && exp_algebraic == Some(true) {
                if exp_rational == Some(true) {
                    return Ok(Some(true));
                }
                let ne = and(vec![
                    applied_unchecked(Ne, vec![b.clone(), Expr::zero()]),
                    applied_unchecked(Ne, vec![b.clone(), Expr::one()]),
                ]);
                if cx.ask(&ne, &tru())? == Some(true) && exp_rational == Some(false) {
                    return Ok(Some(false));
                }
            }
            return Ok(None);
        }
        _ => {}
    }
    if let Some(r) = trig_rational_like(cx, Algebraic, e, a) {
        return r;
    }
    Ok(None)
}

// ---------------------------------------------------------------------------
// ntheory
// ---------------------------------------------------------------------------

fn prime(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_int(e) {
        return Ok(Some(
            e.as_number()
                .and_then(|n| n.to_bigint())
                .is_some_and(|n| crate::ntheory::isprime(&n)),
        ));
    }
    if is_rat(e)
        || matches!(
            cst(e),
            Some(Constant::Infinity | Constant::NegInfinity | Constant::I)
        )
    {
        return Ok(Some(false));
    }
    if is_float(e) || is_number_symbol(e) {
        return Ok(prime_number(e));
    }
    if is_nan(e) {
        return Ok(None);
    }
    match e.kind() {
        Kind::Mul(..) => {
            if is_number(e) {
                return Ok(prime_number(e));
            }
            let args = args_of(e);
            for x in &args {
                if !cx.is(Integer, x, a)? {
                    return Ok(None);
                }
            }
            for x in &args {
                if is_number(x) && crate::assumptions::fact(x, "composite") == Some(true) {
                    return Ok(Some(false));
                }
            }
            return Ok(None);
        }
        Kind::Pow(b, x) => {
            if is_number(e) {
                return Ok(prime_number(e));
            }
            if cx.is(Integer, x, a)? && cx.is(Integer, b, a)? {
                let prime_base = cx.askp(Prime, b, a)?;
                if prime_base == Some(false) {
                    return Ok(Some(false));
                }
                let is_exp_one = cx.ask(&applied_unchecked(Eq, vec![x.clone(), Expr::one()]), a)?;
                if is_exp_one == Some(false) {
                    return Ok(Some(false));
                }
                if prime_base == Some(true) && is_exp_one == Some(true) {
                    return Ok(Some(true));
                }
            }
            return Ok(None);
        }
        _ => {}
    }
    if let Some(v) = old(e, Prime) {
        return Ok(Some(v));
    }
    Ok(if is_number(e) { prime_number(e) } else { None })
}

fn composite(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if let Some(v) = old(e, Composite) {
        return Ok(Some(v));
    }
    let pos = cx.askp(Positive, e, a)?;
    if pos != Some(true) {
        return Ok(pos);
    }
    let int = cx.askp(Integer, e, a)?;
    if int != Some(true) {
        return Ok(int);
    }
    let pr = cx.askp(Prime, e, a)?;
    if pr.is_none() {
        return Ok(None);
    }
    let is_one = cx.ask(&applied_unchecked(Eq, vec![e.clone(), Expr::one()]), a)?;
    match is_one {
        Some(true) => Ok(Some(false)),
        None => Ok(None),
        Some(false) => Ok(fnot(pr)),
    }
}

fn even(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_int(e) {
        return Ok(e
            .as_number()
            .and_then(|n| n.to_bigint())
            .map(|n| n.is_even()));
    }
    if is_rat(e)
        || matches!(
            cst(e),
            Some(Constant::Infinity | Constant::NegInfinity | Constant::I)
        )
    {
        return Ok(Some(false));
    }
    if is_number_symbol(e) {
        return Ok(even_number(e));
    }
    if is_nan(e) {
        return Ok(None);
    }
    match e.kind() {
        Kind::Mul(..) => {
            if is_number(e) {
                return Ok(even_number(e));
            }
            let args = args_of(e);
            let (mut ev, mut odd, mut irrational) = (false, 0usize, false);
            let mut acc = Expr::one();
            for x in &args {
                if cx.is(Integer, x, a)? {
                    if cx.is(Even, x, a)? {
                        ev = true;
                    } else if cx.is(Odd, x, a)? {
                        odd += 1;
                    } else if !ev && !acc.is_one() && cx.is(Odd, &crate::ops::add2(&acc, x), a)? {
                        ev = true;
                    }
                } else if cx.is(Irrational, x, a)? {
                    if irrational {
                        return Ok(None);
                    }
                    irrational = true;
                } else {
                    return Ok(None);
                }
                acc = x.clone();
            }
            if irrational {
                return Ok(Some(false));
            }
            if ev {
                return Ok(Some(true));
            }
            if odd == args.len() {
                return Ok(Some(false));
            }
            return Ok(None);
        }
        Kind::Add(..) => {
            if is_number(e) {
                return Ok(even_number(e));
            }
            let mut result = true;
            for x in &args_of(e) {
                if cx.is(Even, x, a)? {
                } else if cx.is(Odd, x, a)? {
                    result = !result;
                } else {
                    return Ok(None);
                }
            }
            return Ok(Some(result));
        }
        Kind::Pow(b, x) => {
            if is_number(e) {
                return Ok(even_number(e));
            }
            if cx.is(Integer, x, a)? {
                if cx.is(Positive, x, a)? {
                    return cx.askp(Even, b, a);
                } else if cx.ask(&and(vec![not(q(Negative, x)), q(Odd, b)]), a)? == Some(true) {
                    return Ok(Some(false));
                } else if b.is_minus_one() {
                    return Ok(Some(false));
                }
            }
            return Ok(None);
        }
        Kind::Func(Func::Abs | Func::Re, args) => {
            if cx.is(Real, &args[0], a)? {
                return cx.askp(Even, &args[0], a);
            }
            return Ok(None);
        }
        Kind::Func(Func::Im, args) => {
            return Ok(if cx.is(Real, &args[0], a)? {
                Some(true)
            } else {
                None
            })
        }
        _ => {}
    }
    if let Some(v) = old(e, Even) {
        return Ok(Some(v));
    }
    Ok(if is_number(e) { even_number(e) } else { None })
}

fn odd(cx: &Cx, e: &Expr, a: &Expr) -> R {
    if let Some(v) = old(e, Pred::Odd) {
        return Ok(Some(v));
    }
    let int = cx.askp(Pred::Integer, e, a)?;
    if int == Some(true) {
        let ev = cx.askp(Pred::Even, e, a)?;
        return Ok(fnot(ev));
    }
    Ok(int)
}

// ---------------------------------------------------------------------------
// calculus
// ---------------------------------------------------------------------------

fn finite(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if is_nan(e) {
        return Ok(None);
    }
    if matches!(
        cst(e),
        Some(Constant::ComplexInfinity | Constant::Infinity | Constant::NegInfinity)
    ) {
        return Ok(Some(false));
    }
    if is_number_class(e) || matches!(cst(e), Some(Constant::Pi | Constant::E | Constant::I)) {
        return Ok(Some(true));
    }
    match e.kind() {
        Kind::Sym(_) => {
            if let Some(v) = crate::assumptions::fact(e, "finite") {
                return Ok(Some(v));
            }
            return Ok(if in_conjuncts(&q(Finite, e), a) {
                Some(true)
            } else {
                None
            });
        }
        Kind::Add(..) => {
            // sign: None = the -1 sentinel, Some(x) = the fuzzy value x
            let mut sign: Option<O> = None;
            let mut result = Some(true);
            for x in &args_of(e) {
                let bounded = cx.askp(Finite, x, a)?;
                if bounded == Some(true) {
                    continue;
                }
                let s = cx.askp(ExtendedPositive, x, a)?;
                let differs = matches!(sign, Some(prev) if prev != s);
                let none_case = s.is_none() && (bounded.is_none() || matches!(sign, Some(None)));
                if differs || none_case {
                    return Ok(None);
                }
                sign = Some(s);
                if result != Some(false) {
                    result = bounded;
                }
            }
            return Ok(result);
        }
        Kind::Mul(..) => {
            let mut result = Some(true);
            let mut possible_zero = false;
            for x in &args_of(e) {
                let bounded = cx.askp(Finite, x, a)?;
                match bounded {
                    Some(true) => {
                        if cx.askp(Zero, x, a)? != Some(false) {
                            if result == Some(false) {
                                return Ok(None);
                            }
                            possible_zero = true;
                        }
                    }
                    None => {
                        if result.is_none() {
                            return Ok(None);
                        }
                        if cx.askp(ExtendedNonzero, x, a)?.is_none() {
                            return Ok(None);
                        }
                        if result != Some(false) {
                            result = None;
                        }
                    }
                    Some(false) => {
                        if possible_zero {
                            return Ok(None);
                        }
                        result = Some(false);
                    }
                }
            }
            return Ok(result);
        }
        Kind::Pow(b, x) => {
            let base_bounded = cx.askp(Finite, b, a)?;
            let exp_bounded = cx.askp(Finite, x, a)?;
            if base_bounded.is_none() && exp_bounded.is_none() {
                return Ok(None);
            }
            if base_bounded == Some(false) && cx.is(ExtendedNonzero, x, a)? {
                return Ok(Some(false));
            }
            if base_bounded == Some(true) && exp_bounded == Some(true) {
                let is_base_zero = cx.askp(Zero, b, a)?;
                let is_exp_negative = cx.askp(Negative, x, a)?;
                if is_base_zero == Some(true) && is_exp_negative == Some(true) {
                    return Ok(Some(false));
                }
                if is_base_zero != Some(false) && is_exp_negative != Some(false) {
                    return Ok(None);
                }
                return Ok(Some(true));
            }
            let absb = crate::functions::func(Func::Abs, vec![b.clone()]);
            let le1 = crate::relational::le(&absb, &Expr::one()) == tru();
            let ge1 = crate::relational::ge(&absb, &Expr::one()) == tru();
            if le1 && cx.is(ExtendedPositive, x, a)? {
                return Ok(Some(true));
            }
            if ge1 && cx.is(ExtendedNegative, x, a)? {
                return Ok(Some(true));
            }
            if ge1 && exp_bounded == Some(false) {
                return Ok(Some(false));
            }
            return Ok(None);
        }
        Kind::Func(Func::Exp, args) => return cx.askp(Finite, &args[0], a),
        Kind::Func(Func::Log, args) => {
            if cx.is(Infinite, &args[0], a)? {
                return Ok(Some(false));
            }
            return cx.ask(&not(q(Zero, &args[0])), a);
        }
        Kind::Func(Func::Cos | Func::Sin | Func::Sign, _) => return Ok(Some(true)),
        _ => {}
    }
    Ok(None)
}

fn infinite(cx: &Cx, e: &Expr, a: &Expr) -> R {
    if !is_expr(e) {
        return Ok(None);
    }
    Ok(fnot(finite(cx, e, a)?))
}

// ---------------------------------------------------------------------------
// common
// ---------------------------------------------------------------------------

fn commutative(cx: &Cx, e: &Expr, a: &Expr) -> R {
    if let Kind::Sym(_) = e.kind() {
        let c = crate::assumptions::fact(e, "commutative");
        return Ok(match c {
            Some(c) => Some(c && !in_conjuncts(&not(q(Pred::Commutative, e)), a)),
            None => {
                if in_conjuncts(&q(Pred::Commutative, e), a) {
                    Some(true)
                } else if in_conjuncts(&not(q(Pred::Commutative, e)), a) {
                    Some(false)
                } else {
                    Some(true)
                }
            }
        });
    }
    if is_number_class(e) {
        return Ok(Some(true));
    }
    let args: Vec<Expr> = match e.kind() {
        Kind::Func(Func::AppliedPredicate, args) => args[1..].to_vec(),
        Kind::Func(Func::Predicate(_), _) => Vec::new(),
        _ => e.args(),
    };
    for x in &args {
        if !cx.is(Pred::Commutative, x, a)? {
            return Ok(Some(false));
        }
    }
    Ok(Some(true))
}

fn is_true(cx: &Cx, e: &Expr, a: &Expr) -> R {
    use crate::logic::LogicOp;
    match e.kind() {
        Kind::Bool(b) => Ok(Some(*b)),
        Kind::Func(Func::AppliedPredicate, _) => cx.ask(e, a),
        Kind::Func(Func::Logic(op), args) => match op {
            LogicOp::Not => {
                if matches!(args[0].kind(), Kind::Sym(_)) {
                    return Ok(None);
                }
                Ok(fnot(cx.ask(&args[0], a)?))
            }
            LogicOp::Or => {
                let mut result = Some(false);
                for x in args {
                    match cx.ask(x, a)? {
                        Some(true) => return Ok(Some(true)),
                        None => result = None,
                        _ => {}
                    }
                }
                Ok(result)
            }
            LogicOp::And => {
                let mut result = Some(true);
                for x in args {
                    match cx.ask(x, a)? {
                        Some(false) => return Ok(Some(false)),
                        None => result = None,
                        _ => {}
                    }
                }
                Ok(result)
            }
            LogicOp::Implies => cx.ask(&or(vec![crate::logic::not(&args[0]), args[1].clone()]), a),
            LogicOp::Equivalent if args.len() == 2 => {
                let Some(pt) = cx.ask(&args[0], a)? else {
                    return Ok(None);
                };
                let Some(qt) = cx.ask(&args[1], a)? else {
                    return Ok(None);
                };
                Ok(Some(pt == qt))
            }
            _ => Ok(None),
        },
        _ => Ok(None),
    }
}

// ---------------------------------------------------------------------------
// binary relations
// ---------------------------------------------------------------------------

/// A conjunct as a binary relation `(predicate, lhs, rhs, negated)`.
fn as_binrel(c: &Expr) -> Option<(Pred, Expr, Expr, bool)> {
    match c.kind() {
        Kind::Rel(op, l, r) => Some((Pred::from_relop(*op), l.clone(), r.clone(), false)),
        Kind::Func(Func::Logic(crate::logic::LogicOp::Not), args) => {
            let (p, l, r, n) = as_binrel(&args[0])?;
            if matches!(args[0].kind(), Kind::Rel(..)) {
                return None;
            }
            Some((p, l, r, !n))
        }
        _ => {
            let (p, args) = as_applied(c)?;
            if p.is_binrel() && args.len() == 2 {
                Some((p, args[0].clone(), args[1].clone(), false))
            } else {
                None
            }
        }
    }
}

fn binrel(cx: &Cx, p: Pred, l: &Expr, r: &Expr, a: &Expr) -> R {
    let rels: Vec<(Pred, Expr, Expr, bool)> = conjuncts(a).iter().filter_map(as_binrel).collect();
    let has = |p: Pred, x: &Expr, y: &Expr, n: bool| {
        rels.iter()
            .any(|(q, l2, r2, n2)| *q == p && l2 == x && r2 == y && *n2 == n)
    };
    let rp = p.reversed();
    if has(p, l, r, false) || has(rp, r, l, false) {
        return Ok(Some(true));
    }
    let np = p.negated().unwrap_or(p);
    let nrp = rp.negated().unwrap_or(rp);
    if has(np, l, r, false) || has(nrp, r, l, false) || has(p, l, r, true) || has(rp, r, l, true) {
        return Ok(Some(false));
    }
    let v = binrel_eval(cx, p, l, r, a)?;
    if v.is_some() {
        return Ok(v);
    }
    if is_matrix(l) || is_matrix(r) {
        return Ok(None);
    }
    let (ls, rs) = (crate::simplify::simplify(l), crate::simplify::simplify(r));
    if ls == *l && rs == *r {
        return Ok(None);
    }
    binrel_eval(cx, p, &ls, &rs, a)
}

fn binrel_eval(cx: &Cx, p: Pred, l: &Expr, r: &Expr, a: &Expr) -> R {
    match p {
        Pred::Eq => is_eq(cx, l, r, a),
        Pred::Ne => Ok(fnot(is_eq(cx, l, r, a)?)),
        Pred::Ge => is_ge(cx, l, r, a),
        Pred::Le => is_ge(cx, r, l, a),
        Pred::Lt => Ok(fnot(is_ge(cx, l, r, a)?)),
        Pred::Gt => Ok(fnot(is_ge(cx, r, l, a)?)),
        _ => Ok(None),
    }
}

fn is_ge(cx: &Cx, l: &Expr, r: &Expr, a: &Expr) -> R {
    use Pred::*;
    if !is_expr(l) || !is_expr(r) {
        return Err(AskError::Type(
            "Can only compare inequalities with Expr".into(),
        ));
    }
    if is_matrix(l) || is_matrix(r) {
        return Ok(None);
    }
    if let Kind::Bool(b) = crate::relational::ge(l, r).kind() {
        return Ok(Some(*b));
    }
    if cx.is(ExtendedReal, l, a)? && cx.is(ExtendedReal, r, a)? {
        if (cx.is(Infinite, l, a)? && cx.is(ExtendedPositive, l, a)?)
            || (cx.is(Infinite, r, a)? && cx.is(ExtendedNegative, r, a)?)
        {
            return Ok(Some(true));
        }
        let diff = crate::ops::sub(l, r);
        if !is_nan(&diff) {
            let rv = cx.askp(ExtendedNonnegative, &diff, a)?;
            if rv.is_some() {
                return Ok(rv);
            }
        }
    }
    Ok(None)
}

fn fuzzy_xor2(x: O, y: O) -> O {
    Some(x? ^ y?)
}

fn is_eq(cx: &Cx, l: &Expr, r: &Expr, a: &Expr) -> R {
    use Pred::*;
    if l == r {
        return Ok(Some(true));
    }
    if is_matrix(l) || is_matrix(r) {
        return Ok(None);
    }
    let lb = matches!(l.kind(), Kind::Bool(_));
    let rb = matches!(r.kind(), Kind::Bool(_));
    if lb && rb {
        return Ok(Some(false));
    }
    if !(l.is_symbol() || r.is_symbol()) && is_expr(l) != is_expr(r) {
        return Ok(Some(false));
    }
    if !is_expr(l) || !is_expr(r) {
        return Ok(None);
    }
    if let Kind::Bool(b) = crate::relational::eq(l, r).kind() {
        return Ok(Some(*b));
    }
    let li = cx.askp(Infinite, l, a)?;
    let ri = cx.askp(Infinite, r, a)?;
    if li == Some(true) || ri == Some(true) {
        if fuzzy_xor2(li, ri) == Some(true) {
            return Ok(Some(false));
        }
        let lr = cx.askp(ExtendedReal, l, a)?;
        let rr = cx.askp(ExtendedReal, r, a)?;
        if fuzzy_xor2(lr, rr) == Some(true) {
            return Ok(Some(false));
        }
        if lr == Some(true) && rr == Some(true) {
            let lp = cx.askp(ExtendedPositive, l, a)?;
            let rp = cx.askp(ExtendedPositive, r, a)?;
            return Ok(fuzzy_xor2(lp, fnot(rp)));
        }
        return Ok(None);
    }
    let dif = crate::ops::sub(l, r);
    match cx.askp(Zero, &dif, a)? {
        Some(true) => return Ok(Some(true)),
        Some(false) if cx.is(Commutative, &dif, a)? => return Ok(Some(false)),
        _ => {}
    }
    Ok(None)
}

// ---------------------------------------------------------------------------
// matrices
// ---------------------------------------------------------------------------

/// A matrix expression as a scalar `Expr` argument.
pub(crate) fn mx(m: &MatrixExpr) -> Expr {
    Expr::from_kind(Kind::Func(Func::MatrixProxy(Arc::new(m.clone())), vec![]))
}

fn as_mx(e: &Expr) -> Option<&MatrixExpr> {
    match e.kind() {
        Kind::Func(Func::MatrixProxy(m), _) => Some(m),
        _ => None,
    }
}

/// SymPy `MatMul.args`: the scalar coefficient (when not 1) then the factors.
pub(super) fn matmul_args(coeff: &Expr, factors: &[MatrixExpr]) -> Vec<Expr> {
    let mut out = Vec::with_capacity(factors.len() + 1);
    if !coeff.is_one() {
        out.push(coeff.clone());
    }
    out.extend(factors.iter().map(mx));
    out
}

fn is_square(m: &MatrixExpr) -> bool {
    m.rows() == m.cols()
}

fn is_shape(m: &MatrixExpr, n: i64) -> bool {
    *m.rows() == Expr::int(n) && *m.cols() == Expr::int(n)
}

fn empty_or_1x1(m: &MatrixExpr) -> bool {
    is_shape(m, 0) || is_shape(m, 1)
}

/// Unevaluated `MatMul(*factors)`.
fn matmul(factors: &[MatrixExpr]) -> Option<Expr> {
    MatrixExpr::mul(Expr::one(), factors, false)
        .ok()
        .map(|m| mx(&m))
}

fn on_diag(rows: &crate::matrixexpr::SliceSpec, cols: &crate::matrixexpr::SliceSpec) -> bool {
    rows == cols
}

/// `MatPow` handlers sharing the "integer exponent, nonnegative or
/// invertible base" test.
fn matpow_nonneg_or_invertible(cx: &Cx, p: Pred, b: &MatrixExpr, x: &Expr, a: &Expr) -> R {
    if !cx.is(Pred::Integer, x, a)? {
        return Ok(None);
    }
    let nn = cx.ask(&not(q(Pred::Negative, x)), a)?;
    if nn == Some(true) || (nn == Some(false) && cx.is(Pred::Invertible, &mx(b), a)?) {
        return cx.askp(p, &mx(b), a);
    }
    Ok(None)
}

fn all_matrices(cx: &Cx, p: Pred, ms: &[MatrixExpr], a: &Expr) -> Result<bool, AskError> {
    for m in ms {
        if !cx.is(p, &mx(m), a)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn any_false(cx: &Cx, p: Pred, ms: &[MatrixExpr], a: &Expr) -> Result<bool, AskError> {
    for m in ms {
        if cx.askp(p, &mx(m), a)? == Some(false) {
            return Ok(true);
        }
    }
    Ok(false)
}

fn matrix_pred(cx: &Cx, p: Pred, e: &Expr, a: &Expr) -> R {
    use Pred::*;
    if let Kind::Func(Func::MatrixTrace(m) | Func::MatrixDeterminant(m), _) = e.kind() {
        // Determinant/Trace registrations of the *_elements predicates.
        if matches!(p, IntegerElements | RealElements | ComplexElements) {
            return closed_group(cx, p, &[mx(m)], a);
        }
        return Ok(None);
    }
    let Some(m) = as_mx(e) else { return Ok(None) };
    if let MatrixKind::Permutation(images) = m.kind() {
        let identity = images.iter().enumerate().all(|(i,j)|i==*j);
        let result = match p {
            Square | Invertible | Fullrank | Orthogonal | Unitary | Normal |
                IntegerElements | RealElements | ComplexElements => Some(true),
            Singular => Some(false),
            Symmetric | Hermitian => Some(images.iter().enumerate().all(|(i,j)|images[*j]==i)),
            Diagonal | UpperTriangular | LowerTriangular | Triangular |
                UnitTriangular | PositiveDefinite => Some(identity),
            _ => None,
        };
        if result.is_some() { return Ok(result); }
    }
    if let MatrixKind::Permute(parent,_,_) = m.kind() {
        if matches!(p,Invertible|Fullrank|Orthogonal|Unitary|Singular|
                    IntegerElements|RealElements|ComplexElements) {
            return cx.askp(p,&mx(parent),a);
        }
    }
    let a_in = |pp: Pred| in_conjuncts(&q(pp, e), a);
    match p {
        Square => Ok(Some(is_square(m))),
        Symmetric => match m.kind() {
            MatrixKind::Mul(_, fs) => {
                if all_matrices(cx, Symmetric, fs, a)? {
                    return Ok(Some(true));
                }
                if cx.is(Diagonal, e, a)? {
                    return Ok(Some(true));
                }
                if fs.len() >= 2 && fs[0].transpose(true).ok().as_ref() == fs.last() {
                    if fs.len() == 2 {
                        return Ok(Some(true));
                    }
                    let Some(inner) = matmul(&fs[1..fs.len() - 1]) else {
                        return Ok(None);
                    };
                    return cx.askp(Symmetric, &inner, a);
                }
                Ok(None)
            }
            MatrixKind::Pow(b, x) => matpow_nonneg_or_invertible(cx, Symmetric, b, x, a),
            MatrixKind::Add(ts) => Ok(Some(all_matrices(cx, Symmetric, ts, a)?)),
            MatrixKind::Symbol(_) => {
                if !is_square(m) {
                    return Ok(Some(false));
                }
                if cx.is(Diagonal, e, a)? || a_in(Symmetric) {
                    return Ok(Some(true));
                }
                Ok(None)
            }
            MatrixKind::One | MatrixKind::Zero => cx.askp(Square, e, a),
            MatrixKind::Inverse(b) | MatrixKind::Transpose(b) => cx.askp(Symmetric, &mx(b), a),
            MatrixKind::Slice(parent, rs, cs) => {
                if cx.is(Diagonal, e, a)? {
                    return Ok(Some(true));
                }
                if !on_diag(rs, cs) {
                    return Ok(None);
                }
                cx.askp(Symmetric, &mx(parent), a)
            }
            MatrixKind::Identity => Ok(Some(true)),
            _ => Ok(None),
        },
        Invertible => match m.kind() {
            MatrixKind::Mul(_, fs) => {
                if all_matrices(cx, Invertible, fs, a)? {
                    return Ok(Some(true));
                }
                if any_false(cx, Invertible, fs, a)? {
                    return Ok(Some(false));
                }
                Ok(None)
            }
            MatrixKind::Pow(b, x) => {
                if !cx.is(Integer, x, a)? {
                    return Ok(None);
                }
                if crate::assumptions::is_negative(x) == Some(false) {
                    return cx.askp(Invertible, &mx(b), a);
                }
                Ok(None)
            }
            MatrixKind::Add(_) => Ok(None),
            MatrixKind::Symbol(_) => {
                if !is_square(m) {
                    return Ok(Some(false));
                }
                Ok(if a_in(Invertible) { Some(true) } else { None })
            }
            MatrixKind::Identity | MatrixKind::Inverse(_) => Ok(Some(true)),
            MatrixKind::Zero => Ok(Some(false)),
            MatrixKind::One => Ok(Some(is_shape(m, 1))),
            MatrixKind::Transpose(b) => cx.askp(Invertible, &mx(b), a),
            MatrixKind::Slice(parent, rs, cs) => {
                if !on_diag(rs, cs) {
                    return Ok(None);
                }
                cx.askp(Invertible, &mx(parent), a)
            }
            MatrixKind::Explicit(mat) => {
                if !mat.is_square() {
                    return Ok(Some(false));
                }
                Ok(Some(mat.rank() == mat.rows()))
            }
            MatrixKind::BlockDiag(blocks) => {
                let mut rv = Some(true);
                for b in blocks {
                    if b.rows() != b.cols() {
                        return Ok(None);
                    }
                }
                for b in blocks {
                    match cx.askp(Invertible, &mx(b), a)? {
                        Some(false) => return Ok(Some(false)),
                        None => rv = None,
                        _ => {}
                    }
                }
                Ok(rv)
            }
            MatrixKind::Block(grid) => {
                if !is_square(m) {
                    return Ok(Some(false));
                }
                if grid.shape() == (1, 1) {
                    return cx.askp(Invertible, &mx(grid.get(0, 0)), a);
                }
                Ok(None)
            }
            _ => Ok(if !is_square(m) { Some(false) } else { None }),
        },
        Orthogonal | Unitary => {
            let unitary = p == Unitary;
            match m.kind() {
                MatrixKind::Mul(c, fs) => {
                    let unit = if unitary {
                        crate::functions::func(Func::Abs, vec![c.clone()]).is_one()
                    } else {
                        c.is_one()
                    };
                    if all_matrices(cx, p, fs, a)? && unit {
                        return Ok(Some(true));
                    }
                    if any_false(cx, Invertible, fs, a)? {
                        return Ok(Some(false));
                    }
                    Ok(None)
                }
                MatrixKind::Pow(b, x) => {
                    if cx.is(Integer, x, a)? {
                        return cx.askp(p, &mx(b), a);
                    }
                    Ok(None)
                }
                MatrixKind::Add(ts) if !unitary => {
                    Ok(if ts.len() == 1 && cx.is(p, &mx(&ts[0]), a)? {
                        Some(true)
                    } else {
                        None
                    })
                }
                MatrixKind::Symbol(_) => {
                    if !is_square(m) || cx.askp(Invertible, e, a)? == Some(false) {
                        return Ok(Some(false));
                    }
                    Ok(if a_in(p) { Some(true) } else { None })
                }
                MatrixKind::Identity => Ok(Some(true)),
                MatrixKind::Zero => Ok(Some(false)),
                MatrixKind::Inverse(b) | MatrixKind::Transpose(b) => cx.askp(p, &mx(b), a),
                MatrixKind::Slice(parent, rs, cs) => {
                    if !on_diag(rs, cs) {
                        return Ok(None);
                    }
                    cx.askp(p, &mx(parent), a)
                }
                _ => Ok(None),
            }
        }
        Fullrank => match m.kind() {
            MatrixKind::Mul(c, fs) => {
                for x in matmul_args(c, fs) {
                    if !cx.is(Fullrank, &x, a)? {
                        return Ok(None);
                    }
                }
                Ok(Some(true))
            }
            MatrixKind::Pow(b, x) => {
                if cx.is(Integer, x, a)? && cx.ask(&not(q(Negative, x)), a)? == Some(true) {
                    return cx.askp(Fullrank, &mx(b), a);
                }
                Ok(None)
            }
            MatrixKind::Identity => Ok(Some(true)),
            MatrixKind::Zero => Ok(Some(false)),
            MatrixKind::One => Ok(Some(is_shape(m, 1))),
            MatrixKind::Inverse(b) | MatrixKind::Transpose(b) => cx.askp(Fullrank, &mx(b), a),
            MatrixKind::Slice(parent, ..) => Ok(if cx.is(Orthogonal, &mx(parent), a)? {
                Some(true)
            } else {
                None
            }),
            _ => Ok(None),
        },
        PositiveDefinite => match m.kind() {
            MatrixKind::Mul(c, fs) => {
                if all_matrices(cx, PositiveDefinite, fs, a)? {
                    match crate::relational::gt(c, &Expr::zero()).kind() {
                        Kind::Bool(true) => return Ok(Some(true)),
                        Kind::Bool(false) => {}
                        _ => {
                            return Err(AskError::Type(
                                "cannot determine truth value of Relational".into(),
                            ))
                        }
                    }
                }
                if fs.len() >= 2
                    && fs[0].transpose(true).ok().as_ref() == fs.last()
                    && cx.is(Fullrank, &mx(&fs[0]), a)?
                {
                    let Some(inner) = matmul(&fs[1..fs.len() - 1]) else {
                        return Ok(None);
                    };
                    return cx.askp(PositiveDefinite, &inner, a);
                }
                Ok(None)
            }
            MatrixKind::Pow(b, _) => Ok(if cx.is(PositiveDefinite, &mx(b), a)? {
                Some(true)
            } else {
                None
            }),
            MatrixKind::Add(ts) => Ok(if all_matrices(cx, PositiveDefinite, ts, a)? {
                Some(true)
            } else {
                None
            }),
            MatrixKind::Symbol(_) => {
                if !is_square(m) {
                    return Ok(Some(false));
                }
                Ok(if a_in(PositiveDefinite) {
                    Some(true)
                } else {
                    None
                })
            }
            MatrixKind::Identity => Ok(Some(true)),
            MatrixKind::Zero => Ok(Some(false)),
            MatrixKind::One => Ok(Some(is_shape(m, 1))),
            MatrixKind::Inverse(b) | MatrixKind::Transpose(b) => {
                cx.askp(PositiveDefinite, &mx(b), a)
            }
            MatrixKind::Slice(parent, rs, cs) => {
                if !on_diag(rs, cs) {
                    return Ok(None);
                }
                cx.askp(PositiveDefinite, &mx(parent), a)
            }
            _ => Ok(None),
        },
        UpperTriangular | LowerTriangular => {
            let other = if p == UpperTriangular {
                LowerTriangular
            } else {
                UpperTriangular
            };
            match m.kind() {
                MatrixKind::Mul(_, fs) => Ok(if all_matrices(cx, p, fs, a)? {
                    Some(true)
                } else {
                    None
                }),
                MatrixKind::Add(ts) => Ok(if all_matrices(cx, p, ts, a)? {
                    Some(true)
                } else {
                    None
                }),
                MatrixKind::Pow(b, x) => matpow_nonneg_or_invertible(cx, p, b, x, a),
                MatrixKind::Symbol(_) => Ok(if a_in(p) { Some(true) } else { None }),
                MatrixKind::Identity | MatrixKind::Zero => Ok(Some(true)),
                MatrixKind::One => Ok(Some(is_shape(m, 1))),
                MatrixKind::Transpose(b) => cx.askp(other, &mx(b), a),
                MatrixKind::Inverse(b) => cx.askp(p, &mx(b), a),
                MatrixKind::Slice(parent, rs, cs) => {
                    if !on_diag(rs, cs) {
                        return Ok(None);
                    }
                    cx.askp(p, &mx(parent), a)
                }
                _ => Ok(None),
            }
        }
        Diagonal => match m.kind() {
            MatrixKind::Mul(_, fs) => {
                if empty_or_1x1(m) {
                    return Ok(Some(true));
                }
                Ok(if all_matrices(cx, Diagonal, fs, a)? {
                    Some(true)
                } else {
                    None
                })
            }
            MatrixKind::Pow(b, x) => matpow_nonneg_or_invertible(cx, Diagonal, b, x, a),
            MatrixKind::Add(ts) => Ok(if all_matrices(cx, Diagonal, ts, a)? {
                Some(true)
            } else {
                None
            }),
            MatrixKind::Symbol(_) => Ok(if empty_or_1x1(m) || a_in(Diagonal) {
                Some(true)
            } else {
                None
            }),
            MatrixKind::One => Ok(Some(is_shape(m, 1))),
            MatrixKind::Inverse(b) | MatrixKind::Transpose(b) => cx.askp(Diagonal, &mx(b), a),
            MatrixKind::Slice(parent, rs, cs) => {
                if empty_or_1x1(m) {
                    return Ok(Some(true));
                }
                if !on_diag(rs, cs) {
                    return Ok(None);
                }
                cx.askp(Diagonal, &mx(parent), a)
            }
            MatrixKind::Diagonal(_)
            | MatrixKind::Diag(_)
            | MatrixKind::Identity
            | MatrixKind::Zero => Ok(Some(true)),
            _ => Ok(None),
        },
        IntegerElements | RealElements | ComplexElements => {
            let scalar = match p {
                IntegerElements => Integer,
                RealElements => Real,
                _ => Complex,
            };
            match m.kind() {
                MatrixKind::HadamardProduct(ts) | MatrixKind::Add(ts) => {
                    closed_group(cx, p, &ts.iter().map(mx).collect::<Vec<_>>(), a)
                }
                MatrixKind::Transpose(b) => closed_group(cx, p, &[mx(b)], a),
                MatrixKind::Inverse(b) if p == ComplexElements => {
                    closed_group(cx, p, &[mx(b), Expr::minus_one()], a)
                }
                MatrixKind::Inverse(b) if p == IntegerElements => {
                    let _ = b;
                    Ok(None)
                }
                MatrixKind::Inverse(b) => {
                    matpow_nonneg_or_invertible(cx, p, b, &Expr::minus_one(), a)
                }
                MatrixKind::Pow(b, x) => {
                    if p == IntegerElements {
                        if !cx.is(Integer, x, a)? {
                            return Ok(None);
                        }
                        if crate::assumptions::is_negative(x) == Some(false) {
                            return cx.askp(p, &mx(b), a);
                        }
                        return Ok(None);
                    }
                    matpow_nonneg_or_invertible(cx, p, b, x, a)
                }
                MatrixKind::Identity | MatrixKind::One | MatrixKind::Zero
                    if p == IntegerElements =>
                {
                    Ok(Some(true))
                }
                MatrixKind::Mul(c, fs) => {
                    let scalars: Vec<Expr> = if c.is_one() { vec![] } else { vec![c.clone()] };
                    let s = closed_group(cx, scalar, &scalars, a)?;
                    if s == Some(false) {
                        return Ok(Some(false));
                    }
                    let mm = closed_group(cx, p, &fs.iter().map(mx).collect::<Vec<_>>(), a)?;
                    Ok(match (s, mm) {
                        (_, Some(false)) => Some(false),
                        (Some(true), Some(true)) => Some(true),
                        _ => None,
                    })
                }
                MatrixKind::Slice(parent, ..) => cx.askp(p, &mx(parent), a),
                MatrixKind::Block(grid) => {
                    for b in grid.data() {
                        if !cx.is(p, &mx(b), a)? {
                            return Ok(Some(false));
                        }
                    }
                    Ok(Some(true))
                }
                _ => Ok(None),
            }
        }
        _ => Ok(None),
    }
}
