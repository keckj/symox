//! Automatic evaluation of `re`, `im`, `arg`, `Abs`, `sign` and `conjugate`.
//!
//! Each evaluator follows the corresponding SymPy `eval` classmethod (and
//! `_eval_conjugate` dispatch for `conjugate`). Returning `None` keeps the
//! application unevaluated. Advanced splitting goes through the bounded
//! `as_real_imag`; resource-limit errors and deep re-entrance keep the
//! application held instead of failing.
use super::as_real_imag;
use crate::expr::Constant;
use crate::functions as f;
use crate::ops::{add, div, mul, mul2, neg, pow, sub};
use crate::{Expr, Func, Kind, Number};
use std::cell::Cell;

thread_local! {
    static DEPTH: Cell<u32> = const { Cell::new(0) };
}
/// Nested `as_real_imag` calls issued from the evaluators themselves.
const MAX_DEPTH: u32 = 8;

fn raw(function: Func, args: Vec<Expr>) -> Expr {
    Expr::from_kind(Kind::Func(function, args))
}
fn fact(e: &Expr, name: &str) -> Option<bool> {
    use crate::assumptions as a;
    let value = a::fact(e, name);
    if value.is_some() || !f::is_numeric_value(e) {
        return value;
    }
    // Numbers: SymPy decides signs by evaluation, which the real-sign
    // queries certify even when the extended fact is left open.
    let either = |yes: Option<bool>, no: Option<bool>| match (yes, no) {
        (Some(true), _) => Some(true),
        (_, Some(true)) => Some(false),
        _ => None,
    };
    match name {
        "extended_positive" | "positive" => either(a::is_positive(e), a::is_nonpositive(e)),
        "extended_negative" | "negative" => either(a::is_negative(e), a::is_nonnegative(e)),
        "extended_nonnegative" | "nonnegative" => either(a::is_nonnegative(e), a::is_negative(e)),
        "extended_nonpositive" | "nonpositive" => either(a::is_nonpositive(e), a::is_positive(e)),
        "zero" => a::is_zero(e),
        "extended_nonzero" => a::is_zero(e).map(|z| !z),
        _ => None,
    }
}
fn holds(e: &Expr, name: &str) -> bool {
    fact(e, name) == Some(true)
}
fn i_times(e: &Expr) -> Expr {
    mul2(&Expr::i(), e)
}
fn minus_i_times(e: &Expr) -> Expr {
    mul(&[Expr::minus_one(), Expr::i(), e.clone()])
}
/// `is_comparable`: a number that evaluates to a real value.
fn comparable(e: &Expr) -> bool {
    f::is_numeric_value(e) && crate::evalf::eval_f64(e).is_some()
}
fn is_nan(e: &Expr) -> bool {
    e.is_constant(Constant::NaN)
}
fn is_zoo(e: &Expr) -> bool {
    e.is_constant(Constant::ComplexInfinity)
}
fn sign_of_number(n: &Number) -> Expr {
    Expr::int(if n.is_positive() {
        1
    } else if n.is_negative() {
        -1
    } else {
        0
    })
}

/// `e.as_real_imag(ignore=ignore)`: `Some(None)` when SymPy returns `None`
/// (the ignored value), `None` when the split is not available here.
fn parts(e: &Expr, ignore: Option<&Expr>) -> Option<Option<(Expr, Expr)>> {
    DEPTH.with(|depth| {
        let level = depth.get();
        if level >= MAX_DEPTH {
            return None;
        }
        depth.set(level + 1);
        let value = as_real_imag(e, true, ignore).ok();
        depth.set(level);
        value
    })
}

/// SymPy's `term.as_coefficient(I)` for the shapes canonical terms take.
fn coefficient_of_i(term: &Expr) -> Option<Expr> {
    let value = match term.kind() {
        Kind::Const(Constant::I) => Expr::one(),
        Kind::Mul(c, factors) => {
            let index = factors
                .iter()
                .position(|(b, x)| b.is_constant(Constant::I) && x.is_one())?;
            let mut rest: Vec<Expr> = factors
                .iter()
                .enumerate()
                .filter(|(k, _)| *k != index)
                .map(|(_, (b, x))| crate::ops::mul_factor(b, x))
                .collect();
            rest.push(Expr::num(c.clone()));
            mul(&rest)
        }
        _ => return None,
    };
    (!f::contains_i(&value)).then_some(value)
}

fn add_terms(e: &Expr) -> Vec<Expr> {
    match e.kind() {
        Kind::Add(c, terms) => {
            let mut out = Vec::with_capacity(terms.len() + 1);
            if !c.is_zero() {
                out.push(Expr::num(c.clone()));
            }
            out.extend(terms.iter().map(|(t, k)| crate::ops::add_term(t, k)));
            out
        }
        _ => vec![e.clone()],
    }
}

/// `re.eval` (`want_re`) and `im.eval`.
pub(crate) fn re_im(a: &Expr, want_re: bool) -> Option<Expr> {
    if is_nan(a) || is_zoo(a) {
        return Some(Expr::nan());
    }
    if holds(a, "extended_real") {
        return Some(if want_re { a.clone() } else { Expr::zero() });
    }
    if holds(a, "imaginary") || holds(&i_times(a), "extended_real") {
        return Some(if want_re { Expr::zero() } else { minus_i_times(a) });
    }
    if let Kind::Func(Func::Conjugate, args) = a.kind() {
        return Some(if want_re { f::re(&args[0]) } else { neg(&f::im(&args[0])) });
    }
    let terms = add_terms(a);
    let (mut included, mut reverted, mut excluded) = (Vec::new(), Vec::new(), Vec::new());
    for term in &terms {
        if let Some(coefficient) = coefficient_of_i(term) {
            if !holds(&coefficient, "extended_real") {
                reverted.push(coefficient);
            } else if !want_re {
                excluded.push(coefficient);
            }
        } else if !f::contains_i(term) && holds(term, "extended_real") {
            if want_re {
                excluded.push(term.clone());
            }
        } else {
            match parts(term, Some(a)) {
                Some(Some((r, i))) => excluded.push(if want_re { r } else { i }),
                _ => included.push(term.clone()),
            }
        }
    }
    if included.len() == terms.len() {
        return None;
    }
    let (a, b, c) = (add(&included), add(&reverted), add(&excluded));
    Some(if want_re {
        add(&[f::re(&a), neg(&f::im(&b)), c])
    } else {
        add(&[f::im(&a), f::re(&b), c])
    })
}

/// SymPy's `atan2.eval`, without the Piecewise and logarithmic fallbacks for
/// numbers whose signs cannot be decided.
fn atan2(y: &Expr, x: &Expr) -> Option<Expr> {
    if is_nan(x) || is_nan(y) {
        return Some(Expr::nan());
    }
    let pi = Expr::pi();
    if x.is_constant(Constant::NegInfinity) {
        if holds(y, "zero") {
            return Some(pi);
        }
        return Some(sub(
            &mul(&[Expr::int(2), pi, f::heaviside(&f::re(y))]),
            &Expr::pi(),
        ));
    }
    if x.is_constant(Constant::Infinity) {
        return Some(Expr::zero());
    }
    let (x, y) = if holds(x, "imaginary")
        && holds(y, "imaginary")
        && f::is_numeric_value(x)
        && f::is_numeric_value(y)
    {
        (f::im(x), f::im(y))
    } else {
        (x.clone(), y.clone())
    };
    let ratio = || {
        if y.is_zero() {
            Expr::zero()
        } else {
            f::atan(&div(&y, &x))
        }
    };
    if holds(&x, "extended_real") && holds(&y, "extended_real") {
        if holds(&x, "positive") {
            return Some(ratio());
        } else if holds(&x, "negative") {
            if holds(&y, "negative") {
                return Some(sub(&ratio(), &pi));
            } else if holds(&y, "nonnegative") {
                return Some(add(&[ratio(), pi]));
            }
        } else if holds(&x, "zero") {
            if holds(&y, "positive") {
                return Some(div(&pi, &Expr::int(2)));
            } else if holds(&y, "negative") {
                return Some(div(&pi, &Expr::int(-2)));
            } else if holds(&y, "zero") {
                return Some(Expr::nan());
            }
        }
    }
    if holds(&y, "zero") && holds(&x, "extended_nonzero") {
        return Some(mul2(&pi, &sub(&Expr::one(), &f::heaviside(&x))));
    }
    None
}

/// `im(z) mod 2*pi`, reduced into `(-pi, pi]`.
fn principal_angle(value: &Expr) -> Option<Expr> {
    let v = crate::evalf::eval_f64(value)?;
    let turns = (v / (2.0 * std::f64::consts::PI)).floor();
    if !turns.is_finite() || turns.abs() > 1e15 {
        return None;
    }
    let two_pi = mul2(&Expr::int(2), &Expr::pi());
    let reduced = sub(value, &mul2(&Expr::int(turns as i64), &two_pi));
    match crate::assumptions::is_positive(&sub(&reduced, &Expr::pi())) {
        Some(true) => Some(sub(&reduced, &two_pi)),
        Some(false) => Some(reduced),
        None => None,
    }
}

fn applied_undefined(e: &Expr, out: &mut Vec<Expr>) {
    if let Kind::Func(Func::Undefined(head), _) = e.kind() {
        if !head.is_internal("Tuple") && !out.contains(e) {
            out.push(e.clone());
        }
    }
    for value in e.args() {
        applied_undefined(&value, out);
    }
}

/// `arg.eval`.
pub(crate) fn arg(a: &Expr) -> Option<Expr> {
    let mut inner = a.clone();
    let mut peeled = 0;
    loop {
        if peeled == 3 {
            return Some(Expr::nan());
        }
        match inner.kind() {
            Kind::Func(Func::Arg, args) => {
                let next = args[0].clone();
                inner = next;
                peeled += 1;
            }
            _ => {
                if peeled == 2 && holds(&inner, "extended_real") {
                    return Some(Expr::nan());
                }
                break;
            }
        }
    }
    if let Kind::Func(Func::Exp, args) = a.kind() {
        let imaginary = f::im(&args[0]);
        if comparable(&imaginary) {
            if let Some(angle) = principal_angle(&imaginary) {
                return Some(angle);
            }
        }
    }
    let reduced = if matches!(a.kind(), Kind::Sym(_) | Kind::Num(_) | Kind::Const(_)) {
        a.clone()
    } else {
        // factor_terms(arg).as_coeff_Mul(), then factors of known sign drop.
        let (content, primitive) = crate::decompose::as_content_primitive(a, false, false)
            .unwrap_or((Number::ONE, a.clone()));
        let (coefficient, rest) = match primitive.kind() {
            Kind::Mul(c, factors) => (
                content.mul(c),
                mul(&factors
                    .iter()
                    .map(|(b, x)| crate::ops::mul_factor(b, x))
                    .collect::<Vec<_>>()),
            ),
            _ => (content, primitive.clone()),
        };
        let rest = match rest.kind() {
            Kind::Mul(_, factors) => mul(&factors
                .iter()
                .map(|(b, x)| {
                    let factor = crate::ops::mul_factor(b, x);
                    let s = f::sign(&factor);
                    if s.is_one() || s.is_minus_one() {
                        s
                    } else {
                        factor
                    }
                })
                .collect::<Vec<_>>()),
            _ => rest,
        };
        mul2(&sign_of_number(&coefficient), &rest)
    };
    let mut undefined = Vec::new();
    applied_undefined(&reduced, &mut undefined);
    if undefined.iter().any(|u| fact(u, "extended_positive").is_none()) {
        return None;
    }
    if let Some(Some((x, y))) = parts(&reduced, None) {
        if let Some(value) = atan2(&y, &x) {
            if is_nan(&value) || f::is_numeric_value(&value) {
                return Some(value);
            }
        }
    }
    (reduced != *a).then(|| raw(Func::Arg, vec![reduced]))
}

/// `sign.eval`.
pub(crate) fn sign(a: &Expr) -> Option<Expr> {
    if let Kind::Mul(c, factors) = a.kind() {
        let mut s = sign_of_number(c);
        let mut unknown = Vec::new();
        for (b, x) in factors {
            let factor = crate::ops::mul_factor(b, x);
            if holds(&factor, "extended_negative") {
                s = neg(&s);
            } else if holds(&factor, "extended_positive") {
            } else if holds(&factor, "imaginary") {
                let imaginary = f::im(&factor);
                if comparable(&imaginary) {
                    s = mul2(&s, &Expr::i());
                    if holds(&imaginary, "extended_negative") {
                        s = neg(&s);
                    }
                } else {
                    unknown.push(factor);
                }
            } else {
                unknown.push(factor);
            }
        }
        if c.is_one() && unknown.len() == factors.len() {
            return None;
        }
        return Some(mul2(&s, &f::sign(&mul(&unknown))));
    }
    if is_nan(a) {
        return Some(Expr::nan());
    }
    if holds(a, "zero") {
        return Some(Expr::zero());
    }
    if holds(a, "extended_positive") {
        return Some(Expr::one());
    }
    if holds(a, "extended_negative") {
        return Some(Expr::minus_one());
    }
    if a.is_func(&Func::Sign) {
        return Some(a.clone());
    }
    // sqrt of a negative value is imaginary even where the general fact
    // engine cannot tell.
    if matches!(a.kind(), Kind::Pow(b, x) if *x == Expr::rational(1, 2) && holds(b, "extended_negative")) {
        return Some(Expr::i());
    }
    if holds(a, "imaginary") {
        if matches!(a.kind(), Kind::Pow(_, x) if *x == Expr::rational(1, 2)) {
            return Some(Expr::i());
        }
        let rotated = minus_i_times(a);
        if holds(&rotated, "extended_positive") {
            return Some(Expr::i());
        }
        if holds(&rotated, "extended_negative") {
            return Some(neg(&Expr::i()));
        }
    }
    None
}

/// `Abs.eval`.
pub(crate) fn abs(a: &Expr) -> Option<Expr> {
    match a.kind() {
        Kind::Num(n) => return Some(Expr::num(n.abs())),
        // _eval_Abs hooks
        Kind::Func(Func::Conjugate, args) => return Some(f::abs(&args[0])),
        Kind::Func(Func::Sign, args) if fact(&args[0], "zero") == Some(false) => {
            return Some(Expr::one())
        }
        _ => {}
    }
    // signsimp: a sum that prefers a leading minus sign is a negated product.
    if !a.is_mul() && crate::decompose::could_extract_minus_sign(a) {
        return Some(f::abs(&neg(a)));
    }
    let (n, d) = crate::simplify::numer_denom(a);
    if !d.is_one()
        && !crate::traverse::free_symbols(&d).is_empty()
        && crate::traverse::free_symbols(&n).is_empty()
    {
        return Some(div(&f::abs(&n), &f::abs(&d)));
    }
    if let Kind::Mul(c, factors) = a.kind() {
        let mut known = vec![Expr::num(c.abs())];
        let mut unknown = Vec::new();
        for (b, x) in factors {
            let factor = crate::ops::mul_factor(b, x);
            match factor.kind() {
                Kind::Pow(base, exponent)
                    if holds(exponent, "integer") && holds(exponent, "negative") =>
                {
                    let value = f::abs(base);
                    if value.is_func(&Func::Abs) {
                        unknown.push(factor.clone());
                    } else {
                        known.push(pow(&value, exponent));
                    }
                }
                _ => {
                    let value = f::abs(&factor);
                    if value.is_func(&Func::Abs) {
                        unknown.push(factor);
                    } else {
                        known.push(value);
                    }
                }
            }
        }
        let unknown = if unknown.is_empty() {
            Expr::one()
        } else {
            raw(Func::Abs, vec![mul(&unknown)])
        };
        return Some(mul2(&mul(&known), &unknown));
    }
    if is_nan(a) {
        return Some(Expr::nan());
    }
    if is_zoo(a) {
        return Some(Expr::oo());
    }
    if let Kind::Pow(base, exponent) = a.kind() {
        if holds(base, "extended_real") {
            if holds(exponent, "integer") {
                if holds(exponent, "even") {
                    return Some(a.clone());
                }
                if base.is_minus_one() {
                    return Some(Expr::one());
                }
                return Some(pow(&f::abs(base), exponent));
            }
            if holds(base, "extended_nonnegative") {
                return Some(pow(base, &f::re(exponent)));
            }
            if holds(base, "extended_negative") {
                return Some(mul2(
                    &pow(&neg(base), &f::re(exponent)),
                    &f::exp(&neg(&mul2(&Expr::pi(), &f::im(exponent)))),
                ));
            }
            return None;
        } else if crate::traverse::free_symbols(base).is_empty() {
            if let Some(Some((x, y))) = parts(&f::log(base), None) {
                let z = add(&[x, i_times(&y)]);
                return Some(f::exp(&f::re(&mul2(exponent, &z))));
            }
        }
    }
    if let Kind::Func(Func::Exp, args) = a.kind() {
        return Some(f::exp(&f::re(&args[0])));
    }
    if let Kind::Func(Func::Undefined(head), _) = a.kind() {
        if !head.is_internal("Tuple") {
            if holds(a, "positive") {
                return Some(a.clone());
            }
            if holds(a, "negative") {
                return Some(neg(a));
            }
            return None;
        }
    }
    if a.is_add()
        && (crate::traverse::has(a, &Expr::oo()) || crate::traverse::has(a, &Expr::neg_oo()))
    {
        if let Some(Some((x, y))) = parts(a, None) {
            if holds(&x, "infinite") || holds(&y, "infinite") {
                return Some(Expr::oo());
            }
        }
    }
    if holds(a, "zero") {
        return Some(Expr::zero());
    }
    if holds(a, "extended_nonnegative") {
        return Some(a.clone());
    }
    if holds(a, "extended_nonpositive") {
        return Some(neg(a));
    }
    if holds(a, "imaginary") {
        let rotated = minus_i_times(a);
        if holds(&rotated, "extended_nonnegative") {
            return Some(rotated);
        }
    }
    if holds(a, "extended_real") {
        return None;
    }
    super::abs_conjugate_rule(a)
}

/// `conjugate.eval`, i.e. the `_eval_conjugate` dispatch.
pub(crate) fn conjugate(a: &Expr) -> Option<Expr> {
    let c = f::conjugate;
    let default = |a: &Expr| {
        if holds(a, "extended_real") {
            Some(a.clone())
        } else if holds(a, "imaginary") {
            Some(neg(a))
        } else {
            None
        }
    };
    match a.kind() {
        Kind::Num(_) | Kind::Const(Constant::NaN) => Some(a.clone()),
        Kind::Const(Constant::I) => Some(neg(a)),
        Kind::Add(k, terms) => {
            let mut out = vec![Expr::num(k.clone())];
            out.extend(terms.iter().map(|(t, k)| mul2(&Expr::num(k.clone()), &c(t))));
            Some(add(&out))
        }
        Kind::Mul(k, factors) => {
            let mut out = vec![Expr::num(k.clone())];
            out.extend(factors.iter().map(|(b, x)| c(&crate::ops::mul_factor(b, x))));
            Some(mul(&out))
        }
        Kind::Pow(base, exponent) => {
            let integer = fact(exponent, "integer");
            let positive = fact(base, "positive");
            if integer == Some(true) {
                return Some(pow(&c(base), exponent));
            }
            if positive == Some(true) {
                return Some(pow(base, &c(exponent)));
            }
            if integer == Some(false) && positive == Some(false) {
                let options = crate::expand::ExpandOptions {
                    basic: false,
                    mul: false,
                    multinomial: false,
                    power_base: false,
                    power_exp: false,
                    log: false,
                    ..crate::expand::ExpandOptions::default()
                };
                if let Ok(expanded) = super::expand_complex(a, &options) {
                    if expanded != *a {
                        return Some(c(&expanded));
                    }
                }
            }
            holds(a, "extended_real").then(|| a.clone())
        }
        Kind::Func(Func::Conjugate, args) => Some(args[0].clone()),
        Kind::Func(
            function @ (Func::Exp
            | Func::Sin
            | Func::Cos
            | Func::Tan
            | Func::Cot
            | Func::Sec
            | Func::Csc
            | Func::Sinh
            | Func::Cosh
            | Func::Tanh
            | Func::Coth
            | Func::Gamma
            | Func::Erf
            | Func::Erfc
            | Func::Sign),
            args,
        ) if args.len() == 1 => Some(f::func(function.clone(), vec![c(&args[0])])),
        Kind::Func(Func::Special(s), args) if matches!(s,
            crate::special::SpecialFn::MathieuC|crate::special::SpecialFn::MathieuS|
            crate::special::SpecialFn::MathieuCPrime|crate::special::SpecialFn::MathieuSPrime)=>
            Some(f::func(Func::Special(*s),args.iter().map(c).collect())),
        Kind::Func(Func::Special(s),args) if matches!(s,crate::special::SpecialFn::Hn1|crate::special::SpecialFn::Hn2)=>{
            let other=if *s==crate::special::SpecialFn::Hn1{crate::special::SpecialFn::Hn2}else{crate::special::SpecialFn::Hn1};
            Some(f::func(Func::Special(other),args.iter().map(c).collect()))
        },
        Kind::Func(Func::Atan2, args) => Some(f::atan2(&c(&args[0]), &c(&args[1]))),
        Kind::Func(Func::Acos, args) => {
            let z = &args[0];
            let real = fact(z, "extended_real");
            if real == Some(false)
                || (real == Some(true)
                    && holds(&add(&[z.clone(), Expr::one()]), "nonnegative")
                    && holds(&sub(z, &Expr::one()), "nonpositive"))
            {
                Some(f::acos(&c(z)))
            } else {
                None
            }
        }
        Kind::Func(Func::LogGamma, args) => {
            let z = &args[0];
            (!z.is_zero() && !z.is_constant(Constant::NegInfinity))
                .then(|| f::loggamma(&c(z)))
        }
        _ => default(a),
    }
}

#[cfg(test)]
mod tests {
    use crate::functions as f;
    use crate::ops::{add2, mul, mul2, neg, pow};
    use crate::{Assumptions, Expr};

    #[test]
    fn projections_follow_sympy_coefficient_rules() {
        let x = Expr::symbol("x");
        let r = Expr::symbol_with("r", Assumptions::REAL);
        let value = mul(&[x.clone(), Expr::i(), r.clone()]);
        assert_eq!(f::re(&value), neg(&mul2(&r, &f::im(&x))));
        assert_eq!(f::im(&value), mul2(&r, &f::re(&x)));
        assert_eq!(f::re(&f::conjugate(&x)), f::re(&x));
        assert_eq!(f::im(&f::conjugate(&x)), neg(&f::im(&x)));
    }

    #[test]
    fn exact_numbers_have_closed_form_arg_abs_and_sign() {
        let one_plus_i = add2(&Expr::one(), &Expr::i());
        assert_eq!(f::arg(&one_plus_i), mul2(&Expr::rational(1, 4), &Expr::pi()));
        let z = add2(&Expr::int(3), &mul2(&Expr::int(4), &Expr::i()));
        assert_eq!(f::abs(&z), Expr::int(5));
        let x = Expr::symbol("x");
        let value = mul(&[Expr::int(2), Expr::i(), x.clone()]);
        assert_eq!(f::sign(&value), mul2(&Expr::i(), &f::sign(&x)));
        assert_eq!(f::sign(&f::sign(&x)), f::sign(&x));
    }

    #[test]
    fn conjugate_and_abs_power_rules() {
        let x = Expr::symbol("x");
        let r = Expr::symbol_with("r", Assumptions::REAL);
        let value = f::exp(&mul2(&Expr::i(), &x));
        assert_eq!(f::conjugate(&value), f::exp(&neg(&mul2(&Expr::i(), &f::conjugate(&x)))));
        assert_eq!(f::conjugate(&f::conjugate(&x)), x);
        assert_eq!(pow(&f::abs(&r), &Expr::int(3)), mul2(&pow(&r, &Expr::int(2)), &f::abs(&r)));
        assert_eq!(f::abs(&neg(&x)), f::abs(&x));
    }
}
