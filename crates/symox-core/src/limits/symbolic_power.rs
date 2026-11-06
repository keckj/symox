//! The leading-term step of the reference `Limit.doit` for a power of the
//! small variable with a symbolic exponent.
//!
//! At a finite point `z0`, `Limit.doit` shifts `z -> z + z0` and, for an
//! expression that is not meromorphic (a power `z**p` with `p` not an
//! `Integer` literal), takes `coeff, ex = newe.leadterm(z, cdir)`. With a
//! coefficient free of `z` the answer is decided by the sign of `ex`:
//!
//! - `ex > 0` gives 0 and `ex == 0` gives `coeff`;
//! - `ex < 0` gives `oo*sign(coeff)` from the right and
//!   `-oo*sign(coeff)*(-1)**(1 + ex)` from the left;
//! - otherwise it raises `NotImplementedError("Not sure of sign of ex")`.
//!
//! The numeric-exponent engine cannot carry a symbolic exponent, so the
//! product `z**p * R(z)` is recognized here: the leading term `c*z**k` of
//! the rest `R` (analytic elementary factors only, a coefficient free of
//! `z` and of `log(z)`) is computed by the engine and `ex = p + k`. Any
//! other shape is left to the engine.

use super::{unknown, LRes};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, mul, mul2, pow, sub};
use crate::series::engine::depends;

/// The message of the reference for an exponent of undecided sign.
pub const NOT_SURE_OF_SIGN: &str = "Not sure of sign of ";

/// Only analytic elementary functions in the rest of the product: the
/// reference first rewrites `Abs`, `sign` and `arg` (`set_signs`), and
/// rounding or piecewise functions have one-sided leading terms of their own.
fn elementary(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) => true,
        Kind::Const(c) => !c.is_infinite() && *c != Constant::NaN,
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => e.args().iter().all(elementary),
        Kind::Func(f, args) => {
            matches!(
                f,
                Func::Exp | Func::Log | Func::Sin | Func::Cos | Func::Tan | Func::Sinh | Func::Cosh | Func::Tanh | Func::Atan
            ) && args.iter().all(elementary)
        }
        _ => false,
    }
}

/// `Some(result)` when `e` is `(x - x0)**p * R` with a symbolic exponent
/// `p` free of `x` and the leading term of `R` decided; `None` otherwise.
pub fn symbolic_power_limit(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> Option<LRes<Expr>> {
    if crate::traverse::has(x0, &Expr::oo()) || crate::traverse::has(x0, &Expr::neg_oo()) || crate::traverse::has(x0, &Expr::zoo()) || crate::traverse::has(x0, &Expr::nan()) {
        return None;
    }
    let (coeff, facs): (Number, Vec<(Expr, Expr)>) = match e.kind() {
        Kind::Pow(b, p) => (Number::Int(1), vec![(b.clone(), p.clone())]),
        Kind::Mul(c, f) => (c.clone(), f.clone()),
        _ => return None,
    };
    let small = sub(x, x0);
    let mut sym_exp = None;
    let mut rest = vec![Expr::num(coeff)];
    for (b, p) in facs {
        if b == small && !depends(&p, x) && p.as_number().is_none() {
            if sym_exp.is_some() {
                return None;
            }
            sym_exp = Some(p);
        } else {
            rest.push(crate::ops::mul_factor(&b, &p));
        }
    }
    let p = sym_exp?;
    if crate::traverse::has(&p, &Expr::oo()) || crate::traverse::has(&p, &Expr::nan()) {
        return None;
    }
    let r = mul(&rest);
    let (c, k) = if depends(&r, x) {
        if !elementary(&r) {
            return None;
        }
        let lt = crate::series::leading_term(&r, x, x0, minus).ok()?;
        if depends(&lt.coeff_logt, x) || depends(&lt.coeff_logt, &lt.logt) || lt.coeff.is_zero() {
            return None;
        }
        if crate::traverse::has(&lt.coeff, &Expr::oo()) || crate::traverse::has(&lt.coeff, &Expr::nan()) {
            return None;
        }
        // the engine expands in `x0 - x > 0` from the left; the reference
        // leading term is in `z = x - x0 < 0`: equal up to `(-1)**k`, which
        // is only a sign for an integer `k`.
        if minus && !lt.exp.is_integer() {
            return None;
        }
        let c = if minus { mul2(&lt.coeff, &pow(&Expr::minus_one(), &Expr::num(lt.exp.clone()))) } else { lt.coeff.clone() };
        (c, Expr::num(lt.exp))
    } else {
        (r, Expr::zero())
    };
    let ex = add(&[p, k]);
    let sign = crate::functions::sign;
    Some(if crate::assumptions::is_positive(&ex) == Some(true) {
        Ok(Expr::zero())
    } else if crate::assumptions::is_zero(&ex) == Some(true) {
        Ok(c)
    } else if crate::assumptions::is_negative(&ex) == Some(true) {
        if minus {
            let unit = pow(&Expr::minus_one(), &add(&[Expr::one(), ex]));
            Ok(mul(&[Expr::neg_oo(), sign(&c), unit]))
        } else {
            Ok(mul2(&Expr::oo(), &sign(&c)))
        }
    } else {
        unknown(format!("{NOT_SURE_OF_SIGN}{ex}"))
    })
}

/// A power `b**p` of an algebraic or logarithmic base at an infinite point,
/// with an exponent `p` free of `x` whose sign is undecided. The reference
/// gruntz cannot conclude there (`ValueError: -sign(p) could not be
/// evaluated`) and its `heuristics` decide:
///
/// - a lone power is rebuilt from the limits of its arguments, so with
///   `b -> ±oo` the limit is `(±oo)**p` (`x**w` and `(x + 1)**w` at `oo`
///   give `oo**w`, `(-x)**w` gives `(-oo)**w`), held when that power is
///   `nan` (`x**I`);
/// - a product with other algebraic or logarithmic factors stays held
///   (`2*x**w`, `x**w*log(x)`: the limit `oo**w` of a factor is neither
///   finite nor infinite);
/// - for an integer `n`, `(x + c)**n` at `oo` reaches the leading-term step
///   as `z**(-n)` and raises `Not sure of sign of -n` (`2*x**n` too).
///
/// The bases with an exponential of `x` are left to the engine: they enter
/// the reference mrv set, whose rewrite raises `Result depends on the sign
/// of sign(p)` instead (`exp(x)**w`).
pub fn infinite_power_limit(e: &Expr, x: &Expr, x0: &Expr) -> Option<LRes<Expr>> {
    let pos = x0.is_constant(Constant::Infinity);
    if !pos && !x0.is_constant(Constant::NegInfinity) {
        return None;
    }
    let (lone, facs): (bool, Vec<(Expr, Expr)>) = match e.kind() {
        Kind::Pow(b, p) => (true, vec![(b.clone(), p.clone())]),
        Kind::Mul(_, f) => (false, f.clone()),
        _ => return None,
    };
    let symbolic = |b: &Expr, p: &Expr| {
        depends(b, x)
            && !depends(p, x)
            && p.as_number().is_none()
            && !crate::traverse::has(p, &Expr::oo())
            && !crate::traverse::has(p, &Expr::nan())
            && crate::assumptions::is_positive(p) != Some(true)
            && crate::assumptions::is_negative(p) != Some(true)
            && crate::assumptions::is_zero(p) != Some(true)
    };
    let powers: Vec<&(Expr, Expr)> = facs.iter().filter(|(b, p)| symbolic(b, p)).collect();
    let [(b, p)] = powers.as_slice() else { return None };
    if !algebraic_or_log(b, x) {
        return None;
    }
    let others_free = facs.iter().all(|f| std::ptr::eq(f, *powers.first().unwrap()) || !depends(&f.0, x) && !depends(&f.1, x));
    if crate::assumptions::is_integer(p) == Some(true) {
        let shift = sub(b, x);
        if pos && !depends(&shift, x) && others_free {
            return Some(unknown(format!("{NOT_SURE_OF_SIGN}{}", crate::ops::neg(p))));
        }
        return None;
    }
    let held = || Some(unknown("symbolic power at infinity (held by the reference heuristics)"));
    if !lone {
        let algebraic = facs.iter().all(|(fb, fp)| !depends(fp, x) && algebraic_or_log(fb, x));
        return if algebraic { held() } else { None };
    }
    let l = super::limit(b, x, x0, super::Dir::Plus).ok()?;
    if !(l.is_constant(Constant::Infinity) || l.is_constant(Constant::NegInfinity)) {
        return None;
    }
    let r = pow(&l, p);
    if crate::traverse::has(&r, &Expr::nan()) || crate::traverse::has(&r, &Expr::zoo()) {
        return held();
    }
    Some(Ok(r))
}

fn algebraic_or_log(e: &Expr, x: &Expr) -> bool {
    if !depends(e, x) {
        return true;
    }
    match e.kind() {
        Kind::Sym(_) => true,
        Kind::Add(..) | Kind::Mul(..) => e.args().iter().all(|a| algebraic_or_log(a, x)),
        Kind::Pow(b, q) => !depends(q, x) && algebraic_or_log(b, x),
        Kind::Func(Func::Log, a) => algebraic_or_log(&a[0], x),
        _ => false,
    }
}

/// `Limit.pow_heuristics` for a base free of `x` (SymPy 1.14): `b**p` is
/// `exp(limit(p*log(b)))`. The reference reaches it when the leading term
/// of `b**p` fails, before `gruntz`; it is taken here only after the
/// native engine does not decide, so an engine value is kept. Examples:
/// `(I/2)**x` at `oo` is 0, `(-2)**x` and `(1 + I)**x` are `zoo`. A `nan`
/// (`I**x`) is still undecided in `limit_one_sided`.
pub fn constant_base_power(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> Option<LRes<Expr>> {
    let Kind::Pow(b, p) = e.kind() else { return None };
    if depends(b, x) || !depends(p, x) || b.is_constant(Constant::E) {
        return None;
    }
    let log_b = crate::functions::func(Func::Log, vec![b.clone()]);
    let l = super::limit_one_sided(&mul2(p, &log_b), x, x0, minus).ok()?;
    Some(Ok(exp_of_limit(&l)))
}

/// `exp(l)` with the clause of SymPy's `exp.eval` for `l = coeff*terms`,
/// `coeff` in `(oo, -oo)` and `terms` a number: `re(terms) == 0` gives
/// `nan`, `re(terms) > 0` with `im(terms) != 0` gives `zoo`,
/// `re(terms) < 0` gives 0 (after `terms -> -terms` for `-oo`). Otherwise
/// the evaluation of the native `exp`.
fn exp_of_limit(l: &Expr) -> Expr {
    let plain = || crate::functions::func(Func::Exp, vec![l.clone()]);
    let Kind::Mul(c, facs) = l.kind() else { return plain() };
    let infinite = |b: &Expr| b.is_constant(Constant::Infinity) || b.is_constant(Constant::NegInfinity);
    let Some(at) = facs.iter().position(|(b, q)| infinite(b) && q.is_one()) else {
        return plain();
    };
    // `terms` up to a positive factor: each `sign(a)` is `a/|a|`
    let mut rest = vec![Expr::num(c.clone())];
    if facs[at].0.is_constant(Constant::NegInfinity) {
        rest.push(Expr::minus_one());
    }
    for (i, (b, q)) in facs.iter().enumerate() {
        if i == at {
            continue;
        }
        match b.kind() {
            Kind::Func(Func::Sign, a) if q.is_one() && a.len() == 1 => rest.push(a[0].clone()),
            _ => rest.push(pow(b, q)),
        }
    }
    let terms = mul(&rest);
    if !crate::traverse::free_symbols(&terms).is_empty() || terms.is_zero() {
        return plain();
    }
    let Ok(Some((re, im))) = crate::complex_parts::as_real_imag_bounded(&terms, true, None) else {
        return plain();
    };
    use crate::assumptions as asm;
    if asm::is_zero(&re) == Some(true) {
        return Expr::nan();
    }
    if asm::is_positive(&re) == Some(true) && asm::is_zero(&im) == Some(false) {
        return Expr::zoo();
    }
    if asm::is_negative(&re) == Some(true) {
        return Expr::zero();
    }
    plain()
}
