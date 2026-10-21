//! Error functions, Fresnel integrals and exponential integrals.

use super::*;
use crate::functions::{erf, erfc};
use crate::ops::{div, mul2, neg, pow, sub};

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    let z = &a[0];
    match s {
        Erfi => eval_erfi(z),
        Erf2 => eval_erf2(&a[0], &a[1]),
        Erfinv => {
            if z.is_minus_one() {
                return Some(Expr::neg_oo());
            }
            if z.is_one() {
                return Some(Expr::oo());
            }
            if z.is_zero() {
                return Some(Expr::zero());
            }
            if let Kind::Func(Func::Erf, x) = z.kind() {
                if crate::assumptions::is_real(&x[0]) == Some(true) {
                    return Some(x[0].clone());
                }
            }
            None
        }
        Erfcinv => {
            if z.is_zero() {
                return Some(Expr::oo());
            }
            if z.is_one() {
                return Some(Expr::zero());
            }
            if *z == Expr::int(2) {
                return Some(Expr::neg_oo());
            }
            None
        }
        FresnelS | FresnelC => {
            if is_inf(z) {
                return Some(rat(1, 2));
            }
            if z.is_zero() {
                return Some(Expr::zero());
            }
            let sign = if s == FresnelS { -1 } else { 1 };
            let mut pref = Expr::one();
            let mut arg = z.clone();
            let mut changed = false;
            if could_extract_minus_sign(&arg) {
                pref = neg(&pref);
                arg = neg(&arg);
                changed = true;
            }
            if let Some(t) = extract_i(&arg) {
                pref = mul2(&mul2(&Expr::int(sign), &Expr::i()), &pref);
                arg = t;
                changed = true;
            }
            if changed {
                return Some(mul2(&pref, &sf(s, vec![arg])));
            }
            None
        }
        Ei => {
            if z.is_zero() {
                return Some(Expr::neg_oo());
            }
            if is_inf(z) {
                return Some(Expr::oo());
            }
            if is_neg_inf(z) {
                return Some(Expr::zero());
            }
            let (nz, n) = super::hyperexpand::extract_branch_factor(z, false);
            if n.is_zero() {
                return None;
            }
            Some(crate::ops::add(&[super::ei(&nz), crate::ops::mul(&[Expr::int(2), Expr::i(), Expr::pi(), n])]))
        }
        Li => {
            if z.is_zero() {
                return Some(Expr::zero());
            }
            if z.is_one() {
                return Some(Expr::neg_oo());
            }
            if is_inf(z) {
                return Some(Expr::oo());
            }
            None
        }
        LiOffset => {
            if is_inf(z) {
                return Some(Expr::oo());
            }
            if *z == Expr::int(2) {
                return Some(Expr::zero());
            }
            None
        }
        Si | Shi | Ci | Chi => eval_trig_integral(s, z),
        Expint => eval_expint(&a[0], &a[1]),
        E1 => Some(super::expint(&Expr::one(), z)),
        _ => None,
    }
}

/// Core `erf` head (SymPy's `erf.eval`).
pub(super) fn eval_erf(a: &Expr) -> Option<Expr> {
    if a.is_zero() {
        return Some(Expr::zero());
    }
    if is_inf(a) {
        return Some(Expr::one());
    }
    if is_neg_inf(a) {
        return Some(Expr::minus_one());
    }
    match a.kind() {
        Kind::Func(Func::Special(SpecialFn::Erfinv), x) => return Some(x[0].clone()),
        Kind::Func(Func::Special(SpecialFn::Erfcinv), x) => return Some(sub(&Expr::one(), &x[0])),
        _ => {}
    }
    if let Some(t) = extract_i(a) {
        if is_inf(&t) || is_neg_inf(&t) {
            return Some(a.clone());
        }
    }
    if could_extract_minus_sign(a) {
        return Some(neg(&erf(&neg(a))));
    }
    None
}

/// Core `erfc` head (SymPy's `erfc.eval`).
pub(super) fn eval_erfc(a: &Expr) -> Option<Expr> {
    if a.is_zero() {
        return Some(Expr::one());
    }
    if is_inf(a) {
        return Some(Expr::zero());
    }
    match a.kind() {
        Kind::Func(Func::Special(SpecialFn::Erfinv), x) => return Some(sub(&Expr::one(), &x[0])),
        Kind::Func(Func::Special(SpecialFn::Erfcinv), x) => return Some(x[0].clone()),
        _ => {}
    }
    if let Some(t) = extract_i(a) {
        if is_inf(&t) || is_neg_inf(&t) {
            return Some(neg(a));
        }
    }
    if could_extract_minus_sign(a) {
        return Some(sub(&Expr::int(2), &erfc(&neg(a))));
    }
    None
}

fn eval_erfi(z: &Expr) -> Option<Expr> {
    if z.is_zero() {
        return Some(Expr::zero());
    }
    if is_inf(z) {
        return Some(Expr::oo());
    }
    if let Some(nz) = extract_i(z) {
        if is_inf(&nz) {
            return Some(Expr::i());
        }
        match nz.kind() {
            Kind::Func(Func::Special(SpecialFn::Erfinv), x) => return Some(mul2(&Expr::i(), &x[0])),
            Kind::Func(Func::Special(SpecialFn::Erfcinv), x) => {
                return Some(mul2(&Expr::i(), &sub(&Expr::one(), &x[0])))
            }
            _ => {}
        }
    }
    if could_extract_minus_sign(z) {
        return Some(neg(&super::erfi(&neg(z))));
    }
    None
}

fn eval_erf2(x: &Expr, y: &Expr) -> Option<Expr> {
    if x == y {
        return Some(Expr::zero());
    }
    let special = |e: &Expr| is_inf(e) || is_neg_inf(e) || e.is_zero();
    if special(x) || special(y) {
        return Some(sub(&erf(y), &erf(x)));
    }
    let (sx, sy) = (could_extract_minus_sign(x), could_extract_minus_sign(y));
    if sx && sy {
        return Some(neg(&super::erf2(&neg(x), &neg(y))));
    }
    if sx || sy {
        return Some(sub(&erf(y), &erf(x)));
    }
    None
}

fn eval_trig_integral(s: SpecialFn, z: &Expr) -> Option<Expr> {
    use SpecialFn::*;
    let pi = Expr::pi();
    if z.is_zero() {
        return Some(match s {
            Si | Shi => Expr::zero(),
            _ => Expr::zoo(),
        });
    }
    if is_inf(z) {
        return Some(match s {
            Si => div(&pi, &Expr::int(2)),
            Ci => Expr::zero(),
            _ => Expr::oo(),
        });
    }
    if is_neg_inf(z) {
        return Some(match s {
            Si => neg(&div(&pi, &Expr::int(2))),
            Ci => mul2(&Expr::i(), &pi),
            Shi => Expr::neg_oo(),
            _ => Expr::oo(),
        });
    }
    let exp_polar = |x: Expr| crate::functions::func(crate::Func::ExpPolar, vec![x]);
    let half_pi_i = div(&mul2(&Expr::i(), &pi), &Expr::int(2));
    // _Ifactor(nz, sign)
    let ifactor = |nz: Expr, sign: i64| -> Expr {
        let sg = Expr::int(sign);
        match s {
            Si => crate::ops::mul(&[Expr::i(), sf(Shi, vec![nz]), sg]),
            Shi => crate::ops::mul(&[Expr::i(), sf(Si, vec![nz]), sg]),
            Ci => crate::ops::add(&[sf(Chi, vec![nz]), crate::ops::mul(&[half_pi_i.clone(), sg])]),
            _ => crate::ops::add(&[sf(Ci, vec![nz]), crate::ops::mul(&[half_pi_i.clone(), sg])]),
        }
    };
    if let Some(nz) = extract_exact_factor(z, &exp_polar(half_pi_i.clone())) {
        return Some(ifactor(nz, 1));
    }
    if matches!(s, Si | Shi) {
        if let Some(t) = extract_i(z) {
            return Some(ifactor(t, 1));
        }
    }
    if let Some(nz) = extract_exact_factor(z, &exp_polar(neg(&half_pi_i))) {
        return Some(ifactor(nz, -1));
    }
    // _minusfactor
    let minus = |nz: Expr| -> Expr {
        match s {
            Si | Shi => neg(&sf(s, vec![nz])),
            _ => crate::ops::add(&[sf(s, vec![nz]), mul2(&Expr::i(), &pi)]),
        }
    };
    if let Some(nz) = extract_exact_factor(z, &exp_polar(mul2(&Expr::i(), &pi))) {
        return Some(minus(nz));
    }
    if matches!(s, Si | Shi) && could_extract_minus_sign(z) {
        return Some(minus(neg(z)));
    }
    let (nz, n) = super::hyperexpand::extract_branch_factor(z, false);
    if n.is_zero() && nz == *z {
        return None;
    }
    Some(match s {
        Si | Shi => sf(s, vec![nz]),
        _ => crate::ops::add(&[crate::ops::mul(&[Expr::int(2), pi.clone(), Expr::i(), n]), sf(s, vec![nz])]),
    })
}

/// `z.extract_multiplicatively(c)` for a non-numeric atom-like factor `c`:
/// `z == c` or `c` is a literal factor of the product `z`.
fn extract_exact_factor(z: &Expr, c: &Expr) -> Option<Expr> {
    if z == c {
        return Some(Expr::one());
    }
    if let Kind::Mul(..) = z.kind() {
        let mut args = z.args();
        let i = args.iter().position(|a| a == c)?;
        args.remove(i);
        return Some(crate::ops::mul(&args));
    }
    None
}

fn eval_expint(nu: &Expr, z: &Expr) -> Option<Expr> {
    let nu2 = super::polar::unpolarify(nu, false);
    if nu2 != *nu {
        return Some(super::expint(&nu2, z));
    }
    if let Some(nq) = as_exact(nu) {
        let int_nonpos = nq.is_integer() && !nq.is_positive();
        let half = !nq.is_integer() && nq.mul(&Number::int(2)).is_integer();
        if int_nonpos || half {
            let v = mul2(
                &pow(z, &sub(nu, &Expr::one())),
                &super::uppergamma(&sub(&Expr::one(), nu), z),
            );
            return Some(super::polar::unpolarify(&crate::simp::expand_mul(&v), false));
        }
    }
    // the branch information, as explained in `lowergamma.eval`
    let (z, n) = super::hyperexpand::extract_branch_factor(z, false);
    if n.is_zero() {
        return None;
    }
    let two_pi_i_n = crate::ops::mul(&[Expr::int(2), Expr::pi(), Expr::i(), n]);
    let nu1 = sub(nu, &Expr::one());
    if crate::assumptions::fact(nu, "integer") == Some(true) {
        // `not nu > 0` raises for an undetermined sign in the reference
        if crate::assumptions::fact(nu, "positive") != Some(true) {
            return None;
        }
        let c = div(
            &crate::ops::mul(&[two_pi_i_n, pow(&Expr::minus_one(), &nu1), pow(&super::polar::unpolarify(&z, false), &nu1)]),
            &super::factorial(&nu1),
        );
        return Some(sub(&super::expint(nu, &z), &c));
    }
    let e = crate::functions::exp(&mul2(&two_pi_i_n, nu));
    Some(crate::ops::add(&[
        crate::ops::mul(&[sub(&e, &Expr::one()), pow(&z, &nu1), crate::functions::gamma(&sub(&Expr::one(), nu))]),
        super::expint(nu, &z),
    ]))
}
