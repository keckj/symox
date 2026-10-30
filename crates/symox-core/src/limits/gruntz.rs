//! The Gruntz algorithm for limits at `x -> +oo` (D. Gruntz, *On Computing
//! Limits in a Symbolic Manipulation System*, PhD thesis, ETH Zürich 1996).
//!
//! * `mrv(e)`: the set of most rapidly varying subexpressions of `e`
//!   (comparability classes compared through `lim log(f)/log(g)`),
//! * `rewrite`: express the mrv set through a single `w -> 0+`,
//! * `mrv_leadterm`: leading term `c0 * w**e0` by the series engine (with
//!   `log(w)` known in terms of `x`),
//! * `limitinf`: sign of `e0` decides; recursion on `c0` when `e0 == 0`.
//!
//! If `x` itself is in the mrv set, the problem is "moved up" by `x -> exp(x)`.

use std::cmp::{Ordering, Reverse};

use rustc_hash::FxHashMap;

use super::{unknown, LRes, LimitError};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions as fx;
use crate::ops::{add, div, mul, mul2, mul_factor, neg, pow, sub};
use crate::series::coeff::C;
use crate::series::engine::{depends, sign_of};
use crate::series::{leading_term_with, SErr};
use crate::subs::xreplace;
use crate::traverse::node_count;

const MAX_DEPTH: u32 = 48;

mod abs_power;
mod symbolic_lt;

pub struct Gruntz {
    x: Expr,
    cache: FxHashMap<Expr, Expr>,
    prep: FxHashMap<Expr, Expr>,
    cmp: FxHashMap<(Expr, Expr), Ordering>,
    depth: u32,
    fresh: u32,
}

fn is_inf(e: &Expr) -> bool {
    e.is_constant(Constant::Infinity) || e.is_constant(Constant::NegInfinity)
}

/// `c*oo` with a factor of unknown sign (`oo*sign(b)`).
fn is_parametric_inf(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Mul(_, f) if f.iter().any(|(b, _)| is_inf(b)))
}

impl Gruntz {
    pub fn new(x: &Expr) -> Gruntz {
        Gruntz {
            x: x.clone(),
            cache: FxHashMap::default(),
            prep: FxHashMap::default(),
            cmp: FxHashMap::default(),
            depth: 0,
            fresh: 0,
        }
    }

    fn fresh(&mut self, base: &str) -> Expr {
        self.fresh += 1;
        Expr::symbol(&format!("_symox_{base}{}", self.fresh))
    }

    /// `lim_{x -> oo} e`
    pub fn limitinf(&mut self, e: &Expr) -> LRes<Expr> {
        if !depends(e, &self.x) {
            return Ok(e.clone());
        }
        if *e == self.x {
            return Ok(Expr::oo());
        }
        if let Some(r) = self.cache.get(e) {
            return Ok(r.clone());
        }
        if self.depth > MAX_DEPTH {
            return unknown("recursion too deep");
        }
        self.depth += 1;
        let r = self.limitinf_impl(e);
        self.depth -= 1;
        let r = r?;
        self.cache.insert(e.clone(), r.clone());
        Ok(r)
    }

    fn limitinf_impl(&mut self, e: &Expr) -> LRes<Expr> {
        let e = self.prepare(e)?;
        if !depends(&e, &self.x) {
            return Ok(e);
        }
        let (c0, e0) = self.mrv_leadterm(&e)?;
        match self.sign(&e0)? {
            1 => Ok(Expr::zero()),
            -1 => self.infinity_with_sign_of(&c0),
            _ => {
                // The canonical product splits `exp(a)*exp(b + c)` into one
                // exponential per term, which can undo the rewrite in the mrv
                // variable (the coefficient comes back unchanged): continue
                // with the exponential factors merged into a single one.
                if c0 == e {
                    if let Some(m) = merge_exp_factors(&c0) {
                        return self.limitinf(&m);
                    }
                }
                self.limitinf(&c0)
            }
        }
    }

    fn infinity_with_sign_of(&mut self, c0: &Expr) -> LRes<Expr> {
        // SymPy's limitinf: a coefficient matching I*a (a free of I) gives
        // c0*oo, e.g. I*exp(x) -> oo*I.
        if !depends(c0, &self.x) && imaginary_multiple(c0) {
            return Ok(mul2(c0, &Expr::oo()));
        }
        if imaginary_multiple(c0) {
            // `c0 = I*a` with `a` real-valued and depending on x: the
            // reference `c0*oo` with a positive x is `I*sign(a)*oo`.
            let a = div(c0, &Expr::i());
            if !crate::traverse::has(&a, &Expr::i()) {
                let r = self.infinity_with_sign_of(&a)?;
                return Ok(mul2(&Expr::i(), &r));
            }
        }
        match self.sign(c0) {
            Ok(1) => Ok(Expr::oo()),
            Ok(-1) => Ok(Expr::neg_oo()),
            Ok(_) => unknown("zero leading coefficient"),
            Err(err) if depends(c0, &self.x) => {
                // The reference `sign(c0, x)` of a product multiplies the
                // signs of its factors, and an x-free factor of unknown sign
                // contributes the symbolic `sign(a)` (`y/(x*log(x))` at 0+
                // is `-oo*sign(y)`).
                if let Kind::Mul(c, facs) = c0.kind() {
                    let (free, dep): (Vec<_>, Vec<_>) = facs.iter().partition(|(b, p)| !depends(b, &self.x) && !depends(p, &self.x));
                    let free: Vec<Expr> = free.into_iter().map(|(b, p)| mul_factor(b, p)).collect();
                    if !free.is_empty() && !dep.is_empty() {
                        // `sign(a)` is the direction `a/|a|` for a complex
                        // `a`; one symbolic sign per factor, as the reference
                        // recursion through `Mul.as_two_terms`.
                        let mut rest = vec![Expr::num(c.clone())];
                        rest.extend(dep.into_iter().map(|(b, p)| mul_factor(b, p)));
                        let mut s = self.sign(&mul(&rest))?;
                        let mut signs = Vec::new();
                        for a in &free {
                            match self.sign(a) {
                                Ok(t) => s *= t,
                                Err(_) => signs.push(fx::sign(a)),
                            }
                        }
                        if s == 0 {
                            return unknown("zero leading coefficient");
                        }
                        if signs.is_empty() {
                            return Ok(if s > 0 { Expr::oo() } else { Expr::neg_oo() });
                        }
                        signs.push(Expr::int(s as i64));
                        signs.push(Expr::oo());
                        return Ok(mul(&signs));
                    }
                }
                Err(err)
            }
            Err(err) => {
                let nonreal = || crate::evalf::eval_c64(c0).is_some_and(|z| z.im.abs() > 1e-12 * (1.0 + z.re.abs()));
                if !depends(c0, &self.x)
                    && ((!crate::functions::is_numeric_value(c0) && crate::evalf::eval_c64(c0).is_none()) || nonreal())
                {
                    // parameter: oo*sign(c0), like SymPy
                    Ok(mul2(&fx::sign(c0), &Expr::oo()))
                } else {
                    Err(err)
                }
            }
        }
    }

    /// Sign of `e` for large `x`.
    pub fn sign(&mut self, e: &Expr) -> LRes<i32> {
        if !depends(e, &self.x) {
            if super::exact_sign::vanishes_after_log_expansion(e) {
                return Ok(0);
            }
            return match sign_of(&C::from_expr(e.clone())) {
                Some(s) if s != 0 => Ok(s),
                Some(_) => Ok(0),
                // the reference `sign(e, x)` asks the assumptions first
                None if crate::assumptions::is_positive(e) == Some(true) => Ok(1),
                None if crate::assumptions::is_negative(e) == Some(true) => Ok(-1),
                None if crate::assumptions::is_zero(e) == Some(true) => Ok(0),
                None => {
                    // An undecided sum of fractions (`-1 + (n - 5)/n`, the
                    // exponent difference of `x**n - x**(n - 5)`) may have a
                    // decided sign once over a common denominator (`-5/n`).
                    let t = crate::simplify::together(e);
                    if t != *e && matches!(e.kind(), Kind::Add(..)) {
                        if let Some(s) = sign_of(&C::from_expr(t.clone())) {
                            return Ok(s);
                        }
                        if crate::assumptions::is_positive(&t) == Some(true) {
                            return Ok(1);
                        }
                        if crate::assumptions::is_negative(&t) == Some(true) {
                            return Ok(-1);
                        }
                    }
                    if let Some(s) = super::exact_sign::sign_after_log_expansion(e) {
                        return Ok(s);
                    }
                    unknown(format!("sign of {e}"))
                }
            };
        }
        if *e == self.x {
            return Ok(1);
        }
        match e.kind() {
            Kind::Mul(c, facs) => {
                let mut s = if c.is_negative() { -1 } else { 1 };
                for (b, p) in facs {
                    s *= self.sign(&mul_factor(b, p))?;
                }
                return Ok(s);
            }
            Kind::Pow(b, p) => {
                if let Some(k) = p.as_number() {
                    if k.is_integer() {
                        let even = k.to_i64().map(|v| v % 2 == 0).unwrap_or(false);
                        if even {
                            return Ok(1);
                        }
                        return self.sign(b);
                    }
                    if self.sign(b)? == 1 {
                        return Ok(1);
                    }
                }
            }
            Kind::Func(Func::Exp, _) => return Ok(1),
            _ => {}
        }
        // The leading coefficient may still depend on x (and, through the
        // gamma rewrites, reproduce itself): bound the recursion instead of
        // overflowing the stack.
        if self.depth > MAX_DEPTH {
            return unknown("recursion too deep");
        }
        self.depth += 1;
        let r = self.sign_of_leading(e);
        self.depth -= 1;
        r
    }

    fn sign_of_leading(&mut self, e: &Expr) -> LRes<i32> {
        let pe = self.prepare(e)?;
        let (c0, _) = self.mrv_leadterm(&pe)?;
        if c0 == *e || c0 == pe {
            // The canonical product splits `exp(a - b)` into `exp(a)*exp(-b)`,
            // which can undo the rewrite in the mrv variable (as in
            // `limitinf`): retry with the exponential factors of every
            // product merged.
            let merged = merge_exp_deep(&c0);
            if merged != c0 {
                return self.sign(&merged);
            }
            return unknown(format!("sign of {e}"));
        }
        self.sign(&c0)
    }

    // ------------------------------------------------------------------
    // rewriting into tractable form
    // ------------------------------------------------------------------

    pub fn prepare(&mut self, e: &Expr) -> LRes<Expr> {
        if !depends(e, &self.x) {
            return Ok(e.clone());
        }
        if let Some(r) = self.prep.get(e) {
            return Ok(r.clone());
        }
        if let Kind::Piecewise(pieces) = e.kind() {
            // The branch that holds on a whole neighbourhood of `x = oo`
            // (the one-sided neighbourhood of the original point): each
            // condition is decided by the eventual sign of `lhs - rhs`.
            // An undecided condition before the selected branch declines.
            for (value, cond) in pieces {
                match self.eventually_true(cond)? {
                    Some(true) => {
                        let r = self.prepare(value)?;
                        self.prep.insert(e.clone(), r.clone());
                        return Ok(r);
                    }
                    Some(false) => continue,
                    None => return unknown(format!("condition {cond} near the limit point")),
                }
            }
            return unknown("no Piecewise branch holds near the limit point");
        }
        if let Some(r) = self.abs_of_real_power(e) {
            let r = self.prepare(&r)?;
            self.prep.insert(e.clone(), r.clone());
            return Ok(r);
        }
        if let Kind::Mul(c, facs) = e.kind() {
            // factors (not bases and exponents separately) must be prepared
            let mut v = vec![Expr::num(c.clone())];
            for (b, p) in facs {
                v.push(self.prepare(&mul_factor(b, p))?);
            }
            let r = mul(&v);
            self.prep.insert(e.clone(), r.clone());
            return Ok(r);
        }
        let mut err: Option<LimitError> = None;
        let x = self.x.clone();
        let e1 = {
            let mut f = |c: &Expr| -> Expr {
                if err.is_some() || !depends(c, &x) {
                    return c.clone();
                }
                match self.prepare(c) {
                    Ok(v) => v,
                    Err(er) => {
                        err = Some(er);
                        c.clone()
                    }
                }
            };
            crate::traverse::map_children(e, &mut f)
        };
        if let Some(er) = err {
            return Err(er);
        }
        let r = self.prepare_node(&e1)?;
        self.prep.insert(e.clone(), r.clone());
        Ok(r)
    }

    /// Truth of a Piecewise condition for all large `x`: `Some(b)` when
    /// proven, `None` when undecided (a condition on other symbols with an
    /// unknown truth value, a relational between non-real values, ...).
    fn eventually_true(&mut self, c: &Expr) -> LRes<Option<bool>> {
        use crate::expr::RelOp;
        use crate::logic::LogicOp;
        Ok(match c.kind() {
            Kind::Bool(b) => Some(*b),
            Kind::Rel(op, l, r) => {
                if !depends(l, &self.x) && !depends(r, &self.x) {
                    return Ok(None);
                }
                let d = sub(l, r);
                if crate::traverse::has(&d, &Expr::i()) {
                    return Ok(None);
                }
                let s = match self.sign(&d) {
                    Ok(s) if s != 0 => s,
                    Ok(_) => return Ok(None),
                    Err(LimitError::Unknown(_)) => return Ok(None),
                    Err(err) => return Err(err),
                };
                Some(match op {
                    RelOp::Lt | RelOp::Le => s < 0,
                    RelOp::Gt | RelOp::Ge => s > 0,
                    RelOp::Eq => false,
                    RelOp::Ne => true,
                })
            }
            Kind::Func(Func::Logic(op), args) => match op {
                LogicOp::Not if args.len() == 1 => self.eventually_true(&args[0])?.map(|b| !b),
                LogicOp::And | LogicOp::Or => {
                    let and = *op == LogicOp::And;
                    let mut undecided = false;
                    for a in args {
                        match self.eventually_true(a)? {
                            Some(b) if b != and => return Ok(Some(b)),
                            Some(_) => {}
                            None => undecided = true,
                        }
                    }
                    if undecided {
                        None
                    } else {
                        Some(and)
                    }
                }
                _ => None,
            },
            _ => None,
        })
    }

    fn prepare_node(&mut self, e: &Expr) -> LRes<Expr> {
        let x = self.x.clone();
        let two = Expr::int(2);
        if let Some(r) = self.infinite_argument_tractable(e) {
            return self.prepare(&r);
        }
        Ok(match e.kind() {
            // `x -> oo` is a positive variable in the reference gruntz, so a
            // Heaviside or DiracDelta whose argument has an eventual sign
            // evaluates (`Heaviside(1/x)` is `1`); `H0` is never reached.
            Kind::Func(f @ (Func::Heaviside | Func::DiracDelta), args) if !args.is_empty() && depends(&args[0], &x) => {
                match self.sign(&args[0]) {
                    Ok(s) if s != 0 && *f == Func::DiracDelta => Expr::zero(),
                    Ok(1) => Expr::one(),
                    Ok(-1) => Expr::zero(),
                    Ok(_) => return unknown(format!("sign of {} near the limit point", args[0])),
                    Err(err) => return Err(err),
                }
            }
            Kind::Pow(b, p) if depends(p, &x) || (depends(b, &x) && !p.is_number()) => {
                fx::exp(&mul2(p, &simp_log(&fx::log(b))))
            }
            Kind::Func(Func::Gamma, args) if args.len() == 1 && depends(&args[0], &x) => self.gamma_tractable(&args[0], e)?,
            Kind::Func(Func::LogGamma, args) if args.len() == 1 && depends(&args[0], &x) => {
                // Stirling's series needs `+oo`; a finite limit point is regular.
                let l = self.limitinf(&args[0])?;
                if !l.is_constant(Constant::Infinity) && !is_finite_limit(&l) {
                    return unknown(format!("loggamma({}) near the limit point", args[0]));
                }
                e.clone()
            }
            Kind::Func(Func::Special(s), args) if args.iter().any(|a| depends(a, &x)) => {
                match self.special_tractable(*s, args)? {
                    Some(r) => self.prepare(&r)?,
                    None => e.clone(),
                }
            }
            Kind::Func(f, args) if args.len() == 1 && depends(&args[0], &x) => {
                let a = &args[0];
                match f {
                    Func::Sinh => div(&sub(&fx::exp(a), &fx::exp(&neg(a))), &two),
                    Func::Cosh => div(&add(&[fx::exp(a), fx::exp(&neg(a))]), &two),
                    Func::Tanh => {
                        let e2 = fx::exp(&mul2(&two, a));
                        div(&sub(&e2, &Expr::one()), &add(&[e2, Expr::one()]))
                    }
                    Func::Coth => {
                        let e2 = fx::exp(&mul2(&two, a));
                        div(&add(&[e2.clone(), Expr::one()]), &sub(&e2, &Expr::one()))
                    }
                    Func::Log => simp_log(e),
                    Func::Asinh => fx::log(&add(&[a.clone(), crate::ops::sqrt(&add(&[pow(a, &two), Expr::one()]))])),
                    Func::Acosh => fx::log(&add(&[a.clone(), crate::ops::sqrt(&sub(&pow(a, &two), &Expr::one()))])),
                    Func::Atanh => div(
                        &sub(&fx::log(&add(&[Expr::one(), a.clone()])), &fx::log(&sub(&Expr::one(), a))),
                        &two,
                    ),
                    Func::Atan | Func::Acot => {
                        let l = self.limitinf(a)?;
                        if is_inf(&l) {
                            let inv = pow(a, &Expr::minus_one());
                            if *f == Func::Acot {
                                fx::atan(&inv)
                            } else {
                                let s = if l.is_constant(Constant::Infinity) { 1 } else { -1 };
                                sub(&mul(&[Expr::rational(s, 2), Expr::pi()]), &fx::atan(&inv))
                            }
                        } else if is_parametric_inf(&l) {
                            // `a*x` with a parameter of unknown sign (or
                            // complex): the reference tractable form of atan,
                            // `I/2*(log(1 - I*u) - log(1 + I*u))`, an identity
                            // of the principal branches; acot(u) = atan(1/u).
                            let r = if *f == Func::Acot {
                                fx::atan(&pow(a, &Expr::minus_one()))
                            } else {
                                let iu = mul2(&Expr::i(), a);
                                let lg = sub(&fx::log(&sub(&Expr::one(), &iu)), &fx::log(&add(&[Expr::one(), iu])));
                                mul(&[Expr::rational(1, 2), Expr::i(), lg])
                            };
                            return self.prepare(&r);
                        } else {
                            e.clone()
                        }
                    }
                    Func::Abs => {
                        let s = self.sign(a)?;
                        mul2(&Expr::int(s as i64), a)
                    }
                    Func::Sign => Expr::int(self.sign(a)? as i64),
                    _ => e.clone(),
                }
            }
            _ => e.clone(),
        })
    }

    /// Tractable forms of the error functions, exponential, logarithmic,
    /// trigonometric and Fresnel integrals and of the incomplete gamma
    /// functions whose argument tends to `±oo` (`li` also at `0+`): the
    /// reference `_eval_rewrite_as_tractable` (`erf`, `erfc`, `erfi`, `Ei`,
    /// `li`, `Li`, `uppergamma`, with `_erfs` and `_eis`) and the
    /// `_eval_aseries` of `Si`, `Ci`, `expint`, `fresnels` and `fresnelc`,
    /// written with the auxiliary functions of [`crate::series::asym`]
    /// that the series engine expands in powers of `1/z`. `Shi` and `Chi`
    /// are `(Ei(z) -+ Ei(-z))/2` (plus `I*pi` for `Chi` at `-oo`).
    fn infinite_argument_tractable(&mut self, e: &Expr) -> Option<Expr> {
        use crate::series::asym::Asym;
        use crate::special::SpecialFn as S;
        let Kind::Func(f, args) = e.kind() else { return None };
        let x = self.x.clone();
        let z = args.last()?;
        if !depends(z, &x) || args[..args.len() - 1].iter().any(|a| depends(a, &x)) {
            return None;
        }
        let interesting = matches!(f, Func::Erf | Func::Erfc)
            || matches!(f, Func::Special(
                S::Erfi | S::Ei | S::Li | S::LiOffset | S::Si | S::Ci | S::Shi | S::Chi | S::Expint
                    | S::E1 | S::UpperGamma | S::LowerGamma | S::FresnelS | S::FresnelC
                    | S::BesselJ | S::BesselY | S::BesselI | S::BesselK));
        if !interesting {
            return None;
        }
        let l = self.limitinf(z).ok()?;
        let s: i64 = if l.is_constant(Constant::Infinity) {
            1
        } else if l.is_constant(Constant::NegInfinity) {
            -1
        } else {
            0
        };
        let one = Expr::one();
        let two = Expr::int(2);
        let half = Expr::rational(1, 2);
        let ipi = mul2(&Expr::i(), &Expr::pi());
        let eis = |u: &Expr| mul2(&fx::exp(u), &Asym::Eis.apply(vec![u.clone()]));
        let mz = neg(z);
        // `uppergamma(a, z) = z**a*exp(-z)*expints(1 - a, z)`
        let upper = |a: &Expr| {
            let nu = sub(&one, a);
            mul(&[pow(z, a), fx::exp(&mz), Asym::Expints.apply(vec![nu, z.clone()])])
        };
        if s == 0 {
            // li(z) = z*_eis(log(z)) for z -> 0+
            if matches!(f, Func::Special(S::Li | S::LiOffset)) && l.is_zero() && self.sign(z).ok()? == 1 {
                let r = mul2(z, &Asym::Eis.apply(vec![fx::log(z)]));
                return Some(if matches!(f, Func::Special(S::LiOffset)) {
                    sub(&r, &crate::special::apply(S::Li, vec![two]))
                } else {
                    r
                });
            }
            return None;
        }
        if let Func::Special(b @ (S::BesselJ | S::BesselY | S::BesselI | S::BesselK)) = f {
            return self.bessel_tractable(*b, &args[0], z, s);
        }
        let sz = Expr::int(s);
        Some(match f {
            Func::Erf | Func::Erfc => {
                let z2 = pow(z, &two);
                // erf(z) = s - _erfs(s*z)*exp(-z**2)
                let tail = mul2(&Asym::Erfs.apply(vec![mul2(&sz, z)]), &fx::exp(&neg(&z2)));
                let erf = sub(&sz, &mul2(&sz, &tail));
                if *f == Func::Erf {
                    erf
                } else {
                    sub(&one, &erf)
                }
            }
            Func::Special(S::Erfi) => mul2(&Asym::Erfis.apply(vec![z.clone()]), &fx::exp(&pow(z, &two))),
            Func::Special(S::Ei) => eis(z),
            Func::Special(S::Li | S::LiOffset) if s == 1 => {
                let r = mul2(z, &Asym::Eis.apply(vec![fx::log(z)]));
                if matches!(f, Func::Special(S::LiOffset)) {
                    sub(&r, &crate::special::apply(S::Li, vec![two]))
                } else {
                    r
                }
            }
            Func::Special(S::Si) => {
                let fz = Asym::SiF.apply(vec![z.clone()]);
                let gz = Asym::SiG.apply(vec![z.clone()]);
                add(&[
                    mul(&[sz, half, Expr::pi()]),
                    neg(&mul2(&fx::cos(z), &fz)),
                    neg(&mul2(&fx::sin(z), &gz)),
                ])
            }
            Func::Special(S::Ci) => {
                let fz = Asym::SiF.apply(vec![z.clone()]);
                let gz = Asym::SiG.apply(vec![z.clone()]);
                let r = sub(&mul2(&fx::sin(z), &fz), &mul2(&fx::cos(z), &gz));
                if s == 1 {
                    r
                } else {
                    add(&[r, ipi])
                }
            }
            Func::Special(S::Shi) => div(&sub(&eis(z), &eis(&mz)), &two),
            Func::Special(S::Chi) => {
                let r = div(&add(&[eis(z), eis(&mz)]), &two);
                if s == 1 {
                    r
                } else {
                    add(&[r, ipi])
                }
            }
            Func::Special(S::Expint) if s == 1 => {
                mul2(&fx::exp(&mz), &Asym::Expints.apply(vec![args[0].clone(), z.clone()]))
            }
            Func::Special(S::E1) if s == 1 => mul2(&fx::exp(&mz), &Asym::Expints.apply(vec![one, z.clone()])),
            Func::Special(S::UpperGamma) if s == 1 => upper(&args[0]),
            Func::Special(S::LowerGamma) if s == 1 => {
                let g = crate::functions::func(Func::Gamma, vec![args[0].clone()]);
                if !is_finite_limit(&g) {
                    return None;
                }
                sub(&g, &upper(&args[0]))
            }
            Func::Special(S::FresnelS | S::FresnelC) => {
                let ph = mul(&[half.clone(), Expr::pi(), pow(z, &two)]);
                let fz = Asym::FresF.apply(vec![z.clone()]);
                let gz = Asym::FresG.apply(vec![z.clone()]);
                let base = mul2(&sz, &half);
                if matches!(f, Func::Special(S::FresnelS)) {
                    add(&[base, neg(&mul2(&fx::cos(&ph), &fz)), neg(&mul2(&fx::sin(&ph), &gz))])
                } else {
                    add(&[base, mul2(&fx::sin(&ph), &fz), neg(&mul2(&fx::cos(&ph), &gz))])
                }
            }
            _ => return None,
        })
    }

    /// Hankel's asymptotic forms of the Bessel functions at a real infinite
    /// argument (the leading terms of the reference
    /// `_eval_as_leading_term`, with the full asymptotic series in the
    /// auxiliary functions): `J` and `Y` oscillate with an amplitude
    /// `sqrt(2/(pi z))`, `I ~ exp(z)/sqrt(2 pi z)`, `K ~ sqrt(pi/(2 z)) exp(-z)`.
    /// At `-oo` an integer order uses `J(nu, -u) = (-1)**nu J(nu, u)`,
    /// `I(nu, -u) = (-1)**nu I(nu, u)`, `Y(nu, -u) = (-1)**nu (Y(nu, u) +
    /// 2 i J(nu, u))` and `K(nu, -u) = (-1)**nu K(nu, u) - i pi I(nu, u)`
    /// (principal branch, `arg(-u) = pi`).
    fn bessel_tractable(&mut self, b: crate::special::SpecialFn, nu: &Expr, z: &Expr, s: i64) -> Option<Expr> {
        use crate::series::asym::Asym;
        use crate::special::SpecialFn as S;
        let one = Expr::one();
        let half = Expr::rational(1, 2);
        let form = |b: S, u: &Expr| -> Expr {
            let aux = |k: Asym| k.apply(vec![nu.clone(), u.clone()]);
            match b {
                S::BesselJ | S::BesselY => {
                    let w = sub(u, &mul(&[add(&[mul2(&Expr::int(2), nu), one.clone()]), Expr::pi(), Expr::rational(1, 4)]));
                    let pref = mul2(&crate::ops::sqrt(&Expr::int(2)), &pow(&mul2(&Expr::pi(), u), &neg(&half)));
                    let p = add(&[one.clone(), aux(Asym::BesP)]);
                    let q = aux(Asym::BesQ);
                    let r = if b == S::BesselJ {
                        sub(&mul2(&fx::cos(&w), &p), &mul2(&fx::sin(&w), &q))
                    } else {
                        add(&[mul2(&fx::sin(&w), &p), mul2(&fx::cos(&w), &q)])
                    };
                    mul2(&pref, &r)
                }
                S::BesselI => mul(&[
                    fx::exp(u),
                    pow(&mul(&[Expr::int(2), Expr::pi(), u.clone()]), &neg(&half)),
                    add(&[one.clone(), aux(Asym::BesI)]),
                ]),
                _ => mul(&[
                    crate::ops::sqrt(&div(&Expr::pi(), &Expr::int(2))),
                    pow(u, &neg(&half)),
                    fx::exp(&neg(u)),
                    add(&[one.clone(), aux(Asym::BesK)]),
                ]),
            }
        };
        if s == 1 {
            return Some(form(b, z));
        }
        if !nu.as_number().is_some_and(|k| k.is_integer()) {
            return None;
        }
        let u = neg(z);
        let sg = pow(&Expr::minus_one(), nu);
        Some(match b {
            S::BesselJ | S::BesselI => mul2(&sg, &form(b, &u)),
            S::BesselY => mul2(&sg, &add(&[form(S::BesselY, &u), mul(&[Expr::int(2), Expr::i(), form(S::BesselJ, &u)])])),
            _ => sub(&mul2(&sg, &form(S::BesselK, &u)), &mul(&[Expr::i(), Expr::pi(), form(S::BesselI, &u)])),
        })
    }

    /// The reference `gamma._eval_rewrite_as_tractable`: `exp(loggamma(z))`,
    /// whose exponent has Stirling's series at `+oo`. A gamma whose argument
    /// has a finite limit is continuous or has a pole the series engine
    /// expands; any other argument is left undecided.
    fn gamma_tractable(&mut self, a: &Expr, e: &Expr) -> LRes<Expr> {
        let l = self.limitinf(a)?;
        if l.is_constant(Constant::Infinity) {
            return Ok(fx::exp(&crate::functions::func(Func::LogGamma, vec![a.clone()])));
        }
        if is_finite_limit(&l) {
            return Ok(e.clone());
        }
        unknown(format!("gamma({a}) near the limit point"))
    }

    /// Tractable forms of the factorial family and of harmonic numbers
    /// (the reference rewrites `factorial` to `gamma` before gruntz and
    /// `binomial`, `RisingFactorial`, `FallingFactorial` and `harmonic`
    /// through `rewrite('tractable')`). `None` keeps the node.
    fn special_tractable(&mut self, s: crate::special::SpecialFn, args: &[Expr]) -> LRes<Option<Expr>> {
        use crate::special::SpecialFn as S;
        let gamma = |a: Expr| crate::functions::func(Func::Gamma, vec![a]);
        let one = Expr::one();
        Ok(Some(match (s, args) {
            (S::Factorial, [a]) => gamma(add(&[a.clone(), one])),
            (S::Binomial, [n, k]) => div(
                &gamma(add(&[n.clone(), one.clone()])),
                &mul2(&gamma(add(&[k.clone(), one.clone()])), &gamma(add(&[n.clone(), neg(k), one]))),
            ),
            (S::RisingFactorial, [a, k]) => div(&gamma(add(&[a.clone(), k.clone()])), &gamma(a.clone())),
            (S::FallingFactorial, [a, k]) => {
                div(&gamma(add(&[a.clone(), one.clone()])), &gamma(add(&[a.clone(), neg(k), one])))
            }
            (S::LowerGamma | S::UpperGamma, [a, z]) if depends(a, &self.x) && !depends(z, &self.x) => {
                // The reference `lowergamma._eval_aseries` in the order
                // (`z**s*exp(-z)*Σ z**k/rf(s, k + 1)`) and
                // `uppergamma = exp(loggamma(s)) - lowergamma(s, z)`, for an
                // order tending to `+oo` and a positive fixed `z` (so that
                // `z**s` is a real exponential).
                if !self.limitinf(a)?.is_constant(Constant::Infinity) || self.sign(z)? != 1 {
                    return Ok(None);
                }
                let lower = mul(&[
                    pow(z, a),
                    fx::exp(&neg(z)),
                    crate::series::asym::Asym::LowerS.apply(vec![z.clone(), a.clone()]),
                ]);
                if s == S::LowerGamma {
                    lower
                } else {
                    sub(&fx::exp(&crate::functions::func(Func::LogGamma, vec![a.clone()])), &lower)
                }
            }
            (S::Harmonic, [n]) => return self.harmonic_tractable(n, &one),
            (S::Harmonic, [n, m]) => return self.harmonic_tractable(n, m),
            (S::Fibonacci, [n]) => {
                // (phi**n - cos(pi*n)/phi**n)/sqrt(5)
                let phi = Expr::constant(Constant::GoldenRatio);
                let pn = pow(&phi, n);
                let osc = div(&fx::cos(&mul2(&Expr::pi(), n)), &pn);
                div(&sub(&pn, &osc), &crate::ops::sqrt(&Expr::int(5)))
            }
            _ => return Ok(None),
        }))
    }

    /// `harmonic(n, m) = polygamma(0, n + 1) + EulerGamma` for `m = 1` and
    /// `(-1)**m*(polygamma(m - 1, 1) - polygamma(m - 1, n + 1))/gamma(m)`
    /// for an integer `m > 1` (the reference `_eval_rewrite_as_polygamma`).
    fn harmonic_tractable(&mut self, n: &Expr, m: &Expr) -> LRes<Option<Expr>> {
        use crate::special::{apply, SpecialFn as S};
        let l = self.limitinf(n)?;
        if is_finite_limit(&l) && !depends(m, &self.x) {
            return Ok(None);
        }
        let order = match m.as_number().filter(|k| k.is_integer()).and_then(|k| k.to_i64()) {
            Some(k) if k >= 1 && l.is_constant(Constant::Infinity) => k,
            _ => return unknown(format!("harmonic({n}, {m}) near the limit point")),
        };
        let n1 = add(&[n.clone(), Expr::one()]);
        let pg = |k: i64, z: Expr| apply(S::Polygamma, vec![Expr::int(k), z]);
        Ok(Some(if order == 1 {
            add(&[pg(0, n1), Expr::constant(Constant::EulerGamma)])
        } else {
            let sgn = if order % 2 == 0 { 1 } else { -1 };
            let diff = sub(&pg(order - 1, Expr::one()), &pg(order - 1, n1));
            div(&mul2(&Expr::int(sgn), &diff), &crate::functions::func(Func::Gamma, vec![m.clone()]))
        }))
    }

    // ------------------------------------------------------------------
    // mrv sets
    // ------------------------------------------------------------------

    fn mrv(&mut self, e: &Expr) -> LRes<Vec<Expr>> {
        let x = self.x.clone();
        if !depends(e, &x) {
            return Ok(vec![]);
        }
        if *e == x {
            return Ok(vec![x]);
        }
        match e.kind() {
            Kind::Add(_, terms) => {
                let mut acc = vec![];
                for (t, _) in terms {
                    let m = self.mrv(t)?;
                    acc = self.mrv_max(acc, m)?;
                }
                Ok(acc)
            }
            Kind::Mul(_, facs) => {
                let mut acc = vec![];
                for (b, p) in facs {
                    let m = self.mrv(&mul_factor(b, p))?;
                    acc = self.mrv_max(acc, m)?;
                }
                Ok(acc)
            }
            Kind::Pow(b, p) => {
                if depends(p, &x) {
                    let e2 = fx::exp(&mul2(p, &fx::log(b)));
                    return self.mrv(&e2);
                }
                self.mrv(b)
            }
            Kind::Func(Func::Exp, a) => {
                let l = self.limitinf(&a[0])?;
                let inner = self.mrv(&a[0])?;
                // `oo*sign(b)` (a parameter of unknown sign) is infinite too,
                // like the reference `any(_.is_infinite for _ in
                // Mul.make_args(limitinf(e.exp, x)))`.
                if is_inf(&l) || is_parametric_inf(&l) {
                    self.mrv_max(vec![e.clone()], inner)
                } else {
                    Ok(inner)
                }
            }
            // SymPy's mrv knows functions only: a held sum, product,
            // integral, limit or substitution that depends on x raises
            // NotImplementedError (its growth is not that of its arguments).
            Kind::Func(Func::Sum | Func::Product | Func::Integral | Func::Limit | Func::Subs | Func::Order, _) => {
                unknown(format!("mrv of {e}"))
            }
            Kind::Func(f, args) => {
                if grows_beyond_its_argument(f) {
                    // The growth of these functions at an infinite argument is
                    // not that of the argument: without a tractable rewrite
                    // the mrv set is unknown.
                    for a in args {
                        if depends(a, &x) && !is_finite_limit(&self.limitinf(a)?) {
                            return unknown(format!("mrv of {e}"));
                        }
                    }
                }
                let mut acc = vec![];
                for a in args {
                    let m = self.mrv(a)?;
                    acc = self.mrv_max(acc, m)?;
                }
                Ok(acc)
            }
            _ => unknown(format!("mrv of {e}")),
        }
    }

    fn mrv_max(&mut self, f: Vec<Expr>, g: Vec<Expr>) -> LRes<Vec<Expr>> {
        if f.is_empty() {
            return Ok(g);
        }
        if g.is_empty() {
            return Ok(f);
        }
        let union = |mut f: Vec<Expr>, g: Vec<Expr>| {
            for a in g {
                if !f.contains(&a) {
                    f.push(a);
                }
            }
            f
        };
        if f.iter().any(|a| g.contains(a)) {
            return Ok(union(f, g));
        }
        match self.compare(&f[0], &g[0])? {
            Ordering::Greater => Ok(f),
            Ordering::Less => Ok(g),
            Ordering::Equal => Ok(union(f, g)),
        }
    }

    fn log_of(&self, a: &Expr) -> Expr {
        match a.kind() {
            Kind::Func(Func::Exp, args) => args[0].clone(),
            _ => fx::log(a),
        }
    }

    fn compare(&mut self, a: &Expr, b: &Expr) -> LRes<Ordering> {
        if let Some(o) = self.cmp.get(&(a.clone(), b.clone())) {
            return Ok(*o);
        }
        let r = self.compare_impl(a, b)?;
        self.cmp.insert((a.clone(), b.clone()), r);
        self.cmp.insert((b.clone(), a.clone()), r.reverse());
        Ok(r)
    }

    fn compare_impl(&mut self, a: &Expr, b: &Expr) -> LRes<Ordering> {
        let la = self.log_of(a);
        let lb = self.log_of(b);
        // A nonzero factor free of x (a parameter of `exp(b*sqrt(x))`)
        // does not change the order of the ratio: only the x-dependent
        // part is passed to the limit, like the numeric case.
        let q = div(&la, &lb);
        let q = match q.kind() {
            Kind::Mul(_, facs) if facs.iter().any(|(f, p)| !depends(f, &self.x) && !depends(p, &self.x)) => {
                let dep: Vec<Expr> = facs
                    .iter()
                    .filter(|(f, p)| depends(f, &self.x) || depends(p, &self.x))
                    .map(|(f, p)| mul_factor(f, p))
                    .collect();
                mul(&dep)
            }
            _ => q,
        };
        let c = self.limitinf(&q)?;
        if c.is_zero() {
            Ok(Ordering::Less)
        } else if is_inf(&c) || is_parametric_inf(&c) {
            // `oo*sign(b)`: the ratio of the logarithms is unbounded
            // whatever the sign of the parameter.
            Ok(Ordering::Greater)
        } else {
            Ok(Ordering::Equal)
        }
    }

    // ------------------------------------------------------------------
    // leading term
    // ------------------------------------------------------------------

    fn move_up(&self, e: &Expr) -> Expr {
        let mut m = FxHashMap::default();
        m.insert(self.x.clone(), fx::exp(&self.x));
        simp_logs(&xreplace(e, &m))
    }

    /// `(c0, e0)` with `e ~ c0 * w**e0`, `w -> 0+` the mrv variable.
    pub fn mrv_leadterm(&mut self, e: &Expr) -> LRes<(Expr, Expr)> {
        let x = self.x.clone();
        if !depends(e, &x) {
            return Ok((e.clone(), Expr::zero()));
        }
        // Reference Mul.flatten gathers `exp(a)*exp(b)` into `exp(a + b)`;
        // the canonical product here keeps one exponential per factor, so
        // reciprocal factors such as `x**(-x - 1)*x**(x + 1)` (rewritten as
        // `exp((-x - 1)*log(x))*exp((x + 1)*log(x))`) would stay unmerged.
        // Products of exponentials whose expanded exponents cancel are
        // merged first.
        let merged = merge_cancelling_exps(e, &x);
        let e = &merged;
        let mut omega = self.mrv(e)?;
        if omega.is_empty() {
            return Ok((e.clone(), Expr::zero()));
        }
        let mut e = e.clone();
        if omega.len() == 1 && omega[0] == x {
            // moving up gives Omega = {exp(x)}, w = exp(-x): directly
            // substitute x -> 1/w with log(w) = -x
            let w = self.fresh("w");
            let mut m = FxHashMap::default();
            m.insert(x.clone(), pow(&w, &Expr::minus_one()));
            let f = xreplace(&e, &m);
            return self.leadterm_in(&f, &w, &neg(&x));
        }
        if omega.contains(&x) {
            // move the problem up (x -> exp(x)); the limit of the leading
            // coefficient is unchanged, so there is no need to move down.
            // The mrv set is moved up too: recomputing it would recurse.
            let up = self.move_up(&e);
            let up = self.prepare(&up)?;
            omega = self.mrv(&up)?;
            if omega.contains(&x) {
                return unknown("cannot move up");
            }
            e = up;
            if omega.is_empty() {
                return Ok((e, Expr::zero()));
            }
        }
        let w = self.fresh("w");
        let (f, logw) = self.rewrite(&e, omega, &w)?;
        self.leadterm_in(&f, &w, &logw)
    }

    fn leadterm_in(&mut self, f: &Expr, w: &Expr, logw: &Expr) -> LRes<(Expr, Expr)> {
        if let Some(r) = self.symbolic_power_leadterm(f, w, logw) {
            return r;
        }
        if let Some(g) = self.split_symbolic_bases(f, w, logw) {
            if let Some(r) = self.symbolic_power_leadterm(&g, w, logw) {
                if r.is_ok() {
                    return r;
                }
            }
        }
        let l = self.fresh("logw");
        let plain = self.lt_cancel(f, w, &l, logw);
        let x = self.x.clone();
        let has_x_exp = |c: &Expr| {
            fn walk(c: &Expr, x: &Expr) -> bool {
                match c.kind() {
                    Kind::Func(Func::Exp, a) if depends(&a[0], x) => true,
                    _ => c.args().iter().any(|a| walk(a, x)),
                }
            }
            walk(c, &x)
        };
        let plain_ok = matches!(&plain, Ok(Some((c, _))) if !c.is_zero() && !has_x_exp(c));
        // With `log(w) = k*x`, `x` is written as `l/k`, so that the expansion
        // sees the cancellations between `x` and `log(w)` like the reference
        // `logx` expansion (`-x/w + log(1/w)/w`, or `exp(-x)/w = 1`, with
        // `log(w) = -x`); only needed when the plain expansion fails, vanishes
        // or leaves exponentials of `x` in the coefficient.
        let k = div(logw, &self.x);
        if !plain_ok && k.as_number().is_some_and(|k| !k.is_zero()) && depends(f, &self.x) {
            let mut m = FxHashMap::default();
            m.insert(self.x.clone(), div(&l, &k));
            // `exp(alpha*l + beta) = w**alpha*exp(beta)`
            let g = crate::simptools::bottom_up(&xreplace(f, &m), &mut |u| match u.kind() {
                Kind::Func(Func::Exp, a) => match crate::series::engine::linear_in(&a[0], &l) {
                    Some((alpha, beta)) if !alpha.is_zero() => {
                        mul2(&pow(w, &Expr::num(alpha)), &fx::exp(&beta))
                    }
                    _ => u.clone(),
                },
                _ => u.clone(),
            });
            if let Ok(Some((c, e))) = self.lt_cancel(&g, w, &l, logw) {
                if !c.is_zero() && !depends(&c, w) {
                    return Ok((c, e));
                }
            }
        }
        let (c, e) = match plain {
            Ok(Some(ce)) => ce,
            Ok(None) => return unknown("cancellation in the mrv leading term"),
            Err(SErr::Pole(s)) => return unknown(format!("series in mrv variable: {s}")),
            Err(_) => self.lt_fallback(f, w, logw)?,
        };
        if depends(&c, w) {
            return unknown("leading coefficient depends on w");
        }
        Ok((c, e))
    }

    /// Powers `w**s` of the mrv variable with symbolic exponents differing
    /// by numbers (`w**(-k)`, `w**(-k - 1)` from `(n + 1)**k` and
    /// `(n + 1)**(k + 1)`): the reference series keeps `w**s` as a factor of
    /// the terms. Each is written `P*w**d` with `P = w**s0` a positive
    /// symbol; when every sum is homogeneous in `P`, `P` is a constant of
    /// the expansion and `c*P**m*w**e` gives `(c, e + m*s0)`.
    fn symbolic_power_leadterm(&mut self, f: &Expr, w: &Expr, logw: &Expr) -> Option<LRes<(Expr, Expr)>> {
        let mut exps: Vec<Expr> = vec![];
        collect_w_powers(f, w, &mut exps);
        let s0 = exps.first()?.clone();
        for s in &exps {
            crate::expand::expand(&sub(s, &s0)).as_number()?;
        }
        // w approaches zero through positive reals, but w**s0 need not be
        // positive when s0 is nonreal: its phase can wind indefinitely. Test
        // with a neutral scale BEFORE positivity can simplify nested roots.
        // Rational dependence on that scale commutes with complex scaling;
        // fractional powers do not in general preserve principal branches.
        if crate::assumptions::exact_queries(|| crate::assumptions::is_real(&s0)) != Some(true) {
            let neutral = Expr::dummy("Pneutral", crate::symbol::Assumptions::NONE);
            let trial = subst_w_powers(f, w, &s0, &neutral)?;
            if !rational_in_scale(&trial, &neutral) { return None; }
        }
        let pp = Expr::dummy("P", crate::symbol::Assumptions::POSITIVE);
        let g = subst_w_powers(f, w, &s0, &pp)?;
        let mut left = vec![];
        collect_w_powers(&g, w, &mut left);
        if !left.is_empty() || p_degree(&g, &pp).is_none() {
            return None;
        }
        // Keep every previously successful calculation, including its exact
        // grouping. Failure to eliminate P from the old coefficient is also a
        // reason to try homogeneity, not an early return before that fallback.
        let previous = match self.leadterm_in(&g, w, logw) {
            Ok((c, e)) => p_degree(&c, &pp).and_then(|m| {
                let c1 = div(&c, &pow(&pp, &Expr::num(m.clone())));
                (!depends(&c1, &pp)).then(|| Ok((c1, add(&[e, mul2(&Expr::num(m), &s0)]))))
            }),
            Err(error) => Some(Err(error)),
        };
        if matches!(&previous, Some(Ok(_))) {
            return previous;
        }
        // P is fresh and positive. Exact homogeneity proves
        // g(P,w) = P**m*g(1,w), without distributing any fractional power
        // across a sign or a complex branch. Only rational degrees enter this
        // new route; all other cases retain their original calculation.
        let Some(m) = p_degree(&g, &pp).filter(|m| m.is_rational()) else { return previous; };
        let mut map = FxHashMap::default();
        map.insert(pp.clone(), Expr::one());
        let plain = xreplace(&g, &map);
        // A finite final coefficient alone is insufficient: cancellation of
        // (k+1)/(k+1) could hide a denominator that vanishes when k=-1.
        // Certify the bases of reciprocal/nonintegral powers before allowing
        // the series arithmetic to cancel such parameter factors.
        if !self.normalized_power_domain(&plain, w, logw) {
            return previous;
        }
        match self.leadterm_in(&plain, w, logw) {
            Ok((c, e)) if !depends(&c, w) && !depends(&c, &pp)
                && crate::assumptions::is_finite(&c) == Some(true)
                && crate::assumptions::is_zero(&c) == Some(false) =>
                Some(Ok((c, add(&[e, mul2(&Expr::num(m), &s0)])))),
            _ => previous,
        }
    }

    /// Sufficient domain certificate for the homogeneity fallback only. The
    /// established series/fallback paths are unchanged when it is undecidable.
    fn normalized_power_domain(&mut self, e: &Expr, w: &Expr, logw: &Expr) -> bool {
        for child in e.args() {
            if !self.normalized_power_domain(&child, w, logw) { return false; }
        }
        let mut certify = |base: &Expr, exponent: &Expr| {
            let integer = exponent.as_number().is_some_and(|p| p.is_integer());
            if exponent.as_number().is_some_and(|p| p.is_integer() && !p.is_negative()) {
                return true;
            }
            let coefficient = if depends(base, w) {
                let l = self.fresh("logw");
                match self.lt_cancel(base, w, &l, logw) {
                    Ok(Some((c, _))) if !depends(&c, w) => c,
                    _ => return false,
                }
            } else { base.clone() };
            if crate::assumptions::is_finite(&coefficient) != Some(true)
                || crate::assumptions::is_zero(&coefficient) != Some(false) {
                return false;
            }
            if integer || crate::assumptions::is_positive(&coefficient) == Some(true)
                || crate::assumptions::is_real(&coefficient) == Some(false) {
                return true;
            }
            // On the negative real branch cut the principal coefficient uses
            // argument +pi. A real base or a certified upper-half-plane base
            // has that continuation. A lower-half-plane approach must retain
            // its old calculation, rather than silently using argument +pi.
            let mut positive = FxHashMap::default();
            positive.insert(w.clone(), Expr::dummy("wpos", crate::symbol::Assumptions::POSITIVE));
            let near = xreplace(base, &positive);
            crate::assumptions::is_real(&near) == Some(true)
                || crate::assumptions::is_nonnegative(&fx::im(&near)) == Some(true)
        };
        match e.kind() {
            Kind::Pow(base, exponent) => certify(base, exponent),
            Kind::Mul(_, factors) => factors.iter().all(|(base, exponent)| certify(base, exponent)),
            _ => true,
        }
    }

    /// `b**p` with a symbolic `p` and a base `b ~ c*w**m` that is real and
    /// positive near `w = 0+` is `c**p*w**(m*p)*(b/(c*w**m))**p` (all three
    /// factors positive), which exposes the power of the mrv variable the
    /// reference `exp(p*log(b))` expansion with `log(w) = logw` produces
    /// (`(1/w + 1)**(k + 1)`, issue 22334). `(1/w)**p` becomes `w**(-p)`.
    fn split_symbolic_bases(&mut self, f: &Expr, w: &Expr, logw: &Expr) -> Option<Expr> {
        let wp = Expr::dummy("wpos", crate::symbol::Assumptions::POSITIVE);
        let l = self.fresh("logw");
        let mut changed = false;
        let g = crate::simptools::bottom_up(f, &mut |u| {
            let Kind::Pow(b, p) = u.kind() else { return u.clone() };
            if p.is_number() || depends(p, w) || !depends(b, w) {
                return u.clone();
            }
            if let Some(s) = w_power_exponent(b, p, w) {
                if b != w {
                    changed = true;
                    return pow(w, &s);
                }
                return u.clone();
            }
            let mut m = FxHashMap::default();
            m.insert(w.clone(), wp.clone());
            if crate::assumptions::is_real(&xreplace(b, &m)) != Some(true) {
                return u.clone();
            }
            let Ok(lt) = leading_term_with(b, w, &Expr::zero(), false, &l, logw) else { return u.clone() };
            let c = lt.coeff.clone();
            if c.is_zero() || depends(&c, w) || depends(&lt.coeff_logt, &l) || depends(&c, &self.x) {
                return u.clone();
            }
            if crate::assumptions::is_positive(&c) != Some(true) {
                return u.clone();
            }
            let wm = pow(w, &Expr::num(lt.exp.clone()));
            changed = true;
            mul(&[pow(&c, p), pow(w, &mul2(&Expr::num(lt.exp.clone()), p)), pow(&div(b, &mul2(&c, &wm)), p)])
        });
        changed.then_some(g)
    }

    /// Leading term of `f` in `w`, the expansion keeping `log(w)` as the
    /// symbol `l`: a leading coefficient that vanishes only once `l = logw`
    /// is substituted is subtracted and the next term is taken.
    fn lt_cancel(&mut self, f: &Expr, w: &Expr, l: &Expr, logw: &Expr) -> Result<Option<(Expr, Expr)>, SErr> {
        let mut g = f.clone();
        for _ in 0..8 {
            let lt = leading_term_with(&g, w, &Expr::zero(), false, l, logw)?;
            if exp_of_logw(&lt.coeff_logt, l) {
                return Err(SErr::Unsupported("symbolic power of the mrv variable".into()));
            }
            if lt.coeff.is_zero() && !lt.coeff_logt.is_zero() {
                g = sub(&g, &mul2(&lt.coeff_logt, &pow(w, &Expr::num(lt.exp))));
                continue;
            }
            return Ok(Some((lt.coeff, Expr::num(lt.exp))));
        }
        Ok(None)
    }

    /// Leading term algebra `(c, e)` with symbolic exponents (`w**log(2)`),
    /// used when the series engine cannot represent the exponents.
    fn lt_fallback(&mut self, f: &Expr, w: &Expr, logw: &Expr) -> LRes<(Expr, Expr)> {
        if !depends(f, w) {
            return Ok((f.clone(), Expr::zero()));
        }
        if f == w {
            return Ok((Expr::one(), Expr::one()));
        }
        let l = self.fresh("logw");
        if let Ok(lt) = leading_term_with(f, w, &Expr::zero(), false, &l, logw) {
            if !exp_of_logw(&lt.coeff_logt, &l) {
                return Ok((lt.coeff, Expr::num(lt.exp)));
            }
        }
        match f.kind() {
            Kind::Mul(c, facs) => {
                let mut cc = vec![Expr::num(c.clone())];
                let mut ee = vec![];
                for (b, p) in facs {
                    let (c1, e1) = self.lt_fallback(&mul_factor(b, p), w, logw)?;
                    cc.push(c1);
                    ee.push(e1);
                }
                Ok((mul(&cc), add(&ee)))
            }
            Kind::Pow(b, p) if !depends(p, w) => {
                let (c1, e1) = self.lt_fallback(b, w, logw)?;
                Ok((pow(&c1, p), mul2(&e1, p)))
            }
            Kind::Add(c, terms) => {
                let mut best: Vec<(Expr, Expr)> = vec![];
                if !c.is_zero() {
                    best.push((Expr::num(c.clone()), Expr::zero()));
                }
                for (t, k) in terms {
                    let (c1, e1) = self.lt_fallback(&crate::ops::add_term(t, k), w, logw)?;
                    match best.first() {
                        None => best.push((c1, e1)),
                        Some((_, eb)) => {
                            let d = sub(&e1, eb);
                            match self.sign(&d)? {
                                0 => best.push((c1, e1)),
                                -1 => best = vec![(c1, e1)],
                                _ => {}
                            }
                        }
                    }
                }
                let c0 = add(&best.iter().map(|b| b.0.clone()).collect::<Vec<_>>());
                if (best.len() > 1 && crate::series::engine::hidden_zero(&c0) != Some(false)) || c0.is_zero() {
                    return self.cancelled_sum(f, w, logw);
                }
                Ok((c0, best[0].1.clone()))
            }
            Kind::Func(Func::Log, a) => {
                let (c1, e1) = self.lt_fallback(&a[0], w, logw)?;
                if c1 == Expr::one() && e1.is_zero() {
                    return self.log_near_one(&a[0], w, logw);
                }
                let r = add(&[fx::log(&c1), mul2(&e1, logw)]);
                if r.is_zero() {
                    return unknown("cancellation in log leading term");
                }
                Ok((r, Expr::zero()))
            }
            Kind::Func(Func::Special(crate::special::SpecialFn::Acoth), a) if a.len() == 1 => {
                // The reference acoth._eval_as_leading_term: 1/arg at an
                // infinite argument, the log form at the branch point 0
                // (`-sign(c)*I*pi/2` for a real leading coefficient c); a
                // value on the cut [-1, 1] needs the approach direction.
                let (c1, e1) = self.lt_fallback(&a[0], w, logw)?;
                match self.sign(&e1)? {
                    -1 => Ok((pow(&c1, &Expr::minus_one()), neg(&e1))),
                    1 => match self.sign(&c1) {
                        Ok(s @ (1 | -1)) if crate::assumptions::is_real(&c1) == Some(true) || c1.is_number() => {
                            Ok((mul(&[Expr::rational(-(s as i64), 2), Expr::i(), Expr::pi()]), Expr::zero()))
                        }
                        _ => unknown("acoth at its branch point"),
                    },
                    _ if c1 == Expr::one() || c1 == Expr::minus_one() => {
                        // Log singularity at +-1 (the reference log form
                        // `(log(1 + 1/u) - log(1 - 1/u))/2`): with
                        // u = s*(1 + r), acoth(u) = s*(log(2)/2 - log(r)/2) + o(1),
                        // and log(r) = log(c2) + e2*logw + o(1).
                        let s = c1.clone();
                        let r = sub(&mul2(&s, &a[0]), &Expr::one());
                        let (c2, e2) = self.lt_fallback(&r, w, logw)?;
                        if self.sign(&e2)? != 1 {
                            return unknown("acoth at +-1");
                        }
                        let lr = add(&[fx::log(&c2), mul2(&e2, logw)]);
                        let v = mul2(&s, &div(&sub(&fx::log(&Expr::int(2)), &lr), &Expr::int(2)));
                        Ok((v, Expr::zero()))
                    }
                    _ => {
                        let r = crate::special::apply(crate::special::SpecialFn::Acoth, vec![c1.clone()]);
                        let on_cut = crate::evalf::eval_c64(&c1).map_or(true, |z| z.im.abs() < 1e-12 && z.re.abs() <= 1.0);
                        if on_cut || !is_finite_limit(&r) {
                            return unknown("acoth on its branch cut");
                        }
                        Ok((r, Expr::zero()))
                    }
                }
            }
            Kind::Func(f0, a) if a.len() == 1 => {
                let (c1, e1) = self.lt_fallback(&a[0], w, logw)?;
                let s1 = self.sign(&e1)?;
                let v = match s1 {
                    1 => Expr::zero(),
                    0 => c1.clone(),
                    _ => return unknown("pole in symbolic leading term"),
                };
                let r = crate::functions::func(f0.clone(), vec![v]);
                if r.is_zero() {
                    if s1 == 1 {
                        return self.func_at_zero(f0, &c1, &e1);
                    }
                    return unknown("zero of a function in symbolic leading term");
                }
                // A singular value (Ei(0) = -oo, Ci(0) = -oo) is not a
                // leading coefficient: the log singularities of two such
                // terms may cancel (Ei(-x) - Ei(-2*x) -> log(2)), which
                // `-oo + oo = nan` would hide.
                if !is_finite_limit(&r) {
                    return unknown("singular value of a function in symbolic leading term");
                }
                Ok((r, Expr::zero()))
            }
            _ => unknown(format!("leading term of {f}")),
        }
    }
    fn rewrite(&mut self, e: &Expr, mut omega: Vec<Expr>, w: &Expr) -> LRes<(Expr, Expr)> {
        let arg = |f: &Expr| -> LRes<Expr> {
            match f.kind() {
                Kind::Func(Func::Exp, a) => Ok(a[0].clone()),
                _ => unknown(format!("mrv element {f} is not an exponential")),
            }
        };
        omega.sort_by_key(|f| Reverse(node_count(f)));
        let g = omega.last().unwrap().clone();
        omega.reverse(); // smallest first: replacements of nested elements
        let garg = arg(&g)?;
        let sig = self.sign(&garg)?;
        let wsym = match sig {
            1 => pow(w, &Expr::minus_one()),
            -1 => w.clone(),
            _ => return unknown("bad mrv element"),
        };
        // All elements are substituted simultaneously (the canonical product
        // of exponentials would otherwise recombine replaced and pending
        // elements); nested elements are rewritten inside the replacements.
        let mut m: FxHashMap<Expr, Expr> = FxHashMap::default();
        for h in &omega {
            let harg = arg(h)?;
            let c = self.limitinf(&div(&harg, &garg))?;
            if depends(&c, &self.x) || is_inf(&c) || c.is_constant(Constant::NaN) || c.is_constant(Constant::ComplexInfinity) {
                return unknown("incomparable mrv elements");
            }
            let a = xreplace(&sub(&harg, &mul2(&c, &garg)), &m);
            let repl = mul2(&fx::exp(&a), &pow(&wsym, &c));
            m.insert(h.clone(), repl);
        }
        let f = xreplace(e, &m);
        let logw = if sig == -1 { garg } else { neg(&garg) };
        Ok((f, logw))
    }
}

/// Does the leading coefficient hold `exp(s*log(w))`, a power `w**s` with a
/// symbolic exponent? Such a factor is in the mrv class itself (`x**b/x`
/// gives `w**(1 - b)`), so the numeric exponent of the term is not its order
/// and the symbolic exponent algebra (which needs the sign of `1 - b`) must
/// decide instead.
fn exp_of_logw(c: &Expr, l: &Expr) -> bool {
    match c.kind() {
        Kind::Func(Func::Exp, a) if depends(&a[0], l) => true,
        _ => depends(c, l) && c.args().iter().any(|a| exp_of_logw(a, l)),
    }
}

/// `exp(a)*exp(b)*r` as `exp(a + b)*r` (integer powers only), or `None`
/// when the product has fewer than two exponential factors.
fn merge_exp_factors(e: &Expr) -> Option<Expr> {
    let Kind::Mul(c, facs) = e.kind() else { return None };
    let mut args = vec![];
    let mut rest = vec![Expr::num(c.clone())];
    for (b, p) in facs {
        match b.kind() {
            Kind::Func(Func::Exp, a) if p.as_number().is_some_and(|k| k.is_integer()) => {
                args.push(mul2(&a[0], p))
            }
            _ => rest.push(mul_factor(b, p)),
        }
    }
    if args.len() < 2 {
        return None;
    }
    rest.push(fx::exp(&add(&args)));
    let m = mul(&rest);
    (m != *e).then_some(m)
}

/// Merges the exponential factors (integer powers, exponents depending on
/// `x`) of every product of `e` into one exponential when the expanded sum of
/// their exponents is smaller than the sum itself (terms cancel).
fn merge_cancelling_exps(e: &Expr, x: &Expr) -> Expr {
    crate::simptools::bottom_up(e, &mut |u| {
        let Kind::Mul(c, facs) = u.kind() else { return u.clone() };
        let mut args = vec![];
        let mut rest = vec![Expr::num(c.clone())];
        for (b, p) in facs {
            match b.kind() {
                Kind::Func(Func::Exp, a) if depends(&a[0], x) && p.as_number().is_some_and(|k| k.is_integer()) => {
                    args.push(mul2(&a[0], p))
                }
                _ => rest.push(mul_factor(b, p)),
            }
        }
        if args.len() < 2 {
            return u.clone();
        }
        let s = add(&args);
        let ex = crate::simp::expand_mul(&s);
        if node_count(&ex) >= node_count(&s) {
            return u.clone();
        }
        rest.push(fx::exp(&ex));
        mul(&rest)
    })
}

/// [`merge_exp_factors`] applied to every product of `e`.
fn merge_exp_deep(e: &Expr) -> Expr {
    crate::simptools::bottom_up(e, &mut |u| merge_exp_factors(u).unwrap_or_else(|| u.clone()))
}

/// A limit value that is a finite number or a finite parameter expression.
fn is_finite_limit(l: &Expr) -> bool {
    fn bad(e: &Expr) -> bool {
        match e.kind() {
            Kind::Const(c) => c.is_infinite() || *c == Constant::NaN,
            Kind::Func(Func::Limit | Func::Order, _) => true,
            _ => e.args().iter().any(bad),
        }
    }
    !bad(l)
}

/// Functions whose growth or decay at an infinite argument leaves the
/// comparability class of the argument (factorial-like, exponential-like or
/// exponentially decaying), so that the mrv set of their arguments does not
/// bound them.
fn grows_beyond_its_argument(f: &Func) -> bool {
    use crate::special::SpecialFn as S;
    match f {
        Func::Gamma => true,
        Func::Special(s) => matches!(
            s,
            S::Factorial
                | S::Factorial2
                | S::Subfactorial
                | S::Binomial
                | S::RisingFactorial
                | S::FallingFactorial
                | S::Catalan
                | S::Fibonacci
                | S::Lucas
                | S::Bell
                | S::Partition
                | S::Tribonacci
                | S::Bernoulli
                | S::Euler
                | S::Genocchi
                | S::Andre
                | S::Motzkin
                | S::Beta
                | S::LowerGamma
                | S::UpperGamma
                | S::Erfi
                | S::Ei
                | S::Shi
                | S::Chi
                | S::Expint
                | S::E1
                | S::BesselI
                | S::BesselK
                | S::Hankel1
                | S::Hankel2
                | S::Hn1
                | S::Hn2
                | S::AiryAi
                | S::AiryBi
                | S::AiryAiPrime
                | S::AiryBiPrime
                | S::Hyper
                | S::MeijerG
                | S::MarcumQ
                | S::AppellF1
                | S::Lerchphi
        ),
        _ => false,
    }
}

/// `log(exp(a)) = a`, `log(exp(a)**k) = k*a` (real variables).
fn simp_log(e: &Expr) -> Expr {
    if let Kind::Func(Func::Log, a) = e.kind() {
        match a[0].kind() {
            Kind::Func(Func::Exp, b) => return b[0].clone(),
            Kind::Pow(b, k) if b.is_func(&Func::Exp) => {
                if let Kind::Func(_, bb) = b.kind() {
                    return mul2(k, &bb[0]);
                }
            }
            _ => {}
        }
    }
    e.clone()
}

fn simp_logs(e: &Expr) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    let e1 = crate::traverse::map_children(e, &mut |c| simp_logs(c));
    if e1.is_func(&Func::Log) {
        simp_log(&e1)
    } else {
        e1
    }
}

/// Whether `e` matches SymPy's `I*Wild("a", exclude=[I])`.
fn imaginary_multiple(e: &Expr) -> bool {
    let i = Expr::i();
    if *e == i {
        return true;
    }
    if !e.is_mul() {
        return false;
    }
    let factors = e.args();
    let k = factors.iter().filter(|f| **f == i).count();
    k == 1 && factors.iter().filter(|f| **f != i).all(|f| !crate::traverse::has(f, &i))
}

/// The exponent of `w` in `b**p` (`w**p`, `(1/w)**p`) for a symbolic `p`.
fn w_power_exponent(b: &Expr, p: &Expr, w: &Expr) -> Option<Expr> {
    if p.is_number() || depends(p, w) {
        return None;
    }
    if b == w {
        return Some(p.clone());
    }
    match b.kind() {
        Kind::Pow(bb, q) if bb == w && q.is_number() => Some(mul2(q, p)),
        Kind::Mul(c, f) if f.len() == 1 && c.is_one() && f[0].0 == *w && f[0].1.is_number() => Some(mul2(&f[0].1, p)),
        _ => None,
    }
}

fn collect_w_powers(e: &Expr, w: &Expr, out: &mut Vec<Expr>) {
    if !depends(e, w) {
        return;
    }
    match e.kind() {
        Kind::Pow(b, p) => {
            if let Some(s) = w_power_exponent(b, p, w) {
                out.push(s);
                return;
            }
            collect_w_powers(b, w, out);
            collect_w_powers(p, w, out);
        }
        Kind::Mul(_, facs) => {
            for (b, p) in facs {
                if let Some(s) = w_power_exponent(b, p, w) {
                    out.push(s);
                } else {
                    collect_w_powers(b, w, out);
                    collect_w_powers(p, w, out);
                }
            }
        }
        _ => {
            for a in e.args() {
                collect_w_powers(&a, w, out);
            }
        }
    }
}

fn subst_w_power(b: &Expr, p: &Expr, w: &Expr, s0: &Expr, pp: &Expr) -> Option<Option<Expr>> {
    match w_power_exponent(b, p, w) {
        Some(s) => {
            let d = crate::expand::expand(&sub(&s, s0));
            let d = d.as_number()?.clone();
            Some(Some(mul2(pp, &pow(w, &Expr::num(d)))))
        }
        None => Some(None),
    }
}

fn subst_w_powers(e: &Expr, w: &Expr, s0: &Expr, pp: &Expr) -> Option<Expr> {
    if !depends(e, w) {
        return Some(e.clone());
    }
    match e.kind() {
        Kind::Pow(b, p) => match subst_w_power(b, p, w, s0, pp)? {
            Some(r) => Some(r),
            None => Some(pow(&subst_w_powers(b, w, s0, pp)?, &subst_w_powers(p, w, s0, pp)?)),
        },
        Kind::Mul(c, facs) => {
            let mut v = vec![Expr::num(c.clone())];
            for (b, p) in facs {
                match subst_w_power(b, p, w, s0, pp)? {
                    Some(r) => v.push(r),
                    None => v.push(pow(&subst_w_powers(b, w, s0, pp)?, &subst_w_powers(p, w, s0, pp)?)),
                }
            }
            Some(mul(&v))
        }
        _ => {
            let mut ok = true;
            let r = crate::traverse::map_children(e, &mut |c: &Expr| match subst_w_powers(c, w, s0, pp) {
                Some(v) => v,
                None => {
                    ok = false;
                    c.clone()
                }
            });
            ok.then_some(r)
        }
    }
}

/// The scale-dependent operations must be rational, over arbitrary
/// scale-independent coefficients. Together with p_degree homogeneity this
/// proves g(P,w)=P**m*g(1,w) for every nonzero complex P, not only positive P.
/// In particular an integral *total degree* does not make sqrt(P**2) safe.
fn rational_in_scale(e: &Expr, scale: &Expr) -> bool {
    if !depends(e, scale) || e == scale { return true; }
    let power = |base: &Expr, exponent: &Expr| {
        !depends(exponent, scale) && (!depends(base, scale)
            || (exponent.as_number().is_some_and(|n| n.is_integer())
                && rational_in_scale(base, scale)))
    };
    match e.kind() {
        Kind::Pow(base, exponent) => power(base, exponent),
        Kind::Mul(_, factors) => factors.iter().all(|(base, exponent)| power(base, exponent)),
        Kind::Add(_, terms) => terms.iter().all(|(term, _)| rational_in_scale(term, scale)),
        _ => false,
    }
}

/// Degree of `e` in `P` when every sum is homogeneous in `P` and `P`
/// occurs only in products and numeric powers.
fn p_degree(e: &Expr, pp: &Expr) -> Option<crate::number::Number> {
    use crate::number::Number;
    if !depends(e, pp) {
        return Some(Number::ZERO);
    }
    if e == pp {
        return Some(Number::ONE);
    }
    let factor = |b: &Expr, p: &Expr| -> Option<Number> {
        if !depends(b, pp) && !depends(p, pp) {
            return Some(Number::ZERO);
        }
        if depends(p, pp) {
            return None;
        }
        let k = p.as_number()?;
        Some(p_degree(b, pp)?.mul(k))
    };
    match e.kind() {
        Kind::Pow(b, p) => factor(b, p),
        Kind::Mul(_, facs) => {
            let mut d = Number::ZERO;
            for (b, p) in facs {
                d = d.add(&factor(b, p)?);
            }
            Some(d)
        }
        Kind::Add(c, terms) => {
            let mut d: Option<Number> = if c.is_zero() { None } else { Some(Number::ZERO) };
            for (t, _) in terms {
                let dt = p_degree(t, pp)?;
                match &d {
                    None => d = Some(dt),
                    Some(d0) if *d0 == dt => {}
                    _ => return None,
                }
            }
            d
        }
        _ => None,
    }
}

#[cfg(test)]
mod homogeneous_domain_tests {
    use super::*;
    use crate::symbol::Assumptions as A;

    #[test]
    fn normalized_reciprocal_requires_a_nonzero_leading_coefficient() {
        let w=Expr::symbol("w"); let x=Expr::symbol_with("x",A::POSITIVE);
        for (facts,expected) in [(A::POSITIVE,true),(A::NONNEGATIVE,true),(A::REAL,false),(A::NEGATIVE,false),(A::NONE,false)] {
            let k=Expr::symbol_with("k",facts);
            let denominator=sub(&pow(&add(&[Expr::one(),w.clone()]),&add(&[k,Expr::one()])),&Expr::one());
            let expression=pow(&denominator,&Expr::minus_one());
            assert_eq!(Gruntz::new(&x).normalized_power_domain(&expression,&w,&neg(&x)),expected,"{expression}");
        }
    }

    #[test]
    fn normalized_fractional_powers_preserve_branch_sides() {
        let w=Expr::symbol("w"); let x=Expr::symbol_with("x",A::POSITIVE);
        let real_negative=pow(&sub(&Expr::minus_one(),&w),&Expr::rational(1,2));
        assert!(Gruntz::new(&x).normalized_power_domain(&real_negative,&w,&neg(&x)));
        for (correction,allowed) in [(mul2(&Expr::i(),&w),true),(neg(&mul2(&Expr::i(),&w)),false)] {
            let branch=pow(&add(&[Expr::minus_one(),correction]),&Expr::rational(1,2));
            assert_eq!(Gruntz::new(&x).normalized_power_domain(&branch,&w,&neg(&x)),allowed);
        }
    }
}

#[cfg(test)]
mod winding_scale_tests {
    use super::*;
    #[test]
    fn complex_scale_certificate_requires_integer_dependent_powers() {
        let p=Expr::symbol("P");
        for source in ["P","P**3","1/P","(P*w+P)/(P*w-P)","sqrt(w)*P"] {
            let e=crate::parser::parse_str(source).unwrap();
            assert!(rational_in_scale(&e,&p),"{e}");
        }
        for source in ["sqrt(P)","sqrt(P*w)","sqrt(-P)","sqrt(P**2)","sqrt(P)/sqrt(-P)","P**P","P**w","exp(P)"] {
            let e=crate::parser::parse_str(source).unwrap();
            assert!(!rational_in_scale(&e,&p),"{e}");
        }
    }
    #[test]
    fn positive_scale_simplification_is_not_used_as_a_winding_certificate() {
        let w=Expr::symbol("w");let x=Expr::symbol("x");
        let e=crate::parser::parse_str("sqrt(w**I*((1+w)**(1+I)-1))/sqrt(-w**I*((1+w)**(1+I)-1))").unwrap();
        assert!(Gruntz::new(&x).symbolic_power_leadterm(&e,&w,&neg(&x)).is_none());
    }
    #[test]
    fn proven_real_scale_retains_nonintegral_homogeneous_powers() {
        let w=Expr::symbol("w");let x=Expr::symbol("x");
        let k=Expr::symbol_with("k",crate::symbol::Assumptions::REAL);
        let base=pow(&w,&k);
        let e=div(&pow(&base,&Expr::rational(1,2)),&pow(&neg(&base),&Expr::rational(1,2)));
        let actual=Gruntz::new(&x).symbolic_power_leadterm(&e,&w,&neg(&x)).unwrap().unwrap();
        assert_eq!(actual,(neg(&Expr::i()),Expr::zero()));
    }

}
