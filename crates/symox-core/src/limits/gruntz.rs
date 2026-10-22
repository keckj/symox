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
        match self.sign(c0) {
            Ok(1) => Ok(Expr::oo()),
            Ok(-1) => Ok(Expr::neg_oo()),
            Ok(_) => unknown("zero leading coefficient"),
            Err(err) => {
                if !depends(c0, &self.x)
                    && !crate::functions::is_numeric_value(c0)
                    && crate::evalf::eval_c64(c0).is_none()
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
            return match sign_of(&C::from_expr(e.clone())) {
                Some(s) if s != 0 => Ok(s),
                Some(_) => Ok(0),
                // the reference `sign(e, x)` asks the assumptions first
                None if crate::assumptions::is_positive(e) == Some(true) => Ok(1),
                None if crate::assumptions::is_negative(e) == Some(true) => Ok(-1),
                None if crate::assumptions::is_zero(e) == Some(true) => Ok(0),
                None => unknown(format!("sign of {e}")),
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
        let pe = self.prepare(e)?;
        let (c0, _) = self.mrv_leadterm(&pe)?;
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

    fn prepare_node(&mut self, e: &Expr) -> LRes<Expr> {
        let x = self.x.clone();
        let two = Expr::int(2);
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
                if is_inf(&l) {
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
        let c = self.limitinf(&div(&la, &lb))?;
        if c.is_zero() {
            Ok(Ordering::Less)
        } else if is_inf(&c) {
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
        let l = self.fresh("logw");
        // With `log(w) = k*x`, `x` is written as `l/k` first, so that the
        // expansion sees the cancellations between `x` and `log(w)` like the
        // reference `logx` expansion (`-x/w + log(1/w)/w`, or
        // `exp(-x)/w = 1`, with `log(w) = -x`).
        let k = div(logw, &self.x);
        if k.as_number().is_some_and(|k| !k.is_zero()) && depends(f, &self.x) {
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
        let (c, e) = match self.lt_cancel(f, w, &l, logw) {
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

    /// Leading term of `f` in `w`, the expansion keeping `log(w)` as the
    /// symbol `l`: a leading coefficient that vanishes only once `l = logw`
    /// is substituted is subtracted and the next term is taken.
    fn lt_cancel(&mut self, f: &Expr, w: &Expr, l: &Expr, logw: &Expr) -> Result<Option<(Expr, Expr)>, SErr> {
        let mut g = f.clone();
        for _ in 0..8 {
            let lt = leading_term_with(&g, w, &Expr::zero(), false, l, logw)?;
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
            return Ok((lt.coeff, Expr::num(lt.exp)));
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
                if c0.is_zero() {
                    return unknown("cancellation in symbolic leading term");
                }
                Ok((c0, best[0].1.clone()))
            }
            Kind::Func(Func::Log, a) => {
                let (c1, e1) = self.lt_fallback(&a[0], w, logw)?;
                let r = add(&[fx::log(&c1), mul2(&e1, logw)]);
                if r.is_zero() {
                    return unknown("cancellation in log leading term");
                }
                Ok((r, Expr::zero()))
            }
            Kind::Func(f0, a) if a.len() == 1 => {
                let (c1, e1) = self.lt_fallback(&a[0], w, logw)?;
                let v = match self.sign(&e1)? {
                    1 => Expr::zero(),
                    0 => c1,
                    _ => return unknown("pole in symbolic leading term"),
                };
                let r = crate::functions::func(f0.clone(), vec![v]);
                if r.is_zero() {
                    return unknown("zero of a function in symbolic leading term");
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
