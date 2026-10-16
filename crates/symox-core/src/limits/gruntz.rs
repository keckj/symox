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
use crate::number::Number;
use crate::ops::{add, div, mul, mul2, mul_factor, neg, pow, sub};
use crate::series::coeff::C;
use crate::series::engine::{depends, sign_of};
use crate::series::leading_term_with;
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
        if std::env::var("SYMOX_TRACE").is_ok() { eprintln!("{}limitinf {}", " ".repeat(self.depth as usize), e); }
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
        if e0.is_positive() {
            return Ok(Expr::zero());
        }
        if e0.is_negative() {
            return self.infinity_with_sign_of(&c0);
        }
        self.limitinf(&c0)
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
            Kind::Pow(b, p) if depends(p, &x) || (depends(b, &x) && !p.is_number()) => {
                fx::exp(&mul2(p, &simp_log(&fx::log(b))))
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
            Kind::Func(_, args) => {
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
    pub fn mrv_leadterm(&mut self, e: &Expr) -> LRes<(Expr, Number)> {
        let x = self.x.clone();
        if !depends(e, &x) {
            return Ok((e.clone(), Number::ZERO));
        }
        let mut omega = self.mrv(e)?;
        if omega.is_empty() {
            return Ok((e.clone(), Number::ZERO));
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
            let om_up: Vec<Expr> = omega.iter().map(|f| self.move_up(f)).collect();
            if om_up.iter().all(|f| crate::traverse::has(&up, f)) {
                omega = om_up;
            } else {
                omega = self.mrv(&up)?;
                if omega.contains(&x) {
                    return unknown("cannot move up");
                }
            }
            e = up;
            if omega.is_empty() {
                return Ok((e, Number::ZERO));
            }
        }
        let w = self.fresh("w");
        let (f, logw) = self.rewrite(&e, omega, &w)?;
        self.leadterm_in(&f, &w, &logw)
    }

    fn leadterm_in(&mut self, f: &Expr, w: &Expr, logw: &Expr) -> LRes<(Expr, Number)> {
        let l = self.fresh("logw");
        let lt = leading_term_with(f, w, &Expr::zero(), false, &l, logw)
            .map_err(|er| LimitError::Unknown(format!("series in mrv variable: {er}")))?;
        if depends(&lt.coeff, w) {
            return unknown("leading coefficient depends on w");
        }
        Ok((lt.coeff, lt.exp))
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
