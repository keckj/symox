//! `Basic.rewrite`: SymPy's rewriting engine together with the
//! `_eval_rewrite_as_*` rules of the elementary functions (trigonometric,
//! hyperbolic and their inverses, `exp`, `Pow`, `sinc`, `re`/`im`, `sign`,
//! `Abs`, `arg`, `Max`/`Min`, `floor`/`ceiling`/`frac`, `Mod`, `Heaviside`,
//! `DiracDelta`, `Piecewise` and `E`).
//!
//! Ported from SymPy (BSD-3-Clause, see LICENSE-SYMPY): `Basic.rewrite` /
//! `Basic._rewrite` and the rewrite methods of
//! `sympy/functions/elementary/*.py`, `delta_functions.py`,
//! `core/power.py`, `core/mod.py` and `core/numbers.py`.
//!
//! Semantics follow SymPy: with `deep` the arguments are rewritten first
//! (bottom-up), a node is rewritten only when no pattern is given or it is an
//! instance of one of the patterns, the rule method receives the rewritten
//! arguments, and when it declines the node is rebuilt from the rewritten
//! arguments with the evaluating constructor.

use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};

use super::trig_eval::{chebyshev_t, pi_coeff};
use super::SpecialFn as S;
use crate::assumptions as asm;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{self as fx, func};
use crate::logic::{self, LogicOp};
use crate::number::Number;
use crate::ops::{add, add2, div, inv, mul, mul2, neg, pow, sqrt, sub};
use crate::relational::{eq, ge, gt, le, lt, ne, piecewise, rel};

/// `e.rewrite(*patterns, rule, deep=deep)`; `rule` and `patterns` are SymPy
/// class or function names (`"cos"`, `"Piecewise"`, `"sincos"`, ...).
pub fn rewrite(e: &Expr, patterns: &[String], rule: &str, deep: bool) -> Expr {
    let rule = canonical_rule(rule).to_string();
    let pats: Vec<String> =
        patterns.iter().map(|p| canonical_rule(p).to_string()).filter(|p| has_instance(e, p)).collect();
    if !patterns.is_empty() && pats.is_empty() {
        return e.clone();
    }
    Rewriter { pats, rule, deep }.go(e)
}

/// `e.rewrite(rule)` with the default options.
fn rw(e: &Expr, rule: &str) -> Expr {
    rewrite(e, &[], rule, true)
}

fn canonical_rule(r: &str) -> &str {
    match r {
        "abs" => "Abs",
        "ln" => "log",
        _ => r,
    }
}

struct Rewriter {
    pats: Vec<String>,
    rule: String,
    deep: bool,
}

impl Rewriter {
    fn go(&self, e: &Expr) -> Expr {
        if opaque(e) {
            return if self.deep { crate::traverse::map_children(e, &mut |c| self.go(c)) } else { e.clone() };
        }
        let args = e.args();
        let new_args: Vec<Expr> = if self.deep { args.iter().map(|a| self.go(a)).collect() } else { args.clone() };
        if self.pats.is_empty() || self.pats.iter().any(|p| is_instance(e, p)) {
            if let Some(r) = method(e, &new_args, &self.rule) {
                return r;
            }
        }
        if args.is_empty() {
            return e.clone();
        }
        rebuild(e, &args, new_args)
    }
}

/// Nodes whose public arguments are not plain sub-expressions.
fn opaque(e: &Expr) -> bool {
    match e.kind() {
        Kind::Derivative(..) => true,
        Kind::Func(f, _) => matches!(
            f,
            Func::MatrixElement(_)
                | Func::MatrixProxy(_)
                | Func::MatrixTrace(_)
                | Func::MatrixDeterminant(_)
                | Func::MatrixDotProduct(..)
                | Func::Order
                | Func::Subs
                | Func::Limit
                | Func::Integral
                | Func::Sum
                | Func::Predicate(_)
                | Func::AppliedPredicate
        ),
        _ => false,
    }
}

/// `self.func(*args)`.
fn rebuild(e: &Expr, old: &[Expr], new: Vec<Expr>) -> Expr {
    if let Kind::Func(f, _) = e.kind() {
        return func(f.clone(), new);
    }
    if new.iter().zip(old).all(|(a, b)| a.same_presentation(b)) {
        return e.clone();
    }
    match e.kind() {
        Kind::Add(..) => add(&new),
        // `Mul.flatten` rebuilds every factor as `Pow(*f.as_base_exp())`, which
        // evaluates held `exp`/`Pow` factors (`-exp(2*log(x))` -> `-x**2`)
        Kind::Mul(..) => mul(&new.iter().map(reevaluate_factor).collect::<Vec<_>>()),
        Kind::Pow(..) => pow(&new[0], &new[1]),
        Kind::Rel(op, ..) => rel(*op, &new[0], &new[1]),
        Kind::Piecewise(_) => piecewise(new.chunks(2).map(|c| (c[0].clone(), c[1].clone())).collect()),
        _ => e.clone(),
    }
}

fn reevaluate_factor(f: &Expr) -> Expr {
    match f.kind() {
        Kind::Func(Func::Exp, a) if a.len() == 1 => fx::exp(&a[0]),
        Kind::Pow(b, x) => pow(b, x),
        _ => f.clone(),
    }
}

// --------------------------------------------------------------------
// class tests
// --------------------------------------------------------------------

fn const_class(c: Constant) -> &'static str {
    match c {
        Constant::Pi => "Pi",
        Constant::E => "Exp1",
        Constant::I => "ImaginaryUnit",
        Constant::EulerGamma => "EulerGamma",
        Constant::Infinity => "Infinity",
        Constant::NegInfinity => "NegativeInfinity",
        Constant::ComplexInfinity => "ComplexInfinity",
        Constant::NaN => "NaN",
        Constant::Catalan => "Catalan",
        Constant::GoldenRatio => "GoldenRatio",
        Constant::TribonacciConstant => "TribonacciConstant",
    }
}

fn trig(f: &Func) -> bool {
    matches!(f, Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc)
}

fn hyperbolic(f: &Func) -> bool {
    matches!(f, Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth | Func::Special(S::Sech) | Func::Special(S::Csch))
}

/// `isinstance(e, <class named p>)` for the classes SymPy users pass as
/// rewrite patterns.
fn is_instance(e: &Expr, p: &str) -> bool {
    if p == "Basic" {
        return true;
    }
    match e.kind() {
        Kind::Add(..) => matches!(p, "Add" | "Expr" | "AssocOp"),
        Kind::Mul(..) => matches!(p, "Mul" | "Expr" | "AssocOp"),
        Kind::Pow(..) => matches!(p, "Pow" | "Expr"),
        Kind::Piecewise(_) => matches!(p, "Piecewise" | "Function" | "Expr"),
        Kind::Sym(_) => matches!(p, "Symbol" | "Atom" | "Expr" | "AtomicExpr"),
        Kind::Num(_) => matches!(p, "Number" | "Atom" | "Expr" | "AtomicExpr"),
        Kind::Const(c) => {
            p == const_class(*c)
                || matches!(p, "Atom" | "Expr" | "AtomicExpr")
                || (p == "NumberSymbol"
                    && !matches!(
                        c,
                        Constant::I | Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN
                    ))
        }
        Kind::Rel(op, ..) => {
            matches!(p, "Relational" | "Boolean")
                || p == match op {
                    crate::expr::RelOp::Eq => "Equality",
                    crate::expr::RelOp::Ne => "Unequality",
                    crate::expr::RelOp::Lt => "StrictLessThan",
                    crate::expr::RelOp::Le => "LessThan",
                    crate::expr::RelOp::Gt => "StrictGreaterThan",
                    crate::expr::RelOp::Ge => "GreaterThan",
                }
        }
        Kind::Derivative(..) => matches!(p, "Derivative" | "Expr"),
        Kind::Bool(_) => matches!(p, "Boolean" | "BooleanAtom" | "Atom"),
        Kind::Func(f, _) => {
            if f.name() == p {
                return true;
            }
            if let Func::Logic(_) = f {
                return p == "Boolean";
            }
            if matches!(f, Func::Set(_)) {
                return p == "Set";
            }
            match p {
                "Function" | "Application" | "Expr" => true,
                "TrigonometricFunction" => trig(f),
                "ReciprocalTrigonometricFunction" => matches!(f, Func::Sec | Func::Csc),
                "InverseTrigonometricFunction" => matches!(
                    f,
                    Func::Asin | Func::Acos | Func::Atan | Func::Acot | Func::Special(S::Asec) | Func::Special(S::Acsc)
                ),
                "HyperbolicFunction" => hyperbolic(f),
                "ReciprocalHyperbolicFunction" => matches!(f, Func::Special(S::Sech) | Func::Special(S::Csch)),
                "InverseHyperbolicFunction" => matches!(
                    f,
                    Func::Asinh
                        | Func::Atanh
                        | Func::Acosh
                        | Func::Special(S::Acoth)
                        | Func::Special(S::Asech)
                        | Func::Special(S::Acsch)
                ),
                "ExpBase" => matches!(f, Func::Exp),
                "MinMaxBase" => matches!(f, Func::Max | Func::Min),
                "RoundFunction" => matches!(f, Func::Floor | Func::Ceiling),
                _ => false,
            }
        }
    }
}

fn has_instance(e: &Expr, p: &str) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| found = found || is_instance(s, p));
    found
}

fn has_func(e: &Expr, f: &Func) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| found = found || s.is_func(f));
    found
}

fn has_symbol(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| found = found || s.is_symbol());
    found
}

// --------------------------------------------------------------------
// construction helpers
// --------------------------------------------------------------------

fn ii() -> Expr {
    Expr::i()
}

fn pi() -> Expr {
    Expr::pi()
}

fn half() -> Expr {
    Expr::rational(1, 2)
}

fn int(n: i64) -> Expr {
    Expr::int(n)
}

/// Evaluated `f(x)`.
fn f1(f: Func, x: Expr) -> Expr {
    func(f, vec![x])
}

/// `f(x, evaluate=False)`.
fn held(f: Func, x: Expr) -> Expr {
    Expr::from_kind(Kind::Func(f, vec![x]))
}

fn sp(s: S, args: Vec<Expr>) -> Expr {
    func(Func::Special(s), args)
}

fn m2(a: &Expr, b: &Expr) -> Expr {
    mul2(a, b)
}

fn sq(x: &Expr) -> Expr {
    pow(x, &int(2))
}

fn i_times(x: &Expr) -> Expr {
    mul2(&ii(), x)
}

fn half_pi() -> Expr {
    div(&pi(), &int(2))
}

/// `I*pi/2`.
fn half_ipi() -> Expr {
    mul(&[ii(), pi(), half()])
}

fn true_() -> Expr {
    Expr::boolean(true)
}

fn pw(pairs: Vec<(Expr, Expr)>) -> Expr {
    piecewise(pairs)
}

fn ext_real(e: &Expr) -> bool {
    asm::fact(e, "extended_real") == Some(true)
}

fn single_arg<'a>(e: &'a Expr, f: &Func) -> Option<&'a Expr> {
    match e.kind() {
        Kind::Func(g, a) if g == f && a.len() == 1 => Some(&a[0]),
        _ => None,
    }
}

/// `arg.func(arg.args[0]).rewrite(exp)` for a trigonometric or hyperbolic
/// argument (SymPy's `sin/cos/tan/cot._eval_rewrite_as_exp` preamble).
fn exp_ready(a: &Expr) -> Expr {
    if let Kind::Func(f, args) = a.kind() {
        if (trig(f) || hyperbolic(f)) && !args.is_empty() {
            return rw(&f1(f.clone(), args[0].clone()), "exp");
        }
    }
    a.clone()
}

// --------------------------------------------------------------------
// dispatch
// --------------------------------------------------------------------

fn method(e: &Expr, args: &[Expr], rule: &str) -> Option<Expr> {
    match e.kind() {
        Kind::Func(f, orig) => func_method(f, orig, args, rule, e),
        Kind::Pow(..) => pow_method(&args[0], &args[1], rule),
        Kind::Const(Constant::E) => match rule {
            "sin" => Some(sub(&fx::sin(&add2(&ii(), &half_pi())), &i_times(&fx::sin(&ii())))),
            "cos" => Some(add2(&fx::cos(&ii()), &i_times(&fx::cos(&add2(&ii(), &half_pi()))))),
            _ => None,
        },
        Kind::Piecewise(_) if rule == "KroneckerDelta" => {
            let pairs: Vec<(Expr, Expr)> = args.chunks(2).map(|c| (c[0].clone(), c[1].clone())).collect();
            piecewise_as_kronecker(&pairs)
        }
        _ => None,
    }
}

fn func_method(f: &Func, orig: &[Expr], args: &[Expr], rule: &str, this: &Expr) -> Option<Expr> {
    if args.len() == 1 {
        let a = &args[0];
        match f {
            Func::Sin | Func::Cos | Func::Tan | Func::Cot => return trig_method(f, a, rule, this),
            Func::Sec | Func::Csc => return reciprocal_trig_method(f, &orig[0], a, rule),
            Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth => return hyperbolic_method(f, a, rule),
            Func::Special(S::Sech) | Func::Special(S::Csch) => return reciprocal_hyperbolic_method(f, a, rule),
            _ => {}
        }
        if let Some(r) = inverse_method(f, a, rule) {
            return Some(r);
        }
        return unary_method(f, &orig[0], a, rule);
    }
    match (f, rule) {
        (Func::Atan2, _) if args.len() == 2 => atan2_method(&args[0], &args[1], rule),
        (Func::Max | Func::Min, _) => minmax_method(f, args, rule),
        (Func::Heaviside, _) if args.len() == 2 => heaviside_method(&args[0], &args[1], rule),
        (Func::Special(S::Mod), "floor") if args.len() == 2 => {
            let (a, b) = (&args[0], &args[1]);
            Some(sub(a, &m2(b, &fx::floor(&div(a, b)))))
        }
        _ => None,
    }
}

// --------------------------------------------------------------------
// sin, cos, tan, cot
// --------------------------------------------------------------------

/// `_eval_rewrite_as_<rule>` of `sin`, `cos`, `tan` or `cot` applied to the
/// (rewritten) argument `a`; `this` is SymPy's `self`.
fn trig_method(f: &Func, a: &Expr, rule: &str, this: &Expr) -> Option<Expr> {
    let pi_2 = half_pi();
    match (f, rule) {
        (Func::Sin, "exp") => {
            let a = exp_ready(a);
            let ia = i_times(&a);
            Some(div(&sub(&fx::exp(&ia), &fx::exp(&m2(&neg(&a), &ii()))), &m2(&int(2), &ii())))
        }
        (Func::Sin, "Pow") => {
            let x = single_arg(a, &Func::Log)?;
            let (xi, xmi) = (pow(x, &ii()), pow(x, &neg(&ii())));
            Some(sub(&div(&i_times(&xmi), &int(2)), &div(&i_times(&xi), &int(2))))
        }
        (Func::Sin, "cos") => Some(held(Func::Cos, sub(a, &pi_2))),
        (Func::Sin, "tan") => {
            let t = fx::tan(&m2(&half(), a));
            Some(div(&m2(&int(2), &t), &add2(&Expr::one(), &sq(&t))))
        }
        (Func::Sin, "sincos") => Some(div(&m2(&fx::sin(a), &fx::cos(a)), &fx::cos(a))),
        (Func::Sin, "cot") => {
            let c = fx::cot(&m2(&half(), a));
            let cond = logic::and(vec![eq(&fx::im(a), &Expr::zero()), eq(&super::modulo(a, &pi()), &Expr::zero())]);
            Some(pw(vec![(Expr::zero(), cond), (div(&m2(&int(2), &c), &add2(&Expr::one(), &sq(&c))), true_())]))
        }
        (Func::Sin, "pow") => Some(rw(&rw(this, "cos"), "pow")),
        (Func::Sin, "sqrt") => Some(rw(&rw(this, "cos"), "sqrt")),
        (Func::Sin, "csc") => Some(inv(&fx::csc(a))),
        (Func::Sin, "sec") => Some(inv(&held(Func::Sec, sub(a, &pi_2)))),
        (Func::Sin, "sinc") => Some(m2(a, &sp(S::Sinc, vec![a.clone()]))),
        (Func::Sin, "besselj") => {
            Some(m2(&sqrt(&div(&m2(&pi(), a), &int(2))), &sp(S::BesselJ, vec![half(), a.clone()])))
        }

        (Func::Cos, "exp") => {
            let a = exp_ready(a);
            let ia = i_times(&a);
            Some(div(&add2(&fx::exp(&ia), &fx::exp(&m2(&neg(&a), &ii()))), &int(2)))
        }
        (Func::Cos, "Pow") => {
            let x = single_arg(a, &Func::Log)?;
            Some(add2(&div(&pow(x, &ii()), &int(2)), &div(&pow(x, &neg(&ii())), &int(2))))
        }
        (Func::Cos, "sin") => Some(held(Func::Sin, add2(a, &pi_2))),
        (Func::Cos, "tan") => {
            let t = sq(&fx::tan(&m2(&half(), a)));
            Some(div(&sub(&Expr::one(), &t), &add2(&Expr::one(), &t)))
        }
        (Func::Cos, "sincos") => Some(div(&m2(&fx::sin(a), &fx::cos(a)), &fx::sin(a))),
        (Func::Cos, "cot") => {
            let c = sq(&fx::cot(&m2(&half(), a)));
            let two_pi = m2(&int(2), &pi());
            let cond = logic::and(vec![eq(&fx::im(a), &Expr::zero()), eq(&super::modulo(a, &two_pi), &Expr::zero())]);
            Some(pw(vec![(Expr::one(), cond), (div(&sub(&c, &Expr::one()), &add2(&c, &Expr::one())), true_())]))
        }
        (Func::Cos, "pow") | (Func::Cos, "sqrt") => cos_sqrt(a),
        (Func::Cos, "sec") => Some(inv(&fx::sec(a))),
        (Func::Cos, "csc") => Some(inv(&rw(&fx::sec(a), "csc"))),
        (Func::Cos, "besselj") => Some(pw(vec![
            (m2(&sqrt(&div(&m2(&pi(), a), &int(2))), &sp(S::BesselJ, vec![neg(&half()), a.clone()])), ne(a, &Expr::zero())),
            (Expr::one(), true_()),
        ])),

        (Func::Tan, "Pow") => {
            let x = single_arg(a, &Func::Log)?;
            let (xi, xmi) = (pow(x, &ii()), pow(x, &neg(&ii())));
            Some(div(&i_times(&sub(&xmi, &xi)), &add2(&xmi, &xi)))
        }
        (Func::Tan, "exp") => {
            let a = exp_ready(a);
            let ia = i_times(&a);
            let (ne_, pe) = (fx::exp(&m2(&neg(&a), &ii())), fx::exp(&ia));
            Some(div(&i_times(&sub(&ne_, &pe)), &add2(&ne_, &pe)))
        }
        (Func::Tan, "sin") => Some(div(&m2(&int(2), &sq(&fx::sin(a))), &fx::sin(&m2(&int(2), a)))),
        (Func::Tan, "cos") => Some(div(&held(Func::Cos, sub(a, &pi_2)), &fx::cos(a))),
        (Func::Tan, "sincos") => Some(div(&fx::sin(a), &fx::cos(a))),
        (Func::Tan, "cot") => Some(inv(&fx::cot(a))),
        (Func::Tan, "sec") => Some(div(&rw(&fx::sin(a), "sec"), &rw(&fx::cos(a), "sec"))),
        (Func::Tan, "csc") => Some(div(&rw(&fx::sin(a), "csc"), &rw(&fx::cos(a), "csc"))),
        (Func::Tan | Func::Cot, "pow" | "sqrt") => {
            let y = rw(&rw(this, "cos"), rule);
            if has_func(&y, &Func::Cos) {
                None
            } else {
                Some(y)
            }
        }
        (Func::Tan, "besselj") => {
            Some(div(&sp(S::BesselJ, vec![half(), a.clone()]), &sp(S::BesselJ, vec![neg(&half()), a.clone()])))
        }

        (Func::Cot, "exp") => {
            let a = exp_ready(a);
            let ia = i_times(&a);
            let (ne_, pe) = (fx::exp(&m2(&neg(&a), &ii())), fx::exp(&ia));
            Some(div(&i_times(&add2(&pe, &ne_)), &sub(&pe, &ne_)))
        }
        (Func::Cot, "Pow") => {
            let x = single_arg(a, &Func::Log)?;
            let (xi, xmi) = (pow(x, &ii()), pow(x, &neg(&ii())));
            Some(div(&neg(&i_times(&add2(&xmi, &xi))), &sub(&xmi, &xi)))
        }
        (Func::Cot, "sin") => Some(div(&fx::sin(&m2(&int(2), a)), &m2(&int(2), &sq(&fx::sin(a))))),
        (Func::Cot, "cos") => Some(div(&fx::cos(a), &held(Func::Cos, sub(a, &pi_2)))),
        (Func::Cot, "sincos") => Some(div(&fx::cos(a), &fx::sin(a))),
        (Func::Cot, "tan") => Some(inv(&fx::tan(a))),
        (Func::Cot, "sec") => Some(div(&rw(&fx::cos(a), "sec"), &rw(&fx::sin(a), "sec"))),
        (Func::Cot, "csc") => Some(div(&rw(&fx::cos(a), "csc"), &rw(&fx::sin(a), "csc"))),
        (Func::Cot, "besselj") => {
            Some(div(&sp(S::BesselJ, vec![neg(&half()), a.clone()]), &sp(S::BesselJ, vec![half(), a.clone()])))
        }
        _ => None,
    }
}

/// `sec`/`csc` rules, including `ReciprocalTrigonometricFunction`'s
/// `_rewrite_reciprocal` for `exp`, `Pow`, `pow` and `sqrt`.
fn reciprocal_trig_method(f: &Func, orig: &Expr, a: &Expr, rule: &str) -> Option<Expr> {
    let pi_2 = half_pi();
    let base = if *f == Func::Sec { Func::Cos } else { Func::Sin };
    match (f, rule) {
        (_, "exp" | "Pow" | "pow" | "sqrt") => {
            let this = f1(base.clone(), orig.clone());
            let t = trig_method(&base, a, rule, &this)?;
            if t == f1(base, a.clone()) {
                None
            } else {
                Some(inv(&t))
            }
        }
        (Func::Sec, "cot") => {
            let c = sq(&fx::cot(&div(a, &int(2))));
            Some(div(&add2(&c, &Expr::one()), &sub(&c, &Expr::one())))
        }
        (Func::Sec, "cos") => Some(inv(&fx::cos(a))),
        (Func::Sec, "sincos") => Some(div(&fx::sin(a), &m2(&fx::cos(a), &fx::sin(a)))),
        (Func::Sec, "sin") => Some(inv(&rw(&fx::cos(a), "sin"))),
        (Func::Sec, "tan") => Some(inv(&rw(&fx::cos(a), "tan"))),
        (Func::Sec, "csc") => Some(held(Func::Csc, sub(&pi_2, a))),
        (Func::Sec, "besselj") => Some(pw(vec![
            (
                inv(&m2(&div(&sqrt(&m2(&pi(), a)), &sqrt(&int(2))), &sp(S::BesselJ, vec![neg(&half()), a.clone()]))),
                ne(a, &Expr::zero()),
            ),
            (Expr::one(), true_()),
        ])),
        (Func::Csc, "sin") => Some(inv(&fx::sin(a))),
        (Func::Csc, "sincos") => Some(div(&fx::cos(a), &m2(&fx::sin(a), &fx::cos(a)))),
        (Func::Csc, "cot") => {
            let c = fx::cot(&div(a, &int(2)));
            Some(div(&add2(&Expr::one(), &sq(&c)), &m2(&int(2), &c)))
        }
        (Func::Csc, "cos") => Some(inv(&rw(&fx::sin(a), "cos"))),
        (Func::Csc, "sec") => Some(held(Func::Sec, sub(&pi_2, a))),
        (Func::Csc, "tan") => Some(inv(&rw(&fx::sin(a), "tan"))),
        (Func::Csc, "besselj") => Some(m2(
            &sqrt(&div(&int(2), &pi())),
            &inv(&m2(&sqrt(a), &sp(S::BesselJ, vec![half(), a.clone()]))),
        )),
        _ => None,
    }
}

// --------------------------------------------------------------------
// cos(p*pi/q) in radicals
// --------------------------------------------------------------------

/// `cos._eval_rewrite_as_sqrt`.
fn cos_sqrt(a: &Expr) -> Option<Expr> {
    let pc = pi_coeff(a)?;
    let k = pc.as_number()?.clone();
    if !k.is_rational() || k.is_integer() {
        return None;
    }
    let (p, q) = k.numer_denom()?;
    let (p64, q64) = (p.to_i64()?, q.to_i64()?);
    let table = match q64 {
        3 => Some(half()),
        5 => Some(div(&add2(&sqrt(&int(5)), &Expr::one()), &int(4))),
        17 => Some(cos_17()),
        // cos(pi/257) in radicals is not built (see the README gap notes)
        257 => return None,
        _ => None,
    };
    if let Some(c) = table {
        return Some(chebyshev_t(p64, &c));
    }
    if q64 % 2 == 0 {
        let pico2 = k.mul(&Number::int(2));
        let nval = rw(&fx::cos(&m2(&Expr::num(pico2.clone()), &pi())), "sqrt");
        let x = pico2.add(&Number::ONE).div(&Number::int(2))?;
        // Python's int() truncates toward zero, `% 2` is non-negative
        let t = x.abs().floor().to_bigint()?;
        let t = if x.is_negative() { -t } else { t };
        let odd = !(&t % BigInt::from(2)).is_zero();
        let v = sqrt(&div(&add2(&Expr::one(), &nval), &int(2)));
        return Some(if odd { neg(&v) } else { v });
    }
    let fc = fermat_coords(q64);
    let denoms: Vec<i64> = match &fc {
        Some(f) => f.clone(),
        None => crate::ntheory::factorint(&q)
            .into_iter()
            .map(|(b, e)| b.to_i64().map(|b| b.pow(e)))
            .collect::<Option<Vec<i64>>>()?,
    };
    let apart = ipartfrac(&denoms)?;
    let vals: Vec<Expr> = apart
        .iter()
        .zip(&denoms)
        .map(|(n, d)| {
            let c = Number::from_bigint(p.clone()).mul(&Number::rational_i64(*n, *d).unwrap());
            m2(&Expr::num(c), &pi())
        })
        .collect();
    let pcls = cos_of_sum(&vals);
    if fc.as_ref().is_none_or(|f| f.len() == 1) {
        Some(pcls)
    } else {
        Some(rw(&pcls, "sqrt"))
    }
}

/// `cos(z0 + z1 + ...)._eval_expand_trig()` with `z_i -> vals[i]`.
fn cos_of_sum(vals: &[Expr]) -> Expr {
    if vals.len() == 1 {
        return fx::cos(&vals[0]);
    }
    let (x, rest) = (&vals[0], &vals[1..]);
    sub(&m2(&fx::cos(x), &cos_of_sum(rest)), &m2(&fx::sin(x), &sin_of_sum(rest)))
}

fn sin_of_sum(vals: &[Expr]) -> Expr {
    if vals.len() == 1 {
        return fx::sin(&vals[0]);
    }
    let (x, rest) = (&vals[0], &vals[1..]);
    add2(&m2(&fx::sin(x), &cos_of_sum(rest)), &m2(&sin_of_sum(rest), &fx::cos(x)))
}

fn fermat_coords(mut n: i64) -> Option<Vec<i64>> {
    let mut primes = Vec::new();
    for p in [3, 5, 17, 257, 65537] {
        if n % p == 0 {
            n /= p;
            primes.push(p);
            if n == 1 {
                return Some(primes);
            }
        }
    }
    None
}

/// SymPy's `igcdex` (extended Euclid, `x*a + y*b = gcd(a, b)`).
fn igcdex(a: i128, b: i128) -> (i128, i128, i128) {
    if a == 0 && b == 0 {
        return (0, 1, 0);
    }
    if a == 0 {
        return (0, b.signum(), b.abs());
    }
    if b == 0 {
        return (a.signum(), 0, a.abs());
    }
    let (mut a, xs) = if a < 0 { (-a, -1) } else { (a, 1) };
    let (mut b, ys) = if b < 0 { (-b, -1) } else { (b, 1) };
    let (mut x, mut y, mut r, mut s) = (1i128, 0i128, 0i128, 1i128);
    while b != 0 {
        let (c, q) = (a % b, a / b);
        (a, b, r, s, x, y) = (b, c, x - q * r, y - q * s, r, s);
    }
    (x * xs, y * ys, a)
}

fn migcdex(x: &[i128]) -> (Vec<i128>, i128) {
    match x.len() {
        0 => (vec![], 0),
        1 => (vec![1], x[0]),
        2 => {
            let (u, v, h) = igcdex(x[0], x[1]);
            (vec![u, v], h)
        }
        _ => {
            let (y, g) = migcdex(&x[1..]);
            let (u, v, h) = igcdex(x[0], g);
            let mut out = vec![u];
            out.extend(y.iter().map(|i| v * i));
            (out, h)
        }
    }
}

fn ipartfrac(denoms: &[i64]) -> Option<Vec<i64>> {
    let denom: i128 = denoms.iter().map(|&d| d as i128).product();
    let a: Vec<i128> = denoms.iter().map(|&d| denom / d as i128).collect();
    migcdex(&a).0.into_iter().map(|v| i64::try_from(v).ok()).collect()
}

/// SymPy's `cos_17()`.
fn cos_17() -> Expr {
    let s = |x: &Expr| sqrt(x);
    let s17 = s(&int(17));
    let s2 = s(&int(2));
    let a = s(&sub(&int(17), &s17));
    let b = s(&add2(&int(17), &s17));
    let inner = add(&[
        m2(&s2, &sub(&m2(&int(-8), &b), &m2(&sub(&Expr::one(), &s17), &a))),
        m2(&int(6), &s17),
        int(34),
    ]);
    let t = add2(&div(&add2(&int(15), &s17), &int(32)), &div(&m2(&s2, &add2(&a, &s(&inner))), &int(32)));
    s(&t)
}

// --------------------------------------------------------------------
// hyperbolic functions
// --------------------------------------------------------------------

fn hyperbolic_method(f: &Func, a: &Expr, rule: &str) -> Option<Expr> {
    let one = Expr::one();
    let ia = i_times(a);
    let hip = half_ipi();
    match (f, rule) {
        (Func::Sinh, "exp") => Some(div(&sub(&fx::exp(a), &fx::exp(&neg(a))), &int(2))),
        (Func::Sinh, "sin") => Some(neg(&i_times(&fx::sin(&ia)))),
        (Func::Sinh, "csc") => Some(div(&neg(&ii()), &fx::csc(&ia))),
        (Func::Sinh, "cosh") => Some(neg(&i_times(&fx::cosh(&add2(a, &hip))))),
        (Func::Sinh, "tanh") => {
            let t = fx::tanh(&m2(&half(), a));
            Some(div(&m2(&int(2), &t), &sub(&one, &sq(&t))))
        }
        (Func::Sinh, "coth") => {
            let c = fx::coth(&m2(&half(), a));
            Some(div(&m2(&int(2), &c), &sub(&sq(&c), &one)))
        }
        (Func::Sinh, "csch") => Some(inv(&sp(S::Csch, vec![a.clone()]))),

        (Func::Cosh, "exp") => Some(div(&add2(&fx::exp(a), &fx::exp(&neg(a))), &int(2))),
        (Func::Cosh, "cos") => Some(held(Func::Cos, ia)),
        (Func::Cosh, "sec") => Some(inv(&held(Func::Sec, ia))),
        (Func::Cosh, "sinh") => Some(neg(&i_times(&held(Func::Sinh, add2(a, &hip))))),
        (Func::Cosh, "tanh") => {
            let t = sq(&fx::tanh(&m2(&half(), a)));
            Some(div(&add2(&one, &t), &sub(&one, &t)))
        }
        (Func::Cosh, "coth") => {
            let c = sq(&fx::coth(&m2(&half(), a)));
            Some(div(&add2(&c, &one), &sub(&c, &one)))
        }
        (Func::Cosh, "sech") => Some(inv(&sp(S::Sech, vec![a.clone()]))),

        (Func::Tanh, "exp") => {
            let (ne_, pe) = (fx::exp(&neg(a)), fx::exp(a));
            Some(div(&sub(&pe, &ne_), &add2(&pe, &ne_)))
        }
        (Func::Tanh, "tan") => Some(neg(&i_times(&held(Func::Tan, ia)))),
        (Func::Tanh, "cot") => Some(div(&neg(&ii()), &held(Func::Cot, ia))),
        (Func::Tanh, "sinh") => Some(div(&i_times(&fx::sinh(a)), &held(Func::Sinh, sub(&hip, a)))),
        (Func::Tanh, "cosh") => Some(div(&i_times(&held(Func::Cosh, sub(&hip, a))), &fx::cosh(a))),
        (Func::Tanh, "coth") => Some(inv(&fx::coth(a))),

        (Func::Coth, "exp") => {
            let (ne_, pe) = (fx::exp(&neg(a)), fx::exp(a));
            Some(div(&add2(&pe, &ne_), &sub(&pe, &ne_)))
        }
        (Func::Coth, "sinh") => Some(div(&neg(&i_times(&held(Func::Sinh, sub(&hip, a)))), &fx::sinh(a))),
        (Func::Coth, "cosh") => Some(div(&neg(&i_times(&fx::cosh(a))), &held(Func::Cosh, sub(&hip, a)))),
        (Func::Coth, "tanh") => Some(inv(&fx::tanh(a))),
        _ => None,
    }
}

/// `sech`/`csch` rules, including `ReciprocalHyperbolicFunction`'s
/// `_rewrite_reciprocal` for `exp`, `tanh` and `coth`.
fn reciprocal_hyperbolic_method(f: &Func, a: &Expr, rule: &str) -> Option<Expr> {
    let ia = i_times(a);
    let hip = half_ipi();
    let is_sech = *f == Func::Special(S::Sech);
    match rule {
        "exp" | "tanh" | "coth" => {
            let base = if is_sech { Func::Cosh } else { Func::Sinh };
            let t = hyperbolic_method(&base, a, rule)?;
            if t == f1(base, a.clone()) {
                None
            } else {
                Some(inv(&t))
            }
        }
        "sin" if !is_sech => Some(div(&ii(), &held(Func::Sin, ia))),
        "csc" if !is_sech => Some(i_times(&held(Func::Csc, ia))),
        "cosh" if !is_sech => Some(div(&ii(), &held(Func::Cosh, add2(a, &hip)))),
        "sinh" if !is_sech => Some(inv(&fx::sinh(a))),
        "cos" if is_sech => Some(inv(&held(Func::Cos, ia))),
        "sec" if is_sech => Some(held(Func::Sec, ia)),
        "sinh" if is_sech => Some(div(&ii(), &held(Func::Sinh, add2(a, &hip)))),
        "cosh" if is_sech => Some(inv(&fx::cosh(a))),
        _ => None,
    }
}

// --------------------------------------------------------------------
// inverse trigonometric and hyperbolic functions
// --------------------------------------------------------------------

fn inverse_method(f: &Func, x: &Expr, rule: &str) -> Option<Expr> {
    let one = Expr::one();
    let i = ii();
    let pi_2 = half_pi();
    let x2 = sq(x);
    let ix = i_times(x);
    let rx = inv(x);
    // sqrt(1 - x**2), sqrt(x**2)/x, x*sqrt(1/x**2)
    let s1mx2 = || sqrt(&sub(&one, &x2));
    let sx2x = || div(&sqrt(&x2), x);
    let xs1x2 = || m2(x, &sqrt(&inv(&x2)));
    let asec = |z: Expr| sp(S::Asec, vec![z]);
    let acsc = |z: Expr| sp(S::Acsc, vec![z]);
    Some(match (f, rule) {
        (Func::Asin, "acos") => sub(&pi_2, &fx::acos(x)),
        (Func::Asin, "atan") => m2(&int(2), &fx::atan(&div(x, &add2(&one, &s1mx2())))),
        (Func::Asin, "log") => neg(&i_times(&fx::log(&add2(&ix, &s1mx2())))),
        (Func::Asin, "acot") => m2(&int(2), &fx::acot(&div(&add2(&one, &s1mx2()), x))),
        (Func::Asin, "asec") => sub(&pi_2, &asec(rx)),
        (Func::Asin, "acsc") => acsc(rx),

        (Func::Acos, "log") => add2(&pi_2, &i_times(&fx::log(&add2(&ix, &s1mx2())))),
        (Func::Acos, "asin") => sub(&pi_2, &fx::asin(x)),
        (Func::Acos, "atan") => add2(&fx::atan(&div(&s1mx2(), x)), &m2(&pi_2, &sub(&one, &xs1x2()))),
        (Func::Acos, "acot") => sub(&pi_2, &m2(&int(2), &fx::acot(&div(&add2(&one, &s1mx2()), x)))),
        (Func::Acos, "asec") => asec(rx),
        (Func::Acos, "acsc") => sub(&pi_2, &acsc(rx)),

        (Func::Atan, "log") => {
            m2(&div(&i, &int(2)), &sub(&fx::log(&sub(&one, &ix)), &fx::log(&add2(&one, &ix))))
        }
        (Func::Atan, "asin") => m2(&sx2x(), &sub(&pi_2, &fx::asin(&inv(&sqrt(&add2(&one, &x2)))))),
        (Func::Atan, "acos") => m2(&sx2x(), &fx::acos(&inv(&sqrt(&add2(&one, &x2))))),
        (Func::Atan, "acot") => fx::acot(&rx),
        (Func::Atan, "asec") => m2(&sx2x(), &asec(sqrt(&add2(&one, &x2)))),
        (Func::Atan, "acsc") => m2(&sx2x(), &sub(&pi_2, &acsc(sqrt(&add2(&one, &x2))))),

        (Func::Acot, "log") => {
            let ir = div(&i, x);
            m2(&div(&i, &int(2)), &sub(&fx::log(&sub(&one, &ir)), &fx::log(&add2(&one, &ir))))
        }
        (Func::Acot, "asin") | (Func::Acot, "acos") => {
            let mx2 = neg(&x2);
            let z = div(&sqrt(&mx2), &sqrt(&sub(&mx2, &one)));
            let core = if rule == "asin" { sub(&pi_2, &fx::asin(&z)) } else { fx::acos(&z) };
            m2(&xs1x2(), &core)
        }
        (Func::Acot, "atan") => fx::atan(&rx),
        (Func::Acot, "asec") => m2(&xs1x2(), &asec(sqrt(&div(&add2(&one, &x2), &x2)))),
        (Func::Acot, "acsc") => m2(&xs1x2(), &sub(&pi_2, &acsc(sqrt(&div(&add2(&one, &x2), &x2))))),

        (Func::Special(S::Asec), "log") => {
            add2(&pi_2, &i_times(&fx::log(&add2(&div(&i, x), &sqrt(&sub(&one, &inv(&x2)))))))
        }
        (Func::Special(S::Asec), "asin") => sub(&pi_2, &fx::asin(&rx)),
        (Func::Special(S::Asec), "acos") => fx::acos(&rx),
        (Func::Special(S::Asec), "atan") | (Func::Special(S::Asec), "acot") => {
            let s = sx2x();
            let z = sqrt(&sub(&x2, &one));
            let g = if rule == "atan" { fx::atan(&z) } else { fx::acot(&inv(&z)) };
            add2(&m2(&pi_2, &sub(&one, &s)), &m2(&s, &g))
        }
        (Func::Special(S::Asec), "acsc") => sub(&pi_2, &acsc(x.clone())),

        (Func::Special(S::Acsc), "log") => {
            neg(&i_times(&fx::log(&add2(&div(&i, x), &sqrt(&sub(&one, &inv(&x2)))))))
        }
        (Func::Special(S::Acsc), "asin") => fx::asin(&rx),
        (Func::Special(S::Acsc), "acos") => sub(&pi_2, &fx::acos(&rx)),
        (Func::Special(S::Acsc), "atan") => m2(&sx2x(), &sub(&pi_2, &fx::atan(&sqrt(&sub(&x2, &one))))),
        (Func::Special(S::Acsc), "acot") => {
            m2(&sx2x(), &sub(&pi_2, &fx::acot(&inv(&sqrt(&sub(&x2, &one))))))
        }
        (Func::Special(S::Acsc), "asec") => sub(&pi_2, &asec(x.clone())),

        (Func::Asinh, "log") => fx::log(&add2(x, &sqrt(&add2(&x2, &one)))),
        (Func::Asinh, "atanh") => fx::atanh(&div(x, &sqrt(&add2(&one, &x2)))),
        (Func::Asinh, "acosh") => {
            let r = div(&sqrt(&sub(&one, &ix)), &sqrt(&sub(&ix, &one)));
            i_times(&sub(&m2(&r, &fx::acosh(&ix)), &pi_2))
        }
        (Func::Asinh, "asin") => neg(&i_times(&held(Func::Asin, ix))),
        (Func::Asinh, "acos") => sub(&i_times(&held(Func::Acos, ix)), &i_times(&pi_2)),

        (Func::Acosh, "log") => {
            fx::log(&add2(x, &m2(&sqrt(&add2(x, &one)), &sqrt(&sub(x, &one)))))
        }
        (Func::Acosh, "acos" | "asin" | "asinh") => {
            let r = div(&sqrt(&sub(x, &one)), &sqrt(&sub(&one, x)));
            let g = match rule {
                "acos" => fx::acos(x),
                "asin" => sub(&pi_2, &fx::asin(x)),
                _ => add2(&pi_2, &i_times(&held(Func::Asinh, ix))),
            };
            m2(&r, &g)
        }
        (Func::Acosh, "atanh") => {
            let sxm1 = sqrt(&sub(x, &one));
            let s1mx = sqrt(&sub(&one, x));
            let sx2m1 = sqrt(&sub(&x2, &one));
            add2(
                &mul(&[pi_2.clone(), sxm1.clone(), inv(&s1mx), sub(&one, &xs1x2())]),
                &mul(&[sxm1, sqrt(&add2(x, &one)), inv(&sx2m1), fx::atanh(&div(&sx2m1, x))]),
            )
        }

        (Func::Atanh, "log") => div(&sub(&fx::log(&add2(&one, x)), &fx::log(&sub(&one, x))), &int(2)),
        (Func::Atanh, "asinh") => {
            let fq = sqrt(&inv(&sub(&x2, &one)));
            sub(
                &div(&m2(&pi(), x), &m2(&int(2), &sqrt(&neg(&x2)))),
                &mul(&[sqrt(&neg(x)), s1mx2(), inv(&sqrt(x)), fq.clone(), fx::asinh(&fq)]),
            )
        }

        (Func::Special(S::Acoth), "log") => {
            div(&sub(&fx::log(&add2(&one, &rx)), &fx::log(&sub(&one, &rx))), &int(2))
        }
        (Func::Special(S::Acoth), "atanh") => fx::atanh(&rx),
        (Func::Special(S::Acoth), "asinh") => {
            let t1 = m2(&sqrt(&div(&sub(x, &one), x)), &sqrt(&div(x, &sub(x, &one))));
            let t2 = m2(&sqrt(&add2(&one, &rx)), &sqrt(&div(x, &add2(x, &one))));
            add2(
                &m2(&half_ipi(), &sub(&t1, &t2)),
                &mul(&[x.clone(), sqrt(&inv(&x2)), fx::asinh(&sqrt(&inv(&sub(&x2, &one))))]),
            )
        }

        (Func::Special(S::Asech), "log") => {
            fx::log(&add2(&rx, &m2(&sqrt(&sub(&rx, &one)), &sqrt(&add2(&rx, &one)))))
        }
        (Func::Special(S::Asech), "acosh") => fx::acosh(&rx),
        (Func::Special(S::Asech), "asinh") => {
            let r = div(&sqrt(&sub(&rx, &one)), &sqrt(&sub(&one, &rx)));
            m2(&r, &add2(&i_times(&held(Func::Asinh, div(&i, x))), &pi_2))
        }
        (Func::Special(S::Asech), "atanh") => {
            let core = add(&[
                one.clone(),
                neg(&m2(&sqrt(x), &sqrt(&rx))),
                neg(&mul(&[div(&i, &int(2)), sqrt(&neg(x)), inv(&sqrt(x))])),
                neg(&mul(&[div(&i, &int(2)), sqrt(&x2), inv(&sqrt(&neg(&x2)))])),
            ]);
            add2(
                &mul(&[i.clone(), pi(), core]),
                &mul(&[sqrt(&inv(&add2(x, &one))), sqrt(&add2(x, &one)), fx::atanh(&s1mx2())]),
            )
        }
        (Func::Special(S::Asech), "acsch") => {
            let r = div(&sqrt(&sub(&rx, &one)), &sqrt(&sub(&one, &rx)));
            m2(&r, &sub(&pi_2, &i_times(&held(Func::Special(S::Acsch), ix))))
        }

        (Func::Special(S::Acsch), "log") => fx::log(&add2(&rx, &sqrt(&add2(&inv(&x2), &one)))),
        (Func::Special(S::Acsch), "asinh") => fx::asinh(&rx),
        (Func::Special(S::Acsch), "acosh") => {
            let ir = div(&i, x);
            let r = div(&sqrt(&sub(&one, &ir)), &sqrt(&sub(&ir, &one)));
            i_times(&sub(&m2(&r, &held(Func::Acosh, ir)), &pi_2))
        }
        (Func::Special(S::Acsch), "atanh") => {
            let a2p1 = add2(&x2, &one);
            let inner = sub(&pi_2, &mul(&[sqrt(&neg(&sq(&a2p1))), inv(&a2p1), fx::atanh(&sqrt(&a2p1))]));
            m2(&div(&sqrt(&neg(&x2)), x), &inner)
        }
        _ => return None,
    })
}

fn atan2_method(y: &Expr, x: &Expr, rule: &str) -> Option<Expr> {
    let r2 = add2(&sq(x), &sq(y));
    match rule {
        "log" => Some(neg(&i_times(&fx::log(&div(&add2(x, &i_times(y)), &sqrt(&r2)))))),
        "atan" => Some(pw(vec![
            (m2(&int(2), &fx::atan(&div(y, &add2(x, &sqrt(&r2))))), ne(y, &Expr::zero())),
            (pi(), lt(&fx::re(x), &Expr::zero())),
            (Expr::zero(), ne(x, &Expr::zero())),
            (Expr::nan(), true_()),
        ])),
        "arg" => {
            if ext_real(x) && ext_real(y) {
                return Some(fx::arg(&add2(x, &i_times(y))));
            }
            let n = add2(x, &i_times(y));
            Some(sub(&fx::arg(&div(&n, &sqrt(&r2))), &i_times(&fx::log(&div(&fx::abs(&n), &sqrt(&fx::abs(&r2)))))))
        }
        _ => None,
    }
}

// --------------------------------------------------------------------
// exp, Pow and the complex/integer helpers
// --------------------------------------------------------------------

/// `arg.coeff(c)` for a product containing every factor of the product `c`.
fn mul_coeff(a: &Expr, factors: &[Expr]) -> Option<Expr> {
    let Kind::Mul(..) = a.kind() else { return None };
    let mut rest = a.args();
    for f in factors {
        let i = rest.iter().position(|r| r == f)?;
        rest.remove(i);
    }
    Some(mul(&rest))
}

fn unary_method(f: &Func, orig: &Expr, a: &Expr, rule: &str) -> Option<Expr> {
    let one = Expr::one();
    match (f, rule) {
        (Func::Exp, "sin") => {
            let ia = i_times(a);
            Some(sub(&fx::sin(&add2(&ia, &half_pi())), &i_times(&fx::sin(&ia))))
        }
        (Func::Exp, "cos") => {
            let ia = i_times(a);
            Some(add2(&fx::cos(&ia), &i_times(&fx::cos(&add2(&ia, &half_pi())))))
        }
        (Func::Exp, "tanh") => {
            let t = fx::tanh(&div(a, &int(2)));
            Some(div(&add2(&one, &t), &sub(&one, &t)))
        }
        (Func::Exp, "sqrt") => {
            let c = mul_coeff(a, &[pi(), ii()])?;
            if c.is_zero() || !c.is_number() {
                return None;
            }
            let pc = m2(&pi(), &c);
            let (co, si) = (fx::cos(&pc), fx::sin(&pc));
            if co.is_func(&Func::Cos) || si.is_func(&Func::Sin) {
                return None;
            }
            Some(add2(&co, &i_times(&si)))
        }
        (Func::Exp, "Pow") => {
            if !a.is_mul() {
                return None;
            }
            let lg = a.args().into_iter().find(|t| single_arg(t, &Func::Log).is_some())?;
            let base = single_arg(&lg, &Func::Log)?.clone();
            Some(pow(&base, &mul_coeff(a, &[lg])?))
        }
        (Func::Re, "im") => Some(sub(orig, &i_times(&fx::im(orig)))),
        (Func::Im, "re") => Some(neg(&i_times(&sub(orig, &fx::re(orig))))),
        (Func::Sign, "Piecewise") if ext_real(a) => Some(pw(vec![
            (one.clone(), gt(a, &Expr::zero())),
            (Expr::int(-1), lt(a, &Expr::zero())),
            (Expr::zero(), true_()),
        ])),
        (Func::Sign, "Heaviside") if ext_real(a) => Some(sub(&m2(&fx::heaviside(a), &int(2)), &one)),
        (Func::Sign, "Abs") => Some(pw(vec![(Expr::zero(), eq(a, &Expr::zero())), (div(a, &fx::abs(a)), true_())])),
        (Func::Abs, "Heaviside") if ext_real(a) => {
            Some(m2(a, &sub(&fx::heaviside(a), &fx::heaviside(&neg(a)))))
        }
        (Func::Abs, "Piecewise") => {
            if ext_real(a) {
                Some(pw(vec![(a.clone(), ge(a, &Expr::zero())), (neg(a), true_())]))
            } else if asm::fact(a, "imaginary") == Some(true) {
                let ia = i_times(a);
                Some(pw(vec![(ia.clone(), ge(&ia, &Expr::zero())), (neg(&ia), true_())]))
            } else {
                None
            }
        }
        (Func::Abs, "sign") => Some(div(a, &fx::sign(a))),
        (Func::Abs, "conjugate") => Some(sqrt(&m2(a, &fx::conjugate(a)))),
        (Func::Arg, "atan2") => Some(fx::atan2(&fx::im(orig), &fx::re(orig))),
        (Func::Floor, "ceiling") => Some(neg(&fx::ceiling(&neg(a)))),
        (Func::Floor, "frac") => Some(sub(a, &sp(S::Frac, vec![a.clone()]))),
        (Func::Ceiling, "floor") => Some(neg(&fx::floor(&neg(a)))),
        (Func::Ceiling, "frac") => Some(add2(a, &sp(S::Frac, vec![neg(a)]))),
        (Func::Special(S::Frac), "floor") => Some(sub(a, &fx::floor(a))),
        (Func::Special(S::Frac), "ceiling") => Some(add2(a, &fx::ceiling(&neg(a)))),
        (Func::Special(S::Sinc), "jn") => Some(sp(S::Jn, vec![Expr::zero(), a.clone()])),
        (Func::Special(S::Sinc), "sin") => {
            Some(pw(vec![(div(&fx::sin(a), a), ne(a, &Expr::zero())), (one, true_())]))
        }
        (Func::Heaviside, _) => heaviside_method(a, &half(), rule),
        (Func::DiracDelta, "Piecewise") => Some(pw(vec![
            (fx::dirac_delta(&Expr::zero()), eq(a, &Expr::zero())),
            (Expr::zero(), true_()),
        ])),
        (Func::Max | Func::Min, _) => minmax_method(f, std::slice::from_ref(a), rule),
        _ => None,
    }
}

fn pow_method(b: &Expr, x: &Expr, rule: &str) -> Option<Expr> {
    if rule != "exp" {
        return None;
    }
    if asm::is_zero(b) == Some(true) || has_func(b, &Func::Exp) || has_func(x, &Func::Exp) {
        return Some(pow(b, x));
    }
    if has_symbol(b) {
        let m = m2(&fx::log(b), x);
        // delay evaluation for a non-symbolic exponent (exp(2*log(x)) -> x**2)
        return Some(if has_symbol(x) { fx::exp(&m) } else { held(Func::Exp, m) });
    }
    Some(fx::exp(&m2(&add2(&fx::log(&fx::abs(b)), &i_times(&fx::arg(b))), x)))
}

fn heaviside_method(a: &Expr, h0: &Expr, rule: &str) -> Option<Expr> {
    let zero = Expr::zero();
    match rule {
        "Piecewise" => Some(if h0.is_zero() {
            pw(vec![(zero.clone(), le(a, &zero)), (Expr::one(), true_())])
        } else if h0.is_one() {
            pw(vec![(zero.clone(), lt(a, &zero)), (Expr::one(), true_())])
        } else {
            pw(vec![(zero.clone(), lt(a, &zero)), (h0.clone(), eq(a, &zero)), (Expr::one(), true_())])
        }),
        "sign" if ext_real(a) => {
            let h = func(Func::Heaviside, if h0 == &half() { vec![zero.clone()] } else { vec![zero.clone(), h0.clone()] });
            let s = div(&add2(&fx::sign(a), &Expr::one()), &int(2));
            let pw1 = pw(vec![(s.clone(), ne(a, &zero)), (h.clone(), true_())]);
            Some(pw(vec![(s, eq(&h, &half())), (pw1, true_())]))
        }
        _ => None,
    }
}

fn minmax_method(f: &Func, args: &[Expr], rule: &str) -> Option<Expr> {
    if args.is_empty() {
        return None;
    }
    let is_max = *f == Func::Max;
    match rule {
        "Abs" => {
            let rest = func(f.clone(), args[1..].to_vec());
            let s = div(&add2(&args[0], &rest), &int(2));
            let d = div(&fx::abs(&sub(&args[0], &rest)), &int(2));
            Some(rw(&if is_max { add2(&s, &d) } else { sub(&s, &d) }, "Abs"))
        }
        "Heaviside" => {
            let terms: Vec<Expr> = args
                .iter()
                .map(|j| {
                    let hs: Vec<Expr> = args
                        .iter()
                        .filter(|i| *i != j)
                        .map(|i| fx::heaviside(&if is_max { sub(j, i) } else { sub(i, j) }))
                        .collect();
                    m2(j, &mul(&hs))
                })
                .collect();
            Some(add(&terms))
        }
        "Piecewise" => {
            let op = |a: &Expr, b: &Expr| if is_max { ge(a, b) } else { le(a, b) };
            let pairs: Vec<(Expr, Expr)> = args
                .iter()
                .enumerate()
                .map(|(i, a)| (a.clone(), logic::and(args[i + 1..].iter().map(|b| op(a, b)).collect())))
                .collect();
            Some(pw(pairs))
        }
        _ => None,
    }
}

// --------------------------------------------------------------------
// Piecewise -> KroneckerDelta
// --------------------------------------------------------------------

fn kd_cond(c: &Expr) -> Option<Expr> {
    let one = Expr::one();
    match c.kind() {
        Kind::Rel(crate::expr::RelOp::Eq, a, b) => Some(super::kronecker_delta(a, b)),
        Kind::Rel(crate::expr::RelOp::Ne, a, b) => Some(sub(&one, &super::kronecker_delta(a, b))),
        Kind::Func(Func::Logic(op @ (LogicOp::And | LogicOp::Or | LogicOp::Not)), args) => {
            let (b1, b2) = match op {
                LogicOp::And => (false, false),
                LogicOp::Or => (true, true),
                _ => (true, false),
            };
            let parts: Vec<Expr> = args.iter().map(kd_cond).collect::<Option<_>>()?;
            let k = if b1 { mul(&parts.iter().map(|p| sub(&one, p)).collect::<Vec<_>>()) } else { mul(&parts) };
            Some(if b2 { sub(&one, &k) } else { k })
        }
        _ => None,
    }
}

fn kd_rule(c: &Expr) -> bool {
    matches!(c.kind(), Kind::Rel(crate::expr::RelOp::Eq | crate::expr::RelOp::Ne, ..))
        || matches!(c.kind(), Kind::Func(Func::Logic(LogicOp::And | LogicOp::Or | LogicOp::Not), _))
}

fn piecewise_as_kronecker(pairs: &[(Expr, Expr)]) -> Option<Expr> {
    let mut conditions = Vec::new();
    let mut true_value: Option<Expr> = None;
    for (v, c) in pairs {
        if kd_rule(c) {
            conditions.push((v.clone(), c.clone()));
        } else if c == &true_() {
            if true_value.is_none() {
                true_value = Some(v.clone());
            }
        } else {
            return None;
        }
    }
    let mut result = true_value?;
    for (v, c) in conditions.iter().rev() {
        let k = kd_cond(c)?;
        result = add2(&m2(&k, v), &m2(&sub(&Expr::one(), &k), &result));
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn r(s: &str, rule: &str) -> String {
        rw(&parse_str(s).unwrap(), rule).to_string()
    }

    #[test]
    fn basic_rewrites() {
        assert_eq!(r("sin(x)", "exp"), "-I*(exp(I*x) - exp(-I*x))/2");
        assert_eq!(r("cos(x)", "sin"), "sin(x + pi/2)");
        assert_eq!(r("tanh(x)", "exp"), "(exp(x) - exp(-x))/(exp(x) + exp(-x))");
        assert_eq!(r("Abs(x)", "sign"), "x/sign(x)");
        assert_eq!(r("Mod(x, 3)", "floor"), "x - 3*floor(x/3)");
    }

    #[test]
    fn cos_in_radicals() {
        assert_eq!(r("sin(pi/5)", "sqrt"), "sqrt(5/8 - sqrt(5)/8)");
        assert_eq!(igcdex(2, 3), (-1, 1, 1));
        assert_eq!(igcdex(100, 2004), (-20, 1, 4));
        assert_eq!(ipartfrac(&[2, 3, 5]).unwrap(), vec![1, 7, -14]);
    }
}
