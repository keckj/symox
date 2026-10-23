//! Native port of the `sympy.simplify.trigsimp` drivers: `futrig` (the
//! default `trigsimp` method), `exptrigsimp` and `_trigsimp_inverse`.
use super::exprtools::factor_terms;
use super::fu::{self, cos, is_trig, Less};
use super::methods::count_ops;
use super::{args, as_coeff_mul, bottom_up, expand_mul, map_args, mexpand, mul_args, nodes, normal};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::func;
use crate::ops;

/// SymPy's `HyperbolicFunction`, including the reciprocal `sech`, `csch`.
fn is_hyper(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth | Func::Special(crate::special::SpecialFn::Sech) | Func::Special(crate::special::SpecialFn::Csch), a) if a.len() == 1)
}

fn has(e: &Expr, pred: &dyn Fn(&Expr) -> bool) -> bool {
    let mut hit = false;
    super::preorder(e, &mut |a| {
        if !hit && pred(a) {
            hit = true;
        }
    });
    hit
}

fn has_trig(e: &Expr) -> bool {
    has(e, &is_trig)
}

fn has_hyper(e: &Expr) -> bool {
    has(e, &is_hyper)
}

fn factor(e: &Expr) -> Expr {
    fu::factor(e)
}

/// `isinstance(e, Expr)`.
fn is_expr(e: &Expr) -> bool {
    match e.kind() {
        Kind::Rel(..) | Kind::Bool(_) => false,
        Kind::Func(Func::Logic(_) | Func::Set(_), _) => false,
        _ => true,
    }
}

/// `_is_Expr(e)`: an expression tree made of expressions only.
fn is_expr_tree(e: &Expr) -> bool {
    match e.kind() {
        Kind::Piecewise(_) => false,
        Kind::Derivative(..) => is_expr(e),
        _ => is_expr(e) && (e.is_atom() || args(e).iter().all(is_expr_tree)),
    }
}

/// `_eapply(func, e, cond)`.
fn eapply(f: &dyn Fn(&Expr) -> Expr, e: &Expr, cond: bool) -> Expr {
    if !is_expr(e) {
        return e.clone();
    }
    if is_expr_tree(e) || e.is_atom() {
        return f(e);
    }
    if matches!(e.kind(), Kind::Piecewise(_)) {
        // the arguments are ExprCondPairs, which `_eapply` leaves alone
        return e.clone();
    }
    map_args(e, &mut |a| if !cond || has_trig(a) { eapply(f, a, false) } else { a.clone() })
}

#[derive(Clone, Copy)]
enum R {
    Id,
    T1,
    T2,
    T2i,
    T3,
    T5,
    T6,
    T10,
    T10i,
    T11,
    T11r,
    T12,
    T14,
    T111,
    Morrie,
    T8,
    Factor,
    Mexpand,
    FactorNormal,
    T2iT2,
    EmT5,
    EmT15,
    EmT6,
    EmT16,
    EmT22,
    FtT12,
}

enum S {
    R(R),
    Chain(Vec<S>),
    Min(Vec<S>),
}

fn apply(r: R, e: &Expr) -> Expr {
    let four = Expr::int(4);
    match r {
        R::Id => e.clone(),
        R::T1 => fu::tr1(e),
        R::T2 => fu::tr2(e),
        R::T2i => fu::tr2i(e, false),
        R::T3 => fu::tr3(e),
        R::T5 => fu::tr5(e, &four, false),
        R::T6 => fu::tr6(e, &four, false),
        R::T10 => fu::tr10(e, true),
        R::T10i => fu::tr10i(e),
        R::T11 => fu::tr11(e, None),
        R::T11r => fu::tr11_ratio(e),
        R::T12 => fu::tr12(e, true),
        R::T14 => fu::tr14(e, true),
        R::T111 => fu::tr111(e),
        R::Morrie => fu::trmorrie(e),
        R::T8 => fu::tr8(e, true),
        R::Factor => eapply(&factor, e, true),
        R::Mexpand => eapply(&mexpand, e, true),
        R::FactorNormal => eapply(&|i| factor(&normal(i)), e, true),
        R::T2iT2 => fu::tr2i(&fu::tr2(e), false),
        R::EmT5 => eapply(&expand_mul, &fu::tr5(e, &four, false), true),
        R::EmT15 => eapply(&expand_mul, &fu::tr15(e, &four, false), true),
        R::EmT6 => eapply(&expand_mul, &fu::tr6(e, &four, false), true),
        R::EmT16 => eapply(&expand_mul, &fu::tr16(e, &four, false), true),
        R::EmT22 => eapply(&expand_mul, &fu::tr22(e, &four, false), true),
        R::FtT12 => eapply(&|i| factor_terms(i, false, false, false, true), &fu::tr12(e, true), true),
    }
}

/// `(L(x), x.count_ops(), _nodes(x), len(x.args), x.is_Add)`.
fn lops_less(a: &Expr, b: &Expr) -> bool {
    let key = |x: &Expr| (fu::l_count(x), count_ops(x), nodes(x), if x.is_atom() { 0 } else { args(x).len() }, x.is_add());
    key(a) < key(b)
}

fn run(s: &S, e: &Expr, less: Less) -> Expr {
    match s {
        S::R(r) => apply(*r, e),
        S::Chain(v) => {
            let mut x = e.clone();
            for s in v {
                x = run(s, &x, less);
            }
            x
        }
        S::Min(v) => fu::min_by_measure(v.iter().map(|s| run(s, e, less)).collect(), less),
    }
}

fn futrig_tree() -> S {
    use R::*;
    let r = S::R;
    let alt = |a: R, b: R| S::Min(vec![S::R(a), S::R(b)]);
    S::Min(vec![
        r(Id),
        S::Chain(vec![
            r(T3),
            r(T1),
            r(T12),
            r(Factor),
            r(T2),
            alt(Id, Mexpand),
            r(T2i),
            r(FactorNormal),
            r(T14),
            r(T5),
            r(T10),
            r(T11),
            r(T11r),
            r(T6),
            r(Factor),
            r(T14),
            alt(Id, Mexpand),
            r(T10i),
            r(Morrie),
            alt(Id, T8),
            alt(Id, T2iT2),
            alt(EmT5, EmT15),
            alt(EmT6, EmT16),
            r(T111),
            alt(Id, T2i),
            alt(Id, EmT22),
            r(T1),
            r(T2),
            r(T2i),
            alt(Id, FtT12),
        ]),
    ])
}

fn futrig_node(e: &Expr) -> Expr {
    if !has_trig(e) {
        return e.clone();
    }
    let (coeff, e) = if e.is_mul() {
        let (dep, indep): (Vec<Expr>, Vec<Expr>) = mul_args(e).into_iter().partition(has_trig);
        (Some(ops::mul(&indep)), super::mul_from_args(&dep))
    } else {
        (None, e.clone())
    };
    let r = run(&futrig_tree(), &e, &lops_less);
    match coeff {
        Some(c) => ops::mul2(&c, &r),
        None => r,
    }
}

/// `futrig(e, hyper)`: the default `trigsimp` method.
pub fn futrig(e: &Expr, hyper: bool) -> Expr {
    if e.is_atom() {
        return e.clone();
    }
    let old = e.clone();
    let mut e = bottom_up(e, &mut futrig_node);
    if hyper && has_hyper(&e) {
        let (t, d, back) = fu::hyper_as_trig(&e);
        e = fu::hyper_as_trig_back(&bottom_up(&t, &mut futrig_node), &d, &back);
    }
    if e != old && e.is_mul() && args(&e)[0].is_rational_number() {
        let (c, m) = as_coeff_mul(&e, false);
        e = ops::mul2(&c, &m);
    }
    e
}

// ---------------------------------------------------------------------------
// exptrigsimp

fn exp(a: &Expr) -> Expr {
    func(Func::Exp, vec![a.clone()])
}

fn held(f: Func, a: Expr) -> Expr {
    Expr::from_unevaluated_kind(Kind::Func(f, vec![a]))
}

/// `e.rewrite(exp)` for the trigonometric and hyperbolic functions.
pub fn rewrite_exp(e: &Expr) -> Expr {
    bottom_up(e, &mut |rv| {
        let Kind::Func(f, a) = rv.kind() else { return rv.clone() };
        if a.len() != 1 {
            return rv.clone();
        }
        let x = &a[0];
        let i = Expr::i();
        let two = Expr::int(2);
        let ix = ops::mul2(x, &i);
        let (pos, neg) = (exp(&ix), exp(&ops::neg(&ix)));
        let (hp, hn) = (exp(x), exp(&ops::neg(x)));
        match f {
            Func::Sin => ops::div(&ops::sub(&pos, &neg), &ops::mul2(&two, &i)),
            Func::Cos => ops::div(&ops::add2(&pos, &neg), &two),
            Func::Tan => ops::div(&ops::mul2(&i, &ops::sub(&neg, &pos)), &ops::add2(&neg, &pos)),
            Func::Cot => ops::div(&ops::mul2(&i, &ops::add2(&pos, &neg)), &ops::sub(&pos, &neg)),
            Func::Sec => ops::div(&Expr::one(), &ops::div(&ops::add2(&pos, &neg), &two)),
            Func::Csc => ops::div(&Expr::one(), &ops::div(&ops::sub(&pos, &neg), &ops::mul2(&two, &i))),
            Func::Sinh => ops::div(&ops::sub(&hp, &hn), &two),
            Func::Cosh => ops::div(&ops::add2(&hp, &hn), &two),
            Func::Tanh => ops::div(&ops::sub(&hp, &hn), &ops::add2(&hp, &hn)),
            Func::Coth => ops::div(&ops::add2(&hp, &hn), &ops::sub(&hp, &hn)),
            _ => rv.clone(),
        }
    })
}

/// `e.rewrite(cos)`.
pub fn rewrite_cos(e: &Expr) -> Expr {
    bottom_up(e, &mut |rv| {
        let Kind::Func(f, a) = rv.kind() else { return rv.clone() };
        if a.len() != 1 {
            return rv.clone();
        }
        let x = &a[0];
        let i = Expr::i();
        let shifted = || held(Func::Cos, ops::sub(x, &ops::div(&Expr::pi(), &Expr::int(2))));
        match f {
            Func::Exp => {
                let ix = ops::mul2(&i, x);
                ops::add2(&cos(&ix), &ops::mul2(&i, &cos(&ops::add2(&ix, &ops::div(&Expr::pi(), &Expr::int(2))))))
            }
            Func::Sin => shifted(),
            Func::Tan => ops::div(&shifted(), &cos(x)),
            Func::Cot => ops::div(&cos(x), &shifted()),
            Func::Sec => ops::div(&Expr::one(), &cos(x)),
            Func::Csc => ops::div(&Expr::one(), &shifted()),
            Func::Cosh => held(Func::Cos, ops::mul2(&i, x)),
            _ => rv.clone(),
        }
    })
}

fn exp_trig(e: &Expr) -> Expr {
    let mut choices = vec![e.clone()];
    if has(e, &|a| is_trig(a) || is_hyper(a)) {
        choices.push(rewrite_exp(e));
    }
    choices.push(rewrite_cos(e));
    fu::min_by_measure(choices, &|a: &Expr, b: &Expr| count_ops(a) < count_ops(b))
}

/// `signlog(expr, sign)` of `exptrigsimp`.
fn signlog(e: &Expr, sign: i64) -> Option<(i64, Expr)> {
    if e.is_constant(Constant::E) {
        return Some((sign, Expr::one()));
    }
    if let Kind::Func(Func::Exp, a) = e.kind() {
        return Some((sign, a[0].clone()));
    }
    if let Kind::Pow(b, x) = e.kind() {
        if b.is_constant(Constant::E) {
            return Some((sign, x.clone()));
        }
    }
    if sign == 1 {
        return signlog(&ops::neg(e), -1);
    }
    None
}

type DD = Vec<(Expr, Expr)>;

fn dd_pos(d: &mut DD, k: &Expr) -> usize {
    match d.iter().position(|(a, _)| a == k) {
        Some(p) => p,
        None => {
            d.push((k.clone(), Expr::zero()));
            d.len() - 1
        }
    }
}

fn dd_add(d: &mut DD, k: &Expr, v: &Expr) {
    let p = dd_pos(d, k);
    d[p].1 = ops::add2(&d[p].1, v);
}

fn exptrig_f(rv: &Expr) -> Expr {
    if !rv.is_mul() {
        return rv.clone();
    }
    let rvd: DD = crate::decompose::as_powers_dict(rv);
    let mut newd = rvd.clone();
    let e1 = Expr::e();
    let mut ee = rvd.iter().find(|(k, _)| *k == e1).map(|(_, v)| v.clone()).unwrap_or_else(Expr::zero);
    for (k, m) in &rvd {
        if !(k.is_add() && args(k).len() == 2) {
            continue;
        }
        let ka = args(k);
        let c = ka[0].clone();
        let Some((sign, x)) = signlog(&ops::div(&ka[1], &c), 1) else { continue };
        if x.is_zero() {
            continue;
        }
        let p = dd_pos(&mut newd, k);
        newd[p].1 = ops::sub(&newd[p].1, m);
        let two = Expr::int(2);
        let ex = ops::pow(&e1, &x);
        let sgn = Expr::int(sign);
        if ee == ops::div(&ops::mul2(&ops::neg(&x), m), &two) {
            let p = dd_pos(&mut newd, &e1);
            newd[p].1 = ops::sub(&newd[p].1, &ee);
            ee = Expr::zero();
            let half = ops::div(&x, &two);
            if sign == 1 {
                dd_add(&mut newd, &ops::mul(&[two, c.clone(), func(Func::Cosh, vec![half])]), m);
            } else {
                dd_add(&mut newd, &ops::mul(&[Expr::int(-2), c.clone(), func(Func::Sinh, vec![half])]), m);
            }
        } else {
            let key = ops::sub(&Expr::one(), &ops::mul2(&sgn, &ex));
            let p = dd_pos(&mut newd, &key);
            if newd[p].1 == ops::neg(m) {
                newd.remove(p);
                let t = func(Func::Tanh, vec![ops::div(&x, &two)]);
                if sign == 1 {
                    dd_add(&mut newd, &ops::div(&ops::neg(&c), &t), m);
                } else {
                    dd_add(&mut newd, &ops::mul2(&ops::neg(&c), &t), m);
                }
            } else if crate::assumptions::is_integer(m) == Some(true)
                || crate::assumptions::is_positive(&c) == Some(true)
            {
                dd_add(&mut newd, &ops::add2(&Expr::one(), &ops::mul2(&sgn, &ex)), m);
                dd_add(&mut newd, &c, m);
            } else {
                // Intentional correction: the reference also writes
                // (c*(1 + s*E**x))**m as c**m*(1 + s*E**x)**m here, which is
                // wrong for a non-integer m and a c that is not positive
                // ((exp(a) - 1)**(1/3) is not (-1)**(1/3)*(1 - exp(a))**(1/3)
                // for a > 0); the factor is kept.
                dd_add(&mut newd, k, m);
            }
        }
    }
    ops::mul(&newd.iter().map(|(k, v)| ops::pow(k, v)).collect::<Vec<_>>())
}

fn has_i(e: &Expr) -> bool {
    has(e, &|a| a.is_constant(Constant::I))
}

/// `exptrigsimp(expr)`.
pub fn exptrigsimp(expr: &Expr) -> Expr {
    let mut newexpr = bottom_up(expr, &mut exp_trig);
    newexpr = bottom_up(&newexpr, &mut exptrig_f);
    if has_hyper(&newexpr) {
        let (e, d, back) = fu::hyper_as_trig(&newexpr);
        newexpr = fu::hyper_as_trig_back(&fu::tr2i(&e, false), &d, &back);
    }
    if has_trig(&newexpr) {
        newexpr = fu::tr2i(&newexpr, false);
    }
    if !(has_i(&newexpr) && !has_i(expr)) {
        return newexpr;
    }
    expr.clone()
}

// ---------------------------------------------------------------------------
// _trigsimp_inverse

fn coeff_isneg(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_negative(),
        Kind::Mul(c, _) => c.is_negative(),
        _ => false,
    }
}

fn inverse_f(rv: &Expr) -> Expr {
    let Kind::Func(f, a) = rv.kind() else { return rv.clone() };
    let inner = match f {
        Func::Asin => Some(Func::Sin),
        Func::Acos => Some(Func::Cos),
        Func::Atan => Some(Func::Tan),
        Func::Acot => Some(Func::Cot),
        // SymPy: `g = rv.inverse` with `g()(1)` a TrigonometricFunction,
        // which includes the reciprocal functions sec and csc.
        Func::Special(crate::special::SpecialFn::Asec) => Some(Func::Sec),
        Func::Special(crate::special::SpecialFn::Acsc) => Some(Func::Csc),
        _ => None,
    };
    if let Some(g) = inner {
        if let Kind::Func(h, b) = a[0].kind() {
            if *h == g && b.len() == 1 {
                return b[0].clone();
            }
        }
    }
    if *f == Func::Atan2 && a.len() == 2 {
        let (y, x) = (a[0].clone(), a[1].clone());
        if coeff_isneg(&y) {
            return ops::neg(&inverse_f(&func(Func::Atan2, vec![ops::neg(&y), x])));
        }
        if coeff_isneg(&x) {
            return ops::sub(&Expr::pi(), &inverse_f(&func(Func::Atan2, vec![y, ops::neg(&x)])));
        }
        let fx = |e: &Expr, g: Func| matches!(e.kind(), Kind::Func(h, b) if *h == g && b.len() == 1);
        if !x.is_atom() && !y.is_atom() && args(&x).first() == args(&y).first() {
            let xa = args(&x)[0].clone();
            if fx(&y, Func::Sin) && fx(&x, Func::Cos) {
                return xa;
            }
            if fx(&y, Func::Cos) && fx(&x, Func::Sin) {
                return ops::sub(&ops::div(&Expr::pi(), &Expr::int(2)), &xa);
            }
        }
    }
    rv.clone()
}

/// `_trigsimp_inverse(rv)`: asin(sin(x)) -> x and atan2 of sin/cos.
pub fn trigsimp_inverse(rv: &Expr) -> Expr {
    bottom_up(rv, &mut inverse_f)
}

/// Default `trigsimp(expr)` (`method='matching'`).
pub fn trigsimp(e: &Expr) -> Expr {
    futrig(e, true)
}

