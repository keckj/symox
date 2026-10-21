//! Expression methods used by the simplification ports: `count_ops`,
//! `is_constant` and `equals`.
use super::{args, is_num};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::{assumptions, ops, traverse};
use rustc_hash::FxHashMap;

/// `count_ops(expr)` (scalar expressions).
pub fn count_ops(e: &Expr) -> usize {
    fn count(e: &Expr, memo: &mut FxHashMap<usize, (Expr, usize)>) -> usize {
        if let Some((_, v)) = memo.get(&e.dag_id()) {
            return *v;
        }
        let value = match e.kind() {
            Kind::Num(n) => usize::from(n.is_rational() && n.is_negative()) + usize::from(matches!(n, Number::Rational(_))),
            Kind::Const(Constant::E) => 1,
            Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => 0,
            Kind::Func(Func::MatrixProxy(_), _) => 0,
            Kind::Func(Func::MatrixElement(_), a) => a.iter().map(|x| count(x, memo)).sum(),
            Kind::Add(c, terms) => {
                let mut result = (terms.len() + usize::from(!c.is_zero())).saturating_sub(1);
                let all_negative = (c.is_zero() || c.is_negative()) && terms.iter().all(|(_, k)| k.is_negative());
                result += usize::from(all_negative);
                if !c.is_zero() {
                    result += usize::from(matches!(c, Number::Rational(_)));
                }
                for (term, coefficient) in terms {
                    result += count(&ops::add_term(term, &coefficient.abs()), memo);
                }
                result
            }
            Kind::Mul(c, factors) => {
                let abs = c.abs();
                let mut numerator = Vec::new();
                let mut denominator = Vec::new();
                // `fraction(a)`: negative (or sign-extractable product)
                // exponents, `exp(-y)` included, go to the denominator.
                let neg_exp = |x: &Expr| {
                    x.as_number().is_some_and(|n| n.is_negative())
                        || assumptions::is_negative(x) == Some(true)
                        || (matches!(x.kind(), Kind::Mul(..))
                            && assumptions::is_positive(x) != Some(true)
                            && crate::decompose::could_extract_minus_sign(x))
                };
                for (base, exponent) in factors {
                    if neg_exp(exponent) {
                        denominator.push(ops::pow(base, &ops::neg(exponent)));
                    } else if let (true, Kind::Func(Func::Exp, a)) = (exponent.is_one(), base.kind()) {
                        if a.len() == 1 && neg_exp(&a[0]) {
                            denominator.push(crate::functions::func(Func::Exp, vec![ops::neg(&a[0])]));
                        } else {
                            numerator.push(base.clone());
                        }
                    } else {
                        numerator.push(ops::mul_factor(base, exponent));
                    }
                }
                if abs.is_rational() && !abs.is_integer() {
                    let p = abs.numer();
                    if !p.is_one() {
                        numerator.push(Expr::num(p));
                    }
                    denominator.push(Expr::num(abs.denom()));
                } else if !abs.is_one() {
                    numerator.push(Expr::num(abs));
                }
                usize::from(c.is_negative())
                    + numerator.iter().map(|a| count(a, memo)).sum::<usize>()
                    + denominator.iter().map(|a| count(a, memo)).sum::<usize>()
                    + numerator.len().saturating_sub(1)
                    + denominator.len().saturating_sub(1)
                    + usize::from(!denominator.is_empty())
            }
            Kind::Pow(base, exponent) => {
                if base.is_constant(Constant::E) {
                    1 + count(exponent, memo)
                } else {
                    count(base, memo) + if exponent.is_minus_one() { 1 } else { 1 + count(exponent, memo) }
                }
            }
            Kind::Func(Func::Subs | Func::Limit | Func::Order, a) => a.iter().map(|x| count(x, memo)).sum(),
            Kind::Func(Func::Set(crate::sets::SetOp::IdentityFunction), a) => a.iter().map(|x| count(x, memo)).sum(),
            Kind::Func(Func::Undefined(n), a) if n.is_internal("UnevaluatedExpr") || n.is_internal("Tuple") => a.iter().map(|x| count(x, memo)).sum(),
            _ => 1 + e.args().iter().map(|a| count(a, memo)).sum::<usize>(),
        };
        memo.insert(e.dag_id(), (e.clone(), value));
        value
    }
    count(e, &mut FxHashMap::default())
}

fn subs_all(e: &Expr, syms: &[Expr], v: i64) -> Expr {
    let pairs: Vec<(Expr, Expr)> = syms.iter().map(|s| (s.clone(), Expr::int(v))).collect();
    crate::subs::subs(e, &pairs)
}

/// Numerically certain nonzero value of a number expression.
fn numerically_nonzero(e: &Expr) -> Option<bool> {
    if !is_num(e) {
        return None;
    }
    let v = crate::evalf_mp::evalf(e, 120);
    let mut nonzero = None;
    let parts = if v.is_add() { args(&v) } else { vec![v.clone()] };
    let mut has_float = false;
    for p in parts {
        for f in super::mul_args(&p) {
            if let Some(n) = f.as_number() {
                if !n.is_rational() {
                    has_float = true;
                    if n.to_f64().abs() > 1e-25 {
                        nonzero = Some(true);
                    }
                } else if !n.is_zero() {
                    nonzero = Some(true);
                }
            } else if !f.is_constant(Constant::I) {
                return None;
            }
        }
    }
    let _ = has_float;
    nonzero
}

/// `expr.is_constant(*wrt, simplify=simplify)`.
pub fn is_constant(e: &Expr, wrt: &[Expr]) -> Option<bool> {
    is_constant_with(e, wrt, true)
}

pub fn is_constant_with(e: &Expr, wrt: &[Expr], simplify: bool) -> Option<bool> {
    if is_num(e) {
        return Some(true);
    }
    let free = traverse::free_symbols(e);
    if free.is_empty() {
        return Some(true);
    }
    if !wrt.is_empty() && !wrt.iter().any(|w| free.contains(w)) {
        return Some(true);
    }
    let wrt: Vec<Expr> = if wrt.is_empty() { free.clone() } else { wrt.to_vec() };
    let expr = if simplify { crate::simplify::simplify(e) } else { e.clone() };
    if assumptions::is_zero(&expr) == Some(true) {
        return Some(true);
    }
    let all_free = free.iter().all(|f| wrt.contains(f));
    if all_free {
        let a = subs_all(&expr, &free, 0);
        if !a.is_constant(Constant::NaN) && is_num(&a) {
            let b = subs_all(&expr, &free, 1);
            if !b.is_constant(Constant::NaN) && is_num(&b) && equals(&b, &a) == Some(false) {
                return Some(false);
            }
        }
    }
    for w in &wrt {
        let mut deriv = crate::diff::diff(&expr, w);
        if simplify {
            deriv = crate::simplify::simplify(&deriv);
        }
        if !deriv.is_zero() {
            return Some(false);
        }
    }
    Some(true)
}

fn has_add(e: &Expr) -> bool {
    let mut hit = false;
    super::preorder(e, &mut |a| {
        if a.is_add() {
            hit = true;
        }
    });
    hit
}

fn minpoly_is_zero(e: &Expr) -> Option<bool> {
    let model = crate::poly::algebraic_model::PrimitiveModel::build_radical(e).ok()?;
    let (_, coefficients) = crate::poly::dense::clear_denoms(&model.minimal);
    if coefficients.len() == 2 {
        Some(coefficients.iter().any(|c| c.is_zero()))
    } else {
        Some(false)
    }
}

/// `expr.equals(other)`: `Some(true)`/`Some(false)` when decided.
pub fn equals(a: &Expr, b: &Expr) -> Option<bool> {
    if a == b {
        return Some(true);
    }
    let d = crate::simplify::simplify(&ops::sub(a, b));
    let diff = super::exprtools::factor_terms(&d, true, false, true, false);
    if diff.is_zero() {
        return Some(true);
    }
    if !has_add(&diff) {
        return Some(false);
    }
    let factors = super::mul_args(&ops::as_coeff_mul(&diff).1);
    if factors.len() > 1 {
        let z: Vec<Option<bool>> = factors.iter().map(|f| equals(f, &Expr::zero())).collect();
        if z.iter().all(Option::is_some) {
            return Some(z.iter().any(|v| *v == Some(true)));
        }
    }
    let constant = is_constant_with(&diff, &[], false);
    if constant == Some(false) {
        return Some(false);
    }
    if !is_num(&diff) {
        return None;
    }
    if numerically_nonzero(&diff) == Some(true) {
        return Some(false);
    }
    // exact rewrites standing in for SymPy's surd solving and minimal_polynomial
    let denested = super::mexpand(&super::sqrtdenest::sqrtdenest(&diff, 3));
    if denested.is_zero() || super::mexpand(&super::radsimp::radsimp(&denested, true, 4)).is_zero() {
        return Some(true);
    }
    if let Some(z) = minpoly_is_zero(&diff) {
        return Some(z);
    }
    None
}

fn has_func(e: &Expr) -> bool {
    let mut hit = false;
    super::preorder(e, &mut |a| {
        if matches!(a.kind(), Kind::Func(..) | Kind::Piecewise(..) | Kind::Rel(..) | Kind::Bool(_)) {
            hit = true;
        }
    });
    hit
}

/// `together(expr, deep=True)` (gcd_terms based, as in SymPy).
fn together_deep(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => e.clone(),
        Kind::Add(..) => {
            let parts: Vec<Expr> = args(e).iter().map(together_deep).collect();
            super::exprtools::gcd_terms_list(&parts, false, true, true)
        }
        Kind::Pow(b, x) => ops::pow(&together_deep(b), &together_deep(x)),
        _ => super::map_args(e, &mut together_deep),
    }
}

fn shorter(choices: &[Expr]) -> Expr {
    if choices.iter().all(|c| c == &choices[0]) {
        return choices[0].clone();
    }
    let mut best = choices[0].clone();
    let mut cost = count_ops(&best);
    for c in &choices[1..] {
        let k = count_ops(c);
        if k < cost {
            best = c.clone();
            cost = k;
        }
    }
    best
}

fn cancel(e: &Expr) -> Expr {
    crate::poly::api::cancel(e).unwrap_or_else(|_| e.clone())
}

/// Re-evaluate `Number*Add` products held unevaluated (`hollow_mul`).
fn hollow_mul(e: &Expr) -> Expr {
    super::bottom_up(e, &mut |w| {
        if w.is_mul() {
            let a = args(w);
            if a.len() == 2 && super::is_number_atom(&a[0]) && a[1].is_add() {
                return ops::mul(&a);
            }
        }
        w.clone()
    })
}

/// `simplify(g, ratio=oo)` for the rational functions met by `hypersimp`:
/// the steps of SymPy's `simplify` that act on an expression free of
/// functions. Other expressions use symox's general `simplify`.
pub fn simplify_rational(g: &Expr) -> Expr {
    if has_func(g) {
        return crate::simplify::simplify(g);
    }
    use super::powsimp::{powsimp, Combine, PowsimpOptions};
    let mut expr = super::misc::signsimp(g, true);
    if expr.is_atom() {
        return expr;
    }
    expr = super::bottom_up(&expr, &mut super::normal);
    let p = powsimp(&expr, &PowsimpOptions { deep: false, combine: Combine::All, force: false });
    expr = match crate::decompose::as_content_primitive(&p, false, true) {
        Ok((c, prim)) => ops::mul(&[Expr::num(c), prim]),
        Err(_) => p,
    };
    let e1 = cancel(&expr);
    let expr1 = shorter(&[e1.clone(), cancel(&super::mexpand(&e1))]);
    expr = shorter(&[together_deep(&expr), together_deep(&expr1)]);
    expr = super::exprtools::factor_terms(&expr, false, false, false, false);
    let mut short = shorter(&[expr.clone(), powsimp(&expr, &PowsimpOptions { deep: true, combine: Combine::Exp, force: false })]);
    short = shorter(&[short.clone(), cancel(&short)]);
    short = shorter(&[short.clone(), super::exprtools::factor_terms(&short, false, false, false, true), crate::simp::expand_power_exp(&crate::simp::expand_mul(&short))]);
    expr = hollow_mul(&short);
    let (numer, denom) = super::as_numer_denom(&expr);
    if denom.is_add() {
        let (n, d) = super::radsimp::fraction(&super::radsimp::radsimp(&ops::inv(&denom), false, 1), false);
        if !n.is_one() {
            expr = ops::div(&crate::expand::expand(&ops::mul(&[numer, n])), &d);
        }
    }
    if ops::could_extract_minus_sign(&expr) {
        let (n, d) = super::radsimp::fraction(&expr, false);
        if !d.is_zero() {
            expr = super::misc::signsimp(&ops::div(&ops::neg(&n), &ops::neg(&d)), true);
        }
    }
    expr
}

/// `Basic.compare`: -1, 0 or 1 following SymPy's canonical args ordering.
pub fn compare(a: &Expr, b: &Expr) -> i32 {
    match crate::traverse::args_cmp(a, b) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}
