//! Native `nth_linear_*` detections of `classify_ode`.
//!
//! Ports of SymPy 1.14 `SingleODEProblem.get_linear_coefficients` (with
//! `_lin_eq2dict`), the matchers of the constant-coefficient hints
//! (`nth_linear_constant_coeff_homogeneous`, `..._undetermined_coefficients`,
//! `..._variation_of_parameters`) and of the Euler hints
//! (`nth_linear_euler_eq_homogeneous`, `..._nonhomogeneous_undetermined_coefficients`,
//! `..._nonhomogeneous_variation_of_parameters`), including `_test_term` and
//! the `test` entry of `_undetermined_coefficients_match`.
//!
//! Every function answers `Match`, `NoMatch` (the reference provably rejects)
//! or `Decline` (content whose reference behaviour is not modelled here; the
//! reference matcher must run).

use crate::diff::derivative_unevaluated;
use crate::expr::{Expr, Func, Kind};
use crate::ops::{add, div, mul, mul2, pow};
use crate::simptools::powsimp::{powsimp, Combine, PowsimpOptions};
use crate::traverse::{has_structural as has, preorder};

use super::classify::Outcome;

/// The `nth_linear_*` hint whose matcher is evaluated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NthLinear {
    ConstHomogeneous,
    ConstUndetermined,
    ConstVariation,
    EulerHomogeneous,
    EulerUndetermined,
    EulerVariation,
}

impl NthLinear {
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "nth_linear_constant_coeff_homogeneous" => Self::ConstHomogeneous,
            "nth_linear_constant_coeff_undetermined_coefficients" => Self::ConstUndetermined,
            "nth_linear_constant_coeff_variation_of_parameters" => Self::ConstVariation,
            "nth_linear_euler_eq_homogeneous" => Self::EulerHomogeneous,
            "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients" => Self::EulerUndetermined,
            "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters" => Self::EulerVariation,
            _ => return None,
        })
    }
}

/// Content the structural port does not model: binders (whose bound
/// variables `has_xfree` skips), root objects, relational, boolean and
/// piecewise content, matrix and set objects.
pub(crate) fn unmodelled(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if found {
            return;
        }
        if crate::poly::rootof::is_rootof(s) {
            found = true;
            return;
        }
        match s.kind() {
            Kind::Rel(..) | Kind::Piecewise(_) | Kind::Bool(_) => found = true,
            Kind::Func(f, _) => {
                if matches!(
                    f,
                    Func::Subs
                        | Func::SymbolicDerivative
                        | Func::Limit
                        | Func::Integral
                        | Func::Sum
                        | Func::Product
                        | Func::Order
                        | Func::Logic(_)
                        | Func::Set(_)
                ) || is_matrix_func(f)
                {
                    found = true;
                }
            }
            _ => {}
        }
    });
    found
}

fn is_matrix_func(f: &Func) -> bool {
    matches!(
        f,
        Func::MatrixElement(_)
            | Func::MatrixProxy(_)
            | Func::MatrixTrace(_)
            | Func::MatrixDeterminant(_)
            | Func::MatrixDotProduct(..)
            | Func::Predicate(_)
            | Func::AppliedPredicate
            | Func::GenericExpr
    )
}

fn is_float(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(n) if n.is_float())
}

/// Truth value of the reference `not expr` on an expression: only an exact
/// zero is false. `None` for a float, whose truth the port does not model.
fn falsy(e: &Expr) -> Option<bool> {
    if is_float(e) {
        return None;
    }
    Some(e.is_zero())
}

type Terms = Vec<Option<Vec<Expr>>>;

/// `_lin_eq2dict(a, symset)`: the part of `a` free of `symset` and the
/// coefficient of each entry (by index). `Err` is `PolyNonlinearError`.
fn lin_eq2dict(a: &Expr, symset: &[Expr]) -> Result<(Expr, Vec<Option<Expr>>), ()> {
    if let Some(i) = symset.iter().position(|s| s == a) {
        let mut t = vec![None; symset.len()];
        t[i] = Some(Expr::one());
        return Ok((Expr::zero(), t));
    }
    match a.kind() {
        Kind::Add(..) => {
            let mut lists: Terms = vec![None; symset.len()];
            let mut coeffs = Vec::new();
            for ai in a.args() {
                let (ci, ti) = lin_eq2dict(&ai, symset)?;
                coeffs.push(ci);
                for (j, c) in ti.into_iter().enumerate() {
                    if let Some(c) = c {
                        lists[j].get_or_insert_with(Vec::new).push(c);
                    }
                }
            }
            let terms = lists.into_iter().map(|l| l.map(|cs| add(&cs))).collect();
            Ok((add(&coeffs), terms))
        }
        Kind::Mul(..) => {
            let mut terms: Option<(Vec<Option<Expr>>, Expr)> = None;
            let mut coeffs = Vec::new();
            for ai in a.args() {
                let (ci, ti) = lin_eq2dict(&ai, symset)?;
                if ti.iter().all(Option::is_none) {
                    coeffs.push(ci);
                } else if terms.is_none() {
                    terms = Some((ti, ci));
                } else {
                    return Err(());
                }
            }
            let coeff = mul(&coeffs);
            match terms {
                None => Ok((coeff, vec![None; symset.len()])),
                Some((ti, tc)) => {
                    let ti = ti.into_iter().map(|c| c.map(|c| mul2(&coeff, &c))).collect();
                    Ok((mul2(&coeff, &tc), ti))
                }
            }
        }
        _ => {
            if symset.iter().any(|s| has(a, s)) {
                Err(())
            } else {
                Ok((a.clone(), vec![None; symset.len()]))
            }
        }
    }
}

/// `get_linear_coefficients(eq, f(x), order)`: `[rhs, c_0, ..., c_order]`
/// (the reference dict `{0: c_0, ..., order: c_order, -1: rhs}`), or
/// `NoMatch` for its `None`.
pub fn linear_coefficients(eq: &Expr, fx: &Expr, x: &Expr, order: u32) -> Outcome<Vec<Expr>> {
    if unmodelled(eq) {
        return Outcome::Decline;
    }
    let symset: Vec<Expr> = (0..=order).map(|i| derivative_unevaluated(fx, vec![(x.clone(), i)])).collect();
    let Ok((rhs, terms)) = lin_eq2dict(eq, &symset) else {
        return Outcome::NoMatch;
    };
    if has(&rhs, fx) || terms.iter().flatten().any(|c| has(c, fx)) {
        return Outcome::NoMatch;
    }
    let mut out = vec![rhs];
    out.extend(terms.into_iter().map(|c| c.unwrap_or_else(Expr::zero)));
    Outcome::Match(out)
}

/// The Euler `_test_term(coeff, f(x), order)`.
fn euler_term(coeff: &Expr, fx: &Expr, x: &Expr, order: u32) -> Option<bool> {
    if is_float(coeff) {
        return None;
    }
    if coeff.is_zero() {
        return Some(true);
    }
    if order == 0 {
        return Some(!has(coeff, x));
    }
    let xo = pow(x, &Expr::int(order as i64));
    Some(match coeff.kind() {
        Kind::Mul(..) => !has(coeff, fx) && coeff.args().iter().any(|a| a == &xo),
        Kind::Pow(b, e) => b == x && e == &Expr::int(order as i64),
        _ => order == 1 && coeff == x,
    })
}

fn has_trig(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Func(Func::Sin | Func::Cos, _)) {
            found = true;
        }
    });
    found
}

/// `expr.is_number`: `Some(false)` with a free symbol, `Some(true)` for
/// numbers, constants and known functions of them, `None` otherwise.
fn is_number(e: &Expr) -> Option<bool> {
    let mut symbol = false;
    let mut unknown = false;
    preorder(e, &mut |s| match s.kind() {
        Kind::Sym(_) => symbol = true,
        Kind::Num(_) | Kind::Const(_) | Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => {}
        Kind::Func(f, _) if plain_function(f) || matches!(f, Func::Exp | Func::Sin | Func::Cos | Func::Sinh | Func::Cosh) => {
            if matches!(f, Func::Undefined(_)) {
                unknown = true;
            }
        }
        _ => unknown = true,
    });
    if symbol {
        Some(false)
    } else if unknown {
        None
    } else {
        Some(true)
    }
}

/// Plain `Function` classes outside `(sin, cos, exp, sinh, cosh)`.
fn plain_function(f: &Func) -> bool {
    matches!(
        f,
        Func::Log
            | Func::Tan
            | Func::Cot
            | Func::Sec
            | Func::Csc
            | Func::Asin
            | Func::Acos
            | Func::Atan
            | Func::Acot
            | Func::Tanh
            | Func::Coth
            | Func::Asinh
            | Func::Acosh
            | Func::Atanh
            | Func::Undefined(_)
    )
}

/// `bool(arg.match(a*x + b))` with `a`, `b` excluding `x`, for an argument
/// containing `x`. Decided only for the structural forms whose matching is
/// known: a single term `c*x` (true) and a polynomial with a monomial of
/// degree other than 0 and 1 (no linear form exists).
fn linear_match(arg: &Expr, x: &Expr) -> Option<bool> {
    let terms = if matches!(arg.kind(), Kind::Add(..)) { arg.args() } else { vec![arg.clone()] };
    let mut linear = 0usize;
    let mut nonlinear_monomial = false;
    let mut other = false;
    for t in terms.iter().filter(|t| has(t, x)) {
        let fs = if matches!(t.kind(), Kind::Mul(..)) { t.args() } else { vec![t.clone()] };
        let with_x: Vec<&Expr> = fs.iter().filter(|f| has(f, x)).collect();
        if with_x.len() != 1 {
            other = true;
            continue;
        }
        let f = with_x[0];
        if f == x {
            linear += 1;
        } else if matches!(f.kind(), Kind::Pow(b, e) if b == x && e.is_rational_number()) {
            nonlinear_monomial = true;
        } else {
            other = true;
        }
    }
    if nonlinear_monomial && !other {
        return Some(false);
    }
    if !other && linear == 1 {
        return Some(true);
    }
    None
}

fn all_of(items: impl Iterator<Item = Option<bool>>) -> Option<bool> {
    let mut undecided = false;
    for v in items {
        match v {
            Some(false) => return Some(false),
            None => undecided = true,
            Some(true) => {}
        }
    }
    if undecided {
        None
    } else {
        Some(true)
    }
}

/// The inner `_test_term(expr, x)` of `_undetermined_coefficients_match`.
fn uc_term(e: &Expr, x: &Expr) -> Option<bool> {
    if !has(e, x) {
        return Some(true);
    }
    match e.kind() {
        Kind::Add(..) => all_of(e.args().iter().map(|a| uc_term(a, x))),
        Kind::Mul(..) => {
            let args = e.args();
            if args.iter().filter(|a| has_trig(a)).count() > 1 {
                return Some(false);
            }
            all_of(args.iter().map(|a| uc_term(a, x)))
        }
        Kind::Func(f, args) => {
            if matches!(f, Func::Sin | Func::Cos | Func::Exp | Func::Sinh | Func::Cosh) {
                if args.len() != 1 {
                    return None;
                }
                linear_match(&args[0], x)
            } else if plain_function(f) {
                Some(false)
            } else {
                None
            }
        }
        Kind::Pow(b, ex) => {
            if b.is_symbol() && matches!(ex.kind(), Kind::Num(n) if n.is_integer() && !n.is_negative()) {
                return Some(true);
            }
            match is_number(b)? {
                true => linear_match(ex, x),
                false => Some(false),
            }
        }
        Kind::Sym(_) => Some(true),
        Kind::Derivative(..) => Some(false),
        _ => None,
    }
}

/// `_undetermined_coefficients_match(expr, x)['test']`.
pub fn undetermined_test(expr: &Expr, x: &Expr) -> Option<bool> {
    if unmodelled(expr) {
        return None;
    }
    let opts = PowsimpOptions { deep: false, combine: Combine::Exp, force: false };
    uc_term(&powsimp(expr, &opts), x)
}

/// The Euler undetermined-coefficients test
/// `_undetermined_coefficients_match(posify(rhs.subs(x, exp(x)))...)['test']`,
/// decided for right-hand sides whose terms are products of rationals,
/// symbols, integer powers of symbols and non-negative integer powers of
/// `log(x)`: the substitution maps them onto products of exponentials of
/// `n*x`, powers of `x` and factors free of `x`, all of which pass.
fn euler_undetermined_test(rhs: &Expr, x: &Expr) -> Option<bool> {
    let terms = if matches!(rhs.kind(), Kind::Add(..)) { rhs.args() } else { vec![rhs.clone()] };
    let log_x = Expr::from_kind(Kind::Func(Func::Log, vec![x.clone()]));
    let admitted = |f: &Expr| -> bool {
        match f.kind() {
            Kind::Num(n) => n.is_rational(),
            Kind::Sym(_) => true,
            Kind::Func(Func::Log, _) => f == &log_x,
            Kind::Pow(b, e) => {
                let Kind::Num(n) = e.kind() else { return false };
                n.is_integer() && (b.is_symbol() || (b == &log_x && !n.is_negative()))
            }
            _ => false,
        }
    };
    for t in &terms {
        let fs = if matches!(t.kind(), Kind::Mul(..)) { t.args() } else { vec![t.clone()] };
        if !fs.iter().all(admitted) {
            return None;
        }
    }
    Some(true)
}

/// The matcher of `hint` on `eq` (`eq_high_order_free` or `eq_preprocessed`,
/// as the reference solver uses). A match carries `[r[-1], r[0], ..., r[order]]`,
/// the solver's `self.r`.
pub fn nth_linear(eq: &Expr, fx: &Expr, x: &Expr, order: u32, hint: NthLinear) -> Outcome<Vec<Expr>> {
    if order == 0 {
        return Outcome::NoMatch;
    }
    let m = match linear_coefficients(eq, fx, x, order) {
        Outcome::Match(m) => m,
        other => return other,
    };
    let euler = matches!(hint, NthLinear::EulerHomogeneous | NthLinear::EulerUndetermined | NthLinear::EulerVariation);
    let r = if euler {
        let lead = &m[order as usize + 1];
        if lead.is_zero() || is_float(lead) {
            return Outcome::Decline;
        }
        let factor = div(&pow(x, &Expr::int(order as i64)), lead);
        let r: Vec<Expr> = m.iter().map(|c| mul2(&factor, c)).collect();
        for (i, c) in r.iter().enumerate().skip(1) {
            match euler_term(c, fx, x, (i - 1) as u32) {
                Some(true) => {}
                Some(false) => return Outcome::NoMatch,
                None => return Outcome::Decline,
            }
        }
        r
    } else {
        if m[1..].iter().any(|c| has(c, x)) {
            return Outcome::NoMatch;
        }
        m
    };
    let Some(homogeneous) = falsy(&r[0]) else {
        return Outcome::Decline;
    };
    let accepted = match hint {
        NthLinear::ConstHomogeneous | NthLinear::EulerHomogeneous => Some(homogeneous),
        NthLinear::ConstVariation | NthLinear::EulerVariation => Some(!homogeneous),
        NthLinear::ConstUndetermined if homogeneous => Some(false),
        NthLinear::ConstUndetermined => undetermined_test(&r[0], x),
        NthLinear::EulerUndetermined if homogeneous => Some(false),
        NthLinear::EulerUndetermined => euler_undetermined_test(&r[0], x),
    };
    match accepted {
        Some(true) => Outcome::Match(r),
        Some(false) => Outcome::NoMatch,
        None => Outcome::Decline,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(s: &str) -> Expr {
        parse_str(s).unwrap()
    }

    #[test]
    fn coefficients_and_rejections() {
        let fx = p("f(x)");
        let x = p("x");
        let eq = p("Derivative(f(x), x, 2) + 3*Derivative(f(x), x) + 2*f(x) - exp(x)");
        let Outcome::Match(m) = linear_coefficients(&eq, &fx, &x, 2) else { panic!() };
        assert_eq!(m, vec![p("-exp(x)"), p("2"), p("3"), p("1")]);
        let nonlinear = p("Derivative(f(x), x, 2) + f(x)**2");
        assert!(matches!(linear_coefficients(&nonlinear, &fx, &x, 2), Outcome::NoMatch));
        assert!(matches!(nth_linear(&eq, &fx, &x, 2, NthLinear::ConstUndetermined), Outcome::Match(_)));
        assert!(matches!(nth_linear(&eq, &fx, &x, 2, NthLinear::ConstHomogeneous), Outcome::NoMatch));
        let log_rhs = p("Derivative(f(x), x, 2) + f(x) - log(x)");
        assert!(matches!(nth_linear(&log_rhs, &fx, &x, 2, NthLinear::ConstUndetermined), Outcome::NoMatch));
        let euler = p("x**2*Derivative(f(x), x, 2) + x*Derivative(f(x), x) + f(x) - log(x)");
        assert!(matches!(nth_linear(&euler, &fx, &x, 2, NthLinear::EulerUndetermined), Outcome::Match(_)));
        assert!(matches!(nth_linear(&euler, &fx, &x, 2, NthLinear::ConstVariation), Outcome::NoMatch));
    }
}
