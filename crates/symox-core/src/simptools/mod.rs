//! Native ports of SymPy's structural simplification helpers.
//!
//! The submodules follow `sympy.core.exprtools`, `sympy.simplify.radsimp`,
//! `powsimp`, `sqrtdenest`, `gammasimp`, `ratsimp`, `fu` and the generic
//! `Expr` helpers (`coeff`, `extract_multiplicatively`, ...) they rely on.
//! Unevaluated products and sums are built exactly where SymPy builds them
//! (`_keep_coeff`, `_unevaluated_Mul`, `Add._from_args`), so results print
//! like SymPy's.
pub mod exprtools;
pub mod radsimp;
pub mod powsimp;
pub mod misc;
pub mod sqrtdenest;
pub mod gammasimp;
pub mod ratsimp;
pub mod fu;
pub mod trigsimp;
pub mod methods;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::{assumptions, decompose, ops, traverse};

/// SymPy's public `.args`.
pub(crate) fn args(e: &Expr) -> Vec<Expr> {
    traverse::ordered_args(e)
}

/// `Add.make_args`.
pub(crate) fn add_args(e: &Expr) -> Vec<Expr> {
    if e.is_add() { args(e) } else { vec![e.clone()] }
}

/// `Mul.make_args`.
pub(crate) fn mul_args(e: &Expr) -> Vec<Expr> {
    if e.is_mul() { args(e) } else { vec![e.clone()] }
}

/// `expr.is_Number`: numeric atoms including the infinities.
pub(crate) fn is_number_atom(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(c) => matches!(c, Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN),
        _ => false,
    }
}

/// SymPy's `expr.is_number`: no free symbols and no undefined functions.
pub(crate) fn is_num(e: &Expr) -> bool {
    is_number_atom(e) || crate::functions::is_numeric_value(e)
}

fn factor_pair(a: &Expr) -> (Expr, Expr) {
    match a.kind() {
        Kind::Pow(b, e) => (b.clone(), e.clone()),
        _ => (a.clone(), Expr::one()),
    }
}

/// `Mul._from_args(args)`: a product node holding exactly `args`.
pub fn mul_from_args(values: &[Expr]) -> Expr {
    match values.len() {
        0 => return Expr::one(),
        1 => return values[0].clone(),
        _ => {}
    }
    let (coeff, rest) = match values[0].as_number() {
        Some(n) => (n.clone(), &values[1..]),
        None => (Number::ONE, values),
    };
    let canonical = ops::mul(values);
    if traverse::ordered_args(&canonical) == values {
        return canonical;
    }
    let mut factors: Vec<(Expr, Expr)> = rest.iter().map(factor_pair).collect();
    let distinct = {
        let mut bases: Vec<&Expr> = factors.iter().map(|f| &f.0).collect();
        let n = bases.len();
        bases.sort();
        bases.dedup();
        bases.len() == n
    };
    if distinct && rest.iter().all(|a| a.as_number().is_none() && !a.is_mul()) {
        factors.sort_unstable_by(ops::cmp_factor);
        return Expr::from_kind(Kind::Mul(coeff, factors));
    }
    Expr::from_unevaluated_kind(Kind::Mul(coeff, rest.iter().map(|a| (a.clone(), Expr::one())).collect()))
}

/// `Add._from_args(args)`: a sum node holding exactly `args`.
pub fn add_from_args(values: &[Expr]) -> Expr {
    match values.len() {
        0 => return Expr::zero(),
        1 => return values[0].clone(),
        _ => {}
    }
    let canonical = ops::add(values);
    let mut a = traverse::ordered_args(&canonical);
    let mut b = values.to_vec();
    a.sort();
    b.sort();
    if a == b {
        return canonical;
    }
    let mut coeff = Number::ZERO;
    let mut terms: Vec<(Expr, Number)> = Vec::new();
    for v in values {
        if let Some(n) = v.as_number() {
            coeff = coeff.add(n);
            continue;
        }
        let (k, t) = ops::as_coeff_mul(v);
        terms.push((t, k));
    }
    let n = terms.len();
    let mut keys: Vec<&Expr> = terms.iter().map(|t| &t.0).collect();
    keys.sort();
    keys.dedup();
    if keys.len() == n && !terms.iter().any(|(t, _)| t.is_add()) {
        terms.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        return Expr::from_kind(Kind::Add(coeff, terms));
    }
    Expr::from_unevaluated_kind(Kind::Add(Number::ZERO, values.iter().map(|v| (v.clone(), Number::ONE)).collect()))
}

/// `Mul(coeff, factors, evaluate=False)` with a numeric coefficient.
pub fn held_coeff_mul(coeff: &Number, factors: &Expr) -> Expr {
    if coeff.is_one() {
        return factors.clone();
    }
    match factors.kind() {
        Kind::Mul(c, fs) if c.is_one() => Expr::from_kind(Kind::Mul(coeff.clone(), fs.clone())),
        Kind::Mul(..) => {
            let mut v = vec![Expr::num(coeff.clone())];
            v.extend(args(factors));
            Expr::from_unevaluated_kind(Kind::Mul(Number::ONE, v.into_iter().map(|a| (a, Expr::one())).collect()))
        }
        _ => {
            let (b, e) = factor_pair(factors);
            Expr::from_kind(Kind::Mul(coeff.clone(), vec![(b, e)]))
        }
    }
}

/// SymPy's `_keep_coeff(coeff, factors, clear, sign)`.
pub fn keep_coeff(coeff: &Expr, factors: &Expr, clear: bool, sign: bool) -> Expr {
    let (coeff, factors) = if !is_number_atom(coeff) {
        if is_number_atom(factors) { (factors, coeff) } else { return ops::mul(&[coeff.clone(), factors.clone()]); }
    } else {
        (coeff, factors)
    };
    if factors.is_one() {
        return coeff.clone();
    }
    if coeff.is_one() {
        return factors.clone();
    }
    if coeff.is_minus_one() && !sign {
        return ops::neg(factors);
    }
    let Some(c) = coeff.as_number() else {
        // Infinite coefficients: ordinary evaluation keeps them as factors.
        return ops::mul(&[coeff.clone(), factors.clone()]);
    };
    if factors.is_add() {
        if !clear && c.is_rational() && !c.is_integer() {
            let parts: Vec<(Expr, Expr)> = args(factors)
                .iter()
                .map(|t| {
                    let (k, m) = ops::as_coeff_mul(t);
                    (keep_coeff(&Expr::num(k), coeff, true, false), m)
                })
                .collect();
            if parts.iter().any(|(k, _)| k.is_integer_number()) {
                let terms: Vec<Expr> = parts
                    .iter()
                    .map(|(k, m)| if k.is_one() { m.clone() } else { mul_from_args(&[k.clone(), m.clone()]) })
                    .collect();
                return add_from_args(&terms);
            }
        }
        return Expr::from_kind(Kind::Mul(c.clone(), vec![(factors.clone(), Expr::one())]));
    }
    if factors.is_mul() {
        let mut margs = args(factors);
        if let Some(n) = margs[0].as_number() {
            let m = n.mul(c);
            if m.is_one() {
                margs.remove(0);
            } else {
                margs[0] = Expr::num(m);
            }
        } else {
            margs.insert(0, coeff.clone());
        }
        return mul_from_args(&margs);
    }
    let m = ops::mul(&[coeff.clone(), factors.clone()]);
    if m.is_number() && !factors.is_number() && is_number_atom(&m) {
        return mul_from_args(&[coeff.clone(), factors.clone()]);
    }
    m
}

/// SymPy's `_unevaluated_Mul` (commutative factors).
pub fn unevaluated_mul(values: &[Expr]) -> Expr {
    let mut stack: Vec<Expr> = values.to_vec();
    let mut cargs: Vec<Expr> = Vec::new();
    let mut co = Number::ONE;
    let mut i = 0;
    while i < stack.len() {
        let a = stack[i].clone();
        i += 1;
        if a.is_mul() {
            stack.extend(args(&a));
        } else if let Some(n) = a.as_number() {
            co = co.mul(n);
        } else {
            cargs.push(a);
        }
    }
    cargs.sort_by(|a, b| ops::cmp_factor(&factor_pair(a), &factor_pair(b)));
    if !co.is_one() {
        cargs.insert(0, Expr::num(co));
    }
    mul_from_args(&cargs)
}

/// SymPy's `_unevaluated_Add`.
pub fn unevaluated_add(values: &[Expr]) -> Expr {
    let mut stack: Vec<Expr> = values.to_vec();
    let mut newargs = Vec::new();
    let mut co = Number::ZERO;
    while let Some(a) = stack.pop() {
        if a.is_add() {
            stack.extend(args(&a));
        } else if let Some(n) = a.as_number() {
            co = co.add(n);
        } else {
            newargs.push(a);
        }
    }
    if !co.is_zero() {
        newargs.insert(0, Expr::num(co));
    }
    add_from_args(&newargs)
}

/// `as_coeff_Mul(rational)` with the public argument view.
pub fn as_coeff_mul(e: &Expr, rational: bool) -> (Expr, Expr) {
    decompose::as_coeff_mul_head(e, rational)
}

/// `as_coeff_Add(rational)` with the public argument view.
pub fn as_coeff_add(e: &Expr, rational: bool) -> (Expr, Expr) {
    decompose::as_coeff_add_head(e, rational)
}

/// `args_cnc()` commutative part as a set (split_1=True).
fn cset(e: &Expr) -> Vec<Expr> {
    let (mut c, _) = decompose::args_cnc(e, true);
    c.sort();
    c.dedup();
    c
}

/// SymPy's `Expr.coeff(x, n)` restricted to commutative expressions.
pub fn coeff(e: &Expr, x: &Expr, n: i64) -> Expr {
    if x.is_zero() {
        return Expr::zero();
    }
    if x == e {
        return if n == 1 { Expr::one() } else { Expr::zero() };
    }
    if x.is_one() {
        let co: Vec<Expr> = add_args(e).into_iter().filter(|a| as_coeff_mul(a, false).0.is_one()).collect();
        return ops::add(&co);
    }
    if n == 0 {
        if x.is_add() && e.is_add() {
            let c = coeff(e, x, 1);
            if c.is_zero() {
                return Expr::zero();
            }
            let parts: Vec<Expr> = add_args(&c).iter().map(|a| ops::mul(&[a.clone(), x.clone()])).collect();
            return ops::sub(e, &ops::add(&parts));
        }
        return decompose::as_independent(e, std::slice::from_ref(x), Some(true)).0;
    }
    let x = ops::powi(x, n);
    let xargs = cset(&x);
    let mut co2 = Vec::new();
    for a in add_args(e) {
        let margs = cset(&a);
        if xargs.len() > margs.len() {
            continue;
        }
        let resid: Vec<Expr> = margs.iter().filter(|m| !xargs.contains(m)).cloned().collect();
        if resid.len() + xargs.len() == margs.len() {
            co2.push(ops::mul(&resid));
        }
    }
    ops::add(&co2)
}

/// `Add.primitive()`.
pub fn add_primitive(e: &Expr) -> (Number, Expr) {
    let values = args(e);
    let mut terms = Vec::new();
    let mut inf = false;
    for a in &values {
        let (c, m) = ops::as_coeff_mul(a);
        let (c, m) = if c.is_rational() { (c, m) } else { (Number::ONE, a.clone()) };
        inf |= m.is_constant(Constant::ComplexInfinity);
        terms.push((c, m));
    }
    let mut ngcd = Number::ZERO;
    let mut dlcm = num_bigint::BigInt::from(1);
    for (c, _) in &terms {
        let (p, q) = c.numer_denom().unwrap();
        ngcd = ngcd.gcd(&Number::from_bigint(p));
        dlcm = num_integer::Integer::lcm(&dlcm, &q);
    }
    let _ = inf;
    let dlcm_n = Number::from_bigint(dlcm);
    if ngcd.is_one() && dlcm_n.is_one() {
        return (Number::ONE, e.clone());
    }
    let ngcd = ngcd.abs();
    let content = ngcd.div(&dlcm_n).unwrap();
    let parts: Vec<Expr> = terms
        .iter()
        .map(|(c, m)| keep_coeff(&Expr::num(c.div(&content).unwrap()), m, true, false))
        .collect();
    (content, add_from_args(&parts))
}

/// SymPy's `Expr.primitive()`.
pub fn primitive(e: &Expr) -> (Number, Expr) {
    if e.is_zero() {
        return (Number::ONE, Expr::zero());
    }
    if e.is_add() {
        return add_primitive(e);
    }
    let (c, r) = as_coeff_mul(e, true);
    let c = c.as_number().cloned().unwrap_or(Number::ONE);
    if c.is_negative() {
        return (c.neg(), ops::neg(&r));
    }
    (c, r)
}

/// `Mul.as_two_terms()`.
fn as_two_terms(e: &Expr) -> (Expr, Expr) {
    let a = args(e);
    match a.len() {
        0 | 1 => (Expr::one(), e.clone()),
        2 => (a[0].clone(), a[1].clone()),
        _ => (a[0].clone(), mul_from_args(&a[1..])),
    }
}

fn is_integer_value(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(n) if n.is_integer())
}

/// SymPy's `Expr.extract_multiplicatively(c)`.
pub fn extract_multiplicatively(e: &Expr, c: &Expr) -> Option<Expr> {
    if e.is_constant(Constant::NaN) {
        return None;
    }
    if c.is_one() {
        return Some(e.clone());
    }
    if c == e {
        return Some(Expr::one());
    }
    let mut c = c.clone();
    if c.is_add() {
        let (cc, pc) = add_primitive(&c);
        if !cc.is_one() {
            c = held_coeff_mul(&cc, &pc);
        }
    }
    if c.is_mul() {
        let (a, b) = as_two_terms(&c);
        let x = extract_multiplicatively(e, &a)?;
        return extract_multiplicatively(&x, &b);
    }
    let quotient = ops::div(e, &c);
    match e.kind() {
        Kind::Const(Constant::Infinity) => {
            if assumptions::is_positive(&c) == Some(true) {
                return Some(Expr::oo());
            }
        }
        Kind::Const(Constant::NegInfinity) => {
            if assumptions::is_negative(&c) == Some(true) {
                return Some(Expr::oo());
            }
            if assumptions::is_positive(&c) == Some(true) {
                return Some(Expr::neg_oo());
            }
        }
        Kind::Const(Constant::ComplexInfinity) => {
            if assumptions::is_zero(&c) != Some(true) {
                return Some(Expr::zoo());
            }
        }
        Kind::Num(n) => {
            let q = quotient.as_number()?;
            let ok = if n.is_integer() {
                q.is_integer()
            } else if n.is_rational() {
                q.is_rational()
            } else {
                q.is_float()
            };
            if !ok || (n.is_positive() && q.is_negative()) {
                return None;
            }
            return Some(quotient);
        }
        Kind::Const(_) | Kind::Sym(_) => {
            if quotient.is_mul() {
                let qa = args(&quotient);
                if qa.len() == 2 && is_integer_value(&qa[0]) && qa[0].as_number().is_some_and(Number::is_positive) && qa[1] == *e {
                    return Some(quotient);
                }
            } else if is_integer_value(&quotient) && is_number_atom(&c) {
                return Some(quotient);
            }
        }
        Kind::Add(..) => {
            let (cs, ps) = add_primitive(e);
            if is_number_atom(&c) && !c.is_minus_one() {
                if !cs.is_one() {
                    let cs_e = Expr::num(cs.clone());
                    let xc = if assumptions::is_negative(&c) == Some(true) {
                        extract_multiplicatively(&cs_e, &ops::neg(&c)).map(|v| ops::neg(&v))
                    } else {
                        extract_multiplicatively(&cs_e, &c)
                    };
                    if let Some(xc) = xc {
                        return Some(ops::mul(&[xc, ps]));
                    }
                }
                return None;
            }
            if c == ps {
                return Some(Expr::num(cs));
            }
            let mut newargs = Vec::new();
            for arg in args(&ps) {
                newargs.push(extract_multiplicatively(&arg, &c)?);
            }
            if !cs.is_one() {
                let scaled: Vec<Expr> = newargs.iter().map(|t| ops::mul(&[Expr::num(cs.clone()), t.clone()])).collect();
                return Some(unevaluated_add(&scaled));
            }
            return Some(add_from_args(&newargs));
        }
        Kind::Mul(..) => {
            let mut margs = args(e);
            for i in 0..margs.len() {
                if let Some(newarg) = extract_multiplicatively(&margs[i], &c) {
                    margs[i] = newarg;
                    return Some(ops::mul(&margs));
                }
            }
        }
        Kind::Pow(..) | Kind::Func(crate::expr::Func::Exp, _) => {
            let (sb, se) = decompose::as_base_exp(e);
            let (cb, ce) = decompose::as_base_exp(&c);
            if cb == sb {
                if let Some(ne) = extract_additively(&se, &ce) {
                    return Some(ops::pow(&sb, &ne));
                }
            } else if c == sb {
                if let Some(ne) = extract_additively(&se, &Expr::one()) {
                    return Some(ops::pow(&sb, &ne));
                }
            }
        }
        _ => {}
    }
    None
}

/// `_corem(eq, c)`: `(co, diff)` with `eq == co*c + diff`.
fn corem(e: &Expr, c: &Expr) -> (Expr, Expr) {
    let mut co = Vec::new();
    let mut non = Vec::new();
    for i in add_args(e) {
        let ci = coeff(&i, c, 1);
        if ci.is_zero() { non.push(i) } else { co.push(ci) }
    }
    (ops::add(&co), ops::add(&non))
}

fn number_lt(a: &Number, b: &Number) -> bool {
    a.cmp_value(b) == Some(std::cmp::Ordering::Less)
}

fn or_zero(v: Option<Expr>) -> Expr {
    v.unwrap_or_else(Expr::zero)
}

fn nonzero(v: Expr) -> Option<Expr> {
    if v.is_zero() { None } else { Some(v) }
}

/// SymPy's `Expr.extract_additively(c)`.
pub fn extract_additively(e: &Expr, c: &Expr) -> Option<Expr> {
    if e.is_constant(Constant::NaN) {
        return None;
    }
    if assumptions::is_zero(c) == Some(true) {
        return Some(e.clone());
    }
    if c == e {
        return Some(Expr::zero());
    }
    if e.is_zero() {
        return None;
    }
    if is_number_atom(e) {
        if !is_number_atom(c) {
            return None;
        }
        let (co, cn) = (e.as_number()?, c.as_number()?);
        let diff = co.sub(cn);
        let zero = Number::ZERO;
        let ok = (number_lt(&zero, co) && !number_lt(&diff, &zero) && number_lt(&diff, co))
            || (number_lt(co, &zero) && !number_lt(&zero, &diff) && number_lt(co, &diff));
        return if ok { Some(Expr::num(diff)) } else { None };
    }
    if is_number_atom(c) {
        let (co, t) = as_coeff_add(e, false);
        let xa = extract_additively(&co, c)?;
        return Some(ops::add(&[xa, t]));
    }
    if c.is_add() && is_number_atom(&args(c)[0]) {
        let co = coeff(e, c, 1);
        let xa0 = ops::mul(&[or_zero(extract_additively(&co, &Expr::one())), c.clone()]);
        if !xa0.is_zero() {
            let diff = ops::sub(e, &ops::mul(&[co, c.clone()]));
            let rest = extract_additively(&diff, c).unwrap_or(diff);
            return nonzero(ops::add(&[xa0, rest]));
        }
        let (h, t) = as_coeff_add(c, false);
        let (sh, st) = as_coeff_add(e, false);
        let xa = extract_additively(&sh, &h)?;
        let xa2 = extract_additively(&st, &t)?;
        return Some(ops::add(&[xa, xa2]));
    }
    let (co, diff) = corem(e, c);
    let xa0 = ops::mul(&[or_zero(extract_additively(&co, &Expr::one())), c.clone()]);
    if !xa0.is_zero() {
        let rest = extract_additively(&diff, c).unwrap_or(diff);
        return nonzero(ops::add(&[xa0, rest]));
    }
    let mut cur = e.clone();
    let mut coeffs = Vec::new();
    for a in add_args(c) {
        let (ac, at) = as_coeff_mul(&a, false);
        let co = coeff(&cur, &at, 1);
        if co.is_zero() {
            return None;
        }
        let (coc, cot) = as_coeff_add(&co, false);
        let xa = extract_additively(&coc, &ac)?;
        cur = ops::sub(&cur, &ops::mul(&[co, at.clone()]));
        coeffs.push(ops::mul(&[ops::add(&[cot, xa]), at]));
    }
    coeffs.push(cur);
    Some(ops::add(&coeffs))
}

/// SymPy's `Expr.as_coefficient(expr)`.
pub fn as_coefficient(e: &Expr, x: &Expr) -> Option<Expr> {
    let r = extract_multiplicatively(e, x)?;
    if !r.is_zero() && !traverse::has(&r, x) { Some(r) } else { None }
}

/// SymPy's `Expr.normal()`.
pub fn normal(e: &Expr) -> Expr {
    let (n, d) = crate::simptools::as_numer_denom(e);
    if d.is_one() {
        return n;
    }
    if let Some(dn) = d.as_number() {
        return unevaluated_mul(&[n, Expr::num(Number::ONE.div(dn).unwrap_or(Number::ZERO))]);
    }
    ops::div(&n, &d)
}

/// `expand(e, deep, mul=..., multinomial=..., ...)` with only the given
/// hints enabled; resource failures leave the expression unchanged.
pub fn expand_hints(e: &Expr, mul: bool, multinomial: bool, power_base: bool, power_exp: bool, log: bool, basic: bool) -> Expr {
    use crate::expand::{expand_with_options, ExpandBudget, ExpandOptions};
    let options = ExpandOptions { deep: true, basic, mul, multinomial, power_base, power_exp, log, trig: false, force: false };
    expand_with_options(e, &options, &mut ExpandBudget::default()).unwrap_or_else(|_| e.clone())
}

/// `expand_mul(e, deep=False)`.
pub fn expand_hints_shallow_mul(e: &Expr) -> Expr {
    use crate::expand::{expand_with_options, ExpandBudget, ExpandOptions};
    let options = ExpandOptions { deep: false, basic: false, mul: true, multinomial: false, power_base: false, power_exp: false, log: false, trig: false, force: false };
    expand_with_options(e, &options, &mut ExpandBudget::default()).unwrap_or_else(|_| e.clone())
}

/// `expand_mul(e)`.
pub fn expand_mul(e: &Expr) -> Expr {
    expand_hints(e, true, false, false, false, false, false)
}

/// `expand_multinomial(e)`.
pub fn expand_multinomial(e: &Expr) -> Expr {
    expand_hints(e, false, true, false, false, false, false)
}

/// `_mexpand(e)`: `expand_mul(expand_multinomial(e))`.
pub fn mexpand(e: &Expr) -> Expr {
    expand_mul(&expand_multinomial(e))
}

/// `expr.func(*[f(a) for a in expr.args])` over SymPy's public arguments,
/// rebuilt with the evaluating constructors. Unchanged nodes are returned
/// as they are.
pub fn map_args(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => e.clone(),
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Func(..) | Kind::Rel(..) => {
            if let Kind::Func(fun, _) = e.kind() {
                if matches!(fun, Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(..) | Func::Subs | Func::Integral | Func::Sum | Func::Limit) {
                    return traverse::map_children(e, f);
                }
            }
            let old = args(e);
            let new: Vec<Expr> = old.iter().map(|a| f(a)).collect();
            if new.iter().zip(&old).all(|(a, b)| a.same_presentation(b)) {
                return e.clone();
            }
            rebuild(e, &new)
        }
        _ => traverse::map_children(e, f),
    }
}

/// `expr.func(*args)` with evaluation.
pub fn rebuild(e: &Expr, new: &[Expr]) -> Expr {
    match e.kind() {
        Kind::Add(..) => ops::add(new),
        Kind::Mul(..) => ops::mul(new),
        Kind::Pow(..) => ops::pow(&new[0], &new[1]),
        Kind::Func(fun, _) => crate::functions::func(fun.clone(), new.to_vec()),
        Kind::Rel(op, ..) => crate::relational::rel(*op, &new[0], &new[1]),
        _ => e.clone(),
    }
}

/// Pre-order traversal over public arguments.
pub fn preorder(e: &Expr, f: &mut dyn FnMut(&Expr)) {
    f(e);
    if matches!(e.kind(), Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Func(..) | Kind::Rel(..)) {
        for a in args(e) {
            preorder(&a, f);
        }
    }
}

/// Distinct subexpressions satisfying `pred` (`expr.atoms(Type)`).
pub fn atoms_where(e: &Expr, pred: &dyn Fn(&Expr) -> bool) -> Vec<Expr> {
    let mut out: Vec<Expr> = Vec::new();
    preorder(e, &mut |a| {
        if pred(a) && !out.contains(a) {
            out.push(a.clone());
        }
    });
    out
}

/// `_nodes(e)`: the number of nodes of the public tree.
pub fn nodes(e: &Expr) -> usize {
    let mut n = 0;
    preorder(e, &mut |_| n += 1);
    n
}

/// SymPy's `ordered(seq)`: by node count, then `default_sort_key`.
pub fn ordered(mut v: Vec<Expr>) -> Vec<Expr> {
    v.sort_by(|a, b| nodes(a).cmp(&nodes(b)).then_with(|| crate::sortkey::cmp_default(a, b)));
    v
}

/// `bottom_up(rv, F)`: rebuild public arguments first, then apply `f`.
pub fn bottom_up(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
    let rebuilt = map_args(e, &mut |c| bottom_up(c, f));
    f(&rebuilt)
}

/// Sign of a real numeric expression, decided like SymPy's `e < 0`
/// (assumptions first, then a high precision evaluation).
pub(crate) fn num_sign(e: &Expr) -> Option<std::cmp::Ordering> {
    use std::cmp::Ordering;
    if let Some(n) = e.as_number() {
        return n.cmp_value(&Number::ZERO);
    }
    if assumptions::is_zero(e) == Some(true) {
        return Some(Ordering::Equal);
    }
    if assumptions::is_positive(e) == Some(true) {
        return Some(Ordering::Greater);
    }
    if assumptions::is_negative(e) == Some(true) {
        return Some(Ordering::Less);
    }
    if !is_num(e) {
        return None;
    }
    let v = crate::evalf_mp::evalf(e, 200);
    let f = v.as_number()?.to_f64();
    if f.abs() > 1e-40 {
        return Some(if f > 0.0 { Ordering::Greater } else { Ordering::Less });
    }
    if methods::equals(e, &Expr::zero()) == Some(true) {
        return Some(Ordering::Equal);
    }
    (f != 0.0).then(|| if f > 0.0 { Ordering::Greater } else { Ordering::Less })
}

/// `a < b` for real numbers (false when undecided).
pub(crate) fn num_lt(a: &Expr, b: &Expr) -> bool {
    num_sign(&ops::sub(a, b)) == Some(std::cmp::Ordering::Less)
}

/// `a <= b` for real numbers (false when undecided).
pub(crate) fn num_le(a: &Expr, b: &Expr) -> bool {
    matches!(num_sign(&ops::sub(a, b)), Some(std::cmp::Ordering::Less | std::cmp::Ordering::Equal))
}

/// SymPy's `expr.as_numer_denom()` (Add, Mul, Pow, exp and Rational).
pub fn as_numer_denom(e: &Expr) -> (Expr, Expr) {
    let is_exp = |b: &Expr| b.is_constant(Constant::E);
    match e.kind() {
        Kind::Num(n) if n.is_rational() => (Expr::num(n.numer()), Expr::num(n.denom())),
        Kind::Mul(..) => {
            let (n, d): (Vec<Expr>, Vec<Expr>) = args(e).iter().map(as_numer_denom).unzip();
            (ops::mul(&n), ops::mul(&d))
        }
        Kind::Func(Func::Exp, a) => exp_numer_denom(e, &a[0]),
        Kind::Pow(b, x) if is_exp(b) => exp_numer_denom(e, x),
        Kind::Pow(base, x) => {
            let (mut n, mut d) = as_numer_denom(base);
            let mut neg_exp = assumptions::is_negative(x) == Some(true);
            if x.is_mul() && !neg_exp && assumptions::is_positive(x) != Some(true) {
                neg_exp = decompose::could_extract_minus_sign(x);
            }
            let int_exp = assumptions::is_integer(x) == Some(true);
            if !(assumptions::is_real(&d) == Some(true) || int_exp) {
                n = base.clone();
                d = Expr::one();
            }
            let dnonpos = assumptions::is_nonpositive(&d);
            if dnonpos == Some(true) {
                n = ops::neg(&n);
                d = ops::neg(&d);
            } else if dnonpos.is_none() && !int_exp {
                n = base.clone();
                d = Expr::one();
            }
            let mut x = x.clone();
            if neg_exp {
                std::mem::swap(&mut n, &mut d);
                x = ops::neg(&x);
            }
            if assumptions::is_finite(&x) == Some(false) {
                if n.is_one() && !d.is_one() {
                    return (n, ops::pow(&d, &x));
                }
                if !n.is_one() && d.is_one() {
                    return (ops::pow(&n, &x), d);
                }
                return (e.clone(), Expr::one());
            }
            (ops::pow(&n, &x), ops::pow(&d, &x))
        }
        Kind::Add(..) => {
            let (content, expr) = primitive(e);
            if !expr.is_add() {
                return as_numer_denom(&held_coeff_mul(&content, &expr));
            }
            let (ncon, dcon) = (Expr::num(content.numer()), Expr::num(content.denom()));
            let mut nd: Vec<(Expr, Vec<Expr>)> = Vec::new();
            for f in args(&expr) {
                let (ni, di) = as_numer_denom(&f);
                match nd.iter_mut().find(|(d, _)| *d == di) {
                    Some(slot) => slot.1.push(ni),
                    None => nd.push((di, vec![ni])),
                }
            }
            if nd.len() == 1 {
                let (d, n) = nd.pop().unwrap();
                let terms: Vec<Expr> = n.iter().map(|ni| keep_coeff(&ncon, ni, false, true)).collect();
                return (ops::add(&terms), keep_coeff(&dcon, &d, false, true));
            }
            let denoms: Vec<Expr> = nd.iter().map(|(d, _)| d.clone()).collect();
            let numers: Vec<Expr> = nd.iter().map(|(_, n)| if n.len() == 1 { n[0].clone() } else { ops::add(n) }).collect();
            let mut terms = Vec::new();
            for i in 0..numers.len() {
                let mut f: Vec<Expr> = denoms[..i].to_vec();
                f.push(numers[i].clone());
                f.extend_from_slice(&denoms[i + 1..]);
                terms.push(ops::mul(&f));
            }
            (keep_coeff(&ncon, &ops::add(&terms), false, true), keep_coeff(&dcon, &ops::mul(&denoms), false, true))
        }
        _ => (e.clone(), Expr::one()),
    }
}

fn exp_numer_denom(e: &Expr, x: &Expr) -> (Expr, Expr) {
    let mut neg_exp = assumptions::is_negative(x) == Some(true);
    if !neg_exp && assumptions::is_negative(&ops::neg(x)) != Some(true) {
        neg_exp = decompose::could_extract_minus_sign(x);
    }
    if neg_exp {
        return (Expr::one(), crate::functions::exp(&ops::neg(x)));
    }
    (e.clone(), Expr::one())
}
