//! `sympy.simplify.radsimp`: `fraction`, `numer`, `denom`, `collect_const`,
//! `collect_sqrt`, `rcollect`, `radsimp`, `rad_rationalize` and `split_surds`.
use super::exprtools::{gcd_terms, Factors};
use super::{is_num, add_args, add_from_args, args, is_number_atom, keep_coeff, mexpand, mul_args, nodes, ordered, unevaluated_add, unevaluated_mul};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::poly::int::Rat;
use crate::sortkey::cmp_default;
use crate::{assumptions, decompose, ops, traverse};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;

/// An exact decision for rational linear combinations of positive rational
/// square roots. None means that this proof does not apply.
#[derive(Debug, PartialEq, Eq)]
pub enum SquareRootSum {
    Rational(Rat),
    Irrational,
}

fn rational_square_root(r: &Rat) -> Option<Rat> {
    if r.is_negative() { return None; }
    let n = r.n.isqrt();
    let d = r.d.isqrt();
    (n.mul(&n) == r.n && d.mul(&d) == r.d).then(|| Rat::new(n, d))
}

/// Distinct classes of positive rationals modulo rational squares have
/// linearly independent square roots over QQ. Comparing square classes
/// requires only exact integer square roots of a reduced rational ratio,
/// not factorization, numerical recognition or a minimal polynomial.
pub fn square_root_sum(e: &Expr) -> Option<SquareRootSum> {
    let mut rational = Rat::zero();
    let mut classes: Vec<(Rat, Rat)> = Vec::new();
    for term in add_args(e) {
        let (mut coefficient, mut radicand) = (Rat::one(), Rat::one());
        for factor in mul_args(&term) {
            if let Some(r) = factor.as_number().and_then(Rat::from_number) {
                coefficient = coefficient.mul(&r);
            } else if let Kind::Pow(base, exponent) = factor.kind() {
                if *exponent != Expr::rational(1, 2) { return None; }
                let r = base.as_number().and_then(Rat::from_number)?;
                if r.is_negative() { return None; }
                radicand = radicand.mul(&r);
            } else {
                return None;
            }
        }
        if coefficient.is_zero() || radicand.is_zero() { continue; }
        if let Some(root) = rational_square_root(&radicand) {
            rational = rational.add(&coefficient.mul(&root));
            continue;
        }
        let mut matched = false;
        for (representative, sum) in &mut classes {
            if let Some(root) = rational_square_root(&radicand.div(representative)) {
                *sum = sum.add(&coefficient.mul(&root));
                matched = true;
                break;
            }
        }
        if !matched { classes.push((radicand, coefficient)); }
    }
    Some(if classes.iter().all(|(_, c)| c.is_zero()) {
        SquareRootSum::Rational(rational)
    } else {
        SquareRootSum::Irrational
    })
}

/// The numeric-denominator preparation of radsimp. Unrecognized exact
/// forms keep the existing host recognition callback; an error or signal
/// stops the traversal before any further callback or replacement.
pub fn prepare_numeric_denominators<E>(
    e: &Expr,
    check: &mut impl FnMut() -> Result<(), E>,
    fallback: &mut impl FnMut(&Expr) -> Result<Option<Expr>, E>,
) -> Result<Expr, E> {
    let mut pending = vec![e.clone()];
    let mut visited = FxHashSet::default();
    let mut replacements = FxHashMap::default();
    while let Some(node) = pending.pop() {
        check()?;
        if !visited.insert(node.clone()) { continue; }
        pending.extend(node.args());
        let Kind::Pow(base, exponent) = node.kind() else { continue };
        if !exponent.is_integer_number() || !exponent.as_number().is_some_and(Number::is_negative)
            || !base.is_add() || base.float_info().max_precision.is_some() || !is_num(base) {
            continue;
        }
        let expanded = mexpand(base);
        let candidate = match square_root_sum(&expanded) {
            Some(SquareRootSum::Rational(r)) => Some(Expr::num(r.to_number())),
            Some(SquareRootSum::Irrational) => None,
            None => fallback(&expanded)?,
        };
        if let Some(candidate) = candidate {
            if candidate != *base {
                replacements.insert(node.clone(), ops::pow(&candidate, exponent));
            }
        }
    }
    check()?;
    Ok(crate::subs::xreplace(e, &replacements))
}

/// `fraction(expr, exact)`: `(numerator, denominator)` without
/// rewriting the expression.
pub fn fraction(e: &Expr, exact: bool) -> (Expr, Expr) {
    let mut numer = Vec::new();
    let mut denom = Vec::new();
    for term in mul_args(e) {
        let is_powlike = term.is_pow() || matches!(term.kind(), Kind::Func(Func::Exp, _));
        if is_powlike && decompose::commutativity(&term) == Some(true) {
            let (b, ex) = decompose::as_base_exp(&term);
            if assumptions::is_negative(&ex) == Some(true) {
                if ex.is_minus_one() {
                    denom.push(b);
                } else if exact {
                    if super::methods::is_constant(&ex, &[]) == Some(true) {
                        denom.push(ops::pow(&b, &ops::neg(&ex)));
                    } else {
                        numer.push(term);
                    }
                } else {
                    denom.push(ops::pow(&b, &ops::neg(&ex)));
                }
            } else if assumptions::is_positive(&ex) == Some(true) {
                numer.push(term);
            } else if !exact && ex.is_mul() {
                let (n, d) = crate::simptools::as_numer_denom(&term);
                if !n.is_one() {
                    numer.push(n);
                }
                denom.push(d);
            } else {
                numer.push(term);
            }
        } else if let Some(r) = term.as_number().filter(|n| n.is_rational() && !n.is_integer()) {
            let (p, q) = r.numer_denom().unwrap();
            if p != 1.into() {
                numer.push(Expr::num(Number::from_bigint(p)));
            }
            denom.push(Expr::num(Number::from_bigint(q)));
        } else {
            numer.push(term);
        }
    }
    if exact {
        (held_mul(&numer), held_mul(&denom))
    } else {
        (ops::mul(&numer), ops::mul(&denom))
    }
}

/// `Mul(*args, evaluate=False)`.
pub fn held_mul(values: &[Expr]) -> Expr {
    match values.len() {
        0 => Expr::one(),
        1 => values[0].clone(),
        _ => {
            let canonical = ops::mul(values);
            if args(&canonical) == values {
                return canonical;
            }
            let (coeff, rest) = match values[0].as_number() {
                Some(n) if !n.is_one() => (n.clone(), &values[1..]),
                _ => (Number::ONE, values),
            };
            Expr::from_unevaluated_kind(Kind::Mul(coeff, rest.iter().map(|a| (a.clone(), Expr::one())).collect()))
        }
    }
}

pub fn numer(e: &Expr, exact: bool) -> Expr {
    fraction(e, exact).0
}

pub fn denom(e: &Expr, exact: bool) -> Expr {
    fraction(e, exact).1
}

/// `collect_const(expr, *vars, Numbers=numbers)`.
pub fn collect_const(e: &Expr, vars: &[Expr], numbers: bool) -> Expr {
    if !e.is_add() {
        return e.clone();
    }
    let mut expr = e.clone();
    let recurse = vars.is_empty();
    let mut vs: Vec<Expr> = if recurse {
        let mut found: Vec<Expr> = Vec::new();
        for a in args(e) {
            for m in mul_args(&a) {
                if is_num(&m) && !found.contains(&m) {
                    found.push(m);
                }
            }
        }
        found
    } else {
        let mut found: Vec<Expr> = Vec::new();
        for v in vars {
            if !found.contains(v) {
                found.push(v.clone());
            }
        }
        found
    };
    if !numbers {
        vs.retain(|v| !is_number_atom(v));
    }
    let mut vs = ordered(vs);
    let mut idx = 0;
    while idx < vs.len() {
        let v = vs[idx].clone();
        idx += 1;
        let fv = Factors::from_expr(&v);
        let mut collected: Vec<Expr> = Vec::new();
        let mut others: Vec<Expr> = Vec::new();
        for m in add_args(&expr) {
            let f = Factors::from_expr(&m);
            let (q, r) = f.div(&fv);
            if r.is_one() {
                let lost_integer = q.factors.iter().any(|(k, now)| {
                    f.factors.iter().any(|(kw, was)| kw == k && was.is_integer_number() && !now.is_integer_number())
                });
                if !lost_integer {
                    collected.push(q.as_expr());
                    continue;
                }
            }
            others.push(m);
        }
        // `ordered(terms)` on the keys {v, 1}.
        let mut keys: Vec<(Expr, Vec<Expr>)> = Vec::new();
        if !collected.is_empty() {
            keys.push((v.clone(), collected));
        }
        if !others.is_empty() {
            keys.push((Expr::one(), others));
        }
        keys.sort_by(|a, b| nodes(&a.0).cmp(&nodes(&b.0)).then_with(|| cmp_default(&a.0, &b.0)));
        let mut out = Vec::new();
        let mut hit = false;
        let mut uneval = false;
        for (k, terms) in keys {
            if k.is_one() {
                out.extend(terms);
                continue;
            }
            let val = if terms.len() > 1 {
                let s = ops::add(&terms);
                hit = true;
                if recurse && s != expr {
                    vs.push(s.clone());
                }
                s
            } else {
                terms[0].clone()
            };
            if numbers && is_number_atom(&k) && val.is_add() {
                out.push(keep_coeff(&k, &val, true, true));
                uneval = true;
            } else {
                out.push(ops::mul(&[k.clone(), val]));
            }
        }
        if hit {
            expr = if uneval { unevaluated_add(&out) } else { ops::add(&out) };
            if !expr.is_add() {
                break;
            }
        }
    }
    expr
}

fn is_sqrt_like(m: &Expr) -> bool {
    if m.is_constant(Constant::I) {
        return true;
    }
    if let Kind::Pow(_, x) = m.kind() {
        if let Some(n) = x.as_number() {
            if let Some((_, q)) = n.numer_denom() {
                return q == 2.into();
            }
        }
    }
    false
}

/// `collect_sqrt(expr)` (evaluate=True).
pub fn collect_sqrt(e: &Expr) -> Expr {
    let Ok((coeff, expr)) = decompose::as_content_primitive(e, false, true) else { return e.clone() };
    let mut vars: Vec<Expr> = Vec::new();
    for a in add_args(&expr) {
        for m in decompose::args_cnc(&a, true).0 {
            if is_num(&m) && is_sqrt_like(&m) && !vars.contains(&m) {
                vars.push(m);
            }
        }
    }
    let d = collect_const(&expr, &vars, false);
    ops::mul(&[Expr::num(coeff), d])
}

/// `collect_sqrt(expr, evaluate=False)`: `(terms, number of radical terms)`.
pub fn collect_sqrt_terms(e: &Expr) -> (Vec<Expr>, usize) {
    let Ok((coeff, expr)) = decompose::as_content_primitive(e, false, true) else { return (vec![e.clone()], 0) };
    let mut vars: Vec<Expr> = Vec::new();
    for a in add_args(&expr) {
        for m in decompose::args_cnc(&a, true).0 {
            if is_num(&m) && is_sqrt_like(&m) && !vars.contains(&m) {
                vars.push(m);
            }
        }
    }
    let d = collect_const(&expr, &vars, false);
    let hit = expr != d;
    let mut nrad = 0;
    let mut out = ordered(add_args(&d));
    for m in out.iter_mut() {
        let (c, _) = decompose::args_cnc(m, true);
        if c.iter().any(is_sqrt_like) {
            nrad += 1;
        }
        *m = ops::mul(&[m.clone(), Expr::num(coeff.clone())]);
    }
    if !(hit || nrad > 0) {
        out = vec![ops::add(&out)];
    }
    (out, nrad)
}

/// `rcollect(expr, *vars)`.
pub fn rcollect(e: &Expr, vars: &[Expr]) -> Expr {
    if e.is_atom() || !vars.iter().any(|v| traverse::has(e, v)) {
        return e.clone();
    }
    let rebuilt = super::map_args(e, &mut |a| rcollect(a, vars));
    if rebuilt.is_add() {
        crate::poly::api::collect(&rebuilt, vars)
    } else {
        rebuilt
    }
}

fn ispow2(d: &Expr, symbolic: bool, log2: bool) -> bool {
    let Kind::Pow(_, e) = d.kind() else { return false };
    let rational_q = e.as_number().and_then(Number::numer_denom).map(|(_, q)| q);
    if rational_q.as_ref().is_some_and(|q| *q == 2.into()) {
        return true;
    }
    if symbolic && rational_q.is_none() && denom(e, false) == Expr::int(2) {
        return true;
    }
    if log2 {
        let mut q = num_bigint::BigInt::from(1);
        if let Some(rq) = rational_q {
            q = rq;
        } else if symbolic {
            if let Some(n) = denom(e, false).as_number().and_then(Number::to_bigint) {
                q = n;
            }
        }
        if q != 1.into() && q > 0.into() && (&q & (&q - 1u32)) == 0.into() {
            return true;
        }
    }
    false
}

fn has_bad(d: &Expr) -> bool {
    let mut bad = false;
    super::preorder(d, &mut |a| {
        if a.is_zero() || a.is_constant(Constant::NaN) || a.is_constant(Constant::ComplexInfinity) {
            bad = true;
        }
    });
    bad
}

/// The conjugate multiplier of `radsimp._num`.
fn conj_num(rterms: &[(Expr, Expr)]) -> Expr {
    let s = |i: usize| ops::sqrt(&rterms[i].0);
    let (big, small): (Vec<Expr>, Vec<Expr>) = rterms.iter().cloned().unzip();
    let m = |v: &[Expr]| ops::mul(v);
    let i = |n: i64| Expr::int(n);
    match rterms.len() {
        1 => s(0),
        2 => ops::sub(&m(&[s(0), small[0].clone()]), &m(&[s(1), small[1].clone()])),
        3 => {
            let (a, b, c) = (&small[0], &small[1], &small[2]);
            let (aa, bb, cc) = (&big[0], &big[1], &big[2]);
            let f1 = ops::add(&[m(&[s(0), a.clone()]), m(&[s(1), b.clone()]), ops::neg(&m(&[s(2), c.clone()]))]);
            let f2 = ops::add(&[
                m(&[i(2), s(0), s(1), a.clone(), b.clone()]),
                ops::neg(&m(&[aa.clone(), ops::powi(a, 2)])),
                ops::neg(&m(&[bb.clone(), ops::powi(b, 2)])),
                m(&[cc.clone(), ops::powi(c, 2)]),
            ]);
            m(&[f1, f2])
        }
        _ => {
            let (a, b, c, d) = (&small[0], &small[1], &small[2], &small[3]);
            let (aa, bb, cc, dd) = (&big[0], &big[1], &big[2], &big[3]);
            let p2 = |x: &Expr| ops::powi(x, 2);
            let p4 = |x: &Expr| ops::powi(x, 4);
            let f1 = ops::add(&[m(&[s(0), a.clone()]), m(&[s(1), b.clone()]), ops::neg(&m(&[s(2), c.clone()])), ops::neg(&m(&[s(3), d.clone()]))]);
            let f2 = ops::add(&[
                m(&[i(2), s(0), s(1), a.clone(), b.clone()]),
                ops::neg(&m(&[aa.clone(), p2(a)])),
                ops::neg(&m(&[bb.clone(), p2(b)])),
                m(&[i(-2), s(2), s(3), c.clone(), d.clone()]),
                m(&[cc.clone(), p2(c)]),
                m(&[dd.clone(), p2(d)]),
            ]);
            let f3 = ops::add(&[
                m(&[i(-8), s(0), s(1), s(2), s(3), a.clone(), b.clone(), c.clone(), d.clone()]),
                m(&[p2(aa), p4(a)]),
                m(&[i(-2), aa.clone(), bb.clone(), p2(a), p2(b)]),
                m(&[i(-2), aa.clone(), cc.clone(), p2(a), p2(c)]),
                m(&[i(-2), aa.clone(), dd.clone(), p2(a), p2(d)]),
                m(&[p2(bb), p4(b)]),
                m(&[i(-2), bb.clone(), cc.clone(), p2(b), p2(c)]),
                m(&[i(-2), bb.clone(), dd.clone(), p2(b), p2(d)]),
                m(&[p2(cc), p4(c)]),
                m(&[i(-2), cc.clone(), dd.clone(), p2(c), p2(d)]),
                m(&[p2(dd), p4(d)]),
            ]);
            m(&[f1, f2, f3])
        }
    }
}

/// `ordered` applied to `(radicands, terms)` pairs.
fn cmp_rterm(a: &(Vec<Expr>, Vec<Expr>), b: &(Vec<Expr>, Vec<Expr>)) -> Ordering {
    let n = |x: &(Vec<Expr>, Vec<Expr>)| 3 + x.0.iter().map(nodes).sum::<usize>() + x.1.iter().map(nodes).sum::<usize>();
    let seq = |x: &[Expr], y: &[Expr]| {
        x.len().cmp(&y.len()).then_with(|| {
            for (p, q) in x.iter().zip(y) {
                let c = cmp_default(p, q);
                if c != Ordering::Equal {
                    return c;
                }
            }
            Ordering::Equal
        })
    };
    n(a).cmp(&n(b)).then_with(|| seq(&a.0, &b.0)).then_with(|| seq(&a.1, &b.1))
}

struct Radsimp {
    symbolic: bool,
    max_terms: usize,
}

impl Radsimp {
    fn handle(&self, expr: &Expr) -> Expr {
        let result = self.handle_inner(expr);
        undefined_algebraic_product(&result).unwrap_or(result)
    }

    fn handle_inner(&self, expr: &Expr) -> Expr {
        if expr.is_atom() {
            return expr.clone();
        }
        if matches!(expr.kind(), Kind::Rel(..) | Kind::Bool(_) | Kind::Piecewise(_)) {
            return super::map_args(expr, &mut |a| self.handle(a));
        }
        let (n, d) = fraction(expr, false);
        if d.is_atom() && n.is_atom() {
            return expr.clone();
        } else if !n.is_atom() {
            let n = super::map_args(&n, &mut |a| self.handle(a));
            return unevaluated_mul(&[n, self.handle(&ops::inv(&d))]);
        } else if !n.is_one() {
            return unevaluated_mul(&[n, self.handle(&ops::inv(&d))]);
        } else if d.is_mul() {
            let parts: Vec<Expr> = args(&d).iter().map(|f| self.handle(&ops::inv(f))).collect();
            return unevaluated_mul(&parts);
        }
        if !self.symbolic && !traverse::free_symbols(&d).is_empty() {
            return expr.clone();
        }
        let mut n = n;
        if ispow2(&d, self.symbolic, false) {
            let Kind::Pow(base, ex) = d.kind() else { unreachable!() };
            let d2 = ops::pow(&super::sqrtdenest::sqrtdenest(&ops::sqrt(base), 3), &numer(ex, false));
            if d2 != d {
                return self.handle(&ops::inv(&d2));
            }
        } else if let Kind::Pow(base, ex) = d.kind() {
            if assumptions::is_integer(ex) == Some(true) || assumptions::is_positive(base) == Some(true) {
                return ops::pow(&self.handle(&ops::inv(base)), ex);
            }
        }
        if !(d.is_add() || ispow2(&d, self.symbolic, false)) {
            return ops::inv(&super::map_args(&d, &mut |a| self.handle(a)));
        }
        let mut keep = true;
        let mut d = mexpand(&d);
        if d.is_atom() {
            return ops::inv(&d);
        }
        if is_num(&d) {
            let nd = crate::simp::nsimplify(&d, &[]);
            if nd.as_number().is_some() && super::methods::equals(&nd, &d) == Some(true) {
                return ops::mul(&[n, ops::inv(&nd)]);
            }
        }
        loop {
            let mut collected: Vec<(Vec<Expr>, Vec<Expr>)> = Vec::new();
            for m in add_args(&d) {
                let mut p2 = Vec::new();
                let mut other = Vec::new();
                for i in mul_args(&m) {
                    if ispow2(&i, self.symbolic, true) {
                        let Kind::Pow(b, x) = i.kind() else { unreachable!() };
                        p2.push(if *x == Expr::rational(1, 2) { b.clone() } else { ops::pow(b, &ops::mul(&[Expr::int(2), x.clone()])) });
                    } else if i.is_constant(Constant::I) {
                        p2.push(Expr::minus_one());
                    } else {
                        other.push(i);
                    }
                }
                let key = ordered(p2);
                let val = ops::mul(&other);
                if let Some(slot) = collected.iter_mut().find(|(k, _)| *k == key) {
                    slot.1.push(val);
                } else {
                    collected.push((key, vec![val]));
                }
            }
            collected.sort_by(cmp_rterm);
            let rterms: Vec<(Expr, Expr)> = collected.iter().map(|(i, j)| (ops::mul(i), ops::add(j))).collect();
            let nrad = rterms.len() - usize::from(rterms[0].0.is_one());
            if nrad < 1 {
                break;
            } else if nrad > self.max_terms {
                keep = false;
                break;
            }
            if rterms.len() > 4 {
                if rterms.iter().all(|(x, y)| x.is_integer_number() && ops::powi(y, 2).is_rational_number()) {
                    let terms: Vec<Expr> = rterms.iter().map(|(x, y)| ops::mul(&[ops::sqrt(x), y.clone()])).collect();
                    let (nd, dd) = rad_rationalize(&Expr::one(), &add_from_args(&terms));
                    n = ops::mul(&[n, nd]);
                    d = dd;
                } else {
                    keep = false;
                }
                break;
            }
            let num = super::powsimp::powsimp(&conj_num(&rterms), &super::powsimp::PowsimpOptions::default());
            n = ops::mul(&[n, num.clone()]);
            d = ops::mul(&[d, num]);
            d = super::powsimp::powdenest(&mexpand(&d), self.symbolic, false);
            if has_bad(&d) {
                return expr.clone();
            }
            if d.is_atom() {
                break;
            }
        }
        if !keep {
            return expr.clone();
        }
        unevaluated_mul(&[n, ops::inv(&d)])
    }
}

/// `radsimp(expr, symbolic, max_terms)`.
pub fn radsimp(e: &Expr, symbolic: bool, max_terms: usize) -> Expr {
    if matches!(e.kind(), Kind::Rel(..) | Kind::Bool(_) | Kind::Piecewise(_)) {
        return super::map_args(e, &mut |a| radsimp(a, symbolic, max_terms));
    }
    let rs = Radsimp { symbolic, max_terms };
    let (coeff, expr) = super::as_coeff_add(e, false);
    let expr = super::normal(&expr);
    let old = fraction(&expr, false);
    let (mut n, mut d) = fraction(&rs.handle(&expr), false);
    if old != (n.clone(), d.clone()) {
        if !d.is_atom() {
            let was = (n.clone(), d.clone());
            n = super::misc::signsimp(&n, false);
            d = super::misc::signsimp(&d, false);
            let u = Factors::from_expr(&unevaluated_mul(&[n.clone(), ops::inv(&d)]));
            let parts: Vec<Expr> = u.factors.iter().map(|(k, v)| ops::pow(k, v)).collect();
            let u = unevaluated_mul(&parts);
            let nd = fraction(&u, false);
            n = nd.0;
            d = nd.1;
            if old == (n.clone(), d.clone()) {
                n = was.0;
                d = was.1;
            }
        }
        n = super::expand_mul(&n);
        if is_number_atom(&d) || d.is_add() {
            let (n2, d2) = fraction(&gcd_terms(&unevaluated_mul(&[n.clone(), ops::inv(&d)]), false, true, true), false);
            if is_number_atom(&d2) || super::methods::count_ops(&d2) <= super::methods::count_ops(&d) {
                n = super::misc::signsimp(&n2, true);
                d = super::misc::signsimp(&d2, true);
                if n.is_mul() && is_number_atom(&args(&n)[0]) {
                    n = ops::mul(&args(&n));
                }
            }
        }
    }
    let result = ops::add(&[coeff, unevaluated_mul(&[n, ops::inv(&d)])]);
    undefined_algebraic_product(&result).unwrap_or(result)
}

/// Rationalization can expose zoo before recognizing that the remaining
/// radical numerator is exactly zero. Decide only that concrete product:
/// no symbolic/host query and no numerical cancellation proves a zero.
fn undefined_algebraic_product(e: &Expr) -> Option<Expr> {
    let Kind::Mul(coefficient, factors) = e.kind() else { return None };
    if !factors.iter().any(|(base, exponent)|
        base.is_constant(Constant::ComplexInfinity) && exponent.is_one()) {
        return None;
    }
    let mut finite = vec![Expr::num(coefficient.clone())];
    for (base, exponent) in factors {
        if base.is_constant(Constant::ComplexInfinity) && exponent.is_one() { continue; }
        finite.push(ops::pow(base, exponent));
    }
    let finite = ops::mul(&finite);
    if finite.has_host() || !traverse::free_symbols(&finite).is_empty() { return None; }
    let zero = match assumptions::is_zero_exact(&finite) {
        Some(value) => value,
        None => assumptions::algebraic_zero(&finite),
    };
    zero.then(|| Expr::constant(Constant::NaN))
}

/// `rad_rationalize(num, den)`.
pub fn rad_rationalize(num: &Expr, den: &Expr) -> (Expr, Expr) {
    if !den.is_add() {
        return (num.clone(), den.clone());
    }
    let (g, a, b) = split_surds(den);
    let a = ops::mul(&[a, ops::sqrt(&g)]);
    let num = mexpand(&ops::mul(&[ops::sub(&a, &b), num.clone()]));
    let den = mexpand(&ops::sub(&ops::powi(&a, 2), &ops::powi(&b, 2)));
    rad_rationalize(&num, &den)
}

fn gcd2(a: &Expr, b: &Expr) -> Expr {
    match (a.as_number().and_then(Number::to_bigint), b.as_number().and_then(Number::to_bigint)) {
        (Some(x), Some(y)) => {
            use num_integer::Integer;
            Expr::num(Number::from_bigint(x.gcd(&y)))
        }
        _ => crate::poly::api::gcd(a, b, &[]).unwrap_or_else(|_| Expr::one()),
    }
}

/// `_split_gcd(*a)`.
pub fn split_gcd(a: &[Expr]) -> (Expr, Vec<Expr>, Vec<Expr>) {
    let mut g = a[0].clone();
    let mut b1 = vec![g.clone()];
    let mut b2 = Vec::new();
    for x in &a[1..] {
        let g1 = gcd2(&g, x);
        if g1.is_one() {
            b2.push(x.clone());
        } else {
            g = g1;
            b1.push(x.clone());
        }
    }
    (g, b1, b2)
}

/// `split_surds(expr)`: `(g, a, b)` with `expr == sqrt(g)*a + b`.
pub fn split_surds(e: &Expr) -> (Expr, Expr, Expr) {
    let mut terms = args(e);
    terms.sort_by(cmp_default);
    let coeff_muls: Vec<(Expr, Expr)> = terms.iter().map(|x| super::as_coeff_mul(x, false)).collect();
    let mut surds: Vec<Expr> = coeff_muls.iter().filter(|(_, s)| s.is_pow()).map(|(_, s)| ops::powi(s, 2)).collect();
    surds.sort_by(cmp_default);
    let (g, b1, mut b2) = split_gcd(&surds);
    let mut g2 = g.clone();
    if b2.is_empty() && b1.len() >= 2 {
        let b1n: Vec<Expr> = b1.iter().map(|x| ops::div(x, &g)).filter(|x| !x.is_one()).collect();
        let (g1, _, nb2) = split_gcd(&b1n);
        b2 = nb2;
        g2 = ops::mul(&[g.clone(), g1]);
    }
    let _ = b2;
    let mut a1v = Vec::new();
    let mut a2v = Vec::new();
    for (c, s) in &coeff_muls {
        if let Kind::Pow(base, x) = s.kind() {
            if *x == Expr::rational(1, 2) {
                if b1.contains(base) {
                    a1v.push(ops::mul(&[c.clone(), ops::sqrt(&ops::div(base, &g2))]));
                } else {
                    a2v.push(ops::mul(&[c.clone(), s.clone()]));
                }
                continue;
            }
        }
        a2v.push(ops::mul(&[c.clone(), s.clone()]));
    }
    (g2, ops::add(&a1v), ops::add(&a2v))
}

#[cfg(test)]
mod numeric_denominator_tests {
    use super::*;

    fn parse(s: &str) -> Expr { crate::parser::parse_str(s).unwrap() }

    #[test]
    fn rational_square_classes_are_compared_exactly() {
        let sqrt8 = Expr::from_unevaluated_kind(Kind::Pow(Expr::int(8), Expr::rational(1, 2)));
        let sum = ops::add(&[Expr::int(3), sqrt8, ops::mul2(&Expr::int(-2), &ops::sqrt(&Expr::int(2)))]);
        assert_eq!(square_root_sum(&sum), Some(SquareRootSum::Rational(Rat::from_i64(3))));
        assert_eq!(square_root_sum(&parse("sqrt(2)+sqrt(3)+sqrt(5)")), Some(SquareRootSum::Irrational));
        assert_eq!(square_root_sum(&parse("sqrt(2)+sqrt(3)+sqrt(5)+sqrt(7)+sqrt(11)")), Some(SquareRootSum::Irrational));
        assert_eq!(square_root_sum(&parse("sqrt(2)+I")), None);
    }

    #[test]
    fn proven_irrational_denominators_do_not_use_recognition() {
        let expression = parse("1/(sqrt(2)+sqrt(3)+sqrt(5))");
        let result = prepare_numeric_denominators(&expression, &mut || Ok::<_, ()>(()), &mut |_| {
            panic!("a proved irrational square-root sum needs no numerical recognition")
        }).unwrap();
        assert_eq!(result, expression);
    }

    #[test]
    fn unknown_denominators_keep_fallback_and_propagate_interrupts() {
        let expression = parse("1/(exp(1)+sqrt(2))+1/(pi+sqrt(3))");
        let mut calls = 0;
        let result = prepare_numeric_denominators(&expression, &mut || Ok::<_, &str>(()), &mut |_| {
            calls += 1;
            Err("interrupt")
        });
        assert_eq!(result, Err("interrupt"));
        assert_eq!(calls, 1);
    }

    #[test]
    fn undefined_radical_product_requires_an_exact_zero_certificate() {
        let zero = parse("sqrt(5-sqrt(5))*sqrt(40-8*sqrt(5))-2*sqrt(2)*(5-sqrt(5))");
        assert!(!zero.is_zero());
        let held = |factor: Expr| Expr::from_unevaluated_kind(Kind::Mul(
            Number::ONE, vec![(Expr::constant(Constant::ComplexInfinity), Expr::one()),
                (factor, Expr::one())],
        ));
        assert_eq!(undefined_algebraic_product(&held(zero.clone())), Some(Expr::nan()));
        let tiny = parse("1/10000000000000000000000000000000000000000");
        assert_eq!(undefined_algebraic_product(&held(ops::add2(&zero, &tiny))), None);
        assert_eq!(undefined_algebraic_product(&held(parse("x"))), None);
    }
}
