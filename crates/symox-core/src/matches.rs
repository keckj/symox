//! SymPy 1.14's `matches` protocol, natively: `expr.match(pattern)`.
//!
//! This is the reference algorithm of sympy/core (basic.py `Basic.matches`,
//! operations.py `AssocOp._matches_commutative`, add.py and mul.py
//! `_matches_simple`/`_combine_inverse`, mul.py `as_base_exp`, power.py
//! `Pow.matches`, symbol.py `Wild.matches`): the same traversal order, the same elimination of the
//! exact part, the same rewriting of integer powers and coefficients, the
//! same collection on non-Wild symbols and the same `old` flag, so that the
//! first match found, among several valid ones, is the reference's.
//!
//! [`crate::pattern::pattern_match`] is a different (backtracking)
//! algorithm kept for its internal callers.  Here, an input the port does
//! not reproduce exactly is declined ([`Stop::Decline`]) so that the caller
//! falls back to the reference implementation; a decline never stands for
//! "no match".
//!
//! Declined inputs: noncommutative products, matrix, host, binder, limit,
//! Order, logic, set, relational, Piecewise, Max/Min nodes and infinities;
//! a Float facing a side without free symbols in a multiplicative inverse
//! (where the reference compares numerically); a collection whose
//! result depends on the (hash) order of more than one symbol; and a search
//! that exceeds a step budget.

use rustc_hash::FxHashMap;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::{decompose, ops, sortkey, subs, traverse};

pub type MatchError = crate::symbol::WildcardError;

/// Why a match was not decided.
#[derive(Debug)]
pub enum Stop {
    /// The port does not reproduce the reference on this input.
    Decline(&'static str),
    /// A Wild property raised.
    Error(MatchError),
}

type Outcome = Result<Option<Bindings>, Stop>;

/// Steps (calls of `matches`) before the search is declined.
const BUDGET: usize = 200_000;

/// Insertion-ordered bindings; binding a Wild again overwrites it, as
/// `repl_dict[self] = expr` in the reference.
#[derive(Clone, Default, Debug)]
pub struct Bindings(Vec<(Expr, Expr)>);

impl Bindings {
    fn with(&self, key: &Expr, value: &Expr) -> Bindings {
        let mut out = self.clone();
        match out.0.iter_mut().find(|(k, _)| k == key) {
            Some(slot) => slot.1 = value.clone(),
            None => out.0.push((key.clone(), value.clone())),
        }
        out
    }
    fn map(&self) -> FxHashMap<Expr, Expr> {
        self.0.iter().cloned().collect()
    }
    pub fn into_vec(self) -> Vec<(Expr, Expr)> {
        self.0
    }
}

/// `expr.match(pattern, old)` of SymPy 1.14.
pub fn reference_match(expr: &Expr, pattern: &Expr, old: bool) -> Result<Option<FxHashMap<Expr, Expr>>, Stop> {
    Ok(reference_match_ordered(expr, pattern, old)?.map(|b| b.map()))
}

/// [`reference_match`] keeping the order in which the Wilds were bound.
pub fn reference_match_ordered(expr: &Expr, pattern: &Expr, old: bool) -> Outcome {
    admit(pattern)?;
    admit(expr)?;
    Matcher { steps: 0 }.matches(pattern, expr, &Bindings::default(), old)
}

/// Inputs whose reference semantics the port does not cover.
fn admit(e: &Expr) -> Result<(), Stop> {
    if e.has_noncommutative_host() || e.has_host() || e.has_binding_identity() || e.has_bound_generator()
        || e.has_limit() || e.has_symbolic_derivative() || e.has_nonfinite()
    {
        return Err(Stop::Decline("unsupported node class"));
    }
    fn visit(e: &Expr, seen: &mut rustc_hash::FxHashSet<usize>) -> Result<(), Stop> {
        if !seen.insert(e.dag_id()) {
            return Ok(());
        }
        let ok = match e.kind() {
            Kind::Num(_) | Kind::Sym(_) => true,
            Kind::Const(c) => !matches!(c, Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN),
            Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => true,
            Kind::Derivative(body, variables) => {
                visit(body, seen)?;
                for (v, _) in variables {
                    visit(v, seen)?;
                }
                return Ok(());
            }
            Kind::Func(f, _) => match f {
                Func::Undefined(head) => !head.internal,
                Func::Max | Func::Min | Func::Order | Func::Subs | Func::SymbolicDerivative | Func::Limit
                | Func::Integral | Func::Sum | Func::Product | Func::Logic(_) | Func::Set(_)
                | Func::MatrixElement(_) | Func::MatrixProxy(_) | Func::MatrixTrace(_)
                | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(..) | Func::Predicate(_)
                | Func::AppliedPredicate | Func::GenericExpr => false,
                _ => true,
            },
            Kind::Rel(..) | Kind::Piecewise(_) | Kind::Bool(_) => false,
        };
        if !ok {
            return Err(Stop::Decline("unsupported node class"));
        }
        for a in e.args() {
            visit(&a, seen)?;
        }
        Ok(())
    }
    visit(e, &mut Default::default())
}

fn is_wild(e: &Expr) -> bool {
    e.as_sym().is_some_and(|s| s.is_wild())
}

fn has_wild(e: &Expr) -> bool {
    crate::pattern::has_wild(e)
}

fn args(e: &Expr) -> Vec<Expr> {
    traverse::ordered_args(e)
}

/// `Basic.has` of the reference: a sum or product query is also found as
/// a subset of the arguments of a sum or product (`AssocOp._has_matcher`).
pub fn reference_has(e: &Expr, q: &Expr) -> bool {
    if traverse::has_structural(e, q) {
        return true;
    }
    if !q.is_mul() || q.has_noncommutative_host() {
        return false;
    }
    let wanted = args(q);
    let mut found = false;
    let mut seen = rustc_hash::FxHashSet::default();
    fn visit(e: &Expr, wanted: &[Expr], found: &mut bool, seen: &mut rustc_hash::FxHashSet<usize>) {
        if *found || !seen.insert(e.dag_id()) {
            return;
        }
        if e.is_mul() {
            let have = args(e);
            if wanted.iter().all(|w| have.contains(w)) {
                *found = true;
                return;
            }
        }
        for a in args(e) {
            visit(&a, wanted, found, seen);
        }
    }
    visit(e, &wanted, &mut found, &mut seen);
    found
}

fn is_number(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(_))
}

fn is_float(e: &Expr) -> bool {
    e.as_number().is_some_and(Number::is_float)
}

/// `_nodes` of `sympy.core.sorting`, doubled (a Float counts 1/2).
fn nodes2(e: &Expr) -> usize {
    fn count(e: &Expr) -> usize {
        if is_float(e) {
            return 1;
        }
        2 + args(e).iter().map(count).sum::<usize>()
    }
    if let Kind::Derivative(body, variables) = e.kind() {
        return nodes2(body) + variables.iter().map(|(_, n)| 2 * *n as usize).sum::<usize>();
    }
    count(e)
}

/// `ordered(seq)`: stable, by node count then `default_sort_key`.
fn ordered(mut v: Vec<Expr>) -> Vec<Expr> {
    v.sort_by(|a, b| nodes2(a).cmp(&nodes2(b)).then_with(|| sortkey::cmp_default(a, b)));
    v
}

/// The key `term.args[0]` of a product with a Number coefficient, else 0.
fn numeric_first(t: &Expr) -> Number {
    if t.is_mul() {
        if let Some(n) = args(t).first().and_then(|a| a.as_number()) {
            return n.clone();
        }
    }
    Number::ZERO
}

fn sort_numeric_first(v: &mut [Expr]) -> Result<(), Stop> {
    let keys: Vec<Number> = v.iter().map(numeric_first).collect();
    for k in &keys {
        if k.cmp_value(&Number::ZERO).is_none() {
            return Err(Stop::Decline("unordered numeric key"));
        }
    }
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| keys[a].cmp_value(&keys[b]).unwrap_or(std::cmp::Ordering::Equal));
    let sorted: Vec<Expr> = idx.iter().map(|&i| v[i].clone()).collect();
    v.clone_from_slice(&sorted);
    Ok(())
}

/// `_from_args`: a node holding exactly `values` (identity when empty).
fn raw(add: bool, values: &[Expr]) -> Expr {
    if add {
        crate::simptools::add_from_args(values)
    } else {
        crate::simptools::mul_from_args(values)
    }
}

fn make_args(add: bool, e: &Expr) -> Vec<Expr> {
    if (add && e.is_add()) || (!add && e.is_mul()) {
        args(e)
    } else {
        vec![e.clone()]
    }
}

fn count_ops(e: &Expr) -> Result<num_bigint::BigInt, Stop> {
    match crate::inspection::operation_counts::count(e) {
        Ok(Some(c)) => Ok(c.total()),
        _ => Err(Stop::Decline("count_ops")),
    }
}

/// `_coeff_isneg`.
fn coeff_isneg(a: &Expr) -> bool {
    let first = if a.is_mul() { args(a).into_iter().next().unwrap_or_else(|| a.clone()) } else { a.clone() };
    first.as_number().is_some_and(Number::is_negative)
}

/// `Add._combine_inverse` (infinities are declined on admission).
fn add_inverse(lhs: &Expr, rhs: &Expr) -> Expr {
    let rv = ops::sub(lhs, rhs);
    let srv = crate::simptools::misc::signsimp(&rv, true);
    if is_number(&srv) { srv } else { rv }
}

/// `Mul._combine_inverse`.
fn mul_inverse(lhs: &Expr, rhs: &Expr) -> Result<Expr, Stop> {
    if lhs == rhs {
        return Ok(Expr::one());
    }
    // `check`: a Float against a comparable number compares numerically.
    if (is_float(lhs) && traverse::free_symbols(rhs).is_empty()) || (is_float(rhs) && traverse::free_symbols(lhs).is_empty()) {
        return Err(Stop::Decline("Float in a multiplicative inverse"));
    }
    let (mut lhs, mut rhs) = (lhs.clone(), rhs.clone());
    if [&lhs, &rhs].iter().any(|i| i.is_pow() || i.is_mul()) {
        // A literal I does not combine with a power of -1.
        let d = Expr::dummy("I", crate::symbol::Assumptions::NONE);
        let to: FxHashMap<Expr, Expr> = [(Expr::i(), d.clone())].into_iter().collect();
        let back: FxHashMap<Expr, Expr> = [(d, Expr::i())].into_iter().collect();
        let mut a = decompose::as_powers_dict(&subs::xreplace(&lhs, &to));
        let mut b = decompose::as_powers_dict(&subs::xreplace(&rhs, &to));
        let blen = b.len();
        let keys: Vec<Expr> = b.iter().map(|(k, _)| k.clone()).collect();
        for bi in keys {
            if let Some(ia) = a.iter().position(|(k, _)| *k == bi) {
                let ib = b.iter().position(|(k, _)| *k == bi).unwrap();
                let (_, eb) = b.remove(ib);
                a[ia].1 = ops::sub(&a[ia].1, &eb);
                if a[ia].1.is_zero() {
                    a.remove(ia);
                }
            }
        }
        if b.len() != blen {
            let product = |v: &[(Expr, Expr)]| ops::mul(&v.iter().map(|(k, e)| ops::pow(k, e)).collect::<Vec<_>>());
            lhs = subs::xreplace(&product(&a), &back);
            rhs = subs::xreplace(&product(&b), &back);
        }
    }
    let rv = ops::div(&lhs, &rhs);
    let srv = crate::simptools::misc::signsimp(&rv, true);
    Ok(if is_number(&srv) { srv } else { rv })
}

/// `collect(expr, symbols)` with the reference's (hash-ordered) set; a
/// result that depends on the order is declined.
fn collect_set(e: &Expr, symbols: &[Expr]) -> Result<Expr, Stop> {
    use crate::poly::collect::{run, Collection, Options};
    let once = |order: &[Expr]| -> Result<Expr, Stop> {
        let options = Options { exact: false, discover: true, distribute_order_term: true };
        match run::<()>(e, order, options, true, None) {
            Ok(Collection::Expression(v)) => Ok(v),
            _ => Err(Stop::Decline("collect")),
        }
    };
    if symbols.len() > 3 {
        return Err(Stop::Decline("collect on more than three symbols"));
    }
    let first = once(symbols)?;
    let n = symbols.len();
    let mut perm: Vec<usize> = (0..n).collect();
    // Heap's algorithm over the remaining orders.
    let mut c = vec![0usize; n];
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 { perm.swap(0, i) } else { perm.swap(c[i], i) }
            let order: Vec<Expr> = perm.iter().map(|&k| symbols[k].clone()).collect();
            if !once(&order)?.same_presentation(&first) {
                return Err(Stop::Decline("collect depends on the symbol order"));
            }
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    Ok(first)
}

struct Matcher {
    steps: usize,
}

/// `as_base_exp` with `Mul.as_base_exp`: a product whose factors share one
/// integer exponent is the power of the product of their bases.
fn base_exp(e: &Expr) -> (Expr, Expr) {
    if !e.is_mul() {
        return decompose::as_base_exp(e);
    }
    let mut first: Option<Expr> = None;
    let mut bases = Vec::new();
    for factor in args(e) {
        let (b, x) = decompose::as_base_exp(&factor);
        match &first {
            None => first = Some(x),
            Some(f) => {
                if &x != f || !x.is_integer_number() {
                    return (e.clone(), Expr::one());
                }
            }
        }
        bases.push(b);
    }
    (ops::mul(&bases), first.unwrap_or_else(Expr::one))
}

impl Matcher {
    fn matches(&mut self, pattern: &Expr, expr: &Expr, d: &Bindings, old: bool) -> Outcome {
        self.steps += 1;
        if self.steps > BUDGET {
            return Err(Stop::Decline("step budget"));
        }
        if is_wild(pattern) {
            return self.wild(pattern, expr, d);
        }
        match pattern.kind() {
            Kind::Add(..) => self.commutative(pattern, expr, d, old),
            Kind::Mul(..) => self.commutative(pattern, expr, d, old),
            Kind::Pow(..) => self.pow(pattern, expr, d, old),
            Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => {
                Ok((pattern == expr).then(|| d.clone()))
            }
            _ => self.basic(pattern, expr, d, old),
        }
    }

    /// `Wild.matches`.
    fn wild(&mut self, pattern: &Expr, expr: &Expr, d: &Bindings) -> Outcome {
        let s = pattern.as_sym().unwrap();
        if s.exclude().iter().any(|x| reference_has(expr, x)) {
            return Ok(None);
        }
        if let Some(properties) = s.wild_properties() {
            match properties.evaluate(expr) {
                Ok(true) => {}
                Ok(false) => return Ok(None),
                Err(error) => return Err(Stop::Error(error)),
            }
        }
        Ok(Some(d.with(pattern, expr)))
    }

    /// `Basic.matches`.
    fn basic(&mut self, pattern: &Expr, expr: &Expr, d: &Bindings, old: bool) -> Outcome {
        let same = match (pattern.kind(), expr.kind()) {
            (Kind::Func(f, _), Kind::Func(g, _)) => f == g,
            (Kind::Derivative(..), Kind::Derivative(..)) => true,
            (Kind::Pow(..), Kind::Pow(..)) => true,
            _ => false,
        };
        if !same {
            return Ok(None);
        }
        if pattern == expr {
            return Ok(Some(d.clone()));
        }
        let (pargs, eargs) = (args(pattern), args(expr));
        if pargs.len() != eargs.len() {
            return Ok(None);
        }
        let mut d = d.clone();
        for (arg, other) in pargs.iter().zip(&eargs) {
            if arg == other {
                continue;
            }
            match self.matches(&subs::xreplace(arg, &d.map()), other, &d, old)? {
                Some(next) => d = next,
                None => return Ok(None),
            }
        }
        Ok(Some(d))
    }

    /// `Pow.matches`.
    fn pow(&mut self, pattern: &Expr, expr: &Expr, d: &Bindings, old: bool) -> Outcome {
        let Kind::Pow(base, exponent) = pattern.kind() else { unreachable!() };
        if expr.is_one() {
            if let Some(found) = self.matches(exponent, &Expr::zero(), d, false)? {
                return Ok(Some(found));
            }
        }
        if let Kind::Func(Func::Undefined(head), _) = expr.kind() {
            if head.internal {
                return Ok(None);
            }
        }
        let (b, e) = base_exp(expr);
        let (sb, se) = base_exp(pattern);
        let nonzero = !(is_number(expr) && expr.is_zero());
        if sb.is_symbol() && se.is_integer_number() && nonzero {
            let target = if crate::assumptions::is_rational(&e) == Some(true) {
                ops::pow(&b, &ops::div(&e, &se))
            } else {
                ops::pow(expr, &ops::inv(&se))
            };
            return self.matches(&sb, &target, d, false);
        }
        let Some(found) = self.matches(base, &b, d, false)? else { return Ok(None) };
        match self.matches(&subs::xreplace(exponent, &found.map()), &e, &found, false)? {
            Some(found) => Ok(Some(found)),
            None => self.basic(pattern, expr, d, old),
        }
    }

    /// `AssocOp._matches_commutative` with the `_matches_simple` and
    /// `_combine_inverse` of Add and Mul.
    fn commutative(&mut self, pattern: &Expr, expr: &Expr, d: &Bindings, old: bool) -> Outcome {
        let add = pattern.is_add();
        if let Kind::Func(Func::Undefined(head), _) = expr.kind() {
            if head.internal {
                return Ok(None);
            }
        }
        if pattern == expr {
            return Ok(Some(d.clone()));
        }
        // _matches_simple
        let simple = if add {
            let (coeff, terms) = decompose::as_coeff_add(pattern, &[]);
            if terms.len() == 1 { self.matches(&terms[0], &ops::sub(expr, &coeff), d, false)? } else { None }
        } else {
            let (coeff, terms) = crate::simptools::as_coeff_mul(pattern, false);
            let terms = make_args(false, &terms);
            if terms.len() == 1 { self.matches(&terms[0], &mul_inverse(expr, &coeff)?, d, false)? } else { None }
        };
        if simple.is_some() {
            return Ok(simple);
        }
        // eliminate the exact part
        let (mut wild_part, mut exact_part) = (Vec::new(), Vec::new());
        for p in args(pattern) {
            if has_wild(&p) && !reference_has(expr, &p) { wild_part.push(p) } else { exact_part.push(p) }
        }
        if exact_part.is_empty() {
            wild_part = ordered(wild_part);
            if add {
                sort_numeric_first(&mut wild_part)?;
            }
        } else {
            let exact = raw(add, &exact_part);
            let free = traverse::free_symbols(expr);
            if !free.is_empty() && traverse::free_symbols(&exact).iter().any(|s| !free.contains(s)) {
                return Ok(None);
            }
            let newexpr = if add { add_inverse(expr, &exact) } else { mul_inverse(expr, &exact)? };
            if !old && (expr.is_add() || expr.is_mul()) {
                let check = if coeff_isneg(&newexpr) { ops::neg(&newexpr) } else { newexpr.clone() };
                if count_ops(&check)? > count_ops(expr)? {
                    return Ok(None);
                }
            }
            let newpattern = raw(add, &wild_part);
            return self.matches(&newpattern, &newexpr, d, false);
        }
        let identity = if add { Expr::zero() } else { Expr::one() };
        let mut expr = expr.clone();
        let mut i = 0;
        let mut saw: Vec<Expr> = Vec::new();
        while !saw.contains(&expr) {
            saw.push(expr.clone());
            let mut terms = ordered(make_args(add, &expr));
            if add && expr.is_add() {
                sort_numeric_first(&mut terms)?;
            }
            let mut expr_list = vec![identity.clone()];
            expr_list.extend(terms);
            for last_op in expr_list.iter().rev() {
                for w in wild_part.iter().rev() {
                    if let Some(d1) = self.matches(w, last_op, d, false)? {
                        let instance = subs::xreplace(pattern, &d1.map());
                        if let Some(d2) = self.matches(&instance, &expr, &d1, false)? {
                            return Ok(Some(d2));
                        }
                    }
                }
            }
            if i == 0 {
                if !add {
                    // make e**i look like Mul
                    if let Kind::Pow(b, e) = expr.kind() {
                        if let Some(n) = e.as_number().filter(|n| n.is_integer()) {
                            expr = if n.is_positive() {
                                raw(false, &[b.clone(), ops::pow(b, &ops::sub(e, &Expr::one()))])
                            } else {
                                raw(false, &[ops::inv(b), ops::pow(b, &ops::add2(e, &Expr::one()))])
                            };
                            i += 1;
                            continue;
                        }
                    }
                } else {
                    // make i*e look like Add
                    let (c, e) = crate::simptools::as_coeff_mul(&expr, false);
                    if let Some(n) = c.as_number() {
                        if n.abs().cmp_value(&Number::ONE) == Some(std::cmp::Ordering::Greater) {
                            expr = if n.is_positive() {
                                raw(true, &[e.clone(), ops::mul2(&ops::sub(&c, &Expr::one()), &e)])
                            } else {
                                raw(true, &[ops::neg(&e), ops::mul2(&ops::add2(&c, &Expr::one()), &e)])
                            };
                            i += 1;
                            continue;
                        }
                    }
                    // try collection on non-Wild symbols
                    let was = expr.clone();
                    let mut did: Vec<Expr> = Vec::new();
                    for w in wild_part.iter().rev() {
                        let c = ops::mul(&make_args(false, w).into_iter().filter(|f| !has_wild(f)).collect::<Vec<_>>());
                        let mut free: Vec<Expr> = traverse::free_symbols(&c).into_iter().filter(|s| !did.contains(s)).collect();
                        if !free.is_empty() {
                            free.sort_by(sortkey::cmp_default);
                            did.extend(free.iter().cloned());
                            expr = collect_set(&expr, &free)?;
                        }
                    }
                    if !expr.same_presentation(&was) {
                        continue;
                    }
                }
            }
            break;
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sym(n: &str) -> Expr {
        Expr::symbol(n)
    }

    fn check(expr: &Expr, pattern: &Expr, expected: Option<Vec<(Expr, Expr)>>) {
        let got = reference_match(expr, pattern, false).unwrap_or_else(|_| panic!("declined"));
        let expected = expected.map(|v| v.into_iter().collect::<FxHashMap<_, _>>());
        assert_eq!(got, expected, "{expr} against {pattern}");
    }

    #[test]
    fn reference_choices() {
        let (x, y) = (sym("x"), sym("y"));
        let a = Expr::wild("a", vec![x.clone()]);
        let b = Expr::wild("b", vec![x.clone()]);
        // (2*x + 3).match(a*x + b) -> {a: 2, b: 3}
        let e = ops::add2(&ops::mul2(&Expr::int(2), &x), &Expr::int(3));
        let p = ops::add2(&ops::mul2(&a, &x), &b);
        check(&e, &p, Some(vec![(a.clone(), Expr::int(2)), (b.clone(), Expr::int(3))]));
        // (x**2).match(c**2*x**(2*a)) -> {a: 1, c: 1}
        let c = Expr::wild("c", vec![x.clone()]);
        let p = ops::mul2(&ops::pow(&c, &Expr::int(2)), &ops::pow(&x, &ops::mul2(&Expr::int(2), &a)));
        check(&ops::pow(&x, &Expr::int(2)), &p, Some(vec![(a.clone(), Expr::one()), (c.clone(), Expr::one())]));
        // (2*x*y*z).has(x*y) as in the reference
        let z = sym("z");
        let e = ops::mul(&[Expr::int(2), x.clone(), y.clone(), z]);
        assert!(reference_has(&e, &ops::mul2(&x, &y)));
        assert!(!traverse::has_structural(&e, &ops::mul2(&x, &y)) || reference_has(&e, &ops::mul2(&x, &y)));
        // y does not match a*x (a excludes x but 0 is not y)
        check(&y, &ops::mul2(&a, &x), None);
    }
}
