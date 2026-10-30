//! Pattern matching of the first-order `classify_ode` hints.
//!
//! The reference matches each first-order hint by `eq.match(pattern)` on a
//! collected form of the equation. For the patterns used there the match is
//! unique (every wild excludes the derivative or the function), so it is the
//! additive decomposition computed here. Each step mirrors the reference
//! matcher (`_add_inverse`, `_mul_inverse`, the `count_ops` guard, the
//! rebuilt `pattern.xreplace(...)` that must reproduce the matched term). A
//! shape whose reference outcome depends on matcher details that are not
//! modelled declines, and the caller runs the reference matcher instead.
//!
//! * `de`: `d + e*df` (`separable`, `1st_exact`, `almost_linear`,
//!   `1st_homogeneous_coeff_*`, `linear_coefficients`, `separable_reduced`,
//!   the order-one power series match).
//! * `linear`: `df + P*f - Q` on `collect(expand(eq/eq.coeff(df)))`.
//! * `bernoulli`: `df + P*f - Q*f**n` on the same form.
//! * `riccati_special`: `a*df + b*f**2 + c*f/x + d/x**2`, rejected early when
//!   no match is possible, otherwise matched on the collected form.

use num_bigint::BigInt;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::ops::{div, mul, neg, sub};
use crate::poly::collect::{self, Collection, Options};
use crate::traverse::{has, preorder};

/// Result of one pattern match.
#[derive(Clone, Debug)]
pub enum Outcome<T> {
    Match(T),
    /// The reference matcher provably returns `None`.
    NoMatch,
    /// Not decided natively: the reference matcher must run.
    Decline,
}

/// The operations whose public semantics the matcher depends on. The
/// binding supplies the exact public `cancel`, `expand` and `count_ops`.
/// `None` from `cancel`/`expand` declines.
pub struct Hooks<'a> {
    pub cancel: &'a mut dyn FnMut(&Expr) -> Option<Expr>,
    pub expand: &'a mut dyn FnMut(&Expr) -> Option<Expr>,
    pub count_ops: &'a mut dyn FnMut(&Expr) -> BigInt,
}

/// The first-order pattern matches of one equation.
#[derive(Clone, Debug)]
pub struct FirstOrder {
    /// `(d, e)` of `d + e*df` on `collect(eq, [df, f], cancel)`.
    pub de: Outcome<(Expr, Expr)>,
    /// `(P, Q)` of `df + P*f - Q`.
    pub linear: Outcome<(Expr, Expr)>,
    /// `(P, Q, n)` of `df + P*f - Q*f**n`; `n` is unbound (`None`) when
    /// the equation has no term free of `f` and no power of `f`.
    pub bernoulli: Outcome<(Expr, Expr, Option<Expr>)>,
    /// Never `Match`: a match is left to the reference.
    pub riccati_special: Outcome<[Expr; 4]>,
}

fn terms(e: &Expr) -> Vec<Expr> {
    if matches!(e.kind(), Kind::Add(..)) { e.args() } else { vec![e.clone()] }
}

fn factors(e: &Expr) -> Vec<Expr> {
    if matches!(e.kind(), Kind::Mul(..)) { e.args() } else { vec![e.clone()] }
}

fn is_add_or_mul(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Add(..) | Kind::Mul(..))
}

/// `_coeff_isneg`: a negative leading number of a product, or a negative number.
fn coeff_isneg(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_negative(),
        Kind::Mul(c, _) => !c.is_one() && c.is_negative(),
        _ => false,
    }
}

/// Content whose reference matching takes branches not modelled here:
/// infinities (the `oo` dummy of `_add_inverse`), `I` (the dummy swap of
/// `_mul_inverse`), `nan`, floats (whose `cancel` entry rewrites logs),
/// piecewise and relational content.
pub(crate) fn blocked(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |x| {
        if found {
            return;
        }
        match x.kind() {
            Kind::Const(Constant::I | Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN) => found = true,
            Kind::Num(n) if n.is_float() => found = true,
            Kind::Piecewise(_) | Kind::Rel(..) | Kind::Bool(_) => found = true,
            _ => {}
        }
    });
    found
}

/// Factors of a product whose `as_powers_dict` rebuild in `_mul_inverse`
/// could combine differently from the product itself: a repeated base, a
/// number next to a power of a number, an exponential.
fn powers_rebuild_unsafe(fs: &[Expr]) -> bool {
    let mut bases: Vec<Expr> = Vec::new();
    let mut number = false;
    let mut number_base = false;
    for f in fs {
        // `exp(a).as_base_exp()` is `(E, a)`: exponentials share the base `E`.
        let base = match f.kind() {
            Kind::Pow(b, _) => b.clone(),
            Kind::Func(Func::Exp, _) => Expr::constant(Constant::E),
            Kind::Num(_) => {
                number = true;
                f.clone()
            }
            _ => f.clone(),
        };
        if matches!(base.kind(), Kind::Num(_) | Kind::Const(_)) && matches!(f.kind(), Kind::Pow(..)) {
            number_base = true;
        }
        if bases.contains(&base) {
            return true;
        }
        bases.push(base);
    }
    number && number_base
}

/// `_mul_inverse(t, factor)` for a factor `factor` of the product `t`
/// (exponent one), with the matcher's `count_ops` guard. `None` declines.
fn quotient(t: &Expr, factor: &Expr, hooks: &mut Hooks<'_>) -> Option<Expr> {
    if t == factor {
        return Some(Expr::one());
    }
    let fs = factors(t);
    if !matches!(t.kind(), Kind::Mul(..)) || powers_rebuild_unsafe(&fs) {
        return None;
    }
    let position = fs.iter().position(|f| f == factor)?;
    let rest: Vec<Expr> = fs.iter().enumerate().filter(|(i, _)| *i != position).map(|(_, f)| f.clone()).collect();
    let q = mul(&rest);
    let check = if coeff_isneg(&q) { neg(&q) } else { q.clone() };
    if (hooks.count_ops)(&check) > (hooks.count_ops)(t) {
        return None;
    }
    Some(q)
}

/// `_mul_inverse(t, -1)`. `None` declines.
fn negated(t: &Expr) -> Option<Expr> {
    let plain = div(t, &Expr::int(-1));
    if let Kind::Mul(c, _) = t.kind() {
        if c.is_minus_one() {
            let fs = factors(t);
            if powers_rebuild_unsafe(&fs) {
                return None;
            }
            let rebuilt = mul(&fs[1..]);
            if rebuilt != plain {
                return None;
            }
        }
    }
    Some(plain)
}

/// `count_ops` guard of the exact-part elimination `expr - exact`.
fn remainder_admitted(expr: &Expr, remainder: &Expr, hooks: &mut Hooks<'_>) -> bool {
    if !is_add_or_mul(expr) {
        return true;
    }
    let check = if coeff_isneg(remainder) { neg(remainder) } else { remainder.clone() };
    (hooks.count_ops)(&check) <= (hooks.count_ops)(expr)
}

/// Linear occurrence of `target` in the term `t`: `Some(Ok(c))` when
/// `t == c*target` with `c` free of `target`, `Some(Err(()))` when the term
/// is such a product whose reference quotient is not reproduced, and `None`
/// when `target` occurs otherwise.
fn linear_cofactor(t: &Expr, target: &Expr, hooks: &mut Hooks<'_>) -> Option<Result<Expr, ()>> {
    if t == target {
        return Some(Ok(Expr::one()));
    }
    if !matches!(t.kind(), Kind::Mul(..)) {
        return None;
    }
    let fs = factors(t);
    let hits: Vec<usize> = fs.iter().enumerate().filter(|(_, f)| has(f, target)).map(|(i, _)| i).collect();
    if hits.len() != 1 || fs[hits[0]] != *target {
        return None;
    }
    let Some(c) = quotient(t, target, hooks) else { return Some(Err(())) };
    if has(&c, target) {
        return None;
    }
    // `pattern.xreplace({w: c})` rebuilds `c*target`; the exact-part
    // elimination needs it to be the matched term itself.
    if mul(&[c.clone(), target.clone()]) != *t {
        return Some(Err(()));
    }
    Some(Ok(c))
}

/// `d + e*df` with `d` excluding `df` and `d_exclude`, `e` excluding `df`.
pub fn split_de(c: &Expr, df: &Expr, d_exclude: &[Expr], hooks: &mut Hooks<'_>) -> Outcome<(Expr, Expr)> {
    if c.is_zero() || blocked(c) {
        return Outcome::Decline;
    }
    let mut piece: Option<(Expr, Expr)> = None;
    for t in terms(c) {
        if !has(&t, df) {
            continue;
        }
        match linear_cofactor(&t, df, hooks) {
            None => return Outcome::NoMatch,
            Some(Err(())) => return Outcome::Decline,
            Some(Ok(e)) => {
                if piece.is_some() {
                    return Outcome::Decline;
                }
                piece = Some((t, e));
            }
        }
    }
    let (d, e) = match piece {
        None => (c.clone(), Expr::zero()),
        Some((t, e)) => {
            let d = sub(c, &t);
            if !remainder_admitted(c, &d, hooks) {
                return Outcome::Decline;
            }
            (d, e)
        }
    };
    if d_exclude.iter().any(|w| has(&d, w)) {
        return Outcome::NoMatch;
    }
    Outcome::Match((d, e))
}

/// The reference `Expr.coeff(target)` for a commutative `target`: the sum of
/// the cofactors of the terms having `target` as a factor.
pub(crate) fn coeff(e: &Expr, target: &Expr) -> Expr {
    if e == target {
        return Expr::one();
    }
    let mut out = Vec::new();
    for t in terms(e) {
        let fs = factors(&t);
        if let Some(i) = fs.iter().position(|f| f == target) {
            let rest: Vec<Expr> = fs.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, f)| f.clone()).collect();
            out.push(mul(&rest));
        }
    }
    crate::ops::add(&out)
}

/// Classification of the pieces of `R = eq2 - df` as `c*f` (at most one),
/// `c*f**k` and `f`-free pieces.
struct Pieces {
    linear: Option<(Expr, Expr)>,
    power: Vec<(Expr, Expr, Expr)>,
    free: Vec<Expr>,
}

/// `None`: some piece has `f` in another form (no representation exists);
/// `Some(Err)`: a shape not modelled.
fn pieces(r: &Expr, fx: &Expr, df: &Expr, hooks: &mut Hooks<'_>) -> Option<Result<Pieces, ()>> {
    let mut out = Pieces { linear: None, power: Vec::new(), free: Vec::new() };
    if r.is_zero() {
        return Some(Ok(out));
    }
    for t in terms(r) {
        if has(&t, df) {
            return Some(Err(()));
        }
        if !has(&t, fx) {
            out.free.push(t);
            continue;
        }
        match linear_cofactor(&t, fx, hooks) {
            Some(Err(())) => return Some(Err(())),
            Some(Ok(c)) => {
                if out.linear.is_some() {
                    return Some(Err(()));
                }
                out.linear = Some((t, c));
                continue;
            }
            None => {}
        }
        // c*f**k
        let fs = factors(&t);
        let hits: Vec<usize> = fs.iter().enumerate().filter(|(_, f)| has(f, fx)).map(|(i, _)| i).collect();
        if hits.len() != 1 {
            return None;
        }
        let Kind::Pow(b, k) = fs[hits[0]].kind() else { return None };
        if b != fx || has(k, fx) {
            return None;
        }
        out.power.push((t.clone(), fs[hits[0]].clone(), k.clone()));
    }
    Some(Ok(out))
}

/// `collect(e, [df, f], func=cancel)`.
pub(crate) fn collected(e: &Expr, df: &Expr, fx: &Expr, hooks: &mut Hooks<'_>) -> Option<Expr> {
    let cancel = &mut hooks.cancel;
    let mut callback = |value: &Expr| -> Result<Expr, ()> { cancel(value).ok_or(()) };
    let options = Options { exact: false, discover: false, distribute_order_term: true };
    match collect::run(e, &[df.clone(), fx.clone()], options, true, Some(&mut callback)) {
        Ok(Collection::Expression(v)) => Some(v),
        _ => None,
    }
}

/// `df + P*f - Q` and `df + P*f - Q*f**n` on `eq2`.
fn linear_forms(eq2: &Expr, fx: &Expr, x: &Expr, df: &Expr, hooks: &mut Hooks<'_>) -> (Outcome<(Expr, Expr)>, Outcome<(Expr, Expr, Option<Expr>)>) {
    let decline = (Outcome::Decline, Outcome::Decline);
    if blocked(eq2) || !terms(eq2).iter().any(|t| t == df) {
        return decline;
    }
    let r = sub(eq2, df);
    if !remainder_admitted(eq2, &r, hooks) {
        return decline;
    }
    let p = match pieces(&r, fx, df, hooks) {
        None => return (Outcome::NoMatch, Outcome::NoMatch),
        Some(Err(())) => return decline,
        Some(Ok(p)) => p,
    };
    // The exact part `P*f` leaves `F = R - c*f`.
    let (p_value, rest) = match &p.linear {
        Some((t, c)) => {
            let rest = sub(&r, t);
            if !remainder_admitted(&r, &rest, hooks) {
                return decline;
            }
            (c.clone(), rest)
        }
        None => (Expr::zero(), r.clone()),
    };
    let linear = if !p.power.is_empty() {
        Outcome::NoMatch
    } else if rest.is_zero() {
        Outcome::Match((p_value.clone(), Expr::zero()))
    } else if p.linear.is_none() && p.free.len() == 1 {
        // `-Q` against the single term, then `P*f` against zero.
        match negated(&rest) {
            Some(q) if neg(&q) == rest => Outcome::Match((Expr::zero(), q)),
            _ => Outcome::Decline,
        }
    } else {
        // `_mul_simple`: Q = rest/(-1).
        match negated(&rest) {
            Some(q) => Outcome::Match((p_value.clone(), q)),
            None => Outcome::Decline,
        }
    };
    let bernoulli = match (p.power.len(), p.free.len()) {
        (0, 0) => Outcome::Match((p_value.clone(), Expr::zero(), None)),
        (1, 0) => {
            let (t, power, k) = &p.power[0];
            if has(k, x) || has(k, df) || matches!(power.kind(), Kind::Pow(_, k) if k.is_one()) {
                Outcome::Decline
            } else {
                match bernoulli_power(t, power, hooks) {
                    Some(q) => Outcome::Match((p_value.clone(), q, Some(k.clone()))),
                    None => Outcome::Decline,
                }
            }
        }
        (0, 1) => {
            let t = &p.free[0];
            match negated(t) {
                Some(q) if neg(&q) == *t => Outcome::Match((p_value.clone(), q, Some(Expr::zero()))),
                _ => Outcome::Decline,
            }
        }
        (0, _) => Outcome::Decline,
        _ => Outcome::NoMatch,
    };
    (linear, bernoulli)
}

/// `-Q*f**n` against `t = B*f**k`: `Q = -B` through `_mul_inverse(t, -1)`
/// and `_mul_inverse(-t, f**k)`.
fn bernoulli_power(t: &Expr, power: &Expr, hooks: &mut Hooks<'_>) -> Option<Expr> {
    let n = negated(t)?;
    if is_add_or_mul(t) {
        let check = if coeff_isneg(&n) { neg(&n) } else { n.clone() };
        if (hooks.count_ops)(&check) > (hooks.count_ops)(t) {
            return None;
        }
    }
    let q = quotient(&n, power, hooks)?;
    if has(&q, &power.args()[0]) {
        return None;
    }
    if neg(&mul(&[q.clone(), power.clone()])) != *t {
        return None;
    }
    Some(q)
}

/// No match of `a*df + b*f**2 + c*f/x + d/x**2` (`a`, `b` nonzero, all
/// coefficients free of `x`, `f` and `df`) when the collected form has no
/// linear `df` term with such a coefficient, or no `f**2` at all (the term
/// `b*f**2`, `b != 0`, must then equal a term of the collected form).
/// Otherwise the pattern itself is matched on the collected form.
fn riccati_special(c: Option<&Expr>, de: &Outcome<(Expr, Expr)>, fx: &Expr, x: &Expr, df: &Expr) -> Outcome<[Expr; 4]> {
    match de {
        Outcome::NoMatch => Outcome::NoMatch,
        Outcome::Match((_, e)) if e.is_zero() || has(e, x) || has(e, fx) => Outcome::NoMatch,
        Outcome::Match(_) if c.is_some_and(|c| !has(c, &crate::ops::pow(fx, &Expr::int(2)))) => Outcome::NoMatch,
        _ => c.map_or(Outcome::Decline, |c| super::reducible::riccati_special(c, fx, x, df)),
    }
}

/// All first-order matches of `eq_expanded` (the expanded prepared
/// equation) in the function `fx = f(x)` with `df = f(x).diff(x)`.
pub fn first_order(eq_expanded: &Expr, fx: &Expr, x: &Expr, df: &Expr, hooks: &mut Hooks<'_>) -> FirstOrder {
    let c = if blocked(eq_expanded) { None } else { collected(eq_expanded, df, fx, hooks) };
    let de = match &c {
        Some(c) => split_de(c, df, &[], hooks),
        None => Outcome::Decline,
    };
    let riccati = riccati_special(c.as_ref(), &de, fx, x, df);
    let (linear, bernoulli) = {
        let k = coeff(eq_expanded, df);
        if k.is_zero() || blocked(eq_expanded) {
            (Outcome::Decline, Outcome::Decline)
        } else {
            match (hooks.expand)(&div(eq_expanded, &k)).and_then(|e| collected(&e, df, fx, hooks)) {
                Some(eq2) => linear_forms(&eq2, fx, x, df, hooks),
                None => (Outcome::Decline, Outcome::Decline),
            }
        }
    };
    FirstOrder { de, linear, bernoulli, riccati_special: riccati }
}

/// The reference `allhints`, in order.
pub const ALL_HINTS: &[&str] = &[
    "factorable", "nth_algebraic", "separable", "1st_exact", "1st_linear", "Bernoulli",
    "1st_rational_riccati", "Riccati_special_minus2", "1st_homogeneous_coeff_best",
    "1st_homogeneous_coeff_subs_indep_div_dep", "1st_homogeneous_coeff_subs_dep_div_indep",
    "almost_linear", "linear_coefficients", "separable_reduced", "1st_power_series", "lie_group",
    "nth_linear_constant_coeff_homogeneous", "nth_linear_euler_eq_homogeneous",
    "nth_linear_constant_coeff_undetermined_coefficients",
    "nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients",
    "nth_linear_constant_coeff_variation_of_parameters",
    "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters", "Liouville", "2nd_linear_airy",
    "2nd_linear_bessel", "2nd_hypergeometric", "2nd_hypergeometric_Integral", "nth_order_reducible",
    "2nd_power_series_ordinary", "2nd_power_series_regular", "nth_algebraic_Integral",
    "separable_Integral", "1st_exact_Integral", "1st_linear_Integral", "Bernoulli_Integral",
    "1st_homogeneous_coeff_subs_indep_div_dep_Integral",
    "1st_homogeneous_coeff_subs_dep_div_indep_Integral", "almost_linear_Integral",
    "linear_coefficients_Integral", "separable_reduced_Integral",
    "nth_linear_constant_coeff_variation_of_parameters_Integral",
    "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters_Integral", "Liouville_Integral",
    "2nd_nonlinear_autonomous_conserved", "2nd_nonlinear_autonomous_conserved_Integral",
];

/// Solvers whose match also registers `hint + "_Integral"`.
pub const HAS_INTEGRAL: &[&str] = &[
    "nth_algebraic", "separable", "1st_exact", "1st_linear", "Bernoulli",
    "1st_homogeneous_coeff_subs_indep_div_dep", "1st_homogeneous_coeff_subs_dep_div_indep",
    "almost_linear", "linear_coefficients", "separable_reduced",
    "nth_linear_constant_coeff_variation_of_parameters",
    "nth_linear_euler_eq_nonhomogeneous_variation_of_parameters", "Liouville",
    "2nd_hypergeometric", "2nd_nonlinear_autonomous_conserved",
];

/// The matched solver hints (in matching order) with their `_Integral`
/// entries, as `matching_hints` is filled by the reference.
pub fn registered_hints(matched: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    for h in matched {
        out.push(h.clone());
        if HAS_INTEGRAL.contains(&h.as_str()) {
            out.push(format!("{h}_Integral"));
        }
    }
    out
}

/// `[i for i in allhints if i in matching_hints]`.
pub fn ordered_hints(registered: &[String]) -> Vec<String> {
    ALL_HINTS.iter().filter(|h| registered.iter().any(|r| r == *h)).map(|h| (*h).to_string()).collect()
}

/// Convenience for callers mapping an outcome to its parts.
pub fn outcome_parts<T: Clone>(o: &Outcome<T>) -> Option<Option<T>> {
    match o {
        Outcome::Match(v) => Some(Some(v.clone())),
        Outcome::NoMatch => Some(None),
        Outcome::Decline => None,
    }
}

/// The symbols of `e` as `e.atoms(Symbol)` sees them, or `None` when that
/// set may differ from the free symbols (a bound variable).
fn symbol_atoms(e: &Expr) -> Option<usize> {
    let mut seen: Vec<Expr> = Vec::new();
    preorder(e, &mut |x| {
        if x.is_symbol() && !seen.contains(x) {
            seen.push(x.clone());
        }
    });
    let free = crate::traverse::free_symbols(e);
    (free.len() == seen.len() && seen.iter().all(|s| free.contains(s))).then_some(seen.len())
}

/// Whether `e` has an application of the head of `fx` (`e.atoms(f)`).
fn has_head(e: &Expr, fx: &Expr) -> bool {
    let Kind::Func(head, _) = fx.kind() else { return true };
    let mut found = false;
    preorder(e, &mut |x| {
        if let Kind::Func(h, _) = x.kind() {
            found |= h == head;
        }
    });
    found
}

/// `collect(e, fx)`.
fn collected_in(e: &Expr, fx: &Expr) -> Option<Expr> {
    let options = Options { exact: false, discover: false, distribute_order_term: true };
    match collect::run::<()>(e, &[fx.clone()], options, true, None) {
        Ok(Collection::Expression(v)) => Some(v),
        _ => None,
    }
}

/// The reference `match_riccati(eq, f, x)`: `(b0, b1, b2)` of
/// `df = b0 + b1*f + b2*f**2` with rational coefficients in `x`.
pub fn riccati(eq: &Expr, fx: &Expr, x: &Expr, df: &Expr, hooks: &mut Hooks<'_>) -> Outcome<(Expr, Expr, Expr)> {
    if blocked(eq) {
        return Outcome::Decline;
    }
    let Some(eq) = (hooks.expand)(eq).and_then(|e| collected_in(&e, fx)) else { return Outcome::Decline };
    let cf = coeff(&eq, df);
    if cf.is_zero() || !matches!(eq.kind(), Kind::Add(..)) {
        return Outcome::NoMatch;
    }
    let mut parts = Vec::new();
    for t in terms(&eq) {
        match (hooks.cancel)(&div(&t, &cf)) {
            Some(v) => parts.push(v),
            None => return Outcome::Decline,
        }
    }
    let Some(eq) = collected_in(&crate::ops::add(&parts), fx) else { return Outcome::Decline };
    if blocked(&eq) {
        return Outcome::Decline;
    }
    let f2 = crate::ops::pow(fx, &Expr::int(2));
    let b1 = neg(&coeff(&eq, fx));
    let b2 = neg(&coeff(&eq, &f2));
    let rest = sub(&sub(&sub(df, &crate::ops::mul2(&b1, fx)), &crate::ops::mul2(&b2, &f2)), &eq);
    let Some(b0) = (hooks.expand)(&rest) else { return Outcome::Decline };
    let funcs = [b0, b1, b2];
    for b in &funcs {
        if blocked(b) {
            return Outcome::Decline;
        }
        match symbol_atoms(b) {
            None => return Outcome::Decline,
            Some(n) if n > 1 => return Outcome::NoMatch,
            _ => {}
        }
    }
    let [b0, b1, b2] = funcs;
    if has_head(&b0, fx) || b2.is_zero() {
        return Outcome::NoMatch;
    }
    let rational = |b: &Expr| crate::inspection::classify(b, std::slice::from_ref(x), crate::inspection::PolynomialKind::Rational) == Some(true);
    if !(rational(&b0) && rational(&b1) && rational(&b2)) {
        return Outcome::NoMatch;
    }
    Outcome::Match((b0, b1, b2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn setup() -> (Expr, Expr, Expr) {
        let x = Expr::symbol("x");
        let f = parse_str("f(x)").unwrap();
        let df = crate::diff::derivative_unevaluated(&f, vec![(x.clone(), 1)]);
        (x, f, df)
    }

    fn run(eq: &Expr, f: &Expr, x: &Expr, df: &Expr) -> FirstOrder {
        let mut cancel = |e: &Expr| crate::poly::api::cancel(e).ok();
        let mut expand = |e: &Expr| Some(crate::expand::expand(e));
        let mut count = |e: &Expr| BigInt::from(crate::simptools::methods::count_ops(e));
        let mut hooks = Hooks { cancel: &mut cancel, expand: &mut expand, count_ops: &mut count };
        first_order(eq, f, x, df, &mut hooks)
    }

    #[test]
    fn separable_and_linear_decompositions() {
        let (x, f, df) = setup();
        // f' - x*f
        let eq = sub(&df, &mul(&[x.clone(), f.clone()]));
        let r = run(&eq, &f, &x, &df);
        let Outcome::Match((d, e)) = r.de else { panic!("de") };
        assert_eq!(d, neg(&mul(&[x.clone(), f.clone()])));
        assert!(e.is_one());
        let Outcome::Match((p, q)) = r.linear else { panic!("linear") };
        assert_eq!(p, neg(&x));
        assert!(q.is_zero());
        let Outcome::Match((p, q, n)) = r.bernoulli else { panic!("bernoulli") };
        assert_eq!(p, neg(&x));
        assert!(q.is_zero() && n.is_none());
        assert!(matches!(r.riccati_special, Outcome::NoMatch));
    }

    #[test]
    fn bernoulli_power_and_nonlinear_derivative() {
        let (x, f, df) = setup();
        // f' + f/x - x*f**3
        let f3 = crate::ops::pow(&f, &Expr::int(3));
        let eq = crate::ops::add(&[df.clone(), div(&f, &x), neg(&mul(&[x.clone(), f3]))]);
        let r = run(&eq, &f, &x, &df);
        assert!(matches!(r.linear, Outcome::NoMatch));
        let Outcome::Match((p, q, n)) = r.bernoulli else { panic!("bernoulli") };
        assert_eq!(p, div(&Expr::one(), &x));
        assert_eq!(q, x);
        assert_eq!(n, Some(Expr::int(3)));
        // f'**2 - f: the derivative is not linear.
        let eq = sub(&crate::ops::pow(&df, &Expr::int(2)), &f);
        let r = run(&eq, &f, &x, &df);
        assert!(matches!(r.de, Outcome::NoMatch));
    }

    #[test]
    fn unmodelled_content_declines() {
        let (x, f, df) = setup();
        let eq = sub(&df, &mul(&[Expr::constant(Constant::I), f.clone()]));
        let r = run(&eq, &f, &x, &df);
        assert!(matches!(r.de, Outcome::Decline));
        assert!(matches!(r.linear, Outcome::Decline));
        assert!(matches!(r.bernoulli, Outcome::Decline));
    }

    #[test]
    fn hint_order_follows_allhints() {
        let matched = vec!["separable".to_string(), "1st_power_series".to_string(), "factorable".to_string()];
        let ordered = ordered_hints(&registered_hints(&matched));
        assert_eq!(ordered, vec!["factorable", "separable", "1st_power_series", "separable_Integral"]);
    }

    fn riccati_of(eq: &Expr, f: &Expr, x: &Expr, df: &Expr) -> Outcome<(Expr, Expr, Expr)> {
        let mut cancel = |e: &Expr| crate::poly::api::cancel(e).ok();
        let mut expand = |e: &Expr| Some(crate::expand::expand(e));
        let mut count = |e: &Expr| BigInt::from(crate::simptools::methods::count_ops(e));
        let mut hooks = Hooks { cancel: &mut cancel, expand: &mut expand, count_ops: &mut count };
        riccati(eq, f, x, df, &mut hooks)
    }

    #[test]
    fn riccati_coefficients() {
        let (x, f, df) = setup();
        let f2 = crate::ops::pow(&f, &Expr::int(2));
        // 2*f' - 2*x*f**2 - 4/x: b0 = 2/x, b1 = 0, b2 = x
        let eq = crate::ops::add(&[mul(&[Expr::int(2), df.clone()]), neg(&mul(&[Expr::int(2), x.clone(), f2.clone()])), neg(&div(&Expr::int(4), &x))]);
        let Outcome::Match((b0, b1, b2)) = riccati_of(&eq, &f, &x, &df) else { panic!("riccati") };
        assert_eq!(b0, div(&Expr::int(2), &x));
        assert!(b1.is_zero());
        assert_eq!(b2, x);
        // no f**2 term: b2 == 0
        let eq = sub(&df, &mul(&[x.clone(), f.clone()]));
        assert!(matches!(riccati_of(&eq, &f, &x, &df), Outcome::NoMatch));
        // a coefficient that is not rational in x
        let eq = sub(&df, &mul(&[parse_str("exp(x)").unwrap(), f2.clone()]));
        assert!(matches!(riccati_of(&eq, &f, &x, &df), Outcome::NoMatch));
        // two symbols in a coefficient
        let eq = sub(&df, &mul(&[Expr::symbol("a"), x.clone(), f2]));
        assert!(matches!(riccati_of(&eq, &f, &x, &df), Outcome::NoMatch));
    }
}
