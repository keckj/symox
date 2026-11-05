//! Native `2nd_hypergeometric` (2F1) steps of `dsolve`.
//!
//! Ports of SymPy 1.14 `equivalence_hypergeometric`, `equivalence`,
//! `match_2nd_2F1_hypergeometric` and `get_sol_2F1_hypergeometric`
//! (sympy/solvers/ode/hypergeometric.py). The order of operations, the
//! traversal of the invariant's powers and the choice of the singular points
//! and roots follow the reference literally.
//!
//! `cancel`, `factor`, `simplify`, `expand` and `roots` are the public
//! callables supplied by the binding (`Hooks`); a hook that fails, or a
//! step whose reference semantics the port does not model (a float, a
//! symbolic power of `x`, an incomparable root), declines and the caller
//! runs the reference body, which then raises or answers as in SymPy.

use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::diff::diff;
use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, div, mul, mul2, neg, pow, sqrt, sub};
use crate::subs::subs1;
use crate::symbol::Assumptions;
use crate::traverse::{has_structural as has, ordered_args, preorder};

use super::classify::Outcome;

/// The public callables the reference uses; `None` declines.
pub trait Hooks {
    fn cancel(&mut self, e: &Expr) -> Option<Expr>;
    fn factor(&mut self, e: &Expr) -> Option<Expr>;
    fn simplify(&mut self, e: &Expr) -> Option<Expr>;
    fn expand(&mut self, e: &Expr) -> Option<Expr>;
    /// The keys of `roots(e, x)`, in the order of the returned dict.
    fn roots(&mut self, e: &Expr, x: &Expr) -> Option<Vec<Expr>>;
}

/// `match_2nd_2F1_hypergeometric` of the 2F1 equivalence of `(A, B)`.
#[derive(Clone, Debug, PartialEq)]
pub struct Match2F1 {
    pub a: Expr,
    pub b: Expr,
    pub c: Expr,
    pub k: Expr,
    pub mobius: Expr,
}

/// `equivalence_hypergeometric` when it reports the 2F1 type.
#[derive(Clone, Debug, PartialEq)]
pub struct Equivalence {
    pub i0: Expr,
    pub k: Expr,
    pub sing_point: Vec<Expr>,
}

macro_rules! hook {
    ($e:expr) => {
        match $e {
            Some(v) => v,
            None => return Outcome::Decline,
        }
    };
}

fn has_float(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Num(n) if n.is_float()) {
            found = true;
        }
    });
    found
}

/// `e.has(Symbol)`: any symbol atom, dummies and wilds included.
fn has_symbol(e: &Expr) -> bool {
    let mut found = false;
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Sym(_)) {
            found = true;
        }
    });
    found
}

fn ratio(e: &Expr) -> Option<BigRational> {
    match e.kind() {
        Kind::Num(n) => n.to_ratio(),
        _ => None,
    }
}

fn rational(r: BigRational) -> Expr {
    Expr::num(Number::from_ratio(r))
}

/// The nested `_power_counting` of the reference: the exponents of the
/// powers of `x` met by descending into every argument that contains `x`.
fn power_counting(values: &[Expr], x: &Expr, out: &mut Vec<Expr>) {
    for val in values {
        if !has(val, x) {
            continue;
        }
        match val.kind() {
            Kind::Pow(base, exponent) if base == x => push_unique(out, exponent.clone()),
            _ if val == x => push_unique(out, Expr::one()),
            _ => power_counting(&ordered_args(val), x, out),
        }
    }
}

fn push_unique(out: &mut Vec<Expr>, e: Expr) {
    if !out.contains(&e) {
        out.push(e);
    }
}

/// The reference `gcd` of a set of rationals: the gcd of the numerators
/// over the lcm of the denominators; `None` for another exponent.
fn rational_gcd(values: &[Expr]) -> Option<BigRational> {
    let mut numerator = num_bigint::BigInt::zero();
    let mut denominator = num_bigint::BigInt::one();
    for v in values {
        let r = ratio(v)?;
        numerator = numerator.gcd(r.numer());
        denominator = denominator.lcm(r.denom());
    }
    Some(BigRational::new(numerator.abs(), denominator))
}

/// `equivalence(max_num_pow, dem_pow)`: the 2F1 rows of table 1 of Chan and
/// Cheb-Terrab.
fn equivalence(max_num_pow: i64, dem_pow: &[i64]) -> bool {
    let rows: &[&[i64]] = match max_num_pow {
        2 => &[&[2, 2], &[2, 2, 2]],
        1 => &[&[1, 2, 2], &[2, 2, 2], &[1, 2], &[2, 2]],
        0 => &[&[1, 1, 2], &[2, 2], &[1, 2, 2], &[1, 1], &[2], &[1, 2], &[2, 2]],
        _ => &[],
    };
    rows.contains(&dem_pow)
}

fn small_integer(r: &BigRational) -> Option<i64> {
    if !r.is_integer() {
        return None;
    }
    i64::try_from(r.numer()).ok()
}

/// `equivalence_hypergeometric(A, B, f(x))` restricted to its 2F1 answer.
pub fn equivalence_hypergeometric(a_coeff: &Expr, b_coeff: &Expr, x: &Expr, h: &mut dyn Hooks) -> Outcome<Equivalence> {
    if has_float(a_coeff) || has_float(b_coeff) || !x.is_symbol() {
        return Outcome::Decline;
    }
    let quarter = Expr::rational(1, 4);
    // making given equation in normal form
    let normal = add(&[div(&diff(a_coeff, x), &Expr::int(2)), mul2(&pow(a_coeff, &Expr::int(2)), &quarter), neg(b_coeff)]);
    let i1 = hook!(h.cancel(&normal));
    let i1 = hook!(h.factor(&i1));
    // computing shifted invariant(J1) of the equation
    let j1 = hook!(h.cancel(&add(&[mul2(&pow(x, &Expr::int(2)), &i1), quarter.clone()])));
    let j1 = hook!(h.factor(&j1));
    let (num, dem) = crate::simptools::as_numer_denom(&j1);
    let num = crate::simptools::powsimp::powdenest(&hook!(h.expand(&num)), false, false);
    let dem = crate::simptools::powsimp::powdenest(&hook!(h.expand(&dem)), false, false);
    let mut powers = vec![Expr::zero()];
    power_counting(std::slice::from_ref(&num), x, &mut powers);
    power_counting(std::slice::from_ref(&dem), x, &mut powers);
    let Some(k) = rational_gcd(&powers) else {
        return Outcome::Decline;
    };
    if k.is_zero() {
        // J1 is free of x: the reference divides by k**2 = 0 and its I0 is
        // zoo, which is not a rational function of x.
        return Outcome::NoMatch;
    }
    let k = rational(k);
    // computing I0 of the given equation
    let shifted = div(&sub(&div(&j1, &pow(&k, &Expr::int(2))), &quarter), &pow(&pow(x, &k), &Expr::int(2)));
    let i0 = hook!(h.factor(&shifted));
    let i0 = hook!(h.simplify(&i0));
    let i0 = crate::simptools::powsimp::powdenest(&i0, true, false);
    let substituted = subs1(&i0, x, &pow(x, &div(&Expr::one(), &k)));
    let i0 = crate::simptools::powsimp::powdenest(&substituted, true, false);
    let i0 = hook!(h.cancel(&i0));
    let i0 = hook!(h.factor(&i0));
    if crate::inspection::classify(&i0, std::slice::from_ref(x), crate::inspection::PolynomialKind::Rational) != Some(true) {
        return Outcome::NoMatch;
    }
    let (num, dem) = crate::simptools::as_numer_denom(&i0);
    let mut num_powers = vec![Expr::zero()];
    power_counting(std::slice::from_ref(&num), x, &mut num_powers);
    let mut max_num_pow: Option<BigRational> = None;
    for p in &num_powers {
        let Some(r) = ratio(p) else {
            return Outcome::Decline;
        };
        if max_num_pow.as_ref().is_none_or(|m| r > *m) {
            max_num_pow = Some(r);
        }
    }
    let mut sing_point = Vec::new();
    let mut dem_pow = Vec::new();
    // calculating singular point of I0.
    for arg in ordered_args(&dem) {
        if !has(&arg, x) {
            continue;
        }
        // (x - a)**n gives its base, a bare (x - a) itself; both record the
        // exponent of `as_base_exp`.
        let (base, exponent) = crate::decompose::as_base_exp(&arg);
        let base = if matches!(arg.kind(), Kind::Pow(..)) { base } else { arg.clone() };
        let Some(exponent) = ratio(&exponent) else {
            return Outcome::Decline;
        };
        dem_pow.push(exponent);
        let keys = hook!(h.roots(&base, x));
        let Some(first) = keys.into_iter().next() else {
            return Outcome::Decline;
        };
        sing_point.push(first);
    }
    dem_pow.sort();
    let Some(max_num_pow) = max_num_pow.as_ref().and_then(small_integer) else {
        return Outcome::NoMatch;
    };
    let dem_pow: Option<Vec<i64>> = dem_pow.iter().map(small_integer).collect();
    match dem_pow {
        Some(dem_pow) if equivalence(max_num_pow, &dem_pow) => Outcome::Match(Equivalence { i0, k, sing_point }),
        _ => Outcome::NoMatch,
    }
}

/// `min` of the reference over exact real roots; `None` when two roots do
/// not compare (the reference raises).
fn least(values: &[Expr]) -> Option<Expr> {
    // Python's min compares with `<`: a non-real value raises TypeError and
    // an undecided relation cannot be made a bool. Both decline (only real
    // values are compared; a single value is not compared).
    if values.len() > 1 && values.iter().any(|v| crate::assumptions::is_real(v) != Some(true)) {
        return None;
    }
    let mut best: Option<Expr> = None;
    for v in values {
        match &best {
            None => best = Some(v.clone()),
            Some(b) => match crate::assumptions::is_negative(&sub(v, b)) {
                Some(true) => best = Some(v.clone()),
                Some(false) => {}
                None => return None,
            },
        }
    }
    best
}

/// `match_2nd_2F1_hypergeometric(I0, k, sing_point, f(x))`.
pub fn match_2f1(eqv: &Equivalence, x: &Expr, h: &mut dyn Hooks) -> Outcome<Match2F1> {
    let sing_point = &eqv.sing_point;
    if sing_point.iter().any(has_float) {
        return Outcome::Decline;
    }
    let standard = sing_point.len() == 2 && sing_point[0].is_zero() && sing_point[1].is_one();
    let (mob, t) = if !standard {
        if sing_point.len() < 2 {
            // the reference indexes sing_point[1]
            return Outcome::Decline;
        }
        let fresh = |name: &str| Expr::dummy(name, Assumptions::NONE);
        let (alpha, beta, gamma, delta) = (fresh("alpha"), fresh("beta"), fresh("gamma"), fresh("delta"));
        // solving the equations of the mobius transformation
        let beta_value = neg(&mul2(&alpha, &sing_point[0]));
        let delta_value = neg(&mul2(&gamma, &sing_point[1]));
        let gamma_value = if sing_point.len() == 3 {
            div(&add(&[beta_value.clone(), mul2(&sing_point[2], &alpha)]), &sub(&sing_point[2], &sing_point[1]))
        } else {
            alpha.clone()
        };
        let substitute = |e: &Expr| subs1(&subs1(&subs1(e, &beta, &beta_value), &delta, &delta_value), &gamma, &gamma_value);
        let mob = div(&add(&[mul2(&alpha, x), beta.clone()]), &add(&[mul2(&gamma, x), delta.clone()]));
        let mob = hook!(h.cancel(&substitute(&mob)));
        let t = div(&sub(&beta, &mul2(&delta, x)), &sub(&mul2(&gamma, x), &alpha));
        let t = hook!(h.cancel(&substitute(&t)));
        if [&alpha, &beta, &gamma, &delta].iter().any(|w| has(&mob, w) || has(&t, w)) {
            return Outcome::Decline;
        }
        (mob, t)
    } else {
        (x.clone(), x.clone())
    };
    // applying mobius transformation in I to make it into I0.
    let transformed = mul2(&subs1(&eqv.i0, x, &t), &pow(&diff(&t, x), &Expr::int(2)));
    let transformed = hook!(h.factor(&transformed));
    let two = Expr::int(2);
    let i0_dem = mul(&[Expr::int(4), pow(x, &two), pow(&sub(x, &Expr::one()), &two)]);
    let expanded = hook!(h.cancel(&mul2(&transformed, &i0_dem)));
    let expanded = hook!(h.expand(&expanded));
    let x2 = pow(x, &two);
    let options = crate::poly::collect::Options { exact: false, discover: false, distribute_order_term: true };
    let collected = match crate::poly::collect::run::<()>(&expanded, &[x2.clone(), x.clone()], options, false, None) {
        Ok(crate::poly::collect::Collection::Dictionary(v)) => v,
        _ => return Outcome::Decline,
    };
    let coeff = |key: &Expr| collected.iter().find(|(k, _)| k == key).map(|(_, v)| v.clone()).unwrap_or_else(Expr::zero);
    let (l2, l1, l0) = (coeff(&x2), coeff(x), coeff(&Expr::one()));
    let one = Expr::one();
    let mut c = sub(&one, &hook!(h.factor(&sqrt(&add(&[one.clone(), l0.clone()])))));
    if !has_symbol(&c) {
        let wild_c = Expr::dummy("c", Assumptions::NONE);
        let equation = sub(&l0, &mul2(&wild_c, &sub(&wild_c, &two)));
        let keys = hook!(h.roots(&equation, &wild_c));
        c = hook!(least(&keys));
    }
    let s = hook!(h.factor(&sqrt(&add(&[l2, one]))));
    let radicand = add(&[pow(&c, &two), pow(&s, &two), l1, neg(&mul2(&two, &c))]);
    let r = sub(&c, &hook!(h.factor(&sqrt(&radicand))));
    let a = div(&add(&[r.clone(), s.clone()]), &two);
    let b = div(&sub(&r, &s), &two);
    Outcome::Match(Match2F1 {
        a: hook!(h.simplify(&a)),
        b: hook!(h.simplify(&b)),
        c: hook!(h.simplify(&c)),
        k: eqv.k.clone(),
        mobius: mob,
    })
}

fn integral(f: &Expr, x: &Expr) -> Expr {
    crate::binders::held(Func::Integral, f, vec![vec![x.clone()]]).expect("valid integration variable")
}

/// `exp(logcombine(Integral(f, x), force=True))`.
fn exp_integral(f: &Expr, x: &Expr) -> Expr {
    crate::functions::exp(&crate::simp::logcombine(&integral(f, x), true))
}

/// The right-hand side of `get_sol_2F1_hypergeometric` with the numbered
/// constants `c0`, `c1`; `NoMatch` for its `None` (no 2F1 basis for these
/// parameters).
pub fn solution_2f1(m: &Match2F1, a_coeff: &Expr, x: &Expr, c0: &Expr, c1: &Expr, h: &mut dyn Hooks) -> Outcome<Expr> {
    use crate::special::hyper;
    let (a, b, c) = (&m.a, &m.b, &m.c);
    if [a, b, c, &m.k, &m.mobius, a_coeff].iter().any(|e| has_float(e)) {
        return Outcome::Decline;
    }
    let one = Expr::one();
    let two = Expr::int(2);
    let sol = if crate::assumptions::is_integer(c) == Some(false) {
        let first = mul2(c0, &hyper(&[a.clone(), b.clone()], std::slice::from_ref(c), x));
        let shifted = [add(&[a.clone(), neg(c), one.clone()]), add(&[b.clone(), neg(c), one.clone()])];
        let second = mul(&[c1.clone(), hyper(&shifted, &[sub(&two, c)], x), pow(x, &sub(&one, c))]);
        add(&[first, second])
    } else if c.is_one() {
        let h0 = hyper(&[a.clone(), b.clone()], std::slice::from_ref(c), x);
        let expanded = match crate::special::hyperexpand::hyperexpand_options(&h0, false, Some("default"), None) {
            Ok(v) => v,
            Err(_) => return Outcome::Decline,
        };
        let inner = div(&add(&[neg(&mul2(&add(&[a.clone(), b.clone(), one.clone()]), x)), c.clone()]), &sub(&pow(x, &two), x));
        let y2 = mul2(&integral(&div(&crate::functions::exp(&integral(&inner, x)), &pow(&expanded, &two)), x), &h0);
        add(&[mul2(c0, &h0), mul2(c1, &y2)])
    } else if crate::assumptions::is_integer(&sub(&sub(c, a), b)) == Some(false) {
        let z = sub(&one, x);
        let first = mul2(c0, &hyper(&[a.clone(), b.clone()], &[add(&[one.clone(), a.clone(), b.clone(), neg(c)])], &z));
        let cab = sub(&sub(c, a), b);
        let second = mul(&[
            c1.clone(),
            hyper(&[sub(c, a), sub(c, b)], &[add(&[one.clone(), cab.clone()])], &z),
            pow(&z, &cab),
        ]);
        add(&[first, second])
    } else {
        return Outcome::NoMatch;
    };
    // applying transformation in the solution
    let mob = &m.mobius;
    let dtdx = hook!(h.simplify(&div(&one, &diff(mob, x))));
    let quadratic = sub(&pow(x, &two), x);
    let at_mob = |e: &Expr| subs1(e, x, mob);
    let big_b = mul2(&at_mob(&sub(&mul2(&add(&[a.clone(), b.clone(), one.clone()]), x), c)), &dtdx);
    let big_b = hook!(h.factor(&add(&[big_b, mul2(&at_mob(&quadratic), &mul2(&diff(&dtdx, x), &dtdx))])));
    let big_a = hook!(h.factor(&mul2(&at_mob(&quadratic), &pow(&dtdx, &two))));
    let e = exp_integral(&hook!(h.cancel(&div(&big_b, &mul2(&two, &big_a)))), x);
    let xk = pow(x, &m.k);
    let sol = subs1(&subs1(&sol, x, mob), x, &xk);
    let e = subs1(&e, x, &xk);
    let power = pow(x, &div(&sub(&one, &m.k), &two));
    let factor = if crate::assumptions::is_zero(a_coeff) == Some(true) {
        mul2(&e, &power)
    } else {
        let e1 = exp_integral(&div(a_coeff, &two), x);
        mul2(&div(&e, &e1), &power)
    };
    Outcome::Match(mul2(&hook!(h.cancel(&factor)), &sol))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(source: &str) -> Expr {
        crate::parser::parse_str(source).unwrap()
    }

    #[test]
    fn power_counting_descends_into_non_x_bases() {
        let x = p("x");
        let mut out = vec![Expr::zero()];
        power_counting(&[p("3*x**(1/2)*(x+1)**2 + exp(x**3)")], &x, &mut out);
        assert!(out.contains(&p("1/2")) && out.contains(&p("1")) && out.contains(&p("3")));
        assert!(!out.contains(&p("2")));
        assert_eq!(rational_gcd(&out), Some(BigRational::new(1.into(), 2.into())));
    }

    #[test]
    fn rational_gcd_matches_reference() {
        assert_eq!(rational_gcd(&[p("0")]), Some(BigRational::zero()));
        assert_eq!(rational_gcd(&[p("0"), p("-1"), p("2")]), Some(BigRational::one()));
        assert_eq!(rational_gcd(&[p("2/7"), p("4/7")]), Some(BigRational::new(2.into(), 7.into())));
        assert_eq!(rational_gcd(&[p("0"), p("n")]), None);
    }

    #[test]
    fn equivalence_table() {
        assert!(equivalence(2, &[2, 2]));
        assert!(equivalence(0, &[1, 1]));
        assert!(!equivalence(2, &[1, 2]));
        assert!(!equivalence(3, &[2, 2]));
    }

    #[test]
    fn least_requires_comparable_roots() {
        assert_eq!(least(&[p("3"), p("-1")]), Some(p("-1")));
        assert_eq!(least(&[p("1+sqrt(2)"), p("1-sqrt(2)")]), Some(p("1-sqrt(2)")));
        assert_eq!(least(&[p("1+I"), p("1-I")]), None);
        assert_eq!(least(&[p("1+I")]), Some(p("1+I")));
        assert_eq!(least(&[p("n"), p("1")]), None);
    }
}
