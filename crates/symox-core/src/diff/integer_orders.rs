//! Closed forms for exact nonnegative *literal* integer derivative orders.
//!
//! This is an optional certificate, never a replacement for the old fallback.
//! No input order is narrowed to u32. Unsupported variables/functions,
//! nonfinite or inexact coefficients and nonliteral orders return None.
//!
//! Proofs: differentiating sin/cos cycles after four steps; exp is a fixed
//! point; sinh/cosh cycle after two steps (including canonical images of
//! sin(I*x)/cos(I*x)). An affine argument contributes its slope at every step.
//! For a monomial/power of an affine polynomial, induction gives
//! D^n(u^m)=(m)_n*(u')^n*u^(m-n) for 0<=n<=m, and zero for n>m.
//! A structural polynomial degree bound d proves D^n(P)=0 for n>d;
//! D^d(P)=d!*[x^d]P even when the computed leading coefficient cancels.
//! Bounds and exponents are BigInt; powers/falling factorials may remain exact
//! unevaluated expressions instead of allocating an astronomically large int.
use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive, Zero};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, mul, mul2, neg, pow};

const MAX_MATERIALIZED_BITS: usize = 1 << 20;

/// Literal means exact numeric representation, never Float(0)/Float(1),
/// a Boolean or an integer-assumed Symbol.
pub(super) fn literal(e: &Expr) -> Option<BigInt> {
    match e.kind() {
        Kind::Num(Number::Int(n)) if *n >= 0 => Some((*n).into()),
        Kind::Num(Number::BigInt(n)) if !n.is_negative() => Some(n.clone()),
        _ => None,
    }
}
fn integer(n: &BigInt) -> Expr {
    Expr::num(Number::from_bigint(n.clone()))
}
fn ordinary_variable(x: &Expr) -> bool {
    x.as_sym().is_some_and(|s| {
        x.host_handle().is_none()
            && s.scalar_identity().is_none()
            && crate::decompose::commutativity(x) == Some(true)
    })
}
/// Optimization budget only: exhausting it declines the certificate and
/// preserves the established computation, never a mathematical zero.
fn bounded_input(e: &Expr) -> bool {
    let mut pending = vec![(e.clone(), 0usize)];
    let mut seen = FxHashSet::default();
    while let Some((node, depth)) = pending.pop() {
        if depth > 256 {
            return false;
        }
        if !seen.insert(node.dag_id()) {
            continue;
        }
        if seen.len() > 4096 {
            return false;
        }
        node.for_each_owned(&mut |child| pending.push((child.clone(), depth + 1)));
    }
    true
}
fn finite_coefficient(e: &Expr) -> bool {
    !e.has_host()
        && !e.has_nonfinite()
        && e.float_info().max_precision.is_none()
        && crate::functions::is_numeric_value(e)
        && crate::decompose::commutativity(e) == Some(true)
        && crate::assumptions::exact_queries(|| crate::assumptions::fact(e, "finite") == Some(true))
}
fn held_power(base: &Expr, exponent: &BigInt) -> Expr {
    if exponent.is_zero() {
        return Expr::one();
    }
    if exponent.is_one() {
        return base.clone();
    }
    if base.is_zero() {
        return Expr::zero();
    }
    if base.is_one() {
        return Expr::one();
    }
    if base.is_minus_one() {
        return if (exponent % 2u8).is_zero() {
            Expr::one()
        } else {
            Expr::minus_one()
        };
    }
    if base.is_constant(Constant::E) {
        return pow(base, &integer(exponent));
    }
    if base.is_constant(Constant::I) || *base == neg(&Expr::i()) {
        let mut phase = (exponent % 4u8).to_u8().unwrap();
        if !base.is_constant(Constant::I) { phase=(4-phase)%4; }
        return match phase {
            0 => Expr::one(),
            1 => Expr::i(),
            2 => Expr::minus_one(),
            _ => neg(&Expr::i()),
        };
    }
    if let Some((a, b)) = base.as_number().and_then(Number::numer_denom) {
        let bits = a.bits().max(b.bits()).max(1);
        if exponent * BigInt::from(bits) <= BigInt::from(MAX_MATERIALIZED_BITS) {
            return pow(base, &integer(exponent));
        }
        if a.is_negative() {
            let positive = Expr::num(Number::rational(-a, b).unwrap());
            let value = crate::ops::deferred::hold(&positive,&integer(exponent));
            return if (exponent % 2u8).is_zero() {
                value
            } else {
                neg(&value)
            };
        }
    } else if exponent <= &BigInt::from(64u8) {
        return pow(base, &integer(exponent));
    }
    crate::ops::deferred::hold(base,&integer(exponent))
}
fn factorial(n: &BigInt) -> Expr {
    // log2(n!) <= n*ceil(log2(n+1)): this storage estimate never
    // constrains the public factorial evaluator or explicit reconstruction.
    if n * BigInt::from(n.bits().max(1)) <= BigInt::from(MAX_MATERIALIZED_BITS) {
        return crate::special::factorial(&integer(n));
    }
    Expr::from_unevaluated_kind(Kind::Func(
        Func::Special(crate::special::SpecialFn::Factorial),
        vec![integer(n)],
    ))
}
fn falling(m: &BigInt, n: &BigInt) -> Expr {
    if n.is_zero() {
        return Expr::one();
    }
    if n.is_one() {
        return integer(m);
    }
    if m - n <= BigInt::one() {
        return factorial(m);
    }
    if n * BigInt::from(m.bits().max(1)) <= BigInt::from(MAX_MATERIALIZED_BITS) {
        let count = n.to_usize().unwrap();
        let mut value = BigInt::one();
        for k in 0..count {
            value *= m - BigInt::from(k);
        }
        return integer(&value);
    }
    Expr::from_unevaluated_kind(Kind::Func(
        Func::Special(crate::special::SpecialFn::FallingFactorial),
        vec![integer(m), integer(n)],
    ))
}

#[derive(Clone)]
struct Polynomial {
    degree: BigInt,
    leading: Expr,
}
struct Kernel<'a> {
    x: &'a Expr,
    order: &'a BigInt,
    polynomials: FxHashMap<usize, (Expr, Option<Polynomial>)>,
    degrees: FxHashMap<usize, (Expr, Option<BigInt>)>,
}
impl Kernel<'_> {
    fn independent(&self, e: &Expr) -> bool {
        !crate::traverse::free_symbols(e).contains(self.x)
    }
    /// Prove a degree upper bound without constructing any coefficient power.
    /// This is enough for zero beyond the degree, even for huge factored input.
    fn degree(&mut self, e: &Expr) -> Option<BigInt> {
        if let Some((_, degree)) = self.degrees.get(&e.dag_id()) {
            return degree.clone();
        }
        let result = self.degree_inner(e);
        self.degrees.insert(e.dag_id(), (e.clone(), result.clone()));
        result
    }
    fn degree_inner(&mut self, e: &Expr) -> Option<BigInt> {
        if e == self.x {
            return Some(BigInt::one());
        }
        if self.independent(e) {
            return finite_coefficient(e).then(BigInt::zero);
        }
        match e.kind() {
            Kind::Add(_, terms) => {
                let mut degree = BigInt::zero();
                for (term, _) in terms {
                    degree = degree.max(self.degree(term)?);
                }
                Some(degree)
            }
            Kind::Mul(_, factors) => {
                let mut degree = BigInt::zero();
                for (base, exponent) in factors {
                    degree += self.degree(base)? * literal(exponent)?;
                }
                Some(degree)
            }
            Kind::Pow(base, exponent) => Some(self.degree(base)? * literal(exponent)?),
            _ => None,
        }
    }
    fn polynomial(&mut self, e: &Expr) -> Option<Polynomial> {
        if let Some((_, info)) = self.polynomials.get(&e.dag_id()) {
            return info.clone();
        }
        let result = self.polynomial_inner(e);
        self.polynomials
            .insert(e.dag_id(), (e.clone(), result.clone()));
        result
    }
    fn polynomial_inner(&mut self, e: &Expr) -> Option<Polynomial> {
        if e == self.x {
            return Some(Polynomial {
                degree: BigInt::one(),
                leading: Expr::one(),
            });
        }
        if self.independent(e) {
            return finite_coefficient(e).then(|| Polynomial {
                degree: BigInt::zero(),
                leading: e.clone(),
            });
        }
        match e.kind() {
            Kind::Add(c, terms) => {
                let mut items = vec![(
                    Polynomial {
                        degree: BigInt::zero(),
                        leading: Expr::num(c.clone()),
                    },
                    Expr::one(),
                )];
                for (term, coefficient) in terms {
                    items.push((self.polynomial(term)?, Expr::num(coefficient.clone())));
                }
                let degree = items.iter().map(|(p, _)| p.degree.clone()).max()?;
                let leading = add(&items
                    .iter()
                    .filter(|(p, _)| p.degree == degree)
                    .map(|(p, c)| mul2(c, &p.leading))
                    .collect::<Vec<_>>());
                Some(Polynomial { degree, leading })
            }
            Kind::Mul(c, factors) => {
                let mut degree = BigInt::zero();
                let mut leading = vec![Expr::num(c.clone())];
                for (base, exponent) in factors {
                    let exponent = literal(exponent)?;
                    let part = self.polynomial(base)?;
                    degree += part.degree * &exponent;
                    leading.push(held_power(&part.leading, &exponent));
                }
                Some(Polynomial {
                    degree,
                    leading: mul(&leading),
                })
            }
            Kind::Pow(base, exponent) => {
                let exponent = literal(exponent)?;
                let part = self.polynomial(base)?;
                Some(Polynomial {
                    degree: part.degree * &exponent,
                    leading: held_power(&part.leading, &exponent),
                })
            }
            _ => None,
        }
    }
    fn affine(&mut self, e: &Expr) -> Option<Expr> {
        if self.degree(e)? > BigInt::one() {
            return None;
        }
        let info = self.polynomial(e)?;
        if info.degree > BigInt::one() {
            None
        } else {
            Some(if info.degree.is_zero() {
                Expr::zero()
            } else {
                info.leading
            })
        }
    }
    fn value(&mut self, e: &Expr) -> Option<Expr> {
        if self.order.is_zero() {
            return Some(e.clone());
        }
        if self.independent(e) {
            return finite_coefficient(e).then(Expr::zero);
        }
        if let Kind::Func(
            function @ (Func::Sin | Func::Cos | Func::Sinh | Func::Cosh | Func::Exp),
            args,
        ) = e.kind()
        {
            if args.len() != 1 {
                return None;
            }
            let slope = self.affine(&args[0])?;
            let factor = held_power(&slope, self.order);
            let phase = if matches!(function, Func::Exp) {
                e.clone()
            } else if matches!(function, Func::Sinh | Func::Cosh) {
                if (self.order % 2u8).is_zero() {
                    e.clone()
                } else if matches!(function, Func::Sinh) {
                    crate::functions::cosh(&args[0])
                } else {
                    crate::functions::sinh(&args[0])
                }
            } else {
                let sine = matches!(function, Func::Sin);
                let shift = (self.order % 4u8).to_u8().unwrap();
                match (sine, shift) {
                    (true, 0) | (false, 3) => crate::functions::sin(&args[0]),
                    (true, 1) | (false, 0) => crate::functions::cos(&args[0]),
                    (true, 2) | (false, 1) => neg(&crate::functions::sin(&args[0])),
                    _ => neg(&crate::functions::cos(&args[0])),
                }
            };
            return Some(mul2(&factor, &phase));
        }
        let degree = self.degree(e);
        if degree.as_ref().is_some_and(|degree| self.order > degree) {
            return Some(Expr::zero());
        }
        if let Some(poly) = degree.and_then(|_| self.polynomial(e)) {
            if self.order > &poly.degree {
                return Some(Expr::zero());
            }
            if self.order == &poly.degree {
                return Some(mul2(&factorial(self.order), &poly.leading));
            }
            if e == self.x {
                return Some(Expr::one());
            }
            if let Kind::Pow(base, exponent) = e.kind() {
                let exponent = literal(exponent)?;
                if self.order <= &exponent {
                    if let Some(slope) = self.affine(base) {
                        let rest = &exponent - self.order;
                        return Some(mul(&[
                            falling(&exponent, self.order),
                            held_power(&slope, self.order),
                            held_power(base, &rest),
                        ]));
                    }
                }
            }
        }
        match e.kind() {
            Kind::Add(_, terms) => {
                let mut values = Vec::with_capacity(terms.len());
                for (term, coefficient) in terms {
                    values.push(mul2(&Expr::num(coefficient.clone()), &self.value(term)?));
                }
                Some(add(&values))
            }
            Kind::Mul(coefficient, factors) => {
                let mut constants = vec![Expr::num(coefficient.clone())];
                let mut dependent = None;
                for (base, exponent) in factors {
                    let factor = if exponent.is_one() {
                        base.clone()
                    } else if let Some(power) = literal(exponent) {
                        held_power(base, &power)
                    } else {
                        Expr::from_unevaluated_kind(Kind::Pow(base.clone(), exponent.clone()))
                    };
                    if self.independent(&factor) {
                        if !finite_coefficient(&factor) {
                            return None;
                        }
                        constants.push(factor);
                    } else {
                        if dependent.is_some() {
                            return None;
                        }
                        dependent = Some(factor);
                    }
                }
                let derivative = self.value(&dependent?)?;
                constants.push(derivative);
                Some(mul(&constants))
            }
            _ => None,
        }
    }
}
/// None means no certificate; callers retain the complete old computation.
pub(super) fn evaluate(e: &Expr, x: &Expr, order: &BigInt) -> Option<Expr> {
    if order.is_negative() {
        return None;
    }
    if order.is_zero() {
        return Some(e.clone());
    }
    if !ordinary_variable(x)
        || e.has_host()
        || e.has_nonfinite()
        || e.has_matrix()
        || !bounded_input(e)
        || e.float_info().max_precision.is_some()
    {
        return None;
    }
    crate::assumptions::exact_queries(|| {
        if crate::assumptions::fact(x, "finite") == Some(false) {
            return None;
        }
        Kernel {
            x,
            order,
            polynomials: FxHashMap::default(),
            degrees: FxHashMap::default(),
        }
        .value(e)
    })
}
pub(super) fn groups(body: &Expr, variables: &[(Expr, Expr)]) -> Option<Expr> {
    let mut value = body.clone();
    for (variable, order) in variables {
        value = evaluate(&value, variable, &literal(order)?)?;
    }
    Some(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(s: &str) -> Expr {
        crate::parser::parse_str(s).unwrap()
    }
    #[test]
    fn periodic_orders_are_not_narrowed_to_a_machine_integer() {
        let x = Expr::symbol("x");
        let base = (BigInt::one() << 160usize) + BigInt::from(u32::MAX);
        for shift in 0..4u8 {
            let n = &base + shift;
            let r = (&n % 4u8).to_u8().unwrap();
            let sine = [p("sin(x)"), p("cos(x)"), p("-sin(x)"), p("-cos(x)")];
            let cosine = [p("cos(x)"), p("-sin(x)"), p("-cos(x)"), p("sin(x)")];
            assert_eq!(
                evaluate(&p("sin(x)"), &x, &n),
                Some(sine[r as usize].clone())
            );
            assert_eq!(
                evaluate(&p("cos(x)"), &x, &n),
                Some(cosine[r as usize].clone())
            );
            assert_eq!(evaluate(&p("exp(x)"), &x, &n), Some(p("exp(x)")));
            assert_eq!(
                evaluate(&p("sinh(x)"), &x, &n),
                Some(p(if r % 2 == 0 { "sinh(x)" } else { "cosh(x)" }))
            );
            assert_eq!(
                evaluate(&p("cosh(x)"), &x, &n),
                Some(p(if r % 2 == 0 { "cosh(x)" } else { "sinh(x)" }))
            );
        }
    }
    #[test]
    fn literal_classifier_keeps_floats_booleans_and_parameters_distinct() {
        assert_eq!(literal(&Expr::int(7)), Some(7.into()));
        for e in [
            Expr::float(0.0),
            Expr::float(1.0),
            Expr::boolean(false),
            Expr::boolean(true),
            Expr::symbol("n"),
            Expr::minus_one(),
            Expr::rational(1, 2),
        ] {
            assert!(literal(&e).is_none());
        }
    }
    #[test]
    fn closed_forms_match_independent_repeated_first_derivatives() {
        let x = Expr::symbol("x");
        for body in [
            "sin(x)",
            "cos(-2*x+1)",
            "exp(3*x-2)",
            "sinh(2*x+1)",
            "cosh(-2*x+1)",
            "x**12",
            "(x+1)**7",
            "3*x**6+5*x**3-2",
            "sin(2*x)+x**5",
        ] {
            let original = p(body);
            let mut repeated = original.clone();
            for n in 0..20i64 {
                let computed = evaluate(&original, &x, &n.into()).unwrap();
                assert_eq!(
                    crate::expand::expand(&computed),
                    crate::expand::expand(&repeated),
                    "{body}, n={n}"
                );
                repeated = super::super::diff(&repeated, &x);
            }
        }
    }
    #[test]
    fn sparse_degree_bound_and_large_affine_power_do_not_expand() {
        let x = Expr::symbol("x");
        let m = BigInt::one() << 130usize;
        let body = pow(&p("x+1"), &integer(&m));
        assert_eq!(evaluate(&body, &x, &(&m + 1u8)), Some(Expr::zero()));
        let top = evaluate(&body, &x, &m).unwrap();
        assert!(
            matches!(top.kind(),Kind::Func(Func::Special(crate::special::SpecialFn::Factorial),a) if a==&vec![integer(&m)])
        );
        let next = evaluate(&body, &x, &(&m - 1u8)).unwrap();
        assert!(crate::traverse::has_structural(&next, &top));
        assert!(crate::traverse::node_count(&next) < 20);
    }
    #[test]
    fn undecidable_inexact_and_nonfinite_coefficients_have_no_certificate() {
        let x = Expr::symbol("x");
        let n = BigInt::one() << 100usize;
        for body in [
            "f(1)*x",
            "a*x**2",
            "sin(x**2)",
            "exp(x**2)",
            "x**(-1)",
            "oo*x**2",
            "zoo*x",
            "nan",
        ] {
            assert!(evaluate(&p(body), &x, &n).is_none(), "{body}");
        }
        let inexact = mul2(&Expr::float(2.0), &x.powi(2));
        assert!(evaluate(&inexact, &x, &n).is_none());
    }
    #[test]
    fn exhausted_optimization_budget_declines_without_a_zero() {
        let x = Expr::symbol("x");
        let mut body = x.clone();
        for _ in 0..260 {
            body = Expr::from_unevaluated_kind(Kind::Func(Func::Sin, vec![body]));
        }
        assert!(!bounded_input(&body));
        assert!(evaluate(&body, &x, &BigInt::from(2u8)).is_none());
        assert_eq!(evaluate(&body, &x, &BigInt::zero()), Some(body));
    }
    #[test]
    fn huge_coefficient_is_held_exactly_instead_of_materialized() {
        let m = BigInt::one() << 180usize;
        let n = BigInt::one() << 120usize;
        let coefficient = falling(&m, &n);
        assert!(
            matches!(coefficient.kind(),Kind::Func(Func::Special(crate::special::SpecialFn::FallingFactorial),a) if a==&vec![integer(&m),integer(&n)])
        );
        assert_eq!(falling(&8.into(), &5.into()), Expr::int(6720));
        let power = held_power(&Expr::int(2), &n);
        assert!(
            matches!(power.kind(),Kind::Pow(base,exponent) if base==&Expr::int(2)&&exponent==&integer(&n))
        );
    }
}
