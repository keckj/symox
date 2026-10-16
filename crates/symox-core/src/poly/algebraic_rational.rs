//! Exact univariate rational operations over the existing algebraic fields.
//! All inputs are unified before coefficient arithmetic; Euclid and the CRT
//! operate in the genuine quotient field, never on radical spellings.
use super::algebraic::{AlgebraicPoly, AlgebraicScope, NumberField, AN};
use super::convert::{PResult, PolyError};
use super::dense;
use super::int::Coeff;
use super::polyobj::unify_gens;
use crate::expr::Expr;
use crate::ops::{add, div, pow};

pub fn unify(a: &AlgebraicPoly, b: &AlgebraicPoly) -> PResult<(AlgebraicPoly, AlgebraicPoly)> {
    let gens = unify_gens(&a.gens, &b.gens);
    if gens.len() != 1 {
        return Err(PolyError::Unsupported(
            "algebraic rational operations require one polynomial generator".into(),
        ));
    }
    let field = if a.number_field == b.number_field {
        a.number_field.clone()
    } else {
        NumberField::from_generators(&[
            a.number_field.generator.clone(),
            b.number_field.generator.clone(),
        ])?
    };
    let convert = |p: &AlgebraicPoly| {
        if p.number_field == field && p.gens == gens {
            Ok(p.clone())
        } else {
            AlgebraicPoly::from_expr(&p.as_expr(), &gens, field.clone())
        }
    };
    Ok((convert(a)?, convert(b)?))
}
fn dense_poly(template: &AlgebraicPoly, coefficients: &[AN]) -> AlgebraicPoly {
    template.with_rep(dense::to_mpoly(coefficients, 0, 1))
}
fn exact_quotient(a: &[AN], b: &[AN]) -> PResult<Vec<AN>> {
    if b.is_empty() {
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    let (q, r) = dense::divrem(a, b);
    if !r.is_empty() {
        return Err(PolyError::Failed(
            "nonexact algebraic rational quotient".into(),
        ));
    }
    Ok(q)
}
/// `(s,t,h)` with s*a+t*b=h and h monic. Unlike Python's legacy gcdex
/// boundary this native identity also handles a nonzero a with b=0.
pub fn gcdex(
    a: &AlgebraicPoly,
    b: &AlgebraicPoly,
) -> PResult<(AlgebraicPoly, AlgebraicPoly, AlgebraicPoly)> {
    let (a, b) = unify(a, b)?;
    if a.rep.is_zero() && b.rep.is_zero() {
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    let _scope = AlgebraicScope::new(a.number_field.clone());
    let (h, s, t) =
        dense::gcdex_field(&dense::from_mpoly(&a.rep, 0), &dense::from_mpoly(&b.rep, 0));
    Ok((dense_poly(&a, &s), dense_poly(&a, &t), dense_poly(&a, &h)))
}
pub fn invert(a: &AlgebraicPoly, b: &AlgebraicPoly) -> PResult<AlgebraicPoly> {
    let (a, b) = unify(a, b)?;
    if b.rep.is_zero() {
        if a.rep.total_degree() > 0 {
            return Err(PolyError::NotInvertible("zero divisor".into()));
        }
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    let _scope = AlgebraicScope::new(a.number_field.clone());
    let modulus = dense::from_mpoly(&b.rep, 0);
    let (h, s, _) = dense::gcdex_field(&dense::from_mpoly(&a.rep, 0), &modulus);
    if h.len() != 1 {
        return Err(PolyError::NotInvertible("zero divisor".into()));
    }
    Ok(dense_poly(&a, &dense::rem(&s, &modulus)))
}
/// Cancel a monic common gcd, preserving the original scalar units. This is
/// also the pair form of cancel: coefficient=1 and these two exact quotients.
pub fn cancel(a: &AlgebraicPoly, b: &AlgebraicPoly) -> PResult<(AlgebraicPoly, AlgebraicPoly)> {
    let (a, b) = unify(a, b)?;
    if b.rep.is_zero() {
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    let _scope = AlgebraicScope::new(a.number_field.clone());
    if a.rep.is_zero() {
        return Ok((a.clone(), dense_poly(&a, &[AN::one()])));
    }
    let numerator = dense::from_mpoly(&a.rep, 0);
    let denominator = dense::from_mpoly(&b.rep, 0);
    let gcd = dense::gcd_field(&numerator, &denominator);
    Ok((
        dense_poly(&a, &exact_quotient(&numerator, &gcd)?),
        dense_poly(&a, &exact_quotient(&denominator, &gcd)?),
    ))
}
#[derive(Clone, Debug)]
pub struct FractionTerm {
    pub numerator: AlgebraicPoly,
    pub factor: AlgebraicPoly,
    pub power: u32,
}
#[derive(Clone, Debug)]
pub struct PartialFractions {
    pub polynomial: AlgebraicPoly,
    pub terms: Vec<FractionTerm>,
}
impl PartialFractions {
    pub fn as_expr(&self) -> Expr {
        let mut terms = vec![self.polynomial.as_expr()];
        terms.extend(self.terms.iter().map(|t| {
            div(
                &t.numerator.as_expr(),
                &pow(&t.factor.as_expr(), &Expr::int(t.power as i64)),
            )
        }));
        add(&terms)
    }
    /// Independent exact certificate, after multiplication by the original
    /// denominator. Each returned denominator must divide it exactly.
    pub fn verify(&self, numerator: &AlgebraicPoly, denominator: &AlgebraicPoly) -> bool {
        let Ok((numerator, denominator)) = unify(numerator, denominator) else {
            return false;
        };
        if denominator.rep.is_zero() {
            return false;
        }
        let Ok(polynomial) = numerator.convert(&self.polynomial.as_expr()) else {
            return false;
        };
        let mut out = polynomial.mul(&denominator);
        for term in &self.terms {
            let Ok(factor) = numerator.convert(&term.factor.as_expr()) else {
                return false;
            };
            let Ok(value) = numerator.convert(&term.numerator.as_expr()) else {
                return false;
            };
            let Ok((q, r)) = denominator.div(&factor.pow(term.power)) else {
                return false;
            };
            if !r.rep.is_zero() {
                return false;
            }
            out = out.add(&value.mul(&q));
        }
        out.rep == numerator.rep
    }
}
/// Partial fractions with CRT modulo each irreducible factor power, followed
/// by its factor-adic expansion. Repeated factors need no derivative formula.
pub fn apart(a: &AlgebraicPoly, b: &AlgebraicPoly) -> PResult<PartialFractions> {
    let (a, b) = cancel(a, b)?;
    let _scope = AlgebraicScope::new(a.number_field.clone());
    let numerator = dense::from_mpoly(&a.rep, 0);
    let denominator = dense::from_mpoly(&b.rep, 0);
    let (polynomial, rem) = dense::divrem(&numerator, &denominator);
    let mut out = PartialFractions {
        polynomial: dense_poly(&a, &polynomial),
        terms: Vec::new(),
    };
    if rem.is_empty() {
        return Ok(out);
    }
    let (_, factors) = b.factor_list()?;
    for (factor, k) in factors {
        let fd = dense::from_mpoly(&factor.rep, 0);
        if fd.len() < 2 {
            return Err(PolyError::Failed(
                "constant factor in algebraic partial fractions".into(),
            ));
        }
        let mut power = vec![AN::one()];
        for _ in 0..k {
            power = dense::mul(&power, &fd);
        }
        let complement = exact_quotient(&denominator, &power)?;
        let (gcd, inverse, _) = dense::gcdex_field(&dense::rem(&complement, &power), &power);
        if gcd.len() != 1 || !gcd[0].is_one() {
            return Err(PolyError::Failed(
                "noncoprime algebraic partial-fraction factors".into(),
            ));
        }
        let mut digits = dense::rem(&dense::mul(&rem, &inverse), &power);
        for i in 0..k {
            if digits.is_empty() {
                break;
            }
            let (next, digit) = dense::divrem(&digits, &fd);
            if !digit.is_empty() {
                out.terms.push(FractionTerm {
                    numerator: dense_poly(&a, &digit),
                    factor: factor.clone(),
                    power: k - i,
                });
            }
            digits = next;
        }
        if !digits.is_empty() {
            return Err(PolyError::Failed(
                "incomplete algebraic factor-adic expansion".into(),
            ));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    const FIELDS: [&str; 7] = [
        "sqrt(2)",
        "I",
        "2^(1/3)",
        "sqrt(2)+sqrt(3)",
        "2^(1/3)+sqrt(2)",
        "2^(1/3)+I*sqrt(3)",
        "sqrt(2)+sqrt(3)+sqrt(5)",
    ];
    #[test]
    fn exact_bezout_inverse_and_cancel_certificates_in_every_field_family() {
        let x = Expr::symbol("x");
        for source in FIELDS {
            let generator = parse_str(source).unwrap();
            let field = NumberField::from_generator(&generator).unwrap();
            let p = |e: &Expr| {
                AlgebraicPoly::from_expr(e, std::slice::from_ref(&x), field.clone()).unwrap()
            };
            let a = p(&crate::ops::add2(&x, &generator));
            let b = p(&crate::ops::sub(
                &crate::ops::sub(&x.powi(2), &generator.powi(2)),
                &Expr::one(),
            ));
            let (s, t, h) = gcdex(&a, &b).unwrap();
            assert_eq!(s.mul(&a).add(&t.mul(&b)).rep, h.rep, "{source}");
            let inverse = invert(&a, &b).unwrap();
            assert!(a.mul(&inverse).div(&b).unwrap().1.rep.is_one(), "{source}");
            let common = p(&crate::ops::sub(&x, &generator));
            let numerator = common.mul(&a);
            let denominator = common.mul(&b);
            let (n, d) = cancel(&numerator, &denominator).unwrap();
            assert_eq!(n.mul(&denominator).rep, d.mul(&numerator).rep, "{source}");
            assert!(n.gcd(&d).rep.is_one(), "{source}");
            assert!(matches!(
                invert(&numerator, &denominator),
                Err(PolyError::NotInvertible(_))
            ));
            let zero = p(&Expr::zero());
            assert!(matches!(
                invert(&zero, &zero),
                Err(PolyError::Unsupported(_))
            ));
            assert!(matches!(cancel(&a, &zero), Err(PolyError::Unsupported(_))));
        }
    }
    #[test]
    fn crt_and_factor_adic_digits_have_independent_exact_certificates() {
        let x = Expr::symbol("x");
        for source in FIELDS {
            let generator = parse_str(source).unwrap();
            let field = NumberField::from_generator(&generator).unwrap();
            let p = |e: &Expr| {
                AlgebraicPoly::from_expr(e, std::slice::from_ref(&x), field.clone()).unwrap()
            };
            let a = p(&crate::ops::add2(&x, &generator));
            let linear = p(&crate::ops::sub(&x, &generator));
            let other = p(&crate::ops::add2(&x, &Expr::one()));
            let quadratic = p(&crate::ops::add2(&x.powi(2), &Expr::one()));
            for b in [
                linear.mul(&other),
                linear.pow(3).mul(&other.pow(2)),
                linear.mul(&quadratic),
            ] {
                let decomposition = apart(&a, &b).unwrap();
                assert!(decomposition.verify(&a, &b), "{source}");
                assert!(!decomposition.terms.is_empty());
                for term in &decomposition.terms {
                    assert!(term.numerator.rep.degree(0) < term.factor.rep.degree(0));
                    assert!(term.power > 0);
                }
            }
        }
    }
    #[test]
    fn incompatible_fields_and_generators_fail_before_arithmetic() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let field = |s| NumberField::from_generator(&parse_str(s).unwrap()).unwrap();
        let p = AlgebraicPoly::from_expr(&x, std::slice::from_ref(&x), field("2^(1/3)")).unwrap();
        let q = AlgebraicPoly::from_expr(&x, std::slice::from_ref(&x), field("3^(1/3)")).unwrap();
        assert!(gcdex(&p, &q).is_err());
        assert!(invert(&p, &q).is_err());
        assert!(cancel(&p, &q).is_err());
        assert!(apart(&p, &q).is_err());
        let q = AlgebraicPoly::from_expr(&y, std::slice::from_ref(&y), field("sqrt(2)")).unwrap();
        assert!(gcdex(&p, &q).is_err());
        let a: super::super::fraction::Polynomial = p.into();
        let b: super::super::fraction::Polynomial = q.into();
        assert!(a.checked_algebraic_pair(&b).is_ok());
        let c: super::super::fraction::Polynomial =
            super::super::polyobj::Poly::from_expr(&x, std::slice::from_ref(&x))
                .unwrap()
                .into();
        assert!(a.checked_algebraic_pair(&c).unwrap().is_some());
    }
}
