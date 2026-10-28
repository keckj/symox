//! Checked polynomial antiderivatives with truthful coefficient domains.
use super::algebraic::{AlgebraicScope, AN};
use super::convert::{PResult, PolyError};
use super::fraction::Polynomial;
use super::int::{Coeff, Field, Rat};
use super::mpoly::MPoly;
use super::polyobj::Domain;
use super::ratfunc::{Scope, RF};

impl Polynomial {
    /// SymPy promotes the entire coefficient ring before integrating, even
    /// when the result happens to have integral coefficients or is zero.
    pub fn integration_domain(&self, auto: bool) -> Self {
        let mut p = self.clone();
        if auto {
            if let Some(e)=p.expression.as_mut(){
                if let Some(g)=e.ground.as_ref().filter(|g|g.inexact.is_some()){
                    e.ground=Some(g.to_field());return p;
                }
            }
            if let Some(field) = &mut p.field {
                field.promote_field();
            } else if let Some(a) = p.algebraic.as_ref().filter(|a| a.gaussian == Some(true)) {
                // ZZ_I integrates over its field QQ_I.
                let mut a = a.clone();
                a.gaussian = Some(false);
                p = a.into();
            } else if p.algebraic.is_none() && p.modulus.is_none() && p.base.domain == Domain::ZZ {
                p.base.domain = Domain::QQ;
            }
        }
        p
    }
    pub fn integrate_checked(&self, variable: usize, auto: bool) -> PResult<Self> {
        if variable >= self.true_gens().len() {
            return Err(PolyError::Generator(
                "integration generator is out of range".into(),
            ));
        }
        if let Some(e) = &self.expression {
            return Ok(e.integrate_with_options(1, variable, auto)?.into());
        }
        let p = self.integration_domain(auto);
        let degree = |powers: &[u32]| {
            powers[variable].checked_add(1).ok_or_else(|| {
                PolyError::Unsupported("polynomial integration exponent exceeds u32".into())
            })
        };
        if let Some(a) = &p.algebraic {
            let _scope = AlgebraicScope::new(a.number_field.clone());
            let terms = a
                .rep
                .iter()
                .map(|(powers, c)| {
                    let k = degree(powers)?;
                    let mut powers = powers.to_vec();
                    powers[variable] = k;
                    let q = c.div(&AN::from_i64(k as i64));
                    if a.gaussian == Some(true) && !super::gaussian::is_integral(&q) {
                        let show = |z: &AN| {
                            let (x, y) = super::gaussian::parts(z);
                            format!("({} + {}*I)", crate::printer::to_str(&crate::expr::Expr::num(x.to_number())), crate::printer::to_str(&crate::expr::Expr::num(y.to_number())))
                        };
                        return Err(PolyError::ExactQuotientFailed(format!(
                            "{} does not divide {} in ZZ_I",
                            show(&AN::from_i64(k as i64)),
                            show(c)
                        )));
                    }
                    Ok((powers, q))
                })
                .collect::<PResult<Vec<_>>>()?;
            return Ok(a.with_rep(MPoly::from_terms(a.gens.len(), terms)).into());
        }
        if let Some(a) = &p.field {
            let _scope = Scope::new(a.parameters.len());
            let terms = a
                .rep
                .iter()
                .map(|(powers, c)| {
                    let k = degree(powers)?;
                    let mut powers = powers.to_vec();
                    powers[variable] = k;
                    Ok((powers, c.div(&RF::from_i64(k as i64))))
                })
                .collect::<PResult<Vec<_>>>()?;
            let result = a.with_rep(MPoly::from_terms(a.gens.len(), terms));
            if result.display.ends_with(']') && !result.belongs_to_ring() {
                return Err(PolyError::ExactQuotientFailed("integration coefficient does not divide exactly in the selected coefficient ring".into()));
            }
            return Ok(result.into());
        }
        let terms = p
            .base
            .rep
            .iter()
            .map(|(powers, c)| {
                let k = degree(powers)?;
                if p.modulus.is_some_and(|prime| (k as u64) % prime == 0) {
                    return Err(PolyError::NotInvertible(
                        "integration denominator is zero in the coefficient field".into(),
                    ));
                }
                let mut powers = powers.to_vec();
                powers[variable] = k;
                Ok((powers, c.div(&Rat::from_i64(k as i64))))
            })
            .collect::<PResult<Vec<_>>>()?;
        let result = p.base.with_rep(MPoly::from_terms(p.base.gens.len(), terms));
        if result.domain == Domain::ZZ && result.modulus.is_none() && !result.rep.is_integral() {
            return Err(PolyError::ExactQuotientFailed(
                "integration coefficient does not divide exactly in ZZ".into(),
            ));
        }
        Ok(result.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        algebraic::{AlgebraicPoly, NumberField},
        domains::ExactDomain,
        fraction::FractionPoly,
        polyobj::Poly,
    };
    use crate::{parser::parse_str, Expr};
    #[test]
    fn domains_and_derivative_certificates_survive_exact_integration() {
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        for base in ["ZZ", "QQ"] {
            for ring in [true, false] {
                let descriptor =
                    ExactDomain::functions(base == "ZZ", ring, vec![a.clone()], "grevlex").unwrap();
                let mut value = FractionPoly::from_expr(
                    &parse_str("a*x").unwrap(),
                    std::slice::from_ref(&x),
                    std::slice::from_ref(&a),
                    &descriptor.to_string(),
                )
                .unwrap();
                value.coefficient_domain = descriptor.clone();
                let p: Polynomial = value.into();
                let q = p.integrate_checked(0, true).unwrap();
                assert_eq!(
                    q.diff(0).field.as_ref().unwrap().rep,
                    p.field.as_ref().unwrap().rep
                );
                assert_eq!(
                    q.field.as_ref().unwrap().coefficient_domain,
                    descriptor.field()
                );
                if ring && base == "ZZ" {
                    assert!(matches!(
                        p.integrate_checked(0, false),
                        Err(PolyError::ExactQuotientFailed(_))
                    ));
                } else {
                    assert_eq!(
                        p.integrate_checked(0, false)
                            .unwrap()
                            .field
                            .as_ref()
                            .unwrap()
                            .coefficient_domain,
                        descriptor
                    );
                }
            }
        }
        let field = NumberField::from_generator(&parse_str("sqrt(2)").unwrap()).unwrap();
        let p: Polynomial = AlgebraicPoly::from_expr(
            &parse_str("sqrt(2)*x").unwrap(),
            std::slice::from_ref(&x),
            field,
        )
        .unwrap()
        .into();
        assert_eq!(
            p.integrate_checked(0, false)
                .unwrap()
                .diff(0)
                .algebraic
                .unwrap()
                .rep,
            p.algebraic.unwrap().rep
        );
        let p: Polynomial = Poly::from_expr(&parse_str("x").unwrap(), std::slice::from_ref(&x))
            .unwrap()
            .set_modulus(5)
            .into();
        assert_eq!(
            p.integrate_checked(0, true).unwrap().diff(0).base.rep,
            p.base.rep
        );
        let p: Polynomial = Poly::from_expr(&parse_str("x^4").unwrap(), std::slice::from_ref(&x))
            .unwrap()
            .set_modulus(5)
            .into();
        assert!(matches!(
            p.integrate_checked(0, true),
            Err(PolyError::NotInvertible(_))
        ));
    }
}
