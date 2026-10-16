//! Typed polynomial adapters for elimination by subresultant PRS.
use super::algebraic::AlgebraicScope;
use super::convert::{PResult, PolyError};
use super::f5b::{Fp, PrimeScope};
use super::fraction::Polynomial;
use super::int::{Coeff, Rat};
use super::mpoly::MPoly;
use super::pseudo::PseudoCoeff;
use super::ratfunc::Scope;

fn calculate<C: PseudoCoeff>(
    f: &MPoly<C>,
    g: &MPoly<C>,
    discriminant: bool,
    wrap: impl Fn(MPoly<C>) -> PResult<Polynomial>,
) -> PResult<(Polynomial, Vec<Polynomial>)> {
    if discriminant {
        return Ok((wrap(super::prs::discriminant(f, 0)?)?, Vec::new()));
    }
    let result = super::prs::subresultants(f, g, 0)?;
    Ok((
        wrap(result.resultant)?,
        result
            .sequence
            .into_iter()
            .map(wrap)
            .collect::<PResult<_>>()?,
    ))
}

fn remove<C: Coeff>(p: &MPoly<C>) -> PResult<MPoly<C>> {
    if p.nvars < 2 || p.degree(0) > 0 {
        return Err(PolyError::Failed(
            "elimination result still contains its first generator".into(),
        ));
    }
    Ok(MPoly::from_terms(
        p.nvars - 1,
        p.iter()
            .map(|(e, c)| (e[1..].to_vec(), c.clone()))
            .collect(),
    ))
}

impl Polynomial {
    fn subresultant_data(&self, other: &Self, discriminant: bool) -> PResult<(Self, Vec<Self>)> {
        if let Some((a, b)) = self.checked_algebraic_pair(other)? {
            let _scope = AlgebraicScope::new(a.number_field.clone());
            return calculate(&a.rep, &b.rep, discriminant, |r| Ok(a.with_rep(r).into()));
        }
        if let Some((a, b)) = self.checked_fraction_pair(other)? {
            let _scope = Scope::new(a.parameters.len());
            return calculate(&a.rep, &b.rep, discriminant, |r| {
                let result = a.with_rep(r);
                if result.display.ends_with(']') && !result.belongs_to_ring() {
                    return Err(PolyError::Failed(
                        "subresultant normalization escaped its coefficient ring".into(),
                    ));
                }
                Ok(result.into())
            });
        }
        if let Some(prime) = self.modulus.or(other.modulus) {
            if !self.base.domain_gens().is_empty() || !other.base.domain_gens().is_empty() {
                return Err(PolyError::Unsupported(
                    "finite-field subresultants require scalar coefficients".into(),
                ));
            }
            if self.modulus.is_some() && other.modulus.is_some() && self.modulus != other.modulus {
                return Err(PolyError::Unsupported(
                    "subresultants over incompatible finite fields".into(),
                ));
            }
            let (a, b) = self.base.unify(&other.base);
            let _scope = PrimeScope::new(prime)?;
            let convert = |p: &MPoly<Rat>| -> PResult<MPoly<Fp>> {
                Ok(MPoly::from_terms(
                    p.nvars,
                    p.iter()
                        .map(|(e, c)| Ok((e.to_vec(), Fp::from_rat(c, prime)?)))
                        .collect::<PResult<_>>()?,
                ))
            };
            return calculate(&convert(&a.rep)?, &convert(&b.rep)?, discriminant, |r| {
                Ok(a.with_rep(r.map_coeffs(Coeff::to_rat))
                    .set_modulus(prime)
                    .into())
            });
        }
        let (a, b) = self.base.unify(&other.base);
        calculate(&a.rep, &b.rep, discriminant, |r| {
            if a.domain == super::polyobj::Domain::ZZ && !r.is_integral() {
                return Err(PolyError::Failed(
                    "subresultant normalization escaped ZZ".into(),
                ));
            }
            Ok(a.with_rep(r).into())
        })
    }

    /// Subresultant PRS, with every original main generator retained.
    pub fn subresultants(&self, other: &Self) -> PResult<Vec<Self>> {
        Ok(self.subresultant_data(other, false)?.1)
    }
    /// Resultant and PRS, still represented with the eliminated variable at
    /// exponent zero. Public wrappers remove that variable only from the result.
    pub fn resultant_prs(&self, other: &Self) -> PResult<(Self, Vec<Self>)> {
        let (result, mut sequence) = self.subresultant_data(other, false)?;
        // A resultant involving zero uses the empty PRS convention, while
        // subresultants(0,g) independently returns [g].
        if sequence.len() < 2 {
            sequence.clear();
        }
        Ok((result, sequence))
    }
    pub fn discriminant_poly(&self) -> PResult<Self> {
        Ok(self.subresultant_data(self, true)?.0)
    }
    /// Drop precisely the eliminated main generator and preserve typed
    /// coefficient descriptors, including unused parameters and aliases.
    pub fn eliminate_first(&self) -> PResult<Self> {
        if self.true_gens().len() < 2 {
            return Err(PolyError::Failed(
                "an expression represents a univariate elimination result".into(),
            ));
        }
        if let Some(a) = &self.algebraic {
            let mut result = a.with_rep(remove(&a.rep)?);
            result.gens.remove(0);
            return Ok(result.into());
        }
        if let Some(a) = &self.field {
            let mut result = a.with_rep(remove(&a.rep)?);
            result.gens.remove(0);
            return Ok(result.into());
        }
        let mut result = self.base.with_rep(remove(&self.base.rep)?);
        result.gens.remove(0);
        result.ngens -= 1;
        Ok(result.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        algebraic::{AlgebraicPoly, NumberField},
        fraction::FractionPoly,
        polyobj::Poly,
    };
    use crate::{parser::parse_str, Expr};
    #[test]
    fn adapters_preserve_domains_and_drop_only_the_first_generator() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let z = Expr::symbol("z");
        let a = Expr::symbol("a");
        let f = parse_str("x^2+y*x+z").unwrap();
        let g = parse_str("2*x+y").unwrap();
        let gens = vec![x, y.clone(), z.clone()];
        let mut cases: Vec<Polynomial> = vec![
            Poly::from_expr(&f, &gens).unwrap().into(),
            Poly::from_expr(&f, &gens).unwrap().set_modulus(5).into(),
        ];
        let field = NumberField::from_generator(&parse_str("2^(1/3)+I").unwrap()).unwrap();
        cases.push(AlgebraicPoly::from_expr(&f, &gens, field).unwrap().into());
        for domain in ["ZZ[a]", "QQ(a)"] {
            cases.push(
                FractionPoly::from_expr(&f, &gens, std::slice::from_ref(&a), domain)
                    .unwrap()
                    .into(),
            );
        }
        for p in cases {
            let q: Polynomial = if let Some(field) = &p.field {
                field.convert(&g).unwrap().into()
            } else if let Some(field) = &p.algebraic {
                field.convert(&g).unwrap().into()
            } else {
                p.base
                    .with_rep(Poly::from_expr(&g, &gens).unwrap().rep)
                    .into()
            };
            let (result, sequence) = p.resultant_prs(&q).unwrap();
            assert_eq!(sequence.len(), 3);
            for stage in sequence {
                assert_eq!(stage.true_gens(), gens);
                assert_eq!(stage.domain_str(), p.domain_str());
            }
            for r in [result, p.discriminant_poly().unwrap()] {
                let r = r.eliminate_first().unwrap();
                assert_eq!(r.true_gens(), [y.clone(), z.clone()]);
                assert_eq!(r.domain_str(), p.domain_str());
                if let Some(field) = &r.field {
                    assert_eq!(field.parameters, [a.clone()]);
                    assert!(field.rep.coeffs().iter().all(|c| c.n.nvars == 1));
                }
                if let Some(field) = &r.algebraic {
                    assert_eq!(
                        field.number_field,
                        p.algebraic.as_ref().unwrap().number_field
                    );
                }
            }
        }
    }
}
