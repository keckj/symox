//! Immutable, typed exact coefficient domains and closed domain elements.
use super::algebraic::{AlgebraicScope, NumberField, AN};
use super::convert::{PResult, PolyError};
use super::f5b::Fp;
use super::fraction::{coefficient, coefficient_expr, coefficient_in_ring};
use super::int::{Coeff, Field, Rat};
use super::ratfunc::{Scope, RF};
use crate::Expr;
use std::sync::Arc;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ExactDomain {
    ZZ,
    QQ,
    RR,
    GF(u64),
    Algebraic(Arc<NumberField>),
    Functions {
        integral: bool,
        ring: bool,
        parameters: Vec<Expr>,
        order: String,
    },
}
impl std::fmt::Display for ExactDomain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZZ => write!(f, "ZZ"),
            Self::QQ => write!(f, "QQ"),
            Self::RR => write!(f, "RR"),
            Self::GF(p) => write!(f, "GF({p})"),
            Self::Algebraic(k) => write!(f, "{}", k.display()),
            Self::Functions {
                integral,
                ring,
                parameters,
                ..
            } => write!(
                f,
                "{}{}{}{}",
                if *integral { "ZZ" } else { "QQ" },
                if *ring { "[" } else { "(" },
                parameters
                    .iter()
                    .map(crate::printer::to_str)
                    .collect::<Vec<_>>()
                    .join(","),
                if *ring { "]" } else { ")" }
            ),
        }
    }
}
fn failure(s: impl Into<String>) -> PolyError {
    PolyError::Failed(s.into())
}
impl ExactDomain {
    pub fn functions(
        integral: bool,
        ring: bool,
        parameters: Vec<Expr>,
        order: &str,
    ) -> PResult<Self> {
        if parameters.is_empty() {
            return Err(failure(
                "a polynomial ring or fraction field needs generators",
            ));
        }
        if super::groebner::MonomialOrder::parse(order).is_none() {
            return Err(failure("order must be lex, grlex, or grevlex"));
        }
        for (i, p) in parameters.iter().enumerate() {
            if !p.is_symbol() || parameters[..i].contains(p) {
                return Err(failure("domain generators must be distinct symbols"));
            }
        }
        Ok(Self::Functions {
            integral,
            ring,
            parameters,
            order: order.into(),
        })
    }
    pub fn finite(prime: u64) -> PResult<Self> {
        Fp::from_rat(&Rat::zero(), prime)?;
        Ok(Self::GF(prime))
    }
    pub fn is_field(&self) -> bool {
        !matches!(self, Self::ZZ | Self::Functions { ring: true, .. })
    }
    pub fn characteristic(&self) -> u64 {
        if let Self::GF(p) = self {
            *p
        } else {
            0
        }
    }
    pub fn field(&self) -> Self {
        match self {
            Self::ZZ => Self::QQ,
            Self::Functions {
                integral,
                parameters,
                order,
                ..
            } => Self::Functions {
                integral: *integral,
                ring: false,
                parameters: parameters.clone(),
                order: order.clone(),
            },
            _ => self.clone(),
        }
    }
    pub fn ring(&self) -> PResult<Self> {
        match self {
            Self::ZZ | Self::QQ => Ok(Self::ZZ),
            Self::Functions {
                integral,
                parameters,
                order,
                ..
            } => Ok(Self::Functions {
                integral: *integral,
                ring: true,
                parameters: parameters.clone(),
                order: order.clone(),
            }),
            _ => Err(PolyError::Unsupported(
                "this field has no implemented associated ring".into(),
            )),
        }
    }
    pub fn convert(&self, e: &Expr) -> PResult<DomainElement> {
        let value=match self {
            Self::RR=>return Err(PolyError::Unsupported("real coefficient elements require an inexact domain backend that is not implemented".into())),
            Self::ZZ|Self::QQ|Self::GF(_)=>{
                let r=e.as_number().and_then(Rat::from_number).ok_or_else(||failure("expected an exact rational number"))?;
                if matches!(self,Self::ZZ|Self::GF(_))&&!r.is_integer(){return Err(failure("expected an integer"));}
                if let Self::GF(p)=self{ElementValue::Finite(Fp::from_rat(&r,*p)?)}else{ElementValue::Rational(r)}
            },
            Self::Algebraic(field)=>{let _scope=AlgebraicScope::new(field.clone());ElementValue::Algebraic(AN::from_expr(e,field.clone())?)},
            Self::Functions{integral,ring,parameters,..}=>{
                let _scope=Scope::new(parameters.len());let value=coefficient(e,parameters)?;
                if *ring&&!coefficient_in_ring(&value,*integral){return Err(failure("coefficient does not belong to this polynomial ring"));}
                ElementValue::Function(value)
            }
        };
        Ok(DomainElement {
            domain: self.clone(),
            value,
        })
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ElementValue {
    Rational(Rat),
    Finite(Fp),
    Algebraic(AN),
    Function(RF),
}
#[derive(Clone, Debug)]
pub struct DomainElement {
    pub domain: ExactDomain,
    pub value: ElementValue,
}
impl PartialEq for DomainElement {
    fn eq(&self, other: &Self) -> bool {
        match (&self.value, &other.value) {
            (ElementValue::Algebraic(a), ElementValue::Algebraic(b)) => {
                a.field.same_embedding(&b.field) && a.coefficients == b.coefficients
            }
            _ => self.domain == other.domain && self.value == other.value,
        }
    }
}
impl Eq for DomainElement {}
impl std::hash::Hash for DomainElement {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        if let ElementValue::Algebraic(a) = &self.value {
            3u8.hash(state);
            a.field.generator.hash(state);
            a.field.minimal.hash(state);
            a.field.radicand.hash(state);
            a.field.model.hash(state);
            a.coefficients.hash(state);
        } else {
            self.domain.hash(state);
            self.value.hash(state);
        }
    }
}
impl DomainElement {
    pub fn compatible_embedding(&self, other: &Self) -> bool {
        self.domain == other.domain
            || matches!((&self.value,&other.value),(ElementValue::Algebraic(a),ElementValue::Algebraic(b)) if a.field.same_embedding(&b.field))
    }
    pub fn as_expr(&self) -> Expr {
        match &self.value {
            ElementValue::Rational(r) => Expr::num(r.to_number()),
            ElementValue::Finite(f) => Expr::num(f.to_rat().to_number()),
            ElementValue::Algebraic(a) => a.as_expr(),
            ElementValue::Function(r) => {
                let ExactDomain::Functions { parameters, .. } = &self.domain else {
                    unreachable!()
                };
                coefficient_expr(r, parameters)
            }
        }
    }
    pub fn is_zero(&self) -> bool {
        match &self.value {
            ElementValue::Rational(r) => r.is_zero(),
            ElementValue::Finite(r) => r.is_zero(),
            ElementValue::Algebraic(r) => r.is_zero(),
            ElementValue::Function(r) => r.is_zero(),
        }
    }
    pub fn binary(&self, other: &Self, op: char) -> PResult<Self> {
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if op == '/' && other.is_zero() {
            return Err(PolyError::Unsupported("division by zero".into()));
        }
        macro_rules! calc {
            ($a:expr,$b:expr) => {
                match op {
                    '+' => $a.add($b),
                    '-' => $a.sub($b),
                    '*' => $a.mul($b),
                    '/' => $a.div($b),
                    _ => return Err(failure("unknown coefficient operation")),
                }
            };
        }
        let value = match (&self.value, &other.value) {
            (ElementValue::Rational(a), ElementValue::Rational(b)) => {
                ElementValue::Rational(calc!(a, b))
            }
            (ElementValue::Finite(a), ElementValue::Finite(b)) => ElementValue::Finite(calc!(a, b)),
            (ElementValue::Algebraic(a), ElementValue::Algebraic(b)) => {
                let _scope = AlgebraicScope::new(a.field.clone());
                ElementValue::Algebraic(calc!(a, b))
            }
            (ElementValue::Function(a), ElementValue::Function(b)) => {
                let ExactDomain::Functions { parameters, .. } = &self.domain else {
                    unreachable!()
                };
                let _scope = Scope::new(parameters.len());
                ElementValue::Function(calc!(a, b))
            }
            _ => return Err(failure("incompatible coefficient representations")),
        };
        let result = Self {
            domain: self.domain.clone(),
            value,
        };
        match (&self.domain, &result.value) {
            (ExactDomain::ZZ, ElementValue::Rational(r)) if !r.is_integer() => {
                return Err(failure("non-exact division in ZZ"))
            }
            (
                ExactDomain::Functions {
                    integral,
                    ring: true,
                    ..
                },
                ElementValue::Function(r),
            ) if !coefficient_in_ring(r, *integral) => {
                return Err(failure("non-exact division in polynomial ring"))
            }
            _ => {}
        }
        Ok(result)
    }
    pub fn pow(&self, exponent: i64) -> PResult<Self> {
        if exponent < 0 && self.is_zero() {
            return Err(PolyError::Unsupported("division by zero".into()));
        }
        let mut result = self.domain.convert(&Expr::one())?;
        let mut base = self.clone();
        let mut n = exponent.unsigned_abs();
        while n > 0 {
            if n & 1 != 0 {
                result = result.binary(&base, '*')?;
            }
            n >>= 1;
            if n > 0 {
                base = base.binary(&base, '*')?;
            }
        }
        if exponent < 0 {
            self.domain.convert(&Expr::one())?.binary(&result, '/')
        } else {
            Ok(result)
        }
    }
}

impl DomainElement {
    pub fn gcd(&self, other: &Self) -> PResult<Self> {
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        let value = match (&self.value, &other.value) {
            (ElementValue::Rational(a), ElementValue::Rational(b)) => {
                ElementValue::Rational(Rat::new(a.n.gcd(&b.n), a.d.lcm(&b.d)))
            }
            (ElementValue::Function(a), ElementValue::Function(b)) => {
                let ExactDomain::Functions {
                    integral,
                    parameters,
                    ..
                } = &self.domain
                else {
                    unreachable!()
                };
                let _scope = Scope::new(parameters.len());
                if matches!(
                    self.domain,
                    ExactDomain::Functions {
                        integral: false,
                        ring: true,
                        ..
                    }
                ) {
                    let as_poly = |r: &RF| -> PResult<super::MPoly<Rat>> {
                        if !r.d.is_constant() {
                            return Err(failure("coefficient is not in a polynomial ring"));
                        }
                        Ok(r.n.to_rat().scale(&Rat::new(super::Int::one(), r.d.lc())))
                    };
                    let a = as_poly(a)?;
                    let b = as_poly(b)?;
                    let h = if a.is_zero() {
                        if b.lc().is_negative() {
                            b.neg()
                        } else {
                            b
                        }
                    } else if b.is_zero() {
                        if a.lc().is_negative() {
                            a.neg()
                        } else {
                            a
                        }
                    } else if a.len() == 1 || b.len() == 1 {
                        let mut powers = vec![u32::MAX; a.nvars];
                        let mut coefficient = Rat::zero();
                        for (exponents, c) in a.iter().chain(b.iter()) {
                            for (p, e) in powers.iter_mut().zip(exponents) {
                                *p = (*p).min(*e);
                            }
                            coefficient =
                                Rat::new(coefficient.n.gcd(&c.n), coefficient.d.lcm(&c.d));
                        }
                        super::MPoly::monomial(&powers, coefficient)
                    } else {
                        super::gcd::gcd_mq(&a, &b).0
                    };
                    return Ok(Self {
                        domain: self.domain.clone(),
                        value: ElementValue::Function(RF::from_rat_poly(&h)),
                    });
                }
                let (n, _, _) = super::gcd::gcd_mz(&a.n, &b.n);
                let (d, _, _) = super::gcd::gcd_mz(&a.d, &b.d);
                let mut denominator =
                    a.d.div_exact(&d)
                        .ok_or_else(|| failure("denominator gcd certificate failed"))?
                        .mul(&b.d);
                let mut numerator = n;
                if !*integral && !numerator.is_zero() {
                    let denominator_lc = denominator.lc();
                    let numerator_lc =
                        if !a.n.is_zero() && !b.n.is_zero() && a.n.len() > 1 && b.n.len() > 1 {
                            numerator.lc()
                        } else {
                            super::Int::one()
                        };
                    numerator = numerator.scale(&denominator_lc);
                    denominator = denominator.scale(&numerator_lc);
                }
                ElementValue::Function(RF::new(numerator, denominator))
            }
            _ => {
                return self.domain.convert(&if self.is_zero() && other.is_zero() {
                    Expr::zero()
                } else {
                    Expr::one()
                })
            }
        };
        Ok(Self {
            domain: self.domain.clone(),
            value,
        })
    }
    pub fn lcm(&self, other: &Self) -> PResult<Self> {
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if self.is_zero() || other.is_zero() {
            return self.domain.convert(&Expr::zero());
        }
        if let (ElementValue::Rational(a), ElementValue::Rational(b)) = (&self.value, &other.value)
        {
            return Ok(Self {
                domain: self.domain.clone(),
                value: ElementValue::Rational(Rat::new(a.n.lcm(&b.n), a.d.gcd(&b.d))),
            });
        }
        if let (ElementValue::Function(a), ElementValue::Function(b)) = (&self.value, &other.value)
        {
            let ExactDomain::Functions {
                integral,
                parameters,
                ..
            } = &self.domain
            else {
                unreachable!()
            };
            let _scope = Scope::new(parameters.len());
            let (gn, _, _) = super::gcd::gcd_mz(&a.n, &b.n);
            let (gd, _, _) = super::gcd::gcd_mz(&a.d, &b.d);
            let mut n =
                a.n.div_exact(&gn)
                    .ok_or_else(|| failure("numerator lcm certificate failed"))?
                    .mul(&b.n);
            let mut d = gd;
            if !*integral {
                let nlc = n.lc();
                let dlc = if a.d.len() > 1 && b.d.len() > 1 {
                    d.lc()
                } else {
                    super::Int::one()
                };
                n = n.scale(&dlc);
                d = d.scale(&nlc);
            }
            return Ok(Self {
                domain: self.domain.clone(),
                value: ElementValue::Function(RF::new(n, d)),
            });
        }
        if matches!(self.domain, ExactDomain::GF(_)) {
            return self.domain.convert(&Expr::one());
        }
        self.binary(other, '*')?.binary(&self.gcd(other)?, '/')
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn exact_closed_fields_interleaved_scopes() {
        let a = Expr::symbol("a");
        let b = Expr::symbol("b");
        let fields = vec![
            ExactDomain::QQ,
            ExactDomain::finite(5).unwrap(),
            ExactDomain::finite(7).unwrap(),
            ExactDomain::functions(false, false, vec![a.clone(), b], "grevlex").unwrap(),
            ExactDomain::Algebraic(
                NumberField::from_generator(&parse_str("sqrt(2)").unwrap()).unwrap(),
            ),
            ExactDomain::Algebraic(
                NumberField::from_generator(&parse_str("2^(1/3)+sqrt(3)").unwrap()).unwrap(),
            ),
        ];
        let values: Vec<_> = fields
            .iter()
            .map(|k| k.convert(&Expr::int(3)).unwrap())
            .collect();
        for value in &values {
            let inverse = value.pow(-1).unwrap();
            assert_eq!(
                value.binary(&inverse, '*').unwrap(),
                value.domain.convert(&Expr::one()).unwrap()
            );
        }
        for i in 0..values.len() {
            for j in 0..values.len() {
                if i != j {
                    assert!(values[i].binary(&values[j], '+').is_err());
                }
            }
        }
        let ring = ExactDomain::functions(false, true, vec![a.clone()], "lex").unwrap();
        assert!(ring.convert(&a.powi(-1)).is_err());
        assert!(ring.convert(&a).unwrap().pow(-1).is_err());
        assert!(ExactDomain::ZZ
            .convert(&Expr::int(2))
            .unwrap()
            .pow(-1)
            .is_err());
    }
    #[test]
    fn aliases_are_descriptor_metadata_not_element_embeddings() {
        let a = NumberField::from_generator(&parse_str("sqrt(2)").unwrap()).unwrap();
        let mut b = (*a).clone();
        b.alias = Some(Expr::symbol("beta"));
        let k = ExactDomain::Algebraic(a);
        let l = ExactDomain::Algebraic(Arc::new(b));
        assert_ne!(k, l);
        let x = k.convert(&parse_str("sqrt(2)").unwrap()).unwrap();
        let y = l.convert(&parse_str("sqrt(2)").unwrap()).unwrap();
        assert_eq!(x, y);
        assert_eq!(x.binary(&y, '/').unwrap(), k.convert(&Expr::one()).unwrap());
        let mut values = std::collections::HashSet::new();
        values.insert(x);
        values.insert(y);
        assert_eq!(values.len(), 1);
    }
    #[test]
    fn exact_domain_gcd_and_fraction_certificates() {
        let k = ExactDomain::QQ;
        let a = k.convert(&parse_str("2/3").unwrap()).unwrap();
        let b = k.convert(&parse_str("4/5").unwrap()).unwrap();
        assert_eq!(a.gcd(&b).unwrap().as_expr(), parse_str("2/15").unwrap());
        assert_eq!(a.lcm(&b).unwrap().as_expr(), Expr::int(4));
        let x = Expr::symbol("x");
        let k = ExactDomain::functions(false, false, vec![x], "lex").unwrap();
        let a = k.convert(&parse_str("2*x/(x+1)").unwrap()).unwrap();
        let b = k.convert(&parse_str("3*x/(x+2)").unwrap()).unwrap();
        let g = a.gcd(&b).unwrap();
        assert_eq!(
            g,
            k.convert(&parse_str("x/((x+1)*(x+2))").unwrap()).unwrap()
        );
    }
}
