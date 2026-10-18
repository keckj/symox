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
    Expression { raw: bool },
    Inexact { precision: u32, complex: bool },
    Gaussian { integral: bool },
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
            Self::Expression { raw } => write!(f,"{}",if *raw {"EXRAW"} else {"EX"}),
            Self::Inexact {complex,..} => write!(f,"{}",if *complex {"CC"} else {"RR"}),
            Self::Gaussian {integral} => write!(f,"{}",if *integral {"ZZ_I"} else {"QQ_I"}),
        }
    }
}
fn failure(s: impl Into<String>) -> PolyError {
    PolyError::Failed(s.into())
}
impl ExactDomain {
    pub fn inexact_precision(&self)->Option<u32> {
        match self {Self::RR=>Some(53),Self::Inexact{precision,..}=>Some(*precision),_=>None}
    }
    pub fn is_complex(&self)->bool {matches!(self,Self::Inexact{complex:true,..})}
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
        !matches!(self, Self::ZZ | Self::Functions { ring: true, .. } | Self::Gaussian {integral:true})
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
            Self::Gaussian {..} => Self::Gaussian {integral:false},
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
            Self::Gaussian {..} => Ok(Self::Gaussian {integral:true}),
            Self::Expression { raw:false } => Ok(self.clone()),
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
            Self::RR | Self::Inexact{..} => {
                let precision=self.inexact_precision().unwrap();
                let value=crate::evalf_mp::eval_mp(e,precision).ok_or_else(||failure("expected a finite supported numerical value"))?;
                if !self.is_complex()&&!value.is_real(){return Err(failure("expected a real numerical value"));}
                ElementValue::Inexact(crate::evalf_mp::to_expr(value,precision))
            },
            Self::ZZ|Self::QQ|Self::GF(_)=>{
                let r=e.as_number().and_then(Rat::from_number).ok_or_else(||failure("expected an exact rational number"))?;
                if matches!(self,Self::ZZ|Self::GF(_))&&!r.is_integer(){return Err(failure("expected an integer"));}
                if let Self::GF(p)=self{ElementValue::Finite(Fp::from_rat(&r,*p)?)}else{ElementValue::Rational(r)}
            },
            Self::Algebraic(field)=>{let _scope=AlgebraicScope::new(field.clone());ElementValue::Algebraic(AN::from_expr(e,field.clone())?)},
            Self::Gaussian {integral} => {
                let field=NumberField::from_generator(&Expr::i())?;
                let _scope=AlgebraicScope::new(field.clone());
                let value=AN::from_expr(e,field)?;
                if *integral&&value.coefficients.iter().any(|c|!c.is_integer()) {return Err(failure("expected Gaussian integer coordinates"));}
                ElementValue::Algebraic(value)
            },
            Self::Functions{integral,ring,parameters,..}=>{
                let _scope=Scope::new(parameters.len());let value=coefficient(e,parameters)?;
                if *ring&&!coefficient_in_ring(&value,*integral){return Err(failure("coefficient does not belong to this polynomial ring"));}
                ElementValue::Function(value)
            },
            Self::Expression { raw } => {
                if crate::decompose::commutativity(e)==Some(false) {
                    return Err(failure("expression coefficients must be commutative"));
                }
                ElementValue::Expression(if *raw {e.clone()} else {crate::simplify::simplify(e)})
            },
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
    Expression(Expr),
    Inexact(Expr),
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
        if self.domain!=other.domain&&(matches!(self.domain,ExactDomain::Gaussian{..})||matches!(other.domain,ExactDomain::Gaussian{..})) {return false;}
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
            ElementValue::Expression(e) | ElementValue::Inexact(e) => e.clone(),
        }
    }
    pub fn is_zero(&self) -> bool {
        match &self.value {
            ElementValue::Rational(r) => r.is_zero(),
            ElementValue::Finite(r) => r.is_zero(),
            ElementValue::Algebraic(r) => r.is_zero(),
            ElementValue::Function(r) => r.is_zero(),
            ElementValue::Expression(e) | ElementValue::Inexact(e) => e.is_zero(),
        }
    }
    /// Exact multivariate polynomial view of a ring element. Rational function
    /// denominators are constant in QQ[parameters], and one in ZZ[parameters].
    pub fn ring_polynomial(&self) -> PResult<super::MPoly<Rat>> {
        let (
            ExactDomain::Functions {
                integral,
                ring: true,
                parameters,
                ..
            },
            ElementValue::Function(r),
        ) = (&self.domain, &self.value)
        else {
            return Err(failure("a polynomial-ring element is required"));
        };
        if r.n.nvars != parameters.len()
            || r.d.nvars != parameters.len()
            || r.d.is_zero()
            || !coefficient_in_ring(r, *integral)
        {
            return Err(failure(
                "invalid polynomial-ring coefficient representation",
            ));
        }
        Ok(r.n.to_rat().scale(&Rat::new(super::Int::one(), r.d.lc())))
    }
    pub fn div_rem(&self, other: &Self) -> PResult<(Self, Self)> {
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if other.is_zero() {
            return Err(PolyError::Unsupported("division by zero".into()));
        }
        if matches!(self.domain,ExactDomain::Gaussian{integral:true}) {
            let (ElementValue::Algebraic(a),ElementValue::Algebraic(b))=(&self.value,&other.value) else {unreachable!()};
            let _scope=AlgebraicScope::new(a.field.clone());
            let quotient=a.div(b);
            let rounded=quotient.coefficients.iter().map(|q| {
                let shifted=q.add(&Rat::new(super::Int::one(),super::Int::from_i64(2)));
                Rat::from(shifted.n.div_mod_floor(&shifted.d).0)
            }).collect();
            let q=Self{domain:self.domain.clone(),value:ElementValue::Algebraic(AN::new(rounded,a.field.clone()))};
            let r=self.binary(&q.binary(other,'*')?,'-')?;
            return Ok((q,r));
        }
        if let ExactDomain::Functions {
            integral,
            ring: true,
            parameters,
            order,
        } = &self.domain
        {
            let a = self.ring_polynomial()?;
            let b = other.ring_polynomial()?;
            let order = super::groebner::MonomialOrder::parse(order)
                .ok_or_else(|| failure("invalid monomial order"))?;
            let (mut quotients, remainder) = super::groebner::reduce_with_options(
                &a,
                &[b],
                order,
                super::groebner::CoefficientField { modulus: None },
                *integral,
            );
            let _scope = Scope::new(parameters.len());
            let wrap = |p: &super::MPoly<Rat>| Self {
                domain: self.domain.clone(),
                value: ElementValue::Function(RF::from_rat_poly(p)),
            };
            return Ok((wrap(&quotients.remove(0)), wrap(&remainder)));
        }
        if let (ExactDomain::ZZ, ElementValue::Rational(a), ElementValue::Rational(b)) =
            (&self.domain, &self.value, &other.value)
        {
            let (q, r) = a.n.div_mod_floor(&b.n);
            let wrap = |v: super::Int| Self {
                domain: ExactDomain::ZZ,
                value: ElementValue::Rational(Rat::from(v)),
            };
            return Ok((wrap(q), wrap(r)));
        }
        if self.domain.is_field() {
            return Ok((
                self.binary(other, '/')?,
                self.domain.convert(&Expr::zero())?,
            ));
        }
        Err(PolyError::Unsupported(
            "division is not implemented for this coefficient domain".into(),
        ))
    }
    pub fn exquo(&self, other: &Self) -> PResult<Self> {
        let (q, r) = self.div_rem(other)?;
        if !r.is_zero() {
            return Err(PolyError::ExactQuotientFailed(
                "polynomial or integer division is not exact in the selected ring".into(),
            ));
        }
        Ok(q)
    }
    pub fn gcdex(&self, other: &Self) -> PResult<(Self, Self, Self)> {
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if matches!(self.domain,ExactDomain::Gaussian{integral:true}) {
            let zero=self.domain.convert(&Expr::zero())?;
            let one=self.domain.convert(&Expr::one())?;
            let (mut r0,mut r1)=(self.clone(),other.clone());
            let (mut s0,mut s1,mut t0,mut t1)=(one.clone(),zero.clone(),zero,one);
            while !r1.is_zero() {
                let (q,r)=r0.div_rem(&r1)?;
                let next_s=s0.binary(&q.binary(&s1,'*')?,'-')?;
                let next_t=t0.binary(&q.binary(&t1,'*')?,'-')?;
                (r0,r1)=(r1,r);(s0,s1)=(s1,next_s);(t0,t1)=(t1,next_t);
            }
            let unit=r0.gaussian_canonical_unit()?;
            return Ok((s0.binary(&unit,'*')?,t0.binary(&unit,'*')?,r0.binary(&unit,'*')?));
        }
        if let (ExactDomain::ZZ, ElementValue::Rational(a), ElementValue::Rational(b)) =
            (&self.domain, &self.value, &other.value)
        {
            let (g, s, t) = a.n.gcdex(&b.n);
            let wrap = |v: super::Int| Self {
                domain: ExactDomain::ZZ,
                value: ElementValue::Rational(Rat::from(v)),
            };
            return Ok((wrap(s), wrap(t), wrap(g)));
        }
        if !self.domain.is_field() {
            return Err(PolyError::Unsupported(
                "extended gcd requires a supported coefficient field or ZZ".into(),
            ));
        }
        let zero = self.domain.convert(&Expr::zero())?;
        if self.is_zero() && other.is_zero() {
            return Ok((zero.clone(), self.domain.convert(&Expr::one())?, zero));
        }
        let g = self.gcd(other)?;
        if self.is_zero() {
            Ok((zero, g.binary(other, '/')?, g))
        } else {
            Ok((g.binary(self, '/')?, zero, g))
        }
    }
    pub fn binary(&self, other: &Self, op: char) -> PResult<Self> {
        if op == '/' && matches!(self.domain, ExactDomain::Functions { ring: true, .. } | ExactDomain::Gaussian {integral:true}) {
            return self.exquo(other);
        }
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
            (ElementValue::Expression(a),ElementValue::Expression(b)) |
            (ElementValue::Inexact(a),ElementValue::Inexact(b)) => {
                let expression=match op {
                    '+'=>crate::ops::add2(a,b), '-'=>crate::ops::sub(a,b),
                    '*'=>crate::ops::mul2(a,b), '/'=>crate::ops::div(a,b),
                    _=>return Err(failure("unknown coefficient operation")),
                };
                return self.domain.convert(&expression);
            },
            _ => return Err(failure("incompatible coefficient representations")),
        };
        let result = Self {
            domain: self.domain.clone(),
            value,
        };
        match (&self.domain, &result.value) {
            (ExactDomain::Gaussian{integral:true},ElementValue::Algebraic(a)) if a.coefficients.iter().any(|c|!c.is_integer())=>return Err(failure("non-exact division in Gaussian integer ring")),
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
    pub fn gaussian_parts(&self)->Option<(Rat,Rat)> {
        if !matches!(self.domain,ExactDomain::Gaussian{..}) {return None;}
        let ElementValue::Algebraic(a)=&self.value else {return None;};
        Some((a.coefficients.first().cloned().unwrap_or_else(Rat::zero),a.coefficients.get(1).cloned().unwrap_or_else(Rat::zero)))
    }
    pub fn gaussian_canonical_unit(&self)->PResult<Self> {
        let (x,y)=self.gaussian_parts().ok_or_else(||failure("expected a Gaussian element"))?;
        let quadrant=if !y.is_zero()&&!y.is_negative() {if !x.is_zero()&&!x.is_negative(){0}else{1}}
            else if y.is_negative(){if x.is_negative(){2}else{3}}
            else if x.is_negative(){2}else{0};
        self.domain.convert(&[Expr::one(),crate::ops::neg(&Expr::i()),Expr::int(-1),Expr::i()][quadrant])
    }
    pub fn gaussian_fraction(&self)->PResult<(Self,Self)> {
        let (x,y)=self.gaussian_parts().ok_or_else(||failure("expected a Gaussian element"))?;
        let d=x.d.lcm(&y.d);
        let numerator=crate::ops::mul2(&self.as_expr(),&Expr::num(d.to_number()));
        let ring=ExactDomain::Gaussian{integral:true};
        Ok((ring.convert(&crate::expand::expand(&numerator))?,ring.convert(&Expr::num(d.to_number()))?))
    }
}

/// Recover the stored real value exactly, optionally limiting its denominator
/// by the real-field convention. Continued fractions run on native big ints.
pub fn rationalize_real(e:&Expr, precision:u32, limit:bool)->PResult<Rat> {
    use num_bigint::{BigInt,Sign};
    use num_rational::BigRational;
    let value=crate::evalf_mp::eval_mp(e,precision).ok_or_else(||failure("expected a finite numerical value"))?;
    if !value.is_real(){return Err(failure("a nonreal value has no rational conversion"));}
    let exact=value.re.to_ratio_exact();
    let max_denominator=((BigInt::from(1)<<precision as usize)/BigInt::from(200)).max(BigInt::from(99));
    let result=if !limit || exact.denom()<=&max_denominator {exact} else {
        let negative=exact.numer().sign()==Sign::Minus;
        let mut n=if negative {-exact.numer()} else {exact.numer().clone()};
        let mut d=exact.denom().clone();
        let positive=BigRational::new(n.clone(),d.clone());
        let (mut p0,mut q0,mut p1,mut q1)=(BigInt::from(0),BigInt::from(1),BigInt::from(1),BigInt::from(0));
        loop {
            let a=&n/&d;
            let q2=&q0+&a*&q1;
            if q2>max_denominator {break;}
            (p0,p1)=(p1.clone(),&p0+&a*&p1);
            (q0,q1)=(q1,q2);
            (n,d)=(d.clone(),n-&a*&d);
            if d==BigInt::from(0){break;}
        }
        let k=(&max_denominator-&q0)/&q1;
        let b1=BigRational::new(&p0+&k*&p1,&q0+&k*&q1);
        let b2=BigRational::new(p1,q1);
        // A denominator limit must never erase a small nonzero value.
        let mut best=if b1.numer()==&BigInt::from(0)||b2.numer()==&BigInt::from(0) {positive.clone()}
            else if num_traits::Signed::abs(&(&positive-&b2))<=num_traits::Signed::abs(&(&positive-&b1)) {b2} else {b1};
        if negative {best=-best;}
        best
    };
    Ok(Rat::from_number(&crate::number::Number::from_ratio(result)).unwrap())
}

impl DomainElement {
    pub fn gcd(&self, other: &Self) -> PResult<Self> {
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if matches!(self.domain,ExactDomain::Gaussian{integral:true}) {return self.gcdex(other).map(|(_,_,g)|g);}
        if matches!(self.domain,ExactDomain::Gaussian{integral:false}) {
            let (an,ad)=self.gaussian_fraction()?;let (bn,bd)=other.gaussian_fraction()?;
            let p=an.gcd(&bn)?;let q=ad.lcm(&bd)?;
            return self.domain.convert(&p.as_expr())?.binary(&self.domain.convert(&q.as_expr())?,'/');
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
                        let h = super::gcd::gcd_mq(&a, &b).0;
                        let order = if let ExactDomain::Functions { order, .. } = &self.domain {
                            super::groebner::MonomialOrder::parse(order).unwrap()
                        } else {
                            unreachable!()
                        };
                        h.scale(&h.coeff(order.leading_index(&h)).inv())
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
        if matches!(self.domain,ExactDomain::Gaussian{integral:true}) {
            return self.binary(other,'*')?.binary(&self.gcd(other)?,'/');
        }
        if matches!(self.domain,ExactDomain::Gaussian{integral:false}) {
            let (an,ad)=self.gaussian_fraction()?;let (bn,bd)=other.gaussian_fraction()?;
            let p=an.lcm(&bn)?;let q=ad.gcd(&bd)?;
            return self.domain.convert(&p.as_expr())?.binary(&self.domain.convert(&q.as_expr())?,'/');
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
            if matches!(
                self.domain,
                ExactDomain::Functions {
                    integral: false,
                    ring: true,
                    ..
                }
            ) {
                let a = self.ring_polynomial()?;
                let b = other.ring_polynomial()?;
                let gcd = super::gcd::gcd_mq(&a, &b).0;
                let p = a
                    .div_exact(&gcd)
                    .ok_or_else(|| failure("polynomial lcm certificate failed"))?
                    .mul(&b);
                let order = if let ExactDomain::Functions { order, .. } = &self.domain {
                    super::groebner::MonomialOrder::parse(order).unwrap()
                } else {
                    unreachable!()
                };
                let p = p.scale(&p.coeff(order.leading_index(&p)).inv());
                return Ok(Self {
                    domain: self.domain.clone(),
                    value: ElementValue::Function(RF::from_rat_poly(&p)),
                });
            }
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
    fn gaussian_euclidean_division_and_bezout_certificates() {
        let k=ExactDomain::Gaussian{integral:true};
        let b=k.convert(&parse_str("1+2*I").unwrap()).unwrap();
        for re in -7..=7 {for im in -7..=7 {
            let expression=crate::ops::add2(&Expr::int(re),&crate::ops::mul2(&Expr::i(),&Expr::int(im)));
            let a=k.convert(&expression).unwrap();
            let(q,r)=a.div_rem(&b).unwrap();
            assert_eq!(q.binary(&b,'*').unwrap().binary(&r,'+').unwrap(),a);
            let (rx,ry)=r.gaussian_parts().unwrap();
            assert!(rx.mul(&rx).add(&ry.mul(&ry)).sub(&Rat::from_i64(5)).is_negative());
            let(s,t,g)=a.gcdex(&b).unwrap();
            assert_eq!(s.binary(&a,'*').unwrap().binary(&t.binary(&b,'*').unwrap(),'+').unwrap(),g);
            assert!(a.div_rem(&g).unwrap().1.is_zero());
        }}
        let q=ExactDomain::Gaussian{integral:false};
        let a=q.convert(&parse_str("(1+I)/2").unwrap()).unwrap();
        let b=q.convert(&parse_str("(1-I)/4").unwrap()).unwrap();
        let g=a.gcd(&b).unwrap();
        assert_eq!(g.as_expr(),parse_str("(1+I)/4").unwrap());
        assert_eq!(a.binary(&a.pow(-1).unwrap(),'*').unwrap(),q.convert(&Expr::one()).unwrap());
    }
    #[test]
    fn numerical_fields_preserve_decimal_precision_and_rational_conversion() {
        for precision in [53,100,200] {
            let k=ExactDomain::Inexact{precision,complex:false};
            let literal=crate::parser::parse_str_with_precision("0.1",precision).unwrap();
            let a=k.convert(&literal).unwrap();
            assert_eq!(a.as_expr().float_info().max_precision,Some(precision));
            assert_eq!(rationalize_real(&a.as_expr(),precision,true).unwrap(),Rat::new(1.into(),10.into()));
            let exact=rationalize_real(&a.as_expr(),precision,false).unwrap();
            assert_eq!(k.convert(&Expr::num(exact.to_number())).unwrap(),a);
            let tiny=crate::parser::parse_str_with_precision("1e-100",precision).unwrap();
            assert_eq!(rationalize_real(&tiny,precision,true).unwrap(),rationalize_real(&tiny,precision,false).unwrap());
            let complex=ExactDomain::Inexact{precision,complex:true};
            let value=complex.convert(&crate::parser::parse_str_with_precision("0.1+0.2*I",precision).unwrap()).unwrap();
            let numeric=crate::evalf_mp::eval_mp(&value.as_expr(),precision).unwrap();
            assert_eq!(numeric.re,crate::bigfloat::BigFloat::from_decimal("0.1",precision).unwrap());
            assert_eq!(numeric.im,crate::bigfloat::BigFloat::from_decimal("0.2",precision).unwrap());
            assert!(k.convert(&Expr::i()).is_err());
        }
    }
    #[test]
    fn expression_domains_are_closed_and_keep_raw_presentation() {
        for raw in [false,true] {
            let k=ExactDomain::Expression{raw};
            let a=k.convert(&parse_str("sin(x)+y").unwrap()).unwrap();
            assert_eq!(a.binary(&a,'-').unwrap(),k.convert(&Expr::zero()).unwrap());
            assert_eq!(a.binary(&a,'/').unwrap(),k.convert(&Expr::one()).unwrap());
            let b=k.convert(&parse_str("x+1").unwrap()).unwrap();
            assert_eq!(a.binary(&b,'*').unwrap().binary(&b,'/').unwrap(),a);
        }
    }
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
    #[test]
    fn polynomial_ring_division_certificates_all_orders() {
        let params=vec![Expr::symbol("a"),Expr::symbol("b")];
        for order in ["lex","grlex","grevlex"] {for integral in [false,true] {
            let k=ExactDomain::functions(integral,true,params.clone(),order).unwrap();
            for (a,b) in [("a^2+2*a","2*a"),("a*b+1","2*a+2*b"),("(a+2*b^2)*(a+b)+1","a+2*b^2"),("0","a+1")] {
                let a=k.convert(&parse_str(a).unwrap()).unwrap();let b=k.convert(&parse_str(b).unwrap()).unwrap();
                let(q,r)=a.div_rem(&b).unwrap();
                assert_eq!(q.binary(&b,'*').unwrap().binary(&r,'+').unwrap(),a,"{order} integral={integral}");
                assert_eq!(q.domain,k);assert_eq!(r.domain,k);
                let(q2,r2)=r.div_rem(&b).unwrap();assert!(q2.is_zero());assert_eq!(r2,r);
                if r.is_zero(){assert_eq!(a.exquo(&b).unwrap(),q);}else{assert!(matches!(a.exquo(&b),Err(PolyError::ExactQuotientFailed(_))));}
            }
            let a=k.convert(&parse_str("(a+b)*(a-b)").unwrap()).unwrap();let b=k.convert(&parse_str("a-b").unwrap()).unwrap();
            assert_eq!(a.binary(&b,'/').unwrap(),k.convert(&parse_str("a+b").unwrap()).unwrap());
            assert!(a.div_rem(&k.convert(&Expr::zero()).unwrap()).is_err());
        }}
    }
    #[test]
    fn field_gcdex_normalization_matches_its_gcd() {
        for k in [ExactDomain::QQ,ExactDomain::functions(false,false,vec![Expr::symbol("a")],"grevlex").unwrap()] {
            let pairs=if k==ExactDomain::QQ {vec![("2/3","4/5"),("0","4/5"),("2/3","0"),("0","0")]}else{vec![("2*a/(a+1)","3*a/(a+2)"),("0","3*a/(a+2)"),("a","0"),("0","0")]};
            for(a,b)in pairs{let a=k.convert(&parse_str(a).unwrap()).unwrap();let b=k.convert(&parse_str(b).unwrap()).unwrap();let(s,t,h)=a.gcdex(&b).unwrap();assert_eq!(s.binary(&a,'*').unwrap().binary(&t.binary(&b,'*').unwrap(),'+').unwrap(),h);assert_eq!(h,a.gcd(&b).unwrap());}
        }
    }

}
