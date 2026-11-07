//! Immutable, typed exact coefficient domains and closed domain elements.
use super::algebraic::{AlgebraicScope, NumberField, AN};
use super::convert::{PResult, PolyError};
use super::f5b::Fp;
use super::fraction::{coefficient, coefficient_expr, coefficient_in_ring};
use super::int::{Coeff, Field, Rat};
use super::ratfunc::{Scope, RF};
use crate::Expr;
use std::sync::Arc;
mod complex_rounding;

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
    /// A polynomial ring (`ring`) or fraction field over a ground domain
    /// other than ZZ/QQ, as the reference `ZZ[x][y]`, `QQ(x)[y]`, `ZZ_I[x]`
    /// or `QQ<sqrt(2)>[x]`. Elements are stored flattened over
    /// `flat_parameters()`, outermost generators first, so the lex order of
    /// the flat representation is the recursive order of the reference.
    Composite {
        ground: Arc<ExactDomain>,
        ring: bool,
        parameters: Vec<Expr>,
        order: String,
    },
}
/// The innermost coefficients of a composite domain.
#[derive(Clone, Debug)]
pub enum CompositeBase {
    Rational { integral: bool },
    Algebraic { field: Arc<NumberField>, gaussian: Option<bool> },
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
            Self::Composite { ground, ring, parameters, .. } => write!(
                f,
                "{}{}{}{}",
                ground,
                if *ring { "[" } else { "(" },
                parameters.iter().map(crate::printer::to_str).collect::<Vec<_>>().join(","),
                if *ring { "]" } else { ")" }
            ),
        }
    }
}
/// `(len(p), p.terms())` of the reference `PolyElement.sort_key()`, compared
/// lexicographically: `terms()` lists `(monom, coeff)` by decreasing ring
/// order, and the tuples themselves compare monomials as plain tuples.
fn reference_poly_cmp(p: &super::MPoly<Rat>, q: &super::MPoly<Rat>, order: &str) -> Option<std::cmp::Ordering> {
    fn terms<'a>(p: &'a super::MPoly<Rat>, order: &str) -> Option<Vec<(&'a [u32], &'a Rat)>> {
        let key = |m: &[u32]| -> (u64, Vec<i64>) {
            let degree = m.iter().map(|&e| e as u64).sum();
            match order {
                "lex" => (0, m.iter().map(|&e| e as i64).collect()),
                "grlex" => (degree, m.iter().map(|&e| e as i64).collect()),
                _ => (degree, m.iter().rev().map(|&e| -(e as i64)).collect()),
            }
        };
        if !matches!(order, "lex" | "grlex" | "grevlex") {
            return None;
        }
        let mut terms: Vec<_> = p.iter().collect();
        terms.sort_by_key(|term| std::cmp::Reverse(key(term.0)));
        Some(terms)
    }
    let (a, b) = (terms(p, order)?, terms(q, order)?);
    let mut ordering = a.len().cmp(&b.len());
    for (s, t) in a.iter().zip(&b) {
        if ordering.is_ne() {
            break;
        }
        ordering = s.0.cmp(t.0).then_with(|| s.1.cmp(t.1));
    }
    Some(ordering)
}
fn failure(s: impl Into<String>) -> PolyError {
    PolyError::Failed(s.into())
}
/// The reference `FracElement` normal form of `num/den` over `ZZ_I`:
/// integral sides, `PolyElement.cancel` over the Gaussian integers, and a
/// canonical leading coefficient of the denominator.
fn gaussian_fraction_normal(num: &super::MPoly<AN>, den: &super::MPoly<AN>) -> PResult<(super::MPoly<AN>, super::MPoly<AN>)> {
    let field = super::gaussian::field();
    let _scope = AlgebraicScope::new(field.clone());
    let (a, _) = super::gaussian::clear_denoms(num);
    let (b, _) = super::gaussian::clear_denoms(den);
    let common = AN::rational(Rat::from(a.lcm(&b)), field);
    super::gaussian::expr_cancel(&num.scale(&common), &den.scale(&common))
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
    ) -> PResult<Self> {super::mpoly::checked_operation(||{
        if parameters.is_empty() {
            return Err(failure(
                "a polynomial ring or fraction field needs generators",
            ));
        }
        if super::groebner::MonomialOrder::parse(order).is_none() {
            return Err(failure("order must be lex, grlex, or grevlex"));
        }
        for (i, p) in parameters.iter().enumerate() {
            // The reference also accepts generator expressions such as
            // `ZZ[pi]` or `ZZ(sqrt(a))`, but never numbers.
            if p.as_number().is_some() || parameters[..i].contains(p) {
                return Err(failure("domain generators must be distinct nonnumeric expressions"));
            }
        }
        Ok(Self::Functions {
            integral,
            ring,
            parameters,
            order: order.into(),
        })
    })}
    /// The reference `ground.poly_ring(*parameters)` (`ring`) or
    /// `ground.frac_field(*parameters)`.
    pub fn composite(ground: &Self, ring: bool, parameters: Vec<Expr>, order: &str) -> PResult<Self> {super::mpoly::checked_operation(||{
        match ground {
            Self::ZZ | Self::QQ => return Self::functions(matches!(ground, Self::ZZ), ring, parameters, order),
            Self::Functions { .. } | Self::Composite { .. } | Self::Gaussian { .. } | Self::Algebraic(_) | Self::RR | Self::Inexact { .. } => {}
            _ => return Err(PolyError::Unsupported(
                "polynomial rings and fraction fields over this ground domain are not implemented".into(),
            )),
        }
        let Self::Functions { parameters, order, .. } = Self::functions(true, ring, parameters, order)? else {
            unreachable!()
        };
        let inner = ground.flat_parameters();
        if parameters.iter().any(|p| inner.contains(p)) {
            return Err(PolyError::Generator("polynomial ring and it's ground domain share generators".into()));
        }
        Ok(Self::Composite { ground: Arc::new(ground.clone()), ring, parameters, order })
    })}
    /// The variables of the flattened element representation: the
    /// generators of every composite level, outermost first.
    pub fn flat_parameters(&self) -> Vec<Expr> {
        match self {
            Self::Functions { parameters, .. } => parameters.clone(),
            Self::Composite { ground, parameters, .. } => {
                let mut all = parameters.clone();
                all.extend(ground.flat_parameters());
                all
            }
            _ => Vec::new(),
        }
    }
    pub fn composite_base(&self) -> Option<CompositeBase> {
        match self {
            Self::Functions { integral, .. } => Some(CompositeBase::Rational { integral: *integral }),
            Self::Composite { ground, .. } => ground.composite_base(),
            Self::Gaussian { integral } => Some(CompositeBase::Algebraic {
                field: super::gaussian::field(),
                gaussian: Some(*integral),
            }),
            Self::Algebraic(k) => Some(CompositeBase::Algebraic { field: k.clone(), gaussian: None }),
            _ => None,
        }
    }
    /// Whether every level of a composite domain is a polynomial ring.
    pub fn all_rings(&self) -> bool {
        match self {
            Self::Functions { ring, .. } => *ring,
            Self::Composite { ground, ring, .. } => *ring && ground.all_rings(),
            Self::Gaussian { .. } | Self::Algebraic(_) => true,
            _ => false,
        }
    }
    /// The flat domain with the same elements as an all-ring composite over
    /// ZZ or QQ: `ZZ[x][y]` computes in `ZZ[y,x]`.
    pub fn flat_ring(&self) -> Option<Self> {
        match (self, self.composite_base()?) {
            (Self::Composite { .. }, CompositeBase::Rational { integral }) if self.all_rings() => {
                Some(Self::Functions { integral, ring: true, parameters: self.flat_parameters(), order: "lex".into() })
            }
            _ => None,
        }
    }
    /// The flat domain with the same elements as a composite over ZZ or QQ
    /// whose levels below the outermost one are polynomial rings:
    /// `ZZ[x][y]` is `ZZ[y,x]` and `ZZ[x](y)` is `ZZ(y,x)`.
    pub fn flat_domain(&self) -> Option<Self> {
        match (self, self.composite_base()?) {
            (Self::Composite { ground, ring, .. }, CompositeBase::Rational { integral }) if ground.all_rings() => {
                Some(Self::Functions { integral, ring: *ring, parameters: self.flat_parameters(), order: "lex".into() })
            }
            _ => None,
        }
    }
    /// The same composite over QQ instead of ZZ.
    pub fn with_rational_field(&self) -> Self {
        match self {
            Self::ZZ => Self::QQ,
            Self::Functions { ring, parameters, order, .. } => Self::Functions { integral: false, ring: *ring, parameters: parameters.clone(), order: order.clone() },
            Self::Composite { ground, ring, parameters, order } => Self::Composite { ground: Arc::new(ground.with_rational_field()), ring: *ring, parameters: parameters.clone(), order: order.clone() },
            Self::Gaussian { .. } => Self::Gaussian { integral: false },
            other => other.clone(),
        }
    }
    /// Membership of a flattened rational function: each polynomial level
    /// above the innermost fraction field keeps its generators out of the
    /// reduced denominator.
    fn rational_member(&self, r: &RF, offset: usize) -> bool {
        match self {
            Self::Functions { integral, ring, .. } => !ring || coefficient_in_ring(r, *integral),
            Self::Composite { ground, ring, parameters, .. } => {
                if !ring {
                    return true;
                }
                let k = parameters.len();
                if r.d.iter().any(|(e, _)| e[offset..offset + k].iter().any(|&p| p != 0)) {
                    return false;
                }
                ground.rational_member(r, offset + k)
            }
            _ => false,
        }
    }
    /// The `RR[y]`/`CC(y)` composite ground of this domain, if any.
    pub(crate) fn inexact_ground(&self) -> Option<super::gaussian_ground::GaussianGround> {
        if !matches!(self, Self::Composite { .. }) {
            return None;
        }
        super::gaussian_ground::GaussianGround::from_domain(self).filter(|g| g.inexact.is_some())
    }
    fn composite_member(&self, value: &ElementValue) -> bool {
        match (self.composite_base(), value) {
            (Some(CompositeBase::Rational { .. }), ElementValue::Function(r)) => self.rational_member(r, 0),
            (Some(CompositeBase::Algebraic { gaussian, .. }), ElementValue::AlgebraicPoly(p)) => {
                gaussian != Some(true) || super::gaussian::poly_is_integral(p)
            }
            (_, ElementValue::GaussianFraction(..)) => self.gaussian_fraction_ground().is_some(),
            _ => false,
        }
    }
    pub fn finite(prime: u64) -> PResult<Self> {super::mpoly::checked_operation(||{
        Fp::from_rat(&Rat::zero(), prime)?;
        Ok(Self::GF(prime))
    })}
    pub fn is_field(&self) -> bool {
        !matches!(self, Self::ZZ | Self::Functions { ring: true, .. } | Self::Gaussian {integral:true} | Self::Composite { ring: true, .. })
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
            Self::Composite { ground, parameters, order, .. } => Self::Composite {
                ground: ground.clone(),
                ring: false,
                parameters: parameters.clone(),
                order: order.clone(),
            },
            _ => self.clone(),
        }
    }
    pub fn ring(&self) -> PResult<Self> {super::mpoly::checked_operation(||{
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
            Self::Composite { ground, parameters, order, .. } => Ok(Self::Composite {
                ground: ground.clone(),
                ring: true,
                parameters: parameters.clone(),
                order: order.clone(),
            }),
            _ => Err(PolyError::Unsupported(
                "this field has no implemented associated ring".into(),
            )),
        }
    })}
    /// The Gaussian ground of a fraction field directly over `ZZ_I` or
    /// `QQ_I` (`ZZ_I(y)`, `QQ_I(N,pi)`).
    pub fn gaussian_fraction_ground(&self) -> Option<super::gaussian_ground::GaussianGround> {
        if !matches!(self, Self::Composite { ring: false, .. }) {
            return None;
        }
        super::gaussian_ground::GaussianGround::from_domain(self).filter(|g| g.inexact.is_none() && !g.ring)
    }
    /// The reference `FracField` name used in conversion errors.
    pub fn fraction_field_name(&self) -> String {
        let Self::Composite { ground, parameters, order, .. } = self else { return self.to_string() };
        let symbols = parameters.iter().map(crate::printer::to_str).collect::<Vec<_>>().join(", ");
        format!("Rational function field in {symbols} over {ground} with {order} order")
    }
    /// The reference `PolyRing` name used in conversion errors.
    pub fn polynomial_ring_name(&self) -> String {
        let Self::Composite { ground, parameters, order, .. } = self else { return self.to_string() };
        let symbols = parameters.iter().map(crate::printer::to_str).collect::<Vec<_>>().join(", ");
        format!("Polynomial ring in {symbols} over {ground} with {order} order")
    }
    pub fn convert(&self, e: &Expr) -> PResult<DomainElement> {super::mpoly::checked_operation(||{
        let value=match self {
            Self::Composite { .. } if self.gaussian_fraction_ground().is_some() => {
                let g = self.gaussian_fraction_ground().unwrap();
                // FracField._rebuild_expr: an integer power of a
                // non-constant base is a raw FracElement power (over
                // `QQ<a>` the whole reference walk is `g.element`).
                if e.is_pow() && g.algebraic.is_none() {
                    let args = e.args();
                    if let Some(k) = args[1].as_i64().filter(|k| *k != 1) {
                        let base = self.convert(&args[0])?;
                        if let ElementValue::GaussianFraction(p, q) = &base.value {
                            if !(p.is_constant() && q.is_constant()) {
                                return base.pow(k);
                            }
                        }
                    }
                }
                let (p, q) = g.element(e).map_err(|_| PolyError::Reference {
                    kind: "ValueError",
                    message: format!(
                        "expected an expression convertible to a rational function in {}, got {}",
                        self.fraction_field_name(),
                        crate::printer::to_str(e)
                    ),
                })?;
                let (p, q) = if g.algebraic.is_some() { (p, q) } else { gaussian_fraction_normal(&p, &q)? };
                ElementValue::GaussianFraction(p, q)
            }
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
                // `ExpressionDomain.from_sympy`/`EXRAW.from_sympy` keep the
                // expression as it is; only the arithmetic normalizes.
                let _ = raw;
                ElementValue::Expression(e.clone())
            },
            Self::Composite { .. } if self.inexact_ground().is_some() => {
                let g = self.inexact_ground().unwrap();
                ElementValue::Inexact(super::inexact_ground::element_expr(&g, e)?)
            }
            Self::Composite { .. } => {
                let parameters = self.flat_parameters();
                let value = match self.composite_base() {
                    Some(CompositeBase::Rational { .. }) => {
                        let _scope = Scope::new(parameters.len());
                        ElementValue::Function(coefficient(e, &parameters)?)
                    }
                    Some(CompositeBase::Algebraic { field, .. }) => {
                        if !self.all_rings() {
                            return Err(PolyError::Unsupported(
                                "fraction fields over algebraic ground domains are not implemented".into(),
                            ));
                        }
                        // PolyRing.from_expr: a CoercionFailed of the
                        // `_rebuild_expr` walk is the reference ValueError.
                        let rep = super::algebraic::AlgebraicPoly::from_expr(e, &parameters, field).map_err(|_| PolyError::Reference {
                            kind: "ValueError",
                            message: format!(
                                "expected an expression convertible to a polynomial in {}, got {}",
                                self.polynomial_ring_name(),
                                crate::printer::to_str(e)
                            ),
                        })?;
                        ElementValue::AlgebraicPoly(rep.rep)
                    }
                    None => unreachable!(),
                };
                if !self.composite_member(&value) {
                    return Err(failure(format!(
                        "expected an expression convertible to a polynomial in {self}, got {}",
                        crate::printer::to_str(e)
                    )));
                }
                value
            }
        };
        Ok(DomainElement {
            domain: self.clone(),
            value,
        })
    })}
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ElementValue {
    Rational(Rat),
    Finite(Fp),
    Algebraic(AN),
    Function(RF),
    Expression(Expr),
    Inexact(Expr),
    /// A flattened polynomial over a Gaussian or algebraic base field.
    AlgebraicPoly(super::MPoly<AN>),
    /// An element `numer/denom` of `ZZ_I(y...)` or `QQ_I(y...)`, as the
    /// reference `FracElement`: both sides have Gaussian integer
    /// coefficients and are reduced by `PolyElement.cancel`, except after
    /// the operations the reference leaves raw (powers).
    GaussianFraction(super::MPoly<AN>, super::MPoly<AN>),
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
    /// `numer` and `denom` of an element of `RR(y...)` or `CC(y...)`, as
    /// expressions of the associated ring.
    pub fn inexact_fraction_parts(&self) -> Option<PResult<(Expr, Expr)>> {
        if !matches!(self.domain, ExactDomain::Composite { ring: false, .. }) {
            return None;
        }
        let g = self.domain.inexact_ground()?;
        let ElementValue::Inexact(e) = &self.value else { return None };
        Some(super::inexact_ground::element_parts(&g, e))
    }
    /// Element operations of `RR[y]`, `CC(y)`, ... (see `inexact_ground`).
    fn inexact_op(&self, other: &Self, op: &str) -> Option<PResult<Self>> {
        let g = self.domain.inexact_ground()?;
        if self.domain != other.domain {
            return None;
        }
        let (ElementValue::Inexact(a), ElementValue::Inexact(b)) = (&self.value, &other.value) else { return None };
        Some(super::inexact_ground::element_op(&g, a, b, op).map(|e| Self { domain: self.domain.clone(), value: ElementValue::Inexact(e) }))
    }
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
            ElementValue::Function(r) => coefficient_expr(r, &self.domain.flat_parameters()),
            ElementValue::GaussianFraction(p, q) => {
                self.domain.gaussian_fraction_ground().expect("Gaussian fraction field").element_expr(p, q)
            }
            ElementValue::AlgebraicPoly(p) => {
                let parameters = self.domain.flat_parameters();
                crate::ops::add(&p.iter().map(|(e, c)| {
                    let mut factors = vec![c.as_expr()];
                    for (g, &k) in parameters.iter().zip(e) {
                        if k != 0 {
                            factors.push(crate::ops::pow(g, &Expr::int(k as i64)));
                        }
                    }
                    crate::ops::mul(&factors)
                }).collect::<Vec<_>>())
            }
            // The reference real field returns the Float `0.0` for its zero.
            ElementValue::Inexact(e) if e.is_zero() && !self.domain.is_complex() && !matches!(self.domain, ExactDomain::Composite { .. }) => {
                let precision = self.domain.inexact_precision().unwrap_or(53);
                Expr::num(crate::number::Number::mp_float(crate::bigfloat::BigFloat::zero(), precision))
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
            ElementValue::AlgebraicPoly(p) => p.is_zero(),
            ElementValue::GaussianFraction(p, _) => p.is_zero(),
            ElementValue::Expression(e) | ElementValue::Inexact(e) => e.is_zero(),
        }
    }
    /// The ordering of the reference `PolyElement.sort_key()` (ZZ[x], QQ[x])
    /// or `FracElement.sort_key()` (ZZ(x), QQ(x)) between two elements of
    /// the same domain: `(len(p), p.terms())`, terms in decreasing ring
    /// order, and `(denom.sort_key(), numer.sort_key())` for fractions. The
    /// components are totally ordered, so Python's tuple comparison reduces
    /// to this lexicographic order. None outside these domains, for another
    /// domain or for a term order other than lex, grlex and grevlex.
    pub fn reference_sort_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let ExactDomain::Functions { ring, order, .. } = &self.domain else { return None };
        if self.domain != other.domain {
            return None;
        }
        let (ElementValue::Function(a), ElementValue::Function(b)) = (&self.value, &other.value) else {
            return None;
        };
        if *ring {
            let (p, q) = (self.ring_polynomial().ok()?, other.ring_polynomial().ok()?);
            return reference_poly_cmp(&p, &q, order);
        }
        match reference_poly_cmp(&a.d.to_rat(), &b.d.to_rat(), order)? {
            std::cmp::Ordering::Equal => reference_poly_cmp(&a.n.to_rat(), &b.n.to_rat(), order),
            ordering => Some(ordering),
        }
    }
    /// The coefficient that the reference `PolyElement.__eq__` compares with
    /// a foreign truthy scalar, for ZZ/QQ polynomial rings and rational
    /// function fields: `p.get(zero_monom)` when `len(p) <= 1`, and for a
    /// fraction `f.numer.get(zero_monom)` with `f.denom == 1` required as
    /// well. `Ok(None)` means the reference answer is False without any
    /// coefficient comparison, `Err(())` that the domain is not covered.
    /// The flag is the integral ground (ZZ coefficients rather than QQ).
    #[allow(clippy::result_unit_err)]
    pub fn reference_scalar_coefficient(&self) -> Result<Option<(Rat, bool, bool)>, ()> {
        let (ExactDomain::Functions { integral, ring, .. }, ElementValue::Function(r)) = (&self.domain, &self.value) else {
            return Err(());
        };
        if *ring {
            let p = self.ring_polynomial().map_err(|_| ())?;
            return Ok(match p.constant_value() {
                Some(c) if !c.is_zero() && p.len() == 1 => Some((c, *integral, true)),
                _ => None,
            });
        }
        let numer = r.n.to_rat();
        Ok(match numer.constant_value() {
            Some(c) if !c.is_zero() && numer.len() == 1 => Some((c, *integral, r.d.is_one())),
            _ => None,
        })
    }
    /// Exact multivariate polynomial view of a ring element. Rational function
    /// denominators are constant in QQ[parameters], and one in ZZ[parameters].
    pub fn ring_polynomial(&self) -> PResult<super::MPoly<Rat>> {super::mpoly::checked_operation(||{
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
    })}
    pub fn div_rem(&self, other: &Self) -> PResult<(Self, Self)> {super::mpoly::checked_operation(||{
        if self.inexact_op(other, "+").is_some() {
            if other.is_zero() {
                return Err(PolyError::Unsupported("division by zero".into()));
            }
            return Ok((self.inexact_op(other, "quo").unwrap()?, self.inexact_op(other, "rem").unwrap()?));
        }
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if other.is_zero() {
            return Err(PolyError::Unsupported("division by zero".into()));
        }
        if let Some((a, b)) = self.flat_operands(other) {
            let (q, r) = a.div_rem(&b)?;
            return Ok((self.unflatten(q)?, self.unflatten(r)?));
        }
        if let Some((a, b)) = self.algebraic_operands(other) {
            let (q, r) = a.ring_div(&b)?;
            return Ok((self.unflatten_algebraic(q.rep)?, self.unflatten_algebraic(r.rep)?));
        }
        if matches!(self.domain, ExactDomain::Composite { ring: true, .. }) {
            return Err(PolyError::Unsupported(
                "division in a polynomial ring over a fraction field is not implemented".into(),
            ));
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
    })}
    pub fn exquo(&self, other: &Self) -> PResult<Self> {super::mpoly::checked_operation(||{
        let (q, r) = self.div_rem(other)?;
        if !r.is_zero() {
            return Err(PolyError::ExactQuotientFailed(
                "polynomial or integer division is not exact in the selected ring".into(),
            ));
        }
        Ok(q)
    })}
    pub fn gcdex(&self, other: &Self) -> PResult<(Self, Self, Self)> {super::mpoly::checked_operation(||{
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
    })}
    pub fn binary(&self, other: &Self, op: char) -> PResult<Self> {super::mpoly::checked_operation(||{
        if let Some(r) = self.inexact_op(other, &op.to_string()) {
            return r;
        }
        if op == '/' && matches!(self.domain, ExactDomain::Functions { ring: true, .. } | ExactDomain::Gaussian {integral:true} | ExactDomain::Composite { ring: true, .. }) {
            return self.exquo(other);
        }
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if op == '/' && other.is_zero() {
            if matches!(self.value, ElementValue::GaussianFraction(..)) {
                return Err(PolyError::Reference { kind: "ZeroDivisionError", message: String::new() });
            }
            return Err(PolyError::Unsupported("division by zero".into()));
        }
        if let (ElementValue::GaussianFraction(a, b), ElementValue::GaussianFraction(c, d)) = (&self.value, &other.value) {
            let ground = self.domain.gaussian_fraction_ground().ok_or_else(|| failure("unexpected fraction element"))?;
            let _scope = ground.scope();
            // FracElement.__add__/__sub__/__mul__/__truediv__: `f.new(...)`,
            // which cancels; equal denominators are kept as they are, and
            // a zero operand of a sum is returned as it is.
            if matches!(op, '+' | '-') && ground.algebraic.is_some() && (a.is_zero() || c.is_zero()) {
                let value = if c.is_zero() {
                    ElementValue::GaussianFraction(a.clone(), b.clone())
                } else if op == '+' {
                    ElementValue::GaussianFraction(c.clone(), d.clone())
                } else {
                    ElementValue::GaussianFraction(c.neg(), d.clone())
                };
                return Ok(Self { domain: self.domain.clone(), value });
            }
            let (p, q) = match op {
                '+' | '-' => {
                    let sum = |x: &super::MPoly<AN>, y: &super::MPoly<AN>| if op == '+' { x.add(y) } else { x.sub(y) };
                    if b == d { (sum(a, c), b.clone()) } else { (sum(&a.mul(d), &b.mul(c)), b.mul(d)) }
                }
                '*' => (a.mul(c), b.mul(d)),
                '/' => (a.mul(d), b.mul(c)),
                _ => return Err(failure("unknown coefficient operation")),
            };
            let (p, q) = ground.cancel_pair(&p, &q)?;
            return Ok(Self { domain: self.domain.clone(), value: ElementValue::GaussianFraction(p, q) });
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
                let _scope = Scope::new(self.domain.flat_parameters().len());
                ElementValue::Function(calc!(a, b))
            }
            (ElementValue::AlgebraicPoly(a), ElementValue::AlgebraicPoly(b)) => {
                let Some(CompositeBase::Algebraic { field, .. }) = self.domain.composite_base() else {
                    unreachable!()
                };
                let _scope = AlgebraicScope::new(field);
                ElementValue::AlgebraicPoly(match op {
                    '+' => a.add(b),
                    '-' => a.sub(b),
                    '*' => a.mul(b),
                    _ => return Err(failure("unknown coefficient operation")),
                })
            }
            (ElementValue::Expression(a),ElementValue::Expression(b)) if matches!(self.domain, ExactDomain::Expression{raw:false}) => {
                // `Expression.__add__`/`__sub__`/`__mul__`/`__truediv__`:
                // zero shortcuts, the product of two numbers as it is, and
                // otherwise `Expression.simplify` of the plain result.
                let keep = |e: Expr| Ok(Self{domain:self.domain.clone(),value:ElementValue::Expression(e)});
                match op {
                    '+' if b.is_zero() => return keep(a.clone()),
                    '+' if a.is_zero() => return keep(b.clone()),
                    '-' if b.is_zero() => return keep(a.clone()),
                    '-' if a.is_zero() => return keep(crate::ops::neg(b)),
                    '*' if a.is_zero() || b.is_zero() => return keep(Expr::zero()),
                    '*' if a.as_number().is_some() && b.as_number().is_some() => return keep(crate::ops::mul2(a,b)),
                    _ => {}
                }
                let expression=match op {
                    '+'=>crate::ops::add2(a,b), '-'=>crate::ops::sub(a,b),
                    '*'=>crate::ops::mul2(a,b), '/'=>crate::ops::div(a,b),
                    _=>return Err(failure("unknown coefficient operation")),
                };
                return keep(super::expression_domain::simplify_element(&expression)?);
            },
            (ElementValue::Expression(a),ElementValue::Expression(b)) |
            (ElementValue::Inexact(a),ElementValue::Inexact(b)) => {
                if matches!(op,'*'|'/') && self.domain.is_complex() {
                    if let Some(value)=complex_rounding::binary(a,b,op,self.domain.inexact_precision().unwrap()) {
                        return Ok(Self{domain:self.domain.clone(),value:ElementValue::Inexact(value)});
                    }
                }
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
            (ExactDomain::Composite { .. }, value) if !self.domain.composite_member(value) => {
                return Err(failure("non-exact division in polynomial ring"))
            }
            _ => {}
        }
        Ok(result)
    })}
    pub fn pow(&self, exponent: i64) -> PResult<Self> {super::mpoly::checked_operation(||{
        if let ElementValue::GaussianFraction(p, q) = &self.value {
            // FracElement.__pow__ uses raw_new: no cancellation, and a
            // negative power swaps the sides without normalizing the unit.
            if exponent < 0 && p.is_zero() {
                return Err(PolyError::Reference { kind: "ZeroDivisionError", message: String::new() });
            }
            let _scope = self.domain.gaussian_fraction_ground().ok_or_else(|| failure("unexpected fraction element"))?.scope();
            let k = u32::try_from(exponent.unsigned_abs()).map_err(|_| failure("exponent too large"))?;
            let (p, q) = (p.pow(k), q.pow(k));
            let value = if exponent < 0 { ElementValue::GaussianFraction(q, p) } else { ElementValue::GaussianFraction(p, q) };
            return Ok(Self { domain: self.domain.clone(), value });
        }
        if exponent < 0 && self.is_zero() {
            return Err(PolyError::Unsupported("division by zero".into()));
        }
        if let (ExactDomain::Expression { raw: false }, ElementValue::Expression(e)) = (&self.domain, &self.value) {
            // `Expression.__pow__`: `simplify(f.ex**n)`
            let value = super::expression_domain::simplify_element(&crate::ops::pow(e, &Expr::int(exponent)))?;
            return Ok(Self { domain: self.domain.clone(), value: ElementValue::Expression(value) });
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
    })}
}

impl DomainElement {
    pub fn gaussian_parts(&self)->Option<(Rat,Rat)> {
        if !matches!(self.domain,ExactDomain::Gaussian{..}) {return None;}
        let ElementValue::Algebraic(a)=&self.value else {return None;};
        Some((a.coefficients.first().cloned().unwrap_or_else(Rat::zero),a.coefficients.get(1).cloned().unwrap_or_else(Rat::zero)))
    }
    pub fn gaussian_canonical_unit(&self)->PResult<Self> {super::mpoly::checked_operation(||{
        let (x,y)=self.gaussian_parts().ok_or_else(||failure("expected a Gaussian element"))?;
        let quadrant=if !y.is_zero()&&!y.is_negative() {if !x.is_zero()&&!x.is_negative(){0}else{1}}
            else if y.is_negative(){if x.is_negative(){2}else{3}}
            else if x.is_negative(){2}else{0};
        self.domain.convert(&[Expr::one(),crate::ops::neg(&Expr::i()),Expr::int(-1),Expr::i()][quadrant])
    })}
    pub fn gaussian_fraction(&self)->PResult<(Self,Self)> {super::mpoly::checked_operation(||{
        let (x,y)=self.gaussian_parts().ok_or_else(||failure("expected a Gaussian element"))?;
        let d=x.d.lcm(&y.d);
        let numerator=crate::ops::mul2(&self.as_expr(),&Expr::num(d.to_number()));
        let ring=ExactDomain::Gaussian{integral:true};
        Ok((ring.convert(&crate::expand::expand(&numerator))?,ring.convert(&Expr::num(d.to_number()))?))
    })}
}

/// Recover the stored real value exactly, optionally limiting its denominator
/// by the real-field convention. Continued fractions run on native big ints.
pub fn rationalize_real(e:&Expr, precision:u32, limit:bool)->PResult<Rat> {super::mpoly::checked_operation(||{
    use num_bigint::{BigInt,Sign};
    use num_rational::BigRational;
    let value=crate::evalf_mp::eval_mp(e,precision).ok_or_else(||failure("expected a finite numerical value"))?;
    if !value.is_real(){return Err(failure("a nonreal value has no rational conversion"));}
    // Compact literals can now be evaluated beyond the transcendental engine's
    // former magnitude guard. Retain that pre-existing output allocation
    // budget for a materialized integer ratio; widening it is a separate task.
    if !value.re.is_zero()&&value.re.mag().unsigned_abs()>3_321_930{
        return Err(failure("exact rational output exceeds the current conversion resource budget"));
    }
    let exact=value.re.to_ratio_exact();
    if !limit{return Ok(Rat::from_number(&crate::number::Number::from_ratio(exact)).unwrap());}
    let max_denominator=((BigInt::from(1)<<precision as usize)/BigInt::from(200)).max(BigInt::from(99));
    let result=if exact.denom()<=&max_denominator {exact} else {
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
})}

impl DomainElement {
    /// The operands of an all-ring composite over ZZ/QQ in its flat ring.
    fn flat_operands(&self, other: &Self) -> Option<(Self, Self)> {
        let flat = self.domain.flat_ring()?;
        let wrap = |e: &Self| Self { domain: flat.clone(), value: e.value.clone() };
        Some((wrap(self), wrap(other)))
    }
    fn unflatten(&self, flat: Self) -> PResult<Self> {
        if !self.domain.composite_member(&flat.value) {
            return Err(failure("result is outside the selected polynomial ring"));
        }
        Ok(Self { domain: self.domain.clone(), value: flat.value })
    }
    /// The operands of a composite over a Gaussian or algebraic field as
    /// native algebraic polynomials over the flat generators.
    fn algebraic_operands(&self, other: &Self) -> Option<(super::algebraic::AlgebraicPoly, super::algebraic::AlgebraicPoly)> {
        let Some(CompositeBase::Algebraic { field, gaussian }) = self.domain.composite_base() else {
            return None;
        };
        let (ElementValue::AlgebraicPoly(a), ElementValue::AlgebraicPoly(b)) = (&self.value, &other.value) else {
            return None;
        };
        let wrap = |rep: &super::MPoly<AN>| super::algebraic::AlgebraicPoly {
            rep: rep.clone(),
            gens: self.domain.flat_parameters(),
            number_field: field.clone(),
            gaussian,
        };
        Some((wrap(a), wrap(b)))
    }
    fn unflatten_algebraic(&self, rep: super::MPoly<AN>) -> PResult<Self> {
        let value = ElementValue::AlgebraicPoly(rep);
        if !self.domain.composite_member(&value) {
            return Err(failure("result is outside the selected polynomial ring"));
        }
        Ok(Self { domain: self.domain.clone(), value })
    }
    pub fn gcd(&self, other: &Self) -> PResult<Self> {super::mpoly::checked_operation(||{
        if let Some(r) = self.inexact_op(other, "gcd") {
            return r;
        }
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if let Some((a, b)) = self.flat_operands(other) {
            return self.unflatten(a.gcd(&b)?);
        }
        if let Some((a, b)) = self.algebraic_operands(other) {
            return self.unflatten_algebraic(a.gcd(&b).rep);
        }
        if let (ElementValue::GaussianFraction(a, b), ElementValue::GaussianFraction(c, d)) = (&self.value, &other.value) {
            let ground = self.domain.gaussian_fraction_ground().ok_or_else(|| failure("unexpected fraction element"))?;
            let _scope = ground.scope();
            // Field.gcd: ring.gcd of the numerators over ring.lcm of
            // the denominators.
            let ring = self.domain.ring()?;
            let wrap = |p: &super::MPoly<AN>| Self { domain: ring.clone(), value: ElementValue::AlgebraicPoly(p.clone()) };
            let (ElementValue::AlgebraicPoly(p), ElementValue::AlgebraicPoly(q)) =
                (wrap(a).gcd(&wrap(c))?.value, wrap(b).lcm(&wrap(d))?.value)
            else {
                return Err(failure("unexpected Gaussian ring element"));
            };
            let one = super::MPoly::one(p.nvars);
            let normal = |p: &super::MPoly<AN>| if ground.algebraic.is_some() { ground.cancel_pair(p, &one) } else { gaussian_fraction_normal(p, &one) };
            let (p, p1) = normal(&p)?;
            let (q, q1) = normal(&q)?;
            let numerator = Self { domain: self.domain.clone(), value: ElementValue::GaussianFraction(p, p1) };
            let denominator = Self { domain: self.domain.clone(), value: ElementValue::GaussianFraction(q, q1) };
            return numerator.binary(&denominator, '/');
        }
        if let ExactDomain::Composite { ring, .. } = &self.domain {
            if *ring || !matches!(self.domain.composite_base(), Some(CompositeBase::Rational { .. })) {
                return Err(PolyError::Unsupported(
                    "gcd in this composite domain is not implemented".into(),
                ));
            }
            let flat = Self::wrap_flat_field(&self.domain);
            let a = Self { domain: flat.clone(), value: self.value.clone() };
            let b = Self { domain: flat, value: other.value.clone() };
            return Ok(Self { domain: self.domain.clone(), value: a.gcd(&b)?.value });
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
    })}
    fn wrap_flat_field(domain: &ExactDomain) -> ExactDomain {
        let Some(CompositeBase::Rational { integral }) = domain.composite_base() else { unreachable!() };
        ExactDomain::Functions { integral, ring: false, parameters: domain.flat_parameters(), order: "lex".into() }
    }
    pub fn lcm(&self, other: &Self) -> PResult<Self> {super::mpoly::checked_operation(||{
        if let Some(r) = self.inexact_op(other, "lcm") {
            return r;
        }
        if !self.compatible_embedding(other) {
            return Err(failure("incompatible coefficient domains"));
        }
        if let Some((a, b)) = self.flat_operands(other) {
            return self.unflatten(a.lcm(&b)?);
        }
        if let Some((a, b)) = self.algebraic_operands(other) {
            if a.rep.is_zero() || b.rep.is_zero() {
                return self.domain.convert(&Expr::zero());
            }
            return self.unflatten_algebraic(a.lcm(&b).rep);
        }
        if let (ElementValue::GaussianFraction(a, b), ElementValue::GaussianFraction(c, d)) = (&self.value, &other.value) {
            let ground = self.domain.gaussian_fraction_ground().ok_or_else(|| failure("unexpected fraction element"))?;
            let _scope = ground.scope();
            // Field.lcm: ring.lcm of the numerators over ring.gcd of
            // the denominators.
            let ring = self.domain.ring()?;
            let wrap = |p: &super::MPoly<AN>| Self { domain: ring.clone(), value: ElementValue::AlgebraicPoly(p.clone()) };
            let (ElementValue::AlgebraicPoly(p), ElementValue::AlgebraicPoly(q)) =
                (wrap(a).lcm(&wrap(c))?.value, wrap(b).gcd(&wrap(d))?.value)
            else {
                return Err(failure("unexpected Gaussian ring element"));
            };
            let one = super::MPoly::one(p.nvars);
            let normal = |p: &super::MPoly<AN>| if ground.algebraic.is_some() { ground.cancel_pair(p, &one) } else { gaussian_fraction_normal(p, &one) };
            let (p, p1) = normal(&p)?;
            let (q, q1) = normal(&q)?;
            let numerator = Self { domain: self.domain.clone(), value: ElementValue::GaussianFraction(p, p1) };
            let denominator = Self { domain: self.domain.clone(), value: ElementValue::GaussianFraction(q, q1) };
            return numerator.binary(&denominator, '/');
        }
        if let ExactDomain::Composite { ring, .. } = &self.domain {
            if *ring || !matches!(self.domain.composite_base(), Some(CompositeBase::Rational { .. })) {
                return Err(PolyError::Unsupported(
                    "lcm in this composite domain is not implemented".into(),
                ));
            }
            let flat = Self::wrap_flat_field(&self.domain);
            let a = Self { domain: flat.clone(), value: self.value.clone() };
            let b = Self { domain: flat, value: other.value.clone() };
            return Ok(Self { domain: self.domain.clone(), value: a.lcm(&b)?.value });
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
    })}
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
