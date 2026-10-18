//! Exact polynomials over rational function coefficient fields.

use super::convert::{from_exprs_gens, sort_gens, to_expr, to_expr_int, PResult, PolyError};
use super::fieldalg;
use super::groebner::MonomialOrder;
use super::int::{Coeff, Field, Int, Rat};
use super::mpoly::MPoly;
use super::polyobj::Poly;
use super::ratfunc::{Scope, RF};
use crate::expr::{Expr, Kind};
use crate::ops::{add, div, mul, pow};
use std::collections::BTreeMap;

// Dense polynomial division recursively in the declared main generators.
// A coefficient-ring obstruction stops division at that leading coefficient;
// it must not be confused with sparse multivariate ring-element reduction.
pub(crate) fn recursive_division<C: Field>(f: &MPoly<C>, g: &MPoly<C>, v: usize, allow: &impl Fn(&C)->bool) -> (MPoly<C>, MPoly<C>) {
    let n=f.nvars;
    if v>=n || (f.is_constant() && g.is_constant()) {
        if let (Some(a),Some(b))=(f.constant_value(),g.constant_value()) {
            if !b.is_zero() { let c=a.div(&b); if allow(&c) {return (MPoly::constant(c,n),MPoly::zero(n));} }
        }
        return (MPoly::zero(n),f.clone());
    }
    let dg=g.degree(v); let mut q=MPoly::zero(n); let mut r=f.clone();
    if r.is_zero() || r.degree(v)<dg {return (q,r);}
    let lcg=g.lc_in(v);
    loop {
        let dr=r.degree(v);
        let (c,rr)=recursive_division(&r.lc_in(v),&lcg,v+1,allow);
        if !rr.is_zero() || c.is_zero(){break;}
        let mut e=vec![0;n];e[v]=(dr-dg) as u32;
        let t=c.mul_term(&e,&C::one()); q=q.add(&t);r=r.sub(&t.mul(g));
        if r.is_zero()||r.degree(v)<dg||r.degree(v)>=dr{break;}
    }
    (q,r)
}

#[derive(Clone, Debug)]
pub struct FractionPoly {
    pub rep: MPoly<RF>,
    pub gens: Vec<Expr>,
    pub parameters: Vec<Expr>,
    pub display: String,
    pub coefficient_domain: super::domains::ExactDomain,
}

pub fn parse_domain(domain: &str, available: &[Expr]) -> PResult<Option<Vec<Expr>>> {
    if domain.starts_with("QQ<") {return Ok(None);}
    let Some(open) = domain.find(['(', '[']) else {
        return Ok(None);
    };
    let closing = if domain.as_bytes()[open] == b'(' {
        ')'
    } else {
        ']'
    };
    if !matches!(&domain[..open], "QQ" | "ZZ") || !domain.ends_with(closing) {
        return Err(PolyError::Unsupported(format!(
            "unsupported coefficient domain {domain}"
        )));
    }
    let mut parameters = Vec::new();
    for name in domain[open + 1..domain.len() - 1].split(',').map(str::trim) {
        if name.is_empty() {
            return Err(PolyError::Failed(
                "a fraction field needs named parameters".into(),
            ));
        }
        let candidates: Vec<_> = available
            .iter()
            .filter(|s| s.as_symbol().is_some_and(|s| s == name))
            .collect();
        if candidates.len() > 1 {
            return Err(PolyError::Failed(format!(
                "ambiguous parameter assumptions for {name}"
            )));
        }
        let value = candidates
            .first()
            .map(|s| (*s).clone())
            .unwrap_or_else(|| Expr::symbol(name));
        if !matches!(value.kind(), Kind::Sym(_)) || parameters.contains(&value) {
            return Err(PolyError::Failed(
                "fraction field parameters must be distinct symbols".into(),
            ));
        }
        // A parsed name must itself be a symbol, not a hidden expression.
        let parsed = crate::parser::parse_str(name)
            .map_err(|_| PolyError::Failed("invalid fraction field parameter".into()))?;
        if !matches!(parsed.kind(), Kind::Sym(_)) {
            return Err(PolyError::Failed(
                "fraction field parameters must be symbols".into(),
            ));
        }
        parameters.push(value);
    }
    Ok(Some(parameters))
}

pub fn infer_parameters(expressions: &[Expr], gens: &[Expr]) -> Vec<Expr> {
    let mut parameters = Vec::new();
    for e in expressions {
        for s in crate::traverse::free_symbols(e) {
            if !gens.contains(&s) && !parameters.contains(&s) {
                parameters.push(s);
            }
        }
    }
    sort_gens(&mut parameters);
    parameters
}

pub fn domain_name(parameters: &[Expr]) -> String {
    format!(
        "QQ({})",
        parameters
            .iter()
            .map(crate::printer::to_str)
            .collect::<Vec<_>>()
            .join(",")
    )
}

pub fn coefficient_expr(c: &RF, parameters: &[Expr]) -> Expr {
    div(
        &to_expr_int(&c.n, parameters),
        &to_expr_int(&c.d, parameters),
    )
}

pub fn coefficient(value: &Expr, parameters: &[Expr]) -> PResult<RF> {
    let value = crate::simplify::together(value);
    let (num, den) = crate::simplify::numer_denom(&value);
    let converted = from_exprs_gens(&[num, den], parameters)?;
    if converted.float || converted.gens.len() != parameters.len() {
        return Err(PolyError::Unsupported(
            "coefficients must be exact rational functions of the declared parameters".into(),
        ));
    }
    let (nd, n) = converted.polys[0].clear_denoms();
    let (dd, d) = converted.polys[1].clear_denoms();
    if d.is_zero() {
        return Err(PolyError::Failed("zero coefficient denominator".into()));
    }
    Ok(RF::new(n.scale(&dd), d.scale(&nd)))
}

impl FractionPoly {
    pub fn from_expr(e: &Expr, gens: &[Expr], parameters: &[Expr], display: &str) -> PResult<Self> {
        if gens.is_empty() || gens.iter().any(|g| parameters.contains(g)) {
            return Err(PolyError::Failed(
                "polynomial generators and coefficient parameters must be disjoint".into(),
            ));
        }
        if gens
            .iter()
            .enumerate()
            .any(|(i, g)| g.is_number() || gens[..i].contains(g))
        {
            return Err(PolyError::Failed(
                "polynomial generators must be distinct nonnumeric expressions".into(),
            ));
        }
        let _scope = Scope::new(parameters.len());
        let converted = from_exprs_gens(std::slice::from_ref(e), gens)?;
        if converted.float {
            return Err(PolyError::Unsupported(
                "fraction fields require exact coefficients".into(),
            ));
        }
        let mut grouped: BTreeMap<Vec<u32>, Vec<(Vec<u32>, Rat)>> = BTreeMap::new();
        for (exponent, c) in converted.polys[0].iter() {
            grouped
                .entry(exponent[..gens.len()].to_vec())
                .or_default()
                .push((exponent[gens.len()..].to_vec(), c.clone()));
        }
        let coefficient_gens = &converted.gens[gens.len()..];
        let mut terms = Vec::new();
        for (exponent, values) in grouped {
            let value = to_expr(
                &MPoly::from_terms(coefficient_gens.len(), values),
                coefficient_gens,
            );
            let c = coefficient(&value, parameters)?;
            terms.push((exponent, c));
        }
        Ok(Self {
            rep: MPoly::from_terms(gens.len(), terms),
            gens: gens.to_vec(),
            parameters: parameters.to_vec(),
            display: display.into(),
            coefficient_domain: super::domains::ExactDomain::functions(display.starts_with("ZZ"),display.ends_with(']'),parameters.to_vec(),"lex")?,
        })
    }
    pub fn with_rep(&self, rep: MPoly<RF>) -> Self {
        Self {
            rep,
            ..self.clone()
        }
    }
    pub fn as_expr(&self) -> Expr {
        add(&self
            .rep
            .iter()
            .map(|(e, c)| {
                let mut factors = vec![coefficient_expr(c, &self.parameters)];
                for (g, &k) in self.gens.iter().zip(e) {
                    if k != 0 {
                        factors.push(pow(g, &Expr::int(k as i64)));
                    }
                }
                mul(&factors)
            })
            .collect::<Vec<_>>())
    }
    pub fn terms(&self) -> Vec<(Vec<u32>, Expr)> {
        self.rep
            .iter()
            .map(|(e, c)| (e.to_vec(), coefficient_expr(c, &self.parameters)))
            .collect()
    }
    pub fn legacy(&self) -> Poly {
        Poly::from_expr(&self.as_expr(), &self.gens)
            .expect("field polynomial is polynomial in its generators")
    }
    pub fn convert(&self, e: &Expr) -> PResult<Self> {
        Self::from_expr(e, &self.gens, &self.parameters, &self.display).map(|mut p| {p.coefficient_domain=self.coefficient_domain.clone();p})
    }
    pub fn promote_field(&mut self) {
        self.display=self.display.replace('[',"(").replace(']',")");
        self.coefficient_domain=self.coefficient_domain.field();
    }
    pub fn belongs_to_ring(&self) -> bool {
        self.rep
            .coeffs()
            .iter()
            .all(|c| coefficient_in_ring(c, self.display.starts_with("ZZ")))
    }
    pub fn add(&self, other: &Self) -> Self {
        let _s = Scope::new(self.parameters.len());
        self.with_rep(self.rep.add(&other.rep))
    }
    pub fn sub(&self, other: &Self) -> Self {
        let _s = Scope::new(self.parameters.len());
        self.with_rep(self.rep.sub(&other.rep))
    }
    pub fn mul(&self, other: &Self) -> Self {
        let _s = Scope::new(self.parameters.len());
        self.with_rep(self.rep.mul(&other.rep))
    }
    pub fn neg(&self) -> Self {
        self.with_rep(self.rep.neg())
    }
    pub fn pow(&self, k: u32) -> Self {
        let _s = Scope::new(self.parameters.len());
        self.with_rep(self.rep.pow(k))
    }
    pub fn monic(&self) -> Self {
        let _s = Scope::new(self.parameters.len());
        let mut p = self.with_rep(fieldalg::monic(&self.rep, MonomialOrder::Lex));
        if p.display.ends_with(']') && !p.belongs_to_ring() {
            p.promote_field();
        }
        p
    }
    pub fn div(&self,other:&Self)->PResult<(Self,Self)>{self.div_with_options(other,true)}
    pub fn div_with_options(&self, other: &Self, auto:bool) -> PResult<(Self, Self)> {
        if other.rep.is_zero() {
            return Err(PolyError::Unsupported("polynomial division by zero".into()));
        }
        let _s = Scope::new(self.parameters.len());
        let (q,r)=recursive_division(&self.rep,&other.rep,0,&|c|auto||!self.display.ends_with(']')||coefficient_in_ring(c,self.display.starts_with("ZZ")));
        let mut q=self.with_rep(q);
        let mut r = self.with_rep(r);
        if self.display.ends_with(']') && (!q.belongs_to_ring() || !r.belongs_to_ring()) {
            q.promote_field();
            r.promote_field();
        }
        Ok((q, r))
    }
    pub fn diff(&self, i: usize) -> Self {
        let _s = Scope::new(self.parameters.len());
        self.with_rep(self.rep.derivative(i))
    }
    pub fn integrate(&self, i: usize) -> Self {
        let _s = Scope::new(self.parameters.len());
        let terms = self
            .rep
            .iter()
            .map(|(e, c)| {
                let mut e = e.to_vec();
                e[i] += 1;
                let c = c.div(&RF::from_i64(e[i] as i64));
                (e, c)
            })
            .collect();
        let mut p=self.with_rep(MPoly::from_terms(self.gens.len(), terms));
        p.promote_field();
        p
    }
    pub fn gcd(&self, other: &Self) -> Self {
        let _s = Scope::new(self.parameters.len());
        // Gauss' lemma: clear coefficient denominators and compute the gcd
        // in QQ[parameters,generators], then discard coefficient units.
        let a = self.clear_coeff_denominators();
        let b = other.clear_coeff_denominators();
        let (g, _, _) = super::gcd::gcd_mz(&a, &b);
        let mut allgens = self.gens.clone();
        allgens.extend_from_slice(&self.parameters);
        let p = self
            .convert(&to_expr_int(&g, &allgens))
            .expect("gcd polynomial conversion")
            .monic();
        if self.display.ends_with(']') {
            if p.rep.is_zero(){return self.with_rep(p.rep);}
            let content=self.coefficient_content().gcd(&other.coefficient_content()).expect("valid coefficient-ring gcd");
            let mut primitive=p.with_rep(p.rep.scale(&RF::from_poly(p.denominator())));
            primitive.display=self.display.clone();primitive.coefficient_domain=self.coefficient_domain.clone();
            let divisor=primitive.coefficient_content();
            let (super::domains::ElementValue::Function(content),super::domains::ElementValue::Function(divisor))=(content.value,divisor.value)else{unreachable!()};
            self.with_rep(primitive.rep.scale(&content.div(&divisor)))
        } else {
            p
        }
    }
    pub fn lcm(&self,other:&Self)->Self {
        let _scope=Scope::new(self.parameters.len());
        if self.rep.is_zero()||other.rep.is_zero(){return self.with_rep(MPoly::zero(self.gens.len()));}
        if !self.display.ends_with(']') {return self.mul(other).div(&self.gcd(other)).expect("polynomial gcd divides the product").0.monic();}
        // Gauss' lemma over the coefficient ring: remove each input's
        // content first, form the primitive lcm, then restore lcm(contents).
        // Normalizing the product/gcd afterwards loses unit factors when
        // the coefficient domain itself has noncanonical gcd associates.
        let a=self.coefficient_content();let b=other.coefficient_content();
        let target=a.lcm(&b).expect("valid coefficient-ring lcm");
        let(super::domains::ElementValue::Function(a),super::domains::ElementValue::Function(b),super::domains::ElementValue::Function(target))=(a.value,b.value,target.value)else{unreachable!()};
        let a=self.with_rep(self.rep.scale(&a.inv()));let b=other.with_rep(other.rep.scale(&b.inv()));
        let product=a.mul(&b);let gcd=a.gcd(&b);
        let q=product.div_with_options(&gcd,false).expect("primitive gcd divides the product in its ring").0;
        self.with_rep(q.rep.scale(&target))
    }
    fn coefficient_content(&self)->super::domains::DomainElement {
        use super::domains::{DomainElement,ElementValue};
        self.rep.coeffs().iter().fold(DomainElement{domain:self.coefficient_domain.clone(),value:ElementValue::Function(RF::zero())},|g,c|g.gcd(&DomainElement{domain:self.coefficient_domain.clone(),value:ElementValue::Function(c.clone())}).expect("valid polynomial-ring coefficient gcd"))
    }
    pub fn clear_coeff_denominators(&self) -> MPoly<Int> {
        let _s = Scope::new(self.parameters.len());
        let denominator =
            self.rep
                .coeffs()
                .iter()
                .fold(MPoly::one(self.parameters.len()), |a, c| {
                    let (g, _, _) = super::gcd::gcd_mz(&a, &c.d);
                    a.div_exact(&g).unwrap().mul(&c.d)
                });
        let mut terms = Vec::new();
        for (e, c) in self.rep.iter() {
            let numerator = c.n.mul(&denominator.div_exact(&c.d).unwrap());
            for (p, n) in numerator.iter() {
                let mut powers = e.to_vec();
                powers.extend_from_slice(p);
                terms.push((powers, n.clone()));
            }
        }
        MPoly::from_terms(self.gens.len() + self.parameters.len(), terms)
    }
    pub fn ring_normalize(&self) -> Self {
        if self.rep.is_zero() {
            return self.clone();
        }
        let p = self.clear_coeff_denominators();
        let p = p.div_coeff(&p.content()).unwrap();
        let mut allgens = self.gens.clone();
        allgens.extend_from_slice(&self.parameters);
        self.convert(&to_expr_int(&p, &allgens))
            .expect("primitive coefficient polynomial")
    }
    pub fn denominator(&self) -> MPoly<Int> {
        self.rep
            .coeffs()
            .iter()
            .fold(MPoly::one(self.parameters.len()), |a, c| {
                let (g, _, _) = super::gcd::gcd_mz(&a, &c.d);
                a.div_exact(&g).unwrap().mul(&c.d)
            })
    }
    pub fn primitive(&self) -> (RF, Self) {
        let _s = Scope::new(self.parameters.len());
        if self.rep.is_zero() {
            return (RF::zero(), self.clone());
        }
        let denominator = self.denominator();
        let mut content = MPoly::zero(self.parameters.len());
        for c in self.rep.coeffs() {
            let n = c.n.mul(&denominator.div_exact(&c.d).unwrap());
            content = super::gcd::gcd_mz(&content, &n).0;
        }
        let c = RF::new(content, denominator);
        (c.clone(), self.with_rep(self.rep.scale(&c.inv())))
    }
    pub fn factor_list(&self) -> (RF, Vec<(Self, u32)>) {
        let _s = Scope::new(self.parameters.len());
        if self.rep.is_zero() {
            return (RF::zero(), Vec::new());
        }
        let cleared = self.clear_coeff_denominators();
        let (c, factors) = super::factor_multi::factor_list_mz(&cleared);
        let mut unit = RF::new(
            MPoly::constant(c, self.parameters.len()),
            self.denominator(),
        );
        let mut output = Vec::new();
        let mut allgens = self.gens.clone();
        allgens.extend_from_slice(&self.parameters);
        for (f, k) in factors {
            let p = self
                .convert(&to_expr_int(&f, &allgens))
                .expect("fraction field factor");
            if p.rep.is_constant() {
                let c = p.rep.constant_value().unwrap();
                for _ in 0..k {
                    unit = unit.mul(&c);
                }
            } else {
                output.push((p, k));
            }
        }
        (unit, output)
    }
    pub fn sqf_list(&self) -> (RF, Vec<(Self, u32)>) {
        let _s = Scope::new(self.parameters.len());
        let (mut unit, factors) = self.factor_list();
        let mut groups: BTreeMap<u32, Self> = BTreeMap::new();
        for (f, k) in factors {
            let c = f.rep.lc();
            for _ in 0..k {
                unit = unit.mul(&c);
            }
            let f = f.monic();
            groups.entry(k).and_modify(|a| *a = a.mul(&f)).or_insert(f);
        }
        (unit, groups.into_iter().map(|(k, p)| (p, k)).collect())
    }
    pub fn compose(&self, other: &Self) -> Self {
        let _s = Scope::new(self.parameters.len());
        let mut result = MPoly::zero(self.gens.len());
        for (e, c) in self.rep.iter() {
            let mut m = e.to_vec();
            m[0] = 0;
            let term = MPoly::monomial(&m, c.clone()).mul(&other.rep.pow(e[0]));
            result = result.add(&term);
        }
        self.with_rep(result)
    }
}

pub fn coefficient_in_ring(c: &RF, integral: bool) -> bool {
    if integral {
        c.d.is_one()
    } else {
        c.d.is_constant()
    }
}

/// A storage adapter keeps existing QQ/ZZ optimized representations intact.
/// Fraction-domain operations always use the RF polynomial; its ordinary
/// Poly view is used only by legacy metadata consumers.
#[derive(Clone, Debug)]
pub struct Polynomial {
    pub base: Poly,
    pub field: Option<FractionPoly>,
    pub algebraic: Option<super::algebraic::AlgebraicPoly>,
}
impl From<Poly> for Polynomial {
    fn from(base: Poly) -> Self {
        Self { base, field: None, algebraic:None }
    }
}
impl From<FractionPoly> for Polynomial {
    fn from(field: FractionPoly) -> Self {
        Self {
            base: field.legacy(),
            field: Some(field),
            algebraic:None,
        }
    }
}
impl From<super::algebraic::AlgebraicPoly> for Polynomial {
    fn from(p:super::algebraic::AlgebraicPoly)->Self {Self {base:p.legacy(),field:None,algebraic:Some(p)}}
}
impl std::ops::Deref for Polynomial {
    type Target = Poly;
    fn deref(&self) -> &Poly {
        &self.base
    }
}
impl std::ops::DerefMut for Polynomial {
    fn deref_mut(&mut self) -> &mut Poly {
        &mut self.base
    }
}
impl Polynomial {
    pub fn decompose(&self)->PResult<Vec<Self>> {
        if let Some(p)=&self.algebraic {
            super::roots::decompose_preflight(p.gens.len(),p.rep.degree(0))?;
            let _scope=super::algebraic::AlgebraicScope::new(p.number_field.clone());
            let coefficients=super::dense::from_mpoly(&p.rep,0);
            let parts=super::roots::decompose_field(&coefficients,&|_|true)?;
            return Ok(parts.iter().map(|part|p.with_rep(super::dense::to_mpoly(part,0,1)).into()).collect());
        }
        if let Some(p)=&self.field {
            super::roots::decompose_preflight(p.gens.len(),p.rep.degree(0))?;
            let _scope=Scope::new(p.parameters.len());
            let coefficients=super::dense::from_mpoly(&p.rep,0);
            let ring=p.display.ends_with(']');
            let parts=super::roots::decompose_field(&coefficients,&|value|!ring || coefficient_in_ring(value,p.display.starts_with("ZZ")))?;
            return Ok(parts.iter().map(|part|p.with_rep(super::dense::to_mpoly(part,0,1)).into()).collect());
        }
        Ok(self.base.decompose()?.into_iter().map(Into::into).collect())
    }

    pub fn as_expr(&self) -> Expr {
        if let Some(p)=&self.algebraic {return p.as_expr();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.as_expr(), FractionPoly::as_expr)
    }
    pub fn domain_str(&self) -> String {
        if let Some(p)=&self.algebraic {return p.number_field.display();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.domain_str(), |p| p.display.clone())
    }
    pub fn terms(&self) -> Vec<(Vec<u32>, Expr)> {
        if let Some(p)=&self.algebraic {return p.terms();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.terms(), FractionPoly::terms)
    }
    pub fn lc(&self) -> Expr {
        self.terms()
            .into_iter()
            .next()
            .map(|(_, c)| c)
            .unwrap_or_else(Expr::zero)
    }
    pub fn tc(&self) -> Expr {
        self.terms()
            .into_iter()
            .find(|(e, _)| e.iter().all(|&x| x == 0))
            .map(|(_, c)| c)
            .unwrap_or_else(Expr::zero)
    }
    pub fn coeffs(&self) -> Vec<Expr> {
        self.terms().into_iter().map(|(_, c)| c).collect()
    }
    pub fn coeff_monomial(&self, e: &[u32]) -> Expr {
        self.terms()
            .into_iter()
            .find(|(m, _)| m == e)
            .map(|(_, c)| c)
            .unwrap_or_else(Expr::zero)
    }
    pub fn all_coeffs(&self) -> PResult<Vec<Expr>> {
        Ok(super::structure::dense_terms(self)?.into_iter().map(|(_,c)|c).collect())
    }
    /// Checked coercion boundary for public fallible operations. Coefficients
    /// from incompatible fields never enter arithmetic or an expect here.
    pub fn checked_algebraic_pair(&self,o:&Self)->PResult<Option<(super::algebraic::AlgebraicPoly,super::algebraic::AlgebraicPoly)>> {
        let Some(prototype)=self.algebraic.as_ref().or(o.algebraic.as_ref()) else{return Ok(None)};
        if self.field.is_some() || o.field.is_some() || self.modulus.is_some() || o.modulus.is_some() {
            return Err(PolyError::Unsupported("cannot combine algebraic, rational-function, and finite coefficient fields".into()));
        }
        let gens=super::polyobj::unify_gens(self.true_gens(),o.true_gens());
        let field=if let(Some(a),Some(b))=(&self.algebraic,&o.algebraic){
            if a.number_field.same_embedding(&b.number_field) {a.number_field.clone()}else{
                super::algebraic::NumberField::from_generators(&[a.number_field.generator.clone(),b.number_field.generator.clone()])?
            }
        }else{prototype.number_field.clone()};
        let convert=|p:&Self| {
            if let Some(a)=&p.algebraic {if a.number_field==field && a.gens==gens {return Ok(a.clone());}}
            super::algebraic::AlgebraicPoly::from_expr(&p.as_expr(),&gens,field.clone())
        };
        Ok(Some((convert(self)?,convert(o)?)))
    }
    fn algebraic_pair(&self,o:&Self)->Option<(super::algebraic::AlgebraicPoly,super::algebraic::AlgebraicPoly)> {
        self.checked_algebraic_pair(o).expect("compatible algebraic fields")
    }
    pub(crate) fn checked_fraction_pair(&self, o: &Self) -> PResult<Option<(FractionPoly, FractionPoly)>> {
        let Some(prototype)=self.field.as_ref().or(o.field.as_ref()) else{return Ok(None)};
        if self.modulus.is_some()||o.modulus.is_some()||self.algebraic.is_some()||o.algebraic.is_some(){return Err(PolyError::Unsupported("cannot combine rational-function, algebraic, and finite coefficient fields".into()));}
        if let(Some(a),Some(b))=(&self.field,&o.field){if let(super::domains::ExactDomain::Functions{order:a,..},super::domains::ExactDomain::Functions{order:b,..})=(&a.coefficient_domain,&b.coefficient_domain){if a!=b{return Err(PolyError::Unsupported("cannot unify coefficient rings with different monomial orders".into()));}}}
        let mut gens = super::polyobj::unify_gens(self.true_gens(), o.true_gens());
        if gens.is_empty() {
            gens = prototype.gens.clone();
        }
        let mut parameters = prototype.parameters.clone();
        if let Some(f) = &o.field {
            for p in &f.parameters {
                if !parameters.contains(p) {
                    parameters.push(p.clone());
                }
            }
        }
        if let Some(f) = &self.field {
            for p in &f.parameters {
                if !parameters.contains(p) {
                    parameters.push(p.clone());
                }
            }
        }
        for p in infer_parameters(&[self.as_expr(), o.as_expr()], &gens) {
            if !parameters.contains(&p) {
                parameters.push(p);
            }
        }
        if parameters != prototype.parameters {
            sort_gens(&mut parameters);
        }
        let integral=[self,o].iter().all(|p|if let Some(f)=&p.field {matches!(f.coefficient_domain,super::domains::ExactDomain::Functions{integral:true,..})}else{p.domain==super::polyobj::Domain::ZZ});
        let ring=[self,o].iter().all(|p|p.field.as_ref().is_none_or(|f|matches!(f.coefficient_domain,super::domains::ExactDomain::Functions{ring:true,..})));
        let order=if let super::domains::ExactDomain::Functions{order,..}=&prototype.coefficient_domain{order.as_str()}else{"lex"};
        let descriptor=super::domains::ExactDomain::functions(integral,ring,parameters.clone(),order)?;
        let display=descriptor.to_string();
        let mut a = FractionPoly::from_expr(&self.as_expr(), &gens, &parameters, &display)?;
        let mut b = FractionPoly::from_expr(&o.as_expr(), &gens, &parameters, &display)?;
        a.coefficient_domain=descriptor.clone();b.coefficient_domain=descriptor;
        Ok(Some((a, b)))
    }
    fn pair(&self,o:&Self)->Option<(FractionPoly,FractionPoly)>{self.checked_fraction_pair(o).expect("compatible fraction coefficient domains")}

    pub fn add(&self, o: &Self) -> Self {
        if let Some((a,b))=self.algebraic_pair(o){return a.add(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.add(&b).into()
        } else {
            self.base.add(&o.base).into()
        }
    }
    pub fn sub(&self, o: &Self) -> Self {
        if let Some((a,b))=self.algebraic_pair(o){return a.sub(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.sub(&b).into()
        } else {
            self.base.sub(&o.base).into()
        }
    }
    pub fn mul(&self, o: &Self) -> Self {
        if let Some((a,b))=self.algebraic_pair(o){return a.mul(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.mul(&b).into()
        } else {
            self.base.mul(&o.base).into()
        }
    }
    pub fn neg(&self) -> Self {
        if let Some(p)=&self.algebraic {return p.neg().into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.neg().into(), |p| p.neg().into())
    }
    pub fn pow(&self, k: u32) -> Self {
        if let Some(p)=&self.algebraic {return p.pow(k).into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.pow(k).into(), |p| p.pow(k).into())
    }
    pub fn diff(&self, i: usize) -> Self {
        if let Some(p)=&self.algebraic {return p.diff(i).into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.diff(i).into(), |p| p.diff(i).into())
    }
    pub fn integrate(&self, i: usize) -> Self {
        if let Some(p)=&self.algebraic {return p.integrate(i).into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.integrate(i).into(), |p| p.integrate(i).into())
    }
    pub fn div(&self,o:&Self)->PResult<(Self,Self)>{self.div_with_options(o,true)}
    pub fn div_with_options(&self,o:&Self,auto:bool)->PResult<(Self,Self)>{
        if let Some((a,b))=self.checked_algebraic_pair(o)?{let(q,r)=a.div(&b)?;return Ok((q.into(),r.into()));}
        if let Some((a,b))=self.checked_fraction_pair(o)?{let(q,r)=a.div_with_options(&b,auto)?;return Ok((q.into(),r.into()));}
        if let Some(prime)=self.modulus.or(o.modulus) {
            if !self.base.domain_gens().is_empty()||!o.base.domain_gens().is_empty(){return Err(PolyError::Unsupported("finite-field division requires scalar coefficients; finite symbolic coefficient domains are not implemented".into()));}
            if self.modulus.is_some()&&o.modulus.is_some()&&self.modulus!=o.modulus{return Err(PolyError::Unsupported("cannot divide polynomials over different finite fields".into()));}
            let(a,b)=self.base.unify(&o.base);let _scope=super::f5b::PrimeScope::new(prime)?;
            let convert=|p:&MPoly<Rat>|->PResult<MPoly<super::f5b::Fp>>{Ok(MPoly::from_terms(p.nvars,p.iter().map(|(e,c)|Ok((e.to_vec(),super::f5b::Fp::from_rat(c,prime)?))).collect::<PResult<_>>()?))};
            let fa=convert(&a.rep)?;let fb=convert(&b.rep)?;
            if fb.is_zero(){return Err(PolyError::Unsupported("polynomial division by zero in the selected finite field".into()));}
            let(q,r)=recursive_division(&fa,&fb,0,&|_|true);
            return Ok((a.with_rep(q.map_coeffs(Coeff::to_rat)).into(),a.with_rep(r.map_coeffs(Coeff::to_rat)).into()));
        }
        if !auto&&self.modulus.is_none()&&o.modulus.is_none()&&self.domain.max(o.domain)==super::polyobj::Domain::ZZ {
            let(a,b)=self.base.unify(&o.base);if b.rep.is_zero(){return Err(PolyError::Unsupported("polynomial division by zero".into()));}
            let(q,r)=recursive_division(&a.rep,&b.rep,0,&|c|c.is_integer());
            return Ok((a.with_rep(q).into(),a.with_rep(r).into()));
        }
        let(q,r)=self.base.div(&o.base)?;Ok((q.into(),r.into()))
    }
    pub fn pdiv(&self,o:&Self)->PResult<(Self,Self)>{
        if let Some((a,b))=self.checked_algebraic_pair(o)?{
            let _scope=super::algebraic::AlgebraicScope::new(a.number_field.clone());
            let(q,r)=super::pseudo::pdiv(&a.rep,&b.rep,0)?;return Ok((a.with_rep(q).into(),a.with_rep(r).into()));
        }
        if let Some((a,b))=self.checked_fraction_pair(o)?{
            let _scope=Scope::new(a.parameters.len());let(q,r)=super::pseudo::pdiv(&a.rep,&b.rep,0)?;return Ok((a.with_rep(q).into(),a.with_rep(r).into()));
        }
        if let Some(prime)=self.modulus.or(o.modulus){
            if !self.base.domain_gens().is_empty()||!o.base.domain_gens().is_empty(){return Err(PolyError::Unsupported("finite-field pseudo-division requires scalar coefficients".into()));}
            if self.modulus.is_some()&&o.modulus.is_some()&&self.modulus!=o.modulus{return Err(PolyError::Unsupported("cannot pseudo-divide polynomials over different finite fields".into()));}
            let(a,b)=self.base.unify(&o.base);let _scope=super::f5b::PrimeScope::new(prime)?;
            let convert=|p:&MPoly<Rat>|->PResult<MPoly<super::f5b::Fp>>{Ok(MPoly::from_terms(p.nvars,p.iter().map(|(e,c)|Ok((e.to_vec(),super::f5b::Fp::from_rat(c,prime)?))).collect::<PResult<_>>()?))};
            let(q,r)=super::pseudo::pdiv(&convert(&a.rep)?,&convert(&b.rep)?,0)?;
            return Ok((a.with_rep(q.map_coeffs(Coeff::to_rat)).into(),a.with_rep(r.map_coeffs(Coeff::to_rat)).into()));
        }
        let(q,r)=self.base.pdiv(&o.base)?;Ok((q.into(),r.into()))
    }
    pub fn prem(&self,o:&Self)->PResult<Self>{Ok(self.pdiv(o)?.1)}
    pub fn pquo(&self,o:&Self)->PResult<Self>{Ok(self.pdiv(o)?.0)}
    pub fn pexquo(&self,o:&Self)->PResult<Self>{let(q,r)=self.pdiv(o)?;if !r.is_zero(){return Err(PolyError::ExactQuotientFailed("pseudo-division has a nonzero remainder".into()));}Ok(q)}
    pub fn rem(&self,o:&Self)->PResult<Self>{self.rem_with_options(o,true)}
    pub fn quo(&self,o:&Self)->PResult<Self>{self.quo_with_options(o,true)}
    fn retract_division_component(&self,o:&Self,mut result:Self)->Self{
        if let Some((a,_))=self.pair(o){
            if a.display.ends_with(']') {if let Some(p)=result.field.as_mut(){if p.rep.coeffs().iter().all(|c|coefficient_in_ring(c,a.display.starts_with("ZZ"))){p.display=a.display;p.coefficient_domain=a.coefficient_domain;}}}
        }else if self.algebraic.is_none()&&o.algebraic.is_none()&&self.modulus.is_none()&&o.modulus.is_none()&&self.domain.max(o.domain)==super::polyobj::Domain::ZZ&&result.base.rep.is_integral(){result.base.domain=super::polyobj::Domain::ZZ;}
        result
    }
    pub fn rem_with_options(&self,o:&Self,auto:bool)->PResult<Self>{let r=self.div_with_options(o,auto)?.1;Ok(self.retract_division_component(o,r))}
    pub fn quo_with_options(&self,o:&Self,auto:bool)->PResult<Self>{let q=self.div_with_options(o,auto)?.0;Ok(self.retract_division_component(o,q))}
    pub fn exquo(&self,o:&Self)->PResult<Self>{self.exquo_with_options(o,true)}
    pub fn exquo_with_options(&self,o:&Self,auto:bool)->PResult<Self>{
        let(q,r)=self.div_with_options(o,auto)?;
        if !r.is_zero(){return Err(PolyError::ExactQuotientFailed("polynomial division is not exact in the selected coefficient domain".into()));}Ok(q)
    }
    pub fn gcd(&self, o: &Self) -> Self {
        if let Some((a,b))=self.algebraic_pair(o){return a.gcd(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.gcd(&b).into()
        } else {
            self.base.gcd(&o.base).into()
        }
    }
    pub fn lcm(&self, o: &Self) -> Self {
        if self.algebraic.is_some()||o.algebraic.is_some(){if self.is_zero()||o.is_zero(){return self.mul(o);}return self.mul(o).exquo(&self.gcd(o)).unwrap().algebraic.unwrap().monic().into();}
        if let Some((a,b))=self.pair(o){
            a.lcm(&b).into()
        } else {
            self.base.lcm(&o.base).into()
        }
    }
    pub fn cofactors(&self, o: &Self) -> (Self, Self, Self) {
        let g = self.gcd(o);
        if g.is_zero() {
            return (g, self.clone(), o.clone());
        }
        let a = self.exquo_with_options(&g,false).expect("gcd divides first input in its ring");
        let b = o.exquo_with_options(&g,false).expect("gcd divides second input in its ring");
        (g, a, b)
    }
    pub fn to_sympy_string(&self) -> String {
        self.base.to_sympy_string_with_domain(&self.domain_str())
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rational_field_arithmetic_and_basis() {
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        let gens = vec![x.clone()];
        let parameters = vec![a];
        let f = FractionPoly::from_expr(
            &crate::parser::parse_str("x/a+1/(a+1)").unwrap(),
            &gens,
            &parameters,
            "QQ(a)",
        )
        .unwrap();
        let g = f.monic();
        let (q, r) = f.div(&g).unwrap();
        assert!(r.rep.is_zero());
        assert_eq!(q.mul(&g).rep, f.rep);
        let _s = Scope::new(1);
        let b = fieldalg::buchberger(&[f.rep], MonomialOrder::Lex);
        assert_eq!(b, vec![g.rep]);
    }

    #[test]
    fn scoped_coefficient_dimensions_survive_nesting_and_unwinding() {
        let previous = super::super::ratfunc::set_nvars(7);
        {
            let _one = Scope::new(1);
            assert_eq!(RF::one().n.nvars, 1);
            let x = Expr::symbol("x");
            let parameters = vec![Expr::symbol("a"), Expr::symbol("b")];
            let p = FractionPoly::from_expr(
                &crate::parser::parse_str("x/a+1/b").unwrap(),
                &[x],
                &parameters,
                "QQ(a,b)",
            )
            .unwrap();
            let _ = p.gcd(&p.pow(2));
            assert_eq!(RF::zero().n.nvars, 1);
            let caught = std::panic::catch_unwind(|| {
                let _nested = Scope::new(3);
                panic!("test unwinding");
            });
            assert!(caught.is_err());
            assert_eq!(RF::one().n.nvars, 1);
        }
        assert_eq!(RF::one().n.nvars, 7);
        super::super::ratfunc::set_nvars(previous);
        let workers: Vec<_> = (1..=2)
            .map(|n| {
                std::thread::spawn(move || {
                    let _scope = Scope::new(n);
                    for _ in 0..50 {
                        assert_eq!(RF::one().n.nvars, n);
                        std::thread::yield_now();
                    }
                })
            })
            .collect();
        for worker in workers {
            worker.join().unwrap();
        }
    }
}

#[cfg(test)]
mod implicit_ring_gcd_certificates {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn coefficient_first_display_preserves_main_generator_order(){
        let x=Expr::symbol("x");let y=Expr::symbol("y");
        for (source,generator,parameter,domain,expected) in [
            ("x*y+1",x.clone(),y.clone(),"ZZ[y]","Poly(y*x + 1, x, domain='ZZ[y]')"),
            ("x^2+2*x*y+1",y.clone(),x.clone(),"ZZ[x]","Poly(2*x*y + x**2 + 1, y, domain='ZZ[x]')"),
            ("x^2+x*y",y,x,"ZZ[x]","Poly(x*y + x**2, y, domain='ZZ[x]')"),
        ]{
            let p:Polynomial=FractionPoly::from_expr(&parse_str(source).unwrap(),&[generator],&[parameter],domain).unwrap().into();
            assert_eq!(p.to_sympy_string(),expected);
        }
    }
    #[test]
    fn ring_content_survives_fraction_field_primitive_gcd(){
        let x=Expr::symbol("x");let a=Expr::symbol("a");
        for (domain,left,right,expected) in [
            ("ZZ[a]","a*x","a*x^2","a*x"),
            ("ZZ[a]","2*a*x+2*a","4*a*x^2+8*a*x+4*a","2*a*x+2*a"),
            ("QQ[a]","2*a*x/3","4*a*x^2/5","2*a*x/15"),
            ("QQ[a]","2*(a+1)*x/3","4*(a+1)*x^2/5","(a+1)*x"),
            ("QQ[a]","x/2+1","x/3+2/3","x/6+1/3"),
            ("ZZ[a]","0","2*(a+1)*x","2*(a+1)*x"),
        ]{
            let p=|source:&str|FractionPoly::from_expr(&parse_str(source).unwrap(),std::slice::from_ref(&x),std::slice::from_ref(&a),domain).unwrap();
            let left=p(left);let right=p(right);let result=left.gcd(&right);
            assert_eq!(result.rep,p(expected).rep,"{domain}");
            assert!(left.div(&result).unwrap().1.rep.is_zero());assert!(right.div(&result).unwrap().1.rep.is_zero());
            assert!(result.belongs_to_ring());assert_eq!(result.coefficient_domain,left.coefficient_domain);
        }
    }
}

#[cfg(test)]
mod ring_division_certificates {
    use super::*;
    use crate::parser::parse_str;
    fn p(source:&str,domain:&str)->Polynomial{FractionPoly::from_expr(&parse_str(source).unwrap(),&[Expr::symbol("x")],&[Expr::symbol("a")],domain).unwrap().into()}
    #[test]
    fn ring_lcm_and_cofactors_preserve_content_and_membership(){
        for(domain,left,right,expected)in [
            ("ZZ[a]","a*x","a*x^2","a*x^2"),
            ("ZZ[a]","2*a*x","4*a*x^2","4*a*x^2"),
            ("ZZ[a]","2*(a+1)*x","4*(a+1)*x^2","4*(a+1)*x^2"),
            ("QQ[a]","2*a*x","4*a*x^2","a*x^2"),
            ("QQ[a]","2*(a+1)*x","4*(a+1)*x^2","(a+1)*x^2"),
            ("QQ[a]","x/2+1","x/3+2/3","x+2"),
            ("QQ[a]","2*x+1","4*x+2","2*x+1"),
            ("QQ[a]","2*(a+1)*(x+1)","4*(a+1)*(x+1)^2","4*(a+1)*(x+1)^2"),
        ]{
            let a=p(left,domain);let b=p(right,domain);let l=a.lcm(&b);let(g,u,v)=a.cofactors(&b);
            assert_eq!(l.field.as_ref().unwrap().rep,p(expected,domain).field.unwrap().rep,"{domain} {left} {right}");
            for c in [&l,&g,&u,&v]{assert_eq!(c.field.as_ref().unwrap().coefficient_domain,a.field.as_ref().unwrap().coefficient_domain);assert!(c.field.as_ref().unwrap().belongs_to_ring());}
            assert_eq!(g.mul(&u).field.unwrap().rep,a.field.as_ref().unwrap().rep);assert_eq!(g.mul(&v).field.unwrap().rep,b.field.as_ref().unwrap().rep);
            assert!(l.exquo_with_options(&a,false).is_ok());assert!(l.exquo_with_options(&b,false).is_ok());
        }
    }
    #[test]
    fn polynomial_division_is_recursive_and_auto_is_exact(){
        let a=p("x^2+2*x","ZZ[a]");let b=p("2*x","ZZ[a]");
        let(q,r)=a.div_with_options(&b,false).unwrap();assert!(q.is_zero());assert_eq!(r.as_expr(),a.as_expr());
        let a=p("x+1","ZZ[a]");let b=p("2*x+2","ZZ[a]");
        let(q,r)=a.div(&b).unwrap();assert_eq!(q.domain_str(),"ZZ(a)");assert_eq!(r.domain_str(),"ZZ(a)");
        assert_eq!(a.rem(&b).unwrap().domain_str(),"ZZ[a]");
        assert!(matches!(a.exquo_with_options(&b,false),Err(PolyError::ExactQuotientFailed(_))));
        assert_eq!(q.mul(&b).add(&r).as_expr(),a.as_expr());
        let a=p("a*x","ZZ[a]");let b=p("a*x","QQ[a]");
        assert_eq!(a.add(&b).domain_str(),"QQ[a]");
    }
}

#[cfg(test)]
mod typed_dense_division {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn finite_division_is_recursive_and_scopes_restore() {
        let gens=vec![Expr::symbol("x"),Expr::symbol("y")];
        let _outer=super::super::f5b::PrimeScope::new(7).unwrap();
        for prime in [2,5,7,101] {
            for (f,g,expected) in [("x*y+x","y","0"),("(x+y)*(x+2)+1","x+y","x+2")] {
                let p=|s:&str|->Polynomial{Poly::from_expr(&parse_str(s).unwrap(),&gens).unwrap().set_modulus(prime).into()};
                let a=p(f);let b=p(g);let(q,r)=a.div(&b).unwrap();
                assert_eq!(q.base.rep,p(expected).base.rep);assert_eq!(q.mul(&b).add(&r).base.rep,a.base.rep);
                assert_eq!(q.modulus,Some(prime));assert_eq!(r.modulus,Some(prime));
            }
        }
        assert!(super::super::f5b::Fp::from_i64(14).is_zero());
        let a:Polynomial=Poly::from_expr(&parse_str("x+1").unwrap(),&gens).unwrap().set_modulus(5).into();
        let b:Polynomial=Poly::from_expr(&parse_str("x+1").unwrap(),&gens).unwrap().set_modulus(7).into();
        assert!(a.div(&b).is_err());
        let opaque:Polynomial=Poly::from_expr(&parse_str("sqrt(2)*x+1").unwrap(),&[Expr::symbol("x")]).unwrap().set_modulus(5).into();
        assert!(a.div(&opaque).is_err());
    }
    #[test]
    fn algebraic_division_stops_at_coefficient_remainder() {
        let gens=vec![Expr::symbol("x"),Expr::symbol("y")];
        for generator in ["sqrt(2)","2^(1/3)+I"] {
            let field=super::super::algebraic::NumberField::from_generator(&parse_str(generator).unwrap()).unwrap();
            let p=|s:&str|->Polynomial{super::super::algebraic::AlgebraicPoly::from_expr(&parse_str(s).unwrap(),&gens,field.clone()).unwrap().into()};
            for(f,g,expected)in[("x*y+x","y","0"),("(x+sqrt(2))*(y+1)+1","x+sqrt(2)","y+1")] {
                if generator!="sqrt(2)"&&f.contains("sqrt"){continue;}
                let a=p(f);let b=p(g);let(q,r)=a.div(&b).unwrap();
                assert_eq!(q.algebraic.as_ref().unwrap().rep,p(expected).algebraic.unwrap().rep);
                assert_eq!(q.mul(&b).add(&r).algebraic.unwrap().rep,a.algebraic.unwrap().rep);
            }
        }
    }
}

#[cfg(test)]
mod functional_decomposition {
    use super::*;
    use crate::parser::parse_str;

    #[test]
    fn rational_function_parts_reconstruct_and_respect_ring_membership() {
        let x=Expr::symbol("x");let t=Expr::symbol("t");
        for (domain,source,count) in [
            ("QQ(t)","(x^2+x/t)^3+x^2+x/t",2),
            ("QQ[t]","(t*x^2+x)^2",1),
            ("ZZ[t]","4*x^4+4*x^3+x^2",1),
            ("QQ[t]","4*x^4+4*x^3+x^2",2),
            ("ZZ[t]","(x^2+t*x)^3+1",2),
        ] {
            let polynomial:Polynomial=FractionPoly::from_expr(&parse_str(source).unwrap(),std::slice::from_ref(&x),std::slice::from_ref(&t),domain).unwrap().into();
            let parts=polynomial.decompose().unwrap();assert_eq!(parts.len(),count,"{domain}");
            let mut rebuilt=parts[0].field.as_ref().unwrap().clone();
            for inner in &parts[1..] {rebuilt=rebuilt.compose(inner.field.as_ref().unwrap());}
            assert_eq!(rebuilt.rep,polynomial.field.as_ref().unwrap().rep);
            assert!(parts.iter().all(|part|part.domain_str()==domain));
            if domain.ends_with(']') {assert!(parts.iter().all(|part|part.field.as_ref().unwrap().belongs_to_ring()));}
        }
    }

    #[test]
    fn algebraic_parts_use_the_exact_number_field() {
        let x=Expr::symbol("x");
        for generator in ["sqrt(2)","I","sqrt(2)+sqrt(3)"] {
            let field=super::super::algebraic::NumberField::from_generator(&parse_str(generator).unwrap()).unwrap();
            let source=format!("(x^2+({generator})*x)^3+2*(x^2+({generator})*x)+1");
            let polynomial:Polynomial=super::super::algebraic::AlgebraicPoly::from_expr(&parse_str(&source).unwrap(),&[x.clone()],field).unwrap().into();
            let parts=polynomial.decompose().unwrap();assert_eq!(parts.len(),2);
            let rebuilt=parts[0].algebraic.as_ref().unwrap().compose(parts[1].algebraic.as_ref().unwrap());
            assert_eq!(rebuilt.rep,polynomial.algebraic.unwrap().rep);
        }
    }

    #[test]
    fn prime_and_integer_decomposition_keep_their_arithmetic() {
        let x=Expr::symbol("x");
        let polynomial=Poly::from_expr(&parse_str("4*x^4+4*x^3+x^2").unwrap(),&[x.clone()]).unwrap();
        assert_eq!(polynomial.decompose().unwrap().len(),1);
        let _scope=super::super::f5b::PrimeScope::new(7).unwrap();
        for prime in [3,5,11] {
            let polynomial=Poly::from_expr(&parse_str("x^4+2*x^2+1").unwrap(),&[x.clone()]).unwrap().set_modulus(prime);
            let parts=polynomial.decompose().unwrap();assert_eq!(parts.len(),2);
            assert_eq!(parts[0].compose(&parts[1]).rep,polynomial.rep);
        }
        assert!(super::super::f5b::Fp::from_i64(14).is_zero());
        let polynomial=Poly::from_expr(&parse_str("x^4+1").unwrap(),&[x]).unwrap().set_modulus(2);
        assert!(matches!(polynomial.decompose(),Err(PolyError::NotInvertible(_))));
    }
}
