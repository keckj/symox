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

#[derive(Clone, Debug)]
pub struct FractionPoly {
    pub rep: MPoly<RF>,
    pub gens: Vec<Expr>,
    pub parameters: Vec<Expr>,
    pub display: String,
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

fn coefficient(value: &Expr, parameters: &[Expr]) -> PResult<RF> {
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
        Self::from_expr(e, &self.gens, &self.parameters, &self.display)
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
            p.display = p.display.replace('[', "(").replace(']', ")");
        }
        p
    }
    pub fn div(&self, other: &Self) -> PResult<(Self, Self)> {
        if other.rep.is_zero() {
            return Err(PolyError::Unsupported("polynomial division by zero".into()));
        }
        let _s = Scope::new(self.parameters.len());
        let (q, r) = fieldalg::reduce(
            &self.rep,
            std::slice::from_ref(&other.rep),
            MonomialOrder::Lex,
        );
        let mut q = self.with_rep(q.into_iter().next().unwrap());
        let mut r = self.with_rep(r);
        if self.display.ends_with(']') && (!q.belongs_to_ring() || !r.belongs_to_ring()) {
            q.display = q.display.replace('[', "(").replace(']', ")");
            r.display = q.display.clone();
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
        self.with_rep(MPoly::from_terms(self.gens.len(), terms))
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
            let mut p = p.ring_normalize();
            p.display = self.display.clone();
            p
        } else {
            p
        }
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
        if self.field.is_none() && self.algebraic.is_none() {
            return self.base.all_coeffs();
        }
        if self.ngens != 1 {
            return Err(PolyError::Unsupported(
                "univariate polynomial expected".into(),
            ));
        }
        let n = self.base.rep.degree(0).max(0) as u32;
        Ok((0..=n).rev().map(|i| self.coeff_monomial(&[i])).collect())
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
            if a.number_field==b.number_field {a.number_field.clone()}else{
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
    fn pair(&self, o: &Self) -> Option<(FractionPoly, FractionPoly)> {
        let prototype = self.field.as_ref().or(o.field.as_ref())?;
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
        let mut display = domain_name(&parameters);
        if prototype.display.starts_with("ZZ(")
            && self.field.as_ref().is_none_or(|f|!f.display.starts_with("QQ"))
            && o.field.as_ref().is_none_or(|f|!f.display.starts_with("QQ")) {
            display=display.replacen("QQ","ZZ",1);
        }
        if prototype.display.ends_with(']')
            && self.field.as_ref().is_none_or(|f| f.display.ends_with(']'))
            && o.field.as_ref().is_none_or(|f| f.display.ends_with(']'))
        {
            display = display.replace('(', "[").replace(')', "]");
            if prototype.display.starts_with("ZZ")
                && self.domain == super::polyobj::Domain::ZZ
                && o.domain == super::polyobj::Domain::ZZ
            {
                display = display.replacen("QQ", "ZZ", 1);
            }
        }
        let a = FractionPoly::from_expr(&self.as_expr(), &gens, &parameters, &display)
            .expect("compatible fraction fields");
        let b = FractionPoly::from_expr(&o.as_expr(), &gens, &parameters, &display)
            .expect("compatible fraction fields");
        Some((a, b))
    }
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
    pub fn div(&self, o: &Self) -> PResult<(Self, Self)> {
        if let Some((a,b))=self.checked_algebraic_pair(o)?{let(q,r)=a.div(&b)?;return Ok((q.into(),r.into()));}
        if let Some((a, b)) = self.pair(o) {
            let (q, r) = a.div(&b)?;
            Ok((q.into(), r.into()))
        } else {
            let (q, r) = self.base.div(&o.base)?;
            Ok((q.into(), r.into()))
        }
    }
    pub fn rem(&self, o: &Self) -> PResult<Self> {
        if self.field.is_none() && o.field.is_none() && self.algebraic.is_none() && o.algebraic.is_none() {
            return Ok(self.base.rem(&o.base)?.into());
        }
        Ok(self.div(o)?.1)
    }
    pub fn quo(&self, o: &Self) -> PResult<Self> {
        if self.field.is_none() && o.field.is_none() && self.algebraic.is_none() && o.algebraic.is_none() {
            return Ok(self.base.quo(&o.base)?.into());
        }
        Ok(self.div(o)?.0)
    }
    pub fn exquo(&self, o: &Self) -> PResult<Self> {
        let (q, r) = self.div(o)?;
        if r.is_zero() {
            Ok(q)
        } else {
            Err(PolyError::Failed("polynomial division is not exact".into()))
        }
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
        if self.field.is_some() || o.field.is_some() {
            if self.is_zero() || o.is_zero() {
                return self.mul(o);
            }
            self.mul(o)
                .exquo(&self.gcd(o))
                .unwrap()
                .field
                .unwrap()
                .monic()
                .into()
        } else {
            self.base.lcm(&o.base).into()
        }
    }
    pub fn cofactors(&self, o: &Self) -> (Self, Self, Self) {
        let g = self.gcd(o);
        if g.is_zero() {
            return (g, self.clone(), o.clone());
        }
        let a = self.exquo(&g).unwrap();
        let b = o.exquo(&g).unwrap();
        (g, a, b)
    }
    pub fn to_sympy_string(&self) -> String {
        if self.field.is_none() && self.algebraic.is_none() {
            return self.base.to_sympy_string();
        }
        format!(
            "Poly({}, {}, domain='{}')",
            self.as_expr(),
            self.true_gens()
                .iter()
                .map(crate::printer::to_str)
                .collect::<Vec<_>>()
                .join(", "),
            self.domain_str()
        )
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
