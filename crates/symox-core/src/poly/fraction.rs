//! Exact polynomials over rational function coefficient fields.

use super::convert::{from_exprs_gens, sort_gens, to_expr, to_expr_int, PResult, PolyError};
use super::fieldalg;
use super::expression_domain::ExpressionPoly;
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

pub fn parse_domain(domain: &str, available: &[Expr]) -> PResult<Option<Vec<Expr>>> {super::mpoly::checked_operation(||{
    if domain.starts_with("QQ<") {return Ok(None);}
    let Some(open) = domain.find(['(', '[']) else {
        return Ok(None);
    };
    let closing = if domain.as_bytes()[open] == b'(' {
        ')'
    } else {
        ']'
    };
    // The reference fraction-field pattern only admits ZZ and QQ grounds.
    if closing == ')' && matches!(&domain[..open], "ZZ_I" | "QQ_I" | "R" | "RR" | "C" | "CC") {
        return Err(PolyError::Reference {
            kind: "OptionError",
            message: format!("expected a valid domain specification, got {domain}"),
        });
    }
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
        let parsed = crate::parser::parse_str(name)
            .map_err(|_| PolyError::Failed("invalid fraction field parameter".into()))?;
        // As the reference sympifies each name, `ZZ[pi]`, `ZZ[2**pi]` or
        // `QQ(log(x))` name generator expressions; a plain name keeps the
        // assumptions of a matching symbol of the expression.
        let value = if matches!(parsed.kind(), Kind::Sym(_)) {
            candidates
                .first()
                .map(|s| (*s).clone())
                .unwrap_or_else(|| Expr::symbol(name))
        } else {
            parsed
        };
        if value.as_number().is_some() || parameters.contains(&value) {
            return Err(PolyError::Failed(
                "fraction field parameters must be distinct nonnumeric expressions".into(),
            ));
        }
        parameters.push(value);
    }
    Ok(Some(parameters))
})}

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

pub fn coefficient(value: &Expr, parameters: &[Expr]) -> PResult<RF> {super::mpoly::checked_operation(||{
    match coefficient_direct(value, parameters) {
        Ok(c) => Ok(c),
        Err(first) => {
            // The reference `FracField._rebuild_expr` reads a power
            // `b**e` as `g**(e/eg)` for a generator `g = bg**eg` with the
            // same base: over `QQ(exp(1/2))`, `E` is `exp(1/2)**2` and
            // `exp(-3/2)` is `exp(1/2)**-3`.
            match rebuild_powers(value, parameters) {
                Some((rebuilt, dummies)) => coefficient_direct(&rebuilt, &dummies).map_err(|_| first),
                None => Err(first),
            }
        }
    }
})}

/// `value` with the parameters and the integer powers of the power-like
/// parameters (`Pow` or `exp`) replaced by fresh symbols, or `None` when no
/// parameter is a power.
fn rebuild_powers(value: &Expr, parameters: &[Expr]) -> Option<(Expr, Vec<Expr>)> {
    use crate::decompose::as_base_exp;
    let powered = |g: &Expr| matches!(g.kind(), Kind::Pow(..) | Kind::Func(crate::expr::Func::Exp, _));
    if !parameters.iter().any(powered) {
        return None;
    }
    let dummies: Vec<Expr> = parameters.iter().map(|_| Expr::dummy("g", crate::symbol::Assumptions::NONE)).collect();
    let powers: Vec<(usize, Expr, Expr)> = parameters
        .iter()
        .enumerate()
        .filter(|(_, g)| powered(g))
        .map(|(i, g)| {
            let (b, e) = as_base_exp(g);
            (i, b, e)
        })
        .collect();
    fn go(e: &Expr, parameters: &[Expr], dummies: &[Expr], powers: &[(usize, Expr, Expr)]) -> Expr {
        if let Some(i) = parameters.iter().position(|p| p == e) {
            return dummies[i].clone();
        }
        match e.kind() {
            Kind::Num(_) | Kind::Sym(_) | Kind::Bool(_) => e.clone(),
            Kind::Add(..) => add(&e.args().iter().map(|a| go(a, parameters, dummies, powers)).collect::<Vec<_>>()),
            Kind::Mul(c, facs) => {
                let mut fs = vec![Expr::num(c.clone())];
                fs.extend(facs.iter().map(|(b, x)| go(&crate::ops::mul_factor(b, x), parameters, dummies, powers)));
                mul(&fs)
            }
            Kind::Pow(..) | Kind::Func(crate::expr::Func::Exp, _) | Kind::Const(crate::expr::Constant::E) => {
                let (b, x) = as_base_exp(e);
                for (i, bg, eg) in powers {
                    if *bg == b {
                        let q = div(&x, eg);
                        if let Some(n) = q.as_i64() {
                            return pow(&dummies[*i], &Expr::int(n));
                        }
                    }
                }
                if x.as_i64().is_some_and(|n| n != 1) {
                    return pow(&go(&b, parameters, dummies, powers), &x);
                }
                e.clone()
            }
            _ => e.clone(),
        }
    }
    Some((go(value, parameters, &dummies, &powers), dummies))
}

fn coefficient_direct(value: &Expr, parameters: &[Expr]) -> PResult<RF> {
    let together = crate::simplify::together(value);
    let (num, den) = crate::simplify::numer_denom(&together);
    let first = from_exprs_gens(&[num, den], parameters);
    let complete = |c: &super::convert::Converted| !c.float && c.gens.len() == parameters.len();
    let converted = match first {
        Ok(c) if c.float || complete(&c) => c,
        first => {
            // A generator such as `exp(a)` reads `exp(-a)` as `1/exp(a)`
            // only through the reference `as_numer_denom`; without it
            // `exp(-s*log(a))` is an unknown factor containing the
            // parameter `log(a)`.
            let (num, den) = super::api::as_numer_denom(value);
            match from_exprs_gens(&[num, den], parameters) {
                Ok(retry) if complete(&retry) => retry,
                _ => first?,
            }
        }
    };
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
    pub fn from_expr(e: &Expr, gens: &[Expr], parameters: &[Expr], display: &str) -> PResult<Self> {super::mpoly::checked_operation(||{
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
    })}
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
        Poly::from_coefficient_terms(self.terms(),&self.gens)
            .expect("field coefficients have a polynomial metadata view")
    }
    pub fn convert(&self, e: &Expr) -> PResult<Self> {super::mpoly::checked_operation(||{
        Self::from_expr(e, &self.gens, &self.parameters, &self.display).map(|mut p| {p.coefficient_domain=self.coefficient_domain.clone();p})
    })}
    pub fn promote_field(&mut self) {
        self.coefficient_domain=self.coefficient_domain.field();
        // A nested ring `ZZ[x][y]` has the field `ZZ[x](y)`.
        self.display=if matches!(self.coefficient_domain,super::domains::ExactDomain::Composite{..}) {self.coefficient_domain.to_string()} else {self.display.replace('[',"(").replace(']',")")};
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
    pub fn div(&self,other:&Self)->PResult<(Self,Self)>{super::mpoly::checked_operation(||{self.div_with_options(other,true)})}
    pub fn div_with_options(&self, other: &Self, auto:bool) -> PResult<(Self, Self)> {super::mpoly::checked_operation(||{
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
    })}
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
    /// The coefficient-ring `gcd` (`PolynomialRing.gcd`) of two elements.
    fn ring_gcd(&self, a: &RF, b: &RF) -> RF {
        use super::domains::{DomainElement, ElementValue};
        let wrap = |c: &RF| DomainElement { domain: self.coefficient_domain.clone(), value: ElementValue::Function(c.clone()) };
        match wrap(a).gcd(&wrap(b)).expect("valid polynomial-ring coefficient gcd").value {
            ElementValue::Function(c) => c,
            _ => unreachable!(),
        }
    }
    /// `K.is_negative` of a coefficient-ring element: the sign of its
    /// leading ground coefficient.
    fn ring_negative(c: &RF) -> bool {
        !c.n.is_zero() && (c.n.lc().is_negative() != c.d.lc().is_negative())
    }
    /// The reference `dup_rr_prs_gcd` over a polynomial coefficient ring:
    /// the gcd of the contents times the primitive part of the last
    /// subresultant, normalized by `canonical_unit` of its leading term.
    fn reference_ring_gcd(&self, other: &Self) -> PResult<Self> {
        let _s = Scope::new(self.parameters.len());
        let nonnegative = |p: &Self| if Self::ring_negative(&p.rep.lc()) { p.neg() } else { p.clone() };
        match (self.rep.is_zero(), other.rep.is_zero()) {
            (true, true) => return Ok(self.clone()),
            (true, false) => return Ok(self.with_rep(nonnegative(other).rep)),
            (false, true) => return Ok(nonnegative(self)),
            _ => {}
        }
        let (fc, f) = self.reference_primitive();
        let (gc, g) = other.reference_primitive();
        let mut c = self.ring_gcd(&fc, &gc);
        let prs = super::prs::subresultants(&f.rep, &g.rep, 0)?;
        let last = prs.sequence.last().cloned().unwrap_or_else(|| g.rep.clone());
        let (_, h) = self.with_rep(last).reference_primitive();
        if Self::ring_negative(&h.rep.lc()) {
            c = c.neg();
        }
        let _s = Scope::new(self.parameters.len());
        Ok(self.with_rep(h.rep.scale(&c)))
    }
    pub fn gcd(&self, other: &Self) -> Self {
        if self.display.ends_with(']') && self.gens.len() == 1 {
            if let Ok(h) = self.reference_ring_gcd(other) {
                return h;
            }
        }
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
        if self.gens.len()==1 {
            // The reference `dup_rr_lcm`: lcm of the contents times the
            // primitive product over the gcd, with its canonical unit.
            let(fc,f)=self.reference_primitive();let(gc,g)=other.reference_primitive();
            use super::domains::{DomainElement,ElementValue};
            let wrap=|c:&RF|DomainElement{domain:self.coefficient_domain.clone(),value:ElementValue::Function(c.clone())};
            let ElementValue::Function(mut c)=wrap(&fc).lcm(&wrap(&gc)).expect("valid coefficient-ring lcm").value else{unreachable!()};
            let product=f.mul(&g);let gcd=f.gcd(&g);
            let h=product.div_with_options(&gcd,false).expect("primitive gcd divides the product in its ring").0;
            if Self::ring_negative(&h.rep.lc()){c=c.neg();}
            let _scope=Scope::new(self.parameters.len());
            return self.with_rep(h.rep.scale(&c));
        }
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
    /// The polynomial of a cleared numerator in `gens + parameters`.
    fn from_cleared(&self, p: &MPoly<Int>) -> Self {
        let n = self.gens.len();
        let np = self.parameters.len();
        let mut groups: Vec<(Vec<u32>, Vec<(Vec<u32>, Int)>)> = Vec::new();
        for (e, c) in p.iter() {
            let (key, rest) = e.split_at(n);
            match groups.iter_mut().find(|(k, _)| k.as_slice() == key) {
                Some(slot) => slot.1.push((rest.to_vec(), c.clone())),
                None => groups.push((key.to_vec(), vec![(rest.to_vec(), c.clone())])),
            }
        }
        let terms = groups
            .into_iter()
            .map(|(k, t)| (k, RF::from_poly(MPoly::from_terms(np, t))))
            .collect();
        self.with_rep(MPoly::from_terms(n, terms))
    }
    fn lc_is_negative(&self) -> bool {
        !self.rep.is_zero() && self.rep.lc().is_negative()
    }
    /// The reference `dup_half_gcdex` over the fraction field, for one
    /// generator: `s*f = h (mod g)` with `h` monic.
    pub fn half_gcdex(&self, g: &Self) -> PResult<(Self, Self)> {super::mpoly::checked_operation(||{
        if self.gens.len() != 1 {
            return Err(PolyError::Unsupported("univariate polynomial expected".into()));
        }
        let _s = Scope::new(self.parameters.len());
        let (mut f, mut g) = (self.clone(), g.clone());
        let mut a = self.with_rep(MPoly::constant(RF::one(), 1));
        let mut b = self.with_rep(MPoly::zero(1));
        while !g.rep.is_zero() {
            let (q, r) = f.div(&g)?;
            f = g;
            g = r;
            let next = a.sub(&q.mul(&b));
            a = b;
            b = next;
        }
        if f.rep.is_zero() {
            return Err(PolyError::Unsupported("polynomial division by zero".into()));
        }
        let lc = f.rep.lc().clone();
        let a = a.with_rep(a.rep.scale(&lc.inv()));
        Ok((a, f.monic()))
    })}
    /// The reference `dup_gcdex` over the fraction field.
    pub fn gcdex(&self, g: &Self) -> PResult<(Self, Self, Self)> {super::mpoly::checked_operation(||{
        let (s, h) = self.half_gcdex(g)?;
        let (t, r) = h.sub(&s.mul(self)).div(g)?;
        if !r.rep.is_zero() {
            return Err(PolyError::ExactQuotientFailed("extended Euclid cofactor is not exact".into()));
        }
        Ok((s, t, h))
    })}
    /// The reference `dup_invert`: the inverse of `self` modulo `g`.
    pub fn invert(&self, g: &Self) -> PResult<Self> {super::mpoly::checked_operation(||{
        let (s, h) = self.half_gcdex(g)?;
        if h.rep.is_constant() && h.rep.constant_value().is_some_and(|c| c.is_one()) {
            return Ok(s.div(g)?.1);
        }
        Err(PolyError::NotInvertible("zero divisor".into()))
    })}
    /// The reference `dmp_cancel` over a fraction field with an associated
    /// ring: clear both denominators, cancel the gcd over the ring, and
    /// return `(cp, cq, p, q)` with `self/g = (cp*p)/(cq*q)`.
    pub fn reference_cancel(&self, g: &Self) -> (RF, RF, Self, Self) {
        let _s = Scope::new(self.parameters.len());
        let cq = self.denominator();
        let f = self.clear_coeff_denominators();
        let cp = g.denominator();
        let h = g.clear_coeff_denominators();
        let (_, p, q) = super::gcd::gcd_mz(&f, &h);
        let (_, cp, cq) = super::gcd::gcd_mz(&cp, &cq);
        let mut p = self.from_cleared(&p);
        let mut q = self.from_cleared(&q);
        let mut cp = RF::from_poly(cp);
        let cq = RF::from_poly(cq);
        match (p.lc_is_negative(), q.lc_is_negative()) {
            (true, true) => {
                p = p.neg();
                q = q.neg();
            }
            (true, false) => {
                cp = cp.neg();
                p = p.neg();
            }
            (false, true) => {
                cp = cp.neg();
                q = q.neg();
            }
            _ => {}
        }
        (cp, cq, p, q)
    }
    /// The reference `dmp_ground_content`: the coefficient-domain `gcd`
    /// (over a fraction field, numerator gcd over denominator lcm; over a
    /// polynomial ring, `PolyElement.gcd`) folded over the coefficients, recursively by generator and descending degree, which
    /// stops as soon as the running content is one.
    pub fn reference_content(&self) -> RF {
        let _s = Scope::new(self.parameters.len());
        let np = self.parameters.len();
        fn normalized(p: MPoly<Int>) -> MPoly<Int> {
            if !p.is_zero() && p.lc().is_negative() { p.neg() } else { p }
        }
        fn field_gcd(a: &RF, b: &RF, np: usize) -> RF {
            let p = normalized(super::gcd::gcd_mz(&a.n, &b.n).0);
            if p.is_zero() {
                return RF::new(p, MPoly::one(np));
            }
            let g = super::gcd::gcd_mz(&a.d, &b.d).0;
            let q = normalized(a.d.div_exact(&g).unwrap().mul(&b.d));
            RF::new(p, q)
        }
        let ring = self.display.ends_with(']');
        let gcd = |a: &RF, b: &RF| if ring { self.ring_gcd(a, b) } else { field_gcd(a, b, np) };
        fn content(terms: &[(Vec<u32>, RF)], var: usize, nvars: usize, gcd: &dyn Fn(&RF, &RF) -> RF) -> RF {
            let mut cont = RF::zero();
            if var + 1 >= nvars {
                let mut sorted: Vec<&(Vec<u32>, RF)> = terms.iter().collect();
                sorted.sort_by(|a, b| b.0.get(var).cmp(&a.0.get(var)));
                for (_, c) in sorted {
                    cont = gcd(&cont, c);
                    if cont.is_one() {
                        break;
                    }
                }
                return cont;
            }
            let mut degrees: Vec<u32> = terms.iter().map(|t| t.0[var]).collect();
            degrees.sort_unstable_by(|a, b| b.cmp(a));
            degrees.dedup();
            for d in degrees {
                let group: Vec<(Vec<u32>, RF)> = terms.iter().filter(|t| t.0[var] == d).cloned().collect();
                cont = gcd(&cont, &content(&group, var + 1, nvars, gcd));
                if cont.is_one() {
                    break;
                }
            }
            cont
        }
        let terms: Vec<(Vec<u32>, RF)> = self.rep.iter().map(|(e, c)| (e.to_vec(), c.clone())).collect();
        if terms.is_empty() {
            return RF::zero();
        }
        content(&terms, 0, self.gens.len(), &gcd)
    }
    /// The reference `Poly.primitive` over the fraction field, with the
    /// content of `reference_content`.
    pub fn reference_primitive(&self) -> (RF, Self) {
        let c = self.reference_content();
        if c.is_zero() || c.is_one() {
            return (c, self.clone());
        }
        let _s = Scope::new(self.parameters.len());
        (c.clone(), self.with_rep(self.rep.scale(&c.inv())))
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
    /// The reference expression domain `EX`; `base` is then only a
    /// metadata view and is never used for arithmetic.
    pub expression: Option<super::expression_domain::ExpressionPoly>,
}
impl From<Poly> for Polynomial {
    fn from(base: Poly) -> Self {
        Self { base, field: None, algebraic:None, expression:None }
    }
}
impl From<FractionPoly> for Polynomial {
    fn from(field: FractionPoly) -> Self {
        Self {
            base: field.legacy(),
            field: Some(field),
            algebraic:None,
            expression:None,
        }
    }
}
impl From<super::algebraic::AlgebraicPoly> for Polynomial {
    fn from(p:super::algebraic::AlgebraicPoly)->Self {Self {base:p.legacy(),field:None,algebraic:Some(p),expression:None}}
}
impl From<super::expression_domain::ExpressionPoly> for Polynomial {
    fn from(p:super::expression_domain::ExpressionPoly)->Self {
        let p=super::gaussian_ground::normalize_poly(p);
        Self {base:p.legacy(),field:None,algebraic:None,expression:Some(p)}
    }
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
    pub fn decompose(&self)->PResult<Vec<Self>> {super::mpoly::checked_operation(||{
        if let Some(p)=&self.expression {return Ok(p.decompose()?.into_iter().map(Into::into).collect());}
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
    })}

    pub fn as_expr(&self) -> Expr {
        if let Some(p)=&self.expression {return p.as_expr();}
        if let Some(p)=&self.algebraic {return p.as_expr();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.as_expr(), FractionPoly::as_expr)
    }
    pub fn domain_str(&self) -> String {
        if let Some(p)=&self.expression {return p.ground.as_ref().map_or_else(||"EX".into(),|g|g.display());}
        if let Some(p)=&self.algebraic {
            return match p.gaussian {
                Some(true) => "ZZ_I".into(),
                Some(false) => "QQ_I".into(),
                None => p.number_field.display(),
            };
        }
        self.field
            .as_ref()
            .map_or_else(|| self.base.domain_str(), |p| p.display.clone())
    }
    pub fn terms(&self) -> Vec<(Vec<u32>, Expr)> {
        if let Some(p)=&self.expression {return p.terms();}
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
        if let Some(p)=&self.expression{return p.tc();}
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
        if let Some(p)=&self.expression{return p.coefficient(e);}
        self.terms()
            .into_iter()
            .find(|(m, _)| m == e)
            .map(|(_, c)| c)
            .unwrap_or_else(Expr::zero)
    }
    pub fn all_coeffs(&self) -> PResult<Vec<Expr>> {super::mpoly::checked_operation(||{
        let terms=super::structure::dense_terms(self)?;
        if let Some(p)=self.expression.as_ref().filter(|p|!p.is_sparse()){
            return Ok(terms.into_iter().map(|(m,_)|p.coefficient(&m)).collect());
        }
        Ok(terms.into_iter().map(|(_,c)|c).collect())
    })}
    /// Checked coercion boundary for public fallible operations. Coefficients
    /// from incompatible fields never enter arithmetic or an expect here.
    pub fn checked_algebraic_pair(&self,o:&Self)->PResult<Option<(super::algebraic::AlgebraicPoly,super::algebraic::AlgebraicPoly)>> {super::mpoly::checked_operation(||{
        if self.expression.is_some()||o.expression.is_some(){return Ok(None);}
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
        // The reference unification of ZZ_I/QQ_I with each other, with ZZ
        // or QQ, and with a genuine algebraic field.
        let gaussian=match (self.algebraic.as_ref().map(|a|a.gaussian),o.algebraic.as_ref().map(|a|a.gaussian)) {
            (Some(Some(a)),Some(Some(b)))=>Some(a&&b),
            (Some(Some(a)),None)=>Some(a&&o.domain==super::polyobj::Domain::ZZ),
            (None,Some(Some(b)))=>Some(b&&self.domain==super::polyobj::Domain::ZZ),
            _=>None,
        };
        let convert=|p:&Self| {
            if let Some(a)=&p.algebraic {if a.number_field==field && a.gens==gens {let mut a=a.clone();a.gaussian=gaussian;return Ok(a);}}
            let _scope=super::algebraic::AlgebraicScope::new(field.clone());
            let positions:Vec<_>=p.true_gens().iter().map(|g|gens.iter().position(|candidate|candidate==g).unwrap()).collect();
            let terms=p.terms().into_iter().map(|(powers,coefficient)|{
                let mut remapped=vec![0;gens.len()];
                for (index,exponent) in positions.iter().zip(powers){remapped[*index]=exponent;}
                Ok((remapped,super::algebraic::AN::from_expr(&coefficient,field.clone())?))
            }).collect::<PResult<Vec<_>>>()?;
            Ok(super::algebraic::AlgebraicPoly{rep:MPoly::from_terms(gens.len(),terms),gens:gens.clone(),number_field:field.clone(),gaussian})
        };
        Ok(Some((convert(self)?,convert(o)?)))
    })}
    fn algebraic_pair(&self,o:&Self)->Option<(super::algebraic::AlgebraicPoly,super::algebraic::AlgebraicPoly)> {
        self.checked_algebraic_pair(o).expect("compatible algebraic fields")
    }
    pub fn checked_fraction_pair(&self, o: &Self) -> PResult<Option<(FractionPoly, FractionPoly)>> {super::mpoly::checked_operation(||{
        if self.expression.is_some()||o.expression.is_some(){return Ok(None);}
        let Some(prototype)=self.field.as_ref().or(o.field.as_ref()) else{return Ok(None)};
        if self.modulus.is_some()||o.modulus.is_some()||self.algebraic.is_some()||o.algebraic.is_some(){return Err(PolyError::Unsupported("cannot combine rational-function, algebraic, and finite coefficient fields".into()));}
        if let(Some(a),Some(b))=(&self.field,&o.field){if let(super::domains::ExactDomain::Functions{order:a,..},super::domains::ExactDomain::Functions{order:b,..})=(&a.coefficient_domain,&b.coefficient_domain){if a!=b{return Err(PolyError::Unsupported("cannot unify coefficient rings with different monomial orders".into()));}}}
        if let(Some(a),Some(b))=(&self.field,&o.field){
            if a.coefficient_domain==b.coefficient_domain&&a.gens==b.gens{return Ok(Some((a.clone(),b.clone())));}
        }
        let mut gens = super::polyobj::unify_gens(self.true_gens(), o.true_gens());
        if gens.is_empty() {
            gens = prototype.gens.clone();
        }
        if let Some(pair)=self.composite_pair(o,&gens)? {return Ok(Some(pair));}
        // The reference `_unify_gens` of both coefficient domains, the
        // receiver first, never sorted.
        let side = |p: &Self| -> Vec<Expr> {
            match &p.field {
                Some(f) => f.parameters.clone(),
                None => p.base.domain_gens().to_vec(),
            }
        };
        let mut parameters = super::polyobj::unify_gens(&side(self), &side(o));
        let coefficients:Vec<_>=[self,o].iter().flat_map(|p|p.terms().into_iter().map(|(_,c)|c)).collect();
        for p in infer_parameters(&coefficients, &gens) {
            let covered = parameters.iter().any(|q| crate::traverse::has(q, &p));
            if !covered {
                parameters.push(p);
            }
        }
        let integral_side=|p:&Self|if let Some(f)=&p.field {matches!(f.coefficient_domain,super::domains::ExactDomain::Functions{integral:true,..})}else{p.domain==super::polyobj::Domain::ZZ};
        let ring_side=|p:&Self|p.field.as_ref().is_none_or(|f|matches!(f.coefficient_domain,super::domains::ExactDomain::Functions{ring:true,..}));
        let composite=|p:&Self|p.field.is_some()||!p.base.domain_gens().is_empty();
        let ring=ring_side(self)&&ring_side(o);
        // SymPy's `Domain.unify`: a fraction field unified with a polynomial
        // ring over a non-field ground takes the ring of the unified ground
        // (`QQ(a)` and `ZZ[b]` give `ZZ(a,b)`).
        let integral=(integral_side(self)&&integral_side(o))
            ||(composite(self)&&composite(o)&&ring_side(self)!=ring_side(o)&&(integral_side(self)||integral_side(o)));
        let order=if let super::domains::ExactDomain::Functions{order,..}=&prototype.coefficient_domain{order.as_str()}else{"lex"};
        let descriptor=super::domains::ExactDomain::functions(integral,ring,parameters.clone(),order)?;
        let display=descriptor.to_string();
        if parameters.iter().any(|p|gens.contains(p)){return Err(PolyError::Failed("polynomial generators and coefficient parameters must be disjoint".into()));}
        let _scope=Scope::new(parameters.len());
        let convert=|p:&Self|{
            let positions:Vec<_>=p.true_gens().iter().map(|g|gens.iter().position(|candidate|candidate==g).unwrap()).collect();
            let terms=p.terms().into_iter().map(|(powers,c)|{
                let mut remapped=vec![0;gens.len()];
                for (index,exponent) in positions.iter().zip(powers){remapped[*index]=exponent;}
                Ok((remapped,coefficient(&c,&parameters)?))
            }).collect::<PResult<Vec<_>>>()?;
            Ok(FractionPoly{rep:MPoly::from_terms(gens.len(),terms),gens:gens.clone(),parameters:parameters.clone(),display:display.clone(),coefficient_domain:descriptor.clone()})
        };
        Ok(Some((convert(self)?,convert(o)?)))
    })}
    /// Operands over a nested coefficient ring (`ZZ[x][y]`): a plain integer
    /// or rational operand joins it (over `QQ[x][y]` for rationals), as the
    /// reference `Domain.unify` does; other mixtures are not implemented.
    fn composite_pair(&self,o:&Self,gens:&[Expr])->PResult<Option<(FractionPoly,FractionPoly)>> {
        use super::domains::ExactDomain;
        let composite_of=|p:&Self|p.field.as_ref().and_then(|f|matches!(f.coefficient_domain,ExactDomain::Composite{..}).then(||f.coefficient_domain.clone()));
        let Some(mut target)=composite_of(self).or_else(||composite_of(o)) else {return Ok(None)};
        for p in [self,o] {
            if let Some(d)=composite_of(p) {
                if d!=target {
                    if d.with_rational_field()==target.with_rational_field() {target=target.with_rational_field();}
                    else {return Err(PolyError::Unsupported(format!("cannot unify the nested coefficient domains {d} and {target}")));}
                }
            } else if p.field.is_none()&&p.base.domain_gens().is_empty()&&p.modulus.is_none()&&p.algebraic.is_none()&&p.expression.is_none()&&!p.domain.is_rr() {
                if p.domain==super::polyobj::Domain::QQ {target=target.with_rational_field();}
            } else {
                return Err(PolyError::Unsupported(format!("cannot unify the nested coefficient domain {target} with another coefficient domain")));
            }
        }
        let Some(ExactDomain::Functions{parameters,..})=target.flat_domain() else {return Err(PolyError::Unsupported(format!("polynomials over {target} are not implemented")))};
        if parameters.iter().any(|p|gens.contains(p)){return Err(PolyError::Generator("polynomial ring and it's ground domain share generators".into()));}
        let display=target.to_string();
        let _scope=Scope::new(parameters.len());
        let convert=|p:&Self|{
            let positions:Vec<_>=p.true_gens().iter().map(|g|gens.iter().position(|candidate|candidate==g).unwrap()).collect();
            let terms=p.terms().into_iter().map(|(powers,c)|{
                let mut remapped=vec![0;gens.len()];
                for (index,exponent) in positions.iter().zip(powers){remapped[*index]=exponent;}
                Ok((remapped,coefficient(&c,&parameters)?))
            }).collect::<PResult<Vec<_>>>()?;
            Ok::<_,PolyError>(FractionPoly{rep:MPoly::from_terms(gens.len(),terms),gens:gens.to_vec(),parameters:parameters.clone(),display:display.clone(),coefficient_domain:target.clone()})
        };
        Ok(Some((convert(self)?,convert(o)?)))
    }
    fn pair(&self,o:&Self)->Option<(FractionPoly,FractionPoly)>{self.checked_fraction_pair(o).expect("compatible fraction coefficient domains")}

    /// Whether two non-`EX` polynomials unify to the plain field `CC`
    /// (a legacy float `CC` polynomial, or `RR` with `ZZ_I`/`QQ_I`), whose
    /// arithmetic runs on the scalar `CC` ground.
    fn complex_field_pair(&self,o:&Self)->bool{
        if self.modulus.is_some()||o.modulus.is_some(){return false;}
        let (Some(a),Some(b))=(super::gaussian_ground::descriptor(self),super::gaussian_ground::descriptor(o)) else {return false};
        if !a.symbols.is_empty()||!b.symbols.is_empty(){return false;}
        let gens=super::polyobj::unify_gens(self.true_gens(),o.true_gens());
        matches!(super::gaussian_ground::unify(&a,&b,&gens),Ok(Some(g)) if g.is_scalar())
    }
    /// The reference unification with `EX`: every other domain converts
    /// into it coefficient by coefficient, without simplification.
    pub fn checked_expression_pair(&self,o:&Self)->PResult<Option<(ExpressionPoly,ExpressionPoly)>>{super::mpoly::checked_operation(||{
        if self.expression.is_none()&&o.expression.is_none()&&!self.complex_field_pair(o){return Ok(None);}
        if let Some(m)=self.modulus.or(o.modulus){
            let message=if self.expression.is_some(){format!("Cannot unify EX with GF({m})")}else{format!("Cannot unify GF({m}) with EX")};
            return Err(PolyError::Reference{kind:"UnificationFailed",message});
        }
        let gens=super::polyobj::unify_gens(self.true_gens(),o.true_gens());
        // A composite Gaussian ground unifies as the reference domains do.
        if let (Some(a),Some(b))=(super::gaussian_ground::descriptor(self),super::gaussian_ground::descriptor(o)) {
            if let Some(ground)=super::gaussian_ground::unify(&a,&b,&gens)? {
                let convert=|p:&Self|->PResult<ExpressionPoly>{
                    let e=match &p.expression {Some(e)=>e.reorder(&gens)?,None=>ExpressionPoly::from_terms(p.terms(),p.true_gens())?.reorder(&gens)?};
                    ground.convert_poly(&e)
                };
                return Ok(Some((convert(self)?,convert(o)?)));
            }
        }
        let convert=|p:&Self|->PResult<ExpressionPoly>{
            let mut e=match &p.expression {
                Some(e)=>e.reorder(&gens)?,
                None=>ExpressionPoly::from_terms(p.terms(),p.true_gens())?.reorder(&gens)?,
            };
            e.ground=None;
            Ok(e)
        };
        Ok(Some((convert(self)?,convert(o)?)))
    })}
    pub fn to_expression(&self)->PResult<ExpressionPoly>{super::mpoly::checked_operation(||{
        match &self.expression {Some(e)=>Ok(e.clone()),None=>ExpressionPoly::from_terms(self.terms(),self.true_gens())}
    })}
    fn expression_pair(&self,o:&Self)->Option<(ExpressionPoly,ExpressionPoly)>{self.checked_expression_pair(o).expect("EX unification")}
    pub fn add(&self, o: &Self) -> Self {
        if let Some((a,b))=self.expression_pair(o){return a.add(&b).expect("EX arithmetic").into();}
        if let Some((a,b))=self.algebraic_pair(o){return a.add(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.add(&b).into()
        } else {
            self.base.add(&o.base).into()
        }
    }
    pub fn sub(&self, o: &Self) -> Self {
        if let Some((a,b))=self.expression_pair(o){return a.sub(&b).expect("EX arithmetic").into();}
        if let Some((a,b))=self.algebraic_pair(o){return a.sub(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.sub(&b).into()
        } else {
            self.base.sub(&o.base).into()
        }
    }
    pub fn mul(&self, o: &Self) -> Self {
        if let Some((a,b))=self.expression_pair(o){return a.mul(&b).expect("EX arithmetic").into();}
        if let Some((a,b))=self.algebraic_pair(o){return a.mul(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.mul(&b).into()
        } else {
            self.base.mul(&o.base).into()
        }
    }
    pub fn neg(&self) -> Self {
        if let Some(p)=&self.expression {return p.neg().into();}
        if let Some(p)=&self.algebraic {return p.neg().into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.neg().into(), |p| p.neg().into())
    }
    pub fn pow(&self, k: u32) -> Self {
        if let Some(p)=&self.expression {return p.pow(k).expect("EX arithmetic").into();}
        if let Some(p)=&self.algebraic {return p.pow(k).into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.pow(k).into(), |p| p.pow(k).into())
    }
    pub fn diff(&self, i: usize) -> Self {
        if let Some(p)=&self.expression {return p.diff(1,i).expect("EX arithmetic").into();}
        if let Some(p)=&self.algebraic {return p.diff(i).into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.diff(i).into(), |p| p.diff(i).into())
    }
    pub fn integrate(&self, i: usize) -> Self {
        if let Some(p)=&self.expression {return p.integrate(1,i).expect("EX arithmetic").into();}
        if let Some(p)=&self.algebraic {return p.integrate(i).into();}
        self.field
            .as_ref()
            .map_or_else(|| self.base.integrate(i).into(), |p| p.integrate(i).into())
    }
    pub fn div(&self,o:&Self)->PResult<(Self,Self)>{super::mpoly::checked_operation(||{self.div_with_options(o,true)})}
    pub fn div_with_options(&self,o:&Self,auto:bool)->PResult<(Self,Self)>{super::mpoly::checked_operation(||{
        if let Some((a,b))=self.checked_expression_pair(o)?{let(q,r)=a.div(&b)?;return Ok((q.into(),r.into()));}
        if let Some((a,b))=self.checked_algebraic_pair(o)?{let(q,r)=if auto{a.div(&b)?}else{a.ring_div(&b)?};return Ok((q.into(),r.into()));}
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
        match self.base.div(&o.base) {
            Ok((q,r))=>Ok((q.into(),r.into())),
            Err(error) if auto && matches!(&error,PolyError::Unsupported(message) if message.contains("non-polynomial quotient"))=>{
                // Inferred radical coefficients can live in the legacy
                // polynomial view without an explicit extension=True request.
                // Their quotient needs the existing exact number-field engine.
                let mut parameters=self.base.domain_gens().to_vec();
                for parameter in o.base.domain_gens() {
                    if !parameters.contains(parameter) {parameters.push(parameter.clone());}
                }
                if parameters.is_empty() {
                    return Err(error);
                }
                let number_field=if parameters.iter().any(|p|!crate::traverse::free_symbols(p).is_empty()) {None} else {super::algebraic::NumberField::from_generators(&parameters).ok()};
                let Some(field)=number_field else {
                    // `Poly.div(auto=True)` over `ZZ[log(a)]`: the division
                    // runs over the fraction field `ZZ(log(a))`, and both
                    // results return to the ring when they belong to it.
                    let Some((a,b))=self.legacy_ring_pair(o,&parameters)? else {return Err(error)};
                    return a.div_with_options(&b,auto);
                };
                let gens=super::polyobj::unify_gens(self.true_gens(),o.true_gens());
                let a=super::algebraic::AlgebraicPoly::from_expr(&self.as_expr(),&gens,field.clone())?;
                let b=super::algebraic::AlgebraicPoly::from_expr(&o.as_expr(),&gens,field)?;
                let(q,r)=a.div(&b)?;
                Ok((q.into(),r.into()))
            },
            Err(error)=>Err(error),
        }
    })}
    pub fn pdiv(&self,o:&Self)->PResult<(Self,Self)>{super::mpoly::checked_operation(||{
        if let Some((a,b))=self.checked_expression_pair(o)?{let(q,r)=a.pdiv(&b)?;return Ok((q.into(),r.into()));}
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
    })}
    pub fn prem(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{Ok(self.pdiv(o)?.1)})}
    pub fn pquo(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{Ok(self.pdiv(o)?.0)})}
    pub fn pexquo(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{if let Some((a,b))=self.checked_expression_pair(o)?{return Ok(a.pexquo(&b)?.into());}let(q,r)=self.pdiv(o)?;if !r.is_zero(){return Err(PolyError::ExactQuotientFailed("pseudo-division has a nonzero remainder".into()));}Ok(q)})}
    pub fn rem(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{self.rem_with_options(o,true)})}
    pub fn quo(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{self.quo_with_options(o,true)})}
    /// `self` and `o`, two polynomials of the legacy view over coefficient
    /// generators such as `ZZ[log(a)]`, as polynomials over the typed ring
    /// of `parameters`.
    fn legacy_ring_pair(&self,o:&Self,parameters:&[Expr])->PResult<Option<(Self,Self)>>{
        let legacy=|p:&Self|p.field.is_none()&&p.algebraic.is_none()&&p.expression.is_none()&&p.modulus.is_none()&&!p.base.domain.is_rr();
        if !legacy(self)||!legacy(o){return Ok(None);}
        let integral=self.base.domain.max(o.base.domain)==super::polyobj::Domain::ZZ;
        let descriptor=super::domains::ExactDomain::functions(integral,true,parameters.to_vec(),"lex")?;
        let display=descriptor.to_string();
        let gens=super::polyobj::unify_gens(self.true_gens(),o.true_gens());
        let convert=|p:&Self|->PResult<Self>{
            let mut f=FractionPoly::from_expr(&p.as_expr(),&gens,parameters,&display)?;
            f.coefficient_domain=descriptor.clone();
            Ok(f.into())
        };
        Ok(Some((convert(self)?,convert(o)?)))
    }
    fn retract_division_component(&self,o:&Self,mut result:Self)->Self{
        if result.expression.is_some(){return result;}
        if let (Some(p),Ok(Some((a,_))))=(result.algebraic.as_ref(),self.checked_algebraic_pair(o)){
            let p=p.retract_to(&a);
            return p.into();
        }
        if let Some((a,_))=self.pair(o){
            if a.display.ends_with(']') {if let Some(p)=result.field.as_mut(){if p.rep.coeffs().iter().all(|c|coefficient_in_ring(c,a.display.starts_with("ZZ"))){p.display=a.display;p.coefficient_domain=a.coefficient_domain;}}}
        }else if let Some(p)=result.field.as_mut(){
            // legacy `ZZ[log(a)]` operands divided over `ZZ(log(a))`
            if let super::domains::ExactDomain::Functions{integral,ring:false,parameters,order}=&p.coefficient_domain {
                if p.rep.coeffs().iter().all(|c|coefficient_in_ring(c,*integral)) {
                    if let Ok(ring)=super::domains::ExactDomain::functions(*integral,true,parameters.clone(),order) {
                        p.display=ring.to_string();
                        p.coefficient_domain=ring;
                    }
                }
            }
        }else if self.algebraic.is_none()&&o.algebraic.is_none()&&self.modulus.is_none()&&o.modulus.is_none()&&self.domain.max(o.domain)==super::polyobj::Domain::ZZ&&result.base.rep.is_integral(){result.base.domain=super::polyobj::Domain::ZZ;}
        result
    }
    pub fn rem_with_options(&self,o:&Self,auto:bool)->PResult<Self>{super::mpoly::checked_operation(||{let r=self.div_with_options(o,auto)?.1;Ok(self.retract_division_component(o,r))})}
    pub fn quo_with_options(&self,o:&Self,auto:bool)->PResult<Self>{super::mpoly::checked_operation(||{let q=self.div_with_options(o,auto)?.0;Ok(self.retract_division_component(o,q))})}
    pub fn exquo(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{self.exquo_with_options(o,true)})}
    pub fn exquo_with_options(&self,o:&Self,auto:bool)->PResult<Self>{super::mpoly::checked_operation(||{
        if let Some((a,b))=self.checked_expression_pair(o)?{return Ok(a.exquo(&b)?.into());}
        let(q,r)=self.div_with_options(o,auto)?;
        if !r.is_zero(){return Err(PolyError::ExactQuotientFailed("polynomial division is not exact in the selected coefficient domain".into()));}Ok(q)
    })}
    pub fn checked_gcd(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{
        if let Some((a,b))=self.checked_expression_pair(o)?{
            if let Some(g)=gaussian_ground_of(&a,&b){return Ok(g.gcd(&a,&b)?.into());}
            return Ok(a.gcd(&b)?.into());
        }
        Ok(self.gcd(o))
    })}
    pub fn checked_lcm(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{
        if let Some((a,b))=self.checked_expression_pair(o)?{
            if let Some(g)=gaussian_ground_of(&a,&b){return Ok(g.lcm(&a,&b)?.into());}
            return Ok(a.lcm(&b)?.into());
        }
        Ok(self.lcm(o))
    })}
    pub fn checked_cofactors(&self,o:&Self)->PResult<(Self,Self,Self)>{super::mpoly::checked_operation(||{
        if let Some((a,b))=self.checked_expression_pair(o)?{
            let(h,p,q)=match gaussian_ground_of(&a,&b){Some(g)=>g.cofactors(&a,&b)?,None=>a.cofactors(&b)?};
            return Ok((h.into(),p.into(),q.into()));
        }
        Ok(self.cofactors(o))
    })}
    /// Arithmetic that reports `EX` normalization failures instead of
    /// panicking; the other domains are infallible here.
    pub fn checked_add(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{if let Some((a,b))=self.checked_expression_pair(o)?{return Ok(a.add(&b)?.into());}Ok(self.add(o))})}
    pub fn checked_sub(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{if let Some((a,b))=self.checked_expression_pair(o)?{return Ok(a.sub(&b)?.into());}Ok(self.sub(o))})}
    pub fn checked_mul(&self,o:&Self)->PResult<Self>{super::mpoly::checked_operation(||{if let Some((a,b))=self.checked_expression_pair(o)?{return Ok(a.mul(&b)?.into());}Ok(self.mul(o))})}
    pub fn checked_pow(&self,k:u32)->PResult<Self>{super::mpoly::checked_operation(||{if let Some(p)=&self.expression{return Ok(p.pow(k)?.into());}Ok(self.pow(k))})}
    pub fn gcd(&self, o: &Self) -> Self {
        if let Some((a,b))=self.expression_pair(o){
            if let Some(g)=gaussian_ground_of(&a,&b){if let Ok(h)=g.gcd(&a,&b){return h.into();}}
            return a.gcd(&b).expect("EX arithmetic").into();
        }
        if let Some((a,b))=self.algebraic_pair(o){return a.gcd(&b).into();}
        if let Some((a, b)) = self.pair(o) {
            a.gcd(&b).into()
        } else {
            self.base.gcd(&o.base).into()
        }
    }
    pub fn lcm(&self, o: &Self) -> Self {
        if let Some((a,b))=self.expression_pair(o){
            if let Some(g)=gaussian_ground_of(&a,&b){if let Ok(h)=g.lcm(&a,&b){return h.into();}}
            return a.lcm(&b).expect("EX arithmetic").into();
        }
        if let Some((a,b))=self.algebraic_pair(o){return a.lcm(&b).into();}
        if let Some((a,b))=self.pair(o){
            a.lcm(&b).into()
        } else {
            self.base.lcm(&o.base).into()
        }
    }
    pub fn cofactors(&self, o: &Self) -> (Self, Self, Self) {
        if let Some((a,b))=self.expression_pair(o){
            if let Some(g)=gaussian_ground_of(&a,&b){if let Ok((h,p,q))=g.cofactors(&a,&b){return (h.into(),p.into(),q.into());}}
            let(h,p,q)=a.cofactors(&b).expect("EX arithmetic");return (h.into(),p.into(),q.into());
        }
        if let Some((a,b))=self.algebraic_pair(o){let(h,a,b)=a.cofactors(&b);return (h.into(),a.into(),b.into());}
        let g = self.gcd(o);
        if g.is_zero() {
            return (g, self.clone(), o.clone());
        }
        let a = self.exquo_with_options(&g,false).expect("gcd divides first input in its ring");
        let b = o.exquo_with_options(&g,false).expect("gcd divides second input in its ring");
        (g, a, b)
    }
    pub fn to_sympy_string(&self) -> String {
        if let Some(p)=&self.expression {return super::polyobj::format_poly(&p.gens,p.terms(),None,&self.domain_str());}
        // Fraction-field coefficients print as one quotient, `(a**2 + 1)/a`.
        if let Some(f)=&self.field {return super::polyobj::format_poly(&f.gens,f.terms(),self.modulus,&self.domain_str());}
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
    fn coefficient_unification_preserves_formal_composite_generator_slots(){
        let x=Expr::symbol("x");let a=Expr::symbol("a");
        let source:Polynomial=FractionPoly::from_expr(&parse_str("a*x+1").unwrap(),&[x.clone()],&[a.clone()],"ZZ[a]").unwrap().into();
        let generator=parse_str("y+1").unwrap();
        let p=super::super::structure::rename_gens(&source,&[generator.clone()]).unwrap();
        let square=p.mul(&p);
        assert_eq!(square.true_gens(),&[generator]);
        assert_eq!(square.domain_str(),"ZZ[a]");
        assert_eq!(square.terms(),vec![(vec![2],parse_str("a^2").unwrap()),(vec![1],parse_str("2*a").unwrap()),(vec![0],Expr::one())]);
        let other:Polynomial=FractionPoly::from_expr(&parse_str("a*x+2").unwrap(),&[x],&[a],"QQ(a)").unwrap().into();
        let other=super::super::structure::rename_gens(&other,&[parse_str("z+1").unwrap()]).unwrap();
        let product=p.mul(&other);
        // The reference unifies ZZ[a] with QQ(a) into ZZ(a).
        assert_eq!(product.domain_str(),"ZZ(a)");
        assert_eq!(product.terms().len(),4);
        assert_eq!(product.true_gens(),&[parse_str("y+1").unwrap(),parse_str("z+1").unwrap()]);
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

#[cfg(test)]
mod reference_euclid_and_cancel {
    use super::*;
    use crate::parser::parse_str;
    fn poly(source: &str, domain: &str) -> FractionPoly {
        FractionPoly::from_expr(&parse_str(source).unwrap(), &[Expr::symbol("x")], &[Expr::symbol("a")], domain).unwrap()
    }
    fn show(p: &FractionPoly) -> String {
        let p: Polynomial = p.clone().into();
        p.to_sympy_string()
    }
    #[test]
    fn gcdex_half_gcdex_and_invert_match_the_reference() {
        let f = poly("x^2+a", "QQ(a)");
        let g = poly("x+a", "QQ(a)");
        let (s, t, h) = f.gcdex(&g).unwrap();
        assert_eq!(show(&s), "Poly(1/(a**2 + a), x, domain='QQ(a)')");
        assert_eq!(show(&t), "Poly(-1/(a**2 + a)*x + 1/(a + 1), x, domain='QQ(a)')");
        assert_eq!(show(&h), "Poly(1, x, domain='QQ(a)')");
        let (s2, h2) = f.half_gcdex(&g).unwrap();
        assert_eq!(s2.rep, s.rep);
        assert_eq!(h2.rep, h.rep);
        assert_eq!(g.invert(&f).unwrap().rep, t.rep);
        // A common factor is a zero divisor, not an absent inverse.
        let d = poly("x^2-a^2", "QQ(a)");
        assert!(matches!(g.invert(&d), Err(PolyError::NotInvertible(_))));
    }
    #[test]
    fn cancel_and_content_follow_the_reference_normalization() {
        let p = poly("x^2/a-a", "ZZ(a)");
        let q = poly("x/a+1", "ZZ(a)");
        let (cp, cq, r, s) = p.reference_cancel(&q);
        let _scope = Scope::new(1);
        assert_eq!(cp, cq);
        assert_eq!(show(&r), "Poly(x - a, x, domain='ZZ(a)')");
        assert_eq!(show(&s), "Poly(1, x, domain='ZZ(a)')");
        let parameters = [Expr::symbol("a")];
        assert_eq!(coefficient_expr(&p.reference_content(), &parameters), parse_str("1/a").unwrap());
        let (c, primitive) = p.reference_primitive();
        assert_eq!(coefficient_expr(&c, &parameters), parse_str("1/a").unwrap());
        assert_eq!(show(&primitive), "Poly(x**2 - a**2, x, domain='ZZ(a)')");
    }
}

/// The common composite Gaussian ground of a unified `EX` pair.
fn gaussian_ground_of(a:&ExpressionPoly,b:&ExpressionPoly)->Option<super::gaussian_ground::GaussianGround>{
    match (&a.ground,&b.ground) {(Some(x),Some(y)) if x==y=>Some(x.clone()),_=>None}
}
