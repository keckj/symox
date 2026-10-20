//! Expression level polynomial functions with SymPy semantics: `factor`,
//! `factor_list`, `sqf_list`, `gcd`, `lcm`, `cancel`, `div`/`rem`/`quo`,
//! `degree`, `resultant`, `discriminant`, `apart`, `collect`, `together`.

use super::convert::{number_to_rat, PResult, PolyError};
use super::dense;
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use super::polyobj::{sorted_factors, Poly};
use super::resultant::{discriminant_mz, resultant_mz};
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, as_coeff_mul, div, mul, mul_factor, neg, pow};
use crate::traverse::map_children;

// ----------------------------------------------------------------------
// SymPy's together / gcd_terms
// ----------------------------------------------------------------------

/// `Add.primitive()`: `(rational content, primitive sum)`.
pub fn add_primitive(e: &Expr) -> (Number, Expr) {
    let Kind::Add(c, terms) = e.kind() else {
        let (k, r) = as_coeff_mul(e);
        if k.is_rational() {
            return (k.abs(), if k.is_negative() { neg(&r) } else { r });
        }
        return (Number::ONE, e.clone());
    };
    let mut g: Option<(Int, Int)> = None; // (gcd of numerators, lcm of denominators)
    let mut push = |n: &Number| {
        if let Some((p, q)) = n.numer_denom() {
            let (p, q) = (Int::from_big(p), Int::from_big(q));
            g = Some(match g.take() {
                None => (p.abs(), q),
                Some((a, b)) => (a.gcd(&p), b.lcm(&q)),
            });
        } else {
            g = Some(match g.take() {
                None => (Int::one(), Int::one()),
                Some((a, b)) => (a.gcd(&Int::one()), b),
            });
        }
    };
    for (_, k) in terms {
        push(k);
    }
    if !c.is_zero() {
        push(c);
    }
    let (p, q) = g.unwrap();
    if p.is_one() && q.is_one() {
        return (Number::ONE, e.clone());
    }
    let content = Rat::new(p, q).to_number();
    let inv = Number::ONE.div(&content).unwrap();
    (content, mul(&[Expr::num(inv), e.clone()]))
}

#[derive(Clone, Debug, Default)]
struct Factors(Vec<(Expr, i64)>);

impl Factors {
    fn add(&mut self, b: Expr, k: i64) {
        if let Some(x) = self.0.iter_mut().find(|x| x.0 == b) {
            x.1 += k;
        } else {
            self.0.push((b, k));
        }
        self.0.retain(|x| x.1 != 0);
    }
    fn get(&self, b: &Expr) -> i64 {
        self.0.iter().find(|x| &x.0 == b).map(|x| x.1).unwrap_or(0)
    }
    fn gcd(&self, o: &Factors) -> Factors {
        let mut r = Factors::default();
        for (b, k) in &self.0 {
            let m = (*k).min(o.get(b));
            if m > 0 {
                r.0.push((b.clone(), m));
            }
        }
        r
    }
    fn lcm(&self, o: &Factors) -> Factors {
        let mut r = self.clone();
        for (b, k) in &o.0 {
            let cur = r.get(b);
            if *k > cur {
                r.add(b.clone(), k - cur);
            }
        }
        r
    }
    fn mul(&self, o: &Factors) -> Factors {
        let mut r = self.clone();
        for (b, k) in &o.0 {
            r.add(b.clone(), *k);
        }
        r
    }
    fn quo(&self, o: &Factors) -> Factors {
        let mut r = self.clone();
        for (b, k) in &o.0 {
            r.add(b.clone(), -*k);
        }
        r
    }
    fn as_expr(&self) -> Expr {
        let f: Vec<Expr> = self.0.iter().map(|(b, k)| pow(b, &Expr::int(*k))).collect();
        mul(&f)
    }
}

#[derive(Clone, Debug)]
struct Term {
    coeff: Number,
    numer: Factors,
    denom: Factors,
}

impl Term {
    fn new(t: &Expr) -> Term {
        let (mut coeff, rest) = as_coeff_mul(t);
        let mut numer = Factors::default();
        let mut denom = Factors::default();
        let facs: Vec<Expr> = match rest.kind() {
            Kind::Mul(_, fs) => fs.iter().map(|(b, x)| mul_factor(b, x)).collect(),
            Kind::Num(n) if n.is_one() => vec![],
            _ => vec![rest.clone()],
        };
        for f in facs {
            let (mut base, k) = super::convert::decompose_factor_signed(&f);
            if base.is_add() {
                let (cont, prim) = add_primitive(&base);
                if !cont.is_one() {
                    coeff = coeff.mul(&cont.pow_i64(k).unwrap());
                    base = prim;
                }
            }
            if k > 0 {
                numer.add(base, k);
            } else {
                denom.add(base, -k);
            }
        }
        Term { coeff, numer, denom }
    }

    fn gcd(&self, o: &Term) -> Term {
        Term {
            coeff: self.coeff.gcd(&o.coeff).abs(),
            numer: self.numer.gcd(&o.numer),
            denom: self.denom.gcd(&o.denom),
        }
    }

    fn quo(&self, o: &Term) -> Term {
        let coeff = self.coeff.div(&o.coeff).unwrap_or(Number::ONE);
        let n = self.numer.mul(&o.denom);
        let d = self.denom.mul(&o.numer);
        // cancel common bases
        let common = n.gcd(&d);
        Term { coeff, numer: n.quo(&common), denom: d.quo(&common) }
    }

    fn as_expr(&self) -> Expr {
        mul(&[Expr::num(self.coeff.clone()), self.numer.as_expr(), pow(&self.denom.as_expr(), &Expr::minus_one())])
    }
}

/// SymPy's `gcd_terms` on a list of terms (fraction=True).
pub(crate) fn gcd_terms(terms: &[Expr]) -> Expr {
    let terms: Vec<Term> = terms.iter().filter(|t| !t.is_zero()).map(Term::new).collect();
    if terms.is_empty() {
        return Expr::zero();
    }
    let (cont, numer, denom) = if terms.len() == 1 {
        let t = &terms[0];
        (Expr::num(t.coeff.clone()), t.numer.as_expr(), t.denom.as_expr())
    } else {
        let mut c = terms[0].clone();
        for t in &terms[1..] {
            c = c.gcd(t);
        }
        let qs: Vec<Term> = terms.iter().map(|t| t.quo(&c)).collect();
        let mut den = qs[0].denom.clone();
        for t in &qs[1..] {
            den = den.lcm(&t.denom);
        }
        let numers: Vec<Expr> = qs
            .iter()
            .map(|t| mul(&[Expr::num(t.coeff.clone()), t.numer.mul(&den.quo(&t.denom)).as_expr()]))
            .collect();
        (c.as_expr(), add(&numers), den.as_expr())
    };
    let (cont, numer) = if numer.is_add() {
        let (c2, p) = add_primitive(&numer);
        (mul(&[cont, Expr::num(c2)]), p)
    } else {
        (cont, numer)
    };
    let (coeff, factors) = as_coeff_mul(&cont);
    keep_coeff(&Expr::num(coeff), &mul(&[factors, numer, pow(&denom, &Expr::minus_one())]))
}

/// SymPy's `together` (not deep).
pub fn together(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) | Kind::Func(..) => e.clone(),
        Kind::Add(..) => {
            let args: Vec<Expr> = e.args().iter().map(together).collect();
            gcd_terms(&args)
        }
        Kind::Pow(b, x) => pow(&together(b), x),
        Kind::Mul(..) => {
            let args: Vec<Expr> = e.args().iter().map(together).collect();
            mul(&args)
        }
        _ => map_children(e, &mut |c| together(c)),
    }
}

/// SymPy's (recursive) `as_numer_denom`.
pub fn as_numer_denom(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Num(n) => (Expr::num(n.numer()), Expr::num(n.denom())),
        Kind::Mul(c, facs) => {
            let mut num = vec![Expr::num(c.numer())];
            let mut den = vec![Expr::num(c.denom())];
            for (b, x) in facs {
                let (n, d) = as_numer_denom(&mul_factor(b, x));
                num.push(n);
                den.push(d);
            }
            (mul(&num), mul(&den))
        }
        Kind::Pow(b, x) => {
            let int_exp = x.is_integer_number();
            let neg_exp = match x.as_number() {
                Some(n) => n.is_negative(),
                None => crate::ops::could_extract_minus_sign(x),
            };
            let (n, d) = if int_exp { as_numer_denom(b) } else { (b.clone(), Expr::one()) };
            let (n, d) = if !int_exp {
                // only a numeric denominator can be separated safely
                match b.as_number() {
                    Some(r) => (Expr::num(r.numer()), Expr::num(r.denom())),
                    None => (n, d),
                }
            } else {
                (n, d)
            };
            if neg_exp {
                let nx = neg(x);
                (pow(&d, &nx), pow(&n, &nx))
            } else {
                (pow(&n, x), pow(&d, x))
            }
        }
        Kind::Func(Func::Exp, a) if crate::ops::could_extract_minus_sign(&a[0]) => {
            (Expr::one(), crate::functions::exp(&neg(&a[0])))
        }
        Kind::Add(..) => {
            let parts: Vec<(Expr, Expr)> = e.args().iter().map(as_numer_denom).collect();
            if parts.iter().all(|(_, d)| d.is_one()) {
                return (e.clone(), Expr::one());
            }
            // group numerators by denominator
            let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
            for (n, d) in parts {
                if let Some(g) = groups.iter_mut().find(|g| g.0 == d) {
                    g.1.push(n);
                } else {
                    groups.push((d, vec![n]));
                }
            }
            let dens: Vec<Expr> = groups.iter().map(|g| g.0.clone()).collect();
            let mut nums = Vec::with_capacity(groups.len());
            for (i, (_, ns)) in groups.iter().enumerate() {
                let mut f = vec![add(ns)];
                for (j, d) in dens.iter().enumerate() {
                    if j != i {
                        f.push(d.clone());
                    }
                }
                nums.push(mul(&f));
            }
            (add(&nums), mul(&dens))
        }
        _ => (e.clone(), Expr::one()),
    }
}

// ----------------------------------------------------------------------
// factor / factor_list / sqf
// ----------------------------------------------------------------------

/// A factor of a symbolic factorisation.
#[derive(Clone, Debug)]
pub enum FItem {
    P(Poly),
    E(Expr),
}

impl FItem {
    pub fn as_expr(&self) -> Expr {
        match self {
            FItem::P(p) => p.as_expr(),
            FItem::E(e) => e.clone(),
        }
    }
}

fn is_pure_complex(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) => true,
        Kind::Const(Constant::I) => true,
        Kind::Mul(c, facs) => c.is_rational() && facs.len() == 1 && facs[0].0.is_constant(Constant::I) && facs[0].1.is_one(),
        Kind::Add(c, terms) => c.is_rational() && terms.len() == 1 && is_pure_complex(&crate::ops::add_term(&terms[0].0, &terms[0].1)) && terms[0].0.is_constant(Constant::I),
        _ => false,
    }
}

fn rat_expr(c: &Rat, float: bool) -> Expr {
    super::convert::rat_to_expr_float(c, float)
}

/// SymPy's `_symbolic_factor_list`: `(coeff, [(factor, exponent)])`.
fn symbolic_factor_list(e: &Expr, sqf: bool, modulus: Option<u64>) -> PResult<(Expr, Vec<(FItem, Expr)>)> {
    symbolic_factor_list_with(e,sqf,modulus,&[],None)
}
fn symbolic_factor_list_with(e:&Expr,sqf:bool,modulus:Option<u64>,gens:&[Expr],checked:Option<super::factor_budget::Limits>)->PResult<(Expr,Vec<(FItem,Expr)>)> {
    let mut coeff = Expr::one();
    let mut factors: Vec<(FItem, Expr)> = Vec::new();
    let args: Vec<Expr> = match e.kind() {
        Kind::Mul(..) => {
            // SymPy's Mul.args order (ties in the factor sort keep it)
            let mut a = e.args();
            crate::printing::order::sort_factors(&mut a);
            a
        }
        _ => vec![e.clone()],
    };
    for arg in args {
        if is_pure_complex(&arg) {
            coeff = mul(&[coeff, arg]);
            continue;
        }
        let (base, exp) = match arg.kind() {
            Kind::Pow(b, x) => {
                if b.is_number() && x.is_number() {
                    coeff = mul(&[coeff, arg.clone()]);
                    continue;
                }
                if b.is_number() {
                    factors.push((FItem::E(b.clone()), x.clone()));
                    continue;
                }
                (b.clone(), x.clone())
            }
            _ => (arg.clone(), Expr::one()),
        };
        let poly = match Poly::from_expr(&base, gens) {
            Ok(p) => p,
            Err(_) => {
                // Cancellation may remove every polynomial generator. Keep
                // the resulting scalar, rather than the unfactored original
                // sum, when selecting the non-polynomial fallback.
                let scalar = expand(&base);
                if is_pure_complex(&scalar) {
                    coeff = mul(&[coeff, pow(&scalar, &exp)]);
                    continue;
                }
                factors.push((FItem::E(base), exp));
                continue;
            }
        };
        let poly = match modulus {
            Some(m) => {
                if poly.gens.len() > 1 {
                    return Err(PolyError::Unsupported("multivariate polynomials over finite fields".into()));
                }
                poly.set_modulus(m)
            }
            None => poly,
        };
        let(ce,positive,fs)=if let Some(limits)=checked {
            let result=super::checked_factors::compatible_list(&poly.into(),sqf,false,limits)?;
            let coefficient=result.unit.as_expr();
            let positive=coefficient.as_number().is_some_and(|n|n.is_positive());
            (coefficient,positive,result.factors.into_iter().map(|(p,k)|(p.base,k)).collect())
        }else{
            let float = poly.is_float();
            let (c, fs) = if sqf { poly.sqf_list() } else { poly.factor_list() };
            (rat_expr(&c,float),c.n.is_positive(),fs)
        };
        let mut fs: Vec<(FItem, u32)> = fs.into_iter().map(|(p, k)| (FItem::P(p), k)).collect();
        if !ce.is_one() {
            if exp.is_integer_number() {
                coeff = mul(&[coeff, pow(&ce, &exp)]);
            } else if positive {
                factors.push((FItem::E(ce), exp.clone()));
            } else {
                fs.push((FItem::E(ce), 1));
            }
        }
        if exp.is_one() {
            factors.extend(fs.into_iter().map(|(f, k)| (f, Expr::int(k as i64))));
        } else if exp.is_integer_number() {
            factors.extend(fs.into_iter().map(|(f, k)| (f, mul(&[Expr::int(k as i64), exp.clone()]))));
        } else {
            let mut other = Vec::new();
            for (f, k) in fs {
                let fe = f.as_expr();
                if fe.as_number().is_some_and(|n| n.is_positive()) {
                    factors.push((f, mul(&[Expr::int(k as i64), exp.clone()])));
                } else {
                    other.push(pow(&fe, &Expr::int(k as i64)));
                }
            }
            factors.push((FItem::E(mul(&other)), exp.clone()));
        }
    }
    if sqf {
        // multiply factors of equal multiplicity
        let mut ks: Vec<Expr> = Vec::new();
        for (_, k) in &factors {
            if !ks.contains(k) {
                ks.push(k.clone());
            }
        }
        let mut out = Vec::new();
        for k in ks {
            let group: Vec<&FItem> = factors.iter().filter(|(_, kk)| *kk == k).map(|(f, _)| f).collect();
            let item = if group.len() == 1 {
                group[0].clone()
            } else {
                let mut acc: Option<Poly> = None;
                let mut exprs = Vec::new();
                for g in &group {
                    match g {
                        FItem::P(p) => {
                            acc = Some(match acc {
                                None => p.clone(),
                                Some(a) => a.mul(p),
                            })
                        }
                        FItem::E(x) => exprs.push(x.clone()),
                    }
                }
                match (acc, exprs.is_empty()) {
                    (Some(p), true) => FItem::P(p),
                    (Some(p), false) => {
                        exprs.push(p.as_expr());
                        FItem::E(mul(&exprs))
                    }
                    (None, _) => FItem::E(mul(&exprs)),
                }
            };
            out.push((item, k));
        }
        factors = out;
    }
    // collect duplicates
    let mut out: Vec<(FItem, Expr)> = Vec::new();
    for (f, k) in factors {
        let fe = f.as_expr();
        if let Some(x) = out.iter_mut().find(|x| {
            if checked.is_some(){if let(FItem::P(a),FItem::P(b))=(&x.0,&f){if a.gens!=b.gens||a.ngens!=b.ngens||a.domain!=b.domain{return false;}}}
            x.0.as_expr()==fe
        }) {
            x.1 = add(&[x.1.clone(), k]);
        } else {
            out.push((f, k));
        }
    }
    out.retain(|(_, k)| !k.is_zero());
    Ok((coeff, out))
}

/// Natural expression factorization preserves the generators of each source
/// multiplicative component. Re-inferring generators on the *final* factors,
/// or keeping one global generator vector, changes SymPy's stable list order.
pub fn checked_natural_factor_list(e:&Expr,gens:&[Expr],squarefree:bool,prototype:&super::fraction::Polynomial,limits:super::factor_budget::Limits)->PResult<(Expr,Vec<(super::fraction::Polynomial,u32)>)> {
    super::factor_budget::run(limits,||{
        super::factor_preflight::check(e,limits)?;
        // together/as_numer_denom are legacy recursive expression algorithms.
        // Bound their unfolded work, recursion and prospective coefficient
        // products before entering them. The later polynomial engine retains
        // its separate exact arithmetic/storage bounds.
        guard_symbolic_factor_work(e,limits)?;
        let together=together(e);
        super::factor_preflight::check(&together,limits)?;
        guard_symbolic_factor_work(&together,limits)?;
        let(numerator,denominator)=as_numer_denom(&together);
        super::factor_preflight::check(&numerator,limits)?;
        super::factor_preflight::check(&denominator,limits)?;
        let(cp,fp)=symbolic_factor_list_with(&numerator,squarefree,None,gens,Some(limits))?;
        let(cq,fq)=symbolic_factor_list_with(&denominator,squarefree,None,gens,Some(limits))?;
        if !fq.is_empty(){return Err(PolyError::Failed("a polynomial expected: denominator factors require frac=True".into()));}
        let coefficient=div(&cp,&cq);
        let mut factors=Vec::new();
        for(f,k)in fp{
            let k=k.as_i64().filter(|&k|k>0&&k<=limits.degree as i64).ok_or_else(||PolyError::Unsupported("checked factor lists require positive integer multiplicities".into()))? as u32;
            let poly=match f{FItem::P(p)=>p,FItem::E(e)=>Poly::from_expr(&e,gens)?};
            factors.push((poly,k));
        }
        // Exactly the generic factor-list key. Equal keys retain the order
        // of components from canonical Mul.args and the inner factor engine.
        factors.sort_by(|(a,ka),(b,kb)|{
            let degree=a.rep.degree(0).cmp(&b.rep.degree(0));
            let generators=a.ngens.cmp(&b.ngens);
            let main=if squarefree{ka.cmp(kb).then(degree).then(generators)}else{degree.then(generators).then(ka.cmp(kb))};
            main.then_with(||a.domain_str().cmp(&b.domain_str())).then_with(||super::polyobj::rep_cmp(&a.rep,&b.rep,0))
        });
        let factors:Vec<_>=factors.into_iter().map(|(p,k)|(p.into(),k)).collect();
        super::checked_factors::verify_expression_factors(prototype,&coefficient,&factors,limits)?;
        Ok((coefficient,factors))
    })
}
fn guard_symbolic_factor_work(e:&Expr,limits:super::factor_budget::Limits)->PResult<()> {
    let mut stack=vec![(e.clone(),0usize,1u64)];let mut visits=0usize;let mut bits=0u64;
    let size=|n:&Number|->u64{match n{Number::Int(n)=>if n.unsigned_abs()<=1{0}else{64-n.unsigned_abs().leading_zeros() as u64},Number::BigInt(n)=>n.bits(),Number::Rational(n)=>n.numer().bits()+n.denom().bits(),Number::Float(_)=>1075,Number::MpFloat(n,_)=>n.man.bits().saturating_add(u64::try_from(n.exp.unsigned_abs()).unwrap_or(u64::MAX))}};
    while let Some((e,depth,weight))=stack.pop(){
        visits=visits.saturating_add(1);
        if visits>limits.terms||depth>256{return Err(PolyError::Unsupported("factorization resource limit: symbolic component traversal".into()));}
        let(mut coefficient,mut children)=(0u64,Vec::new());
        match e.kind(){
            Kind::Num(n)=>coefficient=size(n),
            Kind::Add(c,xs)=>{coefficient=size(c);for(e,c)in xs{coefficient=coefficient.saturating_add(size(c));children.push((e.clone(),weight));}},
            Kind::Mul(c,xs)=>{coefficient=size(c);for(b,p)in xs{let n=p.as_i64().map(i64::unsigned_abs).unwrap_or(1).max(1);children.push((b.clone(),weight.saturating_mul(n)));}},
            Kind::Pow(b,p)=>{let n=p.as_i64().map(i64::unsigned_abs).unwrap_or(1).max(1);children.push((b.clone(),weight.saturating_mul(n)));},
            Kind::Sym(_)|Kind::Const(_)|Kind::Bool(_)|Kind::Func(..)=>{},
            _=>{for child in e.args(){children.push((child,weight));}},
        }
        bits=bits.saturating_add(coefficient.saturating_mul(weight));
        if bits>limits.coefficient_bits{return Err(PolyError::Unsupported("factorization resource limit: symbolic coefficient products".into()));}
        super::factor_budget::charge(1);
        for(child,weight)in children{stack.push((child,depth+1,weight));}
    }
    Ok(())
}

/// SymPy's `_keep_coeff(coeff, factors)` (clear=True).
pub fn keep_coeff(coeff: &Expr, factors: &Expr) -> Expr {
    let Some(c) = coeff.as_number() else {
        if factors.is_number() {
            return keep_coeff(factors, coeff);
        }
        return mul(&[coeff.clone(), factors.clone()]);
    };
    if factors.is_one() {
        return coeff.clone();
    }
    if c.is_one() {
        return factors.clone();
    }
    if c.is_minus_one() {
        return neg(factors);
    }
    match factors.kind() {
        Kind::Add(..) => Expr::from_kind(Kind::Mul(c.clone(), vec![(factors.clone(), Expr::one())])),
        Kind::Mul(fc, facs) => {
            let nc = fc.mul(c);
            if nc.is_one() {
                if facs.len() == 1 {
                    return mul_factor(&facs[0].0, &facs[0].1);
                }
                return Expr::from_kind(Kind::Mul(Number::ONE, facs.clone()));
            }
            Expr::from_kind(Kind::Mul(nc, facs.clone()))
        }
        _ => {
            let m = mul(&[coeff.clone(), factors.clone()]);
            if m.is_number() && !factors.is_number() {
                return Expr::from_kind(Kind::Mul(c.clone(), vec![(factors.clone(), Expr::one())]));
            }
            m
        }
    }
}

fn factors_product(fs: &[(FItem, Expr)]) -> Expr {
    let parts: Vec<Expr> = fs.iter().map(|(f, k)| pow(&f.as_expr(), k)).collect();
    mul(&parts)
}

fn generic_factor(e: &Expr, sqf: bool, modulus: Option<u64>) -> PResult<Expr> {
    match e.kind() {
        Kind::Rel(..) | Kind::Piecewise(_) => {
            let mut err = None;
            let r = map_children(e, &mut |c| match generic_factor(c, sqf, modulus) {
                Ok(v) => v,
                Err(x) => {
                    err = Some(x);
                    c.clone()
                }
            });
            if let Some(x) = err {
                return Err(x);
            }
            return Ok(r);
        }
        Kind::Bool(_) => return Ok(e.clone()),
        _ => {}
    }
    let t = together(e);
    let (coeff, fs) = symbolic_factor_list(&t, sqf, modulus)?;
    Ok(keep_coeff(&coeff, &factors_product(&fs)))
}

/// `factor(expr)`.
pub fn factor(e: &Expr) -> PResult<Expr> {
    generic_factor(e, false, None)
}

/// `factor(expr, modulus=p)` (univariate).
pub fn factor_mod(e: &Expr, p: u64) -> PResult<Expr> {
    generic_factor(e, false, Some(p))
}

/// `sqf(expr, modulus=p)`.
pub fn sqf_mod(e: &Expr, p: u64) -> PResult<Expr> {
    generic_factor(e, true, Some(p))
}

/// `factor(expr, deep=True)`: factor the arguments of functions and the
/// sub-sums first.
pub fn factor_deep(e: &Expr) -> PResult<Expr> {
    if e.is_atom() {
        return Ok(e.clone());
    }
    let mut err = None;
    let inner = map_children(e, &mut |c| match factor_deep(c) {
        Ok(v) => v,
        Err(x) => {
            err = Some(x);
            c.clone()
        }
    });
    if let Some(x) = err {
        return Err(x);
    }
    match inner.kind() {
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => factor(&inner),
        _ => Ok(inner),
    }
}

/// `sqf(expr)`.
pub fn sqf(e: &Expr) -> PResult<Expr> {
    generic_factor(e, true, None)
}

fn item_poly(f: FItem) -> FItem {
    match f {
        FItem::P(p) => FItem::P(p),
        FItem::E(e) => match Poly::from_expr(&e, &[]) {
            Ok(p) => FItem::P(p),
            Err(_) => FItem::E(e),
        },
    }
}

fn generic_factor_list(e: &Expr, sqf: bool, modulus: Option<u64>) -> PResult<(Expr, Vec<(Expr, Expr)>)> {
    let t = together(e);
    let (n, d) = as_numer_denom(&t);
    let (cp, fp) = symbolic_factor_list(&n, sqf, modulus)?;
    let (cq, fq) = symbolic_factor_list(&d, sqf, modulus)?;
    if !fq.is_empty() {
        return Err(PolyError::Unsupported(format!("a polynomial expected, got {}", crate::printer::to_str(e))));
    }
    let fp: Vec<(FItem, Expr)> = fp.into_iter().map(|(f, k)| (item_poly(f), k)).collect();
    // sort the polynomial factors with SymPy's key, other factors last
    let mut polys: Vec<(Poly, u32, Expr)> = Vec::new();
    let mut others: Vec<(Expr, Expr)> = Vec::new();
    for (f, k) in fp {
        match (f, k.as_i64()) {
            (FItem::P(p), Some(kk)) if kk > 0 => polys.push((p, kk as u32, k)),
            (f, _) => others.push((f.as_expr(), k)),
        }
    }
    let mut keyed: Vec<(Poly, u32)> = polys.iter().map(|(p, k, _)| (p.clone(), *k)).collect();
    match modulus {
        None => sorted_factors(&mut keyed, sqf),
        Some(m) => keyed.sort_by(|(a, ka), (b, kb)| {
            let (ga, gb) = (a.to_gf(), b.to_gf());
            let main = if sqf { ka.cmp(kb).then(ga.len().cmp(&gb.len())) } else { ga.len().cmp(&gb.len()).then(ka.cmp(kb)) };
            let _ = m;
            main.then_with(|| ga.iter().rev().cmp(gb.iter().rev()))
        }),
    }
    let mut out: Vec<(Expr, Expr)> = keyed.into_iter().map(|(p, k)| (p.as_expr(), Expr::int(k as i64))).collect();
    out.extend(others);
    Ok((div(&cp, &cq), out))
}

/// `factor_list(expr)`: `(coeff, [(factor, multiplicity)])`.
pub fn factor_list(e: &Expr) -> PResult<(Expr, Vec<(Expr, Expr)>)> {
    generic_factor_list(e, false, None)
}

/// `factor_list(expr, modulus=p)`.
pub fn factor_list_mod(e: &Expr, p: u64) -> PResult<(Expr, Vec<(Expr, Expr)>)> {
    generic_factor_list(e, false, Some(p))
}

/// `sqf_list(expr, modulus=p)`.
pub fn sqf_list_mod(e: &Expr, p: u64) -> PResult<(Expr, Vec<(Expr, Expr)>)> {
    generic_factor_list(e, true, Some(p))
}

/// `sqf_list(expr)`.
pub fn sqf_list(e: &Expr) -> PResult<(Expr, Vec<(Expr, Expr)>)> {
    generic_factor_list(e, true, None)
}

/// `factor_list(expr, *gens)` / `Poly.factor_list()` with explicit
/// generators: `(coeff, [(poly, multiplicity)])`.
pub fn poly_factor_list(e: &Expr, gens: &[Expr]) -> PResult<(Expr, Vec<(Expr, u32)>)> {
    let p = Poly::from_expr(e, gens)?;
    let (c, fs) = p.factor_list();
    Ok((rat_expr(&c, p.is_float()), fs.into_iter().map(|(f, k)| (f.as_expr(), k)).collect()))
}

// ----------------------------------------------------------------------
// gcd / lcm
// ----------------------------------------------------------------------

fn numbers_of(a: &Expr, b: &Expr) -> Option<(Rat, Rat)> {
    let x = a.as_number().filter(|n| n.is_rational())?;
    let y = b.as_number().filter(|n| n.is_rational())?;
    Some((number_to_rat(x), number_to_rat(y)))
}

fn rat_gcd(a: &Rat, b: &Rat) -> Rat {
    Rat::new(a.n.gcd(&b.n), a.d.lcm(&b.d))
}

/// `gcd(f, g)`.
pub fn gcd(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    if let Some((x, y)) = numbers_of(a, b) {
        return Ok(Expr::num(rat_gcd(&x, &y).to_number()));
    }
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    Ok(ps[0].gcd(&ps[1]).as_expr())
}

/// `lcm(f, g)`.
pub fn lcm(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    if let Some((x, y)) = numbers_of(a, b) {
        if x.is_zero() || y.is_zero() {
            return Ok(Expr::zero());
        }
        let n = x.n.lcm(&y.n);
        let d = x.d.gcd(&y.d);
        return Ok(Expr::num(Rat::new(n, d).to_number()));
    }
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    Ok(ps[0].lcm(&ps[1]).as_expr())
}

fn mod_polys(a: &Expr, b: &Expr, gens: &[Expr], p: u64) -> PResult<(Poly, Poly)> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    if ps[0].gens.len() > 1 {
        return Err(PolyError::Unsupported("multivariate polynomials over finite fields".into()));
    }
    Ok((ps[0].set_modulus(p), ps[1].set_modulus(p)))
}

/// `gcd(f, g, modulus=p)`.
pub fn gcd_mod(a: &Expr, b: &Expr, gens: &[Expr], p: u64) -> PResult<Expr> {
    let (x, y) = mod_polys(a, b, gens, p)?;
    Ok(x.gcd(&y).as_expr())
}

/// `lcm(f, g, modulus=p)`.
pub fn lcm_mod(a: &Expr, b: &Expr, gens: &[Expr], p: u64) -> PResult<Expr> {
    let (x, y) = mod_polys(a, b, gens, p)?;
    Ok(x.lcm(&y).as_expr())
}

/// `div(f, g, modulus=p)`.
pub fn div_mod(a: &Expr, b: &Expr, gens: &[Expr], p: u64) -> PResult<(Expr, Expr)> {
    let (x, y) = mod_polys(a, b, gens, p)?;
    let (q, r) = x.div(&y)?;
    Ok((q.as_expr(), r.as_expr()))
}

/// `cofactors(f, g)`: `(h, f/h, g/h)`.
pub fn cofactors(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<(Expr, Expr, Expr)> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    let (h, x, y) = ps[0].cofactors(&ps[1]);
    Ok((h.as_expr(), x.as_expr(), y.as_expr()))
}

// ----------------------------------------------------------------------
// extended Euclid, terms_gcd, composition
// ----------------------------------------------------------------------

fn uni_pair(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<(Vec<Rat>, Vec<Rat>, Vec<Expr>)> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    if ps[0].gens.len() != 1 {
        return Err(PolyError::Unsupported("univariate polynomials expected".into()));
    }
    Ok((dense::from_mpoly(&ps[0].rep, 0), dense::from_mpoly(&ps[1].rep, 0), ps[0].gens.clone()))
}

fn dense_expr(p: &[Rat], gens: &[Expr]) -> Expr {
    super::convert::to_expr(&dense::to_mpoly(p, 0, 1), gens)
}

/// `gcdex(f, g)`: `(s, t, h)` with `s*f + t*g = h = gcd(f, g)` (monic).
pub fn gcdex(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<(Expr, Expr, Expr)> {
    let (f, g, gens) = uni_pair(a, b, gens)?;
    let (h, s, t) = dense::gcdex_field(&f, &g);
    Ok((dense_expr(&s, &gens), dense_expr(&t, &gens), dense_expr(&h, &gens)))
}

/// `half_gcdex(f, g)`: `(s, h)` with `s*f = h (mod g)`.
pub fn half_gcdex(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<(Expr, Expr)> {
    let (s, _, h) = gcdex(a, b, gens)?;
    Ok((s, h))
}

/// `invert(f, g)`: inverse of `f` modulo `g`.
pub fn invert(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    if let (Some(x), Some(m)) = (a.as_number(), b.as_number()) {
        if let (Some(x), Some(m)) = (Int::from_number(x), Int::from_number(m)) {
            if m.is_zero() { return Err(PolyError::Unsupported("division by zero".into())); }
            let (g, s, _) = x.gcdex(&m);
            if !g.is_one() {
                return Err(PolyError::NotInvertible("inverse does not exist".into()));
            }
            return Ok(Expr::num(s.modulo(&m.abs()).to_number()));
        }
    }
    let (f, g, gens) = uni_pair(a, b, gens)?;
    if g.is_empty() { return Err(PolyError::Unsupported("polynomial division by zero".into())); }
    let (h, s, _) = dense::gcdex_field(&f, &g);
    if h.len() != 1 {
        return Err(PolyError::NotInvertible("zero divisor".into()));
    }
    Ok(dense_expr(&dense::rem(&s, &g), &gens))
}

/// `terms_gcd(f)`: pull out the common monomial and the content.
pub fn terms_gcd(e: &Expr) -> PResult<Expr> {
    let p = match Poly::from_expr(e, &[]) {
        Ok(p) => p,
        Err(_) => return Ok(e.clone()),
    };
    let m = p.rep.min_exps();
    let rest = p.rep.div_monomial(&m);
    let (den, z) = rest.clear_denoms();
    let c = z.content();
    let prim = if c.is_zero() { z.clone() } else { z.div_coeff(&c).unwrap() };
    let coeff = Rat::new(c, den);
    let mut term = vec![];
    for (g, &k) in p.gens.iter().zip(&m) {
        if k > 0 {
            term.push(pow(g, &Expr::int(k as i64)));
        }
    }
    let term = mul(&term);
    if term.is_one() && coeff.is_one() {
        return Ok(e.clone());
    }
    let body = mul(&[term, super::convert::to_expr_int(&prim, &p.gens)]);
    Ok(keep_coeff(&Expr::num(coeff.to_number()), &body))
}

/// `compose(f, g)`: `f(g)` (in the first generator of `f`).
pub fn compose(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    Ok(ps[0].compose(&ps[1]).as_expr())
}

/// `decompose(f)`: functional decomposition `[g1, ..., gn]` with
/// `f = g1(g2(...gn))`.
pub fn decompose(e: &Expr, gens: &[Expr]) -> PResult<Vec<Expr>> {
    let p = Poly::from_expr(e, gens)?;
    Ok(p.decompose()?.into_iter().map(|q| q.as_expr()).collect())
}

/// `sqf_part(f)`.
pub fn sqf_part(e: &Expr, gens: &[Expr]) -> PResult<Expr> {
    let p = Poly::from_expr(e, gens)?;
    Ok(p.sqf_part().as_expr())
}

// ----------------------------------------------------------------------
// cancel
// ----------------------------------------------------------------------

fn contains_bad(e: &Expr) -> bool {
    let mut bad = false;
    crate::traverse::preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Rel(..) | Kind::Piecewise(_) | Kind::Bool(_) | Kind::Derivative(..)) {
            bad = true;
        }
        if let Kind::Const(c) = s.kind() {
            if c.is_infinite() || *c == Constant::NaN {
                bad = true;
            }
        }
    });
    bad
}

/// Cancel common factors of `p/q` (both polynomials over the same ring).
/// Returns `(P, Q)` with `Q` having a positive leading coefficient.
pub fn cancel_polys(p: &MPoly<Rat>, q: &MPoly<Rat>) -> (MPoly<Rat>, MPoly<Rat>) {
    let n = p.nvars;
    if p.is_zero() {
        return (p.clone(), MPoly::one(n));
    }
    let (cq, f) = p.clear_denoms();
    let (cp, g) = q.clear_denoms();
    let (_, pp, qq) = super::gcd::gcd_mz(&f, &g);
    let gg = cp.gcd(&cq);
    let cp = cp.div_exact_unchecked(&gg);
    let cq = cq.div_exact_unchecked(&gg);
    let mut pp = pp.scale(&cp).to_rat();
    let mut qq = qq.scale(&cq).to_rat();
    if qq.lc().is_negative() {
        pp = pp.neg();
        qq = qq.neg();
    }
    (pp, qq)
}

/// `cancel(expr)`.
pub fn cancel(e: &Expr) -> PResult<Expr> {
    match e.kind() {
        Kind::Num(_) | Kind::Rel(..) | Kind::Bool(_) => return Ok(e.clone()),
        _ => {}
    }
    if contains_bad(e) {
        return Ok(map_children(e, &mut |c| cancel(c).unwrap_or_else(|_| c.clone())));
    }
    let (p, q) = as_numer_denom(e);
    let c = match super::convert::from_exprs(&[p.clone(), q.clone()]) {
        Ok(c) => c,
        Err(_) => return Ok(e.clone()),
    };
    if c.ngens == 0 {
        return Ok(expand(e));
    }
    let (pp, qq) = cancel_polys(&c.polys[0], &c.polys[1]);
    let pe = super::convert::to_expr_f(&pp, &c.gens, c.float);
    let qe = super::convert::to_expr_f(&qq, &c.gens, c.float);
    Ok(div(&pe, &qe))
}

// ----------------------------------------------------------------------
// division
// ----------------------------------------------------------------------

/// `div(f, g, *gens)`: `(quotient, remainder)`.
pub fn div_expr(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<(Expr, Expr)> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    let (q, r) = ps[0].div(&ps[1])?;
    Ok((q.as_expr(), r.as_expr()))
}

/// `reduced(f, divisors, *gens)`: ordered multivariate division with
/// remainder in lexicographic monomial order. Automatic generator detection
/// uses the existing deterministic SymPy-compatible generator ordering.
pub fn reduced(a: &Expr, divisors: &[Expr], gens: &[Expr]) -> PResult<(Vec<Expr>, Expr)> {
    use super::convert::{from_exprs, from_exprs_gens, to_expr_f};
    use super::ratfunc::{set_nvars, RF};
    use std::collections::BTreeMap;

    let mut exprs = vec![a.clone()];
    exprs.extend_from_slice(divisors);
    let converted = if gens.is_empty() { from_exprs(&exprs)? } else { from_exprs_gens(&exprs, gens)? };
    if converted.polys.iter().skip(1).any(MPoly::is_zero) {
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    // SymPy's sparse division returns no quotient entries for zero input.
    if converted.polys[0].is_zero() {
        return Ok((Vec::new(), Expr::zero()));
    }
    if converted.ngens == converted.gens.len() {
        let (q, r) = converted.polys[0].reduced(&converted.polys[1..]).unwrap();
        let to_expr = |p: &MPoly<Rat>| to_expr_f(p, &converted.gens, converted.float);
        return Ok((q.iter().map(to_expr).collect(), to_expr(&r)));
    }

    // Explicit generators leave parameters in the coefficient domain. Work
    // over its fraction field, so e.g. reduced(x, [a*x], x) gives 1/a.
    let n = converted.ngens;
    let domain_gens = &converted.gens[n..];
    let old_nvars = set_nvars(domain_gens.len());
    let to_rf_poly = |p: &MPoly<Rat>| {
        let mut grouped: BTreeMap<Vec<u32>, Vec<(Vec<u32>, Rat)>> = BTreeMap::new();
        for (e, c) in p.iter() {
            grouped.entry(e[..n].to_vec()).or_default().push((e[n..].to_vec(), c.clone()));
        }
        MPoly::from_terms(n, grouped.into_iter().map(|(e, terms)| {
            (e, RF::from_rat_poly(&MPoly::from_terms(domain_gens.len(), terms)))
        }).collect())
    };
    let polys: Vec<MPoly<RF>> = converted.polys.iter().map(to_rf_poly).collect();
    let (q, r) = polys[0].reduced(&polys[1..]).unwrap();
    let to_expr = |p: &MPoly<RF>| {
        add(&p.iter().map(|(e, c)| {
            let mut factors = vec![rf_expr(c, domain_gens)];
            for (g, &k) in converted.gens[..n].iter().zip(e) {
                if k != 0 {
                    factors.push(pow(g, &Expr::int(k as i64)));
                }
            }
            mul(&factors)
        }).collect::<Vec<_>>())
    };
    let result = (q.iter().map(to_expr).collect(), to_expr(&r));
    set_nvars(old_nvars);
    Ok(result)
}

/// `rem(f, g, *gens)`.
pub fn rem_expr(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    Ok(div_expr(a, b, gens)?.1)
}

/// `quo(f, g, *gens)`.
pub fn quo_expr(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    Ok(div_expr(a, b, gens)?.0)
}

/// `pdiv(f, g, *gens)`: pseudo-division `(q, r)` with
/// `lc(g)**(deg f - deg g + 1) * f = q*g + r`.
pub fn pdiv_expr(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<(Expr, Expr)> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    let (q, r) = ps[0].pdiv(&ps[1])?;
    Ok((q.as_expr(), r.as_expr()))
}

/// `prem(f, g, *gens)`: pseudo-remainder.
pub fn prem_expr(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    Ok(pdiv_expr(a, b, gens)?.1)
}

/// `pquo(f, g, *gens)`: pseudo-quotient.
pub fn pquo_expr(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    Ok(pdiv_expr(a, b, gens)?.0)
}

/// `pexquo(f, g, *gens)`: exact pseudo-quotient.
pub fn pexquo_expr(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    let (q, r) = pdiv_expr(a, b, gens)?;
    if !r.is_zero() {
        return Err(PolyError::Unsupported(format!(
            "{} does not divide {}",
            crate::printer::to_str(b),
            crate::printer::to_str(a)
        )));
    }
    Ok(q)
}

/// `exquo(f, g)`.
pub fn exquo(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    Ok(ps[0].exquo(&ps[1])?.as_expr())
}

/// `degree(f, gen)`; `None` stands for `-oo`.
pub fn degree(e: &Expr, gen: Option<&Expr>) -> PResult<Option<u32>> {
    if let Some(n) = e.as_number() {
        return Ok(if n.is_zero() { None } else { Some(0) });
    }
    match gen {
        Some(g) => {
            let p = Poly::from_expr(e, std::slice::from_ref(g))?;
            Ok(p.degree(0))
        }
        None => {
            if crate::traverse::free_symbols(e).len() > 1 {
                return Err(PolyError::Unsupported(
                    "A symbolic generator of interest is required for a multivariate expression".into(),
                ));
            }
            let p = Poly::from_expr(e, &[])?;
            Ok(p.degree(0))
        }
    }
}

// ----------------------------------------------------------------------
// resultant / discriminant
// ----------------------------------------------------------------------

fn int_with_denoms(p: &MPoly<Rat>) -> (Int, MPoly<Int>) {
    p.clear_denoms()
}

/// `resultant(f, g, *gens)` (with respect to the first generator).
pub fn resultant(a: &Expr, b: &Expr, gens: &[Expr]) -> PResult<Expr> {
    let ps = Poly::from_exprs(&[a.clone(), b.clone()], gens)?;
    let (da, za) = int_with_denoms(&ps[0].rep);
    let (db, zb) = int_with_denoms(&ps[1].rep);
    let r = resultant_mz(&za, &zb, 0);
    // res(f/da, g/db) = res(f, g) / (da^deg g * db^deg f)
    let m = ps[1].rep.degree(0).max(0) as u32;
    let n = ps[0].rep.degree(0).max(0) as u32;
    let scale = Rat::new(Int::one(), da.pow(m).mul(&db.pow(n)));
    let r = r.to_rat().scale(&scale);
    Ok(super::convert::to_expr_f(&r, &ps[0].gens, ps[0].is_float()))
}

/// `discriminant(f, *gens)`.
pub fn discriminant(f: &Expr, gens: &[Expr]) -> PResult<Expr> {
    let p = Poly::from_expr(f, gens)?;
    let (d, z) = int_with_denoms(&p.rep);
    let r = discriminant_mz(&z, 0);
    // disc(f/d) = disc(f) / d^(2n-2)
    let n = p.rep.degree(0).max(1) as u32;
    let r = r.to_rat().scale(&Rat::new(Int::one(), d.pow(2 * n - 2)));
    Ok(super::convert::to_expr_f(&r, &p.gens, p.is_float()))
}

// ----------------------------------------------------------------------
// apart
// ----------------------------------------------------------------------

/// `apart(f, x)`: univariate partial fraction decomposition over `QQ`.
pub fn apart(f: &Expr, x: Option<&Expr>) -> PResult<Expr> {
    if f.is_atom() {
        return Ok(f.clone());
    }
    let (p, q) = as_numer_denom(f);
    let gens: Vec<Expr> = x.into_iter().cloned().collect();
    let ps = Poly::from_exprs(&[p, q], &gens)?;
    let (pp, qp) = (&ps[0], &ps[1]);
    if pp.ngens == 1 && pp.has_domain_gens() {
        return apart_param(pp, qp);
    }
    if pp.gens.len() > 1 || pp.has_domain_gens() {
        let fc = cancel(f)?;
        if fc != *f {
            return apart(&fc, x);
        }
        return Err(PolyError::Unsupported("multivariate partial fraction decomposition".into()));
    }
    let gen = pp.gens[0].clone();
    // P.cancel(Q)
    let (cp, cq, mut pz, mut qz) = cancel_int(&pp.rep, &qp.rep);
    let common = Rat::new(cp, cq);
    if qz.is_empty() {
        return Err(PolyError::Unsupported("division by zero".into()));
    }
    // polynomial part
    let (poly, rem) = dense::divrem(&dense::to_rat(&pz), &dense::to_rat(&qz));
    let (a, remz) = dense::clear_denoms(&rem);
    pz = remz;
    // P/Q = rem/Q ; rat_clear_denoms: (rem*a) / (Q*a)... keep Q integral
    qz = dense::scale(&qz, &a);
    let poly_e = super::convert::to_expr(&dense::to_mpoly(&poly, 0, 1), std::slice::from_ref(&gen));
    let mut terms: Vec<Expr> = Vec::new();
    if !pz.is_empty() {
        if qz.len() <= 2 {
            let pe = super::convert::to_expr(&dense::to_mpoly(&dense::to_rat(&pz), 0, 1), std::slice::from_ref(&gen));
            let qe = super::convert::to_expr(&dense::to_mpoly(&dense::to_rat(&qz), 0, 1), std::slice::from_ref(&gen));
            terms.push(factor(&div(&pe, &qe))?);
        } else {
            for t in partial_fractions(&pz, &qz) {
                let (h, fct, k) = t;
                let he = super::convert::to_expr(&dense::to_mpoly(&h, 0, 1), std::slice::from_ref(&gen));
                let fe = super::convert::to_expr(&dense::to_mpoly(&dense::to_rat(&fct), 0, 1), std::slice::from_ref(&gen));
                let term = div(&he, &pow(&fe, &Expr::int(k as i64)));
                terms.push(factor(&term)?);
            }
        }
    }
    let mut all = vec![poly_e];
    all.extend(terms);
    Ok(mul(&[Expr::num(common.to_number()), add(&all)]))
}

fn rf_expr(r: &super::ratfunc::RF, gens: &[Expr]) -> Expr {
    let n = super::convert::to_expr_int(&r.n, gens);
    if r.d.is_one() {
        return n;
    }
    div(&n, &super::convert::to_expr_int(&r.d, gens))
}

fn rf_dense_expr(p: &[super::ratfunc::RF], gens: &[Expr]) -> Expr {
    let mut terms = Vec::new();
    for (i, c) in p.iter().enumerate() {
        if c.n.is_zero() {
            continue;
        }
        terms.push(mul(&[rf_expr(c, gens), pow(&gens[0], &Expr::int(i as i64))]));
    }
    add(&terms)
}

fn rf_dense(p: &MPoly<Int>) -> Vec<super::ratfunc::RF> {
    use super::ratfunc::RF;
    let n = p.nvars;
    let d = p.degree(0).max(0) as usize;
    let mut out: Vec<RF> = (0..=d).map(|_| RF::from_poly(MPoly::zero(n))).collect();
    if p.is_zero() {
        return vec![];
    }
    for (k, c) in p.coeffs_in(0) {
        out[k as usize] = RF::from_poly(c);
    }
    out
}

/// `apart` with symbolic parameters: partial fractions over `QQ(params)`.
fn apart_param(pp: &Poly, qp: &Poly) -> PResult<Expr> {
    use super::ratfunc::{set_nvars, RF};
    let n = pp.gens.len();
    let gens = pp.gens.clone();
    let old = set_nvars(n);
    let res = (|| -> PResult<Expr> {
        // P.cancel(Q) over ZZ[params] (dmp_cancel with include=False)
        let (cq0, f) = pp.rep.clear_denoms();
        let (cp0, g) = qp.rep.clear_denoms();
        let (_, mut p, mut q) = super::gcd::gcd_mz(&f, &g);
        let gg = cp0.gcd(&cq0);
        let mut cp = cp0.div_exact_unchecked(&gg);
        let cq = cq0.div_exact_unchecked(&gg);
        let pn = p.lc().is_negative();
        let qn = q.lc().is_negative();
        if pn && qn {
            p = p.neg();
            q = q.neg();
        } else if pn {
            cp = cp.neg();
            p = p.neg();
        } else if qn {
            cp = cp.neg();
            q = q.neg();
        }
        let common = Rat::new(cp, cq);
        if q.degree(0) <= 0 {
            return Err(PolyError::Unsupported("apart: denominator free of the variable".into()));
        }
        let pd = rf_dense(&p);
        let qd = rf_dense(&q);
        let (poly, rem) = if pd.is_empty() { (vec![], vec![]) } else { dense::divrem(&pd, &qd) };
        let mut terms: Vec<Expr> = vec![rf_dense_expr(&poly, &gens)];
        if !rem.is_empty() {
            if q.degree(0) <= 1 {
                let t = div(&rf_dense_expr(&rem, &gens), &rf_dense_expr(&qd, &gens));
                terms.push(factor(&t)?);
            } else {
                let (_, fs) = super::factor_multi::factor_list_mz(&q);
                for (fct, k) in fs {
                    if fct.degree(0) <= 0 {
                        continue;
                    }
                    let fd = rf_dense(&fct);
                    let mut a = vec![<RF as super::int::Coeff>::one()];
                    for _ in 0..k {
                        a = dense::mul(&a, &fd);
                    }
                    let (b, _) = dense::divrem(&qd, &a);
                    let (_, s, _) = dense::gcdex_field(&dense::rem(&b, &a), &a);
                    let mut num = dense::rem(&dense::mul(&rem, &s), &a);
                    let fe = super::convert::to_expr_int(&fct, &gens);
                    let mut i = 0u32;
                    while !num.is_empty() && i < k {
                        let (qt, rr) = dense::divrem(&num, &fd);
                        if !rr.is_empty() {
                            let he = rf_dense_expr(&rr, &gens);
                            let t = div(&he, &pow(&fe, &Expr::int((k - i) as i64)));
                            terms.push(factor(&t)?);
                        }
                        num = qt;
                        i += 1;
                    }
                }
            }
        }
        Ok(mul(&[Expr::num(common.to_number()), add(&terms)]))
    })();
    set_nvars(old);
    res
}

/// `dmp_cancel(include=False)` for univariate polynomials over `QQ`:
/// `(cp, cq, p, q)` with positive leading coefficients.
fn cancel_int(p: &MPoly<Rat>, q: &MPoly<Rat>) -> (Int, Int, Vec<Int>, Vec<Int>) {
    let (cq0, f) = p.clear_denoms();
    let (cp0, g) = q.clear_denoms();
    let (_, pp, qq) = super::gcd::gcd_mz(&f, &g);
    let gg = cp0.gcd(&cq0);
    let mut cp = cp0.div_exact_unchecked(&gg);
    let cq = cq0.div_exact_unchecked(&gg);
    let mut pd = dense::from_mpoly(&pp, 0);
    let mut qd = dense::from_mpoly(&qq, 0);
    let pn = pd.last().is_some_and(|c| c.is_negative());
    let qn = qd.last().is_some_and(|c| c.is_negative());
    if pn && qn {
        pd = dense::neg(&pd);
        qd = dense::neg(&qd);
    } else if pn {
        cp = cp.neg();
        pd = dense::neg(&pd);
    } else if qn {
        cp = cp.neg();
        qd = dense::neg(&qd);
    }
    (cp, cq, pd, qd)
}

/// Partial fractions of `p/q` (`deg p < deg q`, integer coefficients):
/// `[(numerator, factor, power)]` with `p/q = Σ numerator / factor^power`.
fn partial_fractions(p: &[Int], q: &[Int]) -> Vec<(Vec<Rat>, Vec<Int>, u32)> {
    let (_, fs) = super::factor_uni::factor_zx(q);
    // SymPy's order: factor_list order of Q
    let mut fs: Vec<(Vec<Int>, u32)> = fs;
    fs.sort_by(|(a, ka), (b, kb)| {
        a.len().cmp(&b.len()).then(ka.cmp(kb)).then_with(|| {
            for (x, y) in a.iter().rev().zip(b.iter().rev()) {
                let c = x.cmp(y);
                if c != std::cmp::Ordering::Equal {
                    return c;
                }
            }
            std::cmp::Ordering::Equal
        })
    });
    let pq = dense::to_rat(p);
    let qq = dense::to_rat(q);
    let mut out = Vec::new();
    for (f, k) in fs {
        let fq = dense::to_rat(&f);
        let mut a = vec![Rat::one()];
        for _ in 0..k {
            a = dense::mul(&a, &fq);
        }
        let (b, r) = dense::divrem(&qq, &a);
        debug_assert!(r.is_empty());
        // n = p * b^-1 mod a
        let (_, s, _) = dense::gcdex_field(&dense::rem(&b, &a), &a);
        // gcdex returns a monic gcd 1 with s*b + t*a = 1
        let mut nnum = dense::rem(&dense::mul(&pq, &s), &a);
        // f-adic expansion: nnum = Σ r_i f^i
        let mut i = 0u32;
        while !nnum.is_empty() && i < k {
            let (qt, rr) = dense::divrem(&nnum, &fq);
            if !rr.is_empty() {
                out.push((rr, f.clone(), k - i));
            }
            nnum = qt;
            i += 1;
        }
    }
    out
}

// ----------------------------------------------------------------------
// collect
// ----------------------------------------------------------------------

/// `collect(expr, syms)` with default expression-valued options.
pub fn collect(e:&Expr,syms:&[Expr])->Expr{
    match super::collect::run::<std::convert::Infallible>(e,syms,super::collect::Options::default(),true,None){
        Ok(super::collect::Collection::Expression(value))=>value,
        _=>e.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    use crate::printer::to_str;

    fn s(e: &str) -> Expr {
        parse_str(e).unwrap()
    }

    #[test]
    fn checked_natural_factor_components_preserve_generator_provenance_and_order(){
        let limits=super::super::factor_budget::Limits::default();
        for (source,expected) in [
            ("x^4*y-x^2*y^3",vec!["y","x","x-y","x+y"]),
            ("x*y+x",vec!["x","y+1"]),
            ("x^2*y-y",vec!["x-1","y","x+1"]),
            ("x^2*y+x*y^2",vec!["x","y","x+y"]),
            ("-716*x*z",vec!["x","z"]),
            ("432*w^2*x^2*y^6",vec!["w","x","y"]),
        ]{
            let expression=s(source);
            let prototype=Poly::from_expr(&expression,&[]).unwrap().into();
            let(_,factors)=checked_natural_factor_list(&expression,&[],false,&prototype,limits).unwrap();
            assert_eq!(factors.iter().map(|(p,_)|p.as_expr()).collect::<Vec<_>>(),expected.iter().map(|e|s(e)).collect::<Vec<_>>(),"{source}");
        }
        let expression=s("x*y+x");
        let gens=vec![Expr::symbol("x"),Expr::symbol("y")];
        let prototype=Poly::from_expr(&expression,&gens).unwrap().into();
        let(_,factors)=checked_natural_factor_list(&expression,&gens,false,&prototype,limits).unwrap();
        assert_eq!(factors[0].0.as_expr(),s("y+1"));
        assert!(factors.iter().all(|(p,_)|p.true_gens()==gens));
        let mut dag=Expr::symbol("x");
        for _ in 0..40{dag=Expr::from_unevaluated_kind(Kind::Add(Number::ZERO,vec![(dag.clone(),Number::ONE),(dag,Number::ONE)]));}
        assert!(guard_symbolic_factor_work(&dag,limits).is_err());
    }

    #[test]
    fn factor_basic() {
        let cases = [
            ("x**2 - 1", "(x - 1)*(x + 1)"),
            ("2*x**2 + 4*x + 2", "2*(x + 1)**2"),
            ("2*x + 2", "2*(x + 1)"),
            ("-2*x - 2", "-2*(x + 1)"),
            ("-x**2 + 1", "-(x - 1)*(x + 1)"),
            ("x/2 + 1", "(x + 2)/2"),
            ("x**3*y - x*y**3", "x*y*(x - y)*(x + y)"),
            ("(x**2 - 1)/(x**2 + 2*x + 1)", "(x - 1)/(x + 1)"),
            ("sin(x)**2 - 1", "(sin(x) - 1)*(sin(x) + 1)"),
            ("x**2 + 1", "x**2 + 1"),
            ("6*x**2*y + 4*x*y**2", "2*x*y*(3*x + 2*y)"),
        ];
        for (a, b) in cases {
            assert_eq!(to_str(&factor(&s(a)).unwrap()), b, "factor({a})");
        }
    }

    #[test]
    fn cancel_basic() {
        let cases = [
            ("(x**2 - 1)/(x - 1)", "x + 1"),
            ("(2*x + 2)/(4*x**2 - 4)", "1/(2*x - 2)"),
            ("x/2 + 1/(3*x)", "(3*x**2 + 2)/(6*x)"),
            ("(x*y - y)/(x**2 - 1)", "y/(x + 1)"),
            ("1/(-x + 1)", "-1/(x - 1)"),
        ];
        for (a, b) in cases {
            assert_eq!(to_str(&cancel(&s(a)).unwrap()), b, "cancel({a})");
        }
    }

    #[test]
    fn apart_basic() {
        let cases = [
            ("1/(x**2 - 1)", "-1/(2*(x + 1)) + 1/(2*(x - 1))"),
            ("(x**2 + 1)/(x**3 - x)", "1/(x + 1) + 1/(x - 1) - 1/x"),
            ("1/(x*(x + 1)**2)", "-1/(x + 1) - 1/(x + 1)**2 + 1/x"),
            ("x**3/(x**2 - 1)", "x + 1/(2*(x + 1)) + 1/(2*(x - 1))"),
        ];
        for (a, b) in cases {
            assert_eq!(to_str(&apart(&s(a), None).unwrap()), b, "apart({a})");
        }
    }

    #[test]
    fn reduced_multivariate_lex() {
        let x = s("x");
        let y = s("y");
        let cases = [
            ("x**2 + y**2", vec!["x + y"], vec![], vec!["x - y"], "2*y**2"),
            ("x**2 + y**2", vec!["x + y"], vec![y.clone(), x.clone()], vec!["y - x"], "2*x**2"),
            ("x*y**2 + x*y + x", vec!["x*y - 1", "y**2 - 1"], vec![], vec!["y + 1", "0"], "x + y + 1"),
            ("2*x + 1", vec!["3*x"], vec![], vec!["2/3"], "1"),
            ("x", vec!["a*x"], vec![x.clone()], vec!["1/a"], "0"),
            ("x + y", vec![], vec![], vec![], "x + y"),
        ];
        for (expr, divisors, gens, expected_q, expected_r) in cases {
            let ds: Vec<Expr> = divisors.iter().map(|d| s(d)).collect();
            let (q, r) = reduced(&s(expr), &ds, &gens).unwrap();
            assert_eq!(q, expected_q.iter().map(|v| s(v)).collect::<Vec<_>>(), "quotient {expr}");
            assert_eq!(r, s(expected_r), "remainder {expr}");
            let reconstructed = add(&q.iter().zip(&ds).map(|(q, d)| mul(&[q.clone(), d.clone()])).chain(std::iter::once(r)).collect::<Vec<_>>());
            assert_eq!(expand(&reconstructed), expand(&s(expr)));
        }
        assert!(reduced(&x, &[Expr::zero()], &[]).is_err());
    }

    #[test]
    fn misc() {
        assert_eq!(to_str(&gcd(&s("x**2 - 1"), &s("x**2 + 2*x + 1"), &[]).unwrap()), "x + 1");
        assert_eq!(to_str(&gcd(&s("2*x**2 - 2"), &s("4*x + 4"), &[]).unwrap()), "2*x + 2");
        assert_eq!(to_str(&lcm(&s("x**2 - 1"), &s("x + 1"), &[]).unwrap()), "x**2 - 1");
        assert_eq!(to_str(&resultant(&s("x**2 + y"), &s("x*y - 1"), &[s("x")]).unwrap()), "y**3 + 1");
        assert_eq!(to_str(&discriminant(&s("x**2 + b*x + c"), &[s("x")]).unwrap()), "b**2 - 4*c");
        let (q, r) = div_expr(&s("x**2 + 1"), &s("2*x"), &[]).unwrap();
        assert_eq!((to_str(&q), to_str(&r)), ("x/2".to_string(), "1".to_string()));
        assert_eq!(to_str(&collect(&s("a*x**2 + b*x**2 + a*x - b*x + c"), &[s("x")])), "c + x**2*(a + b) + x*(a - b)");
        assert_eq!(to_str(&collect(&s("x**2*y + x*y + x + y"), &[s("x"), s("y")])), "x**2*y + x*(y + 1) + y");
        let (c, fs) = factor_list(&s("x**4*y - x**2*y**3")).unwrap();
        let fs: Vec<String> = fs.iter().map(|(f, k)| format!("({}, {})", to_str(f), to_str(k))).collect();
        assert_eq!(to_str(&c), "1");
        assert_eq!(fs.join(", "), "(y, 1), (x, 2), (x - y, 1), (x + y, 1)");
    }
}
