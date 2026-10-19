//! Cheap exact rational powers and bounded quadratic facts for solver checks.
//! No numerical evaluation and no Groebner basis computation participates.
use num_traits::{One, Signed, ToPrimitive, Zero};

use crate::expr::{Constant, Expr, Kind};
use crate::poly::algebraic::{AN, NumberField};
use crate::poly::int::{Coeff, Rat};

#[derive(Clone, Copy, Default, Debug)]
pub(crate) struct Facts {
    pub integer: Option<bool>,
    pub rational: Option<bool>,
    pub real: Option<bool>,
    pub zero: Option<bool>,
    pub sign: Option<i8>,
    pub imaginary: Option<bool>,
}

/// Principal cube-root conjugate pair: if |z|^2=c^3, c>0, then for
/// u=z^(1/3), c/u=conjugate(u). t=u+c/u is the largest real root of
/// t^3-3*c*t-2*Re(z), in [sqrt(c),2*sqrt(c)]. That cubic is increasing
/// above sqrt(c), so its exact value at a rational threshold certifies
/// the sign of every rational affine expression in t.
pub(super) fn cardano_sign(e:&Expr)->Option<i8>{
    let Kind::Add(constant,terms)=e.kind() else{return None};
    if terms.len()!=2||!constant.is_rational()||terms.iter().any(|(_,c)|!c.is_rational()){return None;}
    let power=|e:&Expr|match e.kind(){
        Kind::Pow(base,exponent) if *exponent==Expr::rational(1,3)=>Some((base.clone(),true)),
        Kind::Pow(base,exponent) if *exponent==Expr::rational(-1,3)=>Some((base.clone(),false)),
        _=>None,
    };
    let (base,positive)=power(&terms[0].0)?;
    let (other,other_positive)=power(&terms[1].0)?;
    if base!=other||positive==other_positive{return None;}
    let (a,b)=if positive{(&terms[0].1,&terms[1].1)}else{(&terms[1].1,&terms[0].1)};
    if a.is_zero(){return None;}
    let c=a.inv()?.mul(b);
    if !c.is_positive(){return None;}
    let (real,imaginary)=super::exact_queries(||crate::complex_parts::as_real_imag(&base,false,None)).ok()??;
    let norm=crate::ops::add2(&real.powi(2),&imaginary.powi(2));
    let expected=Expr::num(c.clone()).powi(3);
    if !norm.is_rational_number()||norm!=expected{return None;}
    let real=real.as_number().filter(|value|value.is_rational())?;
    let threshold=constant.neg().div(a)?;
    let square=threshold.mul(&threshold);
    if !threshold.is_positive()||square.cmp_value(&c)==Some(std::cmp::Ordering::Less){
        return Some(if a.is_positive(){1}else{-1});
    }
    let polynomial=square.mul(&threshold).sub(&c.mul(&threshold).mul(&crate::number::Number::Int(3)))
        .sub(&real.mul(&crate::number::Number::Int(2)));
    let sign=if polynomial.is_positive(){1}else if polynomial.is_negative(){-1}else{0};
    Some(if a.is_positive(){-sign}else{sign})
}

/// Rational powers/products with reciprocal exponent vectors. No expression
/// expansion is needed to look for a matching pair; integer additive factors
/// may have their rational content normalized.
fn reciprocal_factors(term:&Expr)->Option<(Vec<(Expr,Expr)>,crate::number::Number)>{
    let (mut factors,mut scale)=match term.kind(){
        Kind::Pow(base,exponent)=>(vec![(base.clone(),exponent.clone())],crate::number::Number::ONE),
        Kind::Mul(coefficient,factors) if coefficient.is_rational()=>(factors.clone(),coefficient.clone()),
        _=>return None,
    };
    if !factors.iter().all(|(_,exponent)|exponent.is_rational_number()){return None;}
    for (base,exponent) in &mut factors{
        // Canonical input may distribute a rational into an additive base.
        // Pull it back out only for small INTEGER powers, where this is
        // branch-independent and cannot construct an enormous new integer.
        if !matches!(exponent.as_i64(),Some(-2|-1|1|2)){continue;}
        let Kind::Add(constant,terms)=base.kind() else{continue;};
        if constant.is_zero()||!constant.is_rational()||terms.iter().any(|(_,c)|!c.is_rational()){continue;}
        let powered=crate::ops::pow(&Expr::num(constant.clone()),exponent);
        scale=scale.mul(powered.as_number()?);
        *base=Expr::from_kind(Kind::Add(crate::number::Number::ONE,
            terms.iter().map(|(term,c)|(term.clone(),c.div(constant).unwrap())).collect()));
    }
    // Native product display order depends on exponents; reciprocal vectors
    // need a common base order for equality/hash matching.
    factors.sort_by(|a,b|a.0.cmp(&b.0));
    Some((factors,scale))
}

/// If every factor has an exactly rational squared modulus after its rational
/// power, so does their product. Principal powers obey |z^r|²=(|z|²)^r for
/// real r and nonzero finite z, independently of which side of a cut z lies.
fn rational_modulus_squared(factors:&[(Expr,Expr)])->Option<crate::number::Number>{
    let mut product=crate::number::Number::ONE;
    for (base,exponent) in factors{
        let (real,imaginary)=super::exact_queries(||crate::complex_parts::as_real_imag(base,false,None)).ok()??;
        let norm=crate::ops::add2(&real.powi(2),&imaginary.powi(2));
        let norm=norm.as_number().filter(|n|n.is_rational()&&n.is_positive())?;
        let powered=crate::ops::pow(&Expr::num(norm.clone()),exponent);
        let powered=powered.as_number().filter(|n|n.is_rational()&&n.is_positive())?;
        product=product.mul(powered);
    }
    Some(product)
}

/// Identify a*T+b/T pairs with rational a,b and b/a=|T|²>0. They are
/// a*(T+conjugate(T)), hence real. Unlike the narrower Cardano sign theorem,
/// this also covers root-of-unity multiples and arbitrary rational powers.
/// Hash matching avoids an all-pairs term scan. Each product signature needs
/// its factors sorted; other addends still require their own real proofs.
pub(super) fn paired_conjugate_terms(terms:&[(Expr,crate::number::Number)])->Vec<bool>{
    let mut signatures:rustc_hash::FxHashMap<Vec<(Expr,Expr)>,(usize,crate::number::Number)>=Default::default();
    let mut paired=vec![false;terms.len()];
    for (index,(term,coefficient)) in terms.iter().enumerate(){
        if !coefficient.is_rational()||coefficient.is_zero(){continue;}
        let Some((factors,scale))=reciprocal_factors(term) else{continue;};
        let coefficient=coefficient.mul(&scale);
        let inverse:Vec<_>=factors.iter().map(|(base,exponent)|
            (base.clone(),Expr::num(exponent.as_number().unwrap().neg()))).collect();
        if let Some((other,other_coefficient))=signatures.get(&inverse){
            let ratio=coefficient.inv().unwrap().mul(other_coefficient);
            if ratio.is_positive()&&rational_modulus_squared(&factors)==Some(ratio){
                paired[index]=true;paired[*other]=true;
            }
        }
        signatures.insert(factors,(index,coefficient));
    }
    paired
}

pub(super) fn rational_power(base: &Expr, exponent: &Expr) -> Option<Facts> {
    let (n,d)=base.as_number()?.numer_denom()?;
    let (p,q)=exponent.as_number()?.numer_denom()?;
    let mut facts=Facts::default();
    if p.is_zero() { return Some(Facts {integer:Some(true),rational:Some(true),real:Some(true),zero:Some(false),sign:Some(1),imaginary:Some(false)}); }
    if n.is_zero() {
        return Some(if p.is_positive() { Facts {integer:Some(true),rational:Some(true),real:Some(true),zero:Some(true),sign:Some(0),imaginary:Some(false)} }
            else { Facts {integer:Some(false),rational:Some(false),real:Some(false),zero:Some(false),imaginary:Some(false),..facts} });
    }
    facts.zero=Some(false);
    if n.is_negative() && !q.is_one() {
        facts.real=Some(false);facts.rational=Some(false);facts.integer=Some(false);
        facts.imaginary=Some(q==2.into());
        return Some(facts);
    }
    facts.real=Some(true);facts.imaginary=Some(false);
    facts.sign=Some(if n.is_negative() && num_integer::Integer::is_odd(&p) {-1} else {1});
    let a=n.abs();
    let perfect=if q.is_one() || (a.is_one() && d.is_one()) { Some(true) }
        else if q>num_bigint::BigInt::from(a.bits().max(d.bits())) { Some(false) }
        else if a.bits()+d.bits()>16_384 { None }
        else {
            let q=q.to_u32()?;
            Some(crate::number::exact_root(&a,q).is_some() && crate::number::exact_root(&d,q).is_some())
        };
    facts.rational=perfect;
    facts.integer=perfect.map(|rational| rational && if p.is_positive() { d.is_one() } else { a.is_one() });
    Some(facts)
}

/// Restrict the quadratic fallback before AN::from_expr can expand squares.
/// 128 visited nodes, <=32 estimated terms, degree <=16, small integer powers.
fn small_rational(e:&Expr)->bool {
    e.as_number().and_then(|n|n.numer_denom()).is_some_and(|(n,d)|n.bits()+d.bits()<=4096)
}
fn quadratic_shape(e:&Expr, visits:&mut usize, generator:&mut Option<Expr>) -> Option<(usize,usize)> {
    if *visits>=128 {return None;} *visits+=1;
    match e.kind() {
        Kind::Num(_) if small_rational(e) => Some((1,0)),
        Kind::Const(Constant::I) => {if generator.is_none(){*generator=Some(e.clone());}Some((1,1))},
        Kind::Add(..) | Kind::Mul(..) => {
            let is_add=e.is_add();
            let had_generator=generator.is_some();
            let mut terms=if is_add {0usize}else{1usize};let mut degree=0;
            for a in e.args() {
                let (t,d)=quadratic_shape(&a,visits,generator)?;
                terms=if is_add {terms.checked_add(t)?}else{terms.checked_mul(t)?};
                degree=if is_add {degree.max(d)}else{degree+d};
                if terms>32 || degree>16 {return None;}
            }
            // Only square a prospective generator after every child passed
            // the exponent, coefficient and expansion-size guards.
            if !is_add && !had_generator && e.args().iter().all(|a| match a.kind() {
                Kind::Num(_)=>small_rational(a),Kind::Const(Constant::I)=>true,
                Kind::Pow(b,x)=>small_rational(b) && x.as_number().and_then(|n|n.numer_denom()).is_some_and(|(p,q)|q==2.into() && p.abs()<=16.into()),
                _=>false,
            }) && e.powi(2).is_rational_number() { *generator=Some(e.clone()); }
            Some((terms,degree))
        },
        Kind::Pow(base,exponent) => {
            if let Some(k)=exponent.as_i64() {
                if k.unsigned_abs()>8 {return None;}
                let (t,d)=quadratic_shape(base,visits,generator)?;
                let terms=t.checked_pow(k.unsigned_abs() as u32)?;let degree=d*k.unsigned_abs() as usize;
                if terms>32 || degree>16 {None}else{Some((terms,degree))}
            } else {
                let (p,q)=exponent.as_number()?.numer_denom()?;
                if q!=2.into() || p.abs()>16.into() || !small_rational(base) {return None;}
                if generator.is_none(){*generator=Some(crate::ops::sqrt(base));}
                Some((1,1))
            }
        },
        _=>None,
    }
}

pub(crate) fn quadratic_facts(e:&Expr)->Option<Facts> {
    if !crate::traverse::free_symbols(e).is_empty(){return None;}
    let mut generator=None;
    quadratic_shape(e,&mut 0,&mut generator)?;
    let selected=NumberField::quadratic(&generator?).ok()?;
    // NumberField retains its selected embedding, including a negative or
    // imaginary scale. The sign formula below requires the principal sqrt(d)
    // embedding; AN::from_expr then proves the corresponding coefficients.
    let principal=crate::ops::sqrt(&Expr::num(selected.radicand.to_number()));
    let field=NumberField::quadratic(&principal).ok()?;
    let number=AN::from_expr(e,field.clone()).ok()?;
    let a=number.coefficients.first().cloned().unwrap_or_else(Rat::zero);
    let b=number.coefficients.get(1).cloned().unwrap_or_else(Rat::zero);
    let rational=b.is_zero();
    let real=rational || !field.radicand.is_negative();
    let mut facts=Facts {integer:Some(rational && a.d.is_one()),rational:Some(rational),real:Some(real),zero:Some(number.is_zero()),imaginary:Some(!real && a.is_zero()),sign:None};
    if real {
        let sign=|r:&Rat| if r.is_zero(){0}else if r.is_negative(){-1}else{1};
        let (sa,sb)=(sign(&a),sign(&b));
        facts.sign=Some(if sb==0 {sa}else if sa==0 || sa==sb {sb}else {
            match a.mul(&a).cmp(&b.mul(&b).mul(&field.radicand)) {
                std::cmp::Ordering::Equal=>0,std::cmp::Ordering::Greater=>sa,std::cmp::Ordering::Less=>sb,
            }
        });
    }
    Some(facts)
}

#[cfg(test)]
mod cardano_tests{
    use super::*;
    #[test]
    fn exact_norm_and_monotone_cubic_certify_real_signs(){
        let z=crate::parser::parse_str("(3+4*I)/5").unwrap();
        let u=crate::ops::pow(&z,&Expr::rational(1,3));
        let t=crate::ops::add2(&u,&crate::ops::div(&Expr::one(),&u));
        for (shift,scale,expected) in [(0,1,1),(-1,1,1),(-2,1,-1),(1,-1,-1),(2,-1,1)]{
            let value=crate::ops::add2(&Expr::int(shift),&crate::ops::mul2(&Expr::int(scale),&t));
            assert_eq!(cardano_sign(&value),Some(expected));
            assert_eq!(crate::assumptions::exact_queries(||crate::assumptions::is_real(&value)),Some(true));
        }
        let perturbed=crate::ops::add2(&z,&crate::ops::mul2(&Expr::i(),&Expr::rational(1,1000)));
        let u=crate::ops::pow(&perturbed,&Expr::rational(1,3));
        let unproved=crate::ops::add2(&u,&crate::ops::div(&Expr::one(),&u));
        assert_eq!(cardano_sign(&unproved),None);
    }
    #[test]
    fn reciprocal_product_pairs_are_real_only_with_the_exact_modulus(){
        let omega="(-1/2+sqrt(3)*I/2)";
        let z="(-1/2+sqrt(4035)*I/18)";
        for source in [
            format!("{omega}^(4/3)+{omega}^(-4/3)"),
            format!("{omega}^(7/3)+{omega}^(-7/3)+sqrt(2)"),
            format!("{omega}*{z}^(1/3)+7/(3*{omega}*{z}^(1/3))"),
            format!("{omega}^2*{z}^(1/3)+7/(3*{omega}^2*{z}^(1/3))"),
            format!("-{omega}*{z}^(1/3)-7/(3*{omega}*{z}^(1/3))"),
        ]{
            let value=crate::parser::parse_str(&source).unwrap();
            assert_eq!(crate::assumptions::exact_queries(||crate::assumptions::is_real(&value)),Some(true),"{source}");
        }
        for source in [
            format!("{omega}*{z}^(1/3)+1/({omega}*{z}^(1/3))"),
            format!("{omega}*{z}^(1/3)-7/(3*{omega}*{z}^(1/3))"),
            format!("({omega}+I/1000)^(4/3)+({omega}+I/1000)^(-4/3)"),
        ]{
            let value=crate::parser::parse_str(&source).unwrap();
            assert_ne!(crate::assumptions::exact_queries(||crate::assumptions::is_real(&value)),Some(true),"{source}");
        }
    }
}
