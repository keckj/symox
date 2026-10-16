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
