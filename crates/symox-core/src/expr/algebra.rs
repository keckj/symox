//! Size/eligibility metadata for bounded reuse of ordinary exact arithmetic.
//! This metadata never participates in expression identity or serialization.
use super::{Expr,Kind};
use crate::number::Number;

const LIMIT: u16 = 1024;

fn number(value:&Number)->Option<u16>{
    match value{
        Number::Int(_)=>Some(1),
        Number::BigInt(n) if n.bits()<=256=>Some(1+(n.bits()/64) as u16),
        Number::Rational(n) if n.numer().bits()<=256&&n.denom().bits()<=256=>
            Some(1+((n.numer().bits()+n.denom().bits())/64) as u16),
        _=>None,
    }
}

fn child(value:&Expr)->Option<u16>{
    let weight=value.exact_algebra_weight();
    (weight!=0).then_some(weight)
}

fn exponent(value:&Expr)->Option<u16>{
    let n=value.as_number()?;
    if !n.is_integer()||n.is_negative(){return None;}
    number(n)
}

fn add(total:&mut u16,weight:u16)->Option<()>{
    *total=total.checked_add(weight)?;
    (*total<=LIMIT).then_some(())
}

pub(super) fn weight(kind:&Kind,canonical:bool)->u16{
    if !canonical{return 0;}
    fn compute(kind:&Kind)->Option<u16>{
        let mut total=1;
        match kind{
            Kind::Num(n)=>return number(n),
            Kind::Sym(s) if !s.is_dummy()&&!s.is_wild()&&s.name().len()<=256
                &&s.facts().get("commutative")==Some(true)
                &&s.facts().get("finite")!=Some(false)=>
                return Some(1+(s.name().len()/8) as u16),
            Kind::Add(c,terms)=>{
                add(&mut total,number(c)?)?;
                for (term,c) in terms{
                    add(&mut total,child(term)?)?;
                    add(&mut total,number(c)?)?;
                }
            }
            Kind::Mul(c,factors)=>{
                add(&mut total,number(c)?)?;
                for (base,power) in factors{
                    add(&mut total,child(base)?)?;
                    add(&mut total,exponent(power)?)?;
                }
            }
            Kind::Pow(base,power)=>{
                add(&mut total,child(base)?)?;
                add(&mut total,exponent(power)?)?;
            }
            _=>return None,
        }
        Some(total)
    }
    compute(kind).unwrap_or(0)
}
