//! Exact closed-QQ splitting factorization for D(t)=tau(t).
use crate::poly::{dense,Rat};
#[derive(Debug,PartialEq)]
pub struct Split {
    pub normal:Vec<Rat>,pub special:Vec<Rat>,
    pub reuse_input:bool,pub reuse_one:bool,pub special_integer:bool,
}
fn gcd<E>(a:&[Rat],b:&[Rat],check:&mut impl FnMut()->Result<(),E>)->Result<Vec<Rat>,E>{
    let(mut a,mut b)=(a.to_vec(),b.to_vec());
    while !b.is_empty(){check()?;let r=dense::rem(&a,&b);a=b;b=if r.is_empty(){r}else{dense::monic(&r)};}
    check()?;Ok(dense::monic(&a))
}
fn exact(a:&[Rat],b:&[Rat])->Option<Vec<Rat>>{
    if b.is_empty(){return None;}let(q,r)=dense::divrem(a,b);r.is_empty().then_some(q)
}
/// Unsupported zero derivations decline the whole candidate. Coefficients
/// are exact QQ; optional booleans describe the original constant branch.
pub fn split<E>(input:&[Rat],tau:&[Rat],lower_variables:bool,
    check:&mut impl FnMut()->Result<(),E>)->Result<Option<Split>,E>{
    check()?;
    if tau.is_empty()||tau.last().is_some_and(Rat::is_zero)
        ||input.last().is_some_and(Rat::is_zero){return Ok(None);}
    let mut p=input.to_vec();let mut s=vec![Rat::one()];let mut iterations=0;
    if p.is_empty(){return Ok(Some(Split{normal:p,special:s,reuse_input:true,reuse_one:true,special_integer:false}));}
    loop{
        check()?;
        if p.len()==1{
            // Poly.as_poly in different lower generators re-infers the
            // constant domain. Its ZZ gcd absorbs abs(integer content).
            let integer_unit=lower_variables&&p[0].d.is_one();
            if integer_unit{
                let unit=Rat::from(p[0].n.abs());
                p[0]=p[0].div(&unit);s=dense::scale(&s,&unit);
            }
            return Ok(Some(Split{normal:p,special:s,reuse_input:false,reuse_one:false,
                special_integer:lower_variables&&iterations==0}));
        }
        let ordinary=dense::derivative(&p);let dp=dense::mul(&ordinary,tau);
        let h=gcd(&p,&dp,check)?;let g=gcd(&p,&ordinary,check)?;
        let Some(factor)=exact(&h,&g)else{return Ok(None)};
        if factor.len()<=1{
            return Ok(Some(Split{normal:p,special:s,reuse_input:iterations==0,
                reuse_one:iterations==0,special_integer:false}));
        }
        let Some(quotient)=exact(&p,&factor)else{return Ok(None)};
        p=quotient;s=dense::mul(&s,&factor);iterations+=1;
    }
}
#[cfg(test)]mod tests{
    use super::*;
    fn q(a:&[i64])->Vec<Rat>{a.iter().map(|&v|Rat::from_i64(v)).collect()}
    fn checked(p:&[Rat],tau:&[Rat],lower:bool)->Split{
        let s=split(p,tau,lower,&mut||Ok::<_,()>(())).unwrap().unwrap();
        assert_eq!(dense::mul(&s.normal,&s.special),p);
        let ds=dense::mul(&dense::derivative(&s.special),tau);
        assert!(dense::rem(&ds,&s.special).is_empty());
        if !s.normal.is_empty(){
            let sq=dense::divrem(&s.normal,&dense::gcd_field(&s.normal,&dense::derivative(&s.normal))).0;
            let dn=dense::mul(&dense::derivative(&sq),tau);
            assert_eq!(dense::gcd_field(&sq,&dn),q(&[1]));
        }s
    }
    #[test]fn repeated_normal_and_special_factors_have_exact_certificates(){
        for normal_power in 0..=5{for special_power in 0..=5{
            let normal=(0..normal_power).fold(q(&[1]),|a,_|dense::mul(&a,&q(&[1,1])));
            let special=(0..special_power).fold(q(&[1]),|a,_|dense::mul(&a,&q(&[1,0,1])));
            let p=dense::scale(&dense::mul(&normal,&special),&Rat::from_i64(-6));
            for lower in [false,true]{checked(&p,&q(&[1,0,1]),lower);}
        }}
    }
    #[test]fn constant_units_domains_and_passthrough_aliases_follow_the_old_branches(){
        let s=checked(&q(&[6]),&q(&[1]),true);assert_eq!(s.normal,q(&[1]));assert_eq!(s.special,q(&[6]));assert!(s.special_integer);
        let s=checked(&q(&[-6]),&q(&[1]),true);assert_eq!(s.normal,q(&[-1]));assert_eq!(s.special,q(&[6]));
        let s=checked(&q(&[6]),&q(&[1]),false);assert_eq!(s.normal,q(&[6]));assert_eq!(s.special,q(&[1]));assert!(!s.special_integer);
        let s=checked(&q(&[1,1]),&q(&[1]),true);assert!(s.reuse_input&&s.reuse_one);
        let s=checked(&[],&q(&[1]),true);assert!(s.reuse_input&&s.reuse_one);
    }
    #[test]fn zero_derivation_refuses_and_interruptions_propagate(){
        assert!(split(&q(&[1,1]),&[],false,&mut||Ok::<_,()>(())).unwrap().is_none());
        assert_eq!(split(&q(&[1,1]),&q(&[1]),false,&mut||Err("interrupt")),Err("interrupt"));
    }
}
