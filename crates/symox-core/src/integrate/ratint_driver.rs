//! Bounded QQ/ZZ driver for the public rational integration pipeline.
//!
//! Keep cancellation units, polynomial/rational parts, signed LRT groups and
//! their primitive logarithm representatives. Any incomplete subcalculation
//! refuses the whole route; the public caller then runs its full algorithm.
use crate::expr::{Expr,Kind};
use crate::ops::{add,mul2,div};
use crate::poly::{api,dense,fraction::Polynomial,
    polyobj::{Domain,Poly},Rat};

pub enum Input {
    Expression(Expr),
    Fraction(Expr,Expr),
}

fn rational_coefficient(n:&crate::number::Number)->bool {
    Rat::from_number(n)
        .is_some_and(|q|q.n.bits().max(q.d.bits())<=4096)
}

// Integer-polynomial bounds for a representation N/D. Each coefficient has
// absolute value < 2^bits. Product convolution has at most min(degA,degB)+1
// terms; sums add one bit. These deliberately overestimate cancellation.
#[derive(Clone,Copy)]
struct Shape { n:u32,d:u32,nb:u64,db:u64 }
impl Shape {
    fn checked(n:u32,d:u32,nb:u64,db:u64)->Option<Self>{
        (n<=32 && d<=32 && nb<=4096 && db<=4096).then_some(Self{n,d,nb,db})
    }
    fn number(n:&crate::number::Number)->Option<Self>{
        let q=Rat::from_number(n)?;Self::checked(0,0,q.n.bits(),q.d.bits())
    }
    fn convolution(a:u32,b:u32)->u64{
        let terms=u64::from(a.min(b))+1;
        u64::from(64-(terms-1).leading_zeros())
    }
    fn mul(self,other:Self)->Option<Self>{
        Self::checked(self.n.checked_add(other.n)?,self.d.checked_add(other.d)?,
            self.nb+other.nb+Self::convolution(self.n,other.n),
            self.db+other.db+Self::convolution(self.d,other.d))
    }
    fn add(self,other:Self)->Option<Self>{
        let left=self.n.checked_add(other.d)?;let right=other.n.checked_add(self.d)?;
        let bits_left=self.nb+other.db+Self::convolution(self.n,other.d);
        let bits_right=other.nb+self.db+Self::convolution(other.n,self.d);
        Self::checked(left.max(right),self.d.checked_add(other.d)?,bits_left.max(bits_right)+1,
            self.db+other.db+Self::convolution(self.d,other.d))
    }
    fn power(self,exponent:i64)->Option<Self>{
        let power=u32::try_from(exponent.unsigned_abs()).ok().filter(|&n|n<=32)?;
        if power==0{return Self::checked(0,0,1,1);}
        let mut value=Self::checked(self.n.checked_mul(power)?,self.d.checked_mul(power)?,
            u64::from(power)*(self.nb+Self::convolution(self.n,self.n)),
            u64::from(power)*(self.db+Self::convolution(self.d,self.d)))?;
        if exponent<0 {std::mem::swap(&mut value.n,&mut value.d);std::mem::swap(&mut value.nb,&mut value.db);}
        Some(value)
    }
}

// Bound effective degrees and coefficient growth BEFORE as_numer_denom or
// expansion. A bounded exponent at every node alone does not bound nesting.
fn rational_shape(e:&Expr,x:&Expr,budget:&mut usize,depth:u32)->Option<Shape> {
    if depth>64 || *budget==0{return None;}*budget-=1;
    match e.kind() {
        Kind::Num(n)=>Shape::number(n),
        Kind::Sym(_) if e==x=>Shape::checked(1,0,1,1),
        Kind::Add(c,args)=>args.iter().try_fold(Shape::number(c)?,|shape,(a,c)|{
            shape.add(Shape::number(c)?.mul(rational_shape(a,x,budget,depth+1)?)?)
        }),
        Kind::Mul(c,factors)=>factors.iter().try_fold(Shape::number(c)?,|shape,(b,e)|{
            shape.mul(rational_shape(b,x,budget,depth+1)?.power(e.as_i64()?)?)
        }),
        Kind::Pow(b,e)=>rational_shape(b,x,budget,depth+1)?.power(e.as_i64()?),
        _=>None,
    }
}

fn polynomial_degree(e:&Expr,x:&Expr,budget:&mut usize,depth:u32)->Option<u32> {
    if depth>64 || *budget==0{return None;}*budget-=1;
    let power=|b:&Expr,e:&Expr,budget:&mut usize|{
        let n=u32::try_from(e.as_i64()?).ok()?;
        let d=polynomial_degree(b,x,budget,depth+1)?.checked_mul(n)?;
        (d<=32).then_some(d)
    };
    match e.kind() {
        Kind::Num(n) if rational_coefficient(n)=>Some(0),
        Kind::Sym(_) if e==x=>Some(1),
        Kind::Add(c,args) if rational_coefficient(c)=>args.iter().try_fold(0u32,|d,(a,c)|{
            if !rational_coefficient(c){return None;}Some(d.max(polynomial_degree(a,x,budget,depth+1)?))
        }),
        Kind::Mul(c,factors) if rational_coefficient(c)=>{
            factors.iter().try_fold(0u32,|d,(b,e)|{
                let d=d.checked_add(power(b,e,budget)?)?;(d<=32).then_some(d)
            })
        },
        Kind::Pow(b,e)=>power(b,e,budget),
        _=>None,
    }
}

fn dense_poly(p:&[Rat],x:&Expr,domain:Domain)->Poly {
    Poly::new(dense::to_mpoly(p,0,1),vec![x.clone()],1,domain)
}

fn primitive(h:Polynomial)->Option<Polynomial> {
    if let Some(field)=h.field {
        return Some(field.reference_primitive().1.into());
    }
    if h.algebraic.is_some() || h.expression.is_some(){return None;}
    if h.base.has_domain_gens(){Some(h.base.ground_primitive().1.into())}
    else{Some(h.base.primitive().1.into())}
}

pub fn integrate<E>(input:&Input,x:&Expr,symbol:&str,check:&mut impl FnMut()->Result<(),E>)
    -> Result<Option<(Expr,bool)>,E>
{
    check()?;
    let (p,q)=match input {
        Input::Expression(f)=>{
            if rational_shape(f,x,&mut 4096,0).is_none(){return Ok(None);}
            crate::simptools::as_numer_denom(f)
        },
        Input::Fraction(p,q)=>{
            if rational_shape(p,x,&mut 4096,0).is_none() || rational_shape(q,x,&mut 4096,0).is_none(){return Ok(None);}
            (p.clone(),q.clone())
        },
    };
    if rational_shape(&p,x,&mut 4096,0).is_none() || rational_shape(&q,x,&mut 4096,0).is_none()
        || polynomial_degree(&p,x,&mut 4096,0).is_none() || polynomial_degree(&q,x,&mut 4096,0).is_none(){return Ok(None);}
    macro_rules! step {($e:expr)=>{{check()?;let result=$e;check()?;
        let Ok(value)=result else{return Ok(None)};value}}}
    let mut p=step!(Poly::from_expr(&p,std::slice::from_ref(x)));
    let mut q=step!(Poly::from_expr(&q,std::slice::from_ref(x)));
    if p.gens!=[x.clone()] || q.gens!=[x.clone()] || p.ngens!=1 || q.ngens!=1 || q.rep.is_zero(){return Ok(None);}
    p.domain=Domain::QQ;q.domain=Domain::QQ;
    // Reuse the reference dmp_cancel unit convention, already used by apart.
    check()?;
    let (cp,cq,p,q)=api::cancel_int(&p.rep,&q.rep);
    check()?;
    let coefficient=Expr::num(Rat::new(cp,cq).to_number());
    let (p,q)=(dense_poly(&dense::to_rat(&p),x,Domain::QQ),dense_poly(&dense::to_rat(&q),x,Domain::QQ));
    let (polynomial,p)=step!(p.div(&q));
    let mut result=polynomial.integrate(0).as_expr();
    if p.rep.is_zero(){check()?;return Ok(Some((mul2(&coefficient,&result),false)));}
    let Some(reduction)=super::rational_reduction::reduce(&dense::from_mpoly(&p.rep,0),&dense::from_mpoly(&q.rep,0),check)?
        else{return Ok(None)};
    let fraction=|a:&[Rat],b:&[Rat]|api::cancel(&div(&dense_poly(a,x,Domain::QQ).as_expr(),&dense_poly(b,x,Domain::QQ).as_expr()));
    let rational=step!(fraction(&reduction.rational_numerator,&reduction.rational_denominator));
    let logarithmic=step!(fraction(&reduction.logarithmic_numerator,&reduction.logarithmic_denominator));
    if rational_shape(&logarithmic,x,&mut 4096,0).is_none(){return Ok(None);}
    let (p,q)=crate::simptools::as_numer_denom(&logarithmic);
    let p=step!(Poly::from_expr(&p,std::slice::from_ref(x)));
    let q=step!(Poly::from_expr(&q,std::slice::from_ref(x)));
    let (polynomial,remainder)=step!(p.div(&q));
    result=add(&[result,rational,polynomial.integrate(0).as_expr()]);
    let logarithmic=!remainder.rep.is_zero();
    if logarithmic {
        check()?;
        let t=Expr::dummy(symbol,crate::symbol::Assumptions::NONE);
        let Some(groups)=super::logarithmic_reduction::reduce(&remainder,&q,&t,check)? else{return Ok(None)};
        let mut terms=Vec::with_capacity(groups.len());
        for (h,q) in groups {
            check()?;
            // No partial RootSum/atan sum: a missing certificate restarts the
            // complete public driver, including its general RootSum choice.
            if q.field.is_some() || q.algebraic.is_some() || q.expression.is_some()
                || !(1..=16).contains(&q.rep.degree(0)){return Ok(None);}
            let Some(h)=primitive(h) else{return Ok(None)};
            let Some(h)=super::log_to_real::rational_input(&h,x,&t,check)? else{return Ok(None)};
            let Some(value)=super::log_to_real::convert(&h,&dense::from_mpoly(&q.rep,0),x,&t,check)? else{return Ok(None)};
            terms.push(value);
        }
        result=add(&[result,add(&terms)]);
    }
    check()?;
    Ok(Some((mul2(&coefficient,&result),logarithmic)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn signed_units_repeated_poles_and_real_residue_representatives() {
        let x=Expr::symbol("x");
        for (input,expected) in [("1/(x^2+1)","atan(x)"),
            ("-3/(2*x+2)","-3*log(x+1)/2"),
            ("1/x^2","-1/x"),
            ("x^2+3/2","x^3/3+3*x/2")] {
            let result=integrate(&Input::Expression(parse_str(input).unwrap()),&x,"t",&mut ||Ok::<_,()>(())).unwrap().unwrap();
            assert_eq!(result.0,parse_str(expected).unwrap());
        }
    }
    #[test]
    fn unsupported_roots_resources_and_interruptions_never_publish_partial_answers() {
        let x=Expr::symbol("x");let f=parse_str("1/(x^3+x+1)").unwrap();
        assert!(integrate(&Input::Expression(f),&x,"t",&mut ||Ok::<_,()>(())).unwrap().is_none());
        let f=parse_str("x^1000000+1").unwrap();
        assert!(integrate(&Input::Expression(f),&x,"t",&mut ||Ok::<_,()>(())).unwrap().is_none());
        let mut calls=0;
        let result=integrate(&Input::Expression(parse_str("1/(x^2+1)").unwrap()),&x,"t",&mut ||{calls+=1;if calls==12{Err("stop")}else{Ok(())}});
        assert!(matches!(result,Err("stop")));
    }
    #[test]
    fn nested_degrees_and_numeric_powers_are_bounded_before_numerator_extraction() {
        let x=Expr::symbol("x");
        let one=Expr::int(1).as_number().unwrap().clone();
        for mut base in [parse_str("x+1").unwrap(),Expr::int(2)] {
            // Raw nodes keep hidden numeric powers unevaluated. Constructing
            // the witness itself must not execute the dangerous calculation.
            for _ in 0..8 {
                let power=Expr::from_kind(Kind::Pow(base,Expr::int(32)));
                base=Expr::from_kind(Kind::Add(one.clone(),vec![(power,one.clone())]));
            }
            assert!(rational_shape(&base,&x,&mut 4096,0).is_none());
            for input in [Input::Expression(base.clone()),Input::Fraction(base,Expr::one())] {
                assert!(integrate(&input,&x,"t",&mut ||Ok::<_,()>(())).unwrap().is_none());
            }
        }
    }
}
