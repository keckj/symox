//! Exact real logarithm conversion for rational linear/quadratic residues.
//!
//! Enumerate the same real roots or the positive-imaginary conjugate as the
//! general resultant algorithm. Keep each logarithm and each Euclidean atan
//! term: a derivative certificate alone would not fix their branch constants.
use crate::expr::Expr;
use crate::functions::{atan, log};
use crate::ops::{add, mul2, sqrt};
use crate::poly::{algebraic::{AlgebraicScope, NumberField, AN}, convert, dense,
    int::{Coeff, Field, Rat}, mpoly::MPoly, fraction::Polynomial,
    domains::ExactDomain, polyobj::Domain};

trait RealCoefficient: Field {
    fn expression(&self) -> Expr;
    fn rational(q: &Rat) -> Self;
}
impl RealCoefficient for Rat {
    fn expression(&self) -> Expr { Expr::num(self.to_number()) }
    fn rational(q: &Rat) -> Self { q.clone() }
}
impl RealCoefficient for AN {
    fn expression(&self) -> Expr { self.as_expr() }
    fn rational(q: &Rat) -> Self { AN::rational(q.clone(),AN::zero().field) }
}

fn polynomial<C: RealCoefficient>(p: &[C], x: &Expr) -> Expr {
    add(&p.iter().enumerate().filter(|(_,c)|!c.is_zero()).map(|(i,c)|
        mul2(&c.expression(),&x.powi(i as i64))).collect::<Vec<_>>())
}

fn conjugate<C: RealCoefficient,E>(h: &MPoly<Rat>,u: C,v: C,x: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Expr>,E>
{
    let mut powers=vec![(C::one(),C::zero())];
    for j in 0..h.degree(1).max(0) as usize {
        check()?;
        let (a,b)=&powers[j];
        powers.push((a.mul(&u).sub(&b.mul(&v)),a.mul(&v).add(&b.mul(&u))));
    }
    let length=h.degree(0).max(0) as usize+1;
    let (mut a,mut b)=(vec![C::zero();length],vec![C::zero();length]);
    for (e,c) in h.iter() {
        check()?;
        let coefficient=C::rational(c);
        a[e[0] as usize]=a[e[0] as usize].add(&powers[e[1] as usize].0.mul(&coefficient));
        b[e[0] as usize]=b[e[0] as usize].add(&powers[e[1] as usize].1.mul(&coefficient));
    }
    dense::trim(&mut a);dense::trim(&mut b);
    // Preserve the old zero-polynomial/error paths through the full adapter.
    let Some(arguments)=super::log_to_atan::arguments_field(&a,&b,check)? else{return Ok(None)};
    check()?;
    let norm=dense::add(&dense::mul(&a,&a),&dense::mul(&b,&b));
    let terms=arguments.into_iter().map(|p|mul2(&Expr::int(2),&atan(&polynomial(&p,x)))).collect::<Vec<_>>();
    let result=add(&[mul2(&u.expression(),&log(&polynomial(&norm,x))),mul2(&v.expression(),&add(&terms))]);
    check()?;
    Ok(Some(result))
}

/// Read the exact stored polynomial in (x,t), shared by the public binding
/// and the rational integration driver. General fraction/expression domains
/// never pass through their legacy metadata view.
pub fn rational_input<E>(h: &Polynomial,x: &Expr,t: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<MPoly<Rat>>,E>
{
    if h.algebraic.is_some() || h.expression.is_some() || h.base.modulus.is_some()
        || h.base.ngens!=1 || h.base.true_gens()!=[x.clone()] {return Ok(None);}
    let polynomial=if let Some(field)=&h.field {
        if !matches!(&field.coefficient_domain,ExactDomain::Functions{ring:true,parameters,..} if parameters==&[t.clone()])
            || field.parameters!=[t.clone()] || field.gens!=[x.clone()] || field.rep.degree(0)>32 {return Ok(None);}
        let mut terms=Vec::new();
        for (e,c) in field.rep.iter() {
            check()?;
            let Some(d)=c.d.constant_value() else{return Ok(None)};
            for (parameter,n) in c.n.iter() {
                if terms.len()>=1024 || parameter[0]>32{return Ok(None);}
                terms.push((vec![e[0],parameter[0]],Rat::new(n.clone(),d.clone())));
            }
        }
        MPoly::from_terms(2,terms)
    }else{
        if !matches!(h.base.domain,Domain::ZZ|Domain::QQ) || h.base.gens.iter().any(|g|g!=x && g!=t)
            || h.base.rep.len()>1024 || h.base.rep.degrees().iter().any(|&d|d>32){return Ok(None);}
        let parameter=h.base.gens.iter().position(|g|g==t);
        MPoly::from_terms(2,h.base.rep.iter().map(|(e,c)|
            (vec![e[0],parameter.map_or(0,|i|e[i])],c.clone())).collect())
    };
    Ok(Some(polynomial))
}

/// `h` has exact rational coefficients in `(x,t)`, and `q` is ascending
/// in `t`. Degree/bit limits select only this route; refusals retain the
/// complete general algorithm. A signal is always propagated as an error.
pub fn convert<E>(h: &MPoly<Rat>,q: &[Rat],x: &Expr,t: &Expr,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Expr>,E>
{
    check()?;
    if h.nvars!=2 || h.is_zero() || !matches!(q.len(),2|3)
        || q.last().is_some_and(Rat::is_zero) || h.degree(0)>32 || h.degree(1)>32
        || h.len()>1024 || h.coeffs().iter().chain(q).any(|c|c.n.bits().max(c.d.bits())>4096) {
        return Ok(None);
    }
    let expression=convert::to_expr(h,&[x.clone(),t.clone()]);
    let real_term=|root:&Expr|mul2(root,&log(&crate::subs::subs1(&expression,t,root)));
    if q.len()==2 {
        let root=q[0].neg().div(&q[1]).expression();
        check()?;
        let result=real_term(&root);check()?;return Ok(Some(result));
    }
    let four=Rat::from_i64(4);let two=Rat::from_i64(2);
    let discriminant=q[1].mul(&q[1]).sub(&four.mul(&q[2]).mul(&q[0]));
    if discriminant.is_zero(){return Ok(None);}
    let u=q[1].neg().div(&two.mul(&q[2]));
    let radius=discriminant.div(&four.mul(&q[2].mul(&q[2])));
    if !discriminant.is_negative() {
        let v=sqrt(&radius.expression());
        let result=add(&[real_term(&add(&[u.expression(),v.clone()])),
            real_term(&add(&[u.expression(),crate::ops::neg(&v)]))]);
        check()?;return Ok(Some(result));
    }
    let radius=radius.neg();
    let (n,d)=(radius.n.isqrt(),radius.d.isqrt());
    check()?;
    if n.mul(&n)==radius.n && d.mul(&d)==radius.d {
        return conjugate(h,u,Rat::new(n,d),x,check);
    }
    let v=sqrt(&radius.expression());
    let Ok(field)=NumberField::quadratic(&v) else{return Ok(None)};
    let _scope=AlgebraicScope::new(field.clone());
    conjugate(h,AN::rational(u,field.clone()),AN::alpha(field),x,check)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn rational_quadratic_roots_and_conjugates_keep_representatives() {
        let x=Expr::symbol("x");let t=Expr::symbol("t");
        let h=MPoly::from_terms(2,vec![(vec![1,0],Rat::one()),(vec![0,1],Rat::one())]);
        for (q,expected) in [(vec![Rat::one(),Rat::one()],"-log(x-1)"),
            (vec![Rat::one(),Rat::zero(),Rat::one()],"2*atan(x)"),
            (vec![Rat::one(),Rat::zero(),Rat::from_i64(3)],"2*sqrt(3)*atan(sqrt(3)*x)/3")] {
            let result=convert(&h,&q,&x,&t,&mut ||Ok::<_,()>(())).unwrap().unwrap();
            assert_eq!(result,parse_str(expected).unwrap());
        }
    }
    #[test]
    fn repeated_roots_decline_and_signals_never_return_partial_sums() {
        let x=Expr::symbol("x");let t=Expr::symbol("t");
        let h=MPoly::from_terms(2,vec![(vec![1,0],Rat::one()),(vec![0,1],Rat::one())]);
        assert!(convert(&h,&[Rat::one(),Rat::from_i64(2),Rat::one()],&x,&t,&mut ||Ok::<_,()>(())).unwrap().is_none());
        let mut calls=0;
        let result=convert(&h,&[Rat::one(),Rat::zero(),Rat::from_i64(3)],&x,&t,&mut ||{calls+=1;if calls==6{Err("stop")}else{Ok(())}});
        assert!(matches!(result,Err("stop")));
    }
}
