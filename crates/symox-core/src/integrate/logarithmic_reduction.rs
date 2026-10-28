//! Lazard--Rioboo--Trager reduction over the stock rational coefficient rings.
//!
//! Keep the signed square-free groups and the actual PRS representatives.
//! Replacing these by arbitrary associates would preserve a derivative but
//! need not preserve the logarithm branches used by the caller.
use crate::poly::{checked_factors, dense, fraction::Polynomial, mpoly::MPoly,
    polyobj::{Domain, Poly}, pseudo::Budget, Rat};
use crate::poly::factor_budget::Limits;
use crate::Expr;

fn embed(p: &MPoly<Rat>) -> MPoly<Rat> {
    MPoly::from_terms(2,p.iter().map(|(e,c)|(vec![0,e[0]],c.clone())).collect())
}
fn parameter(p: &MPoly<Rat>) -> MPoly<Rat> {
    MPoly::from_terms(1,p.iter().map(|(e,c)|(vec![e[1]],c.clone())).collect())
}
fn natural(p: MPoly<Rat>, gens: Vec<Expr>, ngens: usize) -> Polynomial {
    let domain=if p.is_integral(){Domain::ZZ}else{Domain::QQ};
    Poly::new(p,gens,ngens,domain).into()
}

/// The dispatch is deliberately bounded; a refusal restarts the complete
/// existing algorithm. Signal checks are errors, never converted to refusals.
/// The supplied polynomials have their sole generator equal to `x`.
pub fn reduce<E>(f: &Poly, g: &Poly, t: &Expr, check: &mut impl FnMut()->Result<(),E>)
    -> Result<Option<Vec<(Polynomial,Polynomial)>>,E>
{
    check()?;
    if f.rep.is_zero() || g.rep.is_zero() || f.rep.degree(0)>=g.rep.degree(0)
        || g.rep.degree(0)>32 || f.rep.nvars!=1 || g.rep.nvars!=1
        || f.rep.coeffs().iter().chain(g.rep.coeffs()).any(|c|c.n.bits().max(c.d.bits())>4096) {
        return Ok(None);
    }
    macro_rules! step {($e:expr)=>{{check()?; let value=$e; check()?;
        match value {Ok(v)=>v,Err(_)=>return Ok(None)}}}}
    let lift=|p:&MPoly<Rat>|MPoly::from_terms(2,p.iter().map(|(e,c)|(vec![e[0],0],c.clone())).collect());
    let a=lift(&g.rep);
    let b=lift(&f.rep).sub(&lift(&g.rep.derivative(0)).mul_term(&[0,1],&Rat::one()));
    let prs=step!(crate::poly::prs::subresultants(&a,&b,0));
    if prs.resultant.is_zero() || prs.resultant.degree(0)>0
        || prs.resultant.degree(1)!=g.rep.degree(0)
        || !prs.resultant.iter().any(|(e,_)|e[1]==0) {return Ok(None);}
    // The leading t coefficient is (-1)^deg(g)*Res(g,g'), so its
    // nonvanishing certifies squarefreeness. Its constant coefficient is
    // lc(g)^(deg(g)-1-deg(f))*Res(g,f), certifying coprimality. Reuse it
    // instead of repeating two rational Euclidean chains before the PRS.
    let res=natural(parameter(&prs.resultant),vec![t.clone()],1);
    let mut sqf=step!(checked_factors::compatible_list(&res,true,false,Limits::default()));
    // The complete negative unit is included in the first group before its
    // primitive part is taken, exactly as in rationaltools._include_sign.
    if sqf.unit.rep.lc().is_negative() {
        let Some((first,_))=sqf.factors.first_mut() else{return Ok(None)};
        first.base.rep=first.rep.scale(&sqf.unit.rep.lc());
    }
    let mut budget=Budget::new();
    let mut result=Vec::with_capacity(sqf.factors.len());
    for (q,i) in sqf.factors {
        check()?;
        let (_,q)=q.base.primitive();
        if g.rep.degree(0)==i as i64 {
            result.push((g.clone().into(),q.into()));
            continue;
        }
        let Some(mut h)=prs.sequence.iter().rev().find(|p|p.degree(0)==i as i64).cloned()
            else{return Ok(None)};
        let lc=parameter(&h.lc_in(0));
        let lc_poly=Poly::new(lc.clone(),vec![t.clone()],1,Domain::QQ);
        let factors=step!(checked_factors::compatible_list(&lc_poly.clone().into(),true,true,Limits::default()));
        let qd=dense::from_mpoly(&q.rep,0);
        for (factor,j) in factors.factors {
            check()?;
            // QQ gcd is monic, so multiplying the first LC factor by its
            // negative scalar unit does not change this gcd or its power.
            let common=dense::gcd_field(&dense::from_mpoly(&factor.rep,0),&qd);
            let divisor=embed(&dense::to_mpoly(&common,0,1));
            let divisor=step!(budget.pow(&divisor,u64::from(j)));
            h=step!(budget.exact_div(&h,&divisor));
        }
        check()?;
        let (gcd,inv,_)=dense::gcdex_field(&dense::from_mpoly(&lc,0),&qd);
        if gcd!=vec![Rat::one()] {return Ok(None);}
        let inv=dense::rem(&inv,&qd);
        let leading=h.degree(0);
        let mut terms=vec![(vec![leading as u32,0],Rat::one())];
        // Poly.coeffs/monoms omit zero coefficients. Iterate x degrees in
        // descending order and reduce each retained coefficient modulo q.
        for degree in (0..leading).rev() {
            check()?;
            let coefficient=MPoly::from_terms(1,h.iter().filter(|(e,_)|e[0]==degree as u32)
                .map(|(e,c)|(vec![e[1]],c.clone())).collect());
            if coefficient.is_zero(){continue;}
            let reduced=dense::rem(&dense::mul(&inv,&dense::from_mpoly(&coefficient,0)),&qd);
            terms.extend(reduced.into_iter().enumerate().filter(|(_,c)|!c.is_zero())
                .map(|(j,c)|(vec![degree as u32,j as u32],c)));
        }
        let h=MPoly::from_terms(2,terms);
        let h=if h.degree(1)<=0 {
            natural(MPoly::from_terms(1,h.iter().map(|(e,c)|(vec![e[0]],c.clone())).collect()),g.gens.clone(),1)
        }else{
            let p=natural(h,vec![g.gens[0].clone(),t.clone()],1);
            step!(checked_factors::normalize_natural_factor(p))
        };
        result.push((h,q.into()));
    }
    check()?;
    Ok(Some(result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn signed_groups_and_prs_representatives_are_preserved() {
        let x=Expr::symbol("x"); let t=Expr::symbol("t");
        for (f,g,h,q) in [("1","x^2+x+1","x+3*t/2+1/2","3*t^2+1"),
            ("12","x^2-x-2","x-3*t/8-1/2","-t^2+16")] {
            let poly=|s|Poly::from_expr(&parse_str(s).unwrap(),std::slice::from_ref(&x)).unwrap();
            let r=reduce(&poly(f),&poly(g),&t,&mut ||Ok::<_,()>(())).unwrap().unwrap();
            assert_eq!(r.len(),1);
            assert_eq!(r[0].0.as_expr(),parse_str(h).unwrap());
            assert_eq!(r[0].1.as_expr(),parse_str(q).unwrap());
        }
    }
    #[test]
    fn invalid_inputs_decline_and_interruptions_are_errors() {
        let x=Expr::symbol("x");let t=Expr::symbol("t");
        let p=|s|Poly::from_expr(&parse_str(s).unwrap(),std::slice::from_ref(&x)).unwrap();
        assert!(reduce(&p("1"),&p("x^2"),&t,&mut ||Ok::<_,()>(())).unwrap().is_none());
        let mut calls=0;
        let r=reduce(&p("1"),&p("x^2+x+1"),&t,&mut ||{calls+=1;if calls==7{Err("stop")}else{Ok(())}});
        assert!(matches!(r,Err("stop")));
    }
}
