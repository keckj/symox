//! Certified conditional values for retained rational segments and rays.
//!
//! For finite complex c with -c outside the closed real interval, the path
//! x+c is horizontal and avoids zero. If it lies on the negative real axis,
//! both endpoint principal logarithms have the same +i*pi contribution.
//! Otherwise its imaginary part has constant nonzero sign. Thus their
//! difference is the integral of 1/(x+c), without a logarithmic branch jump.
use crate::expr::{Expr,Kind};
use crate::functions::log;
use crate::logic::{and,not};
use crate::ops::{add2,div,mul2,neg,sqrt,sub};

fn finite_complex(c: &Expr) -> Expr {
    // Retain symbolic membership even when declared assumptions already imply
    // finiteness: later substitutions must still exclude NaN and infinities.
    crate::sets::contains_expr(c,&crate::sets::complexes())
}

fn outside(c: &Expr) -> Expr {
    let interval=crate::sets::interval(&Expr::zero(),&Expr::one(),false,false);
    not(&crate::sets::contains_expr(&neg(c),&interval))
}

/// The principal logarithm and square root are analytic off this closed
/// negative-real ray. Keep symbolic membership through later substitutions.
fn off_cut(a: &Expr) -> Option<Expr> {
    let cut=crate::sets::interval(&Expr::neg_oo(),&Expr::zero(),true,false);
    let excluded=not(&crate::sets::contains_expr(a,&cut));
    // A declared nonreal parameter can make this True immediately. Such a
    // guard would disappear and a later contradictory real substitution
    // could select a false finite continuation. Keep the old retained path.
    if matches!(excluded.kind(),Kind::Bool(true)) && !crate::traverse::free_symbols(a).is_empty() {
        return None;
    }
    let condition=and(vec![finite_complex(a),excluded]);
    (!matches!(condition.kind(),Kind::Bool(false))).then_some(condition)
}

fn conditional(value:Expr,condition:Expr,f:&Expr,x:&Expr,a:&Expr,b:&Expr)->Option<Expr> {
    if matches!(condition.kind(),Kind::Bool(false)){return None;}
    Some(crate::relational::piecewise(vec![
        (value,condition),
        (super::integral_definite_unevaluated(f,x,a,b),Expr::boolean(true)),
    ]))
}

/// Extend only retained results. The logarithm along a->1 stays on one side
/// of its cut. For the tail let s=sqrt(parameter), Re(s)>0; partial fractions
/// along x-is and x+is have no log jumps, with endpoint argument difference pi.
pub(super) fn retained(f:&Expr,x:&Expr,a:&Expr,b:&Expr)->Option<Expr> {
    if !b.is_one() && !(a.is_zero() && b.is_constant(crate::expr::Constant::Infinity)) {
        return None;
    }
    let (num,den)=crate::simplify::numer_denom(f);
    if num.is_one() && x.is_symbol() {
        let dc=super::poly_coeffs(&den,x)?;
        if dc.len()==2 && dc[0].is_zero() && dc[1].is_one() && b.is_one()
            && !super::has(a,x) && crate::assumptions::fact(a,"commutative")==Some(true) {
            let condition=off_cut(a)?;
            return conditional(neg(&log(a)),condition,f,x,a,b);
        }
        if a.is_zero() && b.is_constant(crate::expr::Constant::Infinity)
            && dc.len()==3 && dc[1].is_zero() && dc[2].is_one()
            && crate::assumptions::fact(&dc[0],"commutative")==Some(true) {
            let parameter=&dc[0];
            let condition=off_cut(parameter)?;
            return conditional(div(&Expr::pi(),&mul2(&Expr::int(2),&sqrt(parameter))),
                condition,f,x,a,b);
        }
    }
    unit_interval(f,x,a,b)
}

/// Only the undecided-pole route calls this: already established regular or
/// divergent results are untouched. Conditions describe complex membership,
/// never order a parameter that might be nonreal.
pub(super) fn unit_interval(f: &Expr,x: &Expr,a: &Expr,b: &Expr) -> Option<Expr> {
    if !a.is_zero() || !b.is_one() || !x.is_symbol() { return None; }
    let (numerator,denominator)=crate::simplify::numer_denom(f);
    if !numerator.is_one() { return None; }
    let coefficients=super::poly_coeffs(&denominator,x)?;
    if coefficients.len()!=2 { return None; }
    let held=super::integral_definite_unevaluated(f,x,a,b);
    if coefficients[1].is_one() {
        let c=&coefficients[0];
        if crate::assumptions::fact(c,"commutative")!=Some(true) { return None; }
        let condition=and(vec![finite_complex(c),outside(c)]);
        if matches!(condition.kind(),Kind::Bool(false)) { return None; }
        let value=sub(&log(&add2(&Expr::one(),c)),&log(c));
        return Some(crate::relational::piecewise(vec![
            (value,condition),(held,Expr::boolean(true)),
        ]));
    }
    if !coefficients[0].is_one() { return None; }
    let slope=&coefficients[1];
    if crate::assumptions::fact(slope,"commutative")!=Some(true) { return None; }
    let zero=crate::relational::eq(slope,&Expr::zero());
    // An assumed-nonzero parameter would lose this branch during normal
    // Piecewise construction, and substitution of zero could then select
    // an undefined continuation. Retain the full old Integral in that case.
    if matches!(zero.kind(),Kind::Bool(false)) && !crate::traverse::free_symbols(slope).is_empty() {
        return None;
    }
    let c=div(&Expr::one(),slope);
    let condition=and(vec![finite_complex(slope),crate::relational::ne(slope,&Expr::zero()),outside(&c)]);
    // On this condition the identity x*slope+1=slope*(x+1/slope)
    // proves the same horizontal-segment formula. Do not replace the two
    // endpoint logarithms by log(1+slope), whose branch may differ.
    let value=div(&sub(&log(&add2(&Expr::one(),&c)),&log(&c)),slope);
    Some(crate::relational::piecewise(vec![
        (Expr::one(),zero),
        (value,condition),(held,Expr::boolean(true)),
    ]))
}

/// Recognize the complete certified value/domain construction, not merely
/// a domain guard followed by an arbitrary Integral default. Rebuilding uses
/// only these fixed rules and no integration/general equivalence algorithm.
pub(crate) fn piecewise_original(e:&Expr)->Option<(Expr,Vec<Expr>)> {
    use crate::expr::Func;
    use crate::sets::SetOp;
    if e.has_host() || e.has_held_arithmetic() {return None;}
    let Kind::Piecewise(pieces)=e.kind() else{return None};
    if !(2..=3).contains(&pieces.len()){return None;}
    let (original,last_condition)=pieces.last()?;
    if !matches!(last_condition.kind(),Kind::Bool(true)){return None;}
    let Kind::Func(Func::Integral,args)=original.kind() else{return None};
    if args.len()!=2{return None;}
    let limits=crate::binders::limits(args);
    if limits.len()!=1 || limits[0].len()!=3{return None;}
    let mut complex=false;let mut interval=false;
    for (_,condition) in &pieces[..pieces.len()-1] {
        for (node,_) in crate::inspection::scalar_nodes(condition).iter() {
            if let Kind::Func(Func::Set(SetOp::Contains),values)=node.kind() {
                if values.len()!=2{continue;}
                complex|=values[1]==crate::sets::complexes();
                interval|=matches!(crate::sets::as_set(&values[1]),Some((SetOp::Interval,_)));
            }
        }
    }
    if !complex || !interval{return None;}
    let parameters=crate::traverse::free_symbols(original);
    if parameters.is_empty(){return None;}
    let [x,a,b]=limits[0].as_slice() else{return None};
    let rebuilt=retained(&args[0],x,a,b)?;
    (rebuilt==*e).then(||(original.clone(),parameters))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    fn p(s:&str)->Expr {parse_str(s).unwrap()}
    #[test]
    fn exact_domain_excludes_all_poles_and_keeps_the_logarithms_separate() {
        let f=p("1/(x+a)");let x=Expr::symbol("x");
        let result=unit_interval(&f,&x,&Expr::zero(),&Expr::one()).unwrap();
        let Kind::Piecewise(pieces)=result.kind() else{panic!("undecided domain must remain conditional")};
        assert_eq!(pieces[0].0,p("log(1+a)-log(a)"));
        for c in ["-1","-1/2","0"] {
            assert!(unit_interval(&crate::subs::subs1(&f,&Expr::symbol("a"),&p(c)),&x,&Expr::zero(),&Expr::one()).is_none());
        }
        for c in ["-2","1","I","-2+I"] {
            let value=unit_interval(&crate::subs::subs1(&f,&Expr::symbol("a"),&p(c)),&x,&Expr::zero(),&Expr::one()).unwrap();
            assert!(!matches!(value.kind(),Kind::Piecewise(_)));
        }
    }
    #[test]
    fn reciprocal_slope_preserves_zero_and_both_logarithms() {
        let a=Expr::symbol("a");let x=Expr::symbol("x");
        let result=unit_interval(&p("1/(1+a*x)"),&x,&Expr::zero(),&Expr::one()).unwrap();
        let Kind::Piecewise(pieces)=result.kind() else{panic!("slope domain must remain conditional")};
        assert_eq!(pieces.len(),3);
        assert_eq!(pieces[0].0,Expr::one());
        assert_eq!(pieces[1].0,p("(log(1+1/a)-log(1/a))/a"));
        assert_eq!(crate::subs::subs1(&result,&a,&Expr::zero()),Expr::one());
        for v in ["-2","-1"] {
            let at=crate::subs::subs1(&result,&a,&p(v));
            assert!(matches!(at.kind(),Kind::Func(crate::expr::Func::Integral,_)));
        }
    }

    #[test]
    fn slit_domains_certify_log_segment_and_quadratic_tail() {
        let a=Expr::symbol("a");let x=Expr::symbol("x");
        let lower=retained(&p("1/x"),&x,&a,&Expr::one()).unwrap();
        let tail=retained(&p("1/(x**2+a)"),&x,&Expr::zero(),&Expr::oo()).unwrap();
        for (result,want) in [(&lower,p("-log(a)")),(&tail,p("pi/(2*sqrt(a))"))] {
            let Kind::Piecewise(parts)=result.kind() else{panic!("unknown slit domain must remain conditional")};
            assert_eq!(parts[0].0,want);
            for point in ["-2","-1","0","oo","-oo","zoo","nan"] {
                assert!(matches!(crate::subs::subs1(result,&a,&p(point)).kind(),Kind::Func(crate::expr::Func::Integral,_)));
            }
        }
        for point in ["1","2","I","-2+I"] {
            assert!(!matches!(crate::subs::subs1(&lower,&a,&p(point)).kind(),Kind::Piecewise(_)));
            assert!(!matches!(crate::subs::subs1(&tail,&a,&p(point)).kind(),Kind::Piecewise(_)));
        }
    }

}
