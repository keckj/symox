//! Conservative specialization of certified rational-integral Piecewise rules.
//! A still-symbolic replacement uses the complete original Integral; numeric
//! specializations keep ordinary Piecewise evaluation. Foreign Piecewise and
//! host protocols never select this rule.
use crate::expr::{Expr,Kind};
use rustc_hash::FxHashMap;

fn parameter_pattern(old:&Expr,parameters:&[Expr])->bool {
    if old.has_host() || old.has_calculus_binding(){return false;}
    let free=crate::traverse::free_symbols(old);
    !free.is_empty() && free.iter().all(|symbol|parameters.contains(symbol))
}

pub(super) fn single(e:&Expr,old:&Expr,new:&Expr,transform:&mut impl FnMut(&Expr)->Expr)->Option<Expr> {
    if !matches!(e.kind(),Kind::Piecewise(_)) || new.has_host(){return None;}
    let (original,parameters)=crate::integrate::rational_piecewise_original(e)?;
    if old==new || !parameter_pattern(old,&parameters){return None;}
    let transformed=transform(&original);
    if transformed==original{return None;}
    // A new symbol can be captured by the bound variable, so also retain
    // the original Integral when the RHS is symbolic but its result is closed.
    let symbolic=!crate::traverse::free_symbols(new).is_empty()
        || !crate::traverse::free_symbols(&transformed).is_empty();
    symbolic.then_some(transformed)
}

pub(super) fn exact(e:&Expr,map:&FxHashMap<Expr,Expr>,transform:&mut impl FnMut(&Expr)->Expr)->Option<Expr> {
    if !matches!(e.kind(),Kind::Piecewise(_))
        || map.iter().any(|(old,new)|old.has_host()||new.has_host()){return None;}
    let (original,parameters)=crate::integrate::rational_piecewise_original(e)?;
    let applicable=map.iter().filter(|(old,new)|*old!=*new && parameter_pattern(old,&parameters)).collect::<Vec<_>>();
    if applicable.is_empty(){return None;}
    let limits=crate::binders::limits(&original.args());
    let bound=limits.first()?.first()?;
    // A map may explicitly replace a formula, guard, endpoint or the
    // fallback node itself. Such edits do not specialize its parameters
    // consistently and must keep ordinary exact-replacement semantics.
    for (old,new) in map {
        if old==new || !crate::traverse::has_structural(e,old){continue;}
        if old!=bound && !(parameter_pattern(old,&parameters)
            && crate::traverse::has_structural(&original,old)){return None;}
    }
    let transformed=transform(&original);
    if transformed==original{return None;}
    let symbolic=applicable.iter().any(|(_,new)|!crate::traverse::free_symbols(new).is_empty())
        || !crate::traverse::free_symbols(&transformed).is_empty();
    symbolic.then_some(transformed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    fn p(s:&str)->Expr {parse_str(s).unwrap()}
    fn result()->Expr {crate::integrate::integrate_definite(&p("1/(x+a)"),&p("x"),&Expr::zero(),&Expr::one())}
    #[test]
    fn certified_shape_reverts_but_a_foreign_value_is_unchanged() {
        let value=result();let a=p("a");let b=p("b");
        let (original,_)=crate::integrate::rational_piecewise_original(&value).unwrap();
        assert_eq!(crate::subs::subs1(&value,&a,&b),crate::subs::subs1(&original,&a,&b));
        let Kind::Piecewise(parts)=value.kind() else{panic!()};
        let mut parts=parts.clone();parts[0].0=p("a+1");
        let foreign=crate::relational::piecewise(parts);
        assert!(crate::integrate::rational_piecewise_original(&foreign).is_none());
        assert!(matches!(crate::subs::subs1(&foreign,&a,&b).kind(),Kind::Piecewise(_)));
    }
    #[test]
    fn numeric_specialization_and_bound_substitution_keep_the_formula() {
        let value=result();
        assert_eq!(crate::subs::subs1(&value,&p("a"),&p("2")),p("log(3)-log(2)"));
        assert_eq!(crate::subs::subs1(&value,&p("x"),&p("b")),value);
        assert_eq!(crate::subs::subs1(&value,&p("a"),&p("a")),value);
    }
    #[test]
    fn simultaneous_numeric_masking_retains_at_least_the_original_integral() {
        let value=result();let (original,_)=crate::integrate::rational_piecewise_original(&value).unwrap();
        let pairs=[(p("1/a"),p("2"))];
        let actual=crate::subs::simultaneous(&value,&pairs);
        let expected=crate::subs::simultaneous(&original,&pairs);
        assert_eq!(actual,expected);
        assert_eq!(crate::functions::doit(&actual,true),p("log(3)"));
    }
}
