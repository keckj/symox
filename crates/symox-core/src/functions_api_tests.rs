use super::*;
use crate::expr::UndefinedFunction;
fn parse(source:&str)->Expr{crate::parser::parse_str(source).unwrap()}

#[test]
fn undefined_identity_metadata_and_roundtrip(){
    let x=Expr::symbol("x");
    let spec=UndefinedFunction::new("sin",crate::symbol::Assumptions::REAL,Some(vec![1]));
    let applied=func(Func::Undefined(spec.clone()),vec![x.clone()]);
    assert_ne!(applied,sin(&x));
    assert_eq!(crate::assumptions::is_real(&applied),Some(true));
    assert_eq!(crate::evalf::eval_c64(&func(Func::Undefined(spec),vec![Expr::one()])),None);
    assert_eq!(parse(&crate::printer::srepr(&applied)),applied);
    let tuple_named=func(Func::Undefined(UndefinedFunction::new("Tuple",crate::symbol::Assumptions::NONE,None)),vec![x]);
    assert_eq!(parse(&crate::printer::srepr(&tuple_named)),tuple_named);
    assert_ne!(tuple_named,undefined("Tuple",tuple_named.args()));
}

#[test]
fn held_subs_binding_evaluation_and_undefined_chain_rule(){
    let (x,y)=(Expr::symbol("x"),Expr::symbol("y"));
    let held=crate::subs::held_subs(&parse("x**2+y"),vec![x.clone()],vec![Expr::int(2)]);
    assert_eq!(crate::traverse::free_symbols(&held),vec![y.clone()]);
    assert_eq!(doit(&held,true),parse("y+4"));
    assert_eq!(parse(&crate::printer::srepr(&held)),held);
    let differential=crate::diff::diff(&parse("f(x*x)"),&x);
    assert!(crate::printer::to_str(&differential).contains("Subs(Derivative"));
    assert_eq!(doit(&crate::subs::subs1(&differential,&x,&Expr::zero()),true),Expr::zero());
    let held=crate::subs::held_subs(&parse("f(x)"),vec![x.clone()],vec![y.clone()]);
    assert_eq!(crate::subs::subs1(&held,&x,&Expr::zero()),held);
    assert_eq!(doit(&held,true),parse("f(y)"));
    let derivative=crate::diff::derivative_unevaluated(&parse("f(x,y)"),vec![(x.clone(),1)]);
    assert_eq!(crate::subs::subs1(&derivative,&y,&x),crate::subs::held_subs(&derivative,vec![y],vec![x]));
}

#[test]
fn simultaneous_symbol_substitutions_and_unknown_limits(){
    let (x,y)=(Expr::symbol("x"),Expr::symbol("y"));
    let pairs=vec![(x.clone(),y.clone()),(y.clone(),x.clone())];
    assert_eq!(crate::subs::simultaneous(&parse("x+y"),&pairs),parse("x+y"));
    let held=crate::subs::held_subs(&parse("f(x)"),vec![x.clone()],vec![y.clone()]);
    assert_eq!(crate::subs::simultaneous(&held,&pairs),crate::subs::held_subs(&parse("f(x)"),vec![x.clone()],vec![x.clone()]));
    let derivative=crate::diff::derivative_unevaluated(&parse("f(x)"),vec![(x.clone(),1)]);
    assert_eq!(crate::subs::simultaneous(&derivative,&[(x.clone(),Expr::int(2))]),crate::subs::held_subs(&derivative,vec![x.clone()],vec![Expr::int(2)]));
    let pairs=vec![(x.clone(),Expr::zero()),(y,Expr::zero())];
    assert!(crate::subs::simultaneous(&parse("x/y"),&pairs).is_constant(Constant::NaN));
    assert!(crate::limits::limit(&parse("f(x)"),&x,&Expr::zero(),crate::limits::Dir::Plus).is_err());
    assert_eq!(crate::limits::limit(&parse("sin(x)/x"),&x,&Expr::zero(),crate::limits::Dir::Plus).unwrap(),Expr::one());
}
