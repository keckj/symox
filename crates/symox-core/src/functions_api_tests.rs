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

#[test]
fn exact_sign_exclusions_do_not_imply_real(){
    use crate::symbol::Assumptions as A;
    let a=A::from_kwargs([("nonnegative",false)]).unwrap();
    let q=Expr::symbol_with("q",a);
    assert_eq!(crate::assumptions::is_nonnegative(&q),Some(false));
    assert_eq!(crate::assumptions::is_positive(&q),Some(false));
    assert_eq!(crate::assumptions::is_negative(&q),None);
    assert_eq!(crate::assumptions::is_real(&q),None);
    assert_eq!(crate::assumptions::is_zero(&q),Some(false));
    let b=a.with(A::REAL).closure().unwrap();assert!(b.has(A::NEGATIVE));
    let c=A::from_kwargs([("nonnegative",false),("nonpositive",false)]).unwrap().closure().unwrap();assert!(c.has(A::NOT_REAL));
    assert!(c.with(A::REAL).closure().is_err());
}

#[test]
fn held_native_limits_are_scalar_nodes_with_bound_free_symbols(){
    let (x,y)=(Expr::symbol("x"),Expr::symbol("y"));
    let held=crate::limits::held_limit(&parse("f(x)+y"),&x,&Expr::zero(),crate::limits::Dir::Plus);
    assert!(crate::traverse::has(&held,&x));
    assert_eq!(crate::traverse::free_symbols(&held),vec![y]);
    assert_eq!(parse(&crate::printer::srepr(&held)),held);
    assert_eq!(parse(&crate::printer::to_str(&held)),held);
    assert_eq!(doit(&held,true),held);
    let known=crate::limits::held_limit(&parse("sin(x)/x"),&x,&Expr::zero(),crate::limits::Dir::Both);
    assert_eq!(doit(&known,true),Expr::one());
}

#[test]
fn alpha_subs_identity_preserves_public_args_and_ordered_scopes(){
    let (x,y,z)=(Expr::symbol("x"),Expr::symbol("y"),Expr::symbol("z"));
    let a=crate::subs::held_subs(&parse("f(x)"),vec![x.clone()],vec![Expr::zero()]);
    let b=crate::subs::held_subs(&parse("f(y)"),vec![y.clone()],vec![Expr::zero()]);
    assert_eq!(a,b);assert_eq!(a.hash64(),b.hash64());assert_eq!(a.cmp(&b),std::cmp::Ordering::Equal);
    assert_eq!(a.args()[1].args(),vec![x.clone()]);assert_eq!(b.args()[1].args(),vec![y.clone()]);
    let parent=crate::ops::add(&[a.clone(),z.clone()]);assert!(crate::traverse::has(&parent,&b));
    let mut replacements=rustc_hash::FxHashMap::default();replacements.insert(b.clone(),Expr::one());
    assert_eq!(crate::subs::xreplace(&parent,&replacements),parse("z+1"));
    let ordered=crate::subs::held_subs(&parse("x+y"),vec![x.clone(),y.clone()],vec![y.clone(),z.clone()]);
    let other=crate::subs::held_subs(&parse("x+y"),vec![x.clone(),y.clone()],vec![z.clone(),y.clone()]);assert_ne!(ordered,other);
    assert_ne!(crate::subs::held_subs(&parse("f(x)"),vec![x.clone()],vec![x.clone()]),crate::subs::held_subs(&parse("f(y)"),vec![y.clone()],vec![y.clone()]));
    let nested_a=crate::subs::held_subs(&parse("g(Subs(f(x)+y,y,0),x)"),vec![x.clone()],vec![Expr::one()]);
    let nested_b=crate::subs::held_subs(&parse("g(Subs(f(z)+y,y,0),z)"),vec![z],vec![Expr::one()]);assert_eq!(nested_a,nested_b);
    // Rewriting one branch must not reuse an alpha-equal branch's public tree.
    let pair=undefined("pair",vec![a,b]);let mut replacements=rustc_hash::FxHashMap::default();replacements.insert(x,y.clone());
    let rewritten=crate::subs::xreplace(&pair,&replacements);
    assert_eq!(rewritten.args()[0].args()[1].args(),vec![y]);
}

#[test]
fn alpha_subs_key_keeps_shared_scalar_dag_compact(){
    let (x,y)=(Expr::symbol("x"),Expr::symbol("y"));let(mut a,mut b)=(x.clone(),y.clone());
    for _ in 0..32 {a=Expr::from_kind(Kind::Func(Func::Undefined("pair".into()),vec![a.clone(),a]));b=Expr::from_kind(Kind::Func(Func::Undefined("pair".into()),vec![b.clone(),b]));}
    let a=crate::subs::held_subs(&a,vec![x],vec![Expr::zero()]);let b=crate::subs::held_subs(&b,vec![y],vec![Expr::zero()]);
    assert!(a==b);assert_eq!(a.cmp(&b),std::cmp::Ordering::Equal);
}

#[test]
fn floating_facts_are_cached_on_shared_nodes_and_preserve_precision(){
    let mut dag=Expr::symbol("x");
    for _ in 0..40{dag=undefined("pair",vec![dag.clone(),dag]);}
    assert!(!crate::evalf_mp::has_mp(&dag));assert_eq!(crate::evalf_mp::max_precision(&dag),None);
    let mp=crate::number::Number::mp_float(crate::bigfloat::BigFloat::from_decimal("1.25",300).unwrap(),300);
    let machine=Expr::float(1.25);let precise=Expr::num(mp.clone());assert_eq!(machine.as_number().unwrap().cmp_value(&mp),Some(std::cmp::Ordering::Equal));
    for leaves in [vec![machine.clone(),precise.clone()],vec![precise.clone(),machine.clone()]]{
        let call=undefined("pair",leaves);assert!(crate::evalf_mp::has_mp(&call));assert_eq!(crate::evalf_mp::max_precision(&call),Some(300));
    }
    let x=Expr::symbol("x");
    for raw in [Kind::Add(mp.clone(),vec![(x.clone(),Number::ONE)]),Kind::Add(Number::ZERO,vec![(x.clone(),mp.clone())]),Kind::Mul(mp.clone(),vec![(x.clone(),Expr::one())]),Kind::Mul(Number::ONE,vec![(x,precise.clone())])]{
        let e=Expr::from_kind(raw);assert!(crate::evalf_mp::has_mp(&e));assert_eq!(crate::evalf_mp::max_precision(&e),Some(300));
    }
    let matrix=crate::matrixexpr::MatrixExpr::explicit(crate::matrix::Matrix::from_rows(vec![vec![precise]]).unwrap());
    let trace=Expr::from_kind(Kind::Func(Func::MatrixTrace(std::sync::Arc::new(matrix)),vec![]));
    assert!(crate::evalf_mp::has_mp(&trace));assert_eq!(crate::evalf_mp::max_precision(&trace),Some(300));
}
