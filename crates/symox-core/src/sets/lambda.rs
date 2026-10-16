//! Lambda signatures and exact simultaneous application.
use std::sync::OnceLock;
use rustc_hash::{FxHashMap,FxHashSet};
use super::{as_set,mk,SetOp};
use crate::expr::{Expr,Func,Kind};

fn tuple_args(e:&Expr)->Option<&[Expr]> {
    match e.kind(){Kind::Func(Func::Undefined(head),args) if head.is_internal("Tuple")=>Some(args),_=>None}
}

/// The argument-free identity function is independent of user symbol names.
pub fn identity_function()->Expr {
    static VALUE:OnceLock<Expr>=OnceLock::new();
    VALUE.get_or_init(||mk(SetOp::IdentityFunction,vec![])).clone()
}
fn identity_variable()->Expr {
    static VARIABLE:OnceLock<Expr>=OnceLock::new();
    VARIABLE.get_or_init(||Expr::dummy("x",crate::Assumptions::NONE)).clone()
}

/// Normalized signature. IdentityFunction has a private, stable dummy variable.
pub fn lambda_signature(e:&Expr)->Option<Expr> {
    match as_set(e)? {
        (SetOp::Lambda,args) if args.len()==2=>Some(if args[0].is_symbol(){crate::special::tuple(vec![args[0].clone()])}else{args[0].clone()}),
        (SetOp::IdentityFunction,args) if args.is_empty()=>Some(crate::special::tuple(vec![identity_variable()])),
        _=>None,
    }
}

/// Number of top-level arguments, including destructured tuple arguments.
pub fn lambda_arity(e:&Expr)->Option<usize>{Some(tuple_args(&lambda_signature(e)?)?.len())}

fn variables(signature:&Expr)->Result<Vec<Expr>,String> {
    fn walk(e:&Expr,out:&mut Vec<Expr>,seen:&mut FxHashSet<Expr>)->Result<(),String>{
        if e.is_symbol(){
            if !seen.insert(e.clone()){return Err(format!("Duplicate symbol {e}"));}
            out.push(e.clone());Ok(())
        }else if let Some(args)=tuple_args(e){
            for arg in args{walk(arg,out,seen)?;}Ok(())
        }else{Err(format!("Lambda signature should be only tuples and symbols, not {e}"))}
    }
    if tuple_args(signature).is_none(){return Err("Lambda signature must be a Tuple".into());}
    let mut out=Vec::new();walk(signature,&mut out,&mut FxHashSet::default())?;Ok(out)
}

/// Validate normalized Lambda payloads; constructors also accept a single symbol.
pub fn validate_lambda(op:SetOp,args:&[Expr])->Result<(),String>{
    match op {
        SetOp::IdentityFunction if args.is_empty()=>Ok(()),
        SetOp::IdentityFunction=>Err("IdentityFunction takes no constructor arguments".into()),
        SetOp::Lambda if args.len()==2=>variables(&args[0]).map(|_|()),
        _=>Err("Lambda requires a signature and an expression".into()),
    }
}

/// Construct a held Lambda, retaining nested and empty tuple signatures.
pub fn lambda_checked(signature:&Expr,body:&Expr)->Result<Expr,String>{
    let signature=if signature.is_symbol(){crate::special::tuple(vec![signature.clone()])}else{signature.clone()};
    variables(&signature)?;
    if tuple_args(&signature).is_some_and(|args|args.len()==1 && args[0].is_symbol() && args[0]==*body){return Ok(identity_function());}
    Ok(mk(SetOp::Lambda,vec![signature,body.clone()]))
}

/// Flattened bound variables and body, for scalar and matrix consumers.
pub fn lambda_parts(e:&Expr)->Option<(Vec<Expr>,Expr)>{
    let signature=lambda_signature(e)?;
    let body=match as_set(e)?{(SetOp::IdentityFunction,_)=>identity_variable(),(SetOp::Lambda,args)=>args.get(1)?.clone(),_=>return None};
    Some((variables(&signature).ok()?,body))
}

/// Exact simultaneous replacement, matching SymPy's Lambda.__call__.
/// Capture by a binder in the body is deliberately not prevented.
pub fn lambda_apply(e:&Expr,args:&[Expr])->Result<Expr,String>{
    let signature=lambda_signature(e).ok_or_else(||"not a Lambda".to_string())?;
    let parameters=tuple_args(&signature).ok_or_else(||"invalid Lambda signature".to_string())?;
    if args.len()!=parameters.len(){return Err(format!("{e} takes exactly {} argument{} ({} given)",parameters.len(),if parameters.len()==1{""}else{"s"},args.len()));}
    fn match_args(parameters:&[Expr],args:&[Expr],map:&mut FxHashMap<Expr,Expr>)->Result<(),String>{
        if parameters.len()!=args.len(){return Err("Lambda tuple argument has the wrong length".into());}
        for(parameter,arg)in parameters.iter().zip(args){
            if parameter.is_symbol(){map.insert(parameter.clone(),arg.clone());}
            else if let(Some(pars),Some(values))=(tuple_args(parameter),tuple_args(arg)){match_args(pars,values,map)?;}
            else{return Err(format!("Can't match {arg} and {parameter}"));}
        }Ok(())
    }
    let mut map=FxHashMap::default();match_args(parameters,args,&mut map)?;
    let(_,body)=lambda_parts(e).ok_or_else(||"invalid Lambda signature".to_string())?;
    let result=crate::subs::xreplace(&body,&map);
    validate_lambda_tree(&result)?;
    crate::binders::validate_tree(&result)?;
    crate::limits::validate_held(&result)?;
    Ok(result)
}

/// Check transformed signatures without revisiting shared DAG nodes.
pub fn validate_lambda_tree(e:&Expr)->Result<(),String>{
    fn walk(e:&Expr,seen:&mut FxHashMap<usize,Expr>)->Result<(),String>{
        if !e.has_lambda_binding() || seen.insert(e.dag_id(),e.clone()).is_some(){return Ok(());}
        if let Some((op,args))=as_set(e){if matches!(op,SetOp::Lambda|SetOp::IdentityFunction){validate_lambda(op,args)?;}}
        if let Kind::Func(Func::MatrixDotProduct(a,b),_)=e.kind(){for m in [a,b]{let mut fields=Vec::new();m.scalar_fields(&mut fields);for scalar in fields{walk(&scalar,seen)?;}}}
        if let Kind::Func(Func::MatrixElement(m)|Func::MatrixProxy(m)|Func::MatrixTrace(m)|Func::MatrixDeterminant(m),_)=e.kind(){
            let mut fields=Vec::new();m.scalar_fields(&mut fields);for scalar in fields{walk(&scalar,seen)?;}
        }
        for arg in e.args(){walk(&arg,seen)?;}Ok(())
    }walk(e,&mut FxHashMap::default())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ops,parser::parse_str,traverse};
    fn p(text:&str)->Expr{parse_str(text).unwrap()}

    #[test]
    fn signatures_nested_empty_and_exact_application(){
        let(x,y,z)=(Expr::symbol("x"),Expr::symbol("y"),Expr::symbol("z"));
        let signature=crate::special::tuple(vec![x.clone(),y.clone()]);
        let lam=lambda_checked(&signature,&p("x-y+z")).unwrap();
        assert_eq!(lambda_arity(&lam),Some(2));
        assert_eq!(lambda_parts(&lam).unwrap().0,vec![x.clone(),y.clone()]);
        assert_eq!(lambda_apply(&lam,&[y.clone(),x.clone()]).unwrap(),p("y-x+z"));
        assert!(lambda_apply(&lam,&[x.clone()]).is_err());
        assert_eq!(traverse::free_symbols(&lam),vec![z.clone()]);
        let nested=p("Lambda((x,(y,z)),x+2*y+3*z)");
        assert_eq!(lambda_parts(&nested).unwrap().0,vec![x.clone(),y.clone(),z.clone()]);
        assert_eq!(lambda_apply(&nested,&[Expr::int(1),crate::special::tuple(vec![Expr::int(2),Expr::int(3)])]).unwrap(),Expr::int(14));
        assert!(lambda_apply(&nested,&[Expr::one(),crate::special::tuple(vec![Expr::int(2)])]).is_err());
        assert_eq!(lambda_apply(&p("Lambda((),z)"),&[]).unwrap(),z);
        assert!(lambda_checked(&crate::special::tuple(vec![x.clone(),x.clone()]),&x).is_err());
        assert!(lambda_checked(&ops::add2(&x,&Expr::one()),&x).is_err());
    }

    #[test]
    fn structural_scope_capture_and_transformed_validation(){
        let(x,y,z)=(Expr::symbol("x"),Expr::symbol("y"),Expr::symbol("z"));
        let lam=p("Lambda((x,y),x+y+z)");
        assert_ne!(p("Lambda(x,x+z)"),p("Lambda(y,y+z)"));
        let renamed=crate::subs::subs1(&lam,&x,&Expr::symbol("a"));
        assert_eq!(renamed,p("Lambda((a,y),a+y+z)"));
        let invalid=crate::subs::subs1(&lam,&x,&Expr::int(2));
        assert!(validate_lambda_tree(&invalid).is_err());
        assert_eq!(lambda_apply(&p("Lambda(x,Lambda(y,x+y))"),&[y.clone()]).unwrap(),p("Lambda(y,2*y)"));
        assert_eq!(lambda_apply(&p("Lambda(x,Integral(f(x,y),(y,0,1)))"),&[y]).unwrap(),p("Integral(f(y,y),(y,0,1))"));
        assert_eq!(crate::diff::diff(&lam,&x),Expr::zero());
        assert!(matches!(crate::diff::diff(&lam,&z).kind(),Kind::Derivative(..)));
        assert!(traverse::has(&lam,&x));
    }

    #[test]
    fn identity_and_printer_parser_roundtrips(){
        let x=Expr::symbol("x");let identity=lambda_checked(&x,&x).unwrap();
        assert!(identity.ptr_eq(&identity_function()));
        assert!(identity.args().is_empty());
        assert_eq!(identity.class_name(),"IdentityFunction");
        assert_eq!(identity.to_string(),"Lambda(_x, _x)");
        assert_eq!(crate::printer::srepr(&identity),"IdentityFunction()");
        assert_eq!(p(&crate::printer::srepr(&identity)),identity);
        assert!(traverse::free_symbols(&identity).is_empty());
        assert_eq!(lambda_apply(&identity,&[x.clone()]).unwrap(),x);
        for source in ["Lambda((x,y),x-y)","Lambda((x,(y,z)),f(x,y,z))","Lambda((),3)","Lambda((x,y),(x,y))"] {
            let value=p(source);assert_eq!(p(&crate::printer::srepr(&value)),value);
            assert_eq!(p(&value.to_string()),value);
            let _=crate::sortkey::sort_key(&value);
        }
        let undefined=crate::functions::func(Func::Undefined(crate::expr::UndefinedFunction::new("Lambda",crate::Assumptions::NONE,Some(vec![2]))),vec![Expr::symbol("x"),Expr::symbol("y")]);
        assert!(lambda_parts(&undefined).is_none());
        assert_eq!(p(&crate::printer::srepr(&undefined)),undefined);
    }

    #[test]
    fn application_of_a_shared_body_visits_the_dag_once(){
        let x=Expr::symbol("x");let mut body=x.clone();
        for _ in 0..48{body=crate::functions::undefined("pair",vec![body.clone(),body]);}
        let value=lambda_checked(&x,&body).unwrap();
        let result=lambda_apply(&value,&[Expr::int(7)]).unwrap();
        let mut node=result;
        for _ in 0..48{let args=node.args();assert!(args[0].ptr_eq(&args[1]));node=args[0].clone();}
        assert_eq!(node,Expr::int(7));
    }
}
