//! Held Integral/Sum/Product nodes and their ordered binding scopes.
use rustc_hash::{FxHashMap,FxHashSet};
use crate::expr::{Expr,Func,Kind};
use crate::ops::{add,mul2,neg};

pub fn is_calculus(f:&Func)->bool { matches!(f,Func::Integral|Func::Sum|Func::Product) }
pub fn limits(args:&[Expr])->Vec<Vec<Expr>> { args[1..].iter().map(Expr::args).collect() }

/// SymPy's Integral also binds an applied undefined function such as f(x).
pub fn is_function_variable(v:&Expr)->bool {matches!(v.kind(),Kind::Func(Func::Undefined(head),_) if !head.internal)}

/// Variables which SymPy accepts as a new indefinite integration variable
/// (its `_diff_wrt` rule restricted to representable binding variables).
fn renames_variable(v:&Expr)->bool {v.is_symbol() || is_function_variable(v)}

fn anonymous_dummy()->Expr {Expr::from_kind(Kind::Sym(crate::symbol::Symbol::anonymous_dummy(crate::symbol::Assumptions::NONE)))}

pub fn validate(f:&Func,args:&[Expr])->Result<(),String> {
    if !is_calculus(f) || args.len()<2 { return Err("Integral/Sum/Product requires an integrand and limits".into()); }
    if matches!(args[0].kind(),Kind::Func(Func::MatrixProxy(_),_)) { return Err("matrix-valued calculus binders are not supported".into()); }
    for limit in &args[1..] {
        if !matches!(limit.kind(),Kind::Func(Func::Undefined(head),_) if head.is_internal("Tuple")) {
            return Err("calculus limits must be Tuples".into());
        }
        let values=limit.args();
        if !(1..=3).contains(&values.len()) || (*f!=Func::Integral && values.len()!=3) {
            return Err("Integral expects (x,), (x,b), or (x,a,b); Sum and Product expect (k,a,b)".into());
        }
        if !values[0].is_symbol() && !(*f==Func::Integral && is_function_variable(&values[0])) {
            return Err("calculus binding variables must be symbols".into());
        }
        if values[1..].iter().any(|bound|matches!(bound.kind(),Kind::Func(Func::MatrixProxy(_),_))) {
            return Err("calculus bounds must be scalar expressions".into());
        }
    }
    Ok(())
}

pub fn held(f:Func,body:&Expr,mut limits:Vec<Vec<Expr>>)->Result<Expr,String> {
    let mut body=body.clone();
    if let Kind::Func(inner,args)=body.kind() {
        if *inner==f {
            let mut previous=self::limits(args);previous.append(&mut limits);limits=previous;
            body=args[0].clone();
        }
    }
    let mut args=vec![body];
    args.extend(limits.into_iter().map(crate::special::tuple));
    validate(&f,&args)?;
    Ok(Expr::from_kind(Kind::Func(f,args)))
}

pub fn validate_tree(e:&Expr)->Result<(),String> {
    // args() can synthesize temporaries: retain their owners against ID reuse.
    fn walk(e:&Expr,seen:&mut FxHashMap<usize,Expr>)->Result<(),String> {
        if !e.has_calculus_binding() || seen.insert(e.dag_id(),e.clone()).is_some() { return Ok(()); }
        if let Kind::Func(f,args)=e.kind() { if is_calculus(f) { validate(f,args)?; } }
        if let Kind::Func(Func::MatrixDotProduct(a,b),_)=e.kind(){for parent in [a,b]{let mut fields=Vec::new();parent.scalar_fields(&mut fields);for scalar in fields{walk(&scalar,seen)?;}}}
        if let Kind::Func(Func::MatrixElement(parent)|Func::MatrixProxy(parent)|Func::MatrixTrace(parent)|Func::MatrixDeterminant(parent),_)=e.kind(){
            let mut fields=Vec::new();parent.scalar_fields(&mut fields);
            for scalar in fields{walk(&scalar,seen)?;}
        }
        for arg in crate::inspection::scalar_children(e) { walk(&arg,seen)?; }
        Ok(())
    }
    walk(e,&mut FxHashMap::default())
}

/// Each subsequent (outer) limit binds the body and earlier limits, but
/// never its own endpoints. An indefinite integral introduces a free variable.
pub fn free_symbols(f:&Func,args:&[Expr])->Vec<Expr> {
    // ExprWithLimits.free_symbols masks a non-symbol variable such as f(x)
    // with a Dummy, so the symbols inside a bound f(x) are not free.
    let masks:FxHashMap<Expr,Expr>=limits(args).into_iter()
        .filter(|limit|limit.len()>1 && !limit[0].is_symbol())
        .map(|limit|(limit[0].clone(),anonymous_dummy())).collect();
    let body=if masks.is_empty(){args[0].clone()}else{crate::subs::xreplace(&args[0],&masks)};
    let mut symbols:FxHashSet<Expr>=crate::traverse::free_symbols(&body).into_iter().collect();
    for limit in limits(args) {
        if *f==Func::Integral && limit.len()==1 {
            if limit[0].is_symbol() {symbols.insert(limit[0].clone());}
            else {symbols.extend(crate::traverse::free_symbols(&limit[0]));}
        }
        else {
            symbols.remove(masks.get(&limit[0]).unwrap_or(&limit[0]));
            for endpoint in &limit[1..] { symbols.extend(crate::traverse::free_symbols(endpoint)); }
        }
    }
    let mut result:Vec<_>=symbols.into_iter().collect();result.sort_by(Expr::cmp_struct);result
}

/// Symbol substitutions respect definite scopes. Indefinite variables are
/// free: a symbolic replacement renames them; another expression evaluates
/// the formal antiderivative at an endpoint without replacing its integrand.
pub fn substitute(e:&Expr,pairs:&[(Expr,Expr)],simultaneous:bool)->Expr {
    if pairs.iter().any(|(old,new)|check_substitution(e,old,new).is_err()){return e.clone();}
    if pairs.iter().all(|(old,new)|old==new) {return e.clone();}
    if !simultaneous && pairs.iter().any(|(old,_)|!old.is_symbol()) {
        let mut result=e.clone();
        for (old,new) in pairs {
            if old.is_symbol(){result=substitute(&result,&[(old.clone(),new.clone())],false);continue;}
            let Kind::Func(f,args)=result.kind() else {return crate::subs::subs1(&result,old,new);};
            let old_symbols=crate::traverse::free_symbols(old);
            let mut changed=limits(args);let mut replace_body=true;
            for limit in changed.iter_mut().rev() {
                // SymPy: an indefinite variable equal to old is renamed when
                // new can be differentiated with respect to, else evaluated.
                if *f==Func::Integral && limit.len()==1 && limit[0]==*old {
                    if renames_variable(new) {limit[0]=new.clone();} else {limit.push(new.clone());}
                } else {
                    for endpoint in &mut limit[1..] {*endpoint=crate::subs::subs1(endpoint,old,new);}
                }
                if crate::traverse::free_symbols(&limit[0]).iter().any(|s|old_symbols.contains(s)) {replace_body=false;break;}
            }
            if matches!(old.kind(),Kind::Func(Func::Undefined(head),_) if !head.internal) {
                if check_substitution(&result,old,new).is_err(){return result;}
                replace_body=true;
            }
            let body=if replace_body {crate::subs::subs1(&args[0],old,new)} else {args[0].clone()};
            result=held(f.clone(),&body,changed).unwrap_or(result);
        }
        return result;
    }
    let Kind::Func(f,args)=e.kind() else { return e.clone(); };
    let old_limits=limits(args);
    let transform=|value:&Expr,active:&[(Expr,Expr)]| if simultaneous {crate::subs::simultaneous(value,active)} else {crate::subs::subs(value,active)};
    let mut body_pairs=pairs.to_vec();
    for limit in &old_limits {
        if limit.len()>1 || *f!=Func::Integral {
            body_pairs.retain(|(key,_)|*key!=limit[0]);
        } else {
            body_pairs.retain(|(key,value)|*key!=limit[0] || renames_variable(value));
        }
        if !limit[0].is_symbol() {
            // A symbol inside a bound f(x) is bound with it (SymPy keeps the body).
            let bound=crate::traverse::free_symbols(&limit[0]);
            body_pairs.retain(|(key,_)|!bound.contains(key));
        }
    }
    let body=transform(&args[0],&body_pairs);
    let mut new_limits=Vec::new();
    for (index,limit) in old_limits.iter().enumerate() {
        let active:Vec<_>=pairs.iter().filter(|(key,_)|!old_limits[index+1..].iter().any(|outer|outer[0]==*key && (outer.len()>1 || (*f==Func::Integral && pairs.iter().any(|(old,value)|old==key && !value.is_symbol()))))).cloned().collect();
        let mut changed=vec![limit[0].clone()];
        if *f==Func::Integral && limit.len()==1 {
            if let Some((_,value))=active.iter().find(|(key,_)|*key==limit[0]) {
                if renames_variable(value) { changed[0]=value.clone(); } else { changed.push(value.clone()); }
            }
        } else { changed.extend(limit[1..].iter().map(|endpoint|transform(endpoint,&active))); }
        new_limits.push(changed);
    }
    held(f.clone(),&body,new_limits).unwrap_or_else(|_|e.clone())
}

pub fn doit(f:&Func,args:&[Expr],deep:bool)->Expr {
    if *f==Func::Product {return crate::products::doit(args,deep);}
    let mut value=if deep {crate::functions::doit(&args[0],true)} else {args[0].clone()};
    for mut limit in limits(args) {
        if deep { for endpoint in &mut limit[1..] {*endpoint=crate::functions::doit(endpoint,true);} }
        if !limit[0].is_symbol() {
            // Integrate with respect to f(x) through a fresh symbol, as SymPy does.
            let variable=limit[0].clone();let dummy=anonymous_dummy();
            let body=crate::subs::subs1(&value,&variable,&dummy);limit[0]=dummy.clone();
            let Ok(renamed)=held(f.clone(),&body,vec![limit]) else {return value;};
            value=crate::subs::subs1(&doit(f,&renamed.args(),false),&dummy,&variable);
            continue;
        }
        value=if *f==Func::Sum {
            crate::summation::summation(&value,&limit[0],&limit[1],&limit[2])
        } else { match limit.len() {
            1=>crate::integrate::integrate(&value,&limit[0]),
            2=>match crate::integrate::try_integrate(&value,&limit[0]) {
                Some(antiderivative)=>crate::subs::subs1(&antiderivative,&limit[0],&limit[1]),
                // the undone limit is factored as `Integral.doit` does
                None=>crate::integrate::factor_undone(&held(Func::Integral,&value,vec![limit]).unwrap()),
            },
            _=>crate::integrate::integrate_definite(&value,&limit[0],&limit[1],&limit[2]),
        }};
    }
    value
}

/// FTC/Leibniz use a fresh integration atom so a bound x in the body cannot
/// be confused with a free x in an endpoint. Unknown discrete endpoint
/// derivatives remain held derivatives, rather than a continuous approximation.
pub fn differentiate(e:&Expr,x:&Expr)->Expr {
    if is_function_variable(x) {
        // Differentiate with respect to f(x) as with respect to a symbol.
        let dummy=anonymous_dummy();
        let renamed=crate::subs::subs1(e,x,&dummy);
        if renamed!=*e {return crate::subs::subs1(&differentiate(&renamed,&dummy),&dummy,x);}
    }
    if matches!(x.kind(),Kind::Func(Func::MatrixElement(_),_)) {
        // Sum._eval_derivative by a matrix element: Sum(f.diff(x), limit). The
        // element may match a summand entry with other indices, so no
        // free-symbol shortcut applies.
        if let Kind::Func(Func::Sum,args)=e.kind() {
            let mut limits=limits(args);let outer=limits.pop().unwrap();
            let inner=if limits.is_empty() {args[0].clone()} else {held(Func::Sum,&args[0],limits).unwrap()};
            if outer[1..].iter().any(|endpoint|crate::traverse::has(endpoint,x)) {
                return crate::diff::derivative_unevaluated(e,vec![(x.clone(),1)]);
            }
            let partial=crate::diff::diff(&inner,x);
            return if partial.is_zero() {partial} else {held(Func::Sum,&partial,vec![outer]).unwrap()};
        }
    }
    if !x.is_symbol(){return crate::diff::derivative_unevaluated(e,vec![(x.clone(),1)]);}
    let Kind::Func(f,args)=e.kind() else {return Expr::zero();};
    if !free_symbols(f,args).contains(x) { return Expr::zero(); }
    let mut limits=limits(args);let outer=limits.pop().unwrap();
    let inner=if limits.is_empty() {args[0].clone()} else {held(f.clone(),&args[0],limits).unwrap()};
    if *f==Func::Sum {
        if outer[1..].iter().any(|endpoint|crate::traverse::free_symbols(endpoint).contains(x)) {
            return crate::diff::derivative_unevaluated(e,vec![(x.clone(),1)]);
        }
        let partial=crate::diff::diff(&inner,x);
        return if partial.is_zero() {partial} else {held(Func::Sum,&partial,vec![outer]).unwrap()};
    }
    if *f==Func::Product {
        // SymPy's product rule: sum over h of the factors before and after h
        // times the derivative of the h-th factor, with an anonymous Dummy h.
        if outer[1..].iter().any(|endpoint|crate::traverse::free_symbols(endpoint).contains(x)) {
            return crate::diff::derivative_unevaluated(e,vec![(x.clone(),1)]);
        }
        let h=Expr::from_kind(Kind::Sym(crate::symbol::Symbol::anonymous_dummy(crate::symbol::Assumptions::NONE)));
        let one=Expr::one();
        let before=held(Func::Product,&inner,vec![vec![outer[0].clone(),outer[1].clone(),crate::ops::sub(&h,&one)]]).unwrap();
        let after=held(Func::Product,&inner,vec![vec![outer[0].clone(),crate::ops::add2(&h,&one),outer[2].clone()]]).unwrap();
        let partial=crate::subs::subs1(&crate::diff::diff(&inner,x),&outer[0],&h);
        let term=crate::ops::mul(&[before,after,partial]);
        return held(Func::Sum,&term,vec![vec![h,outer[1].clone(),outer[2].clone()]]).unwrap();
    }
    if outer.len()==1 {
        if outer[0]==*x {return inner;}
        let partial=if outer[0].is_symbol() {crate::diff::diff(&inner,x)} else {
            // SymPy masks a bound f(x) while differentiating its integrand.
            let dummy=anonymous_dummy();
            crate::subs::subs1(&crate::diff::diff(&crate::subs::subs1(&inner,&outer[0],&dummy),x),&dummy,&outer[0])
        };
        return if partial.is_zero() {partial} else {held(Func::Integral,&partial,vec![outer]).unwrap()};
    }
    if inner.has_calculus_binding() && outer[1..].iter().any(|point|matches!(point.kind(),Kind::Func(Func::Undefined(head),_) if !head.internal)) {
        return crate::diff::derivative_unevaluated(e,vec![(x.clone(),1)]);
    }
    let dummy=match outer[0].as_sym() {Some(symbol)=>Expr::from_kind(Kind::Sym(symbol.fresh_dummy_like("integration"))),None=>anonymous_dummy()};
    let renamed=crate::subs::subs1(&inner,&outer[0],&dummy);
    let partial=crate::diff::diff(&renamed,x);
    let mut terms=Vec::new();
    if !partial.is_zero() {
        let partial=crate::subs::subs1(&partial,&dummy,&outer[0]);
        terms.push(held(Func::Integral,&partial,vec![outer.clone()]).unwrap());
    }
    let upper=outer.last().unwrap();let d_upper=crate::diff::diff(upper,x);
    if !d_upper.is_zero() {terms.push(mul2(&crate::subs::subs1(&renamed,&dummy,upper),&d_upper));}
    if outer.len()==3 {
        let d_lower=crate::diff::diff(&outer[1],x);
        if !d_lower.is_zero() {terms.push(neg(&mul2(&crate::subs::subs1(&renamed,&dummy,&outer[1]),&d_lower)));}
    }
    add(&terms)
}

/// Reject unsupported variable shapes and substitutions which would introduce
/// a new bound dependency through an applied undefined function replacement.
pub fn check_substitution(e:&Expr,old:&Expr,new:&Expr)->Result<(),String> {
    // The physical memo must own every transient public argument it has seen.
    fn walk(e:&Expr,old:&Expr,new:&Expr,seen:&mut FxHashMap<usize,Expr>)->Result<(),String>{
        if !e.has_calculus_binding() || e==old || seen.insert(e.dag_id(),e.clone()).is_some(){return Ok(());}
        if let Kind::Func(f,args)=e.kind(){if is_calculus(f){
            let limits=limits(args);
            if let Kind::Func(Func::Undefined(head),old_args)=old.kind(){if !head.internal {
                if limits.iter().any(|limit|crate::traverse::has(new,&limit[0]) && !old_args.contains(&limit[0])) {
                    return Err("substitution cannot create dummy dependencies".into());
                }
            }}
        }}
        if let Kind::Func(Func::MatrixDotProduct(a,b),_)=e.kind(){for parent in [a,b]{let mut fields=Vec::new();parent.scalar_fields(&mut fields);for scalar in fields{walk(&scalar,old,new,seen)?;}}}
        if let Kind::Func(Func::MatrixElement(parent)|Func::MatrixProxy(parent)|Func::MatrixTrace(parent)|Func::MatrixDeterminant(parent),_)=e.kind(){
            let mut fields=Vec::new();parent.scalar_fields(&mut fields);
            for scalar in fields{walk(&scalar,old,new,seen)?;}
        }
        for arg in crate::inspection::scalar_children(e){walk(&arg,old,new,seen)?;}Ok(())
    }
    walk(e,old,new,&mut FxHashMap::default())
}

/// Bound-variable presence is structural; independence is about free symbols.
pub fn depends_on(e:&Expr,x:&Expr)->bool {
    if e.has_calculus_binding() && x.is_symbol() {crate::traverse::free_symbols(e).contains(x)} else {crate::traverse::has(e,x)}
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::sub;
    fn p(text:&str)->Expr {crate::parser::parse_str(text).unwrap()}
    #[test]
    fn held_limits_free_symbols_and_roundtrips() {
        for (source,free) in [("Integral(1,x)",vec![p("x")]),("Integral(x,(x,0,x))",vec![p("x")]),("Integral(x+y,(x,0,y),(y,0,z))",vec![p("z")]),("Sum(x*y,(x,1,n))",vec![p("n"),p("y")])] {
            let e=p(source);let mut actual=crate::traverse::free_symbols(&e);actual.sort();let mut expected=free;expected.sort();assert_eq!(actual,expected);
            assert_eq!(crate::parser::parse_str(&crate::printer::srepr(&e)).unwrap(),e);
            assert_eq!(crate::parser::parse_str(&e.to_string()).unwrap(),e);
        }
        assert_ne!(p("Integral(x,x)"),p("Function('Integral')(x,x)"));
    }
    #[test]
    fn held_evaluation_and_karr_convention() {
        for (source,expected) in [("Integral(x,x)","x**2/2"),("Integral(x,(x,2))","2"),("Integral(x+y,(x,0,y),(y,0,z))","z**3/2"),("Sum(k,(k,2,-1))","-1"),("Sum(sin(k),(k,3,0))","-sin(1)-sin(2)"),("Sum(k,(k,1,0))","0")] {
            let actual=crate::functions::doit(&p(source),true);
            assert!(crate::simplify::simplify(&sub(&actual,&p(expected))).is_zero(),"{source}: {actual}");
        }
    }
    #[test]
    fn validators_preserve_owners_of_synthesized_arguments() {
        // Add.args() synthesizes temporary coefficient Mul nodes. Repeated
        // wrappers and shared powers must not hide a different later branch.
        let valid=p("Integral(f(y),(z,0,1))");
        let forbidden=p("Integral(f(y),(x,0,1))");
        let invalid=Expr::from_kind(Kind::Func(Func::Integral,vec![p("y"),crate::special::tuple(vec![Expr::int(2)])]));
        let branches=|last:Expr|Expr::from_kind(Kind::Add(crate::Number::ZERO,vec![
            (valid.powi(2),crate::Number::Int(2)),
            (valid.powi(3),crate::Number::Int(3)),
            (last.powi(2),crate::Number::Int(4)),
        ]));
        for _ in 0..64 {
            assert!(validate_tree(&branches(invalid.clone())).is_err());
            assert!(check_substitution(&branches(forbidden.clone()),&p("f(y)"),&p("x")).is_err());
        }
    }

    #[test]
    fn composite_substitution_respects_the_first_bound_dependency() {
        let source=p("Integral(f(x)+y,(x,0,y))");
        assert_eq!(crate::subs::subs1(&source,&p("f(x)+y"),&p("z")),source);
        assert_eq!(crate::subs::subs1(&source,&p("f(x)"),&p("z")),p("Integral(z+y,(x,0,y))"));
        assert!(check_substitution(&p("Integral(f(y),(x,0,1))"),&p("f(y)"),&p("x")).is_err());
        assert_eq!(crate::functions::doit(&p("Sum(1,(x,oo,oo))"),true),p("Sum(1,(x,oo,oo))"));
    }
    #[test]
    fn held_scopes_substitution_and_safe_leibniz() {
        let x=p("x");let y=p("y");
        for (source,old,new,expected) in [
            ("Integral(f(x),(x,0,y))","x","z","Integral(f(x),(x,0,y))"),
            ("Integral(f(x),(x,0,x))","x","z","Integral(f(x),(x,0,z))"),
            ("Integral(f(x),x)","x","2","Integral(f(x),(x,2))"),
            ("Integral(f(x),x,x)","x","2","Integral(f(x),x,(x,2))"),
            ("Integral(x+y,(x,0,y),(y,0,z))","y","2","Integral(x+y,(x,0,y),(y,0,z))"),
        ] {assert_eq!(crate::subs::subs1(&p(source),&p(old),&p(new)),p(expected));}
        assert_eq!(crate::subs::simultaneous(&p("Integral(x+y,(x,0,y))"),&[(x.clone(),y.clone()),(y.clone(),x.clone())]),p("Integral(2*x,(x,0,x))"));
        for (source,variable,expected) in [("Integral(f(x),x)","x","f(x)"),("Integral(f(x),(x,0,x))","x","f(x)"),("Integral(y*f(x),(x,0,y))","y","y*f(y)+Integral(f(x),(x,0,y))"),("Sum(k*x,(k,1,n))","x","Sum(k,(k,1,n))")] {
            assert_eq!(crate::diff::diff(&p(source),&p(variable)),p(expected),"{source}");
        }
        assert_eq!(crate::diff::diff(&p("Sum(k,(k,1,x))"),&x),p("Derivative(Sum(k,(k,1,x)),x)"));
    }
}
