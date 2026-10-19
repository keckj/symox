//! Algebraic-atom solving for one non-Symbol target, with exact reconstruction.
use rustc_hash::FxHashMap;

use crate::{Expr,Func,Kind};
use crate::inspection::{classify,scalar_nodes,PolynomialKind};
use crate::ops::{div,mul2,sub};
use crate::symbol::Assumptions;

fn exponential_multiple(argument:&Expr,target_argument:&Expr)->Option<(Expr,Expr)>{
    if target_argument.is_symbol(){
        if let Some(coefficients)=crate::integrate::poly_coeffs(argument,target_argument){
            if coefficients.len()==2&&coefficients[1].is_integer_number()&&!coefficients[1].is_zero(){
                return Some((coefficients[1].clone(),coefficients[0].clone()));
            }
        }
    }
    let multiple=crate::poly::api::cancel(&div(argument,target_argument)).ok()?;
    if multiple.is_integer_number()&&!multiple.is_zero()
        &&crate::expand::expand(&sub(argument,&mul2(&multiple,target_argument))).is_zero(){
        Some((multiple,Expr::zero()))
    }else{None}
}

struct Recast{equations:Vec<Expr>,variable:Expr,restore:FxHashMap<Expr,Expr>}

fn prepare(equations:&[Expr],target:&Expr)->Option<Recast>{
    if target.is_symbol()||target.has_matrix()||equations.iter().any(Expr::has_matrix)||crate::functions::is_numeric_value(target)
        ||crate::assumptions::fact(target,"commutative")==Some(false){return None;}
    if matches!(target.kind(),Kind::Num(_)|Kind::Const(_)|Kind::Bool(_)|Kind::Rel(..)
        |Kind::Func(Func::Set(_)|Func::Logic(_)|Func::AppliedPredicate|Func::Predicate(_)|Func::Order,_)){
        return None;
    }
    let variable=Expr::dummy("_solve_target",Assumptions::NONE);
    // Like reference solve, these are algebraic atoms: do not transfer the
    // analytic range of exp(x), a real-valued f(x), etc. to the fresh Dummy.
    let mut restore:FxHashMap<Expr,Expr>=Default::default();
    restore.insert(variable.clone(),target.clone());
    let opaque_derivatives=matches!(target.kind(),Kind::Func(Func::Undefined(function),_) if !function.internal);
    let mut opaque:FxHashMap<Expr,Expr>=Default::default();
    let originals:Vec<_>=equations.iter().map(super::as_expr).collect();
    if opaque_derivatives{
        for equation in &originals{
            for (node,_) in scalar_nodes(equation){
                let binder=matches!(node.kind(),Kind::Derivative(..)|Kind::Func(
                    Func::Integral|Func::Sum|Func::Product|Func::Limit|Func::Subs,_));
                if binder&&node!=*target&&!opaque.contains_key(&node)
                    &&scalar_nodes(&node).iter().any(|(child,_)|child==target){
                    let dummy=Expr::dummy("_solve_coefficient",Assumptions::NONE);
                    restore.insert(dummy.clone(),node.clone());opaque.insert(node,dummy);
                }
            }
        }
    }
    let mut converted=Vec::with_capacity(equations.len());let mut occurs=false;
    for original in &originals{
        let masked=crate::subs::xreplace(original,&opaque);
        let mut working=crate::subs::subs1(&masked,target,&variable);
        if let Kind::Func(Func::Exp,arguments)=target.kind(){
            let mut powers:FxHashMap<Expr,Expr>=Default::default();
            for (node,_) in scalar_nodes(&working){
                if let Kind::Func(Func::Exp,values)=node.kind(){
                    if let Some((multiple,offset))=exponential_multiple(&values[0],&arguments[0]){
                        let power=crate::ops::pow(&variable,&multiple);
                        powers.insert(node.clone(),mul2(&crate::functions::exp(&offset),&power));
                    }
                }
            }
            working=crate::subs::xreplace(&working,&powers);
        }
        if classify(&working,&[variable.clone()],PolynomialKind::Rational)!=Some(true){return None;}
        // Algebraic substitution may attempt branch-sensitive transforms.
        // Accept the recast only when restoring every atom reconstructs the
        // original exactly, possibly after ordinary ring expansion. No float
        // residual or forced power/log identity is a reconstruction proof.
        let reconstructed=crate::subs::xreplace(&working,&restore);
        if reconstructed!=*original&&!crate::expand::expand(&sub(&reconstructed,original)).is_zero(){return None;}
        occurs|=crate::traverse::has(&working,&variable);
        converted.push(working);
    }
    occurs.then_some(Recast{equations:converted,variable,restore})
}

pub fn solve_atom(equation:&Expr,target:&Expr,check:bool)->Option<Vec<Expr>>{
    let recast=prepare(&[equation.clone()],target)?;
    let mut values=super::solve_with_check(&recast.equations[0],&recast.variable,check)?
        .iter().map(|value|crate::subs::xreplace(value,&recast.restore)).collect();
    super::sort_solutions(&mut values);Some(values)
}

pub fn solve_atom_system(equations:&[Expr],target:&Expr,check:bool)->Option<super::SystemSolution>{
    let recast=prepare(equations,target)?;
    Some(match super::solve_system_with_check(&recast.equations,&[recast.variable],check)?{
        super::SystemSolution::NoSolution=>super::SystemSolution::NoSolution,
        super::SystemSolution::Map(values)=>super::SystemSolution::Map(values.into_iter()
            .map(|(key,value)|(crate::subs::xreplace(&key,&recast.restore),crate::subs::xreplace(&value,&recast.restore))).collect()),
        super::SystemSolution::Tuples(rows)=>super::SystemSolution::Tuples(rows.into_iter()
            .map(|row|row.iter().map(|value|crate::subs::xreplace(value,&recast.restore)).collect()).collect()),
    })
}

#[cfg(test)]
mod tests{
    use super::*;
    fn p(source:&str)->Expr{crate::parser::parse_str(source).unwrap()}

    #[test]
    fn rational_recast_requires_exact_reconstruction(){
        let equation=p("exp(2*x)+x");let target=p("exp(x)");
        let recast=prepare(&[equation.clone()],&target).unwrap();
        let reconstructed=crate::subs::xreplace(&recast.equations[0],&recast.restore);
        assert!(crate::expand::expand(&sub(&reconstructed,&equation)).is_zero());
        assert_eq!(recast.restore[&recast.variable],target);
    }

    #[test]
    fn noninteger_exponential_lifting_is_not_assumed(){
        // exp(x/2) is not the principal square root of exp(x) on all sheets.
        assert!(prepare(&[p("exp(x/2)-1")],&p("exp(x)")).is_none());
        // A fractional power in the new atom is outside the rational scope.
        assert!(prepare(&[p("x**3-1")],&p("x**2")).is_none());
    }
}
