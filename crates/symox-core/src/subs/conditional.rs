//! Substitution of predicates whose parameter domain includes nonfinite scalars.
//!
//! Nested-exponential inversion retains its original equation under
//! `~Contains(parameter, Complexes)`. Until that guard is specialized, ordinary
//! finite-symbol facts (for example exp(z) != 0) cannot discard the equation.
//! This context belongs to the conditional predicate, not to ordinary Eq.subs.

use crate::expr::{Expr, Func, Kind};
use crate::logic::LogicOp;
use crate::sets::SetOp;
use rustc_hash::FxHashSet;

fn any_node(e: &Expr, predicate: &dyn Fn(&Expr) -> bool) -> bool {
    crate::inspection::scalar_nodes(e).iter().any(|(node,_)| predicate(node))
}

// Exact structural substitution without evaluating parameter-dependent
// arithmetic. The public substitution callback still handles atoms and host
// protocols; callers use this path only for symbol subs or exact xreplace.
fn raw(e: &Expr, transform: &mut dyn FnMut(&Expr) -> Expr,
       replacement: &dyn Fn(&Expr) -> Option<Expr>) -> Expr {
    if let Some(value) = replacement(e) { return value; }
    if matches!(e.kind(),Kind::Func(Func::SymbolicDerivative,_)){return transform(e);}
    if e.is_atom() || e.host_handle().is_some_and(|host|host.substitutes()) { return transform(e); }
    let children = crate::inspection::scalar_children(e);
    let mapped = children.iter().map(|child| raw(child, transform, replacement)).collect::<Vec<_>>();
    if mapped.iter().zip(&children).all(|(a,b)| a.same_presentation(b)) { return e.clone(); }
    let kind = match e.kind() {
        Kind::Add(..) => Kind::Add(crate::number::Number::ZERO,
            mapped.into_iter().map(|child| (child, crate::number::Number::ONE)).collect()),
        Kind::Mul(..) => Kind::Mul(crate::number::Number::ONE,
            mapped.into_iter().map(|child| (child, Expr::one())).collect()),
        Kind::Pow(..) => Kind::Pow(mapped[0].clone(), mapped[1].clone()),
        Kind::Func(fun, _) => Kind::Func(fun.clone(), mapped),
        Kind::Rel(op, _, _) => Kind::Rel(*op, mapped[0].clone(), mapped[1].clone()),
        _ => return transform(e),
    };
    Expr::from_unevaluated_kind(kind)
}

// A zero times a factor containing infinity need not be zero: for
// example 0*exp(oo*x) is undefined for x>0 and zero for x<0. Ordinary
// arithmetic currently erases that distinction; this scope must retain it.
fn uncertain_zero_product(e: &Expr) -> bool {
    any_node(e,&|node| {
        if !node.is_mul() { return false; }
        let factors=crate::inspection::scalar_children(node);
        factors.iter().any(Expr::is_zero) && factors.iter().any(|factor| factor.has_nonfinite()
            && !matches!(factor.kind(),Kind::Const(_)))
    })
}

fn nonfinite_variable_calculation(e: &Expr) -> bool {
    e.has_nonfinite() && !crate::traverse::free_symbols(e).is_empty()
}

fn uncertain_arithmetic(e: &Expr) -> bool {
    uncertain_zero_product(e) || any_node(e,&|node| {
        if !node.is_add() && !node.is_mul() {return false;}
        let children=crate::inspection::scalar_children(node);
        children.iter().filter(|child|child.has_nonfinite()).count()>1
            && children.iter().any(nonfinite_variable_calculation)
    })
}

// Canonicalize arithmetic introduced by earlier partial substitutions,
// while retaining every operation still depending on a guarded parameter.
fn settled(e: &Expr, parameters: &FxHashSet<Expr>) -> Expr {
    if matches!(e.kind(),Kind::Func(Func::SymbolicDerivative,_)){return e.clone();}
    if e.is_atom() { return e.clone(); }
    let children=crate::inspection::scalar_children(e);
    let mapped=children.iter().map(|child|settled(child,parameters)).collect::<Vec<_>>();
    let pending=parameters.iter().any(|p|crate::traverse::has(e,p))
        || mapped.iter().any(uncertain_arithmetic)
        || ((e.is_add()||e.is_mul()) && mapped.iter().filter(|child|child.has_nonfinite()).count()>1
            && mapped.iter().any(nonfinite_variable_calculation))
        || (e.is_mul() && mapped.iter().any(Expr::is_zero)
            && mapped.iter().any(|factor|factor.has_nonfinite()&&!matches!(factor.kind(),Kind::Const(_))));
    match e.kind() {
        Kind::Add(..) if !pending=>crate::ops::add(&mapped),
        Kind::Mul(..) if !pending=>crate::ops::mul(&mapped),
        Kind::Pow(..) if !pending=>crate::ops::pow(&mapped[0],&mapped[1]),
        Kind::Func(Func::Exp,_) if !pending=>crate::functions::exp(&mapped[0]),
        _ if mapped.iter().zip(&children).all(|(a,b)|a.same_presentation(b))=>e.clone(),
        Kind::Add(..)=>Expr::from_unevaluated_kind(Kind::Add(crate::number::Number::ZERO,
            mapped.into_iter().map(|e|(e,crate::number::Number::ONE)).collect())),
        Kind::Mul(..)=>Expr::from_unevaluated_kind(Kind::Mul(crate::number::Number::ONE,
            mapped.into_iter().map(|e|(e,Expr::one())).collect())),
        Kind::Pow(..)=>Expr::from_unevaluated_kind(Kind::Pow(mapped[0].clone(),mapped[1].clone())),
        Kind::Func(fun,_) if pending=>Expr::from_unevaluated_kind(Kind::Func(fun.clone(),mapped)),
        Kind::Func(fun,_)=>crate::functions::func(fun.clone(),mapped),
        _=>e.clone(),
    }
}

fn negative_memberships(e: &Expr, negative: bool, out: &mut FxHashSet<Expr>) -> bool {
    match e.kind() {
        Kind::Func(Func::Logic(LogicOp::Not), args) if args.len() == 1 =>
            negative_memberships(&args[0], !negative, out),
        Kind::Func(Func::Logic(LogicOp::And | LogicOp::Or), args) => {
            let mut found = false;
            for child in args { found |= negative_memberships(child, negative, out); }
            found
        }
        Kind::Func(Func::Set(SetOp::Contains), args) if negative && args.len() == 2
            && args[1] == crate::sets::complexes() => {
            let free=crate::traverse::free_symbols(&args[0]).into_iter()
                .filter(|p| crate::assumptions::is_finite(p) != Some(true)).collect::<FxHashSet<_>>();
            out.extend(free.iter().cloned());
            let scalar=settled(&args[0],&free);
            if crate::assumptions::is_finite(&scalar).is_none()
                && !matches!(scalar.kind(),Kind::Const(crate::expr::Constant::NaN)) {
                out.insert(args[0].clone());
            }
            true
        }
        _ => false,
    }
}

// A nonfinite-domain predicate must be a conjunct of the equation. A
// disjunctive alternative is not that scope and keeps its ordinary protocol.
fn guard_parameters(condition: &Expr, out: &mut FxHashSet<Expr>) -> bool {
    let Kind::Func(Func::Logic(LogicOp::And),args)=condition.kind() else {return false;};
    let mut found=false;
    for arg in args {
        if matches!(arg.kind(),Kind::Func(Func::Logic(LogicOp::And),_)) {
            found|=guard_parameters(arg,out);
        } else if !any_node(arg,&|e|matches!(e.kind(),Kind::Rel(..))) {
            found|=negative_memberships(arg,false,out);
        }
    }
    found
}

pub(super) fn substitute(
    condition: &Expr,
    transform: &mut dyn FnMut(&Expr) -> Expr,
    replacement: &dyn Fn(&Expr) -> Option<Expr>,
) -> Option<Expr> {
    crate::assumptions::exact_queries(|| substitute_exact(condition,transform,replacement))
}

fn substitute_exact(
    condition: &Expr,
    transform: &mut dyn FnMut(&Expr) -> Expr,
    replacement: &dyn Fn(&Expr) -> Option<Expr>,
) -> Option<Expr> {
    let mut ignored = FxHashSet::default();
    if !guard_parameters(condition, &mut ignored)
        || !any_node(condition, &|e| matches!(e.kind(), Kind::Rel(..))) {
        return None;
    }
    let mapped = raw(condition, transform, replacement);
    let mut unresolved = FxHashSet::default();
    guard_parameters(&mapped, &mut unresolved);
    let mut protected = false;
    let mut keep_guard = false;
    let mut relation_values = rustc_hash::FxHashMap::default();
    // Relation nodes retain their original tree position in raw(). Walking
    // original/raw pairs lets complete operands use ordinary native subs,
    // including its precision and host behavior, once their guards decide.
    fn relations(e: &Expr, transform: &mut dyn FnMut(&Expr) -> Expr,
        replacement: &dyn Fn(&Expr) -> Option<Expr>, parameters: &FxHashSet<Expr>,
        protected: &mut bool, keep_guard: &mut bool, out: &mut rustc_hash::FxHashMap<Expr, Expr>) {
        if replacement(e).is_some() { return; }
        match e.kind() {
            Kind::Rel(op, lhs, rhs) => {
                let mapped = raw(e, transform, replacement);
                let Kind::Rel(_, left, right) = mapped.kind() else { return; };
                let pending = parameters.iter().any(|p| crate::traverse::has(left,p) || crate::traverse::has(right,p));
                let (left, right) = if pending || nonfinite_variable_calculation(left) || nonfinite_variable_calculation(right) {
                    (settled(left,parameters),settled(right,parameters)) }
                    else { (settled(&transform(lhs),parameters),settled(&transform(rhs),parameters)) };
                // NaN records an undefined specialization, not a proof that
                // its original equation has no solutions. Keep the calculation.
                let undefined = any_node(&left,&|e| *e==Expr::nan())
                    || any_node(&right,&|e| *e==Expr::nan())
                    || nonfinite_variable_calculation(&left) || nonfinite_variable_calculation(&right)
                    || uncertain_arithmetic(&left) || uncertain_arithmetic(&right);
                let result = if pending || undefined {
                    *protected = true;
                    *keep_guard |= pending || !crate::traverse::free_symbols(&left).is_empty()
                        || !crate::traverse::free_symbols(&right).is_empty();
                    Expr::from_unevaluated_kind(Kind::Rel(*op,left,right))
                } else { crate::relational::rel(*op,&left,&right) };
                out.insert(mapped,result);
            }
            Kind::Func(Func::Logic(_),args) => {
                for child in args { relations(child,transform,replacement,parameters,protected,keep_guard,out); }
            }
            _ => {}
        }
    }
    relations(condition,transform,replacement,&unresolved,&mut protected,&mut keep_guard,&mut relation_values);
    fn guard(e: &Expr) -> Expr {
        match e.kind() {
            Kind::Func(Func::Set(SetOp::Contains),args) if args.len()==2 =>
                crate::sets::contains_expr(&args[0],&args[1]),
            Kind::Func(Func::Logic(op),args) => crate::functions::func(Func::Logic(*op),args.iter().map(guard).collect()),
            _ => e.clone(),
        }
    }
    fn rebuild(e: &Expr, protected: bool, keep_guard: bool, relations: &rustc_hash::FxHashMap<Expr,Expr>) -> Expr {
        if let Some(value)=relations.get(e) { return value.clone(); }
        match e.kind() {
            Kind::Func(Func::Logic(LogicOp::Not),_) if {
                let mut parameters=FxHashSet::default();
                negative_memberships(e,false,&mut parameters)
            } => if keep_guard { e.clone() } else { guard(e) },
            Kind::Func(Func::Logic(op),args) => {
                let mapped=args.iter().map(|child|rebuild(child,protected,keep_guard,relations)).collect::<Vec<_>>();
                if !protected { return crate::functions::func(Func::Logic(*op),mapped); }
                // A true nonfinite guard can still carry unresolved parameters.
                // Do not erase that provenance or canonicalize complementary
                // relations with the same finite-symbol shortcut.
                if matches!(op,LogicOp::And|LogicOp::Or) {
                    let absorbing=*op==LogicOp::Or;
                    let mut children=Vec::new();
                    for child in mapped {
                        if let Kind::Bool(value)=child.kind() {
                            if *value==absorbing { return Expr::boolean(absorbing); }
                        } else if !children.contains(&child) { children.push(child); }
                    }
                    crate::logic::ordered(&mut children);
                    match children.len() {
                        0=>Expr::boolean(!absorbing),
                        1=>children.pop().unwrap(),
                        _=>Expr::from_unevaluated_kind(Kind::Func(Func::Logic(*op),children)),
                    }
                } else { Expr::from_unevaluated_kind(Kind::Func(Func::Logic(*op),mapped)) }
            }
            _=>e.clone(),
        }
    }
    Some(rebuild(&mapped,protected,keep_guard,&relation_values))
}

/// Piecewise's normal relation-coefficient canonicalization would evaluate
/// the protected equations again. While a nonfinite guard is unresolved,
/// retain branch order and only eliminate literal boolean conditions.
pub(super) fn piecewise(pieces: Vec<(Expr, Expr)>) -> Expr {
    if !pieces.iter().any(|(_, condition)| substitute(condition, &mut |e| e.clone(), &|_| None).is_some()) {
        return crate::relational::piecewise(pieces);
    }
    let mut kept = Vec::new();
    for (value, condition) in pieces {
        match condition.kind() {
            Kind::Bool(false) => continue,
            Kind::Bool(true) => { kept.push((value, condition)); break; }
            _ => kept.push((value, condition)),
        }
    }
    if kept.is_empty() { Expr::nan() }
    else if kept.len() == 1 && kept[0].1 == Expr::boolean(true) { kept.pop().unwrap().0 }
    else { Expr::from_kind(Kind::Piecewise(kept)) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::RelOp;
    use crate::functions::exp;
    use crate::ops::add2;
    use crate::subs::{subs1, xreplace};
    use rustc_hash::FxHashMap;

    #[test]
    fn a_nonfinite_guard_preserves_the_equation_until_specialized() {
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        let y = Expr::symbol("y");
        let lhs = exp(&add2(&a, &exp(&x)));
        let equation = Expr::from_unevaluated_kind(Kind::Rel(RelOp::Eq, lhs, y.clone()));
        let guard = crate::logic::not(&crate::sets::contains_expr(&a, &crate::sets::complexes()));
        let condition = crate::logic::and(vec![guard, equation.clone()]);
        let partial = subs1(&condition, &y, &Expr::zero());
        assert_ne!(partial, Expr::boolean(false));
        assert_eq!(subs1(&partial, &a, &Expr::neg_oo()), Expr::boolean(true));
        assert_eq!(subs1(&partial, &a, &Expr::zero()), Expr::boolean(false));
        assert_eq!(subs1(&equation, &y, &Expr::zero()), Expr::boolean(false));
        let map = [(y, Expr::zero())].into_iter().collect::<FxHashMap<_, _>>();
        assert_eq!(xreplace(&condition, &map), partial);
    }

    #[test]
    fn finite_base_filtering_uses_the_same_conditional_context() {
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        let y = Expr::symbol("y");
        let equation = Expr::from_unevaluated_kind(Kind::Rel(RelOp::Eq,
            exp(&add2(&a, &exp(&x))), y.clone()));
        let guard = crate::logic::not(&crate::sets::contains_expr(&a, &crate::sets::complexes()));
        let base = crate::sets::finite_set(vec![Expr::zero(), Expr::one()]);
        let set = crate::sets::condition_set(&x, &crate::logic::and(vec![guard, equation]), &base);
        let partial = subs1(&set, &y, &Expr::zero());
        assert_eq!(subs1(&partial, &a, &Expr::neg_oo()), base);
        assert_eq!(subs1(&partial, &a, &Expr::zero()), crate::sets::empty_set());
    }

    #[test]
    fn opposing_infinities_stay_an_undefined_conditional_calculation() {
        let x=Expr::symbol("x"); let a=Expr::symbol("a"); let b=Expr::symbol("b"); let y=Expr::symbol("y");
        let source=crate::parser::parse_str("exp(exp(x)+a+b)-y").unwrap();
        let result=crate::sets::solveset(&source,&x,&crate::sets::complexes()).unwrap();
        let partial=subs1(&result,&y,&Expr::zero());
        let first=subs1(&subs1(&partial,&a,&Expr::neg_oo()),&b,&Expr::oo());
        let second=subs1(&subs1(&partial,&b,&Expr::oo()),&a,&Expr::neg_oo());
        let equation=Expr::from_unevaluated_kind(Kind::Rel(RelOp::Eq,Expr::nan(),Expr::zero()));
        let expected=crate::sets::condition_set(&x,&equation,&crate::sets::complexes());
        assert_eq!(first,expected); assert_eq!(second,expected);
        assert_eq!(subs1(&subs1(&partial,&a,&Expr::neg_oo()),&b,&Expr::one()),crate::sets::complexes());
    }

}
