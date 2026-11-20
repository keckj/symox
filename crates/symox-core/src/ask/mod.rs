//! SymPy's new assumptions system: `Q` predicates, `ask` and `refine`.
//!
//! A predicate is the atom `Func::Predicate(p)`; applying it builds
//! `Func::AppliedPredicate` with arguments `[predicate, args...]`, which
//! mirrors SymPy's `AppliedPredicate.args`. The query engine ports SymPy
//! 1.14's `ask` pipeline: local fact extraction with the generated known
//! facts, predicate handlers, then `satask` with the exact class facts
//! (`sathandlers.py`), decided by the native CDCL solver.

mod cnf;
mod engine;
mod facts_data;
mod handlers;
mod lra;
mod predicates;
mod refine;

pub use engine::{ask, ask_with, AskOptions, ExternalHandler};
pub use predicates::Pred;
pub use refine::{refine, refine_with};

use crate::expr::{Expr, Func, Kind};

/// Errors raised by the assumptions engine (SymPy exception classes).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AskError {
    Value(String),
    Type(String),
}

impl std::fmt::Display for AskError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AskError::Value(m) | AskError::Type(m) => f.write_str(m),
        }
    }
}

/// The predicate atom `Q.<name>`.
pub fn predicate(p: Pred) -> Expr {
    Expr::from_kind(Kind::Func(Func::Predicate(p), vec![]))
}

/// The predicate of a predicate atom.
pub fn as_predicate(e: &Expr) -> Option<Pred> {
    match e.kind() {
        Kind::Func(Func::Predicate(p), _) => Some(*p),
        _ => None,
    }
}

/// `p(*args)`. Binary relations take exactly two arguments.
pub fn applied(p: Pred, args: Vec<Expr>) -> Result<Expr, AskError> {
    if p.is_binrel() && args.len() != 2 {
        return Err(AskError::Value(format!(
            "Binary relation takes two arguments, but got {}.",
            args.len()
        )));
    }
    Ok(applied_unchecked(p, args))
}

pub(crate) fn applied_unchecked(p: Pred, args: Vec<Expr>) -> Expr {
    let mut all = Vec::with_capacity(args.len() + 1);
    all.push(predicate(p));
    all.extend(args);
    Expr::from_kind(Kind::Func(Func::AppliedPredicate, all))
}

/// Unary application `p(arg)`.
pub fn q(p: Pred, arg: &Expr) -> Expr {
    applied_unchecked(p, vec![arg.clone()])
}

/// The predicate and arguments of an applied predicate.
pub fn as_applied(e: &Expr) -> Option<(Pred, &[Expr])> {
    match e.kind() {
        Kind::Func(Func::AppliedPredicate, a) if !a.is_empty() => {
            as_predicate(&a[0]).map(|p| (p, &a[1..]))
        }
        _ => None,
    }
}

/// Rebuild an applied predicate from its SymPy args `(predicate, *arguments)`.
pub fn applied_from_args(args: Vec<Expr>) -> Result<Expr, AskError> {
    let Some(p) = args.first().and_then(as_predicate) else {
        return Err(AskError::Type(
            "AppliedPredicate needs a Predicate as first argument".into(),
        ));
    };
    applied(p, args[1..].to_vec())
}

/// `Q.<name>` printing: `Q.positive`, or the bare name for `Predicate(name)`
/// (SymPy prints both as `Q.name`).
pub fn predicate_str(p: Pred) -> String {
    format!("Q.{}", p.name())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{functions::func, logic, ops};

    #[test]
    fn queries_and_refinement_share_the_native_predicates() {
        let x = Expr::symbol("x");
        let a = q(Pred::Positive, &x);
        assert_eq!(ask(&q(Pred::Real, &x), &a), Ok(Some(true)));
        assert_eq!(ask(&q(Pred::Integer, &x), &a), Ok(None));
        assert_eq!(refine(&func(Func::Abs, vec![x.clone()]), &a), Ok(x.clone()));
        let contradictory = logic::and(vec![a, q(Pred::Negative, &x)]);
        assert!(matches!(
            ask(&q(Pred::Real, &x), &contradictory),
            Err(AskError::Value(_))
        ));
        let square_root = ops::sqrt(&ops::pow(&x, &Expr::int(2)));
        assert_eq!(
            refine(&square_root, &q(Pred::Real, &x)),
            Ok(func(Func::Abs, vec![x]))
        );
    }

    #[test]
    fn shared_scalar_relevance_does_not_expand_the_expression_tree() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let mut e = ops::add2(&x, &y);
        for _ in 0..80 {
            e = Expr::from_unevaluated_kind(Kind::Add(
                crate::number::Number::Int(0),
                vec![
                    (e.clone(), crate::number::Number::Int(1)),
                    (e, crate::number::Number::Int(1)),
                ],
            ));
        }
        let mut found = rustc_hash::FxHashSet::default();
        engine::find_symbols(&q(Pred::Integer, &e), &mut found);
        assert_eq!(found.len(), 2);
        assert!(found.contains(&x) && found.contains(&y));
    }
}
