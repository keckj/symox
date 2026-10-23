//! Structural `Expr.is_number`, distinct from finite numerical evaluation.
//!
//! This property describes the held expression, not its eventual value. In
//! particular, undefined functions and held matrix operations are not numbers,
//! while every numerical constant (including infinities and NaN) is a number.
use rustc_hash::FxHashMap;

use crate::expr::{Expr, Func, Kind};
use crate::matrixexpr::{MatrixExpr, MatrixKind, MatrixOperand};
use crate::sets::SetOp;

/// SymPy-compatible structural number classification, independent of query mode.
///
/// Two bits above the assumption and finite-numeric caches encode the completed
/// Boolean answer. Public argument reconstruction and arithmetic evaluation are
/// deliberately avoided, so shared held arithmetic remains a shared graph.
pub fn is_symbolic_number(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(_) | Kind::Const(_) => return true,
        Kind::Sym(_) | Kind::Bool(_) => return false,
        _ => {}
    }
    match e.cached_assumption(62, 2) {
        2 => return false,
        3 => return true,
        _ => {}
    }
    let answer = match e.kind() {
        Kind::Add(_, terms) => terms.iter().all(|(term, _)| is_symbolic_number(term)),
        Kind::Mul(_, factors) => factors
            .iter()
            .all(|(base, exponent)| is_symbolic_number(base) && is_symbolic_number(exponent)),
        Kind::Pow(base, exponent) => is_symbolic_number(base) && is_symbolic_number(exponent),
        Kind::Func(Func::Integral | Func::Sum | Func::Product, _) => ClosedQuery::default().closed(e, 0),
        Kind::Func(Func::Set(SetOp::IdentityFunction), _) => true,
        // A class-defined function (a Function subclass) is a number when its
        // arguments are, as SymPy's ``Expr.is_number``; an undefined function
        // (AppliedUndef) never is. Other host classes (atoms, Basic and Expr
        // subclasses) keep their own Python answer.
        // A root object is an atom-like number: its polynomial argument is
        // bound by the object, as the reference's ``ComplexRootOf.is_number``.
        Kind::Func(Func::Undefined(head), args) if crate::poly::rootof::is_rootof_head(head, args.len()) => true,
        // A host Number or NumberSymbol subclass atom (hbar) is a number.
        Kind::Func(Func::Undefined(head), args) if args.is_empty() && head.host().is_some_and(|host| host.atomic() && host.number()) => true,
        Kind::Func(Func::Undefined(head), args) => {
            let function = head.host().is_some_and(|host| !host.expression() && !host.atomic());
            (head.is_internal("UnevaluatedExpr") || (function && !args.is_empty()))
                && args.iter().all(is_symbolic_number)
        }
        // Tuple-valued arguments are Basic objects with is_number=False. The
        // typed forms below contain such arguments, or a matrix argument.
        Kind::Func(
            Func::Predicate(_)
            | Func::GenericExpr
            | Func::AppliedPredicate
            | Func::Subs
            | Func::Limit
            | Func::Logic(_)
            | Func::Set(_)
            | Func::MatrixProxy(_)
            | Func::MatrixElement(_)
            | Func::MatrixTrace(_)
            | Func::MatrixDeterminant(_)
            | Func::MatrixDotProduct(_, _),
            _,
        )
        | Kind::Derivative(..)
        | Kind::Rel(..)
        | Kind::Piecewise(_) => false,
        Kind::Func(Func::Order, args) => {
            // Native Order stores limits flat; SymPy exposes each as a Tuple.
            args.len() == 1 && is_symbolic_number(&args[0])
        }
        Kind::Func(_, args) => args.iter().all(is_symbolic_number),
        Kind::Num(_) | Kind::Const(_) | Kind::Sym(_) | Kind::Bool(_) => unreachable!(),
    };
    e.remember_assumption(62, 2, if answer { 3 } else { 2 });
    answer
}

/// Integral/Sum override the generic rule with absence of free variables,
/// including for undefined functions. Check that property directly under each
/// binding scope rather than constructing or sorting free-symbol collections.
/// Physical keys retain owners and distinguish a node reused under two scopes.
struct ClosedQuery {
    scopes: Vec<Vec<Expr>>,
    scalar: FxHashMap<(usize, usize), (Expr, bool)>,
    matrix: FxHashMap<(usize, usize), (MatrixExpr, bool)>,
}
impl Default for ClosedQuery {
    fn default() -> Self {
        Self {
            scopes: vec![Vec::new()],
            scalar: FxHashMap::default(),
            matrix: FxHashMap::default(),
        }
    }
}
impl ClosedQuery {
    fn binds(&self, scope: usize, symbol: &Expr) -> bool {
        self.scopes[scope].contains(symbol)
    }
    fn extend(&mut self, scope: usize, variables: impl IntoIterator<Item = Expr>) -> usize {
        let mut symbols = self.scopes[scope].clone();
        symbols.extend(variables.into_iter().filter(Expr::is_symbol));
        let index = self.scopes.len();
        self.scopes.push(symbols);
        index
    }
    fn all(&mut self, values: &[Expr], scope: usize) -> bool {
        values.iter().all(|value| self.closed(value, scope))
    }
    fn closed(&mut self, e: &Expr, scope: usize) -> bool {
        if let Kind::Sym(_) = e.kind() {
            return self.binds(scope, e);
        }
        if e.mask() == 0 {
            return true;
        }
        let key = (e.dag_id(), scope);
        if let Some((_, value)) = self.scalar.get(&key) {
            return *value;
        }
        let answer = match e.kind() {
            Kind::Add(_, terms) => terms.iter().all(|(term, _)| self.closed(term, scope)),
            Kind::Mul(_, factors) => factors
                .iter()
                .all(|(base, exponent)| self.closed(base, scope) && self.closed(exponent, scope)),
            Kind::Pow(base, exponent) | Kind::Rel(_, base, exponent) => {
                self.closed(base, scope) && self.closed(exponent, scope)
            }
            // Differentiation variables are not free just by being present;
            // native derivative counts are literal integers.
            Kind::Derivative(body, _) => self.closed(body, scope),
            Kind::Piecewise(pieces) => pieces.iter().all(|(value, condition)| {
                self.closed(value, scope) && self.closed(condition, scope)
            }),
            Kind::Func(f @ (Func::Integral | Func::Sum | Func::Product), args) => self.calculus(f, args, scope),
            Kind::Func(Func::Subs, args) if args.len() == 3 => {
                let body_scope = self.extend(scope, args[1].args());
                self.closed(&args[0], body_scope) && self.closed(&args[2], scope)
            }
            Kind::Func(Func::Limit, args) if args.len() == 4 => {
                let body_scope = self.extend(scope, [args[1].clone()]);
                self.closed(&args[0], body_scope) && self.closed(&args[2], scope)
            }
            Kind::Func(Func::Set(SetOp::Lambda | SetOp::IdentityFunction), _) => {
                if let Some((variables, body)) = crate::sets::lambda_parts(e) {
                    let body_scope = self.extend(scope, variables);
                    self.closed(&body, body_scope)
                } else {
                    false
                }
            }
            Kind::Func(Func::Set(SetOp::ConditionSet), args) if args.len() == 3 => {
                let condition_scope = self.extend(scope, [args[0].clone()]);
                self.closed(&args[1], condition_scope) && self.closed(&args[2], scope)
            }
            Kind::Func(
                Func::MatrixProxy(matrix)
                | Func::MatrixElement(matrix)
                | Func::MatrixTrace(matrix)
                | Func::MatrixDeterminant(matrix),
                args,
            ) => self.closed_matrix(matrix, scope) && self.all(args, scope),
            Kind::Func(Func::MatrixDotProduct(left, right), _) => {
                self.closed_matrix(left, scope) && self.closed_matrix(right, scope)
            }
            Kind::Func(_, args) => self.all(args, scope),
            Kind::Num(_) | Kind::Const(_) | Kind::Bool(_) => true,
            Kind::Sym(_) => unreachable!(),
        };
        self.scalar.insert(key, (e.clone(), answer));
        answer
    }
    fn calculus(&mut self, f: &Func, args: &[Expr], mut scope: usize) -> bool {
        if crate::binders::validate(f, args).is_err() {
            return false;
        }
        // Outer limits bind the earlier limits and the body, never their own
        // endpoints. An indefinite integral introduces its variable as free.
        for limit in args[1..].iter().rev() {
            let values = limit.args();
            if values.len() == 1 {
                if !self.binds(scope, &values[0]) {
                    return false;
                }
            } else {
                if !self.all(&values[1..], scope) {
                    return false;
                }
                scope = self.extend(scope, [values[0].clone()]);
            }
        }
        self.closed(&args[0], scope)
    }
    fn closed_matrix(&mut self, matrix: &MatrixExpr, scope: usize) -> bool {
        let key = (matrix.dag_id(), scope);
        if let Some((_, value)) = self.matrix.get(&key) {
            return *value;
        }
        // Visit local metadata once, rather than scalar_fields(), which itself
        // traverses the complete matrix DAG. A MatrixSymbol is a free atom,
        // regardless of its dimensions.
        let local = match matrix.kind() {
            MatrixKind::Symbol(_) => false,
            MatrixKind::Function(lambda) => {
                self.closed(matrix.rows(), scope)
                    && self.closed(matrix.cols(), scope)
                    && self.closed(lambda, scope)
            }
            MatrixKind::ApplyFunction(function, _) => self.closed(function, scope),
            MatrixKind::Identity | MatrixKind::Zero | MatrixKind::One => {
                self.closed(matrix.rows(), scope) && self.closed(matrix.cols(), scope)
            }
            MatrixKind::Explicit(entries) => self.all(entries.data(), scope),
            MatrixKind::Mul(coefficient, _) | MatrixKind::Pow(_, coefficient) => {
                self.closed(coefficient, scope)
            }
            MatrixKind::Slice(_, rows, columns) => rows
                .fields()
                .iter()
                .chain(columns.fields().iter())
                .all(|field| self.closed(field, scope)),
            MatrixKind::HadamardPower(left, right) => [left, right].iter().all(|operand| {
                if let MatrixOperand::Scalar(scalar) = operand {
                    self.closed(scalar, scope)
                } else {
                    true
                }
            }),
            _ => true,
        };
        let answer = local
            && matrix
                .children()
                .iter()
                .all(|child| self.closed_matrix(child, scope));
        self.matrix.insert(key, (matrix.clone(), answer));
        answer
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Assumptions, Number};
    fn p(source: &str) -> Expr {
        crate::parser::parse_str(source).unwrap()
    }
    fn pair(left: Expr, right: Expr) -> Expr {
        Expr::from_kind(Kind::Func(
            Func::Undefined("pair".into()),
            vec![left, right],
        ))
    }
    fn integral(body: &Expr, variable: &Expr) -> Expr {
        crate::binders::held(
            Func::Integral,
            body,
            vec![vec![variable.clone(), Expr::zero(), Expr::one()]],
        )
        .unwrap()
    }
    #[test]
    fn structural_numbers_include_nonfinite_and_held_builtin_functions() {
        for source in [
            "1", "1/2", "pi", "I", "oo", "-oo", "zoo", "nan", "sin(1)", "sqrt(2)", "O(1)",
        ] {
            assert!(is_symbolic_number(&p(source)), "{source}");
        }
        for source in [
            "x",
            "f(1)",
            "f()",
            "sin(x)",
            "Lambda(x,1)",
            "Lambda((),1)",
            "O(x)",
        ] {
            assert!(!is_symbolic_number(&p(source)), "{source}");
        }
        assert!(!is_symbolic_number(&crate::special::tuple(vec![
            Expr::one(),
            Expr::int(2)
        ])));
        assert!(is_symbolic_number(&crate::sets::identity_function()));
        assert!(!crate::functions::is_numeric_value(&Expr::zoo()));
        assert!(!crate::functions::is_numeric_value(&Expr::nan()));
    }
    #[test]
    fn held_arithmetic_looks_at_its_arguments_without_evaluation() {
        let x = Expr::symbol("x");
        let product =
            Expr::from_unevaluated_kind(Kind::Mul(Number::ZERO, vec![(x.clone(), Expr::one())]));
        let power = Expr::from_unevaluated_kind(Kind::Pow(x.clone(), Expr::zero()));
        let sum = Expr::from_unevaluated_kind(Kind::Add(Number::ONE, vec![(x, Number::ZERO)]));
        for e in [product, power, sum] {
            assert!(!is_symbolic_number(&e));
        }
    }
    #[test]
    fn calculus_closure_respects_scopes_and_ignores_function_evaluability() {
        for source in [
            "Integral(f(x),(x,0,1))",
            "Sum(f(),(x,0,1))",
            "Integral(1,(x,0,1))",
            "Integral(Integral(1,x),(x,0,1))",
            "Integral(Sum(x,(x,0,y)),(y,0,1))",
            "Integral(Limit(f(x),x,y),(y,0,1))",
            "Integral(Subs(f(x),x,y),(y,0,1))",
            "Integral(Lambda(x,x+y),(y,0,1))",
        ] {
            assert!(is_symbolic_number(&p(source)), "{source}");
        }
        for source in [
            "Integral(1,x)",
            "Integral(f(y),(x,0,1))",
            "Integral(1,(x,x,1))",
            "Integral(Integral(1,y),(x,0,1))",
            "Sum(f(y),(x,0,1))",
            "Limit(1,x,0)",
            "Subs(1,x,1)",
        ] {
            assert!(!is_symbolic_number(&p(source)), "{source}");
        }
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let body = pair(x.clone(), x.clone());
        assert!(!is_symbolic_number(&body));
        assert!(is_symbolic_number(&integral(&body, &x)));
        assert!(!is_symbolic_number(&integral(&body, &y)));
        let derivative = crate::diff::derivative_unevaluated(&Expr::one(), vec![(x, 1)]);
        assert!(!is_symbolic_number(&derivative));
        // Differentiation variables alone do not become free symbols.
        assert!(is_symbolic_number(&integral(&derivative, &y)));
    }
    #[test]
    fn arbitrary_shared_arithmetic_and_closed_binder_dags_remain_linear() {
        let mut number = Expr::one();
        for _ in 0..80 {
            number = Expr::from_unevaluated_kind(Kind::Add(
                Number::ZERO,
                vec![(number.clone(), Number::ONE), (number, Number::ONE)],
            ));
        }
        assert!(is_symbolic_number(&number));
        assert_eq!(number.cached_assumption(62, 2), 3);
        let x = Expr::dummy("x", Assumptions::NONE);
        let mut body = x.clone();
        for _ in 0..80 {
            body = pair(body.clone(), body);
        }
        assert!(!is_symbolic_number(&body));
        assert!(is_symbolic_number(&integral(&body, &x)));
        assert!(!is_symbolic_number(&body));
    }
    #[test]
    fn concurrent_number_cache_publication_keeps_assumption_fields_separate() {
        let e = p("pi-E");
        let x = Expr::symbol("x");
        let mut body = x.clone();
        for _ in 0..80 {
            body = pair(body.clone(), body);
        }
        let closed = integral(&body, &x);
        let workers: Vec<_> = (0..8)
            .map(|thread| {
                let e = e.clone();
                let body = body.clone();
                let closed = closed.clone();
                std::thread::spawn(move || {
                    for _ in 0..20 {
                        assert!(is_symbolic_number(&e));
                        assert!(is_symbolic_number(&closed));
                        assert!(!is_symbolic_number(&body));
                        if thread % 2 == 0 {
                            assert_eq!(
                                crate::assumptions::exact_queries(|| {
                                    crate::assumptions::is_positive(&e)
                                }),
                                None
                            );
                        } else {
                            assert_eq!(crate::assumptions::is_positive(&e), Some(true));
                        }
                    }
                })
            })
            .collect();
        for worker in workers {
            worker.join().unwrap();
        }
        assert!(is_symbolic_number(&closed));
    }

    #[test]
    fn mode_independent_number_cache_preserves_finite_assumption_banks() {
        for first_exact in [true, false] {
            let e = Expr::from_unevaluated_kind(Kind::Add(
                Number::ZERO,
                vec![(Expr::pi(), Number::ONE), (Expr::e(), Number::MINUS_ONE)],
            ));
            if first_exact {
                assert!(crate::assumptions::exact_queries(|| is_symbolic_number(&e)));
            }
            assert!(is_symbolic_number(&e));
            assert_eq!(crate::assumptions::is_positive(&e), Some(true));
            assert_eq!(
                crate::assumptions::exact_queries(|| crate::assumptions::is_positive(&e)),
                None
            );
            assert!(crate::assumptions::exact_queries(|| is_symbolic_number(&e)));
            assert_eq!(e.cached_assumption(62, 2), 3);
        }
    }
}
