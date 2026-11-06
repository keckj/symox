//! Native `nth_algebraic` and `factorable` solvers.
//!
//! Ports of SymPy 1.14 `NthAlgebraic._matches` and `Factorable._matches` /
//! `Factorable._get_general_solution` (`sympy/solvers/ode/single.py`).
//!
//! The rewriting steps (derivatives replaced by the invertible `diffx` head
//! and restored afterwards, `collect` with a coefficient callback,
//! `fraction`, the factor list, `as_base_exp`, the exponent test, the
//! equation lists and the `Eq` construction) run here. `solve`, `dsolve`,
//! `simplify`, `factor` and `cancel` are the public callables supplied by
//! the binding through [`Hooks`]; an error they raise propagates unchanged,
//! as in the reference. `Outcome::Decline` leaves the decision to the
//! reference body.

use rustc_hash::FxHashMap;

use crate::diff::derivative_unevaluated;
use crate::expr::{Expr, Func, Kind, RelOp};
use crate::functions::func;
use crate::ops::sub;
use crate::poly::collect::{self, Collection, Options};
use crate::relational::rel;
use crate::traverse::{has_structural as has, map_children, ordered_args, preorder};

use super::classify::Outcome;

/// The public callables of the reference bodies. `Ok(None)` from `solve` or
/// `dsolve` stands for a raised `NotImplementedError`.
pub trait Hooks {
    type Error;
    /// `solve(eq, sym)` (with the flags of the calling solver): the list of
    /// solutions.
    fn solve(&mut self, eq: &Expr, sym: &Expr) -> Result<Option<Vec<Expr>>, Self::Error>;
    /// `dsolve(eq, fx)`: its equation, or the list of its equations.
    fn dsolve(&mut self, eq: &Expr, fx: &Expr) -> Result<Option<Vec<Expr>>, Self::Error>;
    fn simplify(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    fn factor(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    fn cancel(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
    fn expand(&mut self, e: &Expr) -> Result<Expr, Self::Error>;
}

/// Content whose rebuilding by `Basic.replace` the port does not model:
/// matrices, symbolic derivatives and host-defined classes.
fn unsupported(e: &Expr) -> bool {
    if e.has_matrix() || e.has_host() {
        return true;
    }
    let mut found = false;
    preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Func(Func::SymbolicDerivative | Func::MatrixProxy(_), _)) {
            found = true;
        }
    });
    found
}

/// Bottom-up `Basic.replace`: `f` sees each node after its children were
/// rebuilt. Shared subtrees are rewritten once.
fn replace_bottom_up(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr) -> Expr {
    fn walk(e: &Expr, f: &mut dyn FnMut(&Expr) -> Expr, memo: &mut FxHashMap<usize, (Expr, Expr)>) -> Expr {
        if let Some((_, v)) = memo.get(&e.dag_id()) {
            return v.clone();
        }
        let rebuilt = match e.kind() {
            // `Derivative(expr, *variable_count)` rebuilt with the new
            // expression merges a nested derivative, as its constructor does.
            Kind::Derivative(inner, vars) => {
                let new = walk(inner, f, memo);
                if new.same_presentation(inner) {
                    e.clone()
                } else {
                    derivative_unevaluated(&new, vars.clone())
                }
            }
            _ => map_children(e, &mut |c| walk(c, f, memo)),
        };
        let value = f(&rebuilt);
        memo.insert(e.dag_id(), (e.clone(), value.clone()));
        value
    }
    walk(e, f, &mut FxHashMap::default())
}

fn diffx_head(template: &Expr) -> Option<Func> {
    match template.kind() {
        Kind::Func(head @ Func::Undefined(_), args) if args.len() == 1 => Some(head.clone()),
        _ => None,
    }
}

/// `replace(eq, var)` of `NthAlgebraic._matches`: each derivative order in
/// `var` becomes one application of `diffx`; other variables stay
/// derivatives.
fn to_diffx(eq: &Expr, x: &Expr, diffx: &Func) -> Expr {
    replace_bottom_up(eq, &mut |node| match node.kind() {
        Kind::Derivative(inner, vars) => {
            let mut out = inner.clone();
            for (v, n) in vars {
                for _ in 0..*n {
                    out = if v == x { func(diffx.clone(), vec![out]) } else { derivative_unevaluated(&out, vec![(v.clone(), 1)]) };
                }
            }
            out
        }
        _ => node.clone(),
    })
}

/// `unreplace(eq, var)`: every `diffx(e)` becomes `Derivative(e, var)`.
fn from_diffx(e: &Expr, x: &Expr, diffx: &Func) -> Expr {
    replace_bottom_up(e, &mut |node| match node.kind() {
        Kind::Func(head, args) if head == diffx && args.len() == 1 => derivative_unevaluated(&args[0], vec![(x.clone(), 1)]),
        _ => node.clone(),
    })
}

/// `NthAlgebraic._matches`: the solutions `Eq(fx, simplify(unreplace(s)))`
/// for `s` in `solve(replace(eq), fx, simplify=False)` (an empty list for
/// no match). `diffx_template` is `diffx(var)`, an application of the
/// invertible head of the reference; its inverse (an `Integral` plus a new
/// constant) is applied by `solve` through the host class.
pub fn nth_algebraic<H: Hooks>(eq: &Expr, fx: &Expr, x: &Expr, diffx_template: &Expr, hooks: &mut H) -> Result<Outcome<Vec<Expr>>, H::Error> {
    let Some(diffx) = diffx_head(diffx_template) else {
        return Ok(Outcome::Decline);
    };
    if unsupported(eq) || !x.is_symbol() {
        return Ok(Outcome::Decline);
    }
    let subs_eqn = to_diffx(eq, x, &diffx);
    let solns = hooks.solve(&subs_eqn, fx)?.unwrap_or_default();
    let mut out = Vec::with_capacity(solns.len());
    for s in &solns {
        let s = hooks.simplify(&from_diffx(s, x, &diffx))?;
        out.push(rel(RelOp::Eq, fx, &s));
    }
    Ok(Outcome::Match(out))
}

/// `Mul.make_args(e)` in `.args` order.
fn mul_args(e: &Expr) -> Vec<Expr> {
    if e.is_mul() {
        ordered_args(e)
    } else {
        vec![e.clone()]
    }
}

/// `Factorable._matches`: `Match(eqs)` (the reference `self.eqs` of a
/// match) or `NoMatch`.
pub fn factorable<H: Hooks>(eq_orig: &Expr, fx: &Expr, x: &Expr, hooks: &mut H) -> Result<Outcome<Vec<Expr>>, H::Error> {
    if unsupported(eq_orig) || crate::decompose::commutativity(eq_orig) != Some(true) || !x.is_symbol() {
        return Ok(Outcome::Decline);
    }
    let df = derivative_unevaluated(fx, vec![(x.clone(), 1)]);
    // eq_orig.collect(f(x), func=cancel)
    let mut failure = None;
    let collected = {
        let mut callback = |e: &Expr| -> Result<Expr, ()> {
            hooks.cancel(e).map_err(|err| {
                failure = Some(err);
            })
        };
        collect::run(eq_orig, std::slice::from_ref(fx), Options { exact: false, discover: false, distribute_order_term: true }, true, Some(&mut callback))
    };
    if let Some(err) = failure {
        return Err(err);
    }
    let eq = match collected {
        Ok(Collection::Expression(v)) => v,
        _ => return Ok(Outcome::Decline),
    };
    let eq = crate::simptools::radsimp::fraction(&hooks.factor(&eq)?, false).0;
    let factors = mul_args(&hooks.factor(&eq)?);
    let roots: Vec<(Expr, Expr)> = factors.iter().filter(|f| !ordered_args(f).is_empty()).map(crate::decompose::as_base_exp).collect();
    let Some(first) = roots.first() else {
        // `roots[0]` raises IndexError in the reference.
        return Ok(Outcome::Decline);
    };
    let several = roots.len() > 1 || {
        match first.1.as_number() {
            Some(n) if n.is_rational() => n.cmp_value(&crate::number::Number::int(1)) == Some(std::cmp::Ordering::Greater),
            _ => return Ok(Outcome::Decline),
        }
    };
    if several {
        let eqs: Vec<Expr> = roots.iter().filter(|(base, _)| has(base, fx)).map(|(base, _)| base.clone()).collect();
        if !eqs.is_empty() {
            return Ok(Outcome::Match(eqs));
        }
    }
    let Some(sols) = hooks.solve(&eq, &df)? else {
        // NotImplementedError propagates from the reference matcher.
        return Ok(Outcome::Decline);
    };
    if !sols.is_empty() {
        let eqs: Vec<Expr> = sols.iter().map(|root| sub(&df, root)).collect();
        // Avoid infinite recursion
        return Ok(if eqs.len() == 1 && &eqs[0] == eq_orig { Outcome::NoMatch } else { Outcome::Match(eqs) });
    }
    let eqs: Vec<Expr> = factors.iter().filter(|f| has(f, fx)).cloned().collect();
    Ok(if !eqs.is_empty() && factors.len() > 1 { Outcome::Match(eqs) } else { Outcome::NoMatch })
}

/// `Factorable._get_general_solution`: `Ok(solutions)`, or `Err(eq)` with
/// the last equation when every `dsolve` raised `NotImplementedError` (the
/// reference then raises `NotImplementedError` naming it).
pub fn factorable_solutions<H: Hooks>(eqs: &[Expr], fx: &Expr, hooks: &mut H) -> Result<Result<Vec<Expr>, Expr>, H::Error> {
    let mut sols = Vec::new();
    for eq in eqs {
        if let Some(sol) = hooks.dsolve(eq, fx)? {
            sols.extend(sol);
        }
    }
    match eqs.last() {
        Some(last) if sols.is_empty() => Ok(Err(last.clone())),
        _ => Ok(Ok(sols)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Func;
    use crate::functions::undefined;
    use crate::ops::{add, mul, mul2, pow};

    struct Mock {
        solves: Vec<(Expr, Expr)>,
        answers: Vec<Option<Vec<Expr>>>,
    }

    impl Hooks for Mock {
        type Error = ();
        fn solve(&mut self, eq: &Expr, sym: &Expr) -> Result<Option<Vec<Expr>>, ()> {
            self.solves.push((eq.clone(), sym.clone()));
            Ok(self.answers.remove(0))
        }
        fn dsolve(&mut self, eq: &Expr, _fx: &Expr) -> Result<Option<Vec<Expr>>, ()> {
            Ok(if has(eq, &Expr::symbol("NI")) { None } else { Some(vec![eq.clone()]) })
        }
        fn simplify(&mut self, e: &Expr) -> Result<Expr, ()> { Ok(e.clone()) }
        fn factor(&mut self, e: &Expr) -> Result<Expr, ()> { crate::poly::api::factor(e).map_err(|_| ()) }
        fn cancel(&mut self, e: &Expr) -> Result<Expr, ()> { crate::poly::api::cancel(e).map_err(|_| ()) }
        fn expand(&mut self, e: &Expr) -> Result<Expr, ()> { Ok(crate::expand::expand(e)) }
    }

    fn setup() -> (Expr, Expr, Expr) {
        let x = Expr::symbol("x");
        let f = undefined("f", vec![x.clone()]);
        let df = derivative_unevaluated(&f, vec![(x.clone(), 1)]);
        (x, f, df)
    }

    #[test]
    fn derivatives_become_diffx_and_back() {
        let (x, f, _) = setup();
        let y = Expr::symbol("y");
        let template = undefined("diffx", vec![x.clone()]);
        let Some(head) = diffx_head(&template) else { panic!() };
        let d2 = derivative_unevaluated(&f, vec![(x.clone(), 2)]);
        let r = to_diffx(&add(&[d2.clone(), mul2(&y, &f)]), &x, &head);
        let dd = func(head.clone(), vec![func(head.clone(), vec![f.clone()])]);
        assert_eq!(r, add(&[dd.clone(), mul2(&y, &f)]));
        assert_eq!(from_diffx(&dd, &x, &head), d2);
        assert!(matches!(head, Func::Undefined(_)));
    }

    #[test]
    fn nth_algebraic_wraps_solutions() {
        let (x, f, df) = setup();
        let template = undefined("diffx", vec![x.clone()]);
        let mut m = Mock { solves: vec![], answers: vec![Some(vec![Expr::zero()])] };
        let eq = mul2(&f, &df);
        let Ok(Outcome::Match(v)) = nth_algebraic(&eq, &f, &x, &template, &mut m) else { panic!() };
        assert_eq!(v, vec![rel(RelOp::Eq, &f, &Expr::zero())]);
        let mut m = Mock { solves: vec![], answers: vec![None] };
        let Ok(Outcome::Match(v)) = nth_algebraic(&eq, &f, &x, &template, &mut m) else { panic!() };
        assert!(v.is_empty());
    }

    #[test]
    fn factorable_splits_powers_and_factors() {
        let (x, f, df) = setup();
        // (f**2 - 4)*(f' + f): two factors with f(x)
        let a = add(&[pow(&f, &Expr::int(2)), Expr::int(-4)]);
        let b = add(&[df.clone(), f.clone()]);
        let eq = crate::expand::expand(&mul(&[a, b.clone()]));
        let mut m = Mock { solves: vec![], answers: vec![] };
        let Ok(Outcome::Match(eqs)) = factorable(&eq, &f, &x, &mut m) else { panic!() };
        assert_eq!(eqs.len(), 3);
        assert!(m.solves.is_empty());
        // A single linear factor goes through solve(eq, df).
        let mut m = Mock { solves: vec![], answers: vec![Some(vec![mul2(&Expr::int(-1), &f)])] };
        let Ok(Outcome::NoMatch) = factorable(&b, &f, &x, &mut m) else { panic!() };
        assert_eq!(m.solves.len(), 1);
    }

    #[test]
    fn factorable_solutions_report_last_equation() {
        let (_, f, _) = setup();
        let ni = Expr::symbol("NI");
        let mut m = Mock { solves: vec![], answers: vec![] };
        let Ok(Ok(sols)) = factorable_solutions(&[f.clone(), ni.clone()], &f, &mut m) else { panic!() };
        assert_eq!(sols, vec![f.clone()]);
        let Ok(Err(last)) = factorable_solutions(&[ni.clone(), mul2(&ni, &f)], &f, &mut m) else { panic!() };
        assert_eq!(last, mul2(&ni, &f));
    }
}
