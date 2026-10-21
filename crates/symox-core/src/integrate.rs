//! Symbolic integration (indefinite and definite).
//!
//! A rule based integrator in the spirit of SymPy's `manualintegrate`:
//! linearity, power rule, a table of elementary antiderivatives composed
//! with linear arguments, rational functions, trigonometric powers,
//! integration by parts (LIATE), cyclic exponential/trigonometric
//! products and u-substitution driven by the integrand's sub-expressions.

use rustc_hash::FxHashSet;

use crate::diff::diff;
use crate::expand::expand;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::*;
use crate::ops::{add, add2, div, mul, mul2, neg, pow, sqrt, sub};
use crate::relational::{ne, piecewise};
use crate::simplify::numer_denom;
use crate::subs::subs1;
use crate::binders::depends_on as has;

const MAX_DEPTH: u32 = 10;

mod rational_definite;
mod dirac;
pub(crate) mod meijerint;
pub(crate) mod piecewise;
pub use piecewise::piecewise_intervals;

/// Stage probe of the Meijer G definite integration (testing only).
/// `meijerint_indefinite(f, x)` with the errors the reference raises, as
/// `(class name, message)`.
pub fn meijerint_indefinite_checked(f: &Expr, x: &Expr) -> Result<Option<Expr>, (&'static str, String)> {
    meijerint::meijerint_indefinite_checked(f, x)
}

pub fn meijerint_debug(op: &str, args: &[Expr]) -> Option<Vec<Expr>> {
    meijerint::debug(op, args)
}

/// Exact absolute value at a distribution support point. None means the
/// arithmetic probe was inconclusive or exceeded its local resource budget.
pub fn distribution_abs_at(expression: &Expr, variable: &Expr, point: &Expr) -> Option<Expr> {
    if !crate::matrix::functions::scalar_substitution_preflight(expression, variable, point) {
        return None;
    }
    crate::assumptions::exact_queries(|| {
        let value = subs1(expression, variable, point);
        if has_bad(&value) { None } else { Some(abs(&value)) }
    })
}

/// Indefinite integral of `f` with respect to `x`. Returns an unevaluated
/// `Integral(f, x)` when no rule applies.
pub fn integrate(f: &Expr, x: &Expr) -> Expr {
    integrate_hint(f, x, None)
}

/// `integrate(f, x, meijerg=...)`: `Some(true)` integrates only by Meijer
/// G-functions, `Some(false)` never, `None` falls back to them when the
/// other rules fail.
pub fn integrate_hint(f: &Expr, x: &Expr, meijerg: Option<bool>) -> Expr {
    let g = piecewise::doit_rewrite(f, x, None);
    doit_indefinite(&g, x, meijerg).unwrap_or_else(|| integral_unevaluated(&g, x))
}

/// One indefinite limit of `Integral.doit` on the rewritten integrand: a
/// Piecewise goes through `Piecewise._eval_integral`, anything else through
/// `_eval_integral`. `None` when the whole integral stays unevaluated.
fn doit_indefinite(g: &Expr, x: &Expr, meijerg: Option<bool>) -> Option<Expr> {
    if piecewise::is_piecewise(g) {
        return match piecewise::eval_integral(g, x, true, meijerg) {
            Ok(r) if !is_held_integral(&r) => Some(r),
            _ => None,
        };
    }
    try_integrate_hint(g, x, meijerg)
}

fn is_held_integral(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Integral, _))
}

/// How `integrate(..., conds=...)` presents the convergence conditions of a
/// Meijer G definite integral.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Conds {
    Piecewise,
    Separate,
    None,
}

/// Definite integral of `f` over `[a, b]` with `conds='piecewise'`.
pub fn integrate_definite(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Expr {
    integrate_definite_conds(f, x, a, b, Conds::Piecewise, None).0
}

fn has_oo(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => true,
        _ => e.args().iter().any(has_oo),
    }
}

fn has_piecewise_or_minmax(e: &Expr) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if matches!(s.kind(), Kind::Piecewise(..) | Kind::Func(Func::Min | Func::Max, _)) {
            found = true;
        }
    });
    found
}

/// One `(x, a, b)` limit of `Integral.doit` with the `conds` and `meijerg`
/// hints. The second component is the separate condition, present only
/// when `conds='separate'` and the Meijer G method gave the value.
///
/// As in the reference, the variable is first replaced by a dummy carrying
/// the sign the bounds imply, and an improper integral over extended real
/// bounds tries `meijerint_definite` before any antiderivative. Its
/// convergence condition is kept, so a parameter range where the integral
/// diverges never receives the analytic continuation as its value.
pub fn integrate_definite_conds(f: &Expr, x: &Expr, a: &Expr, b: &Expr, conds: Conds, meijerg: Option<bool>) -> (Expr, Option<Expr>) {
    use crate::assumptions::{is_nonnegative, is_nonpositive, is_real};
    use crate::symbol::Assumptions;
    if f.is_zero() || a == b {
        return (Expr::zero(), None);
    }
    let both = |p: &dyn Fn(&Expr) -> Option<bool>| p(a) == Some(true) && p(b) == Some(true);
    // A signed dummy would evaluate DiracDelta(x) to zero and lose the half
    // weight of a delta at a finite endpoint, which the native distribution
    // rules keep (a documented difference from the reference, which returns
    // zero there): only realness is attached in that case.
    let mut dirac = false;
    crate::traverse::preorder(f, &mut |s| dirac |= matches!(s.kind(), Kind::Func(Func::DiracDelta, _)));
    let sign = if dirac {
        None
    } else if both(&is_nonnegative) && is_nonnegative(x) != Some(true) {
        Some(Assumptions::POSITIVE)
    } else if both(&is_nonpositive) && is_nonpositive(x) != Some(true) {
        Some(Assumptions::NEGATIVE)
    } else {
        None
    };
    let sign = sign.or_else(|| (both(&is_real) && is_real(x) != Some(true)).then_some(Assumptions::REAL));
    if let Some(sign) = sign {
        if x.is_symbol() {
            let d = Expr::dummy("x", sign);
            let rep = |e: &Expr, from: &Expr, to: &Expr| {
                let mut m = rustc_hash::FxHashMap::default();
                m.insert(from.clone(), to.clone());
                crate::subs::xreplace(e, &m)
            };
            let (v, c) = integrate_definite_conds(&rep(f, x, &d), &d, &rep(a, x, &d), &rep(b, x, &d), conds, meijerg);
            return (rep(&v, &d, x), c.map(|c| rep(&c, &d, x)));
        }
    }
    let rewritten = piecewise::doit_rewrite(f, x, Some((a, b)));
    let f = &rewritten;
    if piecewise::is_piecewise(f) {
        // `piecewise_integrate`, then the Piecewise evaluation; the Meijer G
        // definite method is not tried on a Piecewise
        let big_f = piecewise::piecewise_integrate(f, x, meijerg);
        let value = piecewise::eval_split(&big_f, x, a, b, &|o| eval_generic(o, x, a, b));
        return (value.unwrap_or_else(|| integral_definite_unevaluated(f, x, a, b)), None);
    }
    let extended_real = |e: &Expr| crate::assumptions::fact(e, "extended_real") == Some(true);
    let improper = extended_real(a) && extended_real(b) && (has_oo(a) || has_oo(b));
    let pw = has_piecewise_or_minmax(f);
    let try_meijerg = || -> Option<(Expr, Option<Expr>)> {
        let (value, cond) = meijerint::meijerint_definite(f, x, a, b)?;
        Some(match conds {
            Conds::Piecewise if matches!(cond.kind(), Kind::Bool(true)) => (value, None),
            Conds::Piecewise => (
                Expr::from_kind(Kind::Piecewise(vec![
                    (value, cond),
                    (integral_definite_unevaluated(f, x, a, b), Expr::boolean(true)),
                ])),
                None,
            ),
            Conds::Separate => (value, Some(cond)),
            Conds::None => (value, None),
        })
    };
    if meijerg != Some(false) && !pw && improper {
        if let Some(r) = try_meijerg() {
            return r;
        }
        // the antiderivative by G-functions is likely nonsensical when the
        // definite method failed: `meijerg=True` gives up, otherwise the
        // other methods go on without it
        if meijerg == Some(true) {
            return (integral_definite_unevaluated(f, x, a, b), None);
        }
        return (integrate_definite_native(f, x, a, b, Some(false)), None);
    }
    if meijerg == Some(true) {
        // the antiderivative by `meijerint_indefinite` first, then the
        // definite method
        if !pw {
            if let Some(big_f) = try_integrate_hint(f, x, Some(true)) {
                if let Some(v) = eval_antiderivative(&big_f, f, x, a, b) {
                    return (v, None);
                }
            }
            if let Some(r) = try_meijerg() {
                return r;
            }
        }
        return (integral_definite_unevaluated(f, x, a, b), None);
    }
    (integrate_definite_native(f, x, a, b, meijerg), None)
}

/// Definite integral of `f` over `[a, b]` through the native rules and the
/// antiderivative.
fn integrate_definite_native(f: &Expr, x: &Expr, a: &Expr, b: &Expr, meijerg: Option<bool>) -> Expr {
    if a == b { return Expr::zero(); }
    if rational_definite::compare(a,b) == Some(std::cmp::Ordering::Greater) {
        return neg(&integrate_definite_native(f,x,b,a,meijerg));
    }
    let unevaluated = || integral_definite_unevaluated(f,x,a,b);
    let checked = match rational_definite::analyze(f,x,a,b) {
        rational_definite::Analysis::Value(value) => return value,
        rational_definite::Analysis::Unknown => return unevaluated(),
        rational_definite::Analysis::Regular(value) => value,
        rational_definite::Analysis::NotRational => f.clone(),
    };
    let f = &checked;
    if let Some(value) = gamma_moments(f, x, a, b) {
        return value;
    }
    try_integrate_hint(f, x, meijerg)
        .and_then(|big_f| eval_antiderivative(&big_f, f, x, a, b))
        .unwrap_or_else(|| integral_definite_unevaluated(f, x, a, b))
}

/// The antiderivative `big_f` of `f` between the bounds, through limits
/// when the bound is infinite or the substitution is singular; `None` when
/// a bound value is not settled.
fn eval_antiderivative(big_f: &Expr, f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    let _ = f;
    if piecewise::has_piecewise(big_f) {
        // the reference splits the Piecewise terms off and evaluates them
        // with `Piecewise._eval_interval`
        return piecewise::eval_split(big_f, x, a, b, &|o| eval_generic(o, x, a, b));
    }
    eval_generic(big_f, x, a, b)
}

/// The generic evaluation of an antiderivative between the bounds.
fn eval_generic(big_f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    let at = |p: &Expr, dir: crate::limits::Dir| -> Expr {
        let v = subs1(big_f, x, p);
        let infinite = matches!(p.kind(), Kind::Const(c) if c.is_infinite());
        if infinite || has_bad(&v) {
            if let Ok(l) = crate::limits::limit(big_f, x, p, dir) {
                return l;
            }
        }
        v
    };
    let fb = at(b, crate::limits::Dir::Minus);
    let fa = at(a, crate::limits::Dir::Plus);
    if has_bad(&fb) || has_bad(&fa) || !settled_bound_value(&fb) || !settled_bound_value(&fa) {
        return None;
    }
    interval_singularities(big_f, x, a, b, sub(&fb, &fa)).map(|value| expand_if_numeric(&value))
}

/// `Expr._eval_interval`: a jump of the antiderivative at a pole or a log
/// branch point strictly inside comparable bounds is removed through the
/// one-sided limits there, so `tan` over `(0, pi)` gives nan as in the
/// reference instead of the endpoint difference. `None` when such a limit
/// is undecided.
fn interval_singularities(big_f: &Expr, x: &Expr, a: &Expr, b: &Expr, value: Expr) -> Option<Expr> {
    let comparable = |e: &Expr| {
        matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity))
            || (crate::functions::is_numeric_value(e) && crate::evalf::eval_f64(e).is_some())
    };
    if !comparable(a) || !comparable(b) {
        return Some(value);
    }
    let mut logs = Vec::new();
    crate::traverse::preorder(big_f, &mut |s| {
        if let Kind::Func(Func::Log, args) = s.kind() {
            if has(s, x) && !logs.contains(&args[0]) {
                logs.push(args[0].clone());
            }
        }
    });
    let cancelled = crate::poly::api::cancel(big_f).unwrap_or_else(|_| big_f.clone());
    let den = numer_denom(&cancelled).1;
    if !has(&den, x) && logs.is_empty() {
        return Some(value);
    }
    let increasing = crate::sets::cmp3(a, b) != Some(std::cmp::Ordering::Greater);
    let domain = if increasing { crate::sets::interval(a, b, false, false) } else { crate::sets::interval(b, a, false, false) };
    let mut points: Vec<Expr> = Vec::new();
    for g in std::iter::once(den).chain(logs) {
        if !has(&g, x) {
            continue;
        }
        let Ok(set) = crate::sets::solveset(&g, x, &domain) else { return Some(value) };
        if crate::sets::is_empty_set(&set) {
            continue;
        }
        // An infinite or unsolved set is not iterable: the reference skips
        // the whole correction (its TypeError is caught).
        let Some(found) = crate::sets::finite_elements(&set) else { return Some(value) };
        for p in found {
            if !points.contains(p) {
                points.push(p.clone());
            }
        }
    }
    let mut value = value;
    for p in points {
        if value.is_constant(Constant::NaN) {
            break;
        }
        if !comparable(&p) {
            continue;
        }
        let less = |u: &Expr, v: &Expr| crate::sets::cmp3(u, v) == Some(std::cmp::Ordering::Less);
        let sign = if less(a, &p) && less(&p, b) {
            Expr::one()
        } else if less(b, &p) && less(&p, a) {
            Expr::int(-1)
        } else {
            continue;
        };
        let right = crate::limits::limit(big_f, x, &p, crate::limits::Dir::Plus).ok()?;
        let left = crate::limits::limit(big_f, x, &p, crate::limits::Dir::Minus).ok()?;
        value = add(&[value, mul2(&sign, &sub(&left, &right))]);
    }
    Some(value)
}

/// `∫ P(x) exp(-α x^q) dx` over `(0, oo)` or `(-oo, oo)` through
/// `∫_0^oo x^p exp(-α x^q) dx = Γ((p+1)/q) / (|q| α^((p+1)/q))`
/// for proved finite signed-real q and the appropriate convergence side; the
/// antiderivatives (erf, incomplete gamma) are out of the rule set's reach.
fn gamma_moments(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    let (value,condition)=gamma_moment_conditions(f,x,a,b)?;
    Some(piecewise(vec![(value,condition),(integral_definite_unevaluated(f,x,a,b),Expr::boolean(true))]))
}

/// Gamma moments together with the sufficient convergence condition, before
/// presentation as Piecewise or the public conds='none'/'separate' forms.
/// The exponential rate must lie provably in the open right half-plane.
/// A false convergence condition never licenses analytic continuation as an
/// integral value, even when the caller suppresses displayed conditions.
pub fn gamma_moment_conditions(f: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Option<(Expr,Expr)> {
    if *b != Expr::oo() {
        return None;
    }
    let whole = *a == Expr::neg_oo();
    if !whole && !a.is_zero() {
        return None;
    }
    let ex = expand(f);
    let terms = match ex.kind() {
        Kind::Add(..) => ex.args(),
        _ => vec![ex.clone()],
    };
    let mut out = Vec::with_capacity(terms.len());
    let mut conditions=Vec::with_capacity(terms.len());
    for t in terms {
        let factors: Vec<Expr> = match t.kind() {
            Kind::Mul(c, facs) => std::iter::once(Expr::num(c.clone()))
                .chain(facs.iter().map(|(b, e)| crate::ops::mul_factor(b, e)))
                .collect(),
            _ => vec![t.clone()],
        };
        let (mut coeff, mut p, mut args) = (Vec::new(), Expr::zero(), Vec::new());
        for g in factors {
            if free(&g, x) {
                coeff.push(g);
            } else if g == *x {
                p = add2(&p, &Expr::one());
            } else {
                match g.kind() {
                    Kind::Pow(base, e) if base == x && free(e,x) => p = add2(&p, e),
                    Kind::Func(Func::Exp, u) => args.push(u[0].clone()),
                    _ => return None,
                }
            }
        }
        let phase=expand(&add(&args));
        let (q,c,offset)=if let Some((q,c))=moment_degree(&phase,x){(q,c,Expr::zero())}
        else{
            // A complex coefficient or a free phase offset can distribute
            // into several additive terms. Group their common degree
            // sparsely: a large x^q must not allocate a dense q-entry array.
            let Kind::Add(..)=phase.kind() else{return None;};
            let mut degree=None;let mut coefficients=Vec::new();let mut offsets=Vec::new();
            for term in phase.args(){
                let (q,c)=moment_degree(&term,x)?;
                if q.is_zero(){offsets.push(c);}
                else{
                    if degree.as_ref().is_some_and(|previous|*previous!=q){return None;}
                    degree=Some(q);coefficients.push(c);
                }
            }
            (degree?,add(&coefficients),add(&offsets))
        };
        if !offset.is_zero(){coeff.push(exp(&offset));}
        let alpha = neg(&c);
        if q.is_zero() || crate::assumptions::is_finite(&q)==Some(false)
            || crate::assumptions::is_finite(&alpha)==Some(false)
            || crate::assumptions::exact_queries(||crate::assumptions::is_positive(&re(&alpha))) != Some(true) {
            return None;
        }
        let positive=crate::assumptions::exact_queries(||crate::assumptions::is_positive(&q))==Some(true);
        if !positive&&crate::assumptions::exact_queries(||crate::assumptions::is_negative(&q))!=Some(true){return None;}
        let absolute_q=if positive{q.clone()}else{neg(&q)};
        let condition=if positive{crate::relational::gt(&re(&p),&Expr::int(-1))}
            else{crate::relational::lt(&re(&p),&Expr::int(-1))};
        if matches!(condition.kind(),Kind::Bool(false)){return None;}
        conditions.push(condition);
        let s = div(&add2(&p, &Expr::one()), &q);
        let half = div(&gamma(&s), &mul2(&absolute_q, &pow(&alpha, &s)));
        let value = if !whole {
            half
        } else {
            if crate::assumptions::exact_queries(||crate::assumptions::is_even(&q))!=Some(true){
                return None;
            }
            if crate::assumptions::exact_queries(||crate::assumptions::is_odd(&p))==Some(true){Expr::zero()}
            else if crate::assumptions::exact_queries(||crate::assumptions::is_even(&p))==Some(true){mul2(&Expr::int(2),&half)}
            else{return None;}
        };
        out.push(mul2(&mul(&coeff), &value));
    }
    // Keep one fallback for the ORIGINAL sum. Splitting a convergent sum into
    // individually divergent integrals would give invalid infinity arithmetic.
    Some((add(&out),crate::logic::and(conditions)))
}

/// A sparse monomial view for Mellin moments. Its degree is an expression,
/// so rational, signed-real symbolic and arbitrarily large integer powers do
/// not pass through a machine-sized polynomial-degree counter. Signs and the
/// real-axis domain are proved separately before using the gamma formula.
fn moment_degree(term:&Expr,x:&Expr)->Option<(Expr,Expr)>{
    if free(term,x){return Some((Expr::zero(),term.clone()));}
    if term==x{return Some((Expr::one(),Expr::one()));}
    match term.kind(){
        Kind::Pow(base,exponent) if base==x&&free(exponent,x)=>Some((exponent.clone(),Expr::one())),
        Kind::Mul(coefficient,factors)=>{
            let mut degree=Expr::zero();let mut coefficients=vec![Expr::num(coefficient.clone())];
            for (base,exponent) in factors{
                if base==x&&free(exponent,x){degree=add2(&degree,exponent);}
                else{
                    let factor=crate::ops::mul_factor(base,exponent);
                    if !free(&factor,x){return None;}
                    coefficients.push(factor);
                }
            }
            Some((degree,mul(&coefficients)))
        }
        _=>None,
    }
}

/// A bound value of the antiderivative is usable only when it is either a
/// signed infinity (a proved divergence) or free of every infinity. A
/// residual infinity inside an expression, such as exp(-oo*s), cos(oo) or
/// oo**(s+1), means the limit was not decided for all parameter values and
/// must not be presented as the value of the integral.
fn settled_bound_value(e: &Expr) -> bool {
    if matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity)) {
        return true;
    }
    !has_infinity(e)
}

fn has_infinity(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(c) => c.is_infinite() || matches!(c, Constant::NaN),
        _ => e.args().iter().any(has_infinity),
    }
}

fn has_bad(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(Constant::NaN | Constant::ComplexInfinity) => true,
        _ => e.args().iter().any(has_bad),
    }
}

fn expand_if_numeric(e: &Expr) -> Expr {
    if crate::functions::is_numeric_value(e) {
        let x = expand(e);
        if crate::traverse::node_count(&x) <= crate::traverse::node_count(e) {
            return x;
        }
    }
    e.clone()
}

/// `Integral(f, x)` placeholder.
pub fn integral_unevaluated(f: &Expr, x: &Expr) -> Expr {
    crate::binders::held(Func::Integral,f,vec![vec![x.clone()]]).expect("valid integration variable")
}

pub fn integral_definite_unevaluated(f:&Expr,x:&Expr,a:&Expr,b:&Expr)->Expr {
    crate::binders::held(Func::Integral,f,vec![vec![x.clone(),a.clone(),b.clone()]]).expect("valid integration limits")
}

/// Antiderivative, or `None` if the integrator gives up.
pub fn try_integrate(f: &Expr, x: &Expr) -> Option<Expr> {
    doit_indefinite(&piecewise::doit_rewrite(f, x, None), x, None)
}

fn native_antiderivative(f: &Expr, x: &Expr) -> Option<Expr> {
    dirac::with_input(f, || {
        let mut seen = FxHashSet::default();
        integ(f, x, 0, &mut seen)
    })
}

/// `Integral._eval_integral(f, x, meijerg=...)`. The native rules stand in
/// for the reference stages before the Meijer G method (polynomials,
/// rational functions, Risch and heurisch); a term they leave is then
/// integrated by `meijerint_indefinite` unless `meijerg=False`.
/// `meijerg=True` integrates only by Meijer G-functions.
pub fn try_integrate_hint(f: &Expr, x: &Expr, meijerg: Option<bool>) -> Option<Expr> {
    if meijerg == Some(true) {
        if free(f, x) {
            return Some(mul2(f, x));
        }
        return meijer_terms(f, x, true);
    }
    let native = native_antiderivative(f, x);
    if native.is_some() || meijerg == Some(false) {
        return native;
    }
    meijer_terms(f, x, false)
}

/// The term loop of `_eval_integral` from the Meijer G stage on.
fn meijer_terms(f: &Expr, x: &Expr, only: bool) -> Option<Expr> {
    // `Piecewise._eval_integral` and the order terms have their own paths
    if matches!(f.kind(), Kind::Piecewise(_)) || has_piecewise_or_minmax(f) {
        return None;
    }
    let terms = crate::simptools::add_args(f);
    let mut parts = Vec::new();
    for t in &terms {
        let (coeff, g) = crate::decompose::as_independent(t, std::slice::from_ref(x), None);
        if !only {
            if let Some(h) = native_antiderivative(&g, x) {
                parts.push(mul2(&coeff, &h));
                continue;
            }
        }
        if let Some(h) = meijerint::meijerint_indefinite(&g, x) {
            parts.push(mul2(&coeff, &h));
            continue;
        }
        if terms.len() == 1 {
            let e = crate::simptools::fu::sincos_to_sum(f);
            let options = crate::expand::ExpandOptions { deep: false, ..Default::default() };
            let e = crate::expand::expand_with_options(&e, &options, &mut crate::expand::ExpandBudget::default()).unwrap_or(e);
            if e.is_add() {
                return meijer_terms(&e, x, only);
            }
        }
        return None;
    }
    Some(add(&parts))
}

#[inline]
fn free(e: &Expr, x: &Expr) -> bool {
    !has(e, x)
}

/// `u == a*x + b` with `a != 0`, `a` and `b` free of `x`.
fn linear(u: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    if free(u, x) {
        return None;
    }
    let a = diff(u, x);
    if !free(&a, x) || a.is_zero() {
        return None;
    }
    let b = expand(&sub(u, &mul2(&a, x)));
    if !free(&b, x) {
        return None;
    }
    Some((a, b))
}

/// Dense coefficients `[c0, c1, ...]` of `e` as a polynomial in `x`.
pub fn poly_coeffs(e: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    if free(e, x) {
        return Some(vec![e.clone()]);
    }
    let ex = expand(e);
    let terms: Vec<Expr> = match ex.kind() {
        Kind::Add(..) => ex.args(),
        _ => vec![ex.clone()],
    };
    let mut coeffs: Vec<Vec<Expr>> = Vec::new();
    for t in terms {
        let (deg, c) = monomial_degree(&t, x)?;
        if deg > 10_000 {
            return None;
        }
        if coeffs.len() <= deg {
            coeffs.resize(deg + 1, Vec::new());
        }
        coeffs[deg].push(c);
    }
    Some(coeffs.into_iter().map(|v| add(&v)).collect())
}

fn monomial_degree(t: &Expr, x: &Expr) -> Option<(usize, Expr)> {
    if free(t, x) {
        return Some((0, t.clone()));
    }
    if t == x {
        return Some((1, Expr::one()));
    }
    match t.kind() {
        Kind::Pow(b, e) if b == x => {
            let k = e.as_i64()?;
            if k < 0 {
                return None;
            }
            Some((k as usize, Expr::one()))
        }
        Kind::Mul(c, facs) => {
            let mut deg = None;
            let mut rest = vec![Expr::num(c.clone())];
            for (b, e) in facs {
                if b == x {
                    let k = e.as_i64()?;
                    if k < 0 || deg.is_some() {
                        return None;
                    }
                    deg = Some(k as usize);
                } else {
                    let f = crate::ops::mul_factor(b, e);
                    if !free(&f, x) {
                        return None;
                    }
                    rest.push(f);
                }
            }
            Some((deg.unwrap_or(0), mul(&rest)))
        }
        _ => None,
    }
}

fn from_coeffs(c: &[Expr], x: &Expr) -> Expr {
    let terms: Vec<Expr> = c
        .iter()
        .enumerate()
        .filter(|(_, k)| !k.is_zero())
        .map(|(i, k)| mul2(k, &x.powi(i as i64)))
        .collect();
    add(&terms)
}

/// Polynomial division with exact coefficients: (quotient, remainder).
fn poly_divmod(n: &[Expr], d: &[Expr]) -> Option<(Vec<Expr>, Vec<Expr>)> {
    let dn = d.len().checked_sub(1)?;
    let lead = d.last()?.clone();
    if lead.is_zero() {
        return None;
    }
    let mut r: Vec<Expr> = n.to_vec();
    if r.len() <= dn {
        return Some((vec![Expr::zero()], r));
    }
    let mut q = vec![Expr::zero(); r.len() - dn];
    for i in (dn..r.len()).rev() {
        let c = expand(&div(&r[i], &lead));
        if c.is_zero() {
            continue;
        }
        q[i - dn] = c.clone();
        for j in 0..=dn {
            r[i - dn + j] = expand(&sub(&r[i - dn + j], &mul2(&c, &d[j])));
        }
    }
    r.truncate(dn);
    if r.is_empty() {
        r.push(Expr::zero());
    }
    Some((q, r))
}

fn trim(mut c: Vec<Expr>) -> Vec<Expr> {
    while c.len() > 1 && c.last().is_some_and(|x| x.is_zero()) {
        c.pop();
    }
    c
}

/// Split a product into (factors free of x, the rest).
fn split_const(f: &Expr, x: &Expr) -> (Expr, Expr) {
    match f.kind() {
        Kind::Mul(..) => {
            let (mut c, mut r) = (Vec::new(), Vec::new());
            for a in f.args() {
                if free(&a, x) {
                    c.push(a)
                } else {
                    r.push(a)
                }
            }
            (mul(&c), mul(&r))
        }
        _ => (Expr::one(), f.clone()),
    }
}

fn integ(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    if depth > MAX_DEPTH {
        return None;
    }
    if free(f, x) {
        return Some(mul2(f, x));
    }
    if f == x {
        return Some(div(&x.powi(2), &Expr::int(2)));
    }
    // guard against substitution loops
    if !seen.insert(f.clone()) {
        return None;
    }
    let r = integ_inner(f, x, depth, seen);
    seen.remove(f);
    r
}

fn integ_inner(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    if let Some(value) = dirac::antiderivative(f, x) {
        return Some(value);
    }
    // polynomials: integrate the expanded form term by term
    if let Some(c) = poly_coeffs(f, x) {
        let terms: Vec<Expr> = c
            .iter()
            .enumerate()
            .filter(|(_, k)| !k.is_zero())
            .map(|(i, k)| mul(&[k.clone(), x.powi(i as i64 + 1), Expr::rational(1, i as i64 + 1)]))
            .collect();
        return Some(add(&terms));
    }
    match f.kind() {
        Kind::Add(..) => {
            let mut v = Vec::new();
            for t in f.args() {
                v.push(integ(&t, x, depth + 1, seen)?);
            }
            return Some(add(&v));
        }
        Kind::Mul(..) => {
            let (c, rest) = split_const(f, x);
            if !c.is_one() {
                let r = integ(&rest, x, depth + 1, seen)?;
                // c * Piecewise((F/a, Ne(a, 0)), ...) with c a multiple of a:
                // the a = 0 branch is a constant, drop it (SymPy does too)
                if let Kind::Piecewise(p) = r.kind() {
                    if let Kind::Rel(crate::expr::RelOp::Ne, a, z) = p[0].1.kind() {
                        if z.is_zero() && !has(&div(&c, a), a) {
                            return Some(mul2(&c, &p[0].0));
                        }
                    }
                }
                return Some(mul2(&c, &r));
            }
        }
        _ => {}
    }
    if let Some(r) = table(f, x) {
        return Some(r);
    }
    if let Some(r) = crate::special::integrate::table(f, x) {
        return Some(r);
    }
    if let Some(r) = rational(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = trig_powers(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = cyclic(f, x) {
        return Some(r);
    }
    if let Some(r) = usub(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = poly_exp_trig(f, x) {
        return Some(r);
    }
    if let Some(r) = parts(f, x, depth, seen) {
        return Some(r);
    }
    if let Some(r) = special(f, x) {
        return Some(r);
    }
    // rewrite tan/cot/sec/csc products into sin/cos
    let rw = rewrite_trig(f);
    if rw != *f {
        if let Some(r) = integ(&rw, x, depth + 1, seen) {
            return Some(r);
        }
    }
    let ex = expand(f);
    if ex != *f {
        if let Some(r) = integ(&ex, x, depth + 1, seen) {
            return Some(r);
        }
    }
    // products and powers of sin/cos: turn them into sums of single terms
    if has_head(f, &[Func::Sin, Func::Cos]) {
        let lin = expand(&crate::simptools::fu::tr8(f, false));
        if lin != *f && lin != ex {
            if let Some(r) = integ(&lin, x, depth + 1, seen) {
                return Some(r);
            }
        }
    }
    // exponentials times sin/cos whose cyclic denominator vanishes, as in
    // exp(I*x)*sin(x): integrate the exponential form
    if has_head(f, &[Func::Exp]) && has_head(f, &[Func::Sin, Func::Cos]) {
        let ef = expand(&crate::special::rewrite::rewrite(f, &["sin".to_string(), "cos".to_string()], "exp", true));
        if ef != *f && ef != ex {
            if let Some(r) = integ(&ef, x, depth + 1, seen) {
                return Some(r);
            }
        }
    }
    // hyperbolic functions: integrate the exponential form
    if has_head(f, &[Func::Sinh, Func::Cosh]) {
        let ef = expand(&crate::simptools::trigsimp::rewrite_exp(f));
        if ef != *f {
            return integ(&ef, x, depth + 1, seen);
        }
    }
    None
}

fn has_head(e: &Expr, heads: &[Func]) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(g, _) = s.kind() {
            found |= heads.contains(g);
        }
    });
    found
}

/// `P(x) * exp(a x + c) * {cos, sin}(b x + d)` in closed form, any factor
/// but `P` optional: the antiderivative is `exp(..) * (A cos + B sin)` with
/// polynomials `A`, `B` of the degree of `P`, solved from the top degree down
/// (`a A + b B + A' = P_cos`, `a B - b A + B' = P_sin`).
fn poly_exp_trig(f: &Expr, x: &Expr) -> Option<Expr> {
    let factors = match f.kind() {
        Kind::Mul(..) => f.args(),
        _ => vec![f.clone()],
    };
    let (mut ex, mut tr): (Option<Expr>, Option<(bool, Expr)>) = (None, None);
    let mut poly = Vec::new();
    for t in factors {
        match t.kind() {
            Kind::Func(Func::Exp, a) if ex.is_none() && linear(&a[0], x).is_some() => ex = Some(a[0].clone()),
            Kind::Func(g @ (Func::Sin | Func::Cos), a) if tr.is_none() && linear(&a[0], x).is_some() => {
                tr = Some((*g == Func::Sin, a[0].clone()))
            }
            _ => poly.push(t),
        }
    }
    // a lone exp or trig factor is the table's job, and two factors without
    // a polynomial are the cyclic rule's
    let p = poly_coeffs(&mul(&poly), x)?;
    let p = trim(p);
    if p.len() < 2 && (ex.is_none() || tr.is_none()) {
        return None;
    }
    let (a, ea) = match &ex {
        Some(u) => (linear(u, x)?.0, exp(u)),
        None => (Expr::zero(), Expr::one()),
    };
    let (b, is_sin, c, s) = match &tr {
        Some((is_sin, u)) => (linear(u, x)?.0, *is_sin, cos(u), sin(u)),
        None => (Expr::zero(), false, Expr::one(), Expr::zero()),
    };
    let det = expand(&add2(&a.powi(2), &b.powi(2)));
    if det.is_zero() {
        return None;
    }
    let n = p.len();
    let (mut ac, mut bc) = (vec![Expr::zero(); n + 1], vec![Expr::zero(); n + 1]);
    for k in (0..n).rev() {
        let k1 = Expr::int(k as i64 + 1);
        let (pc, ps) = if is_sin { (Expr::zero(), p[k].clone()) } else { (p[k].clone(), Expr::zero()) };
        let r1 = sub(&pc, &mul2(&k1, &ac[k + 1]));
        let r2 = sub(&ps, &mul2(&k1, &bc[k + 1]));
        ac[k] = expand(&div(&sub(&mul2(&a, &r1), &mul2(&b, &r2)), &det));
        bc[k] = expand(&div(&add2(&mul2(&b, &r1), &mul2(&a, &r2)), &det));
    }
    let pa = from_coeffs(&ac[..n], x);
    let pb = from_coeffs(&bc[..n], x);
    if tr.is_none() {
        // the Risch form: the polynomial over a common denominator times the
        // exponential
        let (mut num, mut den) = crate::poly::api::as_numer_denom(&crate::poly::api::together(&pa));
        if !den.is_number() {
            if let Ok(c) = crate::poly::api::cancel(&pa) {
                (num, den) = crate::poly::api::as_numer_denom(&c);
            }
        }
        let generic = mul(&[num, ea, pow(&den, &Expr::minus_one())]);
        // a coefficient of x that may vanish makes the generic antiderivative
        // undefined there: the value at zero is kept in a Piecewise like the
        // reference Risch integrator
        if den.is_number() || crate::assumptions::fact(&a, "zero") == Some(false) {
            return Some(generic);
        }
        let u = ex.as_ref().expect("exponential factor");
        let f0 = mul2(&mul(&poly), &exp(&expand(&sub(u, &mul2(&a, x)))));
        let special = native_antiderivative(&f0, x).unwrap_or_else(|| integral_unevaluated(&f0, x));
        let (_, monic) = crate::poly::api::add_primitive(&den);
        return Some(piecewise(vec![(generic, ne(&monic, &Expr::zero())), (special, Expr::boolean(true))]));
    }
    let generic = expand(&mul2(&ea, &add2(&mul2(&pa, &c), &mul2(&pb, &s))));
    Some(degenerate_cases(f, x, &det, generic))
}

/// `heurisch_wrapper`: the parameter values that make the denominator `den`
/// of the generic antiderivative vanish are integrated separately and kept
/// in a Piecewise, the generic case last (first, guarded by `Ne`, when there
/// is a single special case).
fn degenerate_cases(f: &Expr, x: &Expr, den: &Expr, generic: Expr) -> Expr {
    // `solve([d], dict=True)`: each factor is solved for its first symbol
    let zeros = |d: &Expr| -> Vec<Vec<(Expr, Expr)>> {
        let factors = match crate::poly::api::factor_list(d) {
            Ok((_, fl)) => fl.into_iter().map(|(g, _)| g).collect(),
            Err(_) => vec![d.clone()],
        };
        let mut out: Vec<Vec<(Expr, Expr)>> = Vec::new();
        for g in factors {
            let mut syms: Vec<Expr> = crate::traverse::free_symbols(&g).into_iter().filter(|s| s != x).collect();
            syms.sort_by(crate::sortkey::cmp_default);
            for s in syms {
                if let Some(v) = crate::solve::solve(&g, &s) {
                    for v in v {
                        let sol = vec![(s.clone(), v)];
                        if !out.contains(&sol) {
                            out.push(sol);
                        }
                    }
                    break;
                }
            }
        }
        out
    };
    if free(den, x) && den.is_number() {
        return generic;
    }
    let mut sols = zeros(den);
    if sols.is_empty() {
        return generic;
    }
    // the poles of the integrand itself are not special cases
    let (_, fden) = numer_denom(f);
    if !fden.is_number() {
        let poles = zeros(&fden);
        sols.retain(|s| !poles.contains(s));
        if sols.is_empty() {
            return generic;
        }
    }
    if sols.len() > 1 {
        let eqs: Vec<Expr> = sols.iter().flatten().map(|(k, v)| sub(k, v)).collect();
        let mut vars: Vec<Expr> = Vec::new();
        for e in &eqs {
            for s in crate::traverse::free_symbols(e) {
                if s != *x && !vars.contains(&s) {
                    vars.push(s);
                }
            }
        }
        vars.sort_by(crate::sortkey::cmp_default);
        match crate::solve::solve_system(&eqs, &vars) {
            Some(crate::solve::SystemSolution::Map(m)) if !m.is_empty() => sols.insert(0, m),
            Some(crate::solve::SystemSolution::NoSolution) => {}
            // the joint case could not be decided
            Some(_) | None => return generic,
        }
    }
    let mut pairs = Vec::new();
    for sol in &sols {
        let fs = sol.iter().fold(f.clone(), |e, (k, v)| subs1(&e, k, v));
        let value = native_antiderivative(&fs, x).unwrap_or_else(|| integral_unevaluated(&fs, x));
        let cond = crate::logic::and(sol.iter().map(|(k, v)| crate::relational::eq(k, v)).collect());
        pairs.push((value, cond));
    }
    if pairs.len() == 1 {
        let sol = &sols[0];
        let cond = crate::logic::or(sol.iter().map(|(k, v)| ne(k, v)).collect());
        let special = pairs.pop().unwrap().0;
        return piecewise(vec![(generic, cond), (special, Expr::boolean(true))]);
    }
    pairs.push((generic, Expr::boolean(true)));
    piecewise(pairs)
}

fn rewrite_trig(e: &Expr) -> Expr {
    match e.kind() {
        Kind::Func(Func::Tan, a) => div(&sin(&a[0]), &cos(&a[0])),
        Kind::Func(Func::Cot, a) => div(&cos(&a[0]), &sin(&a[0])),
        Kind::Func(Func::Sec, a) => div(&Expr::one(), &cos(&a[0])),
        Kind::Func(Func::Csc, a) => div(&Expr::one(), &sin(&a[0])),
        _ if e.is_atom() => e.clone(),
        _ => crate::traverse::map_children(e, &mut rewrite_trig),
    }
}

/// Elementary antiderivatives of `F(a*x + b)`.
fn table(f: &Expr, x: &Expr) -> Option<Expr> {
    // power rule: (a x + b)^n
    if let Kind::Pow(base, n) = f.kind() {
        if free(n, x) {
            if let Some((a, _)) = linear(base, x) {
                if n.is_minus_one() {
                    return Some(div(&log(base), &a));
                }
                if n.is_number() {
                    let n1 = add2(n, &Expr::one());
                    return Some(div(&pow(base, &n1), &mul2(&a, &n1)));
                }
                // symbolic exponent: Piecewise like SymPy
                let n1 = add2(n, &Expr::one());
                return Some(piecewise(vec![
                    (div(&pow(base, &n1), &mul2(&a, &n1)), ne(n, &Expr::minus_one())),
                    (div(&log(base), &a), Expr::boolean(true)),
                ]));
            }
        }
        // c^(a x + b)
        if free(base, x) {
            if let Some((a, _)) = linear(n, x) {
                let lb = log(base);
                let r = div(f, &mul2(&a, &lb));
                if lb.is_number() || crate::functions::is_numeric_value(&lb) {
                    return Some(r);
                }
                return Some(piecewise(vec![(r, ne(&lb, &Expr::zero())), (x.clone(), Expr::boolean(true))]));
            }
        }
        // 1/sqrt(c - d x^2), 1/sqrt(x^2 + c), sqrt(c - d x^2), ...
        if let Some(r) = sqrt_quadratic(base, n, x) {
            return Some(r);
        }
        if let Some(r) = inverse_quadratic(base, n, x) {
            return Some(r);
        }
        // sec^2, csc^2, 1/cos^2, 1/sin^2 (SymPy writes tan as sin/cos)
        if let (Kind::Func(g, a), Some(k)) = (base.kind(), n.as_i64()) {
            if let Some((da, _)) = linear(&a[0], x) {
                let u = &a[0];
                let r = match (g, k) {
                    (Func::Sec, 2) | (Func::Cos, -2) => Some(div(&sin(u), &cos(u))),
                    (Func::Csc, 2) | (Func::Sin, -2) => Some(neg(&div(&cos(u), &sin(u)))),
                    (Func::Tan, 2) => Some(sub(&div(&sin(u), &cos(u)), u)),
                    (Func::Cot, 2) => Some(sub(&neg(&div(&cos(u), &sin(u))), u)),
                    _ => None,
                };
                if let Some(r) = r {
                    return Some(div(&r, &da));
                }
            }
        }
    }
    let Kind::Func(func, args) = f.kind() else {
        return None;
    };
    if args.len() != 1 {
        return None;
    }
    let u = &args[0];
    let (a, _) = linear(u, x)?;
    let sym = numeric_or_piecewise_guard(&a);
    let res = match func {
        Func::Exp => f.clone(),
        Func::Sin => neg(&cos(u)),
        Func::Cos => sin(u),
        Func::Tan => neg(&log(&cos(u))),
        Func::Cot => log(&sin(u)),
        Func::Sec => log(&add2(&tan(u), &sec(u))),
        Func::Csc => neg(&log(&add2(&cot(u), &csc(u)))),
        Func::Sinh => cosh(u),
        Func::Cosh => sinh(u),
        Func::Tanh => sub(u, &log(&add2(&tanh(u), &Expr::one()))),
        Func::Log => sub(&mul2(u, &log(u)), u),
        Func::Atan => sub(&mul2(u, &atan(u)), &div(&log(&add2(&u.powi(2), &Expr::one())), &Expr::int(2))),
        Func::Asin => add2(&mul2(u, &asin(u)), &sqrt(&sub(&Expr::one(), &u.powi(2)))),
        Func::Acos => sub(&mul2(u, &acos(u)), &sqrt(&sub(&Expr::one(), &u.powi(2)))),
        Func::Asinh => sub(&mul2(u, &asinh(u)), &sqrt(&add2(&u.powi(2), &Expr::one()))),
        Func::Atanh => add2(&mul2(u, &atanh(u)), &div(&log(&sub(&u.powi(2), &Expr::one())), &Expr::int(2))),
        Func::Erf => add2(
            &mul2(u, &erf(u)),
            &div(&exp(&neg(&u.powi(2))), &sqrt(&Expr::pi())),
        ),
        _ => return None,
    };
    let r = div(&res, &a);
    if sym {
        // exp(a*x) with symbolic a: SymPy returns a Piecewise
        let at_zero = mul2(&subs1(f, &a, &Expr::zero()), x);
        return Some(piecewise(vec![(r, ne(&a, &Expr::zero())), (at_zero, Expr::boolean(true))]));
    }
    Some(r)
}

fn numeric_or_piecewise_guard(a: &Expr) -> bool {
    !crate::functions::is_numeric_value(a)
}

/// Integrals involving sqrt of a quadratic `c + d x^2` (no linear term).
fn sqrt_quadratic(base: &Expr, n: &Expr, x: &Expr) -> Option<Expr> {
    let half = Expr::rational(1, 2);
    let mhalf = Expr::rational(-1, 2);
    if *n != half && *n != mhalf {
        return None;
    }
    let c = poly_coeffs(base, x)?;
    if c.len() != 3 || !c[1].is_zero() {
        return None;
    }
    let (c0, d) = (&c[0], &c[2]);
    let (c0v, dv) = (crate::evalf::eval_f64(c0)?, crate::evalf::eval_f64(d)?);
    if *n == mhalf {
        if c0v > 0.0 && dv < 0.0 {
            // 1/sqrt(c - k x^2) = asin(sqrt(k/c) x)/sqrt(k)
            let k = neg(d);
            return Some(div(&asin(&mul2(&sqrt(&div(&k, c0)), x)), &sqrt(&k)));
        }
        if c0v > 0.0 && dv > 0.0 {
            return Some(div(&asinh(&mul2(&sqrt(&div(d, c0)), x)), &sqrt(d)));
        }
        if c0v < 0.0 && dv > 0.0 {
            return Some(div(&acosh(&mul2(&sqrt(&div(d, &neg(c0))), x)), &sqrt(d)));
        }
        return None;
    }
    // sqrt(c + d x^2)
    let s = sqrt(base);
    if c0v > 0.0 && dv < 0.0 {
        let k = neg(d);
        let t = mul2(&sqrt(&div(&k, c0)), x);
        return Some(add2(
            &div(&mul2(x, &s), &Expr::int(2)),
            &div(&mul2(c0, &asin(&t)), &mul2(&Expr::int(2), &sqrt(&k))),
        ));
    }
    if c0v > 0.0 && dv > 0.0 {
        let t = mul2(&sqrt(&div(d, c0)), x);
        return Some(add2(
            &div(&mul2(x, &s), &Expr::int(2)),
            &div(&mul2(c0, &asinh(&t)), &mul2(&Expr::int(2), &sqrt(d))),
        ));
    }
    None
}

/// 1/(c + d x^2) -> atan / atanh forms.
fn inverse_quadratic(base: &Expr, n: &Expr, x: &Expr) -> Option<Expr> {
    if !n.is_minus_one() {
        return None;
    }
    let c = poly_coeffs(base, x)?;
    if c.len() != 3 || !c[1].is_zero() {
        return None;
    }
    let (c0, d) = (&c[0], &c[2]);
    let (c0v, dv) = (crate::evalf::eval_f64(c0)?, crate::evalf::eval_f64(d)?);
    if c0v * dv > 0.0 {
        // atan(sqrt(d/c) x) / sqrt(c d)
        let t = mul2(&sqrt(&div(d, c0)), x);
        let r = div(&atan(&t), &sqrt(&mul2(c0, d)));
        return Some(if c0v < 0.0 { neg(&r) } else { r });
    }
    None
}

/// Rational functions with a denominator of degree <= 2, or a power of a
/// linear polynomial.
fn rational(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let (num, den) = numer_denom(f);
    if free(&den, x) {
        return None;
    }
    let nc = poly_coeffs(&num, x)?;
    // denominator a power of a linear polynomial: substitute u = a x + b
    if let Kind::Pow(b, k) = den.kind() {
        if let (Some((a, c)), Some(_)) = (linear(b, x), k.as_i64()) {
            let u = Expr::symbol("_u");
            let xu = div(&sub(&u, &c), &a);
            let numu = expand(&subs1(&num, x, &xu));
            let g = expand(&div(&numu, &pow(&u, k)));
            let gi = integ(&g, &u, depth + 1, seen)?;
            return Some(div(&subs1(&gi, &u, b), &a));
        }
    }
    let dc = trim(poly_coeffs(&den, x)?);
    let nc = trim(nc);
    if dc.len() < 2 {
        return None;
    }
    let (q, r) = poly_divmod(&nc, &dc)?;
    let mut parts = Vec::new();
    let qp = from_coeffs(&q, x);
    if !qp.is_zero() {
        parts.push(integ(&qp, x, depth + 1, seen)?);
    }
    let r = trim(r);
    if r.len() == 1 && r[0].is_zero() {
        return Some(add(&parts));
    }
    match dc.len() {
        2 => {
            // r0 / (d1 x + d0)
            parts.push(div(&mul2(&r[0], &log(&from_coeffs(&dc, x))), &dc[1]));
        }
        3 => {
            let (c, b, a) = (&dc[0], &dc[1], &dc[2]);
            let p = r.get(1).cloned().unwrap_or_else(Expr::zero);
            let qq = r[0].clone();
            let disc = expand(&sub(&b.powi(2), &mul(&[Expr::int(4), a.clone(), c.clone()])));
            let dv = crate::evalf::eval_f64(&disc)?;
            let d_expr = from_coeffs(&dc, x);
            if dv.abs() < 1e-300 {
                return None;
            }
            if dv < 0.0 {
                // (p x + q)/(a x^2 + b x + c)
                //  = p/(2a) * log(den) + (q - p b/(2a)) * 2/sqrt(-disc) * atan((2 a x + b)/sqrt(-disc))
                let s = sqrt(&neg(&disc));
                let two_a = mul2(&Expr::int(2), a);
                let log_part = mul2(&div(&p, &two_a), &log(&d_expr));
                let k = sub(&qq, &div(&mul2(&p, b), &two_a));
                let atan_part = mul(&[
                    k,
                    Expr::int(2),
                    pow(&s, &Expr::minus_one()),
                    atan(&div(&add2(&mul2(&two_a, x), b), &s)),
                ]);
                if !p.is_zero() {
                    parts.push(log_part);
                }
                if !atan_part.is_zero() {
                    parts.push(expand(&atan_part));
                }
            } else {
                // two real roots r1, r2: A/(x - r1) + B/(x - r2)
                let s = sqrt(&disc);
                let two_a = mul2(&Expr::int(2), a);
                let r1 = div(&sub(&neg(b), &s), &two_a);
                let r2 = div(&add2(&neg(b), &s), &two_a);
                let numer = from_coeffs(&r, x);
                // A = num(r1) / (a (r1 - r2)), B = num(r2) / (a (r2 - r1))
                let a1 = expand(&div(&subs1(&numer, x, &r1), &mul2(a, &sub(&r1, &r2))));
                let b1 = expand(&div(&subs1(&numer, x, &r2), &mul2(a, &sub(&r2, &r1))));
                parts.push(mul2(&a1, &log(&sub(x, &r1))));
                parts.push(mul2(&b1, &log(&sub(x, &r2))));
            }
        }
        _ => {
            let rem = from_coeffs(&r, x);
            let piece = match quadratic_power(&rem, &den, x, depth, seen) {
                Some(v) => v,
                None => {
                    // partial fractions over Q, each piece has a linear or
                    // quadratic denominator (or a power of one)
                    let pf = crate::poly::api::apart(&div(&rem, &den), Some(x)).ok()?;
                    if !pf.is_add() {
                        return None;
                    }
                    let mut v = Vec::new();
                    for t in pf.args() {
                        v.push(integ(&t, x, depth + 1, seen)?);
                    }
                    add(&v)
                }
            };
            parts.push(piece);
        }
    }
    Some(add(&parts))
}

/// `(p x + q) / Q^n` with `Q = a x^2 + b x + c`, `n >= 2`, by the reduction
/// `I_n = (2a x + b) / ((n-1) D Q^(n-1)) + 2 (2n-3) a / ((n-1) D) I_(n-1)`,
/// `D = 4ac - b^2`.
fn quadratic_power(num: &Expr, den: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let Kind::Pow(q, k) = den.kind() else { return None };
    let n = k.as_i64().filter(|n| *n >= 2)?;
    let qc = trim(poly_coeffs(q, x)?);
    let pc = trim(poly_coeffs(num, x)?);
    if qc.len() != 3 || pc.len() > 2 {
        return None;
    }
    let (c, b, a) = (&qc[0], &qc[1], &qc[2]);
    let dd = expand(&sub(&mul(&[Expr::int(4), a.clone(), c.clone()]), &b.powi(2)));
    if dd.is_zero() {
        return None;
    }
    let p = pc.get(1).cloned().unwrap_or_else(Expr::zero);
    let two_a = mul2(&Expr::int(2), a);
    // num = p/(2a) Q' + (q0 - p b/(2a))
    let k0 = sub(&pc[0], &div(&mul2(&p, b), &two_a));
    let mut i = integ(&div(&Expr::one(), q), x, depth + 1, seen)?;
    let lin = add2(&mul2(&two_a, x), b);
    for m in 2..=n {
        let m1 = Expr::int(m - 1);
        let head = div(&lin, &mul(&[m1.clone(), dd.clone(), pow(q, &Expr::int(m - 1))]));
        let tail = mul2(&div(&mul(&[Expr::int(2), Expr::int(2 * m - 3), a.clone()]), &mul2(&m1, &dd)), &i);
        i = add2(&head, &tail);
    }
    let log_part = mul2(&div(&p, &mul2(&two_a, &Expr::int(1 - n))), &pow(q, &Expr::int(1 - n)));
    Some(add2(&log_part, &expand(&mul2(&k0, &i))))
}

/// sin(u)^m cos(u)^n with u linear.
fn trig_powers(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let facs: Vec<(Expr, Expr)> = match f.kind() {
        Kind::Mul(c, facs) if c.is_one() => facs.clone(),
        Kind::Pow(b, e) => vec![(b.clone(), e.clone())],
        _ => return None,
    };
    let mut arg: Option<Expr> = None;
    let (mut m, mut n) = (0i64, 0i64);
    for (b, e) in &facs {
        let k = e.as_i64()?;
        match b.kind() {
            Kind::Func(Func::Sin, a) | Kind::Func(Func::Cos, a) => {
                if arg.as_ref().is_some_and(|u| *u != a[0]) {
                    return None;
                }
                arg = Some(a[0].clone());
                if b.is_func(&Func::Sin) {
                    m += k
                } else {
                    n += k
                }
            }
            _ => return None,
        }
    }
    let u = arg?;
    let (a, _) = linear(&u, x)?;
    if m < 0 || n < 0 {
        return None;
    }
    let s = sin(&u);
    let c = cos(&u);
    let t = Expr::symbol("_t");
    if n % 2 == 1 {
        let g = mul(&[pow(&sub(&Expr::one(), &t.powi(2)), &Expr::int((n - 1) / 2)), t.powi(m)]);
        let gi = integ(&expand(&g), &t, depth + 1, seen)?;
        return Some(div(&subs1(&gi, &t, &s), &a));
    }
    if m % 2 == 1 {
        // sin^m cos^n = -(1 - t^2)^((m-1)/2) t^n dt, t = cos u
        let g = mul(&[
            Expr::minus_one(),
            pow(&sub(&Expr::one(), &t.powi(2)), &Expr::int((m - 1) / 2)),
            t.powi(n),
        ]);
        let gi = integ(&expand(&g), &t, depth + 1, seen)?;
        return Some(div(&subs1(&gi, &t, &c), &a));
    }
    // both even: reduction formula on the larger power
    if m >= 2 {
        // ∫ sin^m cos^n = -sin^(m-1) cos^(n+1)/(m+n) + (m-1)/(m+n) ∫ sin^(m-2) cos^n
        let mn = Expr::int(m + n);
        let first = neg(&div(&mul(&[s.powi(m - 1), c.powi(n + 1)]), &mn));
        let rest_f = mul(&[s.powi(m - 2), c.powi(n)]);
        let rest = integ(&rest_f, x, depth + 1, seen)?;
        return Some(add2(&div(&first, &a), &mul2(&Expr::rational(m - 1, m + n), &rest)));
    }
    if n >= 2 {
        let mn = Expr::int(m + n);
        let first = div(&mul(&[s.powi(m + 1), c.powi(n - 1)]), &mn);
        let rest_f = mul(&[s.powi(m), c.powi(n - 2)]);
        let rest = integ(&rest_f, x, depth + 1, seen)?;
        return Some(add2(&div(&first, &a), &mul2(&Expr::rational(n - 1, m + n), &rest)));
    }
    None
}

/// exp(a x) * sin(b x) and exp(a x) * cos(b x).
fn cyclic(f: &Expr, x: &Expr) -> Option<Expr> {
    let Kind::Mul(c, facs) = f.kind() else { return None };
    if !c.is_one() || facs.len() != 2 {
        return None;
    }
    let mut e_arg = None;
    let mut t: Option<(bool, Expr)> = None;
    for (b, e) in facs {
        if !e.is_one() {
            return None;
        }
        match b.kind() {
            Kind::Func(Func::Exp, a) => e_arg = Some(a[0].clone()),
            Kind::Func(Func::Sin, a) => t = Some((true, a[0].clone())),
            Kind::Func(Func::Cos, a) => t = Some((false, a[0].clone())),
            _ => return None,
        }
    }
    let (ea, (is_sin, ta)) = (e_arg?, t?);
    let (a, _) = linear(&ea, x)?;
    let (b, _) = linear(&ta, x)?;
    let den = expand(&add2(&a.powi(2), &b.powi(2)));
    if den.is_zero() {
        return None;
    }
    let ex = exp(&ea);
    let (s, co) = (sin(&ta), cos(&ta));
    let inner = if is_sin {
        sub(&mul2(&a, &s), &mul2(&b, &co))
    } else {
        add2(&mul2(&a, &co), &mul2(&b, &s))
    };
    Some(degenerate_cases(f, x, &den, expand(&div(&mul2(&ex, &inner), &den))))
}

/// Candidate substitutions: function arguments, power bases and exponents,
/// and the functions themselves.
fn candidates(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    if free(e, x) || e == x {
        return;
    }
    match e.kind() {
        Kind::Func(_, args) => {
            out.push(e.clone());
            for a in args {
                if a != x {
                    out.push(a.clone());
                }
                candidates(a, x, out);
            }
        }
        Kind::Pow(b, n) => {
            if b != x {
                out.push(b.clone());
            }
            if !free(n, x) {
                out.push(n.clone());
            }
            if n.is_rational_number() && !n.is_integer_number() {
                out.push(e.clone());
            }
            candidates(b, x, out);
            candidates(n, x, out);
        }
        _ => {
            for a in e.args() {
                candidates(&a, x, out);
            }
        }
    }
}

fn usub(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let mut cands = Vec::new();
    candidates(f, x, &mut cands);
    let mut uniq = FxHashSet::default();
    cands.retain(|c| uniq.insert(c.clone()));
    // prefer bigger sub-expressions first, then SymPy's canonical order
    cands.sort_by_cached_key(|c| (std::cmp::Reverse(crate::traverse::node_count(c)), crate::sortkey::sort_key(c)));
    let usym = Expr::symbol("_u");
    for u in cands {
        if linear(&u, x).is_some() && !matches!(u.kind(), Kind::Func(..)) {
            continue;
        }
        let du = diff(&u, x);
        if du.is_zero() {
            continue;
        }
        let g = div(f, &du);
        let gu = subs1(&g, &u, &usym);
        let gu = if free(&gu, x) {
            gu
        } else {
            // try harder: cancel common factors by expanding the ratio
            let alt = subs1(&expand(&g), &u, &usym);
            if !free(&alt, x) {
                continue;
            }
            alt
        };
        if let Some(r) = integ(&gu, &usym, depth + 1, seen) {
            return Some(subs1(&r, &usym, &u));
        }
    }
    None
}

fn is_poly_factor(e: &Expr, x: &Expr) -> bool {
    poly_coeffs(e, x).is_some() && !free(e, x)
}

/// Integration by parts, choosing `u` by the LIATE rule.
fn parts(f: &Expr, x: &Expr, depth: u32, seen: &mut FxHashSet<Expr>) -> Option<Expr> {
    let factors: Vec<Expr> = match f.kind() {
        Kind::Mul(..) => f.args(),
        Kind::Func(..) => vec![f.clone()],
        Kind::Pow(b, n) if n.is_integer_number() && matches!(b.kind(), Kind::Func(Func::Log, _)) => {
            vec![f.clone()]
        }
        _ => return None,
    };
    // rank: Log, Inverse trig, Algebraic, Trig, Exp
    fn rank(e: &Expr, x: &Expr) -> u8 {
        match e.kind() {
            Kind::Func(Func::Log, _) => 0,
            Kind::Pow(b, _) if matches!(b.kind(), Kind::Func(Func::Log, _)) => 0,
            Kind::Func(Func::Atan | Func::Asin | Func::Acos | Func::Asinh | Func::Atanh | Func::Acot, _) => 1,
            _ if is_poly_factor(e, x) => 2,
            Kind::Func(Func::Sin | Func::Cos | Func::Sinh | Func::Cosh, _) => 3,
            Kind::Func(Func::Exp, _) => 4,
            Kind::Pow(b, _) if !free(&diff(e, x), x) && free(b, x) => 4,
            _ => 5,
        }
    }
    let mut idx: Vec<usize> = (0..factors.len()).collect();
    idx.sort_by_key(|&i| rank(&factors[i], x));
    for &i in idx.iter().take(2) {
        let u = &factors[i];
        if rank(u, x) > 2 {
            continue;
        }
        let dv: Vec<Expr> = factors.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, e)| e.clone()).collect();
        let dv = mul(&dv);
        if rank(u, x) == 2 && dv.is_one() {
            continue;
        }
        let v = integ(&dv, x, depth + 1, seen)?;
        let du = diff(u, x);
        let rest = expand(&mul2(&v, &du));
        let ri = integ(&rest, x, depth + 1, seen)?;
        let r = sub(&mul2(u, &v), &ri);
        // polynomial * exp(..): SymPy factors the exponential out
        if let Kind::Func(Func::Exp, _) = dv.kind() {
            if rank(u, x) == 2 {
                let q = expand(&div(&r, &dv));
                if poly_coeffs(&q, x).is_some() {
                    return Some(mul2(&q, &dv));
                }
            }
        }
        return Some(r);
    }
    None
}

/// Non-elementary special cases.
fn special(f: &Expr, x: &Expr) -> Option<Expr> {
    // exp(-c x^2) -> sqrt(pi) erf(sqrt(c) x) / (2 sqrt(c))
    if let Kind::Func(Func::Exp, a) = f.kind() {
        let c = poly_coeffs(&a[0], x)?;
        if c.len() == 3 && c[0].is_zero() && c[1].is_zero() {
            let k = neg(&c[2]);
            if crate::evalf::eval_f64(&k).is_some_and(|v| v > 0.0) {
                let sk = sqrt(&k);
                return Some(div(&mul2(&sqrt(&Expr::pi()), &erf(&mul2(&sk, x))), &mul2(&Expr::int(2), &sk)));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    fn check(f: &str) {
        let x = Expr::symbol("x");
        let fe = parse_str(f).unwrap();
        let r = try_integrate(&fe, &x).unwrap_or_else(|| panic!("failed to integrate {f}"));
        // verify numerically that d/dx r == f at a few points
        let d = diff(&r, &x);
        for v in [0.3, 0.55, 0.7] {
            let a = eval_c64(&subs(&d, &[(x.clone(), Expr::float(v))]));
            let b = eval_c64(&subs(&fe, &[(x.clone(), Expr::float(v))]));
            match (a, b) {
                (Some(a), Some(b)) => assert!((a - b).norm() < 1e-8 * (1.0 + b.norm()), "{f}: got {r}"),
                _ => panic!("{f}: cannot evaluate {r}"),
            }
        }
    }

    #[test]
    fn basic() {
        for f in [
            "x**3 + 2*x", "1/x", "exp(2*x)", "sin(3*x+1)", "x*exp(x)", "x**2*sin(x)", "log(x)",
            "x*log(x)", "atan(x)", "1/(1+x**2)", "1/sqrt(1-x**2)", "sin(x)**2", "cos(x)**3",
            "sin(x)*cos(x)", "exp(x)*sin(x)", "2*x*exp(x**2)", "x/(x**2+1)", "cos(x)*exp(sin(x))",
            "1/(x**2-1)", "1/(x**2+2*x+5)", "(x+1)/(x**2+x+1)", "tan(x)", "log(x)/x", "1/(x*log(x))",
            "sqrt(1-x**2)", "asin(x)", "x**3*log(x)", "log(x)**2", "1/(2*x+3)**3", "x/(x+1)**2",
            "exp(x)/(1+exp(x))", "sec(x)**2", "sin(x)**3*cos(x)**2", "exp(-x**2)", "x*atan(x)",
        ] {
            check(f);
        }
    }
}
