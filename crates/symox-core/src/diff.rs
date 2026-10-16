//! Symbolic differentiation.

use rustc_hash::FxHashMap;

use crate::expr::{Expr, Func, Kind};
use crate::functions::*;
use crate::ops::{add, div, mul, mul_factor, neg, pow, sqrt, sub};

/// `d e / d x` where `x` is a symbol (or any atom-like expression such as
/// an undefined function application).
pub fn diff(e: &Expr, x: &Expr) -> Expr {
    let mut cache = FxHashMap::default();
    d(e, x, &mut cache)
}

/// n-th derivative.
pub fn diff_n(e: &Expr, x: &Expr, n: u32) -> Expr {
    let mut r = e.clone();
    for _ in 0..n {
        r = diff(&r, x);
        if r.is_zero() {
            break;
        }
    }
    r
}

/// Derivative with respect to several variables in sequence.
pub fn diff_many(e: &Expr, vars: &[(Expr, u32)]) -> Expr {
    let mut r = e.clone();
    for (v, n) in vars {
        r = diff_n(&r, v, *n);
    }
    r
}

/// Build an unevaluated `Derivative` node (merging nested derivatives).
pub fn derivative_unevaluated(f: &Expr, vars: Vec<(Expr, u32)>) -> Expr {
    let vars: Vec<(Expr, u32)> = vars.into_iter().filter(|(_, n)| *n > 0).collect();
    if vars.is_empty() {
        return f.clone();
    }
    let (inner, mut all) = match f.kind() {
        Kind::Derivative(inner, v) => (inner.clone(), v.clone()),
        _ => (f.clone(), Vec::new()),
    };
    for (v, n) in vars {
        match all.last_mut() {
            Some((w, m)) if *w==v => *m += n,
            _ => all.push((v, n)),
        }
    }
    Expr::from_kind(Kind::Derivative(inner, all))
}

fn d(e: &Expr, x: &Expr, cache: &mut FxHashMap<usize, (Expr,Expr)>) -> Expr {
    if e == x {
        return Expr::one();
    }
    if !e.may_contain(x) && x.is_symbol() {
        return Expr::zero();
    }
    match e.kind() {
        Kind::Num(_) | Kind::Const(_) | Kind::Sym(_) | Kind::Bool(_) => return Expr::zero(),
        _ => {}
    }
    if let Some((_,v)) = cache.get(&e.dag_id()) {
        return v.clone();
    }
    let r = match e.kind() {
        Kind::Add(_, terms) => {
            let v: Vec<Expr> = terms
                .iter()
                .map(|(t, k)| {
                    let dt = d(t, x, cache);
                    if dt.is_zero() {
                        dt
                    } else {
                        mul(&[Expr::num(k.clone()), dt])
                    }
                })
                .collect();
            add(&v)
        }
        Kind::Mul(c, facs) => {
            let factors: Vec<Expr> = facs.iter().map(|(b, x)| mul_factor(b, x)).collect();
            let mut terms = Vec::new();
            for i in 0..factors.len() {
                let di = d(&factors[i], x, cache);
                if di.is_zero() {
                    continue;
                }
                let mut v = Vec::with_capacity(factors.len() + 1);
                v.push(Expr::num(c.clone()));
                v.push(di);
                for (j, f) in factors.iter().enumerate() {
                    if j != i {
                        v.push(f.clone());
                    }
                }
                terms.push(mul(&v));
            }
            add(&terms)
        }
        Kind::Pow(b, n) => {
            let db = d(b, x, cache);
            let dn = d(n, x, cache);
            if dn.is_zero() {
                // n * b**(n-1) * b'
                if db.is_zero() {
                    Expr::zero()
                } else {
                    mul(&[n.clone(), pow(b, &sub(n, &Expr::one())), db])
                }
            } else if db.is_zero() {
                // b**n * log(b) * n'
                mul(&[e.clone(), log(b), dn])
            } else {
                // b**n * (n' log b + n b'/b)
                mul(&[e.clone(), add(&[mul(&[dn, log(b)]), mul(&[n.clone(), db, pow(b, &Expr::minus_one())])])])
            }
        }
        Kind::Func(f, args) => d_func(e, f, args, x, cache),
        Kind::Derivative(inner, vars) => {
            let di = d(inner, x, cache);
            match di.kind() {
                Kind::Derivative(..) => {
                    let mut v = vars.clone();
                    match v.iter_mut().find(|(w, _)| w == x) {
                        Some((_, m)) => *m += 1,
                        None => v.push((x.clone(), 1)),
                    }
                    Expr::from_kind(Kind::Derivative(inner.clone(), v))
                }
                _ => diff_many(&di, vars),
            }
        }
        Kind::Piecewise(p) => crate::relational::piecewise(
            p.iter().map(|(ex, c)| (d(ex, x, cache), c.clone())).collect(),
        ),
        _ => derivative_unevaluated(e, vec![(x.clone(), 1)]),
    };
    cache.insert(e.dag_id(), (e.clone(),r.clone()));
    r
}

fn d_func(
    e: &Expr,
    f: &Func,
    args: &[Expr],
    x: &Expr,
    cache: &mut FxHashMap<usize, (Expr,Expr)>,
) -> Expr {
    if matches!(f,Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::IdentityFunction)){
        return if crate::traverse::free_symbols(e).contains(x){derivative_unevaluated(e,vec![(x.clone(),1)])}else{Expr::zero()};
    }
    if let Func::MatrixElement(parent) = f {
        if let Kind::Func(Func::MatrixElement(other), indices) = x.kind() {
            let independent_symbols = matches!(parent.kind(), crate::matrixexpr::MatrixKind::Symbol(_))
                && matches!(other.kind(), crate::matrixexpr::MatrixKind::Symbol(_));
            if independent_symbols && (parent != other || args.iter().zip(indices).any(|(a, b)| {
                a.as_i64().zip(b.as_i64()).is_some_and(|(a, b)| a != b)
            })) {
                return Expr::zero();
            }
            return derivative_unevaluated(e, vec![(x.clone(), 1)]);
        }
        if !crate::traverse::has(e, x) {
            return Expr::zero();
        }
        if args.iter().all(|index| !crate::traverse::has(index, x)) {
            if let Ok(derivative) = parent.diff(x).and_then(|matrix| matrix.element(&args[0], &args[1])) {
                return derivative;
            }
        }
        return derivative_unevaluated(e, vec![(x.clone(), 1)]);
    }
    if let Func::MatrixTrace(parent) = f {
        if let Ok(evaluated) = parent.trace() {
            if evaluated != *e { return d(&evaluated, x, cache); }
        }
        return parent.diff(x).and_then(|matrix| matrix.trace())
            .unwrap_or_else(|_| derivative_unevaluated(e, vec![(x.clone(), 1)]));
    }
    if let Func::MatrixDeterminant(parent) = f {
        if let Ok(evaluated) = parent.determinant() {
            if evaluated != *e { return d(&evaluated, x, cache); }
        }
        return parent.determinant_diff(x)
            .unwrap_or_else(|_| derivative_unevaluated(e, vec![(x.clone(), 1)]));
    }
    if let Func::MatrixProxy(_) = f {
        return Expr::zero();
    }
    if matches!(f,Func::Integral|Func::Sum) {return crate::binders::differentiate(e,x);}
    if matches!(f,Func::Subs | Func::Limit){
        if !crate::traverse::free_symbols(e).contains(x){return Expr::zero();}
        let evaluated=crate::functions::doit(e,true);
        return if evaluated!=*e{d(&evaluated,x,cache)}else{derivative_unevaluated(e,vec![(x.clone(),1)])};
    }
    if let Func::Undefined(function) = f {
        if function.internal || args.iter().all(|arg|arg.is_symbol()) {
            return derivative_unevaluated(e, vec![(x.clone(), 1)]);
        }
        let mut terms=Vec::new();
        for (i,arg) in args.iter().enumerate(){
            let derivative=d(arg,x,cache);
            if derivative.is_zero(){continue;}
            if arg.is_symbol(){terms.push(mul(&[derivative,derivative_unevaluated(e,vec![(arg.clone(),1)])]));continue;}
            // One fresh atom shared by the argument, derivative and substitution
            // binder; user symbols with the same display name cannot be captured.
            let dummy=Expr::dummy(&format!("xi_{}",i+1),crate::Assumptions::NONE);
            let mut arguments=args.to_vec();arguments[i]=dummy.clone();
            let applied=crate::functions::func(f.clone(),arguments);
            let partial=derivative_unevaluated(&applied,vec![(dummy.clone(),1)]);
            let evaluated=crate::subs::held_subs(&partial,vec![dummy],vec![arg.clone()]);
            terms.push(mul(&[derivative,evaluated]));
        }
        return add(&terms);
    }
    if let Func::Atan2 = f {
        let (y, xx) = (&args[0], &args[1]);
        let (dy, dx) = (d(y, x, cache), d(xx, x, cache));
        return div(
            &sub(&mul(&[xx.clone(), dy]), &mul(&[y.clone(), dx])),
            &add(&[xx.powi(2), y.powi(2)]),
        );
    }
    if matches!(f, Func::Max | Func::Min) {
        return derivative_unevaluated(e, vec![(x.clone(), 1)]);
    }
    if let Func::Special(s) = f {
        // chain rule over all arguments
        let mut terms = Vec::new();
        for (i, a) in args.iter().enumerate() {
            if let Some(items) = crate::special::tuple_items(a) {
                // parameter tuples (hyper): only constant ones are supported
                if items.iter().all(|t| d(t, x, cache).is_zero()) {
                    continue;
                }
                return derivative_unevaluated(e, vec![(x.clone(), 1)]);
            }
            let da = d(a, x, cache);
            if da.is_zero() {
                continue;
            }
            match crate::special::diff::partial(*s, args, i) {
                Some(p) => terms.push(mul(&[p, da])),
                None => return derivative_unevaluated(e, vec![(x.clone(), 1)]),
            }
        }
        return add(&terms);
    }
    let u = &args[0];
    let du = d(u, x, cache);
    if du.is_zero() {
        return Expr::zero();
    }
    if matches!(f, Func::Conjugate | Func::Re | Func::Im)
        && crate::assumptions::is_real(x) != Some(true)
    {
        return derivative_unevaluated(e, vec![(x.clone(), 1)]);
    }
    let one = Expr::one();
    let outer = match f {
        Func::Exp => e.clone(),
        Func::Log => pow(u, &Expr::minus_one()),
        Func::Sin => cos(u),
        Func::Cos => neg(&sin(u)),
        Func::Tan => add(&[tan(u).powi(2), one]),
        Func::Cot => neg(&add(&[cot(u).powi(2), one])),
        Func::Sec => mul(&[tan(u), sec(u)]),
        Func::Csc => neg(&mul(&[cot(u), csc(u)])),
        Func::Asin => pow(&sub(&one, &u.powi(2)), &Expr::rational(-1, 2)),
        Func::Acos => neg(&pow(&sub(&one, &u.powi(2)), &Expr::rational(-1, 2))),
        Func::Atan => pow(&add(&[u.powi(2), one]), &Expr::minus_one()),
        Func::Acot => neg(&pow(&add(&[u.powi(2), one]), &Expr::minus_one())),
        Func::Sinh => cosh(u),
        Func::Cosh => sinh(u),
        Func::Tanh => sub(&one, &tanh(u).powi(2)),
        Func::Coth => neg(&pow(&sinh(u), &Expr::int(-2))),
        Func::Asinh => pow(&add(&[u.powi(2), one]), &Expr::rational(-1, 2)),
        Func::Acosh => pow(
            &mul(&[sqrt(&sub(u, &one)), sqrt(&add(&[u.clone(), one]))]),
            &Expr::minus_one(),
        ),
        Func::Atanh => pow(&sub(&one, &u.powi(2)), &Expr::minus_one()),
        Func::Abs => sign(u),
        Func::Sign | Func::Floor | Func::Ceiling => Expr::zero(),
        Func::Erf => mul(&[
            Expr::int(2),
            pow(&Expr::pi(), &Expr::rational(-1, 2)),
            exp(&neg(&u.powi(2))),
        ]),
        Func::Erfc => mul(&[
            Expr::int(-2),
            pow(&Expr::pi(), &Expr::rational(-1, 2)),
            exp(&neg(&u.powi(2))),
        ]),
        Func::Heaviside => dirac_delta(u),
        Func::Gamma => mul(&[e.clone(), crate::special::polygamma(&Expr::zero(), u)]),
        Func::LogGamma => crate::special::polygamma(&Expr::zero(), u),
        Func::Re => return re(&du),
        Func::Im => return im(&du),
        Func::Conjugate => return conjugate(&du),
        _ => return derivative_unevaluated(e, vec![(x.clone(), 1)]),
    };
    mul(&[outer, du])
}
