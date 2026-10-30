//! Symbolic differentiation.

pub mod symbolic;
mod integer_orders;
pub mod order_plan;

use rustc_hash::FxHashMap;

use crate::expr::{Expr, Func, Kind};
use crate::functions::*;
use crate::ops::{add, div, mul, mul2, mul_factor, neg, pow, sqrt, sub};

pub fn is_valid_variable(variable: &Expr) -> bool {
    match variable.kind() {
        Kind::Sym(_) | Kind::Func(Func::MatrixElement(_), _) => true,
        Kind::Func(Func::SymbolicDerivative, _) => symbolic::valid_variable(variable).unwrap_or_else(|error|{symbolic::record(error);false}),
        Kind::Derivative(inner, variables) => {
            let Kind::Func(Func::Undefined(head), arguments) = inner.kind() else { return false; };
            !head.internal && variables.iter().all(|(variable, _)| {
                arguments.iter().filter(|argument| *argument == variable).count() == 1 &&
                arguments.iter().all(|argument| argument == variable || !crate::traverse::has(argument, variable))
            })
        }
        Kind::Func(Func::Undefined(head), _) => !head.internal || head.is_internal("Indexed"),
        _ => false,
    }
}


/// `d e / d x` where `x` is a symbol (or any atom-like expression such as
/// an undefined function application).
pub fn diff(e: &Expr, x: &Expr) -> Expr {
    let mut cache = FxHashMap::default();
    d(e, x, &mut cache)
}

/// n-th derivative.
pub fn diff_n(e: &Expr, x: &Expr, n: u32) -> Expr {
    if n>1 {if let Some(value)=integer_orders::evaluate(e,x,&num_bigint::BigInt::from(n)){return value;}}
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

/// An `AppliedUndef` of the reference: an application of `Function('f')`.
fn applied_undef(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Undefined(head), _) if !head.internal && head.host().is_none())
}

/// The reference's `_nodes`, the first key of `ordered`.
fn sympy_nodes(e: &Expr) -> usize {
    match e.kind() {
        Kind::Derivative(inner, variables) =>
            variables.iter().fold(sympy_nodes(inner), |n, (_, count)| n.saturating_add(*count as usize)),
        _ => crate::traverse::node_count(e),
    }
}

/// `_block` of the reference `_sort_variable_count`: true if `v` may not
/// come before `d`.
fn sort_block(d: &Expr, v: &Expr, wrt: bool) -> bool {
    if d == v {
        return wrt;
    }
    if d.is_symbol() {
        return false;
    }
    if let Some((_,variables))=symbolic::parts(d) {
        return variables.iter().any(|(k,_)|sort_block(k,v,true));
    }
    if let Kind::Derivative(_, variables) = d.kind() {
        return variables.iter().any(|(k, _)| sort_block(k, v, true));
    }
    if !wrt && applied_undef(d) {
        return false;
    }
    if v.is_symbol() {
        return crate::traverse::free_symbols(d).contains(v);
    }
    if applied_undef(v) {
        // _block(d.xreplace({v: D}), D)
        return crate::traverse::has_structural(d, v);
    }
    let free = crate::traverse::free_symbols(v);
    crate::traverse::free_symbols(d).iter().any(|s| free.contains(s))
}

/// SymPy's `Derivative._sort_variable_count`: the canonical order of the
/// differentiation variables, keeping the order of the variables that do
/// not commute during differentiation, with contiguous counts merged.
pub fn sort_variable_count(vc: &[(Expr, u32)]) -> Vec<(Expr, u32)> {
    if vc.len() <= 1 {
        return vc.to_vec();
    }
    let order = variable_order(vc);
    let mut merged: Vec<(Expr, u32)> = Vec::with_capacity(vc.len());
    for i in order {
        let (v, c) = &vc[i];
        if let Some((w, m)) = merged.last_mut() {
            if w == v {
                if let Some(total) = m.checked_add(*c) {
                    *m = total;
                    continue;
                }
            }
        }
        merged.push((v.clone(), *c));
    }
    merged
}

pub(super) fn variable_order<T>(vc: &[(Expr,T)]) -> Vec<usize> {
    let n = vc.len();
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        for j in 0..i {
            if sort_block(&vc[j].0, &vc[i].0, false) {
                edges.push((j, i));
            }
        }
    }
    // O = dict(zip(ordered(uniq(variables)), range(len(vc))))
    let mut unique: Vec<Expr> = Vec::new();
    for (v,_) in vc {
        if !unique.contains(v) {
            unique.push(v.clone());
        }
    }
    let mut keyed: Vec<(crate::number::Number, crate::sortkey::Key, Expr)> =
        unique.into_iter().map(|v| {
            let weight=if symbolic::parts(&v).is_some(){symbolic::node_weight(&v)}else{crate::number::Number::from_bigint(sympy_nodes(&v).into())};
            (weight,crate::sortkey::sort_key(&v),v)
        }).collect();
    keyed.sort_by(|a,b|a.0.cmp_value(&b.0).unwrap_or_else(||a.0.cmp_struct(&b.0)).then_with(||a.1.cmp(&b.1)));
    let rank = |i: usize| keyed.iter().position(|(_, _, v)| *v == vc[i].0).unwrap_or(0);
    // topological_sort((V, E), key=lambda i: O[v(i)])
    let mut start: Vec<usize> = (0..n).filter(|i| !edges.iter().any(|(_, u)| u == i)).collect();
    // sorted(S, key=key, reverse=True). Equal keys are the same variable,
    // whose entries are merged below; they are taken in input order so that
    // counts that cannot be merged (u32 overflow) keep their order.
    start.sort_by(|a, b| rank(*b).cmp(&rank(*a)).then_with(|| b.cmp(a)));
    let mut pending = start;
    let mut order = Vec::with_capacity(n);
    while let Some(node) = pending.pop() {
        order.push(node);
        let targets: Vec<usize> = edges.iter().filter(|(u, _)| *u == node).map(|(_, v)| *v).collect();
        for v in targets {
            if let Some(position) = edges.iter().position(|(u, w)| *u == node && *w == v) {
                edges.remove(position);
            }
            if edges.iter().any(|(_, w)| *w == v) {
                continue;
            }
            let kv = rank(v);
            match pending.iter().position(|s| kv > rank(*s)) {
                Some(i) => pending.insert(i, v),
                None => pending.push(v),
            }
        }
    }
    if order.len() != n {
        return (0..vc.len()).collect();
    }
    order
}

/// `Derivative.canonical` applied to every derivative of `e`, as the
/// reference's evaluated `Derivative` ends with
/// `expr.replace(lambda x: isinstance(x, Derivative), lambda x: x.canonical)`.
pub fn canonical_derivatives(e: &Expr) -> Expr {
    fn walk(e: &Expr, memo: &mut FxHashMap<usize, (Expr, Expr)>) -> Expr {
        if e.is_atom() {
            return e.clone();
        }
        if let Some((_, value)) = memo.get(&e.dag_id()) {
            return value.clone();
        }
        let value = match e.kind() {
            Kind::Num(_) | Kind::Const(_) | Kind::Sym(_) | Kind::Bool(_) => e.clone(),
            Kind::Func(Func::SymbolicDerivative, _) => symbolic::canonical(e),
            Kind::Derivative(inner, variables) => {
                let body = walk(inner, memo);
                let sorted = if variables.len() > 1 { sort_variable_count(variables) } else { variables.clone() };
                if body.same_presentation(inner) && sorted == *variables {
                    e.clone()
                } else {
                    derivative_unevaluated(&body, sorted)
                }
            }
            _ => crate::traverse::map_children(e, &mut |child| walk(child, memo)),
        };
        memo.insert(e.dag_id(), (e.clone(), value.clone()));
        value
    }
    walk(e, &mut FxHashMap::default())
}

/// The reference's `Derivative(e, *vars, evaluate=True)`: the variables in
/// canonical order, then every derivative of the result made canonical.
pub fn diff_evaluated(e: &Expr, vars: &[(Expr, u32)]) -> Expr {
    let sorted = sort_variable_count(vars);
    canonical_derivatives(&diff_many(e, &sorted))
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
        if let Some((w, m)) = all.last_mut() {
            if *w == v {
                if let Some(total) = m.checked_add(n) {
                    *m = total;
                    continue;
                }
            }
        }
        // Keep adjacent groups if their total exceeds u32; the held
        // derivative still represents the exact sum of their orders.
        all.push((v, n));
    }
    Expr::from_kind(Kind::Derivative(inner, all))
}

fn pure_imaginary(value:&Expr)->bool{
    if value.as_sym().is_some_and(|s|s.assumptions().has(crate::Assumptions::IMAGINARY)){return true;}
    crate::assumptions::exact_queries(||re(value).is_zero()&&crate::assumptions::is_zero(value)==Some(false))
}

fn d(e: &Expr, x: &Expr, cache: &mut FxHashMap<usize, (Expr,Expr)>) -> Expr {
    if e == x {
        return Expr::one();
    }
    if let Some(v) = crate::series::order_diff::order_derivative(e, x) {
        return v;
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
                // SymPy's Pow._eval_derivative, self*(b'*n/b), built with
                // the same binary products: 1/b splits a product base, so
                // sqrt(-x)' is sqrt(-x)/(2*x).
                if db.is_zero() {
                    Expr::zero()
                } else {
                    mul2(e, &mul2(&mul2(&db, n), &pow(b, &Expr::minus_one())))
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
                Kind::Derivative(body, more) => {
                    // The reference re-dispatches Derivative(expr,
                    // *(variable_count + [(x, 1)]), evaluate=True), whose
                    // variables are in canonical order.
                    let combined: Vec<(Expr, u32)> = if body == inner {
                        vars.iter().chain(more.iter()).cloned().collect()
                    } else {
                        more.iter().chain(vars.iter()).cloned().collect()
                    };
                    derivative_unevaluated(body, sort_variable_count(&combined))
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

/// SymPy's `ArrayElement._eval_derivative`: an element differentiates only
/// by an element of the same array, to the product of the Kronecker deltas
/// of the index pairs; by anything else (a symbol, an index) it is 0.
fn array_element_derivative(e: &Expr, head: &crate::expr::UndefinedFunction, args: &[Expr], x: &Expr) -> Expr {
    let Kind::Func(Func::Undefined(other), indices) = x.kind() else { return Expr::zero() };
    if other != head || indices.is_empty() || args.is_empty() { return Expr::zero(); }
    if x == e { return Expr::one(); }
    if indices[0] != args[0] || indices.len() != args.len() { return Expr::zero(); }
    let deltas: Vec<Expr> = args[1..].iter().zip(&indices[1..]).map(|(i, j)| {
        crate::functions::func(Func::Special(crate::special::SpecialFn::KroneckerDelta), vec![i.clone(), j.clone()])
    }).collect();
    mul(&deltas)
}

fn d_func(
    e: &Expr,
    f: &Func,
    args: &[Expr],
    x: &Expr,
    cache: &mut FxHashMap<usize, (Expr,Expr)>,
) -> Expr {
    if matches!(f,Func::SymbolicDerivative) {return symbolic::differentiate(e,x);}

    if matches!(f, Func::GenericExpr) {
        return derivative_unevaluated(e, vec![(x.clone(), 1)]);
    }
    if matches!(f,Func::Set(crate::sets::SetOp::Lambda|crate::sets::SetOp::IdentityFunction)){
        return if crate::traverse::free_symbols(e).contains(x){derivative_unevaluated(e,vec![(x.clone(),1)])}else{Expr::zero()};
    }
    if let Func::MatrixElement(parent) = f {
        if let Kind::Func(Func::MatrixElement(other), indices) = x.kind() {
            // MatrixElement._eval_derivative by an element of the same matrix.
            if parent == other {
                return mul2(&crate::matrixexpr::ranged_delta(&args[0], &indices[0], parent.rows()),
                            &crate::matrixexpr::ranged_delta(&args[1], &indices[1], parent.cols()));
            }
            // Inverse: -Sum(M[i, z1]*Y[z1, z2].diff(v)*M[z2, j], (z1, 0, r1 - 1), (z2, 0, r2 - 1)).
            if let crate::matrixexpr::MatrixKind::Inverse(y) = parent.kind() {
                let result = (|| -> crate::matrix::MResult<Expr> {
                    let z1 = Expr::dummy("z1", crate::symbol::Assumptions::NONE);
                    let z2 = Expr::dummy("z2", crate::symbol::Assumptions::NONE);
                    let inner = d(&y.element(&z1, &z2)?, x, cache);
                    let body = mul(&[parent.element(&args[0], &z1)?, inner, parent.element(&z2, &args[1])?]);
                    let one = Expr::one();
                    let limits = vec![vec![z1.clone(), Expr::zero(), sub(y.rows(), &one)],
                                      vec![z2.clone(), Expr::zero(), sub(y.cols(), &one)]];
                    let sum = crate::binders::held(Func::Sum, &body, limits).map_err(crate::matrix::MatrixError::Value)?;
                    Ok(neg(&sum))
                })();
                if let Ok(value) = result { return value; }
            }
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
    if let Func::MatrixDotProduct(a,b)=f {
        let result=(||->crate::matrix::MResult<Expr>{
            let da=a.diff(x)?;let db=b.diff(x)?;
            let zero=|m:&crate::matrixexpr::MatrixExpr|match m.kind(){crate::matrixexpr::MatrixKind::Zero=>true,crate::matrixexpr::MatrixKind::Explicit(m)=>m.data().iter().all(Expr::is_zero),_=>false};
            let mut terms=Vec::new();
            if !zero(&da){terms.push(crate::matrixcontraction::dot_product(&da,b,false)?);}
            if !zero(&db){terms.push(crate::matrixcontraction::dot_product(a,&db,false)?);}
            Ok(add(&terms))
        })();
        return result.unwrap_or_else(|_|derivative_unevaluated(e,vec![(x.clone(),1)]));
    }
    if let Func::MatrixTrace(parent) = f {
        if matches!(x.kind(), Kind::Func(Func::MatrixElement(_), _)) {
            // Trace._eval_derivative by a matrix element: rewrite(Sum).diff(v).
            let index = Expr::dummy("i", crate::symbol::Assumptions::NONE);
            if let Ok(entry) = parent.element(&index, &index) {
                if let Ok(sum) = crate::binders::held(Func::Sum, &entry, vec![vec![index.clone(), Expr::zero(), sub(parent.rows(), &Expr::one())]]) {
                    return d(&sum, x, cache);
                }
            }
            return derivative_unevaluated(e, vec![(x.clone(), 1)]);
        }
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
    if matches!(f,Func::Integral|Func::Sum|Func::Product) {return crate::binders::differentiate(e,x);}
    if matches!(f,Func::Subs) && args.len()==3 {
        if !x.is_symbol() {
            // Derivative.__new__: a non-symbol variable is replaced by a
            // dummy, the derivative taken by it and the dummy substituted back.
            if !crate::traverse::has_structural(e,x){return Expr::zero();}
            let dummy=Expr::dummy("xi",crate::Assumptions::NONE);
            let mut m=FxHashMap::default();m.insert(x.clone(),dummy.clone());
            let masked=crate::subs::xreplace(e,&m);
            if masked==*e{return derivative_unevaluated(e,vec![(x.clone(),1)]);}
            let value=d(&masked,&dummy,&mut FxHashMap::default());
            return crate::subs::subs1(&value,&dummy,x);
        }
        if !crate::traverse::free_symbols(e).contains(x){return Expr::zero();}
        // Subs._eval_derivative: the chain rule over the points, plus the
        // derivative of the expression when it depends on x otherwise.
        let (body,variables,points)=(&args[0],args[1].args(),args[2].args());
        let mut terms=Vec::new();
        for (v,p) in variables.iter().zip(points.iter()){
            let dp=d(p,x,cache);
            if dp.is_zero(){continue;}
            let inner=crate::functions::doit(&crate::subs::held_subs(&diff(body,v),variables.clone(),points.clone()),true);
            terms.push(mul(&[dp,inner]));
        }
        let free=crate::traverse::free_symbols(body);
        if free.contains(x)&&!variables.contains(x){
            terms.push(crate::functions::doit(&crate::subs::held_subs(&diff(body,x),variables.clone(),points.clone()),true));
        }
        return add(&terms);
    }
    if matches!(f,Func::Subs | Func::Limit){
        if !crate::traverse::free_symbols(e).contains(x){return Expr::zero();}
        let evaluated=crate::functions::doit(e,true);
        return if evaluated!=*e{d(&evaluated,x,cache)}else{derivative_unevaluated(e,vec![(x.clone(),1)])};
    }
    if let Func::Undefined(function) = f {
        if function.host().is_some_and(|host| host.array_element()) {
            return array_element_derivative(e, function, args, x);
        }
        if let Some(host) = function.host().filter(|host| host.expression()) {
            // An Expr subclass without a chain rule: its _eval_derivative,
            // otherwise the derivative is held (Expr._eval_derivative_n_times).
            if !crate::traverse::has(e, x) { return Expr::zero(); }
            return host.derivative(e, x).unwrap_or_else(|| derivative_unevaluated(e, vec![(x.clone(), 1)]));
        }
        if let Some(host) = function.host().filter(|host| host.differentiates()) {
            // An overridden Function._eval_derivative replaces the chain rule.
            if !crate::traverse::has(e, x) { return Expr::zero(); }
            return host.derivative(e, x).unwrap_or_else(|| derivative_unevaluated(e, vec![(x.clone(), 1)]));
        }
        if let Some(host) = function.host() {
            // SymPy's Function._eval_derivative: the chain rule over fdiff.
            // An argument without a host partial keeps the default form.
            let mut terms=Vec::new();
            let mut defined=true;
            let mut partials=Vec::with_capacity(args.len());
            for (i,arg) in args.iter().enumerate(){
                let derivative=d(arg,x,cache);
                if derivative.is_zero(){partials.push(None);continue;}
                match host.fdiff(e,i+1){
                    Some(partial)=>partials.push(Some((derivative,partial))),
                    None=>{defined=false;break;}
                }
            }
            if defined {
                for (derivative,partial) in partials.into_iter().flatten(){terms.push(mul(&[partial,derivative]));}
                return add(&terms);
            }
        }
        // A non-symbol variable (g(t), Derivative(g(t), t)) goes through the
        // chain rule below, where symbol arguments contribute zero, as in SymPy.
        if function.internal || (x.is_symbol() && args.iter().all(|arg|arg.is_symbol())) {
            // The bloom mask is only a negative shortcut. Dummy identities
            // can collide with an unrelated differentiation symbol's bit.
            if x.is_symbol() && !crate::traverse::free_symbols(e).contains(x) {
                return Expr::zero();
            }
            return derivative_unevaluated(e, vec![(x.clone(), 1)]);
        }
        let mut terms=Vec::new();
        for (i,arg) in args.iter().enumerate(){
            let derivative=d(arg,x,cache);
            if derivative.is_zero(){continue;}
            // AppliedUndef.fdiff: a symbol argument is held directly only
            // when no other argument depends on it (SymPy issue 8510).
            if arg.is_symbol()&&args.iter().enumerate().all(|(j,other)|j==i||!crate::traverse::free_symbols(other).contains(arg)){terms.push(mul(&[derivative,derivative_unevaluated(e,vec![(arg.clone(),1)])]));continue;}
            // A non-symbol argument the reference can differentiate by
            // (g(t), Derivative(g(t), t)) is held as Derivative(f(g(t)), g(t)).
            // The reference also does this when another argument contains it,
            // f(g(t), g(t)**2), where that derivative by g(t) also acts on the
            // other argument and double counts; those keep the Subs form.
            if !arg.is_symbol()&&is_valid_variable(arg)&&!matches!(arg.kind(),Kind::Func(Func::MatrixElement(_),_))
                &&args.iter().enumerate().all(|(j,other)|j==i||!crate::traverse::has_structural(other,arg)){
                terms.push(mul(&[derivative,derivative_unevaluated(e,vec![(arg.clone(),1)])]));continue;
            }
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
        // SymPy's `MinMaxBase._eval_derivative` with `Max.fdiff` and
        // `Min.fdiff`: each argument contributes a Heaviside factor.
        let is_max = matches!(f, Func::Max);
        let n = args.len();
        let mut terms = Vec::new();
        for (i, a) in args.iter().enumerate() {
            let da = d(a, x, cache);
            if da.is_zero() {
                continue;
            }
            let others = if n == 2 {
                args[1 - i].clone()
            } else {
                let rest: Vec<Expr> = args.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, b)| b.clone()).collect();
                if is_max { max(rest) } else { min(rest) }
            };
            let gap = if is_max { sub(a, &others) } else { sub(&others, a) };
            terms.push(mul(&[heaviside(&gap), da]));
        }
        return add(&terms);
    }
    if let Func::Special(s) = f {
        if s.is_transform() {
            // held integral transforms: the bound variable is not free
            if !crate::traverse::free_symbols(e).contains(x) {
                return Expr::zero();
            }
            return derivative_unevaluated(e, vec![(x.clone(), 1)]);
        }
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
    if matches!(f, Func::Heaviside) && args.len() == 2 {
        let dh0 = d(&args[1], x, cache);
        if !dh0.is_zero() {
            // Heaviside.fdiff(2) raises ArgumentIndexError, so SymPy falls
            // back to Function.fdiff: a symbol H0 absent from the first
            // argument gives Derivative(Heaviside(u, H0), H0), anything else
            // Subs(Derivative(Heaviside(u, xi_2), xi_2), xi_2, H0).
            let h0 = &args[1];
            let partial = if h0.is_symbol() && !crate::traverse::free_symbols(u).contains(h0) {
                derivative_unevaluated(e, vec![(h0.clone(), 1)])
            } else {
                let dummy = Expr::dummy("xi_2", crate::Assumptions::NONE);
                let applied = crate::functions::func(Func::Heaviside, vec![u.clone(), dummy.clone()]);
                let inner = derivative_unevaluated(&applied, vec![(dummy.clone(), 1)]);
                crate::subs::held_subs(&inner, vec![dummy], vec![h0.clone()])
            };
            let held = mul(&[partial, dh0]);
            if du.is_zero() {
                return held;
            }
            return add(&[mul(&[dirac_delta(u), du]), held]);
        }
    }
    if du.is_zero() {
        return Expr::zero();
    }
    if matches!(f, Func::DiracDelta) {
        let k = args.get(1).cloned().unwrap_or_else(Expr::zero);
        if args.len() == 2 && d(&k, x, cache).is_zero() || args.len() == 1 {
            let order = add(&[k, Expr::one()]);
            return mul(&[crate::functions::func(Func::DiracDelta, vec![u.clone(), order]), du]);
        }
    }
    if matches!(f, Func::Conjugate | Func::Re | Func::Im){
        let real_x=crate::assumptions::exact_queries(||crate::assumptions::is_real(x)==Some(true));
        let imaginary_x=pure_imaginary(x);
        if matches!(f,Func::Conjugate){
            return if real_x{conjugate(&du)}else if imaginary_x{neg(&conjugate(&du))}
                else{derivative_unevaluated(e,vec![(x.clone(),1)])};
        }
        if real_x||crate::assumptions::exact_queries(||crate::assumptions::is_real(u)==Some(true)){
            return if matches!(f,Func::Re){re(&du)}else{im(&du)};
        }
        if imaginary_x||pure_imaginary(u){
            let part=if matches!(f,Func::Re){im(&du)}else{re(&du)};
            return mul(&[if matches!(f,Func::Re){Expr::i()}else{neg(&Expr::i())},part]);
        }
        return derivative_unevaluated(e,vec![(x.clone(),1)]);
    }
    if matches!(f,Func::Abs){
        if crate::assumptions::exact_queries(||crate::assumptions::is_real(u)==Some(true))||pure_imaginary(u){
            return mul(&[du,sign(&conjugate(u))]);
        }
        let real=re(u);let imaginary=im(u);
        let numerator=add(&[mul(&[real.clone(),d(&real,x,cache)]),mul(&[imaginary.clone(),d(&imaginary,x,cache)])]);
        // Abs(u).rewrite(sign) = u/sign(u), including held partial
        // derivatives when differentiation is with respect to a complex atom.
        return mul(&[numerator,sign(u),pow(u,&Expr::minus_one())]);
    }
    if matches!(f, Func::Arg) {
        // SymPy's arg._eval_derivative over the real and imaginary parts.
        if let Ok(Some((real, imaginary))) = crate::complex_parts::as_real_imag(u, true, None) {
            let numerator = sub(
                &mul(&[real.clone(), d(&imaginary, x, cache)]),
                &mul(&[imaginary.clone(), d(&real, x, cache)]),
            );
            return mul(&[numerator, pow(&add(&[real.powi(2), imaginary.powi(2)]), &Expr::minus_one())]);
        }
    }
    let one = Expr::one();
    let outer = match f {
        Func::Exp | Func::ExpPolar => e.clone(),
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

#[cfg(test)]
mod dependency_tests {
    use super::*;

    #[test]
    fn held_derivative_order_groups_do_not_overflow() {
        let x = Expr::symbol("x");
        let f = crate::functions::undefined("f", vec![x.clone()]);
        let first = derivative_unevaluated(&f, vec![(x.clone(), u32::MAX)]);
        let second = derivative_unevaluated(&first, vec![(x.clone(), 1)]);
        let Kind::Derivative(inner, orders) = second.kind() else { panic!("held derivative expected"); };
        assert_eq!(*inner, f);
        assert_eq!(orders, &vec![(x.clone(), u32::MAX), (x.clone(), 1)]);
        assert_eq!(diff(&first, &x), second);
        let ordinary = derivative_unevaluated(&derivative_unevaluated(&f, vec![(x.clone(), 2)]),vec![(x.clone(), 3)]);
        assert_eq!(ordinary, derivative_unevaluated(&f, vec![(x, 5)]));
    }

    #[test]
    fn dummy_bloom_collision_cannot_create_a_false_partial_derivative() {
        let y=Expr::symbol("y");
        let dummy=(0..4096).map(|index|Expr::dummy_with_index("integration",crate::Assumptions::NONE,num_bigint::BigInt::from(index)))
            .find(|candidate|candidate.may_contain(&y)).expect("a collision in the 64-bit dependency filter");
        assert_ne!(dummy,y);
        let f=crate::functions::undefined("bloom_dependency",vec![dummy.clone()]);
        assert!(f.may_contain(&y));
        assert_eq!(diff(&f,&y),Expr::zero());
        assert_eq!(diff(&f,&dummy),derivative_unevaluated(&f,vec![(dummy.clone(),1)]));
        let integral=crate::binders::held(Func::Integral,&mul(&[y.clone(),f.clone()]),vec![vec![dummy.clone(),Expr::zero(),y.clone()]]).unwrap();
        let endpoint=crate::functions::undefined("bloom_dependency",vec![y.clone()]);
        let remaining=crate::binders::held(Func::Integral,&f,vec![vec![dummy,Expr::zero(),y.clone()]]).unwrap();
        assert_eq!(diff(&integral,&y),add(&[mul(&[y,endpoint]),remaining]));
    }
}
