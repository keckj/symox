//! Substitution: `subs` (algebraic, sequential) and `xreplace` (exact).

use num_integer::Integer;
use num_traits::ToPrimitive;
use rustc_hash::FxHashMap;

use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, add_from_terms, div, mul, mul2, mul_factor, pow};
use crate::traverse::map_children;
mod memo;

/// Simultaneous exact structural replacement.
pub fn xreplace(e: &Expr, map: &FxHashMap<Expr, Expr>) -> Expr {
    if map.is_empty() {
        return e.clone();
    }
    let mask: u64 = map.keys().fold(0, |m, k| m | k.mask());
    let all_masked = map.keys().all(|k| k.mask() != 0 && !k.has_binding_identity());
    let mut cache: FxHashMap<usize, (Expr,Expr)> = FxHashMap::default();
    xr(e, map, mask, all_masked, false, &mut cache)
}

fn xr(
    e: &Expr,
    map: &FxHashMap<Expr, Expr>,
    mask: u64,
    all_masked: bool,
    held_derivatives:bool,
    cache: &mut FxHashMap<usize, (Expr,Expr)>,
) -> Expr {
    if let Some(v) = map.get(e) {
        return v.clone();
    }
    if e.is_atom() || (all_masked && e.mask() & mask == 0) {
        return e.clone();
    }
    if let Some((_,v)) = cache.get(&e.dag_id()) {
        return v.clone();
    }
    if held_derivatives {
        if matches!(e.kind(),Kind::Func(Func::Integral|Func::Sum|Func::Product,_)) {
            let pairs=map.iter().map(|(key,value)|(key.clone(),value.clone())).collect::<Vec<_>>();
            return crate::binders::substitute(e,&pairs,true);
        }
        if let Kind::Func(Func::Subs,args)=e.kind(){
            // Subs binds variables in its body.  Evaluation points are ordered:
            // a bound variable may occur freely in its own/later points, while
            // earlier points are protected by the later substitution.
            let variables=args[1].args();
            let body_map:FxHashMap<_,_>=map.iter().filter(|(key,_)|!variables.contains(key)).map(|(k,v)|(k.clone(),v.clone())).collect();
            let body=xr(&args[0],&body_map,mask,all_masked,true,&mut FxHashMap::default());
            let points=args[2].args().iter().enumerate().map(|(i,point)|{
                let point_map:FxHashMap<_,_>=map.iter().filter(|(key,_)|!variables[i+1..].contains(key)).map(|(k,v)|(k.clone(),v.clone())).collect();
                xr(point,&point_map,mask,all_masked,true,&mut FxHashMap::default())
            }).collect();
            let result=held_subs(&body,variables,points);
            cache.insert(e.dag_id(),(e.clone(),result.clone()));
            return result;
        }
        if let Kind::Derivative(_,vars)=e.kind(){
            if vars.iter().any(|(v,_)|map.contains_key(v)){
                let mut keys:Vec<_>=map.keys().filter(|key|crate::traverse::has(e,key)).cloned().collect();
                keys.sort_by_cached_key(crate::sortkey::sort_key);
                let points=keys.iter().map(|key|map[key].clone()).collect();
                return held_subs(e,keys,points);
            }
        }
    }
    let r = if e.is_add() || e.is_mul() {
        // Replace public Add terms and Mul powers, including their compact
        // numeric coefficients. Internal map_children skips those nodes.
        let children = crate::inspection::scalar_children(e);
        let replaced = children.iter().map(|c|xr(c,map,mask,all_masked,held_derivatives,cache)).collect::<Vec<_>>();
        if replaced.iter().zip(&children).all(|(a,b)|a.same_presentation(b)) {
            e.clone()
        } else if e.is_add() { add(&replaced) } else { mul(&replaced) }
    } else {
        map_children(e, &mut |c| xr(c, map, mask, all_masked, held_derivatives, cache))
    };
    cache.insert(e.dag_id(), (e.clone(),r.clone()));
    r
}

/// SymPy's `subs`: apply each `(old, new)` pair in turn.
pub fn subs(e: &Expr, pairs: &[(Expr, Expr)]) -> Expr {
    let mut r = e.clone();
    for (old, new) in pairs {
        r = subs1(&r, old, new);
    }
    r
}

pub fn held_subs(e:&Expr, variables:Vec<Expr>,points:Vec<Expr>)->Expr {
    if variables.is_empty(){return e.clone();}
    Expr::from_kind(Kind::Func(Func::Subs,vec![e.clone(),crate::functions::undefined("Tuple",variables),crate::functions::undefined("Tuple",points)]))
}

/// Simultaneous substitution for symbol keys, preserving differentiation
/// variables by routing them through the same held-Subs rules as subs1.
pub fn simultaneous(e:&Expr,pairs:&[(Expr,Expr)])->Expr {
    if pairs.is_empty(){return e.clone();}
    if e.has_host_subs() || pairs.iter().any(|(old,_)|!old.is_symbol()) {
        // Mask algebraic patterns before inserting any replacement value.
        // A common marker protects zero/pole interactions while masking.
        let marker=Expr::dummy("subs_m",Default::default());
        let mut masked=e.clone();
        let mut restore=FxHashMap::default();
        for (old,new) in pairs {
            let temporary=Expr::dummy("subs_d",Default::default());
            masked=subs1(&masked,old,&mul2(&temporary,&marker));
            restore.insert(temporary,new.clone());
        }
        restore.insert(marker,Expr::one());
        return xreplace(&masked,&restore);
    }
    let map:FxHashMap<_,_>=pairs.iter().cloned().collect();
    let mask=map.keys().fold(0,|m,k|m|k.mask());
    let all_masked=map.keys().all(|k|k.mask()!=0 && !k.has_binding_identity());
    xr(e,&map,mask,all_masked,true,&mut FxHashMap::default())
}

pub fn subs_doit(e:&Expr,variables:&[Expr],points:&[Expr])->Expr {
    let active:Vec<_>=variables.iter().cloned().zip(points.iter().cloned()).filter(|(v,p)|v!=p).collect();
    if active.is_empty(){return e.clone();}
    fn bound_derivative(e:&Expr,variables:&[Expr])->bool {
        if let Kind::Derivative(_,vars)=e.kind(){if vars.iter().any(|(v,_)|variables.contains(v)){return true;}}
        e.args().iter().any(|a|bound_derivative(a,variables))
    }
    if bound_derivative(e,variables){return held_subs(e,variables.to_vec(),points.to_vec());}
    // Like SymPy Subs.doit, evaluation follows the specified variable order;
    // this is distinct from expr.subs(..., simultaneous=True).
    subs(e,&active)
}

/// Host solver of a moving Order point: `(old, new, var, point)` gives
/// `limit(d, old, point)` for the solution `d` of `old = new(d)`
/// (`solveset` then `limit`, as `Order._eval_subs` does).
pub type OrderPointHost = fn(&Expr, &Expr, &Expr, &Expr) -> Option<Expr>;

static ORDER_POINT_HOST: std::sync::OnceLock<OrderPointHost> = std::sync::OnceLock::new();

pub fn register_order_point_host(host: OrderPointHost) {
    let _ = ORDER_POINT_HOST.set(host);
}

/// SymPy's `Order._eval_subs` for a substituted order variable: the limit
/// point follows the substitution (`O(x).subs(x, x - 2)` is
/// `O(x - 2, (x, 2))`). An affine replacement is solved here, any other
/// one by the host (`solveset` and `limit`); `None` falls back to
/// substitution in the arguments.
fn order_subs(e: &Expr, old: &Expr, new: &Expr, cache: &mut FxHashMap<usize, (Expr,Expr)>) -> Option<Expr> {
    let (inner, vars) = crate::series::order::order_parts(e)?;
    let i = vars.iter().position(|(v, _)| v == old)?;
    let new_expr = s1(&inner, old, new, cache);
    let mut vars = vars;
    if new.is_symbol() {
        vars[i].0 = new.clone();
    } else {
        let syms = crate::traverse::free_symbols(new);
        let contains_old = syms.contains(old);
        if syms.len() == 1 || contains_old {
            let var = if contains_old { vars[i].0.clone() } else { syms[0].clone() };
            let start = vars[i].1.clone();
            let mut point = subs1(new, &var, &start);
            if point != start {
                // solve old = new(d) for d, then evaluate at old = start
                let slope = crate::diff::diff(new, &var);
                if !slope.is_zero() && !crate::traverse::has(&slope, &var) && start.is_number() {
                    let intercept = subs1(new, &var, &Expr::zero());
                    point = div(&crate::ops::sub(&start, &intercept), &slope);
                } else {
                    point = ORDER_POINT_HOST.get()?(old, new, &var, &start)?;
                }
            }
            vars[i] = (var, point);
        } else {
            vars.remove(i);
        }
    }
    Some(crate::series::order::order(&new_expr, &vars))
}

/// Substitute a single `old -> new` (with algebraic matching of sums,
/// products and powers).
pub fn subs1(e: &Expr, old: &Expr, new: &Expr) -> Expr {
    if let Some(value) = memo::get(e,old,new) { return value; }
    let mut cache = FxHashMap::default();
    let value = s1(e, old, new, &mut cache);
    memo::insert(e,old,new,&value);
    value
}

fn s1(e: &Expr, old: &Expr, new: &Expr, cache: &mut FxHashMap<usize, (Expr,Expr)>) -> Expr {
    if e == old {
        return new.clone();
    }
    // SymPy's Basic._subs consults _eval_subs on every node before falling
    // back to its arguments, so a host hook below defeats the mask filter.
    let hooked = e.has_host_subs();
    // Symbol._eval_subs: a symbol is the power y**1 for an old power of it,
    // so (y).subs(1/y, u) is 1/u.
    let symbol_power = match old.kind() {
        Kind::Pow(b, _) if b.is_symbol() => Some(b),
        _ => None,
    };
    if !hooked && e.is_symbol() {
        if symbol_power == Some(e) {
            if let Some(value) = power_subs(e, &Expr::one(), false, old, new, cache) {
                return value;
            }
        }
        return e.clone();
    }
    if !hooked && (e.is_atom() || (!old.has_binding_identity() && old.mask() & !e.mask() != 0)) {
        return e.clone();
    }
    if let Some((_,v)) = cache.get(&e.dag_id()) {
        return v.clone();
    }
    if hooked {
        if let Some(host) = e.host_handle().filter(|host| host.substitutes()) {
            if let Some(value) = host.subs(e, old, new) {
                cache.insert(e.dag_id(), (e.clone(), value.clone()));
                return value;
            }
        }
        if e.is_atom() { return e.clone(); }
    }
    if matches!(e.kind(),Kind::Func(Func::Integral|Func::Sum|Func::Product,_)) {
        return crate::binders::substitute(e,&[(old.clone(),new.clone())],false);
    }
    if matches!(e.kind(),Kind::Func(Func::Order,_)) {
        if let Some(value)=order_subs(e,old,new,cache) {
            return value;
        }
    }
    if let Kind::Func(Func::Subs,args)=e.kind(){
        let variables=args[1].args();
        if variables.contains(old) && (!new.is_symbol() || crate::traverse::has(&args[2],new) || variables.contains(new)){
            let points=args[2].args().iter().map(|p|s1(p,old,new,cache)).collect();
            return held_subs(&args[0],variables,points);
        }
    }
    // substituting a differentiation variable: keep an unevaluated Subs
    if let Kind::Derivative(_, vars) = e.kind() {
        if vars.iter().any(|(v, _)| v == old) {
            if new.is_symbol() && !crate::traverse::has(e,new) {
                let Kind::Derivative(inner,vars)=e.kind() else{unreachable!()};
                let inner=subs1(inner,old,new);
                let vars=vars.iter().map(|(v,n)|(if v==old{new.clone()}else{v.clone()},*n)).collect();
                return crate::diff::derivative_unevaluated(&inner,vars);
            }
            return held_subs(e,vec![old.clone()],vec![new.clone()]);
        }
        // SymPy: a derivative of the same expression with respect to a
        // sub-multiset of the variables is replaced, leaving the remaining
        // differentiations of new unevaluated.
        if let (Kind::Derivative(inner,vars),Kind::Derivative(old_inner,old_vars))=(e.kind(),old.kind()) {
            if inner==old_inner {
                let count=|list:&[(Expr,u32)],v:&Expr|list.iter().filter(|(w,_)|w==v).map(|(_,n)|*n).sum::<u32>();
                if old_vars.iter().all(|(v,_)|count(old_vars,v)<=count(vars,v)) {
                    let mut rest:Vec<(Expr,u32)>=Vec::new();
                    for (v,_) in vars.iter() {
                        if rest.iter().any(|(w,_)|w==v) {continue;}
                        let remaining=count(vars,v)-count(old_vars,v);
                        if remaining>0 {rest.push((v.clone(),remaining));}
                    }
                    return crate::diff::derivative_unevaluated(new,rest);
                }
            }
        }
        // SymPy's case (1): the substitution may not introduce a
        // differentiation variable into an expression that depends on it; a
        // variable that old already contains is not introduced.
        let Kind::Derivative(inner,_)=e.kind() else{unreachable!()};
        if crate::traverse::has(e,old) && vars.iter().any(|(v,_)|
            crate::traverse::has(new,v) && !crate::traverse::has(old,v) && crate::traverse::has(inner,v)) {
            return held_subs(e,vec![old.clone()],vec![new.clone()]);
        }
        // A variable that is not old can still change, as y in
        // (Derivative(f(y), y)).subs(1/y, u). SymPy's case (2): the new
        // variable may not already be in the expression; case (3): a new
        // variable that is not a valid differentiation variable is held as
        // Subs(Derivative(..., xi_i), xi_i, value).
        let changed:Vec<(Expr,u32)>=vars.iter().map(|(v,n)|(s1(v,old,new,cache),*n)).collect();
        if changed.iter().zip(vars.iter()).any(|(a,b)|a.0!=b.0) {
            if changed.iter().zip(vars.iter()).any(|(a,b)|a.0!=b.0&&a.0.is_symbol()&&crate::traverse::has(inner,&a.0)) {
                return held_subs(e,vec![old.clone()],vec![new.clone()]);
            }
            let mut body=inner.clone();
            let mut wrt=Vec::new();
            let mut dummies=Vec::new();
            let mut points=Vec::new();
            for (i,((value,count),(variable,_))) in changed.iter().zip(vars.iter()).enumerate() {
                if value.is_symbol()||crate::binders::is_function_variable(value)||matches!(value.kind(),Kind::Derivative(..)) {
                    wrt.push((value.clone(),*count));
                    continue;
                }
                let xi=Expr::dummy(&format!("xi_{i}"),crate::symbol::Assumptions::NONE);
                let mut m=FxHashMap::default();
                m.insert(variable.clone(),xi.clone());
                body=xreplace(&body,&m);
                wrt.push((xi.clone(),*count));
                dummies.push(xi);
                points.push(value.clone());
            }
            if !dummies.is_empty() {
                let body=s1(&body,old,new,cache);
                return held_subs(&crate::diff::derivative_unevaluated(&body,wrt),dummies,points);
            }
        }
    }
    let r = match (e.kind(), old.kind()) {
        (Kind::Add(..), Kind::Add(..)) => match add_contains(e, old) {
            Some((k, rest)) => {
                let rest = s1(&rest, old, new, cache);
                add(&[mul2(&Expr::num(k), new), rest])
            }
            None => subs_children(e, old, new, cache),
        },
        (Kind::Mul(..), Kind::Mul(..)) if e.has_noncommutative_host() || old.has_noncommutative_host() =>
            noncommutative_mul_subs(e, old, new, cache).unwrap_or_else(|| subs_children(e, old, new, cache)),
        (Kind::Mul(..), Kind::Mul(..)) => match mul_contains(e, old) {
            Some((count, residual, rest)) => {
                let rest = s1(&rest, old, new, cache);
                // One product, as SymPy's Mul(new, *rest): the residual
                // coefficient is not distributed over a remaining sum.
                mul(&[Expr::num(residual), pow(new, &Expr::int(count)), rest])
            }
            None => {
                // Canonical products store denominator factors with negative
                // powers. Give a positive product in that denominator the
                // same opportunity to match as a product in the numerator.
                // SymPy splits with `fraction`, which leaves sums alone: a
                // substitution in (a/x + b)/7 must not combine it over x.
                let (numerator, denominator) = crate::simptools::radsimp::fraction(e, false);
                if !denominator.is_one() {
                    let n = s1(&numerator, old, new, cache);
                    let d = s1(&denominator, old, new, cache);
                    if n != numerator || d != denominator {
                        div(&n, &d)
                    } else {
                        subs_children(e, old, new, cache)
                    }
                } else {
                    subs_children(e, old, new, cache)
                }
            }
        },
        (Kind::Pow(b, x), _) => power_subs(b, x, false, old, new, cache)
            .unwrap_or_else(|| subs_children(e, old, new, cache)),
        // SymPy's exp._eval_subs views exp(z) and an old exp(w) as E**z and
        // E**w and defers to Pow._eval_subs.
        (Kind::Func(Func::Exp, a), Kind::Const(crate::expr::Constant::E)) => pow(new, &s1(&a[0], old, new, cache)),
        (Kind::Func(Func::Exp, a), Kind::Func(Func::Exp, _)) => power_subs(&Expr::e(), &a[0], true, old, new, cache)
            .unwrap_or_else(|| subs_children(e, old, new, cache)),
        _ => subs_children(e, old, new, cache),
    };
    cache.insert(e.dag_id(), (e.clone(),r.clone()));
    r
}

/// SymPy's Pow._eval_subs for `base**exponent`. `as_exp` marks the
/// exp(z) view, in which an old exp(w) counts as the power E**w.
fn power_subs(base: &Expr, exponent: &Expr, as_exp: bool, old: &Expr, new: &Expr,
    cache: &mut FxHashMap<usize, (Expr,Expr)>) -> Option<Expr> {
    use crate::assumptions::{is_nonnegative, is_positive, is_real, real_signs};
    if old == base {
        return Some(pow(new, &s1(exponent, old, new, cache)));
    }
    // Pow._eval_subs._check for commutative objects: the exponent of old
    // combines with pow when pow is an integer or (b**e)**f == b**(e*f).
    let combine = |own: &Expr, other: &Expr| -> Option<Expr> {
        if old.has_noncommutative_host() { return None; }
        let ratio = div(own, other);
        let integral = ratio.as_number().is_some_and(|n| n.is_integer()
            || n.is_float() && n.to_f64().is_finite() && n.to_f64().fract() == 0.0);
        if integral { return Some(ratio); }
        let (b, e) = crate::decompose::as_base_exp(old);
        let combines = is_positive(&b) == Some(true) && is_real(&e) == Some(true)
            || is_nonnegative(&b) == Some(true) && is_nonnegative(&e) == Some(true);
        combines.then_some(ratio)
    };
    let old_power = match old.kind() {
        Kind::Pow(b, x) => Some((b.clone(), x.clone())),
        Kind::Func(Func::Exp, a) if as_exp => Some((Expr::e(), a[0].clone())),
        _ => None,
    };
    if let Some((old_base, old_exponent)) = &old_power {
        // issue 10829: (4**x).subs(2**x, y) -> y**2 when log(4, 2) is a Number
        if exponent == old_exponent {
            if let Some(l) = rational_log(base, old_base) {
                return Some(pow(new, &Expr::num(l)));
            }
        }
        if base == old_base {
            if !exponent.is_add() {
                let (c1, t1) = symbol_independent(exponent);
                let (c2, t2) = symbol_independent(old_exponent);
                if t1 == t2 {
                    if let Some(k) = combine(&c1, &c2) { return Some(pow(new, &k)); }
                }
            } else {
                // b**(6*x + a).subs(b**(3*x), y) -> y**2 * b**a
                let (c2, t2) = crate::decompose::as_coeff_mul(old_exponent, &[], false);
                let mut factors = Vec::new();
                let mut remaining = Vec::new();
                for term in crate::traverse::ordered_args(exponent) {
                    let term = s1(&term, old, new, cache);
                    let (c1, t1) = crate::decompose::as_coeff_mul(&term, &[], false);
                    if t1 == t2 {
                        if let Some(k) = combine(&c1, &c2) { factors.push(pow(new, &k)); continue; }
                    }
                    remaining.push(term);
                }
                if !factors.is_empty() {
                    let rest = add(&remaining);
                    factors.push(if rest.is_one() { base.clone() } else { pow(base, &rest) });
                    return Some(mul(&factors));
                }
            }
        }
    }
    let old_exp = match old.kind() {
        Kind::Func(Func::Exp, a) => Some(a[0].clone()),
        Kind::Pow(b, x) if b.is_constant(crate::expr::Constant::E) => Some(x.clone()),
        _ => None,
    };
    if let Some(old_exponent) = old_exp {
        if real_signs(exponent).is_some() && is_positive(base) == Some(true) {
            // (2**x).subs(exp(x*log(2)), z) -> z
            let (c1, t1) = symbol_independent(&old_exponent);
            let (c2, t2) = symbol_independent(&mul2(exponent, &crate::functions::log(base)));
            if t1 == t2 {
                if let Some(k) = combine(&c1, &c2) { return Some(pow(new, &k)); }
            }
        }
    }
    None
}

/// `e.as_independent(Symbol, as_Add=False)`.
fn symbol_independent(e: &Expr) -> (Expr, Expr) {
    let dependent = |term: &Expr| !crate::traverse::free_symbols(term).is_empty();
    if !e.is_mul() {
        return if dependent(e) { (Expr::one(), e.clone()) } else { (e.clone(), Expr::one()) };
    }
    let (with, without): (Vec<_>, Vec<_>) = crate::traverse::ordered_args(e).into_iter().partition(dependent);
    (mul(&without), mul(&with))
}

/// `log(arg, base)` when SymPy evaluates it to a Number for positive
/// rationals: multiplicity extracts the powers of the base, and the
/// remaining quotient of logarithms is a Number only for equal integers.
fn rational_log(arg: &Expr, base: &Expr) -> Option<Number> {
    use num_bigint::BigInt;
    let (an, ad) = arg.as_number().filter(|n| n.is_rational() && n.is_positive())?.numer_denom()?;
    let (bn, bd) = base.as_number().filter(|n| n.is_rational() && n.is_positive())?.numer_denom()?;
    if bn == bd { return None; }
    fn count(p: &BigInt, n: &BigInt) -> Option<i64> {
        if *p < BigInt::from(2) { return None; }
        let (mut n, mut k) = (n.clone(), 0);
        while n.is_multiple_of(p) { n /= p; k += 1; }
        Some(k)
    }
    let one = BigInt::from(1);
    let m = if bd == one {
        if an == one { -count(&bn, &ad)? } else { count(&bn, &an)? - count(&bn, &ad)? }
    } else if bn == one {
        count(&bd, &ad)?
    } else {
        count(&bn, &an)?.min(count(&bd, &ad)?) - count(&bd, &an)?.min(count(&bn, &ad)?)
    };
    let a = Number::rational(an, ad)?;
    let b = Number::rational(bn, bd)?;
    let quotient = a.div(&b.pow_i64(m)?)?;
    if quotient.is_one() { return Some(Number::Int(m)); }
    // log(1/q) evaluates to -log(q); any other rational log stays held.
    let normal = |r: &Number| -> (i64, Number) {
        let (p, q) = r.numer_denom().unwrap();
        if p == BigInt::from(1) { (-1, Number::from_bigint(q)) } else { (1, r.clone()) }
    };
    let (sq, iq) = normal(&quotient);
    let (sb, ib) = normal(&b);
    (iq == ib).then(|| Number::Int(m + sq * sb))
}

/// Substitution sees the public expression tree, including coefficients
/// in Add terms and whole powers in Mul factors. map_children deliberately
/// visits the compact internal representation instead, which would miss
/// patterns such as 2*x in 2*x+y, or x**2 in x**2*y.
fn subs_children(e: &Expr, old: &Expr, new: &Expr, cache: &mut FxHashMap<usize, (Expr,Expr)>) -> Expr {
    if e.is_add() || e.is_mul() {
        let args = e.args();
        let replaced: Vec<Expr> = args.iter().map(|c| s1(c, old, new, cache)).collect();
        if replaced.iter().zip(&args).all(|(a,b)|a.same_presentation(b)) {
            return e.clone();
        }
        return if e.is_add() { add(&replaced) } else { mul(&replaced) };
    }
    map_children(e, &mut |c| s1(c, old, new, cache))
}

/// Match a subset of additive terms, including the opposite sign.
/// SymPy matches the coefficients of each term literally: `2*x + 2*y`
/// does not contain `x + y`. Rational constant offsets are independent.
pub fn add_contains(e: &Expr, s: &Expr) -> Option<(Number, Expr)> {
    let (Kind::Add(ce, te), Kind::Add(cs, ts)) = (e.kind(), s.kind()) else {
        return None;
    };
    if ts.len() > te.len() || (!(ce.is_rational() && cs.is_rational()) && ce != cs) {
        return None;
    }
    let lookup: FxHashMap<&Expr, &Number> = te.iter().map(|(t, k)| (t, k)).collect();
    let mut sign: Option<Number> = None;
    for (t, k) in ts {
        let ke = lookup.get(t)?;
        let r = ke.div(k)?;
        if !r.is_one() && !r.is_minus_one() {
            return None;
        }
        match &sign {
            None => sign = Some(r),
            Some(q) if *q == r => {}
            _ => return None,
        }
    }
    let k = sign?;
    let c_rest = ce.sub(&cs.mul(&k));
    let sub_terms: FxHashMap<&Expr, ()> = ts.iter().map(|(t, _)| (t, ())).collect();
    let rest = te.iter().filter(|(t, _)| !sub_terms.contains_key(t)).cloned().collect();
    Some((k, add_from_terms(c_rest, rest)))
}

/// Extract an integer number of copies of a rational exponent, preserving
/// SymPy's extractive divisibility rule for fractional powers.
fn exponent_copies(available: &Expr, wanted: &Expr) -> Option<i64> {
    let a = available.as_number()?.to_ratio()?;
    let b = wanted.as_number()?.to_ratio()?;
    if a.is_integer() && b.is_integer() {
        return (a / b).to_integer().to_i64().filter(|&n| n > 0);
    }
    let aq = a.denom();
    let bq = b.denom();
    if !aq.is_multiple_of(bq) && !bq.is_multiple_of(aq) {
        return None;
    }
    (a / b).to_integer().to_i64().filter(|&n| n > 0)
}

/// Multiplicity of an exact rational factor in a rational coefficient.
/// Numerator and denominator multiplicities are counted independently, so
/// e.g. `3/2` contains one `1/2` and contains no `1/4`.
fn coefficient_copies(available: &Number, wanted: &Number, limit: i64) -> Option<i64> {
    let a = available.to_ratio()?;
    let b = wanted.to_ratio()?;
    let mut numerator = a.numer().clone();
    let mut denominator = a.denom().clone();
    let mut count = 0;
    while count < limit
        && numerator.is_multiple_of(b.numer())
        && denominator.is_multiple_of(b.denom())
    {
        numerator /= b.numer();
        denominator /= b.denom();
        count += 1;
    }
    Some(count)
}

/// SymPy's `Mul._eval_subs` for a product or pattern with factors that may
/// not commute: commutative factors are extracted as from a multiset, a
/// noncommutative pattern replaces consecutive factors whose bases agree
/// in order. `None` falls back to substitution in the arguments.
fn noncommutative_mul_subs(e: &Expr, old: &Expr, new: &Expr, cache: &mut FxHashMap<usize, (Expr,Expr)>) -> Option<Expr> {
    use crate::decompose::{as_base_exp, commutativity};
    let coefficient = |m: &Expr| match m.kind() { Kind::Mul(c, _) if !c.is_one() => Some(c.clone()), _ => None };
    // Keep the replacement literal, so -2*x does not replace 4*x.
    if coefficient(old).is_some_and(|c| c.is_negative()) {
        if let Some(co_self) = coefficient(e) {
            return co_self.is_negative().then(|| subs1(e, &crate::ops::neg(old), &crate::ops::neg(new)));
        }
    }
    enum Nc { Pair(Expr, Number), Done(Expr) }
    fn base_exp(a: &Expr) -> (Expr, Expr) {
        if a.is_pow() || matches!(a.kind(), Kind::Func(Func::Exp, _)) { as_base_exp(a) } else { (a.clone(), Expr::one()) }
    }
    // b**(Rational*e) -> (b**e, Rational); commutative factors are summed
    // by base in first-seen order, noncommutative ones keep their order.
    fn breakup(m: &Expr) -> (Vec<(Expr, Number)>, Vec<(Expr, Number)>) {
        let (mut c, mut nc): (Vec<(Expr, Number)>, Vec<(Expr, Number)>) = (Vec::new(), Vec::new());
        for a in m.args() {
            let (mut b, e) = base_exp(&a);
            let co = match e.kind() { Kind::Num(n) => n.clone(), Kind::Mul(n, _) => n.clone(), _ => Number::ONE };
            if !e.is_one() { b = pow(&b, &div(&e, &Expr::num(co.clone()))); }
            if commutativity(&a) == Some(true) {
                match c.iter_mut().find(|(k, _)| *k == b) { Some(entry) => entry.1 = entry.1.add(&co), None => c.push((b, co)) }
            } else { nc.push((b, co)); }
        }
        (c, nc)
    }
    fn rejoin(b: &Expr, co: &Number) -> Expr {
        let (b, e) = base_exp(b);
        pow(&b, &mul2(&e, &Expr::num(co.clone())))
    }
    // The integer number of times b divides a extractively, else 0.
    fn ndiv(a: &Number, b: &Number) -> i64 {
        let (Some(a), Some(b)) = (a.to_ratio(), b.to_ratio()) else { return 0 };
        if num_traits::Zero::is_zero(&b) || !(a.denom().is_multiple_of(b.denom()) || b.denom().is_multiple_of(a.denom())) { return 0; }
        (a / b).trunc().to_integer().to_i64().unwrap_or(0)
    }
    fn multiplicity(p: &num_bigint::BigInt, n: &num_bigint::BigInt) -> i64 {
        use num_traits::{Signed, Zero};
        let (p, mut n) = (p.abs(), n.abs());
        if p <= num_bigint::BigInt::from(1) || n.is_zero() { return 0; }
        let mut count = 0;
        while (&n % &p).is_zero() { n /= &p; count += 1; }
        count
    }
    fn rational_multiplicity(p: &Number, n: &Number) -> Option<i64> {
        let (pp, pq) = p.numer_denom()?;
        let (np, nq) = n.numer_denom()?;
        let one = num_bigint::BigInt::from(1);
        Some(if pq == one {
            multiplicity(&pp, &np) - multiplicity(&pp, &nq)
        } else if pp == one {
            multiplicity(&pq, &nq)
        } else {
            multiplicity(&pp, &np).min(multiplicity(&pq, &nq)) - multiplicity(&pq, &np).min(multiplicity(&pp, &nq))
        })
    }
    // Give products in the denominator a chance to be changed.
    let mut rv = None;
    let (n, d) = crate::simptools::radsimp::fraction(e, false);
    let mut self2 = e.clone();
    if !d.is_one() {
        self2 = div(&s1(&n, old, new, cache), &s1(&d, old, new, cache));
        if !self2.is_mul() { return Some(s1(&self2, old, new, cache)); }
        if self2 != *e { rv = Some(self2.clone()); }
    }
    let first = |m: &Expr| m.args().into_iter().next();
    let (co_self, co_old) = (first(&self2)?, first(old)?);
    let mut co_xmul = None;
    match (co_self.as_number().filter(|n| n.is_rational()), co_old.as_number().filter(|n| n.is_rational())) {
        (Some(cs), Some(co)) => if cs != co {
            // Expr.extract_multiplicatively of a rational
            co_xmul = cs.div(co).filter(|q| (!cs.is_integer() || q.is_integer()) && !(cs.is_positive() && q.is_negative()));
        },
        (None, Some(_)) => return rv,
        _ => {}
    }
    let (mut c, nc) = breakup(&self2);
    let (old_c, old_nc) = breakup(old);
    let mut co_residual = Number::ONE;
    if let (Some(_), Some(cs), Some(co)) = (&co_xmul, co_self.as_number(), co_old.as_number()) {
        if !co.abs().is_one() {
            let mult = rational_multiplicity(&co.abs(), cs)?;
            let key = Expr::num(cs.clone());
            c.retain(|(b, _)| *b != key);
            let old_key = Expr::num(co.clone());
            match c.iter_mut().find(|(b, _)| *b == old_key) { Some(entry) => entry.1 = entry.1.add(&Number::Int(mult)), None => c.push((old_key, Number::Int(mult))) }
            co_residual = cs.div(&co.pow_i64(mult)?)?;
        }
    }
    let lookup = |list: &[(Expr, Number)], b: &Expr| list.iter().find(|(k, _)| k == b).map(|(_, v)| v.clone());
    let sign = |n: &Number| if n.is_positive() { 1 } else if n.is_negative() { -1 } else { 0 };
    let ok = old_nc.len() <= nc.len() && old_c.len() <= c.len()
        && old_nc.iter().all(|(b, _)| nc.iter().any(|(k, _)| k == b))
        && old_c.iter().all(|(b, _)| lookup(&c, b).is_some())
        && old_c.iter().all(|(b, e)| lookup(&c, b).is_some_and(|ce| sign(&ce) == sign(e)));
    if !ok { return rv; }
    let cdid = if old_c.is_empty() { None } else {
        let mut least = i64::MAX;
        for (b, old_e) in &old_c {
            let r = ndiv(&lookup(&c, b)?, old_e);
            if r == 0 { return rv; }
            least = least.min(r);
        }
        Some(least)
    };
    let mut nc: Vec<Nc> = nc.into_iter().map(|(b, e)| Nc::Pair(b, e)).collect();
    let ncdid = if old_nc.is_empty() {
        for item in nc.iter_mut() { if let Nc::Pair(b, e) = item { *item = Nc::Done(rejoin(b, e)); } }
        None
    } else {
        let mut ncdid = 0i64;
        let take = old_nc.len();
        let mut limit = cdid.unwrap_or(i64::MAX);
        let mut failed = Vec::new();
        let mut i = 0;
        let pair = |item: &Nc| match item { Nc::Pair(b, e) => Some((b.clone(), e.clone())), Nc::Done(_) => None };
        while limit != 0 && i + take <= nc.len() {
            let mut hit = false;
            let mut rat = Vec::new();
            let mut matched = true;
            for j in 0..take {
                let Some((b, e)) = pair(&nc[i + j]) else { matched = false; break };
                if b != old_nc[j].0 { matched = false; break; }
                if j == 0 || j == take - 1 { rat.push(ndiv(&e, &old_nc[j].1)); }
                else if e != old_nc[j].1 { matched = false; break; }
                else { rat.push(1); }
            }
            if matched {
                let mut ndo = rat.iter().copied().min().unwrap_or(0);
                if ndo != 0 {
                    let (b0, e0) = pair(&nc[i]).expect("matched pair");
                    if take == 1 {
                        if let Some(cdid) = cdid { ndo = ndo.min(cdid); }
                        let rest = rejoin(&b0, &e0.sub(&Number::Int(ndo).mul(&old_nc[0].1)));
                        nc[i] = Nc::Done(mul2(&pow(new, &Expr::int(ndo)), &rest));
                    } else {
                        ndo = 1;
                        let left = rejoin(&b0, &e0.sub(&old_nc[0].1));
                        let mid = new.clone();
                        let ir = i + take - 1;
                        let (br, er) = pair(&nc[ir]).expect("matched pair");
                        let right = er.sub(&old_nc[take - 1].1);
                        let replacement = if !right.is_zero() {
                            if i + take < nc.len() { vec![Nc::Done(mul2(&left, &mid)), Nc::Pair(br, right)] }
                            else { vec![Nc::Done(mul(&[left, mid, rejoin(&br, &right)]))] }
                        } else { vec![Nc::Done(mul2(&left, &mid))] };
                        nc.splice(i..i + take, replacement);
                    }
                    limit = limit.saturating_sub(ndo);
                    ncdid += ndo;
                    hit = true;
                }
            }
            if !hit { failed.push(i); }
            i += 1;
        }
        if ncdid == 0 { return rv; }
        failed.extend(i..nc.len());
        for index in failed {
            if let Nc::Pair(b, e) = &nc[index] { let value = s1(&rejoin(b, e), old, new, cache); nc[index] = Nc::Done(value); }
        }
        Some(ncdid)
    };
    let done = match (cdid, ncdid) { (None, Some(n)) => n, (Some(c), None) => c, (Some(c), Some(n)) => c.min(n), (None, None) => 0 };
    let mut margs = Vec::with_capacity(c.len() + 1);
    for (b, e) in &c {
        match lookup(&old_c, b) {
            Some(old_e) => margs.push(rejoin(b, &e.sub(&old_e.mul(&Number::Int(done))))),
            None => margs.push(rejoin(&s1(b, old, new, cache), e)),
        }
    }
    if let Some(cdid) = cdid {
        if ncdid.unwrap_or(0) == 0 { margs.insert(0, pow(new, &Expr::int(cdid))); }
    }
    let nc: Vec<Expr> = nc.into_iter().map(|item| match item { Nc::Done(e) => e, Nc::Pair(b, e) => rejoin(&b, &e) }).collect();
    Some(mul(&[Expr::num(co_residual), mul(&margs), mul(&nc)]))
}

/// Extract as many copies of `s` as occur in `e`, allowing larger powers.
fn mul_contains(e: &Expr, s: &Expr) -> Option<(i64, Number, Expr)> {
    let (Kind::Mul(ce, fe), Kind::Mul(cs, fs)) = (e.kind(), s.kind()) else {
        return None;
    };
    if fs.len() > fe.len() || (cs.is_negative() && !ce.is_negative()) {
        return None;
    }
    let lookup: FxHashMap<&Expr, &Expr> = fe.iter().map(|(b, x)| (b, x)).collect();
    let mut count = i64::MAX;
    for (b, x) in fs {
        count = count.min(exponent_copies(lookup.get(b)?, x)?);
    }
    if count == i64::MAX {
        return None;
    }
    let abs_old = cs.abs();
    if !abs_old.is_one() {
        count = coefficient_copies(ce, &abs_old, count)?;
        if count == 0 {
            return None;
        }
    }
    // Respect the sign when a negative pattern is extracted repeatedly.
    let residual = ce.div(&cs.pow_i64(count)?)?;
    let wanted: FxHashMap<&Expr, &Expr> = fs.iter().map(|(b, x)| (b, x)).collect();
    let mut rest = Vec::with_capacity(fe.len());
    for (b, x) in fe {
        let exponent = match wanted.get(b) {
            Some(old_x) => crate::ops::sub(x, &mul2(&Expr::int(count), old_x)),
            None => x.clone(),
        };
        rest.push(mul_factor(b, &exponent));
    }
    Some((count, residual, mul(&rest)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn s(text: &str) -> Expr {
        parse_str(text).unwrap()
    }

    #[test]
    fn additive_subsets_and_constant_offsets() {
        for (expr, old, expected) in [
            ("x + y + z", "x + y", "w + z"),
            ("x - y + z", "y - x", "z - w"),
            ("x + y", "x + y + 1", "w - 1"),
            ("x + y + 2", "x + y - 1", "w + 3"),
            ("2*x + 2*y + z", "x + y", "2*x + 2*y + z"),
        ] {
            assert_eq!(subs1(&s(expr), &s(old), &s("w")), s(expected), "{expr} subs {old}");
        }
    }

    #[test]
    fn multiplicative_subsets_and_powers() {
        for (expr, old, expected) in [
            ("x**2*y", "x**2", "w*y"),
            ("x**4*y", "x**2", "w**2*y"),
            ("2*x + y", "2*x", "w + y"),
            ("2*x*y", "2*x", "w*y"),
            ("x**2*y*z", "x*y", "w*x*z"),
            ("x**3*y**2*z", "x*y", "w**2*x*z"),
            ("6*x**3*y**2", "2*x*y", "3*w*x**2*y"),
            ("8*x**3*y**3", "2*x*y", "w**3"),
            ("x**2*y**3", "2*x*y", "x**2*y**3"),
            ("3*x**3*y**3/2", "x*y/2", "3*w*x**2*y**2"),
            ("3*x**3*y**3/2", "x*y/4", "3*x**3*y**3/2"),
            ("2*x**3*y**3", "x*y/2", "2*x**3*y**3"),
            ("-x*y*z", "-x*y", "w*z"),
            ("x*y*z", "-x*y", "x*y*z"),
            ("z/(x*y)", "x*y", "z/w"),
            ("z/(x**3*y**2)", "x*y", "z/(w**2*x)"),
            ("x*z/y", "x/y", "w*z"),
            ("(3*a*b + c/x)/7", "3*a", "b*w/7 + c/(7*x)"),
            ("3*a*(x + 5)/7", "3*a", "w*(x + 5)/7"),
        ] {
            assert_eq!(subs1(&s(expr), &s(old), &s("w")), s(expected), "{expr} subs {old}");
        }
    }
}
