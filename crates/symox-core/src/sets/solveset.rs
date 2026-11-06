//! `solveset` (on top of [`crate::solve::solve`], with `ImageSet`s for
//! periodic solutions) and univariate inequalities by sign analysis.

use std::cmp::Ordering;

use super::*;
use crate::expr::{Constant, RelOp};
use crate::functions::is_numeric_value;
use crate::ops::{add2, div, mul2, neg, sub};
use crate::subs::subs1;
use crate::traverse::has;

fn linear(u: &Expr, x: &Expr) -> Option<(Expr, Expr)> {
    if !has(u, x) {
        return None;
    }
    let a = crate::diff::diff(u, x);
    if has(&a, x) || a.is_zero() {
        return None;
    }
    let b = crate::expand::expand(&sub(u, &mul2(&a, x)));
    if has(&b, x) {
        return None;
    }
    Some((a, b))
}

fn proved_nonzero(value:&Expr)->bool{
    if crate::assumptions::exact_queries(||crate::assumptions::is_zero(value))==Some(false){return true;}
    if let Ok(Some((real,imaginary)))=crate::complex_parts::as_real_imag_bounded(value,false,None){
        return [real,imaginary].iter().any(|part|matches!(crate::special::lambert::real_sign(part),Some(-1|1)));
    }
    false
}

fn is_real_domain(domain: &Expr) -> bool {
    is_subset(domain, &reals()) == Some(true)
}

/// The integer dummy of `ImageSet(Lambda(_n, ...), Integers)`.
pub fn dummy_n() -> Expr {
    // Fresh per image; the same atom is reused by its Lambda body and binder.
    // A user parameter Symbol("_n") must remain free.
    Expr::dummy("n", crate::symbol::Assumptions::INTEGER)
}

/// Every periodic branch of one solve uses the same bound integer atom.
/// This keeps presentation independent of branch discovery/allocation order.
#[derive(Default)]
struct SolveContext { periodic_variable: std::sync::OnceLock<Expr>,
    /// The generator occurs once and linearly (the reference `_solve_trig`
    /// inverts only then).
    single: std::cell::Cell<bool> }
impl SolveContext {
    fn periodic_variable(&self)->&Expr { self.periodic_variable.get_or_init(dummy_n) }
}

/// The only periodic elementary generator in `f`, with
/// a linear argument.
fn single_generator(f: &Expr, x: &Expr) -> Option<Expr> {
    fn rec(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
        if !has(e, x) {
            return;
        }
        if let Kind::Func(..) = e.kind() {
            if !out.contains(e) {
                out.push(e.clone());
            }
            return;
        }
        for a in e.args() {
            rec(&a, x, out);
        }
    }
    let mut v = Vec::new();
    rec(f, x, &mut v);
    if v.len() != 1 {
        return None;
    }
    let g = v.pop().unwrap();
    let Kind::Func(fun, a) = g.kind() else { return None };
    if !matches!(fun, Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc | Func::Exp
        | Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth) || a.len() != 1 {
        return None;
    }
    linear(&a[0], x)?;
    Some(g)
}

fn has_periodic_in(e: &Expr, x: &Expr) -> bool {
    match e.kind() {
        Kind::Func(f, a) if matches!(f, Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc
            | Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth) => {
            a.iter().any(|t| has(t, x))
        }
        _ => e.args().iter().any(|a| has_periodic_in(a, x)),
    }
}

/// Prefer a certified representative in `[0,p)`. An unnormalized center
/// defines the same integer lattice, so failure to normalize loses no roots.
fn reduce_mod(c: &Expr, p: &Expr) -> Option<Expr> {
    let ratio=div(c,p);
    let integer=if let Some(number)=ratio.as_number(){
        Some(number.floor())
    }else{crate::special::lambert::integer_round(&ratio,-1).map(crate::number::Number::from_bigint)};
    Some(match integer {
        Some(integer)=>crate::expand::expand(&sub(c,&mul2(&Expr::num(integer),p))),
        None=>c.clone(),
    })
}

#[cfg(test)]
mod periodic_normalization_tests{
    use super::*;
    #[test]
    fn integer_lattice_normalization_has_no_machine_integer_cutoff(){
        let integer=num_bigint::BigInt::from(10u32).pow(100);
        let c=mul2(&add2(&Expr::num(crate::number::Number::from_bigint(integer)),&Expr::rational(1,3)),&Expr::pi());
        assert_eq!(reduce_mod(&c,&mul2(&Expr::int(2),&Expr::pi())),Some(div(&Expr::pi(),&Expr::int(3))));
        assert_eq!(reduce_mod(&Expr::i(),&Expr::pi()),Some(Expr::i()));
    }
}

fn image(period: &Expr, c: &Expr, context: &SolveContext) -> Expr {
    let n = context.periodic_variable();
    image_set(&lambda(n, &add2(&mul2(period, n), c)), &integers())
}

/// Intersect an entire complex lattice with a real domain by its imaginary
/// component, retaining the exact intersection whenever the index is undecided.
fn restricted_image(period:&Expr,c:&Expr,domain:&Expr,context:&SolveContext)->Expr{
    let lattice=image(period,c,context);
    if !is_real_domain(domain){return lattice;}
    let fallback=||intersection(vec![lattice.clone(),domain.clone()]);
    let (Ok(Some((_,step))),Ok(Some((_,offset))))=(
        crate::complex_parts::as_real_imag(period,false,None),
        crate::complex_parts::as_real_imag(c,false,None)) else{return fallback()};
    let zero=|value:&Expr|crate::assumptions::exact_queries(||crate::assumptions::is_zero(value));
    if zero(&step)==Some(true){
        if zero(&offset)==Some(true){
            return if is_op(domain,SetOp::Reals){lattice}else{fallback()};
        }
        if proved_nonzero(&offset){return empty_set();}
    }else if proved_nonzero(&step){
        let index=div(&neg(&offset),&step);
        match crate::assumptions::exact_queries(||crate::assumptions::is_integer(&index)){
            Some(true)=>return filter_domain(vec![crate::expand::expand(&add2(&mul2(period,&index),c))],domain),
            Some(false)=>return empty_set(),
            _=>{},
        }
    }
    fallback()
}

/// The reference `_solve_trig` for a single trig or hyperbolic function of
/// `a*x + b` with a symbolic `a`: the lattice of `trig(u) = k` is inverted
/// through `u = a*x + b` without expansion, as `imageset(Lambda(n, (n -
/// b)/a), ...)`, and guarded by `Ne(numer(a), 0) & Ne(denom(a), 0)`.
fn symbolic_coefficient(f: &Func, a: &Expr, b: &Expr, k: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    if !context.single.get() || crate::traverse::free_symbols(a).is_empty() || matches!(f, Func::Exp) {
        return None;
    }
    let (num, den) = crate::simplify::numer_denom(a);
    let mut cond = crate::logic::and(vec![crate::relational::ne(&num, &Expr::zero()),
        crate::relational::ne(&den, &Expr::zero())]);
    if !is_numeric_value(k) {
        // `_invert_trig_hyp_complex` with a symbolic right-hand side:
        // `2*n*pi + asin(k)`, `2*n*pi + pi - asin(k)`; `2*n*pi + acos(k)`,
        // `2*n*pi - acos(k)`; `n*pi + atan(k)`. Over the reals `sin` and
        // `cos` add the condition `-1 <= k <= 1` (`_invert_trig_hyp_real`).
        use crate::functions::{acos, asin, atan};
        let n = context.periodic_variable().clone();
        let two_pi_n = mul2(&mul2(&Expr::int(2), &Expr::pi()), &n);
        let bodies = match f {
            Func::Sin => vec![add2(&two_pi_n, &asin(k)), crate::ops::add(&[two_pi_n.clone(), Expr::pi(), neg(&asin(k))])],
            Func::Cos => vec![add2(&two_pi_n, &acos(k)), sub(&two_pi_n, &acos(k))],
            Func::Tan if !is_real_domain(domain) => vec![add2(&mul2(&Expr::pi(), &n), &atan(k))],
            _ => return None,
        };
        if is_real_domain(domain) {
            if !is_op(domain, SetOp::Reals) {
                return None;
            }
            cond = crate::logic::and(vec![crate::relational::le(&Expr::minus_one(), k),
                crate::relational::le(k, &Expr::one()), cond]);
        }
        let images = bodies.iter().map(|body| image_set(&lambda(&n, &div(&sub(body, b), a)), &integers())).collect();
        return Some(condition_set(x, &cond, &union(images)));
    }
    if !is_numeric_value(k) {
        return None;
    }
    let u = Expr::dummy("solveset_u", crate::symbol::Assumptions::NONE);
    let unit = periodic(&crate::functions::func(f.clone(), vec![u.clone()]), k, &u, &complexes(), context)?;
    let n = context.periodic_variable().clone();
    let pieces = if is_op(&unit, SetOp::Union) { unit.args() } else { vec![unit.clone()] };
    let mut out = Vec::new();
    for piece in pieces {
        if is_empty_set(&piece) {
            continue;
        }
        if !is_op(&piece, SetOp::ImageSet) {
            return None;
        }
        let parts = piece.args();
        if parts.len() != 2 || !is_op(&parts[1], SetOp::Integers) {
            return None;
        }
        let body = crate::sets::apply_lambda(&parts[0], &n)?;
        if has(&body, &u) {
            return None;
        }
        out.push(image_set(&lambda(&n, &div(&sub(&body, b), a)), &integers()));
    }
    Some(condition_set(x, &cond, &union(out)))
}

/// SymPy's `_invert_trig_hyp_complex` lattices for `sec`, `csc` and `cot`
/// when the reciprocal inverse stays unevaluated: `2*n*pi +- asec(k)`,
/// `2*n*pi + acsc(k)` / `2*n*pi + pi - acsc(k)` and `n*pi + acot(k)`.
/// Over the reals this is used only for a real argument and a real `k`
/// outside `(-1, 1)` (every lattice point is then real); other cases keep
/// the reduction to `cos`, `sin` and `tan`.
#[allow(clippy::too_many_arguments)]
fn reciprocal_trig_images(f: &Func, k: &Expr, a: &Expr, b: &Expr, domain: &Expr, real: bool,
    real_argument: bool, context: &SolveContext) -> Option<Expr> {
    if !is_numeric_value(k) {
        return None;
    }
    let pi = Expr::pi();
    let inverse = match f {
        Func::Sec => crate::special::asec(k),
        Func::Csc => crate::special::acsc(k),
        _ => crate::functions::acot(k),
    };
    let unevaluated = match (f, inverse.kind()) {
        (Func::Sec, Kind::Func(Func::Special(crate::special::SpecialFn::Asec), _))
        | (Func::Csc, Kind::Func(Func::Special(crate::special::SpecialFn::Acsc), _))
        | (Func::Cot, Kind::Func(Func::Acot, _)) => true,
        _ => false,
    };
    if !unevaluated {
        return None;
    }
    if real {
        if !is_op(domain, SetOp::Reals) || !real_argument
            || crate::assumptions::exact_queries(|| crate::assumptions::is_real(k)) != Some(true) {
            return None;
        }
        if !matches!(f, Func::Cot) {
            let above = crate::special::lambert::real_sign(&sub(k, &Expr::one()));
            let below = crate::special::lambert::real_sign(&add2(k, &Expr::one()));
            if above != Some(1) && below != Some(-1) {
                return None;
            }
        }
    }
    let (thetas, base) = match f {
        Func::Sec => (vec![inverse.clone(), neg(&inverse)], mul2(&Expr::int(2), &pi)),
        Func::Csc => (vec![inverse.clone(), sub(&pi, &inverse)], mul2(&Expr::int(2), &pi)),
        _ => (vec![inverse], pi),
    };
    let period = div(&base, a);
    let period = if crate::special::lambert::real_sign(&period) == Some(-1) { neg(&period) } else { period };
    Some(union(thetas.iter().map(|t| image(&period, &crate::expand::expand(&div(&sub(t, b), a)), context)).collect()))
}

/// Periodic solutions of `trig(a*x + b) = k` / `exp(a*x + b) = k`.
fn periodic(g: &Expr, k: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    use crate::functions::{acos, asin, atan, log, asinh, acosh, atanh};
    let Kind::Func(f, args) = g.kind() else { return None };
    if args.len() != 1 {
        return None;
    }
    let (a, b) = linear(&args[0], x)?;
    if has(k, x) {
        return None;
    }
    if !proved_nonzero(&a){return symbolic_coefficient(f, &a, &b, k, x, domain, context);}
    let pi = Expr::pi();
    let two_pi = mul2(&Expr::int(2), &pi);
    let real = is_real_domain(domain);
    let real_argument=real&&crate::assumptions::exact_queries(||crate::assumptions::is_real(&a))==Some(true)
        &&crate::assumptions::exact_queries(||crate::assumptions::is_real(&b))==Some(true);
    if matches!(f,Func::Sec|Func::Csc|Func::Cot){
        if k.is_zero(){
            if !matches!(f,Func::Cot){return Some(empty_set());}
            return periodic(&crate::functions::func(Func::Cos,args.clone()),k,x,domain,context);
        }
        if !proved_nonzero(k){return None;}
        if let Some(images)=reciprocal_trig_images(f,k,&a,&b,domain,real,real_argument,context){
            return Some(images);
        }
        let head=match f{Func::Sec=>Func::Cos,Func::Csc=>Func::Sin,_=>Func::Tan};
        return periodic(&crate::functions::func(head,args.clone()),&div(&Expr::one(),k),x,domain,context);
    }
    if matches!(f, Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth) {
        // A real affine argument admits the real inverse branches. For
        // complex arguments retain the full imaginary period lattice.
        if real_argument {
            if crate::assumptions::exact_queries(||crate::assumptions::is_real(k))==Some(false) {
                return Some(empty_set());
            }
            if matches!(f,Func::Cosh) && crate::special::lambert::real_sign(&sub(k,&Expr::one()))==Some(-1) {
                return Some(empty_set());
            }
            if matches!(f,Func::Tanh) && (matches!(crate::special::lambert::real_sign(&add2(k,&Expr::one())),Some(-1|0))
                || matches!(crate::special::lambert::real_sign(&sub(k,&Expr::one())),Some(0|1))) {
                return Some(empty_set());
            }
            if matches!(f,Func::Coth) && matches!(crate::special::lambert::real_sign(&add2(k,&Expr::one())),Some(0|1))
                && matches!(crate::special::lambert::real_sign(&sub(k,&Expr::one())),Some(-1|0)) {
                return Some(empty_set());
            }
            let inverse=match f {Func::Sinh=>asinh(k),Func::Cosh=>acosh(k),
                Func::Coth=>crate::special::acoth(k),_=>atanh(k)};
            let mut roots=vec![crate::expand::expand(&div(&sub(&inverse,&b),&a))];
            if matches!(f,Func::Cosh) {
                roots.push(crate::expand::expand(&div(&sub(&neg(&inverse),&b),&a)));
            }
            return Some(filter_domain(roots,domain));
        }
        if matches!(f,Func::Tanh | Func::Coth) && (k.is_one() || k.is_minus_one()) {
            return Some(empty_set());
        }
        let imaginary_pi=mul2(&pi,&Expr::i());
        let (centers,base_period)=match f {
            Func::Sinh=>{
                let v=asinh(k);
                (vec![v.clone(),sub(&imaginary_pi,&v)],mul2(&Expr::int(2),&imaginary_pi))
            },
            Func::Cosh=>{
                let v=acosh(k);
                (vec![v.clone(),neg(&v)],mul2(&Expr::int(2),&imaginary_pi))
            },
            Func::Coth=>(vec![crate::special::acoth(k)],imaginary_pi),
            _=>{let v=atanh(k);(vec![v],imaginary_pi)},
        };
        let period=div(&base_period,&a);
        let images=union(centers.iter().map(|center|restricted_image(&period,
            &crate::expand::expand(&div(&sub(center,&b),&a)),domain,context)).collect());
        if matches!(f,Func::Tanh | Func::Coth) && !is_numeric_value(k) {
            let condition=crate::logic::and(vec![crate::relational::ne(k,&Expr::one()),
                crate::relational::ne(k,&Expr::minus_one())]);
            return Some(condition_set(x,&condition,&images));
        }
        return Some(images);
    }
    let thetas: Vec<Expr>;
    let base_period: Expr;
    match f {
        Func::Sin | Func::Cos => {
            if !is_numeric_value(k) {
                if !real {
                    return None;
                }
                // ConditionSet(x, -1 <= k <= 1, general solutions)
                let ts = if matches!(f, Func::Sin) {
                    let s = asin(k);
                    vec![s.clone(), sub(&pi, &s)]
                } else {
                    let c = acos(k);
                    vec![c.clone(), sub(&two_pi, &c)]
                };
                let period = div(&two_pi, &a);
                let imgs: Vec<Expr> = ts
                    .iter()
                    .map(|t| image(&period, &crate::expand::expand(&div(&sub(t, &b), &a)), context))
                    .collect();
                let cond = crate::logic::and(vec![
                    crate::relational::le(&Expr::minus_one(), k),
                    crate::relational::le(k, &Expr::one()),
                ]);
                return Some(condition_set(x, &cond, &union(imgs)));
            }
            if real_argument {
                if crate::assumptions::exact_queries(||crate::assumptions::is_real(k))==Some(false){
                    return Some(empty_set());
                }
                let above=crate::special::lambert::real_sign(&sub(k,&Expr::one()));
                let below=crate::special::lambert::real_sign(&add2(k,&Expr::one()));
                if above==Some(1)||below==Some(-1){return Some(empty_set());}
            }
            let (t1, t2) = if matches!(f, Func::Sin) {
                let s = asin(k);
                (s.clone(), sub(&pi, &s))
            } else {
                let c = acos(k);
                (c.clone(), sub(&two_pi, &c))
            };
            thetas = vec![t1, t2];
            base_period = two_pi.clone();
        }
        Func::Tan => {
            if !is_numeric_value(k){return None;}
            if k==&Expr::i()||k==&neg(&Expr::i()){return Some(empty_set());}
            if real_argument&&crate::assumptions::exact_queries(||crate::assumptions::is_real(k))==Some(false){return Some(empty_set());}
            thetas = vec![atan(k)];
            base_period = pi.clone();
        }
        Func::Exp => {
            if real_argument {
                // a*x + b = log(k) for k > 0
                if is_numeric_value(k) {
                    if crate::assumptions::exact_queries(||crate::assumptions::is_real(k))==Some(false){return Some(empty_set());}
                    match crate::special::lambert::real_sign(k){
                        Some(1)=>return Some(filter_domain(vec![crate::expand::expand(&div(&sub(&log(k),&b),&a))],domain)),
                        Some(0|-1)=>return Some(empty_set()),
                        _=>{},
                    }
                }else{
                    return Some(filter_domain(vec![crate::expand::expand(&div(&sub(&log(k), &b), &a))], domain));
                }
            }
            if k.is_zero() {
                return Some(empty_set());
            }
            // a*x + b = log(k) + 2*n*pi*I
            let period = div(&mul2(&two_pi, &Expr::i()), &a);
            let c = crate::expand::expand(&div(&sub(&log(k), &b), &a));
            return Some(restricted_image(&period,&c,domain,context));
        }
        _ => return None,
    }
    let period = div(&base_period, &a);
    let period = if crate::special::lambert::real_sign(&period)==Some(-1){neg(&period)}else{period};
    let mut cs: Vec<Expr> = Vec::new();
    for t in thetas {
        let t = reduce_mod(&t, &base_period)?;
        let c = crate::expand::expand(&div(&sub(&t, &b), &a));
        let c = reduce_mod(&c, &period)?;
        if !cs.iter().any(|d| cmp3(d, &c) == Some(Ordering::Equal)) {
            cs.push(c);
        }
    }
    Some(union(cs.iter().map(|c|restricted_image(&period,c,domain,context)).collect()))
}

fn any_node(e: &Expr, pred: &dyn Fn(&Expr) -> bool) -> bool {
    let mut found = false;
    crate::traverse::preorder(e, &mut |a| found |= pred(a));
    found
}

/// The exponent of an exponential atom: `u` for `exp(u)` and `u*log(b)` for
/// `b**u` with a positive rational `b`. `Err` for an `x` dependent exponent
/// that is not of this kind (a symbolic or `x` dependent base).
fn exponential_exponent(e: &Expr, x: &Expr) -> Result<Option<(Expr, bool)>, ()> {
    match e.kind() {
        Kind::Func(Func::Exp, a) if a.len() == 1 && has(&a[0], x) => Ok(Some((a[0].clone(), false))),
        Kind::Pow(b, ex) if has(ex, x) => {
            if has(b, x) {
                return Err(());
            }
            if *b == Expr::e() {
                return Ok(Some((ex.clone(), false)));
            }
            if b.as_number().is_some_and(|n| n.is_rational() && n.is_positive()) {
                return Ok(Some((mul2(ex, &crate::functions::log(b)), true)));
            }
            Err(())
        }
        _ => Ok(None),
    }
}

/// The distinct exponential atoms of `e` with their exponents; `false` if
/// an exponent depends on `x` in an unsupported way.
fn collect_exponentials(e: &Expr, x: &Expr, out: &mut Vec<(Expr, Expr, bool)>) -> bool {
    if !has(e, x) {
        return true;
    }
    match exponential_exponent(e, x) {
        Err(()) => return false,
        Ok(Some((u, rational_base))) => {
            if !out.iter().any(|(a, ..)| a == e) {
                out.push((e.clone(), u, rational_base));
            }
            return true;
        }
        Ok(None) => {}
    }
    e.args().iter().all(|a| collect_exponentials(a, x, out))
}

/// `exp(u)**k -> exp(k*u)` for a rational `k`, as SymPy evaluates it.
fn merge_exp_powers(e: &Expr) -> Expr {
    crate::simptools::bottom_up(e, &mut |rv| match rv.kind() {
        Kind::Pow(b, k) if k.as_number().is_some_and(|n| n.is_rational()) => match b.kind() {
            Kind::Func(Func::Exp, a) if a.len() == 1 => crate::functions::func(Func::Exp, vec![crate::expand::expand(&mul2(k, &a[0]))]),
            _ => rv.clone(),
        },
        _ => rv.clone(),
    })
}

/// `x` occurs in an exponent (`exp(u(x))`, `b**u(x)`).
fn has_exponential_in(e: &Expr, x: &Expr) -> bool {
    if !has(e, x) {
        return false;
    }
    match e.kind() {
        Kind::Func(Func::Exp, a) if a.iter().any(|t| has(t, x)) => true,
        Kind::Pow(_, ex) if has(ex, x) => true,
        _ => e.args().iter().any(|a| has_exponential_in(a, x)),
    }
}

/// `x` occurs other than polynomially (inside a function or an exponent).
fn is_transcendental_in(e: &Expr, x: &Expr) -> bool {
    if !has(e, x) {
        return false;
    }
    match e.kind() {
        Kind::Func(..) => true,
        Kind::Pow(b, ex) => has(ex, x) || ex.as_number().is_none_or(|n| !n.is_integer()) || is_transcendental_in(b, x),
        _ => e.args().iter().any(|a| is_transcendental_in(a, x)),
    }
}

/// The canonical form that SymPy's `imageset(Lambda(n, expr), Integers)`
/// gives to an expression linear in `n` (`_set_function(Lambda, Integers)`).
fn canonical_integer_lattice(expr: &Expr, n: &Expr) -> Expr {
    let c = subs1(expr, n, &Expr::zero());
    let fx = sub(expr, &c);
    let f_x = sub(&subs1(expr, n, &neg(n)), &c);
    let neg_count = |e: &Expr| -> usize {
        match e.kind() {
            Kind::Add(..) => e.args().iter().filter(|t| crate::ops::could_extract_minus_sign(t)).count(),
            _ => usize::from(crate::ops::could_extract_minus_sign(e)),
        }
    };
    let expr = if neg_count(&f_x) < neg_count(&fx) { add2(&f_x, &c) } else { expr.clone() };
    // `expr.match(a*n + b)`: every term is free of `n` or has `n` as a
    // plain factor
    let terms = match expr.kind() {
        Kind::Add(..) => expr.args(),
        _ => vec![expr.clone()],
    };
    let mut a_terms = Vec::new();
    let mut b_terms = Vec::new();
    for t in &terms {
        if !has(t, n) {
            b_terms.push(t.clone());
            continue;
        }
        if t == n {
            a_terms.push(Expr::one());
            continue;
        }
        let Kind::Mul(..) = t.kind() else { return expr };
        let factors = t.args();
        let (with, without): (Vec<Expr>, Vec<Expr>) = factors.into_iter().partition(|f| has(f, n));
        if with.len() != 1 || &with[0] != n {
            return expr;
        }
        a_terms.push(crate::ops::mul(&without));
    }
    let a = crate::ops::add(&a_terms);
    if a.is_zero() || any_node(&expr, &|e| e.as_number().is_some_and(|n| n.is_float())) {
        return expr;
    }
    if a.is_one() || a.is_minus_one() {
        b_terms.retain(|t| crate::assumptions::is_integer(t) != Some(true));
    }
    let mut b = crate::ops::add(&b_terms);
    let numeric = |e: &Expr| crate::traverse::free_symbols(e).is_empty();
    if numeric(&b) {
        if let Ok(Some((br, bi))) = crate::complex_parts::as_real_imag(&b, false, None) {
            let real_a = crate::assumptions::exact_queries(|| crate::assumptions::is_real(&a)) == Some(true);
            if real_a && numeric(&a) {
                if !br.is_zero() {
                    if let Some(r) = mod_numeric(&br, &a) {
                        b = add2(&r, &mul2(&Expr::i(), &bi));
                    }
                }
            } else {
                let ai = div(&a, &Expr::i());
                if crate::assumptions::exact_queries(|| crate::assumptions::is_real(&ai)) == Some(true)
                    && numeric(&ai) && !bi.is_zero() {
                    if let Some(r) = mod_numeric(&bi, &ai) {
                        b = add2(&br, &mul2(&Expr::i(), &r));
                    }
                }
            }
        }
    }
    add2(&mul2(&a, n), &b)
}

/// `c % p` for real numeric `c` and `p`, with a certified floor.
fn mod_numeric(c: &Expr, p: &Expr) -> Option<Expr> {
    let ratio = div(c, p);
    let integer = match ratio.as_number() {
        Some(number) => number.floor(),
        None => match crate::special::lambert::integer_round(&ratio, -1) {
            Some(i) => crate::number::Number::from_bigint(i),
            None => {
                let value = crate::evalf_mp::evalf(&ratio, 80);
                let v = value.as_number()?;
                let f = v.floor();
                // refuse a value too close to an integer to be certified
                let low = crate::evalf_mp::evalf(&sub(&ratio, &Expr::num(f.clone())), 80);
                let low = low.as_number()?;
                let eps = crate::number::Number::rational(1.into(), num_bigint::BigInt::from(10u32).pow(30))?;
                let one_minus = crate::number::Number::rational(num_bigint::BigInt::from(10u32).pow(30) - 1, num_bigint::BigInt::from(10u32).pow(30))?;
                if low.cmp_value(&eps) != Some(Ordering::Greater) || low.cmp_value(&one_minus) != Some(Ordering::Less) {
                    return None;
                }
                f
            }
        },
    };
    Some(crate::expand::expand(&sub(c, &mul2(&Expr::num(integer), p))))
}

/// `ImageSet(Lambda(n, expr), Integers)` restricted to `domain`.
fn lattice_in_domain(expr: &Expr, n: &Expr, domain: &Expr, context: &SolveContext) -> Expr {
    let lattice = image_set(&lambda(n, expr), &integers());
    if !is_real_domain(domain) {
        return if is_op(domain, SetOp::Complexes) { lattice } else { intersection(vec![lattice, domain.clone()]) };
    }
    match linear(expr, n) {
        Some((period, c)) => restricted_image(&period, &c, domain, context),
        None => intersection(vec![lattice, domain.clone()]),
    }
}

/// The roots of the polynomial `p` in `t` with multiplicities, when all of
/// them are found in radicals.
fn complete_roots(p: &Expr, t: &Expr) -> Option<Vec<Expr>> {
    let degree = crate::poly::api::degree(p, Some(t)).ok()??;
    if degree == 0 {
        return Some(Vec::new());
    }
    let roots = crate::poly::roots::roots(p, std::slice::from_ref(t)).ok()?;
    if roots.iter().map(|(_, m)| *m).sum::<u32>() != degree
        || roots.iter().any(|(r, _)| has(r, t) || any_node(r, &|e| crate::poly::rootof::is_rootof(e))) {
        return None;
    }
    Some(roots.into_iter().map(|(r, _)| r).collect())
}

/// `_solve_as_poly` for an equation that is a polynomial in a single
/// exponential generator `exp(g*x)`, inverted with `_invert_complex`:
/// `exp(g*x) = k` gives `x = (I*(2*n*pi + arg(k)) + log(Abs(k)))/g`.
/// `num` and `den` are the numerator and denominator of the equation.
/// Rational ratios between the exponents are reduced to the generator
/// with the gcd of the exponents. `None` when the method does not apply or
/// some roots cannot be found.
fn solve_exponential_poly(num: &Expr, den: &Expr, x: &Expr, domain: &Expr, context: &SolveContext,
    allow_symbolic: bool, symbolic_trig: bool) -> Option<Expr> {
    use num_integer::Integer;
    let num = &merge_exp_powers(num);
    let den = &merge_exp_powers(den);
    let mut atoms = Vec::new();
    if !collect_exponentials(num, x, &mut atoms) || atoms.is_empty() {
        return None;
    }
    if !collect_exponentials(den, x, &mut atoms) {
        return None;
    }
    let mut affine = Vec::new();
    for (_, u, rational_base) in &atoms {
        let u = if *rational_base { crate::simp::expand_log(u, true) } else { u.clone() };
        affine.push(linear(&u, x)?);
    }
    let a0 = affine[0].0.clone();
    let mut ratios = Vec::new();
    let (mut g_num, mut l_den) = (num_bigint::BigInt::from(0), num_bigint::BigInt::from(1));
    for (a, _) in &affine {
        let r = crate::expand::expand(&div(a, &a0));
        let (p, q) = r.as_number().filter(|n| n.is_rational())?.numer_denom()?;
        g_num = g_num.gcd(&p);
        l_den = l_den.lcm(&q);
        ratios.push((p, q));
    }
    let mut g = mul2(&a0, &div(&Expr::num(crate::number::Number::from_bigint(g_num.clone())),
        &Expr::num(crate::number::Number::from_bigint(l_den.clone()))));
    let mut flip = false;
    if crate::ops::could_extract_minus_sign(&g) {
        g = neg(&g);
        flip = true;
    }
    // `_solve_trig1` with a symbolic `mu`: the solutions are guarded by
    // the nonvanishing of the coefficient
    let guarded = symbolic_trig && !proved_nonzero(&g) && crate::assumptions::is_zero(&g) != Some(true)
        && !crate::traverse::free_symbols(&g).is_empty();
    if !proved_nonzero(&g) && !(allow_symbolic && crate::assumptions::is_zero(&g) != Some(true)) && !guarded {
        return None;
    }
    // the generator is `exp(g*x + c)`, the first exponent scaled by `g/a0`
    let mut c = crate::expand::expand(&mul2(&div(&g, &a0), &affine[0].1));
    let t = Expr::dummy("solveset_exp", crate::symbol::Assumptions::NONE);
    let mut map = rustc_hash::FxHashMap::default();
    for (((atom, ..), (_, b)), (p, q)) in atoms.iter().zip(&affine).zip(&ratios) {
        let mut k = (p * &l_den / q) / &g_num;
        if flip {
            k = -k;
        }
        let k = i64::try_from(k).ok()?;
        if k.abs() > 64 {
            return None;
        }
        let shift = crate::expand::expand(&sub(b, &mul2(&Expr::int(k), &c)));
        map.insert(atom.clone(), mul2(&crate::functions::func(Func::Exp, vec![shift]), &crate::ops::pow(&t, &Expr::int(k))));
    }
    // the numerator of the equation in `t`; its denominator may only be a
    // power of `t` (`t != 0`)
    let lift = |e: &Expr| -> Option<Expr> {
        let (pn, pd) = crate::simplify::numer_denom(&crate::simplify::together(&crate::subs::xreplace(e, &map)));
        if has(&pd, x) || has(&pd, &t) && !complete_roots(&crate::expand::expand(&pd), &t)?.iter().all(Expr::is_zero) {
            return None;
        }
        let p = crate::expand::expand(&pn);
        if has(&p, x) { None } else { Some(p) }
    };
    let mut p = lift(num)?;
    if allow_symbolic && !has(den, x) {
        // SymPy's generator is that of the cleared numerator (`_solve_trig1`
        // keeps `exp(I*x)` instead): a polynomial
        // in `t**d` is solved in `exp(d*(g*x + c))`
        let terms = match p.kind() { Kind::Add(..) => p.args(), _ => vec![p.clone()] };
        let mut d = 0u32;
        let mut monomials = Vec::new();
        for term in &terms {
            let k = crate::poly::api::degree(term, Some(&t)).ok()??;
            d = d.gcd(&k);
            monomials.push((div(term, &crate::ops::pow(&t, &Expr::int(i64::from(k)))), k));
        }
        if d > 1 && monomials.iter().all(|(coeff, _)| !has(coeff, &t)) {
            p = crate::ops::add(&monomials.iter().map(|(coeff, k)| mul2(coeff, &crate::ops::pow(&t, &Expr::int(i64::from(k / d))))).collect::<Vec<_>>());
            g = mul2(&Expr::int(i64::from(d)), &g);
            c = mul2(&Expr::int(i64::from(d)), &c);
        }
    }
    let roots = complete_roots(&p, &t)?;
    if roots.is_empty() && crate::poly::api::degree(&p, Some(&t)).ok()?? == 0 {
        if p.is_zero() {
            return None;
        }
        if !has(den, x) {
            return Some(empty_set());
        }
    }
    let mut poles = Vec::new();
    if has(den, x) {
        let q = lift(den)?;
        poles = complete_roots(&q, &t)?;
    }
    let rational_base = atoms.iter().any(|(.., r)| *r);
    let n = context.periodic_variable().clone();
    let two_pi = mul2(&Expr::int(2), &Expr::pi());
    let mut parts = Vec::new();
    let mut seen: Vec<Expr> = Vec::new();
    for k in roots {
        if k.is_zero() || seen.iter().any(|s| crate::expand::expand(&sub(s, &k)).is_zero()) {
            continue;
        }
        seen.push(k.clone());
        if poles.iter().any(|s| crate::expand::expand(&sub(s, &k)).is_zero()) {
            continue;
        }
        let modulus = crate::functions::log(&crate::functions::func(Func::Abs, vec![k.clone()]));
        let modulus = if rational_base { crate::simp::expand_log(&modulus, true) } else { modulus };
        let center = add2(&mul2(&Expr::i(), &add2(&mul2(&two_pi, &n), &crate::functions::func(Func::Arg, vec![k.clone()]))), &modulus);
        let mut expr = div(&sub(&center, &c), &g);
        if rational_base {
            expr = crate::expand::expand(&expr);
        }
        if guarded {
            // the composed image `imageset(Lambda(n, n/(I*mu)), ...)` keeps
            // its form, and is intersected with a real domain as a whole
            let lattice = image_set(&lambda(&n, &expr), &integers());
            parts.push(if is_op(domain, SetOp::Complexes) { lattice } else { intersection(vec![lattice, domain.clone()]) });
            continue;
        }
        let expr = canonical_integer_lattice(&expr, &n);
        parts.push(lattice_in_domain(&expr, &n, domain, context));
    }
    if guarded {
        let m = div(&g, &Expr::i());
        let (mn, md) = crate::simplify::numer_denom(&m);
        let strip = |e: &Expr| crate::simptools::primitive(e).1;
        let (mn, md) = (strip(&mn), strip(&md));
        let cond = crate::logic::and(vec![crate::relational::ne(&md, &Expr::zero()),
            crate::relational::ne(&mn, &Expr::zero())]);
        return Some(condition_set(x, &cond, &union(parts)));
    }
    Some(union(parts))
}

/// The right-hand side of the reference `_invert_complex`: a finite set,
/// an image of `Integers` (one integer variable) or of `Integers**k`, the
/// latter possibly under a held affine map `Lambda(x, ...)`.
enum Inverse {
    Finite(Vec<Expr>),
    Lattice(Vec<Expr>, Expr),
    Held(Vec<Expr>, Expr, Expr),
}

/// `imageset(Lambda(v, phi(v)), g_ys)` for an affine `phi`.
fn map_inverse(g: Inverse, phi: &dyn Fn(&Expr) -> Expr, variable: &Expr) -> Inverse {
    match g {
        Inverse::Finite(values) => Inverse::Finite(values.iter().map(phi).collect()),
        Inverse::Lattice(vars, body) if vars.len() == 1 => {
            let composed = phi(&body);
            Inverse::Lattice(vars.clone(), canonical_integer_lattice(&composed, &vars[0]))
        }
        Inverse::Lattice(vars, body) => Inverse::Held(vars, body, phi(variable)),
        Inverse::Held(vars, body, outer) => Inverse::Held(vars, body, phi(&outer)),
    }
}

/// An image domain with every logarithmic preimage constrained to be finite.
/// The tuple binder is the same tuple consumed by the ImageSet Lambda.
fn inverse_image(inverse: Inverse, guards: &[Expr], variable: &Expr) -> Option<Expr> {
    let image = |vars: &[Expr], body: &Expr| -> Option<Expr> {
        let (signature, binder, base) = if vars.len() == 1 {
            (vars[0].clone(), vars[0].clone(), integers())
        } else {
            let tuple = crate::special::tuple(vars.to_vec());
            (crate::special::tuple(vec![tuple.clone()]), tuple, product_set(vec![integers(); vars.len()]))
        };
        let base = condition_set(&binder, &crate::logic::and(guards.to_vec()), &base);
        let lam = lambda_checked(&signature, body).ok()?;
        Some(image_set(&lam, &base))
    };
    match inverse {
        Inverse::Finite(values) => Some(finite_set(values)),
        Inverse::Lattice(vars, body) => image(&vars, &body),
        Inverse::Held(vars, body, outer) => {
            let inner = image(&vars, &body)?;
            Some(image_set(&lambda(variable, &outer), &inner))
        }
    }
}

/// Whether a supported affine/exponential chain takes real values on the
/// real line. An unknown complex coefficient cannot be assumed real: that
/// case uses the complete complex inverse intersected with the real domain.
fn real_exponential_chain(value: &Expr, x: &Expr) -> bool {
    if value == x { return true; }
    if !has(value, x) {
        return crate::assumptions::exact_queries(|| crate::assumptions::is_real(value)) == Some(true);
    }
    match value.kind() {
        Kind::Add(..) | Kind::Mul(..) => value.args().iter().all(|arg| real_exponential_chain(arg, x)),
        Kind::Func(Func::Exp, args) if args.len() == 1 => real_exponential_chain(&args[0], x),
        _ => false,
    }
}

/// Invert nested exponentials with affine arguments. Each complex logarithm
/// adds an integer variable; its argument must be nonzero, including at every
/// previously introduced integer index. Real chains use strictly positive
/// logarithm arguments. Unknown multipliers are split into their nonzero and
/// zero branches rather than silently losing parameter specializations.
fn nested_exponential(f: &Expr, x: &Expr, domain: &Expr) -> Option<Expr> {
    if !any_node(f, &|e| matches!(e.kind(), Kind::Func(Func::Exp, a) if has_exponential_in(&a[0], x))) {
        return None;
    }
    let dependent = if f.is_add() {
        crate::ops::add(&f.args().into_iter().filter(|arg| has(arg, x)).collect::<Vec<_>>())
    } else { f.clone() };
    let real_inverse = is_real_domain(domain) && real_exponential_chain(&dependent, x);
    let ordinary = Expr::symbol("x");
    let variable = if x != &ordinary && has(f, &ordinary) {
        Expr::dummy("x", crate::symbol::Assumptions::NONE)
    } else { ordinary };
    let mut g = Inverse::Finite(vec![Expr::zero()]);
    let mut h = f.clone();
    let two_pi = mul2(&Expr::int(2), &Expr::pi());
    let mut parameters = Vec::new();
    let mut index_guards = Vec::new();
    let mut degenerate = Vec::new();
    let mut domain_guards = Vec::new();
    let mut finite_coefficients = Vec::new();
    let mut guard_scalar = |scalar: &Expr| {
        finite_coefficients.push(contains_expr(scalar, &complexes()));
        // exp(a) is nonzero only for finite a. Guard its parameters too,
        // otherwise a=-oo makes that certified multiplier collapse to zero.
        finite_coefficients.extend(crate::traverse::free_symbols(scalar).iter()
            .map(|parameter| contains_expr(parameter, &complexes())));
    };
    let mut chain = f.clone();
    let mut inverted_exponential = false;
    // Collect the whole chain before an empty finite preimage can stop the
    // inversion: exp(exp(a*x))=1 has no finite-real-a roots, but a=oo still
    // needs its original conditional equation.
    while chain != *x {
        match chain.kind() {
            Kind::Add(..) | Kind::Mul(..) => {
                let (indep, dep): (Vec<Expr>, Vec<Expr>) = chain.args().into_iter().partition(|t| !has(t, x));
                if indep.is_empty() || dep.is_empty() { return None; }
                if chain.is_add() {
                    if inverted_exponential { guard_scalar(&crate::ops::add(&indep)); }
                    chain = crate::ops::add(&dep);
                } else {
                    guard_scalar(&crate::ops::mul(&indep));
                    chain = crate::ops::mul(&dep);
                }
            }
            Kind::Func(Func::Exp, args) if args.len() == 1 => {
                inverted_exponential = true;
                chain = args[0].clone();
            }
            _ => return None,
        }
    }
    crate::traverse::preorder(f, &mut |term| {
        let (_, denominator) = crate::simplify::numer_denom(term);
        if !denominator.is_one() && !has(&denominator, x) {
            domain_guards.push(crate::relational::ne(&denominator, &Expr::zero()));
        }
    });
    // Inversion over C assumes a finite bound value too. Domains containing
    // extended scalars must retain their original-equation complement rather
    // than silently treating x=oo as a finite exponential argument.
    if crate::assumptions::exact_queries(|| is_subset(domain,&complexes())) != Some(true) {
        finite_coefficients.push(contains_expr(x,&complexes()));
    }
    let finish = |branches: Vec<Expr>, finite_coefficients: &[Expr]| {
        let finite = crate::logic::and(finite_coefficients.to_vec());
        let original_equation = Expr::from_unevaluated_kind(Kind::Rel(RelOp::Eq, f.clone(), Expr::zero()));
        // A finite requested domain is exhaustively checked against the
        // original equation. Intersecting it with an indexed family can leave
        // even an exact non-root undecided, and needlessly loses extended
        // bound values. The complementary guards preserve every parameter case.
        let explicit = if finite_elements(domain).is_some() {
            condition_set(x,&crate::logic::and(vec![finite.clone(),original_equation.clone()]),domain)
        } else { condition_set(x, &finite, &union(branches)) };
        // Affine inversion uses subtraction/division in C. Extended scalar
        // specializations can instead collapse exponentials on a half-line;
        // retain the original equation there rather than dividing by infinity
        // or incorrectly declaring the unsupported case empty.
        // Do not let finite-parameter assumptions simplify exp(...) != 0
        // while this branch explicitly allows nonfinite parameters.
        let unsupported = condition_set(x, &crate::logic::and(vec![crate::logic::not(&finite), original_equation]), domain);
        condition_set(x, &crate::logic::and(domain_guards.clone()), &union(vec![explicit, unsupported]))
    };
    // Every successful step removes an outer operation from h. There is no
    // artificial depth cap and no finite prefix of the inverse families.
    while h != *x {
        match h.kind() {
            Kind::Add(..) => {
                let (indep, dep): (Vec<Expr>, Vec<Expr>) = h.args().into_iter().partition(|t| !has(t, x));
                if indep.is_empty() || dep.is_empty() { return None; }
                let c = crate::ops::add(&indep);
                g = map_inverse(g, &|v| sub(v, &c), &variable);
                h = crate::ops::add(&dep);
            }
            Kind::Mul(..) => {
                let (indep, dep): (Vec<Expr>, Vec<Expr>) = h.args().into_iter().partition(|t| !has(t, x));
                if indep.is_empty() || dep.is_empty() { return None; }
                let c = crate::ops::mul(&indep);
                if any_node(&c, &|e| matches!(e.kind(), Kind::Const(Constant::Infinity | Constant::NegInfinity | Constant::ComplexInfinity | Constant::NaN))) {
                    return Some(finish(degenerate, &finite_coefficients));
                }
                let (_, denominator) = crate::simplify::numer_denom(&c);
                parameters.push(crate::relational::ne(&denominator, &Expr::zero()));
                match crate::assumptions::exact_queries(|| crate::assumptions::is_zero(&c)) {
                    Some(false) => {},
                    known => {
                        let constant = subs1(f, &h, &Expr::zero());
                        let condition = crate::logic::and(vec![crate::relational::eq(&c, &Expr::zero()),
                            crate::relational::eq(&constant, &Expr::zero())]);
                        degenerate.push(condition_set(x, &condition, domain));
                        if known == Some(true) { return Some(finish(degenerate, &finite_coefficients)); }
                        parameters.push(crate::relational::ne(&c, &Expr::zero()));
                    }
                }
                g = map_inverse(g, &|v| div(v, &c), &variable);
                h = crate::ops::mul(&dep);
            }
            Kind::Func(Func::Exp, a) if a.len() == 1 => {
                let u = a[0].clone();
                // Affine maps of a multiple-index image are held for the
                // final presentation, but compose them before another log.
                if let Inverse::Held(vars, body, outer) = g {
                    g = Inverse::Lattice(vars, subs1(&outer, &variable, &body));
                }
                g = match g {
                    Inverse::Finite(values) => {
                        let [value] = values.as_slice() else { return None; };
                        let condition = if real_inverse {
                            contains_expr(value, &interval(&Expr::zero(), &Expr::oo(), true, true))
                        } else { crate::logic::and(vec![crate::relational::ne(value, &Expr::zero()),
                            contains_expr(value, &complexes())]) };
                        if matches!(condition.kind(), Kind::Bool(false)) { return Some(finish(degenerate, &finite_coefficients)); }
                        parameters.push(condition);
                        if real_inverse {
                            Inverse::Finite(vec![crate::functions::log(value)])
                        } else {
                            let n = Expr::dummy("n", crate::symbol::Assumptions::INTEGER);
                            let body = add2(&mul2(&Expr::i(), &add2(&mul2(&two_pi, &n), &crate::functions::func(Func::Arg, vec![value.clone()]))),
                                &crate::functions::log(&crate::functions::func(Func::Abs, vec![value.clone()])));
                            Inverse::Lattice(vec![n.clone()], canonical_integer_lattice(&body, &n))
                        }
                    }
                    Inverse::Lattice(vars, body) => {
                        let nonzero_lattice = proved_nonzero(&body) || (vars.len() == 1 &&
                            linear(&body, &vars[0]).is_some_and(|(step, offset)| proved_nonzero(&step) &&
                                crate::assumptions::exact_queries(|| crate::assumptions::is_integer(&div(&neg(&offset), &step))) == Some(false)));
                        if !nonzero_lattice { index_guards.push(crate::relational::ne(&body, &Expr::zero())); }
                        let k = Expr::dummy(&format!("k{}", vars.len()), crate::symbol::Assumptions::INTEGER);
                        let body = add2(&mul2(&Expr::i(), &add2(&mul2(&two_pi, &k), &crate::functions::func(Func::Arg, vec![body.clone()]))),
                            &crate::functions::log(&crate::functions::func(Func::Abs, vec![body])));
                        let mut all = vec![k];
                        all.extend(vars);
                        Inverse::Lattice(all, body)
                    }
                    Inverse::Held(..) => unreachable!("affine image composed before inversion"),
                };
                h = u;
            }
            _ => return None,
        }
    }
    let mut result = inverse_image(g, &index_guards, &variable)?;
    if !is_op(domain, SetOp::Complexes) {
        result = intersection(vec![result, domain.clone()]);
    }
    result = condition_set(x, &crate::logic::and(parameters), &result);
    degenerate.push(result);
    Some(finish(degenerate, &finite_coefficients))
}

/// `_solve_trig1`: an equation in several trigonometric or hyperbolic
/// functions of arguments linear in `x`, rewritten in exponentials.
fn solve_trig_exponential(num: &Expr, den: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    let rewritten = crate::simptools::trigsimp::rewrite_exp(&div(num, den));
    if has_periodic_in(&rewritten, x) {
        return None;
    }
    let (n2, d2) = crate::simplify::numer_denom(&crate::simplify::together(&rewritten));
    let n2 = crate::expand::expand(&n2);
    let d2 = crate::expand::expand(&d2);
    let hyperbolic = any_node(num, &|e| matches!(e.kind(), Kind::Func(Func::Sinh | Func::Cosh | Func::Tanh | Func::Coth, _)));
    solve_exponential_poly(&n2, &d2, x, domain, context, false, !hyperbolic)
}

/// Solutions of `f == 0` restricted to `domain`, as a set.
fn solve_zero(f: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    let real = is_real_domain(domain);
    // product: union of the factors' solutions
    let (num, den) = crate::simplify::numer_denom(&crate::simplify::together(f));
    if !has(&num, x) || super::solveset_steps::vanishes_identically(&num) {
        // `_solve_as_rational`: a zero numerator is solved by the domain
        // minus the zeros of the denominator
        if num.is_zero() || has(&num, x) {
            return Some(exclude_poles(domain.clone(), &den, x, domain, context));
        }
        return Some(empty_set());
    }
    // Preserve independent scalar denominators while inverting a chain.
    // Clearing (exp(exp(x))-y)/a first would lose its a=oo specialization.
    if let Some(result) = nested_exponential(f, x, domain).or_else(|| nested_exponential(&num, x, domain)) {
        let result = if !has(&den, x) {
            condition_set(x, &crate::relational::ne(&den, &Expr::zero()), &result)
        } else { result };
        return Some(exclude_poles(result, &den, x, domain, context));
    }
    if let Kind::Mul(_, facs) = num.kind() {
        if facs.iter().any(|(b, _)| (matches!(b.kind(), Kind::Func(..)) || !real && is_transcendental_in(b, x)) && has(b, x)) && facs.len() > 1 {
            let mut parts = Vec::new();
            for (b, e) in facs {
                if !has(b, x) || e.as_number().is_some_and(|n| n.is_negative()) {
                    continue;
                }
                parts.push(solve_zero(b, x, domain, context)?);
            }
            let u = union(parts);
            return Some(exclude_poles(u, &den, x, domain, context));
        }
    }
    let complex_exponential = !real && has_exponential_in(&num, x);
    if let Some(g) = single_generator(&num, x).filter(|g| !(complex_exponential && matches!(g.kind(), Kind::Func(Func::Exp, _)))) {
        // polynomial in g = trig(a*x + b) (or exp in the complex domain)
        let t = Expr::dummy("solveset_t", crate::symbol::Assumptions::NONE);
        let mut map = rustc_hash::FxHashMap::default();
        map.insert(g.clone(), t.clone());
        let p = crate::subs::xreplace(&num, &map);
        if !has(&p, x) {
            let ks = crate::solve::without_radical_branch(|| crate::solve::solve(&p, &t))?;
            let mut occurrences = 0;
            crate::traverse::preorder(&num, &mut |e| if *e == g { occurrences += 1; });
            context.single.set(occurrences == 1 && crate::integrate::poly_coeffs(&p, &t).is_some_and(|c| c.len() == 2));
            let mut parts = Vec::new();
            for k in ks {
                parts.push(periodic(&g, &k, x, domain, context)?);
            }
            let r = union(parts);
            let r = if is_op(domain, SetOp::Reals) || is_op(domain, SetOp::Complexes) {
                r
            } else {
                intersection(vec![r, domain.clone()])
            };
            return Some(exclude_poles(r, &den, x, domain, context));
        }
    }
    if has_periodic_in(&num, x) {
        // several trigonometric or hyperbolic functions: `_solve_trig1`
        return solve_trig_exponential(&num, &den, x, domain, context);
    }
    if complex_exponential {
        // a polynomial in one exponential generator, inverted over the
        // complex logarithm branches; otherwise unsolved, never the finite
        // principal-branch roots of `solve`
        // `_solve_as_rational`: the solutions of the numerator minus those
        // of the denominator
        let r = match solve_exponential_poly(&num, &Expr::one(), x, domain, context, true, false) {
            Some(r) => r,
            None if !has(&den, x) => return nested_exponential(&num, x, domain),
            None => return None,
        };
        return Some(exclude_poles(r, &den, x, domain, context));
    }
    // absolute values over the reals: split on the sign of the argument
    if real {
        if let Some(r) = solve_abs(&num, x, domain, context) {
            return Some(exclude_poles(r, &den, x, domain, context));
        }
    }
    // solve() intentionally returns finite Lambert representatives. These are
    // not a complete complex solution set, even if W simplifies away (W(E)=1).
    // Symbolic parameter cases also need domain/branch stratification before
    // a finite real set can be certified. Preserve ConditionSet in these cases.
    if (!real || crate::traverse::free_symbols(&num).iter().any(|s| s != x))
        && crate::solve::lambert::solve(&num, x).is_some() {
        return None;
    }
    // the radical branches of `_solve` are not `_solveset`'s: its radical
    // equations go to `_solve_radical` (`radical_step`)
    let sols = crate::solve::without_radical_branch(|| crate::solve::solve(f, x))?;
    let found = filter_domain(sols, domain);
    // `solve` drops a candidate only when it is provably a pole; the
    // reference `_solveset` subtracts the solutions of the denominator, so
    // an undecided pole (`(x - a)/(x - b)` has no zero when `a = b`) stays
    // as a Complement.
    if has(&den, x)
        && !(free_symbols_except(&den, x).is_empty() && free_symbols_except(&found, x).is_empty())
        && !is_empty_set(&found)
    {
        if let Some(poles) = solve_zero(&den, x, domain, context) {
            return Some(complement(&found, &poles));
        }
    }
    Some(found)
}

/// An exact algebraic number (rationals, `I`, sums, products and rational
/// powers): the radical roots of `_solve_as_poly`. Roots with logarithms
/// (`log(4)/log(-2)` of `(-2)**x = 4`) are left undecided: the real
/// inverter they need is not ported, and dropping the complex candidate
/// would turn the missing real root into a false EmptySet.
fn is_exact_algebraic(e: &Expr) -> bool {
    use crate::number::Number;
    let exact = |n: &Number| !matches!(n, Number::Float(_));
    match e.kind() {
        Kind::Num(n) => exact(n),
        Kind::Const(Constant::I) => true,
        Kind::Add(c, terms) => exact(c) && terms.iter().all(|(t, k)| exact(k) && is_exact_algebraic(t)),
        Kind::Mul(c, factors) => exact(c) && factors.iter().all(|(b, x)| {
            is_exact_algebraic(b) && matches!(x.kind(), Kind::Num(n) if exact(n))
        }),
        Kind::Pow(b, x) => is_exact_algebraic(b) && matches!(x.kind(), Kind::Num(n) if exact(n)),
        _ => false,
    }
}

fn free_symbols_except(e: &Expr, x: &Expr) -> Vec<Expr> {
    crate::traverse::free_symbols(e).into_iter().filter(|s| s != x).collect()
}

fn is_empty_set(e: &Expr) -> bool {
    is_op(e, SetOp::EmptySet)
}

fn exclude_poles(s: Expr, den: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Expr {
    if !has(den, x) {
        return s;
    }
    match finite_elements(&s) {
        Some(el) => {
            let keep: Vec<Expr> = el
                .iter()
                .filter(|v| {
                    let d = crate::expand::expand(&subs1(den, x, v));
                    !d.is_zero()
                })
                .cloned()
                .collect();
            let _ = domain;
            finite_set(keep)
        }
        None => match solve_zero(den, x, domain, context) {
            Some(poles) => complement(&s, &poles),
            None => mk(SetOp::Complement, vec![s, condition_set(x, &crate::relational::eq(den, &Expr::zero()), domain)]),
        },
    }
}

fn find_abs(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    if let Kind::Func(Func::Abs, a) = e.kind() {
        if has(&a[0], x) && !out.contains(e) {
            out.push(e.clone());
        }
    }
    for a in e.args() {
        find_abs(&a, x, out);
    }
}

fn solve_abs(f: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    let mut abs = Vec::new();
    find_abs(f, x, &mut abs);
    if abs.len() != 1 {
        return None;
    }
    let Kind::Func(_, a) = abs[0].kind() else { return None };
    let u = &a[0];
    let mut out = Vec::new();
    for sign in [1i64, -1] {
        let repl = if sign == 1 { u.clone() } else { neg(u) };
        let g = subs1(f, &abs[0], &repl);
        let s = solve_zero(&g, x, domain, context)?;
        let el = finite_elements(&s).map(|e| e.to_vec()).unwrap_or_default();
        if !is_empty_set(&s) && finite_elements(&s).is_none() {
            return None;
        }
        for v in el {
            let uv = subs1(u, x, &v);
            let ok = if sign == 1 { cmp3(&uv, &Expr::zero()).map(|o| o != Ordering::Less) } else { lt3(&uv, &Expr::zero()) };
            match ok {
                Some(true) => out.push(v),
                Some(false) => {}
                None => return None,
            }
        }
    }
    Some(intersection(vec![finite_set(out), domain.clone()]))
}

fn filter_domain(sols: Vec<Expr>, domain: &Expr) -> Expr {
    if is_op(domain, SetOp::Complexes) || is_op(domain, SetOp::UniversalSet) {
        return finite_set(sols);
    }
    if is_real_domain(domain) {
        let mut real = Vec::new();
        let mut unknown = Vec::new();
        for s in sols {
            match crate::assumptions::is_real(&s) {
                Some(true) => real.push(s),
                Some(false) => {}
                None => {
                    // The reference `_solve_as_poly` writes the numeric
                    // roots as `expand_complex(s)` before intersecting with
                    // the domain, which decides Cardano roots with a
                    // cube root of unity factor.
                    let expanded = if is_exact_algebraic(&s) {
                        crate::complex_parts::expand_complex(&s, &crate::expand::ExpandOptions::default()).ok()
                    } else {
                        None
                    };
                    let realness = expanded.as_ref().and_then(|e| {
                        crate::assumptions::is_real(e).or_else(|| {
                            // a nonzero imaginary part (decided by the
                            // assumption queries) proves the root not real
                            match crate::complex_parts::as_real_imag_bounded(e, true, None) {
                                Ok(Some((_, im))) if crate::assumptions::is_nonzero(&im) == Some(true) => Some(false),
                                _ => None,
                            }
                        })
                    });
                    match realness {
                        Some(true) => real.push(expanded.unwrap()),
                        Some(false) => {}
                        None => unknown.push(s),
                    }
                }
            }
        }
        let mut parts = vec![intersection(vec![finite_set(real), domain.clone()])];
        if !unknown.is_empty() {
            parts.push(mk(SetOp::Intersection, vec![finite_set(unknown), domain.clone()]));
        }
        return union(parts);
    }
    intersection(vec![finite_set(sols), domain.clone()])
}

/// SymPy's `solveset(f, x, domain)`; `f` may be an expression (`f == 0`),
/// an `Eq`/`Ne` or an inequality (real domains only). Returns an `Err`
/// for unsupported inputs.
pub fn solveset(f: &Expr, x: &Expr, domain: &Expr) -> Result<Expr, String> {
    let result = solveset_with_context(f,x,domain,&SolveContext::default())?;
    if !x.is_symbol() || matches!(f.kind(), Kind::Rel(op, ..) if *op != RelOp::Eq) {
        return Ok(result);
    }
    Ok(super::solveset_steps::check_solutions(f, x, result))
}

fn solveset_with_context(f: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Result<Expr, String> {
    match f.kind() {
        Kind::Bool(true) => return Ok(domain.clone()),
        Kind::Bool(false) => return Ok(empty_set()),
        Kind::Rel(RelOp::Eq, a, b) => return solveset_with_context(&sub(a, b), x, domain, context),
        Kind::Rel(RelOp::Ne, a, b) => {
            let s = solveset_with_context(&sub(a, b), x, domain, context)?;
            return Ok(complement(domain, &s));
        }
        Kind::Rel(..) => {
            if !is_real_domain(domain) {
                return Ok(condition_set(x, f, domain));
            }
            return match solve_univariate_inequality_set(f, x, domain) {
                Ok(s) => Ok(s),
                Err(_) => Ok(condition_set(x, f, domain)),
            };
        }
        _ => {}
    }
    if !x.is_symbol() {
        return Err(format!("{x} is not a symbol"));
    }
    if let Ok(Some(result)) = super::solveset_modular(f, x, domain) {
        return Ok(result);
    }
    if !has(f, x) {
        return Ok(match crate::assumptions::exact_queries(|| crate::assumptions::is_zero(f)) {
            Some(true) => domain.clone(),
            Some(false) => empty_set(),
            None => condition_set(x, &crate::relational::eq(f, &Expr::zero()), domain),
        });
    }
    if super::solveset_steps::vanishes_identically(f) {
        return Ok(domain.clone());
    }
    match solve_zero(f, x, domain, context) {
        Some(s) => Ok(s),
        None => Ok(condition_set(x, &crate::relational::eq(f, &Expr::zero()), domain)),
    }
}

/// `solveset(f, x, S.Reals)`.
pub fn solveset_real(f: &Expr, x: &Expr) -> Result<Expr, String> {
    solveset(f, x, &reals())
}

// ----------------------------------------------------------------------
// inequalities
// ----------------------------------------------------------------------

fn holds(op: RelOp, sign: Ordering) -> bool {
    match op {
        RelOp::Lt => sign == Ordering::Less,
        RelOp::Le => sign != Ordering::Greater,
        RelOp::Gt => sign == Ordering::Greater,
        RelOp::Ge => sign != Ordering::Less,
        RelOp::Eq => sign == Ordering::Equal,
        RelOp::Ne => sign != Ordering::Equal,
    }
}

/// Real domain of definition of `f` (log arguments positive, even roots
/// non-negative), as a set; `None` if it cannot be determined.
fn real_domain(f: &Expr, x: &Expr) -> Option<Expr> {
    let mut conds: Vec<Expr> = Vec::new();
    fn rec(e: &Expr, x: &Expr, conds: &mut Vec<Expr>) -> Option<()> {
        if !has(e, x) {
            return Some(());
        }
        match e.kind() {
            Kind::Func(Func::Log, a) => conds.push(crate::relational::gt(&a[0], &Expr::zero())),
            Kind::Pow(b, ex) => {
                if let Some(n) = ex.as_number() {
                    if !n.is_integer() {
                        let d = n.denom();
                        if d.to_i64().is_some_and(|d| d % 2 == 0) {
                            conds.push(crate::relational::ge(b, &Expr::zero()));
                        }
                    }
                } else if has(ex, x) {
                    return None;
                }
            }
            Kind::Func(f, _) => {
                if !matches!(f, Func::Exp | Func::Abs | Func::Log | Func::Sign) {
                    return None;
                }
            }
            _ => {}
        }
        for a in e.args() {
            rec(&a, x, conds)?;
        }
        Some(())
    }
    rec(f, x, &mut conds)?;
    let mut d = reals();
    for c in conds {
        let s = sign_analysis(&c, x, &reals(), false)?;
        d = intersection(vec![d, s]);
    }
    Some(d)
}

/// Real roots of `g` (a finite set of numbers), via `solveset`.
fn real_roots(g: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    if !has(g, x) {
        return Some(vec![]);
    }
    let s = solve_zero(g, x, &reals(), &SolveContext::default())?;
    if is_empty_set(&s) {
        return Some(vec![]);
    }
    let el = finite_elements(&s)?;
    if el.iter().any(|v| !crate::logic::is_number_like(v)) {
        return None;
    }
    Some(el.to_vec())
}

fn sign_at(f: &Expr, x: &Expr, t: &Expr) -> Option<Ordering> {
    let v = subs1(f, x, t);
    if let Kind::Num(n) = v.kind() {
        return n.cmp_value(&crate::number::Number::ZERO);
    }
    let z = crate::evalf::eval_c64(&v)?;
    if !z.re.is_finite() || z.im.abs() > 1e-9 * (1.0 + z.re.abs()) {
        return None;
    }
    if z.re.abs() < 1e-12 {
        // decide exactly when possible
        return if crate::expand::expand(&v).is_zero() { Some(Ordering::Equal) } else { z.re.partial_cmp(&0.0) };
    }
    z.re.partial_cmp(&0.0)
}

/// Sign analysis of `lhs op rhs` over the reals.
fn sign_analysis(ineq: &Expr, x: &Expr, domain: &Expr, need_domain: bool) -> Option<Expr> {
    let Kind::Rel(op, l, r) = ineq.kind() else {
        return match ineq.kind() {
            Kind::Bool(true) => Some(domain.clone()),
            Kind::Bool(false) => Some(empty_set()),
            _ => None,
        };
    };
    let op = *op;
    let f = sub(l, r);
    if !has(&f, x) {
        return match crate::relational::rel(op, &f, &Expr::zero()).kind() {
            Kind::Bool(true) => Some(domain.clone()),
            Kind::Bool(false) => Some(empty_set()),
            _ => None,
        };
    }
    // A known-sign real coefficient gives an exact half-line even when the
    // root contains parameters, which numeric sign-chart sampling cannot do.
    if let Some((a,b))=linear(&f,x) {
        use crate::assumptions::{is_real,is_positive,is_negative};
        if is_real(&a)==Some(true) && is_real(&b)==Some(true) {
            let direction=if is_positive(&a)==Some(true) {Some(op)}
                else if is_negative(&a)==Some(true) {Some(crate::logic::reversed_op(op))}
                else {None};
            if let Some(direction)=direction {
                let root=div(&neg(&b),&a);
                let solution=match direction {
                    RelOp::Lt=>interval(&Expr::neg_oo(),&root,true,true),
                    RelOp::Le=>interval(&Expr::neg_oo(),&root,true,false),
                    RelOp::Gt=>interval(&root,&Expr::oo(),true,true),
                    RelOp::Ge=>interval(&root,&Expr::oo(),false,true),
                    RelOp::Eq=>finite_set(vec![root.clone()]),
                    RelOp::Ne=>complement(&reals(),&finite_set(vec![root.clone()])),
                };
                return Some(intersection(vec![solution,domain.clone()]));
            }
        }
    }
    if has_periodic(&f) {
        return None;
    }
    let ft = crate::simplify::together(&f);
    let (num, den) = crate::simplify::numer_denom(&ft);
    let dom = if need_domain { real_domain(&f, x)? } else { reals() };
    let mut crit: Vec<Expr> = Vec::new();
    let mut poles: Vec<Expr> = Vec::new();
    // The roots of the numerator are exact zeros. Their radical forms (a
    // Cardano root, for one) need not expand to zero after substitution,
    // and a numeric sign at a zero is rounding noise.
    let zeros = real_roots(&num, x)?;
    crit.extend(zeros.iter().cloned());
    for p in real_roots(&den, x)? {
        poles.push(p.clone());
        crit.push(p);
    }
    // boundaries of the domain of definition
    collect_endpoints(&dom, &mut crit);
    let mut pts: Vec<Expr> = Vec::new();
    for c in crit {
        if !pts.iter().any(|p| cmp3(p, &c) == Some(Ordering::Equal)) {
            pts.push(c);
        }
    }
    let mut bad = false;
    pts.sort_by(|a, b| {
        cmp3(a, b).unwrap_or_else(|| {
            bad = true;
            Ordering::Equal
        })
    });
    if bad {
        return None;
    }
    let in_dom = |t: &Expr| -> bool { is_op(&dom, SetOp::Reals) || contains3(&dom, t) == Some(true) };
    // segments: gap_0, pt_0, gap_1, ..., pt_{k-1}, gap_k
    enum Seg {
        Gap(Expr, Expr),
        Pt(Expr),
    }
    let mut segs: Vec<(Seg, bool)> = Vec::new();
    let k = pts.len();
    for i in 0..=k {
        let lo = if i == 0 { Expr::neg_oo() } else { pts[i - 1].clone() };
        let hi = if i == k { Expr::oo() } else { pts[i].clone() };
        let t = match (i == 0, i == k) {
            (true, true) => Expr::zero(),
            (true, false) => sub(&hi, &Expr::one()),
            (false, true) => add2(&lo, &Expr::one()),
            _ => div(&add2(&lo, &hi), &Expr::int(2)),
        };
        let sat = if in_dom(&t) {holds(op,sign_at(&f,x,&t)?)} else {false};
        segs.push((Seg::Gap(lo, hi), sat));
        if i < k {
            let p = &pts[i];
            let sat = if poles.iter().any(|q| cmp3(q,p)==Some(Ordering::Equal)) || !in_dom(p) {
                false
            } else if zeros.iter().any(|z| cmp3(z,p)==Some(Ordering::Equal)) {
                holds(op,Ordering::Equal)
            } else {
                holds(op,sign_at(&f,x,p)?)
            };
            segs.push((Seg::Pt(p.clone()), sat));
        }
    }
    let mut pieces: Vec<Expr> = Vec::new();
    let mut start: Option<(Expr, bool)> = None;
    let mut last: Option<&Seg> = None;
    let flush = |start: (Expr, bool), last: &Seg, pieces: &mut Vec<Expr>| match last {
        Seg::Gap(_, hi) => pieces.push(interval(&start.0, hi, start.1, true)),
        Seg::Pt(p) => {
            if !start.1 && start.0 == *p {
                pieces.push(finite_set(vec![p.clone()]));
            } else {
                pieces.push(interval(&start.0, p, start.1, false));
            }
        }
    };
    for (seg, sat) in &segs {
        if *sat {
            if start.is_none() {
                start = Some(match seg {
                    Seg::Gap(lo, _) => (lo.clone(), true),
                    Seg::Pt(p) => (p.clone(), false),
                });
            }
            last = Some(seg);
        } else if let Some(s) = start.take() {
            flush(s, last.unwrap(), &mut pieces);
        }
    }
    if let Some(s) = start.take() {
        flush(s, last.unwrap(), &mut pieces);
    }
    let set = union(pieces);
    Some(if is_op(domain, SetOp::Reals) { set } else { intersection(vec![set, domain.clone()]) })
}

fn collect_endpoints(s: &Expr, out: &mut Vec<Expr>) {
    match as_set(s) {
        Some((SetOp::Interval, a)) => {
            for p in &a[..2] {
                if !is_infinite(p) {
                    out.push(p.clone());
                }
            }
        }
        Some((SetOp::Union, a)) => a.iter().for_each(|p| collect_endpoints(p, out)),
        Some((SetOp::FiniteSet, a)) => out.extend(a.iter().cloned()),
        _ => {}
    }
}

fn has_periodic(e: &Expr) -> bool {
    match e.kind() {
        Kind::Func(f, _) if matches!(f, Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc) => true,
        _ => e.args().iter().any(has_periodic),
    }
}

/// Where a periodic inequality is being solved, which decides whether the
/// one-period answer of SymPy's `solve_univariate_inequality` is valid.
#[derive(Clone, Debug)]
pub enum PeriodicWindow {
    /// A direct call of `solveset`, `solve_univariate_inequality` or
    /// `reduce_inequalities`: over `S.Reals`, the reference solves on one
    /// period `[0, T)`.
    Direct,
    /// The answer must hold on the whole real line (an antiderivative, a
    /// condition rewrite): only a range decision (always or never true) is
    /// accepted.
    Whole,
    /// A definite integral over `[lo, hi]`: the one-period answer is used
    /// when `[lo, hi]` lies in `[0, T]`.
    Interval(Expr, Expr),
}

/// Host solver of the periodic branch of `solve_univariate_inequality`
/// (`expr`, `gen`, `domain`, window). `None` keeps the native path.
pub type PeriodicHost = fn(&Expr, &Expr, &Expr, &PeriodicWindow) -> Option<Result<Expr, String>>;

static PERIODIC_HOST: std::sync::OnceLock<PeriodicHost> = std::sync::OnceLock::new();

thread_local! {
    static PERIODIC_WINDOW: std::cell::RefCell<PeriodicWindow> = const { std::cell::RefCell::new(PeriodicWindow::Whole) };
}

/// Install the periodic inequality solver (the Python bindings install the
/// port of the reference branch).
pub fn register_periodic_host(host: PeriodicHost) {
    let _ = PERIODIC_HOST.set(host);
}

/// Run `f` with `window` as the periodic inequality context.
pub fn with_periodic_window<T>(window: PeriodicWindow, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<PeriodicWindow>);
    impl Drop for Restore {
        fn drop(&mut self) {
            if let Some(w) = self.0.take() {
                PERIODIC_WINDOW.with(|cell| *cell.borrow_mut() = w);
            }
        }
    }
    let previous = PERIODIC_WINDOW.with(|cell| std::mem::replace(&mut *cell.borrow_mut(), window));
    let _restore = Restore(Some(previous));
    f()
}

fn periodic_inequality(ineq: &Expr, x: &Expr, domain: &Expr) -> Option<Result<Expr, String>> {
    let Kind::Rel(op, l, r) = ineq.kind() else { return None };
    if matches!(op, RelOp::Eq) || !has_periodic(&sub(l, r)) {
        return None;
    }
    let host = PERIODIC_HOST.get()?;
    let window = PERIODIC_WINDOW.with(|cell| cell.borrow().clone());
    host(ineq, x, domain, &window)
}

fn solve_univariate_inequality_set(ineq: &Expr, x: &Expr, domain: &Expr) -> Result<Expr, String> {
    // SymPy keeps a real generator but its root finding (`solvify` over
    // `S.Reals`) ignores further assumptions: a positive `k` still gives
    // `k**2 > 4` as `(-oo, -2) U (2, oo)`. The sign chart filters roots
    // through the assumptions of `x`, so it runs on a plain real Dummy
    // (otherwise the root `-2` is dropped and the chart can be wrong:
    // `k**2 > 1/(4*pi**2)` became `Ne(k, 1/(2*pi))`).
    if let Some(sym) = x.as_sym() {
        if crate::assumptions::is_real(x) == Some(true) {
            let real_x = Expr::dummy("inequality_gen", crate::symbol::Assumptions::REAL);
            let plain = real_x.as_sym().map(|d| d.assumptions());
            if plain != Some(sym.assumptions()) {
                let result = solve_univariate_inequality_checked(&subs1(ineq, x, &real_x), &real_x, domain)?;
                return Ok(subs1(&result, &real_x, x));
            }
        }
    }
    solve_univariate_inequality_checked(ineq, x, domain)
}

fn solve_univariate_inequality_checked(ineq: &Expr, x: &Expr, domain: &Expr) -> Result<Expr, String> {
    let result = solve_univariate_inequality_inner(ineq, x, domain)?;
    if is_empty_set(&result) {
        return Ok(result);
    }
    let points = super::solveset_steps::relational_singularities(ineq, x, &|d| solve_zero(d, x, &reals(), &SolveContext::default()))
        .map_err(|_| format!("The inequality, {ineq}, cannot be solved using solve_univariate_inequality."))?;
    if points.is_empty() {
        return Ok(result);
    }
    // the reference solves over the reals, excludes the singularities and
    // then intersects with a smaller domain
    if is_op(domain, SetOp::Reals) {
        return Ok(complement(&result, &finite_set(points)));
    }
    let real = solve_univariate_inequality_inner(ineq, x, &reals())?;
    Ok(intersection(vec![complement(&real, &finite_set(points)), domain.clone()]))
}

fn solve_univariate_inequality_inner(ineq: &Expr, x: &Expr, domain: &Expr) -> Result<Expr, String> {
    match ineq.kind() {
        Kind::Rel(RelOp::Eq, ..) => return solveset(ineq, x, domain),
        Kind::Bool(true) => return Ok(domain.clone()),
        Kind::Bool(false) => return Ok(empty_set()),
        Kind::Rel(..) => {}
        _ => return Err(format!("{ineq} is not a relational")),
    }
    // Positive odd integer powers are strictly increasing on the reals.
    // Invert them before asking the finite polynomial-root solver, retaining
    // the real root for negative bounds rather than a complex principal root.
    if is_real_domain(domain) {
        let Kind::Rel(op,left,right)=ineq.kind() else {unreachable!()};
        for (power,bound,direction) in [(left,right,*op),(right,left,crate::logic::reversed_op(*op))] {
            let Kind::Pow(base,exponent)=power.kind() else {continue};
            let Some(order)=exponent.as_i64() else {continue};
            if order<=1 || order%2==0 || has(bound,x) {continue}
            let real_x=Expr::dummy("inequality_x",crate::symbol::Assumptions::REAL);
            if crate::assumptions::is_real(bound)!=Some(true)
                || crate::assumptions::is_real(&subs1(base,x,&real_x))!=Some(true) {continue}
            let root=mul2(&crate::functions::sign(bound),
                &crate::ops::pow(&crate::functions::abs(bound),&Expr::rational(1,order)));
            let reduced=crate::relational::rel(direction,base,&root);
            return solve_univariate_inequality_inner(&reduced,x,domain);
        }
    }
    if let Some(result) = periodic_inequality(ineq, x, domain) {
        return result;
    }
    sign_analysis(ineq, x, domain, true)
        .ok_or_else(|| format!("The inequality, {ineq}, cannot be solved using solve_univariate_inequality."))
}

/// SymPy's `solve_univariate_inequality(expr, gen, relational, domain)`.
pub fn solve_univariate_inequality(ineq: &Expr, x: &Expr, relational: bool, domain: &Expr) -> Result<Expr, String> {
    let s = solve_univariate_inequality_set(ineq, x, domain)?;
    Ok(if relational { as_relational(&s, x) } else { s })
}

/// SymPy's `reduce_inequalities` for univariate polynomial / rational
/// inequalities (and equations) over the reals.
pub fn reduce_inequalities(ineqs: &[Expr], x: Option<&Expr>) -> Result<Expr, String> {
    let gen = match x {
        Some(x) => x.clone(),
        None => {
            let mut syms: Vec<Expr> = Vec::new();
            for e in ineqs {
                for s in crate::traverse::free_symbols(e) {
                    if !syms.contains(&s) {
                        syms.push(s);
                    }
                }
            }
            match syms.len() {
                0 => {
                    let all = crate::logic::and(ineqs.to_vec());
                    return Ok(all);
                }
                1 => syms.pop().unwrap(),
                _ => return Err("only univariate inequalities are supported".into()),
            }
        }
    };
    // like SymPy, decide what the assumptions decide for a real variable
    let real_gen = match gen.as_sym() {
        Some(s) if !s.assumptions().has(crate::symbol::Assumptions::REAL) => {
            Some(Expr::symbol_with(&format!("_real_{}", s.name()), crate::symbol::Assumptions::REAL))
        }
        _ => None,
    };
    let mut sets: Vec<Expr> = Vec::new();
    for e in ineqs {
        let e = &match &real_gen {
            Some(r) => match subs1(e, &gen, r).kind() {
                Kind::Bool(b) => Expr::boolean(*b),
                _ => e.clone(),
            },
            None => e.clone(),
        };
        match e.kind() {
            Kind::Bool(true) => continue,
            Kind::Bool(false) => return Ok(Expr::boolean(false)),
            _ => {}
        }
        if !has(e, &gen) {
            return Err(format!("{e} does not contain {gen}"));
        }
        let s = match solveset(e, &gen, &reals()) {
            Ok(s) if !is_op(&s, SetOp::ConditionSet) => s,
            _ => {
                // linear inequality with symbolic coefficients
                if ineqs.len() == 1 {
                    if let Some(r) = isolate_linear(e, &gen) {
                        return Ok(r);
                    }
                }
                return Err(format!("could not reduce {e}"));
            }
        };
        sets.push(s);
    }
    if sets.is_empty() {
        return Ok(Expr::boolean(true));
    }
    let s = intersection(sets);
    // the reference solves for a real dummy, for which the relational form
    // of the whole real line is True
    if is_op(&s, SetOp::Reals) {
        return Ok(Expr::boolean(true));
    }
    Ok(as_relational(&s, &gen))
}

fn isolate_linear(e: &Expr, x: &Expr) -> Option<Expr> {
    let Kind::Rel(op, l, r) = e.kind() else { return None };
    let (a, b) = linear(&sub(l, r), x)?;
    let rhs = div(&neg(&b), &a);
    let op = if crate::assumptions::is_negative(&a)==Some(true) {crate::logic::reversed_op(*op)}
        else if crate::assumptions::is_positive(&a)==Some(true) {*op}
        else {return None};
    Some(crate::relational::rel(op, x, &rhs))
}

#[cfg(test)]
mod nested_exponential_tests {
    use super::*;
    use crate::parser::parse_str;

    fn p(source: &str) -> Expr { parse_str(source).unwrap() }

    #[test]
    fn symbolic_nested_exp_keeps_parameters_and_excludes_zero_preimages() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let result = solveset(&p("exp(exp(x)) - y"), &x, &complexes()).unwrap();
        assert!(is_op(&result, SetOp::ConditionSet), "{result}");
        assert_eq!(crate::traverse::free_symbols(&result), vec![y.clone()]);
        assert_eq!(subs1(&result, &y, &Expr::zero()), empty_set());
        let one = subs1(&result, &y, &Expr::one());
        let Some((SetOp::ImageSet, image)) = as_set(&one) else { panic!("{one}") };
        let Some((SetOp::ConditionSet, guard)) = as_set(&image[1]) else { panic!("{one}") };
        let excluded = crate::special::tuple(vec![Expr::zero(), Expr::zero()]);
        let admitted = crate::special::tuple(vec![Expr::zero(), Expr::one()]);
        assert_eq!(contains(&image[1], &excluded), Expr::boolean(false));
        assert_eq!(contains(&image[1], &admitted), Expr::boolean(true));
        assert!(has(&guard[1], &guard[0].args()[1]));
    }

    #[test]
    fn real_nested_exp_and_zero_affine_coefficients() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let a = Expr::symbol_with("a", crate::symbol::Assumptions::REAL);
        let real = solveset(&p("exp(exp(x)) - y"), &x, &reals()).unwrap();
        for target in [Expr::zero(), Expr::one(), Expr::minus_one(), Expr::i()] {
            assert_eq!(subs1(&real, &y, &target), empty_set(), "{real} at {target}");
        }
        let two = solveset(&p("exp(exp(x)) - 2"), &x, &reals()).unwrap();
        assert_eq!(two, finite_set(vec![crate::functions::log(&crate::functions::log(&Expr::int(2)))]));
        let source = sub(&crate::functions::exp(&crate::functions::exp(&mul2(&a, &x))), &y);
        let result = solveset(&source, &x, &reals()).unwrap();
        let zero = subs1(&result, &a, &Expr::zero());
        assert_eq!(subs1(&zero, &y, &crate::functions::exp(&Expr::one())), reals(), "{result}");
        assert_eq!(subs1(&zero, &y, &Expr::int(2)), empty_set());
    }

    #[test]
    fn nested_affine_images_compose_and_parameter_poles_stay_excluded() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let a = Expr::symbol("a");
        let result = solveset(&p("exp(2*exp(3*exp(x))) - y"), &x, &complexes()).unwrap();
        assert!(!has(&result, &p("exp(2*exp(3*exp(x)))")), "{result}");
        assert_eq!(crate::traverse::free_symbols(&result), vec![y.clone()]);
        for source in ["exp(exp(x/a)) - y", "(exp(exp(x)) - y)/a"] {
            let result = solveset(&p(source), &x, &complexes()).unwrap();
            assert_eq!(subs1(&subs1(&result, &a, &Expr::zero()), &y, &Expr::int(2)), empty_set(), "{result}");
        }
        let constant = solveset(&a, &x, &complexes()).unwrap();
        assert_eq!(subs1(&constant, &a, &Expr::zero()), complexes());
        assert_eq!(subs1(&constant, &a, &Expr::one()), empty_set());
        let z = Expr::symbol("z");
        let parameter = Expr::symbol("x");
        let result = solveset(&p("exp(exp(2*z)) - x"), &z, &complexes()).unwrap();
        assert_eq!(crate::traverse::free_symbols(&result), vec![parameter]);
    }

    #[test]
    fn extended_finite_bound_values_are_checked_against_the_original_equation() {
        let x=Expr::symbol("x");let y=Expr::symbol("y");
        let domain=finite_set(vec![Expr::neg_oo(),Expr::zero()]);
        let result=solveset(&p("exp(exp(x))-y"),&x,&domain).unwrap();
        assert_eq!(subs1(&result,&y,&Expr::one()),finite_set(vec![Expr::neg_oo()]));
        assert_eq!(subs1(&result,&y,&Expr::e()),finite_set(vec![Expr::zero()]));
    }

    #[test]
    fn nonfinite_affine_parameters_retain_the_original_equation() {
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        let y = Expr::symbol("y");
        let source = p("exp(exp(a*x)) - y");
        for domain in [complexes(), reals()] {
            let result = solveset(&source, &x, &domain).unwrap();
            for infinity in [Expr::oo(), Expr::neg_oo()] {
                for target in [Expr::one(), Expr::int(2)] {
                    let specialized = subs1(&subs1(&result, &a, &infinity), &y, &target);
                    let original = subs1(&subs1(&source, &a, &infinity), &y, &target);
                    let equation=crate::relational::eq(&original,&Expr::zero());
                    assert_eq!(contains_expr(&infinity,&complexes()),Expr::boolean(false));
                    let retained=crate::logic::not(&Expr::from_unevaluated_kind(Kind::Func(Func::Set(SetOp::Contains),vec![infinity.clone(),complexes()])));
                    assert_eq!(specialized,condition_set(&x,&crate::logic::and(vec![retained,equation.clone()]),&domain), "{source} at a={infinity}, y={target}");
                    // The retained guard is exactly true, so this is the same
                    // original equation and complete domain as a direct solve.
                    assert_eq!(solveset(&original,&x,&domain).unwrap(),condition_set(&x,&equation,&domain));
                }
            }
        }
        // Even a finite-target contradiction cannot discard the unsupported
        // infinite inner shift: exp(exp(x)-oo)=0 on finite complex x.
        let source = p("exp(exp(x) + a)");
        let result = solveset(&source, &x, &complexes()).unwrap();
        let original = subs1(&source, &a, &Expr::neg_oo());
        assert_eq!(subs1(&result, &a, &Expr::neg_oo()),
            condition_set(&x, &crate::relational::eq(&original, &Expr::zero()), &complexes()), "source={source}, result={result}");
        let source = p("(exp(exp(x)) - y)/a");
        let result = solveset(&source, &x, &complexes()).unwrap();
        assert_eq!(subs1(&subs1(&result, &a, &Expr::oo()), &y, &Expr::int(2)), complexes());
    }
}
