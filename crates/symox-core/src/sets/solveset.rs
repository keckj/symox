//! `solveset` (on top of [`crate::solve::solve`], with `ImageSet`s for
//! periodic solutions) and univariate inequalities by sign analysis.

use std::cmp::Ordering;

use super::*;
use crate::expr::RelOp;
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
    if let Ok(Some((real,imaginary)))=crate::complex_parts::as_real_imag(value,false,None){
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
struct SolveContext { periodic_variable: std::sync::OnceLock<Expr> }
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

/// Periodic solutions of `trig(a*x + b) = k` / `exp(a*x + b) = k`.
fn periodic(g: &Expr, k: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    use crate::functions::{acos, asin, atan, log, asinh, acosh, atanh};
    let Kind::Func(f, args) = g.kind() else { return None };
    if args.len() != 1 {
        return None;
    }
    let (a, b) = linear(&args[0], x)?;
    if !proved_nonzero(&a){return None;}
    if has(k, x) {
        return None;
    }
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

/// Solutions of `f == 0` restricted to `domain`, as a set.
fn solve_zero(f: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    let real = is_real_domain(domain);
    // product: union of the factors' solutions
    let (num, den) = crate::simplify::numer_denom(&crate::simplify::together(f));
    if !has(&num, x) {
        return Some(empty_set());
    }
    if let Kind::Mul(_, facs) = num.kind() {
        if facs.iter().any(|(b, _)| matches!(b.kind(), Kind::Func(..)) && has(b, x)) && facs.len() > 1 {
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
    if let Some(g) = single_generator(&num, x) {
        // polynomial in g = trig(a*x + b) (or exp in the complex domain)
        let t = Expr::dummy("solveset_t", crate::symbol::Assumptions::NONE);
        let mut map = rustc_hash::FxHashMap::default();
        map.insert(g.clone(), t.clone());
        let p = crate::subs::xreplace(&num, &map);
        if !has(&p, x) {
            let ks = crate::solve::solve(&p, &t)?;
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
        return None;
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
    let sols = crate::solve::solve(f, x)?;
    Some(filter_domain(sols, domain))
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
                None => unknown.push(s),
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
    solveset_with_context(f,x,domain,&SolveContext::default())
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
    if !has(f, x) {
        return Ok(if f.is_zero() { domain.clone() } else { empty_set() });
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
    crit.extend(real_roots(&num, x)?);
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
            let sat = if !poles.iter().any(|q| cmp3(q,p)==Some(Ordering::Equal)) && in_dom(p) {
                holds(op,sign_at(&f,x,p)?)
            } else {false};
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

fn solve_univariate_inequality_set(ineq: &Expr, x: &Expr, domain: &Expr) -> Result<Expr, String> {
    match ineq.kind() {
        Kind::Rel(RelOp::Eq, a, b) => return solveset(&sub(a, b), x, domain),
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
            return solve_univariate_inequality_set(&reduced,x,domain);
        }
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
