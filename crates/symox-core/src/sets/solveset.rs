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

/// The only function of `x` in `f` if it is `sin`/`cos`/`tan`/`exp` of
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
    if !matches!(fun, Func::Sin | Func::Cos | Func::Tan | Func::Exp) || a.len() != 1 {
        return None;
    }
    linear(&a[0], x)?;
    Some(g)
}

fn has_periodic_in(e: &Expr, x: &Expr) -> bool {
    match e.kind() {
        Kind::Func(f, a) if matches!(f, Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc) => {
            a.iter().any(|t| has(t, x))
        }
        _ => e.args().iter().any(|a| has_periodic_in(a, x)),
    }
}

/// Reduce `c` modulo `p` into `[0, p)` (numerically decided).
fn reduce_mod(c: &Expr, p: &Expr) -> Option<Expr> {
    let cv = crate::evalf::eval_f64(c)?;
    let pv = crate::evalf::eval_f64(p)?;
    let k = (cv / pv + 1e-12).floor();
    if k == 0.0 {
        return Some(c.clone());
    }
    Some(crate::expand::expand(&sub(c, &mul2(&Expr::int(k as i64), p))))
}

fn image(period: &Expr, c: &Expr, context: &SolveContext) -> Expr {
    let n = context.periodic_variable();
    image_set(&lambda(n, &add2(&mul2(period, n), c)), &integers())
}

/// Periodic solutions of `trig(a*x + b) = k` / `exp(a*x + b) = k`.
fn periodic(g: &Expr, k: &Expr, x: &Expr, domain: &Expr, context: &SolveContext) -> Option<Expr> {
    use crate::functions::{acos, asin, atan, log};
    let Kind::Func(f, args) = g.kind() else { return None };
    if args.len() != 1 {
        return None;
    }
    let (a, b) = linear(&args[0], x)?;
    if has(k, x) {
        return None;
    }
    let pi = Expr::pi();
    let two_pi = mul2(&Expr::int(2), &pi);
    let real = is_real_domain(domain);
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
            let kv = crate::evalf::eval_c64(k)?;
            if kv.im.abs() > 1e-12 {
                return if real { Some(empty_set()) } else { None };
            }
            if kv.re.abs() > 1.0 + 1e-12 {
                return if real { Some(empty_set()) } else { None };
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
            if !is_numeric_value(k) || crate::assumptions::is_real(k) != Some(true) {
                return None;
            }
            thetas = vec![atan(k)];
            base_period = pi.clone();
        }
        Func::Exp => {
            if real {
                // a*x + b = log(k) for k > 0
                if is_numeric_value(k) {
                    return Some(match cmp3(k, &Expr::zero()) {
                        Some(Ordering::Greater) => finite_set(vec![crate::expand::expand(&div(&sub(&log(k), &b), &a))]),
                        Some(_) => empty_set(),
                        None => return None,
                    });
                }
                return Some(filter_domain(vec![crate::expand::expand(&div(&sub(&log(k), &b), &a))], domain));
            }
            if k.is_zero() {
                return Some(empty_set());
            }
            // a*x + b = log(k) + 2*n*pi*I
            let period = div(&mul2(&two_pi, &Expr::i()), &a);
            let c = crate::expand::expand(&div(&sub(&log(k), &b), &a));
            return Some(image(&period, &c, context));
        }
        _ => return None,
    }
    let period = div(&base_period, &a);
    let period = if crate::evalf::eval_f64(&period)? < 0.0 { neg(&period) } else { period };
    let mut cs: Vec<Expr> = Vec::new();
    for t in thetas {
        let t = reduce_mod(&t, &base_period)?;
        let c = crate::expand::expand(&div(&sub(&t, &b), &a));
        let c = reduce_mod(&c, &period)?;
        if !cs.iter().any(|d| cmp3(d, &c) == Some(Ordering::Equal)) {
            cs.push(c);
        }
    }
    Some(union(cs.iter().map(|c| image(&period, c, context)).collect()))
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
        let t = Expr::symbol("_solveset_t");
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
        let sat = in_dom(&t) && sign_at(&f, x, &t).is_some_and(|s| holds(op, s));
        segs.push((Seg::Gap(lo, hi), sat));
        if i < k {
            let p = &pts[i];
            let sat = !poles.iter().any(|q| cmp3(q, p) == Some(Ordering::Equal))
                && in_dom(p)
                && sign_at(&f, x, p).is_some_and(|s| holds(op, s));
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
    let op = match a.as_number() {
        Some(n) if n.is_negative() => crate::logic::reversed_op(*op),
        Some(_) => *op,
        None => return None,
    };
    Some(crate::relational::rel(op, x, &rhs))
}
