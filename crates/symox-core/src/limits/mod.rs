//! Limits with SymPy semantics.
//!
//! ```
//! use symox::prelude::*;
//! use symox::limits::{limit, Dir};
//!
//! let x = Expr::symbol("x");
//! let e = div(&sin(&x), &x);
//! assert_eq!(limit(&e, &x, &Expr::zero(), Dir::Plus).unwrap(), Expr::one());
//! ```
//!
//! Strategy: a leading term computed by the series engine answers most
//! limits at finite points (and rational functions at infinity) directly;
//! everything else is mapped to `x -> oo` and handed to the Gruntz
//! comparability class (mrv) algorithm in [`gruntz`]. Expressions with
//! bounded oscillating factors (`sin(x)/x` at `oo`) are handled by treating
//! the oscillating parts as bounded parameters.

pub mod exact_sign;
pub mod gruntz;
mod gamma_heuristic;
mod order_elsewhere;
pub mod trace;

use rustc_hash::FxHashMap;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::ops::{add, mul2, neg, pow};
use crate::series::engine::{degree_in, depends};
use crate::series::{leading_term_with, LeadTerm};
use crate::subs::xreplace;

/// Direction of a limit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dir {
    Plus,
    Minus,
    Both,
}

impl Dir {
    pub fn as_str(self)->&'static str {match self{Self::Plus=>"+",Self::Minus=>"-",Self::Both=>"+-"}}
    pub fn parse(s: &str) -> Option<Dir> {
        match s {
            "+" => Some(Dir::Plus),
            "-" => Some(Dir::Minus),
            "+-" | "-+" => Some(Dir::Both),
            _ => None,
        }
    }
}

/// Build the held operator without asserting continuity or evaluating it.
pub fn held_limit(e:&Expr,x:&Expr,x0:&Expr,dir:Dir)->Expr {
    let dir=if x0.is_constant(Constant::Infinity){Dir::Minus}else if x0.is_constant(Constant::NegInfinity){Dir::Plus}else{dir};
    Expr::from_kind(Kind::Func(Func::Limit,vec![e.clone(),x.clone(),x0.clone(),Expr::symbol(dir.as_str())]))
}

/// Validate held nodes after fallible public substitutions/reconstruction.
pub fn validate_held(e:&Expr)->Result<(),String>{
    if !e.has_limit(){return Ok(());}
    fn visit(e:&Expr,seen:&mut rustc_hash::FxHashMap<usize,Expr>)->Result<(),String>{
        if seen.insert(e.dag_id(),e.clone()).is_some(){return Ok(());}
        if let Kind::Func(Func::Limit,args)=e.kind(){
            if args.len()!=4 || !args[1].is_symbol(){return Err("held limits require a symbol variable".into());}
            if crate::traverse::has(&args[2],&args[1]){return Err("limits approaching a variable point are not supported".into());}
        }
        for arg in crate::inspection::scalar_children(e){visit(&arg,seen)?;}
        Ok(())
    }
    visit(e,&mut rustc_hash::FxHashMap::default())
}

/// Surface known nonexistence errors through Python's fallible doit API.
pub fn check_doit(e:&Expr,deep:bool)->LRes<()> {
    if !e.has_limit(){return Ok(());}
    fn visit(e:&Expr,deep:bool,seen:&mut rustc_hash::FxHashMap<usize,Expr>)->LRes<()>{
        if seen.insert(e.dag_id(),e.clone()).is_some(){return Ok(());}
        if deep {for arg in crate::inspection::scalar_children(e){visit(&arg,true,seen)?;}}
        if let Kind::Func(Func::Limit,args)=e.kind(){
            let dir=Dir::parse(args[3].as_sym().map(|s|s.name()).unwrap_or("+"));
            if let Some(dir)=dir {if let Err(error@LimitError::Differ(..))=limit(&args[0],&args[1],&args[2],dir){return Err(error);}}
        }
        Ok(())
    }
    visit(e,deep,&mut rustc_hash::FxHashMap::default())
}

#[derive(Clone, Debug, PartialEq)]
pub enum LimitError {
    /// The algorithm could not decide (return an unevaluated `Limit`).
    Unknown(String),
    /// Two-sided limit with different one-sided values.
    Differ(Expr, Expr),
}

impl std::fmt::Display for LimitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LimitError::Unknown(s) => write!(f, "cannot compute the limit: {s}"),
            LimitError::Differ(l, r) => write!(
                f,
                "The limit does not exist since left hand limit = {l} and right hand limit = {r}"
            ),
        }
    }
}

pub type LRes<T> = Result<T, LimitError>;

pub(crate) fn unknown<T>(s: impl Into<String>) -> LRes<T> {
    Err(LimitError::Unknown(s.into()))
}

fn is_infinite(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(c) if c.is_infinite())
        || matches!(e.kind(), Kind::Mul(_, f) if f.iter().any(|(b, _)| b.is_constant(Constant::Infinity)))
}

/// `limit(e, x, x0, dir)`.
pub fn limit(e: &Expr, x: &Expr, x0: &Expr, dir: Dir) -> LRes<Expr> {
    if !x.is_symbol(){return unknown("limits with non-symbol variables are not implemented");}
    if !depends(e, x) {
        return Ok(e.clone());
    }
    if x0.is_constant(Constant::NaN) {
        // Limit.doit: a nan limit point gives nan.
        return Ok(Expr::nan());
    }
    // `may_have_order` is a mask hint (symbol bits can collide with it):
    // the exact `contains_order` decides.
    if crate::series::contains_order(e) {
        use crate::series::order::{order, order_parts};
        if let Some((inner, vars)) = order_parts(e) {
            // Limit.doit: Order(limit(e.expr, z, z0), *e.args[1:]).
            return Ok(order(&limit(&inner, x, x0, Dir::Plus)?, &vars));
        }
        // An order term `O(g)` at the same variable and point tends to 0
        // when `g` does; the limit is then the limit of the rest
        // (`(1 + O(x)).limit(x, 0)` is 1, as the reference leading term).
        let ex = crate::expand::expand(e);
        let terms = if ex.is_add() { ex.args() } else { vec![ex.clone()] };
        let mut rest = Vec::new();
        let mut orders = Vec::new();
        for t in terms {
            if !crate::series::contains_order(&t) {
                rest.push(t);
                continue;
            }
            // `O(g)` or `f*O(g)` (that is `O(f*g)`) with `f` free of orders
            let parts = match t.kind() {
                Kind::Mul(c, facs) => {
                    let (o, f): (Vec<_>, Vec<_>) = facs.iter().partition(|(b, _)| crate::series::contains_order(b));
                    match o.as_slice() {
                        [(b, p)] if p.is_one() && f.iter().all(|(_, p)| !crate::series::contains_order(p)) => order_parts(b).map(|(g, vars)| {
                            let mut v = vec![Expr::num(c.clone()), g];
                            v.extend(f.iter().map(|(b, p)| crate::ops::mul_factor(b, p)));
                            (crate::ops::mul(&v), vars)
                        }),
                        _ => None,
                    }
                }
                _ => order_parts(&t),
            };
            match parts {
                Some((g, vars)) if !crate::series::contains_order(&g) => orders.push((g, vars)),
                _ => return unknown("order term"),
            }
        }
        if !orders.is_empty() {
            let mut kept = Vec::new();
            for (g, vars) in &orders {
                if order_elsewhere::drops_out(g, vars, x, x0) {
                    continue;
                }
                kept.push((g, vars));
            }
            for (g, vars) in kept {
                if vars.len() != 1 || vars[0].0 != *x || vars[0].1 != *x0 {
                    return unknown("order term at another point");
                }
                if !limit(g, x, x0, dir)?.is_zero() {
                    return unknown("order term that does not vanish at the limit point");
                }
            }
            return limit(&add(&rest), x, x0, dir);
        }
    }
    fn unknown_function(e:&Expr,x:&Expr)->bool {
        if matches!(e.kind(),Kind::Func(Func::Undefined(f),_) if !f.internal) && depends(e,x){return true;}
        e.args().iter().any(|a|unknown_function(a,x))
    }
    if unknown_function(e,x){return unknown("continuity of an undefined function is not known");}
    let infinite_point = x0.is_constant(Constant::Infinity) || x0.is_constant(Constant::NegInfinity);
    match dir {
        _ if infinite_point => limit_one_sided(e, x, x0, false),
        Dir::Plus => limit_one_sided(e, x, x0, false),
        Dir::Minus => limit_one_sided(e, x, x0, true),
        Dir::Both => {
            let r = limit_one_sided(e, x, x0, false)?;
            let l = limit_one_sided(e, x, x0, true)?;
            if r == l {
                return Ok(r);
            }
            if is_infinite(&r) && is_infinite(&l) {
                return Ok(Expr::zoo());
            }
            Err(LimitError::Differ(l, r))
        }
    }
}

fn limit_one_sided(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> LRes<Expr> {
    let r = limit_one_sided_raw(e, x, x0, minus)?;
    // A `nan` computed from a finite expression is an `oo - oo` of two
    // singular leading values (`Ei(-x) - Ei(-2*x)` at 0+, whose log
    // singularities cancel to `log(2)`): undecided, never the limit.
    if has_nan(&r) && !has_nan(e) {
        return unknown("indeterminate combination of singular leading values");
    }
    Ok(r)
}

fn has_nan(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(c) => *c == Constant::NaN,
        Kind::Num(_) | Kind::Sym(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(has_nan),
    }
}

fn limit_one_sided_raw(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> LRes<Expr> {
    if x0.is_constant(Constant::ComplexInfinity) || x0.is_constant(Constant::NaN) {
        return unknown("bad limit point");
    }
    if let Some(r) = fast_path(e, x, x0, minus).filter(|r| !has_nan(r)) {
        return Ok(r);
    }
    // map to x -> oo
    let t = Expr::symbol("_symox_lim_t");
    let inv = pow(&t, &Expr::minus_one());
    let v = if x0.is_constant(Constant::Infinity) {
        t.clone()
    } else if x0.is_constant(Constant::NegInfinity) {
        neg(&t)
    } else if minus {
        add(&[x0.clone(), neg(&inv)])
    } else {
        add(&[x0.clone(), inv])
    };
    let mut m = FxHashMap::default();
    m.insert(x.clone(), v);
    let et = xreplace(e, &m);
    let mut g = gruntz::Gruntz::new(&t);
    let r = match g.limitinf(&et) {
        Ok(r) => Ok(r),
        Err(err) => match bounded_heuristic(&et, &t).or_else(|| function_value(&et, &t)).or_else(|| negative_base_squeeze(&et, &t)).or_else(|| gamma_heuristic::gamma_value(&et, &t)).or_else(|| gamma_heuristic::loggamma_value(&et, &t)) {
            Some(r) => Ok(r),
            None => Err(err),
        },
    }?;
    // The auxiliary functions of the tractable forms never leave the engine.
    if crate::series::asym::contains_asym(&r) {
        return unknown("auxiliary function in the limit value");
    }
    Ok(r)
}

/// The function branch of the reference `heuristics` (`e.func(*[limit(a)
/// for a in e.args])`), restricted to functions continuous at the limit of
/// their argument, where it is the value of the limit:
/// * `harmonic(n, m)` with `m` fixed and `n -> +oo` is `harmonic(oo, m)`
///   (`zeta(m)` for `m > 1`, `oo` for `m <= 1` not proven negative, held
///   when undecided, `nan` and thus undecided for `m < 0`); with `n` tending to a real
///   number `n0 > -1` it is `harmonic(n0, m) = zeta(m) - zeta(m, n0 + 1)`,
///   continuous there.
/// * `subfactorial(a)` with `a -> +oo`: `uppergamma(a + 1, -1)/e`, which is
///   `gamma(a + 1)/e + O(1/a)` for a real order, so `oo`.
///
/// The arguments must be real for a real `t` (no `I`, every other symbol
/// real), so that the argument tends to its limit along the real axis.
fn function_value(e: &Expr, t: &Expr) -> Option<Expr> {
    use crate::special::{apply, SpecialFn as S};
    let Kind::Func(Func::Special(s @ (S::Harmonic | S::Subfactorial)), args) = e.kind() else { return None };
    let real_arg = |a: &Expr| {
        !crate::traverse::has(a, &Expr::i())
            && crate::traverse::free_symbols(a)
                .iter()
                .all(|v| v == t || crate::assumptions::is_real(v) == Some(true))
    };
    let arg = &args[0];
    if !depends(arg, t) || args[1..].iter().any(|a| depends(a, t)) || !real_arg(arg) {
        return None;
    }
    let mut g = gruntz::Gruntz::new(t);
    let l = g.limitinf(arg).ok()?;
    let r = match s {
        S::Subfactorial if l.is_constant(Constant::Infinity) => Expr::oo(),
        S::Harmonic => {
            let mut at = args.clone();
            if l.is_constant(Constant::Infinity) {
                at[0] = Expr::oo();
            } else {
                let n0 = l.as_number()?;
                if !n0.is_rational() || n0.cmp_value(&crate::number::Number::int(-1)) != Some(std::cmp::Ordering::Greater) {
                    return None;
                }
                at[0] = l.clone();
            }
            apply(S::Harmonic, at)
        }
        _ => return None,
    };
    if has_nan(&r) {
        return None;
    }
    Some(r)
}

/// Value of the limit from a leading term `c * t**e` (with `t -> 0+`
/// and `log(t)` written as `lt.logt` in `lt.coeff_logt`).
pub(crate) fn from_leading(lt: &LeadTerm) -> Option<Expr> {
    let c = &lt.coeff_logt;
    // A coefficient that is not a polynomial in log(t) (x**b written as
    // exp(b*log(t)) with a symbolic b) may dominate t**exp: leave such
    // terms to gruntz, which needs the sign of b.
    let k = degree_in(c, &lt.logt)?;
    if lt.exp.is_positive() {
        return Some(Expr::zero());
    }
    if lt.exp.is_zero() && k == 0 {
        return Some(c.clone());
    }
    // infinite: sign of the dominant coefficient (log(t) -> -oo)
    let ck = if k == 0 { c.clone() } else { coeff_of_power(c, &lt.logt, k)? };
    let parity = if k % 2 == 0 { 1 } else { -1 };
    match crate::series::engine::sign_of(&crate::series::coeff::C::from_expr(ck.clone())) {
        Some(s) => Some(if s * parity > 0 { Expr::oo() } else { Expr::neg_oo() }),
        None => {
            if crate::functions::is_numeric_value(&ck) {
                // A complex number: oo*sign(coeff), as SymPy's leading-term
                // limit (`limit(I*x, x, oo)` is `oo*I`).
                let nonreal = crate::evalf::eval_c64(&ck).is_some_and(|z| z.im.abs() > 1e-12 * (1.0 + z.re.abs()));
                if !nonreal {
                    return None;
                }
            }
            if crate::traverse::free_symbols(&ck).is_empty() {
                return None;
            }
            let s = crate::functions::sign(&ck);
            Some(mul2(&mul2(&Expr::int(parity as i64), &s), &Expr::oo()))
        }
    }
}

/// Coefficient of `s**k` in a polynomial in `s`.
pub(crate) fn coeff_of_power(e: &Expr, s: &Expr, k: i64) -> Option<Expr> {
    let terms: Vec<Expr> = if e.is_add() { e.args() } else { vec![e.clone()] };
    let mut acc = Vec::new();
    for t in terms {
        if degree_in(&t, s)? == k {
            let sk = pow(s, &Expr::int(k));
            acc.push(crate::ops::div(&t, &sk));
        }
    }
    Some(add(&acc))
}

fn fast_path(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> Option<Expr> {
    if has_bad(e) {
        return None;
    }
    let logt = Expr::symbol("_symox_lim_logt");
    let lt = match leading_term_with(e, x, x0, minus, &logt, &logt) {
        Ok(lt) => lt,
        Err(_) if crate::series::vanishes_near(e, x, x0, minus) => return Some(Expr::zero()),
        Err(_) => return None,
    };
    if depends(&lt.coeff_logt, x) {
        return None;
    }
    from_leading(&lt)
}

fn has_bad(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(c) => c.is_infinite() || *c == Constant::NaN,
        Kind::Func(Func::Order, _) => true,
        Kind::Num(_) | Kind::Sym(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(has_bad),
    }
}

/// Replace oscillating `sin(u)`/`cos(u)` (with `u -> ±oo`) by bounded
/// parameters; if the limit exists and does not depend on them, it is
/// the limit of the original expression.
/// Oscillating `sin`/`cos` next to a rounding remainder: the reference
/// finds the rounding bounds only through the leading term of the whole
/// expression, which an oscillating `sin(x)` at `oo` defeats, and its
/// term-by-term heuristics then meet `oo - oo` (`sin(x) + floor(x) - x`
/// stays unevaluated there), so the combination is left undecided.
fn mixed_bounded(osc: &FxHashMap<Expr, Expr>, rounding: &Rounding) -> bool {
    !osc.is_empty() && !rounding.map.is_empty()
}

fn bounded_heuristic(e: &Expr, t: &Expr) -> Option<Expr> {
    let mut map: FxHashMap<Expr, Expr> = FxHashMap::default();
    let mut g = gruntz::Gruntz::new(t);
    // The tractable forms of Si, Ci and the Fresnel integrals at an infinite
    // argument hold their oscillating sin and cos.
    let prepared = g.prepare(e).ok()?;
    let e = &prepared;
    collect_oscillating(e, t, &mut g, &mut map);
    let mut params: Vec<Expr> = map.values().cloned().collect();
    let mut rounding = Rounding::new();
    rounding.collect(e, t, &mut g);
    if mixed_bounded(&map, &rounding) {
        return None;
    }
    params.extend(rounding.params.iter().map(|(p, _, _)| p.clone()));
    map.extend(rounding.map);
    if map.is_empty() {
        return None;
    }
    let e2 = xreplace(e, &map);
    // A bounded value is not a constant where it may reach 0: in a log or
    // under a non-polynomial power (sin(x)**(1/x) has no limit, as
    // sin(x) vanishes infinitely often).
    if singular_use(&e2, &params, t) {
        return None;
    }
    let mut g2 = gruntz::Gruntz::new(t);
    let r = g2.limitinf(&e2).ok()?;
    if params.iter().any(|p| depends(&r, p)) {
        return None;
    }
    Some(r)
}

/// Limit of an expression affine in its oscillating parts, as
/// `c0 + sum(c_i * w_i)`: the tractable form is `a0 + sum(a_i * w_i)` where
/// every `w_i` is a `sin(u)` or `cos(u)` with a real infinite limit of `u`
/// (accumulating on `[-1, 1]`), and each `a_i` has a finite limit `c_i`.
/// The bounded remainders of `floor`, `ceiling` and `Mod` at an infinite
/// argument ([`Rounding`]) are further parameters, on `[-1, 0]` or `[0, 1]`.
/// Returns `(c0, [(c_i, lo_i, hi_i)])` with the zero `c_i` left out; an infinite `c0`
/// dominates the bounded part and comes back alone. The caller forms the
/// accumulation bounds, as SymPy's leading terms do (`sin` at an infinite
/// argument is `AccumBounds(-1, 1)`), e.g. `x*Ci(x)` at `oo` is `sin(x)`
/// plus a vanishing part.
pub fn bounded_affine(e: &Expr, x: &Expr, x0: &Expr, minus: bool) -> Option<(Expr, Vec<(Expr, Expr, Expr)>)> {
    if !x.is_symbol() || x0.is_constant(Constant::ComplexInfinity) || x0.is_constant(Constant::NaN) {
        return None;
    }
    if !depends(e, x) || has_nan(e) {
        return None;
    }
    let t = Expr::symbol("_symox_lim_t");
    let et = to_infinity(e, x, x0, minus, &t);
    let mut g = gruntz::Gruntz::new(&t);
    let prepared = g.prepare(&et).ok()?;
    let mut map: FxHashMap<Expr, Expr> = FxHashMap::default();
    collect_real_oscillating(&prepared, &t, &mut g, &mut map);
    let mut rounding = Rounding::new();
    rounding.collect(&prepared, &t, &mut g);
    if map.is_empty() && rounding.map.is_empty() {
        return None;
    }
    if mixed_bounded(&map, &rounding) {
        return bounded_times_limit(&prepared, &t, &map);
    }
    let mut pairs: Vec<(Expr, Expr)> = map.into_iter().collect();
    pairs.sort_by_key(|(_, p)| p.as_sym().map(|s| s.name().to_string()));
    let mut params: Vec<Expr> = pairs.iter().map(|(_, p)| p.clone()).collect();
    let mut ranges: Vec<(Expr, Expr)> = vec![(Expr::minus_one(), Expr::one()); params.len()];
    for (p, lo, hi) in &rounding.params {
        params.push(p.clone());
        ranges.push((lo.clone(), hi.clone()));
    }
    let mut back: FxHashMap<Expr, Expr> = pairs.iter().map(|(w, p)| (w.clone(), p.clone())).collect();
    back.extend(rounding.map);
    let e2 = xreplace(&prepared, &back);
    // An integer power `k >= 2` of a bounded sine or cosine is the bounds
    // power `AccumBounds(-1, 1)**k`: `AccumBounds(0, 1)` for an even `k`,
    // `AccumBounds(-1, 1)` for an odd one (`sin(x)**2` at `oo`); each power
    // is its own bounded value, as in the reference bounds arithmetic.
    let sincos: Vec<Expr> = pairs.iter().map(|(_, p)| p.clone()).collect();
    let mut powers: Vec<(Expr, Expr)> = Vec::new();
    let e2 = bounded_powers(&e2, &sincos, &mut powers);
    for (_, q) in &powers {
        params.push(q.clone());
    }
    for (k, _) in &powers {
        let even = k.as_number().is_some_and(|k| k.to_i64().is_some_and(|k| k % 2 == 0));
        ranges.push(if even { (Expr::zero(), Expr::one()) } else { (Expr::minus_one(), Expr::one()) });
    }
    let (a0, coeffs) = affine_parts(&e2, &params)?;
    let lim = |a: &Expr| -> Option<Expr> {
        if !depends(a, &t) {
            return Some(a.clone());
        }
        let r = gruntz::Gruntz::new(&t).limitinf(a).ok()?;
        if has_nan(&r) || crate::series::asym::contains_asym(&r) || depends(&r, &t) {
            return None;
        }
        Some(r)
    };
    let c0 = lim(&a0)?;
    if params.iter().any(|p| depends(&c0, p)) {
        return None;
    }
    let mut cs = Vec::new();
    for (a, (lo, hi)) in coeffs.iter().zip(ranges) {
        let c = lim(a)?;
        if has_infinity(&c) || params.iter().any(|p| depends(&c, p)) {
            return None;
        }
        if !c.is_zero() {
            cs.push((c, lo, hi));
        }
    }
    if has_infinity(&c0) {
        return Some((c0, Vec::new()));
    }
    if cs.is_empty() {
        return None;
    }
    Some((c0, cs))
}

/// A product of bounded `sin`/`cos` factors (integer powers `k >= 1`) and
/// of a part `f` free of them, next to a rounding function: the
/// reference's term-wise heuristics meet `AccumBounds*0*oo = nan` there and
/// take `limit(f)*AccumBounds` (`cos(x)*floor(x)/x` at `oo` is
/// `AccumBounds(-1, 1)`, `cos(x)*floor(x)/x**2` is 0). Only a real numeric
/// `limit(f)` is accepted; anything else stays undecided.
fn bounded_times_limit(e: &Expr, t: &Expr, osc: &FxHashMap<Expr, Expr>) -> Option<(Expr, Vec<(Expr, Expr, Expr)>)> {
    let Kind::Mul(c, facs) = e.kind() else { return None };
    let mut range: Option<(i64, i64)> = None;
    let mut rest = vec![Expr::num(c.clone())];
    let times = |a: (i64, i64), b: (i64, i64)| {
        let v = [a.0 * b.0, a.0 * b.1, a.1 * b.0, a.1 * b.1];
        (*v.iter().min().unwrap(), *v.iter().max().unwrap())
    };
    for (b, k) in facs {
        if osc.contains_key(b) {
            let k = k.as_number().filter(|n| n.is_integer()).and_then(|n| n.to_i64()).filter(|&k| k >= 1)?;
            let r = if k % 2 == 0 { (0, 1) } else { (-1, 1) };
            range = Some(range.map_or(r, |a| times(a, r)));
        } else {
            let f = crate::ops::mul_factor(b, k);
            if osc.keys().any(|w| crate::traverse::has(&f, w)) {
                return None;
            }
            rest.push(f);
        }
    }
    let (lo, hi) = range?;
    let f = crate::ops::mul(&rest);
    let l = limit(&f, t, &Expr::oo(), Dir::Plus).ok()?;
    if l.is_zero() {
        // 0*AccumBounds is 0
        return Some((Expr::zero(), Vec::new()));
    }
    let numeric = l.is_number() && crate::assumptions::is_finite(&l) == Some(true) && crate::assumptions::is_real(&l) == Some(true);
    numeric.then(|| (Expr::zero(), vec![(l, Expr::int(lo), Expr::int(hi))]))
}

/// Replaces every `p**k` (`p` in `sincos`, `k` an integer `>= 2`) by a new
/// parameter, recorded in `out` as `(k, parameter)`.
fn bounded_powers(e: &Expr, sincos: &[Expr], out: &mut Vec<(Expr, Expr)>) -> Expr {
    if !sincos.iter().any(|p| depends(e, p)) {
        return e.clone();
    }
    let param = |k: &Expr, out: &mut Vec<(Expr, Expr)>| -> Expr {
        let q = Expr::symbol(&format!("_symox_bounded_pow{}", out.len()));
        out.push((k.clone(), q.clone()));
        q
    };
    let power = |k: &Expr| k.as_number().is_some_and(|n| n.is_integer() && n.to_i64().is_some_and(|n| n >= 2));
    match e.kind() {
        Kind::Pow(b, k) if sincos.contains(b) && power(k) => param(k, out),
        Kind::Mul(c, facs) => {
            let mut v = vec![Expr::num(c.clone())];
            for (b, k) in facs {
                if sincos.contains(b) && power(k) {
                    v.push(param(k, out));
                } else {
                    v.push(bounded_powers(&crate::ops::mul_factor(b, k), sincos, out));
                }
            }
            crate::ops::mul(&v)
        }
        Kind::Add(..) => add(&e.args().iter().map(|a| bounded_powers(a, sincos, out)).collect::<Vec<_>>()),
        _ => e.clone(),
    }
}

fn has_infinity(e: &Expr) -> bool {
    match e.kind() {
        Kind::Const(c) => c.is_infinite() || *c == Constant::NaN,
        Kind::Num(_) | Kind::Sym(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(has_infinity),
    }
}

/// `e = a0 + sum(a_i * p_i)` with no `p_j` in the `a`'s, or None.
fn affine_parts(e: &Expr, params: &[Expr]) -> Option<(Expr, Vec<Expr>)> {
    let n = params.len();
    if !params.iter().any(|p| depends(e, p)) {
        return Some((e.clone(), vec![Expr::zero(); n]));
    }
    if let Some(i) = params.iter().position(|p| p == e) {
        let mut v = vec![Expr::zero(); n];
        v[i] = Expr::one();
        return Some((Expr::zero(), v));
    }
    if e.is_add() {
        let mut a0 = Vec::new();
        let mut cs: Vec<Vec<Expr>> = vec![Vec::new(); n];
        for term in e.args() {
            let (b, bs) = affine_parts(&term, params)?;
            a0.push(b);
            for (c, b) in cs.iter_mut().zip(bs) {
                c.push(b);
            }
        }
        return Some((add(&a0), cs.iter().map(|c| add(c)).collect()));
    }
    if let Kind::Mul(..) = e.kind() {
        let factors = e.args();
        let with: Vec<usize> = (0..factors.len()).filter(|&k| params.iter().any(|p| depends(&factors[k], p))).collect();
        if with.len() != 1 {
            return None;
        }
        let k = with[0];
        let rest: Vec<Expr> = factors.iter().enumerate().filter(|&(j, _)| j != k).map(|(_, f)| f.clone()).collect();
        let r = crate::ops::mul(&rest);
        let (b, bs) = affine_parts(&factors[k], params)?;
        return Some((mul2(&r, &b), bs.iter().map(|c| mul2(&r, c)).collect()));
    }
    None
}

/// `collect_oscillating` restricted to arguments with a real infinite limit.
fn collect_real_oscillating(e: &Expr, t: &Expr, g: &mut gruntz::Gruntz, map: &mut FxHashMap<Expr, Expr>) {
    if !depends(e, t) || map.contains_key(e) {
        return;
    }
    if let Kind::Func(Func::Sin | Func::Cos, a) = e.kind() {
        if let Ok(l) = g.limitinf(&a[0]) {
            if l.is_constant(Constant::Infinity) || l.is_constant(Constant::NegInfinity) {
                let p = Expr::symbol(&format!("_symox_bounded{}", map.len()));
                map.insert(e.clone(), p);
                return;
            }
        }
    }
    for a in e.args() {
        collect_real_oscillating(&a, t, g, map);
    }
}

fn to_infinity(e: &Expr, x: &Expr, x0: &Expr, minus: bool, t: &Expr) -> Expr {
    let inv = pow(t, &Expr::minus_one());
    let v = if x0.is_constant(Constant::Infinity) {
        t.clone()
    } else if x0.is_constant(Constant::NegInfinity) {
        neg(t)
    } else if minus {
        add(&[x0.clone(), neg(&inv)])
    } else {
        add(&[x0.clone(), inv])
    };
    let mut m = FxHashMap::default();
    m.insert(x.clone(), v);
    xreplace(e, &m)
}

fn singular_use(e: &Expr, params: &[Expr], t: &Expr) -> bool {
    // a base or log argument made of the bounded values alone (sin(x)**(1/x),
    // log(cos(x) + 1)) may vanish or turn negative; one still depending on
    // the variable ((1 + x/(n + sin(n)))**n) is driven by it
    let has_param = |a: &Expr| !depends(a, t) && params.iter().any(|p| depends(a, p));
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => false,
        Kind::Func(Func::Log, a) if a.iter().any(has_param) => true,
        Kind::Pow(b, p) if has_param(b) && !p.as_number().is_some_and(|k| k.is_integer() && k.is_positive()) => true,
        Kind::Mul(_, f) if f.iter().any(|(b, p)| has_param(b) && !p.as_number().is_some_and(|k| k.is_integer() && k.is_positive())) => true,
        _ => e.args().iter().any(|a| singular_use(a, params, t)),
    }
}

/// `oo`, `-oo`, or `c*oo` with every other factor a real number or the sign
/// of a real expression (`oo*sign(m)` for a real `m`).
fn real_infinity(l: &Expr) -> bool {
    match l.kind() {
        Kind::Const(Constant::Infinity | Constant::NegInfinity) => true,
        Kind::Mul(c, f) => {
            !c.is_zero()
                && f.iter().any(|(b, _)| b.is_constant(Constant::Infinity))
                && f.iter().all(|(b, p)| {
                    b.is_constant(Constant::Infinity)
                        || (p.is_one()
                            && matches!(b.kind(), Kind::Func(Func::Sign, args)
                                if args.len() == 1 && crate::assumptions::is_real(&args[0]) == Some(true)))
                })
        }
        _ => false,
    }
}

/// The bounded remainders of `floor`, `ceiling` and `Mod` of an argument
/// with a real infinite limit, as the reference `RoundFunction._eval_nseries`
/// (`floor(u) = u + AccumBounds(-1, 0)`, `ceiling(u) = u + AccumBounds(0, 1)`
/// for an infinite `u`) and `Mod._eval_nseries` (the rewrite
/// `Mod(u, q) = u - q*floor(u/q)`, so `-q*AccumBounds(-1, 0)`).
///
/// Every occurrence gets its own parameter `p` with its range `(lo, hi)`:
/// two rounding functions of different arguments are independent bounded
/// values, as two accumulation bounds are (`floor(x + 1/2) - floor(x)` at
/// `oo` is `1/2 + AccumBounds(-1, 1)`); they never cancel. The map sends
/// the rounding function to its replacement, `u + p` or `-q*p`.
struct Rounding {
    map: FxHashMap<Expr, Expr>,
    params: Vec<(Expr, Expr, Expr)>,
}

impl Rounding {
    fn new() -> Self {
        Rounding { map: FxHashMap::default(), params: Vec::new() }
    }

    fn param(&mut self, lo: i64, hi: i64) -> Expr {
        let p = Expr::symbol(&format!("_symox_round{}", self.params.len()));
        self.params.push((p.clone(), Expr::int(lo), Expr::int(hi)));
        p
    }

    fn collect(&mut self, e: &Expr, t: &Expr, g: &mut gruntz::Gruntz) {
        if !depends(e, t) || self.map.contains_key(e) {
            return;
        }
        let real_infinite = |u: &Expr, g: &mut gruntz::Gruntz| {
            real_arg(u, t)
                && !has_rounding(u, t)
                && g.limitinf(u).is_ok_and(|l| l.is_constant(Constant::Infinity) || l.is_constant(Constant::NegInfinity))
        };
        match e.kind() {
            Kind::Func(f @ (Func::Floor | Func::Ceiling), a) if a.len() == 1 => {
                if real_infinite(&a[0], g) {
                    let p = if *f == Func::Floor { self.param(-1, 0) } else { self.param(0, 1) };
                    self.map.insert(e.clone(), add(&[a[0].clone(), p]));
                    return;
                }
            }
            Kind::Func(Func::Special(crate::special::SpecialFn::Mod), a) if a.len() == 2 => {
                let q = &a[1];
                if !depends(q, t) && real_arg(q, t) && crate::assumptions::is_zero(q) == Some(false) {
                    let u = crate::ops::div(&a[0], q);
                    if real_infinite(&u, g) {
                        let p = self.param(-1, 0);
                        self.map.insert(e.clone(), neg(&mul2(q, &p)));
                        return;
                    }
                }
            }
            _ => {}
        }
        for a in e.args() {
            self.collect(&a, t, g);
        }
    }
}

/// Squeeze through the modulus: `(-c)**u` with a rational `c > 0` and a real
/// `u` has the modulus `c**u`, so for a product `e` of powers `b**(+-u)`
/// with numeric bases, `|e| = |e'|` where `e'` replaces every `(-c)` by
/// `c`. When `e'` tends to 0, so does `e` (`(-1)**x*2**(-x)`, `(-1/4)**x`
/// at `oo`, issues 22181 and 19739).
fn negative_base_squeeze(e: &Expr, t: &Expr) -> Option<Expr> {
    // The reference reaches these limits only through `pow_heuristics`,
    // after `powsimp` has merged the product into a single power `B**g` of
    // a constant base: every factor is `b**g` or `b**(-g)` for one common
    // real exponent `g` and a numeric base, with no other coefficient
    // (`3*(-1)**x*2**(-x)`, `(-1)**(2*x)*2**(-x)` and `(-1)**x/x` raise
    // there, so they stay undecided here).
    let factors: Vec<(Expr, Expr)> = match e.kind() {
        Kind::Mul(c, f) if c.is_one() => f.clone(),
        Kind::Pow(b, p) => vec![(b.clone(), p.clone())],
        _ => return None,
    };
    // exp(u) is the power E**u
    let factors: Vec<(Expr, Expr)> = factors
        .into_iter()
        .map(|(b, p)| match b.kind() {
            Kind::Func(Func::Exp, a) if p.is_one() && a.len() == 1 => (Expr::e(), a[0].clone()),
            _ => (b, p),
        })
        .collect();
    let g = factors.first().map(|(_, p)| p.clone())?;
    if !depends(&g, t) || !real_arg(&g, t) {
        return None;
    }
    let ng = neg(&g);
    let mut negative = false;
    let mut v = Vec::with_capacity(factors.len());
    for (b, p) in &factors {
        if *p != g && *p != ng {
            return None;
        }
        let modulus = if let Some(c) = b.as_number() {
            if !c.is_rational() || c.is_zero() {
                return None;
            }
            if c.is_negative() {
                negative = true;
                Expr::num(c.neg())
            } else {
                b.clone()
            }
        } else if *b == Expr::e() {
            b.clone()
        } else {
            return None;
        };
        v.push(pow(&modulus, p));
    }
    if !negative {
        return None;
    }
    // |B**g| = |B|**g for a real g: a vanishing modulus squeezes the
    // oscillating sign to 0
    let modulus = crate::ops::mul(&v);
    let r = gruntz::Gruntz::new(t).limitinf(&modulus).ok()?;
    r.is_zero().then(Expr::zero)
}

/// A real value for a real `t`: no `I`, every other symbol real.
fn real_arg(a: &Expr, t: &Expr) -> bool {
    !crate::traverse::has(a, &Expr::i())
        && crate::traverse::free_symbols(a).iter().all(|v| v == t || crate::assumptions::is_real(v) == Some(true))
}

fn has_rounding(e: &Expr, t: &Expr) -> bool {
    if !depends(e, t) {
        return false;
    }
    matches!(e.kind(), Kind::Func(Func::Floor | Func::Ceiling | Func::Special(crate::special::SpecialFn::Mod | crate::special::SpecialFn::Frac), _))
        || e.args().iter().any(|a| has_rounding(a, t))
}

fn collect_oscillating(e: &Expr, t: &Expr, g: &mut gruntz::Gruntz, map: &mut FxHashMap<Expr, Expr>) {
    if !depends(e, t) || map.contains_key(e) {
        return;
    }
    if let Kind::Func(Func::Sin | Func::Cos, a) = e.kind() {
        if let Ok(l) = g.limitinf(&a[0]) {
            // only a real infinity keeps sin and cos bounded: cos(m*x) with a
            // complex m grows like cosh
            if real_infinity(&l) {
                let p = Expr::symbol(&format!("_symox_bounded{}", map.len()));
                map.insert(e.clone(), p);
                return;
            }
        }
    }
    for a in e.args() {
        collect_oscillating(&a, t, g, map);
    }
}

#[cfg(test)]
mod tests;
