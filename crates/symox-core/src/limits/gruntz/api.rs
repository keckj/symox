//! The public functions of the reference `sympy.series.gruntz` module
//! (`compare`, `sign`, `mrv`, `rewrite`, `mrv_leadterm`, `limitinf` and
//! `gruntz`) on top of the native engine.
//!
//! The engine works with its own variable `t -> +oo` (positive by
//! construction); the variable of the caller is mapped to it on entry and
//! back on exit. The mrv set is given, as the reference `SubsSet`, as pairs
//! `(element, dummy)` together with the expression rewritten in the dummies
//! and the rewrites `dummy -> exp(argument in the dummies)` of the
//! exponential elements; `rewrite` takes the same data back.

use std::cmp::Ordering;

use rustc_hash::FxHashMap;

use super::{is_inf, is_parametric_inf, Gruntz};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::limits::LimitError;
use crate::ops::{add, div, mul2, neg, pow, sub};
use crate::series::engine::depends;
use crate::subs::xreplace;

/// The exception class the reference raises.
#[derive(Debug, Clone)]
pub enum ApiError {
    NotImplemented(String),
    Value(String),
}

impl From<LimitError> for ApiError {
    fn from(e: LimitError) -> ApiError {
        match e {
            LimitError::Unknown(s) => ApiError::NotImplemented(s),
            err @ LimitError::Differ(..) => ApiError::Value(err.to_string()),
        }
    }
}

pub type ApiRes<T> = Result<T, ApiError>;

fn engine_var() -> Expr {
    Expr::symbol("_symox_lim_t")
}

fn swap(e: &Expr, from: &Expr, to: &Expr) -> Expr {
    let mut m = FxHashMap::default();
    m.insert(from.clone(), to.clone());
    xreplace(e, &m)
}

/// An applied undefined function of `x`: the reference `mrv` raises
/// `ValueError("MRV set computation for UndefinedFunction is not allowed")`.
fn undefined_function_of(e: &Expr, x: &Expr) -> bool {
    if matches!(e.kind(), Kind::Func(Func::Undefined(f), _) if !f.internal) && depends(e, x) {
        return true;
    }
    e.args().iter().any(|a| undefined_function_of(a, x))
}

fn check_undefined(e: &Expr, x: &Expr) -> ApiRes<()> {
    if undefined_function_of(e, x) {
        return Err(ApiError::Value("MRV set computation for UndefinedFunction is not allowed".into()));
    }
    Ok(())
}

fn exp_arg(e: &Expr) -> Option<&Expr> {
    match e.kind() {
        Kind::Func(Func::Exp, a) => Some(&a[0]),
        _ => None,
    }
}

/// `compare(a, b, x)`: `"<"`, `"="` or `">"` for the comparability classes
/// of `a` and `b` (`lim log(a)/log(b)` zero, finite or infinite).
pub fn compare(a: &Expr, b: &Expr, x: &Expr) -> ApiRes<&'static str> {
    let t = engine_var();
    let (at, bt) = (swap(a, x, &t), swap(b, x, &t));
    let mut g = Gruntz::new(&t);
    let lb = g.log_of(&bt);
    if lb.is_zero() {
        // log(b) = 0: the ratio is zoo (infinite) unless log(a) vanishes too
        return Ok(if g.log_of(&at).is_zero() { "=" } else { ">" });
    }
    let la = g.log_of(&at);
    if !depends(&la, &t) && !depends(&lb, &t) {
        // constants: the ratio is a nonzero number (zero only for a = 1)
        let c = div(&la, &lb);
        return Ok(if c.is_zero() {
            "<"
        } else if is_inf(&c) || c.is_constant(Constant::ComplexInfinity) {
            ">"
        } else {
            "="
        });
    }
    Ok(match g.compare(&at, &bt)? {
        Ordering::Less => "<",
        Ordering::Greater => ">",
        Ordering::Equal => "=",
    })
}

/// `sign(e, x)`: the sign of `e` for all large `x`.
pub fn sign(e: &Expr, x: &Expr) -> ApiRes<i32> {
    let t = engine_var();
    let mut g = Gruntz::new(&t);
    Ok(g.sign(&swap(e, x, &t))?)
}

/// `limitinf(e, x)`: the limit at `x -> +oo`.
pub fn limitinf(e: &Expr, x: &Expr) -> ApiRes<Expr> {
    check_undefined(e, x)?;
    let t = engine_var();
    let mut g = Gruntz::new(&t);
    let r = g.limitinf(&swap(e, x, &t))?;
    if crate::series::asym::contains_asym(&r) {
        return Err(ApiError::NotImplemented("auxiliary function in the limit value".into()));
    }
    Ok(swap(&r, &t, x))
}

/// `mrv_leadterm(e, x)`: `(c0, e0)` with `e ~ c0*w**e0` in the mrv variable.
pub fn mrv_leadterm(e: &Expr, x: &Expr) -> ApiRes<(Expr, Expr)> {
    check_undefined(e, x)?;
    let t = engine_var();
    let mut g = Gruntz::new(&t);
    let (c0, e0) = g.mrv_leadterm(&swap(e, x, &t))?;
    if crate::series::asym::contains_asym(&c0) {
        return Err(ApiError::NotImplemented("auxiliary function in the leading term".into()));
    }
    Ok((swap(&c0, &t, x), swap(&e0, &t, x)))
}

/// The mrv set of `e` as the reference `SubsSet`: `(element, dummy)` pairs,
/// the rewrites `dummy -> exp(argument in the dummies)` of the exponential
/// elements, and `e` with the elements replaced by their dummies.
pub type MrvSet = (Vec<(Expr, Expr)>, Vec<(Expr, Expr)>, Expr);

pub fn mrv(e: &Expr, x: &Expr) -> ApiRes<MrvSet> {
    check_undefined(e, x)?;
    use crate::simptools::powsimp::{powsimp, Combine, PowsimpOptions};
    let e = powsimp(e, &PowsimpOptions { deep: true, combine: Combine::Exp, force: false });
    let t = engine_var();
    let mut g = Gruntz::new(&t);
    let elems: Vec<Expr> = g.mrv(&swap(&e, x, &t))?.iter().map(|f| swap(f, &t, x)).collect();
    let pairs: Vec<(Expr, Expr)> =
        elems.iter().map(|f| (f.clone(), Expr::dummy("Dummy", crate::symbol::Assumptions(0)))).collect();
    let map: FxHashMap<Expr, Expr> = pairs.iter().cloned().collect();
    let rewrites = pairs
        .iter()
        .filter_map(|(f, d)| exp_arg(f).map(|a| (d.clone(), crate::functions::exp(&xreplace(a, &map)))))
        .collect();
    Ok((pairs.clone(), rewrites, xreplace(&e, &map)))
}

/// `rewrite(e, Omega, x, wsym)`: `e` (in the dummies of `omega`) rewritten
/// in the most rapidly varying `w -> 0+` and `log(w)`, as the reference:
/// the elements sorted by the height of their rewrite tree, `w` the last
/// one, every element `f` replaced by `exp(arg(f) - c*arg(g))*w**c` with
/// `c = lim arg(f)/arg(g)`, the powers of `w` cleared of the denominators
/// of the rational `c`, and every node normalized bottom up.
pub fn rewrite(e: &Expr, omega: &[(Expr, Expr)], rewrites: &[(Expr, Expr)], x: &Expr, wsym: &Expr) -> ApiRes<(Expr, Expr)> {
    if omega.is_empty() {
        return Err(ApiError::Value("Length cannot be 0".into()));
    }
    if omega.iter().any(|(f, _)| exp_arg(f).is_none()) {
        return Err(ApiError::Value("Value should be exp".into()));
    }
    let rw: FxHashMap<Expr, Expr> = rewrites.iter().cloned().collect();
    // build_expression_tree: the height of a node is one plus the heights
    // of the elements whose dummy occurs in its rewrite.
    let vars: Vec<Expr> = omega.iter().map(|(_, v)| v.clone()).collect();
    fn height(v: &Expr, vars: &[Expr], rw: &FxHashMap<Expr, Expr>, depth: usize) -> usize {
        let Some(r) = rw.get(v) else { return 1 };
        if depth > vars.len() {
            return 1;
        }
        1 + vars.iter().filter(|v2| crate::traverse::has(r, v2)).map(|v2| height(v2, vars, rw, depth + 1)).sum::<usize>()
    }
    let mut items: Vec<(usize, Expr, Expr)> =
        omega.iter().map(|(f, v)| (height(v, &vars, &rw, 0), f.clone(), v.clone())).collect();
    items.sort_by_key(|it| std::cmp::Reverse(it.0));
    let t = engine_var();
    let mut g = Gruntz::new(&t);
    let mut sig = 0;
    for (_, f, _) in &items {
        let a = swap(exp_arg(f).unwrap(), x, &t);
        sig = match g.sign(&a) {
            Ok(s @ (1 | -1)) => s,
            Ok(s) => return Err(ApiError::NotImplemented(format!("Result depends on the sign of {s}"))),
            Err(_) => {
                return Err(ApiError::NotImplemented(format!("Result depends on the sign of {}", crate::functions::sign(&swap(&a, &t, x)))))
            }
        };
    }
    let garg = exp_arg(&items.last().unwrap().1).unwrap().clone();
    let w = if sig == 1 { pow(wsym, &Expr::minus_one()) } else { wsym.clone() };
    let mut o2 = Vec::new();
    let mut den = num_bigint::BigInt::from(1);
    for (_, f, var) in &items {
        let farg = exp_arg(f).unwrap();
        let c = g.limitinf(&swap(&div(farg, &garg), x, &t))?;
        let c = swap(&c, &t, x);
        if let Some(n) = c.as_number().filter(|n| n.is_rational()) {
            let q = n.denom().to_bigint().unwrap_or_else(|| num_bigint::BigInt::from(1));
            den = num_integer::Integer::lcm(&den, &q);
        }
        let arg = match rw.get(var) {
            Some(r) => match exp_arg(r) {
                Some(a) => a.clone(),
                None => return Err(ApiError::Value("Value should be exp".into())),
            },
            None => farg.clone(),
        };
        o2.push((var.clone(), mul2(&crate::functions::exp(&sub(&arg, &mul2(&c, &garg))), &pow(&w, &c))));
    }
    use crate::simptools::powsimp::{powsimp, Combine, PowsimpOptions};
    let mut f = powsimp(e, &PowsimpOptions { deep: true, combine: Combine::Exp, force: false });
    for (var, b) in &o2 {
        f = swap(&f, var, b);
    }
    if items.iter().any(|(_, _, v)| crate::traverse::has(&f, v)) {
        return Err(ApiError::NotImplemented("an mrv element was not rewritten".into()));
    }
    let mut logw = if sig == 1 { neg(&garg) } else { garg.clone() };
    if den != num_bigint::BigInt::from(1) {
        let k = Expr::num(crate::number::Number::from_bigint(den));
        f = crate::subs::subs1(&f, &w, &pow(&w, &k));
        logw = div(&logw, &k);
    }
    let f = crate::simptools::bottom_up(&f, &mut |n: &Expr| crate::simptools::normal(n));
    Ok((f, logw))
}

/// `gruntz(e, z, z0, dir)`: the reference maps the limit point to `+oo`
/// (`I*oo` is taken as `oo`) and calls `limitinf`, with no heuristics.
pub fn gruntz(e: &Expr, z: &Expr, z0: &Expr, dir: &str) -> ApiRes<Expr> {
    let ioo = mul2(&Expr::i(), &Expr::oo());
    let e0 = if z0.is_constant(Constant::Infinity) || *z0 == ioo {
        e.clone()
    } else if z0.is_constant(Constant::NegInfinity) || *z0 == neg(&ioo) {
        swap(e, z, &neg(z))
    } else {
        let inv = pow(z, &Expr::minus_one());
        match dir {
            "-" => swap(e, z, &sub(z0, &inv)),
            "+" => swap(e, z, &add(&[z0.clone(), inv])),
            _ => return Err(ApiError::NotImplemented("dir must be '+' or '-'".into())),
        }
    };
    let r = limitinf(&e0, z)?;
    if is_parametric_inf(&r) || is_inf(&r) || !depends(&r, z) {
        return Ok(r);
    }
    Err(ApiError::NotImplemented("limit depends on the variable".into()))
}
