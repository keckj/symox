//! `Expr.series` of the held operations `Subs`, `Derivative` and
//! `Integral` at `0` from the right: their `_eval_nseries` (the expansion
//! of the operand, every term wrapped again) under the reference
//! `Expr.series` loop (more terms requested while the order term falls
//! short, `collect(s, x) + O`).

use super::ser::{SErr, SRes};
use super::{get_o, getn, order, remove_o};
use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, div, pow};
use crate::traverse::free_symbols;

/// Whether `e` is a held operation expanded here.
pub fn is_held(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Derivative(..) | Kind::Func(Func::Subs | Func::Integral, _))
}

fn nseries(e: &Expr, x: &Expr, n: i64) -> SRes<Expr> {
    if !free_symbols(e).contains(x) {
        return Ok(e.clone());
    }
    super::series(e, x, &Expr::zero(), n, false)
}

fn terms_of(e: &Expr) -> Vec<Expr> {
    let r = remove_o(e);
    if r.is_add() {
        r.args()
    } else if r.is_zero() {
        vec![]
    } else {
        vec![r]
    }
}

/// `_eval_nseries` of a held operation.
fn eval_nseries(e: &Expr, x: &Expr, n: i64) -> SRes<Expr> {
    match e.kind() {
        Kind::Func(Func::Subs, args) => {
            // Subs._eval_nseries: expand in the variable whose point is x
            let vars = args[1].args();
            let points = args[2].args();
            let other = match points.iter().position(|p| p == x) {
                Some(i) => vars[i].clone(),
                None => x.clone(),
            };
            let arg = nseries(&args[0], &other, n)?;
            let mut rv: Vec<Expr> = terms_of(&arg).into_iter().map(|t| crate::subs::held_subs(&t, vars.clone(), points.clone())).collect();
            if let Some(o) = get_o(&arg) {
                rv.push(crate::subs::subs1(&o, &other, x));
            }
            Ok(add(&rv))
        }
        Kind::Derivative(inner, vars) => {
            // Derivative._eval_nseries: the terms differentiated unevaluated,
            // the order term divided by x
            let arg = nseries(inner, x, n)?;
            let mut rv: Vec<Expr> = terms_of(&arg).into_iter().map(|t| crate::diff::derivative_unevaluated(&t, vars.clone())).collect();
            if let Some(o) = get_o(&arg) {
                rv.push(div(&o, x));
            }
            Ok(add(&rv))
        }
        Kind::Func(Func::Integral, args) => {
            // Integral._eval_nseries: expand in the integration variable
            // whose limits hold x, integrate the terms, the order times x
            let limits = crate::binders::limits(args);
            let symb = limits.iter().find(|l| l[1..].iter().any(|b| b == x)).map_or_else(|| x.clone(), |l| l[0].clone());
            let f = nseries(&args[0], &symb, n)?;
            let terms = remove_o(&f);
            let mut iargs = args.clone();
            iargs[0] = terms;
            let integrated = crate::functions::doit(&Expr::from_kind(Kind::Func(Func::Integral, iargs)), false);
            let mut rv = vec![integrated];
            if let Some(o) = get_o(&f) {
                rv.push(crate::ops::mul2(&crate::subs::subs1(&o, &symb, x), x));
            }
            Ok(add(&rv))
        }
        _ => Err(SErr::Unsupported(e.to_string())),
    }
}

fn getn_int(e: &Expr) -> Option<i64> {
    getn(e).and_then(|k| k.as_number().cloned()).filter(|k| k.is_integer()).and_then(|k| k.to_i64())
}

/// `Expr.series(x, 0, n)` of a held operation (`None` when not handled).
pub fn series_held(e: &Expr, x: &Expr, n: i64) -> Option<SRes<Expr>> {
    if !is_held(e) || !x.is_symbol() {
        return None;
    }
    if !free_symbols(e).contains(x) {
        return Some(Ok(e.clone()));
    }
    Some(series_loop(e, x, n))
}

fn series_loop(e: &Expr, x: &Expr, n: i64) -> SRes<Expr> {
    let ox = |k: &Expr| order(&pow(x, k), &[(x.clone(), Expr::zero())]);
    let nn = n;
    let mut s1 = eval_nseries(e, x, n)?;
    let o;
    if get_o(&s1).is_some() {
        let Some(ngot) = getn_int(&s1) else { return Ok(s1) };
        if ngot > nn {
            let o0 = get_o(&s1).unwrap();
            if n != 0 {
                let r = Expr::num(Number::int(n).div(&Number::int(ngot)).unwrap());
                s1 = add(&[s1.clone(), crate::subs::subs1(&o0, x, &pow(x, &r))]);
            } else {
                s1 = add(&[s1.clone(), order(&Expr::one(), &[(x.clone(), Expr::zero())])]);
            }
        } else if ngot < nn {
            let mut done = false;
            for more in 1..9i64 {
                s1 = eval_nseries(e, x, n + more)?;
                let Some(newn) = getn_int(&s1) else { return Ok(s1) };
                if newn != ngot {
                    // ndo = n + ceiling((n - ngot)*more/(newn - ngot))
                    let (num, den) = if newn > ngot { ((n - ngot) * more, newn - ngot) } else { ((ngot - n) * more, ngot - newn) };
                    let mut ndo = n - (-num).div_euclid(den);
                    s1 = eval_nseries(e, x, ndo)?;
                    let mut guard = 0;
                    while getn_int(&s1).is_some_and(|k| k < nn) {
                        s1 = eval_nseries(e, x, ndo)?;
                        ndo += 1;
                        guard += 1;
                        if guard > 32 {
                            return Err(SErr::Unsupported("series".into()));
                        }
                    }
                    done = true;
                    break;
                }
            }
            if !done {
                return Err(SErr::Unsupported(format!("Could not calculate {n} terms for {e}")));
            }
            s1 = add(&[s1.clone(), ox(&Expr::int(n))]);
        }
        o = get_o(&s1).unwrap_or_else(Expr::zero);
        s1 = remove_o(&s1);
    } else if super::contains_order(&s1) {
        return Ok(s1);
    } else {
        let on = ox(&Expr::int(n));
        let done = crate::functions::doit(&s1, true);
        o = if remove_o(&add(&[done.clone(), on.clone()])) == done { Expr::zero() } else { on };
    }
    Ok(add(&[crate::poly::api::collect(&s1, std::slice::from_ref(x)), o]))
}

/// `Derivative._eval_as_leading_term`: the first term of the operand's
/// expansion (in increasing order) whose derivative does not vanish gives
/// the leading term; `0` when every term vanishes (`Derivative(x**3, y)`).
/// A held operation free of `x` is its own leading term. `None` when `e`
/// is not handled here.
pub fn leading_term_held(e: &Expr, x: &Expr) -> Option<SRes<super::LeadTerm>> {
    if !is_held(e) || !x.is_symbol() {
        return None;
    }
    let constant = |c: Expr| super::LeadTerm { coeff: c.clone(), coeff_logt: c, logt: super::logt_symbol(), exp: Number::ZERO };
    if !free_symbols(e).contains(x) {
        return Some(Ok(constant(e.clone())));
    }
    let Kind::Derivative(inner, vars) = e.kind() else { return None };
    let mut n = 6;
    for _ in 0..4 {
        let s = match nseries(inner, x, n) {
            Ok(s) => s,
            Err(err) => return Some(Err(err)),
        };
        let mut terms: Vec<(Number, Expr)> = Vec::new();
        for t in terms_of(&s) {
            match super::leading_term(&t, x, &Expr::zero(), false) {
                Ok(lt) => terms.push((lt.exp, t)),
                Err(err) => return Some(Err(err)),
            }
        }
        terms.sort_by(|a, b| a.0.cmp_value(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        for (_, t) in &terms {
            let d = crate::functions::doit(&crate::diff::derivative_unevaluated(t, vars.clone()), false);
            if !d.is_zero() {
                return Some(super::leading_term(&d, x, &Expr::zero(), false));
            }
        }
        if get_o(&s).is_none() {
            return Some(Ok(constant(Expr::zero())));
        }
        n *= 2;
    }
    Some(Err(SErr::Unsupported(format!("leading term of {e}"))))
}
