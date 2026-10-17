//! Limit manipulations of held `Sum`/`Product` (and `Integral`) nodes, after
//! SymPy's `ExprWithLimits` and `ExprWithIntLimits`: `index`,
//! `reverse_order`, `change_index`, `reorder_limit` and the
//! `has_finite_limits`, `has_reversed_limits` and `has_empty_sequence`
//! properties.

use crate::assumptions as q;
use crate::binders::{held, limits};
use crate::expr::{Expr, Func, Kind};
use crate::ops::{add, div, inv, mul, neg, sub};

/// Errors named after the Python exception raised by SymPy.
#[derive(Debug, Clone, PartialEq)]
pub enum LimitsError {
    Value(String),
    Reorder(String),
}

fn parts(e: &Expr, int_limits: bool) -> Result<(Func, Vec<Expr>), LimitsError> {
    match e.kind() {
        Kind::Func(f @ (Func::Sum | Func::Product), args) => Ok((f.clone(), args.clone())),
        Kind::Func(Func::Integral, args) if !int_limits => Ok((Func::Integral, args.clone())),
        _ => Err(LimitsError::Value(format!("{e} has no integer limits"))),
    }
}

fn rebuild(f: Func, body: &Expr, all: Vec<Vec<Expr>>) -> Result<Expr, LimitsError> {
    held(f, body, all).map_err(LimitsError::Value)
}

/// `expr.index(x)`: the position of the limit whose variable is `x`.
pub fn index(e: &Expr, x: &Expr) -> Result<usize, LimitsError> {
    let (_, args) = parts(e, true)?;
    let variables: Vec<Expr> = limits(&args).into_iter().map(|l| l[0].clone()).collect();
    match variables.iter().filter(|v| *v == x).count() {
        1 => Ok(variables.iter().position(|v| v == x).unwrap()),
        _ => Err(LimitsError::Value(format!("({e}, 'Number of instances of variable not equal to one')"))),
    }
}

/// `expr.reverse_order(*indices)`: each selected limit `(k, a, b)` becomes
/// `(k, b + 1, a - 1)` (Karr's convention), negating a Sum and inverting
/// the term of a Product once per reversed limit.
pub fn reverse_order(e: &Expr, indices: &[usize]) -> Result<Expr, LimitsError> {
    let (f, args) = parts(e, true)?;
    let mut flips = 0usize;
    let all: Vec<Vec<Expr>> = limits(&args).into_iter().enumerate().map(|(i, l)| {
        if indices.contains(&i) {
            flips += 1;
            vec![l[0].clone(), add(&[l[2].clone(), Expr::one()]), sub(&l[1], &Expr::one())]
        } else {
            l
        }
    }).collect();
    let odd = flips % 2 == 1;
    let body = match (&f, odd) {
        (_, false) => args[0].clone(),
        (Func::Sum, true) => neg(&args[0]),
        (_, true) => inv(&args[0]),
    };
    rebuild(f, &body, all)
}

/// `expr.change_index(var, trafo, newvar)`: substitute the linear index
/// transformation `trafo` (`+-var + beta`) into the limits and the term.
pub fn change_index(e: &Expr, var: &Expr, trafo: &Expr, newvar: Option<&Expr>) -> Result<Expr, LimitsError> {
    let (f, args) = parts(e, true)?;
    let newvar = newvar.unwrap_or(var).clone();
    let poly = crate::poly::polyobj::Poly::from_expr(trafo, std::slice::from_ref(var))
        .map_err(|_| LimitsError::Value(format!("{trafo} is not a polynomial in {var}")))?;
    if poly.degree(0) != Some(1) {
        return Err(LimitsError::Value("Index transformation is not linear".into()));
    }
    let alpha = crate::diff::diff(trafo, var);
    let beta = crate::subs::subs1(trafo, var, &Expr::zero());
    let mut found = false;
    let mut all = Vec::new();
    for l in limits(&args) {
        if l[0] != *var {
            all.push(l);
            continue;
        }
        found = true;
        let image = |x: &Expr| add(&[mul(&[alpha.clone(), x.clone()]), beta.clone()]);
        if crate::traverse::free_symbols(&alpha).is_empty() && !alpha.is_one() && alpha != Expr::minus_one() {
            return Err(LimitsError::Value("Linear transformation results in non-linear summation stepsize".into()));
        }
        if alpha.is_one() {
            all.push(vec![newvar.clone(), image(&l[1]), image(&l[2])]);
        } else {
            all.push(vec![newvar.clone(), image(&l[2]), image(&l[1])]);
        }
    }
    if !found {
        return Err(LimitsError::Value(format!("{var} is not an index of {e}")));
    }
    let body = crate::subs::subs1(&args[0], var, &div(&sub(var, &beta), &alpha));
    let body = crate::subs::subs1(&body, var, &newvar);
    rebuild(f, &body, all)
}

/// `expr.reorder_limit(x, y)`: swap two limits whose bounds do not involve
/// any index variable.
pub fn reorder_limit(e: &Expr, x: usize, y: usize) -> Result<Expr, LimitsError> {
    let (f, args) = parts(e, true)?;
    let all = limits(&args);
    let out_of_range = || LimitsError::Value("tuple index out of range".into());
    let (lx, ly) = (all.get(x).ok_or_else(out_of_range)?.clone(), all.get(y).ok_or_else(out_of_range)?.clone());
    let variables: Vec<Expr> = all.iter().map(|l| l[0].clone()).collect();
    let free = |b: &Expr| crate::traverse::free_symbols(b).iter().any(|s| variables.contains(s));
    if [&lx[1], &lx[2], &ly[1], &ly[2]].into_iter().any(free) {
        return Err(LimitsError::Reorder("could not interchange the two limits specified".into()));
    }
    let swapped = all.into_iter().enumerate().map(|(i, l)| if i == x { ly.clone() } else if i == y { lx.clone() } else { l }).collect();
    rebuild(f, &args[0], swapped)
}

fn is_infinite(e: &Expr) -> Option<bool> {
    q::is_finite(e).map(|finite| !finite)
}

/// `expr.has_finite_limits`.
pub fn has_finite_limits(e: &Expr) -> Result<Option<bool>, LimitsError> {
    let (_, args) = parts(e, false)?;
    let mut unknown = false;
    for l in limits(&args) {
        if l.len() == 3 {
            if l[1..].iter().any(|b| is_infinite(b) == Some(true)) {
                return Ok(Some(false));
            }
            if l[1..].iter().any(|b| is_infinite(b).is_none()) && is_infinite(&l[0]).is_none() {
                unknown = true;
            }
        } else if is_infinite(&l[0]).is_none() {
            unknown = true;
        }
    }
    Ok(if unknown { None } else { Some(true) })
}

/// `expr.has_reversed_limits`.
pub fn has_reversed_limits(e: &Expr) -> Result<Option<bool>, LimitsError> {
    let (_, args) = parts(e, false)?;
    let mut unknown = false;
    for l in limits(&args) {
        if l.len() == 2 {
            unknown = true;
            continue;
        }
        if l.len() != 3 {
            return Ok(None);
        }
        let dif = sub(&l[2], &l[1]);
        if extended_negative(&dif) == Some(true) {
            return Ok(Some(true));
        }
        if extended_negative(&dif) != Some(false) {
            unknown = true;
        }
    }
    Ok(if unknown { None } else { Some(false) })
}

/// `is_extended_negative`, with `-oo` negative and `nan` undecided.
fn extended_negative(e: &Expr) -> Option<bool> {
    use crate::expr::Constant;
    if e.is_constant(Constant::NegInfinity) {
        return Some(true);
    }
    if e.is_constant(Constant::Infinity) {
        return Some(false);
    }
    if e.is_constant(Constant::NaN) || e.is_constant(Constant::ComplexInfinity) {
        return None;
    }
    if q::is_negative(e) == Some(true) {
        return Some(true);
    }
    if q::is_nonnegative(e) == Some(true) {
        return Some(false);
    }
    None
}

/// `expr.has_empty_sequence`: some limit has `a - b == 1`.
pub fn has_empty_sequence(e: &Expr) -> Result<Option<bool>, LimitsError> {
    let (_, args) = parts(e, true)?;
    let mut unknown = false;
    for l in limits(&args) {
        match q::is_zero(&sub(&sub(&l[1], &l[2]), &Expr::one())) {
            Some(true) => return Ok(Some(true)),
            Some(false) => {}
            None => unknown = true,
        }
    }
    Ok(if unknown { None } else { Some(false) })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn e(s: &str) -> Expr {
        parse_str(s).unwrap()
    }

    #[test]
    fn limit_manipulations_match_sympy() {
        let k = Expr::symbol("k");
        assert_eq!(reverse_order(&e("Sum(k, (k, 1, n))"), &[0]).unwrap().to_string(), "Sum(-k, (k, n + 1, 0))");
        assert_eq!(reverse_order(&e("Product(k, (k, 1, n))"), &[0]).unwrap().to_string(), "Product(1/k, (k, n + 1, 0))");
        assert_eq!(change_index(&e("Sum(k, (k, 1, n))"), &k, &e("k + 1"), None).unwrap().to_string(), "Sum(k - 1, (k, 2, n + 1))");
        assert_eq!(index(&e("Sum(k*j, (k, 1, n), (j, 1, m))"), &e("j")).unwrap(), 1);
        assert_eq!(has_empty_sequence(&e("Sum(k, (k, 1, 0))")).unwrap(), Some(true));
        assert_eq!(has_finite_limits(&e("Sum(k, (k, 1, oo))")).unwrap(), Some(false));
        assert!(matches!(reorder_limit(&e("Sum(k*j, (k, 1, j), (j, 1, m))"), 0, 1), Err(LimitsError::Reorder(_))));
    }
}
