//! Complete branch decomposition for factorable/triangular rational systems.
//!
//! Every division by an undecided coefficient is split into zero/nonzero
//! cases. No numerical zero test or generic-coefficient assumption is used.
use super::convert::{from_exprs_gens, PResult, PolyError};
use crate::expr::{Constant, Expr, Kind};
use crate::integrate::poly_coeffs;
use crate::ops::{add, div, mul, neg, pow, sub};
use crate::simplify::numer_denom;
use crate::subs::subs1;
use crate::traverse::{free_symbols, has};
use rustc_hash::FxHashMap;

#[derive(Clone, Debug)]
pub struct Branch {
    pub values: Vec<Expr>,
    /// Excluded values for each surviving free variable. Roots may depend on
    /// other free variables, as in x != y.
    pub excluded: Vec<Vec<Expr>>,
}

#[derive(Clone)]
struct State { equations: Vec<Expr>, nonzero: Vec<Expr>, values: Vec<Expr>, free: Vec<usize> }

fn unsupported(s: &str) -> PolyError { PolyError::Unsupported(s.into()) }
fn normal(e: &Expr) -> Expr {
    let e = crate::expand::expand(e);
    let e = super::api::cancel(&e).unwrap_or(e);
    crate::simplify::simplify(&e)
}
fn zero(e: &Expr) -> bool {
    e.is_zero() || (free_symbols(e).is_empty() && crate::assumptions::exact_queries(|| crate::assumptions::is_zero(e)) == Some(true))
}
fn nonzero(e: &Expr, restrictions: &[Expr]) -> bool {
    // nonlinsolve describes complex varieties, independently of assumptions
    // attached to its generator symbols (as opposed to scalar solve()).
    if free_symbols(e).is_empty() && crate::assumptions::exact_queries(|| crate::assumptions::is_zero(e)) == Some(false) { return true; }
    if restrictions.iter().any(|q| e == q || normal(&neg(e)) == *q) { return true; }
    match e.kind() {
        Kind::Const(c) => *c != Constant::NaN,
        Kind::Mul(c, factors) => !c.is_zero() && factors.iter().all(|(b, _)| nonzero(b, restrictions)),
        Kind::Pow(b, exponent) if exponent.as_i64().is_some_and(|n| n != 0) => nonzero(b, restrictions),
        _ => false,
    }
}
fn factors(e: &Expr) -> Vec<Expr> {
    let f = super::api::factor(e).unwrap_or_else(|_| e.clone());
    match f.kind() {
        Kind::Mul(_, bases) => bases.iter().map(|(b, _)| b.clone()).collect(),
        Kind::Pow(b, exponent) if exponent.as_i64().is_some_and(|n| n > 0) => vec![b.clone()],
        _ => vec![f],
    }
}
fn coefficients(e: &Expr, x: &Expr) -> Option<Vec<Expr>> {
    let mut c: Vec<_> = poly_coeffs(e, x)?.iter().map(normal).collect();
    while c.last().is_some_and(zero) { c.pop(); }
    Some(c)
}
fn from_coefficients(c: &[Expr], x: &Expr) -> Expr {
    add(&c.iter().enumerate().map(|(i, a)| mul(&[a.clone(), pow(x, &Expr::int(i as i64))])).collect::<Vec<_>>())
}
/// Linear and quadratic formulas have no hidden exceptional denominator once
/// their leading coefficient has been certified nonzero.
fn roots(c: &[Expr]) -> PResult<Vec<Expr>> {
    let result = match c.len() {
        2 => vec![normal(&div(&neg(&c[0]), &c[1]))],
        3 => {
            let d = sub(&pow(&c[1], &Expr::int(2)), &mul(&[Expr::int(4), c[2].clone(), c[0].clone()]));
            let radical = pow(&normal(&d), &Expr::rational(1,2));
            let denominator = mul(&[Expr::int(2), c[2].clone()]);
            vec![normal(&div(&sub(&neg(&c[1]), &radical), &denominator)), normal(&div(&add(&[neg(&c[1]), radical]), &denominator))]
        }
        _ => return Err(unsupported("complete parametric elimination supports factorable linear and quadratic equations only")),
    };
    let mut result = result;
    result.sort_by_cached_key(crate::sortkey::sort_key);
    result.dedup();
    Ok(result)
}

fn assign(mut s: State, variables: &[Expr], j: usize, root: &Expr) -> State {
    let substitute = |e: &Expr| normal(&subs1(e, &variables[j], root));
    s.equations = s.equations.iter().map(substitute).collect();
    s.nonzero = s.nonzero.iter().map(substitute).collect();
    s.values = s.values.iter().map(substitute).collect();
    s.values[j] = root.clone();
    s.free.retain(|&k| k != j);
    let (_, den) = numer_denom(root);
    if !den.is_one() { s.nonzero.push(den); }
    s
}

fn descend(mut s: State, variables: &[Expr], budget: &mut usize) -> PResult<Vec<Branch>> {
    if *budget == 0 { return Err(unsupported("parametric branch decomposition exceeded its 4096-node limit")); }
    *budget -= 1;
    let mut equations = Vec::new();
    for e in &s.equations {
        let (n, d) = numer_denom(e);
        if !d.is_one() { s.nonzero.push(d); }
        let n = normal(&n);
        if !zero(&n) { equations.push(n); }
    }
    s.equations = equations;
    let mut restrictions = Vec::new();
    for e in &s.nonzero {
        let (n, d) = numer_denom(e);
        for part in [n, d] {
            for factor in factors(&normal(&part)) {
                let factor = normal(&factor);
                if zero(&factor) { return Ok(vec![]); }
                if !nonzero(&factor, &[]) && !restrictions.contains(&factor) { restrictions.push(factor); }
            }
        }
    }
    s.nonzero = restrictions;
    for e in &s.equations {
        if s.free.iter().all(|&j| !has(e, &variables[j])) {
            if nonzero(e, &s.nonzero) { return Ok(vec![]); }
            return Err(unsupported("a branch imposes an unresolved condition on external parameters; include those parameters among the unknowns"));
        }
        if nonzero(e, &s.nonzero) { return Ok(vec![]); }
    }
    // A product vanishes iff at least one factor vanishes, including at all
    // exceptional coefficient values. Overlapping components are harmless.
    for i in 0..s.equations.len() {
        let fs = factors(&s.equations[i]);
        if fs.len() > 1 || fs.first().is_some_and(|f| *f != s.equations[i]) {
            let mut out = Vec::new();
            for factor in fs {
                let mut branch = s.clone();
                branch.equations[i] = factor;
                out.extend(descend(branch, variables, budget)?);
            }
            return Ok(out);
        }
    }
    if !s.equations.is_empty() {
        let mut choices = Vec::new();
        for (i, e) in s.equations.iter().enumerate() {
            for &j in &s.free {
                if let Some(c) = coefficients(e, &variables[j]) {
                    if (2..=3).contains(&c.len()) {
                        choices.push((!nonzero(c.last().unwrap(), &s.nonzero), c.len(), i, j, c));
                    }
                }
            }
        }
        choices.sort_by_key(|c| (c.0, c.1, c.2, c.3));
        let (_, _, i, j, c) = choices.into_iter().next().ok_or_else(|| unsupported("system is outside factorable linear/quadratic triangular parametric elimination"))?;
        let leading = c.last().unwrap();
        let mut out = Vec::new();
        if !nonzero(leading, &s.nonzero) {
            let mut exceptional = s.clone();
            exceptional.equations[i] = from_coefficients(&c[..c.len()-1], &variables[j]);
            exceptional.equations.push(leading.clone());
            out.extend(descend(exceptional, variables, budget)?);
            s.nonzero.push(leading.clone());
        }
        s.equations.remove(i);
        for root in roots(&c)? {
            out.extend(descend(assign(s.clone(), variables, j, &root), variables, budget)?);
        }
        return Ok(out);
    }
    // Turn each nonzero polynomial into excluded roots of a free parameter.
    // If its leading coefficient can vanish, split that case before division.
    let mut excluded = vec![Vec::new(); variables.len()];
    for i in 0..s.nonzero.len() {
        let e = &s.nonzero[i];
        let mut choices = Vec::new();
        for &j in &s.free {
            if let Some(c) = coefficients(e, &variables[j]) {
                if (2..=3).contains(&c.len()) {
                    choices.push((!nonzero(c.last().unwrap(), &s.nonzero), c.len(), j, c));
                }
            }
        }
        choices.sort_by_key(|c| (c.0, c.1, c.2));
        let (_, _, j, c) = choices.into_iter().next().ok_or_else(|| unsupported("parameter exclusion cannot be represented by complete linear/quadratic root domains"))?;
        if !nonzero(c.last().unwrap(), &s.nonzero) {
            let mut exceptional = s.clone();
            exceptional.equations.push(c.last().unwrap().clone());
            exceptional.nonzero[i] = from_coefficients(&c[..c.len()-1], &variables[j]);
            let mut out = descend(exceptional, variables, budget)?;
            s.nonzero.push(c.last().unwrap().clone());
            out.extend(descend(s, variables, budget)?);
            return Ok(out);
        }
        excluded[j].extend(roots(&c)?);
    }
    for roots in &mut excluded {
        roots.sort_by_cached_key(crate::sortkey::sort_key);
        roots.dedup();
    }
    Ok(vec![Branch { values: s.values, excluded }])
}

/// Decompose the complete complex zero set in the supported triangular class.
/// Unsupported branches invalidate the whole result, never just one component.
pub fn nonlinsolve(equations: &[Expr], variables: &[Expr]) -> PResult<Vec<Branch>> {
    crate::simplify::with_native_selection(|| nonlinsolve_native(equations, variables))
}

fn nonlinsolve_native(equations: &[Expr], variables: &[Expr]) -> PResult<Vec<Branch>> {
    if variables.is_empty() || variables.iter().enumerate().any(|(i,x)| !matches!(x.kind(), Kind::Sym(_)) || variables[..i].contains(x)) {
        return Err(PolyError::Failed("generators must be distinct symbols".into()));
    }
    let equations: Vec<_> = equations.iter().map(|e| match e.kind() {
        Kind::Bool(true) => Expr::zero(), Kind::Bool(false) => Expr::one(), _ => crate::solve::as_expr(e),
    }).collect();
    let parts: Vec<_> = equations.iter().map(numer_denom).collect();
    let numerators: Vec<_> = parts.iter().map(|(n,_)| n.clone()).collect();
    let all_parts: Vec<_> = parts.iter().flat_map(|(n,d)| [n.clone(),d.clone()]).collect();
    let converted = from_exprs_gens(&all_parts, variables)?;
    if converted.float { return Err(unsupported("polynomial system solving requires exact coefficients")); }
    let restrictions: Vec<_> = parts.into_iter().map(|(_,d)| d).collect();
    // Existing finite elimination remains preferable when all branches have
    // exact root/residual certificates. Parameter coefficients stay on the
    // branch-complete path: generic fraction-field pivots are insufficient.
    if converted.gens.len() == converted.ngens {
        match super::systems::solve_poly_system(&numerators,variables) {
            Ok(None) => return Ok(vec![]),
            Ok(Some(solutions)) => {
                let mut finite = Vec::new();
                let mut supported = true;
                for solution in solutions {
                    let mapping: FxHashMap<_,_> = variables.iter().cloned().zip(solution.iter().cloned()).collect();
                    let state = State { equations: vec![], nonzero: restrictions.iter().map(|e| crate::subs::xreplace(e,&mapping)).collect(), values: solution, free: vec![] };
                    match descend(state,variables,&mut 4096) {
                        Ok(rows) => finite.extend(rows),
                        Err(PolyError::Unsupported(_)) => { supported=false; break; },
                        Err(e) => return Err(e),
                    }
                }
                if supported { return Ok(finite); }
            },
            Err(PolyError::Unsupported(_)) => {},
            Err(e) => return Err(e),
        }
    }
    let s = State { equations: numerators, nonzero: restrictions, values: variables.to_vec(), free: (0..variables.len()).collect() };
    let mut branches = descend(s, variables, &mut 4096)?;
    branches.sort_by_cached_key(|b| (b.values.iter().map(crate::sortkey::sort_key).collect::<Vec<_>>(), b.excluded.iter().map(|r| r.iter().map(crate::sortkey::sort_key).collect::<Vec<_>>()).collect::<Vec<_>>()));
    branches.dedup_by(|a,b| a.values == b.values && a.excluded == b.excluded);
    // Factor alternatives can produce intersections already contained in a
    // larger component. Remove them only with exact substitution certificates.
    let mut reduced: Vec<Branch> = Vec::new();
    for branch in branches {
        if reduced.iter().any(|general| contains_branch(general, &branch, variables)) { continue; }
        reduced.retain(|other| !contains_branch(&branch, other, variables));
        reduced.push(branch);
    }
    Ok(reduced)
}

fn contains_branch(general: &Branch, specific: &Branch, variables: &[Expr]) -> bool {
    let mapping: FxHashMap<_,_> = variables.iter().enumerate()
        .filter(|(i,x)| general.values[*i] == **x)
        .map(|(i,x)| (x.clone(),specific.values[i].clone())).collect();
    if general.values.iter().zip(&specific.values).any(|(a,b)| !zero(&normal(&sub(&crate::subs::xreplace(a,&mapping),b)))) { return false; }
    let restrictions: Vec<_> = specific.excluded.iter().enumerate().flat_map(|(i,rs)| rs.iter().map(move |r| normal(&sub(&variables[i],r)))).collect();
    general.excluded.iter().enumerate().all(|(i,rs)| rs.iter().all(|root| {
        let root = crate::subs::xreplace(root,&mapping);
        nonzero(&normal(&sub(&specific.values[i], &root)),&restrictions)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn expr(s: &str) -> Expr { crate::parser::parse_str(s).unwrap() }
    #[test]
    fn coefficient_zero_branches_and_pole_domains() {
        let vars = vec![expr("x"),expr("y"),expr("z")];
        let branches = nonlinsolve(&[expr("x*y-z")], &vars).unwrap();
        // Eliminating z directly covers every value of x and y.
        assert_eq!(branches.len(),1);
        assert_eq!(branches[0].values[2],expr("x*y"));
        let branches = nonlinsolve(&[expr("x*y")], &vars[..2]).unwrap();
        assert_eq!(branches.len(),2);
        let branches = nonlinsolve(&[expr("x*y-1")], &vars[..2]).unwrap();
        assert_eq!(branches.len(),1);
        assert_eq!(branches[0].excluded[1],vec![Expr::zero()]);
        assert!(nonlinsolve(&[expr("a*x")], &vars[..1]).is_err());
        assert!(nonlinsolve(&[expr("x**5-x+1")], &vars[..1]).is_err());
    }
    #[test]
    fn rational_exclusions_and_nongeneric_components() {
        let vars = vec![expr("x"),expr("y")];
        let branches = nonlinsolve(&[expr("x/(y-1)")], &vars).unwrap();
        assert_eq!(branches.len(),1);
        assert_eq!(branches[0].values[0],Expr::zero());
        assert_eq!(branches[0].excluded[1],vec![Expr::one()]);
        let branches = nonlinsolve(&[expr("x*y-x"),expr("x*y*y-x")], &vars).unwrap();
        assert_eq!(branches.len(),2);
    }
}
