//! Finite polynomial systems: lex elimination and complete radical roots.

use super::convert::{from_exprs_gens, to_expr, PResult, PolyError};
use super::groebner::{buchberger, is_zero_dimensional, CoefficientField, MonomialOrder};
use crate::expand::expand;
use crate::expr::{Expr, Kind};
use crate::integrate::poly_coeffs;
use crate::subs::subs1;
use crate::traverse::has;

fn unsupported(message: &str) -> PolyError {
    PolyError::Unsupported(message.into())
}

fn normalize(e: &Expr) -> Expr {
    crate::simplify::simplify(&expand(e))
}

fn scalar_zero_status(e: &Expr) -> Option<bool> {
    if e.is_zero() {
        return Some(true);
    }
    crate::assumptions::exact_queries(|| crate::assumptions::is_zero(e))
        .or_else(|| super::radical_cert::is_zero(e))
}

fn complete_roots(e: &Expr, x: &Expr) -> PResult<Vec<Expr>> {
    let coefficients = poly_coeffs(e, x).ok_or_else(|| {
        unsupported("root back substitution produced a non-polynomial expression")
    })?;
    let degree = coefficients
        .iter()
        .rposition(|c| !normalize(c).is_zero())
        .unwrap_or(0);
    let roots = super::roots::roots(e, std::slice::from_ref(x))?;
    if roots.iter().map(|(_, m)| *m as usize).sum::<usize>() != degree {
        return Err(unsupported("a complete solution set cannot be expressed by the implemented radical root algorithms"));
    }
    if coefficients.iter().all(Expr::is_rational_number) {
        // SymPy's solve_generic keeps the raw `roots` output of a univariate
        // polynomial over QQ (`roots(f.ltrim(gen))`), neither simplified nor
        // re-verified: the exact radical formulas are complete there. A
        // simplify of the substituted residual of cubic or quartic radicals
        // can take unbounded time in the algebraic multivariate gcd.
        let mut values: Vec<_> = roots.into_iter().map(|(root, _)| root).collect();
        values.sort_by_cached_key(crate::sortkey::sort_key);
        values.dedup();
        return Ok(values);
    }
    // The root -c0/c1 of a linear polynomial with a certified nonzero
    // leading coefficient is exact by construction (SymPy's roots_linear),
    // so its residual needs no zero certificate.
    let linear = degree == 1 && scalar_zero_status(&normalize(&coefficients[1])) == Some(false);
    let mut values = Vec::new();
    for (root, _) in roots {
        // Expansion can turn related radicals into independent-looking
        // powers. Preserve the original formula for algebraic certificates.
        let residual = subs1(e, x, &root);
        if !linear
            && super::radical_cert::is_zero(&residual) != Some(true)
            && scalar_zero_status(&normalize(&residual)) != Some(true) {
            return Err(unsupported("an exact zero certificate for a radical root is not implemented"));
        }
        values.push(normalize(&root));
    }
    values.sort_by_cached_key(crate::sortkey::sort_key);
    values.dedup();
    Ok(values)
}

fn descend(
    equations: &[Expr],
    variables: &[Expr],
    remaining: &[usize],
    assignment: &[Option<Expr>],
    lex: bool,
) -> PResult<Vec<Vec<Expr>>> {
    let equations: Vec<_> = equations
        .iter()
        .map(normalize)
        .filter(|e| !e.is_zero())
        .collect();
    for e in &equations {
        if remaining.iter().all(|&j| !has(e, &variables[j])) {
            match scalar_zero_status(e) {
                Some(true) => {},
                Some(false) => return Ok(Vec::new()),
                None => return Err(unsupported("an exact zero certificate for a radical residual or external parameter condition is not implemented")),
            }
        }
    }
    let equations: Vec<_> = equations
        .into_iter()
        .filter(|e| remaining.iter().any(|&j| has(e, &variables[j])))
        .collect();
    if remaining.is_empty() {
        return Ok(vec![assignment
            .iter()
            .map(|e| e.clone().expect("assigned variable"))
            .collect()]);
    }
    if equations.is_empty() {
        return Err(unsupported(
            "only zero-dimensional polynomial systems are supported",
        ));
    }
    let candidates: Vec<_> = if lex {
        vec![*remaining.last().unwrap()]
    } else {
        remaining.to_vec()
    };
    let mut best: Option<(usize, usize, usize)> = None;
    for &j in &candidates {
        for (i, e) in equations.iter().enumerate() {
            let Some(c) = poly_coeffs(e, &variables[j]) else {
                continue;
            };
            let degree = c.iter().rposition(|a| !normalize(a).is_zero()).unwrap_or(0);
            if degree == 0 {
                continue;
            }
            let independent = remaining.iter().all(|&k| k == j || !has(e, &variables[k]));
            let linear_safe =
                !lex && degree == 1 && remaining.iter().all(|&k| !has(&c[1], &variables[k]));
            if !(independent || linear_safe) {
                continue;
            }
            if best.is_none_or(|(_, _, d)| degree < d) {
                best = Some((i, j, degree));
            }
        }
    }
    let (i, j, _) = best.ok_or_else(|| {
        unsupported(
            "coefficient-field elimination is not implemented for this symbolic parameter system",
        )
    })?;
    let roots = complete_roots(&equations[i], &variables[j])?;
    let next: Vec<_> = remaining.iter().copied().filter(|&k| k != j).collect();
    let mut result = Vec::new();
    for root in roots {
        let rest: Vec<_> = equations
            .iter()
            .enumerate()
            .filter(|(k, _)| *k != i)
            .map(|(_, e)| normalize(&subs1(e, &variables[j], &root)))
            .collect();
        let mut values = assignment.to_vec();
        for value in values.iter_mut().flatten() {
            if has(value, &variables[j]) {
                *value = normalize(&subs1(value, &variables[j], &root));
            }
        }
        values[j] = Some(root);
        result.extend(descend(&rest, variables, &next, &values, lex)?);
    }
    Ok(result)
}

/// Solve a polynomial system completely where the native radical algorithms
/// apply. `None` means an inconsistent system, never an unsupported one.
///
/// Rational systems use a lexicographic Groebner basis. Systems with symbolic
/// coefficients support exact linear elimination and triangular univariate
/// equations; more general coefficient-field elimination is explicit.
pub fn solve_poly_system(
    equations: &[Expr],
    variables: &[Expr],
) -> PResult<Option<Vec<Vec<Expr>>>> {
    if variables.is_empty() {
        return Err(unsupported("polynomial system solving requires generators"));
    }
    for (i, x) in variables.iter().enumerate() {
        if !matches!(x.kind(), Kind::Sym(_)) || variables[..i].contains(x) {
            return Err(PolyError::Failed(
                "generators must be distinct symbols".into(),
            ));
        }
    }
    let equations: Vec<_> = equations
        .iter()
        .map(|e| match e.kind() {
            Kind::Bool(true) => Expr::zero(),
            Kind::Bool(false) => Expr::one(),
            _ => crate::solve::as_expr(e),
        })
        .collect();
    let converted = from_exprs_gens(&equations, variables)?;
    if converted.float {
        return Err(unsupported(
            "polynomial system solving requires exact coefficients",
        ));
    }
    let lex = converted.gens.len() == converted.ngens;
    let elimination = if lex {
        let basis = buchberger(
            &converted.polys,
            MonomialOrder::Lex,
            CoefficientField { modulus: None },
        );
        if basis.iter().any(|p| p.is_constant() && !p.is_zero()) {
            return Ok(None);
        }
        if !is_zero_dimensional(&basis, MonomialOrder::Lex, variables.len()) {
            return Err(unsupported(
                "only zero-dimensional polynomial systems are supported",
            ));
        }
        basis.iter().map(|p| to_expr(p, variables)).collect()
    } else {
        equations
    };
    let remaining: Vec<_> = (0..variables.len()).collect();
    let mut solutions = descend(
        &elimination,
        variables,
        &remaining,
        &vec![None; variables.len()],
        lex,
    )?;
    solutions.sort_by_cached_key(|t| t.iter().map(crate::sortkey::sort_key).collect::<Vec<_>>());
    solutions.dedup();
    Ok(if solutions.is_empty() {
        None
    } else {
        Some(solutions)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    #[test]
    fn exact_systems_and_substitution_certificates() {
        let variables = vec![Expr::symbol("x"), Expr::symbol("y")];
        for (input, count) in [
            (vec!["x**2+y**2-5", "x-y-1"], 2),
            (vec!["x**2-1", "y**2-1"], 4),
            (vec!["x**2-a", "y-x"], 2),
        ] {
            let equations: Vec<_> = input.iter().map(|s| parse_str(s).unwrap()).collect();
            let solutions = solve_poly_system(&equations, &variables).unwrap().unwrap();
            assert_eq!(solutions.len(), count);
            for solution in solutions {
                for equation in &equations {
                    let value = variables
                        .iter()
                        .zip(&solution)
                        .fold(equation.clone(), |e, (x, v)| subs1(&e, x, v));
                    assert_eq!(scalar_zero_status(&normalize(&value)), Some(true));
                }
            }
        }
        assert!(solve_poly_system(&[Expr::one()], &variables)
            .unwrap()
            .is_none());
        assert!(solve_poly_system(&[variables[0].clone()], &variables).is_err());
    }
}
