//! Typed checked factor views. Scalar units are never raised accidentally to
//! the multiplicity of a factor when the unit is included in the output.
use super::algebraic::AlgebraicScope;
use super::convert::{PResult, PolyError};
use super::factor_budget::{self, Limits};
use super::fraction::Polynomial;
use super::int::{Coeff, Rat};
use super::mpoly::MPoly;
use super::polyobj::Domain;
use super::pseudo::{Budget, PseudoCoeff};
use super::ratfunc::Scope;
use crate::Expr;

#[derive(Clone, Debug)]
pub struct Factorization {
    pub unit: Polynomial,
    pub factors: Vec<(Polynomial, u32)>,
}
fn failure(message: &str) -> PolyError {
    PolyError::Unsupported(message.into())
}
fn check_rep<C: PseudoCoeff>(rep: &MPoly<C>, limits: Limits) -> PResult<()> {
    if rep.nvars > 8 {
        return Err(failure(
            "checked factorization supports at most eight polynomial/coefficient variables",
        ));
    }
    if rep.degrees().iter().any(|&d| d as usize > limits.degree)
        || rep.total_degree() > limits.degree as i64
    {
        return Err(failure("factorization resource limit: input degree"));
    }
    factor_budget::sparse(rep.len(), rep.nvars);
    let mut storage = 0u64;
    for c in rep.coeffs() {
        c.validate_size()?;
        storage = storage.saturating_add(c.budget_weight() + rep.nvars as u64 * 32);
    }
    factor_budget::storage(storage);
    Ok(())
}
fn preflight(p: &Polynomial, limits: Limits) -> PResult<()> {
    if p.ngens == 0 {
        return Err(failure("polynomial generators are required"));
    }
    if let Some(f) = &p.algebraic {
        if f.gens.len() != 1 {
            return Err(failure(
                "multivariate number-field factorization is not implemented",
            ));
        }
        check_rep(&f.rep, limits)?;
        if f.number_field
            .degree()
            .saturating_mul(f.rep.degree(0).max(0) as usize)
            > 256
        {
            return Err(failure(
                "factorization resource limit: algebraic norm degree exceeds 256",
            ));
        }
        for c in &f.number_field.minimal {
            factor_budget::bits(c.n.bits());
            factor_budget::bits(c.d.bits());
        }
        for c in f.rep.coeffs() {
            for value in &c.coefficients {
                factor_budget::bits(value.n.bits());
                factor_budget::bits(value.d.bits());
            }
        }
    } else if let Some(f) = &p.field {
        if f.gens.len() + f.parameters.len() > 8 {
            return Err(failure(
                "checked factorization supports at most eight polynomial/coefficient variables",
            ));
        }
        check_rep(&f.rep, limits)?;
        for c in f.rep.coeffs() {
            check_rep(&c.n, limits)?;
            check_rep(&c.d, limits)?;
            for c in c.n.coeffs().iter().chain(c.d.coeffs()) {
                factor_budget::bits(c.bits());
            }
        }
    } else {
        if let Some(prime) = p.modulus {
            super::f5b::PrimeScope::new(prime)?;
        }
        if p.domain == Domain::RR || !p.domain_gens().is_empty() {
            return Err(failure(
                "checked factor views require an implemented exact coefficient domain",
            ));
        }
        if p.modulus.is_some() && p.ngens != 1 {
            return Err(failure(
                "multivariate finite-field factorization is not implemented",
            ));
        }
        check_rep(&p.rep, limits)?;
        for c in p.rep.coeffs() {
            factor_budget::bits(c.n.bits());
            factor_budget::bits(c.d.bits());
        }
    }
    Ok(())
}
fn raw(p: &Polynomial, squarefree: bool) -> PResult<Factorization> {
    if let Some(f) = &p.algebraic {
        let _scope = AlgebraicScope::new(f.number_field.clone());
        let (unit, factors) = if squarefree {
            f.sqf_list()
        } else {
            f.factor_list()?
        };
        Ok(Factorization {
            unit: f.with_rep(MPoly::constant(unit, f.gens.len())).into(),
            factors: factors.into_iter().map(|(p, k)| (p.into(), k)).collect(),
        })
    } else if let Some(f) = &p.field {
        let _scope = Scope::new(f.parameters.len());
        let (unit, factors) = if squarefree && !f.display.ends_with(']') {
            f.sqf_list()
        } else {
            let (unit, factors) = f.factor_list();
            if squarefree {
                let mut groups =
                    std::collections::BTreeMap::<u32, super::fraction::FractionPoly>::new();
                let mut budget = Budget::new();
                for (factor, power) in factors {
                    if let Some(group) = groups.get_mut(&power) {
                        group.rep = budget.mul(&group.rep, &factor.rep)?;
                    } else {
                        groups.insert(power, factor);
                    }
                }
                (
                    unit,
                    groups
                        .into_iter()
                        .map(|(power, factor)| (factor, power))
                        .collect(),
                )
            } else {
                (unit, factors)
            }
        };
        Ok(Factorization {
            unit: f.with_rep(MPoly::constant(unit, f.gens.len())).into(),
            factors: factors.into_iter().map(|(p, k)| (p.into(), k)).collect(),
        })
    } else {
        let (unit, factors) = if squarefree {
            p.base.sqf_list()
        } else {
            p.base.factor_list()
        };
        Ok(Factorization {
            unit: p.base.with_rep(MPoly::constant(unit, p.gens.len())).into(),
            factors: factors.into_iter().map(|(p, k)| (p.into(), k)).collect(),
        })
    }
}
fn one(p: &Polynomial) -> Polynomial {
    if let Some(f) = &p.algebraic {
        let _scope = AlgebraicScope::new(f.number_field.clone());
        return f.with_rep(MPoly::one(f.gens.len())).into();
    }
    if let Some(f) = &p.field {
        let _scope = Scope::new(f.parameters.len());
        return f.with_rep(MPoly::one(f.gens.len())).into();
    }
    p.base.with_rep(MPoly::one(p.gens.len())).into()
}
fn is_one(p: &Polynomial) -> bool {
    p.algebraic.as_ref().map_or_else(
        || {
            p.field
                .as_ref()
                .map_or_else(|| p.rep.is_one(), |f| f.rep.is_one())
        },
        |f| f.rep.is_one(),
    )
}
fn multiply(a: &Polynomial, b: &Polynomial, budget: &mut Budget) -> PResult<Polynomial> {
    if let (Some(a), Some(b)) = (&a.algebraic, &b.algebraic) {
        let _scope = AlgebraicScope::new(a.number_field.clone());
        return Ok(a.with_rep(budget.mul(&a.rep, &b.rep)?).into());
    }
    if let (Some(a), Some(b)) = (&a.field, &b.field) {
        let _scope = Scope::new(a.parameters.len());
        return Ok(a.with_rep(budget.mul(&a.rep, &b.rep)?).into());
    }
    Ok(a.base.with_rep(budget.mul(&a.rep, &b.rep)?).into())
}
fn power(p: &Polynomial, mut exponent: u32, budget: &mut Budget) -> PResult<Polynomial> {
    let mut out = one(p);
    let mut base = p.clone();
    while exponent != 0 {
        if exponent & 1 != 0 {
            out = multiply(&out, &base, budget)?;
        }
        exponent >>= 1;
        if exponent != 0 {
            base = multiply(&base, &base, budget)?;
        }
    }
    Ok(out)
}
fn same(a: &Polynomial, b: &Polynomial) -> bool {
    if let (Some(a), Some(b)) = (&a.algebraic, &b.algebraic) {
        return a.rep == b.rep;
    }
    if let (Some(a), Some(b)) = (&a.field, &b.field) {
        return a.rep == b.rep;
    }
    a.rep == b.rep
}
fn verify(p: &Polynomial, factorization: &Factorization) -> PResult<()> {
    let mut budget = Budget::new();
    let mut result = factorization.unit.clone();
    for (factor, multiplicity) in &factorization.factors {
        if *multiplicity == 0 {
            return Err(failure("invalid factor multiplicity"));
        }
        let factor = power(factor, *multiplicity, &mut budget)?;
        result = multiply(&result, &factor, &mut budget)?;
    }
    if !same(&result, p) {
        return Err(failure(
            "factorization failed its exact reconstruction certificate",
        ));
    }
    if let Some(p) = &p.field {
        if p.display.ends_with(']') {
            for factor in std::iter::once(&factorization.unit)
                .chain(factorization.factors.iter().map(|(p, _)| p))
            {
                if !factor.field.as_ref().is_some_and(|f| f.belongs_to_ring()) {
                    return Err(failure("factorization escaped its coefficient ring"));
                }
            }
        }
    }
    Ok(())
}
pub fn factor_list(p: &Polynomial, squarefree: bool, limits: Limits) -> PResult<Factorization> {
    factor_budget::run(limits, || {
        preflight(p, limits)?;
        let result = raw(p, squarefree)?;
        verify(p, &result)?;
        Ok(result)
    })
}
pub fn include(
    p: &Polynomial,
    squarefree: bool,
    all: bool,
    limits: Limits,
) -> PResult<Vec<(Polynomial, u32)>> {
    factor_budget::run(limits, || {
        preflight(p, limits)?;
        let factorization = raw(p, squarefree)?;
        verify(p, &factorization)?;
        let Factorization { unit, mut factors } = factorization;
        if factors.is_empty() {
            return Ok(vec![(unit, 1)]);
        }
        if squarefree && all {
            let highest = factors.iter().map(|(_, k)| *k).max().unwrap_or(1);
            if highest as usize > limits.degree {
                return Err(failure(
                    "factorization resource limit: square-free multiplicity slots",
                ));
            }
            let mut filled = Vec::new();
            for exponent in 1..=highest {
                if let Some((factor, _)) = factors.iter().find(|(_, k)| *k == exponent) {
                    filled.push((factor.clone(), exponent));
                } else {
                    filled.push((one(p), exponent));
                }
            }
            factors = filled;
        }
        if !is_one(&unit) {
            if let Some((factor, _)) = factors.iter_mut().find(|(_, k)| *k == 1) {
                *factor = multiply(&unit, factor, &mut Budget::new())?;
            } else {
                factors.insert(0, (unit, 1));
            }
        }
        // Unlike SymPy 1.14 factor_list_include, a scalar is never blindly
        // absorbed into a repeated factor: that would multiply its exponent.
        verify(
            p,
            &Factorization {
                unit: one(p),
                factors: factors.clone(),
            },
        )?;
        Ok(factors)
    })
}
fn roots<C: PseudoCoeff>(
    original: &MPoly<C>,
    factors: impl Iterator<Item = (MPoly<C>, u32)>,
) -> PResult<Vec<(C, u32)>> {
    let mut roots = Vec::<(C, u32)>::new();
    for (factor, multiplicity) in factors {
        if factor.degree(0) != 1 {
            continue;
        }
        let root = factor.constant_term().neg().checked_div(&factor.lc())?;
        let mut value = C::zero();
        let mut previous = original.degree(0).max(0) as u32;
        for (e, c) in original.iter() {
            let mut n = previous - e[0];
            let mut base = root.clone();
            let mut power = C::one();
            while n != 0 {
                factor_budget::charge(1);
                if n & 1 != 0 {
                    power = power.checked_mul(&base)?;
                }
                n >>= 1;
                if n != 0 {
                    base = base.checked_mul(&base)?;
                }
            }
            value = value.checked_mul(&power)?.checked_add(c)?;
            previous = e[0];
        }
        for _ in 0..previous {
            factor_budget::charge(1);
            value = value.checked_mul(&root)?;
        }
        if !value.is_zero() {
            return Err(failure("ground root failed exact substitution certificate"));
        }
        if let Some((_, power)) = roots.iter_mut().find(|(r, _)| *r == root) {
            *power = power
                .checked_add(multiplicity)
                .ok_or_else(|| failure("root multiplicity overflow"))?;
        } else {
            roots.push((root, multiplicity));
        }
    }
    Ok(roots)
}
pub fn ground_roots(p: &Polynomial, limits: Limits) -> PResult<Vec<(Expr, u32)>> {
    factor_budget::run(limits, || {
        if p.ngens != 1 {
            return Err(failure("ground_roots requires one polynomial generator"));
        }
        preflight(p, limits)?;
        let factorization = raw(p, false)?;
        verify(p, &factorization)?;
        let mut result = if let Some(f) = &p.algebraic {
            let _scope = AlgebraicScope::new(f.number_field.clone());
            roots(
                &f.rep,
                factorization
                    .factors
                    .into_iter()
                    .map(|(p, k)| (p.algebraic.unwrap().rep, k)),
            )?
            .into_iter()
            .map(|(r, k)| (r.as_expr(), k))
            .collect::<Vec<_>>()
        } else if let Some(f) = &p.field {
            let _scope = Scope::new(f.parameters.len());
            roots(
                &f.rep,
                factorization
                    .factors
                    .into_iter()
                    .map(|(p, k)| (p.field.unwrap().rep, k)),
            )?
            .into_iter()
            .map(|(r, k)| (super::fraction::coefficient_expr(&r, &f.parameters), k))
            .collect()
        } else if let Some(prime) = p.modulus {
            let _scope = super::f5b::PrimeScope::new(prime)?;
            let convert = |p: &MPoly<Rat>| {
                p.map_coeffs(|c| {
                    super::f5b::Fp::from_rat(c, prime).expect("finite polynomial coefficients")
                })
            };
            roots(
                &convert(&p.rep),
                factorization
                    .factors
                    .into_iter()
                    .map(|(p, k)| (convert(&p.rep), k)),
            )?
            .into_iter()
            .map(|(r, k)| (Expr::num(r.to_rat().to_number()), k))
            .collect()
        } else {
            roots(
                &p.rep,
                factorization
                    .factors
                    .into_iter()
                    .map(|(p, k)| (p.base.rep, k)),
            )?
            .into_iter()
            .map(|(r, k)| (Expr::num(r.to_number()), k))
            .collect()
        };
        result.sort_by_cached_key(|(root, _)| crate::sortkey::sort_key(root));
        Ok(result)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    use crate::poly::{
        algebraic::{AlgebraicPoly, NumberField},
        fraction::FractionPoly,
        polyobj::Poly,
    };
    fn expr(source: &str) -> Expr {
        parse_str(source).unwrap()
    }
    fn p(source: &str) -> Polynomial {
        Poly::from_expr(&expr(source), &[Expr::symbol("x")])
            .unwrap()
            .into()
    }
    #[test]
    fn exact_views_keep_units_multiplicities_and_associated_field_roots() {
        let f = p("2*(x+1)^2");
        let output = include(&f, false, false, Limits::default()).unwrap();
        assert_eq!(output.len(), 2);
        assert!(output.iter().any(|(factor, k)| factor.rep == p("2").rep && *k == 1));
        let roots = ground_roots(&p("(2*x-1)^3*(x^2-2)"), Limits::default()).unwrap();
        assert_eq!(roots, vec![(Expr::rational(1, 2), 3)]);
        let f = p("(x+1)^3");
        let output = include(&f, true, true, Limits::default()).unwrap();
        assert_eq!(
            output.iter().map(|(_, k)| *k).collect::<Vec<_>>(),
            vec![1, 2, 3]
        );
        assert!(is_one(&output[0].0) && is_one(&output[1].0));
    }
    #[test]
    fn rational_function_rings_and_algebraic_scopes_preserve_exact_values() {
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        let b = Expr::symbol("unused");
        for domain in ["ZZ[a,unused]", "QQ[a,unused]", "QQ(a,unused)"] {
            let f: Polynomial = FractionPoly::from_expr(
                &expr("(a*x-1)^2*(x-2)"),
                std::slice::from_ref(&x),
                &[a.clone(), b.clone()],
                domain,
            )
            .unwrap()
            .into();
            let roots = ground_roots(&f, Limits::default()).unwrap();
            assert!(roots.contains(&(expr("1/a"), 2)));
            assert!(roots.contains(&(Expr::int(2), 1)));
            for squarefree in [false, true] {
                for (q, _) in include(&f, squarefree, true, Limits::default()).unwrap() {
                    assert_eq!(q.domain_str(), domain);
                    assert_eq!(q.field.unwrap().parameters, [a.clone(), b.clone()]);
                }
            }
        }
        let field = NumberField::from_generator(&expr("sqrt(2)")).unwrap();
        let f: Polynomial = AlgebraicPoly::from_expr(&expr("(x^2-2)*(x-1)^2"), &[x], field.clone())
            .unwrap()
            .into();
        let roots = ground_roots(&f, Limits::default()).unwrap();
        assert_eq!(roots.len(), 3);
        assert!(roots.contains(&(expr("sqrt(2)"), 1)));
        let error = ground_roots(
            &f,
            Limits {
                attempts: 0,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("algebraic norm shifts"));
        assert!(!factor_budget::active());
    }
    #[test]
    fn finite_characteristic_and_preflight_failures_are_exact() {
        let f: Polynomial = p("x^5-x").base.set_modulus(5).into();
        let roots = ground_roots(&f, Limits::default()).unwrap();
        assert_eq!(roots.len(), 5);
        let output = include(&f, false, false, Limits::default()).unwrap();
        assert_eq!(output.len(), 5);
        let f = p("x^1000000000+1");
        assert!(ground_roots(&f, Limits::default())
            .unwrap_err()
            .to_string()
            .contains("resource limit"));
        assert!(!factor_budget::active());
        assert!(ground_roots(&p("0"), Limits::default()).unwrap().is_empty());
    }
}
