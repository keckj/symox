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
/// Attach the coefficient-ring descriptor carried by a natural source Poly.
/// This conversion runs in the native factor budget, before Python rendering;
/// it never re-infers main generators or discards source parameters.
pub fn normalize_natural_factor(p:Polynomial)->PResult<Polynomial>{
    if p.field.is_some()||p.algebraic.is_some()||p.domain_gens().is_empty()
        ||p.modulus.is_some()||p.domain.is_rr()||p.domain_gens().iter().any(|g|!g.is_symbol()){
        return Ok(p);
    }
    let parameters=p.domain_gens().to_vec();
    let descriptor=super::domains::ExactDomain::functions(p.domain==Domain::ZZ,true,parameters.clone(),"lex")?;
    let mut field=super::fraction::FractionPoly::from_expr(&p.as_expr(),p.true_gens(),&parameters,&descriptor.to_string())?;
    field.coefficient_domain=descriptor;
    if !field.belongs_to_ring(){return Err(failure("natural factor escaped its coefficient ring"));}
    Ok(field.into())
}
/// Variables (generators plus symbolic coefficients) of a checked
/// factorization; SymPy has no such limit, this only bounds the dense
/// exponent vectors.
const MAX_VARIABLES: usize = 64;
fn failure(message: &str) -> PolyError {
    PolyError::Unsupported(message.into())
}
fn check_rep<C: PseudoCoeff>(rep: &MPoly<C>, limits: Limits) -> PResult<()> {
    if rep.nvars > MAX_VARIABLES {
        return Err(failure(
            "factorization resource limit: more than 64 polynomial/coefficient variables",
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
    preflight_mode(p, limits, false)
}
fn preflight_mode(p: &Polynomial, limits: Limits, legacy: bool) -> PResult<()> {
    if p.ngens == 0 && !p.rep.is_constant() {
        return Err(failure("polynomial generators are required"));
    }
    if let Some(f) = &p.algebraic {
        if f.gens.len() > MAX_VARIABLES {
            return Err(failure(
                "factorization resource limit: more than 64 polynomial/coefficient variables",
            ));
        }
        check_rep(&f.rep, limits)?;
        let degree = (0..f.gens.len())
            .map(|v| f.rep.degree(v).max(0) as usize)
            .fold(0usize, usize::saturating_add);
        if f.number_field.degree().saturating_mul(degree) > 256
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
        if f.gens.len() + f.parameters.len() > MAX_VARIABLES {
            return Err(failure(
                "factorization resource limit: more than 64 polynomial/coefficient variables",
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
        if !legacy && (p.domain.is_rr() || !p.domain_gens().is_empty()) {
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
        let (mut unit, factors) = if squarefree {
            f.sqf_list()
        } else {
            f.factor_list()?
        };
        // The reference `dmp_sqf_list` over ZZ_I loses the unit carried by
        // the normalised gcds (`(I*x + 1)**2` gives `(1, [(-x + I, 2)])`,
        // which is `-(I*x + 1)**2`). Keep the reference factors and take the
        // unit that makes the product exact: an intentional correction.
        if squarefree && f.gaussian.is_some() && !f.rep.is_zero() {
            use super::int::Field;
            let mut budget = Budget::new();
            let mut product = MPoly::one(f.gens.len());
            for (g, k) in &factors {
                for _ in 0..*k {
                    product = budget.mul(&product, &g.rep)?;
                }
            }
            if !product.is_zero() {
                let c = f.rep.lc().mul(&product.lc().inv());
                if c != unit && product.scale(&c) == f.rep {
                    unit = c;
                }
            }
        }
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
        let mut result = Factorization {
            unit: p.base.with_rep(MPoly::constant(unit, p.gens.len())).into(),
            factors: factors.into_iter().map(|(p, k)| (p.into(), k)).collect(),
        };
        // SymPy preserves integer-primitive factors over QQ for factor_list,
        // but makes square-free groups monic. These are distinct conventions.
        if squarefree && p.domain == Domain::QQ && p.modulus.is_none() && p.domain_gens().is_empty()
        {
            let mut budget = Budget::new();
            for (f, k) in &mut result.factors {
                let scalar = f
                    .base
                    .with_rep(MPoly::constant(f.rep.lc(), f.gens.len()))
                    .into();
                let lc = f.rep.lc();
                f.base.rep = budget.exact_div(&f.rep, &MPoly::constant(lc, f.gens.len()))?;
                result.unit =
                    multiply(&result.unit, &power(&scalar, *k, &mut budget)?, &mut budget)?;
            }
        }
        if squarefree {
            result.factors.sort_by_key(|(_, k)| *k);
        }
        Ok(result)
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
    if a.domain.is_rr() || b.domain.is_rr() {
        // RR factors come from the bounded-denominator rational image of
        // the coefficients (RR.to_rational) and are rounded back to floats,
        // so the reconstruction is checked to a relative tolerance.
        let (x, y) = (a.rep.iter().collect::<Vec<_>>(), b.rep.iter().collect::<Vec<_>>());
        let close = |u: &Rat, v: &Rat| {
            let (u, v) = (u.to_f64(), v.to_f64());
            (u - v).abs() <= 1e-12 * u.abs().max(v.abs())
        };
        return x.len() == y.len()
            && x.iter().zip(y.iter()).all(|((m, u), (n, v))| m == n && close(u, v));
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
/// Existing opaque/inexact coefficient polynomials retain their exact stored
/// arithmetic. They do not acquire a fictitious typed Domain; their factors
/// are checked in the existing expanded generator representation instead.
/// The reference views over `EX`: square-free decomposition runs in the
/// domain, while `dmp_factor_list` raises `DomainError`, so `Poly.factor_list`
/// returns a constant as the unit and anything else as its only factor.
fn expression_factorization(
    p: &super::expression_domain::ExpressionPoly,
    squarefree: bool,
    all: bool,
) -> PResult<Factorization> {
    let constant = |c: Expr| {
        super::expression_domain::ExpressionPoly::from_terms(vec![(vec![0; p.gens.len()], c)], &p.gens)
    };
    if squarefree {
        let (c, factors) = p.sqf_list(all)?;
        return Ok(Factorization {
            unit: constant(c)?.into(),
            factors: factors.into_iter().map(|(f, k)| (f.into(), k as u32)).collect(),
        });
    }
    if p.degree(0) == 0 {
        return Ok(Factorization { unit: p.clone().into(), factors: Vec::new() });
    }
    Ok(Factorization { unit: constant(Expr::one())?.into(), factors: vec![(p.clone().into(), 1)] })
}
pub fn compatible_list(
    p: &Polynomial,
    squarefree: bool,
    all: bool,
    limits: Limits,
) -> PResult<Factorization> {
    if let Some(e) = &p.expression {
        return expression_factorization(e, squarefree, all);
    }
    factor_budget::run(limits, || {
        let legacy = p.field.is_none()
            && p.algebraic.is_none()
            && (p.domain.is_rr() || !p.domain_gens().is_empty());
        preflight_mode(p, limits, legacy)?;
        let mut result = raw(p, if legacy && !p.domain.is_rr() { false } else { squarefree })?;
        if legacy {
            let mut kept = Vec::new();
            let mut budget = Budget::new();
            for (factor, k) in result.factors {
                if factor
                    .rep
                    .iter()
                    .all(|(e, _)| e[..p.ngens].iter().all(|&e| e == 0))
                {
                    result.unit =
                        multiply(&result.unit, &power(&factor, k, &mut budget)?, &mut budget)?;
                } else {
                    kept.push((factor, k));
                }
            }
            result.factors = kept;
            if squarefree {
                let mut groups = std::collections::BTreeMap::<u32, Polynomial>::new();
                for (factor, k) in result.factors {
                    if let Some(group) = groups.get_mut(&k) {
                        *group = multiply(group, &factor, &mut budget)?;
                    } else {
                        groups.insert(k, factor);
                    }
                }
                result.factors = groups.into_iter().map(|(k, p)| (p, k)).collect();
            }
        }
        if squarefree {
            result.factors.sort_by_key(|(_, k)| *k);
        }
        if all && !result.factors.is_empty() {
            let last = result.factors.iter().map(|(_, k)| *k).max().unwrap_or(0);
            if last as usize > limits.degree {
                return Err(failure(
                    "factorization resource limit: square-free multiplicity slots",
                ));
            }
            let mut factors = result.factors.into_iter().peekable();
            result.factors = (1..=last)
                .map(|k| {
                    if factors.peek().is_some_and(|(_, j)| *j == k) {
                        factors.next().unwrap()
                    } else {
                        (one(p), k)
                    }
                })
                .collect();
        }
        verify(p, &result)?;
        Ok(result)
    })
}
pub fn squarefree_part(p: &Polynomial, limits: Limits) -> PResult<Polynomial> {
    if let Some(e) = &p.expression {
        return Ok(e.sqf_part()?.into());
    }
    factor_budget::run(limits, || {
        let result = compatible_list(p, true, false, limits)?;
        if p.is_zero() {
            return Ok(p.clone());
        }
        let mut product = one(p);
        let mut budget = Budget::new();
        for (f, _) in result.factors {
            product = multiply(&product, &f, &mut budget)?;
        }
        Ok(product)
    })
}
/// Symbolic list factorization first separates a coefficient denominator.
/// This is deliberately distinct from factoring a Poly over QQ/QQ(a): the
/// factors of an expression with no domain selection are inferred from its
/// integral numerator, while the coefficient can remain outside that ring.
pub fn symbolic_numerator(p: &Polynomial, limits: Limits) -> PResult<(Polynomial, Expr)> {
    factor_budget::run(limits, || {
        let legacy = p.field.is_none()
            && p.algebraic.is_none()
            && (p.domain.is_rr() || !p.domain_gens().is_empty());
        preflight_mode(p, limits, legacy)?;
        if let Some(f) = &p.field {
            let _scope = Scope::new(f.parameters.len());
            let denominator = f.denominator();
            let mut budget = Budget::new();
            let mut terms = Vec::new();
            for (e, c) in f.rep.iter() {
                let quotient = budget.exact_div(&denominator, &c.d)?;
                let numerator = budget.mul(&c.n, &quotient)?;
                terms.push((e.to_vec(), super::ratfunc::RF::from_poly(numerator)));
            }
            let mut numerator = f.with_rep(MPoly::from_terms(f.gens.len(), terms));
            let super::domains::ExactDomain::Functions { order, .. } = &f.coefficient_domain else {
                return Err(failure("unsupported symbolic coefficient descriptor"));
            };
            numerator.coefficient_domain =
                super::domains::ExactDomain::functions(true, true, f.parameters.clone(), order)?;
            numerator.display = numerator.coefficient_domain.to_string();
            let coefficient = super::ratfunc::RF::from_poly(denominator);
            let expected =
                budget.mul(&f.rep, &MPoly::constant(coefficient.clone(), f.gens.len()))?;
            if expected != numerator.rep {
                return Err(failure(
                    "symbolic numerator failed exact denominator-clearing certificate",
                ));
            }
            let denominator = super::fraction::coefficient_expr(&coefficient, &f.parameters);
            return Ok((numerator.into(), denominator));
        }
        if p.algebraic.is_some() || p.modulus.is_some() || legacy {
            return Ok((p.clone(), Expr::one()));
        }
        let (denominator, mut numerator) = p.base.clear_denoms();
        numerator.domain = Domain::ZZ;
        let denominator = Expr::num(denominator.to_number());
        Ok((numerator.into(), denominator))
    })
}
pub fn symbolic_list(
    p: &Polynomial,
    squarefree: bool,
    limits: Limits,
) -> PResult<(Factorization, Expr)> {
    factor_budget::run(limits, || {
        let (numerator, denominator) = symbolic_numerator(p, limits)?;
        let mut result = compatible_list(&numerator, squarefree, false, limits)?;
        // No domain was selected for this symbolic expression. Infer the
        // coefficient parameters of each computed factor; parameters which
        // occurred only in the external unit must not leak into that factor.
        for (factor, _) in &mut result.factors {
            if let Some(f) = &factor.field {
                let used: Vec<_> = (0..f.parameters.len())
                    .filter(|&i| {
                        f.rep
                            .coeffs()
                            .iter()
                            .any(|c| c.n.degree(i) > 0 || c.d.degree(i) > 0)
                    })
                    .collect();
                if used.is_empty() {
                    *factor = factor.base.clone().into();
                    continue;
                }
                if used.len() == f.parameters.len() {
                    continue;
                }
                let mut permutation = vec![0; f.parameters.len()];
                for (j, &i) in used.iter().enumerate() {
                    permutation[i] = j;
                }
                let _scope = Scope::new(used.len());
                let rep = f.rep.map_coeffs(|c| {
                    super::ratfunc::RF::new(
                        c.n.permute(&permutation, used.len()),
                        c.d.permute(&permutation, used.len()),
                    )
                });
                let mut f = f.clone();
                f.parameters = used.iter().map(|&i| f.parameters[i].clone()).collect();
                let super::domains::ExactDomain::Functions {
                    integral,
                    ring,
                    order,
                    ..
                } = &f.coefficient_domain
                else {
                    return Err(failure("unsupported symbolic factor domain"));
                };
                f.coefficient_domain = super::domains::ExactDomain::functions(
                    *integral,
                    *ring,
                    f.parameters.clone(),
                    order,
                )?;
                f.display = f.coefficient_domain.to_string();
                f.rep = rep;
                *factor = f.into();
            }
        }
        // Symbolic sorting uses the *inferred* domains of the final factors,
        // not the common parameter ring in which the numerator was factored.
        result.factors.sort_by(|(a, ka), (b, kb)| {
            let degree = a.rep.degree(0).cmp(&b.rep.degree(0));
            let generators = a.ngens.cmp(&b.ngens);
            let main = if squarefree {
                ka.cmp(kb).then(degree).then(generators)
            } else {
                degree.then(generators).then(ka.cmp(kb))
            };
            main.then_with(|| a.domain_str().cmp(&b.domain_str()))
                .then_with(|| {
                    if a.rep.nvars == b.rep.nvars {
                        super::polyobj::rep_cmp(&a.rep, &b.rep, 0)
                    } else {
                        crate::sortkey::sort_key(&a.as_expr())
                            .cmp(&crate::sortkey::sort_key(&b.as_expr()))
                    }
                })
        });
        Ok((result, denominator))
    })
}
/// Reconstruct heterogeneous source-component factors in one common stored
/// polynomial representation. Individual factor generator vectors remain
/// untouched; the common prototype is used only for the exact certificate.
pub fn verify_expression_factors(prototype:&Polynomial,coefficient:&Expr,factors:&[(Polynomial,u32)],limits:Limits)->PResult<()> {
    factor_budget::run(limits,||{
        let convert=|e:&Expr|->PResult<Polynomial>{
            if let Some(f)=&prototype.field{return f.convert(e).map(Into::into);}
            if let Some(f)=&prototype.algebraic{return f.convert(e).map(Into::into);}
            let p=super::polyobj::Poly::from_expr(e,&prototype.base.gens)?;
            Ok(prototype.base.with_rep(p.rep).into())
        };
        let typed=(||{
            let unit=convert(coefficient)?;
            let factors=factors.iter().map(|(p,k)|Ok((convert(&p.as_expr())?,*k))).collect::<PResult<_>>()?;
            verify(prototype,&Factorization{unit,factors})
        })();
        match typed {
            Ok(())=>Ok(()),
            // The components are polynomials in their own generators. When
            // these are algebraic numbers (`sqrt(2)`) or radicals of other
            // generators, the typed product treats them as independent
            // indeterminates although the expanded source has already used
            // their relations (`sqrt(2)**2 = 2`): `(x - sqrt(2))*(x + sqrt(2))`
            // is `x**2 - 2` in the generators `x`, but `x**2 - s**2` in
            // `x, s`. Certify those factorizations on the canonical expanded
            // expressions instead, where the relations are applied.
            Err(error@(PolyError::Unsupported(_)|PolyError::Failed(_)|PolyError::Generator(_)))
                if !error.to_string().contains("resource limit")
                    && expression_certificate(prototype,coefficient,factors)? => Ok(()),
            Err(error)=>Err(error),
        }
    })
}
/// Whether `coefficient*prod(f**k)` and the prototype expand to the same
/// canonical expression; a structural zero difference is an exact proof.
fn expression_certificate(prototype:&Polynomial,coefficient:&Expr,factors:&[(Polynomial,u32)])->PResult<bool>{
    let mut product=vec![coefficient.clone()];
    for (p,k) in factors {
        if *k==0 { return Ok(false); }
        product.push(crate::ops::pow(&p.as_expr(),&Expr::int(*k as i64)));
    }
    let product=crate::expand::expand(&crate::ops::mul(&product));
    let source=crate::expand::expand(&prototype.as_expr());
    Ok(crate::expand::expand(&crate::ops::sub(&product,&source)).is_zero())
}
pub fn include(
    p: &Polynomial,
    squarefree: bool,
    all: bool,
    limits: Limits,
) -> PResult<Vec<(Polynomial, u32)>> {
    if let Some(e) = &p.expression {
        if !squarefree {
            return Ok(vec![(p.clone(), 1)]);
        }
        return Ok(e.sqf_list_include(all)?.into_iter().map(|(f, k)| (f.into(), k as u32)).collect());
    }
    factor_budget::run(limits, || {
        // Polynomial rings over ZZ[pi], ZZ[log(a)] or RR take the list
        // factorization of compatible_list, which keeps the content in the
        // unit as dmp_sqf_list does, and has already filled `all`.
        let legacy = p.field.is_none()
            && p.algebraic.is_none()
            && (p.domain.is_rr() || !p.domain_gens().is_empty());
        let factorization = if legacy {
            compatible_list(p, squarefree, all, limits)?
        } else {
            preflight(p, limits)?;
            let factorization = raw(p, squarefree)?;
            verify(p, &factorization)?;
            factorization
        };
        let Factorization { unit, mut factors } = factorization;
        if factors.is_empty() {
            return Ok(vec![(unit, 1)]);
        }
        if squarefree && all && !legacy {
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
        if let Some(e) = &p.expression {
            if e.degree(0) != 1 {
                return Ok(Vec::new());
            }
            let c = e.all_coeffs();
            return Ok(vec![(crate::ops::div(&crate::ops::neg(&c[1]), &c[0]), 1)]);
        }
        preflight(p, limits)?;
        let factorization = raw(p, false)?;
        verify(p, &factorization)?;
        let result = if let Some(f) = &p.algebraic {
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
        // The reference keeps the order of the linear factors of factor_list.
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
    #[test]
    fn typed_squarefree_parts_and_distinct_factor_list_conventions() {
        let x = Expr::symbol("x");
        let mut p =
            Poly::from_expr(&parse_str("(2*x+1)^2/3").unwrap(), std::slice::from_ref(&x)).unwrap();
        p.domain = Domain::QQ;
        let p: Polynomial = p.into();
        let factors = compatible_list(&p, false, false, Limits::default()).unwrap();
        assert_eq!(factors.unit.as_expr(), parse_str("1/3").unwrap());
        assert_eq!(factors.factors[0].0.as_expr(), parse_str("2*x+1").unwrap());
        let squarefree = compatible_list(&p, true, true, Limits::default()).unwrap();
        assert_eq!(squarefree.unit.as_expr(), parse_str("4/3").unwrap());
        assert_eq!(squarefree.factors.len(), 2);
        assert!(squarefree.factors[0].0.rep.is_one());
        assert_eq!(
            squarefree.factors[1].0.as_expr(),
            parse_str("x+1/2").unwrap()
        );
        assert_eq!(
            squarefree_part(&p, Limits::default()).unwrap().as_expr(),
            parse_str("x+1/2").unwrap()
        );
        let mut finite =
            Poly::from_expr(&parse_str("(x^2+1)^2").unwrap(), std::slice::from_ref(&x)).unwrap();
        finite.modulus = Some(2);
        finite = finite.with_rep(finite.rep.clone());
        assert_eq!(
            squarefree_part(&finite.into(), Limits::default())
                .unwrap()
                .as_expr(),
            parse_str("x+1").unwrap()
        );
        let zero = Poly::from_expr(&Expr::zero(), std::slice::from_ref(&x)).unwrap();
        let mut zero: Polynomial = zero.into();
        zero.base.modulus = Some(5);
        assert!(squarefree_part(&zero, Limits::default()).unwrap().is_zero());
        let a = Expr::symbol("a");
        for domain in ["QQ(a)", "QQ[a]", "ZZ[a]"] {
            let f = FractionPoly::from_expr(
                &parse_str("a^2*(a*x+1)^2").unwrap(),
                std::slice::from_ref(&x),
                std::slice::from_ref(&a),
                domain,
            )
            .unwrap();
            let part = squarefree_part(&f.clone().into(), Limits::default()).unwrap();
            let expected = if domain.ends_with(')') {
                "x+1/a"
            } else {
                "a*x+1"
            };
            assert_eq!(part.as_expr(), parse_str(expected).unwrap());
            assert_eq!(part.field.unwrap().coefficient_domain, f.coefficient_domain);
        }
    }
    #[test]
    fn symbolic_numerators_have_separate_coefficients_and_natural_factor_domains() {
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        let p =
            Poly::from_expr(&parse_str("(2*x+1)^2/3").unwrap(), std::slice::from_ref(&x)).unwrap();
        let (result, denominator) = symbolic_list(&p.into(), true, Limits::default()).unwrap();
        assert_eq!(denominator, Expr::int(3));
        assert_eq!(result.factors[0].0.domain, Domain::ZZ);
        assert_eq!(result.factors[0].0.as_expr(), parse_str("2*x+1").unwrap());
        for expression in ["x/a", "(a*x+1)^2/a^2"] {
            let p = FractionPoly::from_expr(
                &parse_str(expression).unwrap(),
                std::slice::from_ref(&x),
                std::slice::from_ref(&a),
                "ZZ(a)",
            )
            .unwrap();
            let (result, denominator) =
                symbolic_list(&p.clone().into(), true, Limits::default()).unwrap();
            assert_eq!(result.factors[0].1, if expression == "x/a" { 1 } else { 2 });
            let value = crate::ops::div(&result.unit.as_expr(), &denominator);
            let mut expected = value;
            for (f, k) in result.factors {
                expected = crate::ops::mul2(&expected, &f.as_expr().powi(k as i64));
                if expression == "x/a" {
                    assert!(f.field.is_none());
                }
            }
            let difference = crate::ops::sub(&expected, &p.as_expr());
            assert!(crate::poly::api::cancel(&difference).unwrap().is_zero());
        }
        let b = Expr::symbol("b");
        let p = FractionPoly::from_expr(
            &parse_str("(a*x+1)*(b*x+1)/((a+1)*(b+1))").unwrap(),
            std::slice::from_ref(&x),
            &[a, b],
            "ZZ(a,b)",
        )
        .unwrap();
        let (result, _) = symbolic_list(&p.into(), false, Limits::default()).unwrap();
        assert_eq!(result.factors[0].0.as_expr(), parse_str("a*x+1").unwrap());
        assert_eq!(result.factors[1].0.as_expr(), parse_str("b*x+1").unwrap());
    }
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
        assert!(output
            .iter()
            .any(|(factor, k)| factor.rep == p("2").rep && *k == 1));
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
