//! Checked Brown--Traub subresultant PRS over the supported exact rings.
//!
//! The scalar subresultants are retained: after an abnormal degree drop the
//! final constant polynomial in the PRS need not equal the resultant.
use super::convert::{PResult, PolyError};
use super::mpoly::MPoly;
use super::pseudo::{pdiv_with_budget, Budget, PseudoCoeff};

pub(crate) struct Prs<C> {
    pub sequence: Vec<MPoly<C>>,
    pub resultant: MPoly<C>,
}

pub(crate) fn subresultants<C: PseudoCoeff>(
    first: &MPoly<C>,
    second: &MPoly<C>,
    variable: usize,
) -> PResult<Prs<C>> {
    let mut budget = Budget::new();
    budget.check(first)?;
    budget.check(second)?;
    if first.nvars != second.nvars || variable >= first.nvars {
        return Err(PolyError::Failed("invalid subresultant generators".into()));
    }
    let unit = C::one();
    if first
        .coeffs()
        .iter()
        .chain(second.coeffs())
        .any(|c| !unit.compatible(c))
    {
        return Err(PolyError::Failed(
            "incompatible subresultant coefficient context".into(),
        ));
    }
    let dimensions = first.nvars;
    let zero = MPoly::zero(dimensions);
    let (mut f, mut g) = (first.clone(), second.clone());
    let swapped = f.degree(variable) < g.degree(variable);
    if swapped {
        std::mem::swap(&mut f, &mut g);
    }
    if f.is_zero() {
        return Ok(Prs {
            sequence: Vec::new(),
            resultant: zero,
        });
    }
    if g.is_zero() {
        return Ok(Prs {
            sequence: vec![f],
            resultant: zero,
        });
    }
    let reversed_sign = swapped && f.degree(variable) % 2 == 1 && g.degree(variable) % 2 == 1;
    let mut degree = g.degree(variable);
    let mut gap = (f.degree(variable) - degree) as u64;
    let mut retained = 0;
    budget.retain(&f, &mut retained)?;
    budget.retain(&g, &mut retained)?;
    let mut sequence = vec![f.clone(), g.clone()];
    let mut h = pdiv_with_budget(&f, &g, variable, &mut budget)?.1;
    if (gap + 1) % 2 == 1 {
        h = h.neg();
    }
    let mut leading = g.lc_in(variable);
    let mut scalar = budget.pow(&leading, gap)?;
    let mut c = scalar.neg();
    // This includes constants: res(c,d)=1 under the empty Sylvester
    // determinant convention, while res(f,d)=d^degree(f).
    while !h.is_zero() {
        budget.charge(1)?;
        budget.check(&h)?;
        let next_degree = h.degree(variable);
        if next_degree >= degree {
            return Err(PolyError::Failed(
                "subresultant PRS did not lower the degree".into(),
            ));
        }
        budget.retain(&h, &mut retained)?;
        sequence.push(h.clone());
        f = g;
        g = h;
        gap = (degree - next_degree) as u64;
        degree = next_degree;
        let power = budget.pow(&c, gap)?;
        let denominator = budget.mul(&leading.neg(), &power)?;
        let remainder = pdiv_with_budget(&f, &g, variable, &mut budget)?.1;
        h = budget.exact_div(&remainder, &denominator)?;
        leading = g.lc_in(variable);
        c = if gap > 1 {
            let numerator = budget.pow(&leading.neg(), gap)?;
            let denominator = budget.pow(&c, gap - 1)?;
            budget.exact_div(&numerator, &denominator)?
        } else {
            leading.neg()
        };
        scalar = c.neg();
    }
    let resultant = if g.degree(variable) > 0 {
        zero
    } else if reversed_sign {
        scalar.neg()
    } else {
        scalar
    };
    Ok(Prs {
        sequence,
        resultant,
    })
}

pub(crate) fn discriminant<C: PseudoCoeff>(f: &MPoly<C>, variable: usize) -> PResult<MPoly<C>> {
    let mut budget = Budget::new();
    budget.check(f)?;
    if variable >= f.nvars {
        return Err(PolyError::Failed("invalid discriminant generator".into()));
    }
    let degree = f.degree(variable);
    if degree <= 0 {
        return Ok(MPoly::zero(f.nvars));
    }
    // Derivative coefficients may grow before a normal MPoly::derivative
    // could validate them. Construct through checked multiplication instead.
    let mut derivative = Vec::with_capacity(f.len());
    for (exponent, coefficient) in f.iter() {
        if exponent[variable] == 0 {
            continue;
        }
        budget.charge(1)?;
        let mut e = exponent.to_vec();
        let c = coefficient.checked_mul(&C::from_i64(e[variable] as i64))?;
        e[variable] -= 1;
        derivative.push((e, c));
    }
    let derivative = MPoly::from_terms(f.nvars, derivative);
    let resultant = subresultants(f, &derivative, variable)?.resultant;
    let value = budget.exact_div(&resultant, &f.lc_in(variable))?;
    // n(n-1)/2 is odd exactly for n=2 or 3 modulo four, without overflow.
    Ok(if degree % 4 == 2 || degree % 4 == 3 {
        value.neg()
    } else {
        value
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        int::{Coeff, Rat},
        polyobj::Poly,
    };
    use crate::{parser::parse_str, Expr};
    fn p(source: &str) -> MPoly<Rat> {
        Poly::from_expr(
            &parse_str(source).unwrap(),
            &[Expr::symbol("x"), Expr::symbol("y")],
        )
        .unwrap()
        .rep
    }
    #[test]
    fn abnormal_drops_contents_and_original_order_sign() {
        let f = p("x^4+1");
        let g = p("x^2+1");
        let r = subresultants(&f, &g, 0).unwrap();
        assert_eq!(r.sequence, vec![f, g, p("-2")]);
        assert_eq!(r.resultant, p("4"));
        let f = p("x+1");
        let g = p("x^3+2");
        assert_eq!(subresultants(&f, &g, 0).unwrap().resultant, p("1"));
        assert_eq!(subresultants(&g, &f, 0).unwrap().resultant, p("-1"));
        let f = p("3*x^2*y-y^3-4");
        let g = p("x^2+x*y^3-9");
        let r = subresultants(&f, &g, 0).unwrap();
        assert_eq!(r.sequence[2], p("3*x*y^4+y^3-27*y+4"));
        assert_eq!(
            r.resultant,
            p("-3*y^10-12*y^7+y^6-54*y^4+8*y^3+729*y^2-216*y+16")
        );
    }
    #[test]
    fn degenerate_inputs_and_discriminants() {
        let zero = MPoly::zero(2);
        assert!(subresultants(&zero, &zero, 0).unwrap().sequence.is_empty());
        assert_eq!(
            subresultants(&zero, &p("x"), 0).unwrap().sequence,
            vec![p("x")]
        );
        assert_eq!(
            subresultants(&p("2"), &p("3"), 0).unwrap().resultant,
            p("1")
        );
        assert_eq!(
            subresultants(&p("2*x^3"), &p("3"), 0).unwrap().resultant,
            p("27")
        );
        assert!(subresultants(&p("x^3"), &p("2*x"), 0)
            .unwrap()
            .resultant
            .is_zero());
        assert_eq!(discriminant(&p("y*x^2+x+1"), 0).unwrap(), p("1-4*y"));
        assert_eq!(discriminant(&p("2*x+1"), 0).unwrap(), p("1"));
        assert!(discriminant(&p("2"), 0).unwrap().is_zero());
        assert_eq!(Rat::one(), Rat::from_i64(1));
    }

    fn determinant<C: Coeff>(matrix: &[Vec<MPoly<C>>], variables: usize) -> MPoly<C> {
        if matrix.is_empty() {
            return MPoly::one(variables);
        }
        let mut result = MPoly::zero(variables);
        for column in 0..matrix.len() {
            let minor: Vec<_> = matrix[1..]
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter(|(j, _)| *j != column)
                        .map(|(_, c)| c.clone())
                        .collect()
                })
                .collect();
            let term = matrix[0][column].mul(&determinant(&minor, variables));
            result = if column % 2 == 0 {
                result.add(&term)
            } else {
                result.sub(&term)
            };
        }
        result
    }
    fn sylvester<C: Coeff>(f: &MPoly<C>, g: &MPoly<C>) -> MPoly<C> {
        let n = f.degree(0) as usize;
        let m = g.degree(0) as usize;
        let coefficient = |p: &MPoly<C>, k: usize| {
            MPoly::from_terms(
                p.nvars,
                p.iter()
                    .filter(|(e, _)| e[0] == k as u32)
                    .map(|(e, c)| {
                        let mut e = e.to_vec();
                        e[0] = 0;
                        (e, c.clone())
                    })
                    .collect(),
            )
        };
        let mut matrix = vec![vec![MPoly::zero(f.nvars); n + m]; n + m];
        for (r, row) in matrix.iter_mut().take(m).enumerate() {
            for k in 0..=n {
                row[r + k] = coefficient(f, n - k);
            }
        }
        for (r, row) in matrix.iter_mut().skip(m).enumerate() {
            for k in 0..=m {
                row[r + k] = coefficient(g, m - k);
            }
        }
        determinant(&matrix, f.nvars)
    }
    fn exact_certificates<C: PseudoCoeff>(f: &MPoly<C>, g: &MPoly<C>) {
        let forward = subresultants(f, g, 0).unwrap();
        assert_eq!(forward.resultant, sylvester(f, g));
        let reverse = subresultants(g, f, 0).unwrap().resultant;
        assert_eq!(
            reverse,
            if f.degree(0) % 2 == 1 && g.degree(0) % 2 == 1 {
                forward.resultant.neg()
            } else {
                forward.resultant.clone()
            }
        );
        let derivative = f.derivative(0);
        let value = discriminant(f, 0).unwrap().mul(&f.lc_in(0));
        let expected = sylvester(f, &derivative);
        assert_eq!(
            value,
            if f.degree(0) % 4 == 2 || f.degree(0) % 4 == 3 {
                expected.neg()
            } else {
                expected
            }
        );
        for stages in forward.sequence.windows(3) {
            assert!(stages[2].degree(0) < stages[1].degree(0));
            let remainder = super::super::pseudo::pdiv(&stages[0], &stages[1], 0)
                .unwrap()
                .1;
            assert!(remainder.div_exact(&stages[2]).is_some());
        }
    }
    #[test]
    fn sylvester_certificates_and_scopes_for_all_coefficient_families() {
        use crate::poly::{
            algebraic::{AlgebraicPoly, AlgebraicScope, NumberField, AN},
            f5b::{Fp, PrimeScope},
            fraction::FractionPoly,
            ratfunc::{Scope, RF},
        };
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let a = Expr::symbol("a");
        let gens = vec![x, y];
        exact_certificates(&p("x^3+2*x^2+y*x+1"), &p("2*x^2+x+y+1"));
        let _outer_prime = PrimeScope::new(7).unwrap();
        for prime in [2, 5, 101] {
            let _scope = PrimeScope::new(prime).unwrap();
            let to_fp = |p: MPoly<Rat>| p.map_coeffs(|r| Fp::from_rat(r, prime).unwrap());
            exact_certificates(&to_fp(p("x^3+2*x^2+y*x+1")), &to_fp(p("2*x^2+x+y+1")));
            if prime == 5 {
                assert_eq!(discriminant(&to_fp(p("x^2+1")), 0).unwrap(), to_fp(p("1")));
                assert!(subresultants(&to_fp(p("x^5+1")), &MPoly::zero(2), 0)
                    .unwrap()
                    .resultant
                    .is_zero());
                assert!(discriminant(&to_fp(p("x^5+1")), 0).unwrap().is_zero());
            }
        }
        assert!(Fp::from_i64(14).is_zero());
        let _outer_rf = Scope::new(3);
        for domain in ["ZZ[a]", "QQ[a]", "QQ(a)"] {
            let f = FractionPoly::from_expr(
                &parse_str("a*x^2+a*x+a").unwrap(),
                &gens,
                std::slice::from_ref(&a),
                domain,
            )
            .unwrap();
            let g = f.convert(&parse_str("2*a*x+2*a").unwrap()).unwrap();
            let _scope = Scope::new(1);
            exact_certificates(&f.rep, &g.rep);
            let expected = f.convert(&parse_str("4*a^3").unwrap()).unwrap();
            assert_eq!(
                subresultants(&f.rep, &g.rep, 0).unwrap().resultant,
                expected.rep
            );
        }
        assert_eq!(RF::one().n.nvars, 3);
        let outer = NumberField::from_generator(&parse_str("sqrt(3)").unwrap()).unwrap();
        let _outer_an = AlgebraicScope::new(outer.clone());
        for generator in ["sqrt(2)", "2^(1/3)+I"] {
            let field = NumberField::from_generator(&parse_str(generator).unwrap()).unwrap();
            let f = AlgebraicPoly::from_expr(
                &parse_str(&format!("x^3+({generator})*x^2+y*x+1")).unwrap(),
                &gens,
                field.clone(),
            )
            .unwrap();
            let g = f.convert(&parse_str("x^2+x+y+1").unwrap()).unwrap();
            let _scope = AlgebraicScope::new(field);
            exact_certificates(&f.rep, &g.rep);
        }
        assert_eq!(AN::one().field, outer);
    }
}
