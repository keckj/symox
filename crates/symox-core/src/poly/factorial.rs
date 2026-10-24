//! Certified greatest factorial factors over the represented characteristic-zero
//! domains. A factor (h,k) denotes h(x)h(x-1)...h(x-k+1), not h(x)^k.
//! Paule's gcd/shift descent needs no irreducible factorization. Every operation
//! shares one fallible budget; in particular shifts never densify an unchecked
//! sparse exponent and the implementation does not recurse on the degree.
use super::algebraic::AlgebraicScope;
use super::convert::{PResult, PolyError};
use super::domains::ExactDomain;
use super::fraction::Polynomial;
use super::mpoly::MPoly;
use super::polyobj::Domain;
use super::pseudo::{Budget, PseudoCoeff};
use super::ratfunc::Scope;

const MAX_DEGREE: i64 = 1024;
const MAX_STORAGE: u64 = 64_000_000;
fn failure(s: &str) -> PolyError {
    PolyError::Unsupported(s.into())
}
struct Engine {
    budget: Budget,
    retained: u64,
}
impl Engine {
    fn new() -> Self {
        Self {
            budget: Budget::new(),
            retained: 0,
        }
    }
    fn check<C: PseudoCoeff>(&self, p: &MPoly<C>) -> PResult<()> {
        if p.nvars != 1 {
            return Err(failure(
                "greatest factorial factorization requires one polynomial generator",
            ));
        }
        if p.degree(0) > MAX_DEGREE {
            return Err(failure(
                "greatest factorial factorization resource limit: degree exceeds 1024",
            ));
        }
        self.budget.check(p)
    }
    fn retain<C: PseudoCoeff>(&mut self, p: &MPoly<C>) -> PResult<()> {
        self.check(p)?;
        self.budget.retain(p, &mut self.retained)
    }
    fn monic<C: PseudoCoeff>(&mut self, p: &MPoly<C>) -> PResult<MPoly<C>> {
        if p.is_zero() {
            return Ok(p.clone());
        }
        self.budget.exact_div(p, &MPoly::constant(p.lc(), 1))
    }
    fn remainder<C: PseudoCoeff>(&mut self, a: &MPoly<C>, b: &MPoly<C>) -> PResult<MPoly<C>> {
        if b.is_zero() {
            return Err(failure("invalid factorial-factor gcd divisor"));
        }
        let mut r = a.clone();
        while !r.is_zero() && r.degree(0) >= b.degree(0) {
            self.budget.charge(1)?;
            let previous = r.degree(0);
            let coefficient = r.lc().checked_div(&b.lc())?;
            let term = MPoly::monomial(&[(previous - b.degree(0)) as u32], coefficient);
            let product = self.budget.mul(&term, b)?;
            r = self.budget.merge(&r, &product, true)?;
            if !r.is_zero() && r.degree(0) >= previous {
                return Err(failure("factorial-factor gcd did not lower its degree"));
            }
        }
        Ok(r)
    }
    fn gcd<C: PseudoCoeff>(&mut self, a: &MPoly<C>, b: &MPoly<C>) -> PResult<MPoly<C>> {
        let (mut a, mut b) = (a.clone(), b.clone());
        while !b.is_zero() {
            self.budget.charge(1)?;
            let r = self.remainder(&a, &b)?;
            a = b;
            b = self.monic(&r)?;
        }
        self.monic(&a)
    }
    fn shift<C: PseudoCoeff>(&mut self, p: &MPoly<C>, by: i64) -> PResult<MPoly<C>> {
        self.check(p)?;
        if p.is_zero() || by == 0 {
            return Ok(p.clone());
        }
        let n = p.degree(0) as u64;
        // The Taylor shift can fill every position through the original degree.
        // Bound the container and the *minimum* coefficient storage first;
        // Budget.{mul,merge} check actual coefficient growth before arithmetic.
        let weight = p
            .coeffs()
            .iter()
            .map(PseudoCoeff::weight)
            .max()
            .unwrap_or(0);
        let storage = (n + 1).saturating_mul(weight.saturating_add(32));
        if storage > MAX_STORAGE {
            return Err(failure(
                "greatest factorial factorization resource limit: shift storage",
            ));
        }
        self.budget.charge((n + 1).saturating_mul(n + 1))?;
        let linear = MPoly::from_terms(1, vec![(vec![1], C::one()), (vec![0], C::from_i64(by))]);
        let mut result = MPoly::zero(1);
        let mut terms = p.iter().peekable();
        for power in (0..=n as u32).rev() {
            result = self.budget.mul(&result, &linear)?;
            if terms.peek().is_some_and(|(e, _)| e[0] == power) {
                let (_, coefficient) = terms.next().expect("peeked coefficient");
                result =
                    self.budget
                        .merge(&result, &MPoly::constant(coefficient.clone(), 1), false)?;
            }
        }
        self.check(&result)?;
        Ok(result)
    }
    fn falling<C: PseudoCoeff>(&mut self, h: &MPoly<C>, k: u32) -> PResult<MPoly<C>> {
        let mut result = MPoly::one(1);
        for j in 0..k {
            let shifted = self.shift(h, -i64::from(j))?;
            result = self.budget.mul(&result, &shifted)?;
        }
        Ok(result)
    }
    fn factor<C: PseudoCoeff>(
        &mut self,
        f: &MPoly<C>,
        belongs: impl Fn(&MPoly<C>) -> bool,
    ) -> PResult<Vec<(MPoly<C>, u32)>> {
        self.check(f)?;
        if f.is_zero() {
            return Err(failure(
                "greatest factorial factorization does not exist for a zero polynomial",
            ));
        }
        if f.degree(0) == 0 {
            return Ok(Vec::new());
        }
        let normalized = self.monic(f)?;
        if !belongs(&normalized) {
            return Err(PolyError::ExactQuotientFailed(
                "the leading coefficient does not divide every coefficient in the selected domain"
                    .into(),
            ));
        }
        let mut p = normalized.clone();
        let mut stack = Vec::new();
        while p.degree(0) > 0 {
            let shifted = self.shift(&p, 1)?;
            let g = self.gcd(&p, &shifted)?;
            if g.degree(0) >= p.degree(0) {
                return Err(failure(
                    "factorial-factor gcd/shift descent did not lower its degree",
                ));
            }
            self.retain(&p)?;
            stack.push(p);
            p = g;
        }
        let mut factors: Vec<(MPoly<C>, u32)> = Vec::new();
        while let Some(p) = stack.pop() {
            self.budget.charge(1)?;
            // Fundamental lemma: gcd(p,p(x+1)) = product F(h,k-1).
            let shifted = self.shift(&p, 1)?;
            let mut divisor = self.gcd(&p, &shifted)?;
            for (h, k) in &mut factors {
                let shifted = self.shift(h, -i64::from(*k))?;
                divisor = self.budget.mul(&divisor, &shifted)?;
                *k = k
                    .checked_add(1)
                    .ok_or_else(|| failure("factorial length overflow"))?;
            }
            let remainder = self.budget.exact_div(&p, &divisor)?;
            if remainder.degree(0) > 0 {
                factors.insert(0, (remainder, 1));
            }
        }
        let mut product = MPoly::one(1);
        let mut falling = Vec::new();
        let mut degree = 0u64;
        for (h, k) in &factors {
            if *k == 0 || h.degree(0) <= 0 || !h.lc().is_one() || !belongs(h) {
                return Err(failure(
                    "invalid greatest factorial factor metadata or coefficient domain",
                ));
            }
            degree = degree.saturating_add((*k as u64).saturating_mul(h.degree(0) as u64));
            let chain = self.falling(h, *k)?;
            product = self.budget.mul(&product, &chain)?;
            self.retain(&chain)?;
            falling.push(chain);
        }
        if product != normalized || degree != normalized.degree(0) as u64 {
            return Err(failure(
                "greatest factorial factorization failed its exact reconstruction certificate",
            ));
        }
        // GFF3 prohibits overlaps extending any maximal falling-factor chain.
        for (i, chain) in falling.iter().enumerate() {
            for (h, k) in factors.iter().skip(i) {
                for shift in [1, -i64::from(*k)] {
                    let shifted = self.shift(h, shift)?;
                    if !self.gcd(chain, &shifted)?.is_one() {
                        return Err(failure(
                            "greatest factorial factorization failed its maximality certificate",
                        ));
                    }
                }
            }
        }
        Ok(factors)
    }
}

/// The returned factors reconstruct the monic normalization of `p`. The input
/// coefficient domain must permit that normalization without promotion.
pub fn gff_list(p: &Polynomial) -> PResult<Vec<(Polynomial, u32)>> {
    if p.ngens != 1 {
        return Err(failure(
            "greatest factorial factorization requires one polynomial generator",
        ));
    }
    if p.modulus.is_some() {
        return Err(failure(
            "greatest factorial factorization requires characteristic zero",
        ));
    }
    let mut engine = Engine::new();
    if let Some(f) = &p.algebraic {
        let _scope = AlgebraicScope::new(f.number_field.clone());
        return engine.factor(&f.rep, |_| true).map(|fs| {
            fs.into_iter()
                .map(|(rep, k)| (f.with_rep(rep).into(), k))
                .collect()
        });
    }
    if let Some(f) = &p.field {
        let _scope = Scope::new(f.parameters.len());
        let ring = matches!(
            f.coefficient_domain,
            ExactDomain::Functions { ring: true, .. }
        );
        return engine
            .factor(&f.rep, |rep| {
                !ring || f.with_rep(rep.clone()).belongs_to_ring()
            })
            .map(|fs| {
                fs.into_iter()
                    .map(|(rep, k)| (f.with_rep(rep).into(), k))
                    .collect()
            });
    }
    if p.domain.is_rr() || !p.domain_gens().is_empty() {
        return Err(failure(
            "greatest factorial factorization requires an implemented exact coefficient domain",
        ));
    }
    engine
        .factor(&p.rep, |rep| p.domain != Domain::ZZ || rep.is_integral())
        .map(|fs| {
            fs.into_iter()
                .map(|(rep, k)| (p.base.with_rep(rep).into(), k))
                .collect()
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{polyobj::Poly, Rat};
    use crate::{parser::parse_str, Expr};
    fn p(s: &str) -> Polynomial {
        Poly::from_expr(&parse_str(s).unwrap(), &[Expr::symbol("x")])
            .unwrap()
            .into()
    }
    #[test]
    fn exact_factorial_descent_and_units() {
        let f = p("2*x^5+4*x^4-2*x^3-4*x^2");
        let result = gff_list(&f).unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0.rep, p("x").rep);
        assert_eq!(result[0].1, 1);
        assert_eq!(result[1].0.rep, p("x+2").rep);
        assert_eq!(result[1].1, 4);
        assert!(gff_list(&p("7")).unwrap().is_empty());
        assert!(gff_list(&p("0")).is_err());
        assert!(matches!(
            gff_list(&p("2*x+1")),
            Err(PolyError::ExactQuotientFailed(_))
        ));
        assert!(gff_list(&p(
            "(x^3+7)*( (x-1)^3+7)*(x^2+5*x)*((x-1)^2+5*(x-1))*((x-2)^2+5*(x-2))"
        ))
        .is_ok());
    }
    #[test]
    fn characteristic_zero_strict_degree_and_sparse_budget() {
        let mut f = p("x^2-x");
        f.base.modulus = Some(2);
        assert!(gff_list(&f)
            .unwrap_err()
            .to_string()
            .contains("characteristic zero"));
        let f = p("x^1000000000+1");
        assert!(gff_list(&f)
            .unwrap_err()
            .to_string()
            .contains("resource limit"));
        let mut engine = Engine::new();
        let sparse = MPoly::from_terms(1, vec![(vec![1000000000], Rat::one())]);
        assert!(engine.shift(&sparse, 1).is_err());
        assert!(gff_list(&p("x*(x-1)*(x-2)")).is_ok());
    }
    #[test]
    fn coefficient_fields_rings_and_scopes_have_exact_factorial_chains() {
        use crate::poly::algebraic::{AlgebraicPoly, NumberField};
        use crate::poly::fraction::FractionPoly;
        let x = Expr::symbol("x");
        let a = Expr::symbol("a");
        let b = Expr::symbol("unused");
        for domain in [
            "QQ(a,unused)",
            "ZZ(a,unused)",
            "QQ[a,unused]",
            "ZZ[a,unused]",
        ] {
            let f = FractionPoly::from_expr(
                &parse_str("a*(x+a)*(x+a-1)*(x+a-2)").unwrap(),
                std::slice::from_ref(&x),
                &[a.clone(), b.clone()],
                domain,
            )
            .unwrap();
            let fs = gff_list(&f.clone().into()).unwrap();
            assert_eq!(fs.len(), 1);
            assert_eq!(fs[0].1, 3);
            let got = fs[0].0.field.as_ref().unwrap();
            assert_eq!(got.parameters, f.parameters);
            assert_eq!(got.coefficient_domain, f.coefficient_domain);
            assert_eq!(got.as_expr(), parse_str("x+a").unwrap());
        }
        let f = FractionPoly::from_expr(
            &parse_str("a*x+1").unwrap(),
            std::slice::from_ref(&x),
            &[a, b],
            "QQ[a,unused]",
        )
        .unwrap();
        assert!(matches!(
            gff_list(&f.into()),
            Err(PolyError::ExactQuotientFailed(_))
        ));
        for generator in ["sqrt(2)", "I", "2^(1/3)+I"] {
            let field = NumberField::from_generator(&parse_str(generator).unwrap()).unwrap();
            let expression = parse_str(&format!("(x+({generator}))*(x+({generator})-1)")).unwrap();
            let f = AlgebraicPoly::from_expr(&expression, std::slice::from_ref(&x), field).unwrap();
            let fs = gff_list(&f.clone().into()).unwrap();
            assert_eq!(fs.len(), 1);
            assert_eq!(fs[0].1, 2);
            assert_eq!(
                fs[0].0.algebraic.as_ref().unwrap().number_field,
                f.number_field
            );
        }
        assert_eq!(gff_list(&p("x^2")).unwrap()[0].1, 1);
    }
}
