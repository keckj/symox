//! Pseudo-division over exact coefficient rings, with bounded intermediates.
use super::algebraic::AN;
use super::convert::{PResult, PolyError};
use super::f5b::Fp;
use super::int::{Coeff, Int, Rat};
use super::mpoly::MPoly;
use super::ratfunc::RF;
use std::collections::BTreeMap;

const MAX_TERMS: usize = 100_000;
const MAX_BITS: u64 = 262_144;
const MAX_STORAGE: u64 = 64_000_000;
const MAX_WORK: u64 = 10_000_000;
const MAX_INNER_PAIRS: usize = 20_000;
fn limit(reason: &str) -> PolyError {
    PolyError::Unsupported(format!("pseudo-division resource limit: {reason}"))
}
fn bits(n: u64) -> PResult<()> {
    if n > MAX_BITS {
        Err(limit("coefficient size"))
    } else {
        Ok(())
    }
}
fn rat_weight(a: &Rat) -> u64 {
    128 + a.n.bits() + a.d.bits()
}
fn rat_pair(a: &Rat, b: &Rat, add: bool) -> PResult<()> {
    let n = if add {
        (a.n.bits() + b.d.bits()).max(b.n.bits() + a.d.bits()) + 1
    } else {
        a.n.bits() + b.n.bits()
    };
    bits(n)?;
    bits(a.d.bits() + b.d.bits())
}
fn int_poly_weight(a: &MPoly<Int>) -> u64 {
    a.coeffs()
        .iter()
        .map(|c| 128 + c.bits() + a.nvars as u64 * 32)
        .sum()
}
fn int_product(a: &MPoly<Int>, b: &MPoly<Int>) -> PResult<()> {
    if a.nvars != b.nvars {
        return Err(PolyError::Failed(
            "incompatible rational-function coefficient dimensions".into(),
        ));
    }
    let pairs = a.len().saturating_mul(b.len());
    if pairs > MAX_INNER_PAIRS {
        return Err(limit("rational-function coefficient terms"));
    }
    for (x, y) in a.degrees().iter().zip(b.degrees()) {
        if x.checked_add(y).is_none() {
            return Err(limit("coefficient exponent overflow"));
        }
    }
    let ba = a.coeffs().iter().map(Int::bits).max().unwrap_or(0);
    let bb = b.coeffs().iter().map(Int::bits).max().unwrap_or(0);
    let product_bits = ba + bb + 64 - (a.len().min(b.len()).max(1) as u64).leading_zeros() as u64;
    bits(product_bits)?;
    let storage =
        (pairs as u64).saturating_mul(product_bits.saturating_add(128 + a.nvars as u64 * 32));
    if storage > MAX_STORAGE {
        return Err(limit("rational-function coefficient storage"));
    }
    Ok(())
}
/// Local arithmetic preflights; existing public Coeff implementations remain
/// unrestricted outside pseudo-division. Factory scopes are supplied by each
/// polynomial adapter, as for the other generic exact polynomial algorithms.
pub(crate) trait PseudoCoeff: Coeff {
    fn weight(&self) -> u64;
    fn validate_size(&self) -> PResult<()>;
    fn compatible(&self, other: &Self) -> bool;
    fn preflight(&self, other: &Self, add: bool) -> PResult<()>;
    fn checked_mul(&self, other: &Self) -> PResult<Self> {
        self.validate_size()?;
        other.validate_size()?;
        if !self.compatible(other) {
            return Err(PolyError::Failed(
                "incompatible pseudo-division coefficient domains".into(),
            ));
        }
        if self.is_zero() || other.is_zero() {
            return Ok(Self::zero());
        }
        if self.is_one() {
            return Ok(other.clone());
        }
        if other.is_one() {
            return Ok(self.clone());
        }
        self.preflight(other, false)?;
        let value = self.mul(other);
        value.validate_size()?;
        Ok(value)
    }
    fn checked_add(&self, other: &Self) -> PResult<Self> {
        self.validate_size()?;
        other.validate_size()?;
        if !self.compatible(other) {
            return Err(PolyError::Failed(
                "incompatible pseudo-division coefficient domains".into(),
            ));
        }
        if self.is_zero() {
            return Ok(other.clone());
        }
        if other.is_zero() {
            return Ok(self.clone());
        }
        self.preflight(other, true)?;
        let value = self.add(other);
        value.validate_size()?;
        Ok(value)
    }
}
impl PseudoCoeff for Int {
    fn validate_size(&self) -> PResult<()> {
        bits(self.bits())
    }
    fn weight(&self) -> u64 {
        128 + self.bits()
    }
    fn compatible(&self, _: &Self) -> bool {
        true
    }
    fn preflight(&self, b: &Self, add: bool) -> PResult<()> {
        bits(if add {
            self.bits().max(b.bits()) + 1
        } else {
            self.bits() + b.bits()
        })
    }
}
impl PseudoCoeff for Rat {
    fn validate_size(&self) -> PResult<()> {
        bits(self.n.bits())?;
        bits(self.d.bits())
    }
    fn weight(&self) -> u64 {
        rat_weight(self)
    }
    fn compatible(&self, _: &Self) -> bool {
        true
    }
    fn preflight(&self, b: &Self, add: bool) -> PResult<()> {
        rat_pair(self, b, add)
    }
}
impl PseudoCoeff for Fp {
    fn validate_size(&self) -> PResult<()> {
        Ok(())
    }
    fn weight(&self) -> u64 {
        192
    }
    fn compatible(&self, b: &Self) -> bool {
        self.modulus() == b.modulus()
    }
    fn preflight(&self, _: &Self, _: bool) -> PResult<()> {
        Ok(())
    }
}
impl PseudoCoeff for RF {
    fn validate_size(&self) -> PResult<()> {
        for c in self.n.coeffs().iter().chain(self.d.coeffs()) {
            bits(c.bits())?;
        }
        Ok(())
    }
    fn weight(&self) -> u64 {
        256 + int_poly_weight(&self.n) + int_poly_weight(&self.d)
    }
    fn compatible(&self, b: &Self) -> bool {
        self.n.nvars == b.n.nvars
            && self.d.nvars == self.n.nvars
            && b.d.nvars == b.n.nvars
            && !self.d.is_zero()
            && !b.d.is_zero()
    }
    fn preflight(&self, b: &Self, add: bool) -> PResult<()> {
        if add && self.d == b.d {
            for (e, a) in self.n.iter() {
                let other = b.n.coeff_at(e);
                if !other.is_zero() {
                    bits(
                        a.bits().max(other.bits())
                            + u64::from(a.is_negative() == other.is_negative()),
                    )?;
                }
            }
            return Ok(());
        }
        int_product(&self.d, &b.d)?;
        if add {
            int_product(&self.n, &b.d)?;
            int_product(&b.n, &self.d)
        } else {
            int_product(&self.n, &b.n)
        }
    }
}
impl PseudoCoeff for AN {
    fn validate_size(&self) -> PResult<()> {
        for c in &self.coefficients {
            c.validate_size()?;
        }
        Ok(())
    }
    fn weight(&self) -> u64 {
        256 + self.coefficients.iter().map(rat_weight).sum::<u64>()
    }
    fn compatible(&self, b: &Self) -> bool {
        self.field.same_embedding(&b.field)
    }
    fn preflight(&self, b: &Self, add: bool) -> PResult<()> {
        if add {
            for (a, b) in self.coefficients.iter().zip(&b.coefficients) {
                rat_pair(a, b, true)?;
            }
            return Ok(());
        }
        let a = self.coefficients.iter().map(rat_weight).max().unwrap_or(0);
        let c = b.coefficients.iter().map(rat_weight).max().unwrap_or(0);
        let count = self.coefficients.len().min(b.coefficients.len()) as u64;
        let reduction = (self.coefficients.len() + b.coefficients.len())
            .saturating_sub(self.field.degree() + 1) as u64;
        let minimal = self.field.minimal.iter().map(rat_weight).max().unwrap_or(0);
        bits(
            (a + c).saturating_mul(count).saturating_add(
                reduction
                    .saturating_mul(minimal)
                    .saturating_mul(self.field.degree() as u64),
            ),
        )
    }
}
struct Budget {
    work: u64,
}
impl Budget {
    fn charge(&mut self, n: u64) -> PResult<()> {
        self.work = self
            .work
            .checked_sub(n)
            .ok_or_else(|| limit("polynomial work"))?;
        Ok(())
    }
    fn check<C: PseudoCoeff>(&self, p: &MPoly<C>) -> PResult<()> {
        for c in p.coeffs() {
            c.validate_size()?;
        }
        if p.len() > MAX_TERMS {
            return Err(limit("polynomial terms"));
        }
        let size = p.coeffs().iter().try_fold(0u64, |sum, c| {
            sum.checked_add(c.weight() + p.nvars as u64 * 32)
                .ok_or_else(|| limit("polynomial storage"))
        })?;
        if size > MAX_STORAGE {
            return Err(limit("polynomial storage"));
        }
        Ok(())
    }
    fn merge<C: PseudoCoeff>(
        &mut self,
        a: &MPoly<C>,
        b: &MPoly<C>,
        subtract: bool,
    ) -> PResult<MPoly<C>> {
        let mut map: BTreeMap<Vec<u32>, C> =
            a.iter().map(|(e, c)| (e.to_vec(), c.clone())).collect();
        let mut storage = a
            .coeffs()
            .iter()
            .map(|c| c.weight() + a.nvars as u64 * 32)
            .sum::<u64>();
        for (e, c) in b.iter() {
            self.charge(1 + c.weight() / 64)?;
            let c = if subtract { c.neg() } else { c.clone() };
            let c = if let Some(old) = map.remove(e) {
                storage = storage.saturating_sub(old.weight() + a.nvars as u64 * 32);
                old.checked_add(&c)?
            } else {
                c
            };
            if !c.is_zero() {
                storage = storage
                    .checked_add(c.weight() + a.nvars as u64 * 32)
                    .ok_or_else(|| limit("polynomial storage"))?;
                if storage > MAX_STORAGE || map.len() >= MAX_TERMS {
                    return Err(limit("polynomial storage or terms"));
                }
                map.insert(e.to_vec(), c);
            }
        }
        let p = MPoly::from_terms(a.nvars, map.into_iter().collect());
        self.check(&p)?;
        Ok(p)
    }
    fn mul<C: PseudoCoeff>(&mut self, a: &MPoly<C>, b: &MPoly<C>) -> PResult<MPoly<C>> {
        if a.is_zero() || b.is_zero() {
            return Ok(MPoly::zero(a.nvars));
        }
        for (x, y) in a.degrees().iter().zip(b.degrees()) {
            if x.checked_add(y).is_none() {
                return Err(limit("polynomial exponent overflow"));
            }
        }
        self.charge(a.len().saturating_mul(b.len()) as u64)?;
        let mut map: BTreeMap<Vec<u32>, C> = BTreeMap::new();
        let mut storage = 0u64;
        for (e, ca) in a.iter() {
            for (f, cb) in b.iter() {
                self.charge(1 + (ca.weight() + cb.weight()) / 64)?;
                let c = ca.checked_mul(cb)?;
                if c.is_zero() {
                    continue;
                }
                let powers: Vec<_> = e.iter().zip(f).map(|(x, y)| x + y).collect();
                let c = if let Some(old) = map.remove(&powers) {
                    storage = storage.saturating_sub(old.weight() + a.nvars as u64 * 32);
                    old.checked_add(&c)?
                } else {
                    c
                };
                if !c.is_zero() {
                    storage = storage
                        .checked_add(c.weight() + a.nvars as u64 * 32)
                        .ok_or_else(|| limit("polynomial storage"))?;
                    map.insert(powers, c);
                }
                if storage > MAX_STORAGE || map.len() > MAX_TERMS {
                    return Err(limit("polynomial storage or terms"));
                }
            }
        }
        Ok(MPoly::from_terms(a.nvars, map.into_iter().collect()))
    }
    fn pow<C: PseudoCoeff>(&mut self, a: &MPoly<C>, mut n: u64) -> PResult<MPoly<C>> {
        if n == 0 || a.is_one() {
            return Ok(MPoly::one(a.nvars));
        }
        let mut b = a.clone();
        let mut out = MPoly::one(a.nvars);
        while n != 0 {
            if n & 1 != 0 {
                out = self.mul(&out, &b)?;
            }
            n >>= 1;
            if n != 0 {
                b = self.mul(&b, &b)?;
            }
        }
        Ok(out)
    }
}
/// `LC(g)^max(deg(f)-deg(g)+1,0) f = q*g+r`, with `deg(r)<deg(g)`.
/// The selected variable is a main polynomial generator, never a coefficient
/// parameter. Multiplication and addition suffice over every supported ring.
pub(crate) fn pdiv<C: PseudoCoeff>(
    f: &MPoly<C>,
    g: &MPoly<C>,
    variable: usize,
) -> PResult<(MPoly<C>, MPoly<C>)> {
    if f.nvars != g.nvars || variable >= f.nvars {
        return Err(PolyError::Failed(
            "invalid pseudo-division generators".into(),
        ));
    }
    if g.is_zero() {
        return Err(PolyError::Unsupported("polynomial division by zero".into()));
    }
    let mut budget = Budget { work: MAX_WORK };
    budget.check(f)?;
    budget.check(g)?;
    let unit = C::one();
    if f.coeffs()
        .iter()
        .chain(g.coeffs())
        .any(|c| !unit.compatible(c))
    {
        return Err(PolyError::Failed(
            "incompatible pseudo-division coefficient context".into(),
        ));
    }
    let degree = g.degree(variable);
    let mut remaining = u64::try_from((f.degree(variable) - degree + 1).max(0))
        .map_err(|_| limit("scaling exponent overflow"))?;
    let lc = g.lc_in(variable);
    let mut q = MPoly::zero(f.nvars);
    let mut r = f.clone();
    while !r.is_zero() && r.degree(variable) >= degree {
        budget.charge(1)?;
        let previous = r.degree(variable);
        let shift =
            u32::try_from(previous - degree).map_err(|_| limit("quotient exponent overflow"))?;
        let lead = r.lc_in(variable);
        // lead has exponent zero in this variable, so this shift cannot
        // overflow any other exponent or exceed the input degree.
        let mut exponents = vec![0; f.nvars];
        exponents[variable] = shift;
        let term = lead.mul_term(&exponents, &unit);
        let lq = budget.mul(&lc, &q)?;
        q = budget.merge(&lq, &term, false)?;
        let lr = budget.mul(&lc, &r)?;
        let tg = budget.mul(&term, g)?;
        r = budget.merge(&lr, &tg, true)?;
        remaining = remaining
            .checked_sub(1)
            .ok_or_else(|| PolyError::Failed("pseudo-division degree invariant failed".into()))?;
        if !r.is_zero() && r.degree(variable) >= previous {
            return Err(PolyError::Failed(
                "pseudo-division did not lower the remainder degree".into(),
            ));
        }
    }
    if remaining != 0 {
        let scale = budget.pow(&lc, remaining)?;
        q = budget.mul(&scale, &q)?;
        r = budget.mul(&scale, &r)?;
    }
    Ok((q, r))
}

#[cfg(test)]
mod tests {
    use super::super::{
        algebraic::{AlgebraicPoly, AlgebraicScope, NumberField},
        f5b::PrimeScope,
        fraction::{FractionPoly, Polynomial},
        polyobj::Poly,
        ratfunc::Scope,
    };
    use super::*;
    use crate::{parser::parse_str, Expr};
    fn parse(s: &str) -> Expr {
        parse_str(s).unwrap()
    }
    fn certificate<C: Coeff>(a: &MPoly<C>, b: &MPoly<C>, q: &MPoly<C>, r: &MPoly<C>) {
        let exponent = (a.degree(0) - b.degree(0) + 1).max(0) as u32;
        assert_eq!(b.lc_in(0).pow(exponent).mul(a), q.mul(b).add(r));
        assert!(r.is_zero() || r.degree(0) < b.degree(0));
    }
    #[test]
    fn rational_multivariate_pseudo_division_preserves_dense_identity() {
        let gens = vec![Expr::symbol("x"), Expr::symbol("y")];
        for (f, g, q, r) in [
            ("x^2+1", "2*x+1", "2*x-1", "5"),
            ("x*y+x", "y", "x*y^2+x*y", "0"),
            ("x+1", "x^2+1", "0", "x+1"),
            ("0", "2*x+1", "0", "0"),
            ("x^2+1", "2", "4*x^2+4", "0"),
            (
                "(y+1)*x^2+y*x+1",
                "(y+2)*x+y",
                "(y+1)*(y+2)*x+y",
                "y^2+4*y+4-y^2",
            ),
        ] {
            let a = Poly::from_expr(&parse(f), &gens).unwrap();
            let b = Poly::from_expr(&parse(g), &gens).unwrap();
            let (u, v) = a.pdiv(&b).unwrap();
            certificate(&a.rep, &b.rep, &u.rep, &v.rep);
            if !f.starts_with("(y+") {
                assert_eq!(u.rep, Poly::from_expr(&parse(q), &gens).unwrap().rep);
                assert_eq!(v.rep, Poly::from_expr(&parse(r), &gens).unwrap().rep);
            }
        }
    }
    #[test]
    fn finite_and_algebraic_scopes_and_certificates() {
        let gens = vec![Expr::symbol("x"), Expr::symbol("y")];
        let _outer = PrimeScope::new(7).unwrap();
        for prime in [2, 5, 7, 101] {
            let a: Polynomial = Poly::from_expr(&parse("x^3+y*x+1"), &gens)
                .unwrap()
                .set_modulus(prime)
                .into();
            let b: Polynomial = Poly::from_expr(&parse("(y+1)*x+1"), &gens)
                .unwrap()
                .set_modulus(prime)
                .into();
            let (q, r) = a.pdiv(&b).unwrap();
            let _scope = PrimeScope::new(prime).unwrap();
            let fp = |p: &Polynomial| p.rep.map_coeffs(|c| Fp::from_rat(c, prime).unwrap());
            certificate(&fp(&a), &fp(&b), &fp(&q), &fp(&r));
            assert_eq!(q.modulus, Some(prime));
        }
        assert!(Fp::from_i64(14).is_zero());
        let outer = NumberField::from_generator(&parse("sqrt(3)")).unwrap();
        let _outer = AlgebraicScope::new(outer.clone());
        for generator in ["sqrt(2)", "I", "2^(1/3)"] {
            let field = NumberField::from_generator(&parse(generator)).unwrap();
            let c = parse(generator);
            let a: Polynomial = AlgebraicPoly::from_expr(
                &crate::ops::add(&[
                    parse("x^3+y*x+1"),
                    crate::ops::mul(&[c.clone(), parse("x^2")]),
                ]),
                &gens,
                field.clone(),
            )
            .unwrap()
            .into();
            let b: Polynomial = AlgebraicPoly::from_expr(
                &crate::ops::add(&[crate::ops::mul(&[c, parse("x")]), parse("y+1")]),
                &gens,
                field.clone(),
            )
            .unwrap()
            .into();
            let (q, r) = a.pdiv(&b).unwrap();
            let _scope = AlgebraicScope::new(field.clone());
            certificate(
                &a.algebraic.unwrap().rep,
                &b.algebraic.unwrap().rep,
                &q.algebraic.as_ref().unwrap().rep,
                &r.algebraic.unwrap().rep,
            );
            assert_eq!(q.algebraic.unwrap().number_field, field);
        }
        assert_eq!(AN::one().field, outer);
    }
    #[test]
    fn rational_function_and_ring_contexts_stay_closed() {
        let gens = vec![Expr::symbol("x"), Expr::symbol("y")];
        let parameters = vec![Expr::symbol("a"), Expr::symbol("b")];
        let _outer = Scope::new(3);
        for display in ["ZZ[a,b]", "QQ[a,b]", "ZZ(a,b)", "QQ(a,b)"] {
            let a: Polynomial =
                FractionPoly::from_expr(&parse("a*x^3+b*x*y+1"), &gens, &parameters, display)
                    .unwrap()
                    .into();
            let b: Polynomial =
                FractionPoly::from_expr(&parse("(a+b)*x+y+1"), &gens, &parameters, display)
                    .unwrap()
                    .into();
            let (q, r) = a.pdiv(&b).unwrap();
            let _scope = Scope::new(2);
            certificate(
                &a.field.as_ref().unwrap().rep,
                &b.field.unwrap().rep,
                &q.field.as_ref().unwrap().rep,
                &r.field.as_ref().unwrap().rep,
            );
            assert_eq!(
                q.field.as_ref().unwrap().coefficient_domain,
                a.field.as_ref().unwrap().coefficient_domain
            );
            assert_eq!(r.domain_str(), display);
            if display.ends_with(']') {
                assert!(q.field.unwrap().belongs_to_ring());
                assert!(r.field.unwrap().belongs_to_ring());
            }
        }
        assert_eq!(RF::one().n.nvars, 3);
    }
    #[test]
    fn bounds_reject_growth_before_overflow_but_allow_sparse_unit_scaling() {
        let n = u32::MAX;
        let f = MPoly::monomial(&[n], Rat::one());
        let one = MPoly::one(1);
        let (q, r) = pdiv(&f, &one, 0).unwrap();
        assert_eq!(q, f);
        assert!(r.is_zero());
        let huge = MPoly::monomial(&[1_000_000_000], Rat::one());
        let two_x = MPoly::monomial(&[1], Rat::from_i64(2));
        assert!(pdiv(&huge, &two_x, 0)
            .unwrap_err()
            .to_string()
            .contains("resource limit"));
        let f = MPoly::monomial(&[2, n], Rat::one());
        let g = MPoly::from_terms(2, vec![(vec![1, n], Rat::one()), (vec![0, 0], Rat::one())]);
        assert!(pdiv(&f, &g, 0)
            .unwrap_err()
            .to_string()
            .contains("exponent overflow"));
    }
    #[test]
    fn invalid_fields_and_zero_divisors_are_errors() {
        let gens = vec![Expr::symbol("x")];
        let parameters = vec![Expr::symbol("a")];
        let a: Polynomial = FractionPoly::from_expr(&parse("a*x+1"), &gens, &parameters, "QQ(a)")
            .unwrap()
            .into();
        let b: Polynomial = Poly::from_expr(&parse("x+1"), &gens)
            .unwrap()
            .set_modulus(5)
            .into();
        assert!(a.pdiv(&b).is_err());
        let zero: Polynomial = Poly::from_expr(&Expr::zero(), &gens).unwrap().into();
        assert!(a.pdiv(&zero).is_err());
    }
    #[test]
    fn rational_function_numerators_obey_the_individual_ceiling() {
        let _scope = Scope::new(1);
        let large = Int::from_big(num_bigint::BigInt::from(1u32) << (MAX_BITS - 1) as usize);
        let c = RF::new(MPoly::constant(large.clone(), 1), MPoly::one(1));
        let g = MPoly::from_terms(1, vec![(vec![1], RF::one()), (vec![0], RF::one())]);
        let cancellation = MPoly::from_terms(1, vec![(vec![2], c.clone()), (vec![1], c.clone())]);
        let (q, r) = pdiv(&cancellation, &g, 0).unwrap();
        assert!(r.is_zero());
        assert_eq!(q, MPoly::monomial(&[1], c.clone()));
        let growth = MPoly::from_terms(
            1,
            vec![
                (vec![2], c.clone()),
                (vec![1], c.neg()),
                (vec![0], RF::one()),
            ],
        );
        assert!(pdiv(&growth, &g, 0)
            .unwrap_err()
            .to_string()
            .contains("coefficient size"));
        let oversized = RF::new(MPoly::constant(large.mul(&Int::from(2)), 1), MPoly::one(1));
        let f = MPoly::constant(oversized, 1);
        assert!(pdiv(&f, &MPoly::one(1), 0)
            .unwrap_err()
            .to_string()
            .contains("coefficient size"));
    }
    #[test]
    fn rational_function_product_storage_is_checked_before_multiplication() {
        let _scope = Scope::new(2);
        let c = Int::from_big(num_bigint::BigInt::from(1u32) << 119_999usize);
        let a = RF::new(
            MPoly::from_terms(2, (0..141).map(|j| (vec![j, 0], c.clone())).collect()),
            MPoly::one(2),
        );
        let b = RF::new(
            MPoly::from_terms(2, (0..141).map(|j| (vec![0, j], c.clone())).collect()),
            MPoly::one(2),
        );
        a.validate_size().unwrap();
        b.validate_size().unwrap();
        assert!(a.weight() < MAX_STORAGE && b.weight() < MAX_STORAGE);
        // Distinct variable axes give 141² distinct output terms. Check only
        // the preflight: even a regression must not allocate that huge product.
        assert!(a
            .preflight(&b, false)
            .unwrap_err()
            .to_string()
            .contains("coefficient storage"));
        assert!(a.preflight(&RF::from_i64(2), false).is_ok());
    }
}
