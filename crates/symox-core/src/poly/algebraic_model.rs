//! Certified primitive generators inside small exact rational algebras.
//!
//! Supported ambient fields are a positive rational cube root, one quadratic
//! field, or the compositum of two independent quadratic fields. Expressions
//! are evaluated with their original radical embeddings, then Krylov linear
//! dependence constructs the minimal polynomial and field membership maps.

use super::convert::{PResult, PolyError};
use super::dense;
use super::int::{Coeff, Int, Rat};
use crate::expr::{Expr, Kind};
use crate::ops::{add, mul, pow};

fn unsupported() -> PolyError {
    PolyError::Unsupported("supported algebraic extensions are positive rational cube roots and composita of at most two quadratic fields".into())
}
fn rational(e: &Expr) -> Option<Rat> {
    e.as_number().and_then(Rat::from_number)
}
fn exact_root(q: &Rat, n: u32) -> Option<Rat> {
    if q.is_negative() {
        return None;
    }
    let a = q.n.to_big().nth_root(n);
    let b = q.d.to_big().nth_root(n);
    if a.pow(n) == q.n.to_big() && b.pow(n) == q.d.to_big() {
        Some(Rat::new(Int::from(a), Int::from(b)))
    } else {
        None
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Radical {
    expression: Expr,
    order: u32,
    radicand: Rat,
}
fn radical_power(base: &Expr, exponent: &Expr) -> PResult<Option<(Radical, i64)>> {
    let Some(exponent) = rational(exponent) else {
        return Ok(None);
    };
    if exponent.d.is_one() {
        return Ok(None);
    }
    let Some(n) = exponent.d.as_i64() else {
        return Err(unsupported());
    };
    let Some(q) = rational(base) else {
        return Err(unsupported());
    };
    if q.is_negative() || q.is_zero() || !matches!(n, 2 | 3) {
        return Err(unsupported());
    }
    let m = exponent.n.as_i64().ok_or_else(unsupported)?;
    Ok(Some((
        Radical {
            expression: pow(base, &Expr::rational(1, n)),
            order: n as u32,
            radicand: q,
        },
        m,
    )))
}
/// Principal positive rational cube roots may be canonicalized as products
/// of integer radicals. Their syntax proves reality; an exact rational cube
/// then identifies a single cubic field, including rational radicands.
fn pure_cubic(e: &Expr) -> Option<Radical> {
    fn syntax(e: &Expr) -> Option<bool> {
        if rational(e).is_some() {
            return Some(false);
        }
        match e.kind() {
            Kind::Add(c, terms) => {
                Rat::from_number(c)?;
                let mut has = false;
                for (t, c) in terms {
                    Rat::from_number(c)?;
                    has |= syntax(t)?;
                }
                Some(has)
            }
            Kind::Mul(c, factors) => {
                Rat::from_number(c)?;
                let mut has = false;
                for (b, x) in factors {
                    has |= power_syntax(b, x)?;
                }
                Some(has)
            }
            Kind::Pow(b, x) => power_syntax(b, x),
            _ => None,
        }
    }
    fn power_syntax(b: &Expr, x: &Expr) -> Option<bool> {
        let exponent = rational(x)?;
        if exponent.d.is_one() {
            return syntax(b);
        }
        let base = rational(b)?;
        (exponent.d == Int::from(3) && !base.is_negative() && !base.is_zero()).then_some(true)
    }
    if !syntax(e)? {
        return None;
    }
    let cube = crate::expand::expand(&e.powi(3));
    let q = rational(&cube)?;
    if q.is_zero() || q.is_negative() {
        return None;
    }
    Some(Radical {
        expression: e.clone(),
        order: 3,
        radicand: q,
    })
}
fn collect(e: &Expr, out: &mut Vec<Radical>) -> PResult<()> {
    if rational(e).is_some() {
        return Ok(());
    }
    if let Some(atom) = pure_cubic(e) {
        if !out.contains(&atom) {
            out.push(atom);
        }
        return Ok(());
    }
    if e == &Expr::i() {
        let radical = Radical {
            expression: e.clone(),
            order: 2,
            radicand: Rat::from_i64(-1),
        };
        if !out.contains(&radical) {
            out.push(radical);
        }
        return Ok(());
    }
    match e.kind() {
        Kind::Add(c, terms) => {
            if Rat::from_number(c).is_none() {
                return Err(unsupported());
            }
            for (t, c) in terms {
                if Rat::from_number(c).is_none() {
                    return Err(unsupported());
                }
                collect(t, out)?;
            }
        }
        Kind::Mul(c, factors) => {
            if Rat::from_number(c).is_none() {
                return Err(unsupported());
            }
            for (b, x) in factors {
                collect_power(b, x, out)?;
            }
        }
        Kind::Pow(b, x) => collect_power(b, x, out)?,
        _ => return Err(unsupported()),
    }
    Ok(())
}
fn collect_power(b: &Expr, x: &Expr, out: &mut Vec<Radical>) -> PResult<()> {
    if let Some((radical, _)) = radical_power(b, x)? {
        if !out.contains(&radical) {
            out.push(radical);
        }
        Ok(())
    } else if x.as_i64().is_some() {
        collect(b, out)
    } else {
        Err(unsupported())
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Tensor {
    bases: Vec<Radical>,
    dimension: usize,
}
impl Tensor {
    fn new(expressions: &[Expr]) -> PResult<Self> {
        let mut atoms = Vec::new();
        for expression in expressions {
            collect(expression, &mut atoms)?;
        }
        atoms.sort_by_cached_key(|a| crate::sortkey::sort_key(&a.expression));
        let mut result = Self {
            bases: Vec::new(),
            dimension: 1,
        };
        for atom in atoms {
            if result.radical(&atom).is_ok() {
                continue;
            }
            if result.same_power_class(&atom) {
                return Err(PolyError::Unsupported(
                    "could not certify the radical embedding".into(),
                ));
            }
            if (!result.bases.is_empty() && result.bases[0].order != atom.order)
                || (atom.order == 3 && !result.bases.is_empty())
                || (atom.order == 2 && result.bases.len() >= 2)
            {
                return Err(unsupported());
            }
            result.dimension *= atom.order as usize;
            result.bases.push(atom);
        }
        Ok(result)
    }
    fn unit(&self, i: usize) -> Vec<Rat> {
        let mut value = vec![Rat::zero(); self.dimension];
        value[i] = Rat::one();
        value
    }
    fn constant(&self, c: Rat) -> Vec<Rat> {
        let mut value = vec![Rat::zero(); self.dimension];
        value[0] = c;
        value
    }
    fn powers(&self, mut i: usize) -> Vec<u32> {
        self.bases
            .iter()
            .map(|b| {
                let power = (i % b.order as usize) as u32;
                i /= b.order as usize;
                power
            })
            .collect()
    }
    fn monomial_expr(&self, i: usize) -> Expr {
        mul(&self
            .powers(i)
            .iter()
            .zip(&self.bases)
            .map(|(k, b)| b.expression.powi(*k as i64))
            .collect::<Vec<_>>())
    }
    fn same_power_class(&self, atom: &Radical) -> bool {
        if self.bases.iter().any(|b| b.order != atom.order) {
            return false;
        }
        (0..self.dimension).any(|i| {
            let q = self
                .powers(i)
                .iter()
                .zip(&self.bases)
                .fold(Rat::one(), |q, (p, b)| q.mul(&b.radicand.pow(*p)));
            exact_root(&atom.radicand.div(&q), atom.order).is_some()
        })
    }
    fn radical(&self, atom: &Radical) -> PResult<Vec<Rat>> {
        if self.bases.iter().any(|b| b.order != atom.order) {
            return Err(unsupported());
        }
        for i in 0..self.dimension {
            let mut q = Rat::one();
            for (p, b) in self.powers(i).iter().zip(&self.bases) {
                q = q.mul(&b.radicand.pow(*p));
            }
            let ratio = atom.radicand.div(&q);
            let Some(c) = exact_root(&ratio, atom.order) else {
                continue;
            };
            let candidate = mul(&[Expr::num(c.to_number()), self.monomial_expr(i)]);
            if candidate == atom.expression {
                return Ok(scale(&self.unit(i), &c));
            }
            if atom.order == 2 && crate::ops::neg(&candidate) == atom.expression {
                return Ok(scale(&self.unit(i), &c.neg()));
            }
        }
        Err(unsupported())
    }
    pub fn multiply(&self, a: &[Rat], b: &[Rat]) -> Vec<Rat> {
        let mut result = vec![Rat::zero(); self.dimension];
        for (i, ai) in a.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
            for (j, bj) in b.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
                let mut coefficient = ai.mul(bj);
                let mut index = 0;
                let mut stride = 1;
                for ((u, v), base) in self.powers(i).iter().zip(self.powers(j)).zip(&self.bases) {
                    let exponent = *u + v;
                    coefficient = coefficient.mul(&base.radicand.pow(exponent / base.order));
                    index += (exponent % base.order) as usize * stride;
                    stride *= base.order as usize;
                }
                result[index] = result[index].add(&coefficient);
            }
        }
        result
    }
    fn inverse(&self, a: &[Rat]) -> PResult<Vec<Rat>> {
        if a.iter().all(Coeff::is_zero) {
            return Err(PolyError::Unsupported(
                "division by zero algebraic coefficient".into(),
            ));
        }
        let columns: Vec<_> = (0..self.dimension)
            .map(|i| self.multiply(a, &self.unit(i)))
            .collect();
        solve_columns(&columns, &self.unit(0))
            .ok_or_else(|| PolyError::Failed("noninvertible ambient algebra coefficient".into()))
    }
    fn power(&self, a: &[Rat], k: i64) -> PResult<Vec<Rat>> {
        let mut base = if k < 0 { self.inverse(a)? } else { a.to_vec() };
        let mut k = k.unsigned_abs();
        let mut result = self.unit(0);
        while k != 0 {
            if k & 1 != 0 {
                result = self.multiply(&result, &base);
            }
            k >>= 1;
            if k != 0 {
                base = self.multiply(&base, &base);
            }
        }
        Ok(result)
    }
    fn eval_power(&self, b: &Expr, x: &Expr) -> PResult<Vec<Rat>> {
        if let Some((atom, m)) = radical_power(b, x)? {
            self.power(&self.radical(&atom)?, m)
        } else {
            self.power(&self.evaluate(b)?, x.as_i64().ok_or_else(unsupported)?)
        }
    }
    pub fn evaluate(&self, e: &Expr) -> PResult<Vec<Rat>> {
        if let Some(c) = rational(e) {
            return Ok(self.constant(c));
        }
        if let Some(atom) = pure_cubic(e) {
            return self.radical(&atom);
        }
        if e == &Expr::i() {
            return self.radical(&Radical {
                expression: e.clone(),
                order: 2,
                radicand: Rat::from_i64(-1),
            });
        }
        match e.kind() {
            Kind::Add(c, terms) => {
                let mut result = self.constant(Rat::from_number(c).ok_or_else(unsupported)?);
                for (t, c) in terms {
                    result = add_vectors(
                        &result,
                        &scale(
                            &self.evaluate(t)?,
                            &Rat::from_number(c).ok_or_else(unsupported)?,
                        ),
                    );
                }
                Ok(result)
            }
            Kind::Mul(c, factors) => {
                let mut result = self.constant(Rat::from_number(c).ok_or_else(unsupported)?);
                for (b, x) in factors {
                    result = self.multiply(&result, &self.eval_power(b, x)?);
                }
                Ok(result)
            }
            Kind::Pow(b, x) => self.eval_power(b, x),
            _ => Err(unsupported()),
        }
    }
}
fn scale(a: &[Rat], c: &Rat) -> Vec<Rat> {
    a.iter().map(|a| a.mul(c)).collect()
}
fn add_vectors(a: &[Rat], b: &[Rat]) -> Vec<Rat> {
    a.iter().zip(b).map(|(a, b)| a.add(b)).collect()
}
/// Solve a rectangular, linearly independent rational-column system. The
/// augmented residual is checked; a vector outside the column space fails.
fn solve_columns(columns: &[Vec<Rat>], value: &[Rat]) -> Option<Vec<Rat>> {
    let width = columns.len();
    let height = value.len();
    if width == 0 {
        return value.iter().all(Coeff::is_zero).then(Vec::new);
    }
    let mut rows: Vec<_> = (0..height)
        .map(|i| {
            columns
                .iter()
                .map(|c| c[i].clone())
                .chain(std::iter::once(value[i].clone()))
                .collect::<Vec<_>>()
        })
        .collect();
    let mut pivots = Vec::new();
    let mut next = 0;
    for col in 0..width {
        let pivot = (next..height).find(|&i| !rows[i][col].is_zero())?;
        rows.swap(pivot, next);
        let inverse = rows[next][col].inv();
        for j in col..=width {
            rows[next][j] = rows[next][j].mul(&inverse);
        }
        for i in 0..height {
            if i == next {
                continue;
            }
            let factor = rows[i][col].clone();
            for j in col..=width {
                rows[i][j] = rows[i][j].sub(&factor.mul(&rows[next][j]));
            }
        }
        pivots.push(next);
        next += 1;
    }
    if rows
        .iter()
        .skip(next)
        .any(|r| r[..width].iter().all(Coeff::is_zero) && !r[width].is_zero())
    {
        return None;
    }
    Some(pivots.into_iter().map(|i| rows[i][width].clone()).collect())
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct PrimitiveModel {
    pub algebra: Tensor,
    pub powers: Vec<Vec<Rat>>,
    pub minimal: Vec<Rat>,
}
impl PrimitiveModel {
    fn from_algebra(algebra: Tensor, generator: &Expr) -> PResult<Self> {
        let alpha = algebra.evaluate(generator)?;
        let mut powers = vec![algebra.unit(0)];
        let mut value = alpha.clone();
        for _ in 1..=algebra.dimension {
            if let Some(coefficients) = solve_columns(&powers, &value) {
                let mut minimal: Vec<_> = coefficients.into_iter().map(|c| c.neg()).collect();
                minimal.push(Rat::one());
                // A vanishing polynomial alone does not define a field.
                // Certify irreducibility in QQ[T] with the exact integer
                // factorizer, independently of the Krylov construction.
                let (_, integer) = dense::clear_denoms(&minimal);
                let (_, factors) = super::factor_uni::factor_zx(&integer);
                if factors.len() != 1 || factors[0].1 != 1 || factors[0].0.len() != minimal.len() {
                    return Err(PolyError::Failed(
                        "constructed algebraic modulus is not irreducible".into(),
                    ));
                }
                let model = Self {
                    algebra,
                    powers,
                    minimal,
                };
                if !model.verify(generator) {
                    return Err(PolyError::Failed(
                        "invalid algebraic generator certificate".into(),
                    ));
                }
                return Ok(model);
            }
            powers.push(value.clone());
            value = algebra.multiply(&value, &alpha);
        }
        Err(PolyError::Failed(
            "could not construct a minimal polynomial".into(),
        ))
    }
    pub fn build(generator: &Expr) -> PResult<Self> {
        Self::from_algebra(Tensor::new(std::slice::from_ref(generator))?, generator)
    }
    pub fn build_many(generators: &[Expr]) -> PResult<(Expr, Self)> {
        Self::build_many_with_weights(generators).map(|(g, m, _)| (g, m))
    }
    pub fn build_many_with_weights(generators: &[Expr]) -> PResult<(Expr, Self, Vec<i64>)> {
        let Some(first) = generators.first() else {
            return Err(unsupported());
        };
        let algebra = Tensor::new(generators)?;
        let mut generator = first.clone();
        let mut model = Self::from_algebra(algebra.clone(), &generator)?;
        let mut accepted = vec![first.clone()];
        let mut weights = vec![0; generators.len()];
        weights[0] = 1;
        for (position, next) in generators.iter().enumerate().skip(1) {
            if model.coefficient(next).is_ok() {
                accepted.push(next.clone());
                continue;
            }
            let mut replacement = None;
            for c in 1..=algebra.dimension * algebra.dimension + 1 {
                let candidate =
                    add(&[generator.clone(), mul(&[Expr::int(c as i64), next.clone()])]);
                let candidate_model = Self::from_algebra(algebra.clone(), &candidate)?;
                if accepted
                    .iter()
                    .chain(std::iter::once(next))
                    .all(|e| candidate_model.coefficient(e).is_ok())
                {
                    replacement = Some((candidate, candidate_model, c as i64));
                    break;
                }
            }
            let Some((new_generator, new_model, coefficient)) = replacement else {
                return Err(PolyError::Failed(
                    "could not certify a primitive element".into(),
                ));
            };
            generator = new_generator;
            model = new_model;
            weights[position] = coefficient;
            accepted.push(next.clone());
        }
        Ok((generator, model, weights))
    }
    pub fn coefficient(&self, e: &Expr) -> PResult<Vec<Rat>> {
        let value = self.algebra.evaluate(e)?;
        let mut coefficients = solve_columns(&self.powers, &value).ok_or_else(|| {
            PolyError::Unsupported("coefficient is outside the selected algebraic subfield".into())
        })?;
        dense::trim(&mut coefficients);
        Ok(coefficients)
    }
    pub fn verify(&self, generator: &Expr) -> bool {
        let Ok(alpha) = self.algebra.evaluate(generator) else {
            return false;
        };
        if self.powers.len() + 1 != self.minimal.len()
            || !self.minimal.last().is_some_and(Coeff::is_one)
        {
            return false;
        }
        let mut value = self.algebra.constant(Rat::zero());
        for c in self.minimal.iter().rev() {
            value = add_vectors(
                &self.algebra.multiply(&value, &alpha),
                &self.algebra.constant(c.clone()),
            );
        }
        if value.iter().any(|c| !c.is_zero()) {
            return false;
        }
        let mut power = self.algebra.unit(0);
        for (i, expected) in self.powers.iter().enumerate() {
            if solve_columns(&self.powers[..i], expected).is_some() {
                return false;
            }
            if &power != expected {
                return false;
            }
            power = self.algebra.multiply(&power, &alpha);
        }
        true
    }
}

/// Collect only numeric radical atoms for extension=True, while leaving
/// polynomial generators and transcendental coefficient domains untouched.
pub fn inferred_generators(expressions: &[Expr]) -> PResult<Vec<Expr>> {
    fn scan(e: &Expr, atoms: &mut Vec<Expr>) -> PResult<()> {
        if e == &Expr::i() {
            if !atoms.contains(e) {
                atoms.push(e.clone());
            }
            return Ok(());
        }
        if let Kind::Pow(b, x) = e.kind() {
            if rational(b).is_some() && rational(x).is_some_and(|r| !r.d.is_one()) {
                let Some((atom, _)) = radical_power(b, x)? else {
                    return Err(unsupported());
                };
                if !atoms.contains(&atom.expression) {
                    atoms.push(atom.expression);
                }
                return Ok(());
            }
        }
        for a in e.args() {
            scan(&a, atoms)?;
        }
        Ok(())
    }
    let mut result = Vec::new();
    for expression in expressions {
        scan(expression, &mut result)?;
    }
    result.sort_by_cached_key(crate::sortkey::sort_key);
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::super::algebraic::{AlgebraicPoly, AlgebraicScope, NumberField, AN};
    use super::super::int::Field;
    use super::*;
    fn expr(s: &str) -> Expr {
        crate::parser::parse_str(s).unwrap()
    }
    #[test]
    fn primitive_and_subfield_minimal_polynomials_are_certified() {
        for (generator, degree) in [
            ("2^(1/3)", 3),
            ("1+2^(1/3)", 3),
            ("1/(1+2^(1/3))", 3),
            ("sqrt(2)+I", 4),
            ("sqrt(2)+sqrt(3)", 4),
            ("(sqrt(2)+I)^2", 2),
            ("1+sqrt(2)", 2),
            ("1/(sqrt(2)+I)", 4),
        ] {
            let generator = expr(generator);
            let field = NumberField::from_generator(&generator).unwrap();
            assert_eq!(field.degree(), degree, "{generator}");
            assert!(field.model.as_ref().is_some_and(|m| m.verify(&generator)));
            let (_, z) = dense::clear_denoms(&field.minimal);
            let (_, factors) = super::super::factor_uni::factor_zx(&z);
            assert_eq!(factors.len(), 1);
            assert_eq!(factors[0].1, 1);
            assert_eq!(factors[0].0.len(), degree + 1);
            let _scope = AlgebraicScope::new(field.clone());
            let a = AN::alpha(field.clone()).add(&AN::from_i64(2));
            assert!(a.mul(&a.inv()).is_one());
            assert!(AN::from_expr(&generator, field.clone())
                .unwrap()
                .sub(&AN::alpha(field))
                .is_zero());
        }
    }
    #[test]
    fn embeddings_choose_actual_radical_branches_and_proper_subfields() {
        let field = NumberField::from_generator(&expr("sqrt(2)+I")).unwrap();
        let u = AN::from_expr(&expr("sqrt(2)"), field.clone()).unwrap();
        let v = AN::from_expr(&Expr::i(), field.clone()).unwrap();
        assert_eq!(
            u.coefficients,
            vec![
                Rat::zero(),
                Rat::new(Int::from(5), Int::from(6)),
                Rat::zero(),
                Rat::new(Int::from(-1), Int::from(6))
            ]
        );
        assert_eq!(
            v.coefficients,
            vec![
                Rat::zero(),
                Rat::new(Int::from(1), Int::from(6)),
                Rat::zero(),
                Rat::new(Int::from(1), Int::from(6))
            ]
        );
        assert_eq!(u.mul(&u).to_rat(), Rat::from_i64(2));
        assert_eq!(v.mul(&v).to_rat(), Rat::from_i64(-1));
        assert!(u.add(&v).sub(&AN::alpha(field.clone())).is_zero());
        assert!(AN::from_expr(&expr("sqrt(3)"), field).is_err());
        let field = NumberField::from_generator(&expr("(sqrt(2)+I)^2")).unwrap();
        assert_eq!(field.degree(), 2);
        assert!(AN::from_expr(&Expr::i(), field.clone()).is_err());
        assert!(AN::from_expr(&expr("sqrt(2)"), field.clone()).is_err());
        assert_eq!(
            AN::from_expr(&expr("I*sqrt(2)"), field)
                .unwrap()
                .coefficients,
            vec![
                Rat::new(Int::from(-1), Int::from(2)),
                Rat::new(Int::from(1), Int::from(2))
            ]
        );
    }
    #[test]
    fn dependent_extensions_cubic_embeddings_and_joint_generators() {
        let field = NumberField::from_generators(&[expr("sqrt(2)"), expr("sqrt(8)")]).unwrap();
        assert_eq!(field.degree(), 2);
        assert!(field.model.is_none());
        let field =
            NumberField::from_generators(&[expr("sqrt(2)"), expr("sqrt(3)"), expr("sqrt(6)")])
                .unwrap();
        assert_eq!(field.degree(), 4);
        assert!(field.model.as_ref().unwrap().verify(&field.generator));
        assert!(AN::from_expr(&expr("sqrt(6)"), field).is_ok());
        let field = NumberField::from_generators(&[expr("2^(1/3)"), expr("4^(1/3)")]).unwrap();
        assert_eq!(field.degree(), 3);
        let a = AN::from_expr(&expr("2^(1/3)"), field.clone()).unwrap();
        let b = AN::from_expr(&expr("4^(1/3)"), field).unwrap();
        assert_eq!(a.mul(&a), b);
        for generators in [
            vec![expr("sqrt(2)"), expr("sqrt(3)"), expr("sqrt(5)")],
            vec![expr("2^(1/3)"), expr("sqrt(2)")],
            vec![expr("(-2)^(1/3)")],
            vec![expr("2^(1/5)")],
            vec![expr("sqrt(1+sqrt(2))")],
        ] {
            assert!(
                NumberField::from_generators(&generators).is_err(),
                "{generators:?}"
            );
        }
    }
    #[test]
    fn resultants_give_exact_norms_and_factorization_in_higher_fields() {
        let x = expr("x");
        for generator in [
            expr("2^(1/3)"),
            expr("sqrt(2)+I"),
            expr("sqrt(2)+sqrt(3)"),
            expr("1+sqrt(2)"),
        ] {
            let field = NumberField::from_generator(&generator).unwrap();
            let _scope = AlgebraicScope::new(field.clone());
            let p = AlgebraicPoly::from_expr(
                &crate::ops::sub(&x, &generator),
                std::slice::from_ref(&x),
                field.clone(),
            )
            .unwrap();
            let norm = p.norm().unwrap();
            assert_eq!(dense::from_mpoly(&norm.rep, 0), field.minimal);
            let target =
                AlgebraicPoly::from_expr(&norm.as_expr(), std::slice::from_ref(&x), field.clone())
                    .unwrap();
            let (c, factors) = target.factor_list().unwrap();
            assert_eq!(
                factors.len(),
                if field.degree() == 3 {
                    2
                } else {
                    field.degree()
                }
            );
            let result = factors.iter().fold(
                target.with_rep(super::super::MPoly::constant(c, 1)),
                |p, (f, k)| p.mul(&f.pow(*k)),
            );
            assert_eq!(result.rep, target.rep);
            for (f, _) in factors {
                assert!(target.div(&f).unwrap().1.rep.is_zero());
            }
        }
    }
}

#[cfg(test)]
mod rational_cubic_tests {
    use super::super::algebraic::{AlgebraicScope, NumberField, AN};
    use super::super::int::Field;
    use super::*;
    #[test]
    fn canonical_rational_cube_roots_are_one_field() {
        for source in ["(2/3)^(1/3)", "(5/7)^(1/3)", "(1/2)^(1/3)"] {
            let generator = crate::parser::parse_str(source).unwrap();
            let field = NumberField::from_generator(&generator).unwrap();
            assert_eq!(field.degree(), 3);
            let _scope = AlgebraicScope::new(field.clone());
            let alpha = AN::alpha(field.clone());
            assert!(alpha.mul(&alpha.inv()).is_one());
            assert_eq!(
                AN::from_expr(&generator.powi(2), field.clone()).unwrap(),
                alpha.mul(&alpha)
            );
            assert!(field.model.as_ref().unwrap().verify(&generator));
        }
    }
}

#[cfg(test)]
mod ideal_certificate_tests {
    use super::super::algebraic::{AlgebraicPoly, AlgebraicScope, NumberField};
    use super::super::groebner::MonomialOrder;
    use super::*;
    #[test]
    fn higher_field_groebner_provenance_and_order_conversion() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        for source in [
            "2^(1/3)",
            "(2/3)^(1/3)",
            "sqrt(2)+I",
            "sqrt(2)+sqrt(3)",
            "(sqrt(2)+I)^2",
        ] {
            let generator = crate::parser::parse_str(source).unwrap();
            let field = NumberField::from_generator(&generator).unwrap();
            let _scope = AlgebraicScope::new(field.clone());
            let equations = [
                crate::ops::sub(
                    &crate::ops::mul(&[generator, x.clone(), y.clone()]),
                    &Expr::one(),
                ),
                crate::ops::sub(&x.powi(2), &y),
            ];
            let input: Vec<_> = equations
                .iter()
                .map(|e| {
                    AlgebraicPoly::from_expr(e, &[x.clone(), y.clone()], field.clone())
                        .unwrap()
                        .rep
                })
                .collect();
            for order in [
                MonomialOrder::Lex,
                MonomialOrder::Grlex,
                MonomialOrder::Grevlex,
            ] {
                let certified = super::super::f5b::f5b_certified(&input, order).unwrap();
                assert!(certified.verify(&input, order), "{source} {order:?}");
                assert_eq!(
                    certified.basis,
                    super::super::fieldalg::buchberger(&input, order)
                );
                let converted =
                    super::super::fieldalg::fglm(&certified.basis, order, MonomialOrder::Lex, 2)
                        .unwrap();
                assert_eq!(
                    converted,
                    super::super::fieldalg::buchberger(&input, MonomialOrder::Lex)
                );
            }
        }
    }
}

#[cfg(test)]
mod closure_tests {
    use super::super::algebraic::{AlgebraicPoly, AlgebraicScope, NumberField, AN};
    use super::super::int::Field;
    use super::*;
    #[test]
    fn canonical_products_scaled_subgenerators_and_expression_roundtrips() {
        let x = Expr::symbol("x");
        for (a, b, c, divisor, square) in [(2, 3, 6, 1, 6), (6, 15, 10, 3, 10)] {
            let u = crate::ops::sqrt(&Expr::int(a));
            let v = crate::ops::sqrt(&Expr::int(b));
            let target = crate::ops::sqrt(&Expr::int(c));
            let field = NumberField::from_generators(&[u.clone(), v.clone()]).unwrap();
            assert_eq!(field.degree(), 4);
            let _scope = AlgebraicScope::new(field.clone());
            let product = AN::from_expr(&u, field.clone())
                .unwrap()
                .mul(&AN::from_expr(&v, field.clone()).unwrap())
                .div(&AN::from_i64(divisor));
            assert_eq!(AN::from_expr(&target, field.clone()).unwrap(), product);
            assert_eq!(
                AN::from_expr(&product.as_expr(), field.clone()).unwrap(),
                product
            );
            let p = AlgebraicPoly::from_expr(
                &crate::ops::add2(&x, &target),
                std::slice::from_ref(&x),
                field.clone(),
            )
            .unwrap();
            assert_eq!(p.convert(&p.as_expr()).unwrap().rep, p.rep);
            let divisor = AlgebraicPoly::from_expr(
                &crate::ops::sub(&x, &target),
                std::slice::from_ref(&x),
                field.clone(),
            )
            .unwrap();
            let input = AlgebraicPoly::from_expr(
                &crate::ops::sub(&x.powi(2), &Expr::int(square)),
                std::slice::from_ref(&x),
                field,
            )
            .unwrap();
            assert!(input.div(&divisor).unwrap().1.rep.is_zero());
        }
    }
}
