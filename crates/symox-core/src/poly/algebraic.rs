//! Exact number fields QQ[T]/(m(T)) and polynomials with algebraic coefficients.
//!
//! The initial public constructors select genuine quadratic fields. Every
//! operation reduces modulo the minimal polynomial; inverses use extended
//! Euclid, never floating point evaluation or an external symbolic engine.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::sync::Arc;

use super::convert::{from_exprs_gens, to_expr, PResult, PolyError};
use super::dense;
use super::fieldalg;
use super::groebner::MonomialOrder;
use super::int::{Coeff, Field, Rat};
use super::mpoly::MPoly;
use super::polyobj::Poly;
use crate::expr::{Expr, Kind};
use crate::ops::{add, mul, pow};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct NumberField {
    pub generator: Expr,
    pub minimal: Vec<Rat>,
    pub radicand: Rat,
}

impl NumberField {
    pub fn quadratic(generator: &Expr) -> PResult<Arc<Self>> {
        if !crate::traverse::free_symbols(generator).is_empty() {
            return Err(PolyError::Unsupported(
                "algebraic extensions must be numeric quadratic generators".into(),
            ));
        }
        let square = crate::expand::expand(&generator.powi(2));
        let d = square
            .as_number()
            .and_then(Rat::from_number)
            .ok_or_else(|| {
                PolyError::Unsupported("only a single quadratic extension is implemented".into())
            })?;
        if generator.is_rational_number() || d.is_zero() {
            return Err(PolyError::Unsupported(
                "a nonrational quadratic generator is required".into(),
            ));
        }
        // A degree-two quotient is a field only when T²-d is irreducible
        // over QQ. A held radical can square to a rational perfect square
        // without its generator node having simplified to a rational number.
        if !d.is_negative() {
            let numerator_root = d.n.isqrt();
            let denominator_root = d.d.isqrt();
            if numerator_root.mul(&numerator_root) == d.n
                && denominator_root.mul(&denominator_root) == d.d
            {
                return Err(PolyError::Unsupported(
                    "the quadratic modulus is reducible over QQ".into(),
                ));
            }
        }
        Ok(Arc::new(Self {
            generator: generator.clone(),
            minimal: vec![d.neg(), Rat::zero(), Rat::one()],
            radicand: d,
        }))
    }
    pub fn display(&self) -> String {
        format!("QQ<{}>", self.generator)
    }
    pub fn parse(domain: &str) -> PResult<Option<Arc<Self>>> {
        if !domain.starts_with("QQ<") {
            return Ok(None);
        }
        if !domain.ends_with('>') {
            return Err(PolyError::Failed(
                "malformed algebraic coefficient domain".into(),
            ));
        }
        let e = crate::parser::parse_str(&domain[3..domain.len() - 1])
            .map_err(|_| PolyError::Failed("invalid algebraic generator".into()))?;
        Self::quadratic(&e).map(Some)
    }
}

thread_local! {static CURRENT:RefCell<Option<Arc<NumberField>>>=const {RefCell::new(None)};}
pub struct AlgebraicScope(Option<Arc<NumberField>>);
impl AlgebraicScope {
    pub fn new(field: Arc<NumberField>) -> Self {
        Self(CURRENT.with(|c| c.replace(Some(field))))
    }
}
impl Drop for AlgebraicScope {
    fn drop(&mut self) {
        CURRENT.with(|c| c.replace(self.0.take()));
    }
}
fn current() -> Arc<NumberField> {
    CURRENT
        .with(|c| c.borrow().clone())
        .expect("algebraic coefficient scope")
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AN {
    pub coefficients: Vec<Rat>,
    pub field: Arc<NumberField>,
}
impl AN {
    pub fn new(mut coefficients: Vec<Rat>, field: Arc<NumberField>) -> Self {
        dense::trim(&mut coefficients);
        if coefficients.len() >= field.minimal.len() {
            coefficients = dense::rem(&coefficients, &field.minimal);
        }
        Self {
            coefficients,
            field,
        }
    }
    pub fn rational(c: Rat, field: Arc<NumberField>) -> Self {
        Self::new(vec![c], field)
    }
    pub fn alpha(field: Arc<NumberField>) -> Self {
        Self::new(vec![Rat::zero(), Rat::one()], field)
    }
    pub fn as_expr(&self) -> Expr {
        add(&self
            .coefficients
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.is_zero())
            .map(|(i, c)| {
                mul(&[
                    Expr::num(c.to_number()),
                    pow(&self.field.generator, &Expr::int(i as i64)),
                ])
            })
            .collect::<Vec<_>>())
    }
    pub fn pow(&self, k: i64) -> Self {
        let mut base = if k < 0 { self.inv() } else { self.clone() };
        let mut k = k.unsigned_abs();
        let mut result = Self::rational(Rat::one(), self.field.clone());
        while k != 0 {
            if k & 1 != 0 {
                result = result.mul(&base);
            }
            k >>= 1;
            if k != 0 {
                base = base.mul(&base);
            }
        }
        result
    }
    pub fn conjugate(&self) -> Self {
        Self::new(
            self.coefficients
                .iter()
                .enumerate()
                .map(|(i, c)| if i % 2 == 0 { c.clone() } else { c.neg() })
                .collect(),
            self.field.clone(),
        )
    }
    pub fn from_expr(e: &Expr, field: Arc<NumberField>) -> PResult<Self> {
        if e == &field.generator {
            return Ok(Self::alpha(field));
        }
        if let Some(r) = e.as_number().and_then(Rat::from_number) {
            return Ok(Self::rational(r, field));
        }
        if let Ok(other) = NumberField::quadratic(e) {
            let ratio = other.radicand.div(&field.radicand);
            let root = crate::ops::sqrt(&Expr::num(ratio.to_number()));
            if let Some(c) = root.as_number().and_then(Rat::from_number) {
                // Choose the embedding by exact expression equality, not by
                // signs inferred from approximate numerical roots.
                let candidate = mul(&[Expr::num(c.to_number()), field.generator.clone()]);
                if candidate == *e {
                    return Ok(Self::alpha(field.clone()).mul(&Self::rational(c, field)));
                }
                if crate::ops::neg(&candidate) == *e {
                    return Ok(Self::alpha(field.clone()).mul(&Self::rational(c.neg(), field)));
                }
            }
        }
        let out = match e.kind() {
            Kind::Add(c, terms) => {
                let Some(c) = Rat::from_number(c) else {
                    return Err(unsupported());
                };
                let mut out = Self::rational(c, field.clone());
                for (t, c) in terms {
                    let Some(c) = Rat::from_number(c) else {
                        return Err(unsupported());
                    };
                    out = out.add(
                        &Self::from_expr(t, field.clone())?.mul(&Self::rational(c, field.clone())),
                    );
                }
                out
            }
            Kind::Mul(c, factors) => {
                let Some(c) = Rat::from_number(c) else {
                    return Err(unsupported());
                };
                let mut out = Self::rational(c, field.clone());
                for (b, x) in factors {
                    let value = if x.is_one() {
                        Self::from_expr(b, field.clone())?
                    } else {
                        Self::from_expr(&pow(b, x), field.clone())?
                    };
                    out = out.mul(&value);
                }
                out
            }
            Kind::Pow(b, x) => {
                let k = x.as_i64().ok_or_else(unsupported)?;
                let base = Self::from_expr(b, field.clone())?;
                if k < 0 && base.is_zero() {
                    return Err(PolyError::Failed(
                        "zero algebraic coefficient denominator".into(),
                    ));
                }
                base.pow(k)
            }
            _ => return Err(unsupported()),
        };
        Ok(out)
    }
}
fn unsupported() -> PolyError {
    PolyError::Unsupported("coefficient is not in the selected quadratic number field".into())
}
impl Coeff for AN {
    fn zero() -> Self {
        Self::new(Vec::new(), current())
    }
    fn one() -> Self {
        Self::rational(Rat::one(), current())
    }
    fn from_i64(i: i64) -> Self {
        Self::rational(Rat::from_i64(i), current())
    }
    fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
    }
    fn is_one(&self) -> bool {
        self.coefficients.len() == 1 && self.coefficients[0].is_one()
    }
    fn add(&self, o: &Self) -> Self {
        Self::new(
            dense::add(&self.coefficients, &o.coefficients),
            self.field.clone(),
        )
    }
    fn sub(&self, o: &Self) -> Self {
        Self::new(
            dense::sub(&self.coefficients, &o.coefficients),
            self.field.clone(),
        )
    }
    fn mul(&self, o: &Self) -> Self {
        Self::new(
            dense::mul(&self.coefficients, &o.coefficients),
            self.field.clone(),
        )
    }
    fn neg(&self) -> Self {
        Self::new(dense::neg(&self.coefficients), self.field.clone())
    }
    fn is_negative(&self) -> bool {
        self.coefficients.last().is_some_and(Rat::is_negative)
    }
    fn div_exact(&self, o: &Self) -> Option<Self> {
        if o.is_zero() {
            None
        } else {
            Some(self.mul(&o.inv()))
        }
    }
    fn to_rat(&self) -> Rat {
        if self.coefficients.is_empty() {
            Rat::zero()
        } else {
            assert_eq!(
                self.coefficients.len(),
                1,
                "nonrational algebraic coefficient"
            );
            self.coefficients[0].clone()
        }
    }
}
impl Field for AN {
    fn inv(&self) -> Self {
        assert!(!self.is_zero(), "division by zero algebraic coefficient");
        let (g, s, _) = dense::gcdex_field(&self.coefficients, &self.field.minimal);
        assert_eq!(g, vec![Rat::one()], "irreducible number field modulus");
        Self::new(s, self.field.clone())
    }
}

pub fn infer(expressions: &[Expr]) -> PResult<Option<Arc<NumberField>>> {
    fn scan(e: &Expr, out: &mut Vec<Expr>) {
        let (_, tail) = crate::ops::as_coeff_mul(e);
        if let Ok(f) = NumberField::quadratic(&tail) {
            if !out.contains(&f.generator) {
                out.push(f.generator.clone());
            }
            return;
        }
        for a in e.args() {
            scan(&a, out);
        }
    }
    let mut candidates = Vec::new();
    for e in expressions {
        scan(e, &mut candidates);
    }
    let Some(generator) = candidates.first() else {
        return Ok(None);
    };
    let field = NumberField::quadratic(generator)?;
    let _scope = AlgebraicScope::new(field.clone());
    for e in candidates {
        AN::from_expr(&e, field.clone())?;
    }
    Ok(Some(field))
}

#[derive(Clone, Debug)]
pub struct AlgebraicPoly {
    pub rep: MPoly<AN>,
    pub gens: Vec<Expr>,
    pub number_field: Arc<NumberField>,
}
impl AlgebraicPoly {
    pub fn from_expr(e: &Expr, gens: &[Expr], field: Arc<NumberField>) -> PResult<Self> {
        if gens.is_empty()
            || gens
                .iter()
                .enumerate()
                .any(|(i, g)| g.is_number() || gens[..i].contains(g))
        {
            return Err(PolyError::Failed(
                "polynomial generators must be distinct nonnumeric expressions".into(),
            ));
        }
        let _scope = AlgebraicScope::new(field.clone());
        let c = from_exprs_gens(std::slice::from_ref(e), gens)?;
        if c.float {
            return Err(PolyError::Unsupported(
                "algebraic domains require exact coefficients".into(),
            ));
        }
        let mut grouped: BTreeMap<Vec<u32>, Vec<(Vec<u32>, Rat)>> = BTreeMap::new();
        for (e, c) in c.polys[0].iter() {
            grouped
                .entry(e[..gens.len()].to_vec())
                .or_default()
                .push((e[gens.len()..].to_vec(), c.clone()));
        }
        let mut terms = Vec::new();
        for (e, ts) in grouped {
            let value = to_expr(
                &MPoly::from_terms(c.gens.len() - gens.len(), ts),
                &c.gens[gens.len()..],
            );
            terms.push((e, AN::from_expr(&value, field.clone())?));
        }
        Ok(Self {
            rep: MPoly::from_terms(gens.len(), terms),
            gens: gens.to_vec(),
            number_field: field,
        })
    }
    pub fn convert(&self, e: &Expr) -> PResult<Self> {
        Self::from_expr(e, &self.gens, self.number_field.clone())
    }
    pub fn with_rep(&self, rep: MPoly<AN>) -> Self {
        Self {
            rep,
            ..self.clone()
        }
    }
    pub fn as_expr(&self) -> Expr {
        add(&self
            .rep
            .iter()
            .map(|(e, c)| {
                let mut v = vec![c.as_expr()];
                for (g, &k) in self.gens.iter().zip(e) {
                    if k != 0 {
                        v.push(pow(g, &Expr::int(k as i64)));
                    }
                }
                mul(&v)
            })
            .collect::<Vec<_>>())
    }
    pub fn terms(&self) -> Vec<(Vec<u32>, Expr)> {
        self.rep
            .iter()
            .map(|(e, c)| (e.to_vec(), c.as_expr()))
            .collect()
    }
    pub fn legacy(&self) -> Poly {
        Poly::from_expr(&self.as_expr(), &self.gens).expect("algebraic polynomial view")
    }
    pub fn add(&self, o: &Self) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        self.with_rep(self.rep.add(&o.rep))
    }
    pub fn sub(&self, o: &Self) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        self.with_rep(self.rep.sub(&o.rep))
    }
    pub fn mul(&self, o: &Self) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        self.with_rep(self.rep.mul(&o.rep))
    }
    pub fn neg(&self) -> Self {
        self.with_rep(self.rep.neg())
    }
    pub fn pow(&self, k: u32) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        self.with_rep(self.rep.pow(k))
    }
    pub fn monic(&self) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        self.with_rep(fieldalg::monic(&self.rep, MonomialOrder::Lex))
    }
    pub fn div(&self, o: &Self) -> PResult<(Self, Self)> {
        if o.rep.is_zero() {
            return Err(PolyError::Unsupported("polynomial division by zero".into()));
        }
        let _s = AlgebraicScope::new(self.number_field.clone());
        let (q, r) = fieldalg::reduce(&self.rep, std::slice::from_ref(&o.rep), MonomialOrder::Lex);
        Ok((
            self.with_rep(q.into_iter().next().unwrap()),
            self.with_rep(r),
        ))
    }
    pub fn diff(&self, i: usize) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        self.with_rep(self.rep.derivative(i))
    }
    pub fn integrate(&self, i: usize) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        let ts = self
            .rep
            .iter()
            .map(|(e, c)| {
                let mut e = e.to_vec();
                e[i] += 1;
                let c = c.div(&AN::from_i64(e[i] as i64));
                (e, c)
            })
            .collect();
        self.with_rep(MPoly::from_terms(self.gens.len(), ts))
    }
    pub fn compose(&self, o: &Self) -> Self {
        let _s = AlgebraicScope::new(self.number_field.clone());
        let coefficients = dense::from_mpoly(&self.rep, 0);
        let mut result = self.with_rep(MPoly::zero(1));
        for c in coefficients.iter().rev() {
            result = result
                .mul(o)
                .add(&self.with_rep(MPoly::constant(c.clone(), 1)));
        }
        result
    }
    pub fn gcd(&self, o: &Self) -> Self {
        assert_eq!(self.gens.len(), 1, "algebraic gcd requires one generator");
        let _s = AlgebraicScope::new(self.number_field.clone());
        let a = dense::from_mpoly(&self.rep, 0);
        let b = dense::from_mpoly(&o.rep, 0);
        self.with_rep(dense::to_mpoly(&dense::gcd_field(&a, &b), 0, 1))
    }
    pub fn resultant(&self, other: &Self) -> PResult<AN> {
        if self.gens.len() != 1 {
            return Err(PolyError::Unsupported(
                "algebraic resultants require one generator".into(),
            ));
        }
        let _scope = AlgebraicScope::new(self.number_field.clone());
        if self.rep.is_zero() || other.rep.is_zero() {
            return Ok(AN::zero());
        }
        let a = dense::from_mpoly(&self.rep, 0);
        let b = dense::from_mpoly(&other.rep, 0);
        let n = a.len() - 1;
        let m = b.len() - 1;
        let size = n + m;
        let mut matrix = vec![vec![AN::zero(); size]; size];
        for row in 0..m {
            for (column, c) in a.iter().rev().enumerate() {
                matrix[row][row + column] = c.clone();
            }
        }
        for row in 0..n {
            for (column, c) in b.iter().rev().enumerate() {
                matrix[m + row][row + column] = c.clone();
            }
        }
        let mut determinant = AN::one();
        for col in 0..size {
            let Some(pivot) = (col..size).find(|&row| !matrix[row][col].is_zero()) else {
                return Ok(AN::zero());
            };
            if pivot != col {
                matrix.swap(pivot, col);
                determinant = determinant.neg();
            }
            let diagonal = matrix[col][col].clone();
            determinant = determinant.mul(&diagonal);
            for row in col + 1..size {
                let factor = matrix[row][col].div(&diagonal);
                for j in col + 1..size {
                    matrix[row][j] = matrix[row][j].sub(&factor.mul(&matrix[col][j]));
                }
            }
        }
        Ok(determinant)
    }
    pub fn discriminant(&self) -> PResult<AN> {
        let _scope = AlgebraicScope::new(self.number_field.clone());
        let degree = self.rep.degree(0);
        if degree < 1 {
            return Ok(AN::zero());
        }
        let value = self.resultant(&self.diff(0))?.div(&self.rep.lc());
        Ok(if (degree * (degree - 1) / 2) % 2 == 0 {
            value
        } else {
            value.neg()
        })
    }
    pub fn sqf_list(&self) -> (AN, Vec<(Self, u32)>) {
        let _s = AlgebraicScope::new(self.number_field.clone());
        if self.rep.is_zero() {
            return (AN::zero(), Vec::new());
        }
        let leading = self.rep.lc();
        let f = dense::monic(&dense::from_mpoly(&self.rep, 0));
        let mut c = dense::gcd_field(&f, &dense::derivative(&f));
        let mut w = dense::divrem(&f, &c).0;
        let mut out = Vec::new();
        let mut i = 1;
        while w.len() > 1 {
            let y = dense::gcd_field(&w, &c);
            let z = dense::divrem(&w, &y).0;
            if z.len() > 1 {
                out.push((self.with_rep(dense::to_mpoly(&z, 0, 1)), i));
            }
            w = y;
            c = dense::divrem(&c, &w).0;
            i += 1;
        }
        (leading, out)
    }
    pub fn factor_list(&self) -> PResult<(AN, Vec<(Self, u32)>)> {
        if self.gens.len() != 1 {
            return Err(PolyError::Unsupported(
                "multivariate factorization over number fields is not yet implemented".into(),
            ));
        }
        let _s = AlgebraicScope::new(self.number_field.clone());
        let (c, parts) = self.sqf_list();
        let mut output = Vec::new();
        for (f, k) in parts {
            for factor in f.factor_squarefree()? {
                output.push((factor, k));
            }
        }
        output.sort_by_cached_key(|(p, k)| {
            (p.rep.degree(0), *k, crate::sortkey::sort_key(&p.as_expr()))
        });
        Ok((c, output))
    }
    fn factor_squarefree(&self) -> PResult<Vec<Self>> {
        let field = self.number_field.clone();
        let f = dense::from_mpoly(&self.rep, 0);
        let shift = |p: &[AN], a: &AN| {
            let mut out = Vec::new();
            for c in p.iter().rev() {
                out = dense::mul(&out, &[a.clone(), AN::one()]);
                out = dense::add(&out, std::slice::from_ref(c));
            }
            out
        };
        let alpha = AN::alpha(field.clone());
        for s in 0..64 {
            let delta = alpha.mul(&AN::from_i64(s));
            let shifted = shift(&f, &delta.neg());
            let conjugate: Vec<_> = shifted.iter().map(AN::conjugate).collect();
            let norm = dense::mul(&shifted, &conjugate);
            if norm.iter().any(|a| a.coefficients.len() > 1) {
                return Err(PolyError::Failed("nonrational algebraic norm".into()));
            }
            let norm: Vec<_> = norm.iter().map(Coeff::to_rat).collect();
            if dense::gcd_field(&norm, &dense::derivative(&norm)).len() > 1 {
                continue;
            }
            let (_, z) = dense::clear_denoms(&norm);
            let (_, factors) = super::factor_uni::factor_zx(&z);
            let mut remaining = f.clone();
            let mut out = Vec::new();
            for (g, _) in factors {
                let g: Vec<_> = g
                    .into_iter()
                    .map(|c| AN::rational(Rat::from(c), field.clone()))
                    .collect();
                let lifted = shift(&g, &delta);
                let divisor = dense::gcd_field(&remaining, &lifted);
                if divisor.len() > 1 {
                    remaining = dense::divrem(&remaining, &divisor).0;
                    out.push(self.with_rep(dense::to_mpoly(&divisor, 0, 1)));
                }
            }
            if remaining.len() > 1 {
                return Err(PolyError::Failed(
                    "incomplete algebraic norm factorization".into(),
                ));
            }
            return Ok(out);
        }
        Err(PolyError::Unsupported(
            "could not separate conjugate roots for exact algebraic factorization".into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reject_held_rational_radicals_with_reducible_modulus() {
        for (base, exponent) in [
            (Expr::int(4), Expr::rational(1, 2)),
            (Expr::rational(4, 9), Expr::rational(1, 2)),
        ] {
            let generator = Expr::from_unevaluated_kind(Kind::Pow(base, exponent));
            assert!(!generator.is_rational_number());
            assert!(NumberField::quadratic(&generator).is_err());
        }
        assert!(NumberField::quadratic(&crate::ops::sqrt(&Expr::int(2))).is_ok());
        assert!(NumberField::quadratic(&Expr::i()).is_ok());
    }
    #[test]
    fn quotient_inverses_and_factorization() {
        for generator in [crate::ops::sqrt(&Expr::int(2)), Expr::i()] {
            let field = NumberField::quadratic(&generator).unwrap();
            let _s = AlgebraicScope::new(field.clone());
            let alpha = AN::alpha(field.clone());
            assert_eq!(alpha.mul(&alpha).to_rat(), field.radicand);
            let a = alpha.add(&AN::from_i64(1));
            assert!(a.mul(&a.inv()).is_one());
            let x = Expr::symbol("x");
            let expression = crate::ops::sub(&x.powi(2), &Expr::num(field.radicand.to_number()));
            let p = AlgebraicPoly::from_expr(&expression, &[x], field).unwrap();
            let (c, fs) = p.factor_list().unwrap();
            assert_eq!(fs.len(), 2);
            let result = fs
                .into_iter()
                .fold(p.with_rep(MPoly::constant(c, 1)), |a, (b, k)| {
                    a.mul(&b.pow(k))
                });
            assert_eq!(result.rep, p.rep);
        }
    }
}

#[cfg(test)]
mod scope_tests {
    use super::*;
    #[test]
    fn nested_and_thread_local_number_fields_restore_context() {
        let a = NumberField::quadratic(&crate::ops::sqrt(&Expr::int(2))).unwrap();
        let b = NumberField::quadratic(&Expr::i()).unwrap();
        let _outer = AlgebraicScope::new(a.clone());
        assert_eq!(AN::one().field, a);
        {
            let _inner = AlgebraicScope::new(b.clone());
            assert_eq!(AN::one().field, b);
        }
        assert_eq!(AN::one().field, a);
        let worker = std::thread::spawn(move || {
            let _s = AlgebraicScope::new(b.clone());
            assert_eq!(AN::one().field, b);
        });
        worker.join().unwrap();
        assert_eq!(AN::one().field, a);
        let result = std::panic::catch_unwind(|| {
            let _s = AlgebraicScope::new(NumberField::quadratic(&Expr::i()).unwrap());
            panic!("scope restoration");
        });
        assert!(result.is_err());
        assert_eq!(AN::one().field, a);
    }
    #[test]
    fn exact_inverse_reduction_and_differing_embeddings() {
        let a = NumberField::quadratic(&crate::ops::sqrt(&Expr::int(2))).unwrap();
        let _s = AlgebraicScope::new(a.clone());
        let alpha = AN::alpha(a.clone());
        let one = AN::one();
        let rationalized = AN::from_expr(
            &crate::ops::div(&Expr::one(), &crate::ops::add2(&Expr::one(), &a.generator)),
            a.clone(),
        )
        .unwrap();
        assert_eq!(rationalized, alpha.sub(&one));
        assert_eq!(
            AN::from_expr(&crate::ops::sqrt(&Expr::int(8)), a.clone()).unwrap(),
            alpha.mul(&AN::from_i64(2))
        );
        assert!(AN::from_expr(&crate::ops::sqrt(&Expr::int(3)), a).is_err());
    }
}
