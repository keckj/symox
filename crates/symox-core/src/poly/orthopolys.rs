//! Orthogonal polynomial sequences (`sympy.polys.orthopolys`).
//!
//! The dense three-term recurrences follow SymPy's `dup_*` generators on
//! descending coefficient lists. The ground domain is chosen like SymPy's
//! `construct_domain(params, field=True)`: rationals (`QQ`), machine floats
//! (`RR`) or the fraction field `ZZ(a, b, ...)` of the symbolic parameters.
//! Any other parameter (algebraic numbers, transcendental constants, mixed
//! floats and symbols) is handled by exact expression arithmetic with
//! cancellation, and the domain is left to the caller's inference.
use super::convert::{number_to_rat, PResult, PolyError};
use super::fraction::{coefficient, coefficient_expr, infer_parameters};
use super::int::{Field, Rat};
use super::ratfunc::{Scope, RF};
use crate::expr::{Expr, Kind};
use crate::number::Number;

/// The polynomial families generated here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Family {
    Jacobi,
    Gegenbauer,
    ChebyshevT,
    ChebyshevU,
    Hermite,
    HermiteProb,
    Legendre,
    Laguerre,
    /// `fn(n, z)` of the spherical Bessel functions, in powers of `1/z`.
    SphericalBesselFn,
    /// `fn(-n, z)`.
    SphericalBesselFnMinus,
}

impl Family {
    pub fn parse(name: &str) -> Option<Self> {
        Some(match name {
            "jacobi" => Self::Jacobi,
            "gegenbauer" => Self::Gegenbauer,
            "chebyshevt" => Self::ChebyshevT,
            "chebyshevu" => Self::ChebyshevU,
            "hermite" => Self::Hermite,
            "hermite_prob" => Self::HermiteProb,
            "legendre" => Self::Legendre,
            "laguerre" => Self::Laguerre,
            "spherical_bessel_fn" => Self::SphericalBesselFn,
            "spherical_bessel_fn_minus" => Self::SphericalBesselFnMinus,
            _ => return None,
        })
    }
    fn arity(self) -> usize {
        match self {
            Self::Jacobi => 2,
            Self::Gegenbauer | Self::Laguerre => 1,
            _ => 0,
        }
    }
    /// Fixed ground domain of the parameter-free families (SymPy's `K`).
    fn fixed_domain(self) -> Option<&'static str> {
        match self {
            Self::Legendre => Some("QQ"),
            Self::Jacobi | Self::Gegenbauer | Self::Laguerre => None,
            _ => Some("ZZ"),
        }
    }
}

/// Descending coefficients and the SymPy domain name, when it is known.
#[derive(Clone, Debug)]
pub struct Named {
    pub coeffs: Vec<Expr>,
    pub domain: Option<String>,
}

/// Arithmetic of a coefficient domain.
trait Ground {
    type E: Clone;
    fn int(&self, i: i64) -> Self::E;
    fn add(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn sub(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn mul(&self, a: &Self::E, b: &Self::E) -> Self::E;
    fn div(&self, a: &Self::E, b: &Self::E) -> PResult<Self::E>;
    fn is_zero(&self, a: &Self::E) -> bool;
}

struct Exact<C>(std::marker::PhantomData<C>);
impl<C: Field> Ground for Exact<C> {
    type E = C;
    fn int(&self, i: i64) -> C {
        C::from_i64(i)
    }
    fn add(&self, a: &C, b: &C) -> C {
        a.add(b)
    }
    fn sub(&self, a: &C, b: &C) -> C {
        a.sub(b)
    }
    fn mul(&self, a: &C, b: &C) -> C {
        a.mul(b)
    }
    fn div(&self, a: &C, b: &C) -> PResult<C> {
        if b.is_zero() {
            return Err(PolyError::Unsupported("integer division or modulo by zero".into()));
        }
        Ok(a.div(b))
    }
    fn is_zero(&self, a: &C) -> bool {
        a.is_zero()
    }
}

/// `RR`: IEEE double arithmetic is SymPy's 53-bit `RealField`.
struct Real;
impl Ground for Real {
    type E = f64;
    fn int(&self, i: i64) -> f64 {
        i as f64
    }
    fn add(&self, a: &f64, b: &f64) -> f64 {
        a + b
    }
    fn sub(&self, a: &f64, b: &f64) -> f64 {
        a - b
    }
    fn mul(&self, a: &f64, b: &f64) -> f64 {
        a * b
    }
    fn div(&self, a: &f64, b: &f64) -> PResult<f64> {
        if *b == 0.0 {
            return Err(PolyError::Unsupported("integer division or modulo by zero".into()));
        }
        Ok(a / b)
    }
    fn is_zero(&self, a: &f64) -> bool {
        *a == 0.0
    }
}

/// Exact expression arithmetic, kept canceled (the `EX` fallback).
struct Symbolic;
impl Symbolic {
    fn norm(&self, e: Expr) -> Expr {
        let c = super::api::cancel(&e).unwrap_or(e);
        let (_, den) = crate::simplify::numer_denom(&c);
        if crate::traverse::free_symbols(&den).is_empty() {
            crate::expand::expand(&c)
        } else {
            c
        }
    }
}
impl Ground for Symbolic {
    type E = Expr;
    fn int(&self, i: i64) -> Expr {
        Expr::int(i)
    }
    fn add(&self, a: &Expr, b: &Expr) -> Expr {
        self.norm(crate::ops::add2(a, b))
    }
    fn sub(&self, a: &Expr, b: &Expr) -> Expr {
        self.norm(crate::ops::sub(a, b))
    }
    fn mul(&self, a: &Expr, b: &Expr) -> Expr {
        self.norm(crate::ops::mul2(a, b))
    }
    fn div(&self, a: &Expr, b: &Expr) -> PResult<Expr> {
        if b.is_zero() {
            return Err(PolyError::Unsupported("integer division or modulo by zero".into()));
        }
        Ok(self.norm(crate::ops::div(a, b)))
    }
    fn is_zero(&self, a: &Expr) -> bool {
        a.is_zero()
    }
}

// Dense descending-list helpers mirroring SymPy's densearith.
fn strip<K: Ground>(k: &K, mut f: Vec<K::E>) -> Vec<K::E> {
    let lead = f.iter().take_while(|c| k.is_zero(c)).count();
    f.drain(..lead);
    f
}
fn mul_ground<K: Ground>(k: &K, f: &[K::E], c: &K::E) -> Vec<K::E> {
    if k.is_zero(c) || f.is_empty() {
        return Vec::new();
    }
    f.iter().map(|a| k.mul(a, c)).collect()
}
fn lshift<K: Ground>(k: &K, f: &[K::E]) -> Vec<K::E> {
    if f.is_empty() {
        return Vec::new();
    }
    let mut out = f.to_vec();
    out.push(k.int(0));
    out
}
fn combine<K: Ground>(k: &K, f: &[K::E], g: &[K::E], subtract: bool) -> Vec<K::E> {
    let n = f.len().max(g.len());
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let a = (i + f.len()).checked_sub(n).map(|j| &f[j]);
        let b = (i + g.len()).checked_sub(n).map(|j| &g[j]);
        out.push(match (a, b) {
            (Some(a), Some(b)) if subtract => k.sub(a, b),
            (Some(a), Some(b)) => k.add(a, b),
            (Some(a), None) => a.clone(),
            (None, Some(b)) if subtract => k.sub(&k.int(0), b),
            (None, Some(b)) => b.clone(),
            (None, None) => unreachable!(),
        });
    }
    strip(k, out)
}
fn sub<K: Ground>(k: &K, f: &[K::E], g: &[K::E]) -> Vec<K::E> {
    combine(k, f, g, true)
}
fn add<K: Ground>(k: &K, f: &[K::E], g: &[K::E]) -> Vec<K::E> {
    combine(k, f, g, false)
}
fn mul<K: Ground>(k: &K, f: &[K::E], g: &[K::E]) -> Vec<K::E> {
    if f.is_empty() || g.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<Option<K::E>> = vec![None; f.len() + g.len() - 1];
    for (i, a) in f.iter().enumerate() {
        for (j, b) in g.iter().enumerate() {
            let t = k.mul(a, b);
            out[i + j] = Some(match out[i + j].take() {
                Some(s) => k.add(&s, &t),
                None => t,
            });
        }
    }
    strip(k, out.into_iter().map(|c| c.unwrap_or_else(|| k.int(0))).collect())
}

fn jacobi<K: Ground>(k: &K, n: usize, a: &K::E, b: &K::E) -> PResult<Vec<K::E>> {
    if n < 1 {
        return Ok(vec![k.int(1)]);
    }
    let (one, two) = (k.int(1), k.int(2));
    let ab = k.add(a, b);
    let mut m2 = vec![one.clone()];
    let mut m1 = vec![k.add(&k.div(&ab, &two)?, &one), k.div(&k.sub(a, b), &two)?];
    for i in 2..=n {
        let i = k.int(i as i64);
        let two_i = k.mul(&two, &i);
        let s = k.add(&ab, &two_i);
        let den = k.mul(&k.mul(&i, &k.add(&ab, &i)), &k.sub(&s, &two));
        let sm1 = k.sub(&s, &one);
        let sm2 = k.sub(&s, &two);
        let f0 = k.div(&k.mul(&sm1, &k.sub(&k.mul(a, a), &k.mul(b, b))), &k.mul(&two, &den))?;
        let f1 = k.div(&k.mul(&k.mul(&sm1, &sm2), &s), &k.mul(&two, &den))?;
        let f2 = k.div(
            &k.mul(&k.mul(&k.sub(&k.add(a, &i), &one), &k.sub(&k.add(b, &i), &one)), &s),
            &den,
        )?;
        let p0 = mul_ground(k, &m1, &f0);
        let p1 = mul_ground(k, &lshift(k, &m1), &f1);
        let p2 = mul_ground(k, &m2, &f2);
        let next = sub(k, &add(k, &p0, &p1), &p2);
        m2 = std::mem::replace(&mut m1, next);
    }
    Ok(m1)
}

fn gegenbauer<K: Ground>(k: &K, n: usize, a: &K::E) -> PResult<Vec<K::E>> {
    if n < 1 {
        return Ok(vec![k.int(1)]);
    }
    let (one, two) = (k.int(1), k.int(2));
    let mut m2 = vec![one.clone()];
    let mut m1 = vec![k.mul(&two, a), k.int(0)];
    let am1 = k.mul(&two, &k.sub(a, &one));
    for i in 2..=n {
        let ratio = k.div(&am1, &k.int(i as i64))?;
        let p1 = mul_ground(k, &lshift(k, &m1), &k.add(&ratio, &two));
        let p2 = mul_ground(k, &m2, &k.add(&ratio, &one));
        let next = sub(k, &p1, &p2);
        m2 = std::mem::replace(&mut m1, next);
    }
    Ok(m1)
}

fn laguerre<K: Ground>(k: &K, n: usize, alpha: &K::E) -> PResult<Vec<K::E>> {
    let one = k.int(1);
    let mut m2 = vec![k.int(0)];
    let mut m1 = vec![one.clone()];
    let am1 = k.sub(alpha, &one);
    for i in 1..=n {
        let i = k.int(i as i64);
        let ratio = k.div(&am1, &i)?;
        let factor = [k.div(&k.sub(&k.int(0), &one), &i)?, k.add(&ratio, &k.int(2))];
        let a = mul(k, &m1, &factor);
        let b = mul_ground(k, &m2, &k.add(&ratio, &one));
        let next = sub(k, &a, &b);
        m2 = std::mem::replace(&mut m1, next);
    }
    Ok(m1)
}

/// `T_n`, `U_n`, `H_n`, `He_n`, `P_n` and the spherical Bessel `fn`.
fn integral<K: Ground>(k: &K, family: Family, n: usize) -> PResult<Vec<K::E>> {
    use Family::*;
    let one = || k.int(1);
    Ok(match family {
        ChebyshevT | ChebyshevU => {
            if n < 1 {
                return Ok(vec![one()]);
            }
            let lead = if family == ChebyshevT { 1 } else { 2 };
            let mut m2 = vec![one()];
            let mut m1 = vec![k.int(lead), k.int(0)];
            for _ in 2..=n {
                let next = sub(k, &mul_ground(k, &lshift(k, &m1), &k.int(2)), &m2);
                m2 = std::mem::replace(&mut m1, next);
            }
            m1
        }
        Hermite | HermiteProb => {
            if n < 1 {
                return Ok(vec![one()]);
            }
            let lead = if family == Hermite { 2 } else { 1 };
            let mut m2 = vec![one()];
            let mut m1 = vec![k.int(lead), k.int(0)];
            for i in 2..=n {
                let d = sub(k, &lshift(k, &m1), &mul_ground(k, &m2, &k.int(i as i64 - 1)));
                let next = if family == Hermite { mul_ground(k, &d, &k.int(2)) } else { d };
                m2 = std::mem::replace(&mut m1, next);
            }
            m1
        }
        Legendre => {
            if n < 1 {
                return Ok(vec![one()]);
            }
            let mut m2 = vec![one()];
            let mut m1 = vec![one(), k.int(0)];
            for i in 2..=n {
                let ii = k.int(i as i64);
                let a = mul_ground(k, &lshift(k, &m1), &k.div(&k.int(2 * i as i64 - 1), &ii)?);
                let b = mul_ground(k, &m2, &k.div(&k.int(i as i64 - 1), &ii)?);
                let next = sub(k, &a, &b);
                m2 = std::mem::replace(&mut m1, next);
            }
            m1
        }
        SphericalBesselFn => {
            if n < 1 {
                return Ok(vec![one(), k.int(0)]);
            }
            let mut m2 = vec![one()];
            let mut m1 = vec![one(), k.int(0)];
            for i in 2..=n {
                let next = sub(k, &mul_ground(k, &lshift(k, &m1), &k.int(2 * i as i64 - 1)), &m2);
                m2 = std::mem::replace(&mut m1, next);
            }
            lshift(k, &m1)
        }
        SphericalBesselFnMinus => {
            let mut m2 = vec![one(), k.int(0)];
            let mut m1 = vec![k.int(0)];
            for i in 2..=n {
                let next = sub(k, &mul_ground(k, &lshift(k, &m1), &k.int(3 - 2 * i as i64)), &m2);
                m2 = std::mem::replace(&mut m1, next);
            }
            m1
        }
        Jacobi | Gegenbauer | Laguerre => unreachable!(),
    })
}

fn generate<K: Ground>(k: &K, family: Family, n: usize, params: &[K::E]) -> PResult<Vec<K::E>> {
    match family {
        Family::Jacobi => jacobi(k, n, &params[0], &params[1]),
        Family::Gegenbauer => gegenbauer(k, n, &params[0]),
        Family::Laguerre => laguerre(k, n, &params[0]),
        _ => integral(k, family, n),
    }
}

/// Upper bound on the index, keeping the quadratic recurrences bounded.
pub const MAX_INDEX: usize = 1 << 20;

/// Coefficients (descending) of the `n`-th member of `family` with the given
/// parameters, together with SymPy's ground domain name when it is known.
pub fn named_poly(family: Family, n: usize, params: &[Expr]) -> PResult<Named> {
    if params.len() != family.arity() {
        return Err(PolyError::Failed(format!(
            "expected {} parameters, got {}",
            family.arity(),
            params.len()
        )));
    }
    if n > MAX_INDEX {
        return Err(PolyError::Unsupported(format!(
            "polynomial index {n} exceeds the supported bound {MAX_INDEX}"
        )));
    }
    if let Some(domain) = family.fixed_domain() {
        let k = Exact::<Rat>(Default::default());
        let coeffs = integral(&k, family, n)?;
        return Ok(Named {
            coeffs: coeffs.iter().map(|c| Expr::num(c.to_number())).collect(),
            domain: Some(domain.into()),
        });
    }
    let numbers: Option<Vec<&Number>> = params.iter().map(Expr::as_number).collect();
    if let Some(numbers) = &numbers {
        if numbers.iter().all(|v| !v.is_float()) {
            let k = Exact::<Rat>(Default::default());
            let values: Vec<Rat> = numbers.iter().map(|v| number_to_rat(v)).collect();
            let coeffs = generate(&k, family, n, &values)?;
            return Ok(Named {
                coeffs: coeffs.iter().map(|c| Expr::num(c.to_number())).collect(),
                domain: Some("QQ".into()),
            });
        }
        if numbers.iter().all(|v| v.float_precision().is_none_or(|p| p == 53)) {
            let values: Vec<f64> = numbers
                .iter()
                .map(|v| match v {
                    Number::Float(f) => *f,
                    other => number_to_rat(other).to_f64(),
                })
                .collect();
            let coeffs = generate(&Real, family, n, &values)?;
            return Ok(Named {
                coeffs: coeffs.into_iter().map(|c| Expr::num(Number::float(c))).collect(),
                domain: Some("RR".into()),
            });
        }
    }
    let floating = params.iter().any(|p| p.float_info().max_precision.is_some());
    let symbols = infer_parameters(params, &[]);
    if !floating && !symbols.is_empty() && symbols.iter().all(|s| matches!(s.kind(), Kind::Sym(..))) {
        let _scope = Scope::new(symbols.len());
        if let Ok(values) = params.iter().map(|p| coefficient(p, &symbols)).collect::<PResult<Vec<RF>>>() {
            let k = Exact::<RF>(Default::default());
            let coeffs = generate(&k, family, n, &values)?;
            let names: Vec<String> = symbols.iter().map(crate::printer::to_str).collect();
            return Ok(Named {
                coeffs: coeffs.iter().map(|c| coefficient_expr(c, &symbols)).collect(),
                domain: Some(format!("ZZ({})", names.join(","))),
            });
        }
    }
    let coeffs = generate(&Symbolic, family, n, params)?;
    Ok(Named { coeffs, domain: None })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn show(f: Family, n: usize, params: &[&str]) -> (Vec<String>, Option<String>) {
        let params: Vec<Expr> = params.iter().map(|p| parse_str(p).unwrap()).collect();
        let r = named_poly(f, n, &params).unwrap();
        (r.coeffs.iter().map(crate::printer::to_str).collect(), r.domain)
    }

    #[test]
    fn integral_families_match_sympy() {
        assert_eq!(show(Family::ChebyshevT, 4, &[]).0, ["8", "0", "-8", "0", "1"]);
        assert_eq!(show(Family::ChebyshevU, 3, &[]).0, ["8", "0", "-4", "0"]);
        assert_eq!(show(Family::Hermite, 3, &[]).0, ["8", "0", "-12", "0"]);
        assert_eq!(show(Family::HermiteProb, 4, &[]).0, ["1", "0", "-6", "0", "3"]);
        assert_eq!(show(Family::Legendre, 2, &[]), (vec!["3/2".into(), "0".into(), "-1/2".into()], Some("QQ".into())));
        assert_eq!(show(Family::SphericalBesselFn, 2, &[]).0, ["3", "0", "-1", "0"]);
        assert_eq!(show(Family::SphericalBesselFn, 0, &[]).0, ["1", "0"]);
        assert_eq!(show(Family::SphericalBesselFnMinus, 1, &[]).0, ["0"]);
        assert_eq!(show(Family::SphericalBesselFnMinus, 3, &[]).0, ["3", "0", "0"]);
        assert_eq!(show(Family::SphericalBesselFnMinus, 2, &[]).0, ["-1", "0"]);
    }

    #[test]
    fn parametric_families_choose_sympy_domains() {
        assert_eq!(show(Family::Laguerre, 2, &["0"]).0, ["1/2", "-2", "1"]);
        let (c, d) = show(Family::Laguerre, 2, &["a"]);
        assert_eq!(d.as_deref(), Some("ZZ(a)"));
        assert_eq!(c[0], "1/2");
        assert_eq!(show(Family::Jacobi, 1, &["1", "2"]).0, ["5/2", "-1/2"]);
        let (c, d) = show(Family::Jacobi, 3, &["0.3", "0.7"]);
        assert_eq!(d.as_deref(), Some("RR"));
        assert_eq!(c.len(), 4);
        assert_eq!(show(Family::Gegenbauer, 2, &["3/2"]).0, ["15/2", "0", "-3/2"]);
    }

    #[test]
    fn exotic_parameters_use_canceled_expressions() {
        let (c, d) = show(Family::Laguerre, 1, &["pi"]);
        assert!(d.is_none());
        assert_eq!(c, ["-1", "1 + pi"]);
    }
}
