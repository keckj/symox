//! Exact, conservative zero certificates for closed radical expressions.
//!
//! Each radical is replaced by a generator satisfying its defining power
//! equation. Ideal membership proves identities for every compatible choice
//! of roots, hence also for the actual principal roots. A nonzero remainder
//! is deliberately inconclusive: the ideal can have several components.
//!
//! Limits: 10 radicals, power denominators <= 12, numerator exponents <= 64,
//! product of power degrees <= 4096, <= 8192 estimated uncombined terms per
//! expansion and actual input terms in total, and total degree <= 512.
//! Unsupported expressions return None without approximation.
use num_traits::ToPrimitive;
use rustc_hash::FxHashMap;

use super::convert::from_exprs_gens;
use super::groebner::{buchberger, reduce, CoefficientField, MonomialOrder};
use crate::expr::{Constant, Expr, Kind};
use crate::ops::{add, mul, pow, sub};
use crate::simplify::numer_denom;

#[derive(Default)]
struct Radicals { symbols: Vec<Expr>, relations: Vec<Expr>, denominators: Vec<Expr>, seen: FxHashMap<Expr,Expr>, dimension: usize }

/// Bound distribution before performing it. This intentionally overestimates
/// when terms cancel or combine: rejection is preferable to unbounded work in
/// an implicit zero test.
fn expansion_size(e: &Expr) -> Option<usize> {
    fn product(a: usize,b: usize) -> Option<usize> { a.checked_mul(b).filter(|&n| n<=8192) }
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) => Some(1),
        Kind::Add(..) => e.args().iter().try_fold(0usize,|n,a| n.checked_add(expansion_size(a)?).filter(|&n| n<=8192)),
        Kind::Mul(..) => e.args().iter().try_fold(1usize,|n,a| product(n,expansion_size(a)?)),
        Kind::Pow(base,exponent) => {
            let mut power=exponent.as_i64()?.unsigned_abs();
            if power>512 { return None; }
            let mut factor=expansion_size(base)?;
            let mut result=1;
            while power>0 {
                if power&1 != 0 { result=product(result,factor)?; }
                power >>= 1;
                if power>0 { factor=product(factor,factor)?; }
            }
            Some(result)
        }
        _ => None,
    }
}
fn bounded_expand(e: &Expr) -> Option<Expr> {
    expansion_size(e)?;
    Some(crate::expand::expand(e))
}

impl Radicals {
    fn radical(&mut self, original: &Expr, base: Expr, p: i64, q: u32) -> Option<Expr> {
        if let Some(symbol) = self.seen.get(original) { return Some(symbol.clone()); }
        if self.symbols.len() >= 10 { return None; }
        self.dimension = self.dimension.max(1).checked_mul(q as usize)?;
        if self.dimension > 4096 { return None; }
        let symbol = Expr::symbol(&format!("_symox_radical_{}",self.symbols.len()));
        let (n,d) = numer_denom(&bounded_expand(&base)?);
        if !d.is_one() && !self.denominators.contains(&d) { self.denominators.push(d.clone()); }
        let (left,right) = if p > 0 {
            (mul(&[pow(&symbol,&Expr::int(q as i64)),pow(&d,&Expr::int(p))]),pow(&n,&Expr::int(p)))
        } else {
            (mul(&[pow(&symbol,&Expr::int(q as i64)),pow(&n,&Expr::int(-p))]),pow(&d,&Expr::int(-p)))
        };
        self.symbols.push(symbol.clone());
        self.relations.push(sub(&left,&right));
        self.seen.insert(original.clone(),symbol.clone());
        Some(symbol)
    }

    fn rewrite(&mut self, e: &Expr) -> Option<Expr> {
        match e.kind() {
            Kind::Num(n) if n.is_rational() => Some(e.clone()),
            Kind::Const(Constant::I) => self.radical(e,Expr::minus_one(),1,2),
            Kind::Add(..) => Some(add(&e.args().iter().map(|a| self.rewrite(a)).collect::<Option<Vec<_>>>()?)),
            Kind::Mul(..) => Some(mul(&e.args().iter().map(|a| self.rewrite(a)).collect::<Option<Vec<_>>>()?)),
            Kind::Pow(base, exponent) => {
                let original_base = base.clone();
                let base = self.rewrite(base)?;
                if let Some(n) = exponent.as_i64() {
                    if n.unsigned_abs() > 512 { return None; }
                    if n < 0 {
                        let (num,den) = numer_denom(&bounded_expand(&base)?);
                        for guard in [num,den] {
                            if !guard.is_one() && !self.denominators.contains(&guard) { self.denominators.push(guard); }
                        }
                    }
                    return Some(pow(&base,exponent));
                }
                let (p,q) = exponent.as_number()?.numer_denom()?;
                let (p,q) = (p.to_i64()?,q.to_u32()?);
                if p == 0 || p.unsigned_abs() > 64 || !(2..=12).contains(&q) { return None; }
                // All integral powers of one principal q-th root must share
                // the same generator, including inverse powers. Independent
                // variables for a**(1/3) and a**(-1/3) lose this exact identity.
                let key = pow(&original_base,&Expr::rational(1,q as i64));
                let root = self.radical(&key,base,1,q)?;
                if p < 0 && !self.denominators.contains(&root) { self.denominators.push(root.clone()); }
                Some(pow(&root,&Expr::int(p)))
            }
            _ => None,
        }
    }
}

/// Some(true) is an exact algebraic proof; None means no certificate within
/// this bounded algorithm. This routine never guesses that a value is zero.
pub fn is_zero(e: &Expr) -> Option<bool> {
    if e.is_zero() { return Some(true); }
    if !crate::traverse::free_symbols(e).is_empty() { return None; }
    let mut radicals = Radicals::default();
    let rewritten = radicals.rewrite(e)?;
    if radicals.symbols.is_empty() { return rewritten.is_zero().then_some(true); }
    let (numerator,denominator) = numer_denom(&bounded_expand(&rewritten)?);
    let mut expressions = radicals.relations;
    let count = expressions.len();
    expressions.push(numerator);
    expressions.push(denominator);
    expressions.extend(radicals.denominators);
    for expression in &expressions { expansion_size(expression)?; }
    // Newer radicals before their dependencies generally reduce elimination
    // size, while remaining deterministic for repeated certificates.
    radicals.symbols.reverse();
    let converted = from_exprs_gens(&expressions,&radicals.symbols).ok()?;
    if converted.float || converted.gens.len() != converted.ngens || converted.polys.iter().map(|p| p.len()).sum::<usize>() > 8192
        || converted.polys.iter().any(|p| (0..p.len()).any(|i| p.exp(i).iter().map(|&n| n as u64).sum::<u64>() > 512)) { return None; }
    let field = CoefficientField { modulus: None };
    let order = MonomialOrder::Grevlex;
    let basis = buchberger(&converted.polys[..count],order,field);
    if basis.iter().any(|p| p.is_constant() && !p.is_zero()) { return None; }
    let remainder = reduce(&converted.polys[count],&basis,order,field).1;
    if !remainder.is_zero() { return None; }
    for denominator in &converted.polys[count+1..] {
        if denominator.is_constant() && !denominator.is_zero() { continue; }
        // Include internal radicand denominators as well as the final one:
        // clearing t**q*d-n must not hide an undefined radicand n/d.
        let mut pole_ideal = basis.clone();
        pole_ideal.push(denominator.clone());
        let pole_basis = buchberger(&pole_ideal,order,field);
        if !pole_basis.iter().any(|p| p.is_constant() && !p.is_zero()) { return None; }
    }
    Some(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn expr(s: &str) -> Expr { crate::parser::parse_str(s).unwrap() }
    #[test]
    fn closed_radical_ideal_certificates_are_conservative() {
        let t = expr("(sqrt(2)+sqrt(3))**2-5-2*sqrt(6)");
        // The independent sqrt(6) relation does not encode the principal-root
        // product relation, so an unproved identity must remain unknown.
        assert_eq!(is_zero(&t),None);
        assert_eq!(is_zero(&expr("((1+sqrt(2))**2-3-2*sqrt(2))/(1+sqrt(2))")),Some(true));
        assert_eq!(is_zero(&expr("sqrt(2)-1414213562373095/1000000000000000")),None);
        assert_eq!(is_zero(&expr("sin(1)")),None);
    }
    #[test]
    fn cubic_and_quartic_roots_have_exact_polynomial_certificates() {
        let x=Expr::symbol("x");
        for (source,degree) in [("x**3-x+1",3),("x**4+x+1",4)] {
            let polynomial=expr(source);
            let roots=super::super::roots::roots(&polynomial,&[x.clone()]).unwrap();
            assert_eq!(roots.len(),degree);
            for (root,_) in roots {
                let residual=crate::subs::subs1(&polynomial,&x,&root);
                assert_eq!(is_zero(&residual),Some(true),"{}",crate::printer::to_str(&residual));
            }
            assert_eq!(super::super::systems::solve_poly_system(&[polynomial],&[x.clone()]).unwrap().unwrap().len(),degree);
        }
    }
    #[test]
    fn large_held_radical_power_is_rejected_before_expansion() {
        let terms:Vec<_>=[2,3,5,7,11,13,17,19].iter().map(|&n| pow(&Expr::int(n),&Expr::rational(1,2))).collect();
        let held=Expr::from_unevaluated_kind(Kind::Pow(add(&terms),Expr::int(512)));
        let start=std::time::Instant::now();
        assert_eq!(is_zero(&held),None);
        assert!(start.elapsed()<std::time::Duration::from_secs(1));
    }
}
