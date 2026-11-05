//! The subresultant step of the reference `ratint_logpart` when the
//! coefficients of `f/g` lie in an algebraic field `QQ<a>`.
//!
//! The reference forms `b = f - g'*Poly(t, x)`; unification of `QQ<a>` with
//! `ZZ[t]` gives the composite ring `QQ<a>[t]`, and
//! `resultant(g, b, includePRS=True)` runs `dup_inner_subresultants` over that
//! ring. The subresultants are unique, and the Brown-Traub recurrence of
//! [`crate::poly::prs::subresultants`] with `x` as main variable of
//! `QQ<a>[x, t]` performs the same exact divisions, so the sequence and the
//! resultant are the reference ones. Expressions are produced as the
//! reference `PolynomialRing.to_sympy` does: `Mul(c, t**k)` summed over the
//! terms, with `c` the `AlgebraicField.to_sympy` image of the coefficient.
//!
//! Division of a sequence member by a polynomial of `t` free of `x` follows
//! `Poly.quo(auto=True)` over `QQ<a>[t]`: the quotient is computed in
//! `QQ<a>(t)` and retracted to the ring. Only the retractable case is
//! supported; otherwise `None` is returned and the caller keeps the reference
//! path.
use std::sync::Arc;

use crate::expr::Expr;
use crate::ops::{add, mul, pow};
use crate::poly::algebraic::{AlgebraicPoly, AlgebraicScope, NumberField, AN};
use crate::poly::convert::{PResult, PolyError};
use crate::poly::int::{Coeff, Rat};
use crate::poly::mpoly::MPoly;

/// A polynomial of `x` over `QQ<a>[t]`, stored in `QQ<a>[x, t]`.
#[derive(Clone, Debug)]
pub struct ExtensionRingPoly {
    rep: MPoly<AN>,
    field: Arc<NumberField>,
    t: Expr,
}

/// `PolynomialRing.to_sympy` of an element of `QQ<a>[t]` given in the two
/// variables `(x, t)` with `x`-degree zero.
fn ring_element_expr(c: &MPoly<AN>, t: &Expr) -> Expr {
    add(&c
        .iter()
        .map(|(e, a)| {
            let mut factors = vec![a.as_expr()];
            if e[1] != 0 {
                factors.push(pow(t, &Expr::int(i64::from(e[1]))));
            }
            mul(&factors)
        })
        .collect::<Vec<_>>())
}

impl ExtensionRingPoly {
    /// Degree in `x`, `None` for the zero polynomial (`-oo`).
    pub fn degree(&self) -> Option<u32> {
        let d = self.rep.degree(0);
        u32::try_from(d).ok()
    }

    pub fn is_zero(&self) -> bool {
        self.rep.is_zero()
    }

    /// `Poly.LC()`: the leading coefficient in `x` as an expression.
    pub fn lc_expr(&self) -> Expr {
        let _scope = AlgebraicScope::new(self.field.clone());
        ring_element_expr(&self.rep.lc_in(0), &self.t)
    }

    /// `Poly.terms()`: `(x-degree, coefficient)` by decreasing degree.
    pub fn terms(&self) -> Vec<(u32, Expr)> {
        let _scope = AlgebraicScope::new(self.field.clone());
        self.rep
            .coeffs_in(0)
            .into_iter()
            .map(|(k, c)| (k, ring_element_expr(&c, &self.t)))
            .collect()
    }

    /// `self.quo(Poly(d, x))` for `d` a polynomial of `t` only, with
    /// coefficients in `QQ<a>` (given as an expression). `None` when `d`
    /// involves `x`, is zero, cannot be converted, or when a quotient
    /// coefficient does not retract from `QQ<a>(t)` to `QQ<a>[t]`.
    pub fn quo_by_t_polynomial(&self, d: &Expr, x: &Expr) -> Option<Self> {
        let _scope = AlgebraicScope::new(self.field.clone());
        let gens = [x.clone(), self.t.clone()];
        let divisor = AlgebraicPoly::from_expr(d, &gens, self.field.clone()).ok()?;
        if divisor.rep.is_zero() || divisor.rep.degree(0) != 0
            || !divisor.number_field.same_embedding(&self.field) {
            return None;
        }
        let mut terms: Vec<(Vec<u32>, AN)> = Vec::new();
        for (k, c) in self.rep.coeffs_in(0) {
            let q = c.div_exact(&divisor.rep)?;
            for (e, a) in q.iter() {
                terms.push((vec![k, e[1]], a.clone()));
            }
        }
        Some(Self { rep: MPoly::from_terms(2, terms), field: self.field.clone(), t: self.t.clone() })
    }
}

/// Lifts a univariate polynomial of `x` over `QQ<a>` to `QQ<a>[x, t]`.
fn lift(p: &MPoly<AN>) -> MPoly<AN> {
    MPoly::from_terms(2, p.iter().map(|(e, c)| (vec![e[0], 0], c.clone())).collect())
}

/// `res, R = resultant(g, f - g.diff()*Poly(t, x), includePRS=True)` over
/// `QQ<a>[t]`: the resultant as an expression of `t` and the subresultant
/// sequence. `f` and `g` are univariate in `x` over the same field.
pub fn logpart_prs(f: &AlgebraicPoly, g: &AlgebraicPoly, t: &Expr)
    -> PResult<(Expr, Vec<ExtensionRingPoly>)>
{
    if f.gens.len() != 1 || g.gens.len() != 1 || f.gens != g.gens {
        return Err(PolyError::Unsupported("univariate polynomials of one generator expected".into()));
    }
    if !f.number_field.same_embedding(&g.number_field) {
        return Err(PolyError::Unsupported("coefficients in different algebraic fields".into()));
    }
    let field = g.number_field.clone();
    let _scope = AlgebraicScope::new(field.clone());
    let a = lift(&g.rep);
    // g' * t
    let derivative = MPoly::from_terms(2, g.rep.iter().filter(|(e, _)| e[0] > 0)
        .map(|(e, c)| (vec![e[0] - 1, 1], c.mul(&AN::rational(Rat::from_i64(i64::from(e[0])), field.clone()))))
        .collect());
    let b = lift(&f.rep).sub(&derivative);
    let prs = crate::poly::prs::subresultants(&a, &b, 0)?;
    let resultant = ring_element_expr(&prs.resultant, t);
    let sequence = prs.sequence.into_iter()
        .map(|rep| ExtensionRingPoly { rep, field: field.clone(), t: t.clone() })
        .collect();
    Ok((resultant, sequence))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    fn parse(s: &str) -> Option<Expr> { parse_str(s).ok() }

    #[test]
    fn quadratic_extension_logpart_prs() {
        let x = parse("x").unwrap();
        let t = parse("t").unwrap();
        let s2 = parse("sqrt(2)").unwrap();
        let field = NumberField::from_generator(&s2).unwrap();
        let g = AlgebraicPoly::from_expr(&parse("x**2 + sqrt(2)").unwrap(), &[x.clone()], field.clone()).unwrap();
        let f = AlgebraicPoly::from_expr(&parse("1").unwrap(), &[x.clone()], field).unwrap();
        let (res, prs) = logpart_prs(&f, &g, &t).unwrap();
        assert_eq!(res, parse("4*sqrt(2)*t**2 + 1").unwrap());
        assert_eq!(prs.iter().map(|p| p.degree()).collect::<Vec<_>>(), vec![Some(2), Some(1), Some(0)]);
        assert_eq!(prs[1].lc_expr(), parse("-2*t").unwrap());
    }
}
