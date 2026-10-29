//! The reference multivariate gcd over `QQ[y...]` (`dmp_inner_gcd` of
//! `sympy.polys.euclidtools` on a polynomial ring that is neither `ZZ` nor
//! a field): `dmp_multi_deflate`, `_dmp_rr_trivial_gcd` with
//! `_dmp_simplify_gcd`, `dmp_rr_prs_gcd` and, in the last generator,
//! `dup_rr_prs_gcd`, with `dmp_primitive`/`dmp_content` and the ground gcd
//! of `PolyElement.cofactors` over `QQ` (a monomial operand keeps the
//! rational `QQ.gcd` of the coefficients, two multi-term operands give the
//! monic gcd).
//!
//! `RR[y]` polynomials in several generators take this route through
//! `RR[y].get_exact()`, so their gcd keeps the scale chosen by these steps
//! (`gcd(1.5*x*z + 1.5*y, 3.0*x*z + 3.0*y)` is `1.5*x*z + 1.5*y`), where a
//! monic or primitive gcd would differ by a unit.
//!
//! A polynomial is stored flattened as an `MPoly<Rat>` whose first `n`
//! variables are the polynomial generators and whose remaining variables
//! are the parameters `y...` of the ground ring (in the ring's order, so
//! the lex leading coefficient of a ground element is `PolyElement.LC`).
//! Every division is exact; a failed one is reported, never rounded.
use super::convert::{PResult, PolyError};
use super::int::Rat;
use super::mpoly::MPoly;

type P = MPoly<Rat>;

fn inexact_division() -> PolyError {
    PolyError::ExactQuotientFailed("reference QQ[y] gcd: inexact division".into())
}

fn quo(f: &P, g: &P) -> PResult<P> {
    if g.is_one() {
        return Ok(f.clone());
    }
    f.div_exact(g).ok_or_else(inexact_division)
}

/// `QQ.gcd`: gcd of the numerators over the lcm of the denominators.
fn qq_gcd(a: &Rat, b: &Rat) -> Rat {
    Rat::new(a.n.gcd(&b.n), a.d.lcm(&b.d))
}

fn nonnegative(a: &P) -> P {
    if a.lc().is_negative() { a.neg() } else { a.clone() }
}

fn gcd_u32(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

struct Gcd {
    /// Number of polynomial generators.
    n: usize,
    /// The ground is `RR[y...]` itself rather than its exact `QQ[y...]`
    /// (`dup_inner_gcd` does not convert a polynomial ring over `RR`):
    /// `RR.gcd` is one, so a monomial gcd has the coefficient one.
    rr: bool,
}

impl Gcd {
    fn last(&self, v: usize) -> bool {
        v + 1 == self.n
    }

    /// `PolyElement.gcd` over `QQ` (`cofactors(...)[0]`).
    fn ground_gcd(&self, a: &P, b: &P) -> P {
        if a.is_zero() && b.is_zero() {
            return a.clone();
        }
        if a.is_zero() {
            return nonnegative(b);
        }
        if b.is_zero() {
            return nonnegative(a);
        }
        let monomial = |m: &P, o: &P| -> P {
            let mut exps = m.exp(0).to_vec();
            let mut c = m.coeff(0).clone();
            if self.rr {
                c = Rat::one();
            }
            for (e, d) in o.iter() {
                for (x, y) in exps.iter_mut().zip(e) {
                    *x = (*x).min(*y);
                }
                if !self.rr {
                    c = qq_gcd(&c, d);
                }
            }
            P::monomial(&exps, c)
        };
        if a.len() == 1 {
            return monomial(a, b);
        }
        if b.len() == 1 {
            return monomial(b, a);
        }
        super::gcd::gcd_mq(a, b).0
    }

    /// Dense coefficients in `x_v`, leading first, zeros included.
    fn dense(&self, f: &P, v: usize) -> Vec<P> {
        let d = f.degree(v);
        if d < 0 {
            return Vec::new();
        }
        let mut out = vec![P::zero(f.nvars); d as usize + 1];
        for (k, c) in f.coeffs_in(v) {
            out[d as usize - k as usize] = c;
        }
        out
    }

    /// `dmp_content` (or `dup_content` in the last generator).
    fn content(&self, f: &P, v: usize) -> PResult<P> {
        let cs = self.dense(f, v);
        if self.last(v) {
            let mut cont = P::zero(f.nvars);
            for c in &cs {
                cont = self.ground_gcd(&cont, c);
                if cont.is_one() {
                    break;
                }
            }
            return Ok(cont);
        }
        let Some(first) = cs.first() else { return Ok(P::zero(f.nvars)) };
        let mut cont = first.clone();
        for c in &cs[1..] {
            cont = self.inner_gcd(&cont, c, v + 1)?.0;
            if cont.is_one() {
                break;
            }
        }
        Ok(nonnegative(&cont))
    }

    /// `dmp_primitive` (or `dup_primitive`).
    fn primitive(&self, f: &P, v: usize) -> PResult<(P, P)> {
        if f.is_zero() {
            return Ok((self.content(f, v)?, f.clone()));
        }
        let cont = self.content(f, v)?;
        if cont.is_one() {
            return Ok((cont, f.clone()));
        }
        Ok((cont.clone(), quo(f, &cont)?))
    }

    /// `dmp_prem` in `x_v`: `lc(g)**(df - dg + 1) * f` reduced modulo `g`.
    fn prem(&self, f: &P, g: &P, v: usize) -> PResult<P> {
        let (df, dg) = (f.degree(v), g.degree(v));
        if dg < 0 {
            return Err(PolyError::ExactQuotientFailed("polynomial division by zero".into()));
        }
        if df < dg {
            return Ok(f.clone());
        }
        let mut steps = df - dg + 1;
        let lc_g = g.lc_in(v);
        let mut r = f.clone();
        let mut dr = df;
        loop {
            let lc_r = r.lc_in(v);
            let mut e = vec![0u32; f.nvars];
            e[v] = (dr - dg) as u32;
            steps -= 1;
            r = r.mul(&lc_g).sub(&g.mul(&lc_r).mul_term(&e, &Rat::one()));
            let previous = dr;
            dr = r.degree(v);
            if dr < dg {
                break;
            }
            if dr >= previous {
                return Err(PolyError::ExactQuotientFailed("pseudo-division failed".into()));
            }
        }
        Ok(r.mul(&lc_g.pow(steps as u32)))
    }

    /// The last polynomial of `dmp_subresultants(f, g)`.
    fn last_subresultant(&self, f: &P, g: &P, v: usize) -> PResult<P> {
        let (mut f, mut g) = (f.clone(), g.clone());
        let (mut n, mut m) = (f.degree(v), g.degree(v));
        if n < m {
            std::mem::swap(&mut f, &mut g);
            std::mem::swap(&mut n, &mut m);
        }
        if f.is_zero() || g.is_zero() {
            return Ok(f);
        }
        let nv = f.nvars;
        let minus_one = P::constant(Rat::one().neg(), nv);
        let mut d = n - m;
        let mut h = self.prem(&f, &g, v)?.mul(&minus_one.pow((d + 1) as u32));
        let mut lc = g.lc_in(v);
        let mut c = lc.pow(d as u32).neg();
        let mut last = g.clone();
        while !h.is_zero() {
            let k = h.degree(v);
            last = h.clone();
            f = std::mem::replace(&mut g, h);
            d = m - k;
            m = k;
            let b = lc.neg().mul(&c.pow(d as u32));
            h = quo(&self.prem(&f, &g, v)?, &b)?;
            lc = g.lc_in(v);
            c = if d > 1 { quo(&lc.neg().pow(d as u32), &c.pow((d - 1) as u32))? } else { lc.neg() };
        }
        Ok(last)
    }

    /// `dmp_inner_gcd` in the generators `x_v...` (`dup_inner_gcd`, which
    /// does not deflate, in the last one).
    fn inner_gcd(&self, f: &P, g: &P, v: usize) -> PResult<(P, P, P)> {
        if self.last(v) {
            return self.rr_prs_gcd(f, g, v);
        }
        let nv = f.nvars;
        let mut b = vec![0u32; self.n];
        for p in [f, g] {
            for (e, _) in p.iter() {
                for i in v..self.n {
                    b[i] = gcd_u32(b[i], e[i]);
                }
            }
        }
        for x in b.iter_mut() {
            if *x == 0 {
                *x = 1;
            }
        }
        if b.iter().all(|&x| x == 1) {
            return self.rr_prs_gcd(f, g, v);
        }
        let map = |p: &P, deflate: bool| -> P {
            P::from_terms(nv, p.iter().map(|(e, c)| {
                let mut e = e.to_vec();
                for i in v..self.n {
                    e[i] = if deflate { e[i] / b[i] } else { e[i] * b[i] };
                }
                (e, c.clone())
            }).collect())
        };
        let (h, p, q) = self.rr_prs_gcd(&map(f, true), &map(g, true), v)?;
        Ok((map(&h, false), map(&p, false), map(&q, false)))
    }

    fn rr_prs_gcd(&self, f: &P, g: &P, v: usize) -> PResult<(P, P, P)> {
        let nv = f.nvars;
        let (zero, one) = (P::zero(nv), P::one(nv));
        let minus_one = one.neg();
        // `_dup_rr_trivial_gcd` / `_dmp_rr_trivial_gcd`
        if f.is_zero() && g.is_zero() {
            return Ok((zero.clone(), zero.clone(), zero));
        }
        if f.is_zero() {
            return Ok(if g.lc().is_negative() { (g.neg(), zero, minus_one) } else { (g.clone(), zero, one) });
        }
        if g.is_zero() {
            return Ok(if f.lc().is_negative() { (f.neg(), minus_one, zero) } else { (f.clone(), one, zero) });
        }
        if self.last(v) {
            let (fc, ff) = self.primitive(f, v)?;
            let (gc, gg) = self.primitive(g, v)?;
            let mut c = self.ground_gcd(&fc, &gc);
            let (_, h) = self.primitive(&self.last_subresultant(&ff, &gg, v)?, v)?;
            if h.lc().is_negative() {
                c = c.neg();
            }
            let h = h.mul(&c);
            return Ok((h.clone(), quo(f, &h)?, quo(g, &h)?));
        }
        if f.is_one() || g.is_one() {
            return Ok((one, f.clone(), g.clone()));
        }
        // `_dmp_simplify_gcd` (`USE_SIMPLIFY_GCD` is on by default)
        let (df, dg) = (f.degree(v), g.degree(v));
        if df == 0 || dg == 0 {
            let (a, b) = if df == 0 && dg == 0 {
                (f.lc_in(v), g.lc_in(v))
            } else if df == 0 {
                (f.lc_in(v), self.content(g, v)?)
            } else {
                (self.content(f, v)?, g.lc_in(v))
            };
            let h = self.inner_gcd(&a, &b, v + 1)?.0;
            return Ok((h.clone(), quo(f, &h)?, quo(g, &h)?));
        }
        let (fc, ff) = self.primitive(f, v)?;
        let (gc, gg) = self.primitive(g, v)?;
        let h = self.last_subresultant(&ff, &gg, v)?;
        let c = self.rr_prs_gcd(&fc, &gc, v + 1)?.0;
        let (_, h) = self.primitive(&h, v)?;
        let mut h = h.mul(&c);
        if h.lc().is_negative() {
            h = h.neg();
        }
        Ok((h.clone(), quo(f, &h)?, quo(g, &h)?))
    }
}

/// `dmp_inner_gcd(f, g, n - 1, QQ[y...])` for polynomials in `n >= 2`
/// generators (the first `n` variables) over the parameters that follow.
pub fn inner_gcd(f: &P, g: &P, n: usize) -> PResult<(P, P, P)> {
    debug_assert!(n >= 2 && f.nvars == g.nvars && n <= f.nvars);
    super::mpoly::checked_operation(|| Gcd { n, rr: false }.inner_gcd(f, g, 0))
}

/// `dmp_primitive(f, n - 1, QQ[y...])` in the first generator: the content
/// (free of the first generator) and the primitive part.
///
/// Its content gcds are `dmp_gcd` calls one level down over `RR[y...]`:
/// with two generators that is `dup_inner_gcd`, which works over
/// `RR[y...]` itself, and with more it is the converted `QQ[y...]` gcd.
pub fn primitive_main(f: &P, n: usize) -> PResult<(P, P)> {
    debug_assert!(n >= 2 && n <= f.nvars);
    super::mpoly::checked_operation(|| Gcd { n, rr: n == 2 }.primitive(f, 0))
}

/// `RR[y...].gcd` of two polynomials free of the generators (the first
/// `n` variables): `PolyElement.cofactors(...)[0]` over `RR`.
pub fn rr_ground_gcd(a: &P, b: &P) -> P {
    Gcd { n: 1, rr: true }.ground_gcd(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::int::Int;

    fn p(nv: usize, t: &[(&[u32], i64, i64)]) -> P {
        P::from_terms(nv, t.iter().map(|(e, a, b)| (e.to_vec(), Rat::new(Int::S(*a), Int::S(*b)))).collect())
    }

    #[test]
    fn keeps_the_reference_scale() {
        // x, z over QQ[y]: gcd(3/2*x*z + 3/2*y, 3*x*z + 3*y) = 3/2*x*z + 3/2*y
        let f = p(3, &[(&[1, 1, 0], 3, 2), (&[0, 0, 1], 3, 2)]);
        let g = p(3, &[(&[1, 1, 0], 3, 1), (&[0, 0, 1], 3, 1)]);
        let (h, a, b) = inner_gcd(&f, &g, 2).unwrap();
        assert_eq!(h, f);
        assert_eq!(a, P::one(3));
        assert_eq!(b, P::constant(Rat::from_i64(2), 3));
        // gcd(1/2*x**2*z - 1/2*y**2, 3/2*x*z + 3/2*y) = 1/2
        let f = p(3, &[(&[2, 1, 0], 1, 2), (&[0, 0, 2], -1, 2)]);
        let g = p(3, &[(&[1, 1, 0], 3, 2), (&[0, 0, 1], 3, 2)]);
        let (h, _, _) = inner_gcd(&f, &g, 2).unwrap();
        assert_eq!(h, P::constant(Rat::new(Int::S(1), Int::S(2)), 3));
    }
}
