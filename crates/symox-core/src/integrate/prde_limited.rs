//! Base-case `limited_integrate_reduce` of SymPy `prde.py`.
//!
//! Same `Poly` domain rules as [`crate::integrate::prde`]: products and
//! differences unify ZZ/QQ, `quo` retracts, `gcd`/`lcm` normalize and
//! `cancel(include=True)` keeps the unified domain. `None`
//! means the reference would raise or the input is outside the port.
use crate::integrate::prde::{splitfactor, P};

/// `rde.order_at_oo(a, d, t)`; `None` is `oo` (zero numerator).
fn order_at_oo(a: &P, d: &P) -> Option<i64> {
    if a.is_zero() {
        None
    } else {
        Some(d.deg() - a.deg())
    }
}

pub struct LimitedReduce {
    pub a: P,
    pub b: P,
    pub n: i64,
    pub g: (P, P),
    pub v: Vec<(P, P)>,
}

/// `limited_integrate_reduce(fa, fd, G, DE)` with `D = [Poly(1, x)]`; the
/// reference returns `(a, b, a, N, g, V)`.
pub fn limited_integrate_reduce(fa: &P, fd: &P, g: &[(P, P)]) -> Option<LimitedReduce> {
    // An empty G fails to unpack `zip(*E)`; zero denominators make the
    // reference divide by zero or take `-oo` degrees.
    if g.is_empty() || fd.is_zero() || g.iter().any(|(_, d)| d.is_zero()) {
        return None;
    }
    let (dn, ds) = splitfactor(fd)?;
    let e = g.iter().map(|(_, gd)| splitfactor(gd)).collect::<Option<Vec<_>>>()?;
    let mut c = dn;
    for (en, _) in &e {
        c = c.lcm(en)?;
    }
    let hn = c.gcd(&c.diff());
    let mut b = hn.derivation().neg();
    let mut hs = ds;
    for (_, es) in &e {
        hs = hs.lcm(es)?;
    }
    if hs.is_zero() || hn.is_zero() {
        return None;
    }
    let a = hn.mul(&hs);
    b = b.sub(&hn.mul(&hs.derivation()).quo(&hs)?);
    let mut mu = order_at_oo(fa, fd);
    for (ga, gd) in g {
        // `min` of `oo` and an integer is the integer.
        if let Some(o) = order_at_oo(ga, gd) {
            mu = Some(mu.map_or(o, |m| m.min(o)));
        }
    }
    // `DE.d.degree(DE.t) == 0`; `max(0, 1 - oo) == 0`.
    let extra = mu.map_or(0, |m| (1 - m).max(0));
    let n = hn.deg() + hs.deg() + extra;
    let ahn = a.mul(&hn);
    let v = g.iter().map(|(ga, gd)| ahn.neg().mul(ga).cancel(gd)).collect::<Option<Vec<_>>>()?;
    let gg = ahn.mul(fa).cancel(fd)?;
    Some(LimitedReduce { a, b, n, g: gg, v })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::Rat;
    fn z(c: &[i64]) -> P {
        P::new(c.iter().map(|&x| Rat::from_i64(x)).collect(), false)
    }
    #[test]
    fn reduce_matches_reference_base_case() {
        // SymPy 1.14: limited_integrate_reduce(Poly(x, x), Poly(x**2, x),
        // [(Poly(1, x), Poly(x, x))], DE(x)) ==
        // (Poly(x, x), Poly(-1, x), Poly(x, x), 1, (Poly(x, x), Poly(1, x)),
        //  [(Poly(-x, x), Poly(1, x))]).
        let r = limited_integrate_reduce(&z(&[0, 1]), &z(&[0, 0, 1]), &[(z(&[1]), z(&[0, 1]))]).unwrap();
        assert_eq!(r.a, z(&[0, 1]));
        assert_eq!(r.b, z(&[-1]));
        assert_eq!(r.n, 1);
        assert_eq!(r.g, (z(&[0, 1]), z(&[1])));
        assert_eq!(r.v, vec![(z(&[0, -1]), z(&[1]))]);
    }
}
