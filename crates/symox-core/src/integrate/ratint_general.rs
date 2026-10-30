//! Horowitz-Ostrogradsky reduction over every exact coefficient field of
//! the reference `ratint_ratpart`: algebraic fields `QQ<a>` and the Gaussian
//! rationals (`AN`), rational function fields `QQ(a, ...)` (`RF`) and the
//! expression domain `EX`.
//!
//! The caller supplies `u = gcd(g, g')` and `v = g/u` as computed by the
//! reference `g.cofactors(g.diff())`. The undetermined coefficients of
//! `A` (degree < deg u) and `B` (degree < deg v) form a square linear system
//! `f = A' v - A (u' v / u) + B u` that the reference passes to `solve`.
//! Its solution is unique, so an exact Gauss-Jordan elimination gives the
//! same `A`, `B`. The identity is recomputed before anything is returned:
//! a pivot whose zero test was undecided cannot leak a wrong solution.
use crate::poly::{ddm, dense, Field};

/// Exact `A, B` (ascending, trimmed) with `f = A' v - A (u' v / u) + B u`,
/// or `None` when the system is not provably of full rank, the inputs are
/// outside the reference preconditions, or the identity does not hold.
pub fn solve<C: Field, E>(f: &[C], u: &[C], v: &[C], max_size: usize,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Option<(Vec<C>, Vec<C>)>, E>
{
    check()?;
    if u.is_empty() || v.is_empty() || u.last().is_some_and(C::is_zero)
        || v.last().is_some_and(C::is_zero) || f.last().is_some_and(C::is_zero) {
        return Ok(None);
    }
    let n = u.len() - 1;
    let m = v.len() - 1;
    let size = n + m;
    if size == 0 || size > max_size || f.len() > size { return Ok(None); }
    let (correction, remainder) = dense::divrem(&dense::mul(&dense::derivative(u), v), u);
    if !remainder.is_empty() { return Ok(None); }
    check()?;
    let mut system = vec![vec![C::zero(); size + 1]; size];
    for (row, value) in f.iter().enumerate() { system[row][size] = value.clone(); }
    for column in 0..n {
        check()?;
        // A' v - A (u' v / u), for A = x**column.
        if column != 0 {
            let factor = C::from_i64(column as i64);
            for (j, value) in v.iter().enumerate() {
                let cell = &mut system[j + column - 1][column];
                *cell = cell.add(&value.mul(&factor));
            }
        }
        for (j, value) in correction.iter().enumerate() {
            if j + column >= size {
                if !value.is_zero() { return Ok(None); }
                continue;
            }
            let cell = &mut system[j + column][column];
            *cell = cell.sub(value);
        }
    }
    for column in 0..m {
        check()?;
        for (j, value) in u.iter().enumerate() {
            if j + column >= size { return Ok(None); }
            system[j + column][n + column] = value.clone();
        }
    }
    let pivots = ddm::field_irref(&mut system, false, check, &mut |_| Ok(()))?;
    if pivots != Some((0..size).collect()) { return Ok(None); }
    let mut a: Vec<C> = system.iter().take(n).map(|row| row[size].clone()).collect();
    let mut b: Vec<C> = system.iter().skip(n).map(|row| row[size].clone()).collect();
    dense::trim(&mut a);
    dense::trim(&mut b);
    check()?;
    let mut identity = dense::add(
        &dense::sub(&dense::mul(&dense::derivative(&a), v), &dense::mul(&a, &correction)),
        &dense::mul(&b, u),
    );
    dense::trim(&mut identity);
    let mut expected = f.to_vec();
    dense::trim(&mut expected);
    if identity.len() != expected.len()
        || identity.iter().zip(&expected).any(|(p, q)| p != q && !p.sub(q).is_zero()) {
        return Ok(None);
    }
    check()?;
    Ok(Some((a, b)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{algebraic::{AlgebraicScope, AN}, gaussian, Rat};
    #[test]
    fn gaussian_repeated_pole() {
        // 1/(x**2 + I)**2: u = x**2 + I, v = x**2 + I.
        let field = gaussian::field();
        let _scope = AlgebraicScope::new(field.clone());
        let c = |re: i64, im: i64| gaussian::make(Rat::from_i64(re), Rat::from_i64(im));
        let u = vec![c(0, 1), c(0, 0), c(1, 0)];
        let (a, b) = solve(&[c(1, 0)], &u, &u, 64, &mut || Ok::<_, ()>(())).unwrap().unwrap();
        // A = -I*x/2, B = -I/2.
        let half = Rat::from_i64(1).div(&Rat::from_i64(2));
        assert_eq!(a, vec![c(0, 0), gaussian::make(Rat::zero(), half.neg())]);
        assert_eq!(b, vec![gaussian::make(Rat::zero(), half.neg())]);
        let _ = AN::rational(Rat::one(), field);
    }
    #[test]
    fn declines_outside_preconditions() {
        let q = |v: &[i64]| v.iter().map(|&k| Rat::from_i64(k)).collect::<Vec<_>>();
        assert!(solve(&q(&[1]), &q(&[1]), &q(&[1]), 64, &mut || Ok::<_, ()>(())).unwrap().is_none());
        assert!(solve(&q(&[1, 2, 3]), &q(&[0, 1]), &q(&[0, 1]), 64, &mut || Ok::<_, ()>(())).unwrap().is_none());
        assert!(matches!(solve(&q(&[1]), &q(&[0, 1]), &q(&[0, 1]), 64, &mut || Err("stop")), Err("stop")));
    }
}
