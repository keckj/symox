//! Exact Horowitz–Ostrogradsky reduction over QQ.
//!
//! A refusal keeps the complete general-domain implementation. The bounded
//! linear system is solved exactly and the polynomial identity is checked
//! before publishing its rational and square-free parts.
use crate::poly::{ddm, dense, Rat};

pub struct Reduction {
    pub rational_numerator: Vec<Rat>,
    pub rational_denominator: Vec<Rat>,
    pub logarithmic_numerator: Vec<Rat>,
    pub logarithmic_denominator: Vec<Rat>,
}

fn quotient(a: &[Rat], b: &[Rat]) -> Option<Vec<Rat>> {
    let (q, r) = dense::divrem(a, b);
    r.is_empty().then_some(q)
}

fn cancel(a: &[Rat], b: &[Rat]) -> (Vec<Rat>, Vec<Rat>) {
    if a.is_empty() { return (vec![], vec![Rat::one()]); }
    let gcd = dense::gcd_field(a, b);
    let numerator = quotient(a, &gcd).expect("polynomial gcd divides numerator");
    let denominator = quotient(b, &gcd).expect("polynomial gcd divides denominator");
    let scale = dense::lc(&denominator).inv();
    (dense::scale(&numerator, &scale), dense::scale(&denominator, &scale))
}

/// Ascending, trimmed QQ coefficients. Only proper fractions with denominator
/// degree at most 128 select this kernel; all other inputs keep the caller's
/// existing algorithm. Interruption is an error, never a successful refusal.
pub fn reduce<E>(f: &[Rat], g: &[Rat], check: &mut impl FnMut() -> Result<(), E>)
    -> Result<Option<Reduction>, E>
{
    check()?;
    if g.len() < 2 || g.len() > 129 || f.len() >= g.len()
        || g.last().is_some_and(Rat::is_zero) || f.last().is_some_and(Rat::is_zero) {
        return Ok(None);
    }
    let u = dense::gcd_field(g, &dense::derivative(g));
    check()?;
    let Some(v) = quotient(g, &u) else { return Ok(None) };
    let Some(correction) = quotient(&dense::mul(&dense::derivative(&u), &v), &u)
        else { return Ok(None) };
    let n = u.len() - 1;
    let m = v.len() - 1;
    let size = n + m;
    let mut system = vec![vec![Rat::zero(); size + 1]; size];
    for (row, value) in f.iter().enumerate() { system[row][size] = value.clone(); }
    for column in 0..n {
        check()?;
        // A' v - A (u' v / u), for A = x**column.
        if column != 0 {
            let factor = Rat::from_i64(column as i64);
            for (j, value) in v.iter().enumerate() {
                system[j + column - 1][column] = value.mul(&factor);
            }
        }
        for (j, value) in correction.iter().enumerate() {
            let cell = &mut system[j + column][column];
            *cell = cell.sub(value);
        }
    }
    for column in 0..m {
        check()?;
        for (j, value) in u.iter().enumerate() {
            system[j + column][n + column] = value.clone();
        }
    }
    let pivots = ddm::irref(&mut system, false, check, &mut |_| Ok(()))?;
    if pivots != Some((0..size).collect()) { return Ok(None); }
    let mut a: Vec<_> = system.iter().take(n).map(|row| row[size].clone()).collect();
    let mut b: Vec<_> = system.iter().skip(n).map(|row| row[size].clone()).collect();
    dense::trim(&mut a); dense::trim(&mut b);
    check()?;
    let identity = dense::add(
        &dense::sub(&dense::mul(&dense::derivative(&a), &v), &dense::mul(&a, &correction)),
        &dense::mul(&b, &u),
    );
    if identity != f { return Ok(None); }
    let (rational_numerator, rational_denominator) = cancel(&a, &u);
    let (logarithmic_numerator, logarithmic_denominator) = cancel(&b, &v);
    check()?;
    Ok(Some(Reduction { rational_numerator, rational_denominator,
        logarithmic_numerator, logarithmic_denominator }))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn q(values: &[i64]) -> Vec<Rat> { values.iter().map(|&v| Rat::from_i64(v)).collect() }
    #[test]
    fn repeated_pole_and_square_free_parts_are_exact() {
        let r = reduce(&q(&[1]), &q(&[0, 0, 1]), &mut || Ok::<_, ()>(())).unwrap().unwrap();
        assert_eq!(r.rational_numerator, q(&[-1]));
        assert_eq!(r.rational_denominator, q(&[0, 1]));
        assert!(r.logarithmic_numerator.is_empty());
        let r = reduce(&q(&[36]), &q(&[-2, 1, 4, -2, -2, 1]), &mut || Ok::<_, ()>(())).unwrap().unwrap();
        assert_eq!(r.rational_numerator, q(&[6, 12]));
        assert_eq!(r.rational_denominator, q(&[-1, 0, 1]));
        assert_eq!(r.logarithmic_numerator, q(&[12]));
        assert_eq!(r.logarithmic_denominator, q(&[-2, -1, 1]));
    }
    #[test]
    fn refuses_outside_selection_and_propagates_interruption() {
        assert!(reduce(&q(&[1, 1]), &q(&[1, 1]), &mut || Ok::<_, ()>(())).unwrap().is_none());
        assert!(reduce(&q(&[1]), &q(&[0]), &mut || Ok::<_, ()>(())).unwrap().is_none());
        let error = reduce(&q(&[1]), &q(&[0, 0, 1]), &mut || Err("interrupt"));
        assert!(matches!(error, Err("interrupt")));
    }
}
