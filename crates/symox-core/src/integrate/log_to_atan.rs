//! Exact QQ version of the logarithm-to-arctangent Euclidean recurrence.
//!
//! Return each atan argument in recurrence order. In particular, do not
//! combine atan terms: equal derivatives alone do not fix branch constants.
use crate::poly::{dense, int::Field, Rat};

fn gcdex<C: Field,E>(a: &[C], b: &[C], divisor: &impl Fn(&C)->bool,
    check: &mut impl FnMut() -> Result<(), E>) -> Result<Option<(Vec<C>, Vec<C>, Vec<C>)>, E>
{
    // The reference half-gcdex, followed by t=(h-s*a)/b. This fixes the
    // Bezout representative and its monic normalization, not just its identity.
    let (mut f, mut g) = (a.to_vec(), b.to_vec());
    let (mut s, mut t) = (vec![C::one()], vec![]);
    while !g.is_empty() {
        check()?;
        if !divisor(&dense::lc(&g)) { return Ok(None); }
        let (q, r) = dense::divrem(&f, &g);
        f = std::mem::replace(&mut g, r);
        let next = dense::sub(&s, &dense::mul(&q, &t));
        s = std::mem::replace(&mut t, next);
    }
    if f.is_empty() || b.is_empty() || !divisor(&dense::lc(&f)) || !divisor(&dense::lc(b)) { return Ok(None); }
    check()?;
    let inverse = dense::lc(&f).inv();
    let s = dense::scale(&s, &inverse);
    let h = dense::scale(&f, &inverse);
    let (t, r) = dense::divrem(&dense::sub(&h, &dense::mul(&s, a)), b);
    check()?;
    Ok(r.is_empty().then_some((s, t, h)))
}

/// Ascending, trimmed coefficient vectors. A refusal retains the general
/// implementation; interruption is always an error. The degree bound limits
/// this dispatch, never the caller's mathematical capabilities.
pub fn arguments<E>(f: &[Rat], g: &[Rat], check: &mut impl FnMut() -> Result<(), E>)
    -> Result<Option<Vec<Vec<Rat>>>, E>
{
    arguments_field(f,g,check)
}

/// The identical Euclidean representative over an exact coefficient field.
/// Callers must keep that field's coefficient scope alive for this operation.
pub(crate) fn arguments_field<C: Field,E>(f: &[C],g: &[C],check: &mut impl FnMut()->Result<(),E>)
    -> Result<Option<Vec<Vec<C>>>,E>
{
    arguments_divisors(f,g,&|_:&C|true,check)
}

/// `divisor` admits each leading coefficient used as a divisor. A caller
/// whose general representation does not normalise some quotients (the
/// reference `EX` cancel form of an irrational inverse) declines there.
pub(crate) fn arguments_divisors<C: Field,E>(f: &[C],g: &[C],divisor: &impl Fn(&C)->bool,
    check: &mut impl FnMut()->Result<(),E>) -> Result<Option<Vec<Vec<C>>>,E>
{
    check()?;
    if f.is_empty() || g.is_empty() || f.len().max(g.len()) > 129
        || f.last().is_some_and(C::is_zero) || g.last().is_some_and(C::is_zero) {
        return Ok(None);
    }
    let (mut f, mut g) = (f.to_vec(), g.to_vec());
    let mut result = Vec::new();
    loop {
        check()?;
        if f.len() < g.len() { (f, g) = (dense::neg(&g), f); }
        if g.is_empty() || !divisor(&dense::lc(&g)) { return Ok(None); }
        let (p, q) = dense::divrem(&f, &g);
        if q.is_empty() {
            result.push(p);
            check()?;
            return Ok(Some(result));
        }
        let Some((s, t, h)) = gcdex(&g, &dense::neg(&f), divisor, check)?
            else { return Ok(None) };
        check()?;
        let numerator = dense::add(&dense::mul(&f, &s), &dense::mul(&g, &t));
        let (u, remainder) = dense::divrem(&numerator, &h);
        // Both divisibility and strict degree descent follow from Euclid.
        // Keep a refusal distinct from a partial atan sum if an invariant fails.
        if !remainder.is_empty() || t.is_empty() || s.len().max(t.len()) >= f.len() {
            return Ok(None);
        }
        result.push(u);
        (f, g) = (s, t);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn q(values: &[i64]) -> Vec<Rat> { values.iter().map(|&n| Rat::from_i64(n)).collect() }
    #[test]
    fn euclid_representatives_preserve_the_atan_branch_constants() {
        let actual = arguments(&q(&[0,1]), &q(&[1,1]), &mut || Ok::<_, ()>(())).unwrap().unwrap();
        assert_eq!(actual, vec![q(&[1,2]), q(&[1])]);
        let actual = arguments(&q(&[1,0,1]), &q(&[-1,1]), &mut || Ok::<_, ()>(())).unwrap().unwrap();
        assert_eq!(actual, vec![vec![Rat::zero(), Rat::from_i64(-1), Rat::new((-1).into(),2.into()),Rat::new((-1).into(),2.into())], q(&[1,1])]);
        assert_eq!(arguments(&q(&[0,0,1]), &q(&[0,1]), &mut || Ok::<_, ()>(())).unwrap().unwrap(), vec![q(&[0,1])]);
    }
    #[test]
    fn interruption_never_returns_a_partial_sum_and_zero_inputs_decline() {
        assert!(arguments(&[], &q(&[1]), &mut || Ok::<_, ()>(())).unwrap().is_none());
        let mut calls=0;
        let result=arguments(&q(&[1,0,0,1]), &q(&[1,1,1]), &mut || { calls+=1;if calls==4{Err("stop")}else{Ok(())} });
        assert!(matches!(result,Err("stop")));
    }
}
