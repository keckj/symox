//! Correctly rounded products/quotients of finite literal CC domain elements.
//!
//! Each component is a ratio of sums of exact dyadic terms/products. Keep
//! these sums sparse: exponent gaps must never become giant integer shifts.
//! A fast approximate quotient is accepted only after exact comparisons with
//! both rounding-cell midpoints. Otherwise binary search the p-bit mantissa.
use crate::bigfloat::{BigFloat, Exponent};
use crate::evalf_mp::MpComplex;
use crate::{Expr, Number};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::One;
use std::cmp::Ordering;

fn product(a: &BigFloat, b: &BigFloat) -> BigFloat {
    BigFloat::from_parts(&a.man * &b.man, &a.exp + &b.exp)
}

/// An exact sum's leading survivor. A discarded tail is strictly smaller
/// than half its absolute value, so its sign is certified and its magnitude
/// differs from the exact sum by at most one binary exponent.
fn leading(terms: &[BigFloat]) -> Option<BigFloat> {
    let mut terms: Vec<_> = terms.iter().filter(|t| !t.is_zero()).cloned().collect();
    loop {
        if terms.is_empty() { return Some(BigFloat::zero()); }
        if terms.len() == 1 { return terms.pop(); }
        terms.sort_unstable_by(|a,b| b.mag_exact().cmp(&a.mag_exact()));
        let rest = terms.len()-1;
        let log = usize::BITS - (rest-1).leading_zeros(); // ceil(log2(rest))
        if terms[0].mag_exact() > terms[1].mag_exact() + i128::from(log) + 1 {
            return Some(terms.remove(0));
        }
        let b = terms.remove(1); let a = terms.remove(0);
        // Their leading positions are close; alignment is bounded by their
        // already stored mantissa widths and log(number of terms), even for
        // arbitrary-precision exponents. Never align across a large gap.
        let exponent = a.exp.min_ref(&b.exp).clone();
        let am = &a.man << (&a.exp-&exponent).to_usize()?;
        let bm = &b.man << (&b.exp-&exponent).to_usize()?;
        let sum = BigFloat::from_parts(am+bm,exponent);
        if !sum.is_zero() { terms.push(sum); }
    }
}

fn sign(terms: &[BigFloat]) -> Option<Ordering> {
    Some(leading(terms)?.sign().cmp(&0))
}

fn compare(numerator: &[BigFloat], denominator: &[BigFloat], boundary: &BigFloat) -> Option<Ordering> {
    let mut terms = numerator.to_vec();
    terms.extend(denominator.iter().map(|d| product(d,boundary).neg()));
    sign(&terms)
}

fn power(exponent: Exponent) -> BigFloat { BigFloat::from_parts(BigInt::one(),exponent) }

fn magnitude(terms: &[BigFloat]) -> Option<Exponent> {
    let mut exponent: Exponent = leading(terms)?.abs().mag_exact()-1i128;
    while compare(terms,&[BigFloat::one()],&power(exponent.clone()))? == Ordering::Less { exponent-=1; }
    while compare(terms,&[BigFloat::one()],&power(&exponent+1))? != Ordering::Less { exponent+=1; }
    Some(exponent)
}

fn in_cell(n: &[BigFloat], d: &[BigFloat], candidate: &BigFloat, p: u32) -> Option<bool> {
    if candidate.is_zero() { return Some(false); }
    let candidate = candidate.abs().round_to(p);
    let spacing = power(candidate.mag_exact()-i128::from(p));
    let power_of_two = candidate.man.trailing_zeros().is_some_and(|z| z+1==candidate.man.bits());
    let below = if power_of_two { spacing.ldexp(-1) } else { spacing.clone() };
    let guard = p.checked_add(2)?;
    let lower = candidate.sub(&below.ldexp(-1),guard);
    let upper = candidate.add(&spacing.ldexp(-1),guard);
    // At a binade boundary the lower cell has half the spacing; this also
    // handles p=1, whose lower tie belongs here but upper tie belongs above.
    let even = candidate.man.bits() < u64::from(p) || candidate.man.is_even();
    let lower_even = power_of_two || even;
    let low = compare(n,d,&lower)?; let high = compare(n,d,&upper)?;
    Some((low==Ordering::Greater || low==Ordering::Equal && lower_even)
         && (high==Ordering::Less || high==Ordering::Equal && even))
}

fn round_ratio(n: &[BigFloat], d: &[BigFloat], candidate: &BigFloat, p: u32) -> Option<BigFloat> {
    let direction = sign(n)?;
    if direction==Ordering::Equal { return Some(BigFloat::zero()); }
    let negative = direction==Ordering::Less;
    let n: Vec<_> = n.iter().map(|v| if negative { v.neg() } else { v.clone() }).collect();
    if in_cell(&n,d,candidate,p)? {
        return Some(if negative {candidate.abs().neg()} else {candidate.abs()}.round_to(p));
    }
    // No fixed guard budget or agreement between two approximations serves
    // as a certificate. Search a bounded p-bit significand with exact signs.
    let mut exponent = magnitude(&n)?-magnitude(d)?;
    if compare(&n,d,&power(exponent.clone()))?==Ordering::Less { exponent-=1; }
    let scale = exponent-i128::from(p-1);
    let mut lower = BigInt::one() << (p-1) as usize;
    let mut upper = &lower << 1usize;
    while &upper-&lower > BigInt::one() {
        let middle = (&lower+&upper) >> 1usize;
        let boundary = BigFloat::from_parts(middle.clone(),scale.clone());
        if compare(&n,d,&boundary)?==Ordering::Less { upper=middle; } else { lower=middle; }
    }
    let midpoint = BigFloat::from_parts((&lower<<1usize)+1,&scale-1);
    let cmp = compare(&n,d,&midpoint)?;
    if cmp==Ordering::Greater || cmp==Ordering::Equal && lower.is_odd() { lower+=1; }
    let value = BigFloat::from_parts(if negative {-lower} else {lower},scale).round_to(p);
    Some(value)
}

fn literal(expr: &Expr, p: u32) -> Option<MpComplex> {
    let (re,im) = crate::complex_parts::literal_parts(expr)?;
    let part = |e: &Expr| match e.as_number()? {
        Number::MpFloat(v,_) => Some(v.clone()),
        n if !n.is_nan() && n.to_f64().is_finite() => Some(n.to_bigfloat(p)),
        _ => None,
    };
    Some(MpComplex{re:part(&re)?,im:part(&im)?})
}

pub(super) fn binary(a: &Expr, b: &Expr, op: char, p: u32) -> Option<Expr> {
    if p==0 || crate::evalf_mp::host_interrupted() { return None; }
    let a=literal(a,p)?; let b=literal(b,p)?;
    let (re,im,d,approximate)=match op {
        '/' => {
            if b.re.is_zero() && b.im.is_zero() { return None; } // pole error belongs to the caller
            p.checked_add(10)?;
            ([product(&a.re,&b.re),product(&a.im,&b.im)],
             [product(&a.im,&b.re),product(&a.re,&b.im).neg()],
             vec![product(&b.re,&b.re),product(&b.im,&b.im)],a.div(&b,p)?)
        },
        '*' => ([product(&a.re,&b.re),product(&a.im,&b.im).neg()],
                [product(&a.re,&b.im),product(&a.im,&b.re)],vec![BigFloat::one()],a.mul(&b,p)),
        _ => return None,
    };
    let value=MpComplex{re:round_ratio(&re,&d,&approximate.re,p)?,im:round_ratio(&im,&d,&approximate.im,p)?};
    Some(crate::evalf_mp::to_expr(value,p))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_sign_and_magnitude_do_not_materialize_exponent_gaps() {
        let huge = Exponent::from(BigInt::one()<<190usize);
        let large=power(huge.clone()); let one=BigFloat::one();
        assert_eq!(sign(&[large.clone(),one.neg()]),Some(Ordering::Greater));
        assert_eq!(magnitude(&[large.clone(),one.neg()]),Some(&huge-1));
        assert_eq!(sign(&[large.clone(),one.clone(),large.neg()]),Some(Ordering::Greater));
        let tiny=power(-huge);
        assert_eq!(sign(&[one.clone(),one.neg(),tiny.clone()]),Some(Ordering::Greater));
        assert_eq!(sign(&[tiny.clone(),tiny.neg()]),Some(Ordering::Equal));
    }
    #[test]
    fn midpoint_parity_and_binade_transition_at_one_bit() {
        let one=BigFloat::one(); let two=BigFloat::from_i64(2);
        let three_quarters=BigFloat::from_parts(BigInt::from(3),-2);
        let three_halves=three_quarters.ldexp(1);
        assert_eq!(round_ratio(&[three_quarters],&[one.clone()],&one,1),Some(one.clone()));
        assert_eq!(round_ratio(&[three_halves],&[one.clone()],&one,1),Some(two.clone()));
        assert_eq!(round_ratio(&[BigFloat::from_i64(3)],&[one.clone()],&two,1),Some(BigFloat::from_i64(4)));
        assert_eq!(round_ratio(&[BigFloat::from_i64(3).neg()],&[one],&two,1),Some(BigFloat::from_i64(-4)));
    }
    #[test]
    fn exact_search_survives_an_incorrect_or_zero_candidate() {
        for p in [1,53,139,140,211] {
            for n in [-51,-3,-1,1,3,51] {
                let a=BigFloat::from_i64(n);let b=BigFloat::from_i64(37);
                let expected=a.div(&b,p).unwrap();
                for candidate in [BigFloat::zero(),BigFloat::from_i64(99)] {
                    assert_eq!(round_ratio(&[a.clone()],&[b.clone()],&candidate,p),Some(expected.clone()));
                }
            }
        }
    }
    #[test]
    fn huge_gaps_remain_exactly_decidable_near_rounding_midpoints() {
        let huge=Exponent::from(BigInt::one()<<190usize);
        let one=BigFloat::one(); let tiny=power(-huge);
        let midpoint=BigFloat::from_parts(BigInt::from(3),-1);
        assert_eq!(round_ratio(&[midpoint.clone(),tiny.clone()],&[one.clone()],&one,1),Some(BigFloat::from_i64(2)));
        assert_eq!(round_ratio(&[midpoint,tiny.neg()],&[one.clone()],&one,1),Some(one));
    }
    #[test]
    fn unsupported_payloads_and_poles_leave_the_existing_path_available() {
        assert!(binary(&Expr::symbol("x"),&Expr::one(),'/',53).is_none());
        assert!(binary(&Expr::one(),&Expr::zero(),'/',53).is_none());
        assert!(binary(&Expr::one(),&Expr::one(),'?',53).is_none());
        assert!(binary(&Expr::one(),&Expr::one(),'/',0).is_none());
    }

}
