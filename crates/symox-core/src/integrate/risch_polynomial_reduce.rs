//! Exact polynomial reduction for an admitted rational derivation D(t)=d(t).
//! This is the polynomial_reduce loop, not the whole Risch integrator.
use crate::poly::{dense, Int, Rat};

#[derive(Clone, Debug, PartialEq)]
pub struct Reduced {
    pub quotient: Vec<Rat>,
    pub remainder: Vec<Rat>,
    /// The historical constructor starts q over ZZ and promotes on a
    /// nonintegral monomial; p retains its domain through each subtraction.
    pub quotient_qq: bool,
    pub remainder_qq: bool,
    pub steps: usize,
}

/// None declines a zero derivation; no partial result is returned. Every
/// accepted step cancels the leading term exactly and lowers the degree.
pub fn reduce<E>(input: &[Rat], derivation: &[Rat], input_qq: bool,
    checkpoint: &mut impl FnMut() -> Result<(), E>) -> Result<Option<Reduced>, E>
{
    checkpoint()?;
    let mut d = derivation.to_vec();
    dense::trim(&mut d);
    if d.is_empty() { return Ok(None); }
    let mut p = input.to_vec();
    dense::trim(&mut p);
    let mut q = Vec::<Rat>::new();
    let mut q_qq = false;
    let mut p_qq = input_qq;
    let mut steps = 0;
    while !p.is_empty() && p.len() >= d.len() {
        checkpoint()?;
        // m = degree(p)-degree(d)+1 is always positive in this loop.
        let m = p.len()-d.len()+1;
        let degree = Rat::from(Int::from_u64(m as u64));
        let coefficient = dense::lc(&p).div(&degree.mul(&dense::lc(&d)));
        q_qq |= !coefficient.d.is_one();
        if q.len() <= m { q.resize(m+1, Rat::zero()); }
        q[m] = q[m].add(&coefficient);
        // D(c*t^m) = c*m*t^(m-1)*d; c is a QQ constant.
        let scale = coefficient.mul(&degree);
        let mut term = vec![Rat::zero(); m-1];
        for (i, value) in d.iter().enumerate() {
            if i % 64 == 0 { checkpoint()?; }
            let value = value.mul(&scale);
            p_qq |= !value.d.is_one();
            term.push(value);
        }
        let previous_length = p.len();
        p = dense::sub(&p, &term);
        debug_assert!(p.len() < previous_length);
        steps += 1;
    }
    dense::trim(&mut q);
    checkpoint()?;
    Ok(Some(Reduced { quotient:q, remainder:p, quotient_qq:q_qq,
        remainder_qq:p_qq, steps }))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rational(n:i64,d:i64)->Rat { Rat::new(Int::from(n),Int::from(d)) }
    fn certificate(p:&[Rat],d:&[Rat],input_qq:bool) -> Reduced {
        let result=reduce(p,d,input_qq,&mut ||Ok::<_,()>(())).unwrap().unwrap();
        let recomposed=dense::add(&dense::mul(&dense::derivative(&result.quotient),d),&result.remainder);
        let mut expected=p.to_vec();dense::trim(&mut expected);
        assert_eq!(recomposed,expected);
        assert!(result.remainder.is_empty()||result.remainder.len()<d.len());
        result
    }
    #[test]
    fn exact_qq_reduction_has_complete_derivative_and_degree_certificates() {
        for degree in 0..=35 {
            for order in 0..=7 {
                let p=(0..degree).map(|i|rational((i*17+degree*3)%29-14,(i%7)+1)).collect::<Vec<_>>();
                let mut d=(0..=order).map(|i|rational((i*13+order)%17-8,(i%5)+1)).collect::<Vec<_>>();
                d[order as usize]=rational(order+1,order+2);
                certificate(&p,&d,true);
            }
        }
    }
    #[test]
    fn constructor_domain_promotions_are_separate_for_quotient_and_remainder() {
        let q=certificate(&[Rat::zero(),Rat::zero(),Rat::one()],&[Rat::one(),Rat::one()],false);
        assert!(q.quotient_qq);assert!(!q.remainder_qq);
        let q=certificate(&[Rat::one(),Rat::zero(),Rat::one()],&[rational(1,2),Rat::one()],false);
        assert!(q.quotient_qq);assert!(q.remainder_qq);
        let q=certificate(&[Rat::one()],&[Rat::one(),Rat::one()],true);
        assert!(!q.quotient_qq);assert!(q.remainder_qq);assert_eq!(q.steps,0);
    }
    #[test]
    fn zero_derivation_declines_and_interruptions_do_not_publish_partial_results() {
        assert!(reduce(&[Rat::one()],&[],false,&mut ||Ok::<_,()>(())).unwrap().is_none());
        let mut count=0;
        let output=reduce(&[Rat::one(),Rat::one(),Rat::one()],&[Rat::one()],false,&mut ||{
            count+=1;if count==3 {Err("interrupt")}else{Ok(())}
        });
        assert_eq!(output,Err("interrupt"));
    }
}
