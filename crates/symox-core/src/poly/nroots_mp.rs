//! Precision-preserving numerical roots of rational polynomials.
//!
//! Exact square-free decomposition handles multiplicities before simultaneous
//! Aberth iteration. Machine floats are only optional starting guesses; every
//! correction, convergence decision and returned component uses BigFloat.

use super::{dense, int::{Int, Rat}, roots, sqf};
use super::convert::{PResult, PolyError};
use crate::bigfloat::{BigFloat, dps_to_prec};
use crate::evalf_mp::{MpComplex, to_expr};
use crate::expr::Expr;

pub fn rational_roots(c: &[Rat], digits: u32, maxsteps: u32, cleanup: bool) -> PResult<Vec<Expr>> {
    let precision = dps_to_prec(digits);
    let (_, mut integers) = dense::clear_denoms(c);
    dense::trim(&mut integers);
    if integers.len() <= 1 { return Ok(vec![]); }
    let zeros = integers.iter().take_while(|a| a.is_zero()).count();
    let mut values = vec![MpComplex::real(BigFloat::zero()); zeros];
    for (factor, multiplicity) in sqf::sqf_list_zx(&integers[zeros..]).1 {
        let roots = simple_roots(&factor, precision, maxsteps, cleanup).ok_or_else(||
            PolyError::Failed(format!("convergence to root failed; try n < {digits} or maxsteps > {maxsteps}")))?;
        for root in roots {
            values.extend(std::iter::repeat(root).take(multiplicity as usize));
        }
    }
    values.sort_by(|a,b| {
        (!a.im.is_zero()).cmp(&(!b.im.is_zero()))
            .then(a.re.cmp(&b.re))
            .then(a.im.abs().cmp(&b.im.abs()))
            .then(a.im.cmp(&b.im))
    });
    Ok(values.into_iter().map(|z| {
        // The reference builds roots with Float._new, which maps an exact
        // zero mantissa to S.Zero.
        if z.re.is_zero() && z.im.is_zero() { Expr::zero() } else { to_expr(z, precision) }
    }).collect())
}

/// Numerical roots of a polynomial with Gaussian rational coefficients, given
/// by its exact square-free factors `(coefficients (re, im) ascending,
/// multiplicity)`. Every factor is solved with exact Gaussian integer
/// coefficients at the requested precision, like the rational case.
pub fn gaussian_roots(factors: &[(Vec<(Rat, Rat)>, u32)], digits: u32, maxsteps: u32, cleanup: bool) -> PResult<Vec<Expr>> {
    let precision = dps_to_prec(digits);
    let mut values = Vec::new();
    for (factor, multiplicity) in factors {
        let mut common = Int::one();
        for (x, y) in factor { common = common.lcm(&x.d).lcm(&y.d); }
        let scale = |q: &Rat| q.n.mul(&common).div_exact_unchecked(&q.d);
        let mut c: Vec<(Int, Int)> = factor.iter().map(|(x, y)| (scale(x), scale(y))).collect();
        while c.last().is_some_and(|(x, y)| x.is_zero() && y.is_zero()) { c.pop(); }
        let zeros = c.iter().take_while(|(x, y)| x.is_zero() && y.is_zero()).count();
        values.extend(std::iter::repeat(MpComplex::real(BigFloat::zero())).take(zeros * *multiplicity as usize));
        let c = c.split_off(zeros);
        if c.len() <= 1 { continue; }
        let roots = simple_complex_roots(&c, precision, maxsteps, cleanup).ok_or_else(||
            PolyError::Failed(format!("convergence to root failed; try n < {digits} or maxsteps > {maxsteps}")))?;
        for root in roots {
            values.extend(std::iter::repeat(root).take(*multiplicity as usize));
        }
    }
    values.sort_by(|a,b| {
        (!a.im.is_zero()).cmp(&(!b.im.is_zero()))
            .then(a.re.cmp(&b.re))
            .then(a.im.abs().cmp(&b.im.abs()))
            .then(a.im.cmp(&b.im))
    });
    Ok(values.into_iter().map(|z| {
        // The reference builds roots with Float._new, which maps an exact
        // zero mantissa to S.Zero.
        if z.re.is_zero() && z.im.is_zero() { Expr::zero() } else { to_expr(z, precision) }
    }).collect())
}

fn horner_complex(c: &[MpComplex], z: &MpComplex, precision: u32) -> (MpComplex, MpComplex) {
    let mut value = c.last().unwrap().clone();
    let mut derivative = MpComplex::real(BigFloat::zero());
    for coefficient in c[..c.len()-1].iter().rev() {
        derivative = derivative.mul(z, precision).add(&value, precision);
        value = value.mul(z, precision).add(coefficient, precision);
    }
    (value, derivative)
}

fn simple_complex_roots(c: &[(Int, Int)], precision: u32, maxsteps: u32, cleanup: bool) -> Option<Vec<MpComplex>> {
    if maxsteps == 0 { return None; }
    let degree = c.len()-1;
    let height = c.iter().map(|(x, y)| x.to_big().bits().max(y.to_big().bits())).max().unwrap_or(0);
    let guard = u32::try_from(height).ok()?.checked_mul(2)?.checked_add(64)?;
    let working = precision.checked_add(guard.max(degree as u32 * 8))?;
    for attempt in 0..3 {
        let wp = working.checked_add(attempt * guard)?;
        let coeffs: Vec<MpComplex> = c.iter().map(|(x, y)| MpComplex {
            re: BigFloat::from_int(&x.to_big()).round_to(wp),
            im: BigFloat::from_int(&y.to_big()).round_to(wp),
        }).collect();
        let mut values = if degree == 1 {
            vec![coeffs[0].neg().div(&coeffs[1], wp)?]
        } else {
            // The Newton polygon seeds only use the coefficient magnitudes.
            let magnitudes: Vec<BigFloat> = coeffs.iter().map(|z| z.abs(wp)).collect();
            let mut values = polygon_seeds(&magnitudes, 0.37 + attempt as f64*0.11, wp);
            if !iterate_complex(&coeffs, &mut values, precision, wp, maxsteps) { continue; }
            values
        };
        if cleanup {
            for z in &mut values {
                let cutoff = z.re.mag().max(z.im.mag()) - precision as i128 - 8;
                if z.im.mag() < cutoff { z.im = BigFloat::zero(); }
                if !z.im.is_zero() && z.re.mag() < cutoff { z.re = BigFloat::zero(); }
            }
        }
        return Some(values.into_iter().map(|z| z.round_to(precision)).collect());
    }
    None
}

fn iterate_complex(c: &[MpComplex], roots: &mut [MpComplex], target: u32, wp: u32, maxsteps: u32) -> bool {
    let one = MpComplex::real(BigFloat::one());
    for _ in 0..maxsteps {
        let mut converged = true;
        for i in 0..roots.len() {
            let (value,derivative) = horner_complex(c,&roots[i],wp);
            let Some(ratio) = value.div(&derivative,wp) else { return false; };
            let mut sum = MpComplex::real(BigFloat::zero());
            let mut separation = i128::MAX;
            for j in 0..roots.len() {
                if i == j { continue; }
                let difference = roots[i].sub(&roots[j],wp);
                let Some(inverse) = difference.inv(wp) else { return false; };
                separation = separation.min(difference.re.mag().max(difference.im.mag()));
                sum = sum.add(&inverse,wp);
            }
            let Some(step) = ratio.div(&one.sub(&ratio.mul(&sum,wp),wp),wp) else { return false; };
            roots[i] = roots[i].sub(&step,wp);
            let error = step.re.mag().max(step.im.mag());
            let magnitude = roots[i].re.mag().max(roots[i].im.mag());
            if error > magnitude - target as i128 - 8 || error > separation - 12 { converged = false; }
        }
        if converged { return true; }
    }
    false
}

fn horner(c: &[BigFloat], z: &MpComplex, precision: u32) -> (MpComplex, MpComplex) {
    let mut value = MpComplex::real(c.last().unwrap().clone());
    let mut derivative = MpComplex::real(BigFloat::zero());
    for coefficient in c[..c.len()-1].iter().rev() {
        derivative = derivative.mul(z, precision).add(&value, precision);
        value = value.mul(z, precision).add(&MpComplex::real(coefficient.clone()), precision);
    }
    (value, derivative)
}

/// Upper Newton polygon: the horizontal length of a segment counts roots at
/// approximately 2**(-slope). This keeps seeds representable for coefficients
/// and roots far outside the machine-float exponent range.
fn polygon_seeds(c: &[BigFloat], phase: f64, precision: u32) -> Vec<MpComplex> {
    let mut hull: Vec<(usize, i128)> = Vec::new();
    for (i, coefficient) in c.iter().enumerate().filter(|(_,a)| !a.is_zero()) {
        let point = (i, coefficient.mag());
        while hull.len() >= 2 {
            let (a,b) = (hull[hull.len()-2], hull[hull.len()-1]);
            let left = (b.1-a.1) * (point.0-b.0) as i128;
            let right = (point.1-b.1) * (b.0-a.0) as i128;
            if left > right { break; }
            hull.pop();
        }
        hull.push(point);
    }
    let mut seeds = Vec::with_capacity(c.len()-1);
    for segment in hull.windows(2) {
        let (a,b) = (segment[0],segment[1]);
        let count = b.0-a.0;
        let exponent = (a.1-b.1).div_euclid(count as i128);
        let radius = BigFloat::one().ldexp(exponent);
        for i in 0..count {
            let angle = std::f64::consts::TAU * (i as f64 + phase) / count as f64;
            seeds.push(MpComplex {
                re: radius.mul(&BigFloat::from_f64(angle.cos()),precision),
                im: radius.mul(&BigFloat::from_f64(angle.sin()),precision),
            });
        }
    }
    seeds
}

fn machine_seeds(c: &[Int], maxsteps: u32) -> Option<Vec<MpComplex>> {
    let leading = c.last()?.to_f64();
    if !leading.is_finite() || leading == 0.0 { return None; }
    let coefficients: Vec<_> = c.iter().map(|a| a.to_f64()/leading).collect();
    if coefficients.iter().any(|a| !a.is_finite() || (*a != 0.0 && a.abs() < 1e-200))
        || coefficients.iter().zip(c).any(|(a,b)| *a == 0.0 && !b.is_zero()) { return None; }
    let seeds = roots::aberth_steps(&coefficients, maxsteps.min(16) as usize);
    for (i,z) in seeds.iter().enumerate() {
        if !z.re.is_finite() || !z.im.is_finite() || z.norm() == 0.0 { return None; }
        if seeds[..i].iter().any(|other| (*z-*other).norm() < 1e-8*z.norm().max(other.norm())) { return None; }
    }
    Some(seeds.into_iter().map(|z| MpComplex { re: BigFloat::from_f64(z.re), im: BigFloat::from_f64(z.im) }).collect())
}

fn simple_roots(c: &[Int], precision: u32, maxsteps: u32, cleanup: bool) -> Option<Vec<MpComplex>> {
    if maxsteps == 0 { return None; }
    let degree = c.len()-1;
    // Coefficient height guards cancellation in close roots. Increasing the
    // work precision on retry also avoids declaring stagnation a success.
    let height = c.iter().map(|a| a.to_big().bits()).max().unwrap_or(0);
    let guard = u32::try_from(height).ok()?.checked_mul(2)?.checked_add(64)?;
    let working = precision.checked_add(guard.max(degree as u32 * 8))?;
    for attempt in 0..2 {
        let wp = working.checked_add(attempt * guard)?;
        let coeffs: Vec<_> = c.iter().map(|a| BigFloat::from_int(&a.to_big()).round_to(wp)).collect();
        let initial = if attempt == 0 { machine_seeds(c,maxsteps) } else { None };
        let mut values = initial.unwrap_or_else(|| polygon_seeds(&coeffs,0.37 + attempt as f64*0.11,wp));
        if degree == 1 {
            values = vec![MpComplex::real(coeffs[0].neg().div(&coeffs[1],wp)?)];
        } else if !iterate(&coeffs, &mut values, precision, wp, maxsteps) { continue; }
        if cleanup {
            for z in &mut values {
                // Scale cleanup to the root, so a genuine 1e-1000*I root is
                // never erased merely because the requested precision is 60.
                let cutoff = z.re.mag().max(z.im.mag()) - precision as i128 - 8;
                if z.im.mag() < cutoff { z.im = BigFloat::zero(); }
                if !z.im.is_zero() && z.re.mag() < cutoff { z.re = BigFloat::zero(); }
            }
        }
        return Some(values.into_iter().map(|z| z.round_to(precision)).collect());
    }
    None
}

fn iterate(c: &[BigFloat], roots: &mut [MpComplex], target: u32, wp: u32, maxsteps: u32) -> bool {
    let one = MpComplex::real(BigFloat::one());
    for _ in 0..maxsteps {
        let mut converged = true;
        for i in 0..roots.len() {
            let (value,derivative) = horner(c,&roots[i],wp);
            let Some(ratio) = value.div(&derivative,wp) else { return false; };
            let mut sum = MpComplex::real(BigFloat::zero());
            let mut separation = i128::MAX;
            for j in 0..roots.len() {
                if i == j { continue; }
                let difference = roots[i].sub(&roots[j],wp);
                let Some(inverse) = difference.inv(wp) else { return false; };
                separation = separation.min(difference.re.mag().max(difference.im.mag()));
                sum = sum.add(&inverse,wp);
            }
            let Some(step) = ratio.div(&one.sub(&ratio.mul(&sum,wp),wp),wp) else { return false; };
            roots[i] = roots[i].sub(&step,wp);
            let error = step.re.mag().max(step.im.mag());
            let magnitude = roots[i].re.mag().max(roots[i].im.mag());
            if error > magnitude - target as i128 - 8 || error > separation - 12 { converged = false; }
        }
        if converged { return true; }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{evalf_mp::eval_mp, parser::parse_str, poly::polyobj::Poly};

    fn coefficients(source: &str) -> Vec<Rat> {
        let polynomial = Poly::from_expr(&parse_str(source).unwrap(), &[Expr::symbol("x")]).unwrap();
        dense::from_mpoly(&polynomial.rep,0)
    }

    #[test]
    fn requested_precision_is_computed_not_padded() {
        for digits in [15,30,60,100] {
            let c = coefficients("x^10-3*x^3+x^2-1");
            let roots = rational_roots(&c,digits,200,true).unwrap();
            assert_eq!(roots.len(),10);
            let precision = dps_to_prec(digits);
            let coeffs: Vec<_> = c.iter().map(|c| BigFloat::from_ratio(&c.n.to_big(),&c.d.to_big(),precision+60)).collect();
            for root in roots {
                let z = eval_mp(&root,precision+60).unwrap();
                let (residual,derivative) = horner(&coeffs,&z,precision+60);
                let correction = residual.div(&derivative,precision+40).unwrap();
                assert!(correction.abs(precision+20).mag() <= z.abs(precision+20).mag()-precision as i128+4);
            }
        }
    }

    #[test]
    fn multiplicities_zero_and_extreme_exponents_survive() {
        let roots = rational_roots(&coefficients("x^3*(x-1)^4*(x^2+1)^2"),60,100,true).unwrap();
        assert_eq!(roots.len(),11);
        assert_eq!(roots.iter().filter(|r|r.is_zero()).count(),3);
        assert_eq!(roots.iter().filter(|r|eval_mp(r,210).unwrap().re == BigFloat::one()).count(),4);
        for source in ["x^2-1/10^1000", "x^2+1/10^1000", "x^2-10^1000*x+1"] {
            let roots = rational_roots(&coefficients(source),60,200,true).unwrap();
            assert_eq!(roots.len(),2);
            assert!(roots.iter().all(|r|!r.is_zero()),"{source}");
        }
    }

    #[test]
    fn exhausted_iteration_budget_is_an_error() {
        let c = coefficients("x^10-3*x^3+x^2-1");
        assert!(rational_roots(&c,100,0,true).is_err());
        assert!(rational_roots(&c,100,1,true).is_err());
    }
}
