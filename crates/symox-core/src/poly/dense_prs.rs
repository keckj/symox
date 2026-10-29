//! Dense univariate Brown--Traub subresultant PRS over ZZ.
//!
//! This is the arithmetic of SymPy's `dup_inner_subresultants` (and of the
//! sparse `prs::subresultants`) on dense integer coefficient vectors: every
//! pseudo-remainder is divided exactly by `-lc * c**d` before the next step,
//! so the intermediates are subresultants and no content is ever removed.
//! The checks only select this kernel: dense degree, retained and transient
//! coefficient storage and a limb-product work estimate. When one of them
//! is exceeded the kernel declines and the sparse PRS runs unchanged with
//! its own limits, so nothing is refused that was computed before.
use super::convert::{PResult, PolyError};
use super::int::{Int, Rat};
use super::mpoly::MPoly;

/// Dense vectors (as `modular_resultant` and `rootisolation`).
const MAX_DEGREE: i64 = 1 << 20;
/// Total coefficient bits of the retained PRS (256 MiB of limbs); every
/// dense entry is also charged a word.
const MAX_BITS: u64 = 1 << 31;
/// Schoolbook limb products (an upper bound of the actual cost).
const MAX_WORK: u64 = 1 << 38;

/// The kernel declined; the caller runs the sparse PRS.
struct Declined;
type DResult<T> = Result<T, Declined>;

struct Work(u64);
impl Work {
    fn charge(&mut self, n: u64) -> DResult<()> {
        self.0 = self.0.checked_sub(n).ok_or(Declined)?;
        Ok(())
    }
    fn product(&mut self, a: &Int, b: &Int) -> DResult<()> {
        self.charge((1 + a.bits() / 64).saturating_mul(1 + b.bits() / 64))
    }
}

fn storage(p: &[Int]) -> u64 {
    p.iter().map(|c| 64 + c.bits()).sum()
}

fn trim(p: &mut Vec<Int>) {
    while p.last().is_some_and(Int::is_zero) {
        p.pop();
    }
}

fn max_bits(p: &[Int]) -> u64 {
    p.iter().map(Int::bits).max().unwrap_or(0)
}

/// `lc(g)**(deg f - deg g + 1) * f mod g` (`deg f >= deg g >= 0`, `g != 0`).
fn prem(f: &[Int], g: &[Int], work: &mut Work) -> DResult<Vec<Int>> {
    let (df, dg) = (f.len() - 1, g.len() - 1);
    let lc = &g[dg];
    let steps = (df - dg + 1) as u64;
    // Coefficients grow by at most bits(lc) + bits(g) + 1 per elimination.
    let grow = lc.bits() + max_bits(g) + 1;
    let estimate = (f.len() as u64)
        .saturating_mul(64 + max_bits(f).saturating_add(steps.saturating_mul(grow)));
    if estimate > MAX_BITS {
        return Err(Declined);
    }
    work.charge(steps.saturating_mul(f.len() as u64))?;
    let mut r = f.to_vec();
    for i in (dg..=df).rev() {
        let c = std::mem::replace(&mut r[i], Int::zero());
        for x in r[..i].iter_mut() {
            if !x.is_zero() {
                work.product(x, lc)?;
                *x = x.mul(lc);
            }
        }
        if !c.is_zero() {
            let shift = i - dg;
            for (k, gk) in g[..dg].iter().enumerate() {
                if !gk.is_zero() {
                    work.product(&c, gk)?;
                    r[k + shift].sub_mul(&c, gk);
                }
            }
        }
    }
    trim(&mut r);
    Ok(r)
}

fn checked_pow(base: &Int, exponent: u64, work: &mut Work) -> DResult<Int> {
    let exponent = u32::try_from(exponent).map_err(|_| Declined)?;
    let bits = base.bits().saturating_mul(u64::from(exponent));
    if bits > MAX_BITS {
        return Err(Declined);
    }
    work.charge((1 + bits / 64).saturating_mul(1 + bits / 64))?;
    Ok(base.pow(exponent))
}

fn exact(a: &Int, b: &Int, work: &mut Work) -> PResult<Result<Int, Declined>> {
    if work.product(a, b).is_err() {
        return Ok(Err(Declined));
    }
    a.div_exact(b)
        .map(Ok)
        .ok_or_else(|| PolyError::Failed("inexact subresultant normalization".into()))
}

/// The PRS and the resultant as `prs::subresultants` returns them, for
/// nonzero dense integer operands (low degree first). `Ok(None)`: declined.
fn subresultants_zz(first: &[Int], second: &[Int]) -> PResult<Option<(Vec<Vec<Int>>, Int)>> {
    macro_rules! take {
        ($e:expr) => {
            match $e {
                Ok(v) => v,
                Err(Declined) => return Ok(None),
            }
        };
    }
    if first.last().map_or(true, Int::is_zero) || second.last().map_or(true, Int::is_zero) {
        return Err(PolyError::Failed("dense subresultants need nonzero operands".into()));
    }
    let mut work = Work(MAX_WORK);
    let (mut f, mut g) = (first.to_vec(), second.to_vec());
    let swapped = f.len() < g.len();
    if swapped {
        std::mem::swap(&mut f, &mut g);
    }
    let reversed_sign = swapped && f.len() % 2 == 0 && g.len() % 2 == 0;
    let mut retained = storage(&f).saturating_add(storage(&g));
    if retained > MAX_BITS {
        return Ok(None);
    }
    let mut degree = g.len() - 1;
    let mut gap = (f.len() - g.len()) as u64;
    let mut h = take!(prem(&f, &g, &mut work));
    if gap % 2 == 0 {
        for x in h.iter_mut() {
            *x = x.neg();
        }
    }
    let mut leading = g[degree].clone();
    let mut c = take!(checked_pow(&leading, gap, &mut work)).neg();
    let mut sequence = vec![f, g];
    while !h.is_empty() {
        let next_degree = h.len() - 1;
        if next_degree >= degree {
            return Err(PolyError::Failed(
                "subresultant PRS did not lower the degree".into(),
            ));
        }
        retained = retained.saturating_add(storage(&h));
        if retained > MAX_BITS {
            return Ok(None);
        }
        sequence.push(h);
        gap = (degree - next_degree) as u64;
        degree = next_degree;
        let n = sequence.len();
        let power = take!(checked_pow(&c, gap, &mut work));
        take!(work.product(&power, &leading));
        let denominator = power.mul(&leading).neg();
        let remainder = take!(prem(&sequence[n - 2], &sequence[n - 1], &mut work));
        let mut quotient = Vec::with_capacity(remainder.len());
        for x in &remainder {
            quotient.push(take!(exact(x, &denominator, &mut work)?));
        }
        h = quotient;
        leading = sequence[n - 1][degree].clone();
        c = if gap > 1 {
            let numerator = take!(checked_pow(&leading.neg(), gap, &mut work));
            let divisor = take!(checked_pow(&c, gap - 1, &mut work));
            take!(exact(&numerator, &divisor, &mut work)?)
        } else {
            leading.neg()
        };
    }
    let resultant = if degree > 0 {
        Int::zero()
    } else if reversed_sign {
        c
    } else {
        c.neg()
    };
    Ok(Some((sequence, resultant)))
}

/// `prs::subresultants` for univariate integral `MPoly<Rat>` operands of
/// positive degree; `None` for every other input and whenever the dense
/// kernel declines.
pub(crate) fn subresultants_univariate(
    f: &MPoly<Rat>,
    g: &MPoly<Rat>,
) -> PResult<Option<(Vec<MPoly<Rat>>, MPoly<Rat>)>> {
    if f.nvars != 1 || g.nvars != 1 || f.degree(0) < 1 || g.degree(0) < 1 {
        return Ok(None);
    }
    if f.degree(0) > MAX_DEGREE || g.degree(0) > MAX_DEGREE {
        return Ok(None);
    }
    if f.coeffs().iter().chain(g.coeffs()).any(|c| !c.is_integer()) {
        return Ok(None);
    }
    let dense = |p: &MPoly<Rat>| -> Vec<Int> {
        let mut out = vec![Int::zero(); p.degree(0) as usize + 1];
        for (e, c) in p.iter() {
            out[e[0] as usize] = c.n.clone();
        }
        out
    };
    let Some((sequence, resultant)) = subresultants_zz(&dense(f), &dense(g))? else {
        return Ok(None);
    };
    let sparse = |p: Vec<Int>| -> MPoly<Rat> {
        let terms = p
            .into_iter()
            .enumerate()
            .rev()
            .filter(|(_, c)| !c.is_zero())
            .map(|(k, c)| (vec![k as u32], Rat::from(c)))
            .collect();
        MPoly::from_terms(1, terms)
    };
    let resultant = if resultant.is_zero() {
        MPoly::zero(1)
    } else {
        MPoly::from_terms(1, vec![(vec![0], Rat::from(resultant))])
    };
    Ok(Some((sequence.into_iter().map(sparse).collect(), resultant)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::polyobj::Poly;
    use crate::{parser::parse_str, Expr};
    fn p(source: &str) -> MPoly<Rat> {
        Poly::from_expr(&parse_str(source).unwrap(), &[Expr::symbol("x")])
            .unwrap()
            .rep
    }
    #[test]
    fn same_sequence_as_the_sparse_prs() {
        let pairs = [
            ("x^4+1", "x^2+1"),
            ("x-1", "x^3+2"),
            ("x^3-2*x", "x^5+3"),
            ("x^6-1", "x^3-1"),
            ("7*x^5-3", "2*x^2"),
            ("x^8+x^6-3*x^4-3*x^3+8*x^2+2*x-5", "3*x^6+5*x^4-4*x^2-9*x+21"),
            ("(x+1)^20", "(x-2)^15+1"),
            ("(x-1)^3*(x+2)*(x^2+3)", "(x-1)^2*(x^2+5)"),
        ];
        for (a, b) in pairs {
            let (f, g) = (p(a), p(b));
            let sparse = crate::poly::prs::subresultants(&f, &g, 0).unwrap();
            let (sequence, resultant) = subresultants_univariate(&f, &g).unwrap().unwrap();
            assert_eq!(sequence, sparse.sequence, "{a}, {b}");
            assert_eq!(resultant, sparse.resultant, "{a}, {b}");
        }
        assert!(subresultants_univariate(&p("x/2+1"), &p("x^2")).unwrap().is_none());
        assert!(subresultants_univariate(&p("3"), &p("x^2")).unwrap().is_none());
        // Sparse giant gaps decline before any dense work.
        let lead = MPoly::from_terms(1, vec![(vec![1 << 19], Rat::from(Int::from(3i64)))]);
        let f = lead.add(&p("x+1"));
        assert!(subresultants_univariate(&f, &p("2*x^2+1")).unwrap().is_none());
    }
}
