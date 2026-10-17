//! Bronstein's full partial fraction decomposition over `QQ`
//! (`sympy.polys.partfrac.apart_list_full_decomposition`).
//!
//! For a square-free factor `d` of multiplicity `n` of the denominator `Q`,
//! write `Q = q * d**n`. Near a root `a` of `d`,
//! `P/Q = sum_{k=1..n} c_k(a) / (x - a)**k + (regular part)`, and
//! `c_{n-j}(a)` is the `j`-th Taylor coefficient at `t = 0` of
//! `P(a + t) / (q(a + t) * e(a + t)**n)` with `e(x) = d(x)/(x - a)`.
//! The series is computed in `QQ[a]/(d(a))` (a product of fields since `d` is
//! square-free; the constant term is invertible because `q(a) d'(a) != 0`).
//! SymPy's output keeps only the roots where `c` does not vanish:
//! `D = d / gcd(c, d)` and the numerator is `c mod D`, which is the same
//! unique remainder SymPy obtains through its derivative substitutions.
use super::convert::{PResult, PolyError};
use super::dense;
use super::int::Rat;

/// One `(D, numer, exponent)` entry of the decomposition, with dense
/// ascending coefficients: the term is
/// `RootSum(D, Lambda(a, numer(a) / (x - a)**exponent))`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PartialTerm {
    pub denom: Vec<Rat>,
    pub numer: Vec<Rat>,
    pub exponent: usize,
}

const MAX_MULTIPLICITY: usize = 4096;

fn unsupported(what: &str) -> PolyError {
    PolyError::Unsupported(format!("native partial fractions: {what}"))
}

/// `F^(k)(x) / k!` for `k = 0..count` reduced modulo `d`.
fn taylor_mod(f: &[Rat], count: usize, d: &[Rat]) -> Vec<Vec<Rat>> {
    let mut out = Vec::with_capacity(count);
    let mut h = f.to_vec();
    for k in 0..count {
        if k > 0 {
            h = dense::scale(&dense::derivative(&h), &Rat::from_i64(k as i64).inv());
        }
        out.push(dense::rem(&h, d));
    }
    out
}

fn mul_mod(a: &[Rat], b: &[Rat], d: &[Rat]) -> Vec<Rat> {
    dense::rem(&dense::mul(a, b), d)
}

/// Truncated product of two series with coefficients in `QQ[a]/(d)`.
fn series_mul(a: &[Vec<Rat>], b: &[Vec<Rat>], n: usize, d: &[Rat]) -> Vec<Vec<Rat>> {
    (0..n)
        .map(|k| {
            let mut s = Vec::new();
            for i in 0..=k {
                if i < a.len() && k - i < b.len() {
                    s = dense::add(&s, &mul_mod(&a[i], &b[k - i], d));
                }
            }
            s
        })
        .collect()
}

/// The primitive integer associate (positive leading coefficient) of a
/// non-zero rational polynomial.
fn primitive_rat(p: &[Rat]) -> Vec<Rat> {
    let (_, ints) = dense::clear_denoms(p);
    dense::to_rat(&dense::primitive(&ints).1)
}

/// The decomposition of `P/Q` from `P` and the factors `(d, n)` of
/// `Q.sqf_list_include(all=True)`. `integral` selects SymPy's `ZZ` gcd
/// normalization (primitive with a positive leading coefficient) instead of
/// the monic `QQ` one. Constant factors must be units over `ZZ`.
pub fn full_decomposition(
    p: &[Rat],
    factors: &[(Vec<Rat>, usize)],
    integral: bool,
) -> PResult<Vec<PartialTerm>> {
    let mut out = Vec::new();
    for (i, (d, n)) in factors.iter().enumerate() {
        let n = *n;
        if d.is_empty() || n == 0 {
            return Err(unsupported("zero factor"));
        }
        if n > MAX_MULTIPLICITY {
            return Err(PolyError::Failed(
                "polynomial structure resource limit: multiplicity".into(),
            ));
        }
        if integral && (d.iter().any(|c| !c.is_integer()) || !dense::content(&dense::clear_denoms(d).1).is_one()) {
            return Err(unsupported("non-primitive integer factor"));
        }
        if d.len() == 1 {
            // gcd(Pj, c) is 1 for a unit c: D = d and the numerator is zero.
            if integral && !d[0].abs().is_one() {
                return Err(unsupported("constant factor"));
            }
            for j in 0..n {
                out.push(PartialTerm { denom: d.clone(), numer: vec![], exponent: n - j });
            }
            continue;
        }
        let mut q = vec![Rat::one()];
        for (k, (dk, nk)) in factors.iter().enumerate() {
            if k != i {
                for _ in 0..*nk {
                    q = dense::mul(&q, dk);
                }
            }
        }
        let num = taylor_mod(p, n, d);
        // e(a + t) = sum_{k >= 1} d^(k)(a)/k! t^(k-1)
        let e: Vec<Vec<Rat>> = taylor_mod(d, n + 1, d).into_iter().skip(1).collect();
        let mut den = taylor_mod(&q, n, d);
        let mut power = vec![vec![Rat::one()]];
        let (mut base, mut k) = (e, n);
        while k > 0 {
            if k & 1 == 1 {
                power = series_mul(&power, &base, n, d);
            }
            k >>= 1;
            if k > 0 {
                base = series_mul(&base, &base, n, d);
            }
        }
        den = series_mul(&den, &power, n, d);
        let (g, s, _) = dense::gcdex_field(&den[0], d);
        if g.len() != 1 {
            return Err(PolyError::Failed("singular partial fraction denominator".into()));
        }
        let inv0 = s;
        let mut c: Vec<Vec<Rat>> = Vec::with_capacity(n);
        for k in 0..n {
            let mut acc = num[k].clone();
            for i in 1..=k {
                acc = dense::sub(&acc, &mul_mod(&den[i], &c[k - i], d));
            }
            c.push(mul_mod(&acc, &inv0, d));
        }
        for (j, cj) in c.into_iter().enumerate() {
            let mut g = dense::gcd_field(&cj, d);
            if integral {
                g = primitive_rat(&g);
            }
            let (big, r) = dense::divrem(d, &g);
            debug_assert!(r.is_empty());
            let numer = dense::rem(&cj, &big);
            out.push(PartialTerm { denom: big, numer, exponent: n - j });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(v: &[i64]) -> Vec<Rat> {
        v.iter().map(|&c| Rat::from_i64(c)).collect()
    }

    #[test]
    fn simple_and_repeated_roots() {
        // 1/(x - 2)**2 = 1/(x - 2)**2: one exponent-2 term, a zero exponent-1 term.
        let t = full_decomposition(&r(&[1]), &[(r(&[1]), 1), (r(&[-2, 1]), 2)], true).unwrap();
        assert_eq!(t[0], PartialTerm { denom: r(&[1]), numer: vec![], exponent: 1 });
        assert_eq!(t[1], PartialTerm { denom: r(&[-2, 1]), numer: r(&[1]), exponent: 2 });
        assert_eq!(t[2], PartialTerm { denom: r(&[1]), numer: vec![], exponent: 1 });
        // 2x/(x**2 - 1): RootSum(w**2 - 1, a -> 1/(x - a)).
        let t = full_decomposition(&r(&[0, 2]), &[(r(&[-1, 0, 1]), 1)], true).unwrap();
        assert_eq!(t, vec![PartialTerm { denom: r(&[-1, 0, 1]), numer: r(&[1]), exponent: 1 }]);
    }

    #[test]
    fn vanishing_coefficients_shrink_the_root_set() {
        // x/(x**2 - 1)**2 = 1/4 (1/(x-1)**2 - 1/(x+1)**2): exponent 1 vanishes.
        let t = full_decomposition(&r(&[0, 1]), &[(r(&[1]), 1), (r(&[-1, 0, 1]), 2)], true).unwrap();
        let q = Rat::from_i64(1).div(&Rat::from_i64(4));
        assert_eq!(t[1], PartialTerm { denom: r(&[-1, 0, 1]), numer: vec![Rat::zero(), q], exponent: 2 });
        assert_eq!(t[2], PartialTerm { denom: r(&[1]), numer: vec![], exponent: 1 });
    }
}
