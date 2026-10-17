//! Dispersion sets of univariate rational polynomials
//! (`sympy.polys.dispersion.dispersionset`).
//!
//! `J(p, q)` is the set of integers `a >= 0` with `deg gcd(p(x), q(x + a)) >= 1`.
//! As in SymPy, the irreducible factors `s` of `p` and `t` of `q` are paired:
//! `s(x) = t(x + a)` forces equal degree `n` and leading coefficient (the
//! factors are primitive with positive leading coefficient) and
//! `a = (s_{n-1} - t_{n-1}) / (n t_n)`; for `n > 1` the shift is verified.
use super::dense;
use super::int::{Int, Rat};
use super::polyobj::Poly;

fn irreducible_factors(p: &Poly) -> Vec<Vec<Rat>> {
    p.factor_list().1.into_iter().map(|(f, _)| dense::from_mpoly(&f.rep, 0)).collect()
}

/// The dispersion set of two univariate polynomials over `ZZ`/`QQ` of
/// positive degree, sorted increasingly (`q = None` means `q = p`).
pub fn dispersion_set(p: &Poly, q: Option<&Poly>) -> Vec<Int> {
    let fp = irreducible_factors(p);
    let fq = match q {
        Some(q) => irreducible_factors(q),
        None => fp.clone(),
    };
    let mut out: Vec<Int> = Vec::new();
    for s in &fp {
        for t in &fq {
            let n = s.len() - 1;
            if n != t.len() - 1 || n == 0 || s[n] != t[n] {
                continue;
            }
            let alpha = s[n - 1].sub(&t[n - 1]).div(&Rat::from_i64(n as i64).mul(&t[n]));
            if !alpha.is_integer() || alpha.is_negative() || out.contains(&alpha.n) {
                continue;
            }
            if n > 1 {
                let shifted = dense::from_mpoly(&dense::to_mpoly(t, 0, 1).shift(0, &alpha), 0);
                if &shifted != s {
                    continue;
                }
            }
            out.push(alpha.n.clone());
        }
    }
    out.sort();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn poly(s: &str) -> Poly {
        Poly::from_expr(&parse_str(s).unwrap(), &[crate::Expr::symbol("x")]).unwrap()
    }

    #[test]
    fn sympy_examples() {
        let p = poly("(x - 3)*(x + 3)");
        assert_eq!(dispersion_set(&p, None), vec![Int::from_u64(0), Int::from_u64(6)]);
        let p = poly("x**4 - 3*x**2 + 1");
        let q = poly("x**4 - 3*x**2 + 1").shift(&Rat::from_i64(-3));
        // x**4 - 3*x**2 + 1 = (x**2 - x - 1)*(x**2 + x - 1), shifts of each other by 1.
        assert_eq!(dispersion_set(&p, Some(&q)), [2u64, 3, 4].map(Int::from_u64).to_vec());
        assert_eq!(dispersion_set(&poly("x**2 + 1"), None), vec![Int::from_u64(0)]);
    }
}
