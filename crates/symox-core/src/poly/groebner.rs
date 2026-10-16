//! Deterministic Buchberger bases over QQ and prime finite fields.
//!
//! MPoly retains its lexicographic storage invariant. All leading-term
//! operations here use an explicit admissible monomial order.

use std::cmp::Ordering;
use std::collections::VecDeque;

use super::int::Rat;
use super::mpoly::MPoly;
use super::polyobj::reduce_mod;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MonomialOrder {
    Lex,
    Grlex,
    Grevlex,
}

impl MonomialOrder {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "lex" => Some(Self::Lex),
            "grlex" => Some(Self::Grlex),
            "grevlex" => Some(Self::Grevlex),
            _ => None,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::Lex => "lex",
            Self::Grlex => "grlex",
            Self::Grevlex => "grevlex",
        }
    }
    pub fn compare(self, a: &[u32], b: &[u32]) -> Ordering {
        if self != Self::Lex {
            let total = a
                .iter()
                .map(|&x| x as u64)
                .sum::<u64>()
                .cmp(&b.iter().map(|&x| x as u64).sum::<u64>());
            if total != Ordering::Equal {
                return total;
            }
        }
        if self == Self::Grevlex {
            b.iter().rev().cmp(a.iter().rev())
        } else {
            a.cmp(b)
        }
    }
    pub fn leading_index(self, p: &MPoly<Rat>) -> usize {
        (0..p.len())
            .max_by(|&i, &j| self.compare(p.exp(i), p.exp(j)))
            .expect("nonzero polynomial")
    }
}

fn divides(a: &[u32], b: &[u32]) -> bool {
    a.iter().zip(b).all(|(x, y)| x <= y)
}
fn difference(a: &[u32], b: &[u32]) -> Vec<u32> {
    a.iter().zip(b).map(|(x, y)| x - y).collect()
}

/// Arithmetic over QQ, or GF(p) with symmetric rational-integer storage.
#[derive(Clone, Copy, Debug)]
pub struct CoefficientField {
    pub modulus: Option<u64>,
}
impl CoefficientField {
    pub fn normalize(self, p: &MPoly<Rat>) -> MPoly<Rat> {
        match self.modulus {
            Some(m) => reduce_mod(p, m),
            None => p.clone(),
        }
    }
    fn quotient(self, a: &Rat, b: &Rat) -> Rat {
        match self.modulus {
            None => a.div(b),
            Some(p) => {
                let aa = a.n.rem_u64(p);
                let bb = b.n.rem_u64(p);
                let c = super::gf::mulmod(aa, super::gf::invmod(bb, p), p);
                Rat::from_i64(if c > p / 2 {
                    c as i64 - p as i64
                } else {
                    c as i64
                })
            }
        }
    }
    pub fn monic(self, p: &MPoly<Rat>, order: MonomialOrder) -> MPoly<Rat> {
        if p.is_zero() {
            return p.clone();
        }
        self.normalize(&p.scale(&self.quotient(&Rat::one(), p.coeff(order.leading_index(p)))))
    }
}

/// Multivariate division: `f = sum(q[i]*basis[i]) + remainder`.
pub fn reduce(
    f: &MPoly<Rat>,
    basis: &[MPoly<Rat>],
    order: MonomialOrder,
    field: CoefficientField,
) -> (Vec<MPoly<Rat>>, MPoly<Rat>) {
    reduce_with_options(f, basis, order, field, false)
}

/// As `reduce`, optionally requiring quotient coefficients to remain in ZZ.
pub fn reduce_with_options(
    f: &MPoly<Rat>,
    basis: &[MPoly<Rat>],
    order: MonomialOrder,
    field: CoefficientField,
    integral: bool,
) -> (Vec<MPoly<Rat>>, MPoly<Rat>) {
    let n = f.nvars;
    let mut work = field.normalize(f);
    let mut q = vec![MPoly::zero(n); basis.len()];
    let mut r = MPoly::zero(n);
    let leads: Vec<_> = basis
        .iter()
        .map(|g| {
            if g.is_zero() {
                None
            } else {
                Some(order.leading_index(g))
            }
        })
        .collect();
    while !work.is_zero() {
        let i = order.leading_index(&work);
        let exps = work.exp(i).to_vec();
        let coeff = work.coeff(i).clone();
        let divisor = basis.iter().enumerate().find(|(j, g)| {
            leads[*j].is_some_and(|k| {
                divides(g.exp(k), &exps)
                    && (!integral || field.quotient(&coeff, g.coeff(k)).is_integer())
            })
        });
        if let Some((j, g)) = divisor {
            let k = leads[j].unwrap();
            let t = MPoly::monomial(
                &difference(&exps, g.exp(k)),
                field.quotient(&coeff, g.coeff(k)),
            );
            q[j] = field.normalize(&q[j].add(&t));
            work = field.normalize(&work.sub(&t.mul(g)));
        } else {
            let t = MPoly::monomial(&exps, coeff);
            r = r.add(&t);
            work = work.sub(&t);
        }
    }
    (q, r)
}

pub fn s_polynomial(
    f: &MPoly<Rat>,
    g: &MPoly<Rat>,
    order: MonomialOrder,
    field: CoefficientField,
) -> MPoly<Rat> {
    let i = order.leading_index(f);
    let j = order.leading_index(g);
    let lcm: Vec<_> = f
        .exp(i)
        .iter()
        .zip(g.exp(j))
        .map(|(&a, &b)| a.max(b))
        .collect();
    let a = f.mul_term(
        &difference(&lcm, f.exp(i)),
        &field.quotient(&Rat::one(), f.coeff(i)),
    );
    let b = g.mul_term(
        &difference(&lcm, g.exp(j)),
        &field.quotient(&Rat::one(), g.coeff(j)),
    );
    field.normalize(&a.sub(&b))
}

/// A reduced monic Groebner basis, sorted by decreasing leading monomial.
pub fn buchberger(
    polys: &[MPoly<Rat>],
    order: MonomialOrder,
    field: CoefficientField,
) -> Vec<MPoly<Rat>> {
    let mut basis = Vec::new();
    for f in polys {
        let r = reduce(f, &basis, order, field).1;
        if !r.is_zero() {
            let r = field.monic(&r, order);
            if r.is_constant() {
                return vec![MPoly::one(r.nvars)];
            }
            basis.push(r);
        }
    }
    let mut pairs = VecDeque::new();
    for j in 0..basis.len() {
        for i in 0..j {
            pairs.push_back((i, j));
        }
    }
    while let Some((i, j)) = pairs.pop_front() {
        let a = basis[i].exp(order.leading_index(&basis[i]));
        let b = basis[j].exp(order.leading_index(&basis[j]));
        // Buchberger's product criterion: coprime leading monomials.
        if a.iter().zip(b).all(|(&x, &y)| x == 0 || y == 0) {
            continue;
        }
        let s = s_polynomial(&basis[i], &basis[j], order, field);
        let r = reduce(&s, &basis, order, field).1;
        if r.is_zero() {
            continue;
        }
        let r = field.monic(&r, order);
        if r.is_constant() {
            return vec![MPoly::one(r.nvars)];
        }
        let k = basis.len();
        for i in 0..k {
            pairs.push_back((i, k));
        }
        basis.push(r);
    }
    // Remove redundant leading monomials, retaining the first equal one.
    let minimal: Vec<_> = basis
        .iter()
        .enumerate()
        .filter(|(i, f)| {
            let a = f.exp(order.leading_index(f));
            !basis.iter().enumerate().any(|(j, g)| {
                let b = g.exp(order.leading_index(g));
                *i != j && divides(b, a) && (a != b || j < *i)
            })
        })
        .map(|(_, f)| f.clone())
        .collect();
    let mut reduced = Vec::new();
    for (i, f) in minimal.iter().enumerate() {
        let others: Vec<_> = minimal
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, g)| g.clone())
            .collect();
        reduced.push(field.monic(&reduce(f, &others, order, field).1, order));
    }
    reduced.sort_by(|a, b| {
        order.compare(b.exp(order.leading_index(b)), a.exp(order.leading_index(a)))
    });
    reduced
}

/// Standard monomials are finite iff the leading ideal contains a pure
/// positive power of every generator. The unit ideal follows SymPy's
/// convention and is not classified as zero-dimensional.
pub fn is_zero_dimensional(basis: &[MPoly<Rat>], order: MonomialOrder, nvars: usize) -> bool {
    if basis.iter().any(MPoly::is_constant) {
        return false;
    }
    (0..nvars).all(|v| {
        basis.iter().any(|p| {
            let e = p.exp(order.leading_index(p));
            e[v] > 0 && e.iter().enumerate().all(|(j, &k)| j == v || k == 0)
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    use crate::poly::convert::from_exprs;

    #[test]
    fn buchberger_certificates() {
        let es: Vec<_> = ["x*y-1", "y**2-x"]
            .iter()
            .map(|s| parse_str(s).unwrap())
            .collect();
        let input = from_exprs(&es).unwrap().polys;
        for order in [
            MonomialOrder::Lex,
            MonomialOrder::Grlex,
            MonomialOrder::Grevlex,
        ] {
            for modulus in [None, Some(5)] {
                let field = CoefficientField { modulus };
                let g = buchberger(&input, order, field);
                assert!(is_zero_dimensional(&g, order, 2));
                for f in &input {
                    assert!(reduce(f, &g, order, field).1.is_zero());
                }
                for j in 0..g.len() {
                    for i in 0..j {
                        assert!(reduce(
                            &s_polynomial(&g[i], &g[j], order, field),
                            &g,
                            order,
                            field
                        )
                        .1
                        .is_zero());
                    }
                }
                let (q, r) = reduce(&input[0], &g, order, field);
                let reconstructed = q.iter().zip(&g).fold(r, |a, (q, g)| a.add(&q.mul(g)));
                assert_eq!(field.normalize(&reconstructed), field.normalize(&input[0]));
            }
        }
    }

    #[test]
    fn graded_orders_and_extreme_ideals() {
        assert_eq!(
            MonomialOrder::Lex.compare(&[1, 0], &[0, 3]),
            Ordering::Greater
        );
        assert_eq!(
            MonomialOrder::Grlex.compare(&[1, 0], &[0, 3]),
            Ordering::Less
        );
        assert_eq!(
            MonomialOrder::Grevlex.compare(&[1, 0, 1], &[0, 2, 0]),
            Ordering::Less
        );
        let field = CoefficientField { modulus: None };
        assert!(buchberger(&[MPoly::zero(2)], MonomialOrder::Lex, field).is_empty());
        assert_eq!(
            buchberger(
                &[MPoly::constant(Rat::from_i64(7), 2)],
                MonomialOrder::Lex,
                field
            ),
            vec![MPoly::one(2)]
        );
    }
}
