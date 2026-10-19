//! Packed exponent vectors for the division-free characteristic polynomial.
//!
//! Every intermediate monomial in n-by-n Berkowitz is a product of at most
//! n input monomials. Reserve enough bits per exponent for that bound, then
//! encode variable zero in the most significant field. Multiplication becomes
//! integer addition, and integer ordering is exactly lexicographic ordering.
//! If the bound does not fit u128, the caller uses ordinary MPoly arithmetic.

use rustc_hash::FxHashMap;

use super::ff::{berkowitz, Ring};
use super::int::Int;
use super::mpoly::{MPoly, Mono};

/// Private to this module: only the bounded, division-free Berkowitz kernel
/// may multiply these polynomials. General polynomial arithmetic uses MPoly.
#[derive(Clone, PartialEq)]
struct PackedPoly {
    // Decreasing keys, distinct monomials, no zero coefficients.
    terms: Vec<(u128, Int)>,
}

impl PackedPoly {
    fn merge(&self, other: &Self, negate: bool) -> Self {
        let mut terms = Vec::with_capacity(self.terms.len() + other.terms.len());
        let (mut i, mut j) = (0, 0);
        while i < self.terms.len() && j < other.terms.len() {
            let (a, b) = (&self.terms[i], &other.terms[j]);
            match a.0.cmp(&b.0) {
                std::cmp::Ordering::Greater => {
                    terms.push(a.clone());
                    i += 1;
                }
                std::cmp::Ordering::Less => {
                    terms.push((b.0, if negate { b.1.neg() } else { b.1.clone() }));
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    let c = if negate { a.1.sub(&b.1) } else { a.1.add(&b.1) };
                    if !c.is_zero() {
                        terms.push((a.0, c));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        terms.extend(self.terms[i..].iter().cloned());
        terms.extend(
            other.terms[j..]
                .iter()
                .map(|(m, c)| (*m, if negate { c.neg() } else { c.clone() })),
        );
        Self { terms }
    }

    fn mul_term(&self, monomial: u128, coefficient: &Int) -> Self {
        Self {
            terms: self
                .terms
                .iter()
                .map(|(m, c)| (m + monomial, c.mul(coefficient)))
                .collect(),
        }
    }

    fn products(pairs: &[(&Self, &Self)]) -> Self {
        let size: usize = pairs
            .iter()
            .map(|(a, b)| a.terms.len() * b.terms.len())
            .sum();
        if size <= 256 {
            let mut products = Vec::with_capacity(size);
            for (a, b) in pairs {
                for (ma, ca) in &a.terms {
                    for (mb, cb) in &b.terms {
                        products.push((ma + mb, ca.mul(cb)));
                    }
                }
            }
            products.sort_unstable_by(|a, b| b.0.cmp(&a.0));
            let mut terms: Vec<(u128, Int)> = Vec::with_capacity(products.len());
            for (m, c) in products {
                match terms.last_mut() {
                    Some(last) if last.0 == m => last.1 = last.1.add(&c),
                    _ => terms.push((m, c)),
                }
            }
            terms.retain(|(_, c)| !c.is_zero());
            return Self { terms };
        }
        let mut map: FxHashMap<u128, Int> =
            FxHashMap::with_capacity_and_hasher(size, Default::default());
        for (a, b) in pairs {
            for (ma, ca) in &a.terms {
                for (mb, cb) in &b.terms {
                    let product = ca.mul(cb);
                    match map.entry(ma + mb) {
                        std::collections::hash_map::Entry::Occupied(mut entry) => {
                            let sum = entry.get().add(&product);
                            entry.insert(sum);
                        }
                        std::collections::hash_map::Entry::Vacant(entry) => {
                            entry.insert(product);
                        }
                    }
                }
            }
        }
        let mut terms: Vec<_> = map.into_iter().filter(|(_, c)| !c.is_zero()).collect();
        terms.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        Self { terms }
    }
}

impl Ring for PackedPoly {
    fn zero_like(&self) -> Self {
        Self { terms: Vec::new() }
    }

    fn one_like(&self) -> Self {
        Self {
            terms: vec![(0, Int::ONE)],
        }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn is_one(&self) -> bool {
        self.terms.len() == 1 && self.terms[0].0 == 0 && self.terms[0].1.is_one()
    }

    fn add(&self, other: &Self) -> Self {
        self.merge(other, false)
    }

    fn sub(&self, other: &Self) -> Self {
        self.merge(other, true)
    }

    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return self.zero_like();
        }
        if self.terms.len() == 1 {
            return other.mul_term(self.terms[0].0, &self.terms[0].1);
        }
        if other.terms.len() == 1 {
            return self.mul_term(other.terms[0].0, &other.terms[0].1);
        }
        Self::products(&[(self, other)])
    }

    fn neg(&self) -> Self {
        Self {
            terms: self.terms.iter().map(|(m, c)| (*m, c.neg())).collect(),
        }
    }

    fn exquo(&self, _other: &Self) -> Self {
        unreachable!("packed polynomials are only used by division-free Berkowitz")
    }

    fn sum_products(pairs: &[(&Self, &Self)], zero: &Self) -> Self {
        match pairs {
            [] => zero.clone(),
            [(a, b)] => a.mul(b),
            _ => Self::products(pairs),
        }
    }
}

/// Return None unless all intermediate exponents can be packed without carry.
pub(super) fn charpoly(a: &[MPoly], n: usize) -> Option<Vec<MPoly>> {
    let nvars = a.first()?.nvars;
    let max_exponent = a
        .iter()
        .flat_map(|p| &p.terms)
        .flat_map(|(m, _)| m.iter())
        .copied()
        .max()
        .unwrap_or(0);
    let bound = max_exponent.checked_mul(u32::try_from(n).ok()?)?;
    let bits = (u32::BITS - bound.leading_zeros()).max(1);
    if nvars.checked_mul(bits as usize)? > u128::BITS as usize {
        return None;
    }
    let packed: Vec<_> = a
        .iter()
        .map(|p| PackedPoly {
            terms: p
                .terms
                .iter()
                .map(|(m, c)| {
                    let key = m
                        .iter()
                        .fold(0u128, |key, exponent| (key << bits) | u128::from(*exponent));
                    (key, c.clone())
                })
                .collect(),
        })
        .collect();
    let one = PackedPoly {
        terms: vec![(0, Int::ONE)],
    };
    let mask = (1u128 << bits) - 1;
    Some(
        berkowitz(&packed, n, &one)
            .into_iter()
            .map(|p| MPoly {
                nvars,
                terms: p
                    .terms
                    .into_iter()
                    .map(|(key, c)| {
                        let monomial: Mono = (0..nvars)
                            .map(|i| ((key >> ((nvars - 1 - i) * bits as usize)) & mask) as u32)
                            .collect();
                        (monomial, c)
                    })
                    .collect(),
            })
            .collect(),
    )
}
