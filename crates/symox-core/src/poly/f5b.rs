//! Signature-based F5B over exact fields.
//!
//! The signature order is position-over-term with reversed generator indices.
//! Reductions preserve signatures: an admissible reducer must have strictly
//! smaller multiplied signature. Pair rejection uses the comparable and
//! rewritten criteria of Sun--Wang, arXiv:1004.0084v4, Section 3.
//!
//! QQ/GF(p) Buchberger retains its existing implementation. This module also
//! provides an exact prime-field coefficient adapter for its generic F5B path.

use std::cell::Cell;
use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

use super::convert::{PResult, PolyError};
use super::fieldalg;
use super::groebner::MonomialOrder;
use super::int::{Coeff, Field, Rat};
use super::mpoly::MPoly;

#[derive(Clone, Debug, PartialEq, Eq)]
struct Signature {
    monomial: Vec<u32>,
    position: usize,
}
impl Signature {
    fn compare(&self, other: &Self, order: MonomialOrder) -> Ordering {
        other
            .position
            .cmp(&self.position)
            .then_with(|| order.compare(&self.monomial, &other.monomial))
    }
    fn multiplied(&self, monomial: &[u32]) -> Self {
        Self {
            monomial: self
                .monomial
                .iter()
                .zip(monomial)
                .map(|(a, b)| super::mpoly::exponent_add(*a,*b))
                .collect(),
            position: self.position,
        }
    }
}

#[derive(Clone, Debug)]
struct Label<C: Field> {
    signature: Signature,
    polynomial: MPoly<C>,
    number: usize,
    leading: Vec<u32>,
    witness: Option<Vec<MPoly<C>>>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
struct Part {
    label: usize,
    number: usize,
    signature: Signature,
    multiplier: Vec<u32>,
}
impl Part {
    fn compare(&self, other: &Self, order: MonomialOrder) -> Ordering {
        self.signature
            .compare(&other.signature, order)
            .then_with(|| other.number.cmp(&self.number))
            .then_with(|| self.label.cmp(&other.label))
            .then_with(|| self.multiplier.cmp(&other.multiplier))
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
struct Pair {
    first: Part,
    second: Part,
    serial: usize,
    order: MonomialOrder,
}
impl Ord for Pair {
    fn cmp(&self, other: &Self) -> Ordering {
        let rank = |o| match o {
            MonomialOrder::Lex => 0,
            MonomialOrder::Grlex => 1,
            MonomialOrder::Grevlex => 2,
        };
        rank(self.order)
            .cmp(&rank(other.order))
            .then_with(|| self.first.compare(&other.first, self.order))
            .then_with(|| self.second.compare(&other.second, self.order))
            .then_with(|| self.serial.cmp(&other.serial))
    }
}
impl PartialOrd for Pair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
fn divides(a: &[u32], b: &[u32]) -> bool {
    a.iter().zip(b).all(|(a, b)| a <= b)
}
fn difference(a: &[u32], b: &[u32]) -> Vec<u32> {
    a.iter().zip(b).map(|(a, b)| a - b).collect()
}
fn pair<C: Field>(
    i: usize,
    j: usize,
    labels: &[Label<C>],
    order: MonomialOrder,
    serial: usize,
) -> Pair {
    let a = &labels[i];
    let b = &labels[j];
    let lcm: Vec<_> = a
        .leading
        .iter()
        .zip(&b.leading)
        .map(|(a, b)| (*a).max(*b))
        .collect();
    let part = |label: usize, p: &Label<C>| {
        let multiplier = difference(&lcm, &p.leading);
        Part {
            label,
            number: p.number,
            signature: p.signature.multiplied(&multiplier),
            multiplier,
        }
    };
    let mut first = part(i, a);
    let mut second = part(j, b);
    if first.compare(&second, order) == Ordering::Less {
        std::mem::swap(&mut first, &mut second);
    }
    Pair {
        first,
        second,
        serial,
        order,
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct F5Stats {
    pub processed_pairs: usize,
    pub signature_reductions: usize,
    pub comparable_rejections: usize,
    pub rewritten_rejections: usize,
    pub zero_reductions: usize,
}
#[derive(Clone, Copy, Debug)]
enum Rejection {
    Comparable,
    Rewritten,
}
fn reject<C: Field>(part: &Part, labels: &[Label<C>]) -> Option<Rejection> {
    for h in labels {
        if part.signature.position < h.signature.position
            && !h.polynomial.is_zero()
            && divides(&h.leading, &part.signature.monomial)
        {
            return Some(Rejection::Comparable);
        }
        if part.signature.position == h.signature.position
            && part.number < h.number
            && divides(&h.signature.monomial, &part.signature.monomial)
        {
            return Some(Rejection::Rewritten);
        }
    }
    None
}
fn rejected<C: Field>(pair: &Pair, labels: &[Label<C>], stats: &mut F5Stats) -> bool {
    match reject(&pair.first, labels).or_else(|| reject(&pair.second, labels)) {
        Some(Rejection::Comparable) => {
            stats.comparable_rejections += 1;
            true
        }
        Some(Rejection::Rewritten) => {
            stats.rewritten_rejections += 1;
            true
        }
        None => false,
    }
}
fn subtract_witness<C: Field>(
    target: &mut Option<Vec<MPoly<C>>>,
    source: &Option<Vec<MPoly<C>>>,
    monomial: &[u32],
    coefficient: &C,
) {
    if let (Some(target), Some(source)) = (target, source) {
        for (a, b) in target.iter_mut().zip(source) {
            *a = a.sub(&b.mul_term(monomial, coefficient));
        }
    }
}
fn normalize<C: Field>(label: &mut Label<C>, order: MonomialOrder) {
    if label.polynomial.is_zero() {
        label.leading.clear();
        return;
    }
    let i = order.leading_index(&label.polynomial);
    let inverse = label.polynomial.coeff(i).inv();
    label.leading = label.polynomial.exp(i).to_vec();
    label.polynomial = label.polynomial.scale(&inverse);
    if let Some(witness) = &mut label.witness {
        for p in witness {
            *p = p.scale(&inverse);
        }
    }
}
fn signature_reduce<C: Field>(
    mut p: Label<C>,
    labels: &[Label<C>],
    reducers: &[usize],
    order: MonomialOrder,
    stats: &mut F5Stats,
) -> Label<C> {
    while !p.polynomial.is_zero() {
        let i = order.leading_index(&p.polynomial);
        let exponent = p.polynomial.exp(i).to_vec();
        let Some((id, multiplier)) = reducers.iter().find_map(|&id| {
            let h = &labels[id];
            if !divides(&h.leading, &exponent) {
                return None;
            }
            let multiplier = difference(&exponent, &h.leading);
            (h.signature
                .multiplied(&multiplier)
                .compare(&p.signature, order)
                == Ordering::Less)
                .then_some((id, multiplier))
        }) else {
            break;
        };
        let h = &labels[id];
        let j = order.leading_index(&h.polynomial);
        let coefficient = p.polynomial.coeff(i).div(h.polynomial.coeff(j));
        p.polynomial = p
            .polynomial
            .sub(&h.polynomial.mul_term(&multiplier, &coefficient));
        subtract_witness(&mut p.witness, &h.witness, &multiplier, &coefficient);
        stats.signature_reductions += 1;
    }
    normalize(&mut p, order);
    p
}

struct ResultWithTrace<C: Field> {
    basis: Vec<MPoly<C>>,
    witnesses: Option<Vec<Vec<MPoly<C>>>>,
    stats: F5Stats,
}

/// Compute F5B without retaining polynomial transformation matrices.
pub fn f5b<C: Field>(input: &[MPoly<C>], order: MonomialOrder) -> PResult<Vec<MPoly<C>>> {super::mpoly::checked_operation(||{
    Ok(run(input, order, false)?.basis)
})}

/// An independently checkable equality-of-ideals and Groebner certificate.
#[derive(Clone, Debug)]
pub struct CertifiedBasis<C: Field> {
    pub basis: Vec<MPoly<C>>,
    /// `basis[i] = sum(witnesses[i][j] * input[j])`.
    pub witnesses: Vec<Vec<MPoly<C>>>,
    pub stats: F5Stats,
}
pub fn f5b_certified<C: Field>(
    input: &[MPoly<C>],
    order: MonomialOrder,
) -> PResult<CertifiedBasis<C>> {super::mpoly::checked_operation(||{
    let result = run(input, order, true)?;
    Ok(CertifiedBasis {
        basis: result.basis,
        witnesses: result.witnesses.unwrap(),
        stats: result.stats,
    })
})}
impl<C: Field> CertifiedBasis<C> {
    pub fn verify(&self, input: &[MPoly<C>], order: MonomialOrder) -> bool {
        if self.witnesses.len() != self.basis.len() {
            return false;
        }
        let Some(first) = input.first().or_else(|| self.basis.first()) else {
            return self.basis.is_empty();
        };
        let nvars = first.nvars;
        if input.iter().chain(&self.basis).any(|p| p.nvars != nvars) {
            return false;
        }
        for (p, witness) in self.basis.iter().zip(&self.witnesses) {
            if witness.len() != input.len()
                || witness.iter().any(|w| w.nvars != nvars)
                || p.is_zero()
            {
                return false;
            }
            let reconstructed = witness
                .iter()
                .zip(input)
                .fold(MPoly::zero(nvars), |a, (w, f)| a.add(&w.mul(f)));
            if &reconstructed != p || !p.coeff(order.leading_index(p)).is_one() {
                return false;
            }
        }
        if input
            .iter()
            .any(|f| !fieldalg::reduce(f, &self.basis, order).1.is_zero())
        {
            return false;
        }
        for (i, f) in self.basis.iter().enumerate() {
            for (j, g) in self.basis.iter().enumerate() {
                if i != j
                    && f.iter()
                        .any(|(m, _)| divides(g.exp(order.leading_index(g)), m))
                {
                    return false;
                }
                if j <= i {
                    continue;
                }
                let a = order.leading_index(f);
                let b = order.leading_index(g);
                let lcm: Vec<_> = f
                    .exp(a)
                    .iter()
                    .zip(g.exp(b))
                    .map(|(a, b)| (*a).max(*b))
                    .collect();
                let s = f
                    .mul_term(&difference(&lcm, f.exp(a)), &f.coeff(a).inv())
                    .sub(&g.mul_term(&difference(&lcm, g.exp(b)), &g.coeff(b).inv()));
                if !fieldalg::reduce(&s, &self.basis, order).1.is_zero() {
                    return false;
                }
            }
        }
        true
    }
}

fn run<C: Field>(
    input: &[MPoly<C>],
    order: MonomialOrder,
    trace: bool,
) -> PResult<ResultWithTrace<C>> {
    let mut labels: Vec<Label<C>> = Vec::new();
    let mut reducers = Vec::new();
    let mut stats = F5Stats::default();
    let nvars = input.first().map_or(0, |p| p.nvars);
    if input.iter().any(|p| p.nvars != nvars) {
        return Err(PolyError::Failed(
            "incompatible polynomial generator counts".into(),
        ));
    }
    for (position, f) in input.iter().enumerate() {
        if f.is_zero() {
            continue;
        }
        let witness = trace.then(|| {
            (0..input.len())
                .map(|j| {
                    if j == position {
                        MPoly::one(nvars)
                    } else {
                        MPoly::zero(nvars)
                    }
                })
                .collect()
        });
        let mut label = Label {
            signature: Signature {
                monomial: vec![0; nvars],
                position,
            },
            polynomial: f.clone(),
            number: labels.len(),
            leading: Vec::new(),
            witness,
        };
        normalize(&mut label, order);
        reducers.push(labels.len());
        labels.push(label);
    }
    let mut queue = BinaryHeap::new();
    let mut serial = 0;
    for j in 0..labels.len() {
        for i in 0..j {
            queue.push(Reverse(pair(i, j, &labels, order, serial)));
            serial += 1;
        }
    }
    sort_reducers(&mut reducers, &labels, order);
    while let Some(Reverse(cp)) = queue.pop() {
        if rejected(&cp, &labels, &mut stats) {
            continue;
        }
        stats.processed_pairs += 1;
        let a = &labels[cp.first.label];
        let b = &labels[cp.second.label];
        let ca = a.polynomial.coeff(order.leading_index(&a.polynomial)).inv();
        let cb = b.polynomial.coeff(order.leading_index(&b.polynomial)).inv();
        let polynomial = a
            .polynomial
            .mul_term(&cp.first.multiplier, &ca)
            .sub(&b.polynomial.mul_term(&cp.second.multiplier, &cb));
        let mut witness = a.witness.as_ref().map(|w| {
            w.iter()
                .map(|p| p.mul_term(&cp.first.multiplier, &ca))
                .collect()
        });
        subtract_witness(&mut witness, &b.witness, &cp.second.multiplier, &cb);
        let p = Label {
            signature: cp.first.signature,
            polynomial,
            number: labels.len(),
            leading: Vec::new(),
            witness,
        };
        let p = signature_reduce(p, &labels, &reducers, order, &mut stats);
        let zero = p.polynomial.is_zero();
        let k = labels.len();
        labels.push(p);
        if zero {
            stats.zero_reductions += 1;
        } else {
            for &i in &reducers {
                let next = pair(i, k, &labels, order, serial);
                serial += 1;
                if !rejected(&next, &labels, &mut stats) {
                    queue.push(Reverse(next));
                }
            }
            reducers.push(k);
            sort_reducers(&mut reducers, &labels, order);
        }
        // Recheck all queued pairs against the new rewrite rule. Stable
        // creation numbers and append-only labels remain valid references.
        queue.retain(|Reverse(cp)| !rejected(cp, std::slice::from_ref(&labels[k]), &mut stats));
    }
    let minimal: Vec<_> = reducers
        .iter()
        .copied()
        .filter(|&i| {
            !reducers.iter().any(|&j| {
                j != i
                    && divides(&labels[j].leading, &labels[i].leading)
                    && (labels[j].leading != labels[i].leading || j < i)
            })
        })
        .collect();
    let mut output = Vec::new();
    for &i in &minimal {
        let mut label = labels[i].clone();
        let others: Vec<_> = minimal.iter().copied().filter(|&j| j != i).collect();
        let polynomials: Vec<_> = others
            .iter()
            .map(|&j| labels[j].polynomial.clone())
            .collect();
        let (q, r) = fieldalg::reduce(&label.polynomial, &polynomials, order);
        label.polynomial = r;
        if let Some(witness) = &mut label.witness {
            for (quotient, &j) in q.iter().zip(&others) {
                for (w, g) in witness.iter_mut().zip(labels[j].witness.as_ref().unwrap()) {
                    *w = w.sub(&quotient.mul(g));
                }
            }
        }
        normalize(&mut label, order);
        if !label.polynomial.is_zero() {
            output.push(label);
        }
    }
    output.sort_by(|a, b| order.compare(&b.leading, &a.leading));
    let witnesses = trace.then(|| output.iter().map(|p| p.witness.clone().unwrap()).collect());
    Ok(ResultWithTrace {
        basis: output.into_iter().map(|p| p.polynomial).collect(),
        witnesses,
        stats,
    })
}
fn sort_reducers<C: Field>(reducers: &mut [usize], labels: &[Label<C>], order: MonomialOrder) {
    reducers.sort_by(|&i, &j| {
        order
            .compare(&labels[j].leading, &labels[i].leading)
            .then_with(|| i.cmp(&j))
    });
}

thread_local! {static PRIME:Cell<Option<u64>>=const {Cell::new(None)};}
/// Factory context only: existing coefficients also retain their modulus.
pub struct PrimeScope(Option<u64>);
impl PrimeScope {
    pub fn new(prime: u64) -> PResult<Self> {
        if prime < 2 || prime >= 1u64 << 32 || !super::gf::is_prime(prime) {
            return Err(PolyError::Failed(
                "a prime modulus below 2**32 is required".into(),
            ));
        }
        Ok(Self(PRIME.with(|c| c.replace(Some(prime)))))
    }
}
impl Drop for PrimeScope {
    fn drop(&mut self) {
        PRIME.with(|c| c.set(self.0));
    }
}
fn prime() -> u64 {
    PRIME
        .with(Cell::get)
        .expect("prime coefficient factory scope")
}
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Fp {
    value: u64,
    prime: u64,
}
impl Fp {
    pub(crate) fn modulus(&self)->u64{self.prime}
    pub fn from_rat(c: &Rat, prime: u64) -> PResult<Self> {
        if prime < 2 || prime >= 1u64 << 32 || !super::gf::is_prime(prime) {
            return Err(PolyError::Failed(
                "a prime modulus below 2**32 is required".into(),
            ));
        }
        Self::from_validated_rat(c, prime)
    }
    fn from_validated_rat(c: &Rat, prime: u64) -> PResult<Self> {
        let denominator = c.d.rem_u64(prime);
        if denominator == 0 {
            return Err(PolyError::Failed(
                "coefficient denominator is zero in the finite field".into(),
            ));
        }
        Ok(Self {
            value: super::gf::mulmod(
                c.n.rem_u64(prime),
                super::gf::invmod(denominator, prime),
                prime,
            ),
            prime,
        })
    }
    fn same(&self, other: &Self) {
        assert_eq!(
            self.prime, other.prime,
            "incompatible prime coefficient fields"
        );
    }
}
impl Coeff for Fp {
    fn zero() -> Self {
        Self {
            value: 0,
            prime: prime(),
        }
    }
    fn one() -> Self {
        Self {
            value: 1,
            prime: prime(),
        }
    }
    fn from_i64(i: i64) -> Self {
        let prime = prime();
        Self {
            value: i.rem_euclid(prime as i64) as u64,
            prime,
        }
    }
    fn is_zero(&self) -> bool {
        self.value == 0
    }
    fn is_one(&self) -> bool {
        self.value == 1
    }
    fn add(&self, o: &Self) -> Self {
        self.same(o);
        Self {
            value: (self.value + o.value) % self.prime,
            prime: self.prime,
        }
    }
    fn sub(&self, o: &Self) -> Self {
        self.same(o);
        Self {
            value: (self.value + self.prime - o.value) % self.prime,
            prime: self.prime,
        }
    }
    fn mul(&self, o: &Self) -> Self {
        self.same(o);
        Self {
            value: super::gf::mulmod(self.value, o.value, self.prime),
            prime: self.prime,
        }
    }
    fn neg(&self) -> Self {
        Self {
            value: if self.value == 0 {
                0
            } else {
                self.prime - self.value
            },
            prime: self.prime,
        }
    }
    fn is_negative(&self) -> bool {
        self.value > self.prime / 2
    }
    fn div_exact(&self, o: &Self) -> Option<Self> {
        if o.is_zero() {
            None
        } else {
            Some(self.mul(&o.inv()))
        }
    }
    fn to_rat(&self) -> Rat {
        Rat::from_i64(if self.is_negative() {
            self.value as i64 - self.prime as i64
        } else {
            self.value as i64
        })
    }
}
impl Field for Fp {
    fn inv(&self) -> Self {
        assert!(!self.is_zero(), "division by zero finite-field coefficient");
        Self {
            value: super::gf::invmod(self.value, self.prime),
            prime: self.prime,
        }
    }
}
pub fn f5b_modular(
    input: &[MPoly<Rat>],
    order: MonomialOrder,
    prime: u64,
) -> PResult<Vec<MPoly<Rat>>> {super::mpoly::checked_operation(||{
    let _scope = PrimeScope::new(prime)?;
    let converted = input
        .iter()
        .map(|p| {
            Ok(MPoly::from_terms(
                p.nvars,
                p.iter()
                    .map(|(m, c)| Fp::from_validated_rat(c, prime).map(|c| (m.to_vec(), c)))
                    .collect::<PResult<_>>()?,
            ))
        })
        .collect::<PResult<Vec<_>>>()?;
    Ok(f5b(&converted, order)?
        .into_iter()
        .map(|p| {
            MPoly::from_terms(
                p.nvars,
                p.iter().map(|(m, c)| (m.to_vec(), c.to_rat())).collect(),
            )
        })
        .collect())
})}

#[cfg(test)]
mod tests {
    use super::super::convert::from_exprs_gens;
    use super::*;
    use crate::expr::Expr;
    const ORDERS: [MonomialOrder; 3] = [
        MonomialOrder::Lex,
        MonomialOrder::Grlex,
        MonomialOrder::Grevlex,
    ];
    fn expressions(strings: &[&str]) -> Vec<Expr> {
        strings
            .iter()
            .map(|s| crate::parser::parse_str(s).unwrap())
            .collect()
    }
    fn rational(strings: &[&str], gens: &[Expr]) -> Vec<MPoly<Rat>> {
        from_exprs_gens(&expressions(strings), gens).unwrap().polys
    }
    #[test]
    fn rational_certificates_reference_bases_and_signature_criteria() {
        let gens = expressions(&["x", "y", "z"]);
        let mut stats = F5Stats::default();
        for (system, n) in [
            (vec!["x^2 + 2*x*y^2", "x*y + 2*y^3 - 1"], 2),
            (vec!["x+2*y+2*z-1", "x^2+2*y^2+2*z^2-x", "2*x*y+2*y*z-y"], 3),
            (vec!["x^2", "x*y", "y^2"], 2),
            (vec!["x^2-y", "x^3-z"], 3),
            (vec!["x*y-1", "y*z-1", "z*x-1"], 3),
            (vec!["x^2-y", "y^2-z", "z^2-2"], 3),
            (vec!["0", "x+y", "x+y", "2*x+2*y"], 2),
            (vec!["x*y", "x*y+1"], 2),
            (vec!["0", "0"], 2),
            (vec!["2", "x^2-y"], 2),
            (vec!["x^2/2 + y/3", "x*y-1/5"], 2),
        ] {
            let input = rational(&system, &gens[..n]);
            for order in ORDERS {
                let result = f5b_certified(&input, order).unwrap();
                assert!(
                    result.verify(&input, order),
                    "invalid certificate {system:?} {order:?}"
                );
                assert_eq!(
                    result.basis,
                    super::super::groebner::buchberger(
                        &input,
                        order,
                        super::super::groebner::CoefficientField { modulus: None }
                    ),
                    "basis {system:?} {order:?}"
                );
                stats.signature_reductions += result.stats.signature_reductions;
                stats.comparable_rejections += result.stats.comparable_rejections;
                stats.rewritten_rejections += result.stats.rewritten_rejections;
            }
        }
        assert!(
            stats.signature_reductions > 0,
            "must perform actual signature reductions"
        );
        assert!(
            stats.comparable_rejections > 0,
            "must use comparable criterion"
        );
        assert!(
            stats.rewritten_rejections > 0,
            "must use rewritten criterion"
        );
        assert!(f5b_certified::<Rat>(&[], MonomialOrder::Lex)
            .unwrap()
            .verify(&[], MonomialOrder::Lex));
    }
    #[test]
    fn reject_certificate_that_only_proves_input_membership() {
        let gens = expressions(&["x", "y"]);
        let input = rational(&["x^2", "y^2"], &gens);
        let fake = CertifiedBasis {
            basis: vec![MPoly::<Rat>::one(2)],
            witnesses: vec![vec![MPoly::zero(2); 2]],
            stats: F5Stats::default(),
        };
        // [1] passes all S-pairs and makes every input reduce to zero, but
        // it generates a strictly larger ideal. Provenance rejects it.
        assert!(!fake.verify(&input, MonomialOrder::Lex));
        let mut good = f5b_certified(&input, MonomialOrder::Lex).unwrap();
        good.witnesses[0][0] = MPoly::zero(2);
        assert!(!good.verify(&input, MonomialOrder::Lex));
    }
    #[test]
    fn deterministic_quadratic_systems_have_exact_provenance() {
        let mut seed = 918273u64;
        for trial in 0..16 {
            let mut input = Vec::new();
            for _ in 0..2 {
                let mut terms = Vec::new();
                for exponent in [
                    vec![2, 0],
                    vec![1, 1],
                    vec![0, 2],
                    vec![1, 0],
                    vec![0, 1],
                    vec![0, 0],
                ] {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    terms.push((exponent, Rat::from_i64(((seed >> 32) % 5) as i64 - 2)));
                }
                input.push(MPoly::from_terms(2, terms));
            }
            for order in ORDERS {
                let result = f5b_certified(&input, order).unwrap();
                assert!(result.verify(&input, order), "trial {trial} {order:?}");
                assert_eq!(
                    result.basis,
                    super::super::groebner::buchberger(
                        &input,
                        order,
                        super::super::groebner::CoefficientField { modulus: None }
                    )
                );
            }
        }
    }
    #[test]
    fn prime_fields_exact_certificates_and_rational_denominators() {
        let gens = expressions(&["x", "y"]);
        let original = rational(&["2*x^2+y-1", "x*y+3*x+2", "y^2+x"], &gens);
        for prime in [2, 3, 5, 101, 4294967291] {
            let _scope = PrimeScope::new(prime).unwrap();
            let input: Vec<_> = original
                .iter()
                .map(|p| {
                    MPoly::from_terms(
                        2,
                        p.iter()
                            .map(|(m, c)| (m.to_vec(), Fp::from_rat(c, prime).unwrap()))
                            .collect(),
                    )
                })
                .collect();
            for order in ORDERS {
                let result = f5b_certified(&input, order).unwrap();
                assert!(result.verify(&input, order));
                assert_eq!(
                    f5b_modular(&original, order, prime).unwrap(),
                    super::super::groebner::buchberger(
                        &original,
                        order,
                        super::super::groebner::CoefficientField {
                            modulus: Some(prime)
                        }
                    )
                );
            }
        }
        assert!(f5b_modular(&rational(&["x/5+y"], &gens), MonomialOrder::Lex, 5).is_err());
        for prime in [0, 1, 4, 9, 1u64 << 32] {
            assert!(f5b_modular(&original, MonomialOrder::Lex, prime).is_err());
        }
        let input = rational(&["x^2/2+y/3", "x*y+1/2"], &gens);
        assert_eq!(
            f5b_modular(&input, MonomialOrder::Lex, 5).unwrap(),
            super::super::groebner::buchberger(
                &input,
                MonomialOrder::Lex,
                super::super::groebner::CoefficientField { modulus: Some(5) }
            )
        );
    }
    #[test]
    fn fraction_and_algebraic_coefficients_certify_ideal_equality() {
        let gens = expressions(&["x", "y"]);
        let parameters = expressions(&["a", "b"]);
        let _scope = super::super::ratfunc::Scope::new(parameters.len());
        let input: Vec<_> = expressions(&["a*x*y-1", "x^2-y/(b+1)"])
            .iter()
            .map(|e| {
                super::super::fraction::FractionPoly::from_expr(e, &gens, &parameters, "QQ(a,b)")
                    .unwrap()
                    .rep
            })
            .collect();
        for order in ORDERS {
            let result = f5b_certified(&input, order).unwrap();
            assert!(result.verify(&input, order));
            assert_eq!(result.basis, fieldalg::buchberger(&input, order));
        }
        for generator in [crate::ops::sqrt(&Expr::int(2)), Expr::i()] {
            let field = super::super::algebraic::NumberField::quadratic(&generator).unwrap();
            let _scope = super::super::algebraic::AlgebraicScope::new(field.clone());
            let system = [
                crate::ops::sub(
                    &crate::ops::mul(&[generator, gens[0].clone(), gens[1].clone()]),
                    &Expr::one(),
                ),
                crate::ops::sub(&gens[0].powi(2), &gens[1]),
            ];
            let input: Vec<_> = system
                .iter()
                .map(|e| {
                    super::super::algebraic::AlgebraicPoly::from_expr(e, &gens, field.clone())
                        .unwrap()
                        .rep
                })
                .collect();
            for order in ORDERS {
                let result = f5b_certified(&input, order).unwrap();
                assert!(result.verify(&input, order));
                assert_eq!(result.basis, fieldalg::buchberger(&input, order));
            }
        }
    }
    #[test]
    fn signature_reduction_is_strict_and_not_ordinary_remainder() {
        let polynomial = MPoly::<Rat>::var(0, 1);
        let make = |position, number| Label {
            signature: Signature {
                monomial: vec![0],
                position,
            },
            polynomial: polynomial.clone(),
            number,
            leading: vec![1],
            witness: None,
        };
        let reducer = make(2, 0);
        let mut stats = F5Stats::default();
        let blocked = signature_reduce(
            make(4, 1),
            std::slice::from_ref(&reducer),
            &[0],
            MonomialOrder::Lex,
            &mut stats,
        );
        assert_eq!(blocked.polynomial, polynomial);
        let equal = signature_reduce(
            make(2, 1),
            std::slice::from_ref(&reducer),
            &[0],
            MonomialOrder::Lex,
            &mut stats,
        );
        assert_eq!(equal.polynomial, polynomial);
        let admissible =
            signature_reduce(make(1, 1), &[reducer], &[0], MonomialOrder::Lex, &mut stats);
        assert!(admissible.polynomial.is_zero());
        assert_eq!(stats.signature_reductions, 1);
        let signature = Signature {
            monomial: vec![1],
            position: 2,
        };
        assert_eq!(
            signature.compare(&signature, MonomialOrder::Lex),
            Ordering::Equal
        );
    }
    #[test]
    fn prime_factory_scope_restores_on_nesting_threads_and_panic() {
        let coefficient;
        {
            let _scope = PrimeScope::new(5).unwrap();
            coefficient = Fp::from_i64(3);
            {
                let _nested = PrimeScope::new(7).unwrap();
                assert_eq!(Fp::from_i64(6).prime, 7);
            }
            assert_eq!(Fp::one().prime, 5);
            let worker = std::thread::spawn(|| {
                let _scope = PrimeScope::new(101).unwrap();
                Fp::from_i64(100)
            });
            assert_eq!(worker.join().unwrap().prime, 101);
            assert_eq!(Fp::one().prime, 5);
            assert!(std::panic::catch_unwind(|| {
                let _scope = PrimeScope::new(7).unwrap();
                panic!("restore prime context");
            })
            .is_err());
            assert_eq!(Fp::one().prime, 5);
        }
        // Arithmetic on previously constructed coefficients is independent
        // from the factory context and remains valid after scope exit.
        assert!(coefficient.mul(&coefficient.inv()).is_one());
    }
}
