//! `Coeff` view of exact native domain elements for the fraction-free kernels.
//!
//! Every operation calls exactly the element arithmetic that the Python
//! protocols of `symox._core.DomainElement` call for two elements of the same
//! domain (`__add__`/`__sub__`/`__mul__` -> `DomainElement::binary`,
//! `__neg__` -> `zero - a`, `Domain.exquo` -> `DomainElement::exquo`,
//! `Domain.is_one` -> equality with `Domain.one`). The fraction-free kernels
//! therefore perform the same exact operations in the same order as the
//! Python bodies, without Python objects in between.
//!
//! Arithmetic errors (representation limits, non-exact quotients) are not
//! hidden: the first one is recorded in the active scope, which the binding
//! checks after the kernel to decline before publishing anything, so the
//! complete Python algorithm then raises its own error unchanged.
use super::convert::PolyError;
use super::domains::{CompositeBase, DomainElement, ElementValue, ExactDomain};
use super::int::{Coeff, Rat};
use crate::Expr;
use std::cell::RefCell;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DomainCoeff(pub DomainElement);

impl std::hash::Hash for DomainCoeff {
    // Equality may identify different representations (embeddings); the
    // kernels never hash coefficients, a constant hash keeps Hash/Eq coherent.
    fn hash<H: std::hash::Hasher>(&self, _: &mut H) {}
}

struct Context {
    zero: DomainElement,
    one: DomainElement,
    error: Option<PolyError>,
}

thread_local! {
    static CONTEXT: RefCell<Vec<Context>> = const { RefCell::new(Vec::new()) };
}

/// Integral domains whose native elements carry their own arithmetic: GF(p),
/// Gaussian ZZ_I/QQ_I, algebraic fields QQ<a>, rational-ground polynomial
/// rings/fraction fields ZZ[x], QQ[x], QQ(x) (multivariate included) and the
/// composite rings/fields over Gaussian or algebraic grounds (ZZ_I[x],
/// QQ<a>[x]). ZZ/QQ keep their typed kernels; inexact and expression domains
/// are never admitted (their Python bodies use `quo` or symbolic simplification).
pub fn admitted(domain: &ExactDomain) -> bool {
    match domain {
        // Z/nZ with composite n has zero divisors: a quotient by a non-unit
        // can fail where the Python body would have raised earlier or later.
        ExactDomain::GF(p) => super::gf::is_prime(*p),
        ExactDomain::Algebraic(_) | ExactDomain::Gaussian { .. } => true,
        ExactDomain::Functions { .. } => true,
        ExactDomain::Composite { .. } => matches!(
            domain.composite_base(),
            Some(CompositeBase::Rational { .. } | CompositeBase::Algebraic { .. })
        ),
        _ => false,
    }
}

/// Element representation expected for an element of an admitted domain.
pub fn native_member(domain: &ExactDomain, element: &DomainElement) -> bool {
    if element.domain != *domain {
        return false;
    }
    match (&element.domain, &element.value) {
        (ExactDomain::GF(_), ElementValue::Finite(_)) => true,
        (ExactDomain::Algebraic(_) | ExactDomain::Gaussian { .. }, ElementValue::Algebraic(_)) => true,
        (ExactDomain::Functions { .. }, ElementValue::Function(_)) => true,
        (ExactDomain::Composite { .. }, ElementValue::AlgebraicPoly(_) | ElementValue::GaussianFraction(..) | ElementValue::Function(_)) => true,
        _ => false,
    }
}

/// Active domain for `DomainCoeff::zero/one`; nested scopes are restored.
pub struct DomainScope {
    depth: usize,
}

impl DomainScope {
    pub fn new(domain: &ExactDomain) -> Result<Self, PolyError> {
        let zero = domain.convert(&Expr::zero())?;
        let one = domain.convert(&Expr::one())?;
        let depth = CONTEXT.with(|c| {
            let mut c = c.borrow_mut();
            c.push(Context { zero, one, error: None });
            c.len()
        });
        Ok(Self { depth })
    }
    /// True once an operation of this scope failed; results are then void.
    pub fn failed(&self) -> bool {
        CONTEXT.with(|c| c.borrow().get(self.depth - 1).is_some_and(|c| c.error.is_some()))
    }
    pub fn take_error(&self) -> Option<PolyError> {
        CONTEXT.with(|c| c.borrow_mut().get_mut(self.depth - 1).and_then(|c| c.error.take()))
    }
}

impl Drop for DomainScope {
    fn drop(&mut self) {
        CONTEXT.with(|c| c.borrow_mut().truncate(self.depth - 1));
    }
}

fn constant(one: bool) -> DomainElement {
    CONTEXT.with(|c| {
        let c = c.borrow();
        let context = c.last().expect("DomainCoeff used outside a DomainScope");
        if one { context.one.clone() } else { context.zero.clone() }
    })
}

fn record(error: PolyError) {
    CONTEXT.with(|c| {
        if let Some(context) = c.borrow_mut().last_mut() {
            context.error.get_or_insert(error);
        }
    });
}

fn checked(result: Result<DomainElement, PolyError>) -> DomainCoeff {
    match result {
        Ok(value) => DomainCoeff(value),
        Err(error) => {
            record(error);
            DomainCoeff(constant(false))
        }
    }
}

impl Coeff for DomainCoeff {
    fn zero() -> Self {
        DomainCoeff(constant(false))
    }
    fn one() -> Self {
        DomainCoeff(constant(true))
    }
    fn from_i64(i: i64) -> Self {
        let domain = constant(false).domain;
        checked(domain.convert(&Expr::int(i)))
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn is_one(&self) -> bool {
        CONTEXT.with(|c| c.borrow().last().is_some_and(|context| self.0 == context.one))
    }
    fn add(&self, o: &Self) -> Self {
        checked(self.0.binary(&o.0, '+'))
    }
    fn sub(&self, o: &Self) -> Self {
        checked(self.0.binary(&o.0, '-'))
    }
    fn mul(&self, o: &Self) -> Self {
        checked(self.0.binary(&o.0, '*'))
    }
    fn neg(&self) -> Self {
        checked(constant(false).binary(&self.0, '-'))
    }
    fn is_negative(&self) -> bool {
        false
    }
    /// `Domain.exquo`. A failure is recorded (the binding then declines) and
    /// also reported as `None`, so no kernel can continue with a wrong value.
    fn div_exact(&self, o: &Self) -> Option<Self> {
        match self.0.exquo(&o.0) {
            Ok(value) => Some(DomainCoeff(value)),
            Err(error) => {
                record(error);
                None
            }
        }
    }
    /// Not a rational view: these kernels publish elements, never `to_rat`.
    fn to_rat(&self) -> Rat {
        Rat::zero()
    }
}

/// Principal ideal domains of the native Smith kernel beyond the typed ZZ:
/// QQ, GF(p) with p prime, QQ_I, the rational function fields and the
/// univariate polynomial rings `QQ[x]` (the only polynomial rings the
/// reference `is_PID` admits here). Algebraic fields QQ<a> are left to the
/// Python body: the reference `1 / m[0][0]` raises TypeError on an ANP.
pub fn smith_admitted(domain: &ExactDomain) -> bool {
    match domain {
        ExactDomain::QQ => true,
        ExactDomain::GF(p) => super::gf::is_prime(*p),
        ExactDomain::Gaussian { integral } => !integral,
        ExactDomain::Functions { ring: false, .. } => true,
        ExactDomain::Functions { integral: false, ring: true, parameters, .. } => parameters.len() == 1,
        _ => false,
    }
}

/// Element representation expected in a Smith-admitted domain.
pub fn smith_member(domain: &ExactDomain, element: &DomainElement) -> bool {
    match domain {
        ExactDomain::QQ => element.domain == *domain && matches!(element.value, ElementValue::Rational(_)),
        _ => native_member(domain, element),
    }
}

fn unsupported(message: &str) -> PolyError {
    PolyError::Unsupported(message.into())
}

/// The reference Smith operations, each one the native call behind the
/// corresponding `symox._core.Domain` method (`div` -> `div_rem`, `gcd`,
/// `exquo`; `gcdex` -> the univariate field `gcdex` of `QQ[x]` or the
/// element `gcdex`), with failures recorded in the active scope.
impl super::smith::Pid for DomainCoeff {
    fn zero() -> Self { <Self as Coeff>::zero() }
    fn one() -> Self { <Self as Coeff>::one() }
    fn is_zero(&self) -> bool { self.0.is_zero() }
    fn add(&self, o: &Self) -> Self { Coeff::add(self, o) }
    fn mul(&self, o: &Self) -> Self { Coeff::mul(self, o) }
    fn neg(&self) -> Self { Coeff::neg(self) }
    fn div(&self, o: &Self) -> (Self, Self) {
        match self.0.div_rem(&o.0) {
            Ok((q, r)) => (DomainCoeff(q), DomainCoeff(r)),
            Err(error) => {
                record(error);
                (<Self as Coeff>::zero(), <Self as Coeff>::zero())
            }
        }
    }
    fn gcdex(&self, o: &Self) -> (Self, Self, Self) {
        let zero = || <Self as Coeff>::zero();
        let result = if let ExactDomain::Functions { integral: false, ring: true, parameters, .. } = &self.0.domain {
            super::api::gcdex(&self.0.as_expr(), &o.0.as_expr(), parameters).and_then(|(s, t, g)| {
                let domain = &self.0.domain;
                Ok((domain.convert(&s)?, domain.convert(&t)?, domain.convert(&g)?))
            })
        } else {
            self.0.gcdex(&o.0)
        };
        match result {
            Ok((s, t, g)) => (DomainCoeff(s), DomainCoeff(t), DomainCoeff(g)),
            Err(error) => {
                record(error);
                (zero(), zero(), zero())
            }
        }
    }
    fn gcd(&self, o: &Self) -> Self { checked(self.0.gcd(&o.0)) }
    fn exquo(&self, o: &Self) -> Self { checked(self.0.exquo(&o.0)) }
    /// Over a field `c = 1/m[0][0]`; over `QQ[x]` the reference
    /// `canonical_unit` is `-1` exactly when the leading coefficient is
    /// negative (`is_negative` reads the sign of `Poly(...).LC()`).
    fn unit_factor(&self) -> Option<Self> {
        if self.0.domain.is_field() {
            let one = constant(true);
            let c = checked(one.binary(&self.0, '/'));
            return (c.0 != one).then_some(c);
        }
        if !matches!(self.0.domain, ExactDomain::Functions { integral: false, ring: true, .. }) {
            record(unsupported("canonical unit outside the admitted Smith domains"));
            return None;
        }
        match self.0.ring_polynomial() {
            Ok(p) => p.lc().is_negative().then(|| Coeff::neg(&<Self as Coeff>::one())),
            Err(error) => {
                record(error);
                None
            }
        }
    }
    fn index_is_zero(i: usize) -> bool {
        match constant(false).domain {
            ExactDomain::GF(p) => (i as u64) % p == 0,
            _ => i == 0,
        }
    }
    /// The QQ_I zero compares unequal to the int 0 (the reference
    /// `GaussianElement.__eq__` returns NotImplemented), so the reference
    /// reaches `1 / 0`. The rational function field zero equals 0 as the
    /// reference `FracElement`, and the normalisation is skipped.
    fn zero_pivot_raises() -> bool {
        matches!(constant(false).domain, ExactDomain::Gaussian { integral: false })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{ddm, fraction_free, fflu};

    fn gf7(rows: &[&[i64]]) -> Vec<Vec<DomainCoeff>> {
        let k = ExactDomain::GF(7);
        rows.iter().map(|r| r.iter().map(|&v| DomainCoeff(k.convert(&Expr::int(v)).unwrap())).collect()).collect()
    }

    fn gf(p: u64, rows: &[&[i64]]) -> Vec<Vec<DomainCoeff>> {
        let k = ExactDomain::GF(p);
        rows.iter().map(|r| r.iter().map(|&v| DomainCoeff(k.convert(&Expr::int(v)).unwrap())).collect()).collect()
    }

    // Values of sympy 1.14 `_smith_normal_decomp` (full=True).
    #[test]
    fn gf_smith_matches_reference() {
        let k = ExactDomain::GF(7);
        let scope = DomainScope::new(&k).unwrap();
        let (invs, st) = crate::poly::smith::decompose(gf7(&[&[2, 4], &[6, 1]]), 2, 2, true, &mut || Ok::<_, ()>(())).unwrap().unwrap();
        let (s, t) = st.unwrap();
        assert_eq!(invs, gf(7, &[&[1, 1]])[0]);
        assert_eq!(s, gf(7, &[&[4, 0], &[6, 5]]));
        assert_eq!(t, gf(7, &[&[1, 5], &[0, 1]]));
        assert!(!scope.failed());
        drop(scope);
        // Row index 2 compares equal to GF(2).zero: the column search runs.
        let k = ExactDomain::GF(2);
        let scope = DomainScope::new(&k).unwrap();
        let (invs, st) = crate::poly::smith::decompose(gf(2, &[&[0, 1], &[0, 0], &[1, 0]]), 3, 2, true, &mut || Ok::<_, ()>(())).unwrap().unwrap();
        let (s, t) = st.unwrap();
        assert_eq!(invs, gf(2, &[&[1, 1]])[0]);
        assert_eq!(s, gf(2, &[&[1, 0, 0], &[0, 0, 1], &[0, 1, 0]]));
        assert_eq!(t, gf(2, &[&[0, 1], &[1, 0]]));
        assert!(!scope.failed());
    }

    #[test]
    fn gf_bareiss_and_rref_den() {
        let k = ExactDomain::GF(7);
        let scope = DomainScope::new(&k).unwrap();
        let mut a = gf7(&[&[1, 2], &[3, 4]]);
        let det = ddm::idet(&mut a, &mut || Ok::<_, ()>(()), &mut |_| Ok(())).ok().unwrap().unwrap();
        // 1*4 - 2*3 = -2 = 5 mod 7
        assert_eq!(det.0, k.convert(&Expr::int(5)).unwrap());
        let mut b = gf7(&[&[0, 2, 1], &[3, 4, 1]]);
        let (den, pivots) = fraction_free::dense(&mut b, &mut || Ok::<_, ()>(()), &mut |_| Ok(())).ok().unwrap().unwrap();
        assert_eq!(pivots, vec![0, 1]);
        assert!(!den.is_zero());
        let mut c = gf7(&[&[2, 4], &[1, 2]]);
        let cols = fflu::eliminate(&mut c, 2, &mut || Ok::<_, ()>(()), &mut |_| Ok(())).ok().unwrap().unwrap();
        assert_eq!(cols, vec![0]);
        assert!(!scope.failed());
    }

    #[test]
    fn scopes_nest_and_record_errors() {
        let outer = DomainScope::new(&ExactDomain::GF(5)).unwrap();
        {
            let inner = DomainScope::new(&ExactDomain::Gaussian { integral: true }).unwrap();
            let two = DomainCoeff(ExactDomain::Gaussian { integral: true }.convert(&Expr::int(2)).unwrap());
            let three = DomainCoeff(ExactDomain::Gaussian { integral: true }.convert(&Expr::int(3)).unwrap());
            assert!(three.div_exact(&two).is_none());
            assert!(inner.failed());
            assert!(DomainCoeff::one().is_one());
        }
        assert!(!outer.failed());
        assert_eq!(DomainCoeff::one().0.domain, ExactDomain::GF(5));
        assert!(admitted(&ExactDomain::GF(5)) && !admitted(&ExactDomain::GF(6)));
        assert!(!admitted(&ExactDomain::ZZ) && !admitted(&ExactDomain::QQ));
    }
}
