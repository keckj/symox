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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{ddm, fraction_free, fflu};

    fn gf7(rows: &[&[i64]]) -> Vec<Vec<DomainCoeff>> {
        let k = ExactDomain::GF(7);
        rows.iter().map(|r| r.iter().map(|&v| DomainCoeff(k.convert(&Expr::int(v)).unwrap())).collect()).collect()
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
