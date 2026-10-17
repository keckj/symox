//! Scoped limits for the checked factorization entry points.
//!
//! The legacy factor engines return plain values deep into recursive Hensel,
//! recombination and finite-field code. While that code is being reused, a
//! *private*, typed cancellation payload unwinds only as far as `run`, which
//! immediately converts it into a normal PolyError. `resume_unwind` bypasses
//! the panic hook; no cancellation crosses FFI. Unexpected panics are rethrown.
//! This boundary is enabled only for panic=unwind; abort builds explicitly
//! reject checked factorization. No downstream production catch_unwind exists
//! on these Rust-only factorization paths (the existing catches are tests).
//! TLS state is restored by RAII on every exit, including nested calls.
use super::convert::{PResult, PolyError};
use super::int::{Coeff, Int};
use std::cell::Cell;

#[derive(Clone, Copy, Debug)]
pub struct Limits {
    pub work: u64,
    pub attempts: u64,
    pub degree: usize,
    pub terms: usize,
    pub coefficient_bits: u64,
    pub storage_bits: u64,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            work: 20_000_000,
            attempts: 4096,
            degree: 1024,
            terms: 65_536,
            coefficient_bits: 16_384,
            storage_bits: 64_000_000,
        }
    }
}
#[derive(Clone, Copy)]
struct State {
    limits: Limits,
    work: u64,
    attempts: u64,
}
thread_local! {static STATE:Cell<Option<State>>=const{Cell::new(None)};}
#[derive(Debug)]
struct Cancelled(&'static str);
struct Restore(Option<State>);
impl Drop for Restore {
    fn drop(&mut self) {
        STATE.with(|slot| slot.set(self.0));
    }
}
#[cold]
fn cancel(reason: &'static str) -> ! {
    std::panic::resume_unwind(Box::new(Cancelled(reason)))
}
#[inline]
pub(crate) fn active() -> bool {
    STATE.with(|slot| slot.get().is_some())
}
#[inline]
pub(crate) fn charge(work: u64) {
    let failed = STATE.with(|slot| {
        let Some(mut state) = slot.get() else {
            return false;
        };
        let Some(left) = state.work.checked_sub(work) else {
            return true;
        };
        state.work = left;
        slot.set(Some(state));
        false
    });
    if failed {
        cancel("work");
    }
}
pub(crate) fn event(_phase: &'static str) {
    let failed = STATE.with(|slot| {
        let Some(mut state) = slot.get() else {
            return false;
        };
        let Some(left) = state.attempts.checked_sub(1) else {
            return true;
        };
        state.attempts = left;
        slot.set(Some(state));
        false
    });
    if failed {
        cancel(_phase);
    }
    charge(1);
}
pub(crate) fn require(condition: bool, reason: &'static str) {
    if !condition && active() {
        cancel(reason);
    }
}
pub(crate) fn dense(length: usize) {
    if let Some(state) = STATE.with(Cell::get) {
        if length > state.limits.degree.saturating_add(1) {
            cancel("intermediate polynomial degree");
        }
        if (length as u64).saturating_mul(256) > state.limits.storage_bits {
            cancel("dense storage");
        }
    }
    charge(length as u64);
}
pub(crate) fn sparse(terms: usize, variables: usize) {
    if let Some(state) = STATE.with(Cell::get) {
        if terms > state.limits.terms {
            cancel("sparse terms");
        }
        if (terms as u64).saturating_mul(256 + 32 * variables as u64) > state.limits.storage_bits {
            cancel("sparse storage");
        }
    }
    charge(terms as u64);
}
pub(crate) fn bits(bits: u64) {
    if let Some(state) = STATE.with(Cell::get) {
        if bits > state.limits.coefficient_bits {
            cancel("coefficient size");
        }
    }
}
pub(crate) fn storage(bits: u64) {
    if let Some(state) = STATE.with(Cell::get) {
        if bits > state.limits.storage_bits {
            cancel("coefficient storage");
        }
    }
}
pub(crate) fn product<C: Coeff>(a: &[C], b: &[C], terms: usize, variables: usize, sums: usize) {
    if !active() {
        return;
    }
    let a = a.iter().max_by_key(|c| c.budget_weight());
    let b = b.iter().max_by_key(|c| c.budget_weight());
    if let (Some(a), Some(b)) = (a, b) {
        storage(
            a.budget_product_weight(b, sums)
                .saturating_add(32 * variables as u64)
                .saturating_mul(terms as u64),
        );
    }
}
pub(crate) fn int_pair(a: &Int, b: &Int, multiply: bool, subtract: bool) {
    if !active() {
        return;
    }
    let (x, y) = (a.bits(), b.bits());
    let upper = if multiply {
        if a.is_zero() || b.is_zero() {
            0
        } else {
            x.saturating_add(y)
        }
    } else {
        x.max(y) + u64::from(a.is_negative() == (b.is_negative() != subtract))
    };
    bits(upper);
    charge(if multiply {
        1 + x.saturating_mul(y) / 64
    } else {
        1 + x.max(y) / 64
    });
}
pub(crate) fn int_pow(a: &Int, exponent: u32) {
    if !active() {
        return;
    }
    if !a.is_zero() && !a.is_one() && !a.is_minus_one() {
        bits(a.bits().saturating_mul(exponent as u64));
    }
    charge(
        1 + a
            .bits()
            .saturating_mul((32 - exponent.leading_zeros()) as u64)
            / 64,
    );
}
pub(crate) fn fused(sum: &Int, a: &Int, b: &Int, subtract: bool) {
    if !active() {
        return;
    }
    let product = if a.is_zero() || b.is_zero() {
        0
    } else {
        a.bits().saturating_add(b.bits())
    };
    let negative = (a.is_negative() != b.is_negative()) != subtract;
    bits(
        sum.bits().max(product)
            + u64::from(!sum.is_zero() && product != 0 && sum.is_negative() == negative),
    );
    charge(1 + sum.bits().max(product) / 64);
}
/// A transferable allowance for native stages separated by Python coercion or
/// rendering. Only pure Rust closures run with cancellation active: a callback
/// or finalizer cannot carry the private unwind payload through PyO3.
pub struct Session {
    state: State,
}
impl Session {
    pub fn new(limits: Limits) -> Self {
        Self {
            state: State {
                limits,
                work: limits.work,
                attempts: limits.attempts,
            },
        }
    }
    pub fn run<T>(&mut self, operation: impl FnOnce() -> PResult<T>) -> PResult<T> {
        if active() {
            let result = run(self.state.limits, operation);
            self.state = STATE.with(Cell::get).expect("active factor budget");
            return result;
        }
        let restore = Restore(None);
        STATE.with(|slot| slot.set(Some(self.state)));
        let result = run(self.state.limits, operation);
        self.state = STATE.with(Cell::get).expect("session factor budget");
        drop(restore);
        result
    }
}
pub fn run<T>(limits: Limits, operation: impl FnOnce() -> PResult<T>) -> PResult<T> {
    #[cfg(not(panic = "unwind"))]
    {
        let _ = (limits, operation);
        Err(PolyError::Unsupported(
            "checked factorization requires Rust panic=unwind cancellation support".into(),
        ))
    }
    #[cfg(panic = "unwind")]
    {
        let previous = STATE.with(Cell::get);
        // Nested entry points share the outer allowance; they never reset it.
        if previous.is_none() {
            STATE.with(|slot| {
                slot.set(Some(State {
                    limits,
                    work: limits.work,
                    attempts: limits.attempts,
                }))
            });
        }
        let restore = Restore(if previous.is_none() {
            None
        } else {
            STATE.with(Cell::get)
        });
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(operation));
        // A nested scope must retain charges consumed under its parent.
        if previous.is_some() {
            std::mem::forget(restore);
        } else {
            drop(restore);
        }
        match result {
            Ok(result) => result,
            Err(payload) => match payload.downcast::<Cancelled>() {
                Ok(cancelled) => Err(PolyError::Unsupported(format!(
                    "factorization resource limit: {}",
                    cancelled.0
                ))),
                Err(payload) => std::panic::resume_unwind(payload),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn staged_session_shares_allowance_and_deactivates_between_native_calls() {
        let mut session = super::Session::new(super::Limits {
            work: 3,
            ..Default::default()
        });
        session
            .run(|| {
                super::charge(2);
                Ok(())
            })
            .unwrap();
        assert!(!super::active());
        // A Python callback between stages has no ambient private cancellation.
        super::charge(u64::MAX);
        assert!(session
            .run(|| {
                super::charge(2);
                Ok(())
            })
            .unwrap_err()
            .to_string()
            .contains("resource limit"));
        assert!(!super::active());
        let panic = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            session.run::<()>(|| panic!("unexpected"))
        }));
        assert!(panic.is_err());
        assert!(!super::active());
    }
    use super::*;
    use crate::poly::{Int, MPoly};
    #[test]
    fn limits_restore_nested_threads_and_rethrow_unexpected_panics() {
        assert!(!active());
        let result = run(
            Limits {
                work: 2,
                ..Default::default()
            },
            || {
                charge(1);
                assert!(run(Limits::default(), || {
                    charge(2);
                    Ok(())
                })
                .is_err());
                charge(1);
                Ok(())
            },
        );
        assert!(result.is_ok());
        assert!(!active());
        let panic = std::panic::catch_unwind(|| {
            run::<()>(Limits::default(), || {
                std::panic::resume_unwind(Box::new(77u32))
            })
        });
        assert_eq!(*panic.unwrap_err().downcast::<u32>().unwrap(), 77);
        assert!(!active());
        let workers: Vec<_> = (0..3)
            .map(|_| {
                std::thread::spawn(|| {
                    assert!(run(
                        Limits {
                            work: 0,
                            ..Default::default()
                        },
                        || {
                            Int::from(2).pow(10);
                            Ok(())
                        }
                    )
                    .is_err());
                    assert!(!active());
                    assert_eq!(Int::from(2).pow(10), Int::from(1024));
                })
            })
            .collect();
        for worker in workers {
            worker.join().unwrap();
        }
    }
    #[test]
    fn growth_and_dense_shape_stop_before_expensive_allocation() {
        let f = MPoly::monomial(&[1_000_000_000], Int::one());
        let error = run(Limits::default(), || {
            crate::poly::dense::from_mpoly(&f, 0);
            Ok(())
        })
        .unwrap_err();
        assert!(error.to_string().contains("intermediate polynomial degree"));
        assert!(run(
            Limits {
                coefficient_bits: 64,
                ..Default::default()
            },
            || {
                Int::from(2).pow(1_000_000);
                Ok(())
            }
        )
        .unwrap_err()
        .to_string()
        .contains("coefficient size"));
        let c = Int::from_big(num_bigint::BigInt::from(1u32) << 11999usize);
        let a = MPoly::from_terms(2, (0..141).map(|j| (vec![j, 0], c.clone())).collect());
        let b = MPoly::from_terms(2, (0..141).map(|j| (vec![0, j], c.clone())).collect());
        assert!(run(
            Limits {
                coefficient_bits: 100_000,
                ..Default::default()
            },
            || {
                a.mul(&b);
                Ok(())
            }
        )
        .unwrap_err()
        .to_string()
        .contains("coefficient storage"));
        assert!(!active());
    }
    #[test]
    fn finite_field_splitting_consumes_the_scoped_allowance() {
        let error = run(
            Limits {
                attempts: 0,
                ..Default::default()
            },
            || {
                crate::poly::gf::edf(&[2, 3, 1], 1, 5, &mut crate::poly::gf::Rng::new(7));
                Ok(())
            },
        )
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("finite-field splitting attempts"));
        assert!(!active());
    }
    #[test]
    fn fast_small_integer_accumulators_also_obey_the_coefficient_limit() {
        let value = Int::from_u64(u32::MAX as u64);
        let limits = Limits {
            coefficient_bits: 64,
            ..Default::default()
        };
        let error = run(limits, || {
            let mut sum = Int::from_i64(i64::MAX);
            sum.add_mul(&value, &value);
            Ok(())
        })
        .unwrap_err();
        assert!(error.to_string().contains("coefficient size"));
        let error = run(limits, || {
            crate::poly::dense::mul_int(
                &[value.clone(), value.clone()],
                &[value.clone(), value.clone()],
            );
            Ok(())
        })
        .unwrap_err();
        assert!(error.to_string().contains("coefficient size"));
        assert!(!active());
    }
}
