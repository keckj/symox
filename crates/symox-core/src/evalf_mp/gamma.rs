//! Argument-independent Spouge coefficients, at their full working precision.
//!
//! The coefficients are reused only for the exact `(precision, a)` pair that
//! produced them. In particular, a larger precision is never rounded down to
//! answer another key: that would change the existing sequence of roundings.
use std::{
    collections::VecDeque,
    sync::{Arc, Mutex, OnceLock},
};

use crate::bigfloat::{pi, BigFloat, Exponent};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Key {
    precision: u32,
    a: i64,
}

#[derive(Debug, PartialEq, Eq)]
pub(super) struct Coefficients {
    pub(super) values: Box<[BigFloat]>,
}
impl Coefficients {
    fn weight(&self) -> usize {
        self.values
            .iter()
            .fold(std::mem::size_of_val(self.values.as_ref()), |sum, c| {
                let exponent = match &c.exp {
                    Exponent::Big(e) => e.bits().div_ceil(8) as usize,
                    Exponent::Small(_) => 0,
                };
                sum.saturating_add(c.man.bits().div_ceil(8) as usize)
                    .saturating_add(exponent)
            })
    }
}
struct Entry {
    key: Key,
    coefficients: Arc<Coefficients>,
    weight: usize,
}
#[derive(Default)]
struct Cache {
    entries: VecDeque<Entry>,
    weight: usize,
}
const MAX_ENTRIES: usize = 16;
const MAX_WEIGHT: usize = 2 * 1024 * 1024;
static CACHE: OnceLock<Mutex<Cache>> = OnceLock::new();

impl Cache {
    fn get(&mut self, key: Key) -> Option<Arc<Coefficients>> {
        let position = self.entries.iter().position(|entry| entry.key == key)?;
        let entry = self.entries.remove(position).unwrap();
        let value = Arc::clone(&entry.coefficients);
        self.entries.push_back(entry);
        Some(value)
    }

    // Return evictions to the caller so their mantissas are freed after the
    // mutex is released. A single oversized set remains usable, but unretained.
    fn insert(&mut self, key: Key, coefficients: Arc<Coefficients>, weight: usize) -> Vec<Entry> {
        let mut evicted = Vec::new();
        if weight > MAX_WEIGHT {
            return evicted;
        }
        while self.entries.len() >= MAX_ENTRIES || self.weight + weight > MAX_WEIGHT {
            let entry = self.entries.pop_front().unwrap();
            self.weight -= entry.weight;
            evicted.push(entry);
        }
        self.weight += weight;
        self.entries.push_back(Entry {
            key,
            coefficients,
            weight,
        });
        evicted
    }
}

fn visit(key: Key, mut consume: impl FnMut(i64, BigFloat) -> Option<()>) -> Option<()> {
    if super::host_interrupted() {
        return None;
    }
    let wp = key.precision;
    consume(0, pi(wp).ldexp(1).sqrt(wp)?)?;
    let mut fact = BigFloat::one(); // (k-1)!
    for k in 1..key.a {
        if super::host_interrupted() {
            return None;
        }
        // c_k = (-1)^(k-1) (a-k)^(k-1/2) e^(a-k) / (k-1)!
        let amk = BigFloat::from_i64(key.a - k);
        let e1 = amk.powf(&BigFloat::from_i64(2 * k - 1).ldexp(-1), wp)?;
        let e2 = amk.checked_exp(wp)?;
        let mut coefficient = e1.mul(&e2, wp).div(&fact, wp)?;
        if k % 2 == 0 {
            coefficient = coefficient.neg();
        }
        consume(k, coefficient)?;
        fact = fact.mul(&BigFloat::from_i64(k), wp);
    }
    if super::host_interrupted() {
        return None;
    }
    Some(())
}

fn compute(key: Key) -> Option<Coefficients> {
    let mut values = Vec::new();
    visit(key, |_, value| {
        values.push(value);
        Some(())
    })?;
    Some(Coefficients {
        values: values.into_boxed_slice(),
    })
}

fn cacheable(key: Key) -> bool {
    let count = usize::try_from(key.a.max(1)).unwrap_or(usize::MAX);
    // Mantissas are rounded to wp bits; include scalar/exponent overhead.
    let per = std::mem::size_of::<BigFloat>() + (key.precision as usize).div_ceil(8) + 32;
    count
        .checked_mul(per)
        .is_some_and(|weight| weight <= MAX_WEIGHT)
}

fn sum_uncached(z: &BigFloat, key: Key) -> Option<BigFloat> {
    let mut sum = BigFloat::zero();
    visit(key, |k, coefficient| {
        if k == 0 {
            sum = coefficient;
        } else {
            let denominator = z.add(&BigFloat::from_i64(k), key.precision);
            sum = sum.add(
                &coefficient.div(&denominator, key.precision)?,
                key.precision,
            );
        }
        Some(())
    })?;
    Some(sum)
}

pub(super) fn sum(z: &BigFloat, precision: u32, a: i64) -> Option<BigFloat> {
    let key = Key { precision, a };
    if !cacheable(key) {
        return sum_uncached(z, key);
    }
    let coefficients = coefficients(precision, a)?;
    let mut sum = coefficients.values[0].clone();
    for (k, coefficient) in coefficients.values.iter().enumerate().skip(1) {
        if super::host_interrupted() {
            return None;
        }
        let denominator = z.add(&BigFloat::from_i64(k as i64), precision);
        sum = sum.add(&coefficient.div(&denominator, precision)?, precision);
    }
    Some(sum)
}

fn coefficients(precision: u32, a: i64) -> Option<Arc<Coefficients>> {
    if super::host_interrupted() {
        return None;
    }
    let key = Key { precision, a };
    let cache = CACHE.get_or_init(|| Mutex::new(Cache::default()));
    if let Some(value) = cache.lock().unwrap_or_else(|e| e.into_inner()).get(key) {
        return Some(value);
    }
    // Concurrent misses may build the same key independently. No arithmetic,
    // including pi/exp/pow, runs while holding the shared cache mutex.
    let candidate = Arc::new(compute(key)?);
    let weight = candidate.weight();
    let mut cache = cache.lock().unwrap_or_else(|e| e.into_inner());
    if let Some(value) = cache.get(key) {
        return Some(value);
    }
    let evicted = cache.insert(key, Arc::clone(&candidate), weight);
    drop(cache);
    drop(evicted);
    Some(candidate)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cached_coefficients_keep_full_precision_and_parameter() {
        let mut cache = Cache::default();
        for key in [
            Key {
                precision: 96,
                a: 12,
            },
            Key {
                precision: 97,
                a: 12,
            },
            Key {
                precision: 96,
                a: 13,
            },
        ] {
            let fresh = Arc::new(compute(key).unwrap());
            cache.insert(key, Arc::clone(&fresh), fresh.weight());
            assert!(Arc::ptr_eq(&fresh, &cache.get(key).unwrap()));
            assert_eq!(
                coefficients(key.precision, key.a).unwrap().as_ref(),
                fresh.as_ref()
            );
        }
        assert_eq!(cache.entries.len(), 3);
        assert!(cache
            .get(Key {
                precision: 98,
                a: 12
            })
            .is_none());
    }

    #[test]
    fn retained_coefficients_have_entry_and_storage_bounds() {
        let mut cache = Cache::default();
        for a in 1..=MAX_ENTRIES as i64 + 3 {
            let value = Arc::new(Coefficients {
                values: vec![BigFloat::one()].into_boxed_slice(),
            });
            cache.insert(Key { precision: 32, a }, Arc::clone(&value), value.weight());
            assert!(cache.entries.len() <= MAX_ENTRIES);
            assert!(cache.weight <= MAX_WEIGHT);
        }
        assert!(cache
            .get(Key {
                precision: 32,
                a: 1
            })
            .is_none());
        let huge = Arc::new(Coefficients {
            values: vec![BigFloat::one(); MAX_WEIGHT / std::mem::size_of::<BigFloat>() + 1]
                .into_boxed_slice(),
        });
        let before = cache.weight;
        cache.insert(
            Key {
                precision: 32,
                a: 100,
            },
            Arc::clone(&huge),
            huge.weight(),
        );
        assert_eq!(cache.weight, before);
        assert!(cache
            .get(Key {
                precision: 32,
                a: 100
            })
            .is_none());
    }

    #[test]
    fn cached_and_streamed_sums_have_identical_rounding() {
        for key in [
            Key {
                precision: 96,
                a: 12,
            },
            Key {
                precision: 240,
                a: 40,
            },
        ] {
            for numerator in [2, 7, 49] {
                let z = BigFloat::from_i64(numerator).div_i64(3, key.precision);
                assert_eq!(sum(&z, key.precision, key.a), sum_uncached(&z, key));
            }
        }
        assert!(!cacheable(Key {
            precision: u32::MAX,
            a: i64::MAX
        }));
        assert!(cacheable(Key {
            precision: 240,
            a: 40
        }));
    }

    #[test]
    fn concurrent_misses_preserve_coefficients() {
        let barrier = Arc::new(std::sync::Barrier::new(6));
        let handles = (0..6)
            .map(|i| {
                let barrier = Arc::clone(&barrier);
                std::thread::spawn(move || {
                    barrier.wait();
                    let key = Key {
                        precision: 160 + (i % 2),
                        a: 24,
                    };
                    let value = coefficients(key.precision, key.a).unwrap();
                    assert_eq!(value.as_ref(), &compute(key).unwrap());
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            handle.join().unwrap();
        }
    }

    #[test]
    fn a_warm_coefficient_cache_does_not_hide_a_pending_interrupt() {
        let _scope = super::super::HostErrorScope::enter(true);
        assert!(coefficients(100, 16).is_some());
        super::super::HOST_FATAL.with(|cell| cell.set(true));
        assert!(coefficients(100, 16).is_none());
        assert!(compute(Key {
            precision: 100,
            a: 16
        })
        .is_none());
        assert!(super::super::gamma_mp(&BigFloat::from_i64(3), 80).is_none());
        assert!(super::super::host_interrupted());
        super::super::HOST_FATAL.with(|cell| cell.set(false));
    }
}
