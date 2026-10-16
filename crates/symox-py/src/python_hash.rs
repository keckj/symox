//! Python's numeric hash law: equal numeric values must share dictionary keys.
//! https://docs.python.org/3/library/stdtypes.html#hashing-of-numeric-types
use num_bigint::{BigInt, Sign};
use num_traits::ToPrimitive;
use symox::number::Number;

const MODULUS: u64 = if isize::BITS == 64 { (1 << 61) - 1 } else { (1 << 31) - 1 };

fn modpow(mut base: u64, mut exp: u64) -> u64 {
    let mut value = 1;
    while exp != 0 {
        if exp & 1 != 0 { value = ((value as u128 * base as u128) % MODULUS as u128) as u64; }
        base = ((base as u128 * base as u128) % MODULUS as u128) as u64;
        exp >>= 1;
    }
    value
}

fn remainder(n: &BigInt) -> u64 { (n.magnitude() % MODULUS).to_u64().unwrap() }

fn signed(value: u64, negative: bool) -> isize {
    let h = if negative { -(value as isize) } else { value as isize };
    if h == -1 { -2 } else { h }
}

pub fn exact_hash(n: &Number) -> Option<isize> {
    Some(match n {
        Number::Int(n) => signed(n.unsigned_abs() % MODULUS, *n < 0),
        Number::BigInt(n) => signed(remainder(n), n.sign() == Sign::Minus),
        Number::Rational(q) => {
            let denominator = remainder(q.denom());
            let value = if denominator == 0 { 314159 } else {
                ((remainder(q.numer()) as u128 * modpow(denominator, MODULUS - 2) as u128) % MODULUS as u128) as u64
            };
            signed(value, q.numer().sign() == Sign::Minus)
        }
        Number::MpFloat(value, _) => {
            // Reducing the binary exponent avoids constructing a denominator
            // proportional to its magnitude, including subnormal values.
            let exponent = value.exp.rem_euclid((MODULUS - 1) as i64) as u64;
            let h = ((remainder(&value.man) as u128 * modpow(2, exponent) as u128) % MODULUS as u128) as u64;
            signed(h, value.is_negative())
        }
        Number::Float(_) => return None,
    })
}
