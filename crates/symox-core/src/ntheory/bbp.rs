//! Hexadecimal digits of pi by the BBP formula in fixed point, ported from
//! SymPy's `ntheory.bbp_pi`.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

fn dn(n: i64, prec: u64) -> u64 {
    let n = n + 1;
    let v = (n as u64) + prec;
    (64 - v.leading_zeros() as u64 - 1) / 4 + prec + 3
}

fn series(j: u64, n: i64, prec: u64) -> BigInt {
    let d4 = 4 * dn(n, prec);
    let mut s = BigInt::zero();
    let mut d = j;
    for k in 0..=n {
        let bd = BigInt::from(d);
        let r = BigInt::from(16).modpow(&BigInt::from(n - k), &bd);
        s += (r << d4) / &bd;
        d += 8;
    }
    let mut t = BigInt::zero();
    let k = (n + 1) as u64;
    let mut e = d4 as i64 - 4;
    let mut d = 8 * k + j;
    while e >= 0 {
        let dt = (BigInt::one() << e as u64) / d;
        if dt.is_zero() {
            break;
        }
        t += dt;
        e -= 4;
        d += 8;
    }
    s + t
}

/// SymPy's `pi_hex_digits(n, prec)`: `prec` hexadecimal digits of pi
/// starting at the `n`-th one (`n >= 0`, `prec >= 0`).
pub fn pi_hex_digits(n: u64, prec: u64) -> String {
    if prec == 0 {
        return String::new();
    }
    let n = n as i64 - 1;
    let d = dn(n, prec);
    let x: BigInt = 4 * series(1, n, prec) - 2 * series(4, n, prec) - series(5, n, prec) - series(6, n, prec);
    let modulus = BigInt::one() << (4 * d);
    let x = x.mod_floor(&modulus) >> (4 * (d - prec));
    let s = x.to_str_radix(16);
    let width = prec as usize;
    if s.len() < width {
        format!("{}{}", "0".repeat(width - s.len()), s)
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digits() {
        assert_eq!(pi_hex_digits(0, 14), "3243f6a8885a30");
        assert_eq!(pi_hex_digits(1, 14), "243f6a8885a308");
        assert_eq!(pi_hex_digits(10, 14), "5a308d313198a2");
        assert_eq!(pi_hex_digits(91, 10), "e90c6cc0ac");
        assert_eq!(pi_hex_digits(0, 0), "");
    }
}
