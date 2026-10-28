//! Automatic evaluation of combinatorial functions.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::integers as ints;
use super::*;
use crate::functions::gamma;
use crate::ops::{add, div, mul, mul2, neg, pow, sub};

/// Index bound of the evaluations whose cost is not governed by the size
/// of one exact integer (polynomials, harmonic sums, Bell numbers, ...).
const MAX_EXACT: u64 = 100_000;

/// `log2(F_n)` and `log2(L_n)` are below `n log2(golden ratio) + 1`.
/// Index bound of the Bernoulli and Euler numbers (zeta and Dirichlet beta
/// methods, `summation::bernoulli_zeta`, `integers::euler_beta`); the reference needs 19 s for
/// `bernoulli(50000)` and its cost grows faster than `n^2`.
const MAX_BERNOULLI_NUMBER: i64 = 200_000;

const LOG2_GOLDEN: f64 = 0.694_241_913_630_617_3;

fn num(n: BigInt) -> Expr {
    Expr::num(Number::from_bigint(n))
}

fn numq(r: BigRational) -> Expr {
    Expr::num(Number::from_ratio(r))
}

/// Polynomial `sum c_k x^k` from rational coefficients (constant first).
pub(crate) fn poly_q(c: &[BigRational], x: &Expr) -> Expr {
    let terms: Vec<Expr> = c
        .iter()
        .enumerate()
        .filter(|(_, k)| !k.is_zero())
        .map(|(i, k)| mul2(&numq(k.clone()), &x.powi(i as i64)))
        .collect();
    add(&terms)
}

fn poly_z(c: &[BigInt], x: &Expr) -> Expr {
    let terms: Vec<Expr> = c
        .iter()
        .enumerate()
        .filter(|(_, k)| !k.is_zero())
        .map(|(i, k)| mul2(&num(k.clone()), &x.powi(i as i64)))
        .collect();
    add(&terms)
}

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        Factorial => eval_factorial(&a[0]),
        Factorial2 => eval_factorial2(&a[0]),
        Subfactorial => {
            if is_inf(&a[0]) {
                return Some(Expr::oo());
            }
            let n = as_small(&a[0])?;
            if n < 0 || !ints::fits_exact(ints::log2_factorial(n as f64)) {
                return None;
            }
            Some(num(ints::subfactorial(n as u64)))
        }
        Binomial => eval_binomial(&a[0], &a[1]),
        RisingFactorial => eval_rf(&a[0], &a[1], true),
        FallingFactorial => eval_rf(&a[0], &a[1], false),
        Harmonic => eval_harmonic(a),
        Bernoulli => eval_bernoulli(a),
        Euler => eval_euler(a),
        Catalan => eval_catalan(&a[0]),
        Fibonacci => eval_fibonacci(a),
        Lucas => {
            if is_inf(&a[0]) {
                return Some(Expr::oo());
            }
            let n = as_small(&a[0])?;
            if !ints::fits_exact(n.unsigned_abs() as f64 * LOG2_GOLDEN + 2.0) {
                return None;
            }
            Some(num(ints::lucas_number(n)))
        }
        Bell => eval_bell(a),
        Partition => {
            // SymPy: n.is_negative -> 0, n.is_zero -> 1 also for symbols.
            if as_small(&a[0]).is_none() {
                if is_integer(&a[0]) != Some(false) {
                    if is_negative(&a[0]) == Some(true) {
                        return Some(Expr::zero());
                    }
                    if is_zero(&a[0]) == Some(true) {
                        return Some(Expr::one());
                    }
                }
                return None;
            }
            let n = as_small(&a[0])?;
            if n < 0 {
                return Some(Expr::zero());
            }
            if n as u64 > MAX_EXACT {
                return None;
            }
            Some(num(ints::partition_number(n as u64)))
        }
        _ => None,
    }
}

fn eval_factorial(n: &Expr) -> Option<Expr> {
    if let Some(k) = as_int(n) {
        if k.is_negative() {
            return Some(Expr::zoo());
        }
        let k = k.to_u64()?;
        if !ints::fits_exact(ints::log2_factorial(k as f64)) {
            return None;
        }
        return Some(num(ints::factorial(k)));
    }
    if is_inf(n) {
        return Some(Expr::oo());
    }
    if !n.is_number() && is_integer(n) == Some(true) && is_negative(n) == Some(true) {
        return Some(Expr::zoo());
    }
    None
}

fn eval_factorial2(n: &Expr) -> Option<Expr> {
    let k = as_small(n)?;
    if k >= 0 {
        // n!! <= sqrt(n! (n + 1)) for every n >= 0
        if !ints::fits_exact(ints::log2_factorial(k as f64) / 2.0 + 0.5 * (k as f64 + 1.0).log2() + 2.0) {
            return None;
        }
        return Some(num(ints::factorial2(k as u64)));
    }
    if k % 2 == 0 {
        // SymPy raises ValueError
        return None;
    }
    // n!! = (-1)^((-n-1)/2) / (-n-2)!!
    let m = -k - 2;
    let d = if m < 0 { BigInt::one() } else { ints::factorial2(m as u64) };
    let sign = if ((-k - 1) / 2) % 2 == 0 { 1 } else { -1 };
    Some(Expr::num(Number::rational(BigInt::from(sign), d).unwrap()))
}

fn eval_binomial(n: &Expr, k: &Expr) -> Option<Expr> {
    let d = sub(n, k);
    let n_nonneg = is_nonnegative(n);
    let n_isint = is_integer(n);
    let cond = n_nonneg == Some(true) || n_isint == Some(false);
    if is_zero(k) == Some(true) || (cond && is_zero(&d) == Some(true)) {
        return Some(Expr::one());
    }
    let km1 = sub(k, &Expr::one());
    if is_zero(&km1) == Some(true) || (cond && is_zero(&sub(&d, &Expr::one())) == Some(true)) {
        return Some(n.clone());
    }
    if is_integer(k) == Some(true) {
        if is_negative(k) == Some(true)
            || (n_nonneg == Some(true) && n_isint == Some(true) && is_negative(&d) == Some(true))
        {
            return Some(Expr::zero());
        }
        if crate::functions::is_numeric_value(n) {
            let kk = as_int(k)?;
            if let Some(nn) = as_int(n) {
                let positive_n = if nn.is_negative() { &kk - &nn - 1 } else { nn.clone() };
                if kk <= positive_n {
                    // With min(k,n-k)>u64::MAX the exact output already needs
                    // more than u64::MAX bits. Preserve the symbolic binomial;
                    // modular evaluation can still compute a small exact result.
                    let m = std::cmp::min(kk.clone(), &positive_n - &kk).to_u64()?;
                    let total = positive_n.to_f64()?;
                    if !ints::fits_exact(ints::log2_binomial(total, m as f64)) {
                        return None;
                    }
                }
                return Some(num(ints::binomial(&nn, &kk)));
            }
            // non integer number: falling product / k!
            let kk = kk.to_u64()?;
            if kk > MAX_EXACT {
                return None;
            }
            let mut v = Vec::with_capacity(kk as usize + 1);
            for i in 1..=kk {
                v.push(div(&add(&[n.clone(), Expr::int(i as i64 - kk as i64)]), &Expr::int(i as i64)));
            }
            let p = mul(&v);
            return Some(crate::expand::expand(&p));
        }
        return None;
    }
    if n_nonneg == Some(false) && n_isint == Some(true) {
        return Some(Expr::zoo());
    }
    if crate::functions::is_numeric_value(k) {
        let one = Expr::one();
        return Some(div(
            &gamma(&add(&[n.clone(), one.clone()])),
            &mul2(&gamma(&add(&[k.clone(), one.clone()])), &gamma(&add(&[d, one]))),
        ));
    }
    None
}

/// Rising (`up`) or falling factorial.
fn eval_rf(x: &Expr, k: &Expr, up: bool) -> Option<Expr> {
    if up && x.is_one() {
        return Some(super::factorial(k));
    }
    // FallingFactorial.eval: `k.is_integer and x == k` gives factorial(x).
    if !up && x == k && crate::assumptions::is_integer(k) == Some(true) {
        return Some(super::factorial(x));
    }
    if let Some(kk) = as_small(k) {
        if kk == 0 {
            return Some(Expr::one());
        }
        if let Some(v) = integer_rf(x, kk, up) {
            return Some(v);
        }
        if kk.unsigned_abs() > MAX_EXACT {
            return None;
        }
        if kk > 0 {
            if is_inf(x) {
                return Some(Expr::oo());
            }
            if is_neg_inf(x) {
                // rf and ff alike: (-oo)**k by the parity of k
                return Some(if kk % 2 == 1 { Expr::neg_oo() } else { Expr::oo() });
            }
            let mut r = Expr::one();
            for i in 0..kk {
                let t = if up { add(&[x.clone(), Expr::int(i)]) } else { sub(x, &Expr::int(i)) };
                r = mul2(&r, &t);
            }
            return Some(r);
        }
        if is_inf(x) || is_neg_inf(x) {
            return Some(Expr::oo());
        }
        let mut r = Expr::one();
        for i in 1..=(-kk) {
            let t = if up { sub(x, &Expr::int(i)) } else { add(&[x.clone(), Expr::int(i)]) };
            r = mul2(&r, &t);
        }
        return Some(pow(&r, &Expr::minus_one()));
    }
    // Only RisingFactorial.eval has this rule in the reference.
    if up && is_integer(k) == Some(false) && is_integer(x) == Some(true) && is_negative(x) == Some(true) {
        return Some(Expr::zero());
    }
    None
}

/// `rf(x, k)` / `ff(x, k)` for an exact machine integer `x` as one balanced
/// product of the integer factors (the generic path multiplies them one by
/// one). `None` when `x` is not such an integer or the result is too large.
fn integer_rf(x: &Expr, k: i64, up: bool) -> Option<Expr> {
    let x = as_small(x)? as i128;
    let m = k.unsigned_abs() as i128;
    // rf(x, k) = x (x+1) ... (x+k-1), rf(x, -k) = 1/((x-1) ... (x-k)),
    // ff(x, k) = x (x-1) ... (x-k+1), ff(x, -k) = 1/((x+1) ... (x+k))
    let (lo, hi) = match (up, k > 0) {
        (true, true) => (x, x + m - 1),
        (true, false) => (x - m, x - 1),
        (false, true) => (x - m + 1, x),
        (false, false) => (x + 1, x + m),
    };
    let product = if lo <= 0 && hi >= 0 {
        BigInt::zero()
    } else {
        let (a, b) = if hi < 0 { (-hi, -lo) } else { (lo, hi) };
        let a = u64::try_from(a).ok()?;
        let b = u64::try_from(b).ok()?;
        if !ints::fits_exact(m as f64 * (b as f64 + 1.0).log2()) {
            return None;
        }
        let factors: Vec<u64> = (a..=b).collect();
        let p = crate::ntheory::product_u64(&factors);
        if hi < 0 && m % 2 == 1 { -p } else { p }
    };
    let p = num(product);
    Some(if k > 0 { p } else { pow(&p, &Expr::minus_one()) })
}

fn eval_harmonic(a: &[Expr]) -> Option<Expr> {
    let n = &a[0];
    let m = a.get(1).cloned().unwrap_or_else(Expr::one);
    if a.len() == 2 && m.is_one() {
        return Some(super::harmonic(n));
    }
    if is_zero(&m) == Some(true) {
        return Some(n.clone());
    }
    if is_inf(n) {
        if is_negative(&m) == Some(true) {
            return Some(Expr::nan());
        }
        if let Some(mv) = as_exact(&m) {
            if mv.cmp_value(&Number::ONE) != Some(std::cmp::Ordering::Greater) {
                return Some(Expr::oo());
            }
            return Some(super::zeta(&m));
        }
        if m.is_constant(Constant::Infinity) {
            return Some(Expr::one());
        }
        // SymPy: is_le(m, 1) gives oo, is_gt(m, 1) gives zeta(m).
        let shifted = sub(&m, &Expr::one());
        if is_nonpositive(&shifted) == Some(true) {
            return Some(Expr::oo());
        }
        if is_positive(&shifted) == Some(true) {
            return Some(super::zeta(&m));
        }
        return None;
    }
    if let Some(mm) = as_small(&m) {
        if mm <= 0 {
            // (B_{1-m}(n+1) - B_{1-m}) / (1-m)
            let p = (1 - mm) as usize;
            if p > 1000 {
                return None;
            }
            let np1 = add(&[n.clone(), Expr::one()]);
            let bp = poly_q(&ints::bernoulli_poly_coeffs(p), &np1);
            let b0 = numq(ints::bernoulli_number(p));
            return Some(div(&sub(&bp, &b0), &Expr::int(p as i64)));
        }
    }
    if let Some(nn) = as_small(n) {
        if nn < 0 {
            let m_int = is_integer(&m);
            let m_nonpos = is_nonpositive(&m);
            if m_int == Some(false) || m_nonpos == Some(false) {
                return Some(if m.is_one() { Expr::zoo() } else { Expr::nan() });
            }
            return None;
        }
        if nn as u64 > MAX_EXACT {
            return None;
        }
        if let Some(mm) = as_small(&m) {
            if mm > 0 && mm < 1000 {
                return Some(numq(ints::harmonic_number(nn as u64, mm as u32)));
            }
        }
        // generic exponent: sum k^-m
        let v: Vec<Expr> = (1..=nn).map(|k| pow(&Expr::int(k), &neg(&m))).collect();
        return Some(add(&v));
    }
    None
}

fn eval_bernoulli(a: &[Expr]) -> Option<Expr> {
    let n = &a[0];
    // SymPy: bernoulli(n) = 0 for odd n > 1
    if a.len() == 1
        && is_nonnegative(n) == Some(true)
        && crate::assumptions::is_odd(n) == Some(true)
        && is_positive(&sub(n, &Expr::one())) == Some(true)
    {
        return Some(Expr::zero());
    }
    let nn = as_small(n)?;
    if nn < 0 {
        return None;
    }
    if a.len() == 1 {
        if nn > MAX_BERNOULLI_NUMBER {
            return None;
        }
        return Some(numq(ints::bernoulli_number(nn as usize)));
    }
    if nn > 5000 {
        return None;
    }
    let x = &a[1];
    Some(poly_q(&ints::bernoulli_poly_coeffs(nn as usize), x))
}

fn eval_euler(a: &[Expr]) -> Option<Expr> {
    let n = &a[0];
    if is_zero(n) == Some(true) {
        return Some(Expr::one());
    }
    if as_small(n) == Some(-1) {
        // SymPy: euler(-1) = pi/2, euler(-1, x) = digamma((x+1)/2) - digamma(x/2)
        let Some(x) = a.get(1) else {
            return Some(div(&Expr::pi(), &Expr::int(2)));
        };
        let half = |e: &Expr| div(e, &Expr::int(2));
        let digamma = |e: Expr| sf(SpecialFn::Digamma, vec![e]);
        return Some(sub(&digamma(half(&add(&[x.clone(), Expr::one()]))), &digamma(half(x))));
    }
    if a.len() == 1 && !n.is_number() && is_integer(n) == Some(true)
        && crate::assumptions::is_odd(n) == Some(true) && is_positive(n) == Some(true) {
        return Some(Expr::zero());
    }
    let nn = as_small(n)?;
    if nn < 0 {
        return None;
    }
    if a.len() == 1 {
        if nn % 2 == 1 {
            return Some(Expr::zero());
        }
        if nn > MAX_BERNOULLI_NUMBER {
            return None;
        }
        return Some(num(ints::euler_number(nn as usize)));
    }
    if nn > 5000 {
        return None;
    }
    Some(poly_q(&ints::euler_poly_coeffs(nn as usize), &a[1]))
}

fn eval_catalan(n: &Expr) -> Option<Expr> {
    if let Some(k) = as_small(n) {
        if k >= 0 {
            if !ints::fits_exact(2.0 * k as f64) {
                return None;
            }
            return Some(num(ints::catalan_number(k as u64)));
        }
        if k == -1 {
            return Some(rat(-1, 2));
        }
        return Some(Expr::zero());
    }
    if !n.is_number() && is_integer(n) == Some(true) && is_negative(n) == Some(true) {
        let n1 = add(&[n.clone(), Expr::one()]);
        if is_negative(&n1) == Some(true) {
            return Some(Expr::zero());
        }
        if is_zero(&n1) == Some(true) {
            return Some(rat(-1, 2));
        }
    }
    None
}

fn eval_fibonacci(a: &[Expr]) -> Option<Expr> {
    if is_inf(&a[0]) {
        return Some(Expr::oo());
    }
    let nn = as_small(&a[0])?;
    if a.len() == 1 {
        if !ints::fits_exact(nn.unsigned_abs() as f64 * LOG2_GOLDEN + 2.0) {
            return None;
        }
        return Some(num(ints::fibonacci_number(nn)));
    }
    if nn < 0 || nn as u64 > MAX_EXACT {
        return None;
    }
    Some(poly_z(&ints::fibonacci_poly_coeffs(nn as u64), &a[1]))
}

fn eval_bell(a: &[Expr]) -> Option<Expr> {
    if is_inf(&a[0]) && a.len() == 1 {
        return Some(Expr::oo());
    }
    let nn = as_small(&a[0])?;
    if nn < 0 || nn as u64 > MAX_EXACT {
        return None;
    }
    match a.len() {
        1 => Some(num(ints::bell_number(nn as u64))),
        2 => Some(poly_z(&ints::bell_poly_coeffs(nn as u64), &a[1])),
        _ => {
            let k = as_small(&a[1])?;
            let xs = tuple_items(&a[2])?;
            Some(incomplete_bell(nn, k, xs))
        }
    }
}

/// Incomplete Bell polynomial `B_{n,k}(x_1, x_2, ...)`.
fn incomplete_bell(n: i64, k: i64, xs: &[Expr]) -> Expr {
    if n == 0 && k == 0 {
        return Expr::one();
    }
    if n <= 0 || k <= 0 {
        return Expr::zero();
    }
    // B_{n,k} = sum_{m=1}^{n-k+1} C(n-1, m-1) x_m B_{n-m,k-1}
    let n = n as usize;
    let k = k as usize;
    // table[i][j] = B_{i,j}
    let mut t: Vec<Vec<Expr>> = vec![vec![Expr::zero(); k + 1]; n + 1];
    t[0][0] = Expr::one();
    for i in 1..=n {
        for j in 1..=k.min(i) {
            let mut terms = Vec::new();
            let mut c = BigInt::one(); // C(i-1, m-1)
            for m in 1..=(i - j + 1) {
                if m - 1 < xs.len() && !t[i - m][j - 1].is_zero() {
                    terms.push(mul(&[num(c.clone()), xs[m - 1].clone(), t[i - m][j - 1].clone()]));
                }
                c = c * BigInt::from((i - m) as u64) / BigInt::from(m as u64);
            }
            t[i][j] = crate::expand::expand(&add(&terms));
        }
    }
    t[n][k].clone()
}

#[cfg(test)]
mod bigexact_tests {
    use super::*;

    #[test]
    fn integer_rising_factorial_beyond_old_cap() {
        let three = Expr::int(3);
        let k = Expr::int(100_001);
        let direct = integer_rf(&three, 100_001, true);
        assert!(direct.is_some());
        assert!(eval_rf(&three, &k, true).is_some());
        let zero = integer_rf(&Expr::int(-3), 100_001, true).unwrap();
        assert!(zero.is_zero());
        let via_apply = crate::special::apply(SpecialFn::RisingFactorial, vec![Expr::int(-3), Expr::int(100_001)]);
        assert!(via_apply.is_zero(), "apply: {via_apply}");
        let parsed = crate::parser::parse_str("rf(-3, 100001)").unwrap();
        assert!(parsed.is_zero(), "parsed: {parsed}");
    }
}
