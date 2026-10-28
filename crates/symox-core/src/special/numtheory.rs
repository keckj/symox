//! Number theoretic function heads (`totient`, `mobius`, `primepi`,
//! `legendre_symbol`, ...), ported from SymPy's
//! `functions.combinatorial.numbers`: the same automatic evaluation, the
//! same argument errors and the same class assumptions.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use super::*;
use crate::ntheory as nt;
use crate::ops::{add, add2, mul, mul2, pow, sub};

fn num(n: BigInt) -> Expr {
    Expr::num(Number::from_bigint(n))
}

fn is_prime(e: &Expr) -> Option<bool> {
    if let Some(n) = as_int(e) {
        return Some(nt::isprime(&n));
    }
    crate::assumptions::fact(e, "prime")
}

fn is_one(e: &Expr) -> bool {
    as_int(e).is_some_and(|n| n.is_one())
}

/// `1 < e` as a fuzzy boolean.
fn one_lt(e: &Expr) -> Option<bool> {
    if let Some(n) = as_exact(e) {
        return Some(n.cmp_value(&Number::Int(1)) == Some(std::cmp::Ordering::Greater));
    }
    is_positive(&sub(e, &Expr::one()))
}

/// Is `a % p` known to vanish (integers, or `c*rest` with `p | c`)?
fn mod_is_zero(a: &Expr, p: &Expr) -> bool {
    let Some(p) = as_int(p) else { return false };
    if p.is_zero() {
        return false;
    }
    if let Some(a) = as_int(a) {
        return a.mod_floor(&p).is_zero();
    }
    match a.kind() {
        Kind::Mul(c, facs) => {
            let Some(c) = c.to_bigint() else { return false };
            c.mod_floor(&p).is_zero()
                && facs.iter().all(|(b, x)| is_integer(b) == Some(true) && x.as_i64().is_some_and(|k| k >= 0))
        }
        _ => false,
    }
}

/// The SymPy exception (`true`: TypeError, `false`: ValueError) raised by
/// the constructor for these arguments, if any.
pub fn check_args(s: SpecialFn, a: &[Expr]) -> Result<(), (bool, String)> {
    use SpecialFn::*;
    let te = |m: &str| Err((true, m.to_string()));
    let ve = |m: &str| Err((false, m.to_string()));
    match s {
        Totient | ReducedTotient | Mobius | PrimeNu | PrimeOmega | DivisorSigma | UdivisorSigma => {
            if is_integer(&a[0]) == Some(false) {
                return te("n should be an integer");
            }
            if is_positive(&a[0]) == Some(false) {
                return ve("n should be a positive integer");
            }
            if a.len() > 1 {
                if is_integer(&a[1]) == Some(false) {
                    return te("k should be an integer");
                }
                if is_nonnegative(&a[1]) == Some(false) {
                    return ve("k should be a nonnegative integer");
                }
            }
            Ok(())
        }
        LegendreSymbol => {
            if is_integer(&a[0]) == Some(false) {
                return te("a should be an integer");
            }
            if is_integer(&a[1]) == Some(false) {
                return te("p should be an integer");
            }
            if is_prime(&a[1]) == Some(false) || crate::assumptions::is_odd(&a[1]) == Some(false) {
                return ve("p should be an odd prime integer");
            }
            Ok(())
        }
        JacobiSymbol => {
            if is_integer(&a[0]) == Some(false) {
                return te("m should be an integer");
            }
            if is_integer(&a[1]) == Some(false) {
                return te("n should be an integer");
            }
            if is_positive(&a[1]) == Some(false) || crate::assumptions::is_odd(&a[1]) == Some(false) {
                return ve("n should be an odd positive integer");
            }
            Ok(())
        }
        KroneckerSymbol => {
            if is_integer(&a[0]) == Some(false) {
                return te("a should be an integer");
            }
            if is_integer(&a[1]) == Some(false) {
                return te("n should be an integer");
            }
            Ok(())
        }
        PrimePi => {
            let n = &a[0];
            if is_inf(n) || is_neg_inf(n) {
                return Ok(());
            }
            if n.is_constant(Constant::NaN) {
                return te("Invalid NaN comparison");
            }
            if crate::assumptions::is_real(n) == Some(false) {
                return te("n should be a real");
            }
            Ok(())
        }
        Tribonacci => match as_int(&a[0]) {
            Some(n) if n.is_negative() => ve("Tribonacci polynomials are defined only for non-negative integer indices."),
            _ => Ok(()),
        },
        Motzkin => match as_int(&a[0]) {
            Some(n) if !n.is_negative() => Ok(()),
            _ => ve("The provided number must be a positive integer"),
        },
        // SymPy's laguerre.eval rejects a non-integer degree first.
        Laguerre if is_integer(&a[0]) == Some(false) => ve("Error: n should be an integer."),
        // hermite.eval: a numerical degree is evaluated as a polynomial and
        // must not be negative.
        Hermite => match a[0].as_number() {
            Some(n) if n.is_negative() => Err((false, format!("The index n must be nonnegative integer (got {})", a[0]))),
            _ => Ok(()),
        },
        // assoc_legendre.eval: P^{-m}_n is rewritten through P^m_n, m = 0 is
        // legendre and x = 0 has a closed form; other integer indices must
        // satisfy 0 <= |m| <= n.
        AssocLegendre => {
            let (Some(n), Some(m)) = (as_int(&a[0]), as_int(&a[1])) else { return Ok(()) };
            if m.is_zero() || a[2].is_zero() || !a[0].is_number() || !a[1].is_number() {
                return Ok(());
            }
            let m = m.abs();
            if n.is_negative() {
                return Err((false, format!("assoc_legendre : 1st index must be nonnegative integer (got {n})")));
            }
            if m > n {
                return Err((false, format!("assoc_legendre : abs('2nd index') must be <= '1st index' (got {n}, {m})")));
            }
            Ok(())
        }
        _ => Ok(()),
    }
}

pub(super) fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    if check_args(s, a).is_err() {
        return None;
    }
    match s {
        Totient => {
            let n = &a[0];
            if is_one(n) {
                return Some(Expr::one());
            }
            if let Some(v) = as_int(n) {
                return Some(num(nt::totient(&v)));
            }
            (is_prime(n) == Some(true)).then(|| sub(n, &Expr::one()))
        }
        ReducedTotient => {
            let n = &a[0];
            if is_one(n) {
                return Some(Expr::one());
            }
            if let Some(v) = as_int(n) {
                return Some(num(reduced_totient(&v)));
            }
            (is_prime(n) == Some(true)).then(|| sub(n, &Expr::one()))
        }
        PrimeNu | PrimeOmega => {
            let n = &a[0];
            if let Some(v) = as_int(n) {
                let fs = nt::factorint(&v);
                let c: u64 = if s == PrimeNu { fs.len() as u64 } else { fs.values().map(|&e| e as u64).sum() };
                return Some(num(BigInt::from(c)));
            }
            (is_prime(n) == Some(true)).then(Expr::one)
        }
        Mobius => eval_mobius(&a[0]),
        DivisorSigma | UdivisorSigma => {
            let n = &a[0];
            let k = a.get(1).cloned().unwrap_or_else(Expr::one);
            if is_prime(n) == Some(true) {
                return Some(add2(&Expr::one(), &pow(n, &k)));
            }
            if is_one(n) {
                return Some(Expr::one());
            }
            let v = as_int(n)?;
            let fs = nt::factorint(&v);
            if s == UdivisorSigma {
                let f: Vec<Expr> = fs
                    .iter()
                    .map(|(p, e)| add2(&Expr::one(), &pow(&num(p.clone()), &mul2(&k, &Expr::int(*e as i64)))))
                    .collect();
                return Some(mul(&f));
            }
            if let Some(kv) = as_int(&k) {
                let kv = kv.to_u32()?;
                if kv == 0 {
                    return Some(num(fs.values().map(|&e| BigInt::from(e + 1)).product()));
                }
                return Some(num(nt::divisor_sigma(&v, kv)));
            }
            None
        }
        LegendreSymbol => {
            let (x, p) = (&a[0], &a[1]);
            if mod_is_zero(x, p) {
                return Some(Expr::zero());
            }
            if is_one(x) {
                return Some(Expr::one());
            }
            let (x, p) = (as_int(x)?, as_int(p)?);
            Some(Expr::int(nt::residue::jacobi(&x, &p) as i64))
        }
        JacobiSymbol => {
            let (m, n) = (&a[0], &a[1]);
            if is_one(m) || is_one(n) {
                return Some(Expr::one());
            }
            if mod_is_zero(m, n) {
                return Some(Expr::zero());
            }
            let (m, n) = (as_int(m)?, as_int(n)?);
            Some(Expr::int(nt::residue::jacobi(&m, &n) as i64))
        }
        KroneckerSymbol => {
            let (x, n) = (&a[0], &a[1]);
            if is_one(x) || is_one(n) {
                return Some(Expr::one());
            }
            let (x, n) = (as_int(x)?, as_int(n)?);
            Some(Expr::int(nt::residue::kronecker(&x, &n) as i64))
        }
        PrimePi => eval_primepi(&a[0]),
        Carmichael => None,
        Tribonacci => eval_tribonacci(a),
        Genocchi => eval_genocchi(a),
        Andre => eval_andre(&a[0]),
        Motzkin => {
            let n = as_int(&a[0])?.to_u64()?;
            Some(num(motzkin_number(n)))
        }
        _ => None,
    }
}

/// Carmichael's reduced totient `lambda(n)` (`n >= 1`).
pub fn reduced_totient(n: &BigInt) -> BigInt {
    let (m, t) = nt::factor::remove(n, &BigInt::from(2)).unwrap_or((n.clone(), 0));
    let mut l: BigInt = match t {
        0 => BigInt::one(),
        1 | 2 => BigInt::from(t),
        _ => BigInt::one() << (t - 2),
    };
    for (p, k) in nt::factorint(&m) {
        let v = (&p - 1u32) * p.pow(k - 1);
        l = l.lcm(&v);
    }
    l
}

fn eval_mobius(n: &Expr) -> Option<Expr> {
    if is_prime(n) == Some(true) {
        return Some(Expr::int(-1));
    }
    if let Some(v) = as_int(n) {
        return Some(Expr::int(nt::mobius(&v) as i64));
    }
    let parts: Vec<(Expr, Expr)> = match n.kind() {
        Kind::Mul(c, facs) => {
            let mut v = Vec::new();
            if !c.is_one() {
                v.push((Expr::num(c.clone()), Expr::one()));
            }
            v.extend(facs.iter().cloned());
            v
        }
        Kind::Pow(b, x) => vec![(b.clone(), x.clone())],
        _ => vec![(n.clone(), Expr::one())],
    };
    let mut result: Option<BigInt> = None;
    for (m, e) in parts {
        let ok = is_integer(&m) == Some(true)
            && is_positive(&m) == Some(true)
            && is_integer(&e) == Some(true)
            && is_positive(&e) == Some(true);
        if !ok {
            return None;
        }
        let lt = one_lt(&e);
        if lt == Some(true) {
            result = Some(BigInt::zero());
        } else if let Some(mv) = as_int(&m) {
            let fs = nt::factorint(&mv);
            if fs.values().any(|&v| v > 1) {
                result = Some(BigInt::zero());
            } else if lt == Some(false) {
                let s = if fs.len() % 2 == 1 { -BigInt::one() } else { BigInt::one() };
                result = Some(match result {
                    None => s,
                    Some(r) => r * s,
                });
            }
        }
    }
    result.map(num)
}

fn eval_primepi(n: &Expr) -> Option<Expr> {
    if is_inf(n) {
        return Some(n.clone());
    }
    if is_neg_inf(n) {
        return Some(Expr::zero());
    }
    if one_lt_two(n) == Some(true) {
        return Some(Expr::zero());
    }
    let v = match as_exact(n) {
        Some(x) => x.floor().to_bigint()?,
        None => {
            if let Kind::Num(x) = n.kind() {
                x.floor().to_bigint()?
            } else if crate::functions::is_numeric_value(n) {
                as_int(&crate::functions::func(Func::Floor, vec![n.clone()]))?
            } else {
                return None;
            }
        }
    };
    let v = v.to_u64()?;
    Some(num(BigInt::from(nt::primepi(v))))
}

/// `n < 2` as a fuzzy boolean.
fn one_lt_two(n: &Expr) -> Option<bool> {
    if let Kind::Num(x) = n.kind() {
        return x.cmp_value(&Number::Int(2)).map(|o| o == std::cmp::Ordering::Less);
    }
    is_positive(&sub(&Expr::int(2), n))
}

/// Tribonacci polynomial coefficients (constant term first).
fn tribonacci_poly(n: usize) -> Vec<BigInt> {
    let mut seq: Vec<Vec<BigInt>> = vec![vec![], vec![BigInt::one()], vec![BigInt::zero(), BigInt::zero(), BigInt::one()]];
    while seq.len() <= n {
        let l = seq.len();
        let (a, b, c) = (&seq[l - 3], &seq[l - 2], &seq[l - 1]);
        let len = a.len().max(b.len() + 1).max(c.len() + 2);
        let mut r = vec![BigInt::zero(); len];
        for (i, x) in a.iter().enumerate() {
            r[i] += x;
        }
        for (i, x) in b.iter().enumerate() {
            r[i + 1] += x;
        }
        for (i, x) in c.iter().enumerate() {
            r[i + 2] += x;
        }
        seq.push(r);
    }
    seq.swap_remove(n)
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

fn eval_tribonacci(a: &[Expr]) -> Option<Expr> {
    if is_inf(&a[0]) {
        return Some(a[0].clone());
    }
    let n = as_int(&a[0])?.to_usize()?;
    if a.len() == 1 {
        // log2(T_n) < n log2(1.8393), the tribonacci constant
        if !super::integers::fits_exact(n as f64 * 0.879_146_788_819_652 + 2.0) {
            return None;
        }
        return Some(num(super::integers::tribonacci_number(n as u64)));
    }
    Some(poly_z(&tribonacci_poly(n), &a[1]))
}

/// Genocchi numbers `G_0..=G_n` as `2*(1 - 2**k)*B_k` (SymPy's convention,
/// `B_1 = +1/2`).
fn genocchi_numbers(n: usize) -> Option<Vec<BigInt>> {
    let mut out = Vec::with_capacity(n + 1);
    for k in 0..=n {
        let b = sf(SpecialFn::Bernoulli, vec![Expr::int(k as i64)]);
        let g = mul(&[Expr::int(2), sub(&Expr::one(), &pow(&Expr::int(2), &Expr::int(k as i64))), b]);
        out.push(as_int(&g)?);
    }
    Some(out)
}

fn eval_genocchi(a: &[Expr]) -> Option<Expr> {
    let n = &a[0];
    if a.len() == 2 && is_one(&a[1]) {
        return Some(sf(SpecialFn::Genocchi, vec![n.clone()]));
    }
    if is_integer(n) == Some(false) || is_nonnegative(n) == Some(false) {
        return None;
    }
    if a.len() == 1 {
        if crate::assumptions::is_odd(n) == Some(true) && is_positive(&sub(n, &Expr::one())) == Some(true) {
            return Some(Expr::zero());
        }
        if let Kind::Num(_) = n.kind() {
            let b = sf(SpecialFn::Bernoulli, vec![n.clone()]);
            return Some(mul(&[Expr::int(2), sub(&Expr::one(), &pow(&Expr::int(2), n)), b]));
        }
        return None;
    }
    let nv = as_int(n)?.to_usize()?;
    // SymPy's genocchi_poly is -n*euler_poly(n - 1, x), the Appell sequence
    // G_n(x) = sum_k binomial(n, k) G_k(0) x**(n - k) whose constants are
    // G_k(0) = -G_k except G_1(0) = G_1 = -1 (the number convention takes
    // B_1 = +1/2, the polynomial one B_1 = -1/2).
    let g = genocchi_numbers(nv)?;
    let row = nt::misc::binomial_row(nv as u64);
    let mut c = vec![BigInt::zero(); nv + 1];
    for k in 0..=nv {
        let constant = if k == 1 { g[k].clone() } else { -&g[k] };
        c[nv - k] = &row[k] * constant;
    }
    Some(poly_z(&c, &a[1]))
}

fn eval_andre(n: &Expr) -> Option<Expr> {
    if n.is_constant(Constant::NaN) || is_inf(n) {
        return Some(n.clone());
    }
    let v = as_int(n)?;
    if v.is_zero() {
        return Some(Expr::one());
    }
    if v == BigInt::from(-1) {
        return Some(crate::ops::neg(&crate::functions::log(&Expr::int(2))));
    }
    if v == BigInt::from(-2) {
        // -2*Catalan: the Catalan constant is not available natively
        return None;
    }
    if !v.is_negative() && v.is_even() {
        let e = sf(SpecialFn::Euler, vec![n.clone()]);
        return Some(crate::functions::func(Func::Abs, vec![e]));
    }
    if v.is_odd() {
        let m: BigInt = -&v - 1u32;
        // I**m with m even, times Rational(1 - 2**m, 4**m), times zeta(-n)
        let sign = if (&m / 2u32).is_odd() { -1 } else { 1 };
        let mm = m.to_i64()?;
        let r = if mm >= 0 {
            let two_m = BigInt::one() << (mm as u64);
            Expr::num(Number::from_ratio(num_rational::BigRational::new(BigInt::one() - &two_m, &two_m * &two_m)))
        } else {
            // (1 - 2**m) / 4**m = 4**|m| - 2**|m|
            let k = (-mm) as u64;
            num((BigInt::one() << (2 * k)) - (BigInt::one() << k))
        };
        let z = sf(SpecialFn::Zeta, vec![num(-v)]);
        return Some(mul(&[Expr::int(sign), r, z]));
    }
    None
}

/// SymPy's `motzkin(n)` value (`_motzkin(n - 1)`, with `motzkin(0) == 1`).
pub fn motzkin_number(n: u64) -> BigInt {
    if n <= 2 {
        return BigInt::one();
    }
    let (mut m1, mut m2) = (BigInt::one(), BigInt::one());
    for i in 2..n {
        let i = BigInt::from(i);
        let m = ((2u32 * &i + 1u32) * &m2 + (3u32 * &i - 3u32) * &m1) / (&i + 2u32);
        m1 = m2;
        m2 = m;
    }
    m2
}

fn genocchi_plain(a: &[Expr]) -> Option<&Expr> {
    if a.len() > 1 && !is_one(&a[1]) {
        return None;
    }
    let n = &a[0];
    (is_integer(n) == Some(true) && is_nonnegative(n) == Some(true)).then_some(n)
}

/// Class `is_integer` of a number theoretic head.
pub fn integer_fact(s: SpecialFn, a: &[Expr]) -> Option<bool> {
    use SpecialFn::*;
    match s {
        Totient | ReducedTotient | DivisorSigma | UdivisorSigma | PrimeNu | PrimeOmega | PrimePi
        | LegendreSymbol | JacobiSymbol | KroneckerSymbol | Mobius => Some(true),
        Andre => (is_integer(&a[0]) == Some(true) && is_nonnegative(&a[0]) == Some(true)).then_some(true),
        Genocchi => genocchi_plain(a).map(|_| true),
        _ => None,
    }
}

/// Whether the head is known not prime (the symbol heads and `mobius`).
pub fn not_prime(s: SpecialFn) -> bool {
    use SpecialFn::*;
    matches!(s, LegendreSymbol | JacobiSymbol | KroneckerSymbol | Mobius)
}

/// Class sign facts as `(positive, negative, nonnegative)`, each a fuzzy
/// boolean; `None` when the head is not known to be real.
pub fn sign_facts(s: SpecialFn, a: &[Expr]) -> Option<(Option<bool>, Option<bool>, bool)> {
    use SpecialFn::*;
    Some(match s {
        Totient | ReducedTotient | DivisorSigma | UdivisorSigma => (Some(true), Some(false), true),
        PrimeNu | PrimeOmega | PrimePi => (None, Some(false), true),
        LegendreSymbol | JacobiSymbol | KroneckerSymbol | Mobius => (None, None, false),
        Andre => {
            let n = &a[0];
            if is_nonnegative(n) == Some(true) {
                (Some(true), Some(false), true)
            } else if integer_fact(s, a) == Some(true) {
                (None, None, false)
            } else {
                return None;
            }
        }
        Genocchi => {
            let n = genocchi_plain(a)?;
            let odd = crate::assumptions::is_odd(n);
            let half = crate::ops::div(n, &Expr::int(2));
            let pos = if is_zero(n) == Some(true) || odd == Some(true) {
                Some(false)
            } else {
                crate::assumptions::is_even(&half)
            };
            let neg = if odd == Some(true) {
                is_positive(&sub(n, &Expr::one())).map(|b| !b)
            } else {
                crate::assumptions::is_odd(&half)
            };
            (pos, neg, false)
        }
        _ => return None,
    })
}
