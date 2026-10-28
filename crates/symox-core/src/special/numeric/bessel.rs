//! Bessel, Hankel and Airy functions in double precision.
//!
//! Real order and real argument: power series in double-double arithmetic
//! for moderate arguments, Hankel asymptotic expansions (with stable
//! recurrences in the order) for large arguments, and the integral
//! `K_nu(x) = int_0^oo exp(-x cosh t) cosh(nu t) dt` evaluated by the
//! (exponentially convergent) trapezoidal rule.

use num_complex::Complex64 as C64;

use super::dd::{consts, DD};
use super::gamma::{gamma_real, rgamma_real};
use super::{cospi, cr, is_real_c, sinpi, PI};

const EULER: f64 = crate::expr::EULER_GAMMA_F64;
const ASYM_X: f64 = 25.0;

fn is_int(v: f64) -> bool {
    v == v.floor()
}

/// `sum_k (-x^2/4)^k / (k! (nu+1)_k)` (sign = -1) or with `+` (I).
fn series_sum_dd(nu: f64, x: f64, sign: f64) -> DD {
    let t = (DD::new(x) * DD::new(x)).mul_f64(0.25 * sign);
    let mut term = DD::ONE;
    let mut sum = DD::ONE;
    let mut k = 0.0;
    loop {
        k += 1.0;
        let d = DD::new(nu).add_f64(k).mul_f64(k);
        term = term * t / d;
        sum = sum + term;
        if (term.hi.abs() < 1e-33 * sum.hi.abs() && k > x) || k > 2000.0 {
            break;
        }
    }
    sum
}

/// `(x/2)^nu / Gamma(nu+1)`.
fn prefactor(nu: f64, x: f64) -> f64 {
    if nu + 1.0 > 0.0 && nu > 100.0 {
        // log space
        let l = DD::new(x / 2.0).ln().mul_f64(nu) - super::gamma::lngamma_dd(DD::new(nu + 1.0));
        return l.exp().to_f64();
    }
    (x / 2.0).powf(nu) * rgamma_real(nu + 1.0)
}

fn j_series(nu: f64, x: f64) -> DD {
    series_sum_dd(nu, x, -1.0).mul_f64(prefactor(nu, x))
}

fn i_series(nu: f64, x: f64) -> DD {
    series_sum_dd(nu, x, 1.0).mul_f64(prefactor(nu, x))
}

/// Hankel asymptotic `(P, Q)` for order `nu` at `x`.
fn hankel_pq(nu: f64, x: f64) -> (f64, f64) {
    let mu = 4.0 * nu * nu;
    let mut p = 1.0;
    let mut q = 0.0;
    let mut a = 1.0; // a_k / x^k
    let mut last = f64::INFINITY;
    for k in 1..200 {
        let kf = k as f64;
        a *= (mu - (2.0 * kf - 1.0).powi(2)) / (kf * 8.0 * x);
        if a.abs() > last && k > 2 {
            break;
        }
        last = a.abs();
        match k % 4 {
            1 => q += a,
            2 => p -= a,
            3 => q -= a,
            _ => p += a,
        }
        if a.abs() < 1e-17 {
            break;
        }
    }
    (p, q)
}

/// `cos` and `sin` of `x - (nu/2 + 1/4) pi`, reduced in double-double.
fn chi_cos_sin(nu: f64, x: f64) -> (f64, f64) {
    let pi = consts().pi;
    let chi = DD::new(x) - pi * DD::new(nu / 2.0 + 0.25);
    let half_pi = pi.mul_f64(0.5);
    let k = (chi.hi / half_pi.hi).round();
    let r = (chi - half_pi.mul_f64(k)).to_f64();
    let (s, c) = r.sin_cos();
    match (k as i64).rem_euclid(4) {
        0 => (c, s),
        1 => (-s, c),
        2 => (-c, -s),
        _ => (s, -c),
    }
}

/// `(J_nu(x), Y_nu(x))` from the asymptotic expansion.
fn jy_asym(nu: f64, x: f64) -> (f64, f64) {
    let (p, q) = hankel_pq(nu, x);
    let (c, s) = chi_cos_sin(nu, x);
    let f = (2.0 / (PI * x)).sqrt();
    (f * (p * c - q * s), f * (p * s + q * c))
}

fn asym_ok(nu: f64, x: f64) -> bool {
    x >= ASYM_X && nu * nu < 0.5 * x
}

/// Bessel J, real order, real `x > 0`.
pub fn besselj_real(nu: f64, x: f64) -> f64 {
    if x == 0.0 {
        return if nu == 0.0 { 1.0 } else if nu > 0.0 || is_int(nu) { 0.0 } else { f64::INFINITY };
    }
    if nu < 0.0 && is_int(nu) {
        let v = besselj_real(-nu, x);
        return if (nu as i64) % 2 == 0 { v } else { -v };
    }
    if x >= ASYM_X {
        if asym_ok(nu, x) {
            return jy_asym(nu, x).0;
        }
        if nu > 0.0 && nu < x {
            // forward recurrence from nu0 in [0, 1)
            let nu0 = nu - nu.floor();
            let mut a = jy_asym(nu0, x).0;
            let mut b = jy_asym(nu0 + 1.0, x).0;
            let mut k = nu0 + 1.0;
            while k < nu - 0.5 {
                let c = 2.0 * k / x * b - a;
                a = b;
                b = c;
                k += 1.0;
            }
            return if nu0 == nu { a } else { b };
        }
    }
    j_series(nu, x).to_f64()
}

/// Y_n for integer `n >= 0` by its series (double-double), `x < ASYM_X`.
fn yn_series(n: i64, x: f64) -> f64 {
    let cst = consts();
    let h = x / 2.0;
    // first finite sum: (n-k-1)!/k! (x/2)^(2k-n)
    let mut s1 = DD::ZERO;
    if n > 0 {
        let mut f = gamma_real(n as f64); // (n-1)!/0!
        let mut p = h.powi(-(n as i32));
        for k in 0..n {
            s1 = s1 + DD::new(f).mul_f64(p);
            if k + 1 < n {
                f = f / ((n - k - 1) as f64 * (k + 1) as f64);
            }
            p *= h * h;
        }
    }
    // infinite sum: [psi(k+1) + psi(n+k+1)] (-x^2/4)^k (x/2)^n / (k! (n+k)!)
    let t = -(h * h);
    let mut psi1 = -cst.euler; // psi(k+1)
    let mut psi2 = -cst.euler; // psi(n+k+1)
    for j in 1..=n {
        psi2 = psi2 + DD::new(j as f64).recip();
    }
    let mut term = DD::new(h.powi(n as i32) * rgamma_real(n as f64 + 1.0));
    let mut s2 = DD::ZERO;
    let mut k = 0i64;
    loop {
        s2 = s2 + term * (psi1 + psi2);
        k += 1;
        term = term.mul_f64(t).div_f64((k * (n + k)) as f64);
        psi1 = psi1 + DD::new(k as f64).recip();
        psi2 = psi2 + DD::new((n + k) as f64).recip();
        if (term.hi.abs() < 1e-34 && k as f64 > x) || k > 3000 {
            break;
        }
    }
    let jn = j_series(n as f64, x);
    let r = jn.mul_f64(2.0) * DD::new(h).ln() - s1 - s2;
    (r / cst.pi).to_f64()
}

/// Y for non integer order by reflection (double-double).
fn y_reflect(nu: f64, x: f64) -> f64 {
    let j1 = j_series(nu, x);
    let j2 = j_series(-nu, x);
    ((j1.mul_f64(cospi(nu)) - j2) / DD::new(sinpi(nu))).to_f64()
}

/// Bessel Y, real order, real `x > 0`.
pub fn bessely_real(nu: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if is_int(nu) && nu < 0.0 {
        let v = bessely_real(-nu, x);
        return if (nu as i64) % 2 == 0 { v } else { -v };
    }
    if x >= ASYM_X && asym_ok(nu, x) {
        return jy_asym(nu, x).1;
    }
    if nu < 0.0 {
        return y_reflect(nu, x);
    }
    let nu0 = nu - nu.floor();
    let y0 = |v: f64| -> f64 {
        if x >= ASYM_X {
            jy_asym(v, x).1
        } else if is_int(v) {
            yn_series(v as i64, x)
        } else {
            y_reflect(v, x)
        }
    };
    if nu < 1.5 || (x < ASYM_X && !is_int(nu) && nu < 4.0) {
        return y0(nu);
    }
    if is_int(nu) && x < ASYM_X && nu < 3.0 {
        return y0(nu);
    }
    // forward recurrence (stable for Y)
    let mut a = y0(nu0);
    let mut b = y0(nu0 + 1.0);
    let mut k = nu0 + 1.0;
    while k < nu - 0.5 {
        let c = 2.0 * k / x * b - a;
        a = b;
        b = c;
        k += 1.0;
    }
    b
}

/// Bessel I, real order, real `x >= 0`.
pub fn besseli_real(nu: f64, x: f64) -> f64 {
    if x == 0.0 {
        return if nu == 0.0 { 1.0 } else if nu > 0.0 || is_int(nu) { 0.0 } else { f64::INFINITY };
    }
    if x < 0.0 {
        if is_int(nu) {
            let v = besseli_real(nu, -x);
            return if (nu as i64) % 2 == 0 { v } else { -v };
        }
        return f64::NAN;
    }
    if is_int(nu) && nu < 0.0 {
        return besseli_real(-nu, x);
    }
    if x > 30.0 && nu * nu < 0.5 * x {
        // e^x / sqrt(2 pi x) sum (-1)^k a_k / x^k
        let mu = 4.0 * nu * nu;
        let mut s = 1.0;
        let mut a = 1.0;
        let mut last = f64::INFINITY;
        for k in 1..200 {
            let kf = k as f64;
            a *= -(mu - (2.0 * kf - 1.0).powi(2)) / (kf * 8.0 * x);
            if a.abs() > last && k > 2 {
                break;
            }
            last = a.abs();
            s += a;
            if a.abs() < 1e-17 {
                break;
            }
        }
        let e = DD::new(x).exp().to_f64();
        return e / (2.0 * PI * x).sqrt() * s;
    }
    i_series(nu, x).to_f64()
}

/// Bessel K by the trapezoidal rule on `int_0^oo exp(-z cosh t) cosh(nu t) dt`.
pub fn besselk_c(nu: f64, z: C64) -> Option<C64> {
    if z.re <= 0.0 {
        return None;
    }
    let nu = nu.abs();
    let ax = z.norm();
    let dstrip = (z.re / z.im.abs().max(1e-300)).atan().min(PI / 2.0);
    let h = (0.5 / ((ax * ax + nu * nu).sqrt() + 1.0).sqrt()).min(0.2) * (dstrip / (PI / 2.0));
    let real = z.im == 0.0;
    // peak of the integrand exponent -x cosh t + nu t
    let mut sum = C64::new(0.0, 0.0);
    let mut sum_r = DD::ZERO;
    let f = |t: f64| -> C64 {
        let ch = t.cosh();
        let e1 = -z * ch + nu * t;
        let e2 = -z * ch - nu * t;
        (e1.exp() + e2.exp()) * 0.5
    };
    let fr = |t: f64| -> f64 {
        let ch = t.cosh();
        let x = z.re;
        ((-x * ch + nu * t).exp() + (-x * ch - nu * t).exp()) * 0.5
    };
    let tpeak = if nu > 0.0 { (nu / z.re).asinh() } else { 0.0 };
    let mut k = 0usize;
    loop {
        let t = k as f64 * h;
        let w = if k == 0 { 0.5 } else { 1.0 };
        if real {
            let v = fr(t) * w;
            sum_r = sum_r.add_f64(v);
            if t > tpeak && v < 1e-19 * sum_r.hi {
                break;
            }
        } else {
            let v = f(t) * w;
            sum += v;
            if t > tpeak && v.norm() < 1e-19 * sum.norm() {
                break;
            }
        }
        k += 1;
        if k > 200_000 {
            break;
        }
    }
    Some(if real { cr(sum_r.to_f64() * h) } else { sum * h })
}

pub fn besselk_real(nu: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    if x > 705.0 {
        return 0.0;
    }
    besselk_c(nu, cr(x)).unwrap().re
}

// ----------------------------------------------------------------------
// complex arguments (power series, moderate |z|)
// ----------------------------------------------------------------------

fn cpow(b: C64, e: f64) -> C64 {
    if b.im == 0.0 && b.re > 0.0 {
        return cr(b.re.powf(e));
    }
    (b.ln() * e).exp()
}

fn ji_series_c(nu: f64, z: C64, sign: f64) -> Option<C64> {
    if z.norm() > 40.0 {
        return None;
    }
    let t = z * z * (0.25 * sign);
    let mut term = C64::new(1.0, 0.0);
    let mut sum = term;
    for k in 1..3000 {
        let kf = k as f64;
        term = term * t / (kf * (nu + kf));
        sum += term;
        if term.norm() < 1e-17 * sum.norm() && kf > z.norm() {
            break;
        }
    }
    Some(sum * cpow(z / 2.0, nu) * rgamma_real(nu + 1.0))
}

pub fn besselj_c(nu: C64, z: C64) -> Option<C64> {
    if !is_real_c(nu) {
        return None;
    }
    let n = nu.re;
    if is_real_c(z) && z.re >= 0.0 {
        return Some(cr(besselj_real(n, z.re)));
    }
    if is_real_c(z) && is_int(n) {
        let v = besselj_real(n, -z.re);
        return Some(cr(if (n as i64) % 2 == 0 { v } else { -v }));
    }
    if n < 0.0 && is_int(n) {
        let v = besselj_c(cr(-n), z)?;
        return Some(if (n as i64) % 2 == 0 { v } else { -v });
    }
    if z.norm() > 30.0 && z.re > 0.0 && n * n < 0.5 * z.norm() {
        let (p, q) = hankel_pq_c(n, z);
        let chi = z - (n / 2.0 + 0.25) * PI;
        let f = (2.0 / (PI * z)).sqrt();
        return Some(f * (p * chi.cos() - q * chi.sin()));
    }
    ji_series_c(n, z, -1.0)
}

fn hankel_pq_c(nu: f64, z: C64) -> (C64, C64) {
    let mu = 4.0 * nu * nu;
    let mut p = C64::new(1.0, 0.0);
    let mut q = C64::new(0.0, 0.0);
    let mut a = C64::new(1.0, 0.0);
    let mut last = f64::INFINITY;
    for k in 1..200 {
        let kf = k as f64;
        a = a * (mu - (2.0 * kf - 1.0).powi(2)) / (z * (kf * 8.0));
        if a.norm() > last && k > 2 {
            break;
        }
        last = a.norm();
        match k % 4 {
            1 => q += a,
            2 => p -= a,
            3 => q -= a,
            _ => p += a,
        }
        if a.norm() < 1e-17 {
            break;
        }
    }
    (p, q)
}

pub fn besseli_c(nu: C64, z: C64) -> Option<C64> {
    if !is_real_c(nu) {
        return None;
    }
    let n = nu.re;
    if is_real_c(z) && (z.re >= 0.0 || is_int(n)) {
        return Some(cr(besseli_real(n, z.re)));
    }
    if n < 0.0 && is_int(n) {
        return besseli_c(cr(-n), z);
    }
    ji_series_c(n, z, 1.0)
}

pub fn bessely_c(nu: C64, z: C64) -> Option<C64> {
    if !is_real_c(nu) {
        return None;
    }
    let n = nu.re;
    if is_real_c(z) && z.re > 0.0 {
        return Some(cr(bessely_real(n, z.re)));
    }
    if z == cr(0.0) {
        return None;
    }
    if is_int(n) {
        if n < 0.0 {
            let v = bessely_c(cr(-n), z)?;
            return Some(if (n as i64) % 2 == 0 { v } else { -v });
        }
        // integer order: series with logarithm (complex double)
        let ni = n as i64;
        if z.norm() > 30.0 {
            return None;
        }
        let h = z / 2.0;
        let mut s1 = C64::new(0.0, 0.0);
        if ni > 0 {
            let mut f = gamma_real(n);
            let mut p = h.powi(-(ni as i32));
            for k in 0..ni {
                s1 += p * f;
                if k + 1 < ni {
                    f /= (ni - k - 1) as f64 * (k + 1) as f64;
                }
                p *= h * h;
            }
        }
        let t = -(h * h);
        let mut psi1 = -EULER;
        let mut psi2 = -EULER;
        for j in 1..=ni {
            psi2 += 1.0 / j as f64;
        }
        let mut term = h.powi(ni as i32) * rgamma_real(n + 1.0);
        let mut s2 = C64::new(0.0, 0.0);
        for k in 1..3000i64 {
            s2 += term * (psi1 + psi2);
            term = term * t / (k * (ni + k)) as f64;
            psi1 += 1.0 / k as f64;
            psi2 += 1.0 / (ni + k) as f64;
            if term.norm() < 1e-18 * s2.norm() && k as f64 > z.norm() {
                break;
            }
        }
        let jn = besselj_c(nu, z)?;
        return Some((jn * 2.0 * h.ln() - s1 - s2) / PI);
    }
    let j1 = besselj_c(nu, z)?;
    let j2 = besselj_c(cr(-n), z)?;
    Some((j1 * cospi(n) - j2) / sinpi(n))
}

pub fn besselk_cc(nu: C64, z: C64) -> Option<C64> {
    if !is_real_c(nu) {
        return None;
    }
    if is_real_c(z) && z.re > 0.0 {
        return Some(cr(besselk_real(nu.re, z.re)));
    }
    if z.re > 0.0 && ((z.re >= z.im.abs() && z.norm() > 2.0) || z.norm() > 20.0) {
        return besselk_c(nu.re, z);
    }
    if z == cr(0.0) || z.norm() > 30.0 {
        return None;
    }
    let n = nu.re.abs();
    if !is_int(n) {
        // pi/2 (I_{-nu} - I_nu) / sin(nu pi)
        let a = ji_series_c(-n, z, 1.0)?;
        let b = ji_series_c(n, z, 1.0)?;
        return Some((a - b) * (PI / 2.0) / sinpi(n));
    }
    // integer order: series with logarithm
    let ni = n as i64;
    let h = z / 2.0;
    let mut s1 = C64::new(0.0, 0.0);
    if ni > 0 {
        let mut f = gamma_real(n);
        let mut p = h.powi(-(ni as i32));
        for k in 0..ni {
            let t = p * f;
            s1 += if k % 2 == 0 { t } else { -t };
            if k + 1 < ni {
                f /= (ni - k - 1) as f64 * (k + 1) as f64;
            }
            p *= h * h;
        }
    }
    let t = h * h;
    let mut psi1 = -EULER;
    let mut psi2 = -EULER;
    for j in 1..=ni {
        psi2 += 1.0 / j as f64;
    }
    let mut term = h.powi(ni as i32) * rgamma_real(n + 1.0);
    let mut s2 = C64::new(0.0, 0.0);
    for k in 1..3000i64 {
        s2 += term * (psi1 + psi2);
        term = term * t / (k * (ni + k)) as f64;
        psi1 += 1.0 / k as f64;
        psi2 += 1.0 / (ni + k) as f64;
        if term.norm() < 1e-18 * s2.norm() && k as f64 > z.norm() {
            break;
        }
    }
    let i_n = ji_series_c(n, z, 1.0)?;
    let sgn = if ni % 2 == 0 { 1.0 } else { -1.0 };
    Some(-sgn * h.ln() * i_n + s1 * 0.5 + s2 * (sgn * 0.5))
}

pub fn hankel_c(nu: C64, z: C64, first: bool) -> Option<C64> {
    let j = besselj_c(nu, z)?;
    let y = bessely_c(nu, z)?;
    let iy = C64::new(-y.im, y.re);
    Some(if first { j + iy } else { j - iy })
}

/// Spherical Bessel functions `j_n`, `y_n`.
pub fn spherical_c(n: C64, z: C64, first: bool) -> Option<C64> {
    if z == cr(0.0) {
        return if first && n == cr(0.0) { Some(cr(1.0)) } else { None };
    }
    let f = (PI / (2.0 * z)).sqrt();
    let nu = n + 0.5;
    Some(f * if first { besselj_c(nu, z)? } else { bessely_c(nu, z)? })
}

// ----------------------------------------------------------------------
// Airy functions
// ----------------------------------------------------------------------

/// Maclaurin series of `(Ai, Bi, Ai', Bi')`.
fn airy_series(z: C64) -> [C64; 4] {
    // f = sum 3^k (1/3)_k z^(3k) / (3k)!,  g = sum 3^k (2/3)_k z^(3k+1) / (3k+1)!
    let z3 = z * z * z;
    let mut f = C64::new(1.0, 0.0);
    let mut g = z;
    let mut tf = C64::new(1.0, 0.0);
    let mut tg = z;
    // f' = sum_{k>=1} 3^k (1/3)_k z^(3k-1)/(3k-1)!, g' = sum_k 3^k (2/3)_k z^(3k)/(3k)!
    let mut tfp = z * z / 2.0;
    let mut fp = tfp;
    let mut tgp = C64::new(1.0, 0.0);
    let mut gp = tgp;
    for k in 1..400 {
        let kf = k as f64;
        tf = tf * z3 / ((3.0 * kf - 1.0) * (3.0 * kf));
        tg = tg * z3 / ((3.0 * kf) * (3.0 * kf + 1.0));
        f += tf;
        g += tg;
        if k >= 2 {
            tfp = tfp * z3 / ((3.0 * kf - 1.0) * (3.0 * kf - 3.0));
            fp += tfp;
        }
        tgp = tgp * z3 / ((3.0 * kf) * (3.0 * kf - 2.0));
        gp += tgp;
        if tf.norm() + tg.norm() + tfp.norm() + tgp.norm() < 1e-18 * (f.norm() + g.norm()) && kf > 3.0 {
            break;
        }
    }
    let c1 = 1.0 / (3f64.powf(2.0 / 3.0) * gamma_real(2.0 / 3.0)); // Ai(0)
    let c2 = 1.0 / (3f64.powf(1.0 / 3.0) * gamma_real(1.0 / 3.0)); // -Ai'(0)
    let s3 = 3f64.sqrt();
    [
        f * c1 - g * c2,
        (f * c1 + g * c2) * s3,
        fp * c1 - gp * c2,
        (fp * c1 + gp * c2) * s3,
    ]
}

/// `which`: 0 = Ai, 1 = Bi, 2 = Ai', 3 = Bi'.
pub fn airy_c(which: usize, z: C64) -> Option<C64> {
    if is_real_c(z) {
        return Some(cr(airy_real(which, z.re)));
    }
    if z.norm() > 8.0 {
        return None;
    }
    Some(airy_series(z)[which])
}

fn airy_real(which: usize, x: f64) -> f64 {
    if x.abs() <= 1.0 {
        return airy_series(cr(x))[which].re;
    }
    let s3 = 3f64.sqrt();
    if x > 0.0 {
        let zeta = 2.0 / 3.0 * x * x.sqrt();
        match which {
            0 => (x / 3.0).sqrt() / PI * besselk_real(1.0 / 3.0, zeta),
            2 => -x / (PI * s3) * besselk_real(2.0 / 3.0, zeta),
            1 => (x / 3.0).sqrt() * (besseli_real(-1.0 / 3.0, zeta) + besseli_real(1.0 / 3.0, zeta)),
            _ => x / s3 * (besseli_real(-2.0 / 3.0, zeta) + besseli_real(2.0 / 3.0, zeta)),
        }
    } else {
        let t = -x;
        let zeta = 2.0 / 3.0 * t * t.sqrt();
        match which {
            0 => t.sqrt() / 3.0 * (besselj_real(1.0 / 3.0, zeta) + besselj_real(-1.0 / 3.0, zeta)),
            1 => (t / 3.0).sqrt() * (besselj_real(-1.0 / 3.0, zeta) - besselj_real(1.0 / 3.0, zeta)),
            2 => t / 3.0 * (besselj_real(2.0 / 3.0, zeta) - besselj_real(-2.0 / 3.0, zeta)),
            _ => t / s3 * (besselj_real(-2.0 / 3.0, zeta) + besselj_real(2.0 / 3.0, zeta)),
        }
    }
}
