//! Gamma, polygamma, zeta and related functions in double precision.

use std::sync::OnceLock;

use num_complex::Complex64 as C64;
use num_traits::ToPrimitive;

use super::dd::{consts, DD};
use super::{c, cr, is_real_c, sinpi, cospi, PI};

/// `B_{2k}` for k = 0.. as f64.
pub(super) fn bern2() -> &'static [f64] {
    static B: OnceLock<Vec<f64>> = OnceLock::new();
    B.get_or_init(|| {
        (0..40)
            .map(|k| {
                let b = crate::special::integers::bernoulli_number(2 * k);
                b.numer().to_f64().unwrap() / b.denom().to_f64().unwrap()
            })
            .collect()
    })
}

/// `B_{2k} / (2k (2k-1))` (Stirling series coefficients), k >= 1.
fn stirling_coeffs() -> &'static [f64] {
    static S: OnceLock<Vec<f64>> = OnceLock::new();
    S.get_or_init(|| {
        let b = bern2();
        (1..16).map(|k| b[k] / ((2 * k) as f64 * (2 * k - 1) as f64)).collect()
    })
}

const SHIFT: f64 = 20.0;

/// `ln Gamma(x)` for real `x >= 0.5` in double-double precision.
pub(super) fn lngamma_dd(x: DD) -> DD {
    let mut y = x;
    let mut prod = DD::ONE;
    while y.hi < SHIFT {
        prod = prod * y;
        y = y.add_f64(1.0);
    }
    let cst = consts();
    let ly = y.ln();
    let inv = 1.0 / y.hi;
    let inv2 = inv * inv;
    let mut s = 0.0;
    let mut p = inv;
    for c in stirling_coeffs().iter().take(9) {
        s += c * p;
        p *= inv2;
    }
    let main = y.add_f64(-0.5) * ly - y + cst.half_ln_2pi;
    let mut r = main.add_f64(s);
    if prod.hi != 1.0 || prod.lo != 0.0 {
        r = r - prod.ln();
    }
    r
}

/// Gamma function of a real argument (DD internally).
pub fn gamma_real(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == x.floor() {
        if x <= 0.0 {
            return f64::NAN;
        }
        if x <= 171.0 {
            let mut acc = DD::ONE;
            let mut i = 2.0;
            while i < x {
                acc = acc.mul_f64(i);
                i += 1.0;
            }
            return acc.to_f64();
        }
        return f64::INFINITY;
    }
    if x > 171.7 {
        return f64::INFINITY;
    }
    if x < 0.5 {
        // reflection: pi / (sin(pi x) Gamma(1 - x))
        if x < -190.0 {
            return 0.0;
        }
        let omx = DD::ONE - DD::new(x);
        let lg = lngamma_dd(omx);
        let s = sinpi(x);
        // pi / (s * exp(lg)): do it in log space for the magnitude
        let g = lg.exp();
        if g.hi.is_infinite() {
            return 0.0;
        }
        let d = g.mul_f64(s);
        return (consts().pi / d).to_f64();
    }
    lngamma_dd(DD::new(x)).exp().to_f64()
}

/// `ln Gamma(z)` principal branch for complex `z` (real negative axis
/// gives `ln|Gamma(x)| + i pi floor(x)`).
pub fn loggamma_c(z: C64) -> Option<C64> {
    if z.im == 0.0 {
        let x = z.re;
        if x <= 0.0 && x == x.floor() {
            return None;
        }
        if x > 0.0 {
            return Some(cr(lngamma_dd(DD::new(x)).to_f64()));
        }
        // ln|Gamma(x)| = ln pi - ln|sin pi x| - lnGamma(1-x)
        let omx = DD::ONE - DD::new(x);
        let lg = lngamma_dd(omx);
        let re = (consts().pi.ln() - DD::new(sinpi(x).abs()).ln() - lg).to_f64();
        return Some(c(re, PI * x.floor()));
    }
    // shift into the Stirling region with principal logs
    let mut w = z;
    let mut logs = C64::new(0.0, 0.0);
    while w.re < SHIFT || w.norm() < SHIFT {
        if w.re >= SHIFT && w.norm() >= SHIFT {
            break;
        }
        logs += w.ln();
        w += 1.0;
        if w.re > 1e4 {
            break;
        }
    }
    Some(stirling_c(w) - logs)
}

fn stirling_c(w: C64) -> C64 {
    let inv = w.inv();
    let inv2 = inv * inv;
    let mut s = C64::new(0.0, 0.0);
    let mut p = inv;
    for cf in stirling_coeffs().iter().take(12) {
        s += p * *cf;
        p *= inv2;
    }
    (w - 0.5) * w.ln() - w + consts().half_ln_2pi.to_f64() + s
}

/// Gamma function, complex argument.
pub fn gamma_c(z: C64) -> Option<C64> {
    if is_real_c(z) {
        let v = gamma_real(z.re);
        return if v.is_nan() { None } else { Some(cr(v)) };
    }
    if z.re < 0.5 {
        let g = gamma_c(1.0 - z)?;
        return Some(PI / ((z * PI).sin() * g));
    }
    let mut w = z;
    let mut prod = C64::new(1.0, 0.0);
    while w.norm() < SHIFT {
        prod *= w;
        w += 1.0;
    }
    Some(stirling_c(w).exp() / prod)
}

/// `1 / Gamma(x)` (zero at the poles).
pub fn rgamma_real(x: f64) -> f64 {
    if x <= 0.0 && x == x.floor() {
        return 0.0;
    }
    1.0 / gamma_real(x)
}

/// Digamma of a real argument.
pub fn digamma_real(x: f64) -> f64 {
    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }
    if x < 0.5 {
        // psi(x) = psi(1-x) - pi cot(pi x)
        let r = digamma_dd(DD::ONE - DD::new(x));
        let cot = cospi(x) / sinpi(x);
        return (r - DD::new(PI).mul_f64(cot)).to_f64();
    }
    digamma_dd(DD::new(x)).to_f64()
}

fn digamma_dd(x: DD) -> DD {
    let mut y = x;
    let mut acc = DD::ZERO;
    while y.hi < SHIFT {
        acc = acc + y.recip();
        y = y.add_f64(1.0);
    }
    let b = bern2();
    let inv = 1.0 / y.hi;
    let inv2 = inv * inv;
    let mut s = 0.0;
    let mut p = inv2;
    for k in 1..10 {
        s += b[k] / (2 * k) as f64 * p;
        p *= inv2;
    }
    y.ln() - y.recip().mul_f64(0.5) - DD::new(s) - acc
}

/// Digamma, complex argument.
pub fn digamma_c(z: C64) -> Option<C64> {
    if is_real_c(z) {
        let v = digamma_real(z.re);
        return if v.is_nan() { None } else { Some(cr(v)) };
    }
    if z.re < 0.5 {
        let r = digamma_c(1.0 - z)?;
        let pz = z * PI;
        return Some(r - PI * pz.cos() / pz.sin());
    }
    let mut w = z;
    let mut acc = C64::new(0.0, 0.0);
    while w.norm() < SHIFT {
        acc += w.inv();
        w += 1.0;
    }
    let b = bern2();
    let inv2 = (w * w).inv();
    let mut s = C64::new(0.0, 0.0);
    let mut p = inv2;
    for k in 1..12 {
        s += p * (b[k] / (2 * k) as f64);
        p *= inv2;
    }
    Some(w.ln() - 0.5 / w - s - acc)
}

/// `polygamma(m, z)`.
pub fn polygamma_c(m: C64, z: C64) -> Option<C64> {
    if !is_real_c(m) {
        return None;
    }
    let mr = m.re;
    if mr == 0.0 {
        return digamma_c(z);
    }
    if mr == -1.0 {
        return Some(loggamma_c(z)? - consts().half_ln_2pi.to_f64());
    }
    if mr > 0.0 && mr == mr.floor() {
        // (-1)^(m+1) m! zeta(m+1, z)
        let f = gamma_real(mr + 1.0);
        let zt = hurwitz_zeta(cr(mr + 1.0), z)?;
        let sgn = if (mr as i64) % 2 == 1 { 1.0 } else { -1.0 };
        return Some(zt * (f * sgn));
    }
    None
}

/// Euler–Maclaurin summation for the Hurwitz zeta function `zeta(s, a)`.
pub fn hurwitz_zeta(s: C64, a: C64) -> Option<C64> {
    if s == cr(1.0) {
        return None;
    }
    if is_real_c(a) && a.re <= 0.0 && a.re == a.re.floor() {
        return None;
    }
    let b = bern2();
    let sabs = s.norm();
    let n = ((18.0 + sabs * 0.6 - a.re).ceil()).max(0.0) as usize;
    let mut sum = C64::new(0.0, 0.0);
    let mut sum_dd_re = DD::ZERO;
    let sint = is_real_c(s) && s.re == s.re.floor() && s.re.abs() < 1000.0;
    let real = is_real_c(s) && is_real_c(a) && (a.re > 0.0 || sint);
    for k in 0..n {
        let base = a + k as f64;
        if real {
            let b = DD::new(a.re).add_f64(k as f64);
            let t = if sint { b.powi(-(s.re as i64)) } else { DD::new(base.re.powf(-s.re)) };
            sum_dd_re = sum_dd_re + t;
        } else {
            sum += cpow(base, -s);
        }
    }
    let an = a + n as f64;
    let mut tail = cpow(an, 1.0 - s) / (s - 1.0) + cpow(an, -s) * 0.5;
    // sum_j B_2j/(2j)! (s)_{2j-1} an^{-s-2j+1}
    let mut fact = 1.0; // (2j)!
    let mut poch = s; // (s)_{2j-1}
    let mut pw = cpow(an, -s - 1.0);
    let inv2 = (an * an).inv();
    for j in 1..30 {
        fact *= (2 * j - 1) as f64 * (2 * j) as f64;
        let t = poch * pw * (b[j] / fact);
        tail += t;
        if t.norm() < 1e-18 * tail.norm().max(1e-300) {
            break;
        }
        poch = poch * (s + (2 * j - 1) as f64) * (s + (2 * j) as f64);
        pw *= inv2;
    }
    if real {
        let v = (sum_dd_re + DD::new(tail.re)).to_f64();
        return Some(cr(v));
    }
    Some(sum + tail)
}

#[inline]
fn cpow(b: C64, e: C64) -> C64 {
    if b.im == 0.0 && b.re > 0.0 && e.im == 0.0 {
        return cr(b.re.powf(e.re));
    }
    (b.ln() * e).exp()
}

/// Riemann zeta function.
pub fn zeta_c(s: C64) -> Option<C64> {
    if s == cr(1.0) {
        return None;
    }
    if is_real_c(s) {
        let x = s.re;
        if x == 0.0 {
            return Some(cr(-0.5));
        }
        if x < 0.0 && x == x.floor() && (x as i64) % 2 == 0 {
            return Some(cr(0.0));
        }
    }
    if s.re < 0.5 {
        // zeta(s) = 2^s pi^(s-1) sin(pi s / 2) Gamma(1-s) zeta(1-s)
        let oms = 1.0 - s;
        let z1 = zeta_c(oms)?;
        let g = gamma_c(oms)?;
        let sn = if is_real_c(s) { cr(sinpi(s.re / 2.0)) } else { (s * (PI / 2.0)).sin() };
        let f = cpow(cr(2.0), s) * cpow(cr(PI), s - 1.0);
        return Some(f * sn * g * z1);
    }
    if is_real_c(s) || (s.im.abs() < 30.0) {
        return hurwitz_zeta(s, cr(1.0));
    }
    hurwitz_zeta(s, cr(1.0))
}

/// Dirichlet eta function (Borwein's algorithm for `Re s > 0`).
pub fn eta_c(s: C64) -> Option<C64> {
    if s == cr(1.0) {
        return Some(cr(std::f64::consts::LN_2));
    }
    if s.re < 0.0 {
        let z = zeta_c(s)?;
        return Some((1.0 - cpow(cr(2.0), 1.0 - s)) * z);
    }
    let n = 30usize;
    // d_k = n sum_{i=0}^{k} (n+i-1)! 4^i / ((n-i)! (2i)!)
    let mut d = vec![0.0f64; n + 1];
    let mut term = 1.0 / n as f64; // i = 0: (n-1)!/(n!) = 1/n
    let mut acc = term;
    d[0] = n as f64 * acc;
    for i in 1..=n {
        // ratio term_i / term_{i-1} = (n+i-1) * 4 / ((n-i+1) (2i)(2i-1))
        term *= (n + i - 1) as f64 * 4.0 * (n - i + 1) as f64 / ((2 * i) as f64 * (2 * i - 1) as f64);
        acc += term;
        d[i] = n as f64 * acc;
    }
    let dn = d[n];
    let mut sum = C64::new(0.0, 0.0);
    for k in 0..n {
        let t = cpow(cr((k + 1) as f64), -s) * (d[k] - dn);
        if k % 2 == 0 {
            sum += t;
        } else {
            sum -= t;
        }
    }
    Some(-sum / dn)
}

/// Polylogarithm `Li_s(z)`.
pub fn polylog_c(s: C64, z: C64) -> Option<C64> {
    if !is_real_c(s) {
        return None;
    }
    let sr = s.re;
    if z == cr(0.0) {
        return Some(cr(0.0));
    }
    if z == cr(1.0) {
        return zeta_c(s);
    }
    if z == cr(-1.0) {
        return Some(-eta_c(s)?);
    }
    let sint = sr == sr.floor() && sr.abs() < 1e6;
    if sint && sr <= 0.0 {
        return Some(polylog_neg_int((-sr) as usize, z));
    }
    if sint && sr == 1.0 {
        return Some(-(1.0 - z).ln());
    }
    let az = z.norm();
    if az <= 0.5 {
        return Some(polylog_series(sr, z));
    }
    if is_real_c(z) && z.re < 0.0 && az < 2.0 {
        // Li_s(z) = 2^(1-s) Li_s(z^2) - Li_s(-z)
        let a = polylog_c(s, z * z)?;
        let b = polylog_c(s, -z)?;
        return Some(cr((a * 2f64.powf(1.0 - sr) - b).re));
    }
    if az < 2.0 {
        return polylog_logseries(sr, z);
    }
    if sint {
        // inversion: Li_s(z) = -(-1)^s Li_s(1/z) - (2 pi i)^s / s! B_s(1/2 + ln(-z)/(2 pi i))
        let n = sr as usize;
        let inv = polylog_c(s, z.inv())?;
        let tpi = C64::new(0.0, 2.0 * PI);
        let w = 0.5 + (-z).ln() / tpi;
        let bpoly = bernoulli_poly_c(n, w);
        let f = tpi.powi(n as i32) / gamma_real(n as f64 + 1.0);
        let sgn = if n % 2 == 0 { 1.0 } else { -1.0 };
        let mut r = -inv * sgn - f * bpoly;
        if z.im == 0.0 && z.re > 1.0 {
            // value on the cut: approach from below (mpmath convention)
            r = C64::new(r.re, -r.im.abs());
        }
        return Some(r);
    }
    None
}

fn bernoulli_poly_c(n: usize, w: C64) -> C64 {
    let c = crate::special::integers::bernoulli_poly_coeffs(n);
    let mut acc = C64::new(0.0, 0.0);
    for k in (0..=n).rev() {
        let ck = c[k].numer().to_f64().unwrap() / c[k].denom().to_f64().unwrap();
        acc = acc * w + ck;
    }
    acc
}

fn polylog_series(s: f64, z: C64) -> C64 {
    let mut sum = C64::new(0.0, 0.0);
    let mut zk = z;
    for k in 1..400 {
        let t = zk * (k as f64).powf(-s);
        sum += t;
        if t.norm() < 1e-17 * sum.norm() {
            break;
        }
        zk *= z;
    }
    sum
}

fn polylog_neg_int(n: usize, z: C64) -> C64 {
    // Li_{-n}(z) = z sum_{k=0}^{n-1} A(n,k) z^k / (1-z)^{n+1}  (Eulerian numbers)
    if n == 0 {
        return z / (1.0 - z);
    }
    let mut row = vec![1.0f64];
    for m in 2..=n {
        let mut next = vec![0.0f64; m];
        for k in 0..m {
            let a = if k < row.len() { row[k] * (k + 1) as f64 } else { 0.0 };
            let b = if k >= 1 && k - 1 < row.len() { row[k - 1] * (m - k) as f64 } else { 0.0 };
            next[k] = a + b;
        }
        row = next;
    }
    let mut p = C64::new(0.0, 0.0);
    for k in (0..row.len()).rev() {
        p = p * z + row[k];
    }
    z * p / (1.0 - z).powi(n as i32 + 1)
}

fn polylog_logseries(s: f64, z: C64) -> Option<C64> {
    let mu = z.ln();
    let sint = s == s.floor();
    let mut sum = C64::new(0.0, 0.0);
    let mut mk = C64::new(1.0, 0.0); // mu^k / k!
    for k in 0..200usize {
        let sk = s - k as f64;
        if sint && (sk - 1.0).abs() < 0.5 {
            // k = s - 1: harmonic term
            let n = k; // s - 1
            let mut h = 0.0;
            for j in 1..=n {
                h += 1.0 / j as f64;
            }
            sum += mk * (h - (-mu).ln());
        } else {
            let zv = zeta_c(cr(sk))?;
            sum += mk * zv;
        }
        mk = mk * mu / (k + 1) as f64;
        if k > s as usize + 3 && mk.norm() < 1e-18 * sum.norm().max(1e-300) {
            break;
        }
    }
    if !sint {
        sum += gamma_c(cr(1.0 - s))? * cpow(-mu, cr(s - 1.0));
    }
    Some(sum)
}

/// Lerch transcendent `Phi(z, s, a) = sum z^k / (k + a)^s` for `|z| < 1`.
pub fn lerchphi_c(z: C64, s: C64, a: C64) -> Option<C64> {
    if z.norm() >= 0.98 {
        if z == cr(1.0) {
            return hurwitz_zeta(s, a);
        }
        return None;
    }
    let mut sum = C64::new(0.0, 0.0);
    let mut zk = C64::new(1.0, 0.0);
    for k in 0..5000 {
        let t = zk * cpow(a + k as f64, -s);
        sum += t;
        if t.norm() < 1e-17 * sum.norm().max(1e-300) && k > 3 {
            break;
        }
        zk *= z;
    }
    Some(sum)
}

/// Beta function.
pub fn beta_c(x: C64, y: C64) -> Option<C64> {
    if is_real_c(x) && is_real_c(y) {
        let (a, b) = (x.re, y.re);
        let ga = gamma_real(a);
        let gb = gamma_real(b);
        let gab = gamma_real(a + b);
        if ga.is_finite() && gb.is_finite() && gab.is_finite() && gab != 0.0 {
            return Some(cr(ga * gb / gab));
        }
        if a > 0.0 && b > 0.0 {
            let l = lngamma_dd(DD::new(a)) + lngamma_dd(DD::new(b)) - lngamma_dd(DD::new(a + b));
            return Some(cr(l.exp().to_f64()));
        }
        if (a + b) <= 0.0 && (a + b) == (a + b).floor() {
            return Some(cr(0.0));
        }
        return None;
    }
    Some(gamma_c(x)? * gamma_c(y)? / gamma_c(x + y)?)
}

/// Regularisation helper: lower incomplete gamma series `sum x^k/(a)_{k+1}`.
fn gser(a: f64, x: f64) -> f64 {
    let mut sum = DD::new(1.0 / a);
    let mut del = 1.0 / a;
    let mut ap = a;
    for _ in 0..2000 {
        ap += 1.0;
        del *= x / ap;
        sum = sum.add_f64(del);
        if del.abs() < sum.hi.abs() * 1e-18 {
            break;
        }
    }
    sum.to_f64()
}

/// Continued fraction for `Gamma(a, x) e^x x^-a` (modified Lentz).
fn gcf(a: f64, x: f64) -> f64 {
    let tiny = 1e-300;
    let mut b = x + 1.0 - a;
    let mut cc = 1.0 / tiny;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..5000 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < tiny {
            d = tiny;
        }
        cc = b + an / cc;
        if cc.abs() < tiny {
            cc = tiny;
        }
        d = 1.0 / d;
        let del = d * cc;
        h *= del;
        if (del - 1.0).abs() < 1e-17 {
            break;
        }
    }
    h
}

fn xa_emx(a: f64, x: f64) -> f64 {
    // x^a e^-x computed in log space with double-double
    let l = DD::new(x).ln().mul_f64(a) - DD::new(x);
    l.exp().to_f64()
}

/// Lower incomplete gamma, real arguments with `x >= 0`.
pub fn lowergamma_real(a: f64, x: f64) -> Option<f64> {
    if x < 0.0 {
        return None;
    }
    if x == 0.0 {
        return if a > 0.0 { Some(0.0) } else { None };
    }
    if a <= 0.0 && a == a.floor() {
        return None;
    }
    if x < a + 1.0 || a < 0.0 || x < 1.0 {
        return Some(xa_emx(a, x) * gser(a, x));
    }
    Some(gamma_real(a) - xa_emx(a, x) * gcf(a, x))
}

/// Upper incomplete gamma, real arguments with `x > 0`.
pub fn uppergamma_real(a: f64, x: f64) -> Option<f64> {
    if x < 0.0 {
        return None;
    }
    if x == 0.0 {
        return if a > 0.0 { Some(gamma_real(a)) } else { None };
    }
    if x >= a + 1.0 && x >= 1.0 {
        return Some(xa_emx(a, x) * gcf(a, x));
    }
    if a <= 0.0 && a == a.floor() {
        // Gamma(-n, x) = (-1)^n/n! (E1(x) - e^-x sum_{k<n} (-1)^k k!/x^{k+1})
        let n = (-a) as i64;
        let e1 = super::erf::e1_real(x)?;
        let mut s = 0.0;
        let mut f = 1.0;
        for k in 0..n {
            if k > 0 {
                f *= k as f64;
            }
            let t = f / x.powi(k as i32 + 1);
            s += if k % 2 == 0 { t } else { -t };
        }
        let mut nf = 1.0;
        for k in 2..=n {
            nf *= k as f64;
        }
        let sgn = if n % 2 == 0 { 1.0 } else { -1.0 };
        return Some(sgn / nf * (e1 - (-x).exp() * s));
    }
    Some(gamma_real(a) - xa_emx(a, x) * gser(a, x))
}

/// Generalised harmonic number `H_x^{(m)} = zeta(m) - zeta(m, x+1)`.
pub fn harmonic_c(x: C64, m: C64) -> Option<C64> {
    if m == cr(1.0) {
        return Some(digamma_c(x + 1.0)? + consts().euler.to_f64());
    }
    Some(zeta_c(m)? - hurwitz_zeta(m, x + 1.0)?)
}

/// `n!` for real/complex n.
pub fn factorial_c(n: C64) -> Option<C64> {
    gamma_c(n + 1.0)
}

/// `binomial(n, k)` via gamma functions (limits at integers).
pub fn binomial_c(n: C64, k: C64) -> Option<C64> {
    if is_real_c(n) && is_real_c(k) {
        let (nr, kr) = (n.re, k.re);
        if kr == kr.floor() && kr >= 0.0 && kr < 1e6 {
            // falling product, exact for integer arguments
            let mut r = 1.0;
            let kk = kr as i64;
            for i in 1..=kk {
                r *= (nr - kk as f64 + i as f64) / i as f64;
            }
            return Some(cr(r));
        }
        let a = rgamma_real(kr + 1.0) * rgamma_real(nr - kr + 1.0);
        return Some(cr(gamma_real(nr + 1.0) * a));
    }
    Some(gamma_c(n + 1.0)? / (gamma_c(k + 1.0)? * gamma_c(n - k + 1.0)?))
}

/// Rising factorial `Gamma(x+k)/Gamma(x)`.
pub fn rf_c(x: C64, k: C64) -> Option<C64> {
    if is_real_c(x) && is_real_c(k) {
        let (xr, kr) = (x.re, k.re);
        if kr == kr.floor() && kr.abs() < 1e5 {
            let mut r = 1.0;
            if kr >= 0.0 {
                for i in 0..(kr as i64) {
                    r *= xr + i as f64;
                }
            } else {
                for i in 1..=((-kr) as i64) {
                    r /= xr - i as f64;
                }
            }
            return Some(cr(r));
        }
        if xr > 0.0 && xr + kr > 0.0 {
            let l = lngamma_dd(DD::new(xr + kr)) - lngamma_dd(DD::new(xr));
            return Some(cr(l.exp().to_f64()));
        }
        return Some(cr(gamma_real(xr + kr) * rgamma_real(xr)));
    }
    Some(gamma_c(x + k)? / gamma_c(x)?)
}

/// Falling factorial `Gamma(x+1)/Gamma(x-k+1)`.
pub fn ff_c(x: C64, k: C64) -> Option<C64> {
    if is_real_c(x) && is_real_c(k) {
        let (xr, kr) = (x.re, k.re);
        if kr == kr.floor() && kr.abs() < 1e5 {
            let mut r = 1.0;
            if kr >= 0.0 {
                for i in 0..(kr as i64) {
                    r *= xr - i as f64;
                }
            } else {
                for i in 1..=((-kr) as i64) {
                    r /= xr + i as f64;
                }
            }
            return Some(cr(r));
        }
        return Some(cr(gamma_real(xr + 1.0) * rgamma_real(xr - kr + 1.0)));
    }
    Some(gamma_c(x + 1.0)? / gamma_c(x - k + 1.0)?)
}

/// Catalan number for real/complex `n`: `4^n Gamma(n+1/2) / (sqrt(pi) Gamma(n+2))`.
pub fn catalan_c(n: C64) -> Option<C64> {
    let sp = consts().sqrt_pi.to_f64();
    if is_real_c(n) {
        let x = n.re;
        return Some(cr(4f64.powf(x) * gamma_real(x + 0.5) * rgamma_real(x + 2.0) / sp));
    }
    Some(cpow(cr(4.0), n) * gamma_c(n + 0.5)? / (gamma_c(n + 2.0)? * sp))
}

/// Fibonacci number for real/complex `n` (Binet's formula).
pub fn fibonacci_c(n: C64) -> C64 {
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    let s5 = 5f64.sqrt();
    if is_real_c(n) {
        let x = n.re;
        return cr((phi.powf(x) - cospi(x) * phi.powf(-x)) / s5);
    }
    (cpow(cr(phi), n) - (n * PI).cos() * cpow(cr(phi), -n)) / s5
}

/// Lucas number for real/complex `n`.
pub fn lucas_c(n: C64) -> C64 {
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    if is_real_c(n) {
        let x = n.re;
        return cr(phi.powf(x) + cospi(x) * phi.powf(-x));
    }
    cpow(cr(phi), n) + (n * PI).cos() * cpow(cr(phi), -n)
}

/// Bernoulli number/polynomial for non integer index:
/// `B_s(a) = -s zeta(1-s, a)`.
pub fn bernoulli_c(s: C64, a: Option<C64>) -> Option<C64> {
    let a = a.unwrap_or(cr(1.0));
    if a == cr(1.0) {
        return Some(-s * zeta_c(1.0 - s)?);
    }
    Some(-s * hurwitz_zeta(1.0 - s, a)?)
}

/// Euler number for non integer index.
pub fn euler_c(s: C64, x: Option<C64>) -> Option<C64> {
    match x {
        None => {
            // E(s) = 2^{s+1} (zeta(-s, 1/4) - 2^{-s-1}... ) via Euler polynomial at 1/2
            let e = euler_poly_c(s, cr(0.5))?;
            Some(cpow(cr(2.0), s) * e)
        }
        Some(x) => euler_poly_c(s, x),
    }
}

fn euler_poly_c(s: C64, x: C64) -> Option<C64> {
    // E_s(x) = 2/(s+1) (B_{s+1}(x) - 2^{s+1} B_{s+1}(x/2))
    let s1 = s + 1.0;
    let b1 = bernoulli_c(s1, Some(x))?;
    let b2 = bernoulli_c(s1, Some(x * 0.5))?;
    Some(2.0 / s1 * (b1 - cpow(cr(2.0), s1) * b2))
}
