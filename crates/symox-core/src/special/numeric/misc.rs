//! Elliptic integrals, Lambert W, hypergeometric series and orthogonal
//! polynomials of non integer degree.

use num_complex::Complex64 as C64;

use super::gamma::{gamma_c, gamma_real, rgamma_real};
use super::{cr, is_real_c, PI};

// ----------------------------------------------------------------------
// Carlson symmetric integrals (real arguments)
// ----------------------------------------------------------------------

pub fn rf(mut x: f64, mut y: f64, mut z: f64) -> Option<f64> {
    if x < 0.0 || y < 0.0 || z < 0.0 || (x + y == 0.0) || (x + z == 0.0) || (y + z == 0.0) {
        return None;
    }
    for _ in 0..100 {
        let (sx, sy, sz) = (x.sqrt(), y.sqrt(), z.sqrt());
        let lam = sx * sy + sy * sz + sz * sx;
        x = (x + lam) / 4.0;
        y = (y + lam) / 4.0;
        z = (z + lam) / 4.0;
        let a = (x + y + z) / 3.0;
        let d = ((a - x).abs()).max((a - y).abs()).max((a - z).abs()) / a;
        if d < 1e-4 {
            break;
        }
    }
    let a = (x + y + z) / 3.0;
    let (xx, yy) = (1.0 - x / a, 1.0 - y / a);
    let zz = -(xx + yy);
    let e2 = xx * yy - zz * zz;
    let e3 = xx * yy * zz;
    Some((1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0) / a.sqrt())
}

pub fn rd(mut x: f64, mut y: f64, mut z: f64) -> Option<f64> {
    if x < 0.0 || y < 0.0 || z <= 0.0 || x + y == 0.0 {
        return None;
    }
    let mut sum = 0.0;
    let mut fac = 1.0;
    for _ in 0..100 {
        let (sx, sy, sz) = (x.sqrt(), y.sqrt(), z.sqrt());
        let lam = sx * sy + sy * sz + sz * sx;
        sum += fac / (sz * (z + lam));
        fac /= 4.0;
        x = (x + lam) / 4.0;
        y = (y + lam) / 4.0;
        z = (z + lam) / 4.0;
        let a = (x + y + 3.0 * z) / 5.0;
        let d = ((a - x).abs()).max((a - y).abs()).max((a - z).abs()) / a;
        if d < 1e-4 {
            break;
        }
    }
    let a = (x + y + 3.0 * z) / 5.0;
    let (xx, yy) = (1.0 - x / a, 1.0 - y / a);
    let zz = -(xx + yy) / 3.0;
    let e2 = xx * yy - 6.0 * zz * zz;
    let e3 = (3.0 * xx * yy - 8.0 * zz * zz) * zz;
    let e4 = 3.0 * (xx * yy - zz * zz) * zz * zz;
    let e5 = xx * yy * zz * zz * zz;
    let s = 1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0 - 3.0 * e4 / 22.0 - 9.0 * e2 * e3 / 52.0
        + 3.0 * e5 / 26.0;
    Some(3.0 * sum + fac * s / (a * a.sqrt()))
}

fn rc(x: f64, y: f64) -> Option<f64> {
    if y <= 0.0 || x < 0.0 {
        return None;
    }
    if x < y {
        let d = (y - x).sqrt();
        Some((d / y.sqrt()).atan2((x / y).sqrt()) / d)
    } else if x > y {
        let d = (x - y).sqrt();
        Some(((x.sqrt() + d) / y.sqrt()).ln() / d)
    } else {
        Some(1.0 / x.sqrt())
    }
}

pub fn rj(mut x: f64, mut y: f64, mut z: f64, mut p: f64) -> Option<f64> {
    if x < 0.0 || y < 0.0 || z < 0.0 || p <= 0.0 {
        return None;
    }
    let a0 = (x + y + z + 2.0 * p) / 5.0;
    let delta = (p - x) * (p - y) * (p - z);
    let (x0, y0, z0) = (x, y, z);
    let mut a = a0;
    let mut sum = 0.0;
    let mut f4 = 1.0; // 4^-m
    let mut m = 0;
    loop {
        let (sx, sy, sz, sp) = (x.sqrt(), y.sqrt(), z.sqrt(), p.sqrt());
        let lam = sx * sy + sy * sz + sz * sx;
        let d = (sp + sx) * (sp + sy) * (sp + sz);
        let e = f4 * f4 * f4 * delta / (d * d);
        sum += f4 * rc(1.0, 1.0 + e)? / d;
        x = (x + lam) / 4.0;
        y = (y + lam) / 4.0;
        z = (z + lam) / 4.0;
        p = (p + lam) / 4.0;
        a = (a + lam) / 4.0;
        f4 /= 4.0;
        m += 1;
        let q = (a0 - x0).abs().max((a0 - y0).abs()).max((a0 - z0).abs()).max((a0 - p).abs());
        if f4 * q < 1e-4 * a.abs() * 1e-0 && m > 2 {
            let dev = ((a - x).abs()).max((a - y).abs()).max((a - z).abs()).max((a - p).abs()) / a;
            if dev < 1e-4 {
                break;
            }
        }
        if m > 100 {
            break;
        }
    }
    let xx = (a0 - x0) * f4 / a;
    let yy = (a0 - y0) * f4 / a;
    let zz = (a0 - z0) * f4 / a;
    let pp = -(xx + yy + zz) / 2.0;
    let e2 = xx * yy + xx * zz + yy * zz - 3.0 * pp * pp;
    let e3 = xx * yy * zz + 2.0 * e2 * pp + 4.0 * pp * pp * pp;
    let e4 = (2.0 * xx * yy * zz + e2 * pp + 3.0 * pp * pp * pp) * pp;
    let e5 = xx * yy * zz * pp * pp;
    let s = 1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0 - 3.0 * e4 / 22.0 - 9.0 * e2 * e3 / 52.0
        + 3.0 * e5 / 26.0;
    Some(f4 * s / (a * a.sqrt()) + 6.0 * sum)
}

fn agm(a: f64, b: f64) -> (f64, f64) {
    // returns (agm, sum 2^(n-1) c_n^2) for E
    let (mut a, mut b) = (a, b);
    let mut c2sum = 0.0;
    let mut pw = 0.5;
    let c0 = a * a - b * b;
    c2sum += pw * c0;
    for _ in 0..60 {
        let an = (a + b) / 2.0;
        let bn = (a * b).sqrt();
        let cn = (a - b) / 2.0;
        pw *= 2.0;
        c2sum += pw * cn * cn;
        a = an;
        b = bn;
        if cn.abs() < 1e-17 * a {
            break;
        }
    }
    (a, c2sum)
}

/// Complete elliptic integral of the first kind `K(m)`.
pub fn ellipk_c(m: C64) -> Option<C64> {
    if !is_real_c(m) {
        return None;
    }
    let m = m.re;
    if m == 1.0 {
        return None;
    }
    if m > 1.0 {
        let k1 = ellipk_c(cr(1.0 / m))?.re;
        let k2 = ellipk_c(cr(1.0 - 1.0 / m))?.re;
        return Some(C64::new(k1, -k2) / m.sqrt());
    }
    let (g, _) = agm(1.0, (1.0 - m).sqrt());
    Some(cr(PI / (2.0 * g)))
}

/// Complete elliptic integral of the second kind `E(m)`.
pub fn ellipe_c(m: C64) -> Option<C64> {
    if !is_real_c(m) {
        return None;
    }
    let m = m.re;
    if m == 1.0 {
        return Some(cr(1.0));
    }
    if m > 1.0 {
        return None;
    }
    let (g, s) = agm(1.0, (1.0 - m).sqrt());
    let k = PI / (2.0 * g);
    Some(cr(k * (1.0 - s)))
}

/// Reduce `phi = j pi + r`, `|r| <= pi/2`.
fn reduce_phi(phi: f64) -> (f64, f64) {
    let j = (phi / PI).round();
    (j, phi - j * PI)
}

/// Incomplete elliptic integral of the first kind `F(phi | m)`.
pub fn ellipf_c(phi: C64, m: C64) -> Option<C64> {
    if !is_real_c(phi) || !is_real_c(m) {
        return None;
    }
    let (j, r) = reduce_phi(phi.re);
    let m = m.re;
    let (s, c) = r.sin_cos();
    if m * s * s > 1.0 {
        return None;
    }
    let f = s * rf(c * c, 1.0 - m * s * s, 1.0)?;
    let k = if j != 0.0 { ellipk_c(cr(m))?.re } else { 0.0 };
    Some(cr(2.0 * j * k + f))
}

/// Incomplete elliptic integral of the second kind `E(phi | m)`.
pub fn ellipe2_c(phi: C64, m: C64) -> Option<C64> {
    if !is_real_c(phi) || !is_real_c(m) {
        return None;
    }
    let (j, r) = reduce_phi(phi.re);
    let m = m.re;
    let (s, c) = r.sin_cos();
    let d = 1.0 - m * s * s;
    if d < 0.0 {
        return None;
    }
    let e = s * rf(c * c, d, 1.0)? - m / 3.0 * s * s * s * rd(c * c, d, 1.0)?;
    let ec = if j != 0.0 { ellipe_c(cr(m))?.re } else { 0.0 };
    Some(cr(2.0 * j * ec + e))
}

/// Elliptic integral of the third kind `Pi(n; phi | m)` (complete when
/// `phi` is `None`).
pub fn ellippi_c(n: C64, phi: Option<C64>, m: C64) -> Option<C64> {
    if !is_real_c(n) || !is_real_c(m) {
        return None;
    }
    let (n, m) = (n.re, m.re);
    match phi {
        None => {
            if m >= 1.0 || n >= 1.0 {
                return None;
            }
            Some(cr(rf(0.0, 1.0 - m, 1.0)? + n / 3.0 * rj(0.0, 1.0 - m, 1.0, 1.0 - n)?))
        }
        Some(p) => {
            if !is_real_c(p) {
                return None;
            }
            let (j, r) = reduce_phi(p.re);
            let (s, c) = r.sin_cos();
            let d = 1.0 - m * s * s;
            let q = 1.0 - n * s * s;
            if d < 0.0 || q <= 0.0 {
                return None;
            }
            let v = s * rf(c * c, d, 1.0)? + n / 3.0 * s * s * s * rj(c * c, d, 1.0, q)?;
            let full = if j != 0.0 { ellippi_c(cr(n), None, cr(m))?.re } else { 0.0 };
            Some(cr(2.0 * j * full + v))
        }
    }
}

// ----------------------------------------------------------------------
// Lambert W
// ----------------------------------------------------------------------

pub fn lambertw_c(z: C64, k: i64) -> Option<C64> {
    let em1 = (-1.0f64).exp();
    if z == cr(0.0) {
        return if k == 0 { Some(cr(0.0)) } else { None };
    }
    let real0 = k == 0 && is_real_c(z) && z.re >= -em1;
    let realm1 = k == -1 && is_real_c(z) && z.re >= -em1 && z.re < 0.0;
    if real0 || realm1 {
        let x = z.re;
        let p2 = 2.0 * (std::f64::consts::E * x + 1.0);
        let mut w = if p2 < 0.3 {
            let p = p2.max(0.0).sqrt();
            if k == 0 { -1.0 + p - p * p / 3.0 } else { -1.0 - p - p * p / 3.0 }
        } else if k == 0 {
            if x < 3.0 { (1.0 + x).ln() * 0.8 } else { let l = x.ln(); l - l.ln() }
        } else {
            let l = (-x).ln();
            l - (-l).ln()
        };
        for _ in 0..100 {
            let ew = w.exp();
            let f = w * ew - x;
            let wp1 = w + 1.0;
            if wp1 == 0.0 {
                break;
            }
            let dw = f / (ew * wp1 - (w + 2.0) * f / (2.0 * wp1));
            w -= dw;
            if dw.abs() <= 1e-16 * (1.0 + w.abs()) {
                break;
            }
        }
        return Some(cr(w));
    }
    // complex: initial guess
    let two_pi_ik = C64::new(0.0, 2.0 * PI * k as f64);
    let mut w = if k == 0 && z.norm() < 1.0 {
        let p = (2.0 * (std::f64::consts::E * z + 1.0)).sqrt();
        if (z + em1).norm() < 0.3 {
            -1.0 + p - p * p / 3.0
        } else {
            z * (1.0 - z)
        }
    } else {
        let l1 = z.ln() + two_pi_ik;
        l1 - l1.ln()
    };
    for _ in 0..200 {
        let ew = w.exp();
        let f = w * ew - z;
        let wp1 = w + 1.0;
        let dw = f / (ew * wp1 - (w + 2.0) * f / (2.0 * wp1));
        w -= dw;
        if dw.norm() <= 1e-16 * (1.0 + w.norm()) {
            break;
        }
    }
    Some(w)
}

// ----------------------------------------------------------------------
// hypergeometric series
// ----------------------------------------------------------------------

/// `pFq(a; b; z)` by its power series (`p <= q`, or `p = q + 1` with `|z| < 1`).
pub fn hyper_c(a: &[C64], b: &[C64], z: C64) -> Option<C64> {
    if a.len() > b.len() + 1 {
        return None;
    }
    if a.iter().chain(b.iter()).all(|v| is_real_c(*v)) && is_real_c(z) {
        let ar: Vec<f64> = a.iter().map(|v| v.re).collect();
        let br: Vec<f64> = b.iter().map(|v| v.re).collect();
        return hyper_dd(&ar, &br, z.re).map(cr);
    }
    if a.len() == b.len() + 1 && z.norm() >= 0.95 {
        return None;
    }
    for bb in b {
        if is_real_c(*bb) && bb.re <= 0.0 && bb.re == bb.re.floor() {
            return None;
        }
    }
    let mut term = C64::new(1.0, 0.0);
    let mut sum = term;
    for k in 0..20000 {
        let kf = k as f64;
        let mut r = z / (kf + 1.0);
        for aa in a {
            r *= *aa + kf;
        }
        for bb in b {
            r /= *bb + kf;
        }
        term *= r;
        sum += term;
        if term == C64::new(0.0, 0.0) {
            break;
        }
        if term.norm() < 1e-17 * sum.norm() && kf > z.norm() {
            break;
        }
    }
    Some(sum)
}

/// Real hypergeometric series in double-double arithmetic.
fn hyper_dd(a: &[f64], b: &[f64], z: f64) -> Option<f64> {
    use super::dd::DD;
    if a.len() == b.len() + 1 && z.abs() >= 0.95 {
        return None;
    }
    for bb in b {
        if *bb <= 0.0 && *bb == bb.floor() {
            return None;
        }
    }
    let zd = DD::new(z);
    let mut term = DD::ONE;
    let mut sum = DD::ONE;
    for k in 0..20000 {
        let kf = k as f64;
        let mut num = zd;
        for aa in a {
            num = num * DD::new(*aa).add_f64(kf);
        }
        let mut den = DD::new(kf + 1.0);
        for bb in b {
            den = den * DD::new(*bb).add_f64(kf);
        }
        term = term * num / den;
        sum = sum + term;
        if term.hi == 0.0 {
            break;
        }
        if term.hi.abs() < 1e-33 * sum.hi.abs() && kf > z.abs() {
            break;
        }
    }
    Some(sum.to_f64())
}

fn binom_real(n: f64, k: f64) -> f64 {
    gamma_real(n + 1.0) * rgamma_real(k + 1.0) * rgamma_real(n - k + 1.0)
}

// ----------------------------------------------------------------------
// orthogonal polynomials of general degree
// ----------------------------------------------------------------------

/// `which`: index into the orthopoly heads (see `numeric::apply`).
pub fn legendre_c(n: C64, x: C64) -> Option<C64> {
    if !is_real_c(n) {
        return None;
    }
    let z = (1.0 - x) / 2.0;
    if z.norm() > 0.8 {
        return None;
    }
    hyper_c(&[-n, n + 1.0], &[cr(1.0)], z)
}

pub fn assoc_legendre_c(n: C64, m: C64, x: C64) -> Option<C64> {
    if !is_real_c(n) || !is_real_c(m) || !is_real_c(x) {
        return None;
    }
    let (nu, mu, xr) = (n.re, m.re, x.re);
    if xr.abs() >= 1.0 || mu != mu.floor() {
        return None;
    }
    if (1.0 - xr) / 2.0 > 0.8 {
        return None;
    }
    let mi = mu.abs();
    // P^{-m} = 1/m! ((1-x)/(1+x))^{m/2} 2F1(-nu, nu+1; 1+m; (1-x)/2)
    let f = hyper_c(&[cr(-nu), cr(nu + 1.0)], &[cr(1.0 + mi)], cr((1.0 - xr) / 2.0))?.re;
    let pm = ((1.0 - xr) / (1.0 + xr)).powf(mi / 2.0) * rgamma_real(mi + 1.0) * f;
    if mu < 0.0 {
        return Some(cr(pm));
    }
    // P^m = (-1)^m Gamma(nu+m+1)/Gamma(nu-m+1) P^{-m}
    let sgn = if (mi as i64) % 2 == 0 { 1.0 } else { -1.0 };
    Some(cr(sgn * gamma_real(nu + mi + 1.0) * rgamma_real(nu - mi + 1.0) * pm))
}

pub fn chebyshev_c(first: bool, n: C64, x: C64) -> Option<C64> {
    if !is_real_c(n) || !is_real_c(x) {
        return None;
    }
    let (nu, xr) = (n.re, x.re);
    if xr.abs() <= 1.0 {
        let t = xr.acos();
        if first {
            return Some(cr((nu * t).cos()));
        }
        let s = t.sin();
        if s == 0.0 {
            // limit
            let sign = if xr > 0.0 { 1.0 } else { (nu * PI).cos() };
            return Some(cr((nu + 1.0) * sign));
        }
        return Some(cr(((nu + 1.0) * t).sin() / s));
    }
    if xr > 1.0 {
        let t = xr.acosh();
        if first {
            return Some(cr((nu * t).cosh()));
        }
        return Some(cr(((nu + 1.0) * t).sinh() / t.sinh()));
    }
    None
}

pub fn hermite_c(n: C64, x: C64) -> Option<C64> {
    if !is_real_c(n) {
        return None;
    }
    let nu = n.re;
    let x2 = x * x;
    let sp = PI.sqrt();
    let a = hyper_c(&[cr(-nu / 2.0)], &[cr(0.5)], x2)? * rgamma_real((1.0 - nu) / 2.0);
    let b = hyper_c(&[cr((1.0 - nu) / 2.0)], &[cr(1.5)], x2)? * rgamma_real(-nu / 2.0) * 2.0 * x;
    Some((a - b) * 2f64.powf(nu) * sp)
}

pub fn laguerre_c(n: C64, alpha: C64, x: C64) -> Option<C64> {
    if !is_real_c(n) || !is_real_c(alpha) {
        return None;
    }
    let nu = n.re.trunc();
    let a = alpha.re;
    Some(hyper_c(&[cr(-nu)], &[cr(a + 1.0)], x)? * binom_real(nu + a, nu))
}

pub fn jacobi_c(n: C64, a: C64, b: C64, x: C64) -> Option<C64> {
    if !is_real_c(n) {
        return None;
    }
    let nu = n.re.trunc();
    let z = (1.0 - x) / 2.0;
    let pre = if is_real_c(a) {
        binom_real(nu + a.re, nu)
    } else {
        return None;
    };
    Some(hyper_c(&[cr(-nu), a + b + nu + 1.0], &[a + 1.0], z)? * pre)
}

pub fn gegenbauer_c(n: C64, a: C64, x: C64) -> Option<C64> {
    if !is_real_c(n) || !is_real_c(a) {
        return None;
    }
    let nu = n.re.trunc();
    let ar = a.re;
    let z = (1.0 - x) / 2.0;
    let pre = binom_real(nu + 2.0 * ar - 1.0, nu);
    Some(hyper_c(&[cr(-nu), cr(nu + 2.0 * ar)], &[cr(ar + 0.5)], z)? * pre)
}

/// Spherical harmonic `Y_n^m(theta, phi)` for integer `n, m`.
pub fn ynm_c(n: C64, m: C64, theta: C64, phi: C64) -> Option<C64> {
    if !is_real_c(n) || !is_real_c(m) || !is_real_c(theta) {
        return None;
    }
    let (nn, mm) = (n.re, m.re);
    if nn != nn.floor() || mm != mm.floor() || nn < 0.0 || mm.abs() > nn {
        return None;
    }
    let ma = mm.abs() as i64;
    let l = nn as i64;
    let x = theta.re.cos();
    // P_l^m(x) with Condon-Shortley phase by recurrence
    let s = (1.0 - x * x).max(0.0).sqrt();
    let mut pmm = 1.0;
    let mut f = 1.0;
    for _ in 0..ma {
        pmm *= -f * s;
        f += 2.0;
    }
    let p = if l == ma {
        pmm
    } else {
        let mut p1 = x * (2 * ma + 1) as f64 * pmm;
        let mut p0 = pmm;
        for ll in (ma + 2)..=l {
            let p2 = ((2 * ll - 1) as f64 * x * p1 - (ll + ma - 1) as f64 * p0) / (ll - ma) as f64;
            p0 = p1;
            p1 = p2;
        }
        p1
    };
    let mut ratio = 1.0;
    for k in (l - ma + 1)..=(l + ma) {
        ratio /= k as f64;
    }
    let norm = ((2 * l + 1) as f64 / (4.0 * PI) * ratio).sqrt();
    let e = (C64::new(0.0, ma as f64) * phi).exp();
    let y = e * (norm * p);
    if mm < 0.0 {
        let sgn = if ma % 2 == 0 { 1.0 } else { -1.0 };
        // Y_n^{-m} = (-1)^m exp(-2 i m phi) Y_n^m
        let e2 = (C64::new(0.0, -2.0 * ma as f64) * phi).exp();
        return Some(y * e2 * sgn);
    }
    let _ = gamma_c;
    Some(y)
}
