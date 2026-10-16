//! Error functions, Fresnel integrals and exponential integrals.

use num_complex::Complex64 as C64;

use super::dd::{consts, DD};
use super::{cr, is_real_c, PI};

const EULER: f64 = 0.577_215_664_901_532_9;

fn two_over_sqrt_pi() -> DD {
    DD::new(2.0) / consts().sqrt_pi
}

/// `e^{-x^2}` with the square computed exactly.
fn exp_mx2(x: f64) -> DD {
    let x2 = DD::new(x) * DD::new(x);
    (-x2).exp()
}

/// erf for real `|x| <= 3` in double-double:
/// `2x/sqrt(pi) e^{-x^2} sum (2x^2)^n / (2n+1)!!`.
fn erf_series_dd(x: f64) -> DD {
    let x2 = DD::new(x) * DD::new(x);
    let t2 = x2.mul_f64(2.0);
    let mut term = DD::ONE;
    let mut sum = DD::ONE;
    let mut n = 0.0;
    loop {
        n += 1.0;
        term = (term * t2).div_f64(2.0 * n + 1.0);
        sum = sum + term;
        if term.hi < 1e-34 * sum.hi || n > 500.0 {
            break;
        }
    }
    two_over_sqrt_pi().mul_f64(x) * (-x2).exp() * sum
}

/// erfc for `x >= 2` by the continued fraction (evaluated backwards).
fn erfc_cf(x: f64) -> f64 {
    let depth = (30.0 + 600.0 / (x * x)) as i64;
    let mut f = 0.0;
    for k in (1..=depth).rev() {
        f = (k as f64 / 2.0) / (x + f);
    }
    (exp_mx2(x).to_f64()) / (consts().sqrt_pi.to_f64() * (x + f))
}

pub fn erf_real(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x < 0.0 {
        return -erf_real(-x);
    }
    if x < 3.0 {
        return erf_series_dd(x).to_f64();
    }
    if x > 6.0 {
        return 1.0;
    }
    1.0 - erfc_cf(x)
}

pub fn erfc_real(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x < 0.0 {
        return 2.0 - erfc_real(-x);
    }
    if x < 2.0 {
        return (DD::ONE - erf_series_dd(x)).to_f64();
    }
    if x > 27.3 {
        return 0.0;
    }
    erfc_cf(x)
}

pub fn erfi_real(x: f64) -> f64 {
    if x < 0.0 {
        return -erfi_real(-x);
    }
    if x <= 10.0 {
        // 2/sqrt(pi) sum x^(2n+1) / (n! (2n+1))
        let x2 = DD::new(x) * DD::new(x);
        let mut p = DD::new(x); // x^(2n+1)/n!
        let mut sum = DD::new(x);
        let mut n = 0.0;
        loop {
            n += 1.0;
            p = (p * x2).div_f64(n);
            let t = p.div_f64(2.0 * n + 1.0);
            sum = sum + t;
            if t.hi < 1e-33 * sum.hi || n > 1000.0 {
                break;
            }
        }
        return (two_over_sqrt_pi() * sum).to_f64();
    }
    // asymptotic: e^{x^2}/(x sqrt(pi)) sum (2n-1)!!/(2x^2)^n
    let x2 = DD::new(x) * DD::new(x);
    let mut s = 1.0;
    let mut t = 1.0;
    let inv = 1.0 / (2.0 * x * x);
    for n in 1..60 {
        t *= (2 * n - 1) as f64 * inv;
        s += t;
        if t < 1e-17 {
            break;
        }
    }
    (x2.exp() * DD::new(s) / (consts().sqrt_pi.mul_f64(x))).to_f64()
}

/// Complex error function.
pub fn erf_c(z: C64) -> C64 {
    if is_real_c(z) {
        return cr(erf_real(z.re));
    }
    if z.re < 0.0 {
        return -erf_c(-z);
    }
    let az = z.norm();
    let tsp = 2.0 / consts().sqrt_pi.to_f64();
    if az > 6.0 && z.re > z.im.abs() * 0.5 {
        return 1.0 - erfc_asym(z);
    }
    if z.im.abs() > z.re {
        // Maclaurin series sum (-1)^n z^(2n+1)/(n!(2n+1))
        let z2 = z * z;
        let mut p = z;
        let mut sum = z;
        for n in 1..2000 {
            p = -p * z2 / n as f64;
            let t = p / (2 * n + 1) as f64;
            sum += t;
            if t.norm() < 1e-17 * sum.norm() {
                break;
            }
        }
        return sum * tsp;
    }
    // 2z/sqrt(pi) e^{-z^2} sum (2 z^2)^n / (2n+1)!!
    let z2 = z * z;
    let t2 = z2 * 2.0;
    let mut term = C64::new(1.0, 0.0);
    let mut sum = term;
    for n in 1..2000 {
        term = term * t2 / (2 * n + 1) as f64;
        sum += term;
        if term.norm() < 1e-17 * sum.norm() {
            break;
        }
    }
    z * tsp * (-z2).exp() * sum
}

/// Asymptotic / continued fraction erfc for large `|z|` in the right half plane.
fn erfc_asym(z: C64) -> C64 {
    let depth = 60 + (600.0 / z.norm_sqr()) as i64;
    let mut f = C64::new(0.0, 0.0);
    for k in (1..=depth).rev() {
        f = (k as f64 / 2.0) / (z + f);
    }
    (-z * z).exp() / (consts().sqrt_pi.to_f64() * (z + f))
}

pub fn erfc_c(z: C64) -> C64 {
    if is_real_c(z) {
        return cr(erfc_real(z.re));
    }
    if z.re > 2.0 && z.re > z.im.abs() {
        return erfc_asym(z);
    }
    1.0 - erf_c(z)
}

pub fn erfi_c(z: C64) -> C64 {
    if is_real_c(z) {
        return cr(erfi_real(z.re));
    }
    // erfi(z) = -i erf(iz)
    let iz = C64::new(-z.im, z.re);
    let e = erf_c(iz);
    C64::new(e.im, -e.re)
}

/// Inverse error function on (-1, 1).
pub fn erfinv_real(x: f64) -> Option<f64> {
    if !(-1.0..=1.0).contains(&x) {
        return None;
    }
    if x == 1.0 {
        return Some(f64::INFINITY);
    }
    if x == -1.0 {
        return Some(f64::NEG_INFINITY);
    }
    if x == 0.0 {
        return Some(0.0);
    }
    if x < 0.0 {
        return erfinv_real(-x).map(|v| -v);
    }
    if x > 0.9 {
        return erfcinv_real(1.0 - x);
    }
    // Winitzki initial guess, then Newton on erf
    let a = 0.147;
    let l = (1.0 - x * x).ln();
    let t = 2.0 / (PI * a) + l / 2.0;
    let mut y = ((t * t - l / a).sqrt() - t).sqrt();
    let sp = consts().sqrt_pi.to_f64();
    for _ in 0..6 {
        let err = erf_real(y) - x;
        let d = 2.0 / sp * (-y * y).exp();
        let step = err / d;
        // Halley correction
        y -= step / (1.0 + y * step);
        if step.abs() < 1e-17 * y.abs() {
            break;
        }
    }
    Some(y)
}

/// Inverse complementary error function on (0, 2).
pub fn erfcinv_real(x: f64) -> Option<f64> {
    if !(0.0..=2.0).contains(&x) {
        return None;
    }
    if x == 0.0 {
        return Some(f64::INFINITY);
    }
    if x == 2.0 {
        return Some(f64::NEG_INFINITY);
    }
    if x > 1.0 {
        return erfcinv_real(2.0 - x).map(|v| -v);
    }
    if x > 0.1 {
        return erfinv_real(1.0 - x);
    }
    // initial guess from the asymptotic expansion, Newton on erfc
    let l = -(x * consts().sqrt_pi.to_f64()).ln();
    let mut y = (l - 0.5 * l.ln()).max(0.5).sqrt();
    let sp = consts().sqrt_pi.to_f64();
    for _ in 0..30 {
        let err = erfc_real(y) - x;
        let d = -2.0 / sp * exp_mx2(y).to_f64();
        let step = err / d;
        y -= step / (1.0 + y * step);
        if step.abs() < 1e-17 * y.abs() {
            break;
        }
    }
    Some(y)
}

/// Fresnel integrals `(S(x), C(x))`, real `x`.
pub fn fresnel_real(x: f64) -> (f64, f64) {
    if x < 0.0 {
        let (s, c) = fresnel_real(-x);
        return (-s, -c);
    }
    if x <= 3.0 {
        // S = sum (-1)^k (pi/2)^(2k+1) x^(4k+3) / ((2k+1)! (4k+3))
        // C = sum (-1)^k (pi/2)^(2k) x^(4k+1) / ((2k)! (4k+1))
        let hp = consts().pi.mul_f64(0.5);
        let t = hp * DD::new(x) * DD::new(x); // (pi/2) x^2
        let mut p = DD::new(x); // (pi x^2/2)^m x / m!
        let mut s = DD::ZERO;
        let mut c = DD::ZERO;
        let mut m = 0usize;
        loop {
            let term = p.div_f64((2 * m + 1) as f64);
            let signed = if (m / 2) % 2 == 0 { term } else { -term };
            if m % 2 == 0 {
                c = c + signed;
            } else {
                s = s + signed;
            }
            m += 1;
            p = (p * t).div_f64(m as f64);
            if p.hi.abs() < 1e-34 && m > 4 {
                break;
            }
            if m > 400 {
                break;
            }
        }
        return (s.to_f64(), c.to_f64());
    }
    // C + iS = (1+i)/2 (1 - erfc(w)), w = (1-i) sqrt(pi) x / 2, with
    // erfc(w) = e^{-w^2} CF(w) / sqrt(pi) and e^{-w^2} = e^{i pi x^2/2}
    let sp = consts().sqrt_pi.to_f64();
    let w = C64::new(1.0, -1.0) * (sp * x / 2.0);
    let depth = 80 + (2000.0 / w.norm_sqr()) as i64;
    let mut f = C64::new(0.0, 0.0);
    for k in (1..=depth).rev() {
        f = (k as f64 / 2.0) / (w + f);
    }
    let cf = 1.0 / (sp * (w + f));
    let x2 = DD::new(x) * DD::new(x);
    let q = (x2.hi / 4.0).floor() * 4.0;
    let r = (x2.add_f64(-q)).to_f64();
    let ph = C64::new(super::cospi(r / 2.0), super::sinpi(r / 2.0));
    let v = C64::new(0.5, 0.5) * (1.0 - ph * cf);
    (v.im, v.re)
}

/// Complex Fresnel integrals by power series (moderate `|z|`).
pub fn fresnel_c(z: C64) -> Option<(C64, C64)> {
    if is_real_c(z) {
        let (s, c) = fresnel_real(z.re);
        return Some((cr(s), cr(c)));
    }
    if z.norm() > 6.0 {
        return None;
    }
    let t = z * z * (PI / 2.0);
    let mut p = z;
    let mut s = C64::new(0.0, 0.0);
    let mut c = C64::new(0.0, 0.0);
    for m in 0..1000usize {
        let term = p / (2 * m + 1) as f64;
        let signed = if (m / 2) % 2 == 0 { term } else { -term };
        if m % 2 == 0 {
            c += signed;
        } else {
            s += signed;
        }
        p = p * t / (m + 1) as f64;
        if p.norm() < 1e-18 * (s.norm() + c.norm()) && m > 4 {
            break;
        }
    }
    Some((s, c))
}

/// Exponential integral `E1(x)`, real `x > 0`.
pub fn e1_real(x: f64) -> Option<f64> {
    if x <= 0.0 {
        return None;
    }
    if x <= 1.0 {
        // -gamma - ln x - sum (-x)^k / (k k!)
        let mut sum = DD::ZERO;
        let mut p = 1.0;
        for k in 1..60 {
            p *= -x / k as f64;
            let t = p / k as f64;
            sum = sum.add_f64(t);
            if t.abs() < 1e-18 {
                break;
            }
        }
        return Some((-consts().euler - DD::new(x).ln() - sum).to_f64());
    }
    Some(en_cf(1.0, x))
}

/// Continued fraction for `E_n(x)`, `x > 0` (modified Lentz).
fn en_cf(n: f64, x: f64) -> f64 {
    let tiny = 1e-300;
    let mut b = x + n;
    let mut c = 1.0 / tiny;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..10000 {
        let a = -(i as f64) * (n - 1.0 + i as f64);
        b += 2.0;
        d = 1.0 / (a * d + b);
        c = b + a / c;
        let del = c * d;
        h *= del;
        if (del - 1.0).abs() < 1e-17 {
            break;
        }
    }
    h * (-x).exp()
}

/// `E1(z)` complex.
pub fn e1_c(z: C64) -> Option<C64> {
    if is_real_c(z) && z.re > 0.0 {
        return e1_real(z.re).map(cr);
    }
    if z == cr(0.0) {
        return None;
    }
    if z.norm() <= 2.0 || (z.re < 0.0 && z.im.abs() < 1.0 && z.norm() < 40.0) {
        let mut sum = C64::new(0.0, 0.0);
        let mut p = C64::new(1.0, 0.0);
        for k in 1..400 {
            p = -p * z / k as f64;
            let t = p / k as f64;
            sum += t;
            if t.norm() < 1e-18 * sum.norm().max(1e-300) {
                break;
            }
        }
        let mut r = -EULER - z.ln() - sum;
        if is_real_c(z) && z.re < 0.0 {
            // value on the cut (mpmath: -Ei(-x) - i pi)
            r = C64::new(r.re, -PI);
        }
        return Some(r);
    }
    // continued fraction (Lentz), complex
    let tiny = C64::new(1e-300, 0.0);
    let mut b = z + 1.0;
    let mut c = C64::new(1e300, 0.0);
    let mut d = b.inv();
    let mut h = d;
    for i in 1..20000 {
        let a = -(i as f64) * (i as f64);
        b += 2.0;
        d = a * d + b;
        if d.norm() < 1e-300 {
            d = tiny;
        }
        d = d.inv();
        c = b + a / c;
        if c.norm() < 1e-300 {
            c = tiny;
        }
        let del = c * d;
        h *= del;
        if (del - 1.0).norm() < 1e-16 {
            break;
        }
    }
    Some(h * (-z).exp())
}

/// Generalised exponential integral `E_nu(z)`.
pub fn expint_c(nu: C64, z: C64) -> Option<C64> {
    if !is_real_c(nu) {
        return None;
    }
    let n = nu.re;
    if n == 1.0 {
        return e1_c(z);
    }
    if is_real_c(z) && z.re > 0.0 {
        let x = z.re;
        if n == 0.0 {
            return Some(cr((-x).exp() / x));
        }
        if x > 1.0 {
            return Some(cr(en_cf(n, x)));
        }
        if n == n.floor() && n > 1.0 {
            // series
            let ni = n as i64;
            let mut psi = -EULER;
            for k in 1..ni {
                psi += 1.0 / k as f64;
            }
            let mut fact = 1.0;
            for k in 1..ni {
                fact *= k as f64;
            }
            let mut r = (-x).powi(ni as i32 - 1) / fact * (-x.ln() + psi);
            let mut p = 1.0; // (-x)^k / k!
            for k in 0..80i64 {
                if k > 0 {
                    p *= -x / k as f64;
                }
                if k != ni - 1 {
                    r -= p / (k - ni + 1) as f64;
                }
                if p.abs() < 1e-18 && k > ni {
                    break;
                }
            }
            return Some(cr(r));
        }
        // x^(nu-1) Gamma(1-nu, x)
        let g = super::gamma::uppergamma_real(1.0 - n, x)?;
        return Some(cr(x.powf(n - 1.0) * g));
    }
    None
}

/// `Ei(z)`.
pub fn ei_c(z: C64) -> Option<C64> {
    if is_real_c(z) {
        let x = z.re;
        if x == 0.0 {
            return None;
        }
        if x < 0.0 {
            return e1_real(-x).map(|v| cr(-v));
        }
        if x <= 40.0 {
            let mut sum = DD::ZERO;
            let mut p = 1.0;
            for k in 1..200 {
                p *= x / k as f64;
                let t = p / k as f64;
                sum = sum.add_f64(t);
                if t < 1e-18 * sum.hi {
                    break;
                }
            }
            return Some(cr((consts().euler + DD::new(x).ln() + sum).to_f64()));
        }
        // asymptotic e^x/x sum k!/x^k
        let mut s = 1.0;
        let mut t = 1.0;
        for k in 1..60 {
            let nt = t * k as f64 / x;
            if nt > t {
                break;
            }
            t = nt;
            s += t;
            if t < 1e-17 {
                break;
            }
        }
        return Some(cr((DD::new(x).exp().to_f64() / x) * s));
    }
    if z.re < 0.0 && z.norm() > 2.0 {
        let e = e1_c(-z)?;
        let s = if z.im > 0.0 { PI } else { -PI };
        return Some(-e + C64::new(0.0, s));
    }
    let mut sum = C64::new(0.0, 0.0);
    let mut p = C64::new(1.0, 0.0);
    for k in 1..400 {
        p = p * z / k as f64;
        let t = p / k as f64;
        sum += t;
        if t.norm() < 1e-18 * sum.norm() {
            break;
        }
    }
    Some(EULER + z.ln() + sum)
}

/// `li(z) = Ei(ln z)`.
pub fn li_c(z: C64) -> Option<C64> {
    if z == cr(0.0) {
        return Some(cr(0.0));
    }
    if is_real_c(z) && z.re > 0.0 {
        let l = DD::new(z.re).ln().to_f64();
        return ei_c(cr(l));
    }
    ei_c(z.ln())
}

/// `(Si(x), Ci(x))` for real `x > 0`.
fn sici_pos(x: f64) -> (f64, f64) {
    if x <= 4.0 {
        let x2 = x * x;
        let mut p = DD::new(x); // x^m / m!
        let mut si = DD::ZERO;
        let mut ci = DD::ZERO;
        let mut m = 1usize;
        loop {
            // m odd -> Si term, m even -> Ci term
            let term = p.div_f64(m as f64);
            let k = m / 2;
            let signed = if k % 2 == 0 { term } else { -term };
            if m % 2 == 1 {
                si = si + signed;
            } else {
                ci = ci + signed;
            }
            m += 1;
            p = p.mul_f64(x).div_f64(m as f64);
            if p.hi < 1e-34 {
                break;
            }
        }
        let _ = x2;
        let ci = consts().euler + DD::new(x).ln() + ci;
        return (si.to_f64(), ci.to_f64());
    }
    // E1(ix) = -Ci(x) + i (Si(x) - pi/2)
    let e = e1_c(C64::new(0.0, x)).unwrap();
    (PI / 2.0 + e.im, -e.re)
}

/// `(Shi(x), Chi(x))` for real `x > 0`.
fn shichi_pos(x: f64) -> (f64, f64) {
    if x <= 40.0 {
        let mut p = DD::new(x);
        let mut shi = DD::ZERO;
        let mut chi = DD::ZERO;
        let mut m = 1usize;
        loop {
            let term = p.div_f64(m as f64);
            if m % 2 == 1 {
                shi = shi + term;
            } else {
                chi = chi + term;
            }
            m += 1;
            p = p.mul_f64(x).div_f64(m as f64);
            if p.hi < 1e-34 * (shi.hi + chi.hi.abs()) {
                break;
            }
        }
        let chi = consts().euler + DD::new(x).ln() + chi;
        return (shi.to_f64(), chi.to_f64());
    }
    let ei = ei_c(cr(x)).unwrap().re;
    let e1 = e1_real(x).unwrap();
    ((ei + e1) / 2.0, (ei - e1) / 2.0)
}

/// Trigonometric/hyperbolic integrals: `kind` 0 = Si, 1 = Ci, 2 = Shi, 3 = Chi.
pub fn trig_integral(kind: u8, z: C64) -> Option<C64> {
    if is_real_c(z) {
        let x = z.re;
        if x == 0.0 {
            return if kind % 2 == 0 { Some(cr(0.0)) } else { None };
        }
        let ax = x.abs();
        let (odd, even) = if kind < 2 { sici_pos(ax) } else { shichi_pos(ax) };
        return Some(if kind % 2 == 0 {
            cr(if x < 0.0 { -odd } else { odd })
        } else if x < 0.0 {
            C64::new(even, PI)
        } else {
            cr(even)
        });
    }
    if z.norm() > 10.0 {
        return None;
    }
    // complex power series
    let mut p = z; // z^m/m!
    let mut odd = C64::new(0.0, 0.0);
    let mut even = C64::new(0.0, 0.0);
    let hyper = kind >= 2;
    for m in 1..1000usize {
        let term = p / m as f64;
        let signed = if !hyper && (m / 2) % 2 == 1 { -term } else { term };
        if m % 2 == 1 {
            odd += signed;
        } else {
            even += signed;
        }
        p = p * z / (m + 1) as f64;
        if p.norm() < 1e-18 * (odd.norm() + even.norm()) && m > 3 {
            break;
        }
    }
    Some(if kind % 2 == 0 { odd } else { EULER + z.ln() + even })
}
