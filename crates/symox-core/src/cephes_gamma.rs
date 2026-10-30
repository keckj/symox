//! SciPy's float64 `scipy.special.gamma`: the cephes `Gamma` of xsf
//! (`xsf/cephes/gamma.h`, SciPy 1.15), with the same coefficient tables,
//! branches and operation order, and the C library's `exp`, `pow` and `sin`,
//! so that it gives SciPy's value bit for bit. It is used only where SymPy's
//! `lambdify` calls `scipy.special.gamma`; `gamma_real` stays the float
//! gamma everywhere else.

mod libm {
    extern "C" {
        pub fn exp(x: f64) -> f64;
        pub fn pow(x: f64, y: f64) -> f64;
        pub fn sin(x: f64) -> f64;
    }
}

const P: [f64; 7] = [
    1.60119522476751861407E-4,
    1.19135147006586384913E-3,
    1.04213797561761569935E-2,
    4.76367800457137231464E-2,
    2.07448227648435975150E-1,
    4.94214826801497100753E-1,
    9.99999999999999996796E-1,
];
const Q: [f64; 8] = [
    -2.31581873324120129819E-5,
    5.39605580493303397842E-4,
    -4.45641913851797240494E-3,
    1.18139785222060435552E-2,
    3.58236398605498653373E-2,
    -2.34591795718243348568E-1,
    7.14304917030273074085E-2,
    1.00000000000000000320E0,
];
const STIR: [f64; 5] = [
    7.87311395793093628397E-4,
    -2.29549961613378126380E-4,
    -2.68132617805781232825E-3,
    3.47222221605458667310E-3,
    8.33333333333482257126E-2,
];
const MAXGAM: f64 = 171.624376956302725;
const MAXSTIR: f64 = 143.01608;
const SQRT2PI: f64 = 2.50662827463100050242E0;
const PI: f64 = std::f64::consts::PI;

#[inline]
fn polevl(x: f64, c: &[f64]) -> f64 {
    let mut a = c[0];
    for &v in &c[1..] {
        a = a * x + v;
    }
    a
}

/// Stirling's formula, for x > 33.
fn stirf(x: f64) -> f64 {
    if x >= MAXGAM {
        return f64::INFINITY;
    }
    let w = 1.0 / x;
    let w = 1.0 + w * polevl(w, &STIR);
    let mut y = unsafe { libm::exp(x) };
    if x > MAXSTIR {
        // avoid overflow in pow
        let v = unsafe { libm::pow(x, 0.5 * x - 0.25) };
        y = v * (v / y);
    } else {
        y = unsafe { libm::pow(x, x - 0.5) } / y;
    }
    SQRT2PI * y * w
}

/// `scipy.special.gamma` of a float64, bit for bit: `nan` at the negative
/// integers, -inf and nan; `copysign(inf, x)` at zero; inf on overflow.
// The small-argument branch keeps cephes' literal 0.5772156649015329 (the
// correctly rounded Euler constant) as written in the reference source.
#[allow(clippy::approx_constant)]
pub fn gamma(x: f64) -> f64 {
    if !x.is_finite() {
        return if x > 0.0 { x } else { f64::NAN };
    }
    if x == 0.0 {
        return f64::INFINITY.copysign(x);
    }
    let q = x.abs();
    if q > 33.0 {
        let mut sign = 1.0;
        let z = if x < 0.0 {
            let mut p = q.floor();
            if p == q {
                return f64::NAN;
            }
            // C's `int i = p` on x86-64: out of range values give INT_MIN
            let i = if p < 2147483648.0 { p as i32 } else { i32::MIN };
            if i & 1 == 0 {
                sign = -1.0;
            }
            let mut z = q - p;
            if z > 0.5 {
                p += 1.0;
                z = q - p;
            }
            let z = q * unsafe { libm::sin(PI * z) };
            if z == 0.0 {
                return sign * f64::INFINITY;
            }
            PI / (z.abs() * stirf(q))
        } else {
            stirf(x)
        };
        return sign * z;
    }
    let mut x = x;
    let mut z = 1.0;
    while x >= 3.0 {
        x -= 1.0;
        z *= x;
    }
    let mut small = false;
    while x < 0.0 {
        if x > -1.0E-9 {
            small = true;
            break;
        }
        z /= x;
        x += 1.0;
    }
    if !small {
        while x < 2.0 {
            if x < 1.0E-9 {
                small = true;
                break;
            }
            z /= x;
            x += 1.0;
        }
    }
    if small {
        return if x == 0.0 { f64::NAN } else { z / ((1.0 + 0.5772156649015329 * x) * x) };
    }
    if x == 2.0 {
        return z;
    }
    let x = x - 2.0;
    z * polevl(x, &P) / polevl(x, &Q)
}

#[cfg(test)]
mod tests {
    use super::gamma;

    #[test]
    fn scipy_values() {
        // scipy.special.gamma (SciPy 1.15.3), as bit patterns
        let cases: [(f64, u64); 10] = [
            (0.5, 0x3ffc5bf891b4ef6a),
            (3.0, 0x4000000000000000),
            (-0.5, 0xc00c5bf891b4ef6a),
            (171.5, 0x7fe0e1863dcad78b),
            (-40.5, 0xb5ec357aa9ad5056),
            (1e-10, 0x4202a05f1ffb61dd),
            (0.0, f64::INFINITY.to_bits()),
            (-0.0, f64::NEG_INFINITY.to_bits()),
            (-3.0, f64::NAN.to_bits()),
            (f64::NEG_INFINITY, f64::NAN.to_bits()),
        ];
        for (x, bits) in cases {
            assert_eq!(gamma(x).to_bits(), bits, "{x}");
        }
        assert_eq!(gamma(172.0), f64::INFINITY);
    }
}
