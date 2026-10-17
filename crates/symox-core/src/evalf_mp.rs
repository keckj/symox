//! Arbitrary precision numerical evaluation of expressions.

use num_bigint::BigInt;

use crate::bigfloat::{ln2, pi, BigFloat};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::Number;

/// A complex number with arbitrary precision parts.
#[derive(Clone, Debug)]
pub struct MpComplex {
    pub re: BigFloat,
    pub im: BigFloat,
}

impl MpComplex {
    pub fn real(re: BigFloat) -> MpComplex {
        MpComplex { re, im: BigFloat::zero() }
    }

    pub fn is_real(&self) -> bool {
        self.im.is_zero()
    }

    pub(crate) fn add(&self, o: &MpComplex, p: u32) -> MpComplex {
        MpComplex { re: self.re.add(&o.re, p), im: self.im.add(&o.im, p) }
    }

    pub(crate) fn mul(&self, o: &MpComplex, p: u32) -> MpComplex {
        if self.is_real() && o.is_real() {
            return MpComplex::real(self.re.mul(&o.re, p));
        }
        let re = self.re.mul(&o.re, p).sub(&self.im.mul(&o.im, p), p);
        let im = self.re.mul(&o.im, p).add(&self.im.mul(&o.re, p), p);
        MpComplex { re, im }
    }

    pub(crate) fn scale(&self, k: &BigFloat, p: u32) -> MpComplex {
        MpComplex { re: self.re.mul(k, p), im: self.im.mul(k, p) }
    }

    pub(crate) fn inv(&self, p: u32) -> Option<MpComplex> {
        if self.is_real() {
            return Some(MpComplex::real(BigFloat::one().div(&self.re, p)?));
        }
        let d = self.re.mul(&self.re, p + 10).add(&self.im.mul(&self.im, p + 10), p + 10);
        Some(MpComplex { re: self.re.div(&d, p)?, im: self.im.neg().div(&d, p)? })
    }

    pub(crate) fn abs(&self, p: u32) -> BigFloat {
        if self.is_real() {
            return self.re.abs();
        }
        self.re
            .mul(&self.re, p + 10)
            .add(&self.im.mul(&self.im, p + 10), p + 10)
            .sqrt(p)
            .unwrap()
    }

    pub(crate) fn neg(&self) -> MpComplex {
        MpComplex { re: self.re.neg(), im: self.im.neg() }
    }

    pub(crate) fn sub(&self, other: &MpComplex, precision: u32) -> MpComplex {
        self.add(&other.neg(), precision)
    }

    pub(crate) fn div(&self, other: &MpComplex, precision: u32) -> Option<MpComplex> {
        Some(self.mul(&other.inv(precision + 8)?, precision))
    }

    pub(crate) fn conjugate(&self) -> MpComplex {
        MpComplex { re: self.re.clone(), im: self.im.neg() }
    }

    pub(crate) fn round_to(&self, precision: u32) -> MpComplex {
        MpComplex { re: self.re.round_to(precision), im: self.im.round_to(precision) }
    }

    fn arg(&self, p: u32) -> BigFloat {
        BigFloat::atan2(&self.im, &self.re, p)
    }

    fn exp(&self, p: u32) -> Option<MpComplex> {
        // Keep exponential results inside a decimal exponent magnitude of
        // one million, matching the native decimal constructor's bound.
        if self.re.abs().cmp(&BigFloat::from_i64(2_000_000)) == std::cmp::Ordering::Greater { return None; }
        let r = self.re.checked_exp(p + 10)?;
        if self.im.is_zero() {
            return Some(MpComplex::real(r.round_to(p)));
        }
        let (s, c) = self.im.sin_cos(p + 10);
        Some(MpComplex { re: r.mul(&c, p), im: r.mul(&s, p) })
    }

    fn ln(&self, p: u32) -> Option<MpComplex> {
        if self.is_real() && self.re.sign() > 0 {
            return Some(MpComplex::real(self.re.ln(p)?));
        }
        let a = self.abs(p + 10);
        Some(MpComplex { re: a.ln(p)?, im: self.arg(p) })
    }

    fn pow(&self, e: &MpComplex, p: u32) -> Option<MpComplex> {
        if self.re.is_zero() && self.im.is_zero() {
            return if e.re.sign() > 0 { Some(MpComplex::real(BigFloat::zero())) } else { None };
        }
        if self.is_real() && e.is_real() && self.re.sign() > 0 {
            return Some(MpComplex::real(self.re.powf(&e.re, p)?));
        }
        let wp = p + 20 + e.re.mag().max(0) as u32 + e.im.mag().max(0) as u32;
        let l = self.ln(wp)?;
        l.mul(e, wp).exp(p)
    }

    fn powi(&self, n: i64, p: u32) -> Option<MpComplex> {
        if self.is_real() {
            return Some(MpComplex::real(self.re.powi(n, p)?));
        }
        let wp = p + 10 + 64 - n.unsigned_abs().leading_zeros();
        let mut base = self.clone();
        let mut k = n.unsigned_abs();
        let mut acc = MpComplex::real(BigFloat::one());
        while k > 0 {
            if k & 1 == 1 {
                acc = acc.mul(&base, wp);
                if !supported_magnitude(&acc) { return None; }
            }
            k >>= 1;
            if k > 0 {
                base = base.mul(&base, wp);
                if !supported_magnitude(&base) { return None; }
            }
        }
        if n < 0 {
            acc.inv(p)
        } else {
            Some(acc)
        }
    }
}

fn supported_magnitude(value: &MpComplex) -> bool {
    [&value.re, &value.im].iter().all(|v| v.is_zero() || v.mag().unsigned_abs() <= 3_321_930)
}

fn number(n: &Number, p: u32) -> BigFloat {
    match n {
        Number::Int(i) => BigFloat::from_i64(*i),
        Number::BigInt(b) => BigFloat::from_int(b),
        Number::Rational(r) => BigFloat::from_ratio(r.numer(), r.denom(), p),
        Number::Float(f) => BigFloat::from_f64(*f),
        Number::MpFloat(value, _) => value.round_to(p),
    }
}

/// Highest declared precision among the floating point leaves.
pub fn max_precision(e: &Expr) -> Option<u32> {
    e.float_info().max_precision
}

pub fn has_mp(e: &Expr) -> bool {
    e.float_info().has_mp
}

/// Evaluate a numeric expression to `prec` bits. Returns `None` when the
/// expression has free symbols or uses an unsupported function.
pub fn eval_mp(e: &Expr, prec: u32) -> Option<MpComplex> {
    if let Kind::Func(Func::MatrixDotProduct(a,b),_)=e.kind(){return eval_mp(&crate::matrixcontraction::dot_product(a,b,true).ok()?,prec);}
    if matches!(e.kind(),Kind::Func(Func::Subs | Func::Limit | Func::Integral | Func::Sum,_)){
        let evaluated=crate::functions::doit(e,true);
        if evaluated==*e{return None;}
        return eval_mp(&evaluated,prec);
    }
    let wp = prec + 20;
    let r = ev(e, wp)?;
    if !supported_magnitude(&r) { return None; }
    Some(MpComplex { re: r.re.round_to(prec), im: r.im.round_to(prec) })
}

/// Convert a computed value back to native expression numbers, retaining the
/// declared precision. Exact zero components stay exact zero.
pub fn to_expr(value: MpComplex, prec: u32) -> Expr {
    let re = if value.re.is_zero() { Expr::zero() } else { Expr::num(Number::mp_float(value.re, prec)) };
    if value.im.is_zero() { return re; }
    let im = Expr::from_kind(Kind::Mul(Number::mp_float(value.im, prec), vec![(Expr::i(), Expr::one())]));
    if re.is_zero() { im } else { crate::ops::add2(&re, &im) }
}

/// Numerical evaluation with native arbitrary precision numbers. Unsupported
/// numerical functions remain symbolic; this never silently computes them at
/// machine precision. Integer powers and symbolic function arguments remain
/// exact, as in SymPy's partial evalf.
pub fn evalf(e: &Expr, prec: u32) -> Expr {
    if e.is_zero() { return Expr::zero(); }
    if let Some(value) = eval_mp(e, prec) { return to_expr(value, prec); }
    if crate::functions::is_numeric_value(e) { return e.clone(); }
    match e.kind() {
        Kind::Add(..) => crate::ops::add(&e.args().iter().map(|a| evalf(a, prec)).collect::<Vec<_>>()),
        Kind::Mul(..) => crate::ops::mul(&e.args().iter().map(|a| evalf(a, prec)).collect::<Vec<_>>()),
        Kind::Pow(base, exponent) => crate::ops::pow(
            &evalf(base, prec),
            &if exponent.is_integer_number() { exponent.clone() } else { evalf(exponent, prec) },
        ),
        Kind::Rel(op, a, b) => crate::relational::rel(*op, &evalf(a, prec), &evalf(b, prec)),
        Kind::Piecewise(pairs) => crate::relational::piecewise(
            pairs.iter().map(|(value, condition)| (evalf(value, prec), condition.clone())).collect(),
        ),
        _ => e.clone(),
    }
}

pub(crate) fn ev(e: &Expr, p: u32) -> Option<MpComplex> {
    Some(match e.kind() {
        Kind::Num(n) => MpComplex::real(number(n, p)),
        Kind::Const(c) => match c {
            Constant::Pi => MpComplex::real(pi(p)),
            Constant::E => MpComplex::real(BigFloat::one().exp(p)),
            Constant::I => MpComplex { re: BigFloat::zero(), im: BigFloat::one() },
            Constant::EulerGamma => MpComplex::real(euler_gamma(p)),
            Constant::Catalan => MpComplex::real(catalan(p)),
            Constant::GoldenRatio => MpComplex::real(golden_ratio(p)),
            Constant::TribonacciConstant => MpComplex::real(tribonacci_constant(p)),
            _ => return None,
        },
        Kind::Add(c, terms) => eval_sum(c, terms, p)?,
        Kind::Mul(c, facs) => {
            let mut acc = MpComplex::real(number(c, p));
            for (b, x) in facs {
                let v = pow_of(b, x, p)?;
                acc = acc.mul(&v, p);
            }
            acc
        }
        Kind::Pow(b, x) => pow_of(b, x, p)?,
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>ev(&crate::matrixcontraction::dot_product(a,b,true).ok()?,p)?,
        Kind::Func(Func::Subs | Func::Limit | Func::Integral | Func::Sum,_)=>{let evaluated=crate::functions::doit(e,true);if evaluated==*e{return None;}ev(&evaluated,p)?},
        Kind::Func(f, args) => func(f, args, p)?,
        _ => return None,
    })
}

/// Increase working precision when large terms nearly cancel. A fixed guard
/// budget would incorrectly turn sqrt(10**100 + 1) - 10**50 into exact zero.
fn eval_sum(constant: &Number, terms: &[(Expr, Number)], precision: u32) -> Option<MpComplex> {
    let mut working = precision + 32;
    loop {
        let mut sum = MpComplex::real(number(constant, working));
        let mut largest = sum.re.mag();
        for (term, coefficient) in terms {
            let value = ev(term, working)?;
            let value = if coefficient.is_one() { value } else { value.scale(&number(coefficient, working), working) };
            largest = largest.max(value.re.mag()).max(value.im.mag());
            sum = sum.add(&value, working);
        }
        let magnitude = sum.re.mag().max(sum.im.mag());
        let cancellation = largest.saturating_sub(magnitude).max(0) as u64;
        let required = precision as u64 + cancellation + 24;
        if required <= working as u64 { return Some(sum); }
        if working >= precision + 4096 { return None; }
        working = if magnitude < i64::MIN / 8 {
            (working + (working - precision).max(32)).min(precision + 4096)
        } else {
            (required.min((precision + 4096) as u64)) as u32
        };
    }
}

fn pow_of(b: &Expr, x: &Expr, p: u32) -> Option<MpComplex> {
    let bv = ev(b, p)?;
    if let Some(k) = x.as_i64() {
        return bv.powi(k, p);
    }
    if *x == Expr::rational(1, 2) && bv.is_real() && bv.re.sign() >= 0 {
        return Some(MpComplex::real(bv.re.sqrt(p)?));
    }
    let xv = ev(x, p)?;
    bv.pow(&xv, p)
}

fn func(f: &Func, args: &[Expr], p: u32) -> Option<MpComplex> {
    let wp = p + 10;
    if let Func::Special(s) = f {
        return crate::special::mp::eval(*s, args, p);
    }
    if matches!(f, Func::Undefined(name) if name.is_internal("UnevaluatedExpr")) {
        return ev(args.first()?, p);
    }
    if let Func::Atan2 = f {
        let (y, x) = (ev(&args[0], wp)?, ev(&args[1], wp)?);
        if !(y.is_real() && x.is_real()) {
            return None;
        }
        return Some(MpComplex::real(BigFloat::atan2(&y.re, &x.re, p)));
    }
    let a = ev(args.first()?, wp)?;
    let real_only = |v: Option<BigFloat>| v.map(MpComplex::real);
    Some(match f {
        Func::Exp => a.exp(p)?,
        Func::Log => a.ln(p)?,
        Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc if a.is_real() => {
            let (s, c) = a.re.sin_cos(wp);
            MpComplex::real(match f {
                Func::Sin => s,
                Func::Cos => c,
                Func::Tan => s.div(&c, p)?,
                Func::Cot => c.div(&s, p)?,
                Func::Sec => BigFloat::one().div(&c, p)?,
                _ => BigFloat::one().div(&s, p)?,
            })
        }
        Func::Sin | Func::Cos | Func::Tan => {
            // complex: via exponentials
            let iz = MpComplex { re: a.im.neg(), im: a.re.clone() };
            let e1 = iz.exp(wp)?;
            let e2 = e1.inv(wp)?;
            let half = BigFloat { man: BigInt::from(1), exp: -1 };
            let cosv = e1.add(&e2, wp).scale(&half, wp);
            let diff = e1.add(&e2.scale(&BigFloat::from_i64(-1), wp), wp);
            // sin = (e1 - e2) / (2i)
            let sinv = MpComplex { re: diff.im.mul(&half, wp), im: diff.re.neg().mul(&half, wp) };
            match f {
                Func::Sin => sinv,
                Func::Cos => cosv,
                _ => sinv.mul(&cosv.inv(wp)?, p),
            }
        }
        Func::Asin if a.is_real() => real_only(a.re.asin(p))?,
        Func::Acos if a.is_real() => real_only(a.re.acos(p))?,
        Func::Atan if a.is_real() => MpComplex::real(a.re.atan(p)),
        Func::Acot if a.is_real() => MpComplex::real(BigFloat::one().div(&a.re, wp)?.atan(p)),
        Func::Sinh | Func::Cosh | Func::Tanh if a.is_real() => {
            if a.re.abs().cmp(&BigFloat::from_i64(2_000_000)) == std::cmp::Ordering::Greater { return None; }
            let e1 = a.re.checked_exp(wp)?;
            let e2 = BigFloat::one().div(&e1, wp)?;
            let sh = e1.sub(&e2, wp).ldexp(-1);
            let ch = e1.add(&e2, wp).ldexp(-1);
            MpComplex::real(match f {
                Func::Sinh => sh,
                Func::Cosh => ch,
                _ => sh.div(&ch, p)?,
            })
        }
        Func::Asinh if a.is_real() => {
            // ln(x + sqrt(x^2+1))
            let x = &a.re;
            let s = x.mul(x, wp).add(&BigFloat::one(), wp).sqrt(wp)?;
            MpComplex::real(x.add(&s, wp).ln(p)?)
        }
        Func::Acosh if a.is_real() => {
            let x = &a.re;
            let s = x.mul(x, wp).sub(&BigFloat::one(), wp).sqrt(wp)?;
            MpComplex::real(x.add(&s, wp).ln(p)?)
        }
        Func::Atanh if a.is_real() => {
            let x = &a.re;
            let one = BigFloat::one();
            let q = one.add(x, wp).div(&one.sub(x, wp), wp)?;
            MpComplex::real(q.ln(wp)?.ldexp(-1))
        }
        Func::Abs => MpComplex::real(a.abs(p)),
        Func::Re => MpComplex::real(a.re),
        Func::Im => MpComplex::real(a.im),
        Func::Conjugate => MpComplex { re: a.re, im: a.im.neg() },
        Func::Arg => MpComplex::real(a.arg(p)),
        Func::Sign if a.is_real() => MpComplex::real(BigFloat::from_i64(a.re.sign() as i64)),
        Func::Floor | Func::Ceiling if a.is_real() => {
            let t = a.re.to_int_trunc();
            let tf = BigFloat::from_int(&t);
            let c = a.re.cmp(&tf);
            let v = match (f, c) {
                (Func::Floor, std::cmp::Ordering::Less) => t - 1,
                (Func::Ceiling, std::cmp::Ordering::Greater) => t + 1,
                _ => t,
            };
            MpComplex::real(BigFloat::from_int(&v))
        }
        Func::Erf if a.is_real() => MpComplex::real(erf_mp(&a.re, p)?),
        Func::Erfc if a.is_real() => MpComplex::real(BigFloat::one().sub(&erf_mp(&a.re, wp)?, p)),
        Func::Gamma if a.is_real() => MpComplex::real(gamma_mp(&a.re, p)?),
        _ => return None,
    })
}

/// Euler–Mascheroni constant (Brent–McMillan would be faster; this uses
/// the simple exponential integral formula with fixed point arithmetic).
pub fn euler_gamma(prec: u32) -> BigFloat {
    // gamma = H - ln(n) - corrections; use Brent-McMillan B1:
    // with n ~ prec*ln2/4, A = sum_{k} (n^k/k!)^2 (H_k - ln n), B = sum (n^k/k!)^2
    let wp = prec + 40;
    let n = ((wp as f64) * std::f64::consts::LN_2 / 4.0).ceil() as u64 + 1;
    let one = BigInt::from(1) << (wp as usize);
    let ln_n = BigFloat::from_i64(n as i64).ln(wp).unwrap();
    let ln_n_fx = {
        let s = ln_n.exp + wp as i64;
        if s >= 0 {
            &ln_n.man << (s as usize)
        } else {
            &ln_n.man >> ((-s) as usize)
        }
    };
    let n2 = BigInt::from(n * n);
    let mut a = -&ln_n_fx;
    let mut b = one.clone();
    let mut u = a.clone();
    let mut v = one.clone();
    let mut k = 1u64;
    loop {
        let k2 = BigInt::from(k * k);
        u = (&u * &n2) / &k2;
        v = (&v * &n2) / &k2;
        u = &u + (&v / BigInt::from(k));
        if v.bits() == 0 && u.bits() == 0 {
            break;
        }
        a += &u;
        b += &v;
        k += 1;
        if k > 40 * n + 100 {
            break;
        }
    }
    let r = BigFloat::from_int(&a).div(&BigFloat::from_int(&b), prec).unwrap();
    let _ = ln2;
    r
}

/// Catalan's constant with Ramanujan's series
/// `G = pi/8 log(2 + sqrt(3)) + 3/8 sum_n 1/((2n+1)**2 binomial(2n, n))`.
pub fn catalan(prec: u32) -> BigFloat {
    let wp = prec + 30;
    let one = BigInt::from(1) << (wp as usize);
    // r = 1/binomial(2n, n) in fixed point
    let mut r = one.clone();
    let mut sum = BigInt::from(0);
    let mut n = 0u64;
    while r.bits() > 0 {
        let d = BigInt::from(2 * n + 1);
        sum += &r / (&d * &d);
        r = (&r * BigInt::from(n + 1)) / BigInt::from(2 * (2 * n + 1));
        n += 1;
    }
    let series = BigFloat::round(sum, -(wp as i64), wp).mul_i64(3, wp).ldexp(-3);
    let three = BigFloat::from_i64(3).sqrt(wp).unwrap();
    let lg = BigFloat::from_i64(2).add(&three, wp).ln(wp).unwrap();
    pi(wp).mul(&lg, wp).ldexp(-3).add(&series, wp).round_to(prec)
}

/// `GoldenRatio = (1 + sqrt(5))/2`.
pub fn golden_ratio(prec: u32) -> BigFloat {
    let wp = prec + 10;
    BigFloat::from_i64(5).sqrt(wp).unwrap().add(&BigFloat::one(), wp).ldexp(-1).round_to(prec)
}

/// Real root of `x**3 - x**2 - x - 1` (Newton iteration).
pub fn tribonacci_constant(prec: u32) -> BigFloat {
    let wp = prec + 20;
    let mut x = BigFloat::from_f64(1.839_286_755_214_161);
    for _ in 0..64 {
        let x2 = x.mul(&x, wp);
        let x3 = x2.mul(&x, wp);
        let f = x3.sub(&x2, wp).sub(&x, wp).sub(&BigFloat::one(), wp);
        let df = x2.mul_i64(3, wp).sub(&x.mul_i64(2, wp), wp).sub(&BigFloat::one(), wp);
        let step = f.div(&df, wp).unwrap();
        x = x.sub(&step, wp);
        if step.is_zero() || step.mag() < -(wp as i64) {
            break;
        }
    }
    x.round_to(prec)
}

fn erf_mp(x: &BigFloat, prec: u32) -> Option<BigFloat> {
    // erf x = 2/sqrt(pi) * sum (-1)^n x^(2n+1) / (n! (2n+1)); fine for
    // moderate |x| with enough guard bits (cancellation ~ x^2 / ln 2 bits)
    let x2f = x.to_f64() * x.to_f64();
    if x2f > 1e6 {
        return Some(BigFloat::from_i64(x.sign() as i64));
    }
    let wp = prec + 20 + (x2f * 1.5) as u32;
    let x2 = x.mul(x, wp);
    let mut term = x.clone();
    let mut sum = x.clone();
    let mut n = 1i64;
    loop {
        term = term.mul(&x2, wp).neg().div_i64(n, wp);
        let t = term.div_i64(2 * n + 1, wp);
        if t.is_zero() || t.mag() < sum.mag() - wp as i64 - 4 {
            break;
        }
        sum = sum.add(&t, wp);
        n += 1;
    }
    let c = BigFloat::from_i64(2).div(&pi(wp).sqrt(wp)?, wp)?;
    Some(sum.mul(&c, prec))
}

fn gamma_mp(x: &BigFloat, prec: u32) -> Option<BigFloat> {
    // integers: factorial; otherwise Spouge's approximation with a = prec*ln2/ln(2pi)
    let t = x.to_int_trunc();
    if BigFloat::from_int(&t) == *x {
        let n = num_traits::ToPrimitive::to_i64(&t)?;
        if n <= 0 {
            return None;
        }
        let mut acc = BigInt::from(1);
        for i in 2..n {
            acc *= i;
        }
        return Some(BigFloat::from_int(&acc).round_to(prec));
    }
    let wp = prec.saturating_mul(2) + 80;
    if x.sign() < 0 || x.cmp(&BigFloat { man: BigInt::from(1), exp: -1 }) == std::cmp::Ordering::Less {
        // reflection: gamma(x) = pi / (sin(pi x) gamma(1-x))
        let one_minus = BigFloat::one().sub(x, wp);
        let g = gamma_mp(&one_minus, wp)?;
        let p = pi(wp);
        let s = p.mul(x, wp).sin(wp);
        return p.div(&s.mul(&g, wp), prec);
    }
    let a = (((prec + 40) as f64) * std::f64::consts::LN_2 / (2.0 * std::f64::consts::PI).ln()).ceil() as i64;
    let z = x.sub(&BigFloat::one(), wp);
    let af = BigFloat::from_i64(a);
    // c0 = sqrt(2 pi)
    let mut sum = pi(wp).ldexp(1).sqrt(wp)?;
    let mut fact = BigFloat::one(); // (k-1)!
    for k in 1..a {
        let kf = BigFloat::from_i64(k);
        // c_k = (-1)^(k-1) (a-k)^(k-1/2) e^(a-k) / (k-1)!
        let amk = BigFloat::from_i64(a - k);
        let e1 = amk.powf(&BigFloat::from_i64(2 * k - 1).ldexp(-1), wp)?;
        let e2 = amk.checked_exp(wp)?;
        let mut ck = e1.mul(&e2, wp).div(&fact, wp)?;
        if k % 2 == 0 {
            ck = ck.neg();
        }
        sum = sum.add(&ck.div(&z.add(&kf, wp), wp)?, wp);
        fact = fact.mul(&kf, wp);
    }
    // gamma(z+1) = (z+a)^(z+1/2) e^(-(z+a)) * sum
    let za = z.add(&af, wp);
    let half = BigFloat { man: BigInt::from(1), exp: -1 };
    let pw = za.powf(&z.add(&half, wp), wp)?;
    let ex = za.neg().checked_exp(wp)?;
    Some(pw.mul(&ex, wp).mul(&sum, prec))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bigfloat::dps_to_prec;
    use crate::prelude::*;

    fn s(e: &Expr, dps: u32) -> String {
        let v = eval_mp(e, dps_to_prec(dps)).unwrap();
        v.re.to_string_dps(dps, false)
    }

    #[test]
    fn eval() {
        let e = add(&[Expr::pi(), Expr::e()]);
        assert_eq!(s(&e, 30), "5.85987448204883847382293085463");
        assert_eq!(s(&sqrt(&Expr::int(2)), 40), "1.414213562373095048801688724209698078570");
        assert_eq!(s(&exp(&Expr::pi()), 25), "23.14069263277926900572909");
        assert_eq!(s(&Expr::constant(Constant::EulerGamma), 30), "0.577215664901532860606512090082");
        assert_eq!(s(&gamma(&Expr::rational(1, 3)), 25), "2.678938534707747633655693");
        assert_eq!(s(&erf(&Expr::one()), 25), "0.8427007929497148693412206");
    }

    #[test]
    fn precision_survives_native_expression_arithmetic_and_serialization() {
        let p = dps_to_prec(70);
        let value = crate::bigfloat::BigFloat::from_decimal("0.123456789012345678901234567890123456789", p).unwrap();
        let original = Expr::num(Number::mp_float(value.clone(), p));
        assert_eq!(max_precision(&original), Some(p));
        assert_eq!(crate::parser::parse_str(&crate::printer::srepr(&original)).unwrap(), original);
        let increment = crate::parser::parse_str("Float('1e-60', 70)").unwrap();
        let changed = add(&[original.clone(), increment]);
        assert_ne!(changed, original);
        let delta = sub(&changed, &original);
        let result = eval_mp(&delta, p).unwrap();
        assert!((result.re.to_f64() / 1e-60 - 1.0).abs() < 1e-10);
        let neg = Number::mp_float(value.neg(), p);
        assert_eq!(neg.add(original.as_number().unwrap()), Number::ZERO);
    }

    #[test]
    fn evalf_retains_partial_symbolic_precision() {
        let x = Expr::symbol("x");
        let p = dps_to_prec(80);
        let expression = add(&[x.clone(), Expr::pi()]);
        let evaluated = super::evalf(&expression, p);
        assert_eq!(crate::traverse::free_symbols(&evaluated), vec![x]);
        assert_eq!(max_precision(&evaluated), Some(p));
        assert!(crate::printer::to_str(&evaluated).contains("3.141592653589793238462643383279502884197"));
        let raw = Expr::from_kind(Kind::Func(Func::Undefined("unknown_numeric".into()), vec![Expr::pi()]));
        assert_eq!(super::evalf(&raw, p), raw);
    }

    #[test]
    fn cancellation_increases_working_precision() {
        let expression = crate::parser::parse_str("sqrt(10**100 + 1) - 10**50").unwrap();
        let result = eval_mp(&expression, dps_to_prec(60)).unwrap();
        assert!(result.re.sign() > 0);
        assert!(result.re.to_string_dps(50, false).starts_with("5.0000000000000000000000000000000000000000000000000e-51"));
    }

    #[test]
    fn extreme_exponential_inputs_remain_symbolic() {
        let p = dps_to_prec(60);
        let huge = crate::parser::parse_str("Float('1e100', 60)").unwrap();
        for function in [Func::Exp, Func::Sinh, Func::Cosh, Func::Tanh] {
            for input in [huge.clone(), crate::ops::neg(&huge)] {
                let expression = crate::functions::func(function.clone(), vec![input]);
                assert!(matches!(expression.kind(), Kind::Func(..)));
                assert_eq!(super::evalf(&expression, p), expression);
            }
        }
        let exponent = crate::parser::parse_str("Float('1000001', 60)").unwrap();
        let expression = crate::ops::pow(&Expr::int(10), &exponent);
        assert!(matches!(expression.kind(), Kind::Pow(..)));
        assert_eq!(super::evalf(&expression, p), expression);
        let boundary = Expr::from_kind(Kind::Func(Func::Exp, vec![Expr::int(2_000_000)]));
        assert!(eval_mp(&boundary, p).is_some());
        let outside = Expr::from_kind(Kind::Func(Func::Exp, vec![Expr::int(2_000_001)]));
        assert!(eval_mp(&outside, p).is_none());
    }
}
