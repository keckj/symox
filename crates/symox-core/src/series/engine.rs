//! Expansion engine: expression -> truncated Puiseux series.
//!
//! The engine expands an expression in `t` where the variable is
//! `x = x0 + t^D` (or `x = ±t^-D` at infinity). Sub-expressions free of
//! `x` are coefficients. `log(t^D)` is represented by a dedicated symbol
//! (`logt`) which is treated as a constant ("slower than any power").
//!
//! Precision is managed in absolute terms: `ser(e, n)` returns a series
//! correct up to `O(t^n)`; products and powers query the valuations of
//! their operands and recompute them with the extra precision they need.

use rustc_hash::FxHashMap;

use super::coeff::C;
use super::ser::{exp_kernel, log_kernel, pow_kernel, sincos_kernel, tan_kernel, SErr, SRes, Ser};
use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::{self as fx, func};
use crate::number::Number;
use crate::ops::{add, mul, mul2, mul_factor, neg, pow, sub};

/// How far (in t-units beyond the request) we look for a non-zero
/// leading term before giving up.
const MAX_EXTRA: i64 = 64;
const MAX_DEPTH: u32 = 400;

/// Does `e` depend on the symbol `x`?
pub fn depends(e: &Expr, x: &Expr) -> bool {
    if !e.may_contain(x) {
        return false;
    }
    if e == x {
        return true;
    }
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => false,
        Kind::Add(_, terms) => terms.iter().any(|(t, _)| depends(t, x)),
        Kind::Mul(_, facs) => facs.iter().any(|(b, p)| depends(b, x) || depends(p, x)),
        Kind::Pow(b, p) => depends(b, x) || depends(p, x),
        Kind::Func(_, args) => args.iter().any(|a| depends(a, x)),
        _ => e.args().iter().any(|a| depends(a, x)),
    }
}

/// `e = alpha*s + beta` with `alpha` a number and `beta` free of `s`.
pub fn linear_in(e: &Expr, s: &Expr) -> Option<(Number, Expr)> {
    if !depends(e, s) {
        return Some((Number::ZERO, e.clone()));
    }
    if e == s {
        return Some((Number::ONE, Expr::zero()));
    }
    match e.kind() {
        Kind::Mul(c, facs) if facs.len() == 1 && facs[0].0 == *s && facs[0].1.is_one() => {
            Some((c.clone(), Expr::zero()))
        }
        Kind::Add(c, terms) => {
            let mut alpha = Number::ZERO;
            let mut rest = vec![Expr::num(c.clone())];
            for (t, k) in terms {
                if t == s {
                    alpha = alpha.add(k);
                } else if depends(t, s) {
                    return None;
                } else {
                    rest.push(crate::ops::add_term(t, k));
                }
            }
            Some((alpha, add(&rest)))
        }
        _ => None,
    }
}

/// Polynomial degree of `e` in the symbol `s` (`None` if not polynomial).
pub fn degree_in(e: &Expr, s: &Expr) -> Option<i64> {
    if !depends(e, s) {
        return Some(0);
    }
    if e == s {
        return Some(1);
    }
    match e.kind() {
        Kind::Pow(b, p) if b == s => p.as_i64().filter(|k| *k >= 0),
        Kind::Mul(_, facs) => {
            let mut d = 0;
            for (b, p) in facs {
                d += degree_in(&mul_factor(b, p), s)?;
            }
            Some(d)
        }
        Kind::Add(_, terms) => {
            let mut d = 0;
            for (t, _) in terms {
                d = d.max(degree_in(t, s)?);
            }
            Some(d)
        }
        _ => None,
    }
}

fn is_bad_constant(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(c) if c.is_infinite() || *c == Constant::NaN)
}

fn contains_bad(e: &Expr) -> bool {
    if is_bad_constant(e) {
        return true;
    }
    match e.kind() {
        Kind::Num(_) | Kind::Sym(_) | Kind::Const(_) | Kind::Bool(_) => false,
        _ => e.args().iter().any(contains_bad),
    }
}

/// Numeric sign of a coefficient, if it is a real number.
pub fn sign_of(c: &C) -> Option<i32> {
    match c {
        C::N(n) => Some(if n.is_zero() {
            0
        } else if n.is_negative() {
            -1
        } else {
            1
        }),
        C::E(e) => {
            if !crate::functions::is_numeric_value(e) {
                return None;
            }
            let z = crate::evalf::eval_c64(e)?;
            if z.im.abs() > 1e-12 * (1.0 + z.re.abs()) || !z.re.is_finite() {
                return None;
            }
            // too close to zero to be trusted (e.g. sin(1)**2 + cos(1)**2 - 1)
            if z.re.abs() < 1e-12 {
                return None;
            }
            if z.re > 0.0 {
                Some(1)
            } else if z.re < 0.0 {
                Some(-1)
            } else {
                None
            }
        }
    }
}

pub struct Engine {
    pub x: Expr,
    pub xser: Ser,
    /// ramification: `x - x0 = t^D` (or `1/x = t^D`)
    pub d: i64,
    /// symbol standing for `log(t^D)`
    pub logt: Expr,
    memo: FxHashMap<Expr, Ser>,
    depth: u32,
}

impl Engine {
    /// Engine for `x = xser(t)`.
    pub fn new(x: &Expr, xser: Ser, d: i64, logt: &Expr) -> Engine {
        Engine { x: x.clone(), xser, d, logt: logt.clone(), memo: FxHashMap::default(), depth: 0 }
    }

    fn constant(&self, e: &Expr) -> SRes<Ser> {
        if e.is_zero() {
            return Ok(Ser::zero_exact());
        }
        if contains_bad(e) {
            return Err(SErr::Pole(e.to_string()));
        }
        Ok(Ser::constant(C::from_expr(e.clone())))
    }

    /// Series of `e` correct to `O(t^n)` (when achievable).
    pub fn ser(&mut self, e: &Expr, n: i64) -> SRes<Ser> {
        if !depends(e, &self.x) {
            return self.constant(e);
        }
        if *e == self.x {
            return Ok(self.xser.clone());
        }
        if let Some(s) = self.memo.get(e) {
            if s.exact || s.prec >= n {
                return Ok(s.trunc(n));
            }
        }
        if self.depth > MAX_DEPTH {
            return Err(SErr::Unsupported("expression too deep".into()));
        }
        self.depth += 1;
        let r = self.ser_impl(e, n);
        self.depth -= 1;
        // SymPy's Function._eval_nseries names the innermost function whose
        // argument tends to a logarithm: "Cannot expand sin(log(x)) around 0".
        let r = match r {
            Err(SErr::Pole(message)) if self.xser.val > 0 && matches!(e.kind(), Kind::Func(f, _)
                if message == format!("{} of log", f.name())) =>
                return Err(SErr::Pole(format!("Cannot expand {e} around 0"))),
            r => r?,
        };
        let better = match self.memo.get(e) {
            Some(old) => !old.exact && (r.exact || r.prec > old.prec),
            None => true,
        };
        if better {
            self.memo.insert(e.clone(), r.clone());
        }
        Ok(r.trunc(n))
    }

    /// Like [`ser`](Self::ser) but keeps increasing the precision until a
    /// non-zero term is found.
    pub fn ser_nonzero(&mut self, e: &Expr, n: i64) -> SRes<Ser> {
        let mut m = n;
        let mut step = self.d.max(2);
        loop {
            let s = self.ser(e, m)?;
            if !s.is_zero() {
                return Ok(s);
            }
            if s.exact {
                return Err(SErr::Pole(format!("{e} (identically zero)")));
            }
            if m > n + MAX_EXTRA * self.d {
                return Err(SErr::Precision);
            }
            m = s.prec.max(m) + step;
            step *= 2;
        }
    }

    fn ser_impl(&mut self, e: &Expr, n: i64) -> SRes<Ser> {
        match e.kind() {
            Kind::Add(c, terms) => {
                let mut acc = Ser::constant(C::N(c.clone()));
                for (t, k) in terms {
                    let s = self.ser(t, n)?;
                    acc = acc.add(&s.scale(&C::N(k.clone())));
                }
                Ok(acc)
            }
            Kind::Mul(c, facs) => {
                let mut cst: Vec<Expr> = vec![Expr::num(c.clone())];
                let mut fs: Vec<Expr> = Vec::with_capacity(facs.len());
                for (b, p) in facs {
                    let f = mul_factor(b, p);
                    if depends(&f, &self.x) {
                        fs.push(f);
                    } else {
                        cst.push(f);
                    }
                }
                let k = C::from_expr(mul(&cst));
                Ok(self.ser_product(&fs, n)?.scale(&k))
            }
            Kind::Pow(b, p) => self.ser_pow(b, p, n),
            Kind::Func(f, args) => self.ser_func(f, args, n),
            _ => Err(SErr::Unsupported(e.to_string())),
        }
    }

    fn ser_product(&mut self, fs: &[Expr], n: i64) -> SRes<Ser> {
        if fs.len() == 1 {
            return self.ser(&fs[0], n);
        }
        let mut sers: Vec<Ser> = Vec::with_capacity(fs.len());
        for f in fs {
            sers.push(self.ser(f, n)?);
        }
        for _ in 0..6 {
            let vals: Vec<i64> = sers.iter().map(|s| if s.exact && s.c.is_empty() { 0 } else { s.val }).collect();
            if sers.iter().any(|s| s.is_exact_zero()) {
                return Ok(Ser::zero_exact());
            }
            let total: i64 = vals.iter().sum();
            let mut changed = false;
            for i in 0..fs.len() {
                let need = n - (total - vals[i]);
                if !sers[i].exact && sers[i].prec < need {
                    sers[i] = self.ser(&fs[i], need)?;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        let mut acc = sers[0].clone();
        for s in &sers[1..] {
            acc = acc.mul(s);
        }
        // SymPy only keeps products exact when all factors are polynomials
        if acc.exact && !sers.iter().all(|s| s.val >= 0 && self.is_poly(s)) {
            acc.exact = false;
            acc.prec = n;
            acc = acc.norm();
        }
        Ok(acc)
    }

    fn is_poly(&self, s: &Ser) -> bool {
        if self.xser.val < 0 {
            return false;
        }
        s.c.iter().enumerate().all(|(i, c)| c.is_zero() || (s.val + i as i64) % self.d == 0)
    }

    fn ser_pow(&mut self, b: &Expr, p: &Expr, n: i64) -> SRes<Ser> {
        if depends(p, &self.x) {
            let e = fx::exp(&mul2(p, &fx::log(b)));
            return self.ser(&e, n);
        }
        match p.as_number() {
            Some(k) if k.is_integer() && k.is_positive() => {
                let k = k.to_i64().filter(|k| *k < 1 << 40).ok_or_else(|| SErr::Unsupported("huge power".into()))?;
                let mut s = self.ser(b, n)?;
                for _ in 0..6 {
                    let v = s.val;
                    let need = n - (k - 1) * v;
                    if s.exact || s.prec >= need {
                        break;
                    }
                    s = self.ser(b, need)?;
                }
                if s.is_zero() {
                    return Ok(if s.exact { Ser::zero_exact() } else { Ser::zero_inexact(k.saturating_mul(s.prec)) });
                }
                if s.exact && s.c.len() == 1 {
                    let v = s.val.checked_mul(k).ok_or_else(|| SErr::Unsupported("huge power".into()))?;
                    return Ok(Ser::monomial(s.c[0].pow(&C::int(k)), v).trunc(n));
                }
                if s.exact && k <= 64 {
                    // repeated products keep exact results exact
                    let v = s.val.min(0);
                    let mut result = s.clone();
                    for j in 2..=k {
                        result = result.mul(&s).trunc(n - (k - j) * v);
                    }
                    return Ok(result);
                }
                self.pow_general(&s, &C::N(Number::int(k)), Number::int(k), n)
            }
            Some(k) if k.is_rational() => {
                let mut s = self.ser_nonzero(b, n)?;
                for _ in 0..6 {
                    let v = Number::int(s.val);
                    // need = n - (k - 1) v
                    let need = Number::int(n).sub(&k.sub(&Number::ONE).mul(&v));
                    let need = need.floor().to_i64().unwrap_or(n) + 1;
                    if s.exact || s.prec >= need {
                        break;
                    }
                    s = self.ser_nonzero(b, need)?;
                }
                self.pow_general(&s, &C::N(k.clone()), k.clone(), n)
            }
            _ => {
                // symbolic (or float) exponent
                let s = self.ser_nonzero(b, n)?;
                if s.val != 0 {
                    return Err(SErr::Unsupported(format!("({b})**({p})")));
                }
                self.pow_general(&s, &C::from_expr(p.clone()), Number::ZERO, n)
            }
        }
    }

    /// `s ** r`; `rq` is `r` when it is a rational number (used for the
    /// valuation), zero otherwise (then `s.val` must be 0).
    fn pow_general(&mut self, s: &Ser, r: &C, rq: Number, n: i64) -> SRes<Ser> {
        let a0 = s.lead().ok_or(SErr::Precision)?.clone();
        let vr = if s.val == 0 {
            0
        } else {
            let v = Number::int(s.val).mul(&rq);
            if !v.is_integer() {
                let den = v.denom().to_i64().unwrap_or(2);
                return Err(SErr::NeedD(den));
            }
            v.to_i64().ok_or(SErr::Precision)?
        };
        let a0r = a0.pow(r);
        if s.exact && s.c.len() == 1 {
            return Ok(Ser::monomial(a0r, vr));
        }
        let rel = s.rel().min(n - vr);
        if rel <= 0 {
            return Ok(Ser::zero_inexact(n.min(vr.saturating_add(s.rel()))));
        }
        if matches!(r, C::N(k) if k.is_minus_one()) {
            return s.inv(n);
        }
        let len = rel as usize;
        let inv0 = a0.inv().ok_or(SErr::Precision)?;
        let u: Vec<C> = (0..len).map(|i| if i < s.c.len() { s.c[i].mul(&inv0) } else { C::zero() }).collect();
        let p = pow_kernel(&u, r, len);
        Ok(Ser { val: vr, c: p, prec: vr + rel, exact: false }.norm().scale(&a0r))
    }

    fn exp_ser(&mut self, a: &Ser, n: i64) -> SRes<Ser> {
        if !a.is_zero() && a.val < 0 {
            if sign_of(a.lead().unwrap()) == Some(-1) {
                // exponentially small: beyond all orders
                return Ok(Ser::zero_inexact(n));
            }
            return Err(SErr::Pole("exp".into()));
        }
        let c0 = a.at(0);
        let (shift, beta) = match &c0 {
            C::E(e) if depends(e, &self.logt) => {
                let (alpha, beta) =
                    linear_in(e, &self.logt).ok_or_else(|| SErr::Unsupported("exp of log".into()))?;
                // exp(alpha*log(t^D)) = t^(alpha*D)
                let s = alpha.mul(&Number::int(self.d));
                if !s.is_integer() {
                    return Err(SErr::NeedD(s.denom().to_i64().unwrap_or(2)));
                }
                (s.to_i64().ok_or(SErr::Precision)?, beta)
            }
            _ => (0, c0.to_expr()),
        };
        let ebeta = C::from_expr(fx::exp(&beta));
        let rest_zero = a.c.iter().enumerate().all(|(i, c)| a.val + i as i64 == 0 || c.is_zero());
        if rest_zero {
            if a.exact {
                return Ok(Ser::monomial(ebeta, shift));
            }
            let len = (a.prec.min(n - shift)).max(0);
            let mut c = vec![C::zero(); len as usize];
            if len > 0 {
                c[0] = ebeta;
            }
            return Ok(Ser { val: shift, c, prec: shift + len, exact: false }.norm());
        }
        let len = if a.exact { n - shift } else { a.prec.min(n - shift) };
        if len <= 0 {
            return Ok(Ser::zero_inexact(n));
        }
        let mut g = a.dense0(len as usize);
        g[0] = C::zero();
        let e = exp_kernel(&g, len as usize);
        Ok(Ser { val: shift, c: e, prec: shift + len, exact: false }.norm().scale(&ebeta))
    }

    fn log_ser(&mut self, a: &Ser, n: i64) -> SRes<Ser> {
        let a0 = a.lead().ok_or(SErr::Precision)?.clone();
        let v = a.val;
        let mut logc = fx::log(&a0.to_expr());
        if v != 0 {
            let q = Number::rational_i64(v, self.d).unwrap();
            logc = add(&[logc, mul2(&Expr::num(q), &self.logt)]);
        }
        let logc = C::from_expr(logc);
        if a.exact && a.c.len() == 1 {
            return Ok(Ser::constant(logc));
        }
        let rel = a.rel().min(n);
        if rel <= 0 {
            return Ok(Ser::zero_inexact(rel.max(0)));
        }
        let len = rel as usize;
        let inv0 = a0.inv().ok_or(SErr::Precision)?;
        let u: Vec<C> = (0..len).map(|i| if i < a.c.len() { a.c[i].mul(&inv0) } else { C::zero() }).collect();
        let mut l = log_kernel(&u, len);
        l[0] = logc;
        Ok(Ser { val: 0, c: l, prec: rel, exact: false }.norm())
    }

    fn analytic_arg(&self, a: &Ser, what: &str) -> SRes<C> {
        if !a.is_zero() && a.val < 0 {
            return Err(SErr::Pole(what.to_string()));
        }
        let c0 = a.at(0);
        if let C::E(e) = &c0 {
            if depends(e, &self.logt) {
                return Err(SErr::Pole(format!("{what} of log")));
            }
        }
        Ok(c0)
    }

    fn rest_is_zero(a: &Ser) -> bool {
        a.c.iter().enumerate().all(|(i, c)| a.val + i as i64 == 0 || c.is_zero())
    }

    fn ser_func(&mut self, f: &Func, args: &[Expr], n: i64) -> SRes<Ser> {
        if matches!(f, Func::Heaviside | Func::DiracDelta) && !args.is_empty() {
            // One-sided expansion: the argument keeps the sign of its leading
            // coefficient on the side approached, so `Heaviside` is the
            // constant 0 or 1 (its value `H0` at zero is never reached) and
            // `DiracDelta` (any order) vanishes, as with SymPy's positive
            // dummy (`Heaviside(x).series(x)` is `1`).
            let a = self.ser_nonzero(&args[0], n)?;
            return match sign_of(a.lead().unwrap()) {
                Some(s) if s != 0 && *f == Func::DiracDelta => Ok(Ser::zero_exact()),
                Some(1) => Ok(Ser::constant(C::one())),
                Some(-1) => Ok(Ser::zero_exact()),
                // a nonvanishing symbolic value: continuous there, as the
                // reference `limit(Heaviside(x - y), x, 0)` is `Heaviside(-y)`
                None if a.val == 0 && !matches!(a.lead(), Some(C::E(e)) if depends(e, &self.logt)) => {
                    let mut at = args.to_vec();
                    at[0] = a.lead().unwrap().to_expr();
                    self.constant(&func(f.clone(), at))
                }
                _ => Err(SErr::Unsupported(format!("{}({})", f.name(), args[0]))),
            };
        }
        if let Func::Special(crate::special::SpecialFn::Frac) = f {
            // frac(u) = u - floor(u)
            let fl = func(Func::Floor, vec![args[0].clone()]);
            return self.ser(&sub(&args[0], &fl), n);
        }
        if matches!(f, Func::Floor | Func::Ceiling) && args.len() == 1 {
            return self.ser_floor(*f == Func::Floor, &args[0], n);
        }
        if let Func::Special(s) = f {
            use crate::special::SpecialFn;
            let order_arg = match (s, args.len()) {
                (SpecialFn::Polygamma, 2) => Some((args[0].clone(), &args[1])),
                (SpecialFn::Digamma, 1) => Some((Expr::zero(), &args[0])),
                (SpecialFn::Trigamma, 1) => Some((Expr::one(), &args[0])),
                _ => None,
            };
            if let Some((order, z)) = order_arg {
                let order = order.as_number().filter(|k| k.is_integer()).and_then(|k| k.to_i64()).filter(|k| *k >= 0);
                if let Some(order) = order {
                    let a = self.ser(z, n.max(1))?;
                    if !a.is_zero() && a.val < 0 {
                        return self.ser_polygamma_inf(order, z, &a, n);
                    }
                }
                return Err(SErr::Unsupported(format!("{}(...)", f.name())));
            }
        }
        if let Func::Special(s) = f {
            if let Some(res) = self.ser_log_singular(*s, args, n)? {
                return Ok(res);
            }
            if let Some(res) = self.ser_special_taylor(*s, args, n)? {
                return Ok(res);
            }
        }
        if args.len() != 1 {
            return Err(SErr::Unsupported(format!("{}(...)", f.name())));
        }
        let arg = &args[0];
        match f {
            Func::Exp => {
                let a = self.ser(arg, n)?;
                self.exp_ser(&a, n)
            }
            Func::Log => {
                let mut a = self.ser_nonzero(arg, n)?;
                if !a.exact && a.prec - a.val < n {
                    a = self.ser_nonzero(arg, n + a.val)?;
                }
                self.log_ser(&a, n)
            }
            Func::Sin | Func::Cos | Func::Sinh | Func::Cosh => {
                let a = self.ser(arg, n)?;
                let c0 = self.analytic_arg(&a, f.name())?;
                let hyper = matches!(f, Func::Sinh | Func::Cosh);
                let is_sin = matches!(f, Func::Sin | Func::Sinh);
                if Self::rest_is_zero(&a) && a.exact {
                    return self.constant(&func(f.clone(), vec![c0.to_expr()]));
                }
                let len = if a.exact { n } else { a.prec.min(n) };
                if len <= 0 {
                    return Ok(Ser::zero_inexact(n));
                }
                let mut g = a.dense0(len as usize);
                g[0] = C::zero();
                let (s, c) = sincos_kernel(&g, len as usize, hyper);
                let ss = Ser { val: 0, c: s, prec: len, exact: false }.norm();
                let cs = Ser { val: 0, c, prec: len, exact: false }.norm();
                if c0.is_zero() {
                    return Ok(if is_sin { ss } else { cs });
                }
                let e0 = c0.to_expr();
                let (sf, cf) = if hyper { (Func::Sinh, Func::Cosh) } else { (Func::Sin, Func::Cos) };
                let s0 = C::from_expr(func(sf, vec![e0.clone()]));
                let k0 = C::from_expr(func(cf, vec![e0]));
                Ok(if is_sin {
                    cs.scale(&s0).add(&ss.scale(&k0))
                } else if hyper {
                    cs.scale(&k0).add(&ss.scale(&s0))
                } else {
                    cs.scale(&k0).sub(&ss.scale(&s0))
                })
            }
            Func::Tan | Func::Cot | Func::Tanh | Func::Coth => {
                let a = self.ser(arg, n)?;
                let c0 = self.analytic_arg(&a, f.name())?;
                let t0 = func(f.clone(), vec![c0.to_expr()]);
                if contains_bad(&t0) || (c0.is_zero() && matches!(f, Func::Cot | Func::Coth)) {
                    let (num, den) = match f {
                        Func::Tan => (fx::sin(arg), fx::cos(arg)),
                        Func::Cot => (fx::cos(arg), fx::sin(arg)),
                        Func::Tanh => (fx::sinh(arg), fx::cosh(arg)),
                        _ => (fx::cosh(arg), fx::sinh(arg)),
                    };
                    let e = mul2(&num, &pow(&den, &Expr::minus_one()));
                    return self.ser(&e, n);
                }
                if Self::rest_is_zero(&a) && a.exact {
                    return self.constant(&t0);
                }
                let (alpha, sigma) = match f {
                    Func::Tan => (1, 1),
                    Func::Cot => (-1, 1),
                    Func::Tanh => (1, -1),
                    _ => (1, -1),
                };
                let len = if a.exact { n } else { a.prec.min(n) };
                if len <= 0 {
                    return Ok(Ser::zero_inexact(n));
                }
                let g = a.dense0(len as usize);
                let t = tan_kernel(&g, len as usize, &C::from_expr(t0), alpha, sigma);
                Ok(Ser { val: 0, c: t, prec: len, exact: false }.norm())
            }
            Func::Sec => self.ser(&pow(&fx::cos(arg), &Expr::minus_one()), n),
            Func::Csc => self.ser(&pow(&fx::sin(arg), &Expr::minus_one()), n),
            Func::Asin | Func::Acos | Func::Atan | Func::Acot | Func::Asinh | Func::Acosh | Func::Atanh => {
                self.ser_inverse(f, arg, n)
            }
            Func::Erf | Func::Erfc => {
                let a = self.ser(arg, n)?;
                if !a.is_zero() && a.val < 0 {
                    // erf(z) = ±1 - (exponentially small) for z -> ±oo (real)
                    let s = sign_of(a.lead().unwrap()).ok_or_else(|| SErr::Pole(format!("erf({arg})")))?;
                    let v = match (f, s) {
                        (Func::Erf, s) => s as i64,
                        (_, 1) => 0,
                        _ => 2,
                    };
                    let mut r = Ser::zero_inexact(n.max(1));
                    if v != 0 {
                        r = Ser::constant(C::int(v)).add(&r);
                    }
                    return Ok(r);
                }
                let mut m = n;
                for _ in 0..5 {
                    let a = self.ser(arg, m)?;
                    let c0 = self.analytic_arg(&a, f.name())?;
                    let a2 = a.mul(&a).neg();
                    let ex = self.exp_ser(&a2, m)?;
                    let integ = ex.mul(&a.deriv()).integ()?;
                    let k = C::from_expr(mul2(&Expr::int(2), &pow(&Expr::pi(), &Expr::rational(-1, 2))));
                    let mut r = integ.scale(&k).add(&Ser::constant(C::from_expr(fx::erf(&c0.to_expr()))));
                    if *f == Func::Erfc {
                        r = Ser::constant(C::one()).sub(&r);
                    }
                    if r.exact || r.prec >= n {
                        return Ok(r);
                    }
                    m += n - r.prec;
                }
                Err(SErr::Precision)
            }
            Func::Gamma | Func::LogGamma => self.ser_gamma(f, arg, n),
            Func::Abs => {
                let a = self.ser_nonzero(arg, n)?;
                match sign_of(a.lead().unwrap()) {
                    Some(1) => Ok(a),
                    Some(-1) => Ok(a.neg()),
                    _ => Err(SErr::Unsupported(format!("Abs({arg})"))),
                }
            }
            Func::Sign => {
                let a = self.ser_nonzero(arg, n)?;
                match sign_of(a.lead().unwrap()) {
                    Some(s) => Ok(Ser::constant(C::int(s as i64))),
                    _ => Err(SErr::Unsupported(format!("sign({arg})"))),
                }
            }
            _ => Err(SErr::Unsupported(format!("{}({arg})", f.name()))),
        }
    }

    /// One-sided `floor`/`ceiling`: a constant near the point. A non-integer
    /// value `c0` of the argument gives `floor(c0)`/`ceiling(c0)`; at an
    /// integer the side is the sign of the first non-constant term, as in
    /// `RoundFunction._eval_nseries` (`floor(x).series(x)` is `0`).
    fn ser_floor(&mut self, is_floor: bool, arg: &Expr, n: i64) -> SRes<Ser> {
        let unsupported = || SErr::Unsupported(format!("{}({arg})", if is_floor { "floor" } else { "ceiling" }));
        let a = self.ser(arg, n.max(1))?;
        if !a.is_zero() && a.val < 0 {
            return Err(SErr::Pole(format!("floor({arg})")));
        }
        let c0 = a.at(0);
        if let C::E(e) = &c0 {
            if depends(e, &self.logt) {
                return Err(unsupported());
            }
        }
        let c0e = c0.to_expr();
        let integer = c0e.as_number().map(|k| k.is_integer());
        if integer != Some(true) {
            let fl = func(Func::Floor, vec![c0e.clone()]);
            let ce = func(Func::Ceiling, vec![c0e.clone()]);
            let gap = sub(&ce, &fl);
            if fl.as_number().is_some() && gap.as_number().is_some_and(|g| g.is_one()) {
                return Ok(Ser::constant(C::from_expr(if is_floor { fl } else { ce })));
            }
            return Err(unsupported());
        }
        let rest = sub(arg, &c0e);
        let s = match self.ser_nonzero(&rest, n.max(1)) {
            Ok(r) => sign_of(r.lead().unwrap()).ok_or_else(unsupported)?,
            Err(SErr::Pole(msg)) if msg.ends_with("(identically zero)") => 0,
            Err(err) => return Err(err),
        };
        let k = match (is_floor, s) {
            (_, 0) => c0e,
            (true, 1) | (false, -1) => c0e,
            (true, _) => sub(&c0e, &Expr::one()),
            (false, _) => add(&[c0e, Expr::one()]),
        };
        Ok(Ser::constant(C::from_expr(k)))
    }

    /// Inverse trigonometric/hyperbolic functions through
    /// `F(a) = F(a0) + ∫ F'(a) a' dt`.
    fn ser_inverse(&mut self, f: &Func, arg: &Expr, n: i64) -> SRes<Ser> {
        let a = self.ser(arg, n)?;
        if !a.is_zero() && a.val < 0 {
            // atan/acot at infinity
            return match f {
                Func::Atan | Func::Acot => {
                    let s = sign_of(a.lead().unwrap()).ok_or_else(|| SErr::Pole(format!("atan({arg})")))?;
                    let inv = pow(arg, &Expr::minus_one());
                    let e = if *f == Func::Atan {
                        add(&[mul(&[Expr::rational(s as i64, 2), Expr::pi()]), neg(&fx::atan(&inv))])
                    } else {
                        fx::atan(&inv)
                    };
                    self.ser(&e, n)
                }
                _ => Err(SErr::Pole(format!("{}({arg})", f.name()))),
            };
        }
        // (sigma, rho, sign): F' = sign * (1 + sigma a^2)^rho
        let (sigma, rho, sgn): (i64, Number, i64) = match f {
            Func::Atan => (1, Number::MINUS_ONE, 1),
            Func::Acot => (1, Number::MINUS_ONE, -1),
            Func::Asin => (-1, Number::rational_i64(-1, 2).unwrap(), 1),
            Func::Acos => (-1, Number::rational_i64(-1, 2).unwrap(), -1),
            Func::Asinh => (1, Number::rational_i64(-1, 2).unwrap(), 1),
            Func::Acosh => (-1, Number::rational_i64(-1, 2).unwrap(), 1),
            _ => (-1, Number::MINUS_ONE, 1), // atanh
        };
        let mut m = n;
        let mut a = a;
        for _ in 0..6 {
            let c0 = self.analytic_arg(&a, f.name())?;
            let f0 = func(f.clone(), vec![c0.to_expr()]);
            if contains_bad(&f0) {
                return Err(SErr::Pole(format!("{}({arg})", f.name())));
            }
            if a.exact && Self::rest_is_zero(&a) {
                return self.constant(&f0);
            }
            let a2 = a.mul(&a);
            let mut g = Ser::constant(C::one()).add(&a2.scale(&C::int(sigma)));
            if *f == Func::Acosh {
                g = g.neg(); // a^2 - 1
            }
            if g.is_zero() {
                m += self.d.max(2) * 2;
                a = self.ser(arg, m)?;
                continue;
            }
            let h = self.pow_general(&g, &C::N(rho.clone()), rho.clone(), m)?;
            let integ = h.mul(&a.deriv()).integ()?;
            let r = integ.scale(&C::int(sgn)).add(&Ser::constant(C::from_expr(f0)));
            if r.prec >= n {
                return Ok(r);
            }
            m += (n - r.prec).max(1);
            a = self.ser(arg, m)?;
        }
        Err(SErr::Precision)
    }

    /// `gamma`/`loggamma` at integer points through
    /// `loggamma(1 + r) = -EulerGamma r + Σ (-1)^k zeta(k) r^k / k`.
    fn ser_gamma(&mut self, f: &Func, arg: &Expr, n: i64) -> SRes<Ser> {
        let a = self.ser(arg, n.max(1))?;
        if !a.is_zero() && a.val < 0 {
            return self.ser_stirling(f, arg, &a, n);
        }
        let mut m = n;
        for _ in 0..5 {
            let a = self.ser(arg, m)?;
            let c0 = self.analytic_arg(&a, f.name())?;
            let r = a.sub(&Ser::constant(c0.clone()));
            if r.is_zero() {
                if r.exact {
                    return self.constant(&func(f.clone(), vec![c0.to_expr()]));
                }
                m += n.max(1);
                continue;
            }
            let k0 = match &c0 {
                C::N(Number::Int(k)) => *k,
                _ => {
                    // regular point: log Gamma(c0 + r) = log Gamma(c0)
                    // + sum_{k>=1} polygamma(k - 1, c0) r^k / k!
                    let res = self.ser_gamma_regular(f, &c0, &r, m)?;
                    if res.prec >= n {
                        return Ok(res);
                    }
                    m += (n - res.prec).max(1);
                    continue;
                }
            };
            if *f == Func::LogGamma && k0 == 0 {
                // loggamma._eval_nseries: at a zero of the argument the
                // reference expands the intractable form log(gamma(z)).
                return self.ser(&fx::log(&fx::gamma(arg)), n);
            }
            if *f == Func::LogGamma && k0 < 1 {
                return Err(SErr::Pole(format!("loggamma({arg})")));
            }
            let vr = r.val.max(1);
            let prec = if r.exact { m + vr * (k0.unsigned_abs() as i64 + 2) } else { r.prec };
            let kmax = ((prec + vr - 1) / vr).max(1) as usize;
            // L(r) = Σ_{k>=1} l_k r^k (Horner)
            let mut lk: Vec<C> = vec![C::zero(); kmax + 1];
            lk[1] = C::from_expr(neg(&Expr::constant(Constant::EulerGamma)));
            for (k, l) in lk.iter_mut().enumerate().skip(2) {
                let z = C::from_expr(zeta_int(k as i64));
                let z = z.div_int(k as i64);
                *l = if k % 2 == 0 { z } else { z.neg() };
            }
            let lser = r.compose(&lk, prec);
            let mut res = if *f == Func::LogGamma {
                lser
            } else {
                self.exp_ser(&lser, prec)?
            };
            // shift from 1 + r to k0 + r
            if k0 >= 1 {
                for j in 1..k0 {
                    let fac = r.add(&Ser::constant(C::int(j)));
                    res = if *f == Func::LogGamma {
                        res.add(&self.log_ser(&fac, prec)?)
                    } else {
                        res.mul(&fac)
                    };
                }
            } else {
                for j in 0..=(-k0) {
                    let fac = r.sub(&Ser::constant(C::int(j)));
                    res = res.mul(&fac.inv(prec)?);
                }
            }
            if res.prec >= n {
                return Ok(res);
            }
            m += (n - res.prec).max(1);
        }
        Err(SErr::Precision)
    }

    /// `loggamma(z)` for `z -> +oo`: the Stirling series of the reference
    /// `loggamma._eval_aseries`,
    /// `(z - 1/2)*log(z) - z + log(2*pi)/2 + Σ B_2k/(2k(2k - 1) z**(2k - 1))`,
    /// with the remainder `O(z**-(2K + 1))` kept as the precision. `gamma`
    /// has no asymptotic series there (the reference raises `PoleError`).
    fn ser_stirling(&mut self, f: &Func, arg: &Expr, a: &Ser, n: i64) -> SRes<Ser> {
        if *f != Func::LogGamma || sign_of(a.lead().unwrap()) != Some(1) {
            return Err(SErr::Pole(format!("{}({arg})", f.name())));
        }
        let v = -a.val;
        let k_terms = if n <= v { 0 } else { (n - v + 2 * v - 1) / (2 * v) };
        let half = Expr::rational(1, 2);
        let mut terms = vec![
            mul2(&sub(arg, &half), &fx::log(arg)),
            neg(arg),
            mul2(&half, &add(&[fx::log(&Expr::int(2)), fx::log(&Expr::pi())])),
        ];
        let b = bernoulli(2 * k_terms as usize);
        for k in 1..=k_terms {
            let c = b[2 * k as usize].div(&Number::int(2 * k * (2 * k - 1))).unwrap();
            terms.push(mul2(&Expr::num(c), &pow(arg, &Expr::int(1 - 2 * k))));
        }
        let r = self.ser(&add(&terms), n)?;
        Ok(r.add(&Ser::zero_inexact((2 * k_terms + 1) * v)))
    }

    /// `polygamma(N, z)` for `z -> +oo` (Abramowitz & Stegun 6.3.18 and
    /// 6.4.10, as the reference `polygamma._eval_aseries`), with the
    /// remainder of the truncated series kept as the precision.
    fn ser_polygamma_inf(&mut self, order: i64, z: &Expr, a: &Ser, n: i64) -> SRes<Ser> {
        if sign_of(a.lead().unwrap()) != Some(1) {
            return Err(SErr::Pole(format!("polygamma({order}, {z})")));
        }
        let v = -a.val;
        let lead = if order == 0 { 2 } else { order + 2 };
        let k_terms = if n <= lead * v { 0 } else { (n - lead * v + 2 * v - 1) / (2 * v) };
        let b = bernoulli(2 * k_terms as usize);
        let zpow = |k: i64| pow(z, &Expr::int(-k));
        let mut terms = vec![];
        if order == 0 {
            terms.push(fx::log(z));
            terms.push(mul2(&Expr::rational(-1, 2), &zpow(1)));
            for k in 1..=k_terms {
                let c = b[2 * k as usize].div(&Number::int(-2 * k)).unwrap();
                terms.push(mul2(&Expr::num(c), &zpow(2 * k)));
            }
        } else {
            let fact = |m: i64| (1..=m).fold(Number::ONE, |acc, j| acc.mul(&Number::int(j)));
            let sgn = if order % 2 == 1 { Number::ONE } else { Number::MINUS_ONE };
            terms.push(mul2(&Expr::num(sgn.mul(&fact(order - 1))), &zpow(order)));
            terms.push(mul2(&Expr::num(sgn.mul(&fact(order)).div(&Number::int(2)).unwrap()), &zpow(order + 1)));
            for k in 1..=k_terms {
                let c = b[2 * k as usize].mul(&fact(2 * k + order - 1)).div(&fact(2 * k)).unwrap().mul(&sgn);
                terms.push(mul2(&Expr::num(c), &zpow(2 * k + order)));
            }
        }
        let r = self.ser(&add(&terms), n)?;
        Ok(r.add(&Ser::zero_inexact((2 * k_terms + lead) * v)))
    }

    fn ser_gamma_regular(&mut self, f: &Func, c0: &C, r: &Ser, m: i64) -> SRes<Ser> {
        let c0e = c0.to_expr();
        let g0 = func(f.clone(), vec![c0e.clone()]);
        if contains_bad(&g0) {
            return Err(SErr::Pole(g0.to_string()));
        }
        let vr = r.val.max(1);
        let prec = if r.exact { m + vr } else { r.prec };
        let kmax = ((prec + vr - 1) / vr).max(1) as usize;
        let mut lk: Vec<C> = vec![C::zero(); kmax + 1];
        let mut fact = Number::ONE;
        for (k, l) in lk.iter_mut().enumerate().skip(1) {
            fact = fact.mul(&Number::int(k as i64));
            let pg = crate::special::polygamma(&Expr::int(k as i64 - 1), &c0e);
            if contains_bad(&pg) {
                return Err(SErr::Pole(pg.to_string()));
            }
            *l = C::from_expr(crate::ops::div(&pg, &Expr::num(fact.clone())));
        }
        let lser = r.compose(&lk, prec);
        if *f == Func::LogGamma {
            Ok(Ser::constant(C::from_expr(g0)).add(&lser))
        } else {
            Ok(self.exp_ser(&lser, prec)?.scale(&C::from_expr(g0)))
        }
    }

    /// Expansions of the exponential and trigonometric integrals whose
    /// argument tends to zero, where they have a logarithmic singularity
    /// (SymPy 1.14's `_eval_nseries` through the Si/Ei rewrites and the
    /// `_eval_as_leading_term` of Ei, Ci and Chi):
    /// `Ei(u) = EulerGamma + log(u) + Σ u**k/(k*k!)`, minus `I*pi` when
    /// the leading coefficient of `u` has a negative real part, as SymPy's
    /// `Ei._eval_as_leading_term`;
    /// `Ci(u) = EulerGamma + log(u) + Σ (-1)**k*u**(2*k)/(2*k*(2*k)!)`,
    /// `Chi(u) = EulerGamma + log(u) + Σ u**(2*k)/(2*k*(2*k)!)`,
    /// `expint(1, u) = -EulerGamma - log(u) - Σ (-u)**k/(k*k!)`, and
    /// `li(u) = Ei(log(u))` when `u` tends to 1 (SymPy's `li._eval_nseries`).
    fn ser_log_singular(&mut self, s: crate::special::SpecialFn, args: &[Expr], n: i64) -> SRes<Option<Ser>> {
        use crate::special::SpecialFn::*;
        let u = match (s, args.len()) {
            (Ei | Ci | Chi | Li, 1) => &args[0],
            (Expint, 2) if args[0].is_one() => &args[1],
            _ => return Ok(None),
        };
        // The order of the argument's leading term may exceed `n` (`Ei(2*x)`
        // at O(x), `Ei(x**2)`): expand until it is nonzero.
        let Ok(a) = self.ser_nonzero(u, n.max(1)) else { return Ok(None) };
        if a.is_zero() || a.val < 0 {
            return Ok(None);
        }
        let c0 = a.lead().cloned().ok_or(SErr::Precision)?;
        if let C::E(e) = &c0 {
            if depends(e, &self.logt) {
                return Ok(None);
            }
        }
        if s == Li {
            if a.val == 0 && c0.to_expr().is_one() {
                let ei = crate::special::sf(Ei, vec![fx::log(u)]);
                return self.ser(&ei, n).map(Some);
            }
            return Ok(None);
        }
        if a.val == 0 {
            return Ok(None);
        }
        // Terms u**k contribute from t**(k*val/d): enough of them for O(t**n).
        let kmax = (n.max(1) / a.val + 2).max(2);
        let mut terms = vec![Expr::constant(Constant::EulerGamma), fx::log(u)];
        let mut fact = Expr::one();
        for k in 1..=(if matches!(s, Ci | Chi) { 2 * kmax } else { kmax }) {
            fact = mul2(&fact, &Expr::int(k));
            let ke = Expr::int(k);
            let term = match s {
                Ei => Some(crate::ops::div(&pow(u, &ke), &mul2(&ke, &fact))),
                Expint => Some(crate::ops::div(&pow(&neg(u), &ke), &mul2(&ke, &fact))),
                Ci | Chi if k % 2 == 0 => {
                    let t = crate::ops::div(&pow(u, &ke), &mul2(&ke, &fact));
                    Some(if s == Ci && k % 4 == 2 { neg(&t) } else { t })
                }
                _ => None,
            };
            if let Some(term) = term {
                terms.push(term);
            }
        }
        if s == Ei {
            let real = func(Func::Re, vec![c0.to_expr()]);
            if crate::assumptions::is_negative(&real) == Some(true) {
                terms.push(neg(&mul2(&Expr::i(), &Expr::pi())));
            }
        }
        let value = add(&terms);
        let value = if s == Expint { neg(&value) } else { value };
        self.ser(&value, n).map(Some)
    }

    /// Taylor expansions of the Airy, elliptic and `sinc` functions at a
    /// regular point, following SymPy 1.14: `sinc(u)` is `sin(u)/u`, Airy
    /// `Ai` and `Bi` use their `taylor_term` at 0, `elliptic_k(m)` and
    /// `elliptic_e(m)` their hypergeometric series at 0, and the other
    /// cases the derivatives of `Function._eval_nseries` (each term
    /// expanded). `None` leaves the other functions unsupported.
    fn ser_special_taylor(&mut self, s: crate::special::SpecialFn, args: &[Expr], n: i64) -> SRes<Option<Ser>> {
        use crate::special::SpecialFn::*;
        if !matches!(s, AiryAi | AiryBi | AiryAiPrime | AiryBiPrime | EllipticK | EllipticE | EllipticF | EllipticPi | Sinc)
            || !s.nargs().contains(&args.len())
        {
            return Ok(None);
        }
        if s == Sinc {
            let u = &args[0];
            return self.ser(&crate::ops::div(&fx::sin(u), u), n).map(Some);
        }
        let deps: Vec<usize> = (0..args.len()).filter(|&i| depends(&args[i], &self.x)).collect();
        if deps.len() != 1 {
            return Ok(None);
        }
        let i = deps[0];
        let one_arg_hyper = args.len() == 1 && matches!(s, EllipticK | EllipticE);
        let name = crate::special::sf(s, args.to_vec()).to_string();
        let mut m = n;
        for _ in 0..5 {
            let a = self.ser(&args[i], m.max(1))?;
            let c0 = self.analytic_arg(&a, &name)?;
            let r = a.sub(&Ser::constant(c0.clone()));
            let c0e = c0.to_expr();
            let at = |v: &Expr| {
                let mut w = args.to_vec();
                w[i] = v.clone();
                crate::special::sf(s, w)
            };
            if r.is_zero() {
                if r.exact {
                    let v = at(&c0e);
                    if contains_bad(&v) {
                        return Err(SErr::Pole(name));
                    }
                    return self.constant(&v).map(Some);
                }
                m += n.max(1);
                continue;
            }
            let vr = r.val.max(1);
            let prec = if r.exact { m + vr } else { r.prec };
            let kmax = ((prec + vr - 1) / vr).max(1) as usize;
            let at_zero = c0e.is_zero();
            if one_arg_hyper && !at_zero {
                return Ok(None);
            }
            let mut coeffs: Vec<C> = Vec::with_capacity(kmax + 1);
            if at_zero && (one_arg_hyper || matches!(s, AiryAi | AiryBi)) {
                let mut fact = Expr::one();
                for k in 0..=kmax {
                    if k > 0 {
                        fact = mul2(&fact, &Expr::int(k as i64));
                    }
                    let ke = Expr::int(k as i64);
                    let c = match s {
                        AiryAi | AiryBi => {
                            let g = fx::gamma(&crate::ops::div(&add(&[ke.clone(), Expr::one()]), &Expr::int(3)));
                            let sn = fx::sin(&mul(&[Expr::rational(2, 3), Expr::pi(), add(&[ke.clone(), Expr::one()])]));
                            let (sn, root) = if s == AiryAi {
                                (sn, pow(&Expr::int(3), &Expr::rational(2, 3)))
                            } else {
                                (fx::abs(&sn), pow(&Expr::int(3), &Expr::rational(1, 6)))
                            };
                            mul(&[
                                pow(&mul2(&root, &Expr::pi()), &Expr::minus_one()),
                                g,
                                sn,
                                pow(&fact, &Expr::minus_one()),
                                pow(&pow(&Expr::int(3), &Expr::rational(1, 3)), &ke),
                            ])
                        }
                        _ => {
                            let a0 = if s == EllipticK { Expr::rational(1, 2) } else { Expr::rational(-1, 2) };
                            let rfk = |b: &Expr| mul(&(0..k).map(|j| add(&[b.clone(), Expr::int(j as i64)])).collect::<Vec<_>>());
                            let num = mul2(&rfk(&a0), &rfk(&Expr::rational(1, 2)));
                            let den = mul2(&rfk(&Expr::one()), &fact);
                            mul(&[Expr::pi(), Expr::rational(1, 2), crate::ops::div(&num, &den)])
                        }
                    };
                    if contains_bad(&c) {
                        return Err(SErr::Unsupported(name));
                    }
                    coeffs.push(C::from_expr(c));
                }
            } else {
                let u = Expr::dummy("u", Default::default());
                let mut g = at(&u);
                let mut fact = Expr::one();
                for k in 0..=kmax {
                    if k > 0 {
                        fact = mul2(&fact, &Expr::int(k as i64));
                        g = crate::diff::diff(&g, &u);
                    }
                    let v = crate::subs::subs(&g, &[(u.clone(), c0e.clone())]);
                    if contains_bad(&v) || depends(&v, &u) {
                        return Err(SErr::Unsupported(name));
                    }
                    coeffs.push(C::from_expr(crate::expand::expand(&crate::ops::div(&v, &fact))));
                }
            }
            let res = r.compose(&coeffs, prec);
            if res.prec >= n {
                return Ok(Some(res));
            }
            m += (n - res.prec).max(1);
        }
        Err(SErr::Unsupported(name))
    }

    /// Leading term `(coefficient, exponent in t-units)`.
    pub fn leading(&mut self, e: &Expr) -> SRes<(C, i64)> {
        let s = self.ser_nonzero(e, self.d)?;
        Ok((s.lead().unwrap().clone(), s.val))
    }

    pub fn clear(&mut self) {
        self.memo.clear();
    }
}

/// Series of the variable itself for the expansion point.
pub fn var_series(x0: &Expr, minus: bool, d: i64) -> Ser {
    let sgn = if minus { -1 } else { 1 };
    if x0.is_constant(Constant::Infinity) {
        return Ser::monomial(C::int(sgn), -d);
    }
    if x0.is_constant(Constant::NegInfinity) {
        return Ser::monomial(C::int(-sgn), -d);
    }
    let t = Ser::monomial(C::int(sgn), d);
    if x0.is_zero() {
        t
    } else {
        t.add(&Ser::constant(C::from_expr(x0.clone())))
    }
}

/// Bernoulli numbers `B_0 .. B_n` (with `B_1 = -1/2`).
fn bernoulli(n: usize) -> Vec<Number> {
    let mut b: Vec<Number> = vec![Number::ONE];
    for m in 1..=n {
        // B_m = -1/(m+1) Σ_{k<m} C(m+1, k) B_k
        let mut acc = Number::ZERO;
        let mut binom = Number::ONE; // C(m+1, 0)
        for (k, bk) in b.iter().enumerate() {
            acc = acc.add(&binom.mul(bk));
            // C(m+1, k+1) = C(m+1, k) * (m+1-k)/(k+1)
            binom = binom
                .mul(&Number::int((m + 1 - k) as i64))
                .div(&Number::int((k + 1) as i64))
                .unwrap();
        }
        b.push(acc.div(&Number::int(-(m as i64) - 1)).unwrap());
    }
    b
}

/// `zeta(k)` for an integer `k >= 2`: closed form for even `k`.
pub fn zeta_int(k: i64) -> Expr {
    if k % 2 == 0 {
        // zeta(2n) = (-1)^(n+1) B_2n (2 pi)^(2n) / (2 (2n)!)
        let b = bernoulli(k as usize);
        let mut fact = Number::ONE;
        for i in 2..=k {
            fact = fact.mul(&Number::int(i));
        }
        let two_k = Number::int(2).pow_i64(k).unwrap();
        let mut c = b[k as usize].mul(&two_k).div(&fact.mul(&Number::int(2))).unwrap();
        if (k / 2) % 2 == 0 {
            c = c.neg();
        }
        mul2(&Expr::num(c), &pow(&Expr::pi(), &Expr::int(k)))
    } else {
        crate::special::sf(crate::special::SpecialFn::Zeta, vec![Expr::int(k)])
    }
}
