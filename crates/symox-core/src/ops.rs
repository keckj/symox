//! Canonical constructors: `add`, `mul`, `pow` and friends.
//!
//! These implement SymPy's automatic simplification rules (flattening,
//! collecting like terms and powers, numeric folding, exact radicals,
//! distribution of a numeric coefficient over a sum, ...). Every expression
//! created through them is in canonical form, which is what makes
//! structural equality meaningful.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, ToPrimitive};

use crate::expr::{Constant, Expr, Func, Kind};
use crate::number::{exact_root, trial_factor, Number};

// ----------------------------------------------------------------------
// helpers on the internal representation
// ----------------------------------------------------------------------

/// `k * t` where `t` is a term of an `Add` (no numeric coefficient).
pub fn add_term(t: &Expr, k: &Number) -> Expr {
    if k.is_one() {
        return t.clone();
    }
    match t.kind() {
        Kind::Mul(c, facs) if c.is_one() => Expr::from_kind(Kind::Mul(k.clone(), facs.clone())),
        Kind::Pow(b, e) => Expr::from_kind(Kind::Mul(k.clone(), vec![(b.clone(), e.clone())])),
        _ => Expr::from_kind(Kind::Mul(k.clone(), vec![(t.clone(), Expr::one())])),
    }
}

/// `b ** e` where `(b, e)` is a factor of a `Mul`.
pub fn mul_factor(b: &Expr, e: &Expr) -> Expr {
    if e.is_one() {
        b.clone()
    } else {
        Expr::from_kind(Kind::Pow(b.clone(), e.clone()))
    }
}

/// Split an expression into `(numeric coefficient, rest)`.
pub fn as_coeff_mul(e: &Expr) -> (Number, Expr) {
    match e.kind() {
        Kind::Num(n) => (n.clone(), Expr::one()),
        Kind::Mul(c, facs) if !c.is_one() => (c.clone(), mul_from_factors(facs.clone())),
        _ => (Number::ONE, e.clone()),
    }
}

/// Split an expression into `(numeric constant term, rest)`.
pub fn as_coeff_add(e: &Expr) -> (Number, Expr) {
    match e.kind() {
        Kind::Num(n) => (n.clone(), Expr::zero()),
        Kind::Add(c, terms) if !c.is_zero() => {
            (c.clone(), Expr::from_kind(Kind::Add(Number::ZERO, terms.clone())).normalize_add())
        }
        _ => (Number::ZERO, e.clone()),
    }
}

/// `(base, exponent)` view of an expression.
pub fn as_base_exp(e: &Expr) -> (Expr, Expr) {
    match e.kind() {
        Kind::Pow(b, x) => (b.clone(), x.clone()),
        _ => (e.clone(), Expr::one()),
    }
}

/// Build a Mul with coefficient 1 from already canonical sorted factors.
fn mul_from_factors(facs: Vec<(Expr, Expr)>) -> Expr {
    if facs.len() == 1 {
        let (b, e) = &facs[0];
        return mul_factor(b, e);
    }
    Expr::from_kind(Kind::Mul(Number::ONE, facs))
}

impl Expr {
    /// Fix an `Add` node which may have become degenerate (0 or 1 term).
    fn normalize_add(self) -> Expr {
        match self.kind() {
            Kind::Add(c, terms) if terms.len() == 1 && c.is_zero() => add_term(&terms[0].0, &terms[0].1),
            Kind::Add(c, terms) if terms.is_empty() => Expr::num(c.clone()),
            _ => self,
        }
    }
}

// ----------------------------------------------------------------------
// Add
// ----------------------------------------------------------------------

#[derive(Default)]
struct InfFlags {
    pos: bool,
    neg: bool,
    complex: bool,
}

/// Canonical sum of `args`.
pub fn add(args: &[Expr]) -> Expr {
    match args.len() {
        0 => return Expr::zero(),
        1 => return args[0].clone(),
        _ => {}
    }
    let mut coeff = Number::ZERO;
    let mut terms: Vec<(Expr, Number)> = Vec::with_capacity(args.len() + 4);
    let mut inf = InfFlags::default();
    for a in args {
        match a.kind() {
            Kind::Num(n) => coeff = coeff.add(n),
            Kind::Add(c, ts) => {
                coeff = coeff.add(c);
                terms.extend(ts.iter().cloned());
            }
            Kind::Mul(c, facs) if !c.is_one() => {
                // a term `c * oo` never exists: mul normalises it
                terms.push((mul_from_factors(facs.clone()), c.clone()));
            }
            Kind::Const(Constant::NaN) => return Expr::nan(),
            Kind::Const(Constant::Infinity) => inf.pos = true,
            Kind::Const(Constant::NegInfinity) => inf.neg = true,
            Kind::Const(Constant::ComplexInfinity) => inf.complex = true,
            _ => terms.push((a.clone(), Number::ONE)),
        }
    }
    finish_add(coeff, terms, inf)
}

/// Fast binary addition.
pub fn add2(a: &Expr, b: &Expr) -> Expr {
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }
    if let (Kind::Num(x), Kind::Num(y)) = (a.kind(), b.kind()) {
        return Expr::num(x.add(y));
    }
    add(&[a.clone(), b.clone()])
}

fn finish_add(coeff: Number, mut terms: Vec<(Expr, Number)>, inf: InfFlags) -> Expr {
    if coeff.is_nan() {
        return Expr::nan();
    }
    // sort + merge equal terms
    terms.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let mut out: Vec<(Expr, Number)> = Vec::with_capacity(terms.len());
    for (t, k) in terms {
        match out.last_mut() {
            Some((lt, lk)) if *lt == t => *lk = lk.add(&k),
            _ => out.push((t, k)),
        }
    }
    out.retain(|(_, k)| !k.is_zero());
    let mut coeff = coeff;
    if inf.pos || inf.neg || inf.complex {
        if (inf.pos && inf.neg) || (inf.complex && (inf.pos || inf.neg)) {
            return Expr::nan();
        }
        coeff = Number::ZERO;
        let c = if inf.complex {
            Constant::ComplexInfinity
        } else if inf.pos {
            Constant::Infinity
        } else {
            Constant::NegInfinity
        };
        if out.is_empty() {
            return Expr::constant(c);
        }
        out.push((Expr::constant(c), Number::ONE));
        out.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    }
    if out.is_empty() {
        return Expr::num(coeff);
    }
    if out.len() == 1 && coeff.is_zero() {
        let (t, k) = &out[0];
        return add_term(t, k);
    }
    Expr::from_kind(Kind::Add(coeff, out))
}

/// Build a sum directly from a map `term -> coefficient` (terms must be
/// canonical "coefficient-free" terms, distinct). Used by fast paths.
pub fn add_from_terms(coeff: Number, terms: Vec<(Expr, Number)>) -> Expr {
    finish_add(coeff, terms, InfFlags::default())
}

// ----------------------------------------------------------------------
// Mul
// ----------------------------------------------------------------------

/// Canonical product of `args`.
pub fn mul(args: &[Expr]) -> Expr {
    match args.len() {
        0 => return Expr::one(),
        1 => return args[0].clone(),
        _ => {}
    }
    let mut coeff = Number::ONE;
    let mut facs: Vec<(Expr, Expr)> = Vec::with_capacity(args.len() + 4);
    let mut infs: Vec<Constant> = Vec::new();
    for a in args {
        match a.kind() {
            Kind::Num(n) => coeff = coeff.mul(n),
            Kind::Mul(c, fs) => {
                coeff = coeff.mul(c);
                facs.extend(fs.iter().cloned());
            }
            Kind::Pow(b, e) => facs.push((b.clone(), e.clone())),
            Kind::Const(Constant::NaN) => return Expr::nan(),
            Kind::Const(c) if c.is_infinite() => infs.push(*c),
            _ => facs.push((a.clone(), Expr::one())),
        }
    }
    if !infs.is_empty() {
        return mul_infinite(coeff, facs, &infs);
    }
    if coeff.is_nan() {
        return Expr::nan();
    }
    if coeff.is_zero() {
        // 0 * x = 0 (also for floats, like SymPy)
        return Expr::zero();
    }
    finish_mul(coeff, facs, 0)
}

/// Fast binary multiplication.
pub fn mul2(a: &Expr, b: &Expr) -> Expr {
    if a.is_one() {
        return b.clone();
    }
    if b.is_one() {
        return a.clone();
    }
    match (a.kind(), b.kind()) {
        (Kind::Num(x), Kind::Num(y)) => return Expr::num(x.mul(y)),
        (Kind::Num(x), _) if x.is_zero() && !matches!(b.kind(), Kind::Const(_)) => {
            return Expr::zero()
        }
        (_, Kind::Num(y)) if y.is_zero() && !matches!(a.kind(), Kind::Const(_)) => {
            return Expr::zero()
        }
        _ => {}
    }
    mul(&[a.clone(), b.clone()])
}

fn mul_infinite(coeff: Number, facs: Vec<(Expr, Expr)>, infs: &[Constant]) -> Expr {
    if coeff.is_zero() {
        return Expr::nan();
    }
    let mut sign_neg = coeff.is_negative();
    let mut complex = false;
    for c in infs {
        match c {
            Constant::NegInfinity => sign_neg = !sign_neg,
            Constant::ComplexInfinity => complex = true,
            _ => {}
        }
    }
    let inf = if complex {
        Constant::ComplexInfinity
    } else {
        Constant::Infinity
    };
    let rest = if facs.is_empty() {
        Expr::one()
    } else {
        finish_mul(Number::ONE, facs, 0)
    };
    if rest.is_one() {
        return Expr::constant(match (inf, sign_neg) {
            (Constant::Infinity, true) => Constant::NegInfinity,
            (c, _) => c,
        });
    }
    if rest.is_zero() {
        return Expr::nan();
    }
    let c = if sign_neg && !complex {
        Number::MINUS_ONE
    } else {
        Number::ONE
    };
    let (rc, rest) = as_coeff_mul(&rest);
    let mut sign = c;
    if rc.is_negative() && !complex {
        sign = sign.neg();
    }
    let mut fs = match rest.kind() {
        Kind::Mul(_, fs) => fs.clone(),
        _ => {
            let (b, e) = as_base_exp(&rest);
            vec![(b, e)]
        }
    };
    fs.push((Expr::constant(inf), Expr::one()));
    fs.sort_unstable_by(cmp_factor);
    Expr::from_kind(Kind::Mul(sign, fs))
}

/// Canonical order of `Mul` factors: by base, then by exponent.
#[inline]
pub(crate) fn cmp_factor(a: &(Expr, Expr), b: &(Expr, Expr)) -> std::cmp::Ordering {
    a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1))
}

fn is_positive_int_base(b: &Expr) -> bool {
    matches!(b.kind(), Kind::Num(n) if n.is_integer() && n.is_positive())
}

fn finish_mul(mut coeff: Number, mut facs: Vec<(Expr, Expr)>, depth: u32) -> Expr {
    facs.sort_unstable_by(cmp_factor);
    // merge equal bases
    let mut merged: Vec<(Expr, Expr, bool)> = Vec::with_capacity(facs.len());
    {
        let mut pending_exps: Vec<Expr> = Vec::new();
        let mut it = facs.into_iter().peekable();
        while let Some((b, e)) = it.next() {
            pending_exps.clear();
            pending_exps.push(e);
            while let Some((nb, _)) = it.peek() {
                if *nb == b {
                    pending_exps.push(it.next().unwrap().1);
                } else {
                    break;
                }
            }
            if pending_exps.len() == 1 {
                merged.push((b, pending_exps.pop().unwrap(), false));
            } else {
                // SymPy only combines exponents sharing the same non-numeric
                // part: x**y*x**(2*y) = x**(3*y) but x*x**y stays.
                let mut groups: Vec<(Expr, Vec<Expr>)> = Vec::new();
                for e in pending_exps.drain(..) {
                    let (_, t) = as_coeff_mul(&e);
                    match groups.iter_mut().find(|(gt, _)| *gt == t) {
                        Some((_, v)) => v.push(e),
                        None => groups.push((t, vec![e])),
                    }
                }
                for (_, mut v) in groups {
                    if v.len() == 1 {
                        merged.push((b.clone(), v.pop().unwrap(), false));
                    } else {
                        merged.push((b.clone(), add(&v), true));
                    }
                }
            }
        }
    }

    let mut out: Vec<(Expr, Expr)> = Vec::with_capacity(merged.len());
    let mut pending: Vec<Expr> = Vec::new();
    let mut radicals: Vec<(Expr, Expr)> = Vec::new();

    // exp(a)*exp(b)*E: SymPy combines the arguments, then splits the sum
    // into one exp per term again (so exp(x)*exp(-x) = 1).
    let n_exp = merged
        .iter()
        .filter(|(b, e, _)| {
            e.is_one() && (b.is_func(&Func::Exp) || b.is_constant(Constant::E))
        })
        .count();
    if n_exp >= 2 {
        let mut exp_args: Vec<Expr> = Vec::with_capacity(n_exp);
        merged.retain(|(b, e, _)| {
            if !e.is_one() {
                return true;
            }
            match b.kind() {
                Kind::Func(Func::Exp, a) => {
                    exp_args.push(a[0].clone());
                    false
                }
                Kind::Const(Constant::E) => {
                    exp_args.push(Expr::one());
                    false
                }
                _ => true,
            }
        });
        let s = add(&exp_args);
        let mut parts: Vec<Expr> = Vec::new();
        match s.kind() {
            Kind::Add(c, terms) => {
                for (t, k) in terms {
                    parts.push(add_term(t, k));
                }
                if !c.is_zero() {
                    parts.push(Expr::num(c.clone()));
                }
            }
            _ => parts.push(s.clone()),
        }
        for p in parts {
            let v = crate::functions::exp(&p);
            match v.kind() {
                Kind::Func(Func::Exp, _) | Kind::Const(Constant::E) => out.push((v, Expr::one())),
                Kind::Num(n) => coeff = coeff.mul(n),
                _ => pending.push(v),
            }
        }
    }

    for (b, e, was_merged) in merged {
        if e.is_zero() {
            continue;
        }
        if e.is_one() && !b.is_number() {
            out.push((b, e));
            continue;
        }
        let special_base = matches!(b.kind(), Kind::Num(_) | Kind::Const(_));
        if !(was_merged || special_base) {
            out.push((b, e));
            continue;
        }
        let v = pow(&b, &e);
        match v.kind() {
            Kind::Pow(vb, ve) if *vb == b => {
                if b.is_number() && ve.is_rational_number() {
                    radicals.push((b, ve.clone()));
                } else {
                    out.push((b, ve.clone()));
                }
            }
            Kind::Num(n) => coeff = coeff.mul(n),
            Kind::Mul(c, fs) if fs.iter().all(|(fb, _)| fb.is_number()) => {
                coeff = coeff.mul(c);
                radicals.extend(fs.iter().cloned());
            }
            _ => pending.push(v),
        }
    }

    // sqrt(2)*sqrt(3) -> sqrt(6): merge positive integer bases with equal
    // rational exponents.
    if radicals.len() > 1 {
        radicals.sort_unstable_by(|a, b| a.1.cmp(&b.1));
        let mut i = 0;
        while i < radicals.len() {
            let mut j = i + 1;
            while j < radicals.len() && radicals[j].1 == radicals[i].1 {
                j += 1;
            }
            let group = &radicals[i..j];
            let ints: Vec<&(Expr, Expr)> =
                group.iter().filter(|(b, _)| is_positive_int_base(b)).collect();
            if ints.len() > 1 {
                let mut prod = Number::ONE;
                for (b, _) in &ints {
                    prod = prod.mul(b.as_number().unwrap());
                }
                let v = pow(&Expr::num(prod), &group[0].1);
                match v.kind() {
                    Kind::Num(n) => coeff = coeff.mul(n),
                    Kind::Mul(c, fs) => {
                        coeff = coeff.mul(c);
                        out.extend(fs.iter().cloned());
                    }
                    Kind::Pow(b, e) => out.push((b.clone(), e.clone())),
                    _ => pending.push(v),
                }
                for (b, e) in group.iter().filter(|(b, _)| !is_positive_int_base(b)) {
                    out.push((b.clone(), e.clone()));
                }
            } else {
                out.extend(group.iter().cloned());
            }
            i = j;
        }
    } else {
        out.extend(radicals);
    }

    if !pending.is_empty() {
        if depth < 8 {
            pending.push(Expr::num(coeff));
            for (b, e) in out {
                pending.push(mul_factor(&b, &e));
            }
            return mul_depth(&pending, depth + 1);
        }
        // give up simplifying further but never drop factors
        for p in pending {
            match p.kind() {
                Kind::Num(n) => coeff = coeff.mul(n),
                Kind::Mul(c, fs) => {
                    coeff = coeff.mul(c);
                    out.extend(fs.iter().cloned());
                }
                _ => out.push(as_base_exp(&p)),
            }
        }
    }

    if coeff.is_zero() {
        return Expr::zero();
    }
    if out.is_empty() {
        return Expr::num(coeff);
    }
    out.sort_unstable_by(cmp_factor);
    if out.len() == 1 {
        let (b, e) = &out[0];
        if coeff.is_one() {
            return mul_factor(b, e);
        }
        // Number * (a + b) distributes
        if e.is_one() {
            if let Kind::Add(c, terms) = b.kind() {
                return distribute(&coeff, c, terms);
            }
        }
    }
    Expr::from_kind(Kind::Mul(coeff, out))
}

fn mul_depth(args: &[Expr], depth: u32) -> Expr {
    let mut coeff = Number::ONE;
    let mut facs: Vec<(Expr, Expr)> = Vec::with_capacity(args.len());
    let mut infs = Vec::new();
    for a in args {
        match a.kind() {
            Kind::Num(n) => coeff = coeff.mul(n),
            Kind::Mul(c, fs) => {
                coeff = coeff.mul(c);
                facs.extend(fs.iter().cloned());
            }
            Kind::Pow(b, e) => facs.push((b.clone(), e.clone())),
            Kind::Const(Constant::NaN) => return Expr::nan(),
            Kind::Const(c) if c.is_infinite() => infs.push(*c),
            _ => facs.push((a.clone(), Expr::one())),
        }
    }
    if !infs.is_empty() {
        return mul_infinite(coeff, facs, &infs);
    }
    if coeff.is_zero() {
        return Expr::zero();
    }
    finish_mul(coeff, facs, depth)
}

fn distribute(k: &Number, c: &Number, terms: &[(Expr, Number)]) -> Expr {
    let terms: Vec<(Expr, Number)> = terms.iter().map(|(t, n)| (t.clone(), n.mul(k))).collect();
    let c = c.mul(k);
    // coefficients may vanish only through float underflow; be safe
    if terms.iter().any(|(_, n)| n.is_zero()) {
        return add_from_terms(c, terms);
    }
    Expr::from_kind(Kind::Add(c, terms))
}

// ----------------------------------------------------------------------
// Pow
// ----------------------------------------------------------------------

/// Canonical power `b ** e`.
pub fn pow(b: &Expr, e: &Expr) -> Expr {
    if let Kind::Num(en) = e.kind() {
        if en.is_zero() && en.is_rational() {
            return Expr::one();
        }
        if en.is_one() {
            return b.clone();
        }
        if let Kind::Num(bn) = b.kind() {
            return num_pow(bn, en);
        }
    }
    if e.is_constant(Constant::NaN) || b.is_constant(Constant::NaN) {
        return Expr::nan();
    }
    match b.kind() {
        Kind::Num(n) if n.is_one() => {
            if matches!(e.kind(), Kind::Const(c) if c.is_infinite()) {
                return Expr::nan();
            }
            return Expr::one();
        }
        Kind::Num(n) if n.is_zero() => {
            if let Kind::Num(en) = e.kind() {
                return if en.is_positive() { Expr::zero() } else { Expr::zoo() };
            }
            if e.is_constant(Constant::Infinity) {
                return Expr::zero();
            }
        }
        Kind::Const(Constant::E) => return crate::functions::exp(e),
        Kind::Const(Constant::I) => {
            if let Some(k) = e.as_number().and_then(|n| n.to_i64()) {
                return i_pow(k);
            }
        }
        Kind::Const(Constant::Infinity) => {
            if let Kind::Num(en) = e.kind() {
                return if en.is_positive() { Expr::oo() } else { Expr::zero() };
            }
        }
        Kind::Const(Constant::NegInfinity) => {
            if let Some(k) = e.as_i64() {
                return if k > 0 {
                    if k % 2 == 0 {
                        Expr::oo()
                    } else {
                        Expr::neg_oo()
                    }
                } else {
                    Expr::zero()
                };
            }
        }
        Kind::Const(Constant::ComplexInfinity) => {
            if let Kind::Num(en) = e.kind() {
                return if en.is_positive() { Expr::zoo() } else { Expr::zero() };
            }
        }
        Kind::Pow(bb, be) => {
            if e.is_integer_number() {
                return pow(bb, &mul2(be, e));
            }
        }
        Kind::Mul(c, facs) => {
            if e.is_integer_number() {
                let mut v = Vec::with_capacity(facs.len() + 1);
                v.push(num_pow(c, e.as_number().unwrap()));
                for (fb, fe) in facs {
                    v.push(pow(fb, &mul2(fe, e)));
                }
                return mul(&v);
            }
            if e.is_number() && !c.is_one() && c.is_rational() {
                let rest = mul_from_factors(facs.clone());
                if c.is_positive() {
                    return mul2(&num_pow(c, e.as_number().unwrap()), &pow_raw(&rest, e));
                } else if !c.is_minus_one() {
                    let neg_rest = Expr::from_kind(Kind::Mul(Number::MINUS_ONE, facs.clone()));
                    return mul2(&num_pow(&c.neg(), e.as_number().unwrap()), &pow_raw(&neg_rest, e));
                }
            }
        }
        Kind::Func(Func::Exp, args) => {
            if e.is_integer_number() {
                return crate::functions::exp(&mul2(&args[0], e));
            }
        }
        _ => {}
    }
    pow_raw(b, e)
}

/// `b ** e` without simplification beyond trivial exponents.
fn pow_raw(b: &Expr, e: &Expr) -> Expr {
    if e.is_one() {
        return b.clone();
    }
    if e.is_zero() {
        return Expr::one();
    }
    Expr::from_kind(Kind::Pow(b.clone(), e.clone()))
}

fn i_pow(k: i64) -> Expr {
    match k.rem_euclid(4) {
        0 => Expr::one(),
        1 => Expr::i(),
        2 => Expr::minus_one(),
        _ => Expr::from_kind(Kind::Mul(Number::MINUS_ONE, vec![(Expr::i(), Expr::one())])),
    }
}

/// Power of two numbers (exact radicals are simplified like SymPy).
pub fn num_pow(b: &Number, e: &Number) -> Expr {
    if e.is_zero() && e.is_rational() {
        return Expr::one();
    }
    if e.is_one() {
        return Expr::num(b.clone());
    }
    if let Some(k) = e.to_i64() {
        return match b.pow_i64(k) {
            Some(v) => Expr::num(v),
            None => Expr::zoo(),
        };
    }
    if e.is_integer() {
        // huge integer exponent
        if b.is_one() {
            return Expr::one();
        }
        if b.is_zero() {
            return if e.is_positive() { Expr::zero() } else { Expr::zoo() };
        }
        if b.is_minus_one() {
            let odd = e.to_bigint().unwrap().is_odd();
            return Expr::int(if odd { -1 } else { 1 });
        }
        return Expr::from_kind(Kind::Pow(Expr::num(b.clone()), Expr::num(e.clone())));
    }
    if b.is_float() || e.is_float() {
        let (x, y) = (b.to_f64(), e.to_f64());
        if x >= 0.0 || y.fract() == 0.0 {
            return Expr::float(x.powf(y));
        }
        let r = (-x).powf(y);
        let t = std::f64::consts::PI * y;
        return add2(
            &Expr::float(r * t.cos()),
            &Expr::from_kind(Kind::Mul(Number::float(r * t.sin()), vec![(Expr::i(), Expr::one())])),
        );
    }
    // exact base, non-integer rational exponent
    if b.is_zero() {
        return if e.is_positive() { Expr::zero() } else { Expr::zoo() };
    }
    if b.is_one() {
        return Expr::one();
    }
    if b.is_negative() {
        if b.is_minus_one() {
            return minus_one_pow(e);
        }
        return mul2(&minus_one_pow(e), &num_pow(&b.neg(), e));
    }
    match b {
        Number::Rational(r) => {
            let n = Number::from_bigint(r.numer().clone());
            let d = Number::from_bigint(r.denom().clone());
            mul2(&num_pow(&n, e), &num_pow(&d, &e.neg()))
        }
        _ => int_rational_pow(&b.to_bigint().unwrap(), e),
    }
}

/// `(-1) ** (p/q)` with the exponent reduced into (-1, 1].
fn minus_one_pow(e: &Number) -> Expr {
    // t = e - 2*round_half_down(e/2)  in (-1, 1]
    let two = Number::int(2);
    let half = e.div(&two).unwrap();
    let mut k = half.floor();
    let mut t = e.sub(&k.mul(&two));
    // t in [0, 2): bring into (-1, 1]
    if t.cmp_value(&Number::ONE) == Some(std::cmp::Ordering::Greater) {
        k = k.add(&Number::ONE);
        t = e.sub(&k.mul(&two));
    }
    if t.is_zero() {
        return Expr::one();
    }
    if t.is_one() {
        return Expr::minus_one();
    }
    let half = Number::rational_i64(1, 2).unwrap();
    if t == half {
        return Expr::i();
    }
    if t == half.neg() {
        return i_pow(3);
    }
    // SymPy keeps the exponent in (0, 1): (-1)**(-1/3) = -(-1)**(2/3)
    if t.is_negative() {
        let t2 = t.add(&Number::ONE);
        return Expr::from_kind(Kind::Mul(
            Number::MINUS_ONE,
            vec![(Expr::minus_one(), Expr::num(t2))],
        ));
    }
    Expr::from_kind(Kind::Pow(Expr::minus_one(), Expr::num(t)))
}

const TRIAL_LIMIT: u64 = 1 << 15;

/// `n ** (p/q)` for a positive integer `n` and a non-integer rational.
fn int_rational_pow(n: &BigInt, e: &Number) -> Expr {
    let (p, q) = e.numer_denom().unwrap();
    // 1. perfect power reduction: n = m**g  =>  n**e = m**(e*g)
    let (factors, cof) = trial_factor(n, TRIAL_LIMIT);
    let mut g: u64 = 0;
    for (_, m) in &factors {
        g = g.gcd(&(*m as u64));
    }
    if !cof.is_one() {
        // the unfactored cofactor: find its largest perfect-power exponent
        let bits = cof.bits();
        let mut k = (bits / 16).clamp(1, 64);
        let mut cof_root_exp = 1u64;
        while k >= 2 {
            if exact_root(&cof, k as u32).is_some() {
                cof_root_exp = k;
                break;
            }
            k -= 1;
        }
        g = if g == 0 { cof_root_exp } else { g.gcd(&cof_root_exp) };
    }
    if g > 1 {
        let m = exact_root(n, g as u32).unwrap();
        let e2 = e.mul(&Number::int(g as i64));
        return num_pow(&Number::from_bigint(m), &e2);
    }
    // 2. split integer part: e = k + r/q
    let k = e.floor();
    let r = p - k.to_bigint().unwrap() * &q;
    let mut outside = Number::from_bigint(n.clone()).pow_i64(k.to_i64().unwrap_or(0)).unwrap();
    // 3. extract q-th powers: n = A**q * M
    let q_u = q.to_u64().unwrap_or(u64::MAX);
    let mut a = BigInt::one();
    let mut m = BigInt::one();
    for (prime, mult) in &factors {
        let mult = *mult as u64;
        let out_e = mult / q_u;
        let in_e = mult % q_u;
        if out_e > 0 {
            a *= num_traits::pow::Pow::pow(prime, out_e);
        }
        if in_e > 0 {
            m *= num_traits::pow::Pow::pow(prime, in_e);
        }
    }
    if !cof.is_one() {
        match q.to_u32().and_then(|qq| exact_root(&cof, qq)) {
            Some(rt) => a *= rt,
            None => m *= &cof,
        }
    }
    if !a.is_one() {
        let ar = num_traits::pow::Pow::pow(&a, r.to_u64().unwrap());
        outside = outside.mul(&Number::from_bigint(ar));
    }
    if m.is_one() {
        return Expr::num(outside);
    }
    let rad_e = Number::rational(r, q).unwrap();
    let rad = (Expr::num(Number::from_bigint(m)), Expr::num(rad_e));
    if outside.is_one() {
        Expr::from_kind(Kind::Pow(rad.0, rad.1))
    } else {
        Expr::from_kind(Kind::Mul(outside, vec![rad]))
    }
}

// ----------------------------------------------------------------------
// derived operations
// ----------------------------------------------------------------------

pub fn neg(a: &Expr) -> Expr {
    match a.kind() {
        Kind::Num(n) => Expr::num(n.neg()),
        Kind::Mul(c, facs) => {
            let c = c.neg();
            if c.is_one() {
                mul_from_factors(facs.clone())
            } else {
                Expr::from_kind(Kind::Mul(c, facs.clone()))
            }
        }
        Kind::Add(c, terms) => distribute(&Number::MINUS_ONE, c, terms),
        _ => mul2(&Expr::minus_one(), a),
    }
}

pub fn sub(a: &Expr, b: &Expr) -> Expr {
    add2(a, &neg(b))
}

pub fn div(a: &Expr, b: &Expr) -> Expr {
    mul2(a, &pow(b, &Expr::minus_one()))
}

pub fn inv(a: &Expr) -> Expr {
    pow(a, &Expr::minus_one())
}

pub fn sqrt(a: &Expr) -> Expr {
    pow(a, &Expr::rational(1, 2))
}

pub fn powi(a: &Expr, n: i64) -> Expr {
    pow(a, &Expr::int(n))
}

/// Is `e` (syntactically) "negative looking", i.e. has a negative leading
/// numeric coefficient (SymPy's `could_extract_minus_sign` heuristic)?
pub fn could_extract_minus_sign(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_negative(),
        Kind::Const(Constant::NegInfinity) => true,
        Kind::Mul(c, _) => c.is_negative(),
        Kind::Add(c, terms) => {
            let neg = terms.iter().filter(|(_, k)| k.is_negative()).count();
            let pos = terms.len() - neg;
            if neg != pos {
                return neg > pos;
            }
            if !c.is_zero() {
                return c.is_negative();
            }
            // tie: decide on the structurally first term for determinism
            let first = terms.iter().min_by(|a, b| a.0.cmp_struct(&b.0)).unwrap();
            first.1.is_negative()
        }
        _ => false,
    }
}

// operator sugar ---------------------------------------------------------

macro_rules! impl_binop {
    ($tr:ident, $m:ident, $f:expr) => {
        impl std::ops::$tr<Expr> for Expr {
            type Output = Expr;
            fn $m(self, o: Expr) -> Expr {
                $f(&self, &o)
            }
        }
        impl std::ops::$tr<&Expr> for &Expr {
            type Output = Expr;
            fn $m(self, o: &Expr) -> Expr {
                $f(self, o)
            }
        }
        impl std::ops::$tr<&Expr> for Expr {
            type Output = Expr;
            fn $m(self, o: &Expr) -> Expr {
                $f(&self, o)
            }
        }
        impl std::ops::$tr<Expr> for &Expr {
            type Output = Expr;
            fn $m(self, o: Expr) -> Expr {
                $f(self, &o)
            }
        }
        impl std::ops::$tr<i64> for Expr {
            type Output = Expr;
            fn $m(self, o: i64) -> Expr {
                $f(&self, &Expr::int(o))
            }
        }
        impl std::ops::$tr<i64> for &Expr {
            type Output = Expr;
            fn $m(self, o: i64) -> Expr {
                $f(self, &Expr::int(o))
            }
        }
        impl std::ops::$tr<Expr> for i64 {
            type Output = Expr;
            fn $m(self, o: Expr) -> Expr {
                $f(&Expr::int(self), &o)
            }
        }
        impl std::ops::$tr<&Expr> for i64 {
            type Output = Expr;
            fn $m(self, o: &Expr) -> Expr {
                $f(&Expr::int(self), o)
            }
        }
    };
}

impl_binop!(Add, add, add2);
impl_binop!(Sub, sub, sub);
impl_binop!(Mul, mul, mul2);
impl_binop!(Div, div, div);

impl std::ops::Neg for Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        neg(&self)
    }
}
impl std::ops::Neg for &Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        neg(self)
    }
}

impl Expr {
    pub fn pow(&self, e: &Expr) -> Expr {
        pow(self, e)
    }
    pub fn powi(&self, n: i64) -> Expr {
        powi(self, n)
    }
    pub fn sqrt(&self) -> Expr {
        sqrt(self)
    }
}
