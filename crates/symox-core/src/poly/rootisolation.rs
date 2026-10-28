//! Certified real and complex root isolation over ZZ and QQ.
//!
//! This is a faithful port of the reference `rootisolation` module: real
//! roots are isolated by the continued-fraction method with Mobius
//! transformations, and the complex roots of a square-free polynomial by
//! exact quadrant-sequence winding counts with rational bisection of the
//! upper half plane. Every step follows the reference algorithm (same root
//! bounds, same shifts, same bisection order, same refinement of edge roots),
//! so the isolating intervals and rectangles, their refinements and the root
//! order are the reference's exact rational values and not only some valid
//! enclosure. No floating-point value ever decides a root: the only float in
//! the module is the reference's logarithm used to pick a shift for the
//! lower root bound, which affects speed and the chosen interval, never the
//! correctness of the certificate.
//!
//! Polynomials are dense with the leading coefficient first (the reference
//! convention), the zero polynomial being the empty vector.

use std::sync::Arc;

use super::int::{Int, Rat};
use super::{factor_uni, sqf};

pub type ZPoly = Vec<Int>;
pub type QPoly = Vec<Rat>;
/// `(a, b, c, d)` for the map `x -> (a*x + b)/(c*x + d)`.
pub type Mobius = (Int, Int, Int, Int);

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IsoError {
    Domain(String),
    Refinement(String),
    Value(String),
    NotImplemented(String),
    ZeroDivision,
    Index,
}
pub type IsoResult<T> = Result<T, IsoError>;

/// Dense coefficient vectors beyond this degree are refused before they are
/// allocated: sparse inputs such as `x**1000000000 + 1` would otherwise
/// exhaust memory, where the reference cannot represent them densely either.
pub const MAX_DENSE_DEGREE: usize = 1 << 20;

fn i(v: i64) -> Int {
    Int::from(v)
}

/// The reference prints domain rationals as `p` or `p/q`.
pub fn rat_str(r: &Rat) -> String {
    if r.d.is_one() {
        format!("{}", r.n)
    } else {
        format!("{}/{}", r.n, r.d)
    }
}

fn frac(n: &Int, d: &Int) -> IsoResult<Rat> {
    if d.is_zero() {
        return Err(IsoError::ZeroDivision);
    }
    Ok(Rat::new(n.clone(), d.clone()))
}

// ---------------------------------------------------------------------------
// Dense integer polynomial primitives (leading coefficient first).

fn zstrip(mut f: ZPoly) -> ZPoly {
    let k = f.iter().take_while(|c| c.is_zero()).count();
    if k > 0 {
        f.drain(..k);
    }
    f
}

fn qstrip(mut f: QPoly) -> QPoly {
    let k = f.iter().take_while(|c| c.is_zero()).count();
    if k > 0 {
        f.drain(..k);
    }
    f
}

fn ztc(f: &[Int]) -> Int {
    f.last().cloned().unwrap_or_else(Int::zero)
}

fn zdegree(f: &[Int]) -> isize {
    f.len() as isize - 1
}

/// `f(x + a)` by the reference's in-place Horner scheme.
fn zshift(f: &[Int], a: &Int) -> ZPoly {
    let mut f = f.to_vec();
    let n = f.len().saturating_sub(1);
    if a.is_one() {
        for k in (1..=n).rev() {
            for j in 0..k {
                let (l, r) = f.split_at_mut(j + 1);
                r[0].add_assign(&l[j]);
            }
        }
    } else {
        for k in (1..=n).rev() {
            for j in 0..k {
                let (l, r) = f.split_at_mut(j + 1);
                r[0].add_mul(a, &l[j]);
            }
        }
    }
    f
}

/// `f(a*x)`.
fn zscale(f: &[Int], a: &Int) -> ZPoly {
    let mut f = f.to_vec();
    if f.len() < 2 {
        return f;
    }
    let n = f.len() - 1;
    let mut b = a.clone();
    for k in (0..n).rev() {
        f[k] = b.mul(&f[k]);
        b = b.mul(a);
    }
    f
}

/// `f(-x)`.
fn zmirror(f: &[Int]) -> ZPoly {
    let mut f = f.to_vec();
    let n = f.len();
    let mut k = n as isize - 2;
    while k >= 0 {
        f[k as usize] = f[k as usize].neg();
        k -= 2;
    }
    f
}

fn zreverse(f: &[Int]) -> ZPoly {
    zstrip(f.iter().rev().cloned().collect())
}

fn zrshift1(mut f: ZPoly) -> ZPoly {
    f.pop();
    f
}

fn sign_variations(f: &[Int]) -> i64 {
    let mut prev = 0i32;
    let mut k = 0;
    for c in f {
        let s = c.signum();
        if s * prev < 0 {
            k += 1;
        }
        if s != 0 {
            prev = s;
        }
    }
    k
}

fn zmul(a: &[Int], b: &[Int]) -> ZPoly {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut out = vec![Int::zero(); a.len() + b.len() - 1];
    for (x, ca) in a.iter().enumerate() {
        if ca.is_zero() {
            continue;
        }
        for (y, cb) in b.iter().enumerate() {
            out[x + y].add_mul(ca, cb);
        }
    }
    zstrip(out)
}

fn zadd(a: &[Int], b: &[Int]) -> ZPoly {
    let n = a.len().max(b.len());
    let mut out = vec![Int::zero(); n];
    for (k, c) in a.iter().enumerate() {
        out[n - a.len() + k].add_assign(c);
    }
    for (k, c) in b.iter().enumerate() {
        out[n - b.len() + k].add_assign(c);
    }
    zstrip(out)
}

/// `q**n * f(p/q)` for `n = deg f`.
fn ztransform(f: &[Int], p: &[Int], q: &[Int]) -> ZPoly {
    if f.is_empty() {
        return vec![];
    }
    let n = f.len() - 1;
    let mut h = vec![f[0].clone()];
    let mut powers = vec![vec![Int::one()]];
    for _ in 0..n {
        let next = zmul(powers.last().unwrap(), q);
        powers.push(next);
    }
    for (c, qk) in f[1..].iter().zip(powers[1..].iter()) {
        h = zmul(&h, p);
        let term: ZPoly = if c.is_zero() {
            vec![]
        } else {
            zstrip(qk.iter().map(|x| x.mul(c)).collect())
        };
        h = zadd(&h, &term);
    }
    h
}

fn zterms_gcd(f: &[Int]) -> (u32, ZPoly) {
    if f.is_empty() || !ztc(f).is_zero() {
        return (0, f.to_vec());
    }
    let j = f.iter().rev().take_while(|c| c.is_zero()).count();
    (j as u32, f[..f.len() - j].to_vec())
}

/// `lcm(denominators) * f` with integer coefficients.
pub fn clear_denoms(f: &[Rat]) -> ZPoly {
    let mut common = Int::one();
    for c in f {
        common = common.lcm(&c.d);
    }
    f.iter().map(|c| c.n.mul(&common).div_exact_unchecked(&c.d)).collect()
}

fn to_ascending(f: &[Int]) -> Vec<Int> {
    f.iter().rev().cloned().collect()
}

fn to_descending(f: &[Int]) -> ZPoly {
    f.iter().rev().cloned().collect()
}

/// Square-free factors over ZZ by increasing multiplicity, like the
/// reference's `dup_sqf_list` (contents and signs do not affect isolation).
fn zsqf_list(f: &[Int]) -> Vec<(ZPoly, u32)> {
    let (_, parts) = sqf::sqf_list_zx(&to_ascending(f));
    parts
        .into_iter()
        .map(|(g, k)| {
            let (_, g) = super::dense::primitive(&g);
            let mut g = to_descending(&g);
            if g.first().is_some_and(Int::is_negative) {
                g = g.iter().map(Int::neg).collect();
            }
            (g, k)
        })
        .collect()
}

/// Irreducible factors over ZZ in the reference's canonical order: by
/// length, multiplicity, then coefficient list.
fn zfactor_list(f: &[Int]) -> Vec<(ZPoly, u32)> {
    if f.len() < 2 {
        return vec![];
    }
    let (_, factors) = factor_uni::factor_zx(&to_ascending(f));
    let mut out: Vec<(ZPoly, u32)> = factors
        .into_iter()
        .map(|(g, k)| {
            let mut g = to_descending(&g);
            if g[0].is_negative() {
                g = g.iter().map(Int::neg).collect();
            }
            (g, k)
        })
        .collect();
    out.sort_by(|(f, n), (g, m)| (f.len(), *n, f).cmp(&(g.len(), *m, g)));
    out
}

fn qeval(f: &[Rat], x: &Rat) -> Rat {
    let mut r = Rat::zero();
    for c in f {
        r = r.mul(x).add(c);
    }
    r
}

// ---------------------------------------------------------------------------
// Root bounds.

/// `int(math.log(a, 2))` for an integer `a > 0`, reproducing CPython's
/// float evaluation: correctly rounded conversion when `a < 2**1024`,
/// otherwise the `frexp` decomposition.
fn py_log2(a: &Int) -> i64 {
    let ln2 = 2.0f64.ln();
    let num = match a {
        Int::S(v) => (*v as f64).ln(),
        Int::B(_) => {
            let big = a.to_big();
            let mag = big.magnitude();
            let bits = mag.bits();
            if bits < 64 {
                let v = mag.iter_u64_digits().next().unwrap_or(0);
                return ((v as f64).ln() / ln2) as i64;
            }
            let shift = bits - 54;
            let top = (mag >> shift).iter_u64_digits().next().unwrap_or(0);
            let sticky = mag.trailing_zeros().map_or(false, |z| z < shift);
            let mut m = top >> 1;
            let mut exponent = shift as i64 + 1;
            if top & 1 == 1 && (sticky || m & 1 == 1) {
                m += 1;
                if m == 1u64 << 53 {
                    m >>= 1;
                    exponent += 1;
                }
            }
            // value = m * 2**exponent with m < 2**53: frexp exponent is
            // exponent + 53.
            let e = exponent + 53;
            if e <= 1024 {
                (m as f64 * 2f64.powi(exponent as i32)).ln()
            } else {
                let x = m as f64 / (1u64 << 53) as f64;
                x.ln() + ln2 * e as f64
            }
        }
    };
    (num / ln2) as i64
}

fn pow2(k: i64) -> Rat {
    if k >= 0 {
        Rat::from(i(2).pow(k as u32))
    } else {
        Rat::new(Int::one(), i(2).pow((-k) as u32))
    }
}

/// Upper bound on the positive roots (the reference's local-max bound).
pub fn root_upper_bound(f: &[Int]) -> Option<Rat> {
    let n = f.len();
    let mut t = vec![1i64; n];
    let negate = f.first().is_some_and(|c| c.is_negative());
    let g: Vec<Int> = f.iter().rev().map(|c| if negate { c.neg() } else { c.clone() }).collect();
    let mut p: Vec<i64> = Vec::new();
    for k in 0..n {
        if !g[k].is_negative() {
            continue;
        }
        let a = py_log2(&g[k].neg());
        let mut best: Option<(i64, usize)> = None;
        for j in k + 1..n {
            if !g[j].is_positive() {
                continue;
            }
            let q = t[j] + a - py_log2(&g[j]);
            let cand = (q.div_euclid((j - k) as i64), j);
            if best.is_none_or(|b| cand < b) {
                best = Some(cand);
            }
        }
        if let Some((q, j)) = best {
            t[j] += 1;
            p.push(q);
        }
    }
    p.iter().max().map(|m| pow2(m + 1))
}

pub fn root_lower_bound(f: &[Int]) -> Option<Rat> {
    root_upper_bound(&zreverse(f)).map(|b| b.inv())
}

// ---------------------------------------------------------------------------
// Real roots: continued fractions.

fn mobius_from_interval(s: &Rat, t: &Rat) -> Mobius {
    (s.n.clone(), t.n.clone(), s.d.clone(), t.d.clone())
}

fn mobius_to_interval(m: &Mobius) -> IsoResult<(Rat, Rat)> {
    let s = frac(&m.0, &m.2)?;
    let t = frac(&m.1, &m.3)?;
    Ok(if s <= t { (s, t) } else { (t, s) })
}

fn lower_shift(f: &[Int]) -> Int {
    root_lower_bound(f).map(|r| r.floor()).unwrap_or_else(Int::zero)
}

fn step_refine(f: ZPoly, m: Mobius, fast: bool) -> (ZPoly, Mobius) {
    let (mut a, mut b, mut c, mut d) = m;
    if a == b && c == d {
        return (f, (a, b, c, d));
    }
    let mut f = f;
    let mut shift = lower_shift(&f);
    if fast && shift > i(16) {
        f = zscale(&f, &shift);
        a = shift.mul(&a);
        c = shift.mul(&c);
        shift = Int::one();
    }
    if shift >= Int::one() {
        f = zshift(&f, &shift);
        b = shift.mul(&a).add(&b);
        d = shift.mul(&c).add(&d);
        if ztc(&f).is_zero() {
            return (f, (b.clone(), b, d.clone(), d));
        }
    }
    let g = f;
    let f = zshift(&g, &Int::one());
    let (a1, b1, c1, d1) = (a.clone(), a.add(&b), c.clone(), c.add(&d));
    if ztc(&f).is_zero() {
        return (f, (b1.clone(), b1, d1.clone(), d1));
    }
    if sign_variations(&f) == 1 {
        (f, (a1, b1, c1, d1))
    } else {
        let mut f = zshift(&zreverse(&g), &Int::one());
        if ztc(&f).is_zero() {
            f = zrshift1(f);
        }
        (f, (b.clone(), a.add(&b), d.clone(), c.add(&d)))
    }
}

fn width(m: &Mobius) -> IsoResult<Rat> {
    Ok(frac(&m.0, &m.2)?.sub(&frac(&m.1, &m.3)?).abs())
}

fn inner_refine(
    f: ZPoly,
    m: Mobius,
    eps: Option<&Rat>,
    steps: Option<u64>,
    disjoint: Option<&Rat>,
    fast: bool,
) -> IsoResult<(ZPoly, Mobius)> {
    let (mut f, mut m) = (f, m);
    while m.2.is_zero() {
        (f, m) = step_refine(f, m, fast);
    }
    match (eps, steps) {
        (Some(eps), Some(steps)) => {
            for _ in 0..steps {
                if &width(&m)? >= eps {
                    (f, m) = step_refine(f, m, fast);
                } else {
                    break;
                }
            }
        }
        _ => {
            if let Some(eps) = eps {
                while &width(&m)? >= eps {
                    (f, m) = step_refine(f, m, fast);
                }
            }
            if let Some(steps) = steps {
                for _ in 0..steps {
                    (f, m) = step_refine(f, m, fast);
                }
            }
        }
    }
    if let Some(x) = disjoint {
        loop {
            let (u, v) = mobius_to_interval(&m)?;
            if &v <= x || x <= &u {
                break;
            }
            (f, m) = step_refine(f, m, fast);
        }
    }
    Ok((f, m))
}

fn outer_refine(
    f: &[Int],
    s: &Rat,
    t: &Rat,
    eps: Option<&Rat>,
    steps: Option<u64>,
    disjoint: Option<&Rat>,
    fast: bool,
) -> IsoResult<(Rat, Rat)> {
    let m = mobius_from_interval(s, t);
    let f = ztransform(f, &zstrip(vec![m.0.clone(), m.1.clone()]), &zstrip(vec![m.2.clone(), m.3.clone()]));
    if sign_variations(&f) != 1 {
        return Err(IsoError::Refinement(format!(
            "there should be exactly one root in ({}, {}) interval",
            rat_str(s),
            rat_str(t)
        )));
    }
    let (_, m) = inner_refine(f, m, eps, steps, disjoint, fast)?;
    mobius_to_interval(&m)
}

/// Refine the isolating interval `(s, t)` of a real root of `f` over ZZ
/// (the reference's `dup_refine_real_root`).
pub fn refine_real_root(
    f: &[Int],
    s: &Rat,
    t: &Rat,
    eps: Option<&Rat>,
    steps: Option<u64>,
    disjoint: Option<&Rat>,
    fast: bool,
) -> IsoResult<(Rat, Rat)> {
    if s == t {
        return Ok((s.clone(), t.clone()));
    }
    let (mut s, mut t) = if s > t { (t.clone(), s.clone()) } else { (s.clone(), t.clone()) };
    let mut f = f.to_vec();
    let mut negative = false;
    if s.is_negative() {
        if !t.is_positive_or_zero_rat() {
            f = zmirror(&f);
            (s, t) = (t.neg(), s.neg());
            negative = true;
        } else {
            return Err(IsoError::Value(format!(
                "Cannot refine a real root in ({}, {})",
                rat_str(&s),
                rat_str(&t)
            )));
        }
    }
    let mut disjoint = disjoint.cloned();
    if negative {
        disjoint = match disjoint {
            Some(x) if x.is_negative() => Some(x.neg()),
            _ => None,
        };
    }
    let (s, t) = outer_refine(&f, &s, &t, eps, steps, disjoint.as_ref(), fast)?;
    Ok(if negative { (t.neg(), s.neg()) } else { (s, t) })
}

trait RatSign {
    fn is_positive_or_zero_rat(&self) -> bool;
}
impl RatSign for Rat {
    /// `t > 0` (the reference tests `t <= 0` for the mirrored branch).
    fn is_positive_or_zero_rat(&self) -> bool {
        self.n.is_positive()
    }
}

type Found = (ZPoly, Mobius);

fn inner_isolate_real_roots(f: &[Int], eps: Option<&Rat>, fast: bool) -> IsoResult<Vec<Found>> {
    let one = Int::one;
    let k = sign_variations(f);
    if k == 0 {
        return Ok(vec![]);
    }
    if k == 1 {
        return Ok(vec![inner_refine(f.to_vec(), (one(), Int::zero(), Int::zero(), one()), eps, None, None, fast)?]);
    }
    let mut roots = Vec::new();
    let mut stack: Vec<(Int, Int, Int, Int, ZPoly, i64)> =
        vec![(one(), Int::zero(), Int::zero(), one(), f.to_vec(), k)];
    while let Some((mut a, mut b, mut c, mut d, mut f, mut k)) = stack.pop() {
        let mut shift = lower_shift(&f);
        if fast && shift > i(16) {
            f = zscale(&f, &shift);
            a = shift.mul(&a);
            c = shift.mul(&c);
            shift = one();
        }
        if shift >= one() {
            f = zshift(&f, &shift);
            b = shift.mul(&a).add(&b);
            d = shift.mul(&c).add(&d);
            if ztc(&f).is_zero() {
                roots.push((f.clone(), (b.clone(), b.clone(), d.clone(), d.clone())));
                f = zrshift1(f);
            }
            k = sign_variations(&f);
            if k == 0 {
                continue;
            }
            if k == 1 {
                roots.push(inner_refine(f, (a, b, c, d), eps, None, None, fast)?);
                continue;
            }
        }
        let mut f1 = Some(zshift(&f, &one()));
        let (mut a1, mut b1, mut c1, mut d1, mut r) = (a.clone(), a.add(&b), c.clone(), c.add(&d), 0);
        if ztc(f1.as_ref().unwrap()).is_zero() {
            roots.push((f1.clone().unwrap(), (b1.clone(), b1.clone(), d1.clone(), d1.clone())));
            f1 = Some(zrshift1(f1.unwrap()));
            r = 1;
        }
        let mut k1 = sign_variations(f1.as_ref().unwrap());
        let mut k2 = k - k1 - r;
        let (mut a2, mut b2, mut c2, mut d2) = (b.clone(), a.add(&b), d.clone(), c.add(&d));
        let reversed_shift = |f: &[Int]| {
            let mut g = zshift(&zreverse(f), &Int::one());
            if ztc(&g).is_zero() {
                g = zrshift1(g);
            }
            g
        };
        let mut f2 = None;
        if k2 > 1 {
            let g = reversed_shift(&f);
            k2 = sign_variations(&g);
            f2 = Some(g);
        }
        if k1 < k2 {
            std::mem::swap(&mut a1, &mut a2);
            std::mem::swap(&mut b1, &mut b2);
            std::mem::swap(&mut c1, &mut c2);
            std::mem::swap(&mut d1, &mut d2);
            std::mem::swap(&mut f1, &mut f2);
            std::mem::swap(&mut k1, &mut k2);
        }
        if k1 == 0 {
            continue;
        }
        let f1 = f1.unwrap_or_else(|| reversed_shift(&f));
        if k1 == 1 {
            roots.push(inner_refine(f1, (a1, b1, c1, d1), eps, None, None, fast)?);
        } else {
            stack.push((a1, b1, c1, d1, f1, k1));
        }
        if k2 == 0 {
            continue;
        }
        let f2 = f2.unwrap_or_else(|| reversed_shift(&f));
        if k2 == 1 {
            roots.push(inner_refine(f2, (a2, b2, c2, d2), eps, None, None, fast)?);
        } else {
            stack.push((a2, b2, c2, d2, f2, k2));
        }
    }
    Ok(roots)
}

fn discard_if_outside(
    f: ZPoly,
    m: Mobius,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    negative: bool,
    fast: bool,
) -> IsoResult<Option<Found>> {
    let (mut f, mut m) = (f, m);
    loop {
        let (mut u, mut v) = mobius_to_interval(&m)?;
        if negative {
            (u, v) = (v.neg(), u.neg());
        }
        if inf.is_none_or(|x| &u >= x) && sup.is_none_or(|x| &v <= x) {
            return Ok(Some((f, m)));
        } else if sup.is_some_and(|x| &u > x) || inf.is_some_and(|x| &v < x) {
            return Ok(None);
        }
        (f, m) = step_refine(f, m, fast);
    }
}

fn isolate_positive(
    f: &[Int],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    fast: bool,
) -> IsoResult<Vec<Found>> {
    if sup.is_some_and(|s| s.is_negative()) {
        return Ok(vec![]);
    }
    let roots = inner_isolate_real_roots(f, eps, fast)?;
    if inf.is_none() && sup.is_none() {
        return Ok(roots);
    }
    let mut out = Vec::new();
    for (g, m) in roots {
        if let Some(r) = discard_if_outside(g, m, inf, sup, false, fast)? {
            out.push(r);
        }
    }
    Ok(out)
}

fn isolate_negative(
    f: &[Int],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    fast: bool,
) -> IsoResult<Vec<Found>> {
    if inf.is_some_and(|s| !s.is_negative()) {
        return Ok(vec![]);
    }
    let roots = inner_isolate_real_roots(&zmirror(f), eps, fast)?;
    if inf.is_none() && sup.is_none() {
        return Ok(roots);
    }
    let mut out = Vec::new();
    for (g, m) in roots {
        if let Some(r) = discard_if_outside(g, m, inf, sup, true, fast)? {
            out.push(r);
        }
    }
    Ok(out)
}

fn positive_interval(m: &Mobius) -> IsoResult<(Rat, Rat)> {
    mobius_to_interval(m)
}

fn negative_interval(m: &Mobius) -> IsoResult<(Rat, Rat)> {
    let (u, v) = mobius_to_interval(m)?;
    Ok((v.neg(), u.neg()))
}

fn zero_in(inf: Option<&Rat>, sup: Option<&Rat>) -> bool {
    inf.is_none_or(|x| !x.is_positive_rat()) && sup.is_none_or(|x| !x.is_negative())
}

trait RatPos {
    fn is_positive_rat(&self) -> bool;
}
impl RatPos for Rat {
    fn is_positive_rat(&self) -> bool {
        self.n.is_positive()
    }
}

/// Isolating intervals of the real roots of a square-free `f` over ZZ,
/// sorted, together with `f` without its power of `x` (the polynomial the
/// reference stores in its blackbox intervals).
pub fn isolate_real_roots_sqf(
    f: &[Int],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    fast: bool,
) -> IsoResult<(Vec<(Rat, Rat)>, ZPoly)> {
    if zdegree(f) <= 0 {
        return Ok((vec![], f.to_vec()));
    }
    let (j, f) = zterms_gcd(f);
    let mut roots = Vec::new();
    for (_, m) in isolate_negative(&f, eps, inf, sup, fast)? {
        roots.push(negative_interval(&m)?);
    }
    if j > 0 && zero_in(inf, sup) {
        roots.push((Rat::zero(), Rat::zero()));
    }
    for (_, m) in isolate_positive(&f, eps, inf, sup, fast)? {
        roots.push(positive_interval(&m)?);
    }
    roots.sort();
    Ok((roots, f))
}

/// Isolating intervals with multiplicities of the real roots of `f` over ZZ
/// (the reference's `dup_isolate_real_roots`).
pub fn isolate_real_roots(
    f: &[Int],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    fast: bool,
) -> IsoResult<Vec<((Rat, Rat), u32)>> {
    if zdegree(f) <= 0 {
        return Ok(vec![]);
    }
    let (j, f) = zterms_gcd(f);
    let mut out = Vec::new();
    if j > 0 && zero_in(inf, sup) {
        out.push(((Rat::zero(), Rat::zero()), j));
    }
    let factors = zsqf_list(&f);
    if factors.len() == 1 {
        let (g, k) = &factors[0];
        for (_, m) in isolate_negative(g, eps, inf, sup, fast)? {
            out.push((negative_interval(&m)?, *k));
        }
        for (_, m) in isolate_positive(g, eps, inf, sup, fast)? {
            out.push((positive_interval(&m)?, *k));
        }
    } else {
        let (neg, pos) = real_isolate_and_disjoin(factors, eps, inf, sup, false, fast)?;
        for (u, v, k, _) in neg.into_iter().chain(pos) {
            out.push(((u, v), k));
        }
    }
    out.sort_by(|x, y| x.0.cmp(&y.0));
    Ok(out)
}

fn disjoint_p(m: &Mobius, n: &Mobius, strict: bool) -> bool {
    let (mut a1, mut b1, mut c1, mut d1) = m.clone();
    let (mut a2, mut b2, mut c2, mut d2) = n.clone();
    let (a1d1, b1c1) = (a1.mul(&d1), b1.mul(&c1));
    let (a2d2, b2c2) = (a2.mul(&d2), b2.mul(&c2));
    if a1d1 == b1c1 && a2d2 == b2c2 {
        return true;
    }
    if a1d1 > b1c1 {
        std::mem::swap(&mut a1, &mut b1);
        std::mem::swap(&mut c1, &mut d1);
    }
    if a2d2 > b2c2 {
        std::mem::swap(&mut a2, &mut b2);
        std::mem::swap(&mut c2, &mut d2);
    }
    if !strict {
        a2.mul(&d1) >= c2.mul(&b1) || b2.mul(&c1) <= d2.mul(&a1)
    } else {
        a2.mul(&d1) > c2.mul(&b1) || b2.mul(&c1) < d2.mul(&a1)
    }
}

type Tagged<T> = (ZPoly, Mobius, T, ZPoly);

fn disjoin<T>(list: &mut [Tagged<T>], strict: bool, fast: bool) -> IsoResult<()> {
    for x in 0..list.len() {
        for y in x + 1..list.len() {
            while !disjoint_p(&list[x].1, &list[y].1, strict) {
                let (f, m) = inner_refine(std::mem::take(&mut list[x].0), list[x].1.clone(), None, Some(1), None, fast)?;
                list[x].0 = f;
                list[x].1 = m;
                let (g, n) = inner_refine(std::mem::take(&mut list[y].0), list[y].1.clone(), None, Some(1), None, fast)?;
                list[y].0 = g;
                list[y].1 = n;
            }
        }
    }
    Ok(())
}

fn strict_away_from_zero<T>(list: &mut [Tagged<T>], fast: bool) -> IsoResult<()> {
    for entry in list.iter_mut() {
        if entry.1 .0.is_zero() {
            while entry.1 .0.is_zero() {
                let (f, m) = inner_refine(std::mem::take(&mut entry.0), entry.1.clone(), None, Some(1), None, fast)?;
                entry.0 = f;
                entry.1 = m;
            }
            break;
        }
    }
    Ok(())
}

type Disjoined<T> = Vec<(Rat, Rat, T, ZPoly)>;

fn real_isolate_and_disjoin<T: Clone>(
    factors: Vec<(ZPoly, T)>,
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    strict: bool,
    fast: bool,
) -> IsoResult<(Disjoined<T>, Disjoined<T>)> {
    let mut pos: Vec<Tagged<T>> = Vec::new();
    let mut neg: Vec<Tagged<T>> = Vec::new();
    for (f, k) in &factors {
        for (g, m) in isolate_positive(f, eps, inf, sup, fast)? {
            pos.push((g, m, k.clone(), f.clone()));
        }
        for (g, m) in isolate_negative(f, eps, inf, sup, fast)? {
            neg.push((g, m, k.clone(), f.clone()));
        }
    }
    disjoin(&mut pos, strict, fast)?;
    disjoin(&mut neg, strict, fast)?;
    if strict {
        strict_away_from_zero(&mut neg, fast)?;
        strict_away_from_zero(&mut pos, fast)?;
    }
    let mut out_neg = Vec::new();
    for (_, m, k, f) in neg {
        let (u, v) = negative_interval(&m)?;
        out_neg.push((u, v, k, f));
    }
    let mut out_pos = Vec::new();
    for (_, m, k, f) in pos {
        let (u, v) = positive_interval(&m)?;
        out_pos.push((u, v, k, f));
    }
    Ok((out_neg, out_pos))
}

/// An isolating interval of simultaneous real root isolation: the endpoints
/// in stored order (reversed on the top and left rectangle edges), the
/// multiplicity of the root in each input polynomial and the irreducible
/// factor it is a root of.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EdgeRoot {
    pub s: Rat,
    pub t: Rat,
    pub indices: Vec<(usize, u32)>,
    pub h: ZPoly,
}

impl EdgeRoot {
    fn index(&self, k: usize) -> Option<u32> {
        self.indices.iter().find(|(j, _)| *j == k).map(|(_, m)| *m)
    }
}

/// Simultaneous isolation of the real roots of several polynomials over ZZ
/// (the reference's `dup_isolate_real_roots_list` with `basis=True`).
pub fn isolate_real_roots_list(
    polys: &[ZPoly],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    strict: bool,
    fast: bool,
) -> IsoResult<Vec<EdgeRoot>> {
    let zeros = zero_in(inf, sup);
    let mut zero_indices: Vec<(usize, u32)> = Vec::new();
    let mut factors: Vec<(ZPoly, Vec<(usize, u32)>)> = Vec::new();
    for (k, p) in polys.iter().enumerate() {
        let (j, p) = zterms_gcd(p);
        if zeros && j > 0 {
            zero_indices.push((k, j));
        }
        for (f, m) in zfactor_list(&p) {
            match factors.iter_mut().find(|(g, _)| *g == f) {
                Some((_, idx)) => match idx.iter_mut().find(|(j, _)| *j == k) {
                    Some(e) => e.1 = m,
                    None => idx.push((k, m)),
                },
                None => factors.push((f, vec![(k, m)])),
            }
        }
    }
    let (neg, pos) = real_isolate_and_disjoin(factors, eps, inf, sup, strict, fast)?;
    let mut out: Vec<EdgeRoot> = Vec::new();
    for (s, t, indices, h) in neg {
        out.push(EdgeRoot { s, t, indices, h });
    }
    if zeros && !zero_indices.is_empty() {
        out.push(EdgeRoot { s: Rat::zero(), t: Rat::zero(), indices: zero_indices, h: vec![Int::one(), Int::zero()] });
    }
    for (s, t, indices, h) in pos {
        out.push(EdgeRoot { s, t, indices, h });
    }
    out.sort_by(|x, y| (&x.s, &x.t).cmp(&(&y.s, &y.t)));
    Ok(out)
}

/// The Sturm sequence of the monic square-free part of `f` over QQ, leading
/// coefficients first (the reference's `dup_sturm`): the zero polynomial
/// gives `[0]` and a nonzero constant gives `[1]`.
pub fn sturm_sequence(f: &[Rat]) -> Vec<QPoly> {
    use super::dense;
    let asc: Vec<Rat> = f.iter().rev().cloned().collect();
    let asc = {
        let mut a = asc;
        while a.last().is_some_and(|c| c.is_zero()) {
            a.pop();
        }
        a
    };
    let mut p = if asc.is_empty() {
        asc
    } else {
        let g = dense::gcd_field(&asc, &dense::derivative(&asc));
        dense::divrem(&asc, &g).0
    };
    if let Some(lc) = p.last().cloned() {
        for c in p.iter_mut() {
            *c = c.div(&lc);
        }
    }
    let mut chain = vec![p.clone(), dense::derivative(&p)];
    while !chain.last().unwrap().is_empty() {
        let n = chain.len();
        let r = dense::rem(&chain[n - 2], &chain[n - 1]);
        chain.push(dense::neg(&r));
    }
    chain.pop();
    chain.into_iter().map(|s| s.into_iter().rev().collect()).collect()
}

/// Number of real roots of `f` in `[inf, sup]` by Sturm sequences (the
/// reference's `dup_count_real_roots`).
///
/// The reference builds the Sturm sequence of the monic square-free part
/// over QQ. Here every element is replaced by a positive multiple with
/// integer coefficients (primitive pseudo-remainders, signs tracked), which
/// leaves every sign, hence every variation count, unchanged while avoiding
/// rational normalizations. Signs at rational points come from homogeneous
/// integer Horner evaluation. The result is the reference's exact count.
pub fn count_real_roots(f: &[Rat], inf: Option<&Rat>, sup: Option<&Rat>) -> usize {
    if f.len() < 2 {
        return 0;
    }
    // positive multiple of f with integer coefficients, ascending
    let asc: Vec<Int> = clear_denoms(f).into_iter().rev().collect();
    let chain = integer_sturm_chain(&asc);
    let variations = |signs: Vec<i32>| {
        let mut prev = 0i32;
        let mut k = 0i64;
        for s in signs {
            if s * prev < 0 {
                k += 1;
            }
            if s != 0 {
                prev = s;
            }
        }
        k
    };
    let at = |x: Option<&Rat>, minus: bool| -> Vec<i32> {
        chain
            .iter()
            .map(|s| match x {
                Some(x) => sign_at(s, x),
                None => {
                    let lc = s.last().map_or(0, Int::signum);
                    if minus && (s.len() - 1) % 2 == 1 {
                        -lc
                    } else {
                        lc
                    }
                }
            })
            .collect()
    };
    let mut count = (variations(at(inf, true)) - variations(at(sup, false))).unsigned_abs() as usize;
    if let Some(x) = inf {
        if sign_at(&asc, x) == 0 {
            count += 1;
        }
    }
    count
}

/// Sign of the ascending integer polynomial `f` at the rational `x`, from
/// the homogeneous form `q**n f(p/q)` (`q > 0`).
fn sign_at(f: &[Int], x: &Rat) -> i32 {
    let Some(top) = f.last() else { return 0 };
    if x.n.is_zero() {
        return f[0].signum();
    }
    let (p, q) = (&x.n, &x.d);
    let mut h = top.clone();
    if q.is_one() {
        for c in f[..f.len() - 1].iter().rev() {
            h = h.mul(p);
            if !c.is_zero() {
                h.add_assign(c);
            }
        }
    } else {
        let mut qk = Int::one();
        for c in f[..f.len() - 1].iter().rev() {
            h = h.mul(p);
            qk = qk.mul(q);
            if !c.is_zero() {
                h.add_mul(c, &qk);
            }
        }
    }
    h.signum()
}

/// Ascending integer polynomials `T_k = c_k S_k` (`c_k > 0`) for the
/// reference's Sturm sequence `S_k` of the monic square-free part of `f`.
fn integer_sturm_chain(f: &[Int]) -> Vec<Vec<Int>> {
    use super::dense;
    let df = dense::derivative(f);
    let (_, sqf, _) = super::gcd::gcd_zx(f, &df);
    let mut p = positive_primitive(&sqf);
    if p.last().is_some_and(Int::is_negative) {
        p = dense::neg(&p);
    }
    let dp = positive_primitive(&dense::derivative(&p));
    let mut chain = vec![p, dp];
    while !chain.last().unwrap().is_empty() {
        let n = chain.len();
        let b = &chain[n - 1];
        if b.len() == 1 {
            break;
        }
        let a = &chain[n - 2];
        let r = lazy_prem(a, b);
        // prem(a, b) = lc(b)**(da - db + 1) rem(a, b); the next element is
        // -rem(a, b) up to a positive factor.
        let steps = a.len() + 1 - b.len();
        let negative_scale = b.last().unwrap().is_negative() && steps % 2 == 1;
        let r = positive_primitive(&r);
        chain.push(if negative_scale { r } else { dense::neg(&r) });
    }
    if chain.last().unwrap().is_empty() {
        chain.pop();
    }
    chain
}

/// `p / content(p)` with a positive divisor (signs preserved).
fn positive_primitive(p: &[Int]) -> Vec<Int> {
    let g = super::dense::content(p).abs();
    if g.is_zero() || g.is_one() {
        return p.to_vec();
    }
    p.iter().map(|c| c.div_exact_unchecked(&g)).collect()
}

/// The pseudo-remainder `lc(b)**(da - db + 1) a mod b` of ascending integer
/// polynomials (as the classical algorithm), scaling coefficients lazily:
/// a coefficient that the subtraction does not touch is multiplied by the
/// accumulated power of `lc(b)` only when it is next read. A sparse divisor
/// (such as a linear Sturm element) thus costs one product per step
/// instead of rescaling the whole remainder.
fn lazy_prem(a: &[Int], b: &[Int]) -> Vec<Int> {
    use super::dense;
    let m = b.len() - 1;
    if a.len() < b.len() {
        return a.to_vec();
    }
    let n = a.len() - 1;
    let lb = &b[m];
    let terms: Vec<(usize, &Int)> = b[..m].iter().enumerate().filter(|(_, c)| !c.is_zero()).collect();
    let steps = n - m + 1;
    let mut r = a.to_vec();
    // true value of r[k] after t steps: r[k] * lb**(t - stamp[k])
    let mut stamp = vec![0usize; n + 1];
    let mut cache: Option<(usize, Int)> = None;
    let mut power = |e: usize| -> Int {
        match &cache {
            Some((k, v)) if *k == e => return v.clone(),
            Some((k, v)) if *k + 1 == e => {
                let v = v.mul(lb);
                cache = Some((e, v.clone()));
                return v;
            }
            _ => {}
        }
        let v = lb.pow(e as u32);
        cache = Some((e, v.clone()));
        v
    };
    fn scaled(c: &Int, e: usize, lb: &Int, power: &mut dyn FnMut(usize) -> Int) -> Int {
        match e {
            _ if c.is_zero() => Int::zero(),
            0 => c.clone(),
            1 => c.mul(lb),
            _ => c.mul(&power(e)),
        }
    }
    for t in 0..steps {
        let d = n - t;
        let c = scaled(&r[d], t - stamp[d], lb, &mut power);
        r[d] = Int::zero();
        if c.is_zero() {
            continue;
        }
        for &(j, bj) in &terms {
            let k = d - m + j;
            let mut v = scaled(&r[k], t + 1 - stamp[k], lb, &mut power);
            v.sub_mul(&c, bj);
            r[k] = v;
            stamp[k] = t + 1;
        }
    }
    r.truncate(m);
    for k in 0..m {
        let e = steps - stamp[k];
        r[k] = scaled(&r[k], e, lb, &mut power);
    }
    dense::trim(&mut r);
    r
}

// ---------------------------------------------------------------------------
// Complex roots: quadrant sequences on rectangle edges.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Qd {
    OO,
    Q1,
    Q2,
    Q3,
    Q4,
    A1,
    A2,
    A3,
    A4,
}

impl Qd {
    pub fn from_name(name: &str) -> Option<Qd> {
        Some(match name {
            "OO" => Qd::OO,
            "Q1" => Qd::Q1,
            "Q2" => Qd::Q2,
            "Q3" => Qd::Q3,
            "Q4" => Qd::Q4,
            "A1" => Qd::A1,
            "A2" => Qd::A2,
            "A3" => Qd::A3,
            "A4" => Qd::A4,
            _ => return None,
        })
    }
    pub fn name(self) -> &'static str {
        match self {
            Qd::OO => "OO",
            Qd::Q1 => "Q1",
            Qd::Q2 => "Q2",
            Qd::Q3 => "Q3",
            Qd::Q4 => "Q4",
            Qd::A1 => "A1",
            Qd::A2 => "A2",
            Qd::A3 => "A3",
            Qd::A4 => "A4",
        }
    }
}

fn rule_simple(k: &[Qd]) -> Option<i32> {
    use Qd::*;
    Some(match k {
        [Q1, Q1] | [Q2, Q2] | [Q3, Q3] | [Q4, Q4] => 0,
        [A1, Q1] | [A2, Q2] | [A3, Q3] | [A4, Q4] => 1,
        [A1, Q4] | [A2, Q1] | [A3, Q2] | [A4, Q3] => 2,
        [Q1, A2] | [Q2, A3] | [Q3, A4] | [Q4, A1] => 3,
        [Q1, A1] | [Q2, A2] | [Q3, A3] | [Q4, A4] => 4,
        [Q1, Q2] | [Q2, Q3] | [Q3, Q4] | [Q4, Q1] => 5,
        [Q1, Q4] | [Q2, Q1] | [Q3, Q2] | [Q4, Q3] => -5,
        _ => return None,
    })
}

fn rule_ambiguous(k: &[Qd]) -> Option<i32> {
    use Qd::*;
    Some(match k {
        [A1, A3] | [A2, A4] | [A3, A1] | [A4, A2] => 7,
        [Q1, Q3] | [Q2, Q4] | [Q3, Q1] | [Q4, Q2] => 8,
        [A1, A2] | [A2, A3] | [A3, A4] | [A4, A1] => 9,
        [A1, A4] | [A2, A1] | [A3, A2] | [A4, A3] => 10,
        [Q1, A3] | [Q2, A4] | [Q3, A1] | [Q4, A2] => 11,
        [Q1, A4] | [Q2, A1] | [Q3, A2] | [Q4, A3] => 12,
        [A1, Q3] | [A2, Q4] | [A3, Q1] | [A4, Q2] => 13,
        [A1, Q2] | [A2, Q3] | [A3, Q4] | [A4, Q1] => 14,
        [A1, OO, Q1] | [A2, OO, Q2] | [A3, OO, Q3] | [A4, OO, Q4] => -1,
        [A1, OO, Q4] | [A2, OO, Q1] | [A3, OO, Q2] | [A4, OO, Q3] => -2,
        [Q1, OO, A2] | [Q2, OO, A3] | [Q3, OO, A4] | [Q4, OO, A1] => -3,
        [Q1, OO, A1] | [Q2, OO, A2] | [Q3, OO, A3] | [Q4, OO, A4] => -4,
        [A1, OO, A3] | [A2, OO, A4] | [A3, OO, A1] | [A4, OO, A2] => 7,
        [Q1, OO, Q3] | [Q2, OO, Q4] | [Q3, OO, Q1] | [Q4, OO, Q2] => 8,
        [A1, OO, A2] | [A2, OO, A3] | [A3, OO, A4] | [A4, OO, A1] => 9,
        [A1, OO, A4] | [A2, OO, A1] | [A3, OO, A2] | [A4, OO, A3] => 10,
        [Q1, OO, A3] | [Q2, OO, A4] | [Q3, OO, A1] | [Q4, OO, A2] => 11,
        [Q1, OO, A4] | [Q2, OO, A1] | [Q3, OO, A2] | [Q4, OO, A3] => 12,
        [A1, OO, Q3] | [A2, OO, Q4] | [A3, OO, Q1] | [A4, OO, Q2] => 13,
        [A1, OO, Q2] | [A2, OO, Q3] | [A3, OO, Q4] | [A4, OO, Q1] => 14,
        [Q1, OO, Q2] | [Q2, OO, Q3] | [Q3, OO, Q4] | [Q4, OO, Q1] => 15,
        [Q1, OO, Q4] | [Q2, OO, Q1] | [Q3, OO, Q2] | [Q4, OO, Q3] => 16,
        [A1, OO, A1] | [A2, OO, A2] | [A3, OO, A3] | [A4, OO, A4] => 17,
        [Q1, OO, Q1] | [Q2, OO, Q2] | [Q3, OO, Q3] | [Q4, OO, Q4] => 18,
        _ => return None,
    })
}

fn rule_values(t: i32) -> &'static [(i64, i64)] {
    match t {
        0 => &[(0, 1)],
        1 | 3 => &[(1, 4)],
        2 | 4 => &[(-1, 4)],
        -1 | -3 => &[(9, 4), (1, 4)],
        -2 | -4 => &[(7, 4), (-1, 4)],
        5 => &[(1, 2)],
        -5 => &[(-1, 2)],
        7 | 8 => &[(1, 1), (-1, 1)],
        9 | 15 => &[(1, 2), (-3, 2)],
        10 | 16 => &[(3, 2), (-1, 2)],
        11 | 14 => &[(3, 4), (-5, 4)],
        12 | 13 => &[(5, 4), (-3, 4)],
        17 | 18 => &[(2, 1), (0, 1)],
        _ => &[],
    }
}

fn rule_text(k: &[Qd]) -> String {
    format!("({})", k.iter().map(|q| format!("'{}'", q.name())).collect::<Vec<_>>().join(", "))
}

fn classify_point(re: &Rat, im: &Rat) -> Option<Qd> {
    match (re.is_zero(), im.is_zero()) {
        (true, true) => Some(Qd::OO),
        (true, false) => Some(if im.is_negative() { Qd::A4 } else { Qd::A2 }),
        (false, true) => Some(if re.is_negative() { Qd::A3 } else { Qd::A1 }),
        _ => None,
    }
}

fn positive(r: &Rat) -> bool {
    r.n.is_positive()
}

fn intervals_to_quadrants(intervals: &[EdgeRoot], f1: &[Rat], f2: &[Rat], s: &Rat, t: &Rat) -> IsoResult<Vec<Qd>> {
    if intervals.is_empty() {
        return Ok(vec![]);
    }
    let two = Rat::from_i64(2);
    let mut q = Vec::new();
    // One part identically zero on the edge: follow the sign of the other.
    for (zero, other, plus, minus, which) in [(f1, f2, Qd::A2, Qd::A4, 1usize), (f2, f1, Qd::A1, Qd::A3, 0usize)] {
        if !zero.is_empty() {
            continue;
        }
        let mut list = intervals;
        let first = &intervals[0];
        let mut sgn;
        if first.s == first.t && &first.s == s {
            if intervals.len() == 1 {
                return Ok(if positive(&qeval(other, t)) { vec![Qd::OO, plus] } else { vec![Qd::OO, minus] });
            }
            let a = &intervals[1].s;
            if positive(&qeval(other, &s.add(a).div(&two))) {
                q.extend([Qd::OO, plus]);
                sgn = 1;
            } else {
                q.extend([Qd::OO, minus]);
                sgn = -1;
            }
            list = &intervals[1..];
        } else if positive(&qeval(other, s)) {
            q.push(plus);
            sgn = 1;
        } else {
            q.push(minus);
            sgn = -1;
        }
        for root in list {
            q.push(Qd::OO);
            let m = root.index(which).ok_or(IsoError::Index)?;
            if m % 2 == 1 {
                sgn = -sgn;
            }
            if &root.s != t {
                q.push(if sgn > 0 { plus } else { minus });
            }
        }
        return Ok(q);
    }
    let mut list = intervals;
    let mut re = qeval(f1, s);
    let mut im = qeval(f2, s);
    if re.is_zero() || im.is_zero() {
        q.push(classify_point(&re, &im).ok_or(IsoError::Index)?);
        if intervals.len() == 1 {
            re = qeval(f1, t);
            im = qeval(f2, t);
        } else {
            let x = s.add(&intervals[1].s).div(&two);
            re = qeval(f1, &x);
            im = qeval(f2, &x);
        }
        list = &intervals[1..];
    }
    let mut s1 = if positive(&re) { 1 } else { -1 };
    let mut s2 = if positive(&im) { 1 } else { -1 };
    let sgn = |a: i32, b: i32| match (a, b) {
        (1, 1) => Qd::Q1,
        (-1, 1) => Qd::Q2,
        (-1, -1) => Qd::Q3,
        _ => Qd::Q4,
    };
    q.push(sgn(s1, s2));
    for root in list {
        let (a, b) = (&root.s, &root.t);
        if a == b {
            let re = qeval(f1, a);
            let im = qeval(f2, a);
            if let Some(c) = classify_point(&re, &im) {
                q.push(c);
            }
        }
        if let Some(m) = root.index(0) {
            if m % 2 == 1 {
                s1 = -s1;
            }
        }
        if let Some(m) = root.index(1) {
            if m % 2 == 1 {
                s2 = -s2;
            }
        }
        if !(a == b && b == t) {
            q.push(sgn(s1, s2));
        }
    }
    Ok(q)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Exclude {
    None,
    All,
    Some { edges: [bool; 4], corners: [bool; 4] },
}

impl Exclude {
    /// Edge names `S, E, N, W` and corner names `SW, SE, NE, NW`.
    pub fn from_names(names: &[String]) -> Exclude {
        let mut edges = [false; 4];
        let mut corners = [false; 4];
        for (k, e) in ["S", "E", "N", "W"].iter().enumerate() {
            if names.iter().any(|n| n == e) {
                edges[k] = true;
            }
        }
        for (k, c) in ["SW", "SE", "NE", "NW"].iter().enumerate() {
            if names.iter().any(|n| n == c) {
                corners[k] = true;
            }
        }
        Exclude::Some { edges, corners }
    }
}

/// Rules for the quadrant sequences of the four edges; `corner[k]` is the
/// flag of the corner where edge `k` starts.
fn traverse_quadrants(qs: [&[Qd]; 4], exclude: &Exclude) -> IsoResult<Vec<(i32, usize)>> {
    let (edges, corners) = match exclude {
        Exclude::All => ([1, 1, 0, 0], [1, 1, 1, 0]),
        Exclude::None => ([0; 4], [0; 4]),
        Exclude::Some { edges, corners } => (edges.map(usize::from), corners.map(usize::from)),
    };
    let mut rules = Vec::new();
    for (k, full) in qs.iter().enumerate() {
        let mut q: &[Qd] = full;
        if q.is_empty() {
            continue;
        }
        if q[q.len() - 1] == Qd::OO {
            q = &q[..q.len() - 1];
        }
        if q.first() == Some(&Qd::OO) || q.is_empty() {
            let j = (k + 3) % 4;
            q = q.get(1..).ok_or(IsoError::Index)?;
            let prev = qs[j];
            let key = [
                *prev.get(prev.len().wrapping_sub(2)).ok_or(IsoError::Index)?,
                Qd::OO,
                *q.first().ok_or(IsoError::Index)?,
            ];
            let rule = rule_ambiguous(&key)
                .ok_or_else(|| IsoError::NotImplemented(format!("3 element rule (corner): {}", rule_text(&key))))?;
            rules.push((rule, corners[k]));
        }
        let mut q1 = *q.first().ok_or(IsoError::Index)?;
        let mut idx = 1;
        while idx < q.len() {
            let q2 = q[idx];
            idx += 1;
            if q2 != Qd::OO {
                let key = [q1, q2];
                if let Some(rule) = rule_simple(&key) {
                    rules.push((rule, 0));
                } else if let Some(rule) = rule_ambiguous(&key) {
                    rules.push((rule, edges[k]));
                } else {
                    return Err(IsoError::NotImplemented(format!("2 element rule (inside): {}", rule_text(&key))));
                }
                q1 = q2;
            } else {
                let key = [q1, q2, *q.get(idx).ok_or(IsoError::Index)?];
                idx += 1;
                let rule = rule_ambiguous(&key)
                    .ok_or_else(|| IsoError::NotImplemented(format!("3 element rule (edge): {}", rule_text(&key))))?;
                rules.push((rule, edges[k]));
                q1 = key[2];
            }
        }
    }
    Ok(rules)
}

fn winding_number(rules: &[(i32, usize)]) -> IsoResult<i64> {
    let mut total = Rat::zero();
    for (t, k) in rules {
        let (n, d) = *rule_values(*t).get(*k).ok_or(IsoError::Index)?;
        total = total.add(&Rat::new(i(n), i(d)));
    }
    let half = total.div(&Rat::from_i64(2));
    // `int()` truncates toward zero.
    let (q, _) = half.n.div_rem_trunc(&half.d);
    Ok(q.as_i64().unwrap_or(0))
}

/// Real and imaginary parts of `f(x + i*y0)` as polynomials in `x`.
fn horizontal_edge(f: &[Rat], y0: &Rat) -> (QPoly, QPoly) {
    let mut re: QPoly = vec![];
    let mut im: QPoly = vec![];
    for c in f {
        // (re + i im) * (x + i y0) + c
        let n = re.len().max(im.len()) + 1;
        let mut nre = vec![Rat::zero(); n];
        let mut nim = vec![Rat::zero(); n];
        let off_re = n - 1 - re.len();
        let off_im = n - 1 - im.len();
        for (k, a) in re.iter().enumerate() {
            nre[off_re + k] = nre[off_re + k].add(a);
            nim[off_re + k + 1] = nim[off_re + k + 1].add(&a.mul(y0));
        }
        for (k, b) in im.iter().enumerate() {
            nim[off_im + k] = nim[off_im + k].add(b);
            nre[off_im + k + 1] = nre[off_im + k + 1].sub(&b.mul(y0));
        }
        nre[n - 1] = nre[n - 1].add(c);
        re = qstrip(nre);
        im = qstrip(nim);
    }
    (re, im)
}

/// Real and imaginary parts of `f(x0 + i*y)` as polynomials in `y`.
fn vertical_edge(f: &[Rat], x0: &Rat) -> (QPoly, QPoly) {
    let mut re: QPoly = vec![];
    let mut im: QPoly = vec![];
    for c in f {
        // (re + i im) * (x0 + i y) + c = (x0 re - y im + c) + i (x0 im + y re)
        let n = re.len().max(im.len()) + 1;
        let mut nre = vec![Rat::zero(); n];
        let mut nim = vec![Rat::zero(); n];
        let off_re = n - 1 - re.len();
        let off_im = n - 1 - im.len();
        for (k, a) in re.iter().enumerate() {
            nre[off_re + k + 1] = nre[off_re + k + 1].add(&a.mul(x0));
            nim[off_re + k] = nim[off_re + k].add(a);
        }
        for (k, b) in im.iter().enumerate() {
            nre[off_im + k] = nre[off_im + k].sub(b);
            nim[off_im + k + 1] = nim[off_im + k + 1].add(&b.mul(x0));
        }
        nre[n - 1] = nre[n - 1].add(c);
        re = qstrip(nre);
        im = qstrip(nim);
    }
    (re, im)
}

/// The reference's `dup_real_imag`: `f(x + i*y) = f1(x, y) + i*f2(x, y)`
/// as dense polynomials in `x` whose coefficients are dense in `y`.
pub fn real_imag(f: &[Rat]) -> (Vec<QPoly>, Vec<QPoly>) {
    if f.is_empty() {
        return (vec![], vec![]);
    }
    let n = f.len() - 1;
    // coefficient tables c[p][q] of x^p y^q
    let mut c1 = vec![vec![Rat::zero(); n + 1]; n + 1];
    let mut c2 = vec![vec![Rat::zero(); n + 1]; n + 1];
    for (pos, a) in f.iter().enumerate() {
        let j = n - pos;
        let mut binom = Int::one();
        for q in 0..=j {
            if q > 0 {
                binom = binom.mul(&i((j - q + 1) as i64)).div_exact_unchecked(&i(q as i64));
            }
            let term = a.mul_int(&binom);
            let p = j - q;
            match q % 4 {
                0 => c1[p][q] = c1[p][q].add(&term),
                1 => c2[p][q] = c2[p][q].add(&term),
                2 => c1[p][q] = c1[p][q].sub(&term),
                _ => c2[p][q] = c2[p][q].sub(&term),
            }
        }
    }
    let to_dmp = |c: Vec<Vec<Rat>>| -> Vec<QPoly> {
        let mut out: Vec<QPoly> = (0..=n).rev().map(|p| qstrip(c[p].iter().rev().cloned().collect())).collect();
        let k = out.iter().take_while(|g| g.is_empty()).count();
        out.drain(..k);
        out
    };
    (to_dmp(c1), to_dmp(c2))
}

fn edge_roots(f1: &[Rat], f2: &[Rat], inf: &Rat, sup: &Rat) -> IsoResult<Vec<EdgeRoot>> {
    isolate_real_roots_list(&[clear_denoms(f1), clear_denoms(f2)], None, Some(inf), Some(sup), true, true)
}

fn reverse_intervals(list: &[EdgeRoot]) -> Vec<EdgeRoot> {
    list.iter()
        .rev()
        .map(|r| EdgeRoot { s: r.t.clone(), t: r.s.clone(), indices: r.indices.clone(), h: r.h.clone() })
        .collect()
}

/// A rectangle `[a, b]` of the isolation, with the roots of the real and
/// imaginary parts on its edges (bottom, right, top, left), their quadrant
/// sequences and the edge restrictions of the real and imaginary parts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rect {
    pub n: i64,
    pub a: (Rat, Rat),
    pub b: (Rat, Rat),
    pub i: [Vec<EdgeRoot>; 4],
    pub q: [Vec<Qd>; 4],
    pub f1: [QPoly; 4],
    pub f2: [QPoly; 4],
}

fn split_edge(list: &[EdgeRoot], x: &Rat, reversed: bool) -> IsoResult<(Vec<EdgeRoot>, Vec<EdgeRoot>)> {
    let (mut lo, mut hi) = (Vec::new(), Vec::new());
    for r in list {
        let (a, b) = if reversed { (&r.t, &r.s) } else { (&r.s, &r.t) };
        if a == b {
            if a == x {
                lo.push(r.clone());
                hi.push(r.clone());
            } else if a < x {
                lo.push(r.clone());
            } else {
                hi.push(r.clone());
            }
        } else if b <= x {
            lo.push(r.clone());
        } else if a >= x {
            hi.push(r.clone());
        } else {
            let (a, b) = refine_real_root(&r.h, a, b, None, None, Some(x), true)?;
            let (s, t) = if reversed { (b.clone(), a.clone()) } else { (a.clone(), b.clone()) };
            let fresh = EdgeRoot { s, t, indices: r.indices.clone(), h: r.h.clone() };
            if &b <= x {
                lo.push(fresh.clone());
            }
            if &a >= x {
                hi.push(fresh);
            }
        }
    }
    Ok((lo, hi))
}

fn vertical_bisection(r: &Rect, f: &[Rat]) -> IsoResult<(Rect, Rect)> {
    let ((u, v), (s, t)) = (&r.a, &r.b);
    let x = u.add(s).div(&Rat::from_i64(2));
    let (f1v, f2v) = vertical_edge(f, &x);
    let iv = edge_roots(&f1v, &f2v, v, t)?;
    let (i1l, i1r) = split_edge(&r.i[0], &x, false)?;
    let (i3l, i3r) = split_edge(&r.i[2], &x, true)?;
    let i2l = iv.clone();
    let i2r = r.i[1].clone();
    let i4l = r.i[3].clone();
    let i4r = reverse_intervals(&iv);
    let q1l = intervals_to_quadrants(&i1l, &r.f1[0], &r.f2[0], u, &x)?;
    let q2l = intervals_to_quadrants(&i2l, &f1v, &f2v, v, t)?;
    let q3l = intervals_to_quadrants(&i3l, &r.f1[2], &r.f2[2], &x, u)?;
    let q4l = r.q[3].clone();
    let q1r = intervals_to_quadrants(&i1r, &r.f1[0], &r.f2[0], &x, s)?;
    let q2r = r.q[1].clone();
    let q3r = intervals_to_quadrants(&i3r, &r.f1[2], &r.f2[2], s, &x)?;
    let q4r = intervals_to_quadrants(&i4r, &f1v, &f2v, t, v)?;
    let nl = winding_number(&traverse_quadrants([&q1l, &q2l, &q3l, &q4l], &Exclude::All)?)?;
    let nr = winding_number(&traverse_quadrants([&q1r, &q2r, &q3r, &q4r], &Exclude::All)?)?;
    let left = Rect {
        n: nl,
        a: (u.clone(), v.clone()),
        b: (x.clone(), t.clone()),
        i: [i1l, i2l, i3l, i4l],
        q: [q1l, q2l, q3l, q4l],
        f1: [r.f1[0].clone(), f1v.clone(), r.f1[2].clone(), r.f1[3].clone()],
        f2: [r.f2[0].clone(), f2v.clone(), r.f2[2].clone(), r.f2[3].clone()],
    };
    let right = Rect {
        n: nr,
        a: (x, v.clone()),
        b: (s.clone(), t.clone()),
        i: [i1r, i2r, i3r, i4r],
        q: [q1r, q2r, q3r, q4r],
        f1: [r.f1[0].clone(), r.f1[1].clone(), r.f1[2].clone(), f1v],
        f2: [r.f2[0].clone(), r.f2[1].clone(), r.f2[2].clone(), f2v],
    };
    Ok((left, right))
}

fn horizontal_bisection(r: &Rect, f: &[Rat]) -> IsoResult<(Rect, Rect)> {
    let ((u, v), (s, t)) = (&r.a, &r.b);
    let y = v.add(t).div(&Rat::from_i64(2));
    let (f1h, f2h) = horizontal_edge(f, &y);
    let ih = edge_roots(&f1h, &f2h, u, s)?;
    let i1b = r.i[0].clone();
    let i1u = ih.clone();
    let (i2b, i2u) = split_edge(&r.i[1], &y, false)?;
    let i3b = reverse_intervals(&ih);
    let i3u = r.i[2].clone();
    let (i4b, i4u) = split_edge(&r.i[3], &y, true)?;
    let q1b = r.q[0].clone();
    let q2b = intervals_to_quadrants(&i2b, &r.f1[1], &r.f2[1], v, &y)?;
    let q3b = intervals_to_quadrants(&i3b, &f1h, &f2h, s, u)?;
    let q4b = intervals_to_quadrants(&i4b, &r.f1[3], &r.f2[3], &y, v)?;
    let q1u = intervals_to_quadrants(&i1u, &f1h, &f2h, u, s)?;
    let q2u = intervals_to_quadrants(&i2u, &r.f1[1], &r.f2[1], &y, t)?;
    let q3u = r.q[2].clone();
    let q4u = intervals_to_quadrants(&i4u, &r.f1[3], &r.f2[3], t, &y)?;
    let nb = winding_number(&traverse_quadrants([&q1b, &q2b, &q3b, &q4b], &Exclude::All)?)?;
    let nu = winding_number(&traverse_quadrants([&q1u, &q2u, &q3u, &q4u], &Exclude::All)?)?;
    let bottom = Rect {
        n: nb,
        a: (u.clone(), v.clone()),
        b: (s.clone(), y.clone()),
        i: [i1b, i2b, i3b, i4b],
        q: [q1b, q2b, q3b, q4b],
        f1: [r.f1[0].clone(), r.f1[1].clone(), f1h.clone(), r.f1[3].clone()],
        f2: [r.f2[0].clone(), r.f2[1].clone(), f2h.clone(), r.f2[3].clone()],
    };
    let upper = Rect {
        n: nu,
        a: (u.clone(), y),
        b: (s.clone(), t.clone()),
        i: [i1u, i2u, i3u, i4u],
        q: [q1u, q2u, q3u, q4u],
        f1: [f1h, r.f1[1].clone(), r.f1[2].clone(), r.f1[3].clone()],
        f2: [f2h, r.f2[1].clone(), r.f2[2].clone(), r.f2[3].clone()],
    };
    Ok((bottom, upper))
}

/// An isolating rectangle of a non-real root of a square-free polynomial;
/// `conj` selects the mirrored rectangle in the lower half plane.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ComplexBox {
    pub rect: Rect,
    pub f: Arc<QPoly>,
    pub conj: bool,
}

impl ComplexBox {
    pub fn ax(&self) -> Rat {
        self.rect.a.0.clone()
    }
    pub fn bx(&self) -> Rat {
        self.rect.b.0.clone()
    }
    pub fn ay(&self) -> Rat {
        if self.conj {
            self.rect.b.1.neg()
        } else {
            self.rect.a.1.clone()
        }
    }
    pub fn by(&self) -> Rat {
        if self.conj {
            self.rect.a.1.neg()
        } else {
            self.rect.b.1.clone()
        }
    }
    pub fn dx(&self) -> Rat {
        self.rect.b.0.sub(&self.rect.a.0)
    }
    pub fn dy(&self) -> Rat {
        self.rect.b.1.sub(&self.rect.a.1)
    }
    pub fn conjugate(&self) -> ComplexBox {
        ComplexBox { rect: self.rect.clone(), f: self.f.clone(), conj: true }
    }
    /// One bisection keeping the half that contains the root.
    pub fn refine(&self) -> IsoResult<ComplexBox> {
        let r = &self.rect;
        let ((u, v), (s, t)) = (&r.a, &r.b);
        let rect = if s.sub(u) > t.sub(v) {
            let (l, rr) = vertical_bisection(r, &self.f)?;
            if l.n == 1 {
                l
            } else {
                rr
            }
        } else {
            let (b, up) = horizontal_bisection(r, &self.f)?;
            if b.n == 1 {
                b
            } else {
                up
            }
        };
        Ok(ComplexBox { rect: ComplexBox::with_unit(rect), f: self.f.clone(), conj: self.conj })
    }
    fn with_unit(mut r: Rect) -> Rect {
        // The reference keeps the rectangle data without the count.
        r.n = 1;
        r
    }
    pub fn refine_size(&self, dx: &Rat, dy: &Rat) -> IsoResult<ComplexBox> {
        let mut b = self.clone();
        while !(&b.dx() < dx && &b.dy() < dy) {
            b = b.refine()?;
        }
        Ok(b)
    }
    pub fn is_disjoint(&self, other: &ComplexBox) -> bool {
        if self.conj != other.conj {
            return true;
        }
        if self.bx() < other.ax() || other.bx() < self.ax() {
            return true;
        }
        self.by() < other.ay() || other.by() < self.ay()
    }
}

fn root_bound_2(f: &[Rat]) -> Rat {
    let lc = f[0].abs();
    let mut best = Rat::zero();
    for c in f {
        let q = c.abs().div(&lc);
        if q > best {
            best = q;
        }
    }
    best.mul(&Rat::from_i64(2))
}

/// Isolating rectangles of the non-real roots of a square-free `f` over QQ,
/// each upper-half-plane rectangle preceded by its conjugate (the
/// reference's `dup_isolate_complex_roots_sqf` with `blackbox=True`).
pub fn isolate_complex_roots_sqf(
    f: &[Rat],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
) -> IsoResult<Vec<ComplexBox>> {
    let f = qstrip(f.to_vec());
    if f.len() < 2 {
        return Ok(vec![]);
    }
    let bound = root_bound_2(&f);
    let (mut u, v, mut s, t) = (bound.neg(), Rat::zero(), bound.clone(), bound);
    if let Some(x) = inf {
        u = x.clone();
    }
    if let Some(x) = sup {
        s = x.clone();
    }
    if v.is_negative() || t <= v || s <= u {
        return Err(IsoError::Value("not a valid complex isolation rectangle".into()));
    }
    let (f1l1, f2l1) = horizontal_edge(&f, &v);
    let (f1l2, f2l2) = vertical_edge(&f, &s);
    let (f1l3, f2l3) = horizontal_edge(&f, &t);
    let (f1l4, f2l4) = vertical_edge(&f, &u);
    let i1 = edge_roots(&f1l1, &f2l1, &u, &s)?;
    let i2 = edge_roots(&f1l2, &f2l2, &v, &t)?;
    let i3 = reverse_intervals(&edge_roots(&f1l3, &f2l3, &u, &s)?);
    let i4 = reverse_intervals(&edge_roots(&f1l4, &f2l4, &v, &t)?);
    let q1 = intervals_to_quadrants(&i1, &f1l1, &f2l1, &u, &s)?;
    let q2 = intervals_to_quadrants(&i2, &f1l2, &f2l2, &v, &t)?;
    let q3 = intervals_to_quadrants(&i3, &f1l3, &f2l3, &s, &u)?;
    let q4 = intervals_to_quadrants(&i4, &f1l4, &f2l4, &t, &v)?;
    let n = winding_number(&traverse_quadrants([&q1, &q2, &q3, &q4], &Exclude::None)?)?;
    if n == 0 {
        return Ok(vec![]);
    }
    let fa = Arc::new(f);
    let mut rectangles = vec![Rect {
        n,
        a: (u, v),
        b: (s, t),
        i: [i1, i2, i3, i4],
        q: [q1, q2, q3, q4],
        f1: [f1l1, f1l2, f1l3, f1l4],
        f2: [f2l1, f2l2, f2l3, f2l4],
    }];
    let small = |r: &Rect| match eps {
        Some(e) => &r.b.0.sub(&r.a.0) < e && &r.b.1.sub(&r.a.1) < e,
        None => true,
    };
    let mut roots: Vec<Rect> = Vec::new();
    while !rectangles.is_empty() {
        // depth-first: the smallest area, first one on ties
        let mut j = 0;
        let mut min_area: Option<Rat> = None;
        for (k, r) in rectangles.iter().enumerate() {
            let area = r.b.0.sub(&r.a.0).mul(&r.b.1.sub(&r.a.1));
            if min_area.as_ref().is_none_or(|m| &area < m) {
                min_area = Some(area);
                j = k;
            }
        }
        let r = rectangles.remove(j);
        let (first, second) = if r.b.0.sub(&r.a.0) > r.b.1.sub(&r.a.1) {
            vertical_bisection(&r, &fa)?
        } else {
            horizontal_bisection(&r, &fa)?
        };
        for d in [first, second] {
            if d.n >= 1 {
                if d.n == 1 && small(&d) {
                    roots.push(d);
                } else {
                    rectangles.push(d);
                }
            }
        }
    }
    roots.sort_by(|x, y| (&x.a.0, &x.a.1).cmp(&(&y.a.0, &y.a.1)));
    let mut out = Vec::with_capacity(2 * roots.len());
    for r in roots {
        let b = ComplexBox { rect: r, f: fa.clone(), conj: false };
        out.push(b.conjugate());
        out.push(b);
    }
    Ok(out)
}

/// Number of roots of `f` over QQ in the closed rectangle `[inf, sup]` (the
/// reference's `dup_count_complex_roots`).
pub fn count_complex_roots(
    f: &[Rat],
    inf: Option<(Rat, Rat)>,
    sup: Option<(Rat, Rat)>,
    exclude: &Exclude,
) -> IsoResult<i64> {
    let f = qstrip(f.to_vec());
    if f.is_empty() {
        return Err(IsoError::ZeroDivision);
    }
    let bound = if inf.is_none() || sup.is_none() { Some(root_bound_2(&f)) } else { None };
    let (u, v) = inf.unwrap_or_else(|| {
        let b = bound.clone().unwrap();
        (b.neg(), b.neg())
    });
    let (s, t) = sup.unwrap_or_else(|| {
        let b = bound.clone().unwrap();
        (b.clone(), b)
    });
    let (f1l1, f2l1) = horizontal_edge(&f, &v);
    let (f1l2, f2l2) = vertical_edge(&f, &s);
    let (f1l3, f2l3) = horizontal_edge(&f, &t);
    let (f1l4, f2l4) = vertical_edge(&f, &u);
    let i1 = edge_roots(&f1l1, &f2l1, &u, &s)?;
    let i2 = edge_roots(&f1l2, &f2l2, &v, &t)?;
    let i3 = reverse_intervals(&edge_roots(&f1l3, &f2l3, &u, &s)?);
    let i4 = reverse_intervals(&edge_roots(&f1l4, &f2l4, &v, &t)?);
    let q1 = intervals_to_quadrants(&i1, &f1l1, &f2l1, &u, &s)?;
    let q2 = intervals_to_quadrants(&i2, &f1l2, &f2l2, &v, &t)?;
    let q3 = intervals_to_quadrants(&i3, &f1l3, &f2l3, &s, &u)?;
    let q4 = intervals_to_quadrants(&i4, &f1l4, &f2l4, &t, &v)?;
    winding_number(&traverse_quadrants([&q1, &q2, &q3, &q4], exclude)?)
}

/// Real intervals and complex rectangles of a square-free `f` over QQ (the
/// reference's `dup_isolate_all_roots_sqf`).
pub fn isolate_all_roots_sqf(
    f: &[Rat],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    fast: bool,
) -> IsoResult<(Vec<(Rat, Rat)>, Vec<ComplexBox>)> {
    let (real, _) = isolate_real_roots_sqf(&clear_denoms(f), eps, inf, sup, fast)?;
    let complex = isolate_complex_roots_sqf(f, eps, inf, sup)?;
    Ok((real, complex))
}

/// The reference's `dup_isolate_all_roots`: only a power of one square-free
/// polynomial is supported.
#[allow(clippy::type_complexity)]
pub fn isolate_all_roots(
    f: &[Rat],
    eps: Option<&Rat>,
    inf: Option<&Rat>,
    sup: Option<&Rat>,
    fast: bool,
) -> IsoResult<(Vec<((Rat, Rat), u32)>, Vec<(ComplexBox, u32)>)> {
    let f = qstrip(f.to_vec());
    let factors = if f.is_empty() { vec![] } else { zsqf_list(&clear_denoms(&f)) };
    if factors.len() != 1 {
        return Err(IsoError::NotImplemented("only trivial square-free polynomials are supported".into()));
    }
    let (g, k) = &factors[0];
    let g: QPoly = g.iter().map(|c| Rat::from(c.clone())).collect();
    let (real, complex) = isolate_all_roots_sqf(&g, eps, inf, sup, fast)?;
    Ok((real.into_iter().map(|r| (r, *k)).collect(), complex.into_iter().map(|c| (c, *k)).collect()))
}

/// Whether `f` over QQ is square-free.
pub fn is_squarefree(f: &[Rat]) -> bool {
    use super::dense;
    let asc: Vec<Rat> = f.iter().rev().cloned().collect();
    if asc.len() < 2 {
        return true;
    }
    dense::gcd_field(&asc, &dense::derivative(&asc)).len() <= 1
}

// ---------------------------------------------------------------------------
// Real isolating intervals with their Mobius state.

/// A real isolating interval that can be refined: `f` transformed by the
/// Mobius map, possibly of the mirrored polynomial (`neg`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RealBox {
    pub f: ZPoly,
    pub m: Mobius,
    pub neg: bool,
}

impl RealBox {
    /// From an isolating interval `(s, t)` of a root of `f` over ZZ.
    pub fn from_interval(f: &[Int], s: &Rat, t: &Rat) -> IsoResult<RealBox> {
        let (mut s, mut t, mut f, mut neg) = (s.clone(), t.clone(), f.to_vec(), false);
        if s.is_negative() {
            if !t.is_positive_or_zero_rat() {
                f = zmirror(&f);
                (s, t) = (t.neg(), s.neg());
                neg = true;
            } else {
                return Err(IsoError::Value(format!(
                    "Cannot refine a real root in ({}, {})",
                    rat_str(&s),
                    rat_str(&t)
                )));
            }
        }
        let m = mobius_from_interval(&s, &t);
        let f = ztransform(&f, &zstrip(vec![m.0.clone(), m.1.clone()]), &zstrip(vec![m.2.clone(), m.3.clone()]));
        Ok(RealBox { f, m, neg })
    }
    fn left(&self, neg: bool) -> IsoResult<Rat> {
        let (a, b, c, d) = &self.m;
        let ad_lt = a.mul(d) < b.mul(c);
        let ad_gt = a.mul(d) > b.mul(c);
        Ok(if !neg {
            if ad_lt {
                frac(a, c)?
            } else {
                frac(b, d)?
            }
        } else if ad_gt {
            frac(a, c)?.neg()
        } else {
            frac(b, d)?.neg()
        })
    }
    pub fn a(&self) -> IsoResult<Rat> {
        self.left(self.neg)
    }
    pub fn b(&self) -> IsoResult<Rat> {
        Ok(self.left(!self.neg)?.neg())
    }
    pub fn refine(&self) -> IsoResult<RealBox> {
        let (f, m) = inner_refine(self.f.clone(), self.m.clone(), None, Some(1), None, false)?;
        Ok(RealBox { f, m, neg: self.neg })
    }
    pub fn refine_size(&self, dx: &Rat) -> IsoResult<RealBox> {
        let mut r = self.clone();
        while !(&r.b()?.sub(&r.a()?) < dx) {
            r = r.refine()?;
        }
        Ok(r)
    }
}

/// The reference's `dup_isolate_real_roots_sqf(..., blackbox=True)`.
pub fn isolate_real_boxes_sqf(f: &[Int]) -> IsoResult<Vec<RealBox>> {
    let (intervals, g) = isolate_real_roots_sqf(f, None, None, None, false)?;
    intervals.iter().map(|(s, t)| RealBox::from_interval(&g, s, t)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z(v: &[i64]) -> ZPoly {
        v.iter().map(|&c| Int::from(c)).collect()
    }
    fn q(n: i64, d: i64) -> Rat {
        Rat::new(Int::from(n), Int::from(d))
    }

    /// The reference's count: QQ Sturm sequence of the monic square-free
    /// part and rational evaluation (the previous implementation).
    fn qq_count(f: &[Rat], inf: Option<&Rat>, sup: Option<&Rat>) -> usize {
        let chain: Vec<Vec<Rat>> = sturm_sequence(f).into_iter().map(|s| s.into_iter().rev().collect()).collect();
        let sign = |c: &Rat| if c.is_zero() { 0 } else if c.is_negative() { -1 } else { 1 };
        let var = |v: Vec<i32>| {
            let (mut prev, mut k) = (0, 0i64);
            for s in v {
                if s * prev < 0 {
                    k += 1;
                }
                if s != 0 {
                    prev = s;
                }
            }
            k
        };
        let at = |x: Option<&Rat>, minus: bool| -> Vec<i32> {
            chain
                .iter()
                .map(|s| match x {
                    Some(x) => sign(&super::super::dense::eval(s, x)),
                    None => {
                        let l = sign(s.last().unwrap());
                        if minus && (s.len() - 1) % 2 == 1 { -l } else { l }
                    }
                })
                .collect()
        };
        let mut n = (var(at(inf, true)) - var(at(sup, false))).unsigned_abs() as usize;
        let asc: Vec<Rat> = f.iter().rev().cloned().collect();
        if inf.is_some_and(|x| super::super::dense::eval(&asc, x).is_zero()) {
            n += 1;
        }
        n
    }

    #[test]
    fn integer_sturm_count_matches_rational_sturm_count() {
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let points: Vec<Option<Rat>> = vec![None, Some(q(0, 1)), Some(q(1, 1)), Some(q(-1, 1)), Some(q(-3, 2)), Some(q(7, 5)), Some(q(2, 1)), Some(q(-1, 3))];
        for case in 0..60 {
            let deg = 1 + (next() % 9) as usize;
            let mut f: Vec<Rat> = (0..=deg)
                .map(|_| {
                    let v = next();
                    if v % 3 == 0 { q(0, 1) } else { q((v % 11) as i64 - 5, 1 + (v >> 8) as i64 % 4) }
                })
                .collect();
            if f[0].is_zero() {
                f[0] = q(1 + case % 3, 1);
            }
            // repeated factors: multiply by (x - 1)**2 sometimes
            if case % 4 == 0 {
                let asc: Vec<Rat> = f.iter().rev().cloned().collect();
                let sq = super::super::dense::mul(&asc, &[q(1, 1), q(-2, 1), q(1, 1)]);
                f = sq.into_iter().rev().collect();
            }
            for a in &points {
                for b in &points {
                    assert_eq!(count_real_roots(&f, a.as_ref(), b.as_ref()), qq_count(&f, a.as_ref(), b.as_ref()), "{f:?} {a:?} {b:?}");
                }
            }
        }
        let mut trinomial = vec![q(0, 1); 302];
        trinomial[0] = q(1, 1);
        trinomial[300] = q(-1, 1);
        trinomial[301] = q(-1, 1);
        for (a, b) in [(None, None), (Some(q(0, 1)), Some(q(2, 1))), (Some(q(-3, 2)), Some(q(-1, 2)))] {
            assert_eq!(count_real_roots(&trinomial, a.as_ref(), b.as_ref()), qq_count(&trinomial, a.as_ref(), b.as_ref()));
        }
    }

    #[test]
    fn lazy_pseudo_remainder_matches_classical() {
        let a = z(&[3, 0, -1, 0, 0, 7, 0, 2, -5]);
        for b in [z(&[2, 5]), z(&[-3, 0, 0, 1]), z(&[4, -1, 0, 6, 0]), z(&[-7, 1, 1])] {
            let mut b = b;
            super::super::dense::trim(&mut b);
            let a_asc: Vec<Int> = a.iter().rev().cloned().collect();
            let b_asc: Vec<Int> = b.iter().rev().cloned().collect();
            assert_eq!(lazy_prem(&a_asc, &b_asc), super::super::dense::prem(&a_asc, &b_asc));
        }
    }

    #[test]
    fn python_logarithm_matches_float_semantics() {
        assert_eq!(py_log2(&Int::from(8)), 3);
        assert_eq!(py_log2(&Int::from(9)), 3);
        assert_eq!(py_log2(&Int::from(1)), 0);
        // math.log(2**29, 2) == 29.000000000000004
        assert_eq!(py_log2(&Int::from(1 << 29)), 29);
        let big = Int::from_big(num_bigint::BigInt::from(1u8) << 2000usize);
        assert_eq!(py_log2(&big), 2000);
    }

    #[test]
    fn continued_fraction_intervals_match_the_reference() {
        // intervals(x**3 - 2) == [((1, 2), 1)]
        let r = isolate_real_roots(&z(&[1, 0, 0, -2]), None, None, None, false).unwrap();
        assert_eq!(r, vec![((q(1, 1), q(2, 1)), 1)]);
        // intervals(x**4 - 10*x**2 + 1, eps=1/100)
        let r = isolate_real_roots(&z(&[1, 0, -10, 0, 1]), Some(&q(1, 100)), None, None, false).unwrap();
        let got: Vec<_> = r.into_iter().map(|(i, _)| i).collect();
        assert_eq!(
            got,
            vec![(q(-107, 34), q(-22, 7)), (q(-7, 22), q(-6, 19)), (q(6, 19), q(7, 22)), (q(22, 7), q(107, 34))]
        );
        // (x**2-2)*(x**2-3)**2
        let f = z(&[1, 0, -8, 0, 21, 0, -18]);
        let r = isolate_real_roots(&f, None, None, None, false).unwrap();
        assert_eq!(
            r,
            vec![
                ((q(-2, 1), q(-3, 2)), 2),
                ((q(-3, 2), q(-1, 1)), 1),
                ((q(1, 1), q(3, 2)), 1),
                ((q(3, 2), q(2, 1)), 2)
            ]
        );
    }

    #[test]
    fn refinement_matches_the_reference() {
        let f = z(&[1, 0, -2]);
        assert_eq!(refine_real_root(&f, &q(1, 1), &q(2, 1), Some(&q(1, 100)), None, None, false).unwrap(), (q(24, 17), q(17, 12)));
        assert_eq!(refine_real_root(&f, &q(-2, 1), &q(-1, 1), None, Some(3), None, false).unwrap(), (q(-3, 2), q(-7, 5)));
        assert_eq!(refine_real_root(&f, &q(1, 1), &q(2, 1), None, Some(5), None, true).unwrap(), (q(7, 5), q(17, 12)));
        assert!(matches!(refine_real_root(&f, &q(2, 1), &q(3, 1), None, Some(1), None, false), Err(IsoError::Refinement(_))));
    }

    #[test]
    fn complex_rectangles_match_the_reference() {
        // x**5 - x - 1: ((-2 - 2*I, 0), (-2, 2*I), (-2*I, 2), (0, 2 + 2*I))
        let f: QPoly = [1, 0, 0, 0, -1, -1].iter().map(|&c| Rat::from_i64(c)).collect();
        let boxes = isolate_complex_roots_sqf(&f, None, None, None).unwrap();
        let got: Vec<_> = boxes.iter().map(|b| (b.ax(), b.ay(), b.bx(), b.by())).collect();
        let w = |a, b, c, d| (q(a, 1), q(b, 1), q(c, 1), q(d, 1));
        assert_eq!(got, vec![w(-2, -2, 0, 0), w(-2, 0, 0, 2), w(0, -2, 2, 0), w(0, 0, 2, 2)]);
        // count_roots(x**4 + 1, -2-2I, 2+2I) == 4 and in [0, 2+2I] == 1
        let g: QPoly = [1, 0, 0, 0, 1].iter().map(|&c| Rat::from_i64(c)).collect();
        let c = |inf: (i64, i64), sup: (i64, i64)| {
            count_complex_roots(&g, Some((q(inf.0, 1), q(inf.1, 1))), Some((q(sup.0, 1), q(sup.1, 1))), &Exclude::None).unwrap()
        };
        assert_eq!(c((-2, -2), (2, 2)), 4);
        assert_eq!(c((0, 0), (2, 2)), 1);
    }

    #[test]
    fn real_boxes_refine_like_the_reference() {
        // CRootOf(x**3 - x**2 - x - 1, 0)._get_interval() == (1, 2) ...
        let boxes = isolate_real_boxes_sqf(&z(&[1, -1, -1, -1])).unwrap();
        assert_eq!(boxes.len(), 1);
        assert_eq!((boxes[0].a().unwrap(), boxes[0].b().unwrap()), (q(1, 1), q(2, 1)));
        let r = boxes[0].refine_size(&q(1, 10)).unwrap();
        assert!(r.a().unwrap() < r.b().unwrap());
        let neg = isolate_real_boxes_sqf(&z(&[1, 0, -2])).unwrap();
        assert_eq!((neg[0].a().unwrap(), neg[0].b().unwrap()), (q(-2, 1), q(-1, 1)));
    }
}
