//! Minimal polynomials of algebraic expressions by composition.
//!
//! This follows the reference `_minpoly_compose`: rationals, `I`, the
//! golden ratio and the tribonacci constant have fixed polynomials; sums
//! and products combine minimal polynomials through resultants; rational
//! powers use `resultant(mp(y), x**d - y**n)`; `sin`, `cos`, `tan` of
//! rational multiples of `pi`, `exp(I*pi*r)` and indexed roots are
//! factors of explicit polynomials. After every resultant the irreducible
//! factor vanishing at the value is selected.
//!
//! The selection is certified: each candidate is evaluated on a rigorous
//! enclosure of the value and discarded only when the enclosure of its
//! value excludes zero, with the working precision doubled until a single
//! candidate remains. When the expression has free symbols they are
//! coefficients of the ground field `QQ(syms)`; candidates are then
//! compared at the integer points tried by the reference, in the same
//! order, and the first point that separates them decides.
//!
//! The product of rational radicals `c*prod(b_i**e_i)` is handled as one
//! root of `x**L - R` whose factor is selected like any other, so the
//! result is always irreducible. The reference returns `x**L - R` itself
//! there, which is reducible for values like `(-8)**(1/3)`.

use rustc_hash::FxHashMap;

use super::dense;
use super::enclosure::{self, EncError};
use super::factor_multi::factor_list_mz;
use super::factor_uni::factor_zx;
use super::int::{Int, Rat};
use super::mpoly::MPoly;
use super::resultant::resultant_mz;
use crate::expr::{Constant, Expr, Func, Kind};

/// Why no minimal polynomial was produced.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MinpolyError {
    /// `NotAlgebraic("%s does not seem to be an algebraic element")`.
    NotAlgebraic(Expr),
    /// `ZeroDivisionError("%s is zero")`.
    ZeroDivision(Expr),
    /// `NotImplementedError("multiple candidates for the minimal
    /// polynomial of %s")`: no precision or sample point separated the
    /// candidate factors.
    Multiple(Expr),
}

impl std::fmt::Display for MinpolyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MinpolyError::NotAlgebraic(e) => write!(f, "{e} does not seem to be an algebraic element"),
            MinpolyError::ZeroDivision(e) => write!(f, "{e} is zero"),
            MinpolyError::Multiple(e) => write!(f, "multiple candidates for the minimal polynomial of {e}"),
        }
    }
}

type MResult<T> = Result<T, MinpolyError>;

/// Highest working precision for constant values. Distinct irreducible
/// factors never share a root, so this bound is only reached for values
/// whose enclosure cannot be computed (a radical of a number lying on the
/// branch cut that is not recognised as exactly real).
const MAX_PREC: u32 = 1 << 15;
/// Highest working precision at one sample point of the symbols.
const POINT_PREC: u32 = 1 << 10;

/// Composition context: variable 0 is `x`, variables `1..=k` are the
/// symbols of the ground field.
struct Ctx<'a> {
    syms: &'a [Expr],
}

impl Ctx<'_> {
    fn n(&self) -> usize {
        1 + self.syms.len()
    }
    fn x(&self) -> MPoly<Int> {
        MPoly::var(0, self.n())
    }
    fn constant(&self, c: i64) -> MPoly<Int> {
        MPoly::constant(Int::from(c), self.n())
    }
    /// `q*x - p`.
    fn rational(&self, r: &Rat) -> MPoly<Int> {
        self.x().scale(&r.d).sub(&MPoly::constant(r.n.clone(), self.n()))
    }
    fn dense(&self, f: &[Int]) -> MPoly<Int> {
        dense::to_mpoly(f, 0, self.n())
    }
}

/// `p(x)` re-embedded in `n + 1` variables, `x` moved to the new last
/// variable `y`.
fn to_y(p: &MPoly<Int>) -> MPoly<Int> {
    let n = p.nvars;
    let mut perm: Vec<usize> = (0..n).collect();
    perm[0] = n;
    p.permute(&perm, n + 1)
}

/// `p` re-embedded in `n + 1` variables.
fn lift(p: &MPoly<Int>) -> MPoly<Int> {
    let n = p.nvars;
    let perm: Vec<usize> = (0..n).collect();
    p.permute(&perm, n + 1)
}

/// Drops the last variable of a polynomial free of it.
fn drop_last(p: &MPoly<Int>) -> MPoly<Int> {
    let n = p.nvars - 1;
    let perm: Vec<usize> = (0..=n).map(|v| v.min(n)).collect();
    p.permute(&perm, n)
}

/// `p(x - y)` in `n + 1` variables.
fn shifted(p: &MPoly<Int>) -> MPoly<Int> {
    let n = p.nvars + 1;
    let lin = MPoly::var(0, n).sub(&MPoly::var(n - 1, n));
    horner_in_x(p, &lin)
}

/// `sum c_i(syms) * t**i` for `p = sum c_i x**i`, in `n + 1` variables.
fn horner_in_x(p: &MPoly<Int>, t: &MPoly<Int>) -> MPoly<Int> {
    let n = t.nvars;
    let cs = lift(p).coeffs_in(0);
    let top = cs.first().map(|c| c.0).unwrap_or(0);
    let mut acc = MPoly::zero(n);
    let mut idx = 0;
    for k in (0..=top).rev() {
        acc = acc.mul(t);
        if idx < cs.len() && cs[idx].0 == k {
            acc = acc.add(&cs[idx].1);
            idx += 1;
        }
    }
    acc
}

/// `y**deg * p(x/y)` in `n + 1` variables.
fn homogenized(p: &MPoly<Int>) -> MPoly<Int> {
    let n = p.nvars + 1;
    let d = p.degree(0).max(0) as u32;
    let mut terms = Vec::with_capacity(p.len());
    for (e, c) in p.iter() {
        let mut ne = e.to_vec();
        ne.push(d - e[0]);
        terms.push((ne, c.clone()));
    }
    MPoly::from_terms(n, terms)
}

/// `x**deg * p(1/x)`.
fn inverted(p: &MPoly<Int>) -> MPoly<Int> {
    let d = p.degree(0).max(0) as u32;
    let terms = p
        .iter()
        .map(|(e, c)| {
            let mut ne = e.to_vec();
            ne[0] = d - e[0];
            (ne, c.clone())
        })
        .collect();
    MPoly::from_terms(p.nvars, terms)
}

/// Irreducible factors involving `x`.
fn factors(p: &MPoly<Int>) -> Vec<MPoly<Int>> {
    if p.nvars == 1 {
        let (_, fs) = factor_zx(&dense::from_mpoly(p, 0));
        return fs.into_iter().map(|(f, _)| dense::to_mpoly(&f, 0, 1)).collect();
    }
    let (_, fs) = factor_list_mz(p);
    fs.into_iter().map(|(f, _)| f).filter(|f| f.degree(0) > 0).collect()
}

fn rational_of(e: &Expr) -> Option<Rat> {
    e.as_number().and_then(Rat::from_number)
}

impl Ctx<'_> {
    /// The factor of `p` vanishing at `v`.
    fn choose_of(&self, p: &MPoly<Int>, v: &Expr) -> MResult<MPoly<Int>> {
        self.choose(factors(p), v)
    }

    fn choose(&self, mut candidates: Vec<MPoly<Int>>, v: &Expr) -> MResult<MPoly<Int>> {
        if candidates.len() == 1 {
            return Ok(candidates.pop().unwrap());
        }
        if candidates.is_empty() {
            return Err(MinpolyError::Multiple(v.clone()));
        }
        if self.syms.is_empty() {
            let dense: Vec<Vec<Rat>> = candidates.iter().map(|f| dense::to_rat(&dense::from_mpoly(f, 0))).collect();
            return match select(&dense, v, MAX_PREC) {
                Some(i) => Ok(candidates.swap_remove(i)),
                None => Err(MinpolyError::Multiple(v.clone())),
            };
        }
        // Points `subsets(range(5), k, repetition=True)` in the reference
        // order, assigned to the symbols in order.
        let k = self.syms.len();
        let mut point = vec![0u32; k];
        loop {
            let mut map = FxHashMap::default();
            for (s, &i) in self.syms.iter().zip(&point) {
                map.insert(s.clone(), Expr::int(i as i64));
            }
            let value = crate::subs::xreplace(v, &map);
            let values: Vec<Int> = std::iter::once(Int::zero()).chain(point.iter().map(|&i| Int::from(i as i64))).collect();
            let dense: Vec<Vec<Rat>> = candidates
                .iter()
                .map(|f| {
                    let mut g = f.clone();
                    for var in (1..=k).rev() {
                        g = g.eval(var, &values[var]);
                    }
                    dense::to_rat(&dense::from_mpoly(&g, 0))
                })
                .collect();
            if let Some(i) = select(&dense, &value, POINT_PREC) {
                return Ok(candidates.swap_remove(i));
            }
            // next multiset in lexicographic order
            let mut j = k;
            loop {
                if j == 0 {
                    return Err(MinpolyError::Multiple(v.clone()));
                }
                j -= 1;
                if point[j] < 4 {
                    let next = point[j] + 1;
                    for p in &mut point[j..] {
                        *p = next;
                    }
                    break;
                }
            }
        }
    }
}

/// Index of the only candidate (lowest degree first) not excluded at
/// `v`, if one is isolated before `max_prec`.
fn select(candidates: &[Vec<Rat>], v: &Expr, max_prec: u32) -> Option<usize> {
    let mut alive: Vec<usize> = (0..candidates.len()).collect();
    let mut prec = 64u32;
    loop {
        match enclosure::enclose(v, prec) {
            Ok(b) => {
                alive.retain(|&i| {
                    let f = &candidates[i];
                    f.iter().all(|c| c.is_zero()) || !enclosure::horner_rat(f, &b, prec).excludes_zero()
                });
                match alive.len() {
                    0 => return None,
                    1 => return Some(alive[0]),
                    _ => {}
                }
            }
            Err(EncError::NotAlgebraic(_)) => return None,
            Err(EncError::Precision) => {}
        }
        if prec >= max_prec {
            return None;
        }
        prec = (prec * 2).min(max_prec);
    }
}

/// Chebyshev polynomial of the first kind, lowest degree first.
fn chebyshev_t(n: u64) -> Vec<Int> {
    let mut a = vec![Int::one()];
    if n == 0 {
        return a;
    }
    let mut b = vec![Int::zero(), Int::one()];
    for _ in 1..n {
        let mut c = vec![Int::zero(); b.len() + 1];
        for (i, v) in b.iter().enumerate() {
            c[i + 1] = c[i + 1].add(&v.mul_i64(2));
        }
        for (i, v) in a.iter().enumerate() {
            c[i] = c[i].sub(v);
        }
        a = std::mem::replace(&mut b, c);
    }
    b
}

/// The rational `c` with `e == c*a` for the given atom, if any.
fn coeff_of_atom(e: &Expr, atom: &dyn Fn(&Expr) -> bool) -> Option<Rat> {
    if atom(e) {
        return Some(Rat::one());
    }
    if let Kind::Mul(c, factors) = e.kind() {
        let rest = crate::ops::mul(&factors.iter().map(|(b, x)| crate::ops::mul_factor(b, x)).collect::<Vec<_>>());
        if atom(&rest) {
            return Rat::from_number(c);
        }
    }
    None
}

fn is_pi(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Const(Constant::Pi))
}

fn is_i_pi(e: &Expr) -> bool {
    if let Kind::Mul(c, factors) = e.kind() {
        return c.is_one()
            && factors.len() == 2
            && factors.iter().all(|(_, x)| x.is_one())
            && factors.iter().any(|(b, _)| matches!(b.kind(), Kind::Const(Constant::I)))
            && factors.iter().any(|(b, _)| is_pi(b));
    }
    false
}

impl Ctx<'_> {
    fn compose(&self, ex: &Expr) -> MResult<MPoly<Int>> {
        let not = || MinpolyError::NotAlgebraic(ex.clone());
        match ex.kind() {
            Kind::Num(n) => Rat::from_number(n).map(|r| self.rational(&r)).ok_or_else(not),
            Kind::Const(Constant::I) => {
                // x**2 + 1 stays irreducible over QQ(syms).
                Ok(self.x().pow(2).add(&self.constant(1)))
            }
            Kind::Const(Constant::GoldenRatio) => Ok(self.dense(&[Int::from(-1), Int::from(-1), Int::one()])),
            Kind::Const(Constant::TribonacciConstant) => Ok(self.dense(&[Int::from(-1), Int::from(-1), Int::from(-1), Int::one()])),
            Kind::Sym(_) => match self.syms.iter().position(|s| s == ex) {
                Some(i) => Ok(self.x().sub(&MPoly::var(i + 1, self.n()))),
                None => Err(not()),
            },
            Kind::Add(..) => {
                let args = crate::traverse::ordered_args(ex);
                let mut mp = self.compose(&args[0])?;
                let mut p = args[0].clone();
                for px in &args[1..] {
                    let mp2 = self.compose(px)?;
                    let value = crate::ops::add2(&p, px);
                    mp = self.op_add(&mp, &mp2, &value)?;
                    p = value;
                }
                Ok(mp)
            }
            Kind::Mul(..) => self.compose_mul(ex),
            Kind::Pow(b, x) => self.compose_pow(ex, b, x),
            Kind::Func(f @ (Func::Sin | Func::Cos | Func::Tan), args) if args.len() == 1 => {
                let c = coeff_of_atom(&args[0], &is_pi).ok_or_else(not)?;
                let f = match f {
                    Func::Sin => self.cos_pi(&Rat::new(Int::one(), Int::from(2)).sub(&c)),
                    Func::Cos => self.cos_pi(&c),
                    _ => self.tan_pi(&c),
                };
                self.choose_of(&f, ex)
            }
            Kind::Func(Func::Exp, args) if args.len() == 1 => {
                let c = coeff_of_atom(&args[0], &is_i_pi).ok_or_else(not)?;
                // exp(I*pi*p/q) is a root of x**(2q) - 1.
                let q = c.d.as_i64().filter(|&q| q > 0 && q <= 1 << 20).ok_or_else(not)? as u32;
                let f = self.x().pow(2 * q).sub(&self.constant(1));
                self.choose_of(&f, ex)
            }
            Kind::Func(Func::Undefined(_), _) if super::rootof::is_rootof(ex) => {
                let (expr, _) = super::rootof::parts(ex).ok_or_else(not)?;
                let (_, mut f) = super::rootof::integer_poly(&expr).ok_or_else(not)?;
                f.reverse();
                self.choose_of(&self.dense(&f), ex)
            }
            _ => Err(not()),
        }
    }

    /// `T_q(x) - (-1)**p` for `cos(pi*p/q)`.
    fn cos_pi(&self, c: &Rat) -> MPoly<Int> {
        let q = c.d.as_i64().unwrap_or(1).max(1) as u64;
        let mut t = chebyshev_t(q);
        let sign = if c.n.rem_u64(2) == 0 { 1 } else { -1 };
        t[0] = t[0].sub(&Int::from(sign));
        self.dense(&t)
    }

    /// The reference polynomial with root `tan(pi*c)`.
    fn tan_pi(&self, c: &Rat) -> MPoly<Int> {
        let c = c.mul_int(&Int::from(2));
        let n = c.d.clone();
        let odd = c.n.rem_u64(2) == 1;
        let mut a = if odd { Int::one() } else { n.clone() };
        let mut f = vec![Int::zero(); n.as_i64().unwrap_or(1) as usize + 1];
        let mut k = if odd { 0i64 } else { 1 };
        let top = n.as_i64().unwrap_or(1);
        while k <= top {
            f[k as usize] = a.clone();
            // a = -(a*(n-k-1)*(n-k)) // ((k+1)*(k+2))
            let num = a.mul(&n.sub(&Int::from(k + 1))).mul(&n.sub(&Int::from(k)));
            a = num.neg().div_mod_floor(&Int::from((k + 1) * (k + 2))).0;
            k += 2;
        }
        self.dense(&f)
    }

    fn compose_mul(&self, ex: &Expr) -> MResult<MPoly<Int>> {
        let Kind::Mul(c, factors) = ex.kind() else { unreachable!() };
        if !self.syms.is_empty() {
            // Over QQ(syms) the reference folds over the arguments in order.
            let args = crate::traverse::ordered_args(ex);
            let mut mp = self.compose(&args[0])?;
            let mut p = args[0].clone();
            for px in &args[1..] {
                let mp2 = self.compose(px)?;
                let value = crate::ops::mul2(&p, px);
                mp = self.op_mul(&mp, &mp2, &value)?;
                p = value;
            }
            return Ok(mp);
        }
        // Split rational bases with rational exponents from the rest.
        let mut radical: Vec<(Rat, Rat)> = Vec::new();
        let mut rest: Vec<Expr> = Vec::new();
        for (b, x) in factors {
            match (rational_of(b), rational_of(x)) {
                (Some(rb), Some(rx)) => radical.push((rb, rx)),
                _ => rest.push(crate::ops::mul_factor(b, x)),
            }
        }
        let c = Rat::from_number(c).ok_or_else(|| MinpolyError::NotAlgebraic(ex.clone()))?;
        let has_radical = !radical.is_empty() || !c.is_one();
        let mp1 = match rest.len() {
            0 => None,
            1 => Some((self.compose(&rest[0])?, rest[0].clone())),
            _ => {
                let mut mp = self.compose(&rest[0])?;
                let mut p = rest[0].clone();
                for px in &rest[1..] {
                    let mp2 = self.compose(px)?;
                    let value = crate::ops::mul2(&p, px);
                    mp = self.op_mul(&mp, &mp2, &value)?;
                    p = value;
                }
                Some((mp, p))
            }
        };
        if !has_radical {
            return Ok(mp1.expect("a product has factors").0);
        }
        let mut w_factors = vec![Expr::num(c.to_number())];
        w_factors.extend(radical.iter().map(|(b, e)| crate::ops::pow(&Expr::num(b.to_number()), &Expr::num(e.to_number()))));
        let w = crate::ops::mul(&w_factors);
        let mp2 = self.radical(&c, &radical, &w, ex)?;
        match mp1 {
            None => Ok(mp2),
            Some((mp1, p)) => self.op_mul(&mp1, &mp2, &crate::ops::mul2(&p, &w)),
        }
    }

    /// Minimal polynomial of `w = c*prod(b_i**e_i)`, a root of
    /// `x**L - w**L` with `L` the common denominator of the exponents.
    fn radical(&self, c: &Rat, radical: &[(Rat, Rat)], w: &Expr, ex: &Expr) -> MResult<MPoly<Int>> {
        let too_big = || MinpolyError::NotAlgebraic(ex.clone());
        let lcm = radical.iter().fold(Int::one(), |l, (_, e)| l.lcm(&e.d));
        let l = lcm.as_i64().filter(|&l| l > 0 && l <= 1 << 16).ok_or_else(too_big)? as u32;
        let mut power = c.pow(l);
        for (b, e) in radical {
            let k = e.mul_int(&lcm).n.as_i64().filter(|k| k.unsigned_abs() <= 1 << 20).ok_or_else(too_big)?;
            let base = if k < 0 { b.inv() } else { b.clone() };
            power = power.mul(&base.pow(k.unsigned_abs() as u32));
        }
        let f = self.x().pow(l).scale(&power.d).sub(&MPoly::constant(power.n.clone(), self.n()));
        if l == 1 {
            return Ok(f);
        }
        self.choose_of(&f, w)
    }

    fn compose_pow(&self, ex: &Expr, b: &Expr, x: &Expr) -> MResult<MPoly<Int>> {
        if let (Some(rb), Some(rx)) = (rational_of(b), rational_of(x)) {
            return self.radical(&Rat::one(), &[(rb, rx)], ex, ex);
        }
        let mut mp = self.compose(b)?;
        let r = rational_of(x).ok_or_else(|| MinpolyError::NotAlgebraic(b.clone()))?;
        if r.is_negative() {
            if mp == self.x() {
                return Err(MinpolyError::ZeroDivision(b.clone()));
            }
            mp = inverted(&mp);
            if r.n.is_minus_one() && r.d.is_one() {
                return Ok(mp);
            }
        }
        let (n, d) = (r.n.abs(), r.d.clone());
        let n = n.as_i64().filter(|&k| k <= 1 << 16).ok_or_else(|| MinpolyError::NotAlgebraic(b.clone()))? as u32;
        let d = d.as_i64().filter(|&k| k <= 1 << 16).ok_or_else(|| MinpolyError::NotAlgebraic(b.clone()))? as u32;
        let m = self.n() + 1;
        let y = MPoly::var(m - 1, m);
        let g = MPoly::var(0, m).pow(d).sub(&y.pow(n));
        let res = drop_last(&resultant_mz(&to_y(&mp), &g, m - 1));
        self.choose_of(&res, ex)
    }

    fn op_add(&self, mp1: &MPoly<Int>, mp2: &MPoly<Int>, value: &Expr) -> MResult<MPoly<Int>> {
        let m = self.n() + 1;
        let r = drop_last(&resultant_mz(&shifted(mp1), &to_y(mp2), m - 1));
        if mp2.degree(0) == 1 {
            return Ok(r);
        }
        self.choose_of(&r, value)
    }

    fn op_mul(&self, mp1: &MPoly<Int>, mp2: &MPoly<Int>, value: &Expr) -> MResult<MPoly<Int>> {
        let m = self.n() + 1;
        let r = drop_last(&resultant_mz(&homogenized(mp1), &to_y(mp2), m - 1));
        if mp1.degree(0) == 1 || mp2.degree(0) == 1 {
            return Ok(r);
        }
        self.choose_of(&r, value)
    }
}

/// `_mexpand(ex, recursive=True)`.
pub fn mexpand_fixed(ex: &Expr) -> Expr {
    let mut e = ex.clone();
    for _ in 0..16 {
        let next = crate::simptools::mexpand(&e);
        if next == e {
            break;
        }
        e = next;
    }
    e
}

/// Minimal polynomial of `ex` over `QQ(syms)` as a polynomial in the
/// variables `[x, syms...]`, primitive over the integers with the
/// reference sign normalisation (a leading coefficient that is negative
/// is negated).
pub fn minimal_polynomial(ex: &Expr, syms: &[Expr]) -> MResult<MPoly<Int>> {
    let ctx = Ctx { syms };
    let e = mexpand_fixed(ex);
    let p = ctx.compose(&e)?;
    let (_, mut p) = p.primitive();
    let lc = p.lc_in(0);
    let negative = match lc.constant_value() {
        Some(c) => c.is_negative(),
        None => crate::assumptions::is_negative(&mpoly_expr(&lc, &Expr::int(0), syms)) == Some(true),
    };
    if negative {
        p = p.neg();
    }
    Ok(p)
}

/// The expression of `p` in `[x, syms...]`, collected in `x` like the
/// reference `result.collect(x)`.
pub fn mpoly_expr(p: &MPoly<Int>, x: &Expr, syms: &[Expr]) -> Expr {
    let mut terms = Vec::new();
    for (k, c) in p.coeffs_in(0) {
        let mut parts = Vec::new();
        for (e, v) in c.iter() {
            let mut factors = vec![Expr::num(v.to_number())];
            for (s, &d) in syms.iter().zip(&e[1..]) {
                if d > 0 {
                    factors.push(crate::ops::powi(s, d as i64));
                }
            }
            parts.push(crate::ops::mul(&factors));
        }
        let coefficient = crate::ops::add(&parts);
        terms.push(crate::ops::mul(&[coefficient, crate::ops::powi(x, k as i64)]));
    }
    crate::ops::add(&terms)
}

/// Minimal polynomial over `QQ` as integer coefficients, lowest degree
/// first.
pub fn minimal_polynomial_qq(ex: &Expr) -> MResult<Vec<Int>> {
    let p = minimal_polynomial(ex, &[])?;
    Ok(dense::from_mpoly(&p, 0))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn expr(s: &str) -> Expr {
        crate::parser::parse_str(s).unwrap()
    }
    fn mp(s: &str) -> Vec<i64> {
        minimal_polynomial_qq(&expr(s)).unwrap_or_else(|e| panic!("{s}: {e:?}")).iter().map(|c| c.as_i64().unwrap()).collect()
    }

    #[test]
    fn radicals_and_sums() {
        assert_eq!(mp("2^(1/4)+sqrt(2)"), vec![2, -8, -4, 0, 1]);
        assert_eq!(mp("sqrt(2)+sqrt(3)+sqrt(5)"), vec![576, 0, -960, 0, 352, 0, -40, 0, 1]);
        assert_eq!(mp("(-8)^(1/3)"), vec![4, -2, 1]);
        assert_eq!(mp("(-1)^(2/3)"), vec![1, 1, 1]);
        assert_eq!(mp("1/(2^(1/3)+1)"), vec![-1, 3, -3, 3]);
        assert_eq!(mp("(2^(1/3)+1)^(-2)"), vec![-1, 3, 9, 9]);
        assert_eq!(mp("sqrt(2)*sqrt(3)-sqrt(6)"), vec![0, 1]);
        assert_eq!(mp("(1+I)^(1/3)"), vec![2, 0, 0, -2, 0, 0, 1]);
    }

    #[test]
    fn trigonometric_and_roots() {
        assert_eq!(mp("cos(pi/7)"), vec![1, -4, -4, 8]);
        assert_eq!(mp("sin(pi/5)"), vec![5, 0, -20, 0, 16]);
        assert_eq!(mp("tan(pi/8)"), vec![-1, 2, 1]);
        assert_eq!(mp("exp(2*I*pi/5)"), vec![1, 1, 1, 1, 1]);
        let r = super::super::rootof::complex_root(&expr("x^5-x-1"), 0).unwrap();
        let e = crate::ops::add2(&crate::ops::powi(&r, 2), &Expr::int(1));
        let p: Vec<i64> = minimal_polynomial_qq(&e).unwrap().iter().map(|c| c.as_i64().unwrap()).collect();
        assert_eq!(p, vec![-1, 0, -4, 8, -5, 1]);
    }

    #[test]
    fn errors_name_the_failing_element() {
        assert_eq!(minimal_polynomial_qq(&expr("pi+sqrt(2)")), Err(MinpolyError::NotAlgebraic(expr("pi"))));
        assert_eq!(minimal_polynomial_qq(&expr("sqrt(2)^pi")), Err(MinpolyError::NotAlgebraic(expr("2"))));
    }

    #[test]
    fn symbolic_ground_field() {
        let y = expr("y");
        let p = minimal_polynomial(&expr("sqrt(y)+1/y"), std::slice::from_ref(&y)).unwrap();
        assert_eq!(mpoly_expr(&p, &expr("x"), &[y]), expr("x**2*y**2 - 2*x*y - y**3 + 1"));
    }
}
