//! `Expr.nseries(x, 0, n, cdir=c)` for a negative or nonreal `cdir`, where
//! the reference passes the direction down to the step functions.
//!
//! SymPy 1.14 calls `_eval_nseries(x, n, logx, cdir)`. `Add` and `Mul`
//! forward `cdir` to their arguments; most functions ignore it (`exp`,
//! `Abs`, `sqrt` expand as from the right), but these use the direction
//! `ndir = arg.dir(x, cdir)` (`Expr.dir`: the first nonzero derivative
//! coefficient `a_k` at 0, `ndir = a_k*cdir**k`):
//!
//! * `floor(g)` with an integer `g(0) = r`: `r - 1` when `ndir < 0`, `r`
//!   when `ndir > 0`, `NotImplementedError` otherwise (`cdir = I`);
//! * `ceiling(g)`: `r` when `ndir < 0`, `r + 1` when `ndir > 0`;
//! * `frac(g)` with `frac(g(0)) = 0`: `g - g(0)`, plus 1 when `ndir < 0`;
//! * `sign(g)` with `g(0) = 0`: `-1` when `re(ndir) < 0`, else `1`;
//! * `atan(g)` with `g(0) = zoo`: the expansion minus `pi` when
//!   `re(ndir) > 0`;
//! * `acot(g)` with `g(0) = 0`: the expansion minus `pi` when
//!   `re(ndir) < 0`.
//!
//! * `log(g)` (`g != x`): `log(g(c*t))` expanded at `t -> 0+`, with
//!   `log(t)` read as `log(x) - log(c)` and `t = x/c` (the
//!   `- b*log(cdir)` and `res.subs(t, x/cdir)` steps of
//!   `log._eval_nseries`).
//!
//! The direction is read off the leading term `c*x**k` of `g - g(0)`
//! (`g` itself at a pole): the first nonzero derivative is of order
//! `max(1, ceil(k))` with coefficient `c*k*(k - 1)*...` (`c*k` for `k < 0`).
//!
//! `Mul` keeps the reference bookkeeping of its Order term: when the
//! leading terms of all factors exist (`sign` has none), the truncated
//! product gets `O(x**n)` unless it is exactly zero with a vanishing
//! leading term; otherwise the product of the factor expansions is exact.
//! Any other node holding a step function is expanded as before (its
//! `_eval_nseries` does not forward `cdir`). An expansion the engine
//! cannot do declines (`None`) and the caller keeps its previous path.

use crate::assumptions::{is_integer, is_negative, is_positive, is_zero};
use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add, mul, mul2, pow, sub};
use crate::series::engine::depends;
use crate::series::order::{get_o, order, remove_o};
use crate::series::{leading_term, series};
use crate::special::SpecialFn;
use crate::subs::subs1;

/// The outcome: `Ok(None)` declines, `Err` is the reference
/// `NotImplementedError` message.
pub type CdirRes = Result<Option<Expr>, String>;

enum Step {
    /// the node is replaced by this expression before the expansion
    Replace(Expr),
    /// the node expands as from the right
    Plain,
}

struct Ctx<'a> {
    x: &'a Expr,
    n: i64,
    cdir: &'a Expr,
}

fn is_step(e: &Expr) -> bool {
    matches!(
        e.kind(),
        Kind::Func(Func::Floor | Func::Ceiling | Func::Sign | Func::Atan | Func::Acot | Func::Special(SpecialFn::Frac), a)
            if a.len() == 1
    )
}

/// `log(g)` with `g` depending on `x` and not `x` itself.
fn is_branch_log(e: &Expr, x: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Log, a) if a.len() == 1 && a[0] != *x && depends(&a[0], x))
}

fn has_branch_log(e: &Expr, x: &Expr) -> bool {
    depends(e, x) && (is_branch_log(e, x) || e.args().iter().any(|a| has_branch_log(a, x)))
}

fn has_step(e: &Expr, x: &Expr) -> bool {
    if !depends(e, x) {
        return false;
    }
    if is_step(e) || is_branch_log(e, x) {
        return true;
    }
    e.args().iter().any(|a| has_step(a, x))
}

fn finite_number(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(v) if !v.is_nan())
}

/// `re(v)`, or `None` when undecided.
fn real_part(v: &Expr) -> Option<Expr> {
    match crate::complex_parts::as_real_imag(v, true, None) {
        Ok(Some((re, _))) => Some(re),
        _ => None,
    }
}

impl Ctx<'_> {
    fn plain(&self, e: &Expr) -> Option<Expr> {
        if !depends(e, self.x) {
            return Some(e.clone());
        }
        series(e, self.x, &Expr::zero(), self.n, false).ok()
    }

    /// `Expr.dir(x, cdir)` of `g` with `g(0) = g0` (`None` at a pole).
    fn ndir(&self, g: &Expr, g0: Option<&Expr>) -> Option<Expr> {
        let h = match g0 {
            Some(v) => crate::expand::expand(&sub(g, v)),
            None => g.clone(),
        };
        if is_zero(&h) == Some(true) {
            return None;
        }
        let lt = leading_term(&h, self.x, &Expr::zero(), false).ok()?;
        if depends(&lt.coeff, self.x) || depends(&lt.coeff, &lt.logt) || lt.coeff.is_zero() {
            return None;
        }
        let k = lt.exp.clone();
        let (minexp, factor) = if k.is_negative() {
            (1i64, Expr::num(k))
        } else if k.is_positive() {
            let ceil = k.neg().floor().neg();
            let m = ceil.to_i64()?.max(1);
            if m > 64 {
                return None;
            }
            let mut f = Expr::one();
            for j in 0..m {
                f = mul2(&f, &Expr::num(k.sub(&Number::int(j))));
            }
            (m, f)
        } else {
            return None;
        };
        Some(crate::expand::expand(&mul(&[lt.coeff.clone(), factor, pow(self.cdir, &Expr::int(minexp))])))
    }

    fn step(&self, node: &Expr) -> Result<Step, String> {
        let Kind::Func(f, a) = node.kind() else { return Ok(Step::Plain) };
        let g = &a[0];
        let g0 = subs1(g, self.x, &Expr::zero());
        let not_sure = |d: &Expr| Err(format!("Not sure of sign of {d}"));
        match f {
            Func::Floor | Func::Ceiling => {
                if !finite_number(&g0) || is_integer(&g0) != Some(true) {
                    return Ok(Step::Plain);
                }
                let Some(d) = self.ndir(g, Some(&g0)) else { return Ok(Step::Plain) };
                let floor = *f == Func::Floor;
                if is_negative(&d) == Some(true) {
                    Ok(Step::Replace(if floor { sub(&g0, &Expr::one()) } else { g0 }))
                } else if is_positive(&d) == Some(true) {
                    Ok(Step::Replace(if floor { g0 } else { add(&[g0, Expr::one()]) }))
                } else {
                    not_sure(&d)
                }
            }
            Func::Special(SpecialFn::Frac) => {
                if !finite_number(&g0) || is_integer(&g0) != Some(true) {
                    return Ok(Step::Plain);
                }
                let Some(d) = self.ndir(g, Some(&g0)) else { return Ok(Step::Plain) };
                let one = if is_negative(&d) == Some(true) { Expr::one() } else { Expr::zero() };
                Ok(Step::Replace(add(&[sub(g, &g0), one])))
            }
            Func::Sign => {
                if is_zero(&g0) != Some(true) {
                    return Ok(Step::Plain);
                }
                let Some(d) = self.ndir(g, Some(&g0)) else { return Ok(Step::Plain) };
                let Some(re) = real_part(&d) else { return Ok(Step::Plain) };
                match is_negative(&re) {
                    Some(true) => Ok(Step::Replace(Expr::minus_one())),
                    Some(false) => Ok(Step::Replace(Expr::one())),
                    None => Ok(Step::Plain),
                }
            }
            Func::Atan | Func::Acot => {
                let atan = *f == Func::Atan;
                let at_point = if atan { g0.is_constant(crate::expr::Constant::ComplexInfinity) } else { is_zero(&g0) == Some(true) };
                if !at_point {
                    return Ok(Step::Plain);
                }
                // the expansion from the right is the reference one only for
                // a positive leading coefficient (pi/2 - ...)
                let lt = match leading_term(g, self.x, &Expr::zero(), false) {
                    Ok(lt) => lt,
                    Err(_) => return Ok(Step::Plain),
                };
                if is_positive(&lt.coeff) != Some(true) {
                    return Ok(Step::Plain);
                }
                let Some(d) = self.ndir(g, if atan { None } else { Some(&g0) }) else { return Ok(Step::Plain) };
                let Some(re) = real_part(&d) else { return Ok(Step::Plain) };
                let shift = if atan { is_positive(&re) } else { is_negative(&re) };
                match shift {
                    Some(true) => Ok(Step::Replace(sub(node, &Expr::pi()))),
                    Some(false) => Ok(Step::Plain),
                    None => Ok(Step::Plain),
                }
            }
            _ => Ok(Step::Plain),
        }
    }

    /// The node with its replacement, or the node itself.
    fn rewritten(&self, node: &Expr) -> Result<Expr, String> {
        Ok(match self.step(node)? {
            Step::Replace(r) => r,
            Step::Plain => node.clone(),
        })
    }

    fn ns(&self, e: &Expr) -> Result<Option<Expr>, String> {
        if !has_step(e, self.x) {
            return Ok(self.plain(e));
        }
        if is_step(e) {
            let r = self.rewritten(e)?;
            return Ok(self.plain(&r));
        }
        if is_branch_log(e, self.x) {
            return Ok(self.log_series(e));
        }
        match e.kind() {
            Kind::Add(..) => {
                let mut parts = Vec::new();
                for t in e.args() {
                    match self.ns(&t)? {
                        Some(s) => parts.push(s),
                        None => return Ok(None),
                    }
                }
                Ok(Some(add(&parts)))
            }
            Kind::Mul(..) => self.ns_mul(e),
            // other nodes do not forward cdir
            _ => Ok(self.plain(e)),
        }
    }

    /// `log(g)` along `cdir` (see the module notes).
    fn log_series(&self, e: &Expr) -> Option<Expr> {
        let t = Expr::dummy("t", crate::symbol::Assumptions::POSITIVE);
        let z = subs1(&e.args()[0], self.x, &mul2(self.cdir, &t));
        let s = series(&crate::functions::log(&z), &t, &Expr::zero(), self.n, false).ok()?;
        let logt = crate::functions::log(&t);
        let logx = sub(&crate::functions::log(self.x), &crate::functions::log(self.cdir));
        let r = crate::subs::xreplace(&s, &std::iter::once((logt, logx)).collect());
        // any other logarithm of t would not be rewritten
        if has_branch_log(&r, &t) {
            return None;
        }
        let r = subs1(&r, &t, &crate::ops::div(self.x, self.cdir));
        if depends(&r, &t) {
            return None;
        }
        Some(r)
    }

    fn ns_mul(&self, e: &Expr) -> Result<Option<Expr>, String> {
        let facs = e.args();
        // the Order bookkeeping of a product with a logarithm along cdir
        // is not reproduced: keep the previous path
        if facs.iter().any(|f| has_branch_log(f, self.x)) {
            return Ok(None);
        }
        // Mul._eval_nseries first takes the leading term of every factor:
        // `sign` has none, and a sum holding a step function is not
        // reduced to a single node here.
        let try_path = facs.iter().all(|f| !has_step(f, self.x) || (is_step(f) && !matches!(f.kind(), Kind::Func(Func::Sign, _))));
        if !try_path {
            let mut out = Vec::new();
            for f in &facs {
                match self.ns(f)? {
                    Some(s) => out.push(s),
                    None => return Ok(None),
                }
            }
            return Ok(Some(crate::expand::expand(&mul(&out))));
        }
        let mut rew = Vec::new();
        let mut zero_factor = false;
        for f in &facs {
            let r = if is_step(f) && depends(f, self.x) { self.rewritten(f)? } else { f.clone() };
            if r.is_zero() {
                zero_factor = true;
            }
            rew.push(r);
        }
        let Some(s) = self.plain(&mul(&rew)) else { return Ok(None) };
        let res = remove_o(&s);
        let o = || order(&pow(self.x, &Expr::int(self.n)), &[(self.x.clone(), Expr::zero())]);
        if get_o(&s).is_some() {
            return Ok(Some(add(&[res, o()])));
        }
        // `(self - res).subs(x, 0) == 0` and a vanishing leading term keep
        // the product exact
        let at0 = sub(&subs1(e, self.x, &Expr::zero()), &subs1(&res, self.x, &Expr::zero()));
        if zero_factor && self.n > 0 && is_zero(&at0) == Some(true) {
            return Ok(Some(res));
        }
        Ok(Some(add(&[res, o()])))
    }
}

/// `expr.nseries(x, 0, n, cdir=cdir)` (from the right, `logx=None`), for
/// a negative or nonreal numeric `cdir` and an expression holding a step
/// function of `x`; `Ok(None)` otherwise.
pub fn nseries_cdir(e: &Expr, x: &Expr, n: i64, cdir: &Expr) -> CdirRes {
    if !x.is_symbol() || !has_step(e, x) || n < 0 {
        return Ok(None);
    }
    match crate::series::direction::direction(cdir) {
        Ok(Some(crate::series::direction::Direction::Real(true))) => {}
        Ok(Some(crate::series::direction::Direction::Complex(_))) => {}
        _ => return Ok(None),
    }
    let ctx = Ctx { x, n, cdir };
    ctx.ns(e)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::functions::func;

    fn x() -> Expr {
        Expr::symbol("x")
    }

    #[test]
    fn floor_from_the_left() {
        let e = func(Func::Floor, vec![x()]);
        assert_eq!(nseries_cdir(&e, &x(), 3, &Expr::minus_one()).unwrap(), Some(Expr::minus_one()));
        let e = func(Func::Ceiling, vec![x()]);
        assert_eq!(nseries_cdir(&e, &x(), 3, &Expr::minus_one()).unwrap(), Some(Expr::zero()));
    }

    #[test]
    fn floor_along_i_is_not_sure() {
        let e = func(Func::Floor, vec![x()]);
        assert!(nseries_cdir(&e, &x(), 3, &Expr::i()).is_err());
    }

    #[test]
    fn log_from_the_left() {
        let e = crate::functions::log(&crate::ops::neg(&x()));
        let want = sub(&crate::functions::log(&x()), &mul2(&Expr::i(), &Expr::pi()));
        assert_eq!(nseries_cdir(&e, &x(), 3, &Expr::minus_one()).unwrap(), Some(want));
    }

    #[test]
    fn positive_cdir_declines() {
        let e = func(Func::Floor, vec![x()]);
        assert_eq!(nseries_cdir(&e, &x(), 3, &Expr::one()).unwrap(), None);
    }
}
