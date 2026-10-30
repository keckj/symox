//! The front end of SymPy 1.14's `Expr.series` before the expansion itself:
//! the complex direction `cdir`, a non-symbol expansion variable, and the
//! `logx` replacement of `log(x)`.
//!
//! * `cdir` (reference lines after `cdir = sympify(cdir)`): a nonzero real
//!   `cdir` sets the side, `'+'` when it is positive and `'-'` otherwise;
//!   a nonreal one gives `'+'` when its imaginary part is positive. The
//!   expansion is then that of `e(x0 + u*x)` at `0+` with `u = cdir/|cdir|`,
//!   substituted back with `x -> (x - x0)/u`. For a real `u` this is the
//!   one-sided expansion from that side. At an infinite point the reference
//!   expands `e(u/x)` at `0+` and substitutes `x -> u/x`.
//! * A variable that is not a symbol (`series(exp(a*b), a*b)`) is replaced
//!   by a dummy, expanded, and substituted back (`self.subs(x, xpos)`).
//! * `logx`: `log(x)` is an independent symbol. It is replaced by a fresh
//!   symbol before the expansion and in its result, and that symbol by
//!   `logx` only once at the end, so a `logx` that holds `x` (`logx=-x`) is
//!   never expanded or substituted again (`log(-log(x))` with `logx=-x` is
//!   `log(x)`).

use crate::expr::{Expr, Func, Kind};
use crate::ops::{add, div, mul2, neg};
use crate::series::engine::depends;
use crate::series::{SErr, SRes};
use crate::subs::{subs1, xreplace};

/// The side and the unit direction of a nonzero `cdir`.
#[derive(Clone, Debug, PartialEq)]
pub enum Direction {
    /// A real direction: the one-sided expansion from the left when `true`.
    Real(bool),
    /// A nonreal unit direction `cdir/|cdir|`.
    Complex(Expr),
}

/// `None` for a zero `cdir` (the `dir` argument decides), else the
/// reference normalisation. A `cdir` with free symbols is not supported:
/// the reference would substitute a symbolic unit, which the engine cannot
/// expand with a proven side.
pub fn direction(cdir: &Expr) -> SRes<Option<Direction>> {
    use crate::assumptions::{is_positive, is_real, is_zero};
    match is_zero(cdir) {
        Some(true) => return Ok(None),
        Some(false) => {}
        None => return Err(SErr::Unsupported("series with an undecided direction cdir".into())),
    }
    if !crate::traverse::free_symbols(cdir).is_empty() {
        return Err(SErr::Unsupported("series with a symbolic direction cdir".into()));
    }
    match is_real(cdir) {
        Some(true) => match is_positive(cdir) {
            Some(p) => Ok(Some(Direction::Real(!p))),
            None => Err(SErr::Unsupported("series with a direction cdir of unknown sign".into())),
        },
        Some(false) => {
            let u = div(cdir, &crate::functions::abs(cdir));
            if is_real(&u) == Some(true) {
                return Err(SErr::Unsupported("series with an undecided direction cdir".into()));
            }
            Ok(Some(Direction::Complex(u)))
        }
        None => Err(SErr::Unsupported("series with an undecided direction cdir".into())),
    }
}

fn is_infinite_point(x0: &Expr) -> bool {
    use crate::expr::Constant;
    x0.is_constant(Constant::Infinity) || x0.is_constant(Constant::NegInfinity)
}

/// `Expr.series(x, x0, n, dir, logx, cdir)` for a fixed `n`.
pub fn series_front(e: &Expr, x: &Expr, x0: &Expr, n: i64, minus: bool, logx: Option<&Expr>, cdir: Option<&Expr>) -> SRes<Expr> {
    if !x.is_symbol() {
        let d = Expr::dummy("x", crate::symbol::Assumptions::NONE);
        let ed = subs1(e, x, &d);
        if !depends(&ed, &d) {
            return Ok(e.clone());
        }
        let r = series_front(&ed, &d, x0, n, minus, logx, cdir)?;
        return Ok(subs1(&r, &d, x));
    }
    if !depends(e, x) {
        return Ok(e.clone());
    }
    let mut minus = minus;
    let mut x0 = x0.clone();
    if let Some(c) = cdir {
        match direction(c)? {
            None => {}
            Some(Direction::Real(m)) => {
                // At an infinite point the reference expands e(u/x): the
                // point is u*oo whatever the sign of x0.
                if is_infinite_point(&x0) {
                    x0 = if m { Expr::neg_oo() } else { Expr::oo() };
                    minus = false;
                } else {
                    minus = m;
                }
            }
            Some(Direction::Complex(u)) => {
                if logx.is_some() {
                    return Err(SErr::Unsupported("series with logx and a complex direction".into()));
                }
                // The reference then expands at `0+` in a positive dummy
                // (`x` has no sign assumption): `Abs(1 + I*t)` is
                // `sqrt(t**2 + 1)` there, not `1 + I*t`.
                let t = Expr::dummy("x", crate::symbol::Assumptions::POSITIVE);
                if is_infinite_point(&x0) {
                    let s = crate::series::series(&subs1(e, x, &div(&u, &t)), &t, &Expr::zero(), n, false)?;
                    return Ok(subs1(&s, &t, &div(&u, x)));
                }
                let shifted = subs1(e, x, &add(&[x0.clone(), mul2(&u, &t)]));
                let s = crate::series::series(&shifted, &t, &Expr::zero(), n, false)?;
                let back = div(&add(&[x.clone(), neg(&x0)]), &u);
                return Ok(subs1(&s, &t, &back));
            }
        }
    }
    match logx {
        None => crate::series::series(e, x, &x0, n, minus),
        Some(l) => {
            let lx = Expr::from_kind(Kind::Func(Func::Log, vec![x.clone()]));
            let fresh = Expr::dummy("logx", crate::symbol::Assumptions::NONE);
            let to_fresh = |v: &Expr| {
                let mut m = rustc_hash::FxHashMap::default();
                m.insert(lx.clone(), fresh.clone());
                xreplace(v, &m)
            };
            let r = crate::series::series(&to_fresh(e), x, &x0, n, minus)?;
            let r = to_fresh(&r);
            let mut m = rustc_hash::FxHashMap::default();
            m.insert(fresh.clone(), l.clone());
            Ok(xreplace(&r, &m))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expr::Constant;

    fn x() -> Expr {
        Expr::symbol("x")
    }

    #[test]
    fn real_cdir_selects_the_side() {
        let e = Expr::from_kind(Kind::Func(Func::Abs, vec![x()]));
        let r = series_front(&e, &x(), &Expr::zero(), 6, false, None, Some(&Expr::minus_one())).unwrap();
        assert_eq!(r, neg(&x()));
        let r = series_front(&e, &x(), &Expr::zero(), 6, true, None, Some(&Expr::one())).unwrap();
        assert_eq!(r, x());
    }

    #[test]
    fn complex_cdir_rotates_the_variable() {
        // abs(x).series(x, 0, cdir=I) is -I*x in SymPy 1.14
        let e = Expr::from_kind(Kind::Func(Func::Abs, vec![x()]));
        let r = series_front(&e, &x(), &Expr::zero(), 6, false, None, Some(&Expr::i())).unwrap();
        assert_eq!(r, neg(&mul2(&Expr::i(), &x())));
    }

    #[test]
    fn logx_is_substituted_once() {
        // log(-log(x)).nseries(x, logx=-x) is log(x)
        let lx = Expr::from_kind(Kind::Func(Func::Log, vec![x()]));
        let e = Expr::from_kind(Kind::Func(Func::Log, vec![neg(&lx)]));
        let r = series_front(&e, &x(), &Expr::zero(), 6, false, Some(&neg(&x())), None).unwrap();
        assert_eq!(r, lx);
    }

    #[test]
    fn symbolic_cdir_declines() {
        let e = Expr::from_kind(Kind::Func(Func::Abs, vec![x()]));
        assert!(series_front(&e, &x(), &Expr::zero(), 6, false, None, Some(&Expr::symbol("c"))).is_err());
        let _ = Constant::Infinity;
    }
}
