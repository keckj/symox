//! `Abs(b**p)` with a real base of eventual sign `s` and a real exponent:
//! `(s*b)**p`. The reference gruntz runs on a positive variable, where
//! `Abs.eval` gives `base**re(exponent)` for a nonnegative base and
//! `(-base)**re(exponent)*exp(-pi*im(exponent))` for a negative one;
//! `limit(Abs((-x/(x + 1))**x), x, oo)` is `exp(-1)`. The exponent must be
//! real once the variable is positive (an assumption query, never a guess),
//! and the base must have a proven eventual sign.

use super::Gruntz;
use crate::expr::{Expr, Func, Kind};
use crate::ops::{mul2, pow};
use crate::series::engine::depends;
use crate::subs::subs1;

impl Gruntz {
    pub(super) fn abs_of_real_power(&mut self, e: &Expr) -> Option<Expr> {
        let Kind::Func(Func::Abs, args) = e.kind() else { return None };
        let [a] = args.as_slice() else { return None };
        let Kind::Pow(b, p) = a.kind() else { return None };
        if !depends(b, &self.x) {
            return None;
        }
        let xpos = Expr::dummy("x", crate::symbol::Assumptions::POSITIVE);
        if crate::assumptions::is_real(&subs1(p, &self.x, &xpos)) != Some(true)
            || crate::assumptions::is_real(&subs1(b, &self.x, &xpos)) != Some(true)
        {
            return None;
        }
        match self.sign(b) {
            Ok(1) => Some(pow(b, p)),
            Ok(-1) => Some(pow(&mul2(&Expr::minus_one(), b), p)),
            _ => None,
        }
    }
}
