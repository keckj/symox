//! Arbitrary precision evaluation of the special functions (real
//! arguments). Called from [`crate::evalf_mp`].

use crate::evalf_mp::MpComplex;
use crate::expr::Expr;

use super::SpecialFn;

/// Evaluate `s(args)` to `prec` bits, `None` when unsupported.
pub fn eval(_s: SpecialFn, _args: &[Expr], _prec: u32) -> Option<MpComplex> {
    None
}
