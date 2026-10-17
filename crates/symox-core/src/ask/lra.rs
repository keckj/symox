//! Final `lra_satask` stage of `ask`.
//!
//! SymPy 1.14 calls `lra_satask` only for propositions and assumptions made
//! of `Q.eq`/`Q.gt`/... relations between linear real expressions (plus the
//! `Q.positive`/`Q.negative`/... unary predicates). symox does not port the
//! linear real arithmetic solver yet: this stage reports "unknown", so such
//! queries answer `None` where SymPy may decide them.

use super::engine::{Cx, R};
use crate::expr::Expr;

pub(super) fn lra_satask(_cx: &Cx, _prop: &Expr, _assumptions: &Expr) -> R {
    Ok(None)
}
