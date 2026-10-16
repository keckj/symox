//! Generalised hypergeometric function `hyper(ap, bq, z)`.

use super::*;

/// Cancel parameters common to `ap` and `bq`.
pub(super) fn eval(a: &[Expr]) -> Option<Expr> {
    let ap = tuple_items(&a[0])?;
    let bq = tuple_items(&a[1])?;
    let mut ap: Vec<Expr> = ap.to_vec();
    let mut bq: Vec<Expr> = bq.to_vec();
    let mut changed = false;
    let mut i = 0;
    while i < ap.len() {
        if let Some(j) = bq.iter().position(|b| *b == ap[i]) {
            ap.remove(i);
            bq.remove(j);
            changed = true;
        } else {
            i += 1;
        }
    }
    if changed {
        return Some(raw(SpecialFn::Hyper, vec![tuple(ap), tuple(bq), a[2].clone()]));
    }
    None
}
