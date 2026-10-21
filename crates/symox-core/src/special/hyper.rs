//! Generalised hypergeometric function `hyper(ap, bq, z)`.

use super::*;

/// `hyper.__new__` and `hyper.eval`: cancel parameters common to `ap` and
/// `bq`, order and unpolarify the rest, and unpolarify the argument when
/// the series converges everywhere or on the closed unit disc.
pub(super) fn eval(a: &[Expr]) -> Option<Expr> {
    let ap0 = tuple_items(&a[0])?;
    let bq0 = tuple_items(&a[1])?;
    let mut ap: Vec<Expr> = ap0.to_vec();
    let mut bq: Vec<Expr> = bq0.to_vec();
    let mut i = 0;
    while i < ap.len() {
        if let Some(j) = bq.iter().position(|b| *b == ap[i]) {
            ap.remove(i);
            bq.remove(j);
        } else {
            i += 1;
        }
    }
    for v in [&mut ap, &mut bq] {
        v.sort_by(crate::sortkey::cmp_default);
        for x in v.iter_mut() {
            *x = polar::unpolarify(x, false);
        }
    }
    let mut z = a[2].clone();
    let entire = ap.len() <= bq.len();
    if entire
        || (ap.len() == bq.len() + 1
            && matches!(crate::relational::le(&crate::functions::abs(&z), &Expr::one()).kind(), Kind::Bool(true)))
    {
        z = polar::unpolarify(&z, false);
    }
    if ap.as_slice() != ap0 || bq.as_slice() != bq0 || z != a[2] {
        return Some(raw(SpecialFn::Hyper, vec![tuple(ap), tuple(bq), z]));
    }
    None
}
