//! The trigonometric steps of the reference `solve` (SymPy 1.14
//! `sympy/solvers/solvers.py`, BSD-3-Clause, see LICENSE-SYMPY):
//!
//! - `solve` expands the angles of trig sums (`TR10`) and doubles the angles
//!   that are twice the gcd of the arguments (`TR11`) when more than one
//!   distinct trig atom depends on the unknown;
//! - `_solve` rewrites an equation with several trig generators (and no
//!   other function) to `tan`: `sin(x)**p + c*cos(x)**p` through
//!   `TR2i(sin/cos) + 1`, anything else through `TR1(f).rewrite(tan)`;
//! - `_tsolve` falls back to `rewrite(exp)`.

use crate::expr::{Expr, Func, Kind};
use crate::number::Number;
use crate::ops::{add2, as_coeff_mul, div, mul2};
use crate::simplify::{numer_denom, together};
use crate::simptools::fu::{tr1, tr10, tr11, tr2i};
use crate::special::rewrite::rewrite;
use crate::traverse::has;

fn is_trig(func: &Func) -> bool {
    matches!(func, Func::Sin | Func::Cos | Func::Tan | Func::Cot | Func::Sec | Func::Csc)
}

fn trig_atoms(e: &Expr, x: &Expr, out: &mut Vec<Expr>) {
    crate::traverse::preorder(e, &mut |s| {
        if let Kind::Func(func, a) = s.kind() {
            if is_trig(func) && a.len() == 1 && has(s, x) && !out.contains(s) {
                out.push(s.clone());
            }
        }
    });
}

/// `gcd(*args)` of the arguments `c_i*m`; `None` when the rests differ.
fn argument_gcd(args: &[Expr]) -> Option<Expr> {
    let mut parts = args.iter().map(as_coeff_mul);
    let (mut g, rest) = parts.next()?;
    g = g.abs();
    for (c, m) in parts {
        if m != rest {
            return None;
        }
        g = g.gcd(&c.abs());
    }
    (!g.is_zero()).then(|| mul2(&Expr::num(g), &rest))
}

/// The angle expansion of the reference `solve` for one equation.
pub(super) fn expand_angles(f: &Expr, x: &Expr) -> Expr {
    let mut targs = Vec::new();
    trig_atoms(f, x, &mut targs);
    if targs.len() <= 1 {
        return f.clone();
    }
    let (sums, other): (Vec<Expr>, Vec<Expr>) = targs.into_iter().partition(|t| t.args()[0].is_add());
    let mut rep: Vec<(Expr, Expr)> = Vec::new();
    for t in &sums {
        let a = &t.args()[0];
        let (ind, dep) = crate::decompose::as_independent(a, std::slice::from_ref(x), Some(true));
        if dep == *x || crate::ops::neg(&dep) == *x {
            let Kind::Func(func, _) = t.kind() else { continue };
            let expanded = if ind.as_number().is_some() {
                tr10(&crate::functions::func(func.clone(), vec![add2(&dep, &ind)]), true)
            } else {
                let n = Expr::dummy("_n", crate::Assumptions::NONE);
                let e = tr10(&crate::functions::func(func.clone(), vec![add2(&dep, &n)]), true);
                crate::subs::xreplace(&e, &[(n, ind)].into_iter().collect::<rustc_hash::FxHashMap<_, _>>())
            };
            rep.push((t.clone(), expanded));
        }
    }
    if !other.is_empty() && other.len() <= 2 {
        let args: Vec<Expr> = other.iter().map(|t| t.args()[0].clone()).collect();
        if let Some(base) = if args.len() > 1 { argument_gcd(&args) } else { Some(args[0].clone()) } {
            for t in &other {
                rep.push((t.clone(), tr11(t, Some(&base))));
            }
        }
    }
    rep.retain(|(a, b)| a != b);
    if rep.is_empty() {
        return f.clone();
    }
    crate::subs::xreplace(f, &rep.into_iter().collect::<rustc_hash::FxHashMap<_, _>>())
}

/// `_as_base_q`: the base of a generator `b**(p/q)`.
fn generator_base(g: &Expr) -> (Expr, bool) {
    match g.kind() {
        Kind::Pow(b, e) => {
            if let Some(n) = e.as_number() {
                if n.is_rational() {
                    return (b.clone(), n.is_integer());
                }
            }
            let (c, rest) = as_coeff_mul(e);
            if c.is_rational() && c != Number::ONE {
                let q_one = c.is_integer();
                return (crate::ops::pow(b, &rest), q_one);
            }
            (g.clone(), true)
        }
        Kind::Func(Func::Exp, a) => {
            let (c, rest) = as_coeff_mul(&a[0]);
            if c.is_rational() && c != Number::ONE {
                return (crate::functions::exp(&rest), c.is_integer());
            }
            (g.clone(), true)
        }
        _ => (g.clone(), true),
    }
}

/// The trig branch of the multiple-generator case of the reference
/// `_solve`. `Some(rewritten)` is the equation to solve without a check
/// (the caller checks against the original); `None` when the branch does
/// not apply.
pub(super) fn tan_rewrite(num: &Expr, x: &Expr, gens: &[Expr]) -> Option<Expr> {
    if gens.len() <= 1 {
        return None;
    }
    let mut bases: Vec<Expr> = Vec::new();
    let mut all_q_one = true;
    for g in gens {
        let (b, q_one) = generator_base(g);
        all_q_one &= q_one;
        if !bases.contains(&b) {
            bases.push(b);
        }
    }
    if bases.len() <= 1 && all_q_one {
        return None;
    }
    let funcs: Vec<&Expr> = bases.iter().filter(|b| matches!(b.kind(), Kind::Func(..))).collect();
    if funcs.iter().any(|f| !matches!(f.kind(), Kind::Func(func, _) if is_trig(func))) {
        return None;
    }
    let trig: Vec<Expr> = funcs.into_iter().cloned().collect();
    if trig.len() <= 1 {
        return None;
    }
    let mut newf = None;
    if let Kind::Add(..) = num.kind() {
        let mut terms = num.args();
        if terms.len() == 2 {
            let atoms = |t: &Expr| {
                let mut found = Vec::new();
                crate::traverse::preorder(t, &mut |s| {
                    if trig.contains(s) && !found.contains(s) {
                        found.push(s.clone());
                    }
                });
                found
            };
            let (ta, tb) = (atoms(&terms[0]), atoms(&terms[1]));
            if ta.len() == 1 && tb.len() == 1 {
                let (mut a, mut b) = (ta[0].clone(), tb[0].clone());
                if matches!(a.kind(), Kind::Func(Func::Cos, _)) {
                    std::mem::swap(&mut a, &mut b);
                    terms.swap(0, 1);
                }
                if let (Kind::Func(Func::Sin, sa), Kind::Func(Func::Cos, ca)) = (a.kind(), b.kind()) {
                    if sa[0] == ca[0] {
                        let ratio = add2(&tr2i(&div(&terms[0], &terms[1]), false), &Expr::one());
                        let (n, d) = numer_denom(&together(&ratio));
                        if d.as_number().is_some() {
                            newf = Some(n);
                        }
                    }
                }
            }
        }
    }
    let newf = newf.unwrap_or_else(|| rewrite(&tr1(num), &[], "tan", true));
    (newf != *num).then_some(newf)
}

/// The `rewrite(exp)` fallback of the reference `_tsolve` for an equation
/// with trig functions of the unknown.
pub(super) fn exp_rewrite(num: &Expr, x: &Expr) -> Option<Expr> {
    let mut targs = Vec::new();
    trig_atoms(num, x, &mut targs);
    if targs.is_empty() {
        return None;
    }
    // The exp form is solved only as a polynomial in `exp(I*x/q)`: the
    // unknown must occur only in trig functions of `r*x + c`, `r` rational.
    // Anything else (`y*sin(y/x)`) is left unsolved without the rewrite.
    for t in &targs {
        let a = &t.args()[0];
        let k = crate::diff::diff(a, x);
        if !k.is_rational_number() || has(&crate::expand::expand(&crate::ops::sub(a, &mul2(&k, x))), x) {
            return None;
        }
    }
    let u = Expr::dummy("_t", crate::Assumptions::NONE);
    let map: rustc_hash::FxHashMap<Expr, Expr> = targs.iter().map(|t| (t.clone(), u.clone())).collect();
    if has(&crate::subs::xreplace(num, &map), x) {
        return None;
    }
    let r = rewrite(num, &[], "exp", true);
    (r != *num).then_some(r)
}
