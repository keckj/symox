//! Integral transforms (Laplace, Fourier, sine, cosine, Hankel, Mellin).
//!
//! The engines follow SymPy's conventions:
//!
//! * `laplace_transform` returns `(F, plane, cond)` where `Re(s) > plane`
//!   is the half plane of convergence. It is rule based (Bateman's tables,
//!   time/frequency shifts, differentiation, trigonometric pole pairing).
//! * `inverse_laplace_transform` returns `(f, cond)`; rational functions
//!   go through partial fractions, `exp(-a*s)` factors become time shifts.
//! * Fourier, sine, cosine and Hankel transforms return `(F, cond)` and use
//!   SymPy's normalisations (`F(k) = int f(x) exp(-2*pi*I*x*k) dx`, ...).
//! * `mellin_transform` returns `(F, (a, b), cond)` with the fundamental
//!   strip `a < Re(s) < b`; `inverse_mellin_transform` returns `(f, cond)`.
//!
//! When no rule applies, the corresponding held transform head
//! (`LaplaceTransform(f, t, s)`, ...) is returned, exactly like SymPy.

use crate::expr::{Expr, Func, Kind};
use crate::functions::func;
use crate::special::SpecialFn;

mod util;
pub mod laplace;
mod laplace_integration;
mod integral_backend;
pub mod inverse_laplace;
mod ilt_irrational;
mod ilt_integration;
pub mod fourier;
pub mod mellin;
mod inverse_mellin;
pub mod hankel;

pub use fourier::{FourierKind, fourier_type_transform};
pub use hankel::hankel_type_transform;
pub use inverse_laplace::inverse_laplace_transform;
pub use laplace::laplace_transform;
pub use mellin::{inverse_mellin_transform, mellin_transform};

/// Placeholder printed as `_None` for an unspecified inverse Laplace plane
/// (SymPy uses a `Dummy('None')`).
pub fn none_plane() -> Expr {
    Expr::symbol("_None")
}

/// The message of an error that is not an `IntegralTransformError`: the
/// Python binding raises `kind(msg)` instead of falling back.
pub(crate) fn fatal(kind: &str, msg: &str) -> String {
    format!("\u{1}{kind}\u{1}{msg}")
}

/// Is `msg` a fatal error (see `fatal`)?
pub(crate) fn is_fatal(msg: &str) -> bool {
    msg.starts_with('\u{1}')
}

/// Held transform node `name(args...)`.
pub fn held(s: SpecialFn, args: Vec<Expr>) -> Expr {
    debug_assert!(s.is_transform());
    func(Func::Special(s), args)
}

/// Is `e` a held transform node? Returns its head and arguments.
pub fn as_held(e: &Expr) -> Option<(SpecialFn, &[Expr])> {
    match e.kind() {
        Kind::Func(Func::Special(s), args) if s.is_transform() => Some((*s, args.as_slice())),
        _ => None,
    }
}

/// Evaluate a held transform node `s(args...)` (`doit()` without
/// conditions); `None` keeps it held.
pub fn doit_default(s: SpecialFn, args: &[Expr]) -> Option<Expr> {
    let r = transform_doit(s, args, None, false).ok()?.into_iter().next()?;
    if let Some((h, a)) = as_held(&r) {
        if h == s && a == args {
            return None;
        }
    }
    Some(r)
}

/// Direct computation of a transform: the result followed by the
/// auxiliary items SymPy returns with `noconds=False` (strip bounds and
/// condition for Mellin, the condition otherwise). `None` stands for
/// SymPy's `IntegralTransformError`.
fn direct(s: SpecialFn, args: &[Expr], simplify: bool) -> Result<Vec<Expr>, String> {
    use SpecialFn::*;
    if s == InverseMellinTransform {
        return inverse_mellin::compute(&args[0], &args[1], &args[2], &args[3], &args[4])
            .map(|(f, c)| vec![f, c])
            .map_err(inverse_mellin::ImtError::message);
    }
    if let Some(v) = direct_table(s, args, simplify) {
        return Ok(v);
    }
    let unsupported = || "could not compute integral".to_string();
    // the reference does not try the integral when the function contains an
    // undefined function of the variable
    if laplace_integration::has_applied_undef_in(&args[0], &args[1]) {
        return Err(unsupported());
    }
    match s {
        MellinTransform => {
            let (f, (a, b), c) = integral_backend::mellin_integration(&args[0], &args[1], &args[2], simplify)?;
            Ok(vec![f, a, b, c])
        }
        FourierTransform | InverseFourierTransform | SineTransform | InverseSineTransform | CosineTransform
        | InverseCosineTransform => {
            let kind = FourierKind::from_head(s).ok_or_else(unsupported)?;
            let (a, b) = kind.ab();
            let (f, c) = match kind {
                FourierKind::Fourier | FourierKind::InverseFourier => {
                    integral_backend::fourier_integration(&args[0], &args[1], &args[2], &a, &b, simplify)?
                }
                FourierKind::Sine | FourierKind::InverseSine => {
                    integral_backend::sine_cosine_integration(&args[0], &args[1], &args[2], &a, &b, true, simplify)?
                }
                FourierKind::Cosine | FourierKind::InverseCosine => {
                    integral_backend::sine_cosine_integration(&args[0], &args[1], &args[2], &a, &b, false, simplify)?
                }
            };
            Ok(vec![f, c])
        }
        HankelTransform | InverseHankelTransform => {
            let (f, c) = integral_backend::hankel_integration(&args[0], &args[1], &args[2], &args[3], simplify)?;
            Ok(vec![f, c])
        }
        _ => Err(unsupported()),
    }
}

/// The closed form tables of the transforms.
fn direct_table(s: SpecialFn, args: &[Expr], simplify: bool) -> Option<Vec<Expr>> {
    use SpecialFn::*;
    Some(match s {
        MellinTransform => {
            let (f, (a, b), c) = mellin_transform(&args[0], &args[1], &args[2], simplify)?;
            vec![f, a, b, c]
        }
        InverseMellinTransform => {
            let (f, c) = inverse_mellin_transform(&args[0], &args[1], &args[2], (&args[3], &args[4]), simplify)?;
            vec![f, c]
        }
        FourierTransform | InverseFourierTransform | SineTransform | InverseSineTransform | CosineTransform
        | InverseCosineTransform => {
            let kind = FourierKind::from_head(s)?;
            let (f, c) = fourier_type_transform(kind, &args[0], &args[1], &args[2], simplify)?;
            vec![f, c]
        }
        HankelTransform | InverseHankelTransform => {
            let (f, c) = hankel_type_transform(s == InverseHankelTransform, &args[0], &args[1], &args[2], &args[3], simplify)?;
            vec![f, c]
        }
        _ => return None,
    })
}

/// SymPy's `doit` of a held transform node with the hints `simplify`
/// (`None`: the transform's default) and `needeval`.
///
/// Returns `[F]` when there is no auxiliary data (held or partially held
/// results) and `[F, extras...]` otherwise: `[F, plane, cond]` for Laplace,
/// `[f, cond]` for the inverse Laplace transform, `[F, a, b, cond]` for
/// Mellin and `[F, cond]` for the others. `Err` is SymPy's
/// `IntegralTransformError`.
pub fn transform_doit(s: SpecialFn, args: &[Expr], simplify: Option<bool>, needeval: bool) -> Result<Vec<Expr>, String> {
    use SpecialFn::*;
    let name = transform_name(s);
    match s {
        LaplaceTransform => {
            let (f, p, c) = laplace_transform(&args[0], &args[1], &args[2], simplify.unwrap_or(false));
            return Ok(vec![f, p, c]);
        }
        InverseLaplaceTransform => {
            let plane = if args[3] == none_plane() { None } else { Some(&args[3]) };
            let (f, c) = inverse_laplace_transform(&args[0], &args[1], &args[2], plane, simplify.unwrap_or(false));
            return Ok(vec![f, c]);
        }
        HankelTransform | InverseHankelTransform => {
            return direct(s, args, simplify.unwrap_or(true)).map_err(|m| format!("{name} Transform could not be computed: {m}"));
        }
        _ => {}
    }
    let simp = simplify.unwrap_or(true);
    match direct(s, args, simp) {
        Ok(v) => return Ok(v),
        Err(m) if is_fatal(&m) => return Err(m),
        Err(_) => {}
    }
    let x = &args[1];
    let f = if args[0].is_add() { args[0].clone() } else { crate::simp::expand_mul(&args[0]) };
    if f.is_add() {
        let mut res = Vec::new();
        let mut extra: Vec<Vec<Expr>> = Vec::new();
        for t in f.args() {
            let mut a = args.to_vec();
            a[0] = t;
            let r = transform_doit(s, &a, simplify, needeval)?;
            res.push(r[0].clone());
            if r.len() > 1 {
                extra.push(r[1..].to_vec());
            }
        }
        let sum = crate::ops::add(&res);
        let sum = if simp { crate::simplify::simplify(&sum) } else { sum };
        if extra.is_empty() {
            return Ok(vec![sum]);
        }
        if let Some(ex) = collapse(s, extra) {
            let mut v = vec![sum];
            v.extend(ex);
            return Ok(v);
        }
    }
    if needeval {
        return Err(format!("{name} Transform could not be computed: needeval"));
    }
    let (coeff, rest) = util::split_mul(&f, x);
    let mut a = args.to_vec();
    a[0] = rest;
    Ok(vec![crate::ops::mul(&[coeff, held(s, a)])])
}

/// SymPy's `_collapse_extra`: intersect the strips and conjoin conditions.
fn collapse(s: SpecialFn, extra: Vec<Vec<Expr>>) -> Option<Vec<Expr>> {
    if s == SpecialFn::MellinTransform {
        let mut a = Vec::new();
        let mut b = Vec::new();
        let mut c = Vec::new();
        for e in extra {
            if e.len() != 3 {
                return None;
            }
            a.push(e[0].clone());
            b.push(e[1].clone());
            c.push(e[2].clone());
        }
        let a = crate::functions::max(a);
        let b = crate::functions::min(b);
        let c = crate::logic::and(c);
        if util::nonneg(&crate::ops::sub(&a, &b)) || util::truth(&c) == Some(false) {
            return None;
        }
        return Some(vec![a, b, c]);
    }
    let c = crate::logic::and(extra.into_iter().map(|e| e[e.len() - 1].clone()).collect());
    if util::truth(&c) == Some(false) {
        return None;
    }
    Some(vec![c])
}

/// SymPy's `_name` of a transform class.
pub fn transform_name(s: SpecialFn) -> &'static str {
    use SpecialFn::*;
    match s {
        LaplaceTransform => "Laplace",
        InverseLaplaceTransform => "Inverse Laplace",
        FourierTransform => "Fourier",
        InverseFourierTransform => "Inverse Fourier",
        SineTransform => "Sine",
        InverseSineTransform => "Inverse Sine",
        CosineTransform => "Cosine",
        InverseCosineTransform => "Inverse Cosine",
        HankelTransform => "Hankel",
        InverseHankelTransform => "Inverse Hankel",
        MellinTransform => "Mellin",
        InverseMellinTransform => "Inverse Mellin",
        _ => "",
    }
}

#[cfg(test)]
mod tests;

/// `mellin_transform(f, x, s, integrator=integrator, simplify=False,
/// needeval=True)` of the recursive Meijer G rewrite: `(F, (a, b), aux)`,
/// or `None` where the reference raises `IntegralTransformError`.
pub(crate) fn mellin_needeval(
    f: &Expr,
    x: &Expr,
    s: &Expr,
    integrator: &dyn Fn(&Expr, &Expr) -> Expr,
) -> Option<(Expr, (Expr, Expr), Expr)> {
    // `_try_directly`
    if !laplace_integration::has_applied_undef_in(f, x) {
        if let Ok(r) = integral_backend::mellin_integration_with(f, x, s, false, integrator) {
            return Some(r);
        }
    }
    if !f.is_add() {
        return None;
    }
    // linearity, then `MellinTransform._collapse_extra`
    let mut res = Vec::new();
    let (mut a, mut b, mut cond) = (Vec::new(), Vec::new(), Vec::new());
    for t in f.args() {
        let (r, (sa, sb), c) = mellin_needeval(&t, x, s, integrator)?;
        res.push(r);
        a.push(sa);
        b.push(sb);
        cond.push(c);
    }
    let (a, b) = (crate::functions::max(a), crate::functions::min(b));
    let cond = crate::logic::and(cond);
    let empty = util::truth(&crate::relational::ge(&a, &b)) == Some(true);
    if empty || matches!(cond.kind(), Kind::Bool(false)) {
        return None;
    }
    Some((crate::ops::add(&res), (a, b), cond))
}

/// `my_imt` of the recursive Meijer G rewrite: `inverse_mellin_transform(F,
/// s, x, strip, as_meijerg=True, needeval=True)`, retried on `F`
/// simplified when the strip check fails. `None` for every error; the
/// reference lets the errors other than `IntegralTransformError` of the
/// retry propagate out of the integration.
pub(crate) fn inverse_mellin_meijerg(f: &Expr, s: &Expr, x: &Expr, a: &Expr, b: &Expr) -> Option<Expr> {
    let run = |f: &Expr| inverse_mellin::compute_meijerg(f, s, x, a, b);
    match run(f) {
        Ok((g, _)) => Some(g),
        Err(inverse_mellin::ImtError::Strip(_)) => {
            let e = crate::poly::api::cancel(&crate::expand::expand(f)).unwrap_or_else(|_| crate::expand::expand(f));
            run(&crate::simplify::simplify(&e)).ok().map(|(g, _)| g)
        }
        Err(_) => None,
    }
}
