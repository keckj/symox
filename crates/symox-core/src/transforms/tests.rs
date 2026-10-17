use super::*;
use crate::expr::Expr;
use crate::special::SpecialFn;

fn p(s: &str) -> Expr {
    crate::parser::sympify(s).unwrap()
}

#[test]
fn laplace_power() {
    let r = transform_doit(SpecialFn::LaplaceTransform, &[p("t**2"), p("t"), p("s")], None, false).unwrap();
    assert_eq!(r[0], p("2/s**3"));
    assert_eq!(r[1], Expr::zero());
}

#[test]
fn inverse_laplace_exponential() {
    let r = transform_doit(SpecialFn::InverseLaplaceTransform, &[p("1/(s+1)"), p("s"), p("t"), none_plane()], None, false).unwrap();
    assert_eq!(r[0], p("exp(-t)*Heaviside(t)"));
}

#[test]
fn mellin_exponential() {
    let r = transform_doit(SpecialFn::MellinTransform, &[p("exp(-x)"), p("x"), p("s")], None, false).unwrap();
    assert_eq!(r[0], p("gamma(s)"));
    assert_eq!(r[1], Expr::zero());
}

#[test]
fn fourier_gaussian_and_needeval() {
    let r = transform_doit(SpecialFn::FourierTransform, &[p("exp(-x**2)"), p("x"), p("k")], None, false).unwrap();
    assert_eq!(r[0], p("sqrt(pi)*exp(-pi**2*k**2)"));
    assert!(transform_doit(SpecialFn::FourierTransform, &[p("f(x)"), p("x"), p("k")], None, true).is_err());
}

#[test]
fn hankel_failure_is_error() {
    assert!(transform_doit(SpecialFn::HankelTransform, &[p("f(r)"), p("r"), p("k"), p("0")], None, false).is_err());
}
