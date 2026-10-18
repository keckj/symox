use crate::prelude::*;

fn s(e: &str) -> String {
    parse_str(e).unwrap().to_string()
}

#[test]
fn combinatorial() {
    assert_eq!(s("factorial(5)"), "120");
    assert_eq!(s("factorial(n)"), "factorial(n)");
    assert_eq!(s("binomial(5, 2)"), "10");
    assert_eq!(s("binomial(n, 2)"), "binomial(n, 2)");
    assert_eq!(s("rf(x, 3)"), "x*(x + 1)*(x + 2)");
    assert_eq!(s("gamma(n + 1)"), "gamma(n + 1)");
    assert_eq!(s("zeta(2)"), "pi**2/6");
}

#[test]
fn invalid_arity_is_not_indexed() {
    use crate::special::{self, numeric, SpecialFn};
    for head in SpecialFn::ALL {
        for len in 0..=5 {
            if head.nargs().contains(&len) {
                continue;
            }
            let args = vec![Expr::one(); len];
            assert!(numeric::eval_c64(*head, &args).is_none());
            assert!(special::diff::partial(*head, &args, 0).is_none());
            assert!(numeric::apply(*head, &vec![num_complex::Complex64::new(1.0, 0.0); len]).is_none());
        }
    }
}

#[test]
fn mathieu_and_spherical_hankel_differential_rules() {
    use crate::special::{self, SpecialFn as S};
    use crate::functions::{func,cos};
    use crate::expr::Func;
    use crate::ops::{mul,sub,add};
    let x=Expr::symbol("x");let a=Expr::symbol("a");let q=Expr::symbol("q");
    for (head,prime) in [(S::MathieuC,S::MathieuCPrime),(S::MathieuS,S::MathieuSPrime)] {
        let args=vec![a.clone(),q.clone(),x.clone()];
        let value=func(Func::Special(head),args.clone());
        assert_eq!(crate::diff::diff(&value,&x),func(Func::Special(prime),args.clone()));
        let factor=sub(&mul(&[Expr::int(2),q.clone(),cos(&mul(&[Expr::int(2),x.clone()]))]),&a);
        assert_eq!(crate::diff::diff(&crate::diff::diff(&value,&x),&x),mul(&[factor,value]));
        let at_zero=func(Func::Special(head),vec![a.clone(),Expr::zero(),x.clone()]);
        assert!(!at_zero.to_string().contains("mathieu"));
    }
    let h=func(Func::Special(S::Hn1),vec![Expr::zero(),x.clone()]);
    let expected=add(&[special::jn(&Expr::zero(),&x),mul(&[Expr::i(),special::yn(&Expr::zero(),&x)])]);
    assert_eq!(special::expand_func(&h),special::expand_func(&expected));
}

#[test]
fn advanced_numeric_probability_and_double_series() {
    use crate::special::{numeric,SpecialFn as S};
    use num_complex::Complex64 as C;
    let r=|v|C::new(v,0.0);
    let f=numeric::apply(S::AppellF1,&[r(1.0),r(1.0),r(1.0),r(2.0),r(0.2),r(0.3)]).unwrap();
    let expected=((1.0_f64-0.3).ln()-(1.0_f64-0.2).ln())/(0.2-0.3);
    assert!((f.re-expected).abs()<1e-13);
    let q=numeric::apply(S::MarcumQ,&[r(1.0),r(0.0),r(2.0)]).unwrap();
    assert!((q.re-(-2.0_f64).exp()).abs()<1e-14);
    assert!(numeric::apply(S::AppellF1,&[r(1.0),r(1.0),r(1.0),r(2.0),r(2.0),r(3.0)]).is_none());
}

#[test]
fn meijer_slater_regions_and_boundary_are_explicit() {
    use crate::special::{self,SpecialFn as S};
    use crate::functions::func;
    use crate::expr::Func;
    let args=|z:Expr|vec![
        special::tuple(vec![special::tuple(vec![]),special::tuple(vec![Expr::int(2)])]),
        special::tuple(vec![special::tuple(vec![Expr::zero()]),special::tuple(vec![])]),z];
    let small=args(Expr::rational(1,2));
    let value=special::numeric::eval_c64(S::MeijerG,&small).unwrap();
    assert!((value.re-0.5).abs()<1e-14);
    let big=args(Expr::int(2));
    assert_eq!(special::numeric::eval_c64(S::MeijerG,&big).unwrap().re,0.0);
    let boundary=args(Expr::one());
    assert!(special::numeric::eval_c64(S::MeijerG,&boundary).is_none());
    let held=func(Func::Special(S::MeijerG),boundary);
    assert_eq!(special::hyperexpand(&held),held);
}
