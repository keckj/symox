use crate::expr::Expr;
use crate::limits::{limit, Dir};
use crate::parser::parse_str;

fn lim(e: &str, x0: &str, dir: Dir) -> String {
    let x = Expr::symbol("x");
    match limit(&parse_str(e).unwrap(), &x, &parse_str(x0).unwrap(), dir) {
        Ok(r) => r.to_string(),
        Err(e) => format!("ERR {e}"),
    }
}

#[test]
fn classic_limits() {
    assert_eq!(lim("sin(x)/x", "0", Dir::Plus), "1");
    assert_eq!(lim("x**10*exp(-x)", "oo", Dir::Plus), "0");
    assert_eq!(lim("(1+1/x)**x", "oo", Dir::Plus), "E");
    assert_eq!(lim("x*log(x)", "0", Dir::Plus), "0");
    assert_eq!(lim("(x**2+1)/(2*x**2-3)", "oo", Dir::Plus), "1/2");
    assert_eq!(lim("1/x", "0", Dir::Minus), "-oo");
    assert_eq!(lim("1/x", "0", Dir::Both), "zoo");
    assert_eq!(lim("exp(x)*(exp(1/x - exp(-x)) - exp(1/x))", "oo", Dir::Plus), "-1");
}

#[test]
fn gruntz_limits() {
    assert_eq!(lim("2**x/3**x", "oo", Dir::Plus), "0");
    assert_eq!(lim("(3**x+5**x)**(1/x)", "oo", Dir::Plus), "5");
    assert_eq!(lim("exp(exp(x))/exp(exp(x - exp(-x)))", "oo", Dir::Plus), "E");
    assert_eq!(lim("log(log(x*exp(x*exp(x))+1))-exp(exp(log(log(x))+1/x))", "oo", Dir::Plus), "0");
    assert_eq!(lim("exp(-1/x**2)/x**10", "0", Dir::Plus), "0");
    assert_eq!(lim("x**(1/x)", "oo", Dir::Plus), "1");
    assert_eq!(lim("log(x)**10/x", "oo", Dir::Plus), "0");
}

#[test]
fn one_and_two_sided() {
    assert_eq!(lim("exp(-1/x)", "0", Dir::Minus), "oo");
    assert_eq!(lim("1/x**2", "0", Dir::Both), "oo");
    assert_eq!(lim("tan(x)", "pi/2", Dir::Plus), "-oo");
    assert_eq!(lim("abs(x)/x", "0", Dir::Minus), "-1");
    assert!(matches!(
        limit(&parse_str("abs(x)/x").unwrap(), &Expr::symbol("x"), &Expr::zero(), Dir::Both),
        Err(crate::limits::LimitError::Differ(..))
    ));
    assert_eq!(lim("a*x", "oo", Dir::Plus), "oo*sign(a)");
    assert_eq!(lim("sin(x)/x", "oo", Dir::Plus), "0");
    assert!(lim("sin(x)", "oo", Dir::Plus).starts_with("ERR"));
}

fn symbolic_power_difference(n: &Expr, k: &Expr, a: &Expr, b: &Expr) -> Expr {
    use crate::ops::{add2,div,pow,sub};
    let shifted=add2(n,a); let exponent=add2(k,&Expr::one());
    div(&pow(&shifted,k),&sub(&pow(&shifted,&exponent),&pow(&add2(n,b),&exponent)))
}

#[test]
fn homogeneous_symbolic_power_limits_and_shifted_family() {
    use crate::ops::{add2,div,mul2,sub};
    use crate::symbol::Assumptions as A;
    let n=Expr::symbol_with("n",A::POSITIVE);
    for facts in [A::POSITIVE,A::NONNEGATIVE] {
        let k=Expr::symbol_with("k",facts);
        for (a,b) in [(1,0),(3,1),(-2,-5)] {
            let e=symbolic_power_difference(&n,&k,&Expr::int(a),&Expr::int(b));
            let expected=div(&Expr::one(),&mul2(&Expr::int(a-b),&add2(&k,&Expr::one())));
            let actual=limit(&e,&n,&Expr::oo(),Dir::Plus).unwrap();
            assert!(crate::simplify::simplify(&sub(&actual,&expected)).is_zero(),"{e}: {actual} != {expected}");
        }
    }
}

#[test]
fn homogeneous_normalization_does_not_erase_parameter_poles() {
    use crate::ops::{add2,mul2};
    use crate::symbol::Assumptions as A;
    let n=Expr::symbol_with("n",A::POSITIVE);
    for facts in [A::NONE,A::REAL,A::NEGATIVE,A::COMPLEX] {
        let k=Expr::symbol_with("k",facts);
        let ratio=symbolic_power_difference(&n,&k,&Expr::one(),&Expr::zero());
        // The second expression has generic leading coefficient1, but at
        // k=-1 both its numerator and its denominator vanish.
        for e in [ratio.clone(),mul2(&add2(&k,&Expr::one()),&ratio)] {
            assert!(limit(&e,&n,&Expr::oo(),Dir::Plus).is_err(),"uncertified parameter continuation: {e}");
        }
    }
}

#[test]
fn limrest3_values_from_sympy() {
    for (e, x0, want) in [
        ("(Ei(x - exp(-exp(x))) - Ei(x))*exp(-x)*exp(exp(x))*x", "oo", "-1"),
        ("(erf(x - exp(-exp(x))) - erf(x))*exp(exp(x))*exp(x**2)", "oo", "-2/sqrt(pi)"),
        ("uppergamma(x, 1)/gamma(x)", "oo", "1"),
        ("subfactorial(x)", "oo", "oo"),
        ("((x + 1/2)*log(x) - x + log(2*pi)/2 - log(factorial(x)) + 1/(12*x))*x**3", "oo", "1/360"),
        ("harmonic(x, 3)", "oo", "zeta(3)"),
        ("-x + (x + 1/2)*log(x) - loggamma(x + 1) + log(2*pi)/2", "oo", "0"),
        ("lowergamma(x, 1)*x", "oo", "exp(-1)"),
        ("lowergamma(x, 2)*x/2**x", "oo", "exp(-2)"),
        ("harmonic(x)", "0", "0"),
        ("atan(a*x)/x", "oo", "0"),
        ("a*(log(x) - log(a**2*x**2 + 1)/2) - atan(a*x)/x", "oo", "-a*log(a**2)/2"),
    ] {
        assert_eq!(lim(e, x0, Dir::Plus), want, "{e}");
    }
    for (e, want) in [("acoth(sin(x)/x)", "oo"), ("acoth(x - 1)", "-oo"), ("acoth(1 + x) + log(x)/2", "log(2)/2")] {
        assert_eq!(lim(e, "0", Dir::Plus), want, "{e}");
        assert_eq!(lim(e, "0", Dir::Minus), want, "{e}");
    }
}

#[test]
fn limheld_values_from_sympy() {
    assert_eq!(lim("y/x/log(x)", "0", Dir::Plus), "-oo*sign(y)");
    assert_eq!(lim("-y*z*exp(x)/x", "oo", Dir::Plus), "-oo*sign(y)*sign(z)");
    assert_eq!(lim("Piecewise((x**2/2, x <= 1/2), (x/2 - 1/8, True))", "0", Dir::Plus), "0");
    assert_eq!(lim("Piecewise(((x - 9)/5, x < -1), ((x - 9)/5, x > 4), (sqrt(Abs(x - 3)), True))", "-1", Dir::Plus), "2");
    assert_eq!(lim("Piecewise(((x - 9)/5, x < -1), ((x - 9)/5, x > 4), (sqrt(Abs(x - 3)), True))", "-1", Dir::Minus), "-2");
    assert!(lim("Piecewise((x, y < 1), (2, True))", "0", Dir::Plus).starts_with("ERR"));
    assert_eq!(lim("x", "nan", Dir::Plus), "nan");
}

#[test]
fn limheld_order_mask_collision_is_not_an_order_term() {
    // A Dummy whose mask bit coincides with ORDER_MASK (limit_seq's `_n`
    // after one earlier Dummy) is a plain variable: limit(d, d, oo) is oo.
    let oo = parse_str("oo").unwrap();
    let mut found = false;
    for i in 0..4096u32 {
        let d = Expr::dummy("_n", crate::symbol::Assumptions::default());
        if !d.may_have_order() {
            continue;
        }
        found = true;
        assert_eq!(limit(&d, &d, &oo, Dir::Minus).unwrap().to_string(), "oo", "dummy {i}");
        let e = crate::ops::add(&[d.clone(), Expr::int(1)]);
        assert_eq!(limit(&e, &d, &oo, Dir::Minus).unwrap().to_string(), "oo", "dummy {i}");
        break;
    }
    assert!(found, "no Dummy mask collided with ORDER_MASK");
}
