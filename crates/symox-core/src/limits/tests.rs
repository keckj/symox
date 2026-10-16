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
