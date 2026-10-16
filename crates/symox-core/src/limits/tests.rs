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
