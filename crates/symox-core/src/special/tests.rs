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
