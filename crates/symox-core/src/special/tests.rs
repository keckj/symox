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
