use crate::parser::parse_str;
use crate::series::series;
use crate::expr::Expr;

fn p(s: &str) -> Expr {
    parse_str(s).unwrap()
}

fn ser(e: &str, x0: &str, n: i64) -> String {
    let x = Expr::symbol("x");
    series(&p(e), &x, &p(x0), n, false).map(|r| r.to_string()).unwrap_or_else(|e| format!("ERR {e}"))
}

#[test]
fn basic_series() {
    assert_eq!(ser("sin(x)", "0", 6), "x - x**3/6 + x**5/120 + O(x**6)");
    assert_eq!(ser("exp(x)", "0", 3), "1 + x + x**2/2 + O(x**3)");
    assert_eq!(ser("1/(1-x-x**2)", "0", 6), "1 + x + 2*x**2 + 3*x**3 + 5*x**4 + 8*x**5 + O(x**6)");
    assert_eq!(ser("cot(x)", "0", 6), "1/x - x/3 - x**3/45 - 2*x**5/945 + O(x**6)");
    assert_eq!(ser("log(1+x)/x", "0", 6), "1 - x/2 + x**2/3 - x**3/4 + x**4/5 - x**5/6 + O(x**6)");
}

#[test]
fn order_rules() {
    assert_eq!(p("x + O(x**2)").to_string(), "x + O(x**2)");
    assert_eq!(p("x**3 + O(x**2)").to_string(), "O(x**2)");
    assert_eq!(p("O(x)*x").to_string(), "O(x**2)");
    assert_eq!(p("O(2*x + x**2)").to_string(), "O(x)");
    assert_eq!(p("O(x)**2").to_string(), "O(x**2)");
    assert_eq!(p("y*O(x)").to_string(), "O(x)");
}
