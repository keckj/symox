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
fn points_and_exactness() {
    assert_eq!(ser("log(x)", "1", 4), "-1 - (x - 1)**2/2 + (x - 1)**3/3 + x + O((x - 1)**4, (x, 1))");
    assert_eq!(ser("1/(x+1)", "oo", 3), "-1/x**2 + 1/x + O(x**(-3), (x, oo))");
    assert_eq!(ser("atan(x)", "oo", 4), "1/(3*x**3) - 1/x + pi/2 + O(x**(-4), (x, oo))");
    assert_eq!(ser("1/x + x", "0", 6), "x + 1/x");
    assert_eq!(ser("x**10 + x", "0", 6), "x + O(x**6)");
    assert_eq!(ser("(x+1)/x", "0", 3), "1/x + 1 + O(x**3)");
    assert_eq!(ser("sqrt(sin(x))", "0", 3), "sqrt(x) - x**(5/2)/12 + O(x**3)");
    assert_eq!(ser("asin(1-x)", "0", 2), "pi/2 - sqrt(2)*sqrt(x) - sqrt(2)*x**(3/2)/12 + O(x**2)");
}

#[test]
fn log_terms() {
    assert_eq!(ser("x**x", "0", 3), "1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)");
    assert_eq!(ser("log(x)*sin(x)", "0", 4), "x*log(x) - x**3*log(x)/6 + O(x**4*log(x))");
    assert_eq!(ser("log(sin(x))", "0", 4), "log(x) - x**2/6 + O(x**4)");
    assert_eq!(ser("x*log(x) + x", "0", 6), "x*(log(x) + 1)");
}

#[test]
fn symbolic_coefficients() {
    assert_eq!(ser("exp(a*x)", "0", 3), "1 + a*x + a**2*x**2/2 + O(x**3)");
    assert_eq!(ser("sin(x+y)", "0", 3), "sin(y) + x*cos(y) - x**2*sin(y)/2 + O(x**3)");
    assert_eq!(ser("tan(x)", "1", 2), "tan(1) + (1 + tan(1)**2)*(x - 1) + O((x - 1)**2, (x, 1))");
}

#[test]
fn gamma_series() {
    assert_eq!(ser("gamma(x)", "0", 2), "1/x - EulerGamma + x*(EulerGamma**2/2 + pi**2/12) + O(x**2)");
}

#[test]
fn taylor_term_and_remove_o() {
    let x = Expr::symbol("x");
    assert_eq!(crate::series::taylor_term(&p("sin(x)"), &x, 3).unwrap().to_string(), "-x**3/6");
    let s = series(&p("cos(x)"), &x, &Expr::zero(), 4, false).unwrap();
    assert_eq!(crate::series::remove_o(&s).to_string(), "1 - x**2/2");
    assert_eq!(crate::series::get_o(&s).unwrap().to_string(), "O(x**4)");
    assert_eq!(crate::series::getn(&s).unwrap().to_string(), "4");
}

#[test]
fn order_rules() {
    assert_eq!(p("x + O(x**2)").to_string(), "x + O(x**2)");
    assert_eq!(p("x**3 + O(x**2)").to_string(), "O(x**2)");
    assert_eq!(p("O(x)*x").to_string(), "O(x**2)");
    assert_eq!(p("O(2*x + x**2)").to_string(), "O(x)");
    assert_eq!(p("O(x)**2").to_string(), "O(x**2)");
    assert_eq!(p("y*O(x)").to_string(), "O(x)");
    assert_eq!(p("O(x) + x*log(x)").to_string(), "x*log(x) + O(x)");
    assert_eq!(p("O(x) + x**2*log(x)").to_string(), "O(x)");
    assert_eq!(p("O(x**2) + 1/x").to_string(), "1/x + O(x**2)");
    assert_eq!(p("O(x, (x, 1))").to_string(), "O(1, (x, 1))");
    assert_eq!(p("O(x**2 + x, (x, oo))").to_string(), "O(x**2, (x, oo))");
    assert_eq!(p("O(x)*O(y)").to_string(), "O(x*y, x, y)");
    assert_eq!(p("O(x) - O(x)").to_string(), "O(x)");
    assert_eq!(p("2*(x + O(x**2))").to_string(), "2*x + O(x**2)");
    assert_eq!(crate::printer::srepr(&p("O(x)")), "Order(Symbol('x'), Tuple(Symbol('x'), Integer(0)))");
}

#[test]
#[ignore]
fn timing_series() {
    let x = Expr::symbol("x");
    for (src, n) in [("1/cos(x)", 40), ("exp(sin(x))*cos(x)", 30), ("tan(x)", 30)] {
        let e = p(src);
        let t = std::time::Instant::now();
        for _ in 0..20 {
            let _ = series(&e, &x, &Expr::zero(), n, false).unwrap();
        }
        let dt = t.elapsed() / 20;
        let t2 = std::time::Instant::now();
        for _ in 0..20 {
            let logt = crate::series::logt_symbol();
            let _ = crate::series::with_engine(&x, &Expr::zero(), false, &logt, |eng| eng.ser(&e, n)).unwrap();
        }
        let e2 = e.clone();
        let t3 = std::time::Instant::now();
        let mut eng_s = None;
        for _ in 0..20 {
            let logt = crate::series::logt_symbol();
            let mut eng = crate::series::engine::Engine::new(&x, crate::series::engine::var_series(&Expr::zero(), false, 1), 1, &logt);
            let c = eng.ser(&p("cos(x)"), n).unwrap();
            eng_s = Some(c);
        }
        let tc = t3.elapsed() / 20;
        let c = eng_s.unwrap();
        let t4 = std::time::Instant::now();
        for _ in 0..20 {
            let _ = c.inv(n).unwrap();
        }
        println!("{src}: total {:?}, engine {:?}, cos {:?}, inv {:?} {}", dt, t2.elapsed() / 20, tc, t4.elapsed() / 20, e2);
    }
}
