use super::*;
use crate::parser::parse_str;

fn p(s: &str) -> Expr {
    parse_str(s).unwrap()
}

fn s(e: &Expr) -> String {
    e.to_string()
}

#[test]
fn set_evaluation_and_printing() {
    let cases = [
        ("Interval(0, 1)", "Interval(0, 1)"),
        ("Interval(0, 1, True, False)", "Interval.Lopen(0, 1)"),
        ("Interval(1, 0)", "EmptySet"),
        ("Interval(1, 1)", "{1}"),
        ("Interval(-oo, 0)", "Interval(-oo, 0)"),
        ("FiniteSet(3, 1, 2, 1)", "{1, 2, 3}"),
        ("Union(Interval(0, 2), Interval(1, 3))", "Interval(0, 3)"),
        ("Union(Interval(0, 1), FiniteSet(2))", "Union({2}, Interval(0, 1))"),
        ("Union(Interval(0, 1, True, True), FiniteSet(0, 1))", "Interval(0, 1)"),
        ("Intersection(Interval(0, 2), Interval(1, 3))", "Interval(1, 2)"),
        ("Intersection(Interval(0, 1), Interval(1, 2))", "{1}"),
        ("Complement(Interval(0, 3), Interval(1, 2))", "Union(Interval.Ropen(0, 1), Interval.Lopen(2, 3))"),
        ("Complement(Reals, FiniteSet(0))", "Union(Interval.open(-oo, 0), Interval.open(0, oo))"),
        ("Intersection(Integers, Interval(0, 5))", "Range(0, 6, 1)"),
        ("Union(Integers, Reals)", "Reals"),
        ("ConditionSet(x, x > 0, FiniteSet(-1, 1, y))", "Union(ConditionSet(x, x > 0, {y}), {1})"),
    ];
    for (i, o) in cases {
        assert_eq!(s(&p(i)), o, "{i}");
    }
}

#[test]
fn set_properties() {
    let x = Expr::symbol("x");
    assert_eq!(s(&contains(&p("Interval(0, 1)"), &x)), "(0 <= x) & (x <= 1)");
    assert_eq!(s(&contains(&p("FiniteSet(1, 2)"), &x)), "Eq(1, x) | Eq(2, x)");
    assert_eq!(s(&as_relational(&p("Union(Interval(0, 1), Interval(2, 3))"), &x)), "((0 <= x) & (x <= 1)) | ((2 <= x) & (x <= 3))");
    assert_eq!(s(&as_relational(&p("Complement(Reals, FiniteSet(0))"), &x)), "(x > -oo) & (x < oo) & Ne(x, 0)");
    assert_eq!(s(&boundary(&p("Union(Interval(0, 1, False, True), Interval(1, 2, True, False))")).unwrap()), "{0, 1, 2}");
    assert_eq!(s(&measure(&p("Union(Interval(0, 1), Interval(2, 4))")).unwrap()), "3");
    assert_eq!(s(&inf(&p("FiniteSet(x, 2)")).unwrap()), "Min(2, x)");
    assert_eq!(is_subset(&p("Interval(0, 1)"), &p("Interval(0, 1, True, False)")), Some(false));
}

#[test]
fn solving() {
    let x = Expr::symbol("x");
    let ss = |e: &str, d: Expr| s(&solveset(&p(e), &x, &d).unwrap());
    assert_eq!(ss("x**2 - 1", complexes()), "{-1, 1}");
    assert_eq!(ss("x**2 + 1", reals()), "EmptySet");
    assert_eq!(ss("sin(x) - 1", reals()), "ImageSet(Lambda(_n, 2*_n*pi + pi/2), Integers)");
    assert_eq!(ss("tan(2*x) - 1", reals()), "ImageSet(Lambda(_n, _n*pi/2 + pi/8), Integers)");
    assert_eq!(ss("exp(x) - 1", complexes()), "ImageSet(Lambda(_n, 2*_n*I*pi), Integers)");
    assert_eq!(ss("sin(x) - x", reals()), "ConditionSet(x, Eq(-x + sin(x), 0), Reals)");
    assert_eq!(ss("x**2 > 1", reals()), "Union(Interval.open(-oo, -1), Interval.open(1, oo))");
    let sui = |e: &str| s(&solve_univariate_inequality(&p(e), &x, true, &reals()).unwrap());
    assert_eq!(sui("x**2 - 3*x + 2 < 0"), "(1 < x) & (x < 2)");
    assert_eq!(sui("(x - 1)/(x + 1) >= 0"), "((1 <= x) & (x < oo)) | ((-oo < x) & (x < -1))");
    assert_eq!(sui("x**2 > 0"), "(x > -oo) & (x < oo) & Ne(x, 0)");
    assert_eq!(sui("sqrt(x) < 1"), "(0 <= x) & (x < 1)");
    assert_eq!(s(&reduce_inequalities(&[p("x**2 >= 0")], Some(&x)).unwrap()), "True");
    assert_eq!(s(&reduce_inequalities(&[p("x > 1"), p("x < 3")], Some(&x)).unwrap()), "(1 < x) & (x < 3)");
}
