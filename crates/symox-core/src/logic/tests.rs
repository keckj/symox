use super::sat::{mk_lit, solve_dimacs, Solver};
use super::*;
use crate::parser::parse_str;

fn p(s: &str) -> Expr {
    parse_str(s).unwrap()
}

fn s(e: &Expr) -> String {
    e.to_string()
}

#[test]
fn evaluation_and_printing() {
    let cases = [
        ("And(x, y)", "x & y"),
        ("And(y, x)", "x & y"),
        ("Or(x, And(y, z))", "x | (y & z)"),
        ("And(Or(x, y), z)", "z & (x | y)"),
        ("Not(Not(x))", "x"),
        ("Not(And(x, y))", "~(x & y)"),
        ("Xor(x, x)", "False"),
        ("Xor(x, True)", "~x"),
        ("Xor(x, y, True)", "~(x ^ y)"),
        ("Implies(True, x)", "x"),
        ("Implies(x, x)", "True"),
        ("Equivalent(x, y, True)", "x & y"),
        ("Equivalent(x, False)", "~x"),
        ("ITE(x, True, False)", "x"),
        ("ITE(x, False, True)", "~x"),
        ("And(x < 1, y > 2)", "(y > 2) & (x < 1)"),
        ("Not(x < 1)", "x >= 1"),
        ("Not(Eq(x, 1))", "Ne(x, 1)"),
        ("And(x < 1, x >= 1)", "False"),
        ("Or(x < 1, x >= 1)", "True"),
        ("And(x > -oo, x < 1)", "(x > -oo) & (x < 1)"),
        ("And(x < 1, -oo < x, y)", "(-oo < x) & y & (x < 1)"),
        ("And(Eq(x + y, 1), Or(a, b))", "(a | b) & Eq(x + y, 1)"),
        ("And(Implies(x, y), z)", "z & (Implies(x, y))"),
        ("Not(Implies(x, y))", "~(Implies(x, y))"),
        ("Or(Xor(x, y), z)", "z | (x ^ y)"),
        ("x & y | ~z", "~z | (x & y)"),
    ];
    for (i, o) in cases {
        assert_eq!(s(&p(i)), o, "{i}");
    }
}

#[test]
fn normal_forms() {
    let cnf = |e: &str| s(&to_cnf(&p(e), false, false).unwrap());
    let dnf = |e: &str| s(&to_dnf(&p(e), false, false).unwrap());
    assert_eq!(cnf("Or(And(a, b), c)"), "(a | c) & (b | c)");
    assert_eq!(cnf("Implies(a, b)"), "b | ~a");
    assert_eq!(cnf("Equivalent(a, b)"), "(a | ~b) & (b | ~a)");
    assert_eq!(cnf("Xor(a, b, c)"), "(a | b | c) & (a | ~b | ~c) & (b | ~a | ~c) & (c | ~a | ~b)");
    assert_eq!(cnf("ITE(a, b, c)"), "(a | c) & (b | ~a)");
    assert_eq!(dnf("And(Or(a, b), c)"), "(a & c) | (b & c)");
    assert_eq!(s(&to_nnf(&p("Not(ITE(a, b, c))"), true)), "(a | ~c) & (~a | ~b)");
    assert!(is_cnf(&p("And(a, Or(b, c))")));
    assert!(!is_dnf(&p("And(a, Or(b, c))")));
}

#[test]
fn minimization() {
    let sl = |e: &str| s(&simplify_logic(&p(e), None, true, false));
    assert_eq!(sl("Or(And(a, b), And(a, Not(b)))"), "a");
    assert_eq!(sl("Or(And(x, y), And(Not(x), z), And(y, z))"), "(x & y) | (z & ~x)");
    assert_eq!(sl("Xor(a, b)"), "(a & ~b) | (b & ~a)");
    assert_eq!(sl("Or(x < 1, x >= 1)"), "True");
    let v: Vec<Expr> = ["a", "b", "c", "d"].iter().map(|n| Expr::symbol(n)).collect();
    assert_eq!(s(&sop_form(&v, &[1, 3, 7, 11, 15], &[0, 2, 5])), "(c & d) | (~a & ~b)");
}

/// Brute force satisfiability of a small CNF.
fn brute(n: usize, clauses: &[Vec<i32>]) -> bool {
    (0..1u32 << n).any(|m| {
        clauses
            .iter()
            .all(|c| c.iter().any(|&l| ((m >> (l.unsigned_abs() - 1)) & 1 == 1) == (l > 0)))
    })
}

#[test]
fn sat_random_vs_brute_force() {
    let mut seed = 0x2545F4914F6CDD1Du64;
    let mut rnd = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for _ in 0..300 {
        let n = 3 + (rnd() % 10) as usize;
        let m = (rnd() % (5 * n as u64)) as usize + 1;
        let clauses: Vec<Vec<i32>> = (0..m)
            .map(|_| {
                (0..3)
                    .map(|_| {
                        let v = (rnd() % n as u64) as i32 + 1;
                        if rnd() % 2 == 0 {
                            v
                        } else {
                            -v
                        }
                    })
                    .collect()
            })
            .collect();
        let r = solve_dimacs(n, &clauses);
        assert_eq!(r.is_some(), brute(n, &clauses));
        if let Some(model) = r {
            for c in &clauses {
                assert!(c.iter().any(|&l| model[(l.unsigned_abs() - 1) as usize] == (l > 0)));
            }
        }
    }
}

#[test]
fn sat_pigeonhole_unsat() {
    // 7 pigeons, 6 holes
    let (p, h) = (7u32, 6u32);
    let var = |i: u32, j: u32| i * h + j;
    let mut s = Solver::new();
    for _ in 0..p * h {
        s.new_var();
    }
    for i in 0..p {
        let c: Vec<u32> = (0..h).map(|j| mk_lit(var(i, j), false)).collect();
        s.add_clause(&c);
    }
    for j in 0..h {
        for a in 0..p {
            for b in a + 1..p {
                s.add_clause(&[mk_lit(var(a, j), true), mk_lit(var(b, j), true)]);
            }
        }
    }
    assert!(!s.solve());
}

#[test]
fn satisfiable_models() {
    let e = p("And(Or(x, y), Not(x))");
    let m = satisfiable(&e).unwrap();
    assert_eq!(m, vec![(Expr::symbol("x"), false), (Expr::symbol("y"), true)]);
    assert!(satisfiable(&p("And(x, Not(x))")).is_none());
    assert_eq!(satisfiable_all(&p("Xor(a, b, c)")).count(), 4);
    assert_eq!(satisfiable_all(&p("Or(a, b, c, d)")).count(), 15);
    // non-CNF input goes through the Tseitin encoding
    assert_eq!(satisfiable_all(&p("Equivalent(a, ITE(b, c, d))")).count(), 8);
}

#[test]
fn fast_ordering_agrees_with_sort_key() {
    let exprs: Vec<Expr> = [
        "a", "b", "aa", "B", "Not(a)", "Not(b)", "Or(a, b)", "Or(a, Not(b))", "Or(Not(a), b)", "And(a, b)",
        "And(Not(a), b)", "Or(a, b, c)", "And(a, Not(c))", "Or(b, c)", "Not(B)", "Or(Not(a), Not(b))",
    ]
    .iter()
    .map(|t| p(t))
    .collect();
    for a in &exprs {
        for b in &exprs {
            let (la, lb) = (super::light(a).unwrap(), super::light(b).unwrap());
            let fast = super::light_cmp(&la, &lb);
            let slow = ordered_key(a).cmp(&ordered_key(b));
            assert_eq!(fast, slow, "{a} vs {b}");
        }
    }
}
#[test]
fn complete_boolean_piecewise_conditions_are_logical_nodes() {
    let condition = crate::Expr::symbol("condition");
    let yes = crate::Expr::symbol("yes");
    let no = crate::Expr::symbol("no");
    let piecewise = crate::relational::piecewise(vec![
        (yes.clone(), condition.clone()), (no.clone(), crate::Expr::boolean(true)),
    ]);
    assert!(!super::is_boolean_category(&piecewise));
    let expected = super::ite(&condition, &yes, &no);
    assert_eq!(super::boolean_piecewise_to_ite(&piecewise), Some(expected.clone()));
    let outer = crate::relational::piecewise(vec![
        (crate::Expr::one(), piecewise), (crate::Expr::zero(), crate::Expr::boolean(true)),
    ]);
    let crate::Kind::Piecewise(branches) = outer.kind() else { panic!("held numeric branches expected"); };
    assert_eq!(branches[0].1, expected);
    let incomplete = crate::relational::piecewise(vec![(yes, condition)]);
    assert_eq!(super::boolean_piecewise_to_ite(&incomplete), None);
}

#[test]
fn nodes2_counts_the_expanded_tree_without_walking_it() {
    // SymPy's `_nodes` counts tree occurrences; the count of a shared node is
    // reused instead of recomputed.
    let x = Expr::symbol("x");
    let pair = |a: &Expr| crate::functions::undefined("pair", vec![a.clone(), a.clone()]);
    assert_eq!(nodes2(&pair(&x)), 2 + 2 * 2);
    assert_eq!(nodes2(&pair(&pair(&x))), 2 + 2 * (2 + 2 * 2));
    let mut value = x;
    for depth in 1..=100u32 {
        value = pair(&value);
        // tree size 2^(depth+1) - 1 nodes, each counted twice
        let expected = ((1u128 << (depth + 2)) - 2).min(u64::MAX as u128) as u64;
        assert_eq!(nodes2(&value), expected, "depth {depth}");
    }
}
