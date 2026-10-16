//! Unit tests for the matrix module.

use super::*;
use crate::parser::parse_str;

fn p(s: &str) -> Expr {
    parse_str(s).unwrap()
}

fn m(rows: &[&[&str]]) -> Matrix {
    Matrix::from_rows(rows.iter().map(|r| r.iter().map(|s| p(s)).collect()).collect()).unwrap()
}

#[test]
fn printing() {
    let a = m(&[&["1", "2"], &["3", "4"]]);
    assert_eq!(a.to_str(), "Matrix([[1, 2], [3, 4]])");
    assert_eq!(a.to_repr(), "Matrix([\n[1, 2],\n[3, 4]])");
    assert_eq!(Matrix::zeros(0, 3).to_str(), "Matrix(0, 3, [])");
    let c = m(&[&["-1"], &["1"]]);
    assert_eq!(c.to_repr(), "Matrix([\n[-1],\n[ 1]])");
}

#[test]
fn rational_linalg() {
    let a = m(&[&["1", "2"], &["3", "4"]]);
    assert_eq!(a.det().unwrap().to_string(), "-2");
    assert_eq!(a.inv().unwrap().to_str(), "Matrix([[-2, 1], [3/2, -1/2]])");
    let b = m(&[&["1", "2", "3"], &["4", "5", "6"], &["7", "8", "9"]]);
    assert_eq!(b.det().unwrap().to_string(), "0");
    let (r, piv) = b.rref();
    assert_eq!(r.to_str(), "Matrix([[1, 0, -1], [0, 1, 2], [0, 0, 0]])");
    assert_eq!(piv, vec![0, 1]);
    assert_eq!(b.rank(), 2);
    let ns = b.nullspace();
    assert_eq!(ns.len(), 1);
    assert_eq!(ns[0].to_str(), "Matrix([[1], [-2], [1]])");
    let c = Matrix::from_fn(5, 5, |i, j| Expr::int(((i * 7 + j * 3) % 5) as i64 + if i == j { 3 } else { 0 }));
    let ci = c.inv().unwrap();
    assert!(c.matmul(&ci).unwrap() == Matrix::eye(5));
}

#[test]
fn symbolic_det_inv() {
    let a = m(&[&["a", "b"], &["c", "d"]]);
    assert_eq!(a.det().unwrap().to_string(), "a*d - b*c");
    assert_eq!(
        a.inv().unwrap().to_str(),
        "Matrix([[d/(a*d - b*c), -b/(a*d - b*c)], [-c/(a*d - b*c), a/(a*d - b*c)]])"
    );
    let b = m(&[&["x", "1"], &["2", "x"]]);
    assert_eq!(b.inv().unwrap().to_str(), "Matrix([[x/(x**2 - 2), -1/(x**2 - 2)], [-2/(x**2 - 2), x/(x**2 - 2)]])");
    let c = m(&[&["1", "x"], &["1", "y"]]);
    assert_eq!(c.inv().unwrap().to_str(), "Matrix([[-y/(x - y), x/(x - y)], [1/(x - y), -1/(x - y)]])");
    let d = Matrix::from_fn(4, 4, |i, j| {
        let x = Expr::symbol("x");
        crate::ops::add2(&x.powi((i + j) as i64), &Expr::int(i as i64))
    });
    assert_eq!(d.det().unwrap().to_string(), "0");
    let e = m(&[&["x", "1", "0", "0"], &["1", "x", "0", "0"], &["0", "0", "x", "2"], &["0", "0", "3", "x"]]);
    assert_eq!(e.det().unwrap().to_string(), "(x**2 - 6)*(x**2 - 1)");
}

#[test]
fn charpoly_eigen() {
    let a = m(&[&["1", "2"], &["3", "4"]]);
    assert_eq!(a.charpoly("lambda").unwrap().as_expr().to_string(), "lambda**2 - 5*lambda - 2");
    let ev = a.eigenvals().unwrap();
    assert_eq!(ev[0].0.to_string(), "5/2 - sqrt(33)/2");
    let b = m(&[&["2", "1"], &["1", "2"]]);
    let vs = b.eigenvects().unwrap();
    assert_eq!(vs[0].0.to_string(), "1");
    assert_eq!(vs[0].2[0].to_str(), "Matrix([[-1], [1]])");
    assert_eq!(vs[1].2[0].to_str(), "Matrix([[1], [1]])");
    let vs = a.eigenvects().unwrap();
    assert_eq!(vs[0].2[0].to_str(), "Matrix([[-sqrt(33)/6 - 1/2], [1]])");
}

#[test]
fn rref_symbolic() {
    let a = m(&[&["x", "1", "2"], &["1", "y", "3"]]);
    let (r, piv) = a.rref();
    assert_eq!(r.to_str(), "Matrix([[1, 0, (2*y - 3)/(x*y - 1)], [0, 1, (3*x - 2)/(x*y - 1)]])");
    assert_eq!(piv, vec![0, 1]);
    let b = m(&[&["x", "x + 1"]]);
    assert_eq!(b.rref().0.to_str(), "Matrix([[1, 1 + 1/x]])");
}

#[test]
fn powers_and_products() {
    let a = m(&[&["x + 1", "1"], &["1", "x"]]);
    assert_eq!(a.matmul(&a).unwrap().to_str(), "Matrix([[(x + 1)**2 + 1, 2*x + 1], [2*x + 1, x**2 + 1]])");
    assert_eq!(a.powi(2).unwrap().to_str(), "Matrix([[x**2 + 2*x + 2, 2*x + 1], [2*x + 1, x**2 + 1]])");
    let b = m(&[&["1", "2"], &["3", "4"]]);
    assert_eq!(b.powi(-2).unwrap().to_str(), "Matrix([[11/2, -5/2], [-15/4, 7/4]])");
}

#[test]
fn lu() {
    let a = m(&[&["a", "b"], &["c", "d"]]);
    let (l, u, p) = a.lu_decomposition().unwrap();
    assert_eq!(l.to_str(), "Matrix([[1, 0], [c/a, 1]])");
    assert_eq!(u.to_str(), "Matrix([[a, b], [0, d - b*c/a]])");
    assert!(p.is_empty());
    let b = Matrix::column(vec![p_("1"), p_("2")]);
    assert_eq!(a.solve(&b).unwrap().to_str(), "Matrix([[(-2*b + d)/(a*d - b*c)], [(2*a - c)/(a*d - b*c)]])");
}

fn p_(s: &str) -> Expr {
    p(s)
}

fn lcg(seed: &mut u64) -> i64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    (*seed >> 33) as i64
}

#[test]
fn modular_matches_fraction_free() {
    use super::linalg::RatMat;
    let mut seed = 7u64;
    for trial in 0..40 {
        let rows = 6 + (lcg(&mut seed) % 7) as usize;
        let cols = 6 + (lcg(&mut seed) % 9) as usize;
        let rank_cap = if trial % 3 == 0 { 1 + (lcg(&mut seed) % 5) as usize } else { 100 };
        // low rank via product when requested
        let m = if rank_cap < 100 {
            let a = Matrix::from_fn(rows, rank_cap, |_, _| Expr::int(lcg(&mut seed) % 11 - 5));
            let b = Matrix::from_fn(rank_cap, cols, |_, _| Expr::int(lcg(&mut seed) % 11 - 5));
            a.matmul(&b).unwrap()
        } else if trial % 3 == 1 {
            Matrix::from_fn(rows, cols, |_, _| Expr::rational(lcg(&mut seed) % 21 - 10, 1 + lcg(&mut seed) % 7))
        } else {
            Matrix::from_fn(rows, cols, |_, _| Expr::int(lcg(&mut seed) % 2001 - 1000))
        };
        let r = RatMat::from_matrix(&m).unwrap();
        let (a, pa) = r.rref();
        let (b, pb) = r.rref_ff();
        assert_eq!(pa, pb, "pivots differ in trial {trial}");
        assert!(a == b, "rref differs in trial {trial}");
        if rows == cols {
            let x = r.inv();
            let y = r.inv_ff();
            assert!(x == y, "inverse differs in trial {trial}");
        }
        let sq = Matrix::from_fn(rows, rows, |i, j| m.get(i, j % cols).clone());
        let rs = RatMat::from_matrix(&sq).unwrap();
        assert!(rs.inv() == rs.inv_ff(), "square inverse differs in trial {trial}");
    }
}

#[test]
fn nonfinite_matrix_products_keep_zero_terms() {
    let zero = Matrix::zeros(1, 1);
    for value in [Expr::oo(), Expr::neg_oo(), Expr::zoo(), Expr::nan()] {
        let special = Matrix::new(1, 1, vec![value]).unwrap();
        assert!(zero.matmul(&special).unwrap().get(0, 0).is_constant(crate::expr::Constant::NaN));
        assert!(special.matmul(&zero).unwrap().get(0, 0).is_constant(crate::expr::Constant::NaN));
        assert!(zero.matmul_plain(&special).get(0, 0).is_constant(crate::expr::Constant::NaN));
    }
    let held = Expr::from_unevaluated_kind(crate::expr::Kind::Mul(
        crate::number::Number::ONE, vec![(Expr::oo(), Expr::one()), (Expr::symbol("x"), Expr::one())],
    ));
    let special = Matrix::new(1, 1, vec![held]).unwrap();
    assert!(zero.matmul(&special).unwrap().get(0, 0).is_constant(crate::expr::Constant::NaN));
}
