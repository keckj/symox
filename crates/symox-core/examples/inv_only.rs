use symox::matrix::Matrix;
use symox::prelude::*;

fn main() {
    let which = std::env::args().nth(1).unwrap_or_default();
    let x = Expr::symbol("x");
    let y = Expr::symbol("y");
    let m = Matrix::from_fn(4, 4, |i, j| {
        let e = ((i + j) % 3) as i64;
        add(&[x.powi(e), mul(&[y.clone(), Expr::int(i as i64 - j as i64)]), Expr::int((((i + 1) * (j + 2)) % 3) as i64)])
    });
    let syms: Vec<Expr> = (0..9).map(|i| Expr::symbol(&format!("a{i}"))).collect();
    let g = Matrix::new(3, 3, syms).unwrap();
    for _ in 0..20 {
        if which != "g" { let _ = m.inv().unwrap(); }
        if which != "p" { let _ = g.inv().unwrap(); }
    }
}
