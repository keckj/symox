use std::time::Instant;
use symox::matrix::Matrix;
use symox::prelude::*;

fn main() {
    let x = Expr::symbol("x");
    let y = Expr::symbol("y");
    let m = Matrix::from_fn(4, 4, |i, j| {
        let e = ((i + j) % 3) as i64;
        add(&[x.powi(e), mul(&[y.clone(), Expr::int(i as i64 - j as i64)]), Expr::int((((i + 1) * (j + 2)) % 3) as i64)])
    });
    for _ in 0..3 {
        let t = Instant::now();
        let inv = m.inv().unwrap();
        println!("inv {:?}", t.elapsed());
        let _ = inv;
    }
    use symox::matrix::domain::{matrix_domain, Domain, PolyMat};
    let t = Instant::now();
    for _ in 0..10 { let _ = matrix_domain(&m); }
    println!("domain x10 {:?}", t.elapsed());
    let Domain::Poly { ctx, .. } = matrix_domain(&m) else { panic!() };
    let t = Instant::now();
    let (pm, _) = PolyMat::from_matrix(&m, &ctx).unwrap();
    println!("to poly {:?}", t.elapsed());
    let n = 4;
    let one = symox::matrix::mpoly::MPoly::one(ctx.nvars());
    let zero = symox::matrix::mpoly::MPoly::zero(ctx.nvars());
    let mut aug = Vec::new();
    for i in 0..n { for j in 0..n { aug.push(pm.get(i, j).clone()); } for j in 0..n { aug.push(if i == j { one.clone() } else { zero.clone() }); } }
    let t = Instant::now();
    let (den, _) = symox::matrix::ff::rref_den(&mut aug, n, 2 * n, &one);
    println!("rref_den {:?} den terms {}", t.elapsed(), den.terms.len());
    println!("sizes {:?}", aug.iter().map(|p| p.terms.len()).collect::<Vec<_>>());
    let a = pm.get(1, 1).clone(); let b = pm.get(2, 3).clone();
    let t = Instant::now();
    let mut c = a.clone();
    for _ in 0..1000 { c = a.mul(&b); }
    println!("1000 muls {:?} {} x {} -> {}", t.elapsed(), a.terms.len(), b.terms.len(), c.terms.len());
    let t = Instant::now();
    for _ in 0..1000 { let _ = c.divexact(&b); }
    println!("1000 divs {:?}", t.elapsed());
    let t = Instant::now();
    for i in 0..n { for j in 0..n { let _ = symox::matrix::ratfunc::cancel_polys(&aug[i * 2 * n + n + j], &den); } }
    println!("cancel {:?}", t.elapsed());
    let t = Instant::now();
    let d = m.det().unwrap();
    println!("det {:?} {}", t.elapsed(), d);
    let syms: Vec<Expr> = (0..9).map(|i| Expr::symbol(&format!("a{i}"))).collect();
    let g = Matrix::new(3, 3, syms).unwrap();
    for _ in 0..3 {
        let t = Instant::now();
        let _ = g.inv().unwrap();
        println!("inv3 generic {:?}", t.elapsed());
    }
}

#[allow(dead_code)]
fn parts() {}
