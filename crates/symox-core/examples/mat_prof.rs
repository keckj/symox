use std::time::Instant;
use symox::matrix::Matrix;
use symox::prelude::*;

fn main() {
    let mut seed: u64 = 12345;
    let mut rnd = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 33) % 19) as i64 - 9
    };
    let n = 50;
    let m = Matrix::from_fn(n, n, |_, _| Expr::int(rnd()));
    let t = Instant::now();
    let r = symox::matrix::linalg::RatMat::from_matrix(&m).unwrap();
    println!("convert {:?}", t.elapsed());
    let t = Instant::now();
    let (rr, _) = r.rref();
    println!("rref kernel {:?}", t.elapsed());
    let t = Instant::now();
    let _ = rr.to_matrix();
    println!("to_matrix {:?}", t.elapsed());
    let t = Instant::now();
    let _ = m.rref();
    println!("rref total {:?}", t.elapsed());
    let n = 30;
    let mut rnd2 = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((seed >> 33) % 19) as i64 - 9
    };
    let m = Matrix::from_fn(n, n, |_, _| { let a = rnd2(); let b = (rnd2().abs() % 9) + 1; Expr::rational(a, b) });
    let t = Instant::now();
    let r = symox::matrix::linalg::RatMat::from_matrix(&m).unwrap();
    let (a, _s) = r.to_int_rows();
    println!("to_int_rows {:?}", t.elapsed());
    let t = Instant::now();
    let _ = symox::matrix::modular::inv_int(&a, n);
    println!("inv_int {:?}", t.elapsed());
    let t = Instant::now();
    let inv = r.inv().unwrap();
    println!("inv kernel {:?}", t.elapsed());
    let t = Instant::now();
    let _ = r.inv_ff().unwrap();
    println!("inv ff {:?}", t.elapsed());
    let t = Instant::now();
    let _ = inv.to_matrix();
    println!("inv to_matrix {:?}", t.elapsed());
}

#[allow(dead_code)]
fn unused() {}
