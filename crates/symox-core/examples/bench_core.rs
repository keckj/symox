//! Micro benchmarks of the core constructors: `cargo run --release --example bench_core`
use std::time::Instant;

use symox::prelude::*;

fn time<F: FnMut() -> usize>(name: &str, mut f: F) {
    let _ = f();
    let t = Instant::now();
    let mut n = 0;
    let mut iters = 0;
    while t.elapsed().as_secs_f64() < 0.5 {
        n += f();
        iters += 1;
    }
    let dt = t.elapsed().as_secs_f64() / iters as f64;
    println!("{name:40} {:10.3} ms  ({n} checksum)", dt * 1e3);
}

fn main() {
    let syms: Vec<Expr> = (0..50).map(|i| Expr::symbol(&format!("x{i}"))).collect();
    let x = Expr::symbol("x");
    let y = Expr::symbol("y");
    time("10k binary products of symbols", || {
        let mut c = 0;
        for i in 0..10_000 {
            let e = &syms[i % 50] * &syms[(i * 7 + 3) % 50];
            c += e.args().len();
        }
        c
    });
    time("10k num*sym and sym*num", || {
        let mut c = 0;
        for i in 0..10_000 {
            let e = (i as i64) * &syms[i % 50];
            c += e.args().len();
        }
        c
    });
    time("10k x**i * x**j", || {
        let mut c = 0;
        for i in 0..10_000i64 {
            let e = x.powi(i % 7 + 1) * x.powi(i % 5 + 1);
            c += e.args().len();
        }
        c
    });
    time("sum of 2000 products x_i*x_j", || {
        let mut s = Expr::zero();
        for i in 0..2000 {
            s = s + &syms[i % 50] * &syms[(i / 50) % 50];
        }
        s.args().len()
    });
    time("diff of (x+y)**20 expanded, 10x", || {
        let e = expand(&(&x + &y).powi(20));
        let mut c = 0;
        for _ in 0..10 {
            c += diff(&e, &x).args().len();
        }
        c
    });
    time("expand((x+y+1)**30)", || expand(&(&x + &y + 1).powi(30)).args().len());
    time("subs x->y+1 in (x+y)**20 expanded", || {
        let e = expand(&(&x + &y).powi(20));
        subs(&e, &[(x.clone(), &y + 1)]).args().len()
    });
}
