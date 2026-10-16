use num_bigint::BigUint;
use num_integer::Integer;
use std::time::Instant;

fn main() {
    let mut seed = 3u64;
    let mut next = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        seed
    };
    let mut pairs = Vec::new();
    for _ in 0..1000 {
        let mut a = BigUint::from(1u32);
        let mut b = BigUint::from(1u32);
        for _ in 0..8 {
            a = a * BigUint::from(next()) + BigUint::from(next());
            b = b * BigUint::from(next()) + BigUint::from(next());
        }
        pairs.push((a, b));
    }
    for _ in 0..3 {
        let t = Instant::now();
        let mut s = 0u64;
        for (a, b) in &pairs {
            s += a.gcd(b).bits();
        }
        println!("binary {:?} {s}", t.elapsed());
        let t = Instant::now();
        let mut s = 0u64;
        for (a, b) in &pairs {
            s += symox::matrix::int::lehmer_gcd(a.clone(), b.clone()).bits();
        }
        println!("lehmer {:?} {s}", t.elapsed());
    }
}
