use symox::parser::parse_str;
use symox::poly::api;
fn main() {
    let which = std::env::args().nth(1).unwrap_or_default();
    let n: usize = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(20);
    match which.as_str() {
        "gcd" => {
            let g = symox::expand::expand(&parse_str("(x + y + z + 1)**6").unwrap());
            let a = symox::expand::expand(&symox::ops::mul(&[g.clone(), parse_str("(x - y + 2*z + 3)**4").unwrap()]));
            let b = symox::expand::expand(&symox::ops::mul(&[g.clone(), parse_str("(x*y - z + 2)**3").unwrap()]));
            for _ in 0..n { std::hint::black_box(api::gcd(&a, &b, &[]).unwrap()); }
        }
        "wilk" => {
            let mut s = String::from("1");
            for i in 1..=20 { s.push_str(&format!("*(x - {i})")); }
            let e = symox::expand::expand(&parse_str(&s).unwrap());
            for _ in 0..n { std::hint::black_box(api::factor(&e).unwrap()); }
        }
        _ => {}
    }
}
