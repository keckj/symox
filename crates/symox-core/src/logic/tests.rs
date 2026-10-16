use super::*;
use crate::parser::parse_str;

fn p(s: &str) -> Expr {
    parse_str(s).unwrap()
}

fn s(e: &Expr) -> String {
    e.to_string()
}
