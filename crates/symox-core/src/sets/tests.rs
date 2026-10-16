use super::*;
use crate::parser::parse_str;

fn p(s: &str) -> Expr {
    parse_str(s).unwrap()
}
