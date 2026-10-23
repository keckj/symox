//! Printing support for `Order` terms (SymPy compatible).  The term order
//! of an `Add` holding `Order` terms is `printing::order::ordered_terms`.

use crate::expr::Expr;
use crate::printer::to_str;

/// `str` of an `Order` node given its arguments.
pub fn order_str(args: &[Expr]) -> String {
    let mut s = String::from("O(");
    s.push_str(&to_str(&args[0]));
    let pairs: Vec<&[Expr]> = args[1..].chunks(2).collect();
    let all_zero = pairs.iter().all(|p| p.len() == 2 && p[1].is_zero());
    if all_zero {
        if pairs.len() > 1 {
            for p in &pairs {
                s.push_str(", ");
                s.push_str(&to_str(&p[0]));
            }
        }
    } else {
        for p in &pairs {
            s.push_str(", (");
            s.push_str(&to_str(&p[0]));
            if p.len() == 2 {
                s.push_str(", ");
                s.push_str(&to_str(&p[1]));
            }
            s.push(')');
        }
    }
    s.push(')');
    s
}
