//! Sparse Boolean polynomials over GF(2), with idempotent variables.
//!
//! Symbolic conversion does not allocate a truth table: its cost follows
//! the actual number of output monomials rather than 2^number_of_atoms.

use rustc_hash::FxHashSet;
use super::{and, as_logic, is_literal, mk, ordered, LogicOp};
use crate::expr::{Expr, Kind};

type Polynomial = FxHashSet<Vec<Expr>>;

fn toggle(p: &mut Polynomial, m: Vec<Expr>) {
    if !p.insert(m.clone()) { p.remove(&m); }
}

fn constant(b: bool) -> Polynomial {
    if b { [vec![]].into_iter().collect() } else { Polynomial::default() }
}

fn atom(e: &Expr) -> Polynomial {
    match e.kind() {
        Kind::Bool(b) => constant(*b),
        _ => [vec![e.clone()]].into_iter().collect(),
    }
}

fn sum(mut p: Polynomial, q: Polynomial) -> Polynomial {
    for m in q { toggle(&mut p, m); }
    p
}

fn product(p: Polynomial, q: Polynomial) -> Polynomial {
    let mut out = Polynomial::default();
    for a in &p {
        for b in &q {
            let mut m = a.clone();
            m.extend(b.iter().cloned());
            ordered(&mut m);
            m.dedup();
            toggle(&mut out, m);
        }
    }
    out
}

fn complement(p: Polynomial) -> Polynomial { sum(p, constant(true)) }

fn convert(e: &Expr, deep: bool) -> Polynomial {
    let Some((op, args)) = as_logic(e) else { return atom(e) };
    let sub = |e: &Expr| if deep { convert(e, true) } else { atom(e) };
    match op {
        LogicOp::Not => complement(sub(&args[0])),
        LogicOp::Xor => args.iter().fold(constant(false), |p, a| sum(p, sub(a))),
        LogicOp::And => args.iter().fold(constant(true), |p, a| product(p, sub(a))),
        // OR(a_i) = 1 + product(1+a_i) over GF(2).
        LogicOp::Or => complement(args.iter().fold(constant(true), |p, a| product(p, complement(sub(a))))),
        LogicOp::Implies => {
            let a = sub(&args[0]);
            complement(sum(a.clone(), product(a, sub(&args[1]))))
        }
        LogicOp::Equivalent => {
            let yes = args.iter().fold(constant(true), |p, a| product(p, sub(a)));
            let no = args.iter().fold(constant(true), |p, a| product(p, complement(sub(a))));
            sum(yes, no)
        }
        // if c then t else f = f + c*(t+f).
        LogicOp::ITE => {
            let f = sub(&args[2]);
            sum(f.clone(), product(sub(&args[0]), sum(sub(&args[1]), f)))
        }
    }
}

fn expression(p: Polynomial) -> Expr {
    let mut terms = Vec::with_capacity(p.len());
    // Rebuilding And can itself flatten opaque conjunctions (deep=False),
    // so combine identical terms once more after reconstruction.
    let mut seen = FxHashSet::default();
    for monomial in p {
        let term = and(monomial);
        if matches!(term.kind(), Kind::Bool(false)) { continue; }
        if !seen.insert(term.clone()) { seen.remove(&term); }
    }
    terms.extend(seen);
    ordered(&mut terms);
    match terms.len() {
        0 => Expr::boolean(false),
        1 => terms.pop().unwrap(),
        // Unlike ordinary Xor construction, keep the constant 1 in ANF;
        // turning it into Not would leave the requested normal form.
        _ => mk(LogicOp::Xor, terms),
    }
}

/// Is this a SymPy-compatible algebraic normal form? Relational literals
/// may appear alone or as XOR terms; monomials contain Symbols only.
pub fn is_anf(e: &Expr) -> bool {
    if is_literal(e) && !matches!(as_logic(e), Some((LogicOp::Not, _))) { return true; }
    let monomial = |a: &Expr| matches!(as_logic(a), Some((LogicOp::And, terms)) if terms.iter().all(Expr::is_symbol));
    match as_logic(e) {
        Some((LogicOp::And, _)) => monomial(e),
        Some((LogicOp::Xor, terms)) => terms.iter().all(|a| monomial(a) || (is_literal(a) && !matches!(as_logic(a), Some((LogicOp::Not, _))))),
        _ => false,
    }
}

/// Convert to algebraic normal form without dense truth-table expansion.
/// With deep=False, conjunctions and XORs remain unchanged, and operands
/// of other connectives are treated as opaque predicates.
pub fn to_anf(e: &Expr, deep: bool) -> Expr {
    if is_anf(e) { return e.clone(); }
    if !deep && matches!(as_logic(e), Some((LogicOp::And | LogicOp::Xor | LogicOp::ITE, _))) {
        return e.clone();
    }
    expression(convert(e, deep))
}

/// In-place subset Möbius transform, in truth-table order (first variable
/// most significant). Applying this twice recovers the original table.
pub fn anf_coeffs<T: Clone + std::ops::BitXorAssign>(truthvalues: &[T]) -> Vec<T> {
    let mut values = truthvalues.to_vec();
    let mut bit = 1;
    while bit < values.len() {
        for mask in 0..values.len() {
            if mask & bit != 0 {
                let other = values[mask ^ bit].clone();
                values[mask] ^= other;
            }
        }
        bit <<= 1;
    }
    values
}

/// ANF reconstructed from a full result column, with variables ordered
/// as in SymPy (first variable corresponds to the most significant bit).
pub fn anf_form(vars: &[Expr], truthvalues: &[bool]) -> Expr {
    anf_from_coeffs(vars, &anf_coeffs(truthvalues))
}

/// Reconstruct from an already transformed coefficient vector.
pub fn anf_from_coeffs(vars: &[Expr], coefficients: &[bool]) -> Expr {
    let mut p = Polynomial::default();
    for (mask, &coefficient) in coefficients.iter().enumerate() {
        if coefficient {
            let mut m: Vec<Expr> = vars.iter().enumerate()
                .filter(|(i, _)| mask & (1 << (vars.len() - 1 - i)) != 0)
                .map(|(_, e)| e.clone()).collect();
            ordered(&mut m);
            m.dedup();
            toggle(&mut p, m);
        }
    }
    expression(p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::logic::{eval_bool, ite, not, or};

    #[test]
    fn mobius_is_involution() {
        for bits in 0u32..256 {
            let values: Vec<bool> = (0..8).map(|i| bits & (1 << i) != 0).collect();
            assert_eq!(anf_coeffs(&anf_coeffs(&values)), values);
            let vars = [Expr::symbol("a"), Expr::symbol("b"), Expr::symbol("c")];
            let e = anf_form(&vars, &values);
            assert!(is_anf(&e));
            for (row, expected) in values.iter().enumerate() {
                assert_eq!(eval_bool(&e, &|v| vars.iter().position(|x| x == v).map(|i| row & (1 << (2-i)) != 0)), Some(*expected));
            }
        }
    }

    #[test]
    fn nested_and_ite() {
        let (x, y, z) = (Expr::symbol("x"), Expr::symbol("y"), Expr::symbol("z"));
        let e = ite(&x, &or(vec![y.clone(), z.clone()]), &not(&y));
        let a = to_anf(&e, true);
        assert!(is_anf(&a));
        for row in 0..8 {
            let val = |v: &Expr| [x.clone(), y.clone(), z.clone()].iter().position(|x| x == v).map(|i| row & (1 << i) != 0);
            assert_eq!(eval_bool(&e, &val), eval_bool(&a, &val));
        }
        assert_eq!(to_anf(&e, false), e);
    }
}
