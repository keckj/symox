//! The reference `PolyElement` / `FracElement` string form of elements of
//! composite coefficient domains (`ZZ[x][y]`, `ZZ_I[x]`, `QQ<sqrt(2)>[x]`).
//!
//! The reference prints a ring element term by term in the ring's monomial
//! order, each coefficient being printed by its own ground domain and
//! parenthesized by precedence. Elements of nested rings are flattened here,
//! so the levels are rebuilt before printing.
use super::algebraic::AN;
use super::domains::{DomainElement, ElementValue, ExactDomain};
use super::int::{Coeff, Int, Rat};
use super::ratfunc::RF;
use super::MPoly;
use crate::Expr;

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Prec {
    Add,
    Mul,
    Atom,
}

#[derive(Clone)]
enum Node {
    Number(Rat),
    Gaussian(Rat, Rat),
    Anp(AN),
    Poly { symbols: Vec<Expr>, terms: Vec<(Vec<u32>, Node)>, zero: Box<Node> },
    Frac { numer: Box<Node>, denom: Box<Node> },
}

fn rat_str(r: &Rat) -> String {
    Expr::num(r.to_number()).to_string()
}
fn mpq(r: &Rat) -> String {
    format!("MPQ({},{})", Expr::num(r.n.to_number()), Expr::num(r.d.to_number()))
}

impl Node {
    fn is_one(&self) -> bool {
        match self {
            Node::Number(r) => r.is_one(),
            Node::Gaussian(x, y) => x.is_one() && y.is_zero(),
            Node::Anp(a) => a.coefficients.len() == 1 && a.coefficients[0].is_one(),
            Node::Poly { terms, .. } => terms.len() == 1 && terms[0].0.iter().all(|e| *e == 0) && terms[0].1.is_one(),
            Node::Frac { numer, denom } => numer.is_one() && denom.is_one(),
        }
    }
    /// The reference `is_negative`: the sign of the leading coefficient,
    /// always false for Gaussian elements.
    fn is_negative(&self) -> bool {
        match self {
            Node::Number(r) => r.is_negative(),
            Node::Gaussian(..) => false,
            Node::Anp(a) => a.coefficients.last().is_some_and(|c| c.is_negative()),
            Node::Poly { terms, .. } => terms.first().is_some_and(|(_, c)| c.is_negative()),
            Node::Frac { numer, .. } => numer.is_negative(),
        }
    }
    fn neg(&self) -> Self {
        match self {
            Node::Number(r) => Node::Number(r.neg()),
            Node::Gaussian(x, y) => Node::Gaussian(x.neg(), y.neg()),
            Node::Anp(a) => Node::Anp(AN { coefficients: a.coefficients.iter().map(|c| c.neg()).collect(), field: a.field.clone() }),
            Node::Poly { symbols, terms, zero } => Node::Poly {
                symbols: symbols.clone(),
                terms: terms.iter().map(|(e, c)| (e.clone(), c.neg())).collect(),
                zero: zero.clone(),
            },
            Node::Frac { numer, denom } => Node::Frac { numer: Box::new(numer.neg()), denom: denom.clone() },
        }
    }
    fn precedence(&self) -> Prec {
        match self {
            Node::Number(_) | Node::Gaussian(..) | Node::Anp(_) => Prec::Atom,
            Node::Poly { terms, zero, .. } => {
                if terms.len() == 1 && terms[0].1.is_one() && terms[0].0.iter().sum::<u32>() == 1 {
                    Prec::Atom
                } else if terms.is_empty() {
                    zero.precedence()
                } else if terms.len() == 1 && terms[0].0.iter().all(|e| *e == 0) {
                    terms[0].1.precedence()
                } else if terms.len() == 1 {
                    Prec::Mul
                } else {
                    Prec::Add
                }
            }
            Node::Frac { numer, denom } => {
                if denom.is_one() {
                    numer.precedence()
                } else {
                    Prec::Mul
                }
            }
        }
    }
    fn parenthesize(&self, level: Prec) -> String {
        if self.precedence() < level {
            format!("({})", self.render())
        } else {
            self.render()
        }
    }
    fn render(&self) -> String {
        match self {
            Node::Number(r) => rat_str(r),
            Node::Gaussian(x, y) => format!("({} + {}*I)", rat_str(x), rat_str(y)),
            Node::Anp(a) => {
                let rep = a.coefficients.iter().rev().map(mpq).collect::<Vec<_>>().join(", ");
                let modulus = a.field.minimal.iter().rev().map(mpq).collect::<Vec<_>>().join(", ");
                format!("ANP([{rep}], [{modulus}], QQ)")
            }
            Node::Frac { numer, denom } => {
                if denom.is_one() {
                    numer.render()
                } else {
                    format!("{}/{}", numer.parenthesize(Prec::Mul), denom.parenthesize(Prec::Atom))
                }
            }
            Node::Poly { symbols, terms, zero } => {
                if terms.is_empty() {
                    return zero.render();
                }
                let mut parts: Vec<String> = Vec::new();
                for (exps, coeff) in terms {
                    let negative = coeff.is_negative();
                    parts.push(if negative { " - " } else { " + " }.to_string());
                    let scoeff = if exps.iter().all(|e| *e == 0) {
                        let s = coeff.render();
                        if negative && s.starts_with('-') {
                            s[1..].to_string()
                        } else {
                            s
                        }
                    } else {
                        let c = if negative { coeff.neg() } else { coeff.clone() };
                        if c.is_one() {
                            String::new()
                        } else {
                            c.parenthesize(Prec::Mul)
                        }
                    };
                    let mut factors = Vec::new();
                    if !scoeff.is_empty() {
                        factors.push(scoeff);
                    }
                    for (symbol, exp) in symbols.iter().zip(exps) {
                        if *exp == 0 {
                            continue;
                        }
                        // `parenthesize(symbol, PRECEDENCE["Atom"], strict=True)`:
                        // atoms such as `pi` print bare.
                        let name = if symbol.is_symbol() || symbol.args().is_empty() {
                            symbol.to_string()
                        } else {
                            format!("({symbol})")
                        };
                        factors.push(if *exp == 1 { name } else { format!("{name}**{exp}") });
                    }
                    parts.push(factors.join("*"));
                }
                if parts[0] == " + " {
                    parts.remove(0);
                } else if parts[0] == " - " {
                    parts[0] = "-".into();
                }
                parts.concat()
            }
        }
    }
}

/// Descending monomial order of the reference `PolyElement.terms()`.
fn sort_terms(terms: &mut [(Vec<u32>, Node)], order: &str) {
    let key = |e: &Vec<u32>| -> (u32, Vec<i64>) {
        let total = e.iter().sum::<u32>();
        match order {
            "grlex" => (total, e.iter().map(|&v| v as i64).collect()),
            "grevlex" => (total, e.iter().rev().map(|&v| -(v as i64)).collect()),
            _ => (0, e.iter().map(|&v| v as i64).collect()),
        }
    };
    terms.sort_by(|a, b| key(&b.0).cmp(&key(&a.0)));
}

fn constant(p: &MPoly<Int>) -> Int {
    p.constant_term()
}

type Terms<C> = Vec<(Vec<u32>, C)>;

/// The terms of `n` grouped by the variables `offset..offset+k`, the other
/// exponents being kept in the grouped coefficient.
fn group<C: Coeff>(p: &MPoly<C>, offset: usize, k: usize) -> Vec<(Vec<u32>, MPoly<C>)> {
    let mut groups: Vec<(Vec<u32>, Terms<C>)> = Vec::new();
    for (exps, c) in p.iter() {
        let own = exps[offset..offset + k].to_vec();
        let mut rest = exps.to_vec();
        for e in &mut rest[offset..offset + k] {
            *e = 0;
        }
        match groups.iter_mut().find(|(e, _)| *e == own) {
            Some((_, terms)) => terms.push((rest, c.clone())),
            None => groups.push((own, vec![(rest, c.clone())])),
        }
    }
    groups.into_iter().map(|(own, terms)| (own, MPoly::from_terms(p.nvars, terms))).collect()
}

fn ring_level(domain: &ExactDomain) -> ExactDomain {
    match domain {
        ExactDomain::Functions { integral, parameters, order, .. } => ExactDomain::Functions { integral: *integral, ring: true, parameters: parameters.clone(), order: order.clone() },
        ExactDomain::Composite { ground, parameters, order, .. } => ExactDomain::Composite { ground: ground.clone(), ring: true, parameters: parameters.clone(), order: order.clone() },
        other => other.clone(),
    }
}

fn zero_node(domain: &ExactDomain, nvars: usize, offset: usize) -> Node {
    match domain {
        ExactDomain::Gaussian { .. } => Node::Gaussian(Rat::zero(), Rat::zero()),
        ExactDomain::ZZ | ExactDomain::QQ => Node::Number(Rat::zero()),
        _ => rational_node(domain, &MPoly::zero(nvars), &MPoly::one(nvars), offset),
    }
}

/// `n/d` (flat over all parameters) as an element of `domain`, whose own
/// flat parameters start at `offset`.
fn rational_node(domain: &ExactDomain, n: &MPoly<Int>, d: &MPoly<Int>, offset: usize) -> Node {
    let nvars = n.nvars;
    match domain {
        ExactDomain::Functions { ring: false, .. } | ExactDomain::Composite { ring: false, .. } => {
            let ring = ring_level(domain);
            Node::Frac {
                numer: Box::new(rational_node(&ring, n, &MPoly::one(nvars), offset)),
                denom: Box::new(rational_node(&ring, d, &MPoly::one(nvars), offset)),
            }
        }
        ExactDomain::Functions { parameters, order, .. } => {
            let denominator = constant(d);
            let mut terms: Vec<(Vec<u32>, Node)> = n
                .iter()
                .map(|(exps, c)| (exps[offset..offset + parameters.len()].to_vec(), Node::Number(Rat::new(c.clone(), denominator.clone()))))
                .collect();
            sort_terms(&mut terms, order);
            Node::Poly { symbols: parameters.clone(), terms, zero: Box::new(Node::Number(Rat::zero())) }
        }
        ExactDomain::Composite { ground, parameters, order, .. } => {
            let k = parameters.len();
            let mut terms: Vec<(Vec<u32>, Node)> = group(n, offset, k)
                .into_iter()
                .map(|(own, sub)| {
                    let reduced = RF::new(sub, d.clone());
                    (own, rational_node(ground, &reduced.n, &reduced.d, offset + k))
                })
                .collect();
            sort_terms(&mut terms, order);
            Node::Poly { symbols: parameters.clone(), terms, zero: Box::new(zero_node(ground, nvars, offset + k)) }
        }
        _ => Node::Number(Rat::new(constant(n), constant(d))),
    }
}

fn algebraic_node(domain: &ExactDomain, rep: &MPoly<AN>, offset: usize) -> Node {
    match domain {
        ExactDomain::Composite { ground, parameters, order, .. } => {
            let k = parameters.len();
            let mut terms: Vec<(Vec<u32>, Node)> = group(rep, offset, k).into_iter().map(|(own, sub)| (own, algebraic_node(ground, &sub, offset + k))).collect();
            sort_terms(&mut terms, order);
            let zero = match &**ground {
                ExactDomain::Composite { .. } => algebraic_node(ground, &MPoly::zero(rep.nvars), offset + k),
                ExactDomain::Gaussian { .. } => Node::Gaussian(Rat::zero(), Rat::zero()),
                _ => Node::Anp(AN { coefficients: vec![], field: number_field(ground, rep) }),
            };
            Node::Poly { symbols: parameters.clone(), terms, zero: Box::new(zero) }
        }
        ExactDomain::Gaussian { .. } => {
            let c = rep.constant_term();
            Node::Gaussian(c.coefficients.first().cloned().unwrap_or_else(Rat::zero), c.coefficients.get(1).cloned().unwrap_or_else(Rat::zero))
        }
        _ => {
            let field = number_field(domain, rep);
            let c = rep.iter().find(|(e, _)| e.iter().all(|v| *v == 0)).map(|(_, c)| c.clone());
            Node::Anp(c.unwrap_or(AN { coefficients: vec![], field }))
        }
    }
}

fn number_field(domain: &ExactDomain, rep: &MPoly<AN>) -> std::sync::Arc<super::algebraic::NumberField> {
    match domain {
        ExactDomain::Algebraic(k) => k.clone(),
        _ => rep.coeffs().first().map(|c| c.field.clone()).unwrap_or_else(super::gaussian::field),
    }
}

/// The reference string of an element of a composite domain, or `None`
/// for the other domains.
pub fn composite_str(element: &DomainElement) -> Option<String> {
    if !matches!(&element.domain, ExactDomain::Composite { .. } | ExactDomain::Functions { .. }) {
        return None;
    }
    Some(match &element.value {
        ElementValue::Function(r) => rational_node(&element.domain, &r.n, &r.d, 0).render(),
        ElementValue::AlgebraicPoly(rep) => algebraic_node(&element.domain, rep, 0).render(),
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ops::{add2 as add, mul2 as mul};

    fn gaussian_ring(integral: bool, parameters: Vec<Expr>) -> ExactDomain {
        ExactDomain::composite(&ExactDomain::Gaussian { integral }, true, parameters, "lex").unwrap()
    }

    fn shown(domain: &ExactDomain, e: &Expr) -> String {
        composite_str(&domain.convert(e).unwrap()).unwrap()
    }

    #[test]
    fn gaussian_ring_elements_print_like_poly_elements() {
        let y = Expr::symbol("y");
        let zz_i = gaussian_ring(true, vec![y.clone()]);
        assert_eq!(shown(&zz_i, &mul(&Expr::i(), &y)), "(0 + 1*I)*y");
        assert_eq!(shown(&zz_i, &Expr::int(2)), "(2 + 0*I)");
        assert_eq!(shown(&zz_i, &add(&mul(&Expr::i(), &y), &Expr::one())), "(0 + 1*I)*y + (1 + 0*I)");
        let qq_i = gaussian_ring(false, vec![y.clone()]);
        assert_eq!(shown(&qq_i, &mul(&Expr::rational(1, 2), &mul(&Expr::i(), &y))), "(0 + 1/2*I)*y");
    }

    #[test]
    fn atom_generators_print_without_parentheses() {
        // parenthesize(symbol, PRECEDENCE["Atom"], strict=True)
        let pi = Expr::pi();
        let zz_i = gaussian_ring(true, vec![pi.clone()]);
        assert_eq!(shown(&zz_i, &mul(&Expr::i(), &pi)), "(0 + 1*I)*pi");
        // A power is not an atom and keeps its parentheses.
        let root = crate::ops::sqrt(&Expr::int(2));
        let zz_i = gaussian_ring(true, vec![root.clone()]);
        assert_eq!(shown(&zz_i, &mul(&Expr::i(), &root)), "(0 + 1*I)*(sqrt(2))");
    }
}
