//! String printers: SymPy compatible `str` and `srepr`.

use std::cmp::Ordering;

use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::number::{format_float, Number};
use crate::ops::{add_term, mul_factor};

const PREC_ADD: u8 = 40;
const PREC_MUL: u8 = 50;
const PREC_POW: u8 = 60;
const PREC_ATOM: u8 = 100;

fn num_prec(n: &Number) -> u8 {
    match n {
        Number::Rational(_) => PREC_MUL,
        _ if n.is_negative() => PREC_ADD,
        _ => PREC_ATOM,
    }
}

fn prec(e: &Expr) -> u8 {
    match e.kind() {
        Kind::Num(n) => num_prec(n),
        Kind::Const(Constant::NegInfinity) => PREC_ADD,
        Kind::Add(..) => PREC_ADD,
        Kind::Mul(c, _) => {
            if c.is_negative() {
                PREC_ADD
            } else {
                PREC_MUL
            }
        }
        Kind::Pow(_, x) => {
            if x.is_minus_one() || matches!(x.as_number(), Some(n) if n.is_negative() && n.is_rational()) {
                PREC_MUL
            } else if *x == Expr::rational(1, 2) {
                PREC_ATOM
            } else {
                PREC_POW
            }
        }
        Kind::Rel(..) => 20,
        _ => PREC_ATOM,
    }
}

/// SymPy `str()` of an expression.
pub fn to_str(e: &Expr) -> String {
    if let Kind::Num(Number::Float(f)) = e.kind() {
        // SymPy prints a lone Float with full precision
        return crate::number::format_float_full(*f);
    }
    let mut s = String::new();
    write_expr(e, &mut s);
    s
}

fn paren(e: &Expr, p: u8, strict: bool, out: &mut String) {
    let ep = prec(e);
    if ep < p || (strict && ep <= p) {
        out.push('(');
        write_expr(e, out);
        out.push(')');
    } else {
        write_expr(e, out);
    }
}

fn write_expr(e: &Expr, out: &mut String) {
    match e.kind() {
        Kind::Num(n) => out.push_str(&n.to_string()),
        Kind::Sym(s) => out.push_str(s),
        Kind::Const(c) => out.push_str(c.name()),
        Kind::Bool(b) => out.push_str(if *b { "True" } else { "False" }),
        Kind::Add(..) => write_add(e, out),
        Kind::Mul(..) => write_mul(e, out),
        Kind::Pow(b, x) => write_pow(b, x, out),
        Kind::Func(Func::Order, args) => out.push_str(&crate::series::print::order_str(args)),
        Kind::Func(f, args) => {
            out.push_str(f.name());
            out.push('(');
            for (i, a) in args.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                write_expr(a, out);
            }
            out.push(')');
        }
        Kind::Derivative(f, vars) => {
            out.push_str("Derivative(");
            write_expr(f, out);
            for (v, n) in vars {
                out.push_str(", ");
                if *n == 1 {
                    write_expr(v, out);
                } else {
                    out.push('(');
                    write_expr(v, out);
                    out.push_str(&format!(", {n})"));
                }
            }
            out.push(')');
        }
        Kind::Rel(op, a, b) => match op {
            RelOp::Eq | RelOp::Ne => {
                out.push_str(if *op == RelOp::Eq { "Eq(" } else { "Ne(" });
                write_expr(a, out);
                out.push_str(", ");
                write_expr(b, out);
                out.push(')');
            }
            _ => {
                write_expr(a, out);
                out.push(' ');
                out.push_str(op.symbol());
                out.push(' ');
                write_expr(b, out);
            }
        },
        Kind::Piecewise(p) => {
            out.push_str("Piecewise(");
            for (i, (x, c)) in p.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                out.push('(');
                write_expr(x, out);
                out.push_str(", ");
                write_expr(c, out);
                out.push(')');
            }
            out.push(')');
        }
    }
}

// ----------------------------------------------------------------------
// term ordering (approximation of SymPy's lex ordering)
// ----------------------------------------------------------------------

#[derive(Clone)]
struct GenKey {
    rank: u8,
    name: String,
}

fn gen_rank(b: &Expr) -> u8 {
    match b.kind() {
        Kind::Sym(_) => 0,
        Kind::Const(_) => 1,
        Kind::Func(..) => 2,
        Kind::Num(_) => 4,
        _ => 3,
    }
}

/// Decompose a term into (generator, degree) pairs.
fn monomial(t: &Expr) -> Vec<(GenKey, f64)> {
    let mut v = Vec::new();
    let mut push = |b: &Expr, x: &Expr| {
        let deg = x.as_number().map(|n| n.to_f64());
        let (b, deg) = match deg {
            Some(d) => (b.clone(), d),
            None => (mul_factor(b, x), 1.0),
        };
        if b.is_number() {
            return;
        }
        v.push((
            GenKey {
                rank: gen_rank(&b),
                name: to_str(&b),
            },
            deg,
        ));
    };
    match t.kind() {
        Kind::Num(_) => {}
        Kind::Mul(_, facs) => {
            for (b, x) in facs {
                push(b, x)
            }
        }
        Kind::Pow(b, x) => push(b, x),
        _ => push(t, &Expr::one()),
    }
    v
}

fn cmp_gen(a: &GenKey, b: &GenKey) -> Ordering {
    a.rank.cmp(&b.rank).then_with(|| a.name.cmp(&b.name))
}

/// Lex comparison: `Less` means `a` is printed first.
fn cmp_monomials(a: &[(GenKey, f64)], b: &[(GenKey, f64)]) -> Ordering {
    let mut gens: Vec<&GenKey> = a.iter().map(|x| &x.0).chain(b.iter().map(|x| &x.0)).collect();
    gens.sort_by(|x, y| cmp_gen(x, y));
    gens.dedup_by(|x, y| cmp_gen(x, y) == Ordering::Equal);
    let deg = |m: &[(GenKey, f64)], g: &GenKey| -> f64 {
        m.iter()
            .filter(|(k, _)| cmp_gen(k, g) == Ordering::Equal)
            .map(|(_, d)| *d)
            .sum()
    };
    for g in gens {
        let (da, db) = (deg(a, g), deg(b, g));
        if da != db {
            return db.partial_cmp(&da).unwrap_or(Ordering::Equal);
        }
    }
    Ordering::Equal
}

/// Terms of an `Add` in printing order.
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    let Kind::Add(c, terms) = e.kind() else {
        return vec![e.clone()];
    };
    if e.may_have_order() && terms.iter().any(|(t, _)| t.is_func(&Func::Order)) {
        return crate::series::print::series_ordered_terms(e);
    }
    let mut keyed: Vec<(Vec<(GenKey, f64)>, String, Expr)> = terms
        .iter()
        .map(|(t, k)| {
            let full = add_term(t, k);
            (monomial(t), to_str(t), full)
        })
        .collect();
    keyed.sort_by(|a, b| cmp_monomials(&a.0, &b.0).then_with(|| a.1.cmp(&b.1)));
    let mut v: Vec<Expr> = keyed.into_iter().map(|x| x.2).collect();
    if !c.is_zero() {
        v.push(Expr::num(c.clone()));
    }
    v
}

fn write_add(e: &Expr, out: &mut String) {
    let terms = ordered_terms(e);
    for (i, t) in terms.iter().enumerate() {
        let mut s = String::new();
        paren(t, PREC_ADD, false, &mut s);
        if i == 0 {
            out.push_str(&s);
        } else if let Some(rest) = s.strip_prefix('-') {
            out.push_str(" - ");
            out.push_str(rest);
        } else {
            out.push_str(" + ");
            out.push_str(&s);
        }
    }
}

fn factor_rank(b: &Expr) -> u8 {
    match b.kind() {
        Kind::Num(_) => 0,
        Kind::Const(_) => 1,
        Kind::Sym(_) => 2,
        Kind::Pow(..) => 2,
        Kind::Func(..) => 3,
        Kind::Add(..) => 4,
        _ => 5,
    }
}

fn sort_factors(v: &mut [Expr]) {
    v.sort_by_cached_key(|f| {
        let (b, _) = crate::ops::as_base_exp(f);
        (factor_rank(&b), to_str(&b), to_str(f))
    });
}

fn write_mul(e: &Expr, out: &mut String) {
    let Kind::Mul(c, facs) = e.kind() else { unreachable!() };
    let mut c = c.clone();
    if c.is_negative() {
        out.push('-');
        c = c.neg();
    }
    let mut num: Vec<Expr> = Vec::new();
    let mut den: Vec<Expr> = Vec::new();
    match &c {
        Number::Rational(_) => {
            let n = c.numer();
            if !n.is_one() {
                num.push(Expr::num(n));
            }
            den.push(Expr::num(c.denom()));
        }
        _ => {
            if !c.is_one() {
                num.push(Expr::num(c.clone()));
            }
        }
    }
    let mut num_f: Vec<Expr> = Vec::new();
    let mut den_f: Vec<Expr> = Vec::new();
    for (b, x) in facs {
        match x.as_number() {
            Some(n) if n.is_negative() && n.is_rational() => {
                den_f.push(mul_factor(b, &Expr::num(n.neg())));
            }
            _ => num_f.push(mul_factor(b, x)),
        }
    }
    sort_factors(&mut num_f);
    sort_factors(&mut den_f);
    num.extend(num_f);
    // numeric denominator printed last ("3*x*y/2" but "x/(2*y)")
    let mut den_all = den_f;
    den_all.extend(den);
    if den_all.len() > 1 {
        // put the number first inside the denominator
        den_all.sort_by_key(|f| if f.is_number() { 0 } else { 1 });
    }
    let write_list = |v: &[Expr], out: &mut String| {
        for (i, f) in v.iter().enumerate() {
            if i > 0 {
                out.push('*');
            }
            paren(f, PREC_MUL, false, out);
        }
    };
    if num.is_empty() {
        out.push('1');
    } else {
        write_list(&num, out);
    }
    if !den_all.is_empty() {
        out.push('/');
        if den_all.len() == 1 {
            paren(&den_all[0], PREC_MUL, true, out);
        } else {
            out.push('(');
            write_list(&den_all, out);
            out.push(')');
        }
    }
}

fn write_pow(b: &Expr, x: &Expr, out: &mut String) {
    if *x == Expr::rational(1, 2) {
        out.push_str("sqrt(");
        write_expr(b, out);
        out.push(')');
        return;
    }
    if *x == Expr::rational(-1, 2) {
        out.push_str("1/sqrt(");
        write_expr(b, out);
        out.push(')');
        return;
    }
    if x.is_minus_one() {
        out.push_str("1/");
        paren(b, PREC_MUL, true, out);
        return;
    }
    // base
    let base_needs = match b.kind() {
        Kind::Num(n) => n.is_negative() || matches!(n, Number::Rational(_)),
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) | Kind::Rel(..) => true,
        Kind::Const(Constant::NegInfinity) => true,
        _ => false,
    };
    if base_needs {
        out.push('(');
        write_expr(b, out);
        out.push(')');
    } else {
        write_expr(b, out);
    }
    out.push_str("**");
    let exp_needs = match x.kind() {
        Kind::Num(n) => n.is_negative() || matches!(n, Number::Rational(_)),
        Kind::Add(..) | Kind::Mul(..) | Kind::Pow(..) => true,
        Kind::Const(Constant::NegInfinity) => true,
        _ => false,
    };
    if exp_needs {
        out.push('(');
        write_expr(x, out);
        out.push(')');
    } else {
        write_expr(x, out);
    }
}

// ----------------------------------------------------------------------
// srepr
// ----------------------------------------------------------------------

/// SymPy `srepr()` of an expression.
pub fn srepr(e: &Expr) -> String {
    let mut s = String::new();
    write_srepr(e, &mut s);
    s
}

fn write_srepr(e: &Expr, out: &mut String) {
    let list = |name: &str, args: &[Expr], out: &mut String| {
        out.push_str(name);
        out.push('(');
        for (i, a) in args.iter().enumerate() {
            if i > 0 {
                out.push_str(", ");
            }
            write_srepr(a, out);
        }
        out.push(')');
    };
    match e.kind() {
        Kind::Num(n) => match n {
            Number::Int(_) | Number::BigInt(_) => out.push_str(&format!("Integer({n})")),
            Number::Rational(r) => out.push_str(&format!("Rational({}, {})", r.numer(), r.denom())),
            Number::Float(f) => out.push_str(&format!("Float('{}', precision=53)", format_float(*f))),
        },
        Kind::Sym(s) => out.push_str(&format!("Symbol('{s}')")),
        Kind::Const(c) => out.push_str(c.name()),
        Kind::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Kind::Add(..) => list("Add", &ordered_terms(e), out),
        Kind::Mul(..) => list("Mul", &e.args(), out),
        Kind::Pow(..) => list("Pow", &e.args(), out),
        Kind::Func(Func::Undefined(n), args) => list(&format!("Function('{n}')"), args, out),
        Kind::Func(Func::Order, args) => {
            out.push_str("Order(");
            write_srepr(&args[0], out);
            for p in args[1..].chunks(2) {
                out.push_str(", Tuple(");
                write_srepr(&p[0], out);
                out.push_str(", ");
                write_srepr(&p[1], out);
                out.push(')');
            }
            out.push(')');
        }
        Kind::Func(f, args) => list(f.name(), args, out),
        Kind::Derivative(f, vars) => {
            out.push_str("Derivative(");
            write_srepr(f, out);
            for (v, n) in vars {
                out.push_str(", Tuple(");
                write_srepr(v, out);
                out.push_str(&format!(", Integer({n}))"));
            }
            out.push(')');
        }
        Kind::Rel(..) => list(e.class_name(), &e.args(), out),
        Kind::Piecewise(p) => {
            out.push_str("Piecewise(");
            for (i, (x, c)) in p.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                out.push_str("ExprCondPair(");
                write_srepr(x, out);
                out.push_str(", ");
                write_srepr(c, out);
                out.push(')');
            }
            out.push(')');
        }
    }
}
