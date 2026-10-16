//! String printers: SymPy compatible `str` and `srepr`.

use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::number::{format_float, Number};
use crate::ops::mul_factor;

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
        Kind::Sym(s) => {
            out.push_str(s);
            if s.is_wild() {
                out.push('_');
            }
        }
        Kind::Const(c) => out.push_str(c.name()),
        Kind::Bool(b) => out.push_str(if *b { "True" } else { "False" }),
        Kind::Add(..) => write_add(e, out),
        Kind::Mul(..) => write_mul(e, out),
        Kind::Pow(b, x) => write_pow(b, x, out),
        Kind::Func(Func::Undefined(n), args) if &**n == "Tuple" => {
            out.push('(');
            for (i, a) in args.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                write_expr(a, out);
            }
            if args.len() == 1 {
                out.push(',');
            }
            out.push(')');
        }
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

/// Terms of an `Add` in SymPy's printing order.
pub fn ordered_terms(e: &Expr) -> Vec<Expr> {
    crate::sortkey::ordered_terms(e)
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

/// Negative-exponent test used to move factors to the denominator.
fn negative_exponent(x: &Expr) -> Option<Expr> {
    match x.kind() {
        Kind::Num(n) if n.is_negative() => Some(Expr::num(n.neg())),
        Kind::Mul(c, _) if c.is_negative() => Some(crate::ops::neg(x)),
        _ => None,
    }
}

fn write_mul(e: &Expr, out: &mut String) {
    let Kind::Mul(c, facs) = e.kind() else { unreachable!() };
    let mut factors_src: Vec<(Expr, Expr)> = facs.clone();
    let mut c = c.clone();
    if c.is_negative() {
        out.push('-');
        c = c.neg();
    }
    let body = if c.is_one() && factors_src.len() == 1 {
        let (b, x) = factors_src.pop().unwrap();
        mul_factor(&b, &x)
    } else {
        Expr::from_kind(Kind::Mul(c, facs.clone()))
    };
    let ordered = crate::sortkey::ordered_factors(&body);
    let mut num: Vec<Expr> = Vec::new();
    let mut den: Vec<Expr> = Vec::new();
    for item in ordered {
        match item.kind() {
            Kind::Pow(b, x) => match negative_exponent(x) {
                Some(px) => den.push(mul_factor(b, &px)),
                None => num.push(item.clone()),
            },
            Kind::Num(Number::Rational(r)) => {
                let p = Number::from_bigint(r.numer().clone());
                if !p.is_one() {
                    num.push(Expr::num(p));
                }
                den.push(Expr::num(Number::from_bigint(r.denom().clone())));
            }
            _ => num.push(item.clone()),
        }
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
    if !den.is_empty() {
        out.push('/');
        if den.len() == 1 {
            paren(&den[0], PREC_MUL, false, out);
        } else {
            out.push('(');
            write_list(&den, out);
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
        Kind::Sym(s) => {
            let kw = s.given().keywords();
            if kw.is_empty() {
                out.push_str(&format!("Symbol('{s}')"));
            } else {
                out.push_str(&format!("Symbol('{s}', {})", kw.join(", ")));
            }
        }
        Kind::Const(c) => out.push_str(c.name()),
        Kind::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Kind::Add(..) => list("Add", &ordered_terms(e), out),
        Kind::Mul(..) => list("Mul", &e.args(), out),
        Kind::Pow(..) => list("Pow", &e.args(), out),
        Kind::Func(Func::Undefined(n), args) => list(&format!("Function('{n}')"), args, out),
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
