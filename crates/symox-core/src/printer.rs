//! String printers: SymPy compatible `str` and `srepr`.

use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::number::Number;
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
        _ if is_mod(e) => PREC_MUL,
        Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") && args.len() == 1 => prec(&args[0]),
        _ => PREC_ATOM,
    }
}

/// SymPy `str()` of an expression.
pub fn to_str(e: &Expr) -> String {
    let _session = crate::printing::order::session();
    if let Kind::Num(Number::Float(f)) = e.kind() {
        // SymPy prints a lone Float with full precision
        return crate::number::format_float_full(*f);
    }
    if let Kind::Num(Number::MpFloat(value, precision)) = e.kind() {
        return value.to_string_dps(crate::bigfloat::prec_to_dps(*precision), false);
    }
    let mut s = String::new();
    write_expr(e, &mut s);
    s
}

/// SymPy gives `Mod` the precedence of `Mul`, so it is parenthesized inside
/// products with a positive coefficient and inside powers.
fn is_mod(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Func(Func::Special(crate::special::SpecialFn::Mod), _))
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

pub(crate) fn write_expr(e: &Expr, out: &mut String) {
    if e.has_host() {
        if let Some(text) = crate::host::print(e, crate::host::HostFormat::Str) {
            out.push_str(&text);
            return;
        }
    }
    match e.kind() {
        Kind::Num(n) => out.push_str(&n.to_string()),
        Kind::Sym(s) => {
            if s.is_dummy() { out.push('_'); }
            out.push_str(&crate::host::atom_name(e, s));
            if s.is_wild() {
                out.push('_');
            }
        }
        Kind::Const(c) => out.push_str(c.name()),
        Kind::Bool(b) => out.push_str(if *b { "True" } else { "False" }),
        Kind::Add(..) => write_add(e, out),
        Kind::Mul(..) => write_mul(e, out),
        Kind::Pow(b, x) => write_pow(b, x, out),
        Kind::Func(Func::Undefined(name), args) if name.is_internal("UnevaluatedExpr") && args.len() == 1 => {
            write_expr(&args[0], out);
        }
        Kind::Func(Func::Undefined(n), args) if n.is_internal("Tuple") => {
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
        Kind::Func(Func::Logic(op), args) => crate::logic::write_logic(*op, args, out),
        Kind::Func(Func::Set(op), args) => crate::sets::write_set(*op, args, out),
        Kind::Func(Func::Order, args) => out.push_str(&crate::series::print::order_str(args)),
        Kind::Func(f @ (Func::Integral|Func::Sum|Func::Product),args)=>{
            out.push_str(f.name());out.push('(');write_expr(&args[0],out);
            for limit in &args[1..] {
                out.push_str(", ");let values=limit.args();
                if values.len()==1 {write_expr(&values[0],out);} else {write_expr(limit,out);}
            }
            out.push(')');
        },
        Kind::Func(Func::Limit,args)=>{
            out.push_str("Limit(");
            for value in &args[..3]{write_expr(value,out);out.push_str(", ");}
            out.push_str("dir='");out.push_str(args[3].as_sym().map(|s|s.name()).unwrap_or("+"));out.push_str("')");
        },
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>{out.push_str("DotProduct(");out.push_str(&a.to_str());out.push_str(", ");out.push_str(&b.to_str());out.push(')');},
        Kind::Func(Func::MatrixProxy(parent), _) => out.push_str(&parent.to_str()),
        Kind::Func(f @ (Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent)), _) => {
            out.push_str(f.name()); out.push('('); out.push_str(&parent.to_str()); out.push(')');
        }
        Kind::Func(Func::MatrixElement(parent), args) => {
            out.push_str(&parent.to_str()); out.push('[');
            write_expr(&args[0], out); out.push_str(", ");
            write_expr(&args[1], out); out.push(']');
        }
        Kind::Func(Func::Subs,args)=>{
            out.push_str("Subs(");write_expr(&args[0],out);
            for tuple in &args[1..] {
                out.push_str(", ");let values=tuple.args();
                if values.len()==1{write_expr(&values[0],out);}else{write_expr(tuple,out);}
            }
            out.push(')');
        },
        Kind::Func(Func::Predicate(p), _) => out.push_str(&crate::ask::predicate_str(*p)),
        Kind::Func(Func::AppliedPredicate, args) => {
            if let Some((head, rest)) = args.split_first() {
                write_expr(head, out);
                out.push('(');
                for (i, arg) in rest.iter().enumerate() {
                    if i > 0 { out.push_str(", "); }
                    write_expr(arg, out);
                }
                out.push(')');
            }
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
    if let Kind::Add(_, terms) = e.kind() {
        if e.may_have_order() && terms.iter().any(|(t, _)| t.is_func(&Func::Order)) {
            return crate::series::print::series_ordered_terms(e);
        }
    }
    crate::printing::order::ordered_terms(e)
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

/// Sort `Mul` factors in SymPy's `as_ordered_factors` printing order.
pub fn sort_factors(v: &mut [Expr]) {
    crate::printing::order::sort_factors(v)
}

fn write_mul(e: &Expr, out: &mut String) {
    let Kind::Mul(c, facs) = e.kind() else { unreachable!() };
    let mut factors_src: Vec<(Expr, Expr)> = facs.clone();
    let mut c = c.clone();
    // `-oo*x` keeps `-oo` as a factor: print its sign in front like SymPy
    let neg_inf = factors_src.iter().position(|(b, x)| b.is_constant(Constant::NegInfinity) && x.is_one());
    if let Some(i) = neg_inf {
        factors_src[i].0 = Expr::oo();
        c = c.neg();
    }
    let negative = c.is_negative();
    if negative {
        out.push('-');
        c = c.neg();
    }
    let body = if c.is_one() && factors_src.len() == 1 {
        let (b, x) = factors_src.pop().unwrap();
        mul_factor(&b, &x)
    } else {
        Expr::from_kind(Kind::Mul(c, factors_src.clone()))
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
    let write_factor = |f: &Expr, out: &mut String| {
        paren(f, PREC_MUL, !negative && is_mod(f), out);
    };
    let write_list = |v: &[Expr], out: &mut String| {
        for (i, f) in v.iter().enumerate() {
            if i > 0 {
                out.push('*');
            }
            write_factor(f, out);
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
            write_factor(&den[0], out);
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
        _ if is_mod(b) => true,
        Kind::Func(Func::Undefined(name), _) if name.is_internal("UnevaluatedExpr") => prec(b) <= PREC_POW,
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
        _ if is_mod(x) => true,
        Kind::Func(Func::Undefined(name), _) if name.is_internal("UnevaluatedExpr") => prec(x) <= PREC_POW,
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
            // ``repr_dps(53)`` is seventeen significant digits, which recover
            // every finite f64; mpmath's ``to_str`` strips trailing zeros.
            // The human-readable str formatter intentionally rounds to 15.
            Number::Float(f) => out.push_str(&format!("Float('{}', precision=53)", crate::printing::mpf_to_str(*f, 17, true))),
            Number::MpFloat(value, precision) => out.push_str(&format!(
                "Float('{}', precision={precision})",
                // mpmath's ``repr_dps``
                value.to_string_dps(match crate::bigfloat::prec_to_dps(*precision) { 15 => 17, dps => dps + 3 }, true),
            )),
        },
        Kind::Sym(s) => {
            let name=crate::host::atom_name(e,s).replace('\\',"\\\\").replace('\'',"\\'").replace('\n',"\\n").replace('\r',"\\r");
            let class=crate::host::class_name(e).unwrap_or(if s.is_dummy() {"Dummy"} else if s.is_wild() {"Wild"} else {"Symbol"});
            out.push_str(&format!("{class}('{name}'"));
            for keyword in s.given_facts().keywords() { out.push_str(", "); out.push_str(&keyword); }
            if let Some(index)=s.dummy_index() { out.push_str(&format!(", dummy_index={index}")); }
            out.push(')');
        }
        Kind::Const(c) => out.push_str(c.name()),
        Kind::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Kind::Add(..) => list("Add", &ordered_terms(e), out),
        Kind::Mul(..) => list("Mul", &crate::sortkey::ordered_factors(e), out),
        Kind::Pow(..) => list("Pow", &e.args(), out),
        Kind::Func(Func::Undefined(n), args) if n.is_internal("Tuple") => list("Tuple", args, out),
        Kind::Func(Func::Undefined(n), args) if n.is_internal("UnevaluatedExpr") => list("UnevaluatedExpr", args, out),
        Kind::Func(Func::Undefined(n), args) if n.host().is_some() => list(crate::host::class_name(e).unwrap_or(&n.name), args, out),
        Kind::Func(Func::Undefined(n), args) => {
            let escaped=n.name.replace('\\',"\\\\").replace('\'',"\\'").replace('\n',"\\n").replace('\r',"\\r");
            let mut head=format!("Function('{escaped}'");
            for keyword in n.head_facts().keywords(){head.push_str(", ");head.push_str(&keyword);}
            if let Some(counts)=&n.nargs {
                head.push_str(", nargs=(");
                for count in counts {head.push_str(&format!("{count},"));}
                head.push(')');
            }
            head.push(')');list(&head,args,out);
        },
        Kind::Func(Func::Set(crate::sets::SetOp::IdentityFunction), _) => out.push_str("IdentityFunction()"),
        Kind::Func(Func::Set(op), args) if args.is_empty() => out.push_str(op.name()),
        Kind::Func(Func::MatrixDotProduct(a,b),_)=>{out.push_str("DotProduct(");out.push_str(&a.to_repr());out.push_str(", ");out.push_str(&b.to_repr());out.push(')');},
        Kind::Func(Func::MatrixProxy(parent), _) => out.push_str(&parent.to_repr()),
        Kind::Func(f @ (Func::MatrixTrace(parent) | Func::MatrixDeterminant(parent)), _) => {
            out.push_str(f.name()); out.push('('); out.push_str(&parent.to_repr()); out.push(')');
        }
        Kind::Func(Func::MatrixElement(parent), args) => {
            out.push_str("MatrixElement("); out.push_str(&parent.to_repr()); out.push_str(", ");
            write_srepr(&args[0], out); out.push_str(", ");
            write_srepr(&args[1], out); out.push(')');
        }
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
        Kind::Func(Func::Predicate(p), _) => {
            if p.is_undefined() { out.push_str("Predicate("); out.push_str(&format!("'{}'", p.name().replace('\\', "\\\\").replace('\'', "\\'"))); out.push(')'); }
            else { out.push_str(&crate::ask::predicate_str(*p)); }
        }
        Kind::Func(Func::AppliedPredicate, args) => list(e.class_name(), args, out),
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
