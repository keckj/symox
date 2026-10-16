//! A small parser for SymPy-style expression strings (`sympify`).
//!
//! Supports `+ - * / ** ^`, parentheses, function calls, relationals,
//! numbers (integers, decimals, exponents), the usual constants and the
//! constructors printed by `srepr` (`Integer`, `Rational`, `Symbol`, ...).

use std::fmt;

use num_bigint::BigInt;

use crate::expr::{Constant, Expr, Func, RelOp};
use crate::functions::{func, log_base, undefined};
use crate::number::Number;
use crate::ops::{add, div, mul, neg, pow, sqrt};
use crate::relational::{piecewise, rel};

#[derive(Debug, Clone, PartialEq)]
pub struct ParseError {
    pub msg: String,
    pub pos: usize,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "parse error at {}: {}", self.pos, self.msg)
    }
}

impl std::error::Error for ParseError {}

#[derive(Debug, Clone, PartialEq)]
enum Tok {
    Num(String),
    Ident(String),
    Str(String),
    Op(&'static str),
    LParen,
    RParen,
    Comma,
    End,
}

fn tokenize(s: &str) -> Result<Vec<(Tok, usize)>, ParseError> {
    let b = s.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    while i < b.len() {
        let c = b[i] as char;
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        let start = i;
        if c.is_ascii_digit() || (c == '.' && i + 1 < b.len() && (b[i + 1] as char).is_ascii_digit()) {
            while i < b.len() && ((b[i] as char).is_ascii_digit() || b[i] == b'.' || b[i] == b'_') {
                i += 1;
            }
            if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
                let save = i;
                i += 1;
                if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
                    i += 1;
                }
                if i < b.len() && (b[i] as char).is_ascii_digit() {
                    while i < b.len() && (b[i] as char).is_ascii_digit() {
                        i += 1;
                    }
                } else {
                    i = save;
                }
            }
            out.push((Tok::Num(s[start..i].replace('_', "")), start));
            continue;
        }
        if c.is_alphabetic() || c == '_' {
            while i < b.len() && ((b[i] as char).is_alphanumeric() || b[i] == b'_') {
                i += 1;
            }
            out.push((Tok::Ident(s[start..i].to_string()), start));
            continue;
        }
        if c == '\'' || c == '"' {
            i += 1;
            while i < b.len() && b[i] as char != c {
                i += 1;
            }
            if i >= b.len() {
                return Err(ParseError { msg: "unterminated string".into(), pos: start });
            }
            out.push((Tok::Str(s[start + 1..i].to_string()), start));
            i += 1;
            continue;
        }
        let two = if i + 1 < b.len() { &s[i..i + 2] } else { "" };
        let op: Option<&'static str> = match two {
            "**" => Some("**"),
            "<=" => Some("<="),
            ">=" => Some(">="),
            "==" => Some("=="),
            "!=" => Some("!="),
            _ => None,
        };
        if let Some(op) = op {
            out.push((Tok::Op(op), start));
            i += 2;
            continue;
        }
        let tok = match c {
            '+' => Tok::Op("+"),
            '-' => Tok::Op("-"),
            '*' => Tok::Op("*"),
            '/' => Tok::Op("/"),
            '^' => Tok::Op("**"),
            '<' => Tok::Op("<"),
            '>' => Tok::Op(">"),
            '=' => Tok::Op("="),
            '(' | '[' => Tok::LParen,
            ')' | ']' => Tok::RParen,
            ',' => Tok::Comma,
            _ => {
                return Err(ParseError {
                    msg: format!("unexpected character {c:?}"),
                    pos: start,
                })
            }
        };
        out.push((tok, start));
        i += 1;
    }
    out.push((Tok::End, s.len()));
    Ok(out)
}

/// Parser values: expressions, tuples and strings (for `Symbol('x')`).
#[derive(Debug, Clone)]
enum Val {
    E(Expr),
    Tuple(Vec<Val>),
    Str(String),
}

struct Parser {
    toks: Vec<(Tok, usize)>,
    pos: usize,
}

type R<T> = Result<T, ParseError>;

impl Parser {
    fn peek(&self) -> &Tok {
        &self.toks[self.pos].0
    }
    fn at(&self) -> usize {
        self.toks[self.pos].1
    }
    fn next(&mut self) -> Tok {
        let t = self.toks[self.pos].0.clone();
        if self.pos + 1 < self.toks.len() {
            self.pos += 1;
        }
        t
    }
    fn err<T>(&self, msg: impl Into<String>) -> R<T> {
        Err(ParseError { msg: msg.into(), pos: self.at() })
    }
    fn expect(&mut self, t: Tok) -> R<()> {
        if *self.peek() == t {
            self.next();
            Ok(())
        } else {
            self.err(format!("expected {t:?}, found {:?}", self.peek()))
        }
    }

    fn expr(&mut self) -> R<Expr> {
        match self.value()? {
            Val::E(e) => Ok(e),
            v => self.err(format!("expected an expression, found {v:?}")),
        }
    }

    fn value(&mut self) -> R<Val> {
        let lhs = self.arith()?;
        let op = match self.peek() {
            Tok::Op("<") => Some(RelOp::Lt),
            Tok::Op("<=") => Some(RelOp::Le),
            Tok::Op(">") => Some(RelOp::Gt),
            Tok::Op(">=") => Some(RelOp::Ge),
            Tok::Op("==") => Some(RelOp::Eq),
            Tok::Op("!=") => Some(RelOp::Ne),
            _ => None,
        };
        if let Some(op) = op {
            self.next();
            let rhs = self.arith()?;
            let (Val::E(a), Val::E(b)) = (lhs, rhs) else {
                return self.err("relational between non-expressions");
            };
            return Ok(Val::E(rel(op, &a, &b)));
        }
        Ok(lhs)
    }

    fn arith(&mut self) -> R<Val> {
        let first = self.term()?;
        if !matches!(self.peek(), Tok::Op("+") | Tok::Op("-")) {
            return Ok(first);
        }
        let mut terms = vec![self.as_expr(first)?];
        loop {
            match self.peek() {
                Tok::Op("+") => {
                    self.next();
                    let t = self.term()?;
                    terms.push(self.as_expr(t)?);
                }
                Tok::Op("-") => {
                    self.next();
                    let t = self.term()?;
                    terms.push(neg(&self.as_expr(t)?));
                }
                _ => break,
            }
        }
        Ok(Val::E(add(&terms)))
    }

    fn as_expr(&self, v: Val) -> R<Expr> {
        match v {
            Val::E(e) => Ok(e),
            v => self.err(format!("expected an expression, found {v:?}")),
        }
    }

    fn term(&mut self) -> R<Val> {
        let first = self.unary()?;
        if !matches!(self.peek(), Tok::Op("*") | Tok::Op("/")) {
            return Ok(first);
        }
        let mut acc = self.as_expr(first)?;
        loop {
            match self.peek() {
                Tok::Op("*") => {
                    self.next();
                    let f = self.unary()?;
                    acc = mul(&[acc, self.as_expr(f)?]);
                }
                Tok::Op("/") => {
                    self.next();
                    let f = self.unary()?;
                    acc = div(&acc, &self.as_expr(f)?);
                }
                _ => break,
            }
        }
        Ok(Val::E(acc))
    }

    fn unary(&mut self) -> R<Val> {
        match self.peek() {
            Tok::Op("-") => {
                self.next();
                let v = self.unary()?;
                let e = self.as_expr(v)?;
                Ok(Val::E(neg(&e)))
            }
            Tok::Op("+") => {
                self.next();
                self.unary()
            }
            _ => self.power(),
        }
    }

    fn power(&mut self) -> R<Val> {
        let base = self.postfix()?;
        if matches!(self.peek(), Tok::Op("**")) {
            self.next();
            let ex = self.unary()?; // right associative, allows 2**-1
            let (b, e) = (self.as_expr(base)?, self.as_expr(ex)?);
            return Ok(Val::E(pow(&b, &e)));
        }
        Ok(base)
    }

    fn postfix(&mut self) -> R<Val> {
        let mut v = self.atom()?;
        // Function('f')(x) style application
        while matches!(self.peek(), Tok::LParen) {
            let name = match &v {
                Val::Str(s) => s.clone(),
                _ => break,
            };
            let args = self.call_args()?;
            let args: R<Vec<Expr>> = args.into_iter().map(|a| self.as_expr(a)).collect();
            v = Val::E(undefined(&name, args?));
        }
        Ok(v)
    }

    fn call_args(&mut self) -> R<Vec<Val>> {
        self.expect(Tok::LParen)?;
        let mut args = Vec::new();
        if matches!(self.peek(), Tok::RParen) {
            self.next();
            return Ok(args);
        }
        loop {
            // keyword arguments (evaluate=False, precision=53) are ignored
            if matches!(self.peek(), Tok::Ident(_))
                && matches!(self.toks.get(self.pos + 1), Some((Tok::Op("="), _)))
            {
                self.next();
                self.next();
                self.value()?;
                match self.next() {
                    Tok::Comma => continue,
                    Tok::RParen => break,
                    t => return self.err(format!("expected ',' or ')', found {t:?}")),
                }
            }
            args.push(self.value()?);
            match self.next() {
                Tok::Comma => {
                    if matches!(self.peek(), Tok::RParen) {
                        self.next();
                        break;
                    }
                }
                Tok::RParen => break,
                t => return self.err(format!("expected ',' or ')', found {t:?}")),
            }
        }
        Ok(args)
    }

    fn atom(&mut self) -> R<Val> {
        match self.next() {
            Tok::Num(s) => Ok(Val::E(parse_number(&s).ok_or_else(|| ParseError {
                msg: format!("bad number {s}"),
                pos: self.at(),
            })?)),
            Tok::Str(s) => Ok(Val::Str(s)),
            Tok::LParen => {
                if matches!(self.peek(), Tok::RParen) {
                    self.next();
                    return Ok(Val::Tuple(vec![]));
                }
                let first = self.value()?;
                if matches!(self.peek(), Tok::Comma) {
                    let mut items = vec![first];
                    while matches!(self.peek(), Tok::Comma) {
                        self.next();
                        if matches!(self.peek(), Tok::RParen) {
                            break;
                        }
                        items.push(self.value()?);
                    }
                    self.expect(Tok::RParen)?;
                    return Ok(Val::Tuple(items));
                }
                self.expect(Tok::RParen)?;
                Ok(first)
            }
            Tok::Ident(name) => {
                if matches!(self.peek(), Tok::LParen) {
                    let args = self.call_args()?;
                    return self.call(&name, args);
                }
                Ok(Val::E(ident(&name)))
            }
            t => self.err(format!("unexpected token {t:?}")),
        }
    }

    fn call(&mut self, name: &str, args: Vec<Val>) -> R<Val> {
        let exprs = |p: &Parser, args: Vec<Val>| -> R<Vec<Expr>> {
            args.into_iter().map(|a| p.as_expr(a)).collect()
        };
        let one = |p: &Parser, args: Vec<Val>| -> R<Expr> {
            let mut v = exprs(p, args)?;
            if v.len() != 1 {
                return p.err(format!("{name} expects one argument"));
            }
            Ok(v.pop().unwrap())
        };
        Ok(Val::E(match name {
            "Symbol" | "Dummy" => match args.first() {
                Some(Val::Str(s)) => Expr::symbol(s),
                _ => return self.err("Symbol expects a string"),
            },
            "Function" => match args.first() {
                Some(Val::Str(s)) => return Ok(Val::Str(s.clone())),
                _ => return self.err("Function expects a string"),
            },
            "Integer" | "S" | "sympify" | "Float" | "N" => {
                let v = match args.into_iter().next() {
                    Some(Val::Str(s)) => parse_str(&s)?,
                    Some(Val::E(e)) => e,
                    _ => return self.err(format!("{name} expects one argument")),
                };
                if name == "Float" || name == "N" {
                    crate::evalf::evalf(&v)
                } else {
                    v
                }
            }
            "Rational" => {
                let v = exprs(self, args)?;
                match v.len() {
                    1 => v[0].clone(),
                    2 => div(&v[0], &v[1]),
                    _ => return self.err("Rational expects 1 or 2 arguments"),
                }
            }
            "Add" => add(&exprs(self, args)?),
            "Mul" => mul(&exprs(self, args)?),
            "Pow" => {
                let v = exprs(self, args)?;
                if v.len() != 2 {
                    return self.err("Pow expects 2 arguments");
                }
                pow(&v[0], &v[1])
            }
            "sqrt" => sqrt(&one(self, args)?),
            "cbrt" => pow(&one(self, args)?, &Expr::rational(1, 3)),
            "root" => {
                let v = exprs(self, args)?;
                if v.len() != 2 {
                    return self.err("root expects 2 arguments");
                }
                pow(&v[0], &div(&Expr::one(), &v[1]))
            }
            "log" | "ln" if args.len() == 2 => {
                let v = exprs(self, args)?;
                log_base(&v[0], &v[1])
            }
            "Eq" | "Ne" | "Lt" | "Le" | "Gt" | "Ge" => {
                let v = exprs(self, args)?;
                if v.len() != 2 {
                    return self.err(format!("{name} expects 2 arguments"));
                }
                let op = match name {
                    "Eq" => RelOp::Eq,
                    "Ne" => RelOp::Ne,
                    "Lt" => RelOp::Lt,
                    "Le" => RelOp::Le,
                    "Gt" => RelOp::Gt,
                    _ => RelOp::Ge,
                };
                rel(op, &v[0], &v[1])
            }
            "Piecewise" => {
                let mut pieces = Vec::new();
                for a in args {
                    match a {
                        Val::Tuple(t) if t.len() == 2 => {
                            let mut it = t.into_iter();
                            let e = self.as_expr(it.next().unwrap())?;
                            let c = self.as_expr(it.next().unwrap())?;
                            pieces.push((e, c));
                        }
                        _ => return self.err("Piecewise expects (expr, cond) pairs"),
                    }
                }
                piecewise(pieces)
            }
            "Derivative" | "diff" => {
                let mut it = args.into_iter();
                let f = match it.next() {
                    Some(v) => self.as_expr(v)?,
                    None => return self.err("Derivative expects arguments"),
                };
                let mut vars = Vec::new();
                for a in it {
                    match a {
                        Val::E(v) => vars.push((v, 1u32)),
                        Val::Tuple(t) if t.len() == 2 => {
                            let mut it = t.into_iter();
                            let v = self.as_expr(it.next().unwrap())?;
                            let n = self.as_expr(it.next().unwrap())?;
                            let n = n.as_i64().filter(|n| *n >= 0).ok_or_else(|| ParseError {
                                msg: "derivative order must be a non-negative integer".into(),
                                pos: self.at(),
                            })?;
                            vars.push((v, n as u32));
                        }
                        _ => return self.err("bad Derivative variable"),
                    }
                }
                if name == "diff" {
                    let mut r = f;
                    for (v, n) in vars {
                        r = crate::diff::diff_n(&r, &v, n);
                    }
                    r
                } else {
                    crate::diff::derivative_unevaluated(&f, vars)
                }
            }
            "O" | "Order" => {
                let mut it = args.into_iter();
                let e = match it.next() {
                    Some(v) => self.as_expr(v)?,
                    None => return self.err("Order expects arguments"),
                };
                let mut vars = Vec::new();
                for a in it {
                    match a {
                        Val::Tuple(t) if t.len() == 2 => {
                            let mut it = t.into_iter();
                            let v = self.as_expr(it.next().unwrap())?;
                            let p = self.as_expr(it.next().unwrap())?;
                            vars.push((v, p));
                        }
                        v => vars.push((self.as_expr(v)?, Expr::zero())),
                    }
                }
                crate::series::order::order(&e, &vars)
            }
            "Tuple" => {
                return Ok(Val::Tuple(args));
            }
            "ExprCondPair" => {
                return Ok(Val::Tuple(args));
            }
            "expand" => crate::expand::expand(&one(self, args)?),
            _ => {
                let v = exprs(self, args)?;
                match Func::from_name(name) {
                    Some(f) => {
                        if let Some(n) = f.arity() {
                            if v.len() != n {
                                return self.err(format!("{name} expects {n} argument(s)"));
                            }
                        }
                        func(f, v)
                    }
                    None => undefined(name, v),
                }
            }
        }))
    }
}

fn ident(name: &str) -> Expr {
    match name {
        "pi" => Expr::pi(),
        "E" => Expr::e(),
        "I" => Expr::i(),
        "oo" => Expr::oo(),
        "zoo" => Expr::zoo(),
        "nan" => Expr::nan(),
        "EulerGamma" => Expr::constant(Constant::EulerGamma),
        "True" | "true" => Expr::boolean(true),
        "False" | "false" => Expr::boolean(false),
        _ => Expr::symbol(name),
    }
}

/// Parse a numeric literal. Decimal literals become floats like in SymPy.
pub fn parse_number(s: &str) -> Option<Expr> {
    if s.contains(['.', 'e', 'E']) {
        return s.parse::<f64>().ok().map(Expr::float);
    }
    if let Ok(i) = s.parse::<i64>() {
        return Some(Expr::int(i));
    }
    s.parse::<BigInt>().ok().map(|b| Expr::num(Number::from_bigint(b)))
}

/// Parse a string into an expression.
pub fn parse_str(s: &str) -> Result<Expr, ParseError> {
    let toks = tokenize(s)?;
    let mut p = Parser { toks, pos: 0 };
    let e = p.expr()?;
    if !matches!(p.peek(), Tok::End) {
        return p.err(format!("unexpected trailing token {:?}", p.peek()));
    }
    Ok(e)
}

/// `sympify` for strings.
pub fn sympify(s: &str) -> Result<Expr, ParseError> {
    parse_str(s)
}
