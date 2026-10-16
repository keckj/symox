//! A small parser for SymPy-style expression strings (`sympify`).
//!
//! Supports `+ - * / ** ^`, parentheses, function calls, relationals,
//! numbers (integers, decimals, exponents), the usual constants and the
//! constructors printed by `srepr` (`Integer`, `Rational`, `Symbol`, ...).

use std::fmt;
use std::sync::Arc;

use num_bigint::BigInt;

use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::functions::{func, log_base, undefined};
use crate::number::Number;
use crate::ops::{add, div, mul, neg, pow, sqrt};
use crate::relational::{piecewise, rel};
use crate::matrixexpr::MatrixExpr;

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
    LBracket,
    RBracket,
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
            let mut value=String::new();
            while i < b.len() && b[i] as char != c {
                let ch=s[i..].chars().next().unwrap();
                i+=ch.len_utf8();
                if ch=='\\' {
                    let escaped=s[i..].chars().next().ok_or_else(||ParseError{msg:"unterminated string escape".into(),pos:start})?;
                    i+=escaped.len_utf8();
                    match escaped {'n'=>value.push('\n'),'r'=>value.push('\r'),'t'=>value.push('\t'),'\\'|'\''|'"'=>value.push(escaped),_=>{value.push('\\');value.push(escaped);}}
                }else{value.push(ch);}
            }
            if i >= b.len() {
                return Err(ParseError { msg: "unterminated string".into(), pos: start });
            }
            out.push((Tok::Str(value), start));
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
            '&' => Tok::Op("&"),
            '|' => Tok::Op("|"),
            '~' => Tok::Op("~"),
            '(' => Tok::LParen,
            ')' => Tok::RParen,
            '[' => Tok::LBracket,
            ']' => Tok::RBracket,
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
    /// boolean keyword arguments of the call being parsed
    kwargs: Vec<(String, bool)>,
    precision: Option<i64>,
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

    /// Boolean operators `|` and `&` (looser than relationals, which SymPy
    /// always parenthesizes inside them).
    fn value(&mut self) -> R<Val> {
        let first = self.bool_and()?;
        if !matches!(self.peek(), Tok::Op("|")) {
            return Ok(first);
        }
        let mut v = vec![self.as_expr(first)?];
        while matches!(self.peek(), Tok::Op("|")) {
            self.next();
            let t = self.bool_and()?;
            v.push(self.as_expr(t)?);
        }
        Ok(Val::E(crate::logic::or(v)))
    }

    fn bool_and(&mut self) -> R<Val> {
        let first = self.relational()?;
        if !matches!(self.peek(), Tok::Op("&")) {
            return Ok(first);
        }
        let mut v = vec![self.as_expr(first)?];
        while matches!(self.peek(), Tok::Op("&")) {
            self.next();
            let t = self.relational()?;
            v.push(self.as_expr(t)?);
        }
        Ok(Val::E(crate::logic::and(v)))
    }

    fn relational(&mut self) -> R<Val> {
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

    /// An expression, with tuples converted to `Tuple(...)` (for
    /// `hyper((a,), (b,), z)`).
    fn tuple_expr(&self, v: Val) -> R<Expr> {
        match v {
            Val::Tuple(t) => {
                let items = t.into_iter().map(|x| self.tuple_expr(x)).collect::<R<Vec<Expr>>>()?;
                Ok(crate::special::tuple(items))
            }
            v => self.as_expr(v),
        }
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
            Tok::Op("~") => {
                self.next();
                let v = self.unary()?;
                let e = self.as_expr(v)?;
                Ok(Val::E(crate::logic::not(&e)))
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
            // keyword arguments: booleans are kept (Symbol assumptions),
            // others (precision=53) are ignored
            if matches!(self.peek(), Tok::Ident(_))
                && matches!(self.toks.get(self.pos + 1), Some((Tok::Op("="), _)))
            {
                let Tok::Ident(key) = self.next() else { unreachable!() };
                self.next();
                let v = self.value()?;
                if key == "precision" {
                    self.precision = Some(match &v {
                        Val::E(e) => e.as_i64().ok_or_else(|| ParseError { msg: "precision must be an integer".into(), pos: self.at() })?,
                        _ => return self.err("precision must be an integer"),
                    });
                }
                if let Val::E(e) = &v {
                    if let crate::expr::Kind::Bool(b) = e.kind() {
                        self.kwargs.push((key, *b));
                    }
                }
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
            Tok::LBracket => {
                let mut items=Vec::new();
                while !matches!(self.peek(),Tok::RBracket) {
                    items.push(self.value()?);
                    if !matches!(self.peek(),Tok::Comma){break;}
                    self.next();
                }
                self.expect(Tok::RBracket)?;
                Ok(Val::Tuple(items))
            }
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
                    let saved = std::mem::take(&mut self.kwargs);
                    let saved_precision = self.precision.take();
                    let args = self.call_args()?;
                    let kwargs = std::mem::replace(&mut self.kwargs, saved);
                    let precision = std::mem::replace(&mut self.precision, saved_precision);
                    if matches!(name.as_str(), "MatrixSymbol" | "MatAdd" | "MatMul" | "MatPow" | "Identity" | "ZeroMatrix" | "Inverse" | "Transpose" | "Adjoint" | "MatrixConjugate" | "MatrixElement" | "Trace" | "Determinant" | "ImmutableMatrix" | "ImmutableDenseMatrix" | "MutableDenseMatrix" | "Matrix") {
                        let evaluate=kwargs.iter().find(|(k,_)|k=="evaluate").is_some_and(|(_,v)|*v);
                        return self.matrix_call(&name,args,evaluate);
                    }
                    if name == "Float" {
                        if args.is_empty() || args.len() > 2 { return self.err("Float expects one value and optional decimal precision"); }
                        if precision.is_some() && args.len() == 2 { return self.err("specify decimal digits or bit precision, not both"); }
                        let bits = match (precision, args.get(1)) {
                            (Some(p), _) if (2..=40_000).contains(&p) => p as u32,
                            (Some(_), _) => return self.err("bit precision must lie between 2 and 40000"),
                            (None, Some(Val::E(e))) => {
                                let n = e.as_i64().filter(|n| (1..=10_000).contains(n)).ok_or_else(|| ParseError { msg: "decimal precision must lie between 1 and 10000".into(), pos: self.at() })?;
                                crate::bigfloat::dps_to_prec(n as u32)
                            }
                            (None, Some(_)) => return self.err("decimal precision must be an integer"),
                            _ => 53,
                        };
                        let value = match &args[0] {
                            Val::Str(s) => crate::bigfloat::BigFloat::from_decimal(s, bits)
                                .map_err(|msg| ParseError { msg, pos: self.at() })?,
                            Val::E(e) => match e.as_number() {
                                Some(number) => number.to_bigfloat(bits),
                                None => return self.err("Float expects a real number"),
                            },
                            _ => return self.err("Float expects a real number"),
                        };
                        return Ok(Val::E(Expr::num(Number::mp_float(value, bits))));
                    }
                    if (name == "Symbol" || name == "Dummy") && !kwargs.is_empty() {
                        if let Some(Val::Str(s)) = args.first() {
                            let a = crate::symbol::Assumptions::from_kwargs(
                                kwargs.iter().map(|(k, v)| (k.as_str(), *v)),
                            )
                            .map_err(|m| ParseError { msg: m, pos: self.at() })?;
                            return Ok(Val::E(Expr::symbol_with(s, a)));
                        }
                    }
                    if matches!(name.as_str(), "Add" | "Mul" | "Pow")
                        && kwargs.iter().any(|(key, value)| key == "evaluate" && !value)
                    {
                        let v: Vec<Expr> = args.into_iter().map(|a| self.as_expr(a)).collect::<R<_>>()?;
                        if name == "Pow" && v.len() != 2 { return self.err("Pow expects 2 arguments"); }
                        if v.len() >= 2 {
                            let kind = match name.as_str() {
                                "Add" => Kind::Add(Number::Int(0), v.into_iter().map(|a| (a, Number::Int(1))).collect()),
                                "Mul" => Kind::Mul(Number::Int(1), v.into_iter().map(|a| (a, Expr::one())).collect()),
                                "Pow" if v.len() == 2 => Kind::Pow(v[0].clone(), v[1].clone()),
                                _ => return self.err("Pow expects 2 arguments"),
                            };
                            return Ok(Val::E(Expr::from_unevaluated_kind(kind)));
                        }
                        return Ok(Val::E(if name == "Add" { add(&v) } else { mul(&v) }));
                    }
                    return self.call(&name, args);
                }
                Ok(Val::E(ident(&name)))
            }
            t => self.err(format!("unexpected token {t:?}")),
        }
    }

    fn as_matrix(&self, value: &Val) -> R<MatrixExpr> {
        if let Val::E(e)=value { if let Kind::Func(Func::MatrixProxy(m),_)=e.kind() { return Ok((**m).clone()); } }
        self.err("expected a matrix expression")
    }

    fn matrix_call(&self, name:&str, args:Vec<Val>, evaluate:bool)->R<Val> {
        let fail=|msg:String| ParseError{msg,pos:self.at()};
        let scalar=|i:usize| ->R<Expr>{args.get(i).cloned().ok_or_else(||fail(format!("{name}: missing argument"))).and_then(|a|self.as_expr(a))};
        let matrix=|i:usize| ->R<MatrixExpr>{self.as_matrix(args.get(i).ok_or_else(||fail(format!("{name}: missing matrix argument")))?)};
        let arity=|n:usize|->R<()>{if args.len()==n {Ok(())}else{self.err(format!("{name} expects {n} arguments"))}};
        let result=match name {
            "MatrixSymbol"=>{arity(3)?;let Val::Str(s)=&args[0] else{return self.err("MatrixSymbol expects a string name");};MatrixExpr::symbol(s,scalar(1)?,scalar(2)?)},
            "Identity"=>{arity(1)?;MatrixExpr::identity(scalar(0)?)},
            "ZeroMatrix"=>{arity(2)?;MatrixExpr::zero(scalar(0)?,scalar(1)?)},
            "MatAdd"=>MatrixExpr::add(&args.iter().map(|a|self.as_matrix(a)).collect::<R<Vec<_>>>()?,evaluate),
            "MatMul"=>{
                let mut matrices=Vec::new();let mut scalars=Vec::new();
                for a in &args {if let Ok(m)=self.as_matrix(a){matrices.push(m);}else{scalars.push(self.as_expr(a.clone())?);}}
                MatrixExpr::mul(mul(&scalars),&matrices,evaluate)
            },
            "MatPow"=>{arity(2)?;matrix(0)?.pow(scalar(1)?,evaluate)},
            "Inverse"=>{
                if args.len()!=1 && !(args.len()==2 && scalar(1)?==Expr::minus_one()){return self.err("Inverse expects a matrix and optional exponent -1");}
                matrix(0)?.inverse(evaluate)
            },
            "Transpose"=>{arity(1)?;matrix(0)?.transpose(evaluate)},
            "Adjoint"=>{arity(1)?;matrix(0)?.adjoint(evaluate)},
            "MatrixConjugate"=>{arity(1)?;matrix(0)?.conjugate(evaluate)},
            "MatrixElement"=>{arity(3)?;return Ok(Val::E(matrix(0)?.element(&scalar(1)?,&scalar(2)?).map_err(|e|fail(e.to_string()))?));},
            "Trace"=>{arity(1)?;let parent=matrix(0)?;return Ok(Val::E((if evaluate{parent.trace()}else{parent.trace_unevaluated()}).map_err(|e|fail(e.to_string()))?));},
            "Determinant"=>{arity(1)?;let parent=matrix(0)?;return Ok(Val::E((if evaluate{parent.determinant()}else{parent.determinant_unevaluated()}).map_err(|e|fail(e.to_string()))?));},
            "ImmutableMatrix"|"ImmutableDenseMatrix"|"MutableDenseMatrix"|"Matrix"=>{
                fn flatten(p:&Parser,v:Val,out:&mut Vec<Expr>)->R<()>{match v{Val::Tuple(items)=>{for item in items{flatten(p,item,out)?;}},a=>out.push(p.as_expr(a)?)};Ok(())}
                let (rows,cols,entries)=if args.len()==3 {
                    let rows=scalar(0)?.as_i64().and_then(|n|usize::try_from(n).ok()).ok_or_else(||fail("invalid matrix row count".into()))?;
                    let cols=scalar(1)?.as_i64().and_then(|n|usize::try_from(n).ok()).ok_or_else(||fail("invalid matrix column count".into()))?;
                    let mut entries=Vec::new();flatten(self,args[2].clone(),&mut entries)?;(rows,cols,entries)
                }else if args.len()==1 {
                    let Val::Tuple(rows)=&args[0] else{return self.err("matrix literal expects nested rows");};
                    let cols=match rows.first(){Some(Val::Tuple(row))=>row.len(),None=>0,_=>return self.err("matrix literal expects nested rows")};
                    if rows.iter().any(|row|!matches!(row,Val::Tuple(items) if items.len()==cols)){return self.err("matrix rows have unequal lengths");}
                    let mut entries=Vec::new();flatten(self,args[0].clone(),&mut entries)?;(rows.len(),cols,entries)
                }else{return self.err("matrix literal expects rows, columns, and entries");};
                crate::matrix::Matrix::new(rows,cols,entries).map(MatrixExpr::explicit)
            },
            _=>return self.err("unknown matrix constructor"),
        }.map_err(|e|fail(e.to_string()))?;
        Ok(Val::E(Expr::from_kind(Kind::Func(Func::MatrixProxy(Arc::new(result)),vec![]))))
    }

    fn call(&mut self, name: &str, args: Vec<Val>) -> R<Val> {
        if name == "Str" {
            return match args.as_slice() { [Val::Str(s)]=>Ok(Val::Str(s.clone())), _=>self.err("Str expects one string") };
        }
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
            "UnevaluatedExpr" => Expr::from_kind(Kind::Func(
                Func::Undefined("UnevaluatedExpr".into()), vec![one(self, args)?],
            )),
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
            "Eq" | "Ne" | "Lt" | "Le" | "Gt" | "Ge"
            | "Equality" | "Unequality" | "StrictLessThan" | "LessThan"
            | "StrictGreaterThan" | "GreaterThan" => {
                let v = exprs(self, args)?;
                if v.len() != 2 {
                    return self.err(format!("{name} expects 2 arguments"));
                }
                let op = match name {
                    "Eq" | "Equality" => RelOp::Eq,
                    "Ne" | "Unequality" => RelOp::Ne,
                    "Lt" | "StrictLessThan" => RelOp::Lt,
                    "Le" | "LessThan" => RelOp::Le,
                    "Gt" | "StrictGreaterThan" => RelOp::Gt,
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
                let v = if matches!(name, "hyper" | "bell") {
                    args.into_iter().map(|a| self.tuple_expr(a)).collect::<R<Vec<Expr>>>()?
                } else {
                    exprs(self, args)?
                };
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
        _ => crate::sets::named_set(name).unwrap_or_else(|| Expr::symbol(name)),
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
    let mut p = Parser { toks, pos: 0, kwargs: Vec::new(), precision: None };
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

#[cfg(test)]
mod matrix_parser_tests {
    use super::*;
    #[test]
    fn typed_matrix_srepr_roundtrip() {
        let a=MatrixExpr::symbol("A",Expr::int(2),Expr::int(2)).unwrap();
        let b=MatrixExpr::symbol("B",Expr::int(2),Expr::int(2)).unwrap();
        let explicit=MatrixExpr::explicit(crate::matrix::Matrix::new(2,2,vec![Expr::one(),Expr::int(2),Expr::int(3),Expr::int(4)]).unwrap());
        let values=vec![a.clone(),MatrixExpr::add(&[a.clone(),b.clone()],true).unwrap(),MatrixExpr::mul(Expr::int(2),&[a.clone(),b],true).unwrap(),a.transpose(false).unwrap(),a.inverse(false).unwrap(),explicit];
        for value in values {
            let parsed=parse_str(&value.to_repr()).unwrap();
            assert!(matches!(parsed.kind(),Kind::Func(Func::MatrixProxy(m),_) if **m==value));
        }
        let element=a.element(&Expr::zero(),&Expr::one()).unwrap();
        assert_eq!(parse_str(&crate::printer::srepr(&element)).unwrap(),element);
        let matrix=parse_str("MatrixSymbol(Str('A'), Integer(2), Integer(2))").unwrap();
        assert!(matches!(matrix.kind(),Kind::Func(Func::MatrixProxy(m),_) if **m==a));
        assert!(parse_str("MatMul(MatrixSymbol('A',2,3),MatrixSymbol('B',2,2))").is_err());
        assert!(parse_str("ImmutableDenseMatrix([[Integer(1)]])").is_ok());
        let escaped=MatrixExpr::symbol("A\\\"x\n",Expr::int(1),Expr::int(1)).unwrap();
        assert!(matches!(parse_str(&escaped.to_repr()).unwrap().kind(),Kind::Func(Func::MatrixProxy(m),_) if **m==escaped));
    }
}
