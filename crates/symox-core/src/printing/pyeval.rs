//! Bit-exact native evaluation of a printed lambdify return expression.
//!
//! SymPy's generated function evaluates its printed source with Python float
//! semantics: the printed (15 digit) literals, the printed operation order,
//! `**` as C `pow`, and the module functions of the namespace. A native value
//! of the compiled expression can differ from it in the last digits (folded
//! constants, reordered sums, `powi`), so scalar calls evaluate the printed
//! source itself instead: [`parse`] reads the expression, the caller binds its
//! global names ([`Binding`]) and [`Program::eval`] reproduces the Python
//! value, or gives `None` whenever an operation is not finite (the generated
//! function then raises or warns as SymPy's does).
//!
//! Only operations whose Python result is known bit for bit are accepted:
//! IEEE arithmetic and comparisons, C `pow`, the `math` module functions that
//! call the C library, and the NumPy functions that are exact (`sqrt`,
//! `floor`, comparisons, `select`, ...) or call the C library (`sin`, `cos`).
//! NumPy functions with their own SIMD implementations (`exp`, `log`, `tanh`,
//! `arctan`, ...) differ from the C library on some processors: they are not
//! accepted. Arguments are exact integers below 2**53 or floats.

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Num(f64),
    Int(f64),
    Name(String),
    Op(&'static str),
}

fn tokenize(s: &str) -> Option<Vec<Tok>> {
    const OPS: &[&str] = &["**", "<=", ">=", "==", "!=", "+", "-", "*", "/", "(", ")", "[", "]", ",", ".", "<", ">", "="];
    let b = s.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    while i < b.len() {
        let c = b[i];
        if c == b' ' {
            i += 1;
        } else if c.is_ascii_digit() || (c == b'.' && b.get(i + 1).is_some_and(u8::is_ascii_digit)) {
            let start = i;
            let mut float = false;
            while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
                float |= b[i] == b'.';
                i += 1;
            }
            if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
                float = true;
                i += 1;
                if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
                    i += 1;
                }
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
            }
            if i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == b'_') {
                return None;
            }
            let text = &s[start..i];
            let v = if float {
                text.parse::<f64>().ok()?
            } else {
                // Python integers beyond 2**53 do not keep their value as floats.
                let v = text.parse::<u64>().ok()?;
                if v > 1 << 53 {
                    return None;
                }
                v as f64
            };
            out.push(if float { Tok::Num(v) } else { Tok::Int(v) });
        } else if c.is_ascii_alphabetic() || c == b'_' {
            let start = i;
            while i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == b'_') {
                i += 1;
            }
            out.push(Tok::Name(s[start..i].to_string()));
        } else {
            let op = OPS.iter().find(|op| s[i..].starts_with(**op))?;
            if (*op == "*" || *op == "/") && b.get(i + 1) == Some(&c) {
                return None; // `//` (and `**` was matched first)
            }
            i += op.len();
            out.push(Tok::Op(op));
        }
    }
    Some(out)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Bin {
    Add,
    Sub,
    Mul,
    Div,
    Pow,
}

/// The functions with a known bit-exact Python result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Fun {
    Sin,
    Cos,
    Tan,
    Asin,
    Acos,
    Atan,
    Atan2,
    Sinh,
    Cosh,
    Tanh,
    Asinh,
    Acosh,
    Atanh,
    Exp,
    Expm1,
    Log,
    Log1p,
    Log2,
    Log10,
    Sqrt,
    Abs,
    Copysign,
    Floor,
    Ceil,
    Sign,
    Compare(Cmp),
    LogicalAnd,
    LogicalOr,
    LogicalNot,
    AllOf,
    AnyOf,
    Select,
}

impl Fun {
    fn from_key(key: &str) -> Option<Fun> {
        use Fun::*;
        let (module, name) = key.split_once('.')?;
        let common = match name {
            "sqrt" => Some(Sqrt),
            "sin" => Some(Sin),
            "cos" => Some(Cos),
            _ => None,
        };
        match module {
            "math" => common.or(match name {
                "tan" => Some(Tan),
                "asin" => Some(Asin),
                "acos" => Some(Acos),
                "atan" => Some(Atan),
                "atan2" => Some(Atan2),
                "sinh" => Some(Sinh),
                "cosh" => Some(Cosh),
                "tanh" => Some(Tanh),
                "asinh" => Some(Asinh),
                "acosh" => Some(Acosh),
                "atanh" => Some(Atanh),
                "exp" => Some(Exp),
                "expm1" => Some(Expm1),
                "log" => Some(Log),
                "log1p" => Some(Log1p),
                "log2" => Some(Log2),
                "log10" => Some(Log10),
                "fabs" => Some(Abs),
                "copysign" => Some(Copysign),
                _ => None,
            }),
            "numpy" => common.or(match name {
                "absolute" => Some(Abs),
                "sign" => Some(Sign),
                // `math.floor` gives an int: only NumPy's float result
                "floor" => Some(Floor),
                "ceil" => Some(Ceil),
                "less" => Some(Compare(Cmp::Lt)),
                "less_equal" => Some(Compare(Cmp::Le)),
                "greater" => Some(Compare(Cmp::Gt)),
                "greater_equal" => Some(Compare(Cmp::Ge)),
                "equal" => Some(Compare(Cmp::Eq)),
                "not_equal" => Some(Compare(Cmp::Ne)),
                "logical_and" => Some(LogicalAnd),
                "logical_or" => Some(LogicalOr),
                "logical_not" => Some(LogicalNot),
                "select" => Some(Select),
                _ => None,
            }),
            "builtins" if name == "abs" => Some(Abs),
            _ => None,
        }
    }

    fn arity(self) -> usize {
        use Fun::*;
        match self {
            Atan2 | Copysign | Compare(_) | LogicalAnd | LogicalOr => 2,
            AllOf | AnyOf | Select => 0,
            _ => 1,
        }
    }

    fn boolean(self) -> bool {
        matches!(self, Fun::Compare(_) | Fun::LogicalAnd | Fun::LogicalOr | Fun::LogicalNot | Fun::AllOf | Fun::AnyOf)
    }
}

#[derive(Clone, Debug)]
pub(crate) enum Node {
    Num(f64),
    Int(f64),
    Bool(bool),
    Arg(usize),
    Global(usize),
    Attr(usize, String),
    Neg(Box<Node>),
    Pos(Box<Node>),
    Bin(Bin, Box<Node>, Box<Node>),
    Cmp(Cmp, Box<Node>, Box<Node>),
    And(Vec<Node>),
    Or(Vec<Node>),
    Not(Box<Node>),
    IfElse(Box<Node>, Box<Node>, Box<Node>),
    Seq(Vec<Node>),
    Call(Box<Node>, Vec<Node>, Vec<(String, Node)>),
    Apply(Fun, Vec<Node>),
}

struct Parser<'a> {
    toks: Vec<Tok>,
    pos: usize,
    args: &'a [String],
    globals: Vec<String>,
}

impl Parser<'_> {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos)
    }
    fn is_op(&self, op: &str) -> bool {
        matches!(self.peek(), Some(Tok::Op(o)) if *o == op)
    }
    fn is_name(&self, name: &str) -> bool {
        matches!(self.peek(), Some(Tok::Name(n)) if n == name)
    }
    fn eat_op(&mut self, op: &str) -> bool {
        let found = self.is_op(op);
        self.pos += found as usize;
        found
    }
    fn eat_name(&mut self, name: &str) -> bool {
        let found = self.is_name(name);
        self.pos += found as usize;
        found
    }
    fn expect(&mut self, op: &str) -> Option<()> {
        self.eat_op(op).then_some(())
    }

    fn expr(&mut self) -> Option<Node> {
        let body = self.or_test()?;
        if self.eat_name("if") {
            let cond = self.or_test()?;
            if !self.eat_name("else") {
                return None;
            }
            let other = self.expr()?;
            return Some(Node::IfElse(Box::new(cond), Box::new(body), Box::new(other)));
        }
        Some(body)
    }

    fn or_test(&mut self) -> Option<Node> {
        let mut items = vec![self.and_test()?];
        while self.eat_name("or") {
            items.push(self.and_test()?);
        }
        Some(if items.len() == 1 { items.pop()? } else { Node::Or(items) })
    }

    fn and_test(&mut self) -> Option<Node> {
        let mut items = vec![self.not_test()?];
        while self.eat_name("and") {
            items.push(self.not_test()?);
        }
        Some(if items.len() == 1 { items.pop()? } else { Node::And(items) })
    }

    fn not_test(&mut self) -> Option<Node> {
        if self.eat_name("not") {
            return Some(Node::Not(Box::new(self.not_test()?)));
        }
        self.comparison()
    }

    fn comparison(&mut self) -> Option<Node> {
        let left = self.arith()?;
        let op = match self.peek() {
            Some(Tok::Op("<")) => Cmp::Lt,
            Some(Tok::Op("<=")) => Cmp::Le,
            Some(Tok::Op(">")) => Cmp::Gt,
            Some(Tok::Op(">=")) => Cmp::Ge,
            Some(Tok::Op("==")) => Cmp::Eq,
            Some(Tok::Op("!=")) => Cmp::Ne,
            Some(Tok::Name(n)) if n == "is" || n == "in" => return None,
            _ => return Some(left),
        };
        self.pos += 1;
        let right = self.arith()?;
        // chained comparisons are not printed
        if matches!(self.peek(), Some(Tok::Op("<" | "<=" | ">" | ">=" | "==" | "!="))) {
            return None;
        }
        Some(Node::Cmp(op, Box::new(left), Box::new(right)))
    }

    fn arith(&mut self) -> Option<Node> {
        let mut left = self.term()?;
        loop {
            let op = if self.eat_op("+") { Bin::Add } else if self.eat_op("-") { Bin::Sub } else { return Some(left) };
            left = Node::Bin(op, Box::new(left), Box::new(self.term()?));
        }
    }

    fn term(&mut self) -> Option<Node> {
        let mut left = self.factor()?;
        loop {
            let op = if self.eat_op("*") { Bin::Mul } else if self.eat_op("/") { Bin::Div } else { return Some(left) };
            left = Node::Bin(op, Box::new(left), Box::new(self.factor()?));
        }
    }

    fn factor(&mut self) -> Option<Node> {
        if self.eat_op("-") {
            return Some(Node::Neg(Box::new(self.factor()?)));
        }
        if self.eat_op("+") {
            return Some(Node::Pos(Box::new(self.factor()?)));
        }
        let base = self.trailers()?;
        if self.eat_op("**") {
            return Some(Node::Bin(Bin::Pow, Box::new(base), Box::new(self.factor()?)));
        }
        Some(base)
    }

    fn trailers(&mut self) -> Option<Node> {
        let mut node = self.atom()?;
        loop {
            if self.eat_op(".") {
                let Some(Tok::Name(attr)) = self.peek().cloned() else { return None };
                self.pos += 1;
                node = match node {
                    Node::Global(g) => Node::Attr(g, attr),
                    _ => return None,
                };
            } else if self.eat_op("(") {
                let mut args = Vec::new();
                let mut kwargs = Vec::new();
                while !self.eat_op(")") {
                    if let (Some(Tok::Name(k)), Some(Tok::Op("="))) = (self.peek().cloned(), self.toks.get(self.pos + 1).cloned()) {
                        self.pos += 2;
                        kwargs.push((k, self.expr()?));
                    } else if kwargs.is_empty() {
                        args.push(self.expr()?);
                    } else {
                        return None;
                    }
                    if !self.eat_op(",") {
                        self.expect(")")?;
                        break;
                    }
                }
                node = Node::Call(Box::new(node), args, kwargs);
            } else {
                return Some(node);
            }
        }
    }

    fn sequence(&mut self, close: &str) -> Option<(Vec<Node>, bool)> {
        let mut items = Vec::new();
        let mut comma = false;
        while !self.eat_op(close) {
            items.push(self.expr()?);
            if self.eat_op(",") {
                comma = true;
            } else {
                self.expect(close)?;
                break;
            }
        }
        Some((items, comma))
    }

    fn atom(&mut self) -> Option<Node> {
        match self.peek().cloned()? {
            Tok::Num(v) => {
                self.pos += 1;
                Some(Node::Num(v))
            }
            Tok::Int(v) => {
                self.pos += 1;
                Some(Node::Int(v))
            }
            Tok::Name(n) => {
                self.pos += 1;
                const KEYWORDS: &[&str] = &["if", "else", "and", "or", "not", "is", "in", "lambda", "None", "for"];
                if KEYWORDS.contains(&n.as_str()) {
                    return None;
                }
                if n == "True" || n == "False" {
                    return Some(Node::Bool(n == "True"));
                }
                if let Some(i) = self.args.iter().position(|a| *a == n) {
                    return Some(Node::Arg(i));
                }
                let g = match self.globals.iter().position(|a| *a == n) {
                    Some(g) => g,
                    None => {
                        self.globals.push(n);
                        self.globals.len() - 1
                    }
                };
                Some(Node::Global(g))
            }
            Tok::Op("(") => {
                self.pos += 1;
                let (mut items, comma) = self.sequence(")")?;
                if items.len() == 1 && !comma {
                    items.pop()
                } else {
                    Some(Node::Seq(items))
                }
            }
            Tok::Op("[") => {
                self.pos += 1;
                Some(Node::Seq(self.sequence("]")?.0))
            }
            Tok::Op(_) => None,
        }
    }
}

/// A printed return expression with its unresolved global names.
#[derive(Clone, Debug)]
pub struct Parsed {
    root: Node,
    globals: Vec<String>,
}

/// What a global name of the generated function's namespace is.
#[derive(Clone, Debug)]
pub enum Binding {
    /// A Python `float` value (`pi`, `e`, `nan`, ...).
    Float(f64),
    /// A known module function: `"math.sin"`, `"numpy.select"`, `"builtins.abs"`.
    Function(String),
    /// Anything else: the expression is not evaluated natively.
    Other,
}

/// Parse a printed return expression whose argument names are `args`.
pub fn parse(body: &str, args: &[String]) -> Option<Parsed> {
    let mut p = Parser { toks: tokenize(body)?, pos: 0, args, globals: Vec::new() };
    let root = p.expr()?;
    if p.pos != p.toks.len() {
        return None;
    }
    Some(Parsed { root, globals: p.globals })
}

impl Parsed {
    /// The global names, in the order [`Parsed::resolve`] expects their bindings.
    pub fn globals(&self) -> &[String] {
        &self.globals
    }

    /// The parsed tree (read by `crate::lambdify_vec`).
    pub(crate) fn root(&self) -> &Node {
        &self.root
    }

    /// The program once every global name is bound, or `None` when an
    /// operation has no bit-exact native evaluation.
    pub fn resolve(&self, bindings: &[Binding]) -> Option<Program> {
        if bindings.len() != self.globals.len() {
            return None;
        }
        let root = lower(&self.root, bindings)?;
        Some(Program { root })
    }
}

fn function(node: &Node, bindings: &[Binding]) -> Option<Fun> {
    match node {
        Node::Global(g) => match &bindings[*g] {
            Binding::Function(key) => Fun::from_key(key),
            _ => None,
        },
        Node::Attr(g, attr) => match (&bindings[*g], attr.as_str()) {
            (Binding::Function(key), "reduce") => match Fun::from_key(key)? {
                Fun::LogicalAnd => Some(Fun::AllOf),
                Fun::LogicalOr => Some(Fun::AnyOf),
                _ => None,
            },
            _ => None,
        },
        _ => None,
    }
}

fn boolean(node: &Node) -> bool {
    match node {
        Node::Bool(_) | Node::Cmp(..) | Node::Not(_) => true,
        Node::And(items) | Node::Or(items) => items.iter().all(boolean),
        Node::Apply(f, _) => f.boolean(),
        _ => false,
    }
}

fn lower(node: &Node, bindings: &[Binding]) -> Option<Node> {
    let all = |items: &[Node]| items.iter().map(|n| lower(n, bindings)).collect::<Option<Vec<_>>>();
    let one = |n: &Node| lower(n, bindings).map(Box::new);
    Some(match node {
        Node::Num(v) => Node::Num(*v),
        Node::Int(v) => Node::Int(*v),
        Node::Bool(b) => Node::Bool(*b),
        Node::Arg(i) => Node::Arg(*i),
        Node::Global(g) => match &bindings[*g] {
            Binding::Float(v) => Node::Num(*v),
            _ => return None,
        },
        Node::Neg(a) => Node::Neg(one(a)?),
        Node::Pos(a) => Node::Pos(one(a)?),
        Node::Bin(op, a, b) => Node::Bin(*op, one(a)?, one(b)?),
        Node::Cmp(op, a, b) => Node::Cmp(*op, one(a)?, one(b)?),
        // Python's `and`/`or` return an operand: only boolean operands keep
        // a boolean value
        Node::And(items) | Node::Or(items) => {
            let items = all(items)?;
            if !items.iter().all(boolean) {
                return None;
            }
            if matches!(node, Node::And(_)) { Node::And(items) } else { Node::Or(items) }
        }
        Node::Not(a) => Node::Not(one(a)?),
        Node::IfElse(c, a, b) => Node::IfElse(one(c)?, one(a)?, one(b)?),
        Node::Call(callee, args, kwargs) => {
            let f = function(callee, bindings)?;
            let args = all(args)?;
            match f {
                Fun::Select => {
                    // select(conditions, choices, default=value)
                    let [Node::Seq(conds), Node::Seq(choices)] = args.as_slice() else { return None };
                    if conds.len() != choices.len() || conds.is_empty() || !conds.iter().all(boolean)
                        || choices.iter().any(|c| matches!(c, Node::Seq(_))) {
                        return None;
                    }
                    let default = match kwargs.as_slice() {
                        [(k, v)] if k == "default" => lower(v, bindings)?,
                        [] => Node::Int(0.0),
                        _ => return None,
                    };
                    let mut items = conds.clone();
                    items.extend(choices.iter().cloned());
                    items.push(default);
                    Node::Apply(Fun::Select, items)
                }
                Fun::AllOf | Fun::AnyOf => {
                    let [Node::Seq(items)] = args.as_slice() else { return None };
                    if !kwargs.is_empty() || items.is_empty() || !items.iter().all(boolean) {
                        return None;
                    }
                    Node::Apply(f, items.clone())
                }
                _ => {
                    if !kwargs.is_empty() || args.len() != f.arity() || args.iter().any(|a| matches!(a, Node::Seq(_))) {
                        return None;
                    }
                    if matches!(f, Fun::LogicalAnd | Fun::LogicalOr | Fun::LogicalNot) && !args.iter().all(boolean) {
                        return None;
                    }
                    Node::Apply(f, args)
                }
            }
        }
        // a list or tuple argument (select, reduce)
        Node::Seq(items) => Node::Seq(all(items)?),
        Node::Attr(..) | Node::Apply(..) => return None,
    })
}

/// A printed return expression with a bit-exact native evaluation.
#[derive(Clone, Debug)]
pub struct Program {
    root: Node,
}

mod libm {
    extern "C" {
        pub fn sin(x: f64) -> f64;
        pub fn cos(x: f64) -> f64;
        pub fn tan(x: f64) -> f64;
        pub fn asin(x: f64) -> f64;
        pub fn acos(x: f64) -> f64;
        pub fn atan(x: f64) -> f64;
        pub fn atan2(y: f64, x: f64) -> f64;
        pub fn sinh(x: f64) -> f64;
        pub fn cosh(x: f64) -> f64;
        pub fn tanh(x: f64) -> f64;
        pub fn asinh(x: f64) -> f64;
        pub fn acosh(x: f64) -> f64;
        pub fn atanh(x: f64) -> f64;
        pub fn exp(x: f64) -> f64;
        pub fn expm1(x: f64) -> f64;
        pub fn log(x: f64) -> f64;
        pub fn log1p(x: f64) -> f64;
        pub fn log2(x: f64) -> f64;
        pub fn log10(x: f64) -> f64;
        pub fn pow(x: f64, y: f64) -> f64;
    }
}

/// Not finite: the generated function decides (it raises or warns).
struct Stop;

/// A Python value: a float, or an exact integer (or bool) below 2**53.
#[derive(Clone, Copy)]
struct V {
    v: f64,
    int: bool,
}

const EXACT: f64 = 9007199254740992.0;

fn float(v: f64) -> Result<V, Stop> {
    if v.is_finite() { Ok(V { v, int: false }) } else { Err(Stop) }
}

fn int(v: f64) -> Result<V, Stop> {
    // Python integers have no signed zero
    if v.abs() <= EXACT { Ok(V { v: v + 0.0, int: true }) } else { Err(Stop) }
}

fn truth(b: bool) -> V {
    V { v: if b { 1.0 } else { 0.0 }, int: true }
}

fn compare(op: Cmp, a: f64, b: f64) -> V {
    truth(match op {
        Cmp::Lt => a < b,
        Cmp::Le => a <= b,
        Cmp::Gt => a > b,
        Cmp::Ge => a >= b,
        Cmp::Eq => a == b,
        Cmp::Ne => a != b,
    })
}

impl Program {
    /// The value of the generated function at finite float arguments (or
    /// exact integers below 2**53, flagged in `ints`), `None` when it is not
    /// finite somewhere.
    pub fn eval(&self, args: &[f64], ints: &[bool]) -> Option<f64> {
        if args.len() != ints.len() || !args.iter().zip(ints).all(|(a, i)| a.is_finite() && (!i || a.abs() <= EXACT)) {
            return None;
        }
        Self::go(&self.root, args, ints).ok().map(|v| v.v)
    }

    fn go(node: &Node, args: &[f64], ints: &[bool]) -> Result<V, Stop> {
        let go = |n: &Node| Self::go(n, args, ints);
        match node {
            Node::Num(v) => Ok(V { v: *v, int: false }),
            Node::Int(v) => Ok(V { v: *v, int: true }),
            Node::Bool(b) => Ok(truth(*b)),
            Node::Arg(i) => Ok(V { v: *args.get(*i).ok_or(Stop)?, int: ints[*i] }),
            Node::Neg(a) => {
                let a = go(a)?;
                if a.int { int(-a.v) } else { Ok(V { v: -a.v, int: false }) }
            }
            Node::Pos(a) => go(a),
            Node::Bin(op, a, b) => {
                let (a, b) = (go(a)?, go(b)?);
                let both = a.int && b.int;
                match op {
                    Bin::Add if both => int(a.v + b.v),
                    Bin::Sub if both => int(a.v - b.v),
                    Bin::Mul if both => int(a.v * b.v),
                    Bin::Add => float(a.v + b.v),
                    Bin::Sub => float(a.v - b.v),
                    Bin::Mul => float(a.v * b.v),
                    // int / int is correctly rounded, as the float division
                    Bin::Div => float(a.v / b.v),
                    // Python's int ** non-negative int is exact; the float
                    // `**` (and NumPy's scalar power) is C pow
                    Bin::Pow if both && b.v >= 0.0 => int(unsafe { libm::pow(a.v, b.v) }),
                    Bin::Pow => float(unsafe { libm::pow(a.v, b.v) }),
                }
            }
            Node::Cmp(op, a, b) => Ok(compare(*op, go(a)?.v, go(b)?.v)),
            Node::And(items) => {
                for item in items {
                    if go(item)?.v == 0.0 {
                        return Ok(truth(false));
                    }
                }
                Ok(truth(true))
            }
            Node::Or(items) => {
                for item in items {
                    if go(item)?.v != 0.0 {
                        return Ok(truth(true));
                    }
                }
                Ok(truth(false))
            }
            Node::Not(a) => Ok(truth(go(a)?.v == 0.0)),
            Node::IfElse(c, a, b) => {
                if go(c)?.v != 0.0 { go(a) } else { go(b) }
            }
            Node::Apply(f, items) => Self::apply(*f, items, args, ints),
            Node::Global(_) | Node::Attr(..) | Node::Seq(_) | Node::Call(..) => Err(Stop),
        }
    }

    fn apply(f: Fun, items: &[Node], args: &[f64], ints: &[bool]) -> Result<V, Stop> {
        let go = |n: &Node| Self::go(n, args, ints);
        match f {
            Fun::AllOf | Fun::AnyOf => {
                // the tuple is evaluated before the reduction
                let (mut all, mut any) = (true, false);
                for item in items {
                    let v = go(item)?.v != 0.0;
                    all &= v;
                    any |= v;
                }
                return Ok(truth(if f == Fun::AllOf { all } else { any }));
            }
            Fun::Select => {
                // every condition and choice is evaluated
                let n = (items.len() - 1) / 2;
                let mut values = Vec::with_capacity(items.len());
                for item in &items[..2 * n] {
                    values.push(go(item)?);
                }
                let default = go(&items[2 * n])?;
                let chosen = (0..n).find(|&i| values[i].v != 0.0).map_or(default, |i| values[n + i]);
                // a NaN default is the printed `nan`, not an error
                return float(chosen.v);
            }
            _ => {}
        }
        // NumPy functions evaluate every argument
        let a = go(&items[0])?;
        let b = if f.arity() == 2 { go(&items[1])? } else { a };
        let (x, y) = (a.v, b.v);
        Ok(match f {
            Fun::Compare(op) => compare(op, x, y),
            Fun::LogicalAnd => truth(x != 0.0 && y != 0.0),
            Fun::LogicalOr => truth(x != 0.0 || y != 0.0),
            Fun::LogicalNot => truth(x == 0.0),
            // abs and numpy.sign of an int are ints
            Fun::Abs if a.int => int(x.abs())?,
            Fun::Sign if a.int => int(x.signum() * (x != 0.0) as u8 as f64)?,
            _ => float(unsafe {
                match f {
                    Fun::Sin => libm::sin(x),
                    Fun::Cos => libm::cos(x),
                    Fun::Tan => libm::tan(x),
                    Fun::Asin => libm::asin(x),
                    Fun::Acos => libm::acos(x),
                    Fun::Atan => libm::atan(x),
                    Fun::Atan2 => libm::atan2(x, y),
                    Fun::Sinh => libm::sinh(x),
                    Fun::Cosh => libm::cosh(x),
                    Fun::Tanh => libm::tanh(x),
                    Fun::Asinh => libm::asinh(x),
                    Fun::Acosh => libm::acosh(x),
                    Fun::Atanh => libm::atanh(x),
                    Fun::Exp => libm::exp(x),
                    Fun::Expm1 => libm::expm1(x),
                    // math.log of a positive finite value; others raise
                    Fun::Log => if x > 0.0 { libm::log(x) } else { f64::NAN },
                    Fun::Log1p => if x > -1.0 { libm::log1p(x) } else { f64::NAN },
                    Fun::Log2 => if x > 0.0 { libm::log2(x) } else { f64::NAN },
                    Fun::Log10 => if x > 0.0 { libm::log10(x) } else { f64::NAN },
                    Fun::Sqrt => if x >= 0.0 { x.sqrt() } else { f64::NAN },
                    Fun::Abs => x.abs(),
                    Fun::Copysign => x.copysign(y),
                    Fun::Floor => x.floor(),
                    Fun::Ceil => x.ceil(),
                    Fun::Sign => if x > 0.0 { 1.0 } else if x < 0.0 { -1.0 } else { 0.0 },
                    _ => f64::NAN,
                }
            })?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn program(body: &str, args: &[&str], bindings: &[(&str, Binding)]) -> Option<Program> {
        let args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        let parsed = parse(body, &args)?;
        let b: Vec<Binding> = parsed.globals().iter()
            .map(|g| bindings.iter().find(|(n, _)| n == g).map_or(Binding::Other, |(_, b)| b.clone()))
            .collect();
        parsed.resolve(&b)
    }

    #[test]
    fn python_float_semantics() {
        let f = |k: &str| Binding::Function(k.to_string());
        let p = program("15*x + 1.21609067478396", &["x"], &[]).unwrap();
        assert_eq!(p.eval(&[0.3], &[false]), Some(15.0 * 0.3 + 1.21609067478396));
        let p = program("-x**2 + sin(x)", &["x"], &[("sin", f("math.sin"))]).unwrap();
        assert_eq!(p.eval(&[0.5], &[false]), Some(-unsafe { libm::pow(0.5, 2.0) } + 0.5f64.sin()));
        let p = program("((x**2) if (x < 1) else (x))", &["x"], &[]).unwrap();
        assert_eq!(p.eval(&[3.0], &[false]), Some(3.0));
        let p = program("select([less(x, 1),True], [x**2,x], default=nan)", &["x"],
            &[("select", f("numpy.select")), ("less", f("numpy.less")), ("nan", Binding::Float(f64::NAN))]).unwrap();
        assert_eq!(p.eval(&[0.5], &[false]), Some(0.25));
        assert_eq!(p.eval(&[2.0], &[false]), Some(2.0));
        let p = program("logical_and.reduce((greater(x, 1),less(y, 2)))", &["x", "y"],
            &[("logical_and", f("numpy.logical_and")), ("greater", f("numpy.greater")), ("less", f("numpy.less"))]).unwrap();
        assert_eq!(p.eval(&[2.0, 1.0], &[false, false]), Some(1.0));
        assert_eq!(p.eval(&[0.0, 1.0], &[false, false]), Some(0.0));
        // errors and non-finite values stop
        let p = program("1/x", &["x"], &[]).unwrap();
        assert_eq!(p.eval(&[0.0], &[false]), None);
        assert!(program("exp(x)", &["x"], &[("exp", f("numpy.exp"))]).is_none());
        assert!(program("gamma(x)", &["x"], &[("gamma", f("math.gamma"))]).is_none());
        assert!(program("x // 2", &["x"], &[]).is_none());
        assert!(program("x and 1", &["x"], &[]).is_none());
        assert!(program("100000000000000000000*x", &["x"], &[]).is_none());
        // Python integers: exact, without a signed zero
        let p = program("-x*1.5", &["x"], &[]).unwrap();
        assert_eq!(p.eval(&[0.0], &[true]).map(f64::to_bits), Some(0.0f64.to_bits()));
        assert_eq!(p.eval(&[0.0], &[false]).map(f64::to_bits), Some((-0.0f64).to_bits()));
        let p = program("x**3 + 1/x", &["x"], &[]).unwrap();
        assert_eq!(p.eval(&[3.0], &[true]), Some(27.0 + 1.0 / 3.0));
    }
}
