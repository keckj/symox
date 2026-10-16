//! Code printers: C99, C++11, CUDA, Python (`math`), NumPy, Rust and a
//! basic Fortran printer.
//!
//! The C, C++, Python, NumPy and Rust printers reproduce the output of
//! SymPy's `ccode`, `cxxcode`, `pycode`, `NumPyPrinter` and `rust_code`.
//! `cudacode` is a C99 dialect for CUDA device code (`float`/`double`,
//! optional fast-math intrinsics, inlined small integer powers).

use std::cell::RefCell;
use std::fmt;

use rustc_hash::FxHashMap;

use super::{
    float_str, mpf_to_str, ordered_factors, precedence, PREC_ADD, PREC_ATOM, PREC_MUL, PREC_POW,
};
use crate::expr::{Constant, Expr, Func, Kind, RelOp};
use crate::number::Number;
use crate::ops::{div, mul_factor, neg};

/// Target language of a [`CodePrinter`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Lang {
    C,
    Cxx,
    Cuda,
    Python,
    NumPy,
    Rust,
    Fortran,
}

/// Floating point type of the generated code (C family).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum FloatType {
    F32,
    #[default]
    F64,
}

impl FloatType {
    pub fn c_name(self) -> &'static str {
        match self {
            FloatType::F32 => "float",
            FloatType::F64 => "double",
        }
    }
    fn digits(self) -> usize {
        match self {
            FloatType::F32 => 9,
            FloatType::F64 => 17,
        }
    }
    fn literal_suffix(self) -> &'static str {
        match self {
            FloatType::F32 => "F",
            FloatType::F64 => "",
        }
    }
    fn func_suffix(self) -> &'static str {
        match self {
            FloatType::F32 => "f",
            FloatType::F64 => "",
        }
    }
}

/// Error raised for expressions a language cannot express.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CodeError {
    /// Function or object without a translation in the target language.
    Unsupported(String),
    /// A `Piecewise` without a final `True` condition.
    PiecewiseNoDefault,
}

impl fmt::Display for CodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CodeError::Unsupported(s) => write!(f, "Unsupported by the code printer: {s}"),
            CodeError::PiecewiseNoDefault => write!(
                f,
                "All Piecewise expressions must contain an (expr, True) statement to be used as a default condition."
            ),
        }
    }
}

impl std::error::Error for CodeError {}

pub type CodeResult = Result<String, CodeError>;

/// Settings of the generic code printer.
#[derive(Clone, Debug)]
pub struct CodeSettings {
    pub lang: Lang,
    /// Floating point type (C, C++ and CUDA).
    pub real: FloatType,
    /// CUDA: use `__sinf`-style intrinsics (single precision only).
    pub fast_math: bool,
    /// CUDA: print `x**(-1/2)` as `rsqrt(x)`.
    pub use_rsqrt: bool,
    /// CUDA: integer powers of symbols up to this order are written as
    /// products (`x*x*x`). 0 disables.
    pub inline_pow: u32,
    /// C/C++: print `M_PI`, `M_SQRT2`, ... for the known constants.
    pub math_macros: bool,
    /// Extra function name translations (`"f" -> "my_f"`).
    pub user_functions: FxHashMap<String, String>,
    /// Print unknown functions as `name(args)` instead of failing.
    pub allow_unknown_functions: bool,
    /// Fortran standard (66, 77, 90, 95, 2003, 2008); SymPy's default 77.
    pub fortran_standard: u32,
    /// Fortran free source form (default: fixed form, 6 column lead).
    pub fortran_free_form: bool,
}

impl CodeSettings {
    pub fn new(lang: Lang) -> Self {
        CodeSettings {
            lang,
            real: FloatType::F64,
            fast_math: false,
            use_rsqrt: false,
            inline_pow: if lang == Lang::Cuda { 4 } else { 0 },
            math_macros: matches!(lang, Lang::C | Lang::Cxx | Lang::Cuda),
            user_functions: FxHashMap::default(),
            allow_unknown_functions: false,
            fortran_standard: 77,
            fortran_free_form: false,
        }
    }
}

/// Settings of [`cudacode_with`].
#[derive(Clone, Debug)]
pub struct CudaSettings {
    pub dtype: FloatType,
    pub fast_math: bool,
    pub use_rsqrt: bool,
    pub inline_pow: u32,
}

impl Default for CudaSettings {
    fn default() -> Self {
        CudaSettings { dtype: FloatType::F64, fast_math: false, use_rsqrt: true, inline_pow: 4 }
    }
}

impl CudaSettings {
    pub fn to_code_settings(&self) -> CodeSettings {
        let mut s = CodeSettings::new(Lang::Cuda);
        s.real = self.dtype;
        s.fast_math = self.fast_math && self.dtype == FloatType::F32;
        s.use_rsqrt = self.use_rsqrt;
        s.inline_pow = self.inline_pow;
        s.math_macros = self.dtype == FloatType::F64;
        s
    }
}

/// SymPy's `ccode(expr)` (C99).
pub fn ccode(e: &Expr) -> CodeResult {
    code_with(e, &CodeSettings::new(Lang::C))
}

/// SymPy's `cxxcode(expr)` (C++11).
pub fn cxxcode(e: &Expr) -> CodeResult {
    code_with(e, &CodeSettings::new(Lang::Cxx))
}

/// SymPy's `pycode(expr)` (Python with the `math` module).
pub fn pycode(e: &Expr) -> CodeResult {
    code_with(e, &CodeSettings::new(Lang::Python))
}

/// SymPy's `NumPyPrinter().doprint(expr)`.
pub fn numpycode(e: &Expr) -> CodeResult {
    code_with(e, &CodeSettings::new(Lang::NumPy))
}

/// SymPy's `rust_code(expr)`.
pub fn rust_code(e: &Expr) -> CodeResult {
    code_with(e, &CodeSettings::new(Lang::Rust))
}

/// SymPy's `fcode(expr)` (Fortran 77, fixed source form).
pub fn fcode(e: &Expr) -> CodeResult {
    code_with(e, &CodeSettings::new(Lang::Fortran))
}

/// CUDA device code with default settings (double precision).
pub fn cudacode(e: &Expr) -> CodeResult {
    cudacode_with(e, &CudaSettings::default())
}

/// CUDA device code.
pub fn cudacode_with(e: &Expr, s: &CudaSettings) -> CodeResult {
    code_with(e, &s.to_code_settings())
}

/// Print `e` with arbitrary settings.
pub fn code_with(e: &Expr, s: &CodeSettings) -> CodeResult {
    let _session = super::order::session();
    let p = CodePrinter::new(s);
    let body = p.print_top(e)?;
    Ok(p.finish(body))
}

// ----------------------------------------------------------------------
// reserved words
// ----------------------------------------------------------------------

const C_RESERVED: &[&str] = &[
    "auto", "break", "case", "char", "const", "continue", "default", "do", "double", "else",
    "enum", "extern", "float", "for", "goto", "if", "int", "long", "register", "return", "short",
    "signed", "sizeof", "static", "struct", "entry", "switch", "typedef", "union", "unsigned",
    "void", "volatile", "while", "inline", "restrict",
];

const CXX_RESERVED: &[&str] = &[
    "and", "and_eq", "asm", "auto", "bitand", "bitor", "bool", "break", "case", "catch,", "char",
    "class", "compl", "const", "const_cast", "continue", "default", "delete", "do", "double",
    "dynamic_cast", "else", "enum", "explicit", "export", "extern", "false", "float", "for",
    "friend", "goto", "if", "inline", "int", "long", "mutable", "namespace", "new", "not",
    "not_eq", "operator", "or", "or_eq", "private", "protected", "public", "register",
    "reinterpret_cast", "return", "short", "signed", "sizeof", "static", "static_cast", "struct",
    "switch", "template", "this", "throw", "true", "try", "typedef", "typeid", "typename", "union",
    "unsigned", "using", "virtual", "void", "volatile", "wchar_t", "while", "xor", "xor_eq",
    "alignas", "alignof", "char16_t", "char32_t", "constexpr", "decltype", "noexcept", "nullptr",
    "static_assert", "thread_local",
];

const CUDA_EXTRA_RESERVED: &[&str] = &[
    "__global__", "__device__", "__host__", "__shared__", "__constant__", "threadIdx", "blockIdx",
    "blockDim", "gridDim", "warpSize",
];

const PY_RESERVED: &[&str] = &[
    "and", "as", "assert", "break", "class", "continue", "def", "del", "elif", "else", "except",
    "finally", "for", "from", "global", "if", "import", "in", "is", "lambda", "not", "or", "pass",
    "raise", "return", "try", "while", "with", "yield", "None", "False", "nonlocal", "True",
];

const RUST_RESERVED: &[&str] = &[
    "abstract", "as", "async", "await", "become", "box", "break", "const", "continue", "crate",
    "do", "dyn", "else", "enum", "extern", "false", "final", "fn", "for", "gen", "if", "impl", "in",
    "let", "loop", "macro", "match", "mod", "move", "mut", "override", "priv", "pub", "ref",
    "return", "Self", "self", "static", "struct", "super", "trait", "true", "try", "type",
    "typeof", "unsafe", "unsized", "use", "virtual", "where", "while", "yield",
];

/// Is `name` a reserved word of `lang`?
pub fn is_reserved(lang: Lang, name: &str) -> bool {
    match lang {
        Lang::C => C_RESERVED.contains(&name),
        Lang::Cxx => CXX_RESERVED.contains(&name),
        Lang::Cuda => {
            C_RESERVED.contains(&name) || CUDA_EXTRA_RESERVED.contains(&name)
        }
        Lang::Python | Lang::NumPy => PY_RESERVED.contains(&name),
        Lang::Rust => RUST_RESERVED.contains(&name),
        Lang::Fortran => false,
    }
}

// ----------------------------------------------------------------------
// the printer
// ----------------------------------------------------------------------

/// Generic code printer (StrPrinter/CodePrinter semantics).
pub struct CodePrinter<'a> {
    s: &'a CodeSettings,
    macros: Vec<(Expr, &'static str)>,
    number_symbols: RefCell<Vec<(&'static str, f64)>>,
    /// Fortran: bit patterns of the floats produced by `N(..., 17)`
    n17: RefCell<rustc_hash::FxHashSet<u64>>,
    /// Fortran: memo of `N(e)`
    n_memo: RefCell<FxHashMap<Expr, Expr>>,
}

fn is_half(x: &Expr) -> bool {
    *x == Expr::rational(1, 2) || matches!(x.kind(), Kind::Num(Number::Float(f)) if *f == 0.5)
}

fn is_minus_half(x: &Expr) -> bool {
    *x == Expr::rational(-1, 2) || matches!(x.kind(), Kind::Num(Number::Float(f)) if *f == -0.5)
}

fn is_minus_one_valued(x: &Expr) -> bool {
    x.is_minus_one() || matches!(x.kind(), Kind::Num(Number::Float(f)) if *f == -1.0)
}

fn is_third(x: &Expr) -> bool {
    *x == Expr::rational(1, 3)
}

/// Negative rational exponent (moves the factor into a denominator).
fn neg_rational(x: &Expr) -> bool {
    matches!(x.kind(), Kind::Num(n) if n.is_rational() && n.is_negative())
}

/// Real valued numeric expression that is not an integer (SymPy's
/// `is_real and not is_integer` without assumptions).
fn real_non_integer_number(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => !n.is_integer(),
        Kind::Const(c) => matches!(c, Constant::Pi | Constant::E | Constant::EulerGamma),
        _ => {
            super::is_numeric(e)
                && !crate::functions::contains_i(e)
                && !e.is_integer_number()
        }
    }
}

enum Den {
    /// `base ** exp` with a positive exponent (unevaluated)
    Pow(Expr, Expr),
    /// `base` (from `base**-1`)
    Base(Expr),
}

impl<'a> CodePrinter<'a> {
    pub fn new(s: &'a CodeSettings) -> Self {
        let macros = if s.math_macros { math_macros() } else { vec![] };
        CodePrinter { s, macros, number_symbols: RefCell::new(Vec::new()), n17: RefCell::new(Default::default()), n_memo: RefCell::new(Default::default()) }
    }

    fn lang(&self) -> Lang {
        self.s.lang
    }

    /// Final assembly: constant declarations and indentation.
    pub fn finish(&self, body: String) -> String {
        let mut lines = self.declarations();
        lines.push(body);
        if self.lang() == Lang::Fortran {
            return fortran_wrap(&lines, self.s.fortran_free_form).join("\n");
        }
        let code = lines.join("\n");
        match self.lang() {
            Lang::C | Lang::Cxx | Lang::Cuda => indent_code(&code, "   "),
            Lang::Rust => indent_code(&code, "    "),
            _ => code,
        }
    }

    /// Constant declarations (`const double EulerGamma = ...;`) needed by
    /// the code printed so far, one per line.
    pub fn declarations(&self) -> Vec<String> {
        let mut syms = self.number_symbols.borrow().clone();
        syms.sort_by(|a, b| a.0.cmp(b.0));
        syms.dedup_by(|a, b| a.0 == b.0);
        syms.into_iter()
            .filter_map(|(name, v)| match self.lang() {
                Lang::Rust => Some(format!(
                    "const {name}: f64 = {};",
                    self.constant_literal(name, 17).unwrap_or_else(|| mpf_to_str(v, 17, true))
                )),
                Lang::Fortran => Some(format!(
                    "parameter ({name} = {}d0)",
                    self.constant_literal(name, 17).unwrap_or_else(|| mpf_to_str(v, 17, true))
                )),
                Lang::C | Lang::Cxx | Lang::Cuda => Some(format!(
                    "const {} {name} = {};",
                    self.s.real.c_name(),
                    match self.constant_literal(name, self.s.real.digits()) {
                        Some(c) => format!("{c}{}", self.s.real.literal_suffix()),
                        None => self.float_literal(v),
                    }
                )),
                _ => None,
            })
            .collect()
    }

    /// Print a symbol name (reserved words get a `_` suffix).
    pub fn symbol_name(&self, name: &str) -> String {
        self.print_symbol(name)
    }

    fn macro_for(&self, e: &Expr) -> Option<&'static str> {
        if self.macros.is_empty() {
            return None;
        }
        self.macros.iter().find(|(m, _)| m == e).map(|(_, n)| *n)
    }

    /// Precedence of the printed form of `e`.
    fn prec(&self, e: &Expr) -> i32 {
        if let Kind::Pow(b, x) = e.kind() {
            if self.inline_pow_n(b, x).is_some() {
                return PREC_MUL;
            }
        }
        precedence(e)
    }

    fn parenthesize(&self, e: &Expr, level: i32, strict: bool) -> CodeResult {
        let p = self.prec(e);
        let s = self.print(e)?;
        if p < level || (!strict && p <= level) {
            Ok(format!("({s})"))
        } else {
            Ok(s)
        }
    }

    /// Print a whole expression: a lone `Float` gets SymPy's full
    /// precision `str` form in the StrPrinter based languages.
    pub fn print_top(&self, e: &Expr) -> CodeResult {
        let _session = super::order::session();
        if let Kind::Num(Number::Float(f)) = e.kind() {
            if matches!(self.lang(), Lang::Python | Lang::NumPy | Lang::Rust) {
                return Ok(crate::number::format_float_full(*f));
            }
        }
        self.print(e)
    }

    pub fn print(&self, e: &Expr) -> CodeResult {
        match e.kind() {
            Kind::Num(n) => Ok(self.print_number(n)),
            Kind::Sym(name) => Ok(self.print_symbol(name)),
            Kind::Const(c) => self.print_constant(e, *c),
            Kind::Bool(b) => Ok(match self.lang() {
                Lang::Python | Lang::NumPy => if *b { "True" } else { "False" }.into(),
                Lang::Fortran => if *b { ".true." } else { ".false." }.into(),
                _ => if *b { "true" } else { "false" }.into(),
            }),
            Kind::Add(..) => self.print_add(e),
            Kind::Mul(..) => {
                if let Some(m) = self.macro_for(e) {
                    return Ok(m.to_string());
                }
                self.print_mul(e)
            }
            Kind::Pow(b, x) => {
                if let Some(m) = self.macro_for(e) {
                    return Ok(m.to_string());
                }
                self.print_pow(b, x)
            }
            Kind::Func(..) => {
                if let Some(m) = self.macro_for(e) {
                    return Ok(m.to_string());
                }
                self.print_func(e)
            }
            Kind::Rel(op, a, b) => self.print_rel(*op, a, b),
            Kind::Piecewise(p) => self.print_piecewise(p),
            Kind::Derivative(..) => Err(CodeError::Unsupported("Derivative".into())),
        }
    }

    // ------------------------------------------------------------------
    // atoms
    // ------------------------------------------------------------------

    fn print_symbol(&self, name: &str) -> String {
        if is_reserved(self.lang(), name) {
            format!("{name}_")
        } else {
            name.to_string()
        }
    }

    fn float_literal(&self, f: f64) -> String {
        let num = mpf_to_str(f, self.s.real.digits(), false);
        let num = num.strip_prefix('+').unwrap_or(&num).to_string();
        let mut parts: Vec<String> = num.split('e').map(|s| s.to_string()).collect();
        if parts.len() == 1 && !parts[0].contains('.') {
            parts[0].push_str(".0");
        }
        let m = parts[0].trim_end_matches('0').to_string();
        parts[0] = if m.ends_with('.') { format!("{m}0") } else { m };
        format!("{}{}", parts.join("e"), self.s.real.literal_suffix())
    }

    fn print_number(&self, n: &Number) -> String {
        match n {
            Number::Int(_) | Number::BigInt(_) => {
                if self.lang() == Lang::Fortran {
                    format!("{n}")
                } else {
                    n.to_string()
                }
            }
            Number::Rational(r) => {
                let (p, q) = (r.numer(), r.denom());
                match self.lang() {
                    Lang::C | Lang::Cxx | Lang::Cuda => {
                        let sfx = self.s.real.literal_suffix();
                        format!("{p}.0{sfx}/{q}.0{sfx}")
                    }
                    Lang::Python | Lang::NumPy => format!("{p}/{q}"),
                    Lang::Rust => format!("{p}_f64/{q}.0"),
                    Lang::Fortran => format!("{p}.0d0/{q}.0d0"),
                }
            }
            Number::Float(f) => match self.lang() {
                Lang::C | Lang::Cxx | Lang::Cuda => self.float_literal(*f),
                Lang::Fortran => {
                    if self.n17.borrow().contains(&f.to_bits()) {
                        return fortran_n17(*f);
                    }
                    let s = float_str(*f);
                    match s.split_once('e') {
                        Some((m, x)) => format!("{m}d{x}"),
                        None => format!("{s}d0"),
                    }
                }
                _ => float_str(*f),
            },
            Number::MpFloat(..) => {
                let s = n.to_string();
                match self.lang() {
                    Lang::C | Lang::Cxx | Lang::Cuda => format!("{s}{}", self.s.real.literal_suffix()),
                    Lang::Rust => format!("{s}_f64"),
                    Lang::Fortran => match s.split_once('e') {
                        Some((mantissa, exponent)) => format!("{mantissa}d{exponent}"),
                        None => format!("{s}d0"),
                    },
                    _ => s,
                }
            }
        }
    }

    /// Exact decimal expansion of a named constant (SymPy evaluates the
    /// constant to the target precision, not the rounded double).
    fn constant_literal(&self, name: &str, digits: usize) -> Option<String> {
        let exact = match name {
            "pi" => "3.14159265358979323846264338327950288",
            "E" => "2.71828182845904523536028747135266250",
            "EulerGamma" => "0.577215664901532860606512090082402431",
            _ => return None,
        };
        Some(round_decimal(exact, digits))
    }

    fn number_symbol(&self, name: &'static str, v: f64) -> String {
        match self.lang() {
            Lang::Python => self.constant_literal(name, 17).unwrap_or_else(|| mpf_to_str(v, 17, true)),
            Lang::Cuda => match self.constant_literal(name, self.s.real.digits()) {
                Some(s) => format!("{s}{}", self.s.real.literal_suffix()),
                None => self.float_literal(v),
            },
            _ => {
                self.number_symbols.borrow_mut().push((name, v));
                name.to_string()
            }
        }
    }

    fn print_constant(&self, e: &Expr, c: Constant) -> CodeResult {
        if let Some(m) = self.macro_for(e) {
            return Ok(m.to_string());
        }
        let l = self.lang();
        let unsupported = || Err(CodeError::Unsupported(c.name().to_string()));
        Ok(match c {
            Constant::Pi => match l {
                Lang::Python => "math.pi".into(),
                Lang::NumPy => "numpy.pi".into(),
                Lang::Rust => "PI".into(),
                Lang::Fortran => self.number_symbol("pi", std::f64::consts::PI),
                Lang::Cuda => self.number_symbol("pi", std::f64::consts::PI),
                _ => self.number_symbol("pi", std::f64::consts::PI),
            },
            Constant::E => match l {
                Lang::Python => "math.e".into(),
                Lang::NumPy => "numpy.e".into(),
                Lang::Rust => "E".into(),
                Lang::Fortran => self.number_symbol("E", std::f64::consts::E),
                Lang::Cuda => self.number_symbol("E", std::f64::consts::E),
                _ => self.number_symbol("E", std::f64::consts::E),
            },
            Constant::EulerGamma => match l {
                Lang::NumPy => "numpy.euler_gamma".into(),
                _ => self.number_symbol("EulerGamma", 0.577_215_664_901_532_9),
            },
            Constant::I => match l {
                Lang::Python | Lang::NumPy => "1j".into(),
                Lang::Fortran => "cmplx(0,1)".into(),
                _ => "I".into(),
            },
            Constant::Infinity => match l {
                Lang::Python => "math.inf".into(),
                Lang::NumPy => "numpy.inf".into(),
                Lang::Fortran => "(huge(0d0) + 1)".into(),
                _ => "INFINITY".into(),
            },
            Constant::NegInfinity => match l {
                Lang::Python => "float('-inf')".into(),
                Lang::NumPy => "-numpy.inf".into(),
                Lang::Rust => "NEG_INFINITY".into(),
                Lang::Fortran => "-(huge(0d0) + 1)".into(),
                _ => "-INFINITY".into(),
            },
            Constant::NaN => match l {
                Lang::Python => "math.nan".into(),
                Lang::NumPy => "numpy.nan".into(),
                Lang::Fortran => return unsupported(),
                _ => "NAN".into(),
            },
            Constant::ComplexInfinity => match l {
                Lang::Python => "math.nan".into(),
                Lang::NumPy => "numpy.nan".into(),
                _ => return unsupported(),
            },
        })
    }

    // ------------------------------------------------------------------
    // Add / Mul
    // ------------------------------------------------------------------

    /// Rust: integers are written as floats next to float-like values.
    fn rust_float_cast(&self, e: &Expr) -> String {
        match e.kind() {
            Kind::Num(n) if n.is_integer() => float_str(n.to_f64()),
            _ => String::new(),
        }
    }

    fn print_add(&self, e: &Expr) -> CodeResult {
        let terms = crate::printer::ordered_terms(e);
        let rust_floats = self.lang() == Lang::Rust && terms.iter().any(real_non_integer_number);
        let mut out = String::new();
        for (i, term) in terms.iter().enumerate() {
            let mut t = if rust_floats && term.is_integer_number() {
                self.rust_float_cast(term)
            } else {
                self.print(term)?
            };
            let sign = if t.starts_with('-') {
                t.remove(0);
                "-"
            } else {
                "+"
            };
            if self.prec(term) < PREC_ADD {
                t = format!("({t})");
            }
            if i == 0 {
                if sign == "-" {
                    out.push('-');
                }
            } else {
                out.push(' ');
                out.push_str(sign);
                out.push(' ');
            }
            out.push_str(&t);
        }
        Ok(out)
    }

    fn print_mul(&self, e: &Expr) -> CodeResult {
        let Kind::Mul(c, facs) = e.kind() else { unreachable!() };
        let prec = precedence(e);
        let (sign, pos) = if c.is_negative() {
            ("-", Expr::from_kind(Kind::Mul(c.neg(), facs.clone())))
        } else {
            ("", e.clone())
        };
        let args = ordered_factors(&pos);

        // Rust: float-like factors turn every other factor into an opaque
        // `TypeCast` (no numerator/denominator split).
        if self.lang() == Lang::Rust
            && (!c.is_integer() || facs.iter().any(|(b, x)| real_non_integer_number(&mul_factor(b, x))))
        {
            let mut a_str = Vec::with_capacity(args.len());
            let single = args.len() == 1 && sign == "-";
            for a in &args {
                let level = if single { (PREC_POW + PREC_MUL) / 2 } else { prec };
                let s = if a.is_integer_number() { self.rust_float_cast(a) } else { self.print(a)? };
                let p = if a.is_integer_number() { PREC_ATOM } else { self.prec(a) };
                a_str.push(if p <= level { format!("({s})") } else { s });
            }
            return Ok(format!("{sign}{}", a_str.join("*")));
        }

        let mut a: Vec<Expr> = Vec::new();
        let mut b: Vec<Den> = Vec::new();
        for item in args {
            if let Kind::Pow(base, x) = item.kind() {
                if neg_rational(x) {
                    if x.is_minus_one() {
                        b.push(Den::Base(base.clone()));
                    } else {
                        b.push(Den::Pow(base.clone(), neg(x)));
                    }
                    continue;
                }
            }
            a.push(item);
        }
        let a_str: Vec<String> = if a.is_empty() {
            vec![self.print_number(&Number::ONE)]
        } else if a.len() == 1 && sign == "-" {
            vec![self.parenthesize(&a[0], (PREC_POW + PREC_MUL) / 2, false)?]
        } else {
            a.iter()
                .map(|x| match x.kind() {
                    // inline products need no parentheses inside a product
                    Kind::Pow(pb, px) if self.inline_pow_n(pb, px).is_some() => self.print(x),
                    _ => self.parenthesize(x, prec, false),
                })
                .collect::<Result<_, _>>()?
        };
        let mut b_str = Vec::with_capacity(b.len());
        for d in &b {
            b_str.push(match d {
                Den::Base(base) => self.parenthesize(base, prec, false)?,
                Den::Pow(base, x) => {
                    let s = self.print_pow(base, x)?;
                    let p = if self.inline_pow_n(base, x).is_some() { PREC_MUL } else { PREC_POW };
                    if p <= prec {
                        format!("({s})")
                    } else {
                        s
                    }
                }
            });
        }
        let num = a_str.join("*");
        if self.s.fast_math && self.lang() == Lang::Cuda && !b_str.is_empty() {
            // fast single precision division intrinsic
            let den = if b_str.len() == 1 {
                b_str[0].trim_start_matches('(').trim_end_matches(')').to_string()
            } else {
                b_str.join("*")
            };
            let den = if b_str.len() == 1 && !balanced(&den) { b_str[0].clone() } else { den };
            return Ok(format!("{sign}__fdividef({num}, {den})"));
        }
        Ok(match b_str.len() {
            0 => format!("{sign}{num}"),
            1 => format!("{sign}{num}/{}", b_str[0]),
            _ => format!("{sign}{num}/({})", b_str.join("*")),
        })
    }

    // ------------------------------------------------------------------
    // Pow
    // ------------------------------------------------------------------

    /// CUDA: `Some(n)` when `b**x` is printed as an inline product.
    fn inline_pow_n(&self, b: &Expr, x: &Expr) -> Option<i64> {
        if self.lang() != Lang::Cuda || self.s.inline_pow < 2 || !b.is_symbol() {
            return None;
        }
        let n = x.as_i64()?;
        if (2..=self.s.inline_pow as i64).contains(&n) {
            Some(n)
        } else {
            None
        }
    }

    fn c_func(&self, name: &str) -> String {
        match self.lang() {
            Lang::Cxx => format!("std::{name}"),
            _ => format!("{name}{}", self.s.real.func_suffix()),
        }
    }

    fn print_pow(&self, b: &Expr, x: &Expr) -> CodeResult {
        match self.lang() {
            Lang::C | Lang::Cxx | Lang::Cuda => {
                if let Some(n) = self.inline_pow_n(b, x) {
                    let s = self.print(b)?;
                    return Ok(vec![s; n as usize].join("*"));
                }
                if is_minus_one_valued(x) {
                    if self.s.fast_math && self.lang() == Lang::Cuda {
                        return Ok(format!("__fdividef({}, {})", self.float_literal(1.0), self.print(b)?));
                    }
                    return Ok(format!(
                        "{}/{}",
                        self.float_literal(1.0),
                        self.parenthesize(b, PREC_POW, false)?
                    ));
                }
                if is_half(x) {
                    return Ok(format!("{}({})", self.c_func("sqrt"), self.print(b)?));
                }
                if self.lang() == Lang::Cuda && self.s.use_rsqrt && is_minus_half(x) {
                    return Ok(format!("{}({})", self.c_func("rsqrt"), self.print(b)?));
                }
                if is_third(x) {
                    return Ok(format!("{}({})", self.c_func("cbrt"), self.print(b)?));
                }
                if self.lang() == Lang::Cuda {
                    if let Some(n) = x.as_i64() {
                        if n < 0 && b.is_symbol() && (2..=self.s.inline_pow as i64).contains(&-n) {
                            let s = self.print(b)?;
                            return Ok(format!(
                                "{}/({})",
                                self.float_literal(1.0),
                                vec![s; (-n) as usize].join("*")
                            ));
                        }
                    }
                }
                // __powf is exp2(y*log2(x)): wrong for negative bases, so it
                // is only used for non-integer exponents
                let pw = if self.s.fast_math && !x.is_integer_number() {
                    "__powf".to_string()
                } else {
                    self.c_func("pow")
                };
                Ok(format!("{pw}({}, {})", self.print(b)?, self.print(x)?))
            }
            Lang::Python | Lang::NumPy => {
                let sqrt = if self.lang() == Lang::NumPy { "numpy.sqrt" } else { "math.sqrt" };
                if *x == Expr::rational(1, 2) {
                    return Ok(format!("{sqrt}({})", self.print(b)?));
                }
                if *x == Expr::rational(-1, 2) {
                    return Ok(format!("1/{sqrt}({})", self.print(b)?));
                }
                if self.lang() == Lang::NumPy {
                    if let Some(n) = x.as_number() {
                        if n.is_integer() && n.is_negative() {
                            let xs = format!("({})", float_str(n.to_f64()));
                            return Ok(format!("{}**{xs}", self.parenthesize(b, PREC_POW, false)?));
                        }
                    }
                } else if x.is_minus_one() {
                    return Ok(format!("1/{}", self.parenthesize(b, PREC_POW, false)?));
                }
                Ok(format!(
                    "{}**{}",
                    self.parenthesize(b, PREC_POW, false)?,
                    self.parenthesize(x, PREC_POW, false)?
                ))
            }
            Lang::Rust => self.rust_pow(b, x),
            Lang::Fortran => {
                if is_minus_one_valued(x) {
                    return Ok(format!("1d0/{}", self.parenthesize(b, PREC_POW, false)?));
                }
                if is_half(x) {
                    if let Some(n) = b.as_number() {
                        if n.is_integer() {
                            return Ok(format!("sqrt({n}.0d0)"));
                        }
                    }
                    return Ok(format!("sqrt({})", self.print(b)?));
                }
                if *x == Expr::rational(-1, 2) {
                    return Ok(format!("1/sqrt({})", self.print(b)?));
                }
                Ok(format!(
                    "{}**{}",
                    self.parenthesize(b, PREC_POW, false)?,
                    self.parenthesize(x, PREC_POW, false)?
                ))
            }
        }
    }

    /// `_print_caller_var` of the Rust printer.
    fn rust_caller(&self, e: &Expr) -> CodeResult {
        let nargs = match e.kind() {
            Kind::Add(c, t) => t.len() + usize::from(!c.is_zero()),
            Kind::Mul(c, f) => f.len() + usize::from(!c.is_one()),
            Kind::Pow(..) => 2,
            Kind::Func(_, a) => a.len(),
            Kind::Rel(..) => 2,
            Kind::Piecewise(p) => p.len(),
            _ => 0,
        };
        if nargs > 1 {
            return Ok(format!("({})", self.print(e)?));
        }
        match e.kind() {
            Kind::Num(Number::Int(i)) => Ok(format!("{i}_i32")),
            Kind::Num(Number::Float(f)) => Ok(format!("{}_f64", float_str(*f))),
            Kind::Num(Number::Rational(_)) => Ok(format!("({})", self.print(e)?)),
            _ => self.print(e),
        }
    }

    fn rust_pow(&self, b: &Expr, x: &Expr) -> CodeResult {
        let mut base_s: Option<String> = None;
        if b.is_integer_number() && !x.is_integer_number() {
            if let Some(xn) = x.as_number() {
                let v = b.as_number().unwrap().to_f64().powf(xn.to_f64());
                return Ok(float_str(v));
            }
            if b.as_i64() != Some(2) {
                base_s = Some(format!("{}_f64", float_str(b.as_number().unwrap().to_f64())));
            }
        }
        let caller = |e: &Expr| -> CodeResult {
            match &base_s {
                Some(s) if e.ptr_eq(b) => Ok(s.clone()),
                _ => self.rust_caller(e),
            }
        };
        if is_minus_one_valued(x) {
            return Ok(format!("{}.recip()", caller(b)?));
        }
        if is_half(x) {
            return Ok(format!("{}.sqrt()", caller(b)?));
        }
        if is_minus_half(x) {
            return Ok(format!("{}.sqrt().recip()", caller(b)?));
        }
        if is_third(x) {
            return Ok(format!("{}.cbrt()", caller(b)?));
        }
        if matches!(b.as_number(), Some(n) if n.to_f64() == 2.0) {
            return Ok(format!("{}.exp2()", self.rust_caller(x)?));
        }
        if x.is_integer_number() {
            return Ok(format!("{}.powi({})", caller(b)?, self.print(x)?));
        }
        Ok(format!("{}.powf({})", caller(b)?, self.print(x)?))
    }

    // ------------------------------------------------------------------
    // functions
    // ------------------------------------------------------------------

    fn sorted_minmax_args(args: &[Expr]) -> Vec<Expr> {
        let mut v = args.to_vec();
        v.sort_by_cached_key(|a| (!a.is_number(), crate::printer::to_str(a)));
        v
    }

    fn join_args(&self, args: &[Expr]) -> CodeResult {
        let v: Vec<String> = args.iter().map(|a| self.print(a)).collect::<Result<_, _>>()?;
        Ok(v.join(", "))
    }

    /// Rewrite of a function in terms of supported ones (SymPy's
    /// `_rewriteable_functions`), printed between parentheses.
    fn rewrite(&self, f: &Func, args: &[Expr]) -> Option<Expr> {
        let x = args.first()?;
        use crate::functions as fx;
        Some(match f {
            Func::Cot => div(&Expr::one(), &fx::tan(x)),
            Func::Sec => div(&Expr::one(), &fx::cos(x)),
            Func::Csc => div(&Expr::one(), &fx::sin(x)),
            Func::Acot => fx::atan(&div(&Expr::one(), x)),
            Func::Coth => {
                let (a, b) = (fx::exp(x), fx::exp(&neg(x)));
                div(&crate::ops::add2(&a, &b), &crate::ops::sub(&a, &b))
            }
            Func::Heaviside => Expr::from_kind(Kind::Piecewise(vec![
                (Expr::zero(), Expr::from_kind(Kind::Rel(RelOp::Lt, x.clone(), Expr::zero()))),
                (Expr::rational(1, 2), Expr::from_kind(Kind::Rel(RelOp::Eq, x.clone(), Expr::zero()))),
                (Expr::one(), Expr::boolean(true)),
            ])),
            _ => return None,
        })
    }

    fn unsupported(&self, f: &Func, args: &[Expr]) -> CodeResult {
        if let Some(n) = self.s.user_functions.get(f.name()) {
            return Ok(format!("{n}({})", self.join_args(args)?));
        }
        if self.s.allow_unknown_functions {
            return Ok(format!("{}({})", f.name(), self.join_args(args)?));
        }
        Err(CodeError::Unsupported(f.name().to_string()))
    }

    fn print_func(&self, e: &Expr) -> CodeResult {
        let Kind::Func(f, args) = e.kind() else { unreachable!() };
        if matches!(f, Func::MatrixTrace(_) | Func::MatrixDeterminant(_)) {
            let explicit = crate::functions::materialize_matrix_scalars(e)
                .map_err(|error| CodeError::Unsupported(error.to_string()))?;
            return self.print(&explicit);
        }
        if let Func::MatrixElement(parent) = f {
            let crate::matrixexpr::MatrixKind::Symbol(name) = parent.kind() else {
                return Err(CodeError::Unsupported("an unevaluated matrix operation in a scalar entry".into()));
            };
            let name = self.symbol_name(name);
            return match self.lang() {
                Lang::Python | Lang::NumPy => Ok(format!("{name}[{}, {}]", self.print(&args[0])?, self.print(&args[1])?)),
                Lang::Fortran => Ok(format!("{name}({}, {})", self.print(&crate::ops::add2(&args[0], &Expr::one()))?, self.print(&crate::ops::add2(&args[1], &Expr::one()))?)),
                _ => {
                    let index = crate::ops::add2(&crate::ops::mul2(&args[0], parent.cols()), &args[1]);
                    Ok(format!("{name}[{}]", self.print(&index)?))
                }
            };
        }
        if let Func::MatrixProxy(_) = f {
            return Err(CodeError::Unsupported("a matrix expression used as a scalar".into()));
        }
        if let Some(n) = self.s.user_functions.get(f.name()) {
            return Ok(format!("{n}({})", self.join_args(args)?));
        }
        match self.lang() {
            Lang::C | Lang::Cxx | Lang::Cuda => self.c_function(f, args),
            Lang::Python | Lang::NumPy => self.py_function(f, args),
            Lang::Rust => self.rust_function(f, args),
            Lang::Fortran => self.fortran_function(f, args),
        }
    }

    fn c_function(&self, f: &Func, args: &[Expr]) -> CodeResult {
        let name = match f {
            Func::Abs => "fabs",
            Func::Sin => "sin",
            Func::Cos => "cos",
            Func::Tan => "tan",
            Func::Asin => "asin",
            Func::Acos => "acos",
            Func::Atan => "atan",
            Func::Atan2 => "atan2",
            Func::Exp => "exp",
            Func::Log => "log",
            Func::Sinh => "sinh",
            Func::Cosh => "cosh",
            Func::Tanh => "tanh",
            Func::Floor => "floor",
            Func::Ceiling => "ceil",
            Func::LogGamma => "lgamma",
            Func::Erfc => "erfc",
            Func::Erf => "erf",
            Func::Asinh => "asinh",
            Func::Acosh => "acosh",
            Func::Atanh => "atanh",
            Func::Gamma => "tgamma",
            Func::Max | Func::Min => {
                let v = Self::sorted_minmax_args(args);
                let name = match (self.lang(), f) {
                    (Lang::Cxx, Func::Max) => "std::max".to_string(),
                    (Lang::Cxx, _) => "std::min".to_string(),
                    (_, Func::Max) => self.c_func("fmax"),
                    _ => self.c_func("fmin"),
                };
                let mut s = self.print(&v[v.len() - 1])?;
                for a in v[..v.len() - 1].iter().rev() {
                    s = format!("{name}({}, {s})", self.print(a)?);
                }
                return Ok(s);
            }
            Func::Sign => {
                let a = self.print(&args[0])?;
                return Ok(format!("((({a}) > 0) - (({a}) < 0))"));
            }
            _ => {
                if let Some(r) = self.rewrite(f, args) {
                    return Ok(format!("({})", self.print(&r)?));
                }
                return self.unsupported(f, args);
            }
        };
        let fname = if self.s.fast_math {
            match name {
                "sin" | "cos" | "tan" | "exp" | "log" => format!("__{name}f"),
                _ => self.c_func(name),
            }
        } else {
            self.c_func(name)
        };
        Ok(format!("{fname}({})", self.join_args(args)?))
    }

    fn py_function(&self, f: &Func, args: &[Expr]) -> CodeResult {
        let np = self.lang() == Lang::NumPy;
        let math_name = match f {
            Func::Acos => Some("acos"),
            Func::Acosh => Some("acosh"),
            Func::Asin => Some("asin"),
            Func::Asinh => Some("asinh"),
            Func::Atan => Some("atan"),
            Func::Atan2 => Some("atan2"),
            Func::Atanh => Some("atanh"),
            Func::Ceiling => Some("ceil"),
            Func::Cos => Some("cos"),
            Func::Cosh => Some("cosh"),
            Func::Erf => Some("erf"),
            Func::Erfc => Some("erfc"),
            Func::Exp => Some("exp"),
            Func::Floor => Some("floor"),
            Func::Gamma => Some("gamma"),
            Func::LogGamma => Some("lgamma"),
            Func::Log => Some("log"),
            Func::Sin => Some("sin"),
            Func::Sinh => Some("sinh"),
            Func::Tan => Some("tan"),
            Func::Tanh => Some("tanh"),
            _ => None,
        };
        if let Some(m) = math_name {
            let full = if np && !matches!(f, Func::Erf | Func::Erfc | Func::Gamma | Func::LogGamma) {
                let n = match m {
                    "acos" => "arccos",
                    "acosh" => "arccosh",
                    "asin" => "arcsin",
                    "asinh" => "arcsinh",
                    "atan" => "arctan",
                    "atan2" => "arctan2",
                    "atanh" => "arctanh",
                    o => o,
                };
                format!("numpy.{n}")
            } else {
                format!("math.{m}")
            };
            return Ok(format!("{full}({})", self.join_args(args)?));
        }
        match f {
            Func::Abs => Ok(format!("abs({})", self.join_args(args)?)),
            Func::Max | Func::Min => {
                let v = Self::sorted_minmax_args(args);
                if np {
                    let red = if *f == Func::Max { "numpy.maximum" } else { "numpy.minimum" };
                    Ok(format!("functools.reduce({red}, [{}])", self.join_args(&v)?))
                } else {
                    let n = if *f == Func::Max { "max" } else { "min" };
                    Ok(format!("{n}({})", self.join_args(&v)?))
                }
            }
            Func::Sign => {
                let a = self.print(&args[0])?;
                if np {
                    Ok(format!("numpy.sign({a})"))
                } else {
                    Ok(format!("(0.0 if {a} == 0 else math.copysign(1, {a}))"))
                }
            }
            Func::Re if np => Ok(format!("numpy.real({})", self.join_args(args)?)),
            Func::Im if np => Ok(format!("numpy.imag({})", self.join_args(args)?)),
            Func::Arg if np => Ok(format!("numpy.angle({})", self.join_args(args)?)),
            Func::Conjugate if np => Ok(format!("numpy.conjugate({})", self.join_args(args)?)),
            _ => {
                if let Some(r) = self.rewrite(f, args) {
                    return Ok(format!("({})", self.print(&r)?));
                }
                self.unsupported(f, args)
            }
        }
    }

    fn rust_function(&self, f: &Func, args: &[Expr]) -> CodeResult {
        let method = match f {
            Func::Abs => "abs",
            Func::Floor => "floor",
            Func::Ceiling => "ceil",
            Func::Exp => "exp",
            Func::Log => "ln",
            Func::Sin => "sin",
            Func::Cos => "cos",
            Func::Tan => "tan",
            Func::Asin => "asin",
            Func::Acos => "acos",
            Func::Atan => "atan",
            Func::Sinh => "sinh",
            Func::Cosh => "cosh",
            Func::Tanh => "tanh",
            Func::Asinh => "asinh",
            Func::Acosh => "acosh",
            Func::Atanh => "atanh",
            Func::Atan2 => {
                return Ok(format!("{}.atan2({})", self.rust_caller(&args[0])?, self.print(&args[1])?));
            }
            Func::Max | Func::Min => {
                // deviation: SymPy rewrites to a Piecewise; method chains are
                // valid Rust and much shorter
                let v = Self::sorted_minmax_args(args);
                let m = if *f == Func::Max { "max" } else { "min" };
                let mut s = self.rust_caller(&v[0])?;
                for a in &v[1..] {
                    s = format!("{s}.{m}({})", self.print(a)?);
                }
                return Ok(s);
            }
            Func::Sign => {
                let a = self.print(&args[0])?;
                return Ok(format!("(if ({a} == 0.0) {{ 0.0 }} else {{ ({a}).signum() }})"));
            }
            _ => {
                if let Some(r) = self.rewrite(f, args) {
                    return self.print(&r);
                }
                return self.unsupported(f, args);
            }
        };
        Ok(format!("{}.{method}()", self.rust_caller(&args[0])?))
    }

    fn fortran_function(&self, f: &Func, args: &[Expr]) -> CodeResult {
        // the arguments are evaluated with N(): numbers become floats and a
        // function of numbers is evaluated
        let v: Vec<Expr> = args.iter().map(|a| self.fortran_n(a)).collect();
        if !matches!(f, Func::Max | Func::Min | Func::Undefined(_))
            && args.iter().all(|a| super::is_numeric(a) && !crate::functions::contains_i(a))
        {
            let call = crate::functions::func(f.clone(), args.to_vec());
            if let Some(x) = crate::evalf::eval_f64(&call) {
                return Ok(fortran_n17(x));
            }
        }
        let name = match f {
            Func::Sin => "sin",
            Func::Cos => "cos",
            Func::Tan => "tan",
            Func::Asin => "asin",
            Func::Acos => "acos",
            Func::Atan => "atan",
            Func::Atan2 => "atan2",
            Func::Sinh => "sinh",
            Func::Cosh => "cosh",
            Func::Tanh => "tanh",
            Func::Log => "log",
            Func::Exp => "exp",
            Func::Erf => "erf",
            Func::Abs => "abs",
            Func::Conjugate => "conjg",
            Func::Max | Func::Min => {
                // not Functions in SymPy: no float conversion of the args
                let v = Self::sorted_minmax_args(args);
                let n = if *f == Func::Max { "max" } else { "min" };
                return Ok(format!("{n}({})", self.join_args(&v)?));
            }
            Func::Erfc => {
                return Ok(format!("(1 - {})", self.fortran_function(&Func::Erf, args)?));
            }
            Func::Heaviside if self.s.fortran_standard >= 95 => {
                let x = &args[0];
                let pw = vec![
                    (Expr::zero(), Expr::from_kind(Kind::Rel(RelOp::Lt, x.clone(), Expr::zero()))),
                    (Expr::float(0.5), Expr::from_kind(Kind::Rel(RelOp::Eq, x.clone(), Expr::zero()))),
                    (Expr::one(), Expr::boolean(true)),
                ];
                return Ok(format!("({})", self.print_piecewise(&pw)?));
            }
            Func::Sign => {
                if self.s.fortran_standard < 95 {
                    return Err(CodeError::Unsupported("sign (merge requires Fortran 95 or newer)".into()));
                }
                let a = self.print(&args[0])?;
                return Ok(format!("merge(0d0, dsign(1d0, {a}), {a} == 0d0)"));
            }
            _ => {
                if let Some(r) = self.rewrite(f, args) {
                    if !matches!(f, Func::Heaviside) {
                        return Ok(format!("({})", self.print(&r)?));
                    }
                }
                return self.unsupported(f, args);
            }
        };
        let parts: Vec<String> = v.iter().map(|a| self.print(a)).collect::<Result<_, _>>()?;
        Ok(format!("{name}({})", parts.join(", ")))
    }

    /// A float evaluated with `N(..., 17)` (printed with 17 digits).
    fn n17_float(&self, x: f64) -> Expr {
        let e = Expr::float(x);
        if let Kind::Num(Number::Float(f)) = e.kind() {
            self.n17.borrow_mut().insert(f.to_bits());
        }
        e
    }

    /// SymPy's `N(e)` for the Fortran printer: integers become floats,
    /// numeric sub-expressions are evaluated to 17 digit floats.
    fn fortran_n(&self, e: &Expr) -> Expr {
        if e.is_atom() && !e.is_number() {
            return e.clone();
        }
        if let Some(r) = self.n_memo.borrow().get(e) {
            return r.clone();
        }
        let r = self.fortran_n_impl(e);
        self.n_memo.borrow_mut().insert(e.clone(), r.clone());
        r
    }

    fn fortran_n_impl(&self, e: &Expr) -> Expr {
        match e.kind() {
            Kind::Num(Number::Int(_) | Number::BigInt(_)) => Expr::float(e.as_number().unwrap().to_f64()),
            Kind::Num(Number::Float(_)) => e.clone(),
            Kind::Num(n) => self.n17_float(n.to_f64()),
            Kind::Sym(_) | Kind::Bool(_) => e.clone(),
            _ => {
                if super::is_numeric(e) && !crate::functions::contains_i(e) {
                    if let Some(x) = crate::evalf::eval_f64(e) {
                        return self.n17_float(x);
                    }
                }
                match e.kind() {
                    Kind::Mul(c, facs) => {
                        let mut v: Vec<Expr> = Vec::with_capacity(facs.len() + 1);
                        if c.is_minus_one() {
                            v.push(Expr::minus_one());
                        } else if !c.is_one() {
                            v.push(self.fortran_n(&Expr::num(c.clone())));
                        }
                        v.extend(facs.iter().map(|(b, x)| self.fortran_n(&mul_factor(b, x))));
                        crate::ops::mul(&v)
                    }
                    Kind::Add(c, terms) => {
                        let mut v: Vec<Expr> =
                            terms.iter().map(|(t, k)| self.fortran_n(&crate::ops::add_term(t, k))).collect();
                        if !c.is_zero() {
                            v.push(self.fortran_n(&Expr::num(c.clone())));
                        }
                        crate::ops::add(&v)
                    }
                    Kind::Pow(b, x) => {
                        let nx = if x.is_integer_number() { x.clone() } else { self.fortran_n(x) };
                        crate::ops::pow(&self.fortran_n(b), &nx)
                    }
                    _ if e.is_atom() => e.clone(),
                    _ => crate::traverse::map_children(e, &mut |c| self.fortran_n(c)),
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // relationals and piecewise
    // ------------------------------------------------------------------

    fn print_rel(&self, op: RelOp, a: &Expr, b: &Expr) -> CodeResult {
        match self.lang() {
            Lang::Python => Ok(format!("({} {} {})", self.print(a)?, op.symbol(), self.print(b)?)),
            Lang::NumPy => {
                let f = match op {
                    RelOp::Eq => "equal",
                    RelOp::Ne => "not_equal",
                    RelOp::Lt => "less",
                    RelOp::Le => "less_equal",
                    RelOp::Gt => "greater",
                    RelOp::Ge => "greater_equal",
                };
                Ok(format!("numpy.{f}({}, {})", self.print(a)?, self.print(b)?))
            }
            Lang::Rust => {
                let (ai, bi) = (a.is_integer_number(), b.is_integer_number());
                let cast = ai != bi;
                let pa = if cast && ai { self.rust_float_cast(a) } else { self.print(a)? };
                let pb = if cast && bi { self.rust_float_cast(b) } else { self.print(b)? };
                Ok(format!("{pa} {} {pb}", op.symbol()))
            }
            Lang::Fortran => {
                let o = match op {
                    RelOp::Ne => "/=",
                    o => o.symbol(),
                };
                Ok(format!("{} {o} {}", self.print(a)?, self.print(b)?))
            }
            _ => Ok(format!("{} {} {}", self.print(a)?, op.symbol(), self.print(b)?)),
        }
    }

    fn print_piecewise(&self, p: &[(Expr, Expr)]) -> CodeResult {
        let has_default = matches!(p.last().map(|x| x.1.kind()), Some(Kind::Bool(true)));
        match self.lang() {
            Lang::C | Lang::Cxx | Lang::Cuda => {
                if !has_default {
                    return Err(CodeError::PiecewiseNoDefault);
                }
                let mut pairs = Vec::new();
                for (e, c) in &p[..p.len() - 1] {
                    pairs.push(format!("(({}) ? (\n{}\n)\n", self.print(c)?, self.print(e)?));
                }
                let last = format!(": (\n{}\n)", self.print(&p[p.len() - 1].0)?);
                Ok(format!("{}{last}{}", pairs.join(": "), ")".repeat(pairs.len())))
            }
            Lang::Python => {
                let mut out = String::from("(");
                for (i, (e, c)) in p.iter().enumerate() {
                    let last = i + 1 == p.len();
                    out.push_str(&format!("({})", self.print(e)?));
                    let cs = self.print(c)?;
                    if last && cs == "True" {
                        break;
                    }
                    out.push_str(" if ");
                    out.push_str(&cs);
                    out.push_str(if last { " else None" } else { " else " });
                }
                out.push(')');
                Ok(out)
            }
            Lang::NumPy => {
                let conds: Vec<String> = p.iter().map(|(_, c)| self.print(c)).collect::<Result<_, _>>()?;
                let exprs: Vec<String> = p.iter().map(|(e, _)| self.print(e)).collect::<Result<_, _>>()?;
                Ok(format!(
                    "numpy.select([{}], [{}], default=numpy.nan)",
                    conds.join(","),
                    exprs.join(",")
                ))
            }
            Lang::Rust => {
                if !has_default {
                    return Err(CodeError::PiecewiseNoDefault);
                }
                let mut lines: Vec<String> = Vec::new();
                for (i, (e, c)) in p.iter().enumerate() {
                    if i == 0 {
                        lines.push(format!("if ({}) {{", self.print(c)?));
                    } else if i + 1 == p.len() && matches!(c.kind(), Kind::Bool(true)) {
                        lines.last_mut().unwrap().push_str(" else {");
                    } else {
                        let cs = self.print(c)?;
                        lines.last_mut().unwrap().push_str(&format!(" else if ({cs}) {{"));
                    }
                    lines.push(self.print(e)?);
                    lines.push("}".into());
                }
                Ok(lines.join("\n"))
            }
            Lang::Fortran => {
                if self.s.fortran_standard < 95 {
                    return Err(CodeError::Unsupported(
                        "Piecewise as an inline expression requires Fortran 95 or newer".into(),
                    ));
                }
                if !has_default {
                    return Err(CodeError::PiecewiseNoDefault);
                }
                let mut s = self.print(&p[p.len() - 1].0)?;
                for (e, c) in p[..p.len() - 1].iter().rev() {
                    s = format!("merge({}, {s}, {})", self.print(e)?, self.print(c)?);
                }
                Ok(s)
            }
        }
    }
}

/// The `math.h` macros SymPy's C printers use for well known constants.
fn math_macros() -> Vec<(Expr, &'static str)> {
    use crate::functions::log;
    let pi = Expr::pi();
    let two = Expr::int(2);
    vec![
        (Expr::e(), "M_E"),
        (div(&Expr::one(), &log(&two)), "M_LOG2E"),
        (log(&two), "M_LN2"),
        (log(&Expr::int(10)), "M_LN10"),
        (pi.clone(), "M_PI"),
        (div(&pi, &two), "M_PI_2"),
        (div(&pi, &Expr::int(4)), "M_PI_4"),
        (div(&Expr::one(), &pi), "M_1_PI"),
        (div(&two, &pi), "M_2_PI"),
        (div(&two, &crate::ops::sqrt(&pi)), "M_2_SQRTPI"),
        (crate::ops::sqrt(&two), "M_SQRT2"),
        (div(&Expr::one(), &crate::ops::sqrt(&two)), "M_SQRT1_2"),
    ]
}

/// Round a positive decimal string `d.ddd...` to `digits` significant
/// digits (round half up on the exact expansion), trailing zeros stripped.
fn round_decimal(exact: &str, digits: usize) -> String {
    let (ip, fp) = exact.split_once('.').unwrap();
    let all: Vec<u8> = ip.bytes().chain(fp.bytes()).map(|c| c - b'0').collect();
    // position of the first significant digit
    let first = all.iter().position(|d| *d != 0).unwrap_or(0);
    let end = (first + digits).min(all.len());
    let mut kept: Vec<u8> = all[..end].to_vec();
    if end < all.len() && all[end] >= 5 {
        let mut i = kept.len();
        loop {
            i -= 1;
            if kept[i] == 9 {
                kept[i] = 0;
                if i == 0 {
                    kept.insert(0, 1);
                    break;
                }
            } else {
                kept[i] += 1;
                break;
            }
        }
    }
    let int_len = ip.len() + (kept.len() - end);
    let mut s: String = kept[..int_len].iter().map(|d| (d + b'0') as char).collect();
    let mut f: String = kept[int_len..].iter().map(|d| (d + b'0') as char).collect();
    while f.ends_with('0') {
        f.pop();
    }
    if f.is_empty() {
        f.push('0');
    }
    s.push('.');
    s.push_str(&f);
    s
}

/// A float evaluated with `N(..., 17)`, printed with 17 digits.
fn fortran_n17(x: f64) -> String {
    let s = mpf_to_str(x, 17, true);
    match s.split_once('e') {
        Some((m, e)) => format!("{m}d{}", e.strip_prefix('+').unwrap_or(e)),
        None => format!("{s}d0"),
    }
}

/// Are the parentheses of `s` balanced (so that it can stand alone)?
fn balanced(s: &str) -> bool {
    let mut depth = 0i64;
    for c in s.chars() {
        match c {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth < 0 {
                    return false;
                }
            }
            _ => {}
        }
    }
    depth == 0
}

/// SymPy's Fortran line layout: 6 column lead (fixed form) and lines
/// wrapped at column 72 with `@` continuation lines (`&` in free form).
fn fortran_wrap(lines: &[String], free: bool) -> Vec<String> {
    let (lead, cont, trailing) = if free { ("", "      ", " &") } else { ("      ", "     @ ", "") };
    let alnum = |c: char| c == '_' || c == '+' || c == '-' || c == '.' || c.is_ascii_alphanumeric();
    let white = |c: char| c == ' ' || c == '\t' || c == '(' || c == ')';
    let split_pos = |line: &[char], endpos: usize| -> usize {
        if line.len() <= endpos {
            return line.len();
        }
        let split = |p: usize| {
            let (c, d) = (line[p], line[p - 1]);
            (alnum(c) && !alnum(d)) || (!alnum(c) && alnum(d)) || (white(c) && !white(d)) || (!white(c) && white(d))
        };
        let mut pos = endpos;
        while !split(pos) {
            pos -= 1;
            if pos == 0 {
                return endpos;
            }
        }
        pos
    };
    let mut out = Vec::new();
    for raw in lines {
        for l in raw.split('\n') {
            let line: Vec<char> = format!("{lead}{l}").chars().collect();
            let pos = split_pos(&line, 72);
            let mut hunk: String = line[..pos].iter().collect::<String>().trim_end().to_string();
            let mut rest: Vec<char> = line[pos..].iter().copied().skip_while(|c| *c == ' ').collect();
            if !rest.is_empty() {
                hunk.push_str(trailing);
            }
            out.push(hunk);
            while !rest.is_empty() {
                let pos = split_pos(&rest, 65);
                let mut hunk: String = rest[..pos].iter().collect::<String>().trim_end().to_string();
                rest = rest[pos..].iter().copied().skip_while(|c| *c == ' ').collect();
                if !rest.is_empty() {
                    hunk.push_str(trailing);
                }
                out.push(format!("{cont}{hunk}"));
            }
        }
    }
    out
}

/// SymPy's `indent_code` for brace/paren languages.
pub(crate) fn indent_code(code: &str, tab: &str) -> String {
    if !code.contains('\n') {
        return code.to_string();
    }
    let mut out = String::with_capacity(code.len() + 64);
    let mut level: i64 = 0;
    for raw in code.split_inclusive('\n') {
        let line = raw.trim_start_matches([' ', '\t']);
        if line.is_empty() || line == "\n" {
            out.push_str(line);
            continue;
        }
        let body = line.strip_suffix('\n').unwrap_or(line);
        let inc = body.ends_with('{') || body.ends_with('(');
        let dec = line.starts_with('}') || line.starts_with(')');
        if dec {
            level -= 1;
        }
        for _ in 0..level.max(0) {
            out.push_str(tab);
        }
        out.push_str(line);
        if inc {
            level += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn c(s: &str) -> String {
        ccode(&parse_str(s).unwrap()).unwrap()
    }
    fn py(s: &str) -> String {
        pycode(&parse_str(s).unwrap()).unwrap()
    }
    fn rs(s: &str) -> String {
        rust_code(&parse_str(s).unwrap()).unwrap()
    }

    #[test]
    fn constants() {
        assert_eq!(round_decimal("0.577215664901532860606512090082402431", 17), "0.57721566490153286");
        assert_eq!(round_decimal("3.14159265358979323846264338327950288", 9), "3.14159265");
        assert_eq!(round_decimal("2.71828182845904523536028747135266250", 9), "2.71828183");
        assert_eq!(c("EulerGamma"), "const double EulerGamma = 0.57721566490153286;\nEulerGamma");
    }

    #[test]
    fn c_basics() {
        assert_eq!(c("x**2"), "pow(x, 2)");
        assert_eq!(c("1/x"), "1.0/x");
        assert_eq!(c("x/2"), "(1.0/2.0)*x");
        assert_eq!(c("-x/2"), "-1.0/2.0*x");
        assert_eq!(c("sqrt(x)"), "sqrt(x)");
        assert_eq!(c("1/sqrt(x)"), "pow(x, -1.0/2.0)");
        assert_eq!(c("x**(1/3)"), "cbrt(x)");
        assert_eq!(c("pi*x"), "M_PI*x");
        assert_eq!(c("pi/2"), "M_PI_2");
        assert_eq!(c("sec(x)"), "(1.0/cos(x))");
        assert_eq!(c("x + 0.1"), "x + 0.10000000000000001");
        assert_eq!(c("1/(2*x)"), "(1.0/2.0)/x");
        assert_eq!(c("Max(x, y, z)"), "fmax(x, fmax(y, z))");
        assert_eq!(c("Piecewise((x, x > 0), (0, True))"), "((x > 0) ? (\n   x\n)\n: (\n   0\n))");
    }

    #[test]
    fn fortran() {
        let f = |s: &str| fcode(&parse_str(s).unwrap()).unwrap();
        assert_eq!(f("x**2"), "      x**2");
        assert_eq!(f("x/2"), "      (1.0d0/2.0d0)*x");
        assert_eq!(f("1/x"), "      1d0/x");
        assert_eq!(f("pi*x"), "      parameter (pi = 3.1415926535897932d0)\n      pi*x");
        assert_eq!(f("x + 0.1"), "      x + 0.1d0");
        assert_eq!(f("sec(x)"), "      (1d0/cos(x))");
    }

    #[test]
    fn py_and_rust() {
        assert_eq!(py("x**-2"), "x**(-2)");
        assert_eq!(py("x/2"), "(1/2)*x");
        assert_eq!(py("sqrt(x)"), "math.sqrt(x)");
        assert_eq!(py("Piecewise((x, x > 0), (0, True))"), "((x) if (x > 0) else (0))");
        assert_eq!(rs("x**2"), "x.powi(2)");
        assert_eq!(rs("sin(x)**2"), "x.sin().powi(2)");
        assert_eq!(rs("2**x"), "x.exp2()");
        assert_eq!(rs("3**x"), "3.0_f64.powf(x)");
        assert_eq!(rs("x/(2*y)"), "(1_f64/2.0)*x*y.recip()");
        assert_eq!(rs("pi*(x+y)"), "PI*(x + y)");
    }

    #[test]
    fn cuda() {
        let e = parse_str("x**2*sin(y) + 1/sqrt(z)").unwrap();
        let d = cudacode(&e).unwrap();
        assert_eq!(d, "x*x*sin(y) + rsqrt(z)");
        let s = CudaSettings { dtype: FloatType::F32, fast_math: true, ..Default::default() };
        assert_eq!(cudacode_with(&e, &s).unwrap(), "x*x*__sinf(y) + rsqrtf(z)");
        let e = parse_str("y/x**2 + 0.1").unwrap();
        assert_eq!(cudacode_with(&e, &s).unwrap(), "0.1F + __fdividef(y, x*x)");
        let s2 = CudaSettings { dtype: FloatType::F32, ..Default::default() };
        assert_eq!(cudacode_with(&e, &s2).unwrap(), "0.1F + y/(x*x)");
        assert_eq!(cudacode_with(&parse_str("1/(x + y)").unwrap(), &s).unwrap(), "__fdividef(1.0F, x + y)");
    }
}
