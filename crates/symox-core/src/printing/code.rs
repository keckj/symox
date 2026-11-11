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
    /// CUDA: print `x**(-1/2)` as `1.0/sqrt(x)` (two correctly rounded IEEE
    /// operations) instead of `pow(x, -1.0/2.0)`. Opt-in: the value differs
    /// from `pow` at `x = -0.0` (`-inf` instead of `+inf`). Ignored when
    /// `use_rsqrt` applies.
    pub reciprocal_sqrt: bool,
    /// CUDA: integer powers up to this order are written as products
    /// (`x*x*x`, `(x - y)*(x - y)`) for symbols and for small arithmetic bases
    /// (see `cheap_inline_base`). 0 disables.
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
    /// NumPy: SymPy's `SciPyPrinter` names of the `scipy.special` functions.
    pub scipy: bool,
    /// Python/NumPy: the plain forms of SymPy's printers, as lambdify prints
    /// its generated source (no NumPy 2 complex-sign guard).
    pub lambdify: bool,
}

impl CodeSettings {
    pub fn new(lang: Lang) -> Self {
        CodeSettings {
            lang,
            real: FloatType::F64,
            fast_math: false,
            use_rsqrt: false,
            reciprocal_sqrt: false,
            inline_pow: if lang == Lang::Cuda { 4 } else { 0 },
            math_macros: matches!(lang, Lang::C | Lang::Cxx | Lang::Cuda),
            user_functions: FxHashMap::default(),
            allow_unknown_functions: false,
            fortran_standard: 77,
            fortran_free_form: false,
            scipy: false,
            lambdify: false,
        }
    }
}

/// Settings of [`cudacode_with`].
#[derive(Clone, Debug)]
pub struct CudaSettings {
    pub dtype: FloatType,
    pub fast_math: bool,
    pub use_rsqrt: bool,
    /// See [`CodeSettings::reciprocal_sqrt`] (default `false`).
    pub reciprocal_sqrt: bool,
    pub inline_pow: u32,
    /// `M_PI`-style macros; `None` keeps the historical default (double
    /// precision only). NVRTC does not define these macros: use `Some(false)`.
    pub math_macros: Option<bool>,
}

impl Default for CudaSettings {
    fn default() -> Self {
        CudaSettings {
            dtype: FloatType::F64,
            fast_math: false,
            use_rsqrt: true,
            reciprocal_sqrt: false,
            inline_pow: 4,
            math_macros: None,
        }
    }
}

impl CudaSettings {
    pub fn to_code_settings(&self) -> CodeSettings {
        let mut s = CodeSettings::new(Lang::Cuda);
        s.real = self.dtype;
        s.fast_math = self.fast_math && self.dtype == FloatType::F32;
        s.use_rsqrt = self.use_rsqrt;
        s.reciprocal_sqrt = self.reciprocal_sqrt;
        s.inline_pow = self.inline_pow;
        s.math_macros = self.math_macros.unwrap_or(self.dtype == FloatType::F64);
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

/// NumPy 2 complex-sign semantics using APIs available since NumPy 1.21.
/// Bind the input once and dispatch on its runtime dtype. Masked componentwise
/// division avoids both zero warnings and complex subnormal-division overflow.
fn numpy_complex_sign(argument: &str) -> String {
    format!(concat!(
        "(lambda _sgn_value: numpy.sign(_sgn_value) if numpy.__version__.split('.')[0] not in ('0', '1') or not numpy.iscomplexobj(_sgn_value) else ",
        "(lambda _sgn_real, _sgn_imag, _sgn_abs: ",
        "numpy.where(numpy.isnan(_sgn_abs), numpy.nan*(1+1j), ",
        "numpy.where(numpy.isinf(_sgn_abs), ",
        "numpy.where(numpy.isinf(_sgn_real), ",
        "numpy.where(numpy.isinf(_sgn_imag), numpy.nan*(1+1j), numpy.sign(_sgn_real)+0j), ",
        "1j*numpy.sign(_sgn_imag)), ",
        "numpy.divide(_sgn_real, _sgn_abs, out=numpy.zeros_like(_sgn_real), ",
        "where=numpy.isfinite(_sgn_abs)&(_sgn_abs!=0)) + ",
        "1j*numpy.divide(_sgn_imag, _sgn_abs, out=numpy.zeros_like(_sgn_imag), ",
        "where=numpy.isfinite(_sgn_abs)&(_sgn_abs!=0))))",
        ".astype(numpy.asarray(_sgn_value).dtype, copy=False)[()]",
        ")(numpy.real(_sgn_value), numpy.imag(_sgn_value), numpy.abs(_sgn_value)))({})"
    ), argument)
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
        Kind::Const(c) => c.is_number_symbol(),
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
            if self.inline_pow_n(b, x).is_some() || self.reciprocal_sqrt_applies(x) {
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
            Kind::Sym(symbol) => Ok(self.print_symbol(&match symbol.dummy_index() {
                Some(index) => format!("{}_{index}",symbol.name()),
                None => symbol.name().to_owned(),
            })),
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
            Kind::Func(Func::SymbolicDerivative,_) => Err(CodeError::Unsupported("Derivative".into())),
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
            "Catalan" => "0.915965594177219015054603514932384111",
            "GoldenRatio" => "1.61803398874989484820458683436563812",
            "TribonacciConstant" => "1.83928675521416113255185256465328660",
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
                _ => self.number_symbol("EulerGamma", crate::expr::EULER_GAMMA_F64),
            },
            Constant::Catalan | Constant::GoldenRatio | Constant::TribonacciConstant => {
                self.number_symbol(c.name(), c.approx_f64().unwrap())
            }
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
        // A denominator is the right operand of `/`: anything printed at
        // multiplication precedence or below needs parentheses, whatever the
        // sign of the product (`prec` is PREC_ADD for a negative coefficient,
        // which printed `-x/s**2` as `-x/s*s`).
        let den_level = PREC_MUL;
        let mut b_str = Vec::with_capacity(b.len());
        for d in &b {
            b_str.push(match d {
                Den::Base(base) => self.parenthesize(base, den_level, false)?,
                Den::Pow(base, x) => {
                    let s = self.print_pow(base, x)?;
                    let inline = self.inline_pow_n(base, x).is_some();
                    let p = if inline { PREC_MUL } else { PREC_POW };
                    // several denominators are joined with `*` inside one
                    // pair of parentheses: inline products need none there
                    if p <= den_level && !(inline && b.len() > 1) {
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
        if self.lang() != Lang::Cuda || self.s.inline_pow < 2 || !self.cheap_inline_base(b) {
            return None;
        }
        let n = x.as_i64()?;
        if (2..=self.s.inline_pow as i64).contains(&n) {
            Some(n)
        } else {
            None
        }
    }

    /// CUDA: bases repeated textually by inline integer powers: symbols and
    /// small `Add`/`Mul` trees of symbols, numbers and inline powers of
    /// symbols (at most 16 nodes). Function calls are never repeated, so
    /// `sin(x)**2` keeps `pow`. Products are correctly rounded IEEE
    /// operations, evaluated left to right.
    fn cheap_inline_base(&self, b: &Expr) -> bool {
        fn walk(e: &Expr, budget: &mut i32, limit: i64) -> bool {
            *budget -= 1;
            if *budget < 0 {
                return false;
            }
            match e.kind() {
                Kind::Sym(_) => true,
                Kind::Num(n) => !matches!(n, Number::Float(f) if !f.is_finite()),
                Kind::Add(..) | Kind::Mul(..) => {
                    crate::traverse::ordered_args(e).iter().all(|a| walk(a, budget, limit))
                }
                Kind::Pow(pb, px) => {
                    pb.is_symbol() && px.as_i64().is_some_and(|n| (2..=limit).contains(&n))
                }
                _ => false,
            }
        }
        if b.is_symbol() {
            return true;
        }
        let mut budget = 16;
        walk(b, &mut budget, self.s.inline_pow as i64)
    }

    /// CUDA: `x**(-1/2)` printed as `1.0/sqrt(b)`.
    fn reciprocal_sqrt_applies(&self, x: &Expr) -> bool {
        self.lang() == Lang::Cuda && self.s.reciprocal_sqrt && !self.s.use_rsqrt && is_minus_half(x)
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
                    let s = self.parenthesize(b, PREC_MUL, false)?;
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
                if self.reciprocal_sqrt_applies(x) {
                    return Ok(format!(
                        "{}/{}({})",
                        self.float_literal(1.0),
                        self.c_func("sqrt"),
                        self.print(b)?
                    ));
                }
                if is_third(x) {
                    return Ok(format!("{}({})", self.c_func("cbrt"), self.print(b)?));
                }
                if self.lang() == Lang::Cuda {
                    if let Some(n) = x.as_i64() {
                        if n < 0 && self.cheap_inline_base(b) && (2..=self.s.inline_pow as i64).contains(&-n) {
                            let s = self.parenthesize(b, PREC_MUL, false)?;
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
            // SymPy's Heaviside(x, H0).rewrite(Piecewise)
            Func::Heaviside => {
                let h0 = args.get(1).cloned().unwrap_or_else(|| Expr::rational(1, 2));
                let rel = |op, x: &Expr| Expr::from_kind(Kind::Rel(op, x.clone(), Expr::zero()));
                Expr::from_kind(Kind::Piecewise(if h0.is_zero() {
                    vec![(Expr::zero(), rel(RelOp::Le, x)), (Expr::one(), Expr::boolean(true))]
                } else if h0.is_one() {
                    vec![(Expr::zero(), rel(RelOp::Lt, x)), (Expr::one(), Expr::boolean(true))]
                } else {
                    vec![
                        (Expr::zero(), rel(RelOp::Lt, x)),
                        (h0, rel(RelOp::Eq, x)),
                        (Expr::one(), Expr::boolean(true)),
                    ]
                }))
            }
            Func::Special(sp) => return Self::rewrite_elementary(*sp, args),
            _ => return None,
        })
    }

    /// Rewrites of the elementary special heads (SymPy's `_rewriteable_functions`
    /// targets and `Mod` as `p - q*floor(p/q)`).
    fn rewrite_elementary(sp: crate::special::SpecialFn, args: &[Expr]) -> Option<Expr> {
        use crate::functions as fx;
        use crate::ops::{add2, mul2, sub};
        use crate::special::SpecialFn as S;
        let x = &args[0];
        let half = Expr::rational(1, 2);
        let r = || div(&Expr::one(), x);
        Some(match sp {
            S::Sech | S::Csch => {
                let a = mul2(&half, &fx::exp(x));
                let b = mul2(&half, &fx::exp(&neg(x)));
                div(&Expr::one(), &if sp == S::Sech { add2(&a, &b) } else { sub(&a, &b) })
            }
            S::Asec => fx::acos(&r()),
            S::Acsc => fx::asin(&r()),
            S::Asech => fx::log(&add2(
                &mul2(&crate::ops::sqrt(&sub(&r(), &Expr::one())), &crate::ops::sqrt(&add2(&r(), &Expr::one()))),
                &r(),
            )),
            S::Acsch => fx::log(&add2(&crate::ops::sqrt(&add2(&Expr::one(), &crate::ops::pow(x, &Expr::int(-2)))), &r())),
            S::Acoth => mul2(&half, &sub(&fx::log(&add2(&Expr::one(), &r())), &fx::log(&sub(&Expr::one(), &r())))),
            S::Sinc => Expr::from_kind(Kind::Piecewise(vec![
                (div(&fx::sin(x), x), Expr::from_kind(Kind::Rel(RelOp::Ne, x.clone(), Expr::zero()))),
                (Expr::one(), Expr::boolean(true)),
            ])),
            S::Frac => sub(x, &fx::floor(x)),
            S::Mod => sub(x, &mul2(&args[1], &fx::floor(&div(x, &args[1])))),
            S::KroneckerDelta if args.len() == 2 => Expr::from_kind(Kind::Piecewise(vec![
                (Expr::one(), Expr::from_kind(Kind::Rel(RelOp::Eq, x.clone(), args[1].clone()))),
                (Expr::zero(), Expr::boolean(true)),
            ])),
            _ => return None,
        })
    }

    /// `Mod`/`Rem` in C-like languages: `%` for integers, `fmod` otherwise
    /// (like SymPy, `fmod` keeps the sign of the dividend).
    fn c_mod(&self, sp: crate::special::SpecialFn, args: &[Expr]) -> CodeResult {
        use crate::assumptions as a;
        let (p, q) = (&args[0], &args[1]);
        if a::is_integer(p) == Some(true) && a::is_integer(q) == Some(true) {
            let (sp_, sq) = (self.parenthesize(p, PREC_MUL, false)?, self.parenthesize(q, PREC_MUL, false)?);
            let same_sign = (a::is_nonnegative(p) == Some(true) && a::is_nonnegative(q) == Some(true))
                || (a::is_nonpositive(p) == Some(true) && a::is_nonpositive(q) == Some(true));
            if sp == crate::special::SpecialFn::Rem || same_sign {
                return Ok(format!("{sp_} % {sq}"));
            }
            return Ok(format!("(({sp_} % {sq}) + {sq}) % {sq}"));
        }
        Ok(format!("{}({}, {})", self.c_func("fmod"), self.print(p)?, self.print(q)?))
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
        if matches!(f, Func::MatrixTrace(_) | Func::MatrixDeterminant(_) | Func::MatrixDotProduct(_, _)) {
            let explicit = crate::functions::materialize_matrix_scalars(e)
                .map_err(|error| CodeError::Unsupported(error.to_string()))?;
            return self.print(&explicit);
        }
        if let Func::MatrixElement(parent) = f {
            if matches!(parent.kind(), crate::matrixexpr::MatrixKind::Function(_)|crate::matrixexpr::MatrixKind::ApplyFunction(..)) {
                let value = parent.element(&args[0], &args[1]).map_err(|error| CodeError::Unsupported(error.to_string()))?;
                return self.print(&value);
            }
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
        if let Func::Logic(op) = f {
            if matches!(self.lang(), Lang::Python | Lang::NumPy) {
                return self.py_logic(*op, args);
            }
            if self.lang() == Lang::Cuda {
                return self.cuda_logic(*op, args);
            }
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
            Func::Special(sp @ (crate::special::SpecialFn::Mod | crate::special::SpecialFn::Rem)) => {
                return self.c_mod(*sp, args);
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

    /// CUDA C++ boolean operators. Every operand is parenthesized and every
    /// result is wrapped, so no C precedence rule is involved. `Xor`,
    /// `Equivalent` and `Implies` compare the C truth values `!(a)`, which
    /// is exact for `bool` operands and for scalars used as truth values.
    /// `ITE` is a conditional expression: only the selected branch is
    /// evaluated. `&&`/`||` short-circuit; operands are side-effect free.
    fn cuda_logic(&self, op: crate::logic::LogicOp, args: &[Expr]) -> CodeResult {
        use crate::logic::LogicOp;
        let wrap = |e: &Expr| -> CodeResult { Ok(format!("({})", self.print(e)?)) };
        let truth = |e: &Expr| -> CodeResult { Ok(format!("(!{})", wrap(e)?)) };
        let mut sorted = args.to_vec();
        if matches!(op, LogicOp::And | LogicOp::Or | LogicOp::Xor | LogicOp::Equivalent) {
            sorted.sort_by_cached_key(crate::sortkey::sort_key);
        }
        match op {
            LogicOp::And | LogicOp::Or => {
                let sep = if op == LogicOp::And { " && " } else { " || " };
                let parts: Vec<String> = sorted.iter().map(wrap).collect::<Result<_, _>>()?;
                Ok(format!("({})", parts.join(sep)))
            }
            LogicOp::Not => Ok(format!("(!{})", wrap(&args[0])?)),
            LogicOp::Xor => {
                let mut acc = truth(&sorted[0])?;
                for a in &sorted[1..] {
                    acc = format!("({acc} != {})", truth(a)?);
                }
                Ok(acc)
            }
            LogicOp::Equivalent => {
                let truths: Vec<String> = sorted.iter().map(truth).collect::<Result<_, _>>()?;
                let pairs: Vec<String> =
                    truths.windows(2).map(|w| format!("({} == {})", w[0], w[1])).collect();
                Ok(if pairs.is_empty() { "true".into() } else { format!("({})", pairs.join(" && ")) })
            }
            LogicOp::Implies => Ok(format!("({} || {})", truth(&args[0])?, wrap(&args[1])?)),
            LogicOp::ITE => Ok(format!("({} ? {} : {})", wrap(&args[0])?, wrap(&args[1])?, wrap(&args[2])?)),
        }
    }

    /// SymPy's `CodePrinter._print_And/_Or/_Not` (Python keywords) and
    /// `NumPyPrinter._print_And/_Or/_Not` (`numpy.logical_*`).
    fn py_logic(&self, op: crate::logic::LogicOp, args: &[Expr]) -> CodeResult {
        use crate::logic::LogicOp;
        // sympy.printing.precedence: And 30, Or 20, Not 100.
        fn level(e: &Expr) -> i32 {
            match e.kind() {
                Kind::Func(Func::Logic(LogicOp::And), _) => 30,
                Kind::Func(Func::Logic(LogicOp::Or), _) => 20,
                Kind::Func(Func::Logic(LogicOp::Not), _) => 100,
                Kind::Func(Func::Logic(_), _) => 10,
                _ => precedence(e),
            }
        }
        let wrap = |e: &Expr, prec: i32| -> CodeResult {
            let s = self.print(e)?;
            Ok(if level(e) <= prec { format!("({s})") } else { s })
        };
        if self.lang() == Lang::NumPy {
            let name = match op {
                LogicOp::And => "logical_and",
                LogicOp::Or => "logical_or",
                LogicOp::Not => {
                    return Ok(format!("numpy.logical_not({})", self.join_args(args)?));
                }
                other => return Err(CodeError::Unsupported(other.name().into())),
            };
            let parts: Vec<String> = args.iter().map(|a| self.print(a)).collect::<Result<_, _>>()?;
            return Ok(format!("numpy.{name}.reduce(({}))", parts.join(",")));
        }
        match op {
            LogicOp::And | LogicOp::Or => {
                let (word, prec) = if op == LogicOp::And { (" and ", 30) } else { (" or ", 20) };
                let mut sorted = args.to_vec();
                sorted.sort_by_cached_key(crate::sortkey::sort_key);
                let parts: Vec<String> = sorted.iter().map(|a| wrap(a, prec)).collect::<Result<_, _>>()?;
                Ok(parts.join(word))
            }
            LogicOp::Not => Ok(format!("not{}", wrap(&args[0], 100)?)),
            other => Err(CodeError::Unsupported(other.name().into())),
        }
    }

    /// SymPy's `SciPyPrinter` names of `scipy.special` functions.
    fn scipy_special(name: &str) -> Option<&'static str> {
        Some(match name {
            "Ei" => "expi",
            "erf" => "erf",
            "erfc" => "erfc",
            "besselj" => "jv",
            "bessely" => "yv",
            "besseli" => "iv",
            "besselk" => "kv",
            "cosm1" => "cosm1",
            "powm1" => "powm1",
            "factorial" => "factorial",
            "gamma" => "gamma",
            "loggamma" => "gammaln",
            "digamma" => "psi",
            "polygamma" => "polygamma",
            "RisingFactorial" => "poch",
            "jacobi" => "eval_jacobi",
            "gegenbauer" => "eval_gegenbauer",
            "chebyshevt" => "eval_chebyt",
            "chebyshevu" => "eval_chebyu",
            "legendre" => "eval_legendre",
            "hermite" => "eval_hermite",
            "laguerre" => "eval_laguerre",
            "assoc_laguerre" => "eval_genlaguerre",
            "beta" => "beta",
            "LambertW" => "lambertw",
            _ => return None,
        })
    }

    fn py_function(&self, f: &Func, args: &[Expr]) -> CodeResult {
        let np = self.lang() == Lang::NumPy;
        if np && self.s.scipy {
            if let Some(name) = Self::scipy_special(f.name()) {
                return Ok(format!("scipy.special.{name}({})", self.join_args(args)?));
            }
            // SciPyPrinter's own forms of these are not reproduced here.
            if matches!(f.name(), "assoc_legendre" | "lowergamma" | "uppergamma" | "betainc"
                | "betainc_regularized" | "fresnels" | "fresnelc" | "airyai" | "airyaiprime"
                | "airybi" | "airybiprime" | "bernoulli" | "harmonic" | "Si" | "Ci" | "Integral")
            {
                return Err(CodeError::Unsupported(f.name().into()));
            }
        }
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
                    if self.s.lambdify || crate::assumptions::exact_queries(||crate::assumptions::is_real(&args[0])) == Some(true) {
                        Ok(format!("numpy.sign({a})"))
                    } else { Ok(numpy_complex_sign(&a)) }
                } else {
                    Ok(format!("(0.0 if {a} == 0 else math.copysign(1, {a}))"))
                }
            }
            Func::Re if np => Ok(format!("numpy.real({})", self.join_args(args)?)),
            Func::Im if np => Ok(format!("numpy.imag({})", self.join_args(args)?)),
            Func::Arg if np => Ok(format!("numpy.angle({})", self.join_args(args)?)),
            Func::Conjugate if np => Ok(format!("numpy.conjugate({})", self.join_args(args)?)),
            // NumPyPrinter._print_sinc: numpy's normalized sinc of x/pi.
            Func::Special(crate::special::SpecialFn::Sinc) if np && args.len() == 1 => Ok(format!(
                "numpy.sinc({})",
                self.print(&div(&args[0], &Expr::from_kind(Kind::Const(Constant::Pi))))?
            )),
            Func::Special(crate::special::SpecialFn::Mod) if np => Ok(format!("numpy.mod({})", self.join_args(args)?)),
            Func::Special(crate::special::SpecialFn::Mod) => Ok(format!(
                "{} % {}",
                self.parenthesize(&args[0], PREC_MUL, false)?,
                self.parenthesize(&args[1], PREC_MUL, false)?
            )),
            Func::Special(crate::special::SpecialFn::Frac) if !np => {
                Ok(format!("{} % 1", self.parenthesize(&args[0], PREC_MUL, false)?))
            }
            Func::Special(crate::special::SpecialFn::Rem) => {
                Ok(format!("{}({})", if np { "numpy.fmod" } else { "math.fmod" }, self.join_args(args)?))
            }
            Func::Special(crate::special::SpecialFn::KroneckerDelta) if args.len() == 2 && !np => {
                Ok(format!("(1 if {} == {} else 0)", self.print(&args[0])?, self.print(&args[1])?))
            }
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

    #[test]
    fn cuda_negative_quotients_keep_denominators_grouped() {
        let c = |src: &str| cudacode(&parse_str(src).unwrap()).unwrap();
        // regression: a negative coefficient lowered the product precedence
        // and printed -x/s**2 as -x/s*s
        assert_eq!(c("-x/s**2"), "-x/(s*s)");
        assert_eq!(c("x/s**2"), "x/(s*s)");
        assert_eq!(c("-1/s**3"), "-1/(s*s*s)");
        assert_eq!(c("-x/y**3"), "-x/(y*y*y)");
        assert_eq!(c("-a*b/(c*d**2)"), "-a*b/(c*d*d)");
        assert_eq!(c("a*b/(c*d**2)"), "a*b/(c*d*d)");
        let f = CudaSettings { dtype: FloatType::F32, fast_math: true, ..Default::default() };
        assert_eq!(cudacode_with(&parse_str("-x/s**2").unwrap(), &f).unwrap(), "-__fdividef(x, s*s)");
    }

    #[test]
    fn cuda_logic_operators() {
        let c = |src: &str| cudacode(&parse_str(src).unwrap()).unwrap();
        assert_eq!(c("And(x < 1, y > 2)"), "((y > 2) && (x < 1))");
        assert_eq!(c("Or(x < 1, y > 2)"), "((y > 2) || (x < 1))");
        assert_eq!(c("Not(x)"), "(!(x))");
        assert_eq!(c("Xor(x < 1, y > 2)"), "((!(y > 2)) != (!(x < 1)))");
        assert_eq!(c("Implies(x < 1, y > 2)"), "((!(x < 1)) || (y > 2))");
        let pw = c("Piecewise((x, And(x < 1, y > 2)), (y, True))");
        assert!(pw.contains("((y > 2) && (x < 1))"), "{pw}");
    }

    #[test]
    fn cuda_inline_compound_powers() {
        let c = |src: &str| cudacode(&parse_str(src).unwrap()).unwrap();
        assert_eq!(c("(x - y)**2"), "(x - y)*(x - y)");
        assert_eq!(c("(x - y)**3"), "(x - y)*(x - y)*(x - y)");
        assert_eq!(c("(x - y)**5"), "pow(x - y, 5)");
        assert_eq!(c("1/(x + y)**2"), "1.0/((x + y)*(x + y))");
        assert_eq!(c("exp(-(x - y)**2)"), "exp(-((x - y)*(x - y)))");
        // function bases are never repeated
        assert_eq!(c("sin(x)**2"), "pow(sin(x), 2)");
        // inline_pow = 0 disables every inline power
        let s = CudaSettings { inline_pow: 0, ..Default::default() };
        assert_eq!(cudacode_with(&parse_str("(x - y)**2").unwrap(), &s).unwrap(), "pow(x - y, 2)");
        // C is unchanged
        assert_eq!(ccode(&parse_str("(x - y)**2").unwrap()).unwrap(), "pow(x - y, 2)");
    }

    #[test]
    fn cuda_reciprocal_sqrt_and_macros() {
        let e = parse_str("1/sqrt(x + y) + 1").unwrap();
        let plain = CudaSettings { use_rsqrt: false, ..Default::default() };
        assert_eq!(cudacode_with(&e, &plain).unwrap(), "1 + pow(x + y, -1.0/2.0)");
        let r = CudaSettings { use_rsqrt: false, reciprocal_sqrt: true, ..Default::default() };
        assert_eq!(cudacode_with(&e, &r).unwrap(), "1 + 1.0/sqrt(x + y)");
        let r32 = CudaSettings { dtype: FloatType::F32, ..r.clone() };
        assert_eq!(cudacode_with(&e, &r32).unwrap(), "1 + 1.0F/sqrtf(x + y)");
        // use_rsqrt keeps priority (historical default unchanged)
        let both = CudaSettings { reciprocal_sqrt: true, ..Default::default() };
        assert_eq!(cudacode_with(&e, &both).unwrap(), "1 + rsqrt(x + y)");
        let pi = parse_str("pi*x").unwrap();
        assert_eq!(cudacode(&pi).unwrap(), "M_PI*x");
        let nm = CudaSettings { math_macros: Some(false), ..Default::default() };
        assert_eq!(cudacode_with(&pi, &nm).unwrap(), "3.1415926535897932*x");
    }
}


#[cfg(test)]
mod numpy_sign_tests {
    use super::*;
    #[test]
    fn real_sign_stays_short_while_unknown_sign_dispatches_runtime_dtype() {
        let real = Expr::symbol_with("r",crate::Assumptions::REAL);
        assert_eq!(numpycode(&crate::functions::sign(&real)).unwrap(),"numpy.sign(r)");
        let generic = numpycode(&crate::functions::sign(&Expr::symbol("z"))).unwrap();
        assert!(generic.starts_with("(lambda _sgn_value:"));
        assert!(generic.ends_with(")(z)"));
        assert!(generic.contains("numpy.iscomplexobj"));
        assert!(generic.contains("out=numpy.zeros_like"));
        assert!(!generic.contains("symox"));
    }
}

// ----------------------------------------------------------------------
// lambdify source
// ----------------------------------------------------------------------

/// Printer of SymPy's `lambdify` generated source.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LambdaFlavor {
    /// `PythonCodePrinter` (`modules="math"`).
    Math,
    /// `NumPyPrinter` (`modules="numpy"`).
    NumPy,
    /// `SciPyPrinter` (the default modules when SciPy is installed).
    SciPy,
}

/// Module prefixes `_module_format` shortens with `fully_qualified_modules=False`,
/// longest first.
const LAMBDA_MODULES: &[&str] = &[
    "numpy.linalg", "scipy.special", "scipy.constants", "functools", "mpmath", "cmath", "numpy",
    "scipy", "math",
];

fn ident_char(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

/// Remove the module prefixes of `code`, recording the `(module, name)`
/// imports in first-appearance order (as `module_imports` does).
fn strip_modules(code: &str, imports: &mut Vec<(String, String)>) -> String {
    let b: Vec<char> = code.chars().collect();
    let mut out = String::with_capacity(code.len());
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c == '\'' || c == '"' {
            let start = i;
            i += 1;
            while i < b.len() && b[i] != c {
                if b[i] == '\\' {
                    i += 1;
                }
                i += 1;
            }
            i = (i + 1).min(b.len());
            out.extend(&b[start..i]);
            continue;
        }
        let boundary = i == 0 || !(ident_char(b[i - 1]) || b[i - 1] == '.');
        if boundary && (c.is_ascii_alphabetic() || c == '_') {
            let mut j = i;
            while j < b.len() && ident_char(b[j]) {
                j += 1;
            }
            // dotted path starting here
            let mut k = j;
            let mut path: String = b[i..j].iter().collect();
            let mut matched: Option<(usize, String)> = None;
            loop {
                if LAMBDA_MODULES.contains(&path.as_str()) && k < b.len() && b[k] == '.' {
                    matched = Some((k, path.clone()));
                }
                if k < b.len() && b[k] == '.' && k + 1 < b.len() && (b[k + 1].is_ascii_alphabetic() || b[k + 1] == '_') {
                    let mut m = k + 1;
                    while m < b.len() && ident_char(b[m]) {
                        m += 1;
                    }
                    path.push('.');
                    path.extend(&b[k + 1..m]);
                    k = m;
                } else {
                    break;
                }
            }
            if let Some((dot, module)) = matched {
                let mut m = dot + 1;
                while m < b.len() && ident_char(b[m]) {
                    m += 1;
                }
                let name: String = b[dot + 1..m].iter().collect();
                let entry = (module, name.clone());
                if !imports.contains(&entry) {
                    imports.push(entry);
                }
                out.push_str(&name);
                i = m;
                continue;
            }
            out.extend(&b[i..j]);
            i = j;
            continue;
        }
        out.push(c);
        i += 1;
    }
    out
}

/// The expression part of `lambdify`'s generated `return` line and the
/// `from module import name` imports it needs (`fully_qualified_modules=False`,
/// `allow_unknown_functions=True`). Declines (`Err`) for forms whose SymPy
/// printer output is not reproduced.
pub fn lambdify_source(
    e: &Expr,
    flavor: LambdaFlavor,
    user_functions: &FxHashMap<String, String>,
) -> Result<(String, Vec<(String, String)>), CodeError> {
    let mut s = CodeSettings::new(if flavor == LambdaFlavor::Math { Lang::Python } else { Lang::NumPy });
    s.allow_unknown_functions = true;
    s.lambdify = true;
    s.scipy = flavor == LambdaFlavor::SciPy;
    s.user_functions = user_functions.clone();
    let code = code_with(e, &s)?;
    if code.contains('\n') {
        return Err(CodeError::Unsupported("multi-line code".into()));
    }
    let mut imports = Vec::new();
    let body = strip_modules(&code, &mut imports);
    Ok((body, imports))
}

/// Type of a `lambdify` generated function's result when every argument is
/// a Python `float`, as far as it follows from the printed source alone.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LambdaKind {
    /// Not determined: run the generated function itself.
    Unknown = 0,
    /// Python `float`.
    Float = 1,
    /// `numpy.float64`.
    NpFloat = 2,
    /// `numpy.bool_`.
    NpBool = 3,
    /// Python `bool`.
    Bool = 4,
    /// 0-d `numpy.ndarray` of float64 (`select`).
    Array0 = 5,
    /// Python `int` (internal: constant parts).
    Int = 6,
}

/// `head(...)` spanning all of `s`: the dotted head.
fn whole_call(s: &str) -> Option<&str> {
    let open = s.find('(')?;
    let head = &s[..open];
    if head.is_empty() || !head.chars().all(|c| ident_char(c) || c == '.') || !s.ends_with(')') {
        return None;
    }
    let mut depth = 0i32;
    let mut quote: Option<char> = None;
    for (i, c) in s.char_indices().skip(open) {
        if let Some(q) = quote {
            if c == q {
                quote = None;
            }
            continue;
        }
        match c {
            '\'' | '"' => quote = Some(c),
            '(' | '[' => depth += 1,
            ')' | ']' => {
                depth -= 1;
                if depth == 0 && i + 1 != s.len() {
                    return None;
                }
            }
            _ => {}
        }
    }
    Some(head)
}

/// Result type of `lambdify(args, e, modules)`'s generated function for
/// Python `float` arguments (`Unknown` whenever it is not certain).
pub fn lambdify_result_kind(e: &Expr, flavor: LambdaFlavor) -> LambdaKind {
    let mut s = CodeSettings::new(if flavor == LambdaFlavor::Math { Lang::Python } else { Lang::NumPy });
    s.allow_unknown_functions = true;
    s.lambdify = true;
    s.scipy = flavor == LambdaFlavor::SciPy;
    let _session = super::order::session();
    let p = CodePrinter::new(&s);
    let math = flavor == LambdaFlavor::Math;
    fn arith(kinds: impl Iterator<Item = LambdaKind>) -> LambdaKind {
        let mut out = LambdaKind::Int;
        for k in kinds {
            out = match (out, k) {
                (_, LambdaKind::Unknown | LambdaKind::Bool | LambdaKind::NpBool | LambdaKind::Array0) => {
                    return LambdaKind::Unknown
                }
                (LambdaKind::NpFloat, _) | (_, LambdaKind::NpFloat) => LambdaKind::NpFloat,
                (LambdaKind::Float, _) | (_, LambdaKind::Float) => LambdaKind::Float,
                _ => LambdaKind::Int,
            };
        }
        out
    }
    fn num(n: &Number) -> LambdaKind {
        match n {
            Number::Float(_) => LambdaKind::Float,
            n if n.is_integer() => LambdaKind::Int,
            n if n.is_rational() => LambdaKind::Float,
            _ => LambdaKind::Unknown,
        }
    }
    fn go(p: &CodePrinter, e: &Expr, math: bool, top: bool) -> LambdaKind {
        use crate::logic::LogicOp;
        match e.kind() {
            Kind::Sym(_) => LambdaKind::Float,
            Kind::Num(n) => num(n),
            Kind::Const(Constant::Pi | Constant::E) => LambdaKind::Float,
            Kind::Add(c, terms) => {
                let k = arith(terms.iter().flat_map(|(t, c)| [go(p, t, math, false), num(c)]));
                if c.is_zero() { k } else { arith([k, num(c)].into_iter()) }
            }
            Kind::Pow(b, x) => {
                let Ok(s) = p.print(e) else { return LambdaKind::Unknown };
                if s.contains("numpy.") || s.contains("math.") || s.contains("scipy.") {
                    // function forms (sqrt, exp, 1/sqrt, ...): only a single call
                    return match whole_call(&s) {
                        Some(head) if go(p, b, math, false) != LambdaKind::Unknown => call_kind(head, math),
                        _ => LambdaKind::Unknown,
                    };
                }
                if whole_call(&s).is_some() {
                    return LambdaKind::Unknown;
                }
                let k = arith([go(p, b, math, false), go(p, x, math, false)].into_iter());
                // Python `int ** negative int` is a float; keep it unknown
                if k == LambdaKind::Int { LambdaKind::Unknown } else { k }
            }
            Kind::Mul(c, factors) => {
                let k = arith(factors.iter().map(|(b, x)| {
                    if x.is_one() {
                        go(p, b, math, false)
                    } else {
                        go(p, &Expr::from_kind(Kind::Pow(b.clone(), x.clone())), math, false)
                    }
                }));
                let k = if k == LambdaKind::Int { LambdaKind::Unknown } else { k };
                arith([k, num(c)].into_iter())
            }
            Kind::Rel(_, a, b) => {
                let k = arith([go(p, a, math, false), go(p, b, math, false)].into_iter());
                if k == LambdaKind::Unknown {
                    LambdaKind::Unknown
                } else if math {
                    LambdaKind::Bool
                } else {
                    LambdaKind::NpBool
                }
            }
            Kind::Func(Func::Logic(op @ (LogicOp::And | LogicOp::Or | LogicOp::Not)), args) => {
                let want = if math { LambdaKind::Bool } else { LambdaKind::NpBool };
                if !args.iter().all(|a| go(p, a, math, false) == want) {
                    return LambdaKind::Unknown;
                }
                let _ = op;
                want
            }
            Kind::Func(..) => {
                let Ok(s) = p.print(e) else { return LambdaKind::Unknown };
                match whole_call(&s) {
                    Some("numpy.real" | "numpy.imag") if !math => {
                        let Kind::Func(_, args) = e.kind() else { unreachable!() };
                        let k = go(p, &args[0], math, false);
                        if matches!(k, LambdaKind::Float | LambdaKind::NpFloat) { k } else { LambdaKind::Unknown }
                    }
                    Some(head) => {
                        let Kind::Func(_, args) = e.kind() else { unreachable!() };
                        if args.iter().all(|a| matches!(go(p, a, math, false), LambdaKind::Float | LambdaKind::NpFloat | LambdaKind::Int)) {
                            call_kind(head, math)
                        } else {
                            LambdaKind::Unknown
                        }
                    }
                    None => LambdaKind::Unknown,
                }
            }
            Kind::Piecewise(_) if top && !math => {
                let Ok(s) = p.print(e) else { return LambdaKind::Unknown };
                if whole_call(&s) == Some("numpy.select") { LambdaKind::Array0 } else { LambdaKind::Unknown }
            }
            _ => LambdaKind::Unknown,
        }
    }
    /// Result type of a module call on float arguments; `math.*` results
    /// (the C library's last digits) are left to the generated function.
    fn call_kind(head: &str, math: bool) -> LambdaKind {
        if math {
            return match head {
                "math.sin" | "math.cos" | "math.tan" | "math.asin" | "math.acos" | "math.atan"
                | "math.atan2" | "math.sinh" | "math.cosh" | "math.tanh" | "math.asinh"
                | "math.acosh" | "math.atanh" | "math.exp" | "math.log" | "math.sqrt"
                | "math.log10" | "math.log2" | "math.expm1" | "math.log1p" | "math.hypot" => LambdaKind::Float,
                _ => LambdaKind::Unknown,
            };
        }
        if head == "numpy.select" || head == "numpy.sinc" || head.starts_with("numpy.logical") || head.starts_with("numpy.is") {
            return LambdaKind::Unknown;
        }
        // scipy.special values (gamma, erf, Bessel, ...) can differ from the
        // native ones in the last digits: the generated function computes them.
        // except gamma, which `crate::cephes_gamma` gives bit for bit.
        if head.starts_with("numpy.") || head == "scipy.special.gamma" {
            LambdaKind::NpFloat
        } else {
            LambdaKind::Unknown
        }
    }
    if crate::traverse::free_symbols(e).is_empty() {
        // constants keep their Python type (and are not broadcast)
        return LambdaKind::Unknown;
    }
    match go(&p, e, math, true) {
        LambdaKind::Int => LambdaKind::Unknown,
        k => k,
    }
}

/// Result type of `lambdify`'s generated function when some arguments are
/// Python `int`s (the others `float`s) of magnitude at most the returned
/// bound, and when the native double evaluation then gives its value. The
/// generated code keeps `int` arithmetic exact: the bound keeps every integer
/// intermediate below 2**53, where it is a double; a module function then
/// converts it exactly. Functions that keep integer types (`numpy.sign`,
/// `abs`, ...), negative or symbolic powers of integers and `select` give
/// `Unknown`.
pub fn lambdify_int_result_kind(e: &Expr, flavor: LambdaFlavor) -> (LambdaKind, u64) {
    const EXACT: f64 = 9007199254740992.0; // 2**53
    let mut s = CodeSettings::new(if flavor == LambdaFlavor::Math { Lang::Python } else { Lang::NumPy });
    s.allow_unknown_functions = true;
    s.lambdify = true;
    s.scipy = flavor == LambdaFlavor::SciPy;
    let _session = super::order::session();
    let p = CodePrinter::new(&s);
    let math = flavor == LambdaFlavor::Math;
    fn converts(head: &str) -> bool {
        let Some((module, name)) = head.split_once('.') else { return false };
        (module == "numpy" || module == "math") && matches!(name,
            "sin" | "cos" | "tan" | "sinh" | "cosh" | "tanh" | "exp" | "log" | "sqrt" | "log10"
            | "log2" | "expm1" | "log1p" | "arcsin" | "arccos" | "arctan" | "arcsinh" | "arccosh"
            | "arctanh" | "asin" | "acos" | "atan" | "asinh" | "acosh" | "atanh")
    }
    fn join(a: LambdaKind, b: LambdaKind) -> LambdaKind {
        use LambdaKind::*;
        match (a, b) {
            (Unknown | Bool | NpBool | Array0, _) | (_, Unknown | Bool | NpBool | Array0) => Unknown,
            (NpFloat, _) | (_, NpFloat) => NpFloat,
            (Float, _) | (_, Float) => Float,
            _ => Int,
        }
    }
    fn num(n: &Number) -> (LambdaKind, f64) {
        match n {
            Number::Float(_) => (LambdaKind::Float, 0.0),
            n if n.is_integer() => (LambdaKind::Int, n.to_f64().abs()),
            n if n.is_rational() => (LambdaKind::Float, 0.0),
            _ => (LambdaKind::Unknown, 0.0),
        }
    }
    struct Walk<'a, 'b> {
        p: &'a CodePrinter<'b>,
        math: bool,
        bound: f64,
        exact: bool,
        budget: usize,
    }
    impl Walk<'_, '_> {
        // (kind, magnitude bound of an `Int` value)
        fn go(&mut self, e: &Expr) -> (LambdaKind, f64) {
            use LambdaKind::*;
            if self.budget == 0 {
                return (Unknown, 0.0);
            }
            self.budget -= 1;
            let (k, m) = self.kind(e);
            if k == Int && !(m < EXACT) {
                self.exact = false;
            }
            (k, m)
        }
        fn call(&mut self, head: &str, args: &[Expr]) -> LambdaKind {
            let mut int = false;
            for a in args {
                match self.go(a).0 {
                    LambdaKind::Int => int = true,
                    LambdaKind::Float | LambdaKind::NpFloat => {}
                    _ => return LambdaKind::Unknown,
                }
            }
            if int && !converts(head) {
                return LambdaKind::Unknown;
            }
            match head {
                "math.sin" | "math.cos" | "math.tan" | "math.asin" | "math.acos" | "math.atan"
                | "math.sinh" | "math.cosh" | "math.tanh" | "math.asinh" | "math.acosh"
                | "math.atanh" | "math.exp" | "math.log" | "math.sqrt" | "math.log10" | "math.log2"
                | "math.expm1" | "math.log1p" if self.math => LambdaKind::Float,
                _ if self.math => LambdaKind::Unknown,
                "numpy.select" | "numpy.sinc" | "numpy.real" | "numpy.imag" => LambdaKind::Unknown,
                h if h.starts_with("numpy.logical") || h.starts_with("numpy.is") => LambdaKind::Unknown,
                h if h.starts_with("numpy.") => LambdaKind::NpFloat,
                _ => LambdaKind::Unknown,
            }
        }
        fn kind(&mut self, e: &Expr) -> (LambdaKind, f64) {
            use LambdaKind::*;
            match e.kind() {
                Kind::Sym(_) => (Int, self.bound),
                Kind::Num(n) => num(n),
                Kind::Const(Constant::Pi | Constant::E) => (Float, 0.0),
                Kind::Add(c, terms) => {
                    let (mut k, mut m) = num(c);
                    if c.is_zero() {
                        (k, m) = (Int, 0.0);
                    }
                    // partial sums of the integer terms are integers too
                    for (t, c) in terms {
                        let (kt, mt) = self.go(t);
                        let (kc, mc) = num(c);
                        let kt = join(kt, kc);
                        if kt == Int {
                            m += mt * mc;
                        }
                        k = join(k, kt);
                    }
                    if !(m < EXACT) {
                        self.exact = false;
                    }
                    (k, m)
                }
                Kind::Mul(c, factors) => {
                    let (mut k, mut m) = num(c);
                    m = m.max(1.0);
                    for (b, x) in factors {
                        let (kf, mf) = if x.is_one() {
                            self.go(b)
                        } else {
                            self.go(&Expr::from_kind(Kind::Pow(b.clone(), x.clone())))
                        };
                        if kf == Int {
                            m *= mf.max(1.0);
                        }
                        k = join(k, kf);
                    }
                    if !(m < EXACT) {
                        self.exact = false;
                    }
                    (k, m)
                }
                Kind::Pow(b, x) => {
                    let Ok(s) = self.p.print(e) else { return (Unknown, 0.0) };
                    if s.contains("numpy.") || s.contains("math.") || s.contains("scipy.") {
                        return match whole_call(&s) {
                            Some(head) => (self.call(head, &[b.clone(), x.clone()]), 0.0),
                            None => (Unknown, 0.0),
                        };
                    }
                    if whole_call(&s).is_some() {
                        return (Unknown, 0.0);
                    }
                    let (kb, mb) = self.go(b);
                    let (kx, _) = self.go(x);
                    if kb == Int {
                        // `int ** int` is a float for negative exponents
                        return match x.kind() {
                            Kind::Num(n) if n.is_integer() && !n.is_negative() => match n.to_i64() {
                                Some(n) => (Int, mb.max(1.0).powi(n.min(1 << 20) as i32)),
                                None => (Unknown, 0.0),
                            },
                            _ => (Unknown, 0.0),
                        };
                    }
                    match join(kb, kx) {
                        Int => (Unknown, 0.0),
                        k => (k, 0.0),
                    }
                }
                Kind::Rel(_, a, b) => {
                    let k = join(self.go(a).0, self.go(b).0);
                    if k == Unknown {
                        (Unknown, 0.0)
                    } else if self.math {
                        (Bool, 0.0)
                    } else {
                        (NpBool, 0.0)
                    }
                }
                Kind::Func(Func::Logic(crate::logic::LogicOp::And | crate::logic::LogicOp::Or | crate::logic::LogicOp::Not), args) => {
                    let want = if self.math { Bool } else { NpBool };
                    for a in args {
                        if self.go(a).0 != want {
                            return (Unknown, 0.0);
                        }
                    }
                    (want, 0.0)
                }
                Kind::Func(_, args) => {
                    let Ok(s) = self.p.print(e) else { return (Unknown, 0.0) };
                    match whole_call(&s) {
                        Some(head) => {
                            let head = head.to_string();
                            (self.call(&head, args), 0.0)
                        }
                        None => (Unknown, 0.0),
                    }
                }
                _ => (Unknown, 0.0),
            }
        }
    }
    if crate::traverse::free_symbols(e).is_empty() {
        return (LambdaKind::Unknown, 0);
    }
    let run = |bound: f64| {
        let mut w = Walk { p: &p, math, bound, exact: true, budget: 2000 };
        let (k, _) = w.go(e);
        (if w.budget == 0 { LambdaKind::Unknown } else { k }, w.exact)
    };
    let (kind, exact) = run(1.0);
    if matches!(kind, LambdaKind::Unknown | LambdaKind::Int) || !exact {
        return (LambdaKind::Unknown, 0);
    }
    // the largest power of two bound that keeps the integers exact
    let (mut lo, mut hi) = (0u32, 53u32);
    while lo < hi {
        let mid = (lo + hi).div_ceil(2);
        if run(2f64.powi(mid as i32)).1 { lo = mid } else { hi = mid - 1 }
    }
    (kind, 1u64 << lo)
}

/// Size of the tree printed for `e` (with repetitions), counted up to
/// `limit + 1`, and whether it holds an undefined function application.
pub fn lambdify_tree_size(e: &Expr, limit: usize) -> (usize, bool) {
    let mut stack = vec![e.clone()];
    let (mut count, mut undefined) = (0usize, false);
    while let Some(node) = stack.pop() {
        count += 1;
        if count > limit {
            break;
        }
        undefined |= matches!(node.kind(), Kind::Func(Func::Undefined(_), _));
        stack.extend(node.args());
    }
    (count, undefined)
}

/// The global names of a printed `lambdify` return expression (Python's
/// `(?<![\w.])[A-Za-z_]\w*(?!\w|\s*=(?!=))` outside string literals),
/// without `known` names and Python keywords, in first-occurrence order.
pub fn lambdify_global_names(body: &str, known: &[&str]) -> Vec<String> {
    const KEYWORDS: &[&str] = &["False", "None", "True", "and", "as", "assert", "async", "await",
        "break", "class", "continue", "def", "del", "elif", "else", "except", "finally", "for",
        "from", "global", "if", "import", "in", "is", "lambda", "nonlocal", "not", "or", "pass",
        "raise", "return", "try", "while", "with", "yield"];
    let word = |c: char| c.is_alphanumeric() || c == '_';
    // `'[^']*'|"[^"]*"` removed, as `re.sub` does
    let chars: Vec<char> = body.chars().collect();
    let mut text = Vec::with_capacity(chars.len());
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c == '\'' || c == '"' {
            if let Some(end) = chars[i + 1..].iter().position(|&d| d == c) {
                i += end + 2;
                continue;
            }
        }
        text.push(c);
        i += 1;
    }
    let mut out: Vec<String> = Vec::new();
    let mut i = 0;
    while i < text.len() {
        if !word(text[i]) {
            i += 1;
            continue;
        }
        let start = i;
        while i < text.len() && word(text[i]) {
            i += 1;
        }
        let first = text[start];
        if start > 0 && text[start - 1] == '.' || !(first.is_ascii_alphabetic() || first == '_') {
            continue;
        }
        let mut j = i;
        while j < text.len() && text[j].is_whitespace() {
            j += 1;
        }
        if j < text.len() && text[j] == '=' && text.get(j + 1) != Some(&'=') {
            continue;
        }
        let name: String = text[start..i].iter().collect();
        if !known.contains(&name.as_str()) && !KEYWORDS.contains(&name.as_str()) && !out.contains(&name) {
            out.push(name);
        }
    }
    out
}
