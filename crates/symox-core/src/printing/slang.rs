//! Scalar, real-valued Slang expressions, independent of runtime libraries.
//!
//! Scalar conditional expressions short circuit in Slang. Keeping Piecewise
//! branches inside those conditionals avoids evaluating inactive domains.

use std::collections::BTreeMap;

use super::CodeError;
use crate::expr::{Constant, Expr, Func, Kind};
use crate::logic::LogicOp;
use crate::number::Number;

type Result = std::result::Result<String, CodeError>;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SlangDtype {
    #[default]
    Float32,
    Float64,
}

#[derive(Clone, Debug, Default)]
pub struct SlangSettings {
    pub dtype: SlangDtype,
    pub symbol_names: BTreeMap<String, String>,
    pub exact_symbol_names: BTreeMap<Expr, String>,
    pub user_functions: BTreeMap<String, String>,
    pub expression_overrides: BTreeMap<Expr, String>,
}

/// Conservative identifier validation, excluding language and builtin names.
pub fn is_slang_identifier(name: &str) -> bool {
    let mut chars = name.chars();
    let valid = matches!(chars.next(), Some(c) if c == '_' || c.is_ascii_alphabetic())
        && chars.all(|c| c == '_' || c.is_ascii_alphanumeric());
    if !valid || name.starts_with("__") {
        return false;
    }
    let reserved = [
        "true",
        "false",
        "if",
        "else",
        "for",
        "while",
        "do",
        "switch",
        "case",
        "default",
        "break",
        "continue",
        "return",
        "discard",
        "struct",
        "class",
        "interface",
        "extension",
        "enum",
        "typedef",
        "typealias",
        "typename",
        "let",
        "var",
        "func",
        "import",
        "module",
        "namespace",
        "using",
        "public",
        "private",
        "internal",
        "static",
        "const",
        "constexpr",
        "uniform",
        "extern",
        "export",
        "in",
        "out",
        "inout",
        "ref",
        "this",
        "This",
        "Self",
        "void",
        "vector",
        "matrix",
        "sizeof",
        "alignof",
        "operator",
        "property",
        "get",
        "set",
        "generic",
        "where",
        "associatedtype",
        "__init",
        "__subscript",
        "mutating",
        "groupshared",
        "register",
        "packoffset",
        "cbuffer",
        "tbuffer",
        "texture",
        "sampler",
        "RWStructuredBuffer",
        "StructuredBuffer",
        "RWByteAddressBuffer",
        "ByteAddressBuffer",
        "Texture2D",
        "SamplerState",
        "sin",
        "cos",
        "tan",
        "asin",
        "acos",
        "atan",
        "atan2",
        "sinh",
        "cosh",
        "tanh",
        "asinh",
        "acosh",
        "atanh",
        "exp",
        "exp2",
        "log",
        "log2",
        "sqrt",
        "rsqrt",
        "pow",
        "abs",
        "sign",
        "ceil",
        "floor",
        "min",
        "max",
        "clamp",
        "select",
        "fmod",
        "asfloat",
        "asdouble",
        "uint",
        "int",
        "int64_t",
        "uint64_t",
        "int32_t",
        "uint32_t",
        "int16_t",
        "uint16_t",
        "int8_t",
        "uint8_t",
        "float",
        "half",
        "double",
        "bool",
        "float16_t",
        "float32_t",
        "float64_t",
    ];
    if reserved.contains(&name) {
        return false;
    }
    // Builtin vector and matrix types, including generic-size aliases.
    for prefix in ["float", "double", "half", "int", "uint", "bool"] {
        if let Some(suffix) = name.strip_prefix(prefix) {
            if !suffix.is_empty() && suffix.chars().all(|c| c.is_ascii_digit() || c == 'x') {
                return false;
            }
        }
    }
    true
}

fn qualified_function(name: &str) -> bool {
    !name.is_empty()
        && name.split("::").all(|part| {
            let mut chars = part.chars();
            matches!(chars.next(), Some(c) if c == '_' || c.is_ascii_alphabetic())
                && chars.all(|c| c == '_' || c.is_ascii_alphanumeric())
        })
}

pub fn slangcode(expr: &Expr) -> Result {
    slangcode_with(expr, &SlangSettings::default())
}

pub fn slangcode_with(expr: &Expr, settings: &SlangSettings) -> Result {
    // Distinct assumed symbols must not silently alias after textual printing.
    let mut seen = BTreeMap::new();
    for symbol in crate::traverse::free_symbols(expr) {
        let source = symbol.as_symbol().unwrap();
        let target = settings
            .exact_symbol_names
            .get(&symbol)
            .or_else(|| settings.symbol_names.get(source))
            .map(String::as_str)
            .unwrap_or(source);
        if seen
            .insert(target.to_string(), symbol.clone())
            .is_some_and(|other| other != symbol)
        {
            return Err(CodeError::Unsupported(format!(
                "ambiguous Slang symbol {target:?}; supply exact symbol_names"
            )));
        }
    }
    SlangPrinter::new(settings).print(expr)
}

/// The scalar result type used for function outputs and CSE declarations.
pub fn slang_type(expr: &Expr) -> std::result::Result<&'static str, CodeError> {
    match expr.kind() {
        Kind::Bool(_) | Kind::Rel(..) | Kind::Func(Func::Logic(_), _) => Ok("bool"),
        Kind::Piecewise(pairs) => {
            let mut ty = None;
            for (value, _) in pairs {
                let current = slang_type(value)?;
                if ty.is_some_and(|previous| previous != current) {
                    return Err(CodeError::Unsupported(
                        "mixed boolean/numeric Piecewise branches in Slang".into(),
                    ));
                }
                ty = Some(current);
            }
            Ok(ty.unwrap_or("real"))
        }
        Kind::Func(Func::Undefined(name), args)
            if name.as_ref() == "UnevaluatedExpr" && args.len() == 1 =>
        {
            slang_type(&args[0])
        }
        _ => Ok("real"),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ScalarBackend {
    Slang,
    Sycl,
}

pub struct SlangPrinter<'a> {
    settings: &'a SlangSettings,
    backend: ScalarBackend,
}

impl<'a> SlangPrinter<'a> {
    pub fn new(settings: &'a SlangSettings) -> Self {
        Self {
            settings,
            backend: ScalarBackend::Slang,
        }
    }
    pub(crate) fn sycl(settings: &'a SlangSettings) -> Self {
        Self {
            settings,
            backend: ScalarBackend::Sycl,
        }
    }
    fn unsupported(&self, item: &str) -> CodeError {
        let language = if self.backend == ScalarBackend::Slang {
            "Slang"
        } else {
            "SYCL"
        };
        CodeError::Unsupported(format!("{item} in real scalar {language} code generation"))
    }
    fn math(&self, name: &str) -> String {
        if self.backend == ScalarBackend::Slang {
            return name.into();
        }
        format!(
            "sycl::{}",
            match name {
                "abs" => "fabs",
                "min" => "fmin",
                "max" => "fmax",
                _ => name,
            }
        )
    }
    fn call(&self, name: &str, args: &[&str]) -> String {
        format!("{}({})", self.math(name), args.join(", "))
    }
    fn literal(&self, text: &str) -> String {
        let point = if text.contains('.') || text.contains('e') || text.contains('E') {
            ""
        } else {
            ".0"
        };
        let suffix = match (self.backend, self.settings.dtype) {
            (ScalarBackend::Slang, SlangDtype::Float64) => "L",
            (ScalarBackend::Sycl, SlangDtype::Float32) => "F",
            _ => "",
        };
        format!("{text}{point}{suffix}")
    }
    fn nan(&self) -> String {
        self.cast_bits("0x7fc00000u")
    }
    fn cast_bits(&self, bits: &str) -> String {
        if self.backend == ScalarBackend::Sycl {
            let ty = if self.settings.dtype == SlangDtype::Float64 {
                "double"
            } else {
                "float"
            };
            let value = if bits == "0x7fc00000u" {
                "quiet_NaN"
            } else {
                "infinity"
            };
            let sign = if bits == "0xff800000u" { "-" } else { "" };
            return format!("{sign}std::numeric_limits<{ty}>::{value}()");
        }
        let code = format!("asfloat({bits})");
        if self.settings.dtype == SlangDtype::Float64 {
            format!("double({code})")
        } else {
            code
        }
    }
    fn number(&self, n: &Number) -> String {
        if let Some(code) = self
            .settings
            .expression_overrides
            .get(&Expr::num(n.clone()))
        {
            return format!("({code})");
        }
        match n {
            Number::Rational(r) => format!(
                "({}/{})",
                self.literal(&r.numer().to_string()),
                self.literal(&r.denom().to_string())
            ),
            Number::Float(f) if f.is_nan() => self.nan(),
            Number::Float(f) if f.is_infinite() => self.cast_bits(if *f > 0.0 {
                "0x7f800000u"
            } else {
                "0xff800000u"
            }),
            Number::Float(f) => self.literal(&format!("{f:?}")),
            _ => self.literal(&n.to_string()),
        }
    }
    pub fn print(&self, expr: &Expr) -> Result {
        if let Some(code) = self.settings.expression_overrides.get(expr) {
            return Ok(format!("({code})"));
        }
        match expr.kind() {
            Kind::Num(n) => Ok(self.number(n)),
            Kind::Bool(value) => Ok(if *value { "true" } else { "false" }.into()),
            Kind::Const(c) => Ok(match c {
                Constant::Pi => self.literal("3.141592653589793238462643383279502884"),
                Constant::E => self.literal("2.718281828459045235360287471352662498"),
                Constant::EulerGamma => self.literal("0.577215664901532860606512090082402431"),
                Constant::Infinity => self.cast_bits("0x7f800000u"),
                Constant::NegInfinity => self.cast_bits("0xff800000u"),
                Constant::NaN => self.nan(),
                _ => return Err(self.unsupported("complex constant")),
            }),
            Kind::Sym(s) => {
                let name = self
                    .settings
                    .exact_symbol_names
                    .get(expr)
                    .or_else(|| self.settings.symbol_names.get(s.name()))
                    .map(String::as_str)
                    .unwrap_or(s.name());
                if !(match self.backend {
                    ScalarBackend::Slang => is_slang_identifier(name),
                    ScalarBackend::Sycl => super::sycl::is_sycl_identifier(name),
                }) {
                    return Err(self.unsupported(&format!(
                        "unsafe symbol name {name:?}; supply symbol_names"
                    )));
                }
                Ok(name.into())
            }
            Kind::Add(c, terms) => {
                let mut parts = Vec::new();
                if !c.is_zero() {
                    parts.push(self.number(c));
                }
                for (term, coeff) in terms {
                    let term = self.print(term)?;
                    parts.push(if coeff.is_one() {
                        term
                    } else {
                        format!("({}*{term})", self.number(coeff))
                    });
                }
                Ok(if parts.is_empty() {
                    self.literal("0")
                } else {
                    format!("({})", parts.join(" + "))
                })
            }
            Kind::Mul(c, factors) => {
                let mut parts = Vec::new();
                if !c.is_one() {
                    parts.push(self.number(c));
                }
                for (base, exponent) in factors {
                    parts.push(if exponent.is_one() {
                        self.print(base)?
                    } else {
                        self.power(base, exponent)?
                    });
                }
                Ok(if parts.is_empty() {
                    self.literal("1")
                } else {
                    format!("({})", parts.join("*"))
                })
            }
            Kind::Pow(base, exponent) => self.power(base, exponent),
            Kind::Rel(op, a, b) => Ok(format!(
                "({} {} {})",
                self.print(a)?,
                op.symbol(),
                self.print(b)?
            )),
            Kind::Piecewise(pairs) => {
                let boolean = slang_type(expr)? == "bool";
                let mut fallback = if boolean {
                    return self.boolean_piecewise(pairs);
                } else {
                    self.nan()
                };
                for (value, cond) in pairs.iter().rev() {
                    let value = self.print(value)?;
                    if matches!(cond.kind(), Kind::Bool(true)) {
                        fallback = value;
                    } else {
                        fallback = format!("({} ? {value} : {fallback})", self.print(cond)?);
                    }
                }
                Ok(fallback)
            }
            Kind::Func(f, args) => self.function(f, args),
            Kind::Derivative(..) => Err(self.unsupported("unevaluated Derivative")),
        }
    }
    fn boolean_piecewise(&self, pairs: &[(Expr, Expr)]) -> Result {
        if !pairs
            .iter()
            .any(|(_, c)| matches!(c.kind(), Kind::Bool(true)))
        {
            return Err(CodeError::PiecewiseNoDefault);
        }
        let mut out = "false".to_string();
        for (value, cond) in pairs.iter().rev() {
            let value = self.print(value)?;
            if matches!(cond.kind(), Kind::Bool(true)) {
                out = value;
            } else {
                out = format!("({} ? {value} : {out})", self.print(cond)?);
            }
        }
        Ok(out)
    }
    fn power(&self, base: &Expr, exponent: &Expr) -> Result {
        let b = self.print(base)?;
        let one = self.literal("1");
        if let Some(n) = exponent.as_number().and_then(Number::to_i64) {
            let magnitude = n.unsigned_abs();
            let positive = if magnitude == 0 {
                one.clone()
            } else if magnitude <= 16 {
                format!(
                    "({})",
                    std::iter::repeat(format!("({b})"))
                        .take(magnitude as usize)
                        .collect::<Vec<_>>()
                        .join("*")
                )
            } else {
                let absolute = self.call("abs", &[&b]);
                let absolute_power =
                    self.call("pow", &[&absolute, &self.literal(&magnitude.to_string())]);
                if magnitude % 2 == 1 {
                    format!(
                        "(({b}) < {} ? -{absolute_power} : {absolute_power})",
                        self.literal("0")
                    )
                } else {
                    absolute_power
                }
            };
            return Ok(if n < 0 {
                format!("({one}/({positive}))")
            } else {
                positive
            });
        }
        if exponent
            .as_number()
            .is_some_and(|n| n == &Number::rational_i64(1, 2).unwrap())
        {
            return Ok(self.call("sqrt", &[&b]));
        }
        let e = self.print(exponent)?;
        // Restore integer negative-base powers explicitly for portable targets.
        let floor = self.call("floor", &[&e]);
        let abs = self.call("abs", &[&e]);
        let modulo = self.call("fmod", &[&abs, &self.literal("2")]);
        let negative_power = self.call("pow", &[&format!("-({b})"), &e]);
        let positive_power = self.call("pow", &[&b, &e]);
        let zero = self.literal("0");
        Ok(format!("(({b}) < {zero} ? (({e}) == {floor} ? (({modulo} == {zero} ? {one} : -{one})*{negative_power}) : {}) : {positive_power})", self.nan()))
    }
    fn function(&self, f: &Func, args: &[Expr]) -> Result {
        let values = args
            .iter()
            .map(|a| self.print(a))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        if let Some(target) = self.settings.user_functions.get(f.name()) {
            if !qualified_function(target) {
                return Err(self.unsupported("unsafe user function target"));
            }
            return Ok(format!("{target}({})", values.join(", ")));
        }
        let unary = match f {
            Func::Sin => Some("sin"),
            Func::Cos => Some("cos"),
            Func::Tan => Some("tan"),
            Func::Asin => Some("asin"),
            Func::Acos => Some("acos"),
            Func::Atan => Some("atan"),
            Func::Sinh => Some("sinh"),
            Func::Cosh => Some("cosh"),
            Func::Tanh => Some("tanh"),
            Func::Asinh => Some("asinh"),
            Func::Acosh => Some("acosh"),
            Func::Atanh => Some("atanh"),
            Func::Exp => Some("exp"),
            Func::Log => Some("log"),
            Func::Abs => Some("abs"),
            Func::Floor => Some("floor"),
            Func::Ceiling => Some("ceil"),
            _ => None,
        };
        if let Some(name) = unary {
            return Ok(self.call(name, &[&values[0]]));
        }
        if self.backend == ScalarBackend::Sycl {
            match f {
                Func::Erf => return Ok(self.call("erf", &[&values[0]])),
                Func::Erfc => return Ok(self.call("erfc", &[&values[0]])),
                Func::Gamma => return Ok(self.call("tgamma", &[&values[0]])),
                Func::LogGamma => {
                    if crate::assumptions::is_positive(&args[0]) != Some(true) {
                        return Err(
                            self.unsupported("loggamma requires a provably positive real argument")
                        );
                    }
                    return Ok(self.call("lgamma", &[&values[0]]));
                }
                _ => {}
            }
        }
        let zero = self.literal("0");
        let one = self.literal("1");
        match f {
            Func::Undefined(name) if name.as_ref() == "UnevaluatedExpr" && values.len() == 1 => {
                Ok(values[0].clone())
            }
            Func::Atan2 => Ok(self.call("atan2", &[&values[0], &values[1]])),
            Func::Min | Func::Max => {
                let name = if *f == Func::Min { "min" } else { "max" };
                let mut it = values.into_iter();
                let mut result = it.next().ok_or_else(|| self.unsupported("empty Min/Max"))?;
                for value in it {
                    result = self.call(name, &[&result, &value]);
                }
                Ok(result)
            }
            Func::Cot | Func::Sec | Func::Csc | Func::Coth => {
                let target = match f {
                    Func::Cot => "tan",
                    Func::Sec => "cos",
                    Func::Csc => "sin",
                    _ => "tanh",
                };
                Ok(format!("({one}/{})", self.call(target, &[&values[0]])))
            }
            Func::Acot => Ok(format!(
                "({} == {zero} ? {} : {})",
                values[0],
                self.literal("1.570796326794896619231321691639751442"),
                self.call("atan", &[&format!("{one}/({})", values[0])])
            )),
            Func::Heaviside => {
                let x = &values[0];
                let h = values
                    .get(1)
                    .cloned()
                    .unwrap_or_else(|| self.literal("0.5"));
                Ok(format!("(({x}) > {zero} ? {one} : (({x}) < {zero} ? {zero} : (({x}) == {zero} ? {h} : {})))",self.nan()))
            }
            Func::Re | Func::Conjugate => Ok(values[0].clone()),
            Func::Im => Ok(zero),
            Func::Sign => {
                let x = &values[0];
                Ok(format!("(({x}) > {zero} ? {one} : (({x}) < {zero} ? -{one} : (({x}) == {zero} ? {zero} : {})))",self.nan()))
            }
            Func::Arg => Ok(format!(
                "(({}) == {zero} ? {} : {})",
                values[0],
                self.nan(),
                self.call("atan2", &[&zero, &values[0]])
            )),
            Func::Logic(op) => self.logic(*op, &values),
            _ => Err(self.unsupported(f.name())),
        }
    }
    fn logic(&self, op: LogicOp, args: &[String]) -> Result {
        Ok(match op {
            LogicOp::And | LogicOp::Or | LogicOp::Xor => {
                if args.is_empty() {
                    if op == LogicOp::And {
                        "true".into()
                    } else {
                        "false".into()
                    }
                } else {
                    format!(
                        "({})",
                        args.join(match op {
                            LogicOp::And => " && ",
                            LogicOp::Or => " || ",
                            _ => " != ",
                        })
                    )
                }
            }
            LogicOp::Not => format!("(!({}))", args[0]),
            LogicOp::Implies => format!("((!({})) || ({}))", args[0], args[1]),
            LogicOp::ITE => format!("({} ? {} : {})", args[0], args[1], args[2]),
            LogicOp::Equivalent => {
                if args.len() < 2 {
                    "true".into()
                } else {
                    format!(
                        "({})",
                        args[1..]
                            .iter()
                            .map(|a| format!("({} == {a})", args[0]))
                            .collect::<Vec<_>>()
                            .join(" && ")
                    )
                }
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn real_literals_and_signed_powers() {
        let x = Expr::symbol("x");
        let half = Expr::num(Number::rational_i64(1, 2).unwrap());
        assert_eq!(slangcode(&half).unwrap(), "(1.0/2.0)");
        let settings = SlangSettings {
            dtype: SlangDtype::Float64,
            ..Default::default()
        };
        assert_eq!(slangcode_with(&half, &settings).unwrap(), "(1.0L/2.0L)");
        let p = crate::ops::pow(&x, &Expr::num(Number::int(-3)));
        let text = slangcode(&p).unwrap();
        assert!(text.contains("1.0/"));
        assert!(!text.contains("pow("));
        assert!(slangcode(&Expr::symbol("float3")).is_err());
        assert!(slangcode(&Expr::symbol("x;return 0")).is_err());
    }
    #[test]
    fn scalar_piecewise_and_boolean_types() {
        let x = Expr::symbol("x");
        let expression = crate::relational::piecewise(vec![(x.clone(), Expr::boolean(true))]);
        assert_eq!(slang_type(&expression).unwrap(), "real");
        assert_eq!(slangcode(&expression).unwrap(), "x");
        assert_eq!(slang_type(&Expr::boolean(false)).unwrap(), "bool");
    }
    #[test]
    fn domains_boolean_logic_and_explicit_failures() {
        use crate::parser::parse_str;
        let expression = parse_str("Piecewise((log(x), x > 0), (sqrt(-x), True))").unwrap();
        let code = slangcode(&expression).unwrap();
        assert!(code.contains(" ? log(x) : sqrt("));
        assert!(!code.contains("select("));
        let boolean = parse_str("And(x > 0, y < 1)").unwrap();
        assert_eq!(slang_type(&boolean).unwrap(), "bool");
        assert!(slangcode(&boolean).unwrap().contains(" && "));
        assert!(slangcode(&parse_str("gamma(x)").unwrap()).is_err());
        assert!(slangcode(&Expr::from_kind(Kind::Const(Constant::I))).is_err());
        let generic = slangcode(&parse_str("x**y").unwrap()).unwrap();
        assert!(generic.contains("floor(y)"));
        assert!(generic.contains("fmod(abs(y), 2.0)"));
    }
}
