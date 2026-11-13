//! SYCL 2020 scalar expressions for AdaptiveCpp and other SYCL compilers.
//!
//! Expressions use `sycl::` math and `<limits>` for nonfinite constants.
//! Their caller must include `<sycl/sycl.hpp>` and `<limits>`.

use super::slang::{SlangPrinter, SlangSettings};
use super::CodeError;
use crate::expr::{Expr, Kind};
use std::collections::BTreeMap;

pub use super::slang::SlangDtype as SyclDtype;

#[derive(Clone, Debug)]
pub struct SyclSettings {
    pub dtype: SyclDtype,
    pub symbol_names: BTreeMap<String, String>,
    pub exact_symbol_names: BTreeMap<Expr, String>,
    pub user_functions: BTreeMap<String, String>,
    pub expression_overrides: BTreeMap<Expr, String>,
}

impl Default for SyclSettings {
    fn default() -> Self {
        Self {
            dtype: SyclDtype::Float64,
            symbol_names: BTreeMap::new(),
            exact_symbol_names: BTreeMap::new(),
            user_functions: BTreeMap::new(),
            expression_overrides: BTreeMap::new(),
        }
    }
}

/// Conservative C++ identifier validation, excluding implementation names.
pub fn is_sycl_identifier(name: &str) -> bool {
    let mut chars = name.chars();
    if !matches!(chars.next(), Some(c) if c == '_' || c.is_ascii_alphabetic())
        || !chars.all(|c| c == '_' || c.is_ascii_alphanumeric())
        || name.contains("__")
        || (name.starts_with('_') && name.as_bytes().get(1).is_some_and(u8::is_ascii_uppercase))
    {
        return false;
    }
    ![
        "alignas",
        "alignof",
        "and",
        "and_eq",
        "asm",
        "atomic_cancel",
        "atomic_commit",
        "atomic_noexcept",
        "auto",
        "bitand",
        "bitor",
        "bool",
        "break",
        "case",
        "catch",
        "char",
        "char8_t",
        "char16_t",
        "char32_t",
        "class",
        "compl",
        "concept",
        "const",
        "consteval",
        "constexpr",
        "constinit",
        "const_cast",
        "continue",
        "co_await",
        "co_return",
        "co_yield",
        "decltype",
        "default",
        "delete",
        "do",
        "double",
        "dynamic_cast",
        "else",
        "enum",
        "explicit",
        "export",
        "extern",
        "false",
        "float",
        "for",
        "friend",
        "goto",
        "if",
        "inline",
        "int",
        "long",
        "mutable",
        "namespace",
        "new",
        "noexcept",
        "not",
        "not_eq",
        "nullptr",
        "operator",
        "or",
        "or_eq",
        "private",
        "protected",
        "public",
        "reflexpr",
        "register",
        "reinterpret_cast",
        "requires",
        "return",
        "short",
        "signed",
        "sizeof",
        "static",
        "static_assert",
        "static_cast",
        "struct",
        "switch",
        "synchronized",
        "template",
        "this",
        "thread_local",
        "throw",
        "true",
        "try",
        "typedef",
        "typeid",
        "typename",
        "union",
        "unsigned",
        "using",
        "virtual",
        "void",
        "volatile",
        "wchar_t",
        "while",
        "xor",
        "xor_eq",
        "std",
        "sycl",
        "acpp",
    ]
    .contains(&name)
}

pub fn sycl_type(expr: &Expr) -> std::result::Result<&'static str, CodeError> {
    super::slang::slang_type(expr)
}

pub fn syclcode(expr: &Expr) -> std::result::Result<String, CodeError> {
    syclcode_with(expr, &SyclSettings::default())
}

pub fn syclcode_with(
    expr: &Expr,
    settings: &SyclSettings,
) -> std::result::Result<String, CodeError> {
    let mut emitted = BTreeMap::new();
    for symbol in crate::traverse::free_symbols(expr) {
        if settings.expression_overrides.contains_key(&symbol) { continue; }
        let Some(name) = symbol.as_symbol() else {
            return Err(CodeError::Unsupported("matrix-valued SYCL inputs; materialize entries and substitute scalar inputs first".into()));
        };
        let default_name = match symbol.kind() {
            Kind::Sym(s) if s.is_dummy() => format!("{}_{index}", s.name(), index=s.dummy_index().unwrap()),
            _ => name.to_owned(),
        };
        let target = settings
            .exact_symbol_names
            .get(&symbol)
            .or_else(|| settings.symbol_names.get(name))
            .map(String::as_str)
            .unwrap_or(&default_name);
        if emitted
            .insert(target.to_string(), symbol.clone())
            .is_some_and(|previous| previous != symbol)
        {
            return Err(CodeError::Unsupported(format!(
                "ambiguous SYCL symbol {target:?}; supply exact symbol_names"
            )));
        }
    }
    let common = SlangSettings {
        dtype: settings.dtype,
        symbol_names: settings.symbol_names.clone(),
        exact_symbol_names: settings.exact_symbol_names.clone(),
        user_functions: settings.user_functions.clone(),
        expression_overrides: settings.expression_overrides.clone(),
    };
    SlangPrinter::sycl(&common).print(expr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    #[test]
    fn typed_math_and_real_domains() {
        let expr = parse_str("sin(x) + x/3 + sqrt(y)").unwrap();
        let code = syclcode(&expr).unwrap();
        assert!(code.contains("sycl::sin(x)"));
        assert!(code.contains("sycl::sqrt(y)"));
        assert!(code.contains("1.0/3.0"));
        assert!(!code.contains("1.0L"));
        let settings = SyclSettings {
            dtype: SyclDtype::Float32,
            ..Default::default()
        };
        let code = syclcode_with(&expr, &settings).unwrap();
        assert!(code.contains("1.0F/3.0F"));
        let half = SyclSettings {
            dtype: SyclDtype::Float16,
            ..Default::default()
        };
        let code = syclcode_with(&expr, &half).unwrap();
        assert!(code.contains("sycl::half(1.0F)/sycl::half(3.0F)"), "{code}");
        let code = syclcode_with(&parse_str("Piecewise((log(x), x > 0))").unwrap(), &half).unwrap();
        assert!(code.contains("std::numeric_limits<sycl::half>::quiet_NaN()"), "{code}");
        let generic = syclcode(&parse_str("x**y").unwrap()).unwrap();
        assert!(generic.contains("sycl::floor(y)"));
        assert!(generic.contains("sycl::fmod(sycl::fabs(y), 2.0)"));
        assert!(syclcode(&parse_str("gamma(x)").unwrap())
            .unwrap()
            .contains("sycl::tgamma"));
        assert!(syclcode(&parse_str("loggamma(x)").unwrap()).is_err());
    }
    #[test]
    fn limits_lazy_piecewise_and_identifiers() {
        let expr = parse_str("Piecewise((log(x), x > 0))").unwrap();
        let code = syclcode(&expr).unwrap();
        assert!(code.contains(" ? sycl::log(x) : std::numeric_limits<double>::quiet_NaN()"));
        assert!(!code.contains("select("));
        assert!(is_sycl_identifier("sin"));
        assert!(!is_sycl_identifier("class"));
        assert!(!is_sycl_identifier("name__reserved"));
        assert!(!is_sycl_identifier("_Reserved"));
    }
}
