//! Native PyTorch expression printer.
//!
//! Generated expressions use `torch` and `_sx_tensor(value, reference)`.
//! The latter preserves tensor inputs and creates scalar tensors on the
//! reference tensor's device and with its dtype. Complete function generation
//! supplies that helper. Piecewise masks are propagated into operations before
//! evaluation, keeping inactive branches finite for reverse-mode autograd.

use std::collections::BTreeMap;

use crate::expr::{Constant, Expr, Func, Kind};
use crate::logic::LogicOp;
use crate::number::Number;
use crate::special::SpecialFn;

use super::CodeError;

type Result = std::result::Result<String, CodeError>;

#[derive(Clone, Debug, Default)]
pub struct TorchSettings {
    pub symbol_names: BTreeMap<String, String>,
    /// Exact symbol identities, including assumptions, take precedence.
    pub exact_symbol_names: BTreeMap<Expr, String>,
    pub user_functions: BTreeMap<String, String>,
    /// Explicit printer hooks; never reparsed as symbolic expressions.
    pub expression_overrides: BTreeMap<Expr, String>,
    pub reference: Option<String>,
    pub mask: Option<String>,
}

/// A Python identifier that is safe to place in generated source.
pub fn is_python_identifier(name: &str) -> bool {
    let mut chars = name.chars();
    matches!(chars.next(), Some(c) if c == '_' || c.is_ascii_alphabetic())
        && chars.all(|c| c == '_' || c.is_ascii_alphanumeric())
        && !matches!(
            name,
            "False"
                | "None"
                | "True"
                | "and"
                | "as"
                | "assert"
                | "async"
                | "await"
                | "break"
                | "class"
                | "continue"
                | "def"
                | "del"
                | "elif"
                | "else"
                | "except"
                | "finally"
                | "for"
                | "from"
                | "global"
                | "if"
                | "import"
                | "in"
                | "is"
                | "lambda"
                | "nonlocal"
                | "not"
                | "or"
                | "pass"
                | "raise"
                | "return"
                | "try"
                | "while"
                | "with"
                | "yield"
        )
}

pub fn torchcode(expr: &Expr) -> Result {
    torchcode_with(expr, &TorchSettings::default())
}

pub fn torchcode_with(expr: &Expr, settings: &TorchSettings) -> Result {
    TorchPrinter::new(settings).print(expr)
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ArrayBackend {
    Torch,
    Jax,
}

pub struct TorchPrinter<'a> {
    settings: &'a TorchSettings,
    backend: ArrayBackend,
}

impl<'a> TorchPrinter<'a> {
    pub fn new(settings: &'a TorchSettings) -> Self {
        Self {
            settings,
            backend: ArrayBackend::Torch,
        }
    }

    pub(crate) fn jax(settings: &'a TorchSettings) -> Self {
        Self {
            settings,
            backend: ArrayBackend::Jax,
        }
    }

    fn module(&self) -> &'static str {
        match self.backend {
            ArrayBackend::Torch => "torch",
            ArrayBackend::Jax => "jnp",
        }
    }

    fn special_module(&self) -> &'static str {
        match self.backend {
            ArrayBackend::Torch => "torch.special",
            ArrayBackend::Jax => "jsp",
        }
    }

    pub fn print(&self, expr: &Expr) -> Result {
        self.expr(expr, self.settings.mask.as_deref())
    }

    fn unsupported(&self, name: &str) -> CodeError {
        CodeError::Unsupported(format!(
            "{name} in {} code generation",
            if self.backend == ArrayBackend::Torch {
                "PyTorch"
            } else {
                "JAX"
            }
        ))
    }

    fn tensor(&self, code: &str) -> String {
        let reference = self.settings.reference.as_deref().unwrap_or("None");
        format!("_sx_tensor({code}, {reference})")
    }

    fn boolean(&self, code: &str) -> String {
        match self.backend {
            ArrayBackend::Torch => format!("{}.to(dtype=torch.bool)", self.tensor(code)),
            ArrayBackend::Jax => format!("jnp.asarray({}, dtype=jnp.bool_)", self.tensor(code)),
        }
    }

    fn safe(&self, code: String, mask: Option<&str>, neutral: &str) -> String {
        let module = self.module();
        match mask {
            Some(mask) => format!(
                "{module}.where({}, {}, {})",
                self.boolean(mask),
                self.tensor(&code),
                self.tensor(neutral)
            ),
            None => code,
        }
    }

    fn conjunction(&self, lhs: &str, rhs: &str) -> String {
        let module = self.module();
        format!(
            "{module}.logical_and({}, {})",
            self.boolean(lhs),
            self.boolean(rhs)
        )
    }

    fn negate(&self, value: &str) -> String {
        let module = self.module();
        format!("{module}.logical_not({})", self.boolean(value))
    }

    fn number(&self, number: &Number) -> String {
        if let Some(code) = self
            .settings
            .expression_overrides
            .get(&Expr::num(number.clone()))
        {
            return format!("({code})");
        }
        match number {
            Number::Rational(value) => format!("({}/{})", value.numer(), value.denom()),
            Number::Float(value) if value.is_nan() => "float('nan')".into(),
            Number::Float(value) if value.is_infinite() => {
                if *value > 0.0 {
                    "float('inf')".into()
                } else {
                    "-float('inf')".into()
                }
            }
            Number::Float(value) => format!("{value:?}"),
            _ => number.to_string(),
        }
    }

    fn expr(&self, expr: &Expr, mask: Option<&str>) -> Result {
        if let Some(code) = self.settings.expression_overrides.get(expr) {
            return Ok(self.safe(format!("({code})"), mask, "0"));
        }
        let module = self.module();
        match expr.kind() {
            Kind::Num(number) => Ok(self.number(number)),
            Kind::Bool(value) => Ok(if *value { "True" } else { "False" }.into()),
            Kind::Const(constant) => Ok(match constant {
                Constant::Pi => format!("{module}.pi"),
                Constant::E => "2.718281828459045".into(),
                Constant::I => "1j".into(),
                Constant::EulerGamma => "0.5772156649015329".into(),
                Constant::Infinity => "float('inf')".into(),
                Constant::NegInfinity => "-float('inf')".into(),
                Constant::NaN => "float('nan')".into(),
                Constant::ComplexInfinity => return Err(self.unsupported("complex infinity")),
            }),
            Kind::Sym(symbol) => {
                let name = self
                    .settings
                    .exact_symbol_names
                    .get(expr)
                    .or_else(|| self.settings.symbol_names.get(symbol.name()))
                    .map(String::as_str)
                    .unwrap_or(symbol.name());
                if !is_python_identifier(name)
                    || matches!(name, "torch" | "jax" | "jnp" | "jsp" | "_sx_tensor")
                {
                    return Err(self.unsupported(&format!(
                        "unsafe symbol name {name:?}; supply symbol_names"
                    )));
                }
                // Masking leaves prevents overflow in polynomial intermediates
                // even when all nonlinear operators are independently guarded.
                Ok(self.safe(name.to_string(), mask, "0"))
            }
            Kind::Add(coefficient, terms) => {
                let mut parts = Vec::new();
                if !coefficient.is_zero() {
                    parts.push(self.number(coefficient));
                }
                for (term, coefficient) in terms {
                    let term = self.expr(term, mask)?;
                    parts.push(if coefficient.is_one() {
                        term
                    } else {
                        format!("({}*{term})", self.number(coefficient))
                    });
                }
                Ok(if parts.is_empty() {
                    "0".into()
                } else {
                    format!("({})", parts.join(" + "))
                })
            }
            Kind::Mul(coefficient, factors) => {
                let mut parts = Vec::new();
                if !coefficient.is_one() {
                    parts.push(self.number(coefficient));
                }
                for (base, exponent) in factors {
                    parts.push(if exponent.is_one() {
                        self.expr(base, mask)?
                    } else {
                        self.power(base, exponent, mask)?
                    });
                }
                Ok(if parts.is_empty() {
                    "1".into()
                } else {
                    format!("({})", parts.join("*"))
                })
            }
            Kind::Pow(base, exponent) => self.power(base, exponent, mask),
            Kind::Rel(op, lhs, rhs) => Ok(format!(
                "({} {} {})",
                self.expr(lhs, mask)?,
                op.symbol(),
                self.expr(rhs, mask)?
            )),
            Kind::Piecewise(pairs) => self.piecewise(pairs, mask),
            Kind::Func(Func::MatrixElement(parent), args) => {
                let crate::matrixexpr::MatrixKind::Symbol(name) = parent.kind() else {
                    return Err(self.unsupported("an unevaluated matrix operation in a scalar entry"));
                };
                let proxy = Expr::from_kind(Kind::Func(Func::MatrixProxy(parent.clone()), vec![]));
                let name = self.settings.exact_symbol_names.get(&proxy).or_else(|| self.settings.symbol_names.get(name)).map(String::as_str).unwrap_or(name);
                if !is_python_identifier(name) { return Err(self.unsupported("unsafe matrix symbol name")); }
                Ok(format!("{name}[{}, {}]", self.expr(&args[0], None)?, self.expr(&args[1], None)?))
            }
            Kind::Func(Func::MatrixTrace(_) | Func::MatrixDeterminant(_), _) => {
                let explicit = crate::functions::materialize_matrix_scalars(expr)
                    .map_err(|error| self.unsupported(&error.to_string()))?;
                self.expr(&explicit, mask)
            }
            Kind::Func(Func::MatrixProxy(_), _) => Err(self.unsupported("a matrix expression used as a scalar")),
            Kind::Func(function, args) => self.function(function, args, mask),
            Kind::Derivative(..) => Err(self.unsupported("unevaluated Derivative")),
        }
    }

    fn power(&self, base: &Expr, exponent: &Expr, mask: Option<&str>) -> Result {
        let module = self.module();
        let base_code = self.safe(self.expr(base, mask)?, mask, "1");
        let exponent_code = self.expr(exponent, mask)?;
        if exponent
            .as_number()
            .is_some_and(|n| n == &Number::rational_i64(1, 2).unwrap())
        {
            return Ok(format!("{module}.sqrt({})", self.tensor(&base_code)));
        }
        if exponent.is_minus_one() {
            return Ok(format!("(1/({base_code}))"));
        }
        if exponent.is_integer_number() {
            return Ok(format!("(({base_code})**({exponent_code}))"));
        }
        // The exponent is masked as well: enormous inactive exponents can
        // otherwise yield 0*inf in backward, even at the neutral base 1.
        let exponent_code = self.safe(exponent_code, mask, "1");
        Ok(format!(
            "{module}.{}({}, {})",
            if self.backend == ArrayBackend::Torch {
                "pow"
            } else {
                "power"
            },
            self.tensor(&base_code),
            self.tensor(&exponent_code)
        ))
    }

    fn piecewise(&self, pairs: &[(Expr, Expr)], mask: Option<&str>) -> Result {
        let module = self.module();
        let mut remaining = mask.map(str::to_string);
        let mut branches = Vec::new();
        let mut fallback = self.tensor("float('nan')");
        for (value, condition) in pairs {
            if matches!(condition.kind(), Kind::Bool(true)) {
                fallback = self.tensor(&self.expr(value, remaining.as_deref())?);
                break;
            }
            let condition = self.expr(condition, remaining.as_deref())?;
            let active = match remaining.as_deref() {
                Some(parent) => self.conjunction(parent, &condition),
                None => condition.clone(),
            };
            branches.push((
                condition.clone(),
                self.tensor(&self.expr(value, Some(&active))?),
            ));
            let inverse = self.negate(&condition);
            remaining = Some(match remaining.as_deref() {
                Some(parent) => self.conjunction(parent, &inverse),
                None => inverse,
            });
        }
        for (condition, value) in branches.into_iter().rev() {
            fallback = format!(
                "{module}.where({}, {value}, {fallback})",
                self.boolean(&condition)
            );
        }
        Ok(fallback)
    }

    fn unary(&self, name: &str, arg: &Expr, mask: Option<&str>, neutral: &str) -> Result {
        let arg = self.safe(self.expr(arg, mask)?, mask, neutral);
        let target = if let Some(name) = name.strip_prefix("special.") {
            format!("{}.{name}", self.special_module())
        } else if self.backend == ArrayBackend::Jax {
            match name {
                "erf" | "erfc" => format!("jsp.{name}"),
                "lgamma" => "jsp.gammaln".into(),
                "sgn" => "jnp.sign".into(),
                "asin" => "jnp.arcsin".into(),
                "acos" => "jnp.arccos".into(),
                "atan" => "jnp.arctan".into(),
                "asinh" => "jnp.arcsinh".into(),
                "acosh" => "jnp.arccosh".into(),
                "atanh" => "jnp.arctanh".into(),
                _ => format!("jnp.{name}"),
            }
        } else {
            format!("torch.{name}")
        };
        Ok(format!("{target}({})", self.tensor(&arg)))
    }

    fn function(&self, function: &Func, args: &[Expr], mask: Option<&str>) -> Result {
        let module = self.module();
        if let Some(target) = self.settings.user_functions.get(function.name()) {
            if !target.split('.').all(is_python_identifier) {
                return Err(self.unsupported("invalid user function target"));
            }
            let args = args
                .iter()
                .map(|arg| self.expr(arg, mask))
                .collect::<std::result::Result<Vec<_>, _>>()?;
            return Ok(format!("{target}({})", args.join(", ")));
        }
        if *function == Func::LogGamma && crate::assumptions::is_positive(&args[0]) != Some(true) {
            return Err(self.unsupported("loggamma requires a provably positive real argument"));
        }
        let unary = match function {
            Func::Sin => Some(("sin", "0")),
            Func::Cos => Some(("cos", "0")),
            Func::Tan => Some(("tan", "0")),
            Func::Asin => Some(("asin", "0")),
            Func::Acos => Some(("acos", "0")),
            Func::Atan => Some(("atan", "0")),
            Func::Sinh => Some(("sinh", "0")),
            Func::Cosh => Some(("cosh", "0")),
            Func::Tanh => Some(("tanh", "0")),
            Func::Asinh => Some(("asinh", "0")),
            Func::Acosh => Some(("acosh", "2")),
            Func::Atanh => Some(("atanh", "0")),
            Func::Exp => Some(("exp", "0")),
            Func::Log => Some(("log", "1")),
            Func::Abs => Some(("abs", "1")),
            Func::Sign => Some(("sgn", "1")),
            Func::Floor => Some(("floor", "0")),
            Func::Ceiling => Some(("ceil", "0")),
            Func::Erf => Some(("erf", "0")),
            Func::Erfc => Some(("erfc", "0")),
            Func::LogGamma => Some(("lgamma", "1")),
            Func::Re => Some(("real", "0")),
            Func::Conjugate => Some(("conj", "0")),
            Func::Arg => Some(("angle", "1")),
            Func::Special(SpecialFn::Digamma) => Some(("special.digamma", "1")),
            Func::Special(SpecialFn::Erfinv) => Some(("special.erfinv", "0")),
            _ => None,
        };
        if let Some((name, neutral)) = unary {
            return self.unary(name, &args[0], mask, neutral);
        }
        match function {
            Func::Undefined(name) if name.as_ref() == "UnevaluatedExpr" && args.len() == 1 => {
                self.expr(&args[0], mask)
            }
            Func::Im => {
                let arg = self.expr(&args[0], mask)?;
                // Adding a complex scalar preserves float32/64 precision and
                // handles both real and complex inputs without dtype branches.
                Ok(format!("{module}.imag({} + 0j)", self.tensor(&arg)))
            }
            Func::Atan2 => {
                let y = self.safe(self.expr(&args[0], mask)?, mask, "0");
                let x = self.safe(self.expr(&args[1], mask)?, mask, "1");
                Ok(format!(
                    "{module}.{}({}, {})",
                    if self.backend == ArrayBackend::Torch {
                        "atan2"
                    } else {
                        "arctan2"
                    },
                    self.tensor(&y),
                    self.tensor(&x)
                ))
            }
            Func::Max | Func::Min => {
                let target = if *function == Func::Max {
                    "maximum"
                } else {
                    "minimum"
                };
                let mut values = args
                    .iter()
                    .map(|arg| self.expr(arg, mask).map(|s| self.tensor(&s)));
                let mut result = values
                    .next()
                    .ok_or_else(|| self.unsupported("empty Min/Max"))??;
                for value in values {
                    result = format!("{module}.{target}({result}, {})", value?);
                }
                Ok(result)
            }
            Func::Heaviside => {
                let x = self.expr(&args[0], mask)?;
                let middle = if args.len() > 1 {
                    self.expr(&args[1], mask)?
                } else {
                    "(1/2)".into()
                };
                let zero = self.tensor("0");
                let one = self.tensor("1");
                Ok(format!("{module}.where({}, {one}, {module}.where({}, {zero}, {module}.where({}, {}, {})))",
                    self.boolean(&format!("({x} > 0)")), self.boolean(&format!("({x} < 0)")),
                    self.boolean(&format!("({x} == 0)")), self.tensor(&middle), self.tensor("float('nan')")))
            }
            Func::Cot | Func::Sec | Func::Csc | Func::Coth => {
                let (name, neutral) = match function {
                    Func::Cot => ("tan", "1"),
                    Func::Sec => ("cos", "0"),
                    Func::Csc => ("sin", "1"),
                    _ => ("tanh", "1"),
                };
                Ok(format!(
                    "(1/({}))",
                    self.unary(name, &args[0], mask, neutral)?
                ))
            }
            Func::Acot => {
                // SymPy's principal acot is atan(1/x), including pi/2 at 0.
                let x = self.safe(self.expr(&args[0], mask)?, mask, "1");
                let nonzero = format!("({x} != 0)");
                let safe_x = self.safe(x.clone(), Some(&nonzero), "1");
                Ok(format!(
                    "{module}.where({}, {module}.{}(1/({})), {})",
                    self.boolean(&nonzero),
                    if self.backend == ArrayBackend::Torch {
                        "atan"
                    } else {
                        "arctan"
                    },
                    self.tensor(&safe_x),
                    self.tensor(&format!("({module}.pi/2)"))
                ))
            }
            Func::Gamma if self.backend == ArrayBackend::Jax => {
                self.unary("special.gamma", &args[0], mask, "1")
            }
            Func::Gamma => {
                let x = self.safe(self.expr(&args[0], mask)?, mask, "1");
                let x = self.tensor(&x);
                Ok(format!("({module}.exp({module}.lgamma({x}))*{module}.where(({x} < 0), {module}.sign({module}.sin({module}.pi*{x})), {module}.ones_like({x})))"))
            }
            Func::Special(special) => self.special(*special, args, mask),
            Func::Logic(op) => self.logic(*op, args, mask),
            _ => Err(self.unsupported(function.name())),
        }
    }

    fn special(&self, function: SpecialFn, args: &[Expr], mask: Option<&str>) -> Result {
        let module = self.special_module();
        match function {
            SpecialFn::Polygamma => {
                let n = args[0]
                    .as_number()
                    .and_then(Number::to_i64)
                    .filter(|n| *n >= 0)
                    .ok_or_else(|| {
                        self.unsupported("polygamma with nonconstant/noninteger order")
                    })?;
                let x = self.safe(self.expr(&args[1], mask)?, mask, "1");
                Ok(format!("{module}.polygamma({n}, {})", self.tensor(&x)))
            }
            SpecialFn::Trigamma => {
                let x = self.safe(self.expr(&args[0], mask)?, mask, "1");
                Ok(format!("{module}.polygamma(1, {})", self.tensor(&x)))
            }
            SpecialFn::Zeta => {
                let x = self.safe(self.expr(&args[0], mask)?, mask, "2");
                let q = if args.len() == 2 {
                    self.safe(self.expr(&args[1], mask)?, mask, "1")
                } else {
                    "1".into()
                };
                Ok(format!(
                    "{module}.zeta({}, {})",
                    self.tensor(&x),
                    self.tensor(&q)
                ))
            }
            SpecialFn::Erfcinv => {
                let x = self.safe(self.expr(&args[0], mask)?, mask, "1");
                Ok(format!("{module}.erfinv(1 - {})", self.tensor(&x)))
            }
            SpecialFn::BesselJ | SpecialFn::BesselY | SpecialFn::BesselI | SpecialFn::BesselK => {
                let n = args[0].as_number().and_then(Number::to_i64);
                if self.backend == ArrayBackend::Jax && function != SpecialFn::BesselI {
                    return Err(self.unsupported(function.name()));
                }
                let target = match (function, n) {
                    (SpecialFn::BesselJ, Some(0)) => "bessel_j0",
                    (SpecialFn::BesselJ, Some(1)) => "bessel_j1",
                    (SpecialFn::BesselY, Some(0)) => "bessel_y0",
                    (SpecialFn::BesselY, Some(1)) => "bessel_y1",
                    (SpecialFn::BesselI, Some(0)) => "i0",
                    (SpecialFn::BesselI, Some(1)) => "i1",
                    (SpecialFn::BesselK, Some(0)) => "modified_bessel_k0",
                    (SpecialFn::BesselK, Some(1)) => "modified_bessel_k1",
                    _ => {
                        return Err(self.unsupported(&format!(
                            "{} with order other than 0 or 1",
                            function.name()
                        )))
                    }
                };
                self.unary(&format!("special.{target}"), &args[1], mask, "1")
            }
            _ => Err(self.unsupported(function.name())),
        }
    }

    fn logic(&self, op: LogicOp, args: &[Expr], mask: Option<&str>) -> Result {
        let module = self.module();
        let values = args
            .iter()
            .map(|arg| self.expr(arg, mask))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let fold = |target: &str, values: &[String], identity: &str| {
            let mut result = self.boolean(values.first().map(String::as_str).unwrap_or(identity));
            for value in values.iter().skip(1) {
                result = format!("{module}.{target}({result}, {})", self.boolean(value));
            }
            result
        };
        Ok(match op {
            LogicOp::And => fold("logical_and", &values, "True"),
            LogicOp::Or => fold("logical_or", &values, "False"),
            LogicOp::Xor => fold("logical_xor", &values, "False"),
            LogicOp::Not => self.negate(&values[0]),
            LogicOp::Implies => format!(
                "{module}.logical_or({}, {})",
                self.negate(&values[0]),
                self.boolean(&values[1])
            ),
            LogicOp::Equivalent => {
                let pairs = values
                    .windows(2)
                    .map(|v| {
                        format!(
                            "{module}.{}({}, {})",
                            if self.backend == ArrayBackend::Torch {
                                "eq"
                            } else {
                                "equal"
                            },
                            self.boolean(&v[0]),
                            self.boolean(&v[1])
                        )
                    })
                    .collect::<Vec<_>>();
                fold("logical_and", &pairs, "True")
            }
            LogicOp::ITE => format!(
                "{module}.where({}, {}, {})",
                self.boolean(&values[0]),
                self.boolean(&values[1]),
                self.boolean(&values[2])
            ),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    fn code(source: &str) -> String {
        torchcode(&parse_str(source).unwrap()).unwrap()
    }

    #[test]
    fn arithmetic_and_tensors() {
        assert_eq!(
            code("x/3 + sin(y)"),
            "(((1/3)*x) + torch.sin(_sx_tensor(y, None)))"
        );
        assert_eq!(code("sqrt(x)"), "torch.sqrt(_sx_tensor(x, None))");
        assert_eq!(code("x**-1"), "(1/(x))");
        assert!(code("Max(x, y, 3)").contains("torch.maximum"));
    }

    #[test]
    fn piecewise_guards_domains_before_evaluation() {
        let output = code("Piecewise((log(x), x > 0), (sqrt(-x), True))");
        assert!(output.starts_with("torch.where("));
        assert!(output.contains("torch.log(_sx_tensor(torch.where("));
        assert!(output.contains("torch.sqrt(_sx_tensor(torch.where("));
        assert!(output.contains("torch.logical_not"));
        assert!(!output.contains(" if "));
        let output = code("Piecewise((1/x, x > 0), (exp(x), True))");
        assert!(output.contains("torch.where"));
        assert!(output.contains("torch.exp(_sx_tensor(torch.where("));
    }

    #[test]
    fn exact_printer_overrides_and_empty_arithmetic() {
        let expression = parse_str("sin(x) + x**2").unwrap();
        let mut settings = TorchSettings::default();
        settings
            .expression_overrides
            .insert(Expr::symbol("x"), "inputs[0]".into());
        let code = torchcode_with(&expression, &settings).unwrap();
        assert!(code.contains("torch.sin(_sx_tensor((inputs[0]), None))"));
        assert!(code.contains("(((inputs[0]))**(2))"));
        let empty_mul = Expr::from_unevaluated_kind(Kind::Mul(Number::ONE, vec![]));
        assert_eq!(torchcode(&empty_mul).unwrap(), "1");
        let empty_add = Expr::from_unevaluated_kind(Kind::Add(Number::ZERO, vec![]));
        assert_eq!(torchcode(&empty_add).unwrap(), "0");
    }

    #[test]
    fn booleans_specials_and_unsupported() {
        assert!(code("And(x > 0, y < 1)").starts_with("torch.logical_and("));
        assert!(code("polygamma(2, x)").starts_with("torch.special.polygamma(2,"));
        assert!(code("zeta(x)").starts_with("torch.special.zeta("));
        assert!(torchcode(&parse_str("f(x)").unwrap()).is_err());
        assert!(torchcode(&Expr::symbol("class")).is_err());
        assert!(torchcode(&parse_str("besselj(3, x)").unwrap()).is_err());
    }
}
