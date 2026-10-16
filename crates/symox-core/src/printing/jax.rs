//! Native JAX expression printer, using `jnp` (`jax.numpy`), `jsp`
//! (`jax.scipy.special`) and `_sx_tensor(value, reference)`.
//!
//! The shared native traversal selects backend-specific operation mappings
//! before code construction. Piecewise branches receive masks before every
//! potentially singular operation, making `jit` and `grad` work together.

use super::{CodeError, TorchSettings};
use crate::expr::Expr;

pub type JaxSettings = TorchSettings;

pub struct JaxPrinter<'a>(super::pytorch::TorchPrinter<'a>);

impl<'a> JaxPrinter<'a> {
    pub fn new(settings: &'a JaxSettings) -> Self {
        Self(super::pytorch::TorchPrinter::jax(settings))
    }

    pub fn print(&self, expr: &Expr) -> Result<String, CodeError> {
        self.0.print(expr)
    }
}

pub fn jaxcode(expr: &Expr) -> Result<String, CodeError> {
    jaxcode_with(expr, &JaxSettings::default())
}

pub fn jaxcode_with(expr: &Expr, settings: &JaxSettings) -> Result<String, CodeError> {
    JaxPrinter::new(settings).print(expr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;
    fn code(source: &str) -> String {
        jaxcode(&parse_str(source).unwrap()).unwrap()
    }

    #[test]
    fn native_jax_mappings() {
        assert_eq!(code("sin(x)"), "jnp.sin(_sx_tensor(x, None))");
        assert_eq!(code("erf(x)"), "jsp.erf(_sx_tensor(x, None))");
        assert!(code("gamma(x)").starts_with("jsp.gamma("));
        assert!(code("x**y").starts_with("jnp.power("));
        assert!(code("polygamma(2,x)").starts_with("jsp.polygamma(2,"));
        assert!(jaxcode(&parse_str("besselj(0,x)").unwrap()).is_err());
        assert!(jaxcode(&parse_str("loggamma(x)").unwrap()).is_err());
    }

    #[test]
    fn jax_piecewise_uses_masked_domains() {
        let code = code("Piecewise((log(x), x > 0), (sqrt(-x), True))");
        assert!(code.starts_with("jnp.where("));
        assert!(code.contains("jnp.sqrt(_sx_tensor(jnp.where("));
        assert!(code.contains("jnp.asarray("));
        assert!(code.contains("jnp.logical_not("));
        assert!(!code.contains("torch"));
        assert!(!code.contains(" if "));
    }
}
