//! Printers and code generation.
//!
//! * [`latex`] — SymPy compatible LaTeX output,
//! * [`pretty`] — SymPy style 2D unicode pretty printer,
//! * [`cse`] — common subexpression elimination,
//! * [`code`] — C99/C++/CUDA/Python/NumPy/Rust code printers,
//! * [`codegen`] — complete C functions and CUDA kernels.
//!
//! The plain `str`/`srepr` printers live in [`crate::printer`].

pub mod code;
pub mod codegen;
pub mod cse;
pub mod latex;
pub mod jax;
pub mod order;
pub mod pretty;
pub mod pytorch;
pub mod slang;
pub mod sycl;

pub use code::{
    ccode, code_with, cudacode, cudacode_with, cxxcode, fcode, numpycode, pycode, rust_code,
    CodeError, CodePrinter, CodeSettings, CudaSettings, FloatType, Lang,
};
pub use codegen::{ccode_function, cuda_device_function, cuda_kernel, CodegenOptions};
pub use cse::{cse, cse_with, CseOptions};
pub use latex::{latex, latex_with, InvTrigStyle, LatexMode, LatexSettings};
pub use jax::{jaxcode, jaxcode_with, JaxPrinter, JaxSettings};
pub use pretty::{pprint, pretty, pretty_with, PrettySettings};
pub use pytorch::{torchcode, torchcode_with, TorchPrinter, TorchSettings};
pub use slang::{slangcode, slangcode_with, slang_type, SlangDtype, SlangPrinter, SlangSettings};
pub use sycl::{syclcode, syclcode_with, sycl_type, SyclDtype, SyclSettings};

use crate::expr::{Constant, Expr, Kind};
use crate::number::Number;
use crate::ops::mul_factor;

// ----------------------------------------------------------------------
// precedence (sympy.printing.precedence)
// ----------------------------------------------------------------------

pub(crate) const PREC_REL: i32 = 35;
pub(crate) const PREC_ADD: i32 = 40;
pub(crate) const PREC_MUL: i32 = 50;
pub(crate) const PREC_POW: i32 = 60;
pub(crate) const PREC_FUNC: i32 = 70;
pub(crate) const PREC_ATOM: i32 = 1000;

/// `could_extract_minus_sign` for printing purposes.
pub(crate) fn minus_sign(e: &Expr) -> bool {
    crate::ops::could_extract_minus_sign_syntactic(e)
}

/// SymPy's `precedence()` (the one used by the string and code printers).
pub(crate) fn precedence(e: &Expr) -> i32 {
    match e.kind() {
        Kind::Num(n) => match n {
            _ if n.is_negative() => PREC_ADD,
            Number::Rational(_) => PREC_MUL,
            _ => PREC_ATOM,
        },
        Kind::Const(Constant::NegInfinity) => PREC_ADD,
        Kind::Add(..) => PREC_ADD,
        Kind::Mul(c, _) => {
            if c.is_negative() {
                PREC_ADD
            } else {
                PREC_MUL
            }
        }
        Kind::Pow(..) => PREC_POW,
        Kind::Func(crate::expr::Func::Special(crate::special::SpecialFn::Mod), _) => PREC_MUL,
        // A root object is an Expr, not a Function: atom precedence.
        Kind::Func(crate::expr::Func::Undefined(head), args) if crate::poly::rootof::is_rootof_head(head, args.len()) => PREC_ATOM,
        Kind::Func(..) | Kind::Piecewise(_) => PREC_FUNC,
        Kind::Rel(op, ..) => match op {
            crate::expr::RelOp::Eq | crate::expr::RelOp::Ne => PREC_MUL,
            _ => PREC_REL,
        },
        _ => PREC_ATOM,
    }
}

/// SymPy's `precedence_traditional()` (LaTeX and pretty printers).
pub(crate) fn precedence_traditional(e: &Expr) -> i32 {
    match e.kind() {
        Kind::Derivative(..) => PREC_MUL,
        _ => precedence(e),
    }
}

// ----------------------------------------------------------------------
// factor helpers
// ----------------------------------------------------------------------

/// `Mul.as_ordered_factors()`: numeric coefficient first (when not 1), then
/// the factors in printing order. A non-`Mul` gives itself.
pub(crate) fn ordered_factors(e: &Expr) -> Vec<Expr> {
    match e.kind() {
        Kind::Mul(c, facs) => {
            let mut v: Vec<Expr> = facs.iter().map(|(b, x)| mul_factor(b, x)).collect();
            order::sort_commutative_factors(&mut v);
            if !c.is_one() {
                v.insert(0, Expr::num(c.clone()));
            }
            v
        }
        _ => vec![e.clone()],
    }
}

/// Free of symbols and undefined functions. The symbol Bloom mask of the
/// node makes the (common) negative answer O(1) instead of a traversal.
pub(crate) fn is_numeric(e: &Expr) -> bool {
    e.mask() == 0 && crate::functions::is_numeric_value(e)
}

/// SymPy's `is_negative` restricted to what can be decided without
/// assumptions: numbers, and real numeric expressions (decided by
/// evaluating them).
pub(crate) fn is_known_negative(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(n) => n.is_negative(),
        Kind::Const(Constant::NegInfinity) => true,
        Kind::Const(_) | Kind::Sym(_) => false,
        _ => {
            is_numeric(e)
                && !crate::functions::contains_i(e)
                && matches!(crate::evalf::eval_f64(e), Some(v) if v < 0.0)
        }
    }
}

// ----------------------------------------------------------------------
// floats (mpmath's to_str)
// ----------------------------------------------------------------------

/// Format `x` with `dps` significant digits like `mpmath.libmp.to_str`
/// (fixed notation when `-5 < exponent < dps`, otherwise `d.ddde+N`).
pub(crate) fn mpf_to_str(x: f64, dps: usize, strip_zeros: bool) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "+inf".into() } else { "-inf".into() };
    }
    if x == 0.0 {
        return if strip_zeros { "0.0".into() } else { "0.0".into() };
    }
    let dps = dps.max(1);
    let s = format!("{:.*e}", dps - 1, x.abs());
    let (mant, exp) = s.split_once('e').unwrap();
    let exp: i64 = exp.parse().unwrap();
    let digits: String = mant.chars().filter(|c| c.is_ascii_digit()).collect();
    let sign = if x < 0.0 { "-" } else { "" };
    let min_fixed = (-(dps as i64 / 3)).min(-5);
    let max_fixed = dps as i64;
    if min_fixed < exp && exp < max_fixed {
        let (int_part, frac_part) = if exp >= 0 {
            let k = (exp + 1) as usize;
            if k >= digits.len() {
                let mut ip = digits.clone();
                ip.extend(std::iter::repeat('0').take(k - digits.len()));
                (ip, String::new())
            } else {
                (digits[..k].to_string(), digits[k..].to_string())
            }
        } else {
            let mut f = "0".repeat((-exp - 1) as usize);
            f.push_str(&digits);
            ("0".to_string(), f)
        };
        let mut frac = frac_part;
        if strip_zeros {
            while frac.ends_with('0') {
                frac.pop();
            }
            if frac.is_empty() {
                frac.push('0');
            }
        }
        format!("{sign}{int_part}.{frac}")
    } else {
        let mut frac = digits[1..].to_string();
        if strip_zeros {
            while frac.ends_with('0') {
                frac.pop();
            }
            if frac.is_empty() {
                frac.push('0');
            }
        }
        let esign = if exp >= 0 { "+" } else { "-" };
        format!("{sign}{}.{frac}e{esign}{}", &digits[..1], exp.abs())
    }
}

/// SymPy's `str(Float)` inside an expression (15 digits, zeros stripped).
pub(crate) fn float_str(x: f64) -> String {
    let s = mpf_to_str(x, 15, true);
    s.strip_prefix('+').map(|t| t.to_string()).unwrap_or(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mpf_formatting() {
        assert_eq!(mpf_to_str(1e15, 15, true), "1.0e+15");
        assert_eq!(mpf_to_str(1e14, 15, true), "100000000000000.0");
        assert_eq!(mpf_to_str(1e14, 15, false), "100000000000000.");
        assert_eq!(mpf_to_str(1e-5, 15, true), "1.0e-5");
        assert_eq!(mpf_to_str(1e-4, 15, true), "0.0001");
        assert_eq!(mpf_to_str(0.1, 17, false), "0.10000000000000001");
        assert_eq!(mpf_to_str(1e-5, 17, false), "1.0000000000000001e-5");
        assert_eq!(mpf_to_str(-3.25, 15, true), "-3.25");
        assert_eq!(mpf_to_str(1e16, 17, false), "10000000000000000.");
        assert_eq!(mpf_to_str(1.0 / 3.0, 15, true), "0.333333333333333");
    }
}
