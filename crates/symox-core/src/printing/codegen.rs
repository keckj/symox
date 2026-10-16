//! Complete C functions and CUDA kernels evaluating expressions,
//! optionally sharing common subexpressions.
//!
//! ```
//! use symox::prelude::*;
//! use symox::printing::{ccode_function, CodegenOptions};
//!
//! let (x, y) = (Expr::symbol("x"), Expr::symbol("y"));
//! let e = sin(&(&x + &y)) * cos(&(&x + &y));
//! let src = ccode_function("f", &[x, y], &[e], &CodegenOptions::scalar()).unwrap();
//! assert!(src.contains("const double x0 = x + y;"));
//! ```

use super::code::{CodeError, CodePrinter, CodeSettings, CudaSettings, FloatType, Lang};
use super::cse::{cse_with, CseOptions};
use crate::expr::Expr;

/// Options of [`ccode_function`] and [`cuda_kernel`].
#[derive(Clone, Debug)]
pub struct CodegenOptions {
    /// Extract common subexpressions into local constants.
    pub cse: bool,
    /// Floating point type.
    pub dtype: FloatType,
    /// C: loop over arrays (`n` elements) instead of a scalar function.
    pub vectorized: bool,
    /// CUDA: fast-math intrinsics (`__sinf`, ...) in single precision.
    pub fast_math: bool,
    /// CUDA: `rsqrt` for inverse square roots.
    pub use_rsqrt: bool,
    /// CUDA: write `x**n` (`2 <= n <= inline_pow`) as products.
    pub inline_pow: u32,
    /// Emit `#include <math.h>` (C only).
    pub include_header: bool,
}

impl Default for CodegenOptions {
    fn default() -> Self {
        CodegenOptions {
            cse: true,
            dtype: FloatType::F64,
            vectorized: true,
            fast_math: false,
            use_rsqrt: true,
            inline_pow: 4,
            include_header: true,
        }
    }
}

impl CodegenOptions {
    /// Scalar C function (`double f(double x, ...)`).
    pub fn scalar() -> Self {
        CodegenOptions { vectorized: false, ..Default::default() }
    }
}

struct Body {
    /// `(name, value)` local constants in evaluation order.
    temps: Vec<(String, String)>,
    outputs: Vec<String>,
    declarations: Vec<String>,
}

fn arg_names(p: &CodePrinter, args: &[Expr]) -> Result<Vec<String>, CodeError> {
    args.iter()
        .map(|a| match a.as_symbol() {
            Some(s) => Ok(p.symbol_name(s)),
            None => Err(CodeError::Unsupported(format!("argument {a} is not a symbol"))),
        })
        .collect()
}

fn body(p: &CodePrinter, args: &[Expr], exprs: &[Expr], use_cse: bool) -> Result<Body, CodeError> {
    let _session = super::order::session();
    let (reps, reduced) = if use_cse {
        let opts = CseOptions {
            excluded_names: args.iter().filter_map(|a| a.as_symbol().map(|s| s.to_string())).collect(),
            ..Default::default()
        };
        cse_with(exprs, &opts)
    } else {
        (vec![], exprs.to_vec())
    };
    let mut temps = Vec::with_capacity(reps.len());
    for (s, v) in &reps {
        temps.push((p.print(s)?, p.print(v)?));
    }
    let outputs = reduced.iter().map(|e| p.print(e)).collect::<Result<Vec<_>, _>>()?;
    Ok(Body { temps, outputs, declarations: p.declarations() })
}

/// Indent every line of a (possibly multi-line) statement.
fn indent(s: &str, pad: &str) -> String {
    s.lines().map(|l| format!("{pad}{l}")).collect::<Vec<_>>().join("\n")
}

/// A complete C99 function evaluating `exprs` (one output each).
///
/// * vectorized (default):
///   `void name(long n, const double *restrict x_in, ..., double *restrict out0, ...)`
///   evaluates the expressions element-wise over the input arrays;
/// * scalar: `double name(double x, ...)` for one expression, otherwise
///   `void name(double x, ..., double *out0, ...)`.
pub fn ccode_function(name: &str, args: &[Expr], exprs: &[Expr], opts: &CodegenOptions) -> Result<String, CodeError> {
    let mut s = CodeSettings::new(Lang::C);
    s.real = opts.dtype;
    let p = CodePrinter::new(&s);
    let names = arg_names(&p, args)?;
    let b = body(&p, args, exprs, opts.cse)?;
    let t = opts.dtype.c_name();
    let mut out = String::new();
    if opts.include_header {
        out.push_str("#include <math.h>\n\n");
    }
    if opts.vectorized {
        let mut params = vec!["long n".to_string()];
        params.extend(names.iter().map(|a| format!("const {t} *restrict {a}_in")));
        params.extend((0..exprs.len()).map(|k| format!("{t} *restrict out{k}")));
        out.push_str(&format!("void {name}({}) {{\n", params.join(", ")));
        for d in &b.declarations {
            out.push_str(&format!("   {d}\n"));
        }
        out.push_str("   for (long i = 0; i < n; i++) {\n");
        for a in &names {
            out.push_str(&format!("      const {t} {a} = {a}_in[i];\n"));
        }
        for (v, e) in &b.temps {
            out.push_str(&indent(&format!("const {t} {v} = {e};"), "      "));
            out.push('\n');
        }
        for (k, e) in b.outputs.iter().enumerate() {
            out.push_str(&indent(&format!("out{k}[i] = {e};"), "      "));
            out.push('\n');
        }
        out.push_str("   }\n}\n");
    } else {
        let single = exprs.len() == 1;
        let mut params: Vec<String> = names.iter().map(|a| format!("{t} {a}")).collect();
        if !single {
            params.extend((0..exprs.len()).map(|k| format!("{t} *out{k}")));
        }
        let ret = if single { t } else { "void" };
        out.push_str(&format!("{ret} {name}({}) {{\n", params.join(", ")));
        for d in &b.declarations {
            out.push_str(&format!("   {d}\n"));
        }
        for (v, e) in &b.temps {
            out.push_str(&indent(&format!("const {t} {v} = {e};"), "   "));
            out.push('\n');
        }
        if single {
            out.push_str(&indent(&format!("return {};", b.outputs[0]), "   "));
            out.push('\n');
        } else {
            for (k, e) in b.outputs.iter().enumerate() {
                out.push_str(&indent(&format!("*out{k} = {e};"), "   "));
                out.push('\n');
            }
        }
        out.push_str("}\n");
    }
    Ok(out)
}

fn cuda_printer_settings(opts: &CodegenOptions) -> CodeSettings {
    CudaSettings {
        dtype: opts.dtype,
        fast_math: opts.fast_math,
        use_rsqrt: opts.use_rsqrt,
        inline_pow: opts.inline_pow,
    }
    .to_code_settings()
}

/// A CUDA `__device__` function returning the value of one expression, or
/// writing several outputs through pointers.
pub fn cuda_device_function(
    name: &str,
    args: &[Expr],
    exprs: &[Expr],
    opts: &CodegenOptions,
) -> Result<String, CodeError> {
    let s = cuda_printer_settings(opts);
    let p = CodePrinter::new(&s);
    let names = arg_names(&p, args)?;
    let b = body(&p, args, exprs, opts.cse)?;
    let t = opts.dtype.c_name();
    let single = exprs.len() == 1;
    let mut params: Vec<String> = names.iter().map(|a| format!("const {t} {a}")).collect();
    if !single {
        params.extend((0..exprs.len()).map(|k| format!("{t} *__restrict__ out{k}")));
    }
    let ret = if single { t } else { "void" };
    let mut out = format!("__device__ __forceinline__ {ret} {name}({}) {{\n", params.join(", "));
    for d in &b.declarations {
        out.push_str(&format!("   {d}\n"));
    }
    for (v, e) in &b.temps {
        out.push_str(&indent(&format!("const {t} {v} = {e};"), "   "));
        out.push('\n');
    }
    if single {
        out.push_str(&indent(&format!("return {};", b.outputs[0]), "   "));
        out.push('\n');
    } else {
        for (k, e) in b.outputs.iter().enumerate() {
            out.push_str(&indent(&format!("*out{k} = {e};"), "   "));
            out.push('\n');
        }
    }
    out.push_str("}\n");
    Ok(out)
}

/// A CUDA `__global__` kernel evaluating `exprs` element-wise over arrays
/// with a grid-stride loop:
///
/// ```text
/// __global__ void name(const long n, const double *__restrict__ x_in, ...,
///                      double *__restrict__ out0, ...)
/// ```
pub fn cuda_kernel(name: &str, args: &[Expr], exprs: &[Expr], opts: &CodegenOptions) -> Result<String, CodeError> {
    let s = cuda_printer_settings(opts);
    let p = CodePrinter::new(&s);
    let names = arg_names(&p, args)?;
    let b = body(&p, args, exprs, opts.cse)?;
    let t = opts.dtype.c_name();
    let mut params = vec!["const long n".to_string()];
    params.extend(names.iter().map(|a| format!("const {t} *__restrict__ {a}_in")));
    params.extend((0..exprs.len()).map(|k| format!("{t} *__restrict__ out{k}")));
    let mut out = format!("__global__ void {name}({}) {{\n", params.join(", "));
    for d in &b.declarations {
        out.push_str(&format!("   {d}\n"));
    }
    out.push_str(
        "   for (long i = (long)blockIdx.x * blockDim.x + threadIdx.x; i < n;\n        i += (long)blockDim.x * gridDim.x) {\n",
    );
    for a in &names {
        out.push_str(&format!("      const {t} {a} = {a}_in[i];\n"));
    }
    for (v, e) in &b.temps {
        out.push_str(&indent(&format!("const {t} {v} = {e};"), "      "));
        out.push('\n');
    }
    for (k, e) in b.outputs.iter().enumerate() {
        out.push_str(&indent(&format!("out{k}[i] = {e};"), "      "));
        out.push('\n');
    }
    out.push_str("   }\n}\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_str;

    #[test]
    fn c_function_shapes() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let e = parse_str("sin(x + y)**2 + cos(x + y)").unwrap();
        let src = ccode_function("f", &[x.clone(), y.clone()], &[e.clone()], &CodegenOptions::scalar()).unwrap();
        assert!(src.contains("double f(double x, double y) {"), "{src}");
        assert!(src.contains("const double x0 = x + y;"));
        assert!(src.contains("return pow(sin(x0), 2) + cos(x0);"), "{src}");
        let src = ccode_function("g", &[x.clone(), y.clone()], &[e.clone(), x.clone()], &CodegenOptions::default())
            .unwrap();
        assert!(src.contains("out1[i] = x;"));
        let k = cuda_kernel("k", &[x, y], &[e], &CodegenOptions::default()).unwrap();
        assert!(k.starts_with("__global__ void k(const long n"));
        assert!(k.contains("sin(x0)*sin(x0)") || k.contains("pow(sin(x0), 2)"));
    }

    #[test]
    fn excluded_arg_names() {
        let x0 = Expr::symbol("x0");
        let y = Expr::symbol("y");
        let e = parse_str("sin(y + 1) + cos(y + 1)").unwrap();
        let src = ccode_function("f", &[x0, y], &[e], &CodegenOptions::scalar()).unwrap();
        assert!(src.contains("const double x1 = y + 1;"), "{src}");
    }
}
