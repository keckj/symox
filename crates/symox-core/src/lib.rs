//! # symox
//!
//! Fast symbolic mathematics in pure Rust, with SymPy compatible semantics.
//!
//! ```
//! use symox::prelude::*;
//!
//! let x = Expr::symbol("x");
//! let y = Expr::symbol("y");
//! let e = (&x + &y).powi(2);
//! assert_eq!(expand(&e).to_string(), "x**2 + 2*x*y + y**2");
//! assert_eq!(diff(&sin(&x), &x).to_string(), "cos(x)");
//! ```

pub mod assumptions;
pub mod bigfloat;
pub mod diff;
pub mod evalf;
pub mod evalf_mp;
pub mod expand;
pub mod expr;
pub mod functions;
pub mod integrate;
pub mod lambdify;
pub mod ntheory;
pub mod limits;
pub mod number;
pub mod ode;
pub mod ops;
pub mod parser;
pub mod poly;
pub mod pattern;
pub mod printer;
pub mod printing;
pub mod relational;
pub mod simp;
pub mod series;
pub mod simplify;
pub mod solve;
pub mod sortkey;
pub mod subs;
pub mod summation;
pub mod symbol;
pub mod traverse;

pub use expr::{Constant, Expr, Func, Kind, RelOp};
pub use symbol::{Assumptions, Symbol};
pub use number::Number;

/// Everything needed for day to day use.
pub mod prelude {
    pub use crate::diff::{diff, diff_n};
    pub use crate::evalf::{eval_c64, eval_f64, evalf};
    pub use crate::expand::expand;
    pub use crate::expr::{Constant, Expr, Func, Kind, RelOp};
    pub use crate::functions::*;
    pub use crate::integrate::{integrate, integrate_definite, try_integrate};
    pub use crate::number::Number;
    pub use crate::ops::{add, div, mul, neg, pow, sqrt, sub};
    pub use crate::parser::{parse_str, sympify};
    pub use crate::printing::{
        ccode, ccode_function, cse, cuda_kernel, cudacode, cxxcode, latex, pprint, pretty, pycode,
        rust_code,
    };
    pub use crate::relational::{eq, ge, gt, le, lt, ne, piecewise};
    pub use crate::limits::{limit, Dir};
    pub use crate::series::{order, series};
    pub use crate::simplify::{radsimp, simplify, together};
    pub use crate::subs::{subs, xreplace};
    pub use crate::traverse::{free_symbols, has};
}
