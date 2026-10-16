//! Special functions: combinatorial numbers, gamma and zeta families,
//! error functions and exponential integrals, Bessel and Airy functions,
//! elliptic integrals, Lambert W, generalised hypergeometric functions and
//! orthogonal polynomials.
//!
//! Every function is a head [`SpecialFn`] wrapped in
//! [`crate::expr::Func::Special`]. Construction goes through
//! [`crate::functions::func`] which calls [`apply`]: SymPy's automatic
//! evaluation rules are applied (`factorial(5) -> 120`, `zeta(2) ->
//! pi**2/6`, ...) and, like SymPy, the function is evaluated numerically
//! when *all* its arguments are floats.
//!
//! Sub-modules:
//! * `combinatorial`, `gamma`, `zeta`, `error`, `bessel`, `elliptic`,
//!   `hyper`, `orthopoly`: symbolic evaluation rules,
//! * [`diff`]: partial derivatives,
//! * [`numeric`]: double precision (real and complex) evaluation,
//! * [`mp`]: arbitrary precision evaluation,
//! * [`expand`]: `expand_func`, `besselsimp`, `hyperexpand`,
//! * [`integers`]: exact integer sequences (Stirling, Bell, partitions...).

use crate::expr::{Constant, Expr, Func, Kind};
use crate::functions::func;
use crate::number::Number;

mod bessel;
mod combinatorial;
pub mod diff;
mod elliptic;
mod error;
pub mod expand;
mod gamma;
mod hyper;
pub mod integers;
pub mod integrate;
pub mod mp;
pub mod numeric;
mod orthopoly;
mod zeta;

pub use expand::{besselsimp, expand_func, hyperexpand};
pub use integers::{bell_number, multinomial_coefficients, partition_number, stirling};

macro_rules! special_fns {
    ($($v:ident => $name:literal [$($n:literal),*] $auto:literal;)*) => {
        /// Special function heads (see the module documentation).
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
        pub enum SpecialFn {
            $($v),*
        }

        impl SpecialFn {
            /// All heads.
            pub const ALL: &'static [SpecialFn] = &[$(SpecialFn::$v),*];

            /// SymPy name (as printed).
            pub fn name(self) -> &'static str {
                match self {
                    $(SpecialFn::$v => $name),*
                }
            }

            /// Accepted numbers of arguments.
            pub fn nargs(self) -> &'static [usize] {
                match self {
                    $(SpecialFn::$v => &[$($n),*]),*
                }
            }

            /// Whether SymPy evaluates the function numerically when all
            /// arguments are floats.
            pub fn auto_evalf(self) -> bool {
                match self {
                    $(SpecialFn::$v => $auto),*
                }
            }

            /// Look up a head by SymPy name (with the usual aliases).
            pub fn from_name(name: &str) -> Option<SpecialFn> {
                Some(match name {
                    $($name => SpecialFn::$v,)*
                    "rf" => SpecialFn::RisingFactorial,
                    "ff" => SpecialFn::FallingFactorial,
                    _ => return None,
                })
            }
        }
    };
}

special_fns! {
    // combinatorial
    Factorial => "factorial" [1] true;
    Factorial2 => "factorial2" [1] true;
    Subfactorial => "subfactorial" [1] false;
    Binomial => "binomial" [2] true;
    RisingFactorial => "RisingFactorial" [2] true;
    FallingFactorial => "FallingFactorial" [2] true;
    Harmonic => "harmonic" [1, 2] true;
    Bernoulli => "bernoulli" [1, 2] true;
    Euler => "euler" [1, 2] true;
    Catalan => "catalan" [1] true;
    Fibonacci => "fibonacci" [1, 2] true;
    Lucas => "lucas" [1] false;
    Bell => "bell" [1, 2, 3] false;
    Partition => "partition" [1] false;
    // gamma family
    Polygamma => "polygamma" [2] true;
    Digamma => "digamma" [1] true;
    Trigamma => "trigamma" [1] true;
    Beta => "beta" [1, 2] true;
    LowerGamma => "lowergamma" [2] true;
    UpperGamma => "uppergamma" [2] true;
    // zeta family
    Zeta => "zeta" [1, 2] true;
    DirichletEta => "dirichlet_eta" [1] true;
    Polylog => "polylog" [2] true;
    Lerchphi => "lerchphi" [3] true;
    // error functions and exponential integrals
    Erfi => "erfi" [1] true;
    Erf2 => "erf2" [2] false;
    Erfinv => "erfinv" [1] true;
    Erfcinv => "erfcinv" [1] true;
    FresnelS => "fresnels" [1] true;
    FresnelC => "fresnelc" [1] true;
    Ei => "Ei" [1] true;
    Li => "li" [1] true;
    LiOffset => "Li" [1] true;
    Si => "Si" [1] true;
    Ci => "Ci" [1] true;
    Shi => "Shi" [1] true;
    Chi => "Chi" [1] true;
    Expint => "expint" [2] true;
    E1 => "E1" [1] true;
    // Bessel type
    BesselJ => "besselj" [2] true;
    BesselY => "bessely" [2] true;
    BesselI => "besseli" [2] true;
    BesselK => "besselk" [2] true;
    Hankel1 => "hankel1" [2] true;
    Hankel2 => "hankel2" [2] true;
    Jn => "jn" [2] true;
    Yn => "yn" [2] true;
    AiryAi => "airyai" [1] true;
    AiryBi => "airybi" [1] true;
    AiryAiPrime => "airyaiprime" [1] true;
    AiryBiPrime => "airybiprime" [1] true;
    // elliptic integrals, Lambert W, hypergeometric
    EllipticK => "elliptic_k" [1] true;
    EllipticE => "elliptic_e" [1, 2] true;
    EllipticF => "elliptic_f" [2] true;
    EllipticPi => "elliptic_pi" [2, 3] true;
    LambertW => "LambertW" [1, 2] true;
    Hyper => "hyper" [3] true;
    // orthogonal polynomials
    Legendre => "legendre" [2] true;
    AssocLegendre => "assoc_legendre" [3] true;
    ChebyshevT => "chebyshevt" [2] true;
    ChebyshevU => "chebyshevu" [2] true;
    Hermite => "hermite" [2] true;
    Laguerre => "laguerre" [2] true;
    AssocLaguerre => "assoc_laguerre" [3] true;
    Jacobi => "jacobi" [4] true;
    Gegenbauer => "gegenbauer" [3] true;
    Ynm => "Ynm" [4] true;
}

impl SpecialFn {
    /// Fixed number of arguments, `None` when several are accepted.
    pub fn arity(self) -> Option<usize> {
        match self.nargs() {
            [n] => Some(*n),
            _ => None,
        }
    }
}

/// Build the unevaluated node.
pub(crate) fn raw(s: SpecialFn, args: Vec<Expr>) -> Expr {
    Expr::from_kind(Kind::Func(Func::Special(s), args))
}

/// Apply special function `s` to `args` with automatic evaluation.
pub fn apply(s: SpecialFn, args: Vec<Expr>) -> Expr {
    if !s.nargs().contains(&args.len()) {
        return raw(s, args);
    }
    if let Some(v) = eval(s, &args) {
        return v;
    }
    if s.auto_evalf() && should_evalf(&args) {
        if let Some(v) = numeric::eval_exprs(s, &args) {
            return v;
        }
    }
    raw(s, args)
}

fn eval(s: SpecialFn, a: &[Expr]) -> Option<Expr> {
    use SpecialFn::*;
    match s {
        Factorial | Factorial2 | Subfactorial | Binomial | RisingFactorial | FallingFactorial
        | Harmonic | Bernoulli | Euler | Catalan | Fibonacci | Lucas | Bell | Partition => {
            combinatorial::eval(s, a)
        }
        Polygamma | Digamma | Trigamma | Beta | LowerGamma | UpperGamma => gamma::eval(s, a),
        Zeta | DirichletEta | Polylog | Lerchphi => zeta::eval(s, a),
        Erfi | Erf2 | Erfinv | Erfcinv | FresnelS | FresnelC | Ei | Li | LiOffset | Si | Ci | Shi
        | Chi | Expint | E1 => error::eval(s, a),
        BesselJ | BesselY | BesselI | BesselK | Hankel1 | Hankel2 | Jn | Yn | AiryAi | AiryBi
        | AiryAiPrime | AiryBiPrime => bessel::eval(s, a),
        EllipticK | EllipticE | EllipticF | EllipticPi | LambertW => elliptic::eval(s, a),
        Hyper => hyper::eval(a),
        Legendre | AssocLegendre | ChebyshevT | ChebyshevU | Hermite | Laguerre | AssocLaguerre
        | Jacobi | Gegenbauer | Ynm => orthopoly::eval(s, a),
    }
}

/// Automatic evaluation of the core `gamma` head (exact values).
pub(crate) fn eval_gamma(a: &Expr) -> Option<Expr> {
    gamma::eval_gamma(a)
}

/// Automatic evaluation of the core `loggamma` head.
pub(crate) fn eval_loggamma(a: &Expr) -> Option<Expr> {
    gamma::eval_loggamma(a)
}

/// Automatic evaluation of the core `erfc` head.
pub(crate) fn eval_erfc(a: &Expr) -> Option<Expr> {
    error::eval_erfc(a)
}

/// Automatic evaluation of the core `erf` head.
pub(crate) fn eval_erf(a: &Expr) -> Option<Expr> {
    error::eval_erf(a)
}

/// SymPy's `Function._should_evalf`: every argument is a float or a
/// complex number with float parts.
pub(crate) fn should_evalf(args: &[Expr]) -> bool {
    !args.is_empty() && args.iter().all(is_float_like)
}

fn is_float_like(e: &Expr) -> bool {
    match e.kind() {
        Kind::Num(Number::Float(_)) => true,
        Kind::Mul(c, f) => c.is_float() && f.len() == 1 && f[0].0.is_constant(Constant::I) && f[0].1.is_one(),
        Kind::Add(c, t) => {
            c.is_float()
                && t.len() == 1
                && t[0].0.is_constant(Constant::I)
                && t[0].1.is_float()
        }
        _ => false,
    }
}

// ----------------------------------------------------------------------
// small helpers shared by the sub-modules
// ----------------------------------------------------------------------

pub(crate) fn int(i: i64) -> Expr {
    Expr::int(i)
}

pub(crate) fn rat(n: i64, d: i64) -> Expr {
    Expr::rational(n, d)
}

/// Exact integer value.
pub(crate) fn as_int(e: &Expr) -> Option<num_bigint::BigInt> {
    match e.kind() {
        Kind::Num(n) if n.is_integer() && !n.is_float() => n.to_bigint(),
        _ => None,
    }
}

/// Exact machine integer value.
pub(crate) fn as_small(e: &Expr) -> Option<i64> {
    match e.kind() {
        Kind::Num(Number::Int(i)) => Some(*i),
        _ => None,
    }
}

/// Exact rational (non float) number.
pub(crate) fn as_exact(e: &Expr) -> Option<&Number> {
    match e.kind() {
        Kind::Num(n) if !n.is_float() => Some(n),
        _ => None,
    }
}

pub(crate) fn is_int_number(e: &Expr) -> bool {
    matches!(e.kind(), Kind::Num(n) if !n.is_float() && n.is_integer())
}

pub(crate) fn is_integer(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(!n.is_float() && n.is_integer()),
        _ => crate::assumptions::is_integer(e),
    }
}

pub(crate) fn is_positive(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(n.is_positive()),
        _ => crate::assumptions::is_positive(e),
    }
}

pub(crate) fn is_negative(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(n.is_negative()),
        _ => crate::assumptions::is_negative(e),
    }
}

pub(crate) fn is_nonnegative(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(!n.is_negative()),
        _ => crate::assumptions::is_nonnegative(e),
    }
}

pub(crate) fn is_nonpositive(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(!n.is_positive()),
        _ => crate::assumptions::is_nonpositive(e),
    }
}

pub(crate) fn is_zero(e: &Expr) -> Option<bool> {
    match e.kind() {
        Kind::Num(n) => Some(n.is_zero()),
        _ => crate::assumptions::is_zero(e),
    }
}

pub(crate) fn is_inf(e: &Expr) -> bool {
    e.is_constant(Constant::Infinity)
}

pub(crate) fn is_neg_inf(e: &Expr) -> bool {
    e.is_constant(Constant::NegInfinity)
}

/// `e / I` when `e` is syntactically `I * something` (SymPy's
/// `extract_multiplicatively(I)` for the common cases).
pub(crate) fn extract_i(e: &Expr) -> Option<Expr> {
    match e.kind() {
        Kind::Const(Constant::I) => Some(Expr::one()),
        Kind::Mul(c, facs) => {
            let idx = facs.iter().position(|(b, x)| b.is_constant(Constant::I) && x.is_one())?;
            let mut v: Vec<Expr> = Vec::with_capacity(facs.len());
            v.push(Expr::num(c.clone()));
            for (j, (b, x)) in facs.iter().enumerate() {
                if j != idx {
                    v.push(crate::ops::mul_factor(b, x));
                }
            }
            Some(crate::ops::mul(&v))
        }
        _ => None,
    }
}

pub(crate) fn could_extract_minus_sign(e: &Expr) -> bool {
    crate::ops::could_extract_minus_sign(e)
}

macro_rules! ctor {
    ($($fname:ident($($a:ident),*) => $v:ident;)*) => {
        $(
            #[doc = concat!("`", stringify!($fname), "(", stringify!($($a),*), ")`")]
            pub fn $fname($($a: &Expr),*) -> Expr {
                func(Func::Special(SpecialFn::$v), vec![$($a.clone()),*])
            }
        )*
    };
}

ctor! {
    factorial(n) => Factorial;
    factorial2(n) => Factorial2;
    subfactorial(n) => Subfactorial;
    binomial(n, k) => Binomial;
    rf(x, k) => RisingFactorial;
    ff(x, k) => FallingFactorial;
    harmonic(n) => Harmonic;
    harmonic2(n, m) => Harmonic;
    bernoulli(n) => Bernoulli;
    bernoulli_poly(n, x) => Bernoulli;
    euler(n) => Euler;
    euler_poly(n, x) => Euler;
    catalan(n) => Catalan;
    fibonacci(n) => Fibonacci;
    fibonacci_poly(n, x) => Fibonacci;
    lucas(n) => Lucas;
    bell(n) => Bell;
    bell_poly(n, x) => Bell;
    partition(n) => Partition;
    polygamma(n, z) => Polygamma;
    digamma(z) => Digamma;
    trigamma(z) => Trigamma;
    beta(x, y) => Beta;
    lowergamma(a, x) => LowerGamma;
    uppergamma(a, x) => UpperGamma;
    zeta(s) => Zeta;
    hurwitz_zeta(s, a) => Zeta;
    dirichlet_eta(s) => DirichletEta;
    polylog(s, z) => Polylog;
    lerchphi(z, s, a) => Lerchphi;
    erfi(z) => Erfi;
    erf2(x, y) => Erf2;
    erfinv(z) => Erfinv;
    erfcinv(z) => Erfcinv;
    fresnels(z) => FresnelS;
    fresnelc(z) => FresnelC;
    ei(z) => Ei;
    li(z) => Li;
    li_offset(z) => LiOffset;
    si(z) => Si;
    ci(z) => Ci;
    shi(z) => Shi;
    chi(z) => Chi;
    expint(nu, z) => Expint;
    e1(z) => E1;
    besselj(nu, z) => BesselJ;
    bessely(nu, z) => BesselY;
    besseli(nu, z) => BesselI;
    besselk(nu, z) => BesselK;
    hankel1(nu, z) => Hankel1;
    hankel2(nu, z) => Hankel2;
    jn(n, z) => Jn;
    yn(n, z) => Yn;
    airyai(z) => AiryAi;
    airybi(z) => AiryBi;
    airyaiprime(z) => AiryAiPrime;
    airybiprime(z) => AiryBiPrime;
    elliptic_k(m) => EllipticK;
    elliptic_e(m) => EllipticE;
    elliptic_e2(z, m) => EllipticE;
    elliptic_f(z, m) => EllipticF;
    elliptic_pi(n, m) => EllipticPi;
    elliptic_pi3(n, z, m) => EllipticPi;
    lambertw(z) => LambertW;
    lambertw_k(z, k) => LambertW;
    legendre(n, x) => Legendre;
    assoc_legendre(n, m, x) => AssocLegendre;
    chebyshevt(n, x) => ChebyshevT;
    chebyshevu(n, x) => ChebyshevU;
    hermite(n, x) => Hermite;
    laguerre(n, x) => Laguerre;
    assoc_laguerre(n, a, x) => AssocLaguerre;
    jacobi(n, a, b, x) => Jacobi;
    gegenbauer(n, a, x) => Gegenbauer;
    ynm(n, m, theta, phi) => Ynm;
}

/// `bell(n, k, (x_1, ..., x_m))`: incomplete Bell polynomial.
pub fn bell_incomplete(n: &Expr, k: &Expr, xs: &[Expr]) -> Expr {
    func(Func::Special(SpecialFn::Bell), vec![n.clone(), k.clone(), tuple(xs.to_vec())])
}

/// `hyper(ap, bq, z)`: generalised hypergeometric function.
pub fn hyper(ap: &[Expr], bq: &[Expr], z: &Expr) -> Expr {
    func(Func::Special(SpecialFn::Hyper), vec![tuple(ap.to_vec()), tuple(bq.to_vec()), z.clone()])
}

/// A `Tuple(...)` argument (printed `(a, b)`).
pub fn tuple(v: Vec<Expr>) -> Expr {
    crate::functions::undefined("Tuple", v)
}

/// Elements of a `Tuple(...)` argument.
pub fn tuple_items(e: &Expr) -> Option<&[Expr]> {
    match e.kind() {
        Kind::Func(Func::Undefined(n), a) if &**n == "Tuple" => Some(a),
        _ => None,
    }
}

/// Does `e` contain an application of `s`?
pub(crate) fn has_head(e: &Expr, s: SpecialFn) -> bool {
    match e.kind() {
        Kind::Func(Func::Special(g), _) if *g == s => true,
        _ if e.is_atom() => false,
        _ => e.args().iter().any(|a| has_head(a, s)),
    }
}

/// Build `s(args)` (shorthand used across the sub-modules).
pub(crate) fn sf(s: SpecialFn, args: Vec<Expr>) -> Expr {
    func(Func::Special(s), args)
}

#[cfg(test)]
mod tests;
