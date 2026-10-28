//! Compilation of expressions to a fast numerical evaluator.
//!
//! The expression DAG is lowered to a register based instruction list with
//! common sub-expressions shared (thanks to hash-consed structural
//! equality). Evaluation is done column-wise over chunks of points, which
//! amortises dispatch and lets the compiler vectorise the inner loops.

use rustc_hash::FxHashMap;

use crate::special::numeric::{erf::{erf_real, erfc_real}, gamma::gamma_real};
use crate::special::SpecialFn;
use crate::expr::{Constant, Expr, Func, Kind};

#[derive(Clone, Debug)]
enum Exp {
    Int(i32),
    Half,
    MinusHalf,
    Float(f64),
    Reg(u32),
}

#[derive(Clone, Debug)]
enum Ins {
    Const(f64),
    Arg(u32),
    Sum(f64, Vec<(u32, f64)>),
    Prod(f64, Vec<(u32, Exp)>),
    Pow(u32, Exp),
    Call1(fn(f64) -> f64, u32),
    Call2(fn(f64, f64) -> f64, u32, u32),
    /// special function of any arity (real arguments and result)
    CallS(SpecialFn, Vec<u32>),
    Fold(bool, Vec<u32>), // max (true) / min (false)
    Select(Vec<(u32, u32)>),
    Cmp(crate::expr::RelOp, u32, u32),
}

#[derive(Debug)]
pub enum LambdifyError {
    UnknownSymbol(String),
    Unsupported(String),
}

impl std::error::Error for LambdifyError {}

impl std::fmt::Display for LambdifyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LambdifyError::UnknownSymbol(s) => write!(f, "free symbol '{s}' is not an argument"),
            LambdifyError::Unsupported(s) => write!(f, "cannot compile '{s}' to a numerical function"),
        }
    }
}

/// A compiled numerical function of `n_args` real arguments returning
/// `n_out` values.
#[derive(Clone, Debug)]
pub struct Compiled {
    code: Vec<Ins>,
    outputs: Vec<u32>,
    n_args: usize,
}

struct Builder<'a> {
    code: Vec<Ins>,
    memo: FxHashMap<Expr, u32>,
    args: &'a FxHashMap<Expr, u32>,
}

fn f1(f: &Func) -> Option<fn(f64) -> f64> {
    Some(match f {
        Func::Exp => f64::exp,
        Func::Log => f64::ln,
        Func::Sin => f64::sin,
        Func::Cos => f64::cos,
        Func::Tan => f64::tan,
        Func::Cot => |x| 1.0 / x.tan(),
        Func::Sec => |x| 1.0 / x.cos(),
        Func::Csc => |x| 1.0 / x.sin(),
        Func::Asin => f64::asin,
        Func::Acos => f64::acos,
        Func::Atan => f64::atan,
        Func::Acot => |x| (1.0 / x).atan(),
        Func::Sinh => f64::sinh,
        Func::Cosh => f64::cosh,
        Func::Tanh => f64::tanh,
        Func::Coth => |x| 1.0 / x.tanh(),
        Func::Asinh => f64::asinh,
        Func::Acosh => f64::acosh,
        Func::Atanh => f64::atanh,
        Func::Abs => f64::abs,
        Func::Sign => |x| {
            if x > 0.0 {
                1.0
            } else if x < 0.0 {
                -1.0
            } else {
                0.0
            }
        },
        Func::Floor => f64::floor,
        Func::Ceiling => f64::ceil,
        Func::Gamma => gamma_real,
        Func::LogGamma => |x| {
            crate::special::numeric::gamma::loggamma_c(num_complex::Complex64::new(x, 0.0))
                .filter(|z| z.im == 0.0)
                .map(|z| z.re)
                .unwrap_or(f64::NAN)
        },
        Func::Erf => erf_real,
        Func::Erfc => erfc_real,
        Func::Re | Func::Conjugate => |x| x,
        Func::Im => |_| 0.0,
        Func::Arg => |x| if x < 0.0 { std::f64::consts::PI } else { 0.0 },
        Func::Heaviside => |x| {
            if x > 0.0 {
                1.0
            } else if x < 0.0 {
                0.0
            } else {
                0.5
            }
        },
        _ => return None,
    })
}

impl Builder<'_> {
    fn push(&mut self, i: Ins) -> u32 {
        self.code.push(i);
        (self.code.len() - 1) as u32
    }

    fn exp_of(&mut self, x: &Expr) -> Result<Exp, LambdifyError> {
        if let Some(k) = x.as_i64() {
            if let Ok(k) = i32::try_from(k) {
                return Ok(Exp::Int(k));
            }
        }
        if *x == Expr::rational(1, 2) {
            return Ok(Exp::Half);
        }
        if *x == Expr::rational(-1, 2) {
            return Ok(Exp::MinusHalf);
        }
        if let Some(n) = x.as_number() {
            return Ok(Exp::Float(n.to_f64()));
        }
        Ok(Exp::Reg(self.lower(x)?))
    }

    fn lower(&mut self, e: &Expr) -> Result<u32, LambdifyError> {
        if let Some(r) = self.args.get(e) {
            return Ok(*r);
        }
        if let Some(r) = self.memo.get(e) {
            return Ok(*r);
        }
        let ins = match e.kind() {
            Kind::Num(n) => Ins::Const(n.to_f64()),
            Kind::Const(c) => Ins::Const(match c {
                Constant::Pi => std::f64::consts::PI,
                Constant::E => std::f64::consts::E,
                Constant::EulerGamma => crate::expr::EULER_GAMMA_F64,
                Constant::Catalan | Constant::GoldenRatio | Constant::TribonacciConstant => c.approx_f64().unwrap(),
                Constant::Infinity => f64::INFINITY,
                Constant::NegInfinity => f64::NEG_INFINITY,
                Constant::NaN | Constant::ComplexInfinity => f64::NAN,
                Constant::I => return Err(LambdifyError::Unsupported("I (complex)".into())),
            }),
            Kind::Bool(b) => Ins::Const(if *b { 1.0 } else { 0.0 }),
            Kind::Sym(s) => return Err(LambdifyError::UnknownSymbol(s.to_string())),
            Kind::Add(c, terms) => {
                let mut v = Vec::with_capacity(terms.len());
                for (t, k) in terms {
                    v.push((self.lower(t)?, k.to_f64()));
                }
                Ins::Sum(c.to_f64(), v)
            }
            Kind::Mul(c, facs) => {
                let mut v = Vec::with_capacity(facs.len());
                for (b, x) in facs {
                    let rb = self.lower(b)?;
                    let ex = self.exp_of(x)?;
                    v.push((rb, ex));
                }
                Ins::Prod(c.to_f64(), v)
            }
            Kind::Pow(b, x) => {
                let rb = self.lower(b)?;
                let ex = self.exp_of(x)?;
                Ins::Pow(rb, ex)
            }
            Kind::Func(f, args) => match f {
                Func::Atan2 => {
                    let (a, b) = (self.lower(&args[0])?, self.lower(&args[1])?);
                    Ins::Call2(f64::atan2, a, b)
                }
                Func::Max | Func::Min => {
                    let regs: Result<Vec<u32>, _> = args.iter().map(|a| self.lower(a)).collect();
                    Ins::Fold(matches!(f, Func::Max), regs?)
                }
                // booleans are 0.0 / 1.0: And = min, Or = max, Not = 1 - a
                Func::Logic(op @ (crate::logic::LogicOp::And | crate::logic::LogicOp::Or)) => {
                    let regs: Result<Vec<u32>, _> = args.iter().map(|a| self.lower(a)).collect();
                    Ins::Fold(*op == crate::logic::LogicOp::Or, regs?)
                }
                Func::Logic(crate::logic::LogicOp::Not) => Ins::Sum(1.0, vec![(self.lower(&args[0])?, -1.0)]),
                Func::Logic(_) => return self.lower(&crate::logic::eliminate_implications(e)),
                // Heaviside(x, H0): the value at 0 is H0
                Func::Heaviside if args.len() == 2 => {
                    let piecewise = crate::special::rewrite::rewrite(e, &[], "Piecewise", false);
                    if &piecewise == e {
                        return Err(LambdifyError::Unsupported(e.to_string()));
                    }
                    return self.lower(&piecewise);
                }
                Func::Special(s) => {
                    if matches!(s, SpecialFn::Hyper | SpecialFn::Bell) && args.len() == 3 {
                        return Err(LambdifyError::Unsupported(e.to_string()));
                    }
                    let regs: Result<Vec<u32>, _> = args.iter().map(|a| self.lower(a)).collect();
                    Ins::CallS(*s, regs?)
                }
                _ => match f1(f) {
                    Some(g) => {
                        let a = self.lower(&args[0])?;
                        Ins::Call1(g, a)
                    }
                    None => return Err(LambdifyError::Unsupported(e.to_string())),
                },
            },
            Kind::Rel(op, a, b) => {
                let (ra, rb) = (self.lower(a)?, self.lower(b)?);
                Ins::Cmp(*op, ra, rb)
            }
            Kind::Piecewise(p) => {
                let mut v = Vec::with_capacity(p.len());
                for (x, c) in p {
                    v.push((self.lower(x)?, self.lower(c)?));
                }
                Ins::Select(v)
            }
            Kind::Derivative(..) => return Err(LambdifyError::Unsupported(e.to_string())),
        };
        let r = self.push(ins);
        self.memo.insert(e.clone(), r);
        Ok(r)
    }
}

impl Compiled {
    /// Compile `exprs` as functions of `args` (symbols, or any expression
    /// treated as an opaque input such as `f(x)`).
    pub fn new(args: &[Expr], exprs: &[Expr]) -> Result<Compiled, LambdifyError> {
        let mut amap = FxHashMap::default();
        let mut code = Vec::new();
        for (i, a) in args.iter().enumerate() {
            amap.insert(a.clone(), i as u32);
            code.push(Ins::Arg(i as u32));
        }
        let mut b = Builder {
            code,
            memo: FxHashMap::default(),
            args: &amap,
        };
        let mut outputs = Vec::with_capacity(exprs.len());
        for e in exprs {
            outputs.push(b.lower(e)?);
        }
        Ok(Compiled {
            code: b.code,
            outputs,
            n_args: args.len(),
        })
    }

    pub fn n_args(&self) -> usize {
        self.n_args
    }

    pub fn n_outputs(&self) -> usize {
        self.outputs.len()
    }

    /// Evaluate at a single point.
    pub fn call(&self, args: &[f64], out: &mut [f64]) {
        let mut regs = vec![0.0f64; self.code.len()];
        for (i, ins) in self.code.iter().enumerate() {
            regs[i] = match ins {
                Ins::Arg(a) => args[*a as usize],
                Ins::Const(c) => *c,
                Ins::Sum(c, t) => t.iter().fold(*c, |s, (r, k)| s + k * regs[*r as usize]),
                Ins::Prod(c, f) => f
                    .iter()
                    .fold(*c, |p, (r, x)| p * pow_scalar(regs[*r as usize], x, &regs)),
                Ins::Pow(r, x) => pow_scalar(regs[*r as usize], x, &regs),
                Ins::Call1(g, r) => g(regs[*r as usize]),
                Ins::Call2(g, a, b) => g(regs[*a as usize], regs[*b as usize]),
                Ins::CallS(s, rs) => call_special(*s, rs.iter().map(|r| regs[*r as usize])),
                Ins::Fold(is_max, rs) => fold(*is_max, rs.iter().map(|r| regs[*r as usize])),
                Ins::Cmp(op, a, b) => cmp(*op, regs[*a as usize], regs[*b as usize]),
                Ins::Select(p) => select(p.iter().map(|(x, c)| (regs[*x as usize], regs[*c as usize]))),
            };
        }
        for (o, r) in out.iter_mut().zip(&self.outputs) {
            *o = regs[*r as usize];
        }
    }

    /// Evaluate at `n` points. `args[j]` holds the `n` values of argument
    /// `j`; `out[k]` receives the `n` values of output `k`.
    pub fn call_many(&self, args: &[&[f64]], out: &mut [&mut [f64]]) {
        const CHUNK: usize = 256;
        let n = args.first().map(|a| a.len()).unwrap_or(0);
        let mut regs: Vec<[f64; CHUNK]> = vec![[0.0; CHUNK]; self.code.len()];
        let mut start = 0;
        while start < n {
            let len = CHUNK.min(n - start);
            for i in 0..self.code.len() {
                let (prev, cur) = regs.split_at_mut(i);
                let dst = &mut cur[0][..len];
                match &self.code[i] {
                    Ins::Arg(a) => dst.copy_from_slice(&args[*a as usize][start..start + len]),
                    Ins::Const(c) => dst.fill(*c),
                    Ins::Sum(c, t) => {
                        dst.fill(*c);
                        for (r, k) in t {
                            let src = &prev[*r as usize][..len];
                            for (d, s) in dst.iter_mut().zip(src) {
                                *d += k * s;
                            }
                        }
                    }
                    Ins::Prod(c, f) => {
                        dst.fill(*c);
                        for (r, x) in f {
                            let src = &prev[*r as usize][..len];
                            mul_pow_into(dst, src, x, prev, len);
                        }
                    }
                    Ins::Pow(r, x) => {
                        dst.fill(1.0);
                        let src = &prev[*r as usize][..len];
                        mul_pow_into(dst, src, x, prev, len);
                    }
                    Ins::Call1(g, r) => {
                        let src = &prev[*r as usize][..len];
                        for (d, s) in dst.iter_mut().zip(src) {
                            *d = g(*s);
                        }
                    }
                    Ins::Call2(g, a, b) => {
                        for j in 0..len {
                            dst[j] = g(prev[*a as usize][j], prev[*b as usize][j]);
                        }
                    }
                    Ins::CallS(s, rs) => {
                        for j in 0..len {
                            dst[j] = call_special(*s, rs.iter().map(|r| prev[*r as usize][j]));
                        }
                    }
                    Ins::Fold(is_max, rs) => {
                        for j in 0..len {
                            dst[j] = fold(*is_max, rs.iter().map(|r| prev[*r as usize][j]));
                        }
                    }
                    Ins::Cmp(op, a, b) => {
                        for j in 0..len {
                            dst[j] = cmp(*op, prev[*a as usize][j], prev[*b as usize][j]);
                        }
                    }
                    Ins::Select(p) => {
                        for j in 0..len {
                            dst[j] = select(
                                p.iter().map(|(x, c)| (prev[*x as usize][j], prev[*c as usize][j])),
                            );
                        }
                    }
                }
            }
            for (k, r) in self.outputs.iter().enumerate() {
                out[k][start..start + len].copy_from_slice(&regs[*r as usize][..len]);
            }
            start += len;
        }
    }
}

#[inline]
fn powi(x: f64, k: i32) -> f64 {
    match k {
        1 => x,
        2 => x * x,
        3 => x * x * x,
        -1 => 1.0 / x,
        -2 => 1.0 / (x * x),
        _ => x.powi(k),
    }
}

#[inline]
fn pow_scalar(b: f64, x: &Exp, regs: &[f64]) -> f64 {
    match x {
        Exp::Int(k) => powi(b, *k),
        Exp::Half => b.sqrt(),
        Exp::MinusHalf => 1.0 / b.sqrt(),
        Exp::Float(f) => b.powf(*f),
        Exp::Reg(r) => b.powf(regs[*r as usize]),
    }
}

#[inline]
fn mul_pow_into(dst: &mut [f64], src: &[f64], x: &Exp, prev: &[[f64; 256]], len: usize) {
    match x {
        Exp::Int(1) => dst.iter_mut().zip(src).for_each(|(d, s)| *d *= s),
        Exp::Int(2) => dst.iter_mut().zip(src).for_each(|(d, s)| *d *= s * s),
        Exp::Int(-1) => dst.iter_mut().zip(src).for_each(|(d, s)| *d /= s),
        Exp::Int(k) => dst.iter_mut().zip(src).for_each(|(d, s)| *d *= s.powi(*k)),
        Exp::Half => dst.iter_mut().zip(src).for_each(|(d, s)| *d *= s.sqrt()),
        Exp::MinusHalf => dst.iter_mut().zip(src).for_each(|(d, s)| *d /= s.sqrt()),
        Exp::Float(f) => dst.iter_mut().zip(src).for_each(|(d, s)| *d *= s.powf(*f)),
        Exp::Reg(r) => {
            let e = &prev[*r as usize][..len];
            for j in 0..len {
                dst[j] *= src[j].powf(e[j]);
            }
        }
    }
}

fn call_special(s: SpecialFn, args: impl Iterator<Item = f64>) -> f64 {
    let mut buf = [0.0f64; 4];
    let mut n = 0;
    for a in args {
        if n < 4 {
            buf[n] = a;
        }
        n += 1;
    }
    match n {
        1 => crate::special::numeric::real1(s, buf[0]),
        2 => crate::special::numeric::real2(s, buf[0], buf[1]),
        3 | 4 => crate::special::numeric::realn(s, &buf[..n]),
        _ => f64::NAN,
    }
}

#[inline]
fn fold(is_max: bool, it: impl Iterator<Item = f64>) -> f64 {
    let mut acc = if is_max { f64::NEG_INFINITY } else { f64::INFINITY };
    for v in it {
        acc = if is_max { acc.max(v) } else { acc.min(v) };
    }
    acc
}

#[inline]
fn cmp(op: crate::expr::RelOp, a: f64, b: f64) -> f64 {
    use crate::expr::RelOp::*;
    let r = match op {
        Eq => a == b,
        Ne => a != b,
        Lt => a < b,
        Le => a <= b,
        Gt => a > b,
        Ge => a >= b,
    };
    if r {
        1.0
    } else {
        0.0
    }
}

#[inline]
fn select(it: impl Iterator<Item = (f64, f64)>) -> f64 {
    for (x, c) in it {
        if c != 0.0 {
            return x;
        }
    }
    f64::NAN
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;

    #[test]
    fn compile_and_eval() {
        let x = Expr::symbol("x");
        let y = Expr::symbol("y");
        let e = add(&[sin(&x).powi(2), mul(&[Expr::int(3), x.clone(), y.clone()]), sqrt(&y)]);
        let c = Compiled::new(&[x.clone(), y.clone()], &[e]).unwrap();
        let mut out = [0.0];
        c.call(&[0.5, 4.0], &mut out);
        let expect = 0.5f64.sin().powi(2) + 6.0 + 2.0;
        assert!((out[0] - expect).abs() < 1e-14);
        let xs: Vec<f64> = (0..1000).map(|i| i as f64 * 0.001).collect();
        let ys: Vec<f64> = (0..1000).map(|i| 1.0 + i as f64).collect();
        let mut o = vec![0.0; 1000];
        c.call_many(&[&xs, &ys], &mut [&mut o]);
        for i in 0..1000 {
            let v = xs[i].sin().powi(2) + 3.0 * xs[i] * ys[i] + ys[i].sqrt();
            assert!((o[i] - v).abs() < 1e-12);
        }
    }
}
