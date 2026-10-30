//! NumPy float64 array evaluation of a printed lambdify expression.
//!
//! `lambdify(..., 'numpy')` evaluates elementary expressions on float64
//! arrays with the NumPy code printed by the `NumPyPrinter`
//! (`numpy.sin(x) + 3*x**2`). This module evaluates that printed text on the
//! array buffers with NumPy's own semantics, operation by operation over
//! chunks of the arrays, without the temporaries of NumPy's ufunc calls:
//!
//! - Python literals and their combinations follow Python's arithmetic (exact
//!   integers below 2**53, true division, C `pow`), as when the generated
//!   code runs;
//! - scalar arguments (passed as `numpy.float64` by the caller) follow NumPy's
//!   scalar arithmetic, where `**` is C `pow`;
//! - array operations follow the NumPy ufunc loops: IEEE `+ - * /`, the
//!   `**` fast paths of `ndarray.__pow__` (`x**2` is `square`, `x**-1`
//!   `reciprocal`, `x**0.5` `sqrt`), C `pow` otherwise, and the C library
//!   `sin`, `cos`, `exp`, `log` that NumPy's float64 loops call on processors
//!   without AVX-512 (on those, NumPy uses its own SIMD kernels: the caller
//!   then passes `libm = false` and only exact operations are accepted).
//!
//! The result is NumPy's bit for bit. Anything else (other functions or
//! names, comparisons, Python errors in literal arithmetic, a non-finite
//! literal combination) is declined, and the caller keeps NumPy's evaluation.

use crate::printing::pyeval::{self, Bin, Node};

mod libm {
    extern "C" {
        pub fn sin(x: f64) -> f64;
        pub fn cos(x: f64) -> f64;
        pub fn exp(x: f64) -> f64;
        pub fn log(x: f64) -> f64;
        pub fn pow(x: f64, y: f64) -> f64;
    }
}

/// The NumPy functions of the printed text that are evaluated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unary {
    Neg,
    Sin,
    Cos,
    Exp,
    Log,
    Sqrt,
    Square,
    Reciprocal,
}

impl Unary {
    fn from_name(name: &str) -> Option<Unary> {
        Some(match name {
            "sin" => Unary::Sin,
            "cos" => Unary::Cos,
            "exp" => Unary::Exp,
            "log" => Unary::Log,
            "sqrt" => Unary::Sqrt,
            _ => return None,
        })
    }

    /// Whether the value comes from the C library (not exactly rounded).
    fn libm(self) -> bool {
        matches!(self, Unary::Sin | Unary::Cos | Unary::Exp | Unary::Log)
    }

    #[inline(always)]
    fn apply(self, x: f64) -> f64 {
        match self {
            Unary::Neg => -x,
            Unary::Sin => unsafe { libm::sin(x) },
            Unary::Cos => unsafe { libm::cos(x) },
            Unary::Exp => unsafe { libm::exp(x) },
            Unary::Log => unsafe { libm::log(x) },
            Unary::Sqrt => x.sqrt(),
            Unary::Square => x * x,
            Unary::Reciprocal => 1.0 / x,
        }
    }
}

/// The NumPy names a printed text may use: functions (`numpy.sin`) and
/// float constants (`numpy.pi`).
#[derive(Clone, Debug)]
pub struct NumpyText {
    root: Node,
    names: Vec<String>,
}

/// Parse `lambda a, b: body` printed with `numpy.` prefixes; `None` when
/// the text is not of that form.
pub fn parse_lambda(source: &str) -> Option<NumpyText> {
    let rest = source.trim().strip_prefix("lambda")?;
    let (head, body) = rest.split_once(':')?;
    let args: Vec<String> = head.split(',').map(|a| a.trim().to_string()).filter(|a| !a.is_empty()).collect();
    if args.iter().any(|a| !a.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')) {
        return None;
    }
    let parsed = pyeval::parse(body.trim(), &args)?;
    // the only global name is the `numpy` module
    if parsed.globals().iter().any(|g| g != "numpy") {
        return None;
    }
    let mut names = Vec::new();
    collect_attrs(parsed.root(), &mut names);
    Some(NumpyText { root: parsed.root().clone(), names })
}

fn collect_attrs(node: &Node, out: &mut Vec<String>) {
    match node {
        Node::Attr(_, name) => {
            if !out.contains(name) {
                out.push(name.clone());
            }
        }
        Node::Neg(a) | Node::Pos(a) | Node::Not(a) => collect_attrs(a, out),
        Node::Bin(_, a, b) | Node::Cmp(_, a, b) => {
            collect_attrs(a, out);
            collect_attrs(b, out);
        }
        Node::IfElse(a, b, c) => {
            collect_attrs(a, out);
            collect_attrs(b, out);
            collect_attrs(c, out);
        }
        Node::And(v) | Node::Or(v) | Node::Seq(v) | Node::Apply(_, v) => v.iter().for_each(|n| collect_attrs(n, out)),
        Node::Call(f, args, kw) => {
            collect_attrs(f, out);
            args.iter().for_each(|n| collect_attrs(n, out));
            kw.iter().for_each(|(_, n)| collect_attrs(n, out));
        }
        _ => {}
    }
}

/// What a `numpy.<name>` attribute is, as checked by the caller.
#[derive(Clone, Copy, Debug)]
pub enum NumpyName {
    /// The original ufunc of that name.
    Ufunc,
    /// A Python float.
    Float(f64),
    /// Anything else.
    Other,
}

impl NumpyText {
    /// The `numpy.<name>` attributes the text uses.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// The evaluation plan when `scalar[i]` tells which arguments are
    /// scalars (the others are arrays of one common shape), with the
    /// bindings of [`NumpyText::names`]; `libm` allows the C library
    /// functions (`sin`, `cos`, `exp`, `log`, `pow`).
    pub fn plan(&self, scalar: &[bool], bindings: &[NumpyName], libm: bool) -> Option<Plan> {
        if bindings.len() != self.names.len() {
            return None;
        }
        let mut b = PlanBuilder { names: &self.names, bindings, scalar, libm, plan: Plan::default() };
        let root = b.lower(&self.root)?;
        let mut plan = b.plan;
        plan.result = match root {
            Val::Py(v) => {
                plan.scalars.push(SIns::Const(v.v));
                Opd::S((plan.scalars.len() - 1) as u32)
            }
            Val::S(s) => Opd::S(s),
            Val::A(a) => a,
        };
        Some(plan)
    }
}

/// A Python literal value: a float or an exact integer below 2**53.
#[derive(Clone, Copy, Debug)]
struct V {
    v: f64,
    int: bool,
}

const EXACT: f64 = 9007199254740992.0;

/// An operand of an array operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Opd {
    /// Array argument.
    In(u32),
    /// Array register.
    R(u32),
    /// Scalar slot.
    S(u32),
}

#[derive(Clone, Copy, Debug)]
enum Val {
    Py(V),
    S(u32),
    A(Opd),
}

/// Scalar instructions, evaluated once per call with NumPy's float64
/// scalar semantics.
#[derive(Clone, Debug)]
enum SIns {
    Const(f64),
    Arg(u32),
    Bin(Bin, u32, u32),
    Un(Unary, u32),
}

/// Array instructions; register `i` is the result of instruction `i`.
#[derive(Clone, Debug)]
enum AIns {
    Bin(Bin, Opd, Opd),
    Un(Unary, Opd),
    /// `array ** scalar slot`: the NumPy fast paths depend on the value.
    PowS(Opd, u32),
}

/// An evaluation plan of a [`NumpyText`] for one argument pattern.
#[derive(Clone, Debug, Default)]
pub struct Plan {
    scalars: Vec<SIns>,
    arrays: Vec<AIns>,
    result: Opd,
}

impl Default for Opd {
    fn default() -> Opd {
        Opd::S(0)
    }
}

struct PlanBuilder<'a> {
    names: &'a [String],
    bindings: &'a [NumpyName],
    scalar: &'a [bool],
    libm: bool,
    plan: Plan,
}

fn py_float(v: f64) -> Option<V> {
    v.is_finite().then_some(V { v, int: false })
}

fn py_int(v: f64) -> Option<V> {
    (v.abs() <= EXACT).then_some(V { v: v + 0.0, int: true })
}

/// Python arithmetic on literals (the generated code computes them so).
fn py_bin(op: Bin, a: V, b: V) -> Option<V> {
    let both = a.int && b.int;
    match op {
        Bin::Add if both => py_int(a.v + b.v),
        Bin::Sub if both => py_int(a.v - b.v),
        Bin::Mul if both => py_int(a.v * b.v),
        Bin::Add => py_float(a.v + b.v),
        Bin::Sub => py_float(a.v - b.v),
        Bin::Mul => py_float(a.v * b.v),
        // ZeroDivisionError
        Bin::Div if b.v == 0.0 => None,
        Bin::Div => py_float(a.v / b.v),
        Bin::Pow if both && b.v >= 0.0 => py_int(unsafe { libm::pow(a.v, b.v) }),
        // 0 ** negative raises; a negative base with a fractional exponent
        // gives a Python complex
        Bin::Pow if a.v == 0.0 && b.v < 0.0 => None,
        Bin::Pow if a.v < 0.0 && b.v.fract() != 0.0 => None,
        Bin::Pow => py_float(unsafe { libm::pow(a.v, b.v) }),
    }
}

impl PlanBuilder<'_> {
    fn scalar(&mut self, i: SIns) -> u32 {
        self.plan.scalars.push(i);
        (self.plan.scalars.len() - 1) as u32
    }

    fn array(&mut self, i: AIns) -> Opd {
        self.plan.arrays.push(i);
        Opd::R((self.plan.arrays.len() - 1) as u32)
    }

    /// A scalar slot holding `v` (Python literals become float64 operands).
    fn slot(&mut self, v: Val) -> Option<u32> {
        match v {
            Val::Py(p) => Some(self.scalar(SIns::Const(p.v))),
            Val::S(s) => Some(s),
            Val::A(_) => None,
        }
    }

    fn opd(&mut self, v: Val) -> Opd {
        match v {
            Val::A(a) => a,
            other => Opd::S(self.slot(other).unwrap()),
        }
    }

    fn unary(&mut self, f: Unary, a: Val) -> Option<Val> {
        if f.libm() && !self.libm {
            return None;
        }
        Some(match a {
            // a NumPy function of a Python scalar is a numpy.float64
            Val::Py(_) | Val::S(_) => {
                let s = self.slot(a)?;
                Val::S(self.scalar(SIns::Un(f, s)))
            }
            Val::A(x) => Val::A(self.array(AIns::Un(f, x))),
        })
    }

    fn lower(&mut self, node: &Node) -> Option<Val> {
        Some(match node {
            Node::Num(v) => Val::Py(py_float(*v)?),
            Node::Int(v) => Val::Py(py_int(*v)?),
            Node::Arg(i) => {
                if *self.scalar.get(*i)? {
                    Val::S(self.scalar(SIns::Arg(*i as u32)))
                } else {
                    Val::A(Opd::In(*i as u32))
                }
            }
            Node::Attr(_, name) => {
                let k = self.names.iter().position(|n| n == name)?;
                match self.bindings[k] {
                    NumpyName::Float(v) => Val::Py(py_float(v)?),
                    _ => return None,
                }
            }
            Node::Pos(a) => {
                // unary plus copies (an array) or keeps the value
                self.lower(a)?
            }
            Node::Neg(a) => match self.lower(a)? {
                Val::Py(v) => Val::Py(if v.int { py_int(-v.v)? } else { V { v: -v.v, int: false } }),
                other => self.unary(Unary::Neg, other)?,
            },
            Node::Call(f, args, kwargs) => {
                let Node::Attr(_, name) = &**f else { return None };
                let k = self.names.iter().position(|n| n == name)?;
                if !matches!(self.bindings[k], NumpyName::Ufunc) || !kwargs.is_empty() || args.len() != 1 {
                    return None;
                }
                let u = Unary::from_name(name)?;
                let a = self.lower(&args[0])?;
                self.unary(u, a)?
            }
            Node::Bin(op, a, b) => {
                let (a, b) = (self.lower(a)?, self.lower(b)?);
                self.binary(*op, a, b)?
            }
            _ => return None,
        })
    }

    fn binary(&mut self, op: Bin, a: Val, b: Val) -> Option<Val> {
        if let (Val::Py(x), Val::Py(y)) = (a, b) {
            return Some(Val::Py(py_bin(op, x, y)?));
        }
        let pow = op == Bin::Pow;
        if pow && !self.libm {
            // only the exact fast paths of an array base
            let exact = matches!((a, b), (Val::A(_), Val::Py(e)) if e.v == 2.0 || e.v == -1.0 || e.v == 0.5);
            if !exact {
                return None;
            }
        }
        Some(match (a, b) {
            (Val::A(x), Val::Py(e)) if pow => {
                // ndarray.__pow__ fast paths for a Python scalar exponent
                if e.v == 2.0 {
                    Val::A(self.array(AIns::Un(Unary::Square, x)))
                } else if e.v == -1.0 {
                    Val::A(self.array(AIns::Un(Unary::Reciprocal, x)))
                } else if e.v == 0.5 {
                    Val::A(self.array(AIns::Un(Unary::Sqrt, x)))
                } else {
                    let s = self.scalar(SIns::Const(e.v));
                    Val::A(self.array(AIns::Bin(Bin::Pow, x, Opd::S(s))))
                }
            }
            (Val::A(x), Val::S(s)) if pow => Val::A(self.array(AIns::PowS(x, s))),
            (Val::A(_), _) | (_, Val::A(_)) => {
                let (x, y) = (self.opd(a), self.opd(b));
                Val::A(self.array(AIns::Bin(op, x, y)))
            }
            // NumPy scalar arithmetic (`**` is C pow)
            _ => {
                let (x, y) = (self.slot(a)?, self.slot(b)?);
                Val::S(self.scalar(SIns::Bin(op, x, y)))
            }
        })
    }
}

#[inline(always)]
fn bin(op: Bin, a: f64, b: f64) -> f64 {
    match op {
        Bin::Add => a + b,
        Bin::Sub => a - b,
        Bin::Mul => a * b,
        Bin::Div => a / b,
        Bin::Pow => unsafe { libm::pow(a, b) },
    }
}

#[derive(Clone, Copy)]
enum Src<'a> {
    V(&'a [f64]),
    C(f64),
}

#[inline(always)]
fn map2(dst: &mut [f64], a: Src<'_>, b: Src<'_>, f: impl Fn(f64, f64) -> f64) {
    match (a, b) {
        (Src::V(x), Src::V(y)) => {
            for ((d, x), y) in dst.iter_mut().zip(x).zip(y) {
                *d = f(*x, *y);
            }
        }
        (Src::V(x), Src::C(c)) => {
            for (d, x) in dst.iter_mut().zip(x) {
                *d = f(*x, c);
            }
        }
        (Src::C(c), Src::V(y)) => {
            for (d, y) in dst.iter_mut().zip(y) {
                *d = f(c, *y);
            }
        }
        (Src::C(c), Src::C(e)) => dst.fill(f(c, e)),
    }
}

#[inline(always)]
fn map1(dst: &mut [f64], a: Src<'_>, f: impl Fn(f64) -> f64) {
    match a {
        Src::V(x) => {
            for (d, x) in dst.iter_mut().zip(x) {
                *d = f(*x);
            }
        }
        Src::C(c) => dst.fill(f(c)),
    }
}

fn run_bin(dst: &mut [f64], op: Bin, a: Src<'_>, b: Src<'_>) {
    match op {
        Bin::Add => map2(dst, a, b, |x, y| x + y),
        Bin::Sub => map2(dst, a, b, |x, y| x - y),
        Bin::Mul => map2(dst, a, b, |x, y| x * y),
        Bin::Div => map2(dst, a, b, |x, y| x / y),
        Bin::Pow => map2(dst, a, b, |x, y| unsafe { libm::pow(x, y) }),
    }
}

fn run_un(dst: &mut [f64], f: Unary, a: Src<'_>) {
    match f {
        Unary::Neg => map1(dst, a, |x| -x),
        Unary::Sqrt => map1(dst, a, f64::sqrt),
        Unary::Square => map1(dst, a, |x| x * x),
        Unary::Reciprocal => map1(dst, a, |x| 1.0 / x),
        Unary::Sin => map1(dst, a, |x| unsafe { libm::sin(x) }),
        Unary::Cos => map1(dst, a, |x| unsafe { libm::cos(x) }),
        Unary::Exp => map1(dst, a, |x| unsafe { libm::exp(x) }),
        Unary::Log => map1(dst, a, |x| unsafe { libm::log(x) }),
    }
}

const CHUNK: usize = 512;

impl Plan {
    /// Whether the result is a scalar (it does not depend on the arrays).
    pub fn scalar_result(&self) -> bool {
        matches!(self.result, Opd::S(_))
    }

    /// Evaluate with `scalars[i]` the value of scalar argument `i` and
    /// `arrays[i]` the `n` values of array argument `i` (the slices of the
    /// other kind are ignored), writing the `n` results to `out`.
    pub fn eval(&self, scalars: &[f64], arrays: &[&[f64]], out: &mut [f64]) {
        let mut slots = Vec::with_capacity(self.scalars.len());
        for ins in &self.scalars {
            let v = match ins {
                SIns::Const(c) => *c,
                SIns::Arg(i) => scalars[*i as usize],
                SIns::Bin(op, a, b) => bin(*op, slots[*a as usize], slots[*b as usize]),
                SIns::Un(f, a) => f.apply(slots[*a as usize]),
            };
            slots.push(v);
        }
        let n = out.len();
        if let Opd::S(s) = self.result {
            out.fill(slots[s as usize]);
            return;
        }
        let width = CHUNK.min(n.max(1));
        let mut regs = vec![0.0f64; width * self.arrays.len()];
        let mut start = 0;
        while start < n {
            let len = width.min(n - start);
            for (i, ins) in self.arrays.iter().enumerate() {
                let (prev, cur) = regs.split_at_mut(i * width);
                let dst = &mut cur[..len];
                let src = |o: Opd| -> Src<'_> {
                    match o {
                        Opd::In(a) => Src::V(&arrays[a as usize][start..start + len]),
                        Opd::R(r) => Src::V(&prev[r as usize * width..r as usize * width + len]),
                        Opd::S(s) => Src::C(slots[s as usize]),
                    }
                };
                match ins {
                    AIns::Bin(op, a, b) => run_bin(dst, *op, src(*a), src(*b)),
                    AIns::Un(f, a) => run_un(dst, *f, src(*a)),
                    AIns::PowS(a, s) => {
                        let e = slots[*s as usize];
                        let f = if e == 2.0 {
                            Some(Unary::Square)
                        } else if e == -1.0 {
                            Some(Unary::Reciprocal)
                        } else if e == 0.5 {
                            Some(Unary::Sqrt)
                        } else {
                            None
                        };
                        match f {
                            Some(f) => run_un(dst, f, src(*a)),
                            None => run_bin(dst, Bin::Pow, src(*a), Src::C(e)),
                        }
                    }
                }
            }
            let o = &mut out[start..start + len];
            match self.result {
                Opd::In(a) => o.copy_from_slice(&arrays[a as usize][start..start + len]),
                Opd::R(r) => o.copy_from_slice(&regs[r as usize * width..r as usize * width + len]),
                Opd::S(_) => unreachable!(),
            }
            start += len;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plan(src: &str, scalar: &[bool], libm: bool) -> Option<Plan> {
        let t = parse_lambda(src)?;
        let b: Vec<NumpyName> = t
            .names()
            .iter()
            .map(|n| match n.as_str() {
                "pi" => NumpyName::Float(std::f64::consts::PI),
                "e" => NumpyName::Float(std::f64::consts::E),
                _ => NumpyName::Ufunc,
            })
            .collect();
        t.plan(scalar, &b, libm)
    }

    #[test]
    fn numpy_semantics() {
        let xs: Vec<f64> = (0..1500).map(|i| i as f64 * 0.01 - 3.0).collect();
        let p = plan("lambda a: ((a**2) + 3*a + 1)", &[false], true).unwrap();
        let mut out = vec![0.0; xs.len()];
        p.eval(&[0.0], &[&xs], &mut out);
        for (o, x) in out.iter().zip(&xs) {
            assert_eq!(o.to_bits(), ((x * x + 3.0 * x) + 1.0).to_bits());
        }
        let p = plan("lambda a, b: (numpy.exp(-(a**2))*numpy.cos(b) + a**3 + ((1/3))*numpy.pi)", &[false, true], true).unwrap();
        // the scalar of argument k is scalars[k]
        p.eval(&[0.0, 0.25], &[&xs, &[]], &mut out);
        for (o, x) in out.iter().zip(&xs) {
            let v = unsafe { libm::exp(-(x * x)) * libm::cos(0.25) + libm::pow(*x, 3.0) } + (1.0 / 3.0) * std::f64::consts::PI;
            assert_eq!(o.to_bits(), v.to_bits());
        }
        // a scalar argument squared is C pow; constants fill the result
        let p = plan("lambda a, b: (b**2)", &[false, true], true).unwrap();
        assert!(p.scalar_result());
        p.eval(&[0.0, 1.1], &[&xs, &[]], &mut out);
        assert!(out.iter().all(|o| o.to_bits() == unsafe { libm::pow(1.1, 2.0) }.to_bits()));
        // only exact operations without the C library
        assert!(plan("lambda a: (numpy.sin(a))", &[false], false).is_none());
        assert!(plan("lambda a: ((a**2)**((1/2)) + 1/a)", &[false], false).is_some());
        assert!(plan("lambda a: (a**3)", &[false], false).is_none());
        // Python errors and other names are declined
        assert!(plan("lambda a: (a + 1/0)", &[false], true).is_none());
        assert!(plan("lambda a: (numpy.tan(a))", &[false], true).is_none());
        assert!(plan("lambda a: (math.sin(a))", &[false], true).is_none());
        assert!(plan("lambda a: (a > 1)", &[false], true).is_none());
    }
}
