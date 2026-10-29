//! Derivatives of matrix expressions by a matrix symbol: SymPy 1.14's
//! matrix-line algorithm (`_LeftRightArgs`, the
//! `_eval_derivative_matrix_lines` methods and `_matrix_derivative`).
//!
//! An expression derived by `X[m, n]` carries two lines of matrix products:
//! the one contracted to `m` (first line) and the one contracted to `n`
//! (second line). A transposition swaps the lines; a trace joins their ends
//! into the scalar factor `higher`. Results of rank at most two are
//! contracted back to matrix expressions as SymPy does.
//!
//! Every construction mirrors the evaluated products of the reference
//! (`MatMul(a, b).doit()`, `Transpose(a).doit()`, `-expr.T`), so the forms
//! are SymPy's. The kernel declines (`None`) on every case it does not
//! cover: element-wise operations (SymPy's array-expression path), rank
//! four results (array expressions), powers of 1x1 matrices, nested traces,
//! symbolic matrix powers, scalar heads other than Add/Mul/Pow/Trace that
//! contain the variable, and any construction error. The caller then runs
//! its own algorithm unchanged.
use rustc_hash::FxHashSet;
use crate::expr::{Expr, Func, Kind};
use crate::ops::{add2, mul2, pow};
use super::{MatrixExpr, MatrixKind, MatrixOperand, expression_has_matrix};

/// A slot of a line: a scalar or a matrix expression.
#[derive(Clone, Debug, PartialEq)]
pub enum Value { Scalar(Expr), Matrix(MatrixExpr) }

struct Decline;
type R<T> = Result<T, Decline>;
impl From<crate::matrix::MatrixError> for Decline { fn from(_: crate::matrix::MatrixError) -> Self { Decline } }

impl Value {
    fn one() -> Self { Value::Scalar(Expr::one()) }
    fn zero() -> Self { Value::Scalar(Expr::zero()) }
    fn is_one(&self) -> bool { matches!(self, Value::Scalar(e) if e.is_one()) }
    fn is_matrix(&self) -> bool { matches!(self, Value::Matrix(_)) }
    /// `_is_zero`: a ZeroMatrix, or a scalar equal to 0.
    fn is_zero(&self) -> bool {
        match self { Value::Matrix(m) => matches!(m.kind(), MatrixKind::Zero), Value::Scalar(e) => e.is_zero() }
    }
    fn has(&self, x: &MatrixExpr) -> bool {
        match self { Value::Matrix(m) => m.has_matrix(x), Value::Scalar(e) => expression_has_matrix(e, x) }
    }
    fn is_identity_one(&self) -> bool {
        matches!(self, Value::Matrix(m) if matches!(m.kind(), MatrixKind::Identity) && m.rows().is_one())
    }
}

/// `MatMul(*values)` as the Python constructor builds it: the scalars
/// multiply into the coefficient, the matrices keep their order.
fn matmul(values: &[Value]) -> R<MatrixExpr> {
    let mut coefficient = Expr::one();
    let mut factors = Vec::new();
    for v in values {
        match v { Value::Scalar(e) => coefficient = mul2(&coefficient, e), Value::Matrix(m) => factors.push(m.clone()) }
    }
    Ok(MatrixExpr::mul(coefficient, &factors, false)?)
}

/// `_times(a, b)`: `a *= b` on a line slot.
fn times(a: &Value, b: &Value) -> R<Value> {
    if a.is_one() { return Ok(b.clone()); }
    if b.is_one() { return Ok(a.clone()); }
    if let (Value::Scalar(a), Value::Scalar(b)) = (a, b) { return Ok(Value::Scalar(mul2(a, b))); }
    Ok(Value::Matrix(matmul(&[a.clone(), b.clone()])?.doit(true)?))
}

/// `_product(factors)`: `MatMul(*factors).doit(deep=False)` or `Mul`.
fn product(factors: &[Value]) -> R<Value> {
    if factors.len() == 1 { return Ok(factors[0].clone()); }
    if factors.iter().any(Value::is_matrix) { return Ok(Value::Matrix(matmul(factors)?.doit(false)?)); }
    let scalars: Vec<Expr> = factors.iter().filter_map(|v| if let Value::Scalar(e) = v { Some(e.clone()) } else { None }).collect();
    Ok(Value::Scalar(crate::ops::mul(&scalars)))
}

fn sum(terms: Vec<Value>) -> R<Value> {
    let mut iter = terms.into_iter();
    let Some(mut result) = iter.next() else { return Ok(Value::zero()) };
    for t in iter {
        result = match (&result, &t) {
            (Value::Scalar(a), Value::Scalar(b)) => Value::Scalar(add2(a, b)),
            (Value::Matrix(a), Value::Matrix(b)) => Value::Matrix(MatrixExpr::add(&[a.clone(), b.clone()], true)?),
            // MatAdd's validate raises on a scalar term: the caller decides.
            _ => return Err(Decline),
        };
    }
    Ok(result)
}

#[derive(Clone)]
struct Line {
    lines: Vec<Value>,
    first: usize,
    second: usize,
    first_line: usize,
    second_line: usize,
    higher: Value,
}

impl Line {
    fn new(first: Value, second: Value, higher: Value) -> Self {
        Line { lines: vec![first, second], first: 0, second: 1, first_line: 0, second_line: 1, higher }
    }
    fn transpose(mut self) -> Self {
        std::mem::swap(&mut self.first, &mut self.second);
        std::mem::swap(&mut self.first_line, &mut self.second_line);
        self
    }
    fn append_first(&mut self, other: &Value) -> R<()> { self.lines[self.first] = times(&self.lines[self.first], other)?; Ok(()) }
    fn append_second(&mut self, other: &Value) -> R<()> { self.lines[self.second] = times(&self.lines[self.second], other)?; Ok(()) }
    fn build(&self) -> Vec<Value> {
        let mut data = self.lines.clone();
        if !self.higher.is_one() { data.push(self.higher.clone()); }
        data
    }
}

fn identity_or_one(n: &Expr) -> R<Value> {
    Ok(if n.is_one() { Value::one() } else { Value::Matrix(MatrixExpr::identity(n.clone())?) })
}

fn matrix_symbol_lines(e: &MatrixExpr, x: &MatrixExpr) -> R<Vec<Line>> {
    if e != x {
        let first = if e.rows().is_one() { Value::zero() } else { Value::Matrix(MatrixExpr::zero(x.rows().clone(), e.rows().clone())?) };
        let second = if e.cols().is_one() { Value::zero() } else { Value::Matrix(MatrixExpr::zero(x.cols().clone(), e.cols().clone())?) };
        return Ok(vec![Line::new(first, second, Value::one())]);
    }
    Ok(vec![Line::new(identity_or_one(e.rows())?, identity_or_one(e.cols())?, Value::one())])
}

/// `MatMul._eval_derivative_matrix_lines` over the Python arguments of the
/// product (the coefficient first when it is not one).
fn matmul_lines(args: &[Value], rows: &Expr, cols: &Expr, x: &MatrixExpr) -> R<Vec<Line>> {
    let mut out = Vec::new();
    for (ind, arg) in args.iter().enumerate() {
        if !arg.has(x) { continue; }
        let right = &args[ind + 1..];
        let right_mat = if right.is_empty() { Value::Matrix(MatrixExpr::identity(cols.clone())?) } else { product(right)? };
        let left = &args[..ind];
        let left_rev = if left.is_empty() {
            Value::Matrix(MatrixExpr::identity(rows.clone())?)
        } else {
            let reversed = left.iter().rev().map(|v| match v {
                Value::Matrix(m) => Ok(Value::Matrix(m.transpose(false)?.doit(true)?)),
                scalar => Ok(scalar.clone()),
            }).collect::<R<Vec<_>>>()?;
            product(&reversed)?
        };
        for mut line in lines(arg, x)? {
            line.append_first(&left_rev)?;
            line.append_second(&right_mat)?;
            out.push(line);
        }
    }
    Ok(out)
}

fn product_args(e: &MatrixExpr) -> Option<Vec<Value>> {
    let MatrixKind::Mul(coefficient, factors) = e.kind() else { return None };
    let mut args = Vec::with_capacity(factors.len() + 1);
    if !coefficient.is_one() { args.push(Value::Scalar(coefficient.clone())); }
    args.extend(factors.iter().cloned().map(Value::Matrix));
    Some(args)
}

/// `Trace._eval_derivative_matrix_lines` for lines without a scalar factor:
/// `ArrayContraction(ArrayTensorProduct(l0, l1), (1, 3))`, which the array
/// to matrix conversion gives as `l0*l1.T`.
fn trace_lines(m: &MatrixExpr, x: &MatrixExpr) -> R<Vec<Line>> {
    let mut out = lines(&Value::Matrix(m.clone()), x)?;
    for line in &mut out {
        if !line.higher.is_one() || line.lines.len() != 2 { return Err(Decline); }
        let (Value::Matrix(l0), Value::Matrix(l1)) = (&line.lines[0], &line.lines[1]) else { return Err(Decline) };
        let joined = MatrixExpr::mul(Expr::one(), &[l0.clone(), l1.transpose(true)?], true)?;
        line.higher = Value::Matrix(joined);
        line.lines = vec![Value::one(), Value::one()];
        line.first = 0;
        line.second = 1;
    }
    Ok(out)
}

fn inverse_lines(inverse: &MatrixExpr, base: &MatrixExpr, x: &MatrixExpr) -> R<Vec<Line>> {
    let mut out = lines(&Value::Matrix(base.clone()), x)?;
    let left = Value::Matrix(MatrixExpr::mul(Expr::minus_one(), &[inverse.transpose(true)?], true)?);
    let right = Value::Matrix(inverse.clone());
    for line in &mut out {
        line.append_first(&left)?;
        line.append_second(&right)?;
    }
    Ok(out)
}

fn matpow_lines(e: &MatrixExpr, base: &MatrixExpr, exponent: &Expr, x: &MatrixExpr) -> R<Vec<Line>> {
    // A 1x1 base is one contracted chain (`_power_1x1_lines`): not covered.
    if base.rows().is_one() && base.cols().is_one() { return Err(Decline); }
    let Some(n) = exponent.as_i64() else { return Err(Decline) };
    let factor = match n {
        0 => return lines(&Value::Matrix(e.doit(true)?), x),
        -1 => { let inverse = base.inverse(false)?; return inverse_lines(&inverse, base, x); },
        n if n >= 2 => base.clone(),
        n if n <= -2 => base.inverse(false)?,
        _ => return Err(Decline),
    };
    let count = n.unsigned_abs() as usize;
    if count > 64 { return Err(Decline); }
    let product = MatrixExpr::mul(Expr::one(), &vec![factor; count], false)?;
    let Some(args) = product_args(&product) else { return Err(Decline) };
    matmul_lines(&args, product.rows(), product.cols(), x)
}

fn lines(e: &Value, x: &MatrixExpr) -> R<Vec<Line>> {
    match e {
        Value::Matrix(m) => match m.kind() {
            MatrixKind::Symbol(_) => matrix_symbol_lines(m, x),
            MatrixKind::Add(terms) => {
                let mut out = Vec::new();
                for t in terms { out.extend(lines(&Value::Matrix(t.clone()), x)?); }
                Ok(out)
            },
            MatrixKind::Mul(..) => {
                let args = product_args(m).ok_or(Decline)?;
                matmul_lines(&args, m.rows(), m.cols(), x)
            },
            MatrixKind::Transpose(b) => Ok(lines(&Value::Matrix(b.clone()), x)?.into_iter().map(Line::transpose).collect()),
            MatrixKind::Inverse(b) => inverse_lines(m, b, x),
            MatrixKind::Pow(b, exponent) => matpow_lines(m, b, exponent, x),
            _ if !m.has_matrix(x) => Ok(vec![Line::new(Value::one(), Value::one(), Value::zero())]),
            _ => Err(Decline),
        },
        Value::Scalar(s) => match s.kind() {
            Kind::Func(Func::MatrixTrace(m), _) => trace_lines(m, x),
            // A power of a 1x1 matrix expression (`_power_1x1_lines`).
            Kind::Pow(b, _) if matches!(b.kind(), Kind::Func(Func::MatrixProxy(_), _)) => Err(Decline),
            _ => Ok(vec![Line::new(Value::one(), Value::one(), scalar_derivative(s, x)?)]),
        },
    }
}

/// `expr.diff(x)` for a scalar factor of a product rule.
fn derive(e: &Expr, x: &MatrixExpr) -> R<Value> {
    if !expression_has_matrix(e, x) { return Ok(Value::zero()); }
    derivative(&Value::Scalar(e.clone()), x)
}

/// `Expr._eval_derivative` for the scalar heads that hold matrices.
fn scalar_derivative(e: &Expr, x: &MatrixExpr) -> R<Value> {
    if !expression_has_matrix(e, x) { return Ok(Value::zero()); }
    match e.kind() {
        Kind::Add(..) => {
            let mut terms = Vec::new();
            for a in e.args() { let d = derive(&a, x)?; if !d.is_zero() { terms.push(d); } }
            sum(terms)
        },
        Kind::Mul(..) => {
            let args: Vec<Value> = e.args().into_iter().map(Value::Scalar).collect();
            let mut terms = Vec::new();
            for i in 0..args.len() {
                let Value::Scalar(a) = &args[i] else { unreachable!() };
                let d = derive(a, x)?;
                if d.is_zero() { continue; }
                let mut term = Value::one();
                for (j, f) in args.iter().enumerate() { term = times(&term, if j == i { &d } else { f })?; }
                terms.push(term);
            }
            sum(terms)
        },
        Kind::Pow(base, exponent) => {
            if expression_has_matrix(exponent, x) { return Err(Decline); }
            let dbase = derive(base, x)?;
            if dbase.is_zero() { return Ok(Value::zero()); }
            let inverse = pow(base, &Expr::minus_one());
            Ok(match dbase {
                Value::Scalar(d) => Value::Scalar(mul2(e, &mul2(&mul2(&d, exponent), &inverse))),
                Value::Matrix(d) => {
                    let scaled = MatrixExpr::mul(exponent.clone(), &[d], true)?;
                    let divided = MatrixExpr::mul(inverse, &[scaled], true)?;
                    Value::Matrix(MatrixExpr::mul(e.clone(), &[divided], true)?)
                },
            })
        },
        _ => Err(Decline),
    }
}

fn shape_of(v: &Value) -> [Expr; 2] {
    match v { Value::Matrix(m) => [m.rows().clone(), m.cols().clone()], Value::Scalar(_) => [Expr::one(), Expr::one()] }
}

fn contract_one_dims(mut parts: Vec<Value>) -> R<Value> {
    if parts.len() == 1 { return Ok(parts.pop().unwrap()); }
    let (p1, mut p2) = (parts[0].clone(), parts[1].clone());
    if parts.len() == 2 {
        if let (Value::Matrix(a), Value::Matrix(b)) = (&p1, &p2) {
            if a.cols() != b.cols() {
                // A vector derived by a vector.
                if b.rows().is_one() && b.cols().is_one() {
                    return if p2.is_identity_one() { Ok(p1) } else { Err(Decline) };
                }
                if a.cols() == b.rows() { return times(&p1, &p2); }
                return Err(Decline);
            }
        }
    }
    if let Value::Matrix(m) = &p2 { p2 = Value::Matrix(m.transpose(true)?); }
    let pbase = if p1.is_identity_one() { p2 } else if p2.is_identity_one() { p1 } else { times(&p1, &p2)? };
    if parts.len() == 2 { return Ok(pbase); }
    if pbase.is_matrix() { return Err(Decline); }
    let mut rest = Value::one();
    for p in &parts[2..] { rest = times(&rest, p)?; }
    times(&pbase, &rest)
}

/// Whether an element-wise matrix operation occurs anywhere in `v`.
fn has_elementwise(v: &Value) -> bool {
    enum Node { M(MatrixExpr), S(Expr) }
    let mut pending = vec![match v { Value::Matrix(m) => Node::M(m.clone()), Value::Scalar(e) => Node::S(e.clone()) }];
    let mut seen_m = FxHashSet::default();
    let mut seen_s = FxHashSet::default();
    while let Some(node) = pending.pop() {
        match node {
            Node::M(m) => {
                if !seen_m.insert(m.dag_id()) { continue; }
                match m.kind() {
                    MatrixKind::HadamardProduct(_) | MatrixKind::HadamardPower(..) | MatrixKind::ApplyFunction(..) => return true,
                    MatrixKind::Mul(c, _) => pending.push(Node::S(c.clone())),
                    MatrixKind::Pow(_, e) => pending.push(Node::S(e.clone())),
                    MatrixKind::Explicit(data) => pending.extend(data.data().iter().cloned().map(Node::S)),
                    _ => {},
                }
                if let MatrixKind::HadamardPower(a, b) = m.kind() {
                    for o in [a, b] { if let MatrixOperand::Scalar(e) = o { pending.push(Node::S(e.clone())); } }
                }
                pending.extend(m.children().into_iter().map(Node::M));
            },
            Node::S(e) => {
                if !seen_s.insert(e.dag_id()) { continue; }
                if let Kind::Func(f, _) = e.kind() {
                    match f {
                        Func::MatrixElement(p) | Func::MatrixProxy(p) | Func::MatrixTrace(p) | Func::MatrixDeterminant(p) => pending.push(Node::M(p.as_ref().clone())),
                        Func::MatrixDotProduct(a, b) => { pending.push(Node::M(a.as_ref().clone())); pending.push(Node::M(b.as_ref().clone())); },
                        _ => {},
                    }
                }
                pending.extend(e.args().into_iter().map(Node::S));
            },
        }
    }
    false
}

fn derivative(e: &Value, x: &MatrixExpr) -> R<Value> {
    if let Value::Scalar(s) = e {
        if !expression_has_matrix(s, x) { return Ok(Value::Matrix(MatrixExpr::zero(x.rows().clone(), x.cols().clone())?)); }
    }
    if has_elementwise(e) { return Err(Decline); }
    let lines = lines(e, x)?;
    let parts: Vec<Vec<Value>> = lines.iter().map(Line::build).collect();
    let rank = |part: &[Value]| part.iter().flat_map(shape_of).filter(|n| !n.is_one()).count();
    if parts.is_empty() || rank(&parts[0]) > 2 { return Err(Decline); }
    sum(parts.into_iter().map(contract_one_dims).collect::<R<Vec<_>>>()?)
}

/// `expr.diff(x)` for a matrix expression or a scalar and a matrix symbol
/// `x`, when the result has rank at most two; `None` when the kernel does
/// not cover the expression.
pub fn matrix_lines_derivative(e: &Value, x: &MatrixExpr) -> Option<Value> {
    if !matches!(x.kind(), MatrixKind::Symbol(_)) { return None; }
    // A matrix by a matrix with no unit dimension is a rank-four array: only
    // dimensions equal to 1 are contracted, so decline before building lines.
    if let Value::Matrix(m) = e {
        if ![m.rows(), m.cols(), x.rows(), x.cols()].iter().any(|d| d.is_one()) { return None; }
    }
    derivative(e, x).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    fn sym(name: &str, r: i64, c: i64) -> MatrixExpr { MatrixExpr::symbol(name, Expr::int(r), Expr::int(c)).unwrap() }
    fn m(v: Value) -> MatrixExpr { match v { Value::Matrix(m) => m, _ => panic!("scalar") } }

    #[test]
    fn product_and_trace_forms() {
        let x = sym("X", 3, 3);
        let a = sym("A", 3, 3);
        // (A*X).diff(X) is rank four: declined.
        let ax = MatrixExpr::mul(Expr::one(), &[a.clone(), x.clone()], true).unwrap();
        assert!(matrix_lines_derivative(&Value::Matrix(ax.clone()), &x).is_none());
        // Trace(A*X).diff(X) = A.T
        let t = ax.trace().unwrap();
        let d = m(matrix_lines_derivative(&Value::Scalar(t), &x).unwrap());
        assert_eq!(d, a.transpose(true).unwrap());
        // Trace(X).diff(X) = I
        let d = m(matrix_lines_derivative(&Value::Scalar(x.trace().unwrap()), &x).unwrap());
        assert_eq!(d, MatrixExpr::identity(Expr::int(3)).unwrap());
    }

    #[test]
    fn vector_forms() {
        let x = sym("x", 3, 1);
        // (x.T*x).diff(x) = 2*x
        let xtx = MatrixExpr::mul(Expr::one(), &[x.transpose(true).unwrap(), x.clone()], true).unwrap();
        let d = m(matrix_lines_derivative(&Value::Matrix(xtx), &x).unwrap());
        assert_eq!(d, MatrixExpr::mul(Expr::int(2), &[x.clone()], true).unwrap());
    }
}
