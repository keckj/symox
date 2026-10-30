//! SymPy 1.14 `convert_matrix_to_array`
//! (`sympy/tensor/array/expressions/from_matrix_to_array.py`) over the
//! native matrix expressions.
//!
//! The conversion builds the tree of SymPy's canonicalizing constructors
//! (`_array_tensor_product`, `_array_contraction`, `_array_add`,
//! `_permute_dims`, `_array_diagonal`, plus the raw `Reshape` and
//! `ArrayElementwiseApplyFunc`) and evaluates it bottom-up with the
//! canonical forms of [`crate::arrayexpr`], which is the order in which
//! SymPy applies them. Matrices and scalars SymPy keeps as they are become
//! leaves.
use crate::arrayexpr::{AErr, Ctx, LeafInfo, Node, R};
use crate::expr::{Expr, Func, Kind};
use crate::matrixexpr::{MatrixExpr, MatrixKind, MatrixOperand};

#[derive(Clone, Debug)]
pub enum Source { Matrix(MatrixExpr), Scalar(Expr) }

#[derive(Default)]
pub struct Conversion {
    pub leaves: Vec<LeafInfo>,
    pub sources: Vec<Source>,
    pub dims: Vec<Expr>,
}

/// `n` copies of `base`, the `[base for i in range(n)]` of SymPy; an
/// allocation that cannot be made is a MemoryError, as in Python.
fn repeat(base: Node, n: usize) -> R<Vec<Node>> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| AErr::Memory(format!("cannot build a product of {n} factors")))?;
    v.resize(n, base);
    Ok(v)
}

/// A positive integer exponent as a factor count; one that does not fit a
/// machine word cannot be expanded (SymPy's `range` loop exhausts memory).
fn count(exp: &Expr) -> R<usize> {
    exp.as_i64().and_then(|n| usize::try_from(n).ok())
        .ok_or_else(|| AErr::Memory(format!("cannot build a product of {exp} factors")))
}

impl Conversion {
    fn dim(&mut self, d: &Expr) -> usize {
        if let Some(i) = self.dims.iter().position(|x| x == d) { return i; }
        self.dims.push(d.clone());
        self.dims.len() - 1
    }

    fn matrix_leaf(&mut self, m: &MatrixExpr) -> Node {
        let shape = vec![self.dim(m.rows()), self.dim(m.cols())];
        self.leaves.push(LeafInfo { rank: 2, shape: Some(shape), zero: matches!(m.kind(), MatrixKind::Zero) });
        self.sources.push(Source::Matrix(m.clone()));
        Node::Leaf(self.leaves.len() - 1)
    }

    fn scalar_leaf(&mut self, e: &Expr) -> Node {
        self.leaves.push(LeafInfo { rank: 0, shape: Some(vec![]), zero: false });
        self.sources.push(Source::Scalar(e.clone()));
        Node::Leaf(self.leaves.len() - 1)
    }

    fn rank(&self, n: &Node) -> usize {
        let ones: Vec<bool> = vec![];
        let mut order = |_: &[&Node]| Ok(vec![]);
        Ctx { leaves: &self.leaves, ones: &ones, order: &mut order }.rank(n)
    }

    fn power(exp: &Expr) -> R<Option<usize>> {
        if crate::assumptions::is_positive(exp) != Some(true) { return Ok(None); }
        if !exp.is_integer_number() {
            return Err(AErr::Type("object cannot be interpreted as an integer".into()));
        }
        Ok(Some(count(exp)?))
    }

    /// The raw constructor tree of `convert_matrix_to_array(expr)` for a scalar.
    pub fn scalar(&mut self, e: &Expr) -> R<Node> {
        match e.kind() {
            Kind::Func(Func::MatrixTrace(m), _) => {
                let inner = self.matrix(m)?;
                let r = self.rank(&inner);
                Ok(Node::Contract(Box::new(inner), vec![vec![0, r.saturating_sub(1)]]))
            },
            Kind::Mul(..) => Ok(Node::Tensor(e.args().iter().map(|a| self.scalar(a)).collect::<R<_>>()?)),
            Kind::Pow(b, x) => {
                let base = self.scalar(b)?;
                match Self::power(x)? {
                    Some(n) => Ok(Node::Tensor(repeat(base, n)?)),
                    None => Ok(self.scalar_leaf(e)),
                }
            },
            _ => Ok(self.scalar_leaf(e)),
        }
    }

    fn chain(&mut self, scalar: Option<Node>, mats: Vec<Node>) -> Node {
        let groups: Vec<Vec<usize>> = (0..mats.len().saturating_sub(1)).map(|i| vec![2 * i + 1, 2 * i + 2]).collect();
        let mut args = vec![];
        args.extend(scalar);
        args.extend(mats);
        Node::Contract(Box::new(Node::Tensor(args)), groups)
    }

    /// The raw constructor tree of `convert_matrix_to_array(expr)` for a matrix.
    pub fn matrix(&mut self, m: &MatrixExpr) -> R<Node> {
        match m.kind() {
            MatrixKind::Mul(c, factors) => {
                let scalar = if c.is_one() { None } else { Some(Node::Tensor(vec![self.scalar(c)?])) };
                let mats = factors.iter().map(|f| self.matrix(f)).collect::<R<Vec<_>>>()?;
                Ok(self.chain(scalar, mats))
            },
            MatrixKind::Add(a) => Ok(Node::Add(a.iter().map(|x| self.matrix(x)).collect::<R<_>>()?)),
            MatrixKind::Transpose(x) => Ok(Node::Permute(Box::new(self.matrix(x)?), vec![1, 0])),
            MatrixKind::Pow(b, x) => {
                let base = self.matrix(b)?;
                if !x.is_integer_number() {
                    let exponent = self.scalar_leaf(x);
                    let Node::Leaf(i) = exponent else { unreachable!() };
                    return Ok(Node::PowApply(Box::new(base), i));
                }
                match Self::power(x)? {
                    Some(n) => {
                        // `MatMul.fromiter(base, ...)`: SymPy fails on a base
                        // that converted to an array expression
                        if !matches!(base, Node::Leaf(_)) {
                            let class = match &base {
                                Node::Contract(..) => "ArrayContraction",
                                Node::Permute(..) => "PermuteDims",
                                Node::Tensor(..) => "ArrayTensorProduct",
                                Node::Add(..) => "ArrayAdd",
                                Node::Diagonal(..) => "ArrayDiagonal",
                                Node::Reshape(..) => "Reshape",
                                Node::PowApply(..) => "ArrayElementwiseApplyFunc",
                                Node::Zero(..) => "ZeroArray",
                                Node::Leaf(_) => "MatrixExpr",
                            };
                            return Err(AErr::Attribute(format!("'{class}' object has no attribute 'as_coeff_matrices'")));
                        }
                        Ok(self.chain(None, repeat(base, n)?))
                    },
                    None => Ok(self.matrix_leaf(m)),
                }
            },
            MatrixKind::HadamardProduct(a) => Ok(Self::hadamard(a.iter().map(|x| self.matrix(x)).collect::<R<_>>()?)),
            MatrixKind::HadamardPower(MatrixOperand::Matrix(b), MatrixOperand::Scalar(x)) => {
                if x.is_integer_number() && crate::assumptions::is_positive(x) == Some(true) {
                    let n = count(x)?;
                    let base = self.matrix(b)?;
                    return Ok(Self::hadamard(repeat(base, n)?));
                }
                let base = self.matrix_leaf(b);
                let exponent = self.scalar_leaf(x);
                let Node::Leaf(i) = exponent else { unreachable!() };
                Ok(Node::PowApply(Box::new(base), i))
            },
            MatrixKind::Kronecker(a) => {
                let args = a.iter().map(|x| self.matrix(x)).collect::<R<Vec<_>>>()?;
                let n = args.len();
                let perm: Vec<usize> = (0..n).map(|i| 2 * i).chain((0..n).map(|i| 2 * i + 1)).collect();
                let shape = vec![self.dim(m.rows()), self.dim(m.cols())];
                Ok(Node::Reshape(Box::new(Node::Permute(Box::new(Node::Tensor(args)), perm)), shape))
            },
            _ => Ok(self.matrix_leaf(m)),
        }
    }

    fn hadamard(args: Vec<Node>) -> Node {
        let n = args.len();
        let groups = vec![(0..n).map(|i| 2 * i).collect(), (0..n).map(|i| 2 * i + 1).collect()];
        Node::Diagonal(Box::new(Node::Tensor(args)), groups)
    }

    /// Evaluates the constructor tree with SymPy's canonical forms.
    pub fn evaluate(&self, tree: Node, order: &mut dyn FnMut(&[&Node]) -> R<Vec<usize>>) -> R<Node> {
        let ones: Vec<bool> = self.dims.iter().map(Expr::is_one).collect();
        Ctx { leaves: &self.leaves, ones: &ones, order }.doit(tree, true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matmul_chain_and_transpose() {
        let k = Expr::symbol("k");
        let m = MatrixExpr::symbol("M", k.clone(), k.clone()).unwrap();
        let n = MatrixExpr::symbol("N", k.clone(), k.clone()).unwrap();
        let product = MatrixExpr::mul(Expr::one(), &[m.clone(), n.clone()], true).unwrap();
        let mut c = Conversion::default();
        let tree = c.matrix(&product).unwrap();
        let mut order = |n: &[&Node]| Ok((0..n.len()).collect());
        let out = c.evaluate(tree, &mut order).unwrap();
        assert_eq!(out, Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(0), Node::Leaf(1)])), vec![vec![1, 2]]));
        let mut c = Conversion::default();
        let t = c.matrix(&product.transpose(false).unwrap()).unwrap();
        let out = c.evaluate(t, &mut order).unwrap();
        // (M*N).T: the permutation is absorbed by swapping the factors
        assert_eq!(out, Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(1), Node::Leaf(0)])), vec![vec![0, 3]]));
    }
}
