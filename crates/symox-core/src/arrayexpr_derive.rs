//! SymPy 1.14 `array_derive`
//! (`sympy/tensor/array/expressions/arrayexpr_derivatives.py`).
//!
//! The expression is given as a tree of `DNode`s: the structural array
//! nodes (tensor products, sums, permutations, contractions, diagonals,
//! reshapes) and one node per `singledispatch` rule of the reference
//! (array or matrix symbol, constant, transpose, inverse, element-wise
//! function, converted matrix expression, failure). Every object the rules
//! build (`Identity(i)`, `-expr`, `ElementwiseApplyFunction(fdiff, expr)`,
//! ...) is a leaf the caller made beforehand; the derivative itself, i.e.
//! every `_array_tensor_product`, `_permute_dims`, `_array_contraction`,
//! `_array_diagonal` and `_array_add` call of the reference with its index
//! arithmetic and its checks, runs here on the canonical forms of
//! `arrayexpr.rs`.
//!
//! Errors the reference raises from Python objects (a missing `shape`
//! attribute, a failing `fdiff`, `NotImplementedError` for unsupported
//! types, the recursion of an unconvertible matrix expression) are
//! captured by the caller and stored as indices (`DErr::Py`); they are
//! raised at the point where the reference evaluates them.
use crate::arrayexpr::{invert, AErr, Ctx, LeafInfo, Node};

/// A derivation node: one `array_derive` dispatch case.
#[derive(Clone, Debug)]
pub enum DNode {
    /// `ArrayTensorProduct(*args)`: the derivation nodes of the arguments
    /// and their array forms (the factors kept in each addend).
    Tensor(Vec<DNode>, Vec<Node>),
    Add(Vec<DNode>),
    Permute(Box<DNode>, Vec<usize>),
    Contract(Box<DNode>, Vec<Vec<usize>>),
    Diagonal(Box<DNode>, Vec<Vec<usize>>),
    /// `Reshape(expr, shape)` (shape as dimension ids).
    Reshape(Box<DNode>, Vec<usize>),
    /// The `Expr` rule: `ZeroArray(*x.shape)`.
    Scalar,
    /// `ZeroArray(*(x.shape + expr.shape))`: an array or matrix symbol other
    /// than `x`, `Identity`, `OneMatrix` (dimension ids of `expr.shape`).
    ZeroCat(Vec<usize>),
    /// `expr == x` for an ArraySymbol: the leaves `Identity(i)` of its shape.
    SelfArray(Vec<usize>),
    /// `expr == x` for a MatrixSymbol: the leaves `Identity(m)`, `Identity(n)`.
    SelfMatrix(usize, usize),
    /// `Transpose(arg)`.
    Transpose(Box<DNode>),
    /// `Inverse`: the derivation node of `expr.I`, the leaves `-expr` and `expr`.
    Inverse(Box<DNode>, usize, usize),
    /// `ElementwiseApplyFunction`: a failing rank assertion, then the leaf
    /// `ElementwiseApplyFunction(fdiff, expr.expr)` (or the error of
    /// `_get_function_fdiff`), then the node of `expr.expr`.
    EwMatrix { assert_err: Option<usize>, applied: Result<usize, usize>, sub: Box<DNode> },
    /// `ArrayElementwiseApplyFunc`: the leaf `ArrayElementwiseApplyFunc(fdiff,
    /// expr.expr)` (or the error of the fdiff), the node of `expr.expr` and
    /// `get_rank(expr)` (or its error).
    EwArray { applied: Result<usize, usize>, sub: Box<DNode>, rank: Result<usize, usize> },
    /// A matrix expression, derived through `convert_matrix_to_array(expr)`.
    Converted(Box<DNode>),
    /// A dispatch case that raises (the stored Python exception).
    Fail(usize),
}

/// `x`: `x.shape` (or the AttributeError of its absence), `get_shape(x)`
/// (or its error) and `get_rank(x)` (or its error).
#[derive(Clone, Debug)]
pub struct XInfo {
    pub shape: Result<Vec<usize>, usize>,
    pub get_shape: Result<Vec<usize>, usize>,
    pub rank: Result<i64, usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DErr {
    A(AErr),
    /// The stored Python exception of this index.
    Py(usize),
    /// `AttributeError`: the leaf object has no `shape` attribute.
    NoShape(usize),
}

impl From<AErr> for DErr {
    fn from(e: AErr) -> Self { DErr::A(e) }
}

pub type DR<T> = Result<T, DErr>;

pub struct Derive<'a, 'b> {
    pub ctx: Ctx<'b>,
    /// Leaves that have a `shape` attribute (others are scalars).
    pub has_shape: &'a [bool],
    /// The leaf `S.Zero` (what `ZeroArray()` of an empty shape returns).
    pub szero: usize,
    pub x: &'a XInfo,
}

impl Derive<'_, '_> {
    /// `ZeroArray()` with no shape is `S.Zero` in the reference.
    fn fix(&self, n: Node) -> Node {
        match n {
            Node::Zero(s) if s.is_empty() => Node::Leaf(self.szero),
            o => o,
        }
    }

    fn x_shape(&self) -> DR<Vec<usize>> { self.x.shape.clone().map_err(DErr::Py) }

    fn x_get_shape(&self) -> DR<Vec<usize>> { self.x.get_shape.clone().map_err(DErr::Py) }

    fn get_shape(&self, n: &Node) -> DR<Vec<usize>> {
        self.ctx.shape(n).ok_or_else(|| DErr::A(AErr::Type("can only concatenate tuple (not \"NoneType\") to tuple".into())))
    }

    /// `get_rank` of an array node.
    fn get_rank(&self, n: &Node) -> DR<i64> {
        match n {
            Node::Leaf(i) => Ok(self.ctx.leaves[*i].rank as i64),
            _ => self.ctx.shape(n).map(|s| s.len() as i64).ok_or_else(|| DErr::A(AErr::Type("object of type 'NoneType' has no len()".into()))),
        }
    }

    /// `_array_add(*args)`: `ArrayAdd.__new__` checks the ranks, then reads
    /// `arg.shape` of every argument and compares the known shapes.
    fn add(&mut self, args: Vec<Node>) -> DR<Node> {
        let mut ranks = vec![];
        for a in &args {
            let r = self.get_rank(a)?;
            if !ranks.contains(&r) { ranks.push(r); }
        }
        if ranks.len() != 1 { return Err(AErr::Value("summing arrays of different ranks".into()).into()); }
        let mut shapes: Vec<Vec<usize>> = vec![];
        for a in &args {
            if let Node::Leaf(i) = a {
                if !self.has_shape[*i] { return Err(DErr::NoShape(*i)); }
            }
            if let Some(s) = self.ctx.shape(a) {
                if !shapes.contains(&s) { shapes.push(s); }
            }
        }
        if shapes.len() > 1 { return Err(AErr::Value("mismatching shapes in addition".into()).into()); }
        let out = self.ctx.add(args)?;
        Ok(self.fix(out))
    }

    fn permute(&mut self, n: Node, perm: Vec<usize>) -> DR<Node> {
        let out = self.ctx.permute(n, perm)?;
        Ok(self.fix(out))
    }

    fn tensor(&mut self, args: Vec<Node>) -> DR<Node> {
        let out = self.ctx.tensor(args)?;
        Ok(self.fix(out))
    }

    fn contract(&mut self, n: Node, groups: Vec<Vec<usize>>) -> DR<Node> {
        let out = self.ctx.contract(n, groups)?;
        Ok(self.fix(out))
    }

    fn diagonal(&mut self, n: Node, groups: Vec<Vec<usize>>) -> DR<Node> {
        let out = self.ctx.diagonal(n, groups, false)?;
        Ok(self.fix(out))
    }

    fn zeros(&self, dims: Vec<usize>) -> Node { self.fix(Node::Zero(dims)) }

    /// `array_derive(expr, x)`.
    pub fn derive(&mut self, e: &DNode) -> DR<Node> {
        match e {
            DNode::Scalar => Ok(self.zeros(self.x_shape()?)),
            DNode::Fail(i) => Err(DErr::Py(*i)),
            DNode::ZeroCat(dims) => {
                let mut s = self.x_shape()?;
                s.extend(dims.iter().copied());
                Ok(self.zeros(s))
            },
            DNode::SelfArray(ids) => {
                let r = ids.len();
                let perm: Vec<usize> = (0..r).map(|i| 2 * i).chain((0..r).map(|i| 2 * i + 1)).collect();
                // ArrayTensorProduct.fromiter(...) is not canonicalized.
                let raw = Node::Tensor(ids.iter().map(|&i| Node::Leaf(i)).collect());
                self.permute(raw, perm)
            },
            DNode::SelfMatrix(m, n) => {
                let t = self.tensor(vec![Node::Leaf(*m), Node::Leaf(*n)])?;
                self.permute(t, vec![0, 2, 1, 3])
            },
            DNode::Transpose(arg) => {
                let fd = self.derive(arg)?;
                self.permute(fd, vec![0, 1, 3, 2])
            },
            DNode::Inverse(mat, neg, me) => {
                let dexpr = self.derive(mat)?;
                let tp = self.tensor(vec![Node::Leaf(*neg), dexpr, Node::Leaf(*me)])?;
                let mp = self.contract(tp, vec![vec![1, 4], vec![5, 6]])?;
                self.permute(mp, vec![1, 2, 0, 3])
            },
            DNode::EwMatrix { assert_err, applied, sub } => {
                if let Some(i) = assert_err { return Err(DErr::Py(*i)); }
                let leaf = applied.map_err(DErr::Py)?;
                let dexpr = self.derive(sub)?;
                let tp = self.tensor(vec![Node::Leaf(leaf), dexpr])?;
                self.diagonal(tp, vec![vec![0, 4], vec![1, 5]])
            },
            DNode::EwArray { applied, sub, rank } => {
                let leaf = applied.map_err(DErr::Py)?;
                let dsub = self.derive(sub)?;
                let tp = self.tensor(vec![dsub, Node::Leaf(leaf)])?;
                let b = self.x.rank.map_err(DErr::Py)?;
                let c = rank.map_err(DErr::Py)?;
                if b < 0 {
                    return Err(AErr::NotImplemented("array_derive of an elementwise function with respect to an object of undefined rank".into()).into());
                }
                let b = b as usize;
                let groups = (0..c).map(|i| vec![b + i, b + c + i]).collect();
                self.diagonal(tp, groups)
            },
            DNode::Converted(cg) => self.derive(cg),
            DNode::Contract(sub, groups) => {
                let fd = self.derive(sub)?;
                let rx = self.x_get_shape()?.len();
                self.contract(fd, groups.iter().map(|g| g.iter().map(|j| j + rx).collect()).collect())
            },
            DNode::Diagonal(sub, groups) => {
                let ds = self.derive(sub)?;
                let rx = self.x_get_shape()?.len();
                self.diagonal(ds, groups.iter().map(|g| g.iter().map(|j| j + rx).collect()).collect())
            },
            DNode::Add(args) => {
                let mut out = Vec::with_capacity(args.len());
                for a in args { out.push(self.derive(a)?); }
                self.add(out)
            },
            DNode::Permute(sub, perm) => {
                let de = self.derive(sub)?;
                let p = [0, 1].into_iter().chain(perm.iter().map(|i| i + 2)).collect();
                self.permute(de, p)
            },
            DNode::Reshape(sub, dims) => {
                let de = self.derive(sub)?;
                let mut s = self.x_get_shape()?;
                s.extend(dims.iter().copied());
                Ok(Node::Reshape(Box::new(de), s))
            },
            DNode::Tensor(dargs, args) => self.derive_tensor(dargs, args),
        }
    }

    fn derive_tensor(&mut self, dargs: &[DNode], args: &[Node]) -> DR<Node> {
        let mut addends = vec![];
        for (i, d) in dargs.iter().enumerate() {
            let darg = self.derive(d)?;
            if darg == Node::Leaf(self.szero) { continue; }
            let mut prev = 0;
            for a in &args[..i] { prev += self.get_shape(a)?.len(); }
            let mut succ = 0;
            for a in &args[i + 1..] { succ += self.get_shape(a)?.len(); }
            let mut factors: Vec<Node> = args[..i].to_vec();
            factors.push(darg);
            factors.extend(args[i + 1..].iter().cloned());
            let addend = self.tensor(factors)?;
            let tot1 = self.x_get_shape()?.len();
            let tot2 = tot1 + prev;
            let tot3 = tot2 + self.get_shape(&args[i])?.len();
            let tot4 = tot3 + succ;
            let perm: Vec<usize> = (tot1..tot2).chain(0..tot1).chain(tot2..tot3).chain(tot3..tot4).collect();
            addends.push(self.permute(addend, invert(&perm))?);
        }
        match addends.len() {
            0 => Ok(Node::Leaf(self.szero)),
            1 => Ok(addends.pop().unwrap()),
            _ => self.add(addends),
        }
    }
}

/// Runs `array_derive` with a fresh canonical-form context.
pub fn array_derive(
    e: &DNode, leaves: &[LeafInfo], has_shape: &[bool], ones: &[bool], szero: usize, x: &XInfo,
    order: &mut dyn FnMut(&[&Node]) -> Result<Vec<usize>, AErr>,
) -> DR<Node> {
    let ctx = Ctx { leaves, ones, order };
    let mut d = Derive { ctx, has_shape, szero, x };
    d.derive(e)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(e: &DNode, leaves: &[LeafInfo], x: &XInfo) -> DR<Node> {
        let has_shape: Vec<bool> = leaves.iter().map(|l| l.rank > 0).collect();
        let ones = vec![false; 4];
        let mut order = |n: &[&Node]| Ok((0..n.len()).collect());
        array_derive(e, leaves, &has_shape, &ones, 0, x, &mut order)
    }

    #[test]
    fn matrix_symbol_and_transpose() {
        // leaves: 0 S.Zero, 1 Identity(k), 2 Identity(k)
        let leaves = vec![
            LeafInfo { rank: 0, shape: Some(vec![]), zero: false },
            LeafInfo { rank: 2, shape: Some(vec![0, 0]), zero: false },
            LeafInfo { rank: 2, shape: Some(vec![0, 0]), zero: false },
        ];
        let x = XInfo { shape: Ok(vec![0, 0]), get_shape: Ok(vec![0, 0]), rank: Ok(2) };
        let out = run(&DNode::SelfMatrix(1, 2), &leaves, &x).unwrap();
        assert_eq!(out, Node::Permute(Box::new(Node::Tensor(vec![Node::Leaf(1), Node::Leaf(2)])), vec![0, 2, 1, 3]));
        let out = run(&DNode::Transpose(Box::new(DNode::SelfMatrix(1, 2))), &leaves, &x).unwrap();
        // SymPy 1.14: array_derive(X.T, X) is PermuteDims(ArrayTensorProduct(I, I), (1 2 3)), array form [0, 2, 3, 1].
        assert_eq!(out, Node::Permute(Box::new(Node::Tensor(vec![Node::Leaf(1), Node::Leaf(2)])), vec![0, 2, 3, 1]));
        assert_eq!(run(&DNode::ZeroCat(vec![1]), &leaves, &x).unwrap(), Node::Zero(vec![0, 0, 1]));
    }

    #[test]
    fn scalar_x_gives_szero() {
        let leaves = vec![LeafInfo { rank: 0, shape: Some(vec![]), zero: false }];
        let x = XInfo { shape: Ok(vec![]), get_shape: Ok(vec![]), rank: Ok(0) };
        assert_eq!(run(&DNode::Scalar, &leaves, &x).unwrap(), Node::Leaf(0));
        let t = DNode::Tensor(vec![DNode::Scalar, DNode::Scalar], vec![Node::Leaf(0), Node::Leaf(0)]);
        assert_eq!(run(&t, &leaves, &x).unwrap(), Node::Leaf(0));
    }
}
