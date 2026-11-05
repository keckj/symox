//! SymPy 1.14 `convert_array_to_matrix`
//! (`sympy/tensor/array/expressions/from_array_to_matrix.py`) and the
//! `_EditArrayContraction` editor (`array_expressions.py`) over the index
//! trees of [`crate::arrayexpr`].
//!
//! Every pass is a literal port: `_array2matrix`, `_remove_trivial_dims`,
//! `identify_hadamard_products`, `identify_removable_identity_matrices`,
//! `remove_identity_matrices`, `_array_contraction_to_diagonal_multiple_identity`,
//! `_array_diag2contr_diagmatrix`, `_support_function_tp1_recognize`,
//! `ArrayContraction.split_multiple_contractions` and
//! `ArrayContraction.flatten_contraction_of_diagonal`, including the order in
//! which SymPy mutates the editor and the exceptions SymPy raises. The matrix
//! operations SymPy performs on the recognized factors (`MatMul(...).doit()`,
//! `Transpose`, `Trace(...)._normalize()`, `hadamard_product`,
//! `HadamardPower`, `KroneckerProduct`, `DiagMatrix`,
//! `ElementwiseApplyFunction`) are the native matrix expressions.
//!
//! Leaves are a table of sources: native matrices and scalars (interned, so
//! that structural equality of trees is SymPy's `==`), `OneArray`s, opaque
//! array leaves (array symbols, kept by the caller) and the functions of
//! element-wise applications. An input this port does not model returns
//! [`AErr::Unsupported`]; the caller then runs the stock SymPy algorithm.
use crate::arrayexpr::{self as ax, AErr, Ctx, LeafInfo, Node, R};
use crate::expr::{Expr, Func, Kind};
use crate::matrixexpr::{MatrixExpr, MatrixKind, MatrixOperand};
use crate::sortkey::{matrix_sort_key, sort_key, Key};
use rustc_hash::FxHashMap;

#[derive(Clone, Debug)]
pub enum Src {
    Matrix(MatrixExpr),
    Scalar(Expr),
    /// `OneArray(*shape)` (dimension ids).
    One(Vec<usize>),
    /// An array leaf the caller keeps (an `ArraySymbol`).
    Opaque,
    /// The function (a Lambda) of an `ArrayElementwiseApplyFunc`.
    Fn(Expr),
    /// SymPy's `Mul(coeff, *factors)` of a commutative coefficient and
    /// non-commutative array leaves (`OneArray`s, array symbols), in order:
    /// what `OneArray(1, k)*2` or `A*B` evaluate to. It is a rank-zero
    /// operand that is no array (`_a2m_tensor_product` takes it for a
    /// scalar, `MatMul` refuses it).
    NcMul(Expr, Vec<usize>),
    /// SymPy's `MatMul(x, *matrices)` with an array leaf `x` as its
    /// coefficient (`M*OneArray(1, k)`, `A*M`): a matrix expression of the
    /// shape of its matrices that SymPy cannot evaluate or print
    /// (`as_coeff_matrices` raises NotImplementedError).
    NcMatMul(Vec<usize>, Vec<MatrixExpr>),
    /// `Trace` of an [`Src::NcMatMul`] leaf (a scalar that is no `Expr`).
    NcTrace(usize),
}

impl Src {
    /// A leaf SymPy builds but cannot print or evaluate further.
    pub fn is_held(&self) -> bool { matches!(self, Src::NcMatMul(..) | Src::NcTrace(_)) }
}

type Order<'a> = &'a mut dyn FnMut(&[&Node], &[Src], &[Expr]) -> R<Vec<usize>>;

pub struct A2m<'a> {
    pub srcs: Vec<Src>,
    pub info: Vec<LeafInfo>,
    pub dims: Vec<Expr>,
    pub ones: Vec<bool>,
    order: Order<'a>,
    alias: Vec<usize>,
    mats: FxHashMap<MatrixExpr, usize>,
    scals: FxHashMap<Expr, usize>,
    onel: FxHashMap<Vec<usize>, usize>,
    ncm: FxHashMap<(Expr, Vec<usize>), usize>,
    dimid: FxHashMap<Expr, usize>,
}

const NC_MATMUL: &str = "noncommutative scalars in MatMul are not supported.";

fn unsupported<T>(what: &str) -> R<T> { Err(AErr::Unsupported(what.into())) }

/// Python's `seq[i]` (negative indices count from the end).
fn pyidx(len: usize, i: i64, what: &str) -> R<usize> {
    let j = if i < 0 { i + len as i64 } else { i };
    if j < 0 || j >= len as i64 { return Err(AErr::Index(format!("{what} index out of range"))); }
    Ok(j as usize)
}

/// `_af_invert`.
fn af_invert(a: &[i64]) -> R<Vec<usize>> {
    let mut inv = vec![0usize; a.len()];
    for (i, &ai) in a.iter().enumerate() {
        let j = pyidx(a.len(), ai, "list assignment")?;
        inv[j] = i;
    }
    Ok(inv)
}

/// `Permutation(form)` validation of an array form.
fn check_perm(p: &[usize]) -> R<()> {
    let mut seen = vec![false; p.iter().copied().max().map_or(0, |m| m + 1)];
    for &i in p {
        if seen[i] { return Err(AErr::Value("there were repeated elements.".into())); }
        seen[i] = true;
    }
    if seen.len() != p.len() {
        return Err(AErr::Value(format!("Integers 0 through {} must be present.", seen.len() - 1)));
    }
    Ok(())
}

/// `_get_mapping_from_subranks`.
fn mapping(subranks: &[usize]) -> Vec<(usize, usize)> {
    let mut out = vec![];
    for (i, &r) in subranks.iter().enumerate() { for j in 0..r { out.push((i, j)); } }
    out
}

/// `ArrayContraction._get_contraction_tuples`: each contracted axis as
/// (argument, axis within the argument), through `_mapping`.
pub fn contraction_tuples(subranks: &[usize], groups: &[Vec<usize>]) -> R<Vec<Vec<(usize, usize)>>> {
    let map = mapping(subranks);
    groups.iter().map(|g| g.iter().map(|&j| map.get(j).copied()
        .ok_or_else(|| AErr::Index("list index out of range".into()))).collect()).collect()
}

/// `ArrayContraction._contraction_tuples_to_contraction_indices`: the
/// inverse of [`contraction_tuples`], on the cumulative ranks.
pub fn contraction_tuples_to_indices(subranks: &[usize], tuples: &[Vec<(usize, usize)>]) -> R<Vec<Vec<usize>>> {
    let mut cumulative = vec![0];
    cumulative.extend(accumulate(subranks.iter().copied()));
    tuples.iter().map(|g| g.iter().map(|&(j, k)| cumulative.get(j).map(|c| c + k)
        .ok_or_else(|| AErr::Index("list index out of range".into()))).collect()).collect()
}

fn accumulate(v: impl IntoIterator<Item = usize>) -> Vec<usize> {
    let mut acc = 0;
    v.into_iter().map(|x| { acc += x; acc }).collect()
}

/// `_combine_removed`.
pub fn combine_removed(removed1: &[usize], removed2: &[usize]) -> Vec<usize> {
    let mut r1 = removed1.to_vec();
    r1.sort_unstable();
    let mut r2 = removed2.to_vec();
    r2.sort_unstable();
    let (mut i, mut j) = (0, 0);
    let mut out = vec![];
    loop {
        if j >= r2.len() {
            out.extend_from_slice(&r1[i..]);
            break;
        } else if i < r1.len() && r1[i] <= i + r2[j] {
            out.push(r1[i]);
            i += 1;
        } else {
            out.push(i + r2[j]);
            j += 1;
        }
    }
    out
}

fn codegen(n: &Node) -> bool {
    matches!(n, Node::Tensor(_) | Node::Add(_) | Node::Permute(..) | Node::Contract(..) | Node::Diagonal(..)
        | Node::Reshape(..) | Node::PowApply(..) | Node::Apply(..))
}

/// `Trace(m)._normalize()`'s argument: the cyclic rotation of a matrix
/// product starting at its smallest factor (by `default_sort_key`), the
/// product transposed first when that factor is a transpose.
pub fn trace_normalized_arg(m: &MatrixExpr) -> R<MatrixExpr> {
    let MatrixKind::Mul(c, fs) = m.kind() else { return Ok(m.clone()) };
    fn argmin(keys: &[Key]) -> usize {
        let mut best = 0;
        for i in 1..keys.len() { if keys[i] < keys[best] { best = i; } }
        best
    }
    let off = usize::from(!c.is_one());
    let mut keys: Vec<Key> = vec![];
    if off == 1 { keys.push(sort_key(c)); }
    keys.extend(fs.iter().map(|f| match f.kind() { MatrixKind::Transpose(x) => matrix_sort_key(x), _ => matrix_sort_key(f) }));
    if keys.is_empty() { return Ok(m.clone()); }
    let mut indmin = argmin(&keys);
    if indmin >= off && matches!(fs[indmin - off].kind(), MatrixKind::Transpose(_)) {
        let t = m.transpose(false)?.doit(true)?;
        let MatrixKind::Mul(c2, fs2) = t.kind() else { return unsupported("trace of a transposed product that is not a product") };
        let (c, fs) = (c2.clone(), fs2.clone());
        let off2 = usize::from(!c.is_one());
        let mut keys: Vec<Key> = vec![];
        if off2 == 1 { keys.push(sort_key(&c)); }
        keys.extend(fs.iter().map(matrix_sort_key));
        if keys.is_empty() { return unsupported("empty product"); }
        indmin = argmin(&keys);
        return rotate(&c, &fs, indmin, off2);
    }
    rotate(c, fs, indmin, off)
}

/// `MatMul.fromiter(args[indmin:] + args[:indmin])` over SymPy's arguments
/// `[c] + fs` (`c` present when `off == 1`). The native product keeps its
/// scalar coefficient in front: a rotation that would move the coefficient
/// behind matrix factors gives the same product with the coefficient in
/// front (SymPy's unevaluated MatMul prints it in place).
fn rotate(c: &Expr, fs: &[MatrixExpr], indmin: usize, off: usize) -> R<MatrixExpr> {
    let j = if indmin >= off { indmin - off } else { 0 };
    let rotated: Vec<MatrixExpr> = fs[j..].iter().chain(fs[..j].iter()).cloned().collect();
    Ok(MatrixExpr::mul(c.clone(), &rotated, false)?)
}

/// `expr.function.expr.match(w*d**p)` with `w`, `p` excluding `d`: `(w, p)`.
pub(crate) fn match_power(body: &Expr, d: &Expr) -> Option<(Expr, Expr)> {
    let has = |e: &Expr| crate::traverse::has(e, d);
    if !has(body) { return Some((body.clone(), Expr::zero())); }
    if body == d { return Some((Expr::one(), Expr::one())); }
    match body.kind() {
        Kind::Pow(b, e) if b == d && !has(e) => Some((Expr::one(), e.clone())),
        Kind::Mul(..) => {
            let args = body.args();
            let dep: Vec<usize> = (0..args.len()).filter(|&i| has(&args[i])).collect();
            if dep.len() != 1 { return None; }
            let (_, p) = match_power(&args[dep[0]], d).filter(|(w, _)| w.is_one())?;
            let rest: Vec<Expr> = args.iter().enumerate().filter(|(i, _)| *i != dep[0]).map(|(_, a)| a.clone()).collect();
            Some((crate::ops::mul(&rest), p))
        },
        _ => None,
    }
}

// ------------------------------------------------------------------ editor

/// `_ArgE`: an argument of the edited tensor product and the contraction
/// number of each of its axes (negative: diagonal, None: free).
#[derive(Clone, Debug)]
pub struct ArgE {
    pub id: usize,
    pub element: Option<Node>,
    pub indices: Vec<Option<i64>>,
}

/// `_EditArrayContraction`.
#[derive(Clone, Debug, Default)]
pub struct Editor {
    pub args: Vec<ArgE>,
    pub n_contr: i64,
    pub track: Option<Vec<Vec<i64>>>,
    pub next_id: usize,
}

impl Editor {
    pub fn arg(&mut self, element: Option<Node>, indices: Vec<Option<i64>>) -> ArgE {
        let id = self.next_id;
        self.next_id += 1;
        ArgE { id, element, indices }
    }

    pub fn pos(&self, id: usize) -> R<usize> {
        self.args.iter().position(|a| a.id == id).ok_or_else(|| AErr::Value("_ArgE is not in list".into()))
    }

    fn get(&mut self, id: usize) -> R<&mut ArgE> {
        let p = self.pos(id)?;
        Ok(&mut self.args[p])
    }

    pub fn insert_after(&mut self, id: usize, new: ArgE) -> R<()> {
        let p = self.pos(id)?;
        self.args.insert(p + 1, new);
        Ok(())
    }

    pub fn remove(&mut self, id: usize) -> R<()> {
        let p = self.args.iter().position(|a| a.id == id).ok_or_else(|| AErr::Value("list.remove(x): x not in list".into()))?;
        self.args.remove(p);
        Ok(())
    }

    pub fn new_index(&mut self) -> i64 {
        self.n_contr += 1;
        self.n_contr - 1
    }

    pub fn refresh_indices(&mut self) {
        let mut keys: Vec<i64> = vec![];
        for a in &self.args { for i in a.indices.iter().flatten() { if !keys.contains(i) { keys.push(*i); } } }
        keys.sort_unstable();
        self.n_contr = keys.len() as i64;
        for a in &mut self.args {
            a.indices = a.indices.iter().map(|i| i.and_then(|i| keys.iter().position(|&k| k == i).map(|p| p as i64))).collect();
        }
    }

    pub fn contraction_indices(&self) -> R<Vec<Vec<usize>>> {
        let n = self.n_contr.max(0) as usize;
        let mut out: Vec<Vec<usize>> = vec![vec![]; n];
        let mut pos = 0;
        for a in &self.args {
            for j in a.indices.iter() {
                if let Some(j) = j { out[pyidx(n, *j, "list")?].push(pos); }
                pos += 1;
            }
        }
        Ok(out)
    }

    pub fn mapping_for_index(&self, ind: i64) -> R<Vec<(usize, usize)>> {
        if ind >= self.n_contr { return Err(AErr::Value("index value exceeding the index range".into())); }
        let mut out = vec![];
        for (i, a) in self.args.iter().enumerate() {
            for (j, x) in a.indices.iter().enumerate() { if *x == Some(ind) { out.push((i, j)); } }
        }
        Ok(out)
    }

    pub fn ind_rel_pos(&self) -> R<Vec<Vec<(usize, usize)>>> {
        let n = self.n_contr.max(0) as usize;
        let mut out: Vec<Vec<(usize, usize)>> = vec![vec![]; n];
        for (i, a) in self.args.iter().enumerate() {
            for (j, x) in a.indices.iter().enumerate() {
                if let Some(x) = x { out[pyidx(n, *x, "list")?].push((i, j)); }
            }
        }
        Ok(out)
    }

    pub fn count_with(&self, ind: Option<i64>) -> usize { self.args.iter().filter(|a| a.indices.contains(&ind)).count() }

    pub fn args_with(&self, ind: Option<i64>) -> Vec<usize> {
        self.args.iter().filter(|a| a.indices.contains(&ind)).map(|a| a.id).collect()
    }

    pub fn n_diag(&self) -> usize {
        let mut seen: Vec<i64> = vec![];
        for a in &self.args { for i in a.indices.iter().flatten() { if *i < 0 && !seen.contains(i) { seen.push(*i); } } }
        seen.len()
    }

    pub fn track_start(&mut self) {
        let mut perm: Vec<Vec<i64>> = vec![];
        let mut diag: Vec<i64> = vec![];
        let (mut counter, mut counter2) = (0i64, -1i64);
        for a in &self.args {
            let mut p = vec![];
            for i in &a.indices {
                if let Some(i) = i {
                    if *i < 0 { diag.push(counter2); counter2 -= 1; }
                    continue;
                }
                p.push(counter);
                counter += 1;
            }
            perm.push(p);
        }
        let max_ind = perm.iter().map(|p| p.iter().copied().max().unwrap_or(-1)).max().unwrap_or(-1);
        perm.push(diag.iter().map(|i| max_ind - i).collect());
        self.track = Some(perm);
    }

    pub fn track_merge(&mut self, dest: usize, from: usize) -> R<()> {
        let (d, f) = (self.pos(dest)?, self.pos(from)?);
        let track = self.track.as_mut().ok_or_else(|| AErr::Type("'NoneType' object is not subscriptable".into()))?;
        let n = track.len();
        let (d, f) = (pyidx(n, d as i64, "list")?, pyidx(n, f as i64, "pop")?);
        let moved = track[f].clone();
        track[d].extend(moved);
        track.remove(f);
        Ok(())
    }

    pub fn abs_free_range(&self, id: usize) -> R<(usize, usize)> {
        let mut c = 0;
        for a in &self.args {
            let n = a.indices.iter().filter(|i| i.is_none()).count();
            if a.id == id { return Ok((c, c + n)); }
            c += n;
        }
        Err(AErr::Index("argument not found".into()))
    }

    pub fn abs_range(&self, id: usize) -> R<(usize, usize)> {
        let mut c = 0;
        for a in &self.args {
            let n = a.indices.len();
            if a.id == id { return Ok((c, c + n)); }
            c += n;
        }
        Err(AErr::Index("argument not found".into()))
    }
}

// ------------------------------------------------------------------ context

impl<'a> A2m<'a> {
    pub fn new(srcs: Vec<Src>, info: Vec<LeafInfo>, dims: Vec<Expr>, order: Order<'a>) -> Self {
        let ones = dims.iter().map(Expr::is_one).collect();
        let alias = (0..srcs.len()).collect();
        let mut s = A2m { srcs, info, dims, ones, order, alias, mats: FxHashMap::default(), scals: FxHashMap::default(),
            onel: FxHashMap::default(), ncm: FxHashMap::default(), dimid: FxHashMap::default() };
        for (i, d) in s.dims.iter().enumerate() { s.dimid.entry(d.clone()).or_insert(i); }
        for i in 0..s.srcs.len() {
            s.alias[i] = match &s.srcs[i] {
                Src::Matrix(m) => *s.mats.entry(m.clone()).or_insert(i),
                Src::Scalar(e) => *s.scals.entry(e.clone()).or_insert(i),
                Src::One(sh) => *s.onel.entry(sh.clone()).or_insert(i),
                _ => i,
            };
        }
        // The factors of a product refer to the unified leaves.
        for i in 0..s.srcs.len() {
            if let Src::NcMul(c, f) = &s.srcs[i] {
                let (c, f) = (c.clone(), f.iter().map(|&j| s.alias.get(j).copied().unwrap_or(j)).collect::<Vec<_>>());
                s.srcs[i] = Src::NcMul(c.clone(), f.clone());
                s.alias[i] = *s.ncm.entry((c.clone(), f)).or_insert(i);
            };
        }
        s
    }

    /// The tree with each leaf replaced by the first equal leaf of the
    /// table, so that structural comparisons of trees are SymPy's `==`.
    pub fn unify(&self, n: Node) -> Node {
        let rec = |v: Vec<Node>| v.into_iter().map(|x| self.unify(x)).collect();
        match n {
            Node::Leaf(i) => Node::Leaf(self.alias.get(i).copied().unwrap_or(i)),
            Node::Tensor(a) => Node::Tensor(rec(a)),
            Node::Add(a) => Node::Add(rec(a)),
            Node::Permute(e, p) => Node::Permute(Box::new(self.unify(*e)), p),
            Node::Contract(e, g) => Node::Contract(Box::new(self.unify(*e)), g),
            Node::Diagonal(e, g) => Node::Diagonal(Box::new(self.unify(*e)), g),
            Node::Reshape(e, s) => Node::Reshape(Box::new(self.unify(*e)), s),
            Node::PowApply(e, x) => Node::PowApply(Box::new(self.unify(*e)), self.alias.get(x).copied().unwrap_or(x)),
            Node::Apply(e, f) => Node::Apply(Box::new(self.unify(*e)), f),
            z @ Node::Zero(_) => z,
        }
    }

    fn with<T>(&mut self, f: impl FnOnce(&mut Ctx<'_>) -> R<T>) -> R<T> {
        let A2m { info, ones, order, srcs, dims, .. } = self;
        let (srcs, dims): (&[Src], &[Expr]) = (srcs, dims);
        let mut o = |n: &[&Node]| (*order)(n, srcs, dims);
        let mut c = Ctx { leaves: info, ones, order: &mut o };
        f(&mut c)
    }

    fn ro<T>(&self, f: impl FnOnce(&Ctx<'_>) -> T) -> T {
        let mut o = |_: &[&Node]| -> R<Vec<usize>> { Ok(vec![]) };
        let c = Ctx { leaves: &self.info, ones: &self.ones, order: &mut o };
        f(&c)
    }

    pub fn rank(&self, n: &Node) -> usize { self.ro(|c| c.rank(n)) }
    pub fn shape(&self, n: &Node) -> Option<Vec<usize>> { self.ro(|c| c.shape(n)) }
    fn shape_r(&self, n: &Node) -> R<Vec<usize>> {
        self.shape(n).ok_or_else(|| AErr::Unsupported("array of undefined shape".into()))
    }
    pub fn subranks(&self, n: &Node) -> Vec<usize> { self.ro(|c| c.subranks(n)) }
    /// `_get_subrank`.
    fn subrank_of(&self, n: &Node) -> usize { if codegen(n) { self.subranks(n).iter().sum() } else { self.rank(n) } }

    fn has_one(&self, s: &[usize]) -> bool { s.iter().any(|&d| self.ones[d]) }
    fn is_11(&self, s: &[usize]) -> bool { s.len() == 2 && self.ones[s[0]] && self.ones[s[1]] }

    fn dim(&mut self, e: &Expr) -> usize {
        if let Some(&i) = self.dimid.get(e) { return i; }
        self.dims.push(e.clone());
        self.ones.push(e.is_one());
        self.dimid.insert(e.clone(), self.dims.len() - 1);
        self.dims.len() - 1
    }

    fn push(&mut self, src: Src, info: LeafInfo) -> Node {
        self.srcs.push(src);
        self.info.push(info);
        Node::Leaf(self.srcs.len() - 1)
    }

    pub fn mat(&mut self, m: MatrixExpr) -> Node {
        if let Some(&i) = self.mats.get(&m) { return Node::Leaf(i); }
        let shape = vec![self.dim(m.rows()), self.dim(m.cols())];
        let zero = matches!(m.kind(), MatrixKind::Zero);
        self.mats.insert(m.clone(), self.srcs.len());
        self.push(Src::Matrix(m), LeafInfo { rank: 2, shape: Some(shape), zero })
    }

    pub fn scal(&mut self, e: Expr) -> Node {
        if let Some(&i) = self.scals.get(&e) { return Node::Leaf(i); }
        self.scals.insert(e.clone(), self.srcs.len());
        self.push(Src::Scalar(e), LeafInfo { rank: 0, shape: Some(vec![]), zero: false })
    }

    fn one(&mut self, shape: Vec<usize>) -> Node {
        if let Some(&i) = self.onel.get(&shape) { return Node::Leaf(i); }
        self.onel.insert(shape.clone(), self.srcs.len());
        self.push(Src::One(shape.clone()), LeafInfo { rank: shape.len(), shape: Some(shape), zero: false })
    }

    pub fn as_mat(&self, n: &Node) -> Option<&MatrixExpr> {
        match n { Node::Leaf(i) => match &self.srcs[*i] { Src::Matrix(m) => Some(m), _ => None }, _ => None }
    }
    fn as_scal(&self, n: &Node) -> Option<&Expr> {
        match n { Node::Leaf(i) => match &self.srcs[*i] { Src::Scalar(e) => Some(e), _ => None }, _ => None }
    }
    fn is_mat(&self, n: &Node) -> bool { self.as_mat(n).is_some() }
    /// `isinstance(n, MatrixExpr)`.
    fn is_matexpr(&self, n: &Node) -> bool {
        self.is_mat(n) || matches!(n, Node::Leaf(i) if matches!(self.srcs[*i], Src::NcMatMul(..)))
    }
    /// Whether a tree holds a leaf SymPy cannot print ([`Src::is_held`]).
    pub fn holds_unprintable(&self, n: &Node) -> bool {
        match n {
            Node::Leaf(i) => self.srcs.get(*i).is_some_and(Src::is_held),
            Node::Zero(_) => false,
            Node::Tensor(a) | Node::Add(a) => a.iter().any(|x| self.holds_unprintable(x)),
            Node::Permute(e, _) | Node::Contract(e, _) | Node::Diagonal(e, _) | Node::Reshape(e, _) | Node::PowApply(e, _)
            | Node::Apply(e, _) => self.holds_unprintable(e),
        }
    }
    fn is_mat_opt(&self, n: &Option<Node>) -> bool { n.as_ref().is_some_and(|n| self.is_mat(n)) }
    fn is_one_array(&self, n: &Node) -> bool { matches!(n, Node::Leaf(i) if matches!(self.srcs[*i], Src::One(_))) }
    /// `isinstance(x, Identity)`.
    fn is_identity(&self, n: &Option<Node>) -> bool {
        n.as_ref().and_then(|n| self.as_mat(n)).is_some_and(|m| matches!(m.kind(), MatrixKind::Identity))
    }
    /// `getattr(x, "is_Identity", False)` for a 1x1 matrix.
    fn is_identity_attr(&self, m: &MatrixExpr) -> bool {
        match m.kind() {
            MatrixKind::Identity => true,
            MatrixKind::One => m.rows().is_one() && m.cols().is_one(),
            MatrixKind::Explicit(x) => x.shape() == (1, 1) && x.get(0, 0).is_one(),
            _ => false,
        }
    }
    /// `isinstance(x, (MatrixExpr, _ArrayExpr, _CodegenArrayAbstract))`.
    fn arraylike(&self, n: &Node) -> bool {
        match n {
            Node::Leaf(i) => matches!(self.srcs[*i], Src::Matrix(_) | Src::One(_) | Src::Opaque | Src::NcMatMul(..)),
            _ => true,
        }
    }
    fn class(&self, n: &Node) -> R<&'static str> {
        Ok(match n {
            Node::Tensor(_) => "ArrayTensorProduct",
            Node::Add(_) => "ArrayAdd",
            Node::Permute(..) => "PermuteDims",
            Node::Contract(..) => "ArrayContraction",
            Node::Diagonal(..) => "ArrayDiagonal",
            Node::Zero(_) => "ZeroArray",
            Node::Reshape(..) => "Reshape",
            Node::Apply(..) | Node::PowApply(..) => "ArrayElementwiseApplyFunc",
            Node::Leaf(i) => match self.srcs[*i] {
                Src::One(_) => "OneArray",
                Src::Opaque => "ArraySymbol",
                Src::NcMul(..) => "Mul",
                Src::NcMatMul(..) => "MatMul",
                Src::NcTrace(_) => "Trace",
                _ => return unsupported("matrix or scalar class name"),
            },
        })
    }

    // -------------------------------------------------------------- matrix ops

    /// `_mul(a, b)`: `MatMul(a, b).doit()` or the scalar product.
    fn mm(&mut self, a: &Node, b: &Node) -> R<Node> {
        let (ma, mb, sa, sb) = (self.as_mat(a).cloned(), self.as_mat(b).cloned(), self.as_scal(a).cloned(), self.as_scal(b).cloned());
        let out = match (ma, mb, sa, sb) {
            (Some(x), Some(y), ..) => MatrixExpr::mul(Expr::one(), &[x, y], false)?.doit(true)?,
            (Some(x), None, _, Some(s)) | (None, Some(x), Some(s), _) => MatrixExpr::mul(s, &[x], false)?.doit(true)?,
            (None, None, Some(x), Some(y)) => return Ok(self.scal(crate::ops::mul2(&x, &y))),
            (Some(_), None, ..) | (None, Some(_), ..) => {
                // `MatMul(...).doit()` with a non-commutative coefficient.
                let (m, other) = if self.is_mat(a) { (self.as_mat(a).cloned().unwrap(), b) } else { (self.as_mat(b).cloned().unwrap(), a) };
                let Node::Leaf(i) = other else { return unsupported("product of a matrix and an array expression") };
                match self.srcs[*i] {
                    Src::NcMul(..) | Src::NcTrace(_) => return Err(AErr::NotImplemented(NC_MATMUL.into())),
                    Src::One(_) | Src::Opaque => {
                        // `MatMul(x, m).doit()` keeps the array leaf as a
                        // coefficient (its commutativity is unknown).
                        let shape = vec![self.dim(&m.rows().clone()), self.dim(&m.cols().clone())];
                        let i = *i;
                        return Ok(self.push(Src::NcMatMul(vec![i], vec![m]), LeafInfo { rank: 2, shape: Some(shape), zero: false }));
                    },
                    _ => return unsupported("product of a matrix and an array expression"),
                }
            },
            _ => return self.mul_scalars(&[a.clone(), b.clone()]),
        };
        Ok(self.mat(out))
    }

    /// The commutative coefficient and the non-commutative factors of an
    /// operand of SymPy's `Mul`.
    fn mul_parts(&self, n: &Node) -> Option<(Expr, Vec<usize>)> {
        let Node::Leaf(i) = n else { return None };
        match &self.srcs[*i] {
            Src::Scalar(e) => Some((e.clone(), vec![])),
            Src::NcMul(c, f) => Some((c.clone(), f.clone())),
            Src::One(_) | Src::Opaque | Src::NcTrace(_) => Some((Expr::one(), vec![*i])),
            _ => None,
        }
    }

    /// `Mul(*factors)` (`Mul.fromiter`): the commutative factors multiply,
    /// the array factors keep their order.
    fn mul_scalars(&mut self, fs: &[Node]) -> R<Node> {
        let (mut cs, mut nc) = (vec![], vec![]);
        for f in fs {
            let Some((c, n)) = self.mul_parts(f) else { return unsupported("product of an array expression") };
            cs.push(c);
            nc.extend(n);
        }
        if nc.is_empty() { return Ok(self.scal(crate::ops::mul(&cs))); }
        // `Mul` distributes a rational over a sum only for two operands and
        // collects equal adjacent factors into a power; neither is modeled.
        if cs.iter().any(Expr::is_add) { return unsupported("product of an array expression and a sum"); }
        if nc.windows(2).any(|w| w[0] == w[1]) { return unsupported("power of an array expression"); }
        let coeff = crate::ops::mul(&cs);
        if coeff.is_add() { return unsupported("product of an array expression and a sum"); }
        if coeff.is_zero() { return unsupported("product of an array expression and zero"); }
        if coeff.is_one() && nc.len() == 1 { return Ok(Node::Leaf(nc[0])); }
        let key = (coeff, nc);
        if let Some(&i) = self.ncm.get(&key) { return Ok(Node::Leaf(i)); }
        self.ncm.insert(key.clone(), self.srcs.len());
        Ok(self.push(Src::NcMul(key.0, key.1), LeafInfo { rank: 0, shape: Some(vec![]), zero: false }))
    }

    fn mul_matrices(&mut self, coef: Expr, fs: &[MatrixExpr]) -> R<Node> {
        let m = MatrixExpr::mul(coef, fs, false)?.doit(true)?;
        Ok(self.mat(m))
    }

    /// `_a2m_transpose`.
    fn tr(&mut self, n: Node) -> R<Node> {
        if codegen(&n) { return self.permute_checked(n, vec![1, 0]); }
        let Some(m) = self.as_mat(&n).cloned() else { return unsupported("transpose of a non-matrix") };
        let t = m.transpose(false)?.doit(true)?;
        Ok(self.mat(t))
    }

    fn mat_t(&mut self, n: &Node, evaluate: bool) -> R<MatrixExpr> {
        let Some(m) = self.as_mat(n) else { return unsupported("transpose of a non-matrix") };
        Ok(m.transpose(evaluate)?)
    }

    /// The attribute `x.T`, which only matrices have.
    fn attr_t(&mut self, n: &Node) -> R<MatrixExpr> {
        if let Some(m) = self.as_mat(n) { return Ok(m.transpose(true)?); }
        let c = match n {
            Node::Leaf(i) => match &self.srcs[*i] {
                Src::Scalar(e) => e.class_name().to_string(),
                Src::One(_) | Src::Opaque | Src::NcMul(..) | Src::NcTrace(_) => self.class(n)?.to_string(),
                _ => return unsupported("transpose of a non-matrix"),
            },
            _ => return unsupported("transpose of a non-matrix"),
        };
        Err(AErr::Attribute(format!("'{c}' object has no attribute 'T'")))
    }

    /// `Trace` of a square [`Src::NcMatMul`] leaf: a held scalar.
    fn nc_trace(&mut self, n: &Node) -> Option<Node> {
        let Node::Leaf(i) = *n else { return None };
        let Some(LeafInfo { shape: Some(sh), .. }) = self.info.get(i) else { return None };
        if !matches!(self.srcs[i], Src::NcMatMul(..)) || sh[0] != sh[1] { return None; }
        Some(self.push(Src::NcTrace(i), LeafInfo { rank: 0, shape: Some(vec![]), zero: false }))
    }

    /// `_a2m_trace`.
    fn trace(&mut self, n: Node) -> R<Node> {
        if codegen(&n) { return self.with(|c| c.contract(n, vec![vec![0, 1]])); }
        if let Some(t) = self.nc_trace(&n) { return Ok(t); }
        let Some(m) = self.as_mat(&n).cloned() else { return unsupported("trace of a non-matrix") };
        let t = m.trace_unevaluated()?;
        Ok(self.scal(t))
    }

    /// `Trace(m)._normalize()`.
    fn trace_norm(&mut self, n: &Node) -> R<Node> {
        // `_normalize` of a product holding an array coefficient either
        // rotates it (the trace stays held) or transposes it, where
        // `as_coeff_matrices` raises the same NotImplementedError the held
        // trace raises.
        if let Some(t) = self.nc_trace(n) { return Ok(t); }
        let Some(m) = self.as_mat(n).cloned() else { return unsupported("trace of a non-matrix") };
        let t = trace_normalized_arg(&m)?.trace_unevaluated()?;
        Ok(self.scal(t))
    }

    /// `_a2m_add`.
    fn madd(&mut self, args: Vec<Node>) -> R<Node> {
        if !args.iter().any(codegen) {
            let ms: Option<Vec<MatrixExpr>> = args.iter().map(|a| self.as_mat(a).cloned()).collect();
            // `MatAdd` accepts matrix expressions only.
            let Some(ms) = ms else {
                if args.iter().all(|a| self.is_matexpr(a)) { return unsupported("MatAdd of an unevaluated product"); }
                return Err(AErr::Type("Mix of Matrix and Scalar symbols".into()));
            };
            let m = MatrixExpr::add(&ms, false)?.doit(true)?;
            return Ok(self.mat(m));
        }
        self.array_add(args)
    }

    /// `_array_add` with `ArrayAdd`'s rank and shape checks.
    fn array_add(&mut self, args: Vec<Node>) -> R<Node> {
        if args.iter().any(|a| self.as_scal(a).is_some()) { return unsupported("ArrayAdd of a scalar"); }
        let mut ranks: Vec<usize> = args.iter().map(|a| self.rank(a)).collect();
        ranks.sort_unstable();
        ranks.dedup();
        if ranks.len() != 1 { return Err(AErr::Value("summing arrays of different ranks".into())); }
        let mut shapes: Vec<Vec<usize>> = args.iter().filter_map(|a| self.shape(a)).collect();
        shapes.sort_unstable();
        shapes.dedup();
        if shapes.len() > 1 { return Err(AErr::Value("mismatching shapes in addition".into())); }
        self.with(|c| c.add(args))
    }

    /// `_a2m_tensor_product`.
    fn tp(&mut self, args: Vec<Node>) -> R<Node> {
        let mut scalars = vec![];
        let mut arrays = vec![];
        for a in args {
            if self.arraylike(&a) { arrays.push(a); } else {
                if self.mul_parts(&a).is_none() { return unsupported("tensor factor"); }
                scalars.push(a);
            }
        }
        let s = self.mul_scalars(&scalars)?;
        if arrays.is_empty() { return Ok(s); }
        if !self.as_scal(&s).is_some_and(Expr::is_one) {
            if codegen(&arrays[0]) { arrays.insert(0, s); } else { arrays[0] = self.mm(&arrays[0].clone(), &s)?; }
        }
        self.with(|c| c.tensor(arrays))
    }

    /// `hadamard_product(*matrices)`.
    fn hadamard(&mut self, fs: Vec<MatrixExpr>) -> R<MatrixExpr> {
        if fs.is_empty() { return Err(AErr::Type("Empty Hadamard product is undefined".into())); }
        if fs.len() == 1 { return Ok(fs.into_iter().next().unwrap()); }
        Ok(MatrixExpr::hadamard_product(&fs, false)?.doit(true)?)
    }

    /// `_permute_dims` of a computed array form.
    fn permute_checked(&mut self, e: Node, p: Vec<usize>) -> R<Node> {
        check_perm(&p)?;
        self.with(|c| c.permute(e, p))
    }
    fn contract(&mut self, e: Node, g: Vec<Vec<usize>>) -> R<Node> { self.with(|c| c.contract(e, g)) }
    fn diagonal(&mut self, e: Node, g: Vec<Vec<usize>>, trivial: bool) -> R<Node> { self.with(|c| c.diagonal(e, g, trivial)) }
    fn tensor(&mut self, a: Vec<Node>) -> R<Node> { self.with(|c| c.tensor(a)) }

    fn diag_matrix(&mut self, n: &Node, doit: bool) -> R<Node> {
        let Some(m) = self.as_mat(n).cloned() else { return unsupported("DiagMatrix of a non-matrix") };
        let d = MatrixExpr::diag(&m)?;
        let d = if doit { d.doit(true)? } else { d };
        Ok(self.mat(d))
    }

    fn ask_diagonal(&self, n: &Option<Node>) -> bool {
        let Some(m) = n.as_ref().and_then(|n| self.as_mat(n)) else { return false };
        let proxy = Expr::from_kind(Kind::Func(crate::expr::Func::MatrixProxy(std::sync::Arc::new(m.clone())), vec![]));
        let prop = crate::ask::q(crate::ask::Pred::Diagonal, &proxy);
        matches!(crate::ask::ask(&prop, &Expr::boolean(true)), Ok(Some(true)))
    }

    // -------------------------------------------------------------- editor

    /// `_EditArrayContraction(base_array)`.
    pub fn editor(&self, base: &Node) -> R<Editor> {
        let (map, expr, ci, diag): (Vec<(usize, usize)>, &Node, Vec<Vec<usize>>, Vec<Vec<usize>>) = match base {
            Node::Contract(e, cg) => (mapping(&self.subranks(base)), &**e, cg.clone(), vec![]),
            Node::Diagonal(e, dg) => match &**e {
                Node::Contract(e2, cg) => (mapping(&self.subranks(e)), &**e2, cg.clone(), ax::contraction_push_down(cg, dg)),
                _ => (vec![], &**e, vec![], dg.clone()),
            },
            Node::Tensor(_) => (vec![], base, vec![], vec![]),
            _ => return Err(AErr::NotImplemented(String::new())),
        };
        let elems: Vec<Node> = match expr { Node::Tensor(a) => a.clone(), o => vec![o.clone()] };
        let mut ed = Editor::default();
        for e in elems {
            let r = self.rank(&e);
            let a = ed.arg(Some(e), vec![None; r]);
            ed.args.push(a);
        }
        let set = |ed: &mut Editor, map: &[(usize, usize)], j: usize, v: i64| -> R<()> {
            let &(p, r) = map.get(j).ok_or_else(|| AErr::Key(j.to_string()))?;
            let slot = ed.args.get_mut(p).and_then(|a| a.indices.get_mut(r))
                .ok_or_else(|| AErr::Index("list assignment index out of range".into()))?;
            *slot = Some(v);
            Ok(())
        };
        for (i, t) in ci.iter().enumerate() { for &j in t { set(&mut ed, &map, j, i as i64)?; } }
        ed.n_contr = ci.len() as i64;
        let map2 = mapping(&self.subranks(base));
        for (i, t) in diag.iter().enumerate() { for &j in t { set(&mut ed, &map2, j, -1 - i as i64)?; } }
        Ok(ed)
    }

    /// `_EditArrayContraction.merge_scalars`.
    pub fn merge_scalars(&mut self, ed: &mut Editor) -> R<()> {
        let mut scalars = vec![];
        let mut kept = vec![];
        for a in std::mem::take(&mut ed.args) { if a.indices.is_empty() { scalars.push(a); } else { kept.push(a); } }
        ed.args = kept;
        let mut factors = vec![];
        for s in &scalars {
            match &s.element { Some(e) if self.mul_parts(e).is_some() => factors.push(e.clone()), _ => return unsupported("rank zero array factor") }
        }
        let scalar = self.mul_scalars(&factors)?;
        if ed.args.is_empty() {
            let a = ed.arg(Some(scalar), vec![]);
            ed.args.push(a);
        } else {
            let Some(el) = ed.args[0].element.clone() else { return unsupported("removed editor element") };
            ed.args[0].element = Some(self.tp(vec![scalar, el])?);
        }
        Ok(())
    }

    /// `_EditArrayContraction.to_array_contraction`.
    pub fn to_array(&mut self, ed: &mut Editor) -> R<Node> {
        let mut counter = 0usize;
        let mut diag: Vec<(i64, Vec<usize>)> = vec![];
        let mut freq: FxHashMap<Option<i64>, usize> = FxHashMap::default();
        for a in &ed.args { for i in &a.indices { *freq.entry(*i).or_default() += 1; } }
        let free = *freq.get(&None).unwrap_or(&0) as i64;
        let (mut inv1, mut inv2): (Vec<i64>, Vec<i64>) = (vec![], vec![]);
        let mut done: Vec<i64> = vec![];
        let mut counter4 = 0i64;
        for a in &mut ed.args {
            let mut counter2 = 0usize;
            for i in &a.indices {
                let Some(i) = *i else {
                    inv1.push(counter4);
                    counter2 += 1;
                    counter4 += 1;
                    continue;
                };
                if i >= 0 { continue; }
                let key = -1 - i;
                match diag.iter_mut().find(|(k, _)| *k == key) {
                    Some((_, v)) => v.push(counter + counter2),
                    None => diag.push((key, vec![counter + counter2])),
                }
                if freq.get(&Some(i)) == Some(&1) && !done.contains(&i) {
                    inv1.push(free - 1 - i);
                    done.push(i);
                } else if !done.contains(&i) {
                    inv2.push(free - 1 - i);
                    done.push(i);
                }
                counter2 += 1;
            }
            a.indices = a.indices.iter().map(|i| i.filter(|&i| i >= 0)).collect();
            counter += a.indices.iter().filter(|i| i.is_none()).count();
        }
        inv1.extend(inv2);
        let permutation = af_invert(&inv1)?;
        let filtered: Vec<Vec<usize>> = diag.into_iter().map(|(_, v)| v).filter(|v| v.len() > 1).collect();
        self.merge_scalars(ed)?;
        ed.refresh_indices();
        let mut elems = vec![];
        for a in &ed.args {
            match &a.element { Some(e) => elems.push(e.clone()), None => return unsupported("removed editor element") }
        }
        let ci = ed.contraction_indices()?;
        let t = self.tensor(elems)?;
        let e = self.contract(t, ci)?;
        let mut e2 = self.diagonal(e, filtered, false)?;
        if let Some(track) = &ed.track {
            let flat: Vec<i64> = track.iter().flatten().copied().collect();
            let p2 = af_invert(&flat)?;
            e2 = self.permute_checked(e2, p2)?;
        }
        self.permute_checked(e2, permutation)
    }

    // -------------------------------------------------------------- passes

    /// `_get_candidate_for_matmul_from_contraction`: (candidate id, transpose, index).
    fn candidate(&self, scan: &[Option<i64>], remaining: &[ArgE]) -> R<(Option<usize>, bool, i64)> {
        let scan: Vec<i64> = scan.iter().flatten().copied().collect();
        if scan.is_empty() { return Ok((None, false, -1)); }
        let mut transpose = false;
        let mut cand: Option<usize> = None;
        let mut cidx = -1i64;
        for a2 in remaining {
            if !self.is_mat_opt(&a2.element) { continue; }
            for &index in &scan {
                if cidx != -1 && cidx != index { continue; }
                if a2.indices.contains(&Some(index)) {
                    if a2.indices.iter().all(|i| *i == Some(index)) {
                        cand = None;
                        break;
                    }
                    if cand.is_none() {
                        cand = Some(a2.id);
                        cidx = index;
                        transpose = Some(index) == *a2.indices.get(1).ok_or_else(|| AErr::Index("list index out of range".into()))?;
                    } else {
                        cand = None;
                        break;
                    }
                }
            }
        }
        Ok((cand, transpose, cidx))
    }

    /// `_support_function_tp1_recognize(contraction_indices, args)`.
    pub fn tp1(&mut self, ci: Vec<Vec<usize>>, args: Vec<Node>) -> R<Node> {
        if ci.is_empty() { return self.tp(args); }
        let t = self.tensor(args)?;
        let ac = self.contract(t, ci)?;
        let mut ed = self.editor(&ac)?;
        ed.track_start();
        loop {
            let mut stop = true;
            let mut i = 0;
            while i < ed.args.len() {
                if !self.is_mat_opt(&ed.args[i].element) { i += 1; continue; }
                let idx = &ed.args[i].indices;
                let first = *idx.first().ok_or_else(|| AErr::Index("list index out of range".into()))?;
                let second = *idx.get(1).ok_or_else(|| AErr::Index("list index out of range".into()))?;
                let ff = ed.count_with(first);
                let sf = ed.count_with(second);
                if first.is_some() && ff == 1 && first == second {
                    stop = false;
                    let el = ed.args[i].element.clone().unwrap();
                    ed.args[i].element = Some(self.trace_norm(&el)?);
                    ed.args[i].indices = vec![];
                    break;
                }
                let mut scan = vec![];
                if ff == 2 { scan.push(first); }
                if sf == 2 { scan.push(second); }
                let (cand, transpose, found) = self.candidate(&scan, &ed.args[i + 1..])?;
                if let Some(cid) = cand {
                    stop = false;
                    let aid = ed.args[i].id;
                    ed.track_merge(aid, cid)?;
                    let transpose1 = Some(found) == first;
                    // _insert_candidate_into_editor
                    let cpos = ed.pos(cid)?;
                    let cel = ed.args[cpos].element.clone().unwrap();
                    let (other, other_index) = if transpose {
                        let t = self.mat_t(&cel, false)?;
                        (self.mat(t), ed.args[cpos].indices[0])
                    } else {
                        (cel, *ed.args[cpos].indices.get(1).ok_or_else(|| AErr::Index("list index out of range".into()))?)
                    };
                    let ael = ed.args[i].element.clone().unwrap();
                    let left = if transpose1 { let t = self.mat_t(&ael, false)?; self.mat(t) } else { ael };
                    let new_el = self.mm(&left, &other)?;
                    ed.remove(cid)?;
                    let r = self.rank(&new_el);
                    let mut new = ed.arg(Some(new_el), vec![None; r]);
                    new.indices = if Some(found) == first { vec![second, other_index] } else { vec![first, other_index] };
                    let mut set = new.indices.clone();
                    set.sort_unstable();
                    set.dedup();
                    if set.len() == 1 && set[0].is_some() {
                        let el = new.element.clone().unwrap();
                        new.element = Some(self.trace_norm(&el)?);
                        new.indices = vec![];
                    }
                    ed.args[i] = new;
                    break;
                }
                i += 1;
            }
            if stop { break; }
        }
        ed.refresh_indices();
        self.to_array(&mut ed)
    }

    /// `_find_trivial_matrices_rewrite`.
    fn trivial_matrices_rewrite(&mut self, args: &[Node]) -> R<(Node, Vec<usize>)> {
        let mut trivial: Vec<MatrixExpr> = vec![];
        let mut pos: Option<(usize, Expr, Vec<MatrixExpr>, usize)> = None;
        let mut removed = vec![];
        let mut counter = 0;
        let mut kept: Vec<Option<Node>> = args.iter().cloned().map(Some).collect();
        for (i, arg) in args.iter().enumerate() {
            if let Some(m) = self.as_mat(arg) {
                let s = self.shape_r(arg)?;
                if self.is_11(&s) {
                    trivial.push(m.clone());
                    kept[i] = None;
                    removed.extend([counter, counter + 1]);
                } else if pos.is_none() {
                    if let MatrixKind::Mul(c, fs) = m.kind() {
                        for (j, e) in fs.iter().enumerate() {
                            if e.cols().is_one() {
                                pos = Some((i, c.clone(), fs.clone(), j));
                                break;
                            }
                        }
                    }
                }
            }
            counter += self.rank(arg);
        }
        let Some((p, c, fs, j)) = pos else { return Ok((Node::Tensor(args.to_vec()), vec![])) };
        let mut all: Vec<MatrixExpr> = fs[..=j].to_vec();
        all.extend(trivial);
        all.extend_from_slice(&fs[j + 1..]);
        kept[p] = Some(self.mul_matrices(c, &all)?);
        let t = self.tensor(kept.into_iter().flatten().collect())?;
        Ok((t, removed))
    }

    /// `_find_trivial_kronecker_products_broadcast`.
    fn kronecker_broadcast(&mut self, args: &[Node]) -> R<(Node, Vec<usize>)> {
        let mut newargs: Vec<Node> = vec![];
        let mut removed: Vec<usize> = vec![];
        let mut count = 0usize;
        for arg in args {
            count += self.rank(arg);
            let shape = self.shape_r(arg)?;
            let current: Vec<usize> = (1..=shape.len()).rev().map(|i| count - i).collect();
            let last_shape = match newargs.last() { Some(l) => Some(self.shape_r(l)?), None => None };
            if self.is_11(&shape) && last_shape.as_ref().is_some_and(|s| !self.has_one(s))
                && self.is_mat(newargs.last().unwrap()) && self.is_mat(arg) {
                let (a, b) = (self.as_mat(newargs.last().unwrap()).unwrap().clone(), self.as_mat(arg).unwrap().clone());
                let k = MatrixExpr::kronecker(&[a, b], false)?;
                *newargs.last_mut().unwrap() = self.mat(k);
                removed.extend(current);
            } else if !self.has_one(&shape) && last_shape.as_ref().is_some_and(|s| self.is_11(s)) {
                let (Some(a), Some(b)) = (self.as_mat(newargs.last().unwrap()).cloned(), self.as_mat(arg).cloned()) else {
                    return unsupported("KroneckerProduct of an array expression");
                };
                let k = MatrixExpr::kronecker(&[a, b], false)?;
                *newargs.last_mut().unwrap() = self.mat(k);
                let lo = *current.iter().min().ok_or_else(|| AErr::Value("min() arg is an empty sequence".into()))?;
                let prev: Vec<usize> = (0..lo).filter(|i| !removed.contains(i)).collect();
                removed.extend_from_slice(&prev[prev.len().saturating_sub(2)..]);
            } else {
                newargs.push(arg.clone());
            }
        }
        Ok((self.tensor(newargs)?, removed))
    }

    /// `ArrayContraction.flatten_contraction_of_diagonal`.
    pub fn flatten_contraction_of_diagonal(&mut self, expr: Node) -> R<Node> {
        match &expr {
            Node::Contract(e, _) if matches!(**e, Node::Diagonal(..)) => {
                Err(AErr::Type("ArrayDiagonal._push_indices_down() missing 1 required positional argument: 'rank'".into()))
            },
            Node::Contract(..) => Ok(expr),
            _ => Err(AErr::Attribute(format!("'{}' object has no attribute 'flatten_contraction_of_diagonal'", self.class(&expr)?))),
        }
    }

    /// `ArrayContraction.split_multiple_contractions`.
    pub fn split_multiple_contractions(&mut self, expr: Node) -> R<Node> {
        let Node::Contract(inner, ci) = &expr else {
            return Err(AErr::Attribute(format!("'{}' object has no attribute 'split_multiple_contractions'", self.class(&expr)?)));
        };
        let ci = ci.clone();
        let inner_shape = self.shape_r(inner)?;
        let mut ed = self.editor(&expr)?;
        let mut onearray_insert: Vec<usize> = vec![];
        for (indl, links) in ci.iter().enumerate() {
            if links.len() <= 2 { continue; }
            let positions = ed.mapping_for_index(indl as i64)?;
            let current_dim = inner_shape[links[0]];
            let mut not_vectors: Vec<(usize, usize)> = vec![];
            let mut vectors: Vec<(usize, usize)> = vec![];
            for (arg_ind, rel) in positions {
                let a = &ed.args[arg_ind];
                let Some(mat) = &a.element else { return unsupported("removed editor element") };
                let ms = self.shape_r(mat)?;
                let (start, _) = ed.abs_range(a.id)?;
                let other_abs = start as i64 + 1 - rel as i64;
                let elsewhere = ci.iter().enumerate().any(|(li, l)| li != indl && l.iter().any(|&x| x as i64 == other_abs));
                if !self.has_one(&ms) || (self.ones[current_dim] && !self.is_11(&ms)) || elsewhere {
                    not_vectors.push((a.id, rel));
                } else {
                    vectors.push((a.id, rel));
                }
            }
            if not_vectors.len() > 2 { continue; }
            for &(v, _) in &vectors {
                let p = ed.pos(v)?;
                let el = ed.args[p].element.clone().unwrap();
                ed.args[p].element = Some(self.diag_matrix(&el, true)?);
            }
            let mut to_loop: Vec<(usize, usize)> = not_vectors.iter().take(1).copied().collect();
            to_loop.extend(vectors.iter().copied());
            to_loop.extend(not_vectors.iter().skip(1).copied());
            let (first, rel0) = to_loop[0];
            let mut new_index = ed.get(first)?.indices[rel0];
            if to_loop.len() > 2 {
                for &(v, rel) in &to_loop[1..to_loop.len() - 1] {
                    ed.get(v)?.indices[rel] = new_index;
                    new_index = Some(ed.new_index());
                    let a = ed.get(v)?;
                    let none = a.indices.iter().position(Option::is_none).ok_or_else(|| AErr::Value("None is not in list".into()))?;
                    if none as i64 != 1 - rel as i64 { return Err(AErr::Assertion(String::new())); }
                    a.indices[none] = new_index;
                    onearray_insert.push(v);
                }
            }
            let (last, rel) = *to_loop.last().unwrap();
            ed.get(last)?.indices[rel] = new_index;
        }
        for v in onearray_insert {
            let d1 = self.dim(&Expr::one());
            let one = self.one(vec![d1]);
            let a = ed.arg(Some(one), vec![None]);
            ed.insert_after(v, a)?;
        }
        self.to_array(&mut ed)
    }

    /// `identify_hadamard_products`.
    pub fn identify_hadamard_products(&mut self, expr: &Node) -> R<Node> {
        let mut ed = self.editor(expr)?;
        let mut contr_to_args: Vec<(Vec<i64>, Vec<usize>)> = vec![];
        let mut ind_count: FxHashMap<Option<i64>, usize> = FxHashMap::default();
        for a in &ed.args {
            for i in &a.indices { *ind_count.entry(*i).or_default() += 1; }
            if a.indices.contains(&None) { continue; }
            let mut k: Vec<i64> = a.indices.iter().flatten().copied().collect();
            k.sort_unstable();
            k.dedup();
            match contr_to_args.iter_mut().find(|(x, _)| *x == k) { Some((_, v)) => v.push(a.id), None => contr_to_args.push((k, vec![a.id])) }
        }
        let check_transpose = |x: &[Option<i64>]| -> bool {
            let y: Vec<i64> = x.iter().map(|i| { let i = i.unwrap_or(0); if i >= 0 { i } else { -1 - i } }).collect();
            y.windows(2).all(|w| w[0] <= w[1])
        };
        let count = |i: i64| *ind_count.get(&Some(i)).unwrap_or(&0);
        for (k, v) in &contr_to_args {
            let mut make_trace = false;
            let mut first_element: Option<Node> = None;
            if k.len() == 1 && k[0] >= 0 && contr_to_args.iter().filter(|(x, _)| x.contains(&k[0])).count() == 1 {
                make_trace = true;
                first_element = None;
            } else if k.len() != 2 {
                continue;
            }
            if v.len() == 1 { continue; }
            let elems: Vec<(Node, Vec<Option<i64>>)> = v.iter().map(|&id| {
                let a = &ed.args[ed.pos(id)?];
                Ok((a.element.clone().ok_or_else(|| AErr::Unsupported("removed editor element".into()))?, a.indices.clone()))
            }).collect::<R<_>>()?;
            let factors: Vec<(Node, Vec<Option<i64>>)>;
            if k.iter().all(|&j| count(j) == v.len() && j >= 0) {
                make_trace = true;
                let mut fe = elems[0].0.clone();
                if !check_transpose(&elems[0].1) { let t = self.attr_t(&fe)?; fe = self.mat(t); }
                first_element = Some(fe);
                factors = elems[1..].to_vec();
            } else {
                factors = elems.clone();
            }
            if factors.len() > 1 {
                // `Transpose(x)` of an array leaf is a matrix expression; an
                // operand that is no matrix fails `HadamardProduct`'s check.
                let leafy = |e: &Node| matches!(e, Node::Leaf(i) if matches!(self.srcs[*i], Src::One(_) | Src::Opaque));
                let built = factors.iter().all(|(e, idx)| self.is_matexpr(e) || check_transpose(idx) || leafy(e));
                if built && factors.iter().any(|(e, idx)| !self.is_matexpr(e) && check_transpose(idx)) {
                    return Err(AErr::Type("Mix of Matrix and Scalar symbols".into()));
                }
            }
            let mut ms = vec![];
            for (e, idx) in &factors {
                let Some(m) = self.as_mat(e).cloned() else {
                    if factors.len() == 1 && check_transpose(idx) { ms.push(None); continue; }
                    return unsupported("Hadamard product of an array expression");
                };
                ms.push(Some(if check_transpose(idx) { m } else { m.transpose(false)? }));
            }
            let mut hp = if ms.len() == 1 && ms[0].is_none() { factors[0].0.clone() } else {
                let h = self.hadamard(ms.into_iter().flatten().collect())?;
                self.mat(h)
            };
            let mut hp_indices = elems[0].1.clone();
            if !check_transpose(&factors[0].1) { hp_indices.reverse(); }
            if make_trace {
                let hpt = self.attr_t(&hp)?;
                let hpt = self.mat(hpt);
                let prod = match &first_element {
                    Some(fe) => self.mm(fe, &hpt)?,
                    None => { let one = self.scal(Expr::one()); self.mm(&one, &hpt)? },
                };
                hp = self.trace_norm(&prod)?;
                hp_indices = vec![];
            }
            let new = ed.arg(Some(hp), hp_indices);
            ed.insert_after(v[0], new)?;
            for &id in v { ed.remove(id)?; }
        }
        self.to_array(&mut ed)
    }

    /// `identify_removable_identity_matrices`.
    pub fn identify_removable_identity_matrices(&mut self, expr: &Node) -> R<Node> {
        let mut ed = self.editor(expr)?;
        let mut flag = true;
        while flag {
            flag = false;
            let mut i = 0;
            while i < ed.args.len() {
                let a = ed.args[i].clone();
                if self.is_identity(&a.element) {
                    if a.indices == [None, None] { i += 1; continue; }
                    if a.indices.contains(&None) {
                        let ind = a.indices.iter().flatten().next().copied();
                        let counted = ed.count_with(ind);
                        if counted == 1 {
                            // SymPy inserts a bare OneArray into the editor
                            // and fails on it at the next scan.
                            return Err(AErr::Attribute("'OneArray' object has no attribute 'element'".into()));
                        }
                    } else if a.indices[0] == a.indices[1] {
                        let counted = ed.count_with(a.indices[0]);
                        if counted > 1 {
                            ed.remove(a.id)?;
                            flag = true;
                            break;
                        }
                    }
                } else if self.ask_diagonal(&a.element) {
                    if a.indices == [None, None] { i += 1; continue; }
                    if !a.indices.contains(&None) && a.indices.first() == a.indices.get(1) {
                        let ind = a.indices[0];
                        if ed.count_with(ind) == 3 {
                            let ind_new = ed.new_index();
                            let others: Vec<usize> = (0..ed.args.len()).filter(|&j| ed.args[j].id != a.id).collect();
                            let o = *others.get(1).ok_or_else(|| AErr::Index("list index out of range".into()))?;
                            ed.args[o].indices = ed.args[o].indices.iter().map(|&j| if j == ind { Some(ind_new) } else { j }).collect();
                            ed.args[i].indices = vec![ind, Some(ind_new)];
                            flag = true;
                            break;
                        }
                    }
                }
                i += 1;
            }
        }
        self.to_array(&mut ed)
    }

    /// `remove_identity_matrices`.
    pub fn remove_identity_matrices(&mut self, expr: &Node) -> R<(Node, Vec<usize>)> {
        let mut ed = self.editor(expr)?;
        let mut removed: Vec<usize> = vec![];
        let mut permutation_map: Vec<(usize, usize)> = vec![];
        let mut free_map: FxHashMap<usize, usize> = FxHashMap::default();
        let mut acc = 0;
        for a in &ed.args {
            free_map.insert(a.id, acc);
            acc += a.indices.iter().filter(|i| i.is_none()).count();
        }
        let mut update_pairs: Vec<((usize, i64), Vec<Option<i64>>)> = vec![];
        let mut non_identity_leak: Option<usize> = None;
        for ind in 0..ed.n_contr {
            let args = ed.args_with(Some(ind));
            let ids: Vec<usize> = args.iter().copied().filter(|&id| self.is_identity(&ed.args[ed.pos(id).unwrap()].element)).collect();
            let n = ids.len();
            if n + 1 != args.len() || n == 0 { continue; }
            let non_identity = *args.iter().find(|id| !ids.contains(id)).unwrap();
            non_identity_leak = Some(non_identity);
            if ids.iter().any(|&id| !ed.args[ed.pos(id).unwrap()].indices.contains(&None)) { continue; }
            for &id in &ids {
                let p = ed.pos(id)?;
                ed.args[p].element = None;
                let nf = ed.args[p].indices.iter().filter(|i| i.is_none()).count();
                let s = free_map[&id];
                removed.extend(s..s + nf);
            }
            let last_removed = removed.pop().ok_or_else(|| AErr::Index("pop from empty list".into()))?;
            let p = ed.pos(non_identity)?;
            let key = (last_removed, ind);
            let val = ed.args[p].indices.clone();
            match update_pairs.iter_mut().find(|(k, _)| *k == key) { Some((_, v)) => *v = val, None => update_pairs.push((key, val)) }
            ed.args[p].indices = ed.args[p].indices.iter().map(|&i| if i == Some(ind) { None } else { i }).collect();
        }
        removed.sort_unstable();
        let rank = self.rank(expr);
        let shifts = accumulate((0..rank).map(|i| usize::from(removed.contains(&i))));
        for ((last_removed, ind), nii) in &update_pairs {
            let base = free_map[&non_identity_leak.unwrap()];
            let pos: Vec<usize> = nii.iter().enumerate().filter(|(_, e)| **e == Some(*ind)).map(|(i, _)| base + i).collect();
            if pos.len() != 1 { return Err(AErr::Assertion(String::new())); }
            for j in pos {
                match permutation_map.iter_mut().find(|(k, _)| *k == j) { Some((_, v)) => *v = *last_removed, None => permutation_map.push((j, *last_removed)) }
            }
        }
        ed.args.retain(|a| a.element.is_some());
        let ret = self.to_array(&mut ed)?;
        let mut permutation: Vec<i64> = vec![];
        let (mut counter, mut counter2) = (0usize, 0usize);
        for j in 0..rank {
            if removed.contains(&j) { continue; }
            if let Some(&(_, target)) = permutation_map.iter().find(|(k, _)| *k == counter2) {
                let s = *shifts.get(target).ok_or_else(|| AErr::Index("list index out of range".into()))?;
                permutation.push(target as i64 - s as i64);
                counter2 += 1;
            } else {
                while permutation_map.iter().any(|(_, v)| *v == counter) { counter += 1; }
                permutation.push(counter as i64);
                counter += 1;
                counter2 += 1;
            }
        }
        let p = af_invert(&permutation)?;
        Ok((self.permute_checked(ret, p)?, removed))
    }

    /// `_array_contraction_to_diagonal_multiple_identity`.
    pub fn contraction_to_diagonal_multiple_identity(&mut self, expr: &Node) -> R<(Node, Vec<usize>)> {
        let mut ed = self.editor(expr)?;
        ed.track_start();
        let mut removed: Vec<usize> = vec![];
        let mut diag_counter = 0i64;
        for i in 0..ed.n_contr {
            let mut identities = vec![];
            let mut args = vec![];
            for a in &ed.args {
                if !a.indices.contains(&Some(i)) { continue; }
                if self.is_identity(&a.element) { identities.push(a.id); } else { args.push(a.id); }
            }
            if identities.is_empty() || args.len() + identities.len() < 3 { continue; }
            let new_diag = -1 - diag_counter;
            diag_counter += 1;
            let id1 = identities[0];
            if !ed.get(id1)?.indices.contains(&None) { continue; }
            let (lo, hi) = ed.abs_free_range(id1)?;
            if lo >= hi { return Err(AErr::Index("list index out of range".into())); }
            ed.track.as_mut().unwrap().last_mut().unwrap().push(lo as i64);
            ed.get(id1)?.element = None;
            for &id in &identities[1..] {
                ed.get(id)?.element = None;
                let (lo, hi) = ed.abs_free_range(id)?;
                removed.extend(lo..hi);
            }
            for &id in &args {
                let a = ed.get(id)?;
                a.indices = a.indices.iter().map(|&j| if j == Some(i) { Some(new_diag) } else { j }).collect();
            }
        }
        let track = ed.track.take().unwrap();
        let mut kept: Vec<Vec<i64>> = vec![];
        for (j, t) in track.into_iter().enumerate() {
            if ed.args.get(j).is_some_and(|a| a.element.is_none()) { continue; }
            kept.push(t);
        }
        let mut all: Vec<i64> = kept.iter().flatten().copied().collect();
        all.sort_unstable();
        all.dedup();
        ed.track = Some(kept.iter().map(|t| t.iter().map(|j| all.binary_search(j).unwrap() as i64).collect()).collect());
        ed.args.retain(|a| a.element.is_some());
        Ok((self.to_array(&mut ed)?, removed))
    }

    /// `_array_diag2contr_diagmatrix`.
    /// The public `_array_diag2contr_diagmatrix(expr)`: SymPy reads
    /// `expr.expr` and, when it is a tensor product,
    /// `expr.diagonal_indices`, so an argument without these attributes is
    /// an AttributeError.
    pub fn diag2contr_entry(&mut self, expr: Node) -> R<Node> {
        match &expr {
            Node::Diagonal(..) => self.diag2contr_diagmatrix(expr),
            Node::Contract(inner, ..) | Node::Permute(inner, ..) | Node::Reshape(inner, ..) | Node::Apply(inner, _) | Node::PowApply(inner, _) => {
                if matches!(**inner, Node::Tensor(_)) {
                    let c = self.class(&expr)?;
                    return Err(AErr::Attribute(format!("'{c}' object has no attribute 'diagonal_indices'")));
                }
                Ok(expr)
            },
            _ => {
                // Python's AttributeError names the class of the object.
                let c: String = match &expr {
                    Node::Leaf(i) => match &self.srcs[*i] {
                        Src::Matrix(m) => m.kind_name().to_string(),
                        Src::Scalar(e) => e.class_name().to_string(),
                        _ => self.class(&expr)?.to_string(),
                    },
                    _ => self.class(&expr)?.to_string(),
                };
                Err(AErr::Attribute(format!("'{c}' object has no attribute 'expr'")))
            },
        }
    }

    pub fn diag2contr_diagmatrix(&mut self, expr: Node) -> R<Node> {
        let Node::Diagonal(inner, groups) = &expr else { return Ok(expr) };
        let Node::Tensor(targs) = &**inner else { return Ok(expr) };
        let mut args = targs.clone();
        let mut diag: Vec<Option<Vec<usize>>> = groups.iter().cloned().map(Some).collect();
        let map = mapping(&args.iter().map(|a| self.subrank_of(a)).collect::<Vec<_>>());
        let links: Vec<Vec<(usize, usize)>> = groups.iter().map(|g| g.iter().map(|&j| map.get(j).copied().ok_or_else(|| AErr::Key(j.to_string()))).collect::<R<_>>()).collect::<R<_>>()?;
        let mut contr: Vec<[(usize, usize); 2]> = vec![];
        let mut replaced = vec![false; args.len()];
        for (i, rel) in links.iter().enumerate() {
            if groups[i].len() != 2 { continue; }
            let ((p1o, p1i), (p2o, p2i)) = (rel[0], rel[1]);
            let (arg1, arg2) = (args[p1o].clone(), args[p2o].clone());
            if self.rank(&arg1) != 2 || self.rank(&arg2) != 2 {
                if replaced[p1o] || replaced[p2o] { diag[i] = None; }
                continue;
            }
            let (s1, s2) = (self.shape_r(&arg1)?, self.shape_r(&arg2)?);
            let pos1_in2 = pyidx(2, 1 - p1i as i64, "tuple")?;
            let pos2_in2 = pyidx(2, 1 - p2i as i64, "tuple")?;
            let (p1i_, p2i_) = (pyidx(2, p1i as i64, "tuple")?, pyidx(2, p2i as i64, "tuple")?);
            if self.ones[s1[pos1_in2]] {
                let d = if !self.ones[s1[p1i_]] { self.diag_matrix(&arg1, false)? } else { arg1.clone() };
                args.push(d);
                contr.push([(p2o, p2i), (args.len() - 1, p1i)]);
                diag[i] = None;
                args[p1o] = self.one(vec![s1[pos1_in2]]);
                replaced[p1o] = true;
            } else if self.ones[s2[pos2_in2]] {
                let d = if !self.ones[s2[p2i_]] { self.diag_matrix(&arg2, false)? } else { arg2.clone() };
                args.push(d);
                contr.push([(p1o, p1i), (args.len() - 1, p2i)]);
                diag[i] = None;
                args[p2o] = self.one(vec![s2[pos2_in2]]);
                replaced[p2o] = true;
            }
        }
        let diag_new: Vec<Vec<usize>> = diag.into_iter().flatten().collect();
        let mut cumul = vec![0];
        for a in &args { cumul.push(cumul[cumul.len() - 1] + self.rank(a)); }
        let contr2: Vec<Vec<usize>> = contr.iter().map(|c| c.iter().map(|&(a, b)| cumul[a] + b).collect()).collect();
        let t = self.tensor(args)?;
        let tc = self.contract(t, contr2)?;
        self.diagonal(tc, diag_new, false)
    }

    /// `_array2matrix`.
    pub fn a2m(&mut self, expr: Node) -> R<Node> {
        match expr {
            Node::Zero(ref s) if s.len() == 2 => {
                let m = MatrixExpr::zero(self.dims[s[0]].clone(), self.dims[s[1]].clone())?;
                Ok(self.mat(m))
            },
            Node::Tensor(args) => {
                let conv = args.into_iter().map(|a| self.a2m(a)).collect::<R<Vec<_>>>()?;
                self.tp(conv)
            },
            Node::Add(args) => {
                let conv = args.into_iter().map(|a| self.a2m(a)).collect::<R<Vec<_>>>()?;
                self.madd(conv)
            },
            Node::Contract(..) => self.a2m_contraction(expr),
            Node::Diagonal(ref inner, ref groups) => {
                let sub = self.a2m((**inner).clone())?;
                let p = self.diagonal(sub, groups.clone(), false)?;
                let mut p = self.identify_hadamard_products(&p)?;
                if matches!(p, Node::Diagonal(..)) { p = self.diag2contr_diagmatrix(p)?; }
                if p == expr { return Ok(expr); }
                self.a2m(p)
            },
            Node::Permute(..) => self.a2m_permute(expr),
            Node::Apply(inner, f) => {
                let sub = self.a2m(*inner)?;
                let Some(m) = self.as_mat(&sub).cloned() else { return Ok(Node::Apply(Box::new(sub), f)) };
                let Src::Fn(func) = &self.srcs[f] else { return unsupported("element-wise function") };
                let func = func.clone();
                let s = self.shape_r(&sub)?;
                if !self.is_11(&s) {
                    let Some((vars, body)) = crate::sets::lambda_parts(&func) else { return unsupported("element-wise function") };
                    let Some(d) = vars.first() else { return unsupported("element-wise function") };
                    if let Some((w, p)) = match_power(&body, d) {
                        let hp = MatrixExpr::hadamard_power(MatrixOperand::Matrix(m), MatrixOperand::Scalar(p), false)?;
                        let hp = self.mat(hp);
                        let w = self.scal(w);
                        return self.mm(&w, &hp);
                    }
                }
                let a = MatrixExpr::apply_function(func, &m, false)?;
                Ok(self.mat(a))
            },
            Node::PowApply(..) => unsupported("element-wise power"),
            other => Ok(other),
        }
    }

    fn a2m_contraction(&mut self, expr: Node) -> R<Node> {
        let expr = self.flatten_contraction_of_diagonal(expr)?;
        let expr = self.identify_removable_identity_matrices(&expr)?;
        let expr = self.split_multiple_contractions(expr)?;
        let expr = self.identify_hadamard_products(&expr)?;
        let Node::Contract(orig_sub, ci) = &expr else { return self.a2m(expr) };
        let mut subexpr = (**orig_sub).clone();
        let mut ci = ci.clone();
        let is = |c: &[Vec<usize>], v: &[&[usize]]| c.len() == v.len() && c.iter().zip(v).all(|(a, b)| a.as_slice() == *b);
        let special = is(&ci, &[&[0], &[1]]) || {
            if is(&ci, &[&[0]]) { let s = self.shape_r(&subexpr)?; self.ones[s[pyidx(s.len(), 1, "tuple")?]] }
            else if is(&ci, &[&[1]]) { let s = self.shape_r(&subexpr)?; self.ones[s[pyidx(s.len(), 0, "tuple")?]] }
            else { false }
        };
        if special {
            let shape = self.shape_r(&subexpr)?;
            let (r0, r1) = (pyidx(shape.len(), 0, "tuple")?, pyidx(shape.len(), 1, "tuple")?);
            subexpr = self.a2m(subexpr)?;
            if self.is_mat(&subexpr) {
                let one = self.dims.iter().position(Expr::is_one).map(|i| self.dims[i].clone()).unwrap_or_else(Expr::one);
                let left = MatrixExpr::one(one.clone(), self.dims[shape[r0]].clone())?;
                let right = MatrixExpr::one(self.dims[shape[r1]].clone(), one)?;
                let (left, right) = (self.mat(left), self.mat(right));
                let t = self.mm(&left, &subexpr)?;
                return self.mm(&t, &right);
            }
        }
        if matches!(subexpr, Node::Tensor(_)) {
            let conv = self.a2m(subexpr)?;
            let mut newexpr = self.contract(conv, ci.clone())?;
            let Node::Contract(_, g) = &newexpr else {
                return Err(AErr::Attribute(format!("'{}' object has no attribute 'contraction_indices'", self.class(&newexpr)?)));
            };
            ci = g.clone();
            if self.subranks(&newexpr).iter().any(|&r| r > 2) {
                let Node::Tensor(oargs) = &**orig_sub else { return unsupported("contraction of a non tensor product") };
                let choices: Vec<Vec<Node>> = oargs.iter().map(|a| match a { Node::Add(x) => x.clone(), o => vec![o.clone()] }).collect();
                let mut combos: Vec<Vec<Node>> = vec![vec![]];
                for c in &choices {
                    let mut next = vec![];
                    for prefix in &combos { for x in c { let mut p = prefix.clone(); p.push(x.clone()); next.push(p); } }
                    combos = next;
                }
                let mut addends = vec![];
                for j in combos { addends.push(self.tp(j)?); }
                let added = self.array_add(addends)?;
                newexpr = self.contract(added, ci.clone())?;
            }
            if matches!(newexpr, Node::Add(_)) { return self.a2m(newexpr); }
            let Node::Contract(inner, _) = &newexpr else { return Err(AErr::Assertion(String::new())) };
            let args = match &**inner {
                Node::Tensor(a) | Node::Add(a) => a.clone(),
                n => match self.as_mat(n).map(|m| m.kind().clone()) {
                    Some(MatrixKind::Mul(c, fs)) => {
                        let mut v = vec![];
                        if !c.is_one() { v.push(self.scal(c)); }
                        for f in fs { v.push(self.mat(f)); }
                        v
                    },
                    _ => return unsupported("contraction of a non tensor product"),
                },
            };
            return self.tp1(ci, args);
        }
        if !codegen(&subexpr) {
            let ret = self.a2m(subexpr)?;
            let Node::Contract(_, oci) = &expr else { unreachable!() };
            if self.is_mat(&ret) {
                if !is(oci, &[&[0, 1]]) { return Err(AErr::Assertion(String::new())); }
                return self.trace(ret);
            }
            return self.contract(ret, oci.clone());
        }
        unsupported("contraction recognized as None")
    }

    fn a2m_permute(&mut self, expr: Node) -> R<Node> {
        let Node::Permute(inner, perm) = &expr else { unreachable!() };
        if perm.as_slice() == [1, 0] {
            let sub = self.a2m((**inner).clone())?;
            return self.tr(sub);
        }
        match &**inner {
            Node::Tensor(args) => {
                let ranks = self.subranks(inner);
                let inv = ax::invert(perm);
                let total: usize = ranks.iter().sum();
                let newrange: Vec<usize> = (0..total).map(|i| inv.get(i).copied().ok_or_else(|| AErr::Index("list index out of range".into()))).collect::<R<_>>()?;
                let mut newpos = vec![];
                let mut c = 0;
                for &r in &ranks { newpos.push(newrange[c..c + r].to_vec()); c += r; }
                let (mut newargs, mut newperm, mut scalars) = (vec![], vec![], vec![]);
                for (pos, arg) in newpos.iter().zip(args.iter()) {
                    if pos.is_empty() {
                        scalars.push(self.a2m(arg.clone())?);
                    } else if pos.windows(2).all(|w| w[0] <= w[1]) {
                        newargs.push(self.a2m(arg.clone())?);
                        newperm.extend(pos.iter().map(|&x| x as i64));
                    } else if pos.len() == 2 {
                        let a = self.a2m(arg.clone())?;
                        newargs.push(self.tr(a)?);
                        newperm.extend(pos.iter().rev().map(|&x| x as i64));
                    } else {
                        return Err(AErr::NotImplemented(String::new()));
                    }
                }
                scalars.extend(newargs);
                let t = self.tp(scalars)?;
                let p = af_invert(&newperm)?;
                self.permute_checked(t, p)
            },
            Node::Contract(..) => {
                let mml = self.a2m((**inner).clone())?;
                let Node::Tensor(margs) = &mml else { return self.permute_checked(mml, perm.clone()) };
                let margs = margs.clone();
                let n = 2 * margs.len();
                let size = n.max(perm.len());
                let full: Vec<usize> = (0..size).map(|i| if i < perm.len() { perm[i] } else { i }).collect();
                let mut out = vec![];
                for i in 0..margs.len() {
                    let (p1, p2) = (full[2 * i], full[2 * i + 1]);
                    if p1 / 2 != p2 / 2 { return self.permute_checked(mml, full); }
                    let a = margs.get(p1 / 2).cloned().ok_or_else(|| AErr::Index("tuple index out of range".into()))?;
                    out.push(if p1 > p2 { self.tr(a)? } else { a });
                }
                self.tp(out)
            },
            _ => Ok(expr),
        }
    }

    /// `_remove_trivial_dims`.
    pub fn rtd(&mut self, expr: Node) -> R<(Node, Vec<usize>)> {
        match expr {
            Node::Tensor(ref args) => self.rtd_tensor(args.clone()),
            Node::Add(ref args) => {
                let mut newargs = vec![];
                let mut removed = vec![];
                for a in args { let (n, r) = self.rtd(a.clone())?; newargs.push(n); removed.push(r); }
                let mut shapes: Vec<Option<Vec<usize>>> = newargs.iter().map(|a| self.shape(a)).collect();
                shapes.sort();
                shapes.dedup();
                if shapes.len() > 1 { return Ok((expr, vec![])); }
                if removed.is_empty() { return Ok((expr, vec![])); }
                let r0 = removed.swap_remove(0);
                Ok((self.madd(newargs)?, r0))
            },
            Node::Permute(ref inner, ref p) => {
                let (sub, subremoved) = self.rtd((**inner).clone())?;
                let pinv = ax::invert(p);
                let shift = accumulate((0..p.len()).map(|i| usize::from(subremoved.contains(&i))));
                let mut premoved: Vec<usize> = subremoved.iter().map(|&i| pinv.get(i).copied().ok_or_else(|| AErr::Index("list index out of range".into()))).collect::<R<_>>()?;
                let p2: Vec<usize> = p.iter().filter(|e| !subremoved.contains(e)).map(|&e| e - shift[e]).collect();
                let mut newexpr = self.permute_checked(sub, p2)?;
                premoved.sort_unstable();
                if newexpr != expr {
                    let m = self.a2m(newexpr)?;
                    let (n, r2) = self.rtd(m)?;
                    newexpr = n;
                    premoved = combine_removed(&premoved, &r2);
                }
                Ok((newexpr, premoved))
            },
            Node::Contract(..) => self.rtd_contraction(expr),
            Node::Diagonal(ref inner, ref groups) => {
                let (newexpr, mut removed) = self.rtd((**inner).clone())?;
                let rank = self.rank(inner);
                let mut shifts = vec![0];
                shifts.extend(accumulate((0..rank).map(|i| usize::from(removed.contains(&i)))));
                let new_map: Vec<(Vec<usize>, Vec<usize>)> = groups.iter().map(|g| (g.clone(), g.iter().copied().filter(|j| !removed.contains(j)).collect())).collect();
                for (old, new) in &new_map { if new.len() == 1 { removed.retain(|i| !old.contains(i)); } }
                let mut new_diag: Vec<Vec<usize>> = new_map.iter().map(|(_, n)| n.iter().map(|&j| j - shifts[j]).collect()).collect();
                let up = ax::diagonal_push_up(groups, &removed, rank);
                let mut removed: Vec<usize> = up.into_iter().map(|x| x.ok_or_else(|| AErr::Unsupported("diagonalized axis removed".into()))).collect::<R<_>>()?;
                removed.sort_unstable();
                removed.dedup();
                new_diag.retain(|g| !g.is_empty());
                let newexpr2 = if !new_diag.is_empty() { self.diagonal(newexpr, new_diag, true)? } else { newexpr };
                if matches!(newexpr2, Node::Diagonal(..)) {
                    let (n3, r2) = self.remove_diagonalized_identity_matrices(&newexpr2)?;
                    return Ok((n3, combine_removed(&removed, &r2)));
                }
                Ok((newexpr2, removed))
            },
            Node::Apply(inner, f) => {
                let (sub, r) = self.rtd(*inner)?;
                Ok((Node::Apply(Box::new(sub), f), r))
            },
            Node::Leaf(i) => {
                if let Src::Matrix(m) = &self.srcs[i] {
                    if let MatrixKind::ApplyFunction(f, sub) = m.kind() {
                        let (f, sub) = (f.clone(), sub.clone());
                        let subn = self.mat(sub);
                        let (subn, sub_removed) = self.rtd(subn)?;
                        let s = self.shape_r(&subn)?;
                        let Some(sm) = self.as_mat(&subn).cloned() else { return unsupported("element-wise function") };
                        if self.is_11(&s) {
                            // `expr.function(subexpr), removed + [0, 1]`: the
                            // function applied to the 1x1 matrix as a scalar.
                            let proxy = Expr::from_kind(Kind::Func(Func::MatrixProxy(std::sync::Arc::new(sm)), vec![]));
                            let Some(v) = crate::sets::apply_lambda(&f, &proxy) else { return unsupported("element-wise function of a 1x1 matrix") };
                            let mut removed = sub_removed;
                            removed.extend([0, 1]);
                            return Ok((self.scal(v), removed));
                        }
                        let a = MatrixExpr::apply_function(f, &sm, false)?;
                        return Ok((self.mat(a), vec![]));
                    }
                }
                Ok((expr, vec![]))
            },
            other => Ok((other, vec![])),
        }
    }

    fn rtd_tensor(&mut self, args: Vec<Node>) -> R<(Node, Vec<usize>)> {
        let mut removed: Vec<usize> = vec![];
        let mut newargs: Vec<Node> = vec![];
        let mut cumul = vec![0];
        for a in &args { cumul.push(cumul[cumul.len() - 1] + self.rank(a)); }
        let mut pending: Option<usize> = None;
        let mut prev_i = 0usize;
        for (i, arg) in args.iter().enumerate() {
            let current = cumul[i]..cumul[i + 1];
            if self.is_one_array(arg) { removed.extend(current); continue; }
            if !self.is_mat(arg) {
                let (r, rem) = self.rtd(arg.clone())?;
                removed.extend(rem);
                newargs.push(r);
                continue;
            }
            let shape = self.shape_r(arg)?;
            let m = self.as_mat(arg).unwrap().clone();
            if self.is_identity_attr(&m) && self.is_11(&shape) {
                removed.extend(current);
                continue;
            } else if self.is_11(&shape) {
                let (a, _) = self.rtd(arg.clone())?;
                match newargs.last().cloned() {
                    None => newargs.push(a),
                    Some(last) => {
                        let ls = self.shape(&last).unwrap_or_default();
                        if self.has_one(&ls) {
                            let l1 = ls[pyidx(ls.len(), 1, "tuple")?];
                            *newargs.last_mut().unwrap() = if self.ones[l1] { self.mm(&last, &a)? } else { self.mm(&a, &last)? };
                            removed.extend(current);
                        } else {
                            newargs.push(a);
                        }
                    },
                }
            } else if self.has_one(&shape) {
                let k = *shape.iter().find(|&&d| !self.ones[d]).ok_or_else(|| AErr::Index("list index out of range".into()))?;
                match pending {
                    None => { pending = Some(k); prev_i = i; newargs.push(arg.clone()); },
                    Some(pk) if pk == k => {
                        let mut prev = newargs.last().cloned().unwrap();
                        let ps = self.shape_r(&prev)?;
                        let d1;
                        if self.ones[ps[pyidx(ps.len(), 0, "tuple")?]] { d1 = cumul[prev_i]; prev = self.tr(prev)?; } else { d1 = cumul[prev_i] + 1; }
                        let (d2, a) = if self.ones[shape[1]] { (cumul[i] + 1, self.tr(arg.clone())?) } else { (cumul[i], arg.clone()) };
                        *newargs.last_mut().unwrap() = self.mm(&prev, &a)?;
                        pending = None;
                        removed.extend([d1, d2]);
                    },
                    Some(_) => { newargs.push(arg.clone()); pending = Some(k); prev_i = i; },
                }
            } else {
                newargs.push(arg.clone());
                pending = None;
            }
        }
        let mut newexpr = self.tp(newargs)?;
        removed.sort_unstable();
        let mut newremoved = removed;
        if let Node::Tensor(a) = &newexpr {
            let (n, r2) = self.trivial_matrices_rewrite(&a.clone())?;
            newexpr = n;
            newremoved = combine_removed(&newremoved, &r2);
        }
        if let Node::Tensor(a) = &newexpr {
            let (n, r2) = self.kronecker_broadcast(&a.clone())?;
            newexpr = n;
            newremoved = combine_removed(&newremoved, &r2);
        }
        Ok((newexpr, newremoved))
    }

    fn rtd_contraction(&mut self, expr: Node) -> R<(Node, Vec<usize>)> {
        let (new_expr, removed0) = self.contraction_to_diagonal_multiple_identity(&expr)?;
        if new_expr != expr {
            let m = self.a2m(new_expr)?;
            let (n2, r1) = self.rtd(m)?;
            return Ok((n2, combine_removed(&removed0, &r1)));
        }
        let (expr, removed1) = self.remove_identity_matrices(&expr)?;
        let Node::Contract(inner, ci) = &expr else {
            let (e2, r2) = self.rtd(expr)?;
            return Ok((e2, combine_removed(&removed1, &r2)));
        };
        let (newexpr, removed2) = self.rtd((**inner).clone())?;
        let rank = self.rank(inner);
        let shifts = accumulate((0..rank).map(|i| usize::from(removed2.contains(&i))));
        let mut new_ci: Vec<Vec<usize>> = ci.iter().map(|g| g.iter().copied().filter(|j| !removed2.contains(j)).collect()).collect();
        new_ci.retain(|g: &Vec<usize>| !g.is_empty());
        let flat: Vec<usize> = ci.iter().flatten().copied().collect();
        let removed2: Vec<usize> = removed2.into_iter().filter(|i| !flat.contains(i)).collect();
        let new_ci: Vec<Vec<usize>> = new_ci.iter().map(|g| g.iter().map(|&j| shifts.get(j).map(|s| j - s).ok_or_else(|| AErr::Index("list index out of range".into()))).collect::<R<_>>()).collect::<R<_>>()?;
        let up = ax::contraction_push_up(ci, &[removed2]);
        let removed2: Vec<usize> = up.into_iter().flatten().map(|x| x.ok_or_else(|| AErr::Unsupported("contracted axis removed".into()))).collect::<R<_>>()?;
        let removed = combine_removed(&removed1, &removed2);
        Ok((self.contract(newexpr, new_ci)?, removed))
    }

    /// `_remove_diagonalized_identity_matrices`.
    fn remove_diagonalized_identity_matrices(&mut self, expr: &Node) -> R<(Node, Vec<usize>)> {
        let Node::Diagonal(inner, groups) = expr else { return Err(AErr::Assertion(String::new())) };
        let mut ed = self.editor(expr)?;
        let ndiag = ed.n_diag() as i64;
        let mut removed: Vec<usize> = vec![];
        for i in 0..ed.args.len() {
            let a = ed.args[i].clone();
            if !self.is_identity(&a.element) { continue; }
            if !(a.indices.contains(&None) && a.indices.iter().any(|x| x.is_some_and(|x| x < 0))) { continue; }
            let diag_ind = a.indices.iter().flatten().next().copied().unwrap();
            if !(diag_ind < 0 && diag_ind >= -ndiag) { return Err(AErr::Key(diag_ind.to_string())); }
            let others: Vec<usize> = ed.args.iter().filter(|j| j.indices.contains(&Some(diag_ind)) && j.id != a.id).map(|j| j.id).collect();
            let other = *others.first().ok_or_else(|| AErr::Index("list index out of range".into()))?;
            let op = ed.pos(other)?;
            let Some(oel) = ed.args[op].element.clone() else { continue };
            if !self.is_mat(&oel) { continue; }
            let os = self.shape_r(&oel)?;
            if !self.has_one(&os) { continue; }
            let Some(none_index) = ed.args[op].indices.iter().position(Option::is_none) else { continue };
            ed.args[i].element = None;
            ed.args[op].element = Some(self.diag_matrix(&oel, false)?);
            let (start, _) = ed.abs_range(other)?;
            removed.push(start + none_index);
        }
        ed.args.retain(|a| a.element.is_some());
        let rank = self.rank(inner);
        let up = ax::diagonal_push_up(groups, &removed, rank);
        let removed: Vec<usize> = up.into_iter().map(|x| x.ok_or_else(|| AErr::Unsupported("diagonalized axis removed".into()))).collect::<R<_>>()?;
        Ok((self.to_array(&mut ed)?, removed))
    }

    /// `convert_array_to_matrix`.
    pub fn convert(&mut self, expr: Node) -> R<Node> {
        let m = self.a2m(expr)?;
        Ok(self.rtd(m)?.0)
    }

    // ------------------------------------------------- matrix-derivative chains

    /// Python's `a * b` of two matrix expressions (`MatMul(a, b).doit()`).
    fn times(a: &MatrixExpr, b: &MatrixExpr) -> R<MatrixExpr> {
        Ok(MatrixExpr::mul(Expr::one(), &[a.clone(), b.clone()], true)?.doit(true)?)
    }

    /// `_chain_array_to_matrix` (the matrix form of an array expression for
    /// the matrix-derivative lines, not `convert_array_to_matrix`).
    /// `symbols` holds the `MatrixSymbol` of each rank-two `ArraySymbol`
    /// leaf. [`Chain::FunctionMatrix`] is the `FunctionMatrix` of the input,
    /// which only the caller can build: below the top level it declines.
    pub fn chain(&mut self, expr: &Node, symbols: &FxHashMap<usize, MatrixExpr>) -> R<Chain> {
        if matches!(expr, Node::Contract(..) | Node::Tensor(_)) {
            if let Some(n) = self.matrix_chains(expr, symbols)? { return Ok(Chain::Done(n)); }
        }
        if self.rank(expr) != 2 { return Ok(Chain::Keep); }
        match expr {
            Node::Leaf(i) if matches!(self.srcs[*i], Src::Opaque) => {
                let Some(m) = symbols.get(i) else { return unsupported("MatrixSymbol of an array symbol") };
                Ok(Chain::Done(self.mat(m.clone())))
            },
            // A matrix leaf is no array expression: it is returned as it is.
            Node::Leaf(i) if !matches!(self.srcs[*i], Src::One(_)) => Ok(Chain::Keep),
            Node::Add(args) => {
                let mut ms = vec![];
                for a in args { ms.push(self.chain_matrix(a, symbols)?); }
                let m = MatrixExpr::add(&ms, false)?;
                Ok(Chain::Done(self.mat(m)))
            },
            Node::Permute(inner, p) if p.as_slice() == [1, 0] => {
                let m = self.chain_matrix(inner, symbols)?.transpose(true)?;
                Ok(Chain::Done(self.mat(m)))
            },
            Node::Contract(inner, groups) if matches!(**inner, Node::Tensor(_)) => {
                let Node::Tensor(factors) = &**inner else { unreachable!() };
                let n = factors.iter().filter(|a| self.rank(a) == 2).count();
                if groups.iter().cloned().eq((0..n.saturating_sub(1)).map(|i| vec![2 * i + 1, 2 * i + 2])) {
                    // `MatMul` of the converted factors, reached only when a
                    // factor is no matrix (`_matrix_chains` declined).
                    return unsupported("MatMul of array factors");
                }
                Ok(Chain::FunctionMatrix)
            },
            _ => Ok(Chain::FunctionMatrix),
        }
    }

    /// `_chain_array_to_matrix` of a rank-two operand that must be a matrix
    /// (a `FunctionMatrix` below the top level declines).
    fn chain_matrix(&mut self, e: &Node, symbols: &FxHashMap<usize, MatrixExpr>) -> R<MatrixExpr> {
        match self.chain(e, symbols)? {
            Chain::Done(n) => match self.as_mat(&n) { Some(m) => Ok(m.clone()), None => unsupported("matrix chain of a non-matrix") },
            Chain::Keep => match self.as_mat(e) { Some(m) => Ok(m.clone()), None => unsupported("matrix chain of a non-matrix") },
            Chain::FunctionMatrix => unsupported("FunctionMatrix of an array operand"),
        }
    }

    /// `_as_matrix_factor`: the matrix of a tensor factor, or `None`.
    fn as_matrix_factor(&mut self, f: &Node, symbols: &FxHashMap<usize, MatrixExpr>) -> R<Option<MatrixExpr>> {
        if let Some(m) = self.as_mat(f) { return Ok(Some(m.clone())); }
        let arraylike = match f {
            Node::Leaf(i) => matches!(self.srcs[*i], Src::One(_) | Src::Opaque),
            _ => true,
        };
        if !arraylike { return unsupported("matrix factor of an array product"); }
        Ok(match self.chain(f, symbols)? {
            Chain::Done(n) => self.as_mat(&n).cloned(),
            Chain::Keep => None,
            Chain::FunctionMatrix => return unsupported("FunctionMatrix of an array operand"),
        })
    }

    /// `_matrix_chains`: each connected set of contracted matrix axes is a
    /// chain (a matrix product whose ends are free axes) or a cycle (a
    /// trace). `None` when an operand is no matrix or the free axes of a
    /// chain are not consecutive.
    fn matrix_chains(&mut self, expr: &Node, symbols: &FxHashMap<usize, MatrixExpr>) -> R<Option<Node>> {
        let (inner, groups): (&Node, &[Vec<usize>]) = match expr { Node::Contract(e, g) => (&**e, g.as_slice()), _ => (expr, &[]) };
        let factors: Vec<Node> = match inner { Node::Tensor(a) => a.clone(), _ => vec![inner.clone()] };
        if groups.iter().any(|g| g.len() != 2) { return Ok(None); }
        let mut owner: Vec<(usize, usize)> = vec![];
        let mut matrices: Vec<(usize, MatrixExpr)> = vec![];
        let mut scalars: Vec<Expr> = vec![];
        for (k, f) in factors.iter().enumerate() {
            if self.rank(f) == 0 {
                let Some(s) = self.as_scal(f) else { return unsupported("scalar factor of an array product") };
                scalars.push(s.clone());
                continue;
            }
            let Some(m) = self.as_matrix_factor(f, symbols)? else { return Ok(None) };
            matrices.push((k, m));
            owner.push((k, 0));
            owner.push((k, 1));
        }
        let axes = owner.len();
        let position = |k: usize, side: usize| owner.iter().position(|&o| o == (k, side)).unwrap();
        let mut partner: FxHashMap<usize, usize> = FxHashMap::default();
        for g in groups { partner.insert(g[0], g[1]); partner.insert(g[1], g[0]); }
        let free: Vec<usize> = (0..axes).filter(|a| !partner.contains_key(a)).collect();
        let matrix = |k: usize| &matrices.iter().find(|(j, _)| *j == k).unwrap().1;
        let mut seen: Vec<usize> = vec![];
        // The factors from (k, side) on, and the free axis that ends them.
        let walk = |k: usize, side: usize, seen: &mut Vec<usize>| -> R<(Vec<MatrixExpr>, Option<usize>)> {
            let (mut k, mut side) = (k, side);
            let start = (k, side);
            let mut chain = vec![];
            loop {
                chain.push(if side == 0 { matrix(k).clone() } else { matrix(k).transpose(true)? });
                seen.push(k);
                let out = position(k, 1 - side);
                let Some(&p) = partner.get(&out) else { return Ok((chain, Some(out))) };
                let Some(&o) = owner.get(p) else { return Err(AErr::Key(p.to_string())) };
                (k, side) = o;
                if (k, side) == start { return Ok((chain, None)); }
            }
        };
        let product = |chain: Vec<MatrixExpr>| -> R<MatrixExpr> {
            let mut it = chain.into_iter();
            let first = it.next().expect("a chain has a factor");
            it.try_fold(first, |a, b| Self::times(&a, &b))
        };
        let (mut chains, mut ends, mut pending) = (vec![], vec![], free.clone());
        while !pending.is_empty() {
            let start = pending.remove(0);
            let (k, side) = owner[start];
            let (chain, end) = walk(k, side, &mut seen)?;
            let Some(end) = end else { return Ok(None) };
            let Some(at) = pending.iter().position(|&a| a == end) else { return Ok(None) };
            pending.remove(at);
            ends.push(start);
            ends.push(end);
            chains.push(product(chain)?);
        }
        let ks: Vec<usize> = matrices.iter().map(|(k, _)| *k).collect();
        for k in ks {
            if seen.contains(&k) { continue; }
            let (chain, _) = walk(k, 0, &mut seen)?;
            // The trace of a cycle is evaluated (`Trace(...).doit(deep=False)`).
            let t = product(chain)?.trace_unevaluated()?;
            scalars.push(crate::functions::doit(&t, false));
        }
        let scalar = crate::ops::mul(&scalars);
        if chains.is_empty() { return Ok(Some(self.scal(scalar))); }
        if !scalar.is_one() { chains[0] = MatrixExpr::mul(scalar, &[chains[0].clone()], true)?.doit(true)?; }
        if chains.len() == 1 { return Ok(Some(self.mat(chains.pop().unwrap()))); }
        let form: Vec<i64> = ends.iter().map(|e| free.iter().position(|f| f == e).unwrap() as i64).collect();
        let form = af_invert(&form)?;
        let nodes: Vec<Node> = chains.into_iter().map(|m| self.mat(m)).collect();
        Ok(Some(self.permute_checked(Node::Tensor(nodes), form)?))
    }
}

/// The result of [`A2m::chain`].
pub enum Chain {
    Done(Node),
    /// The input itself (`return expr`).
    Keep,
    /// `FunctionMatrix(*shape, Lambda((i, j), expr[i, j]))` of the input.
    FunctionMatrix,
}
