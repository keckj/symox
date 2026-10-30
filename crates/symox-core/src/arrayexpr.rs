//! SymPy 1.14 array-expression canonicalization
//! (`sympy/tensor/array/expressions/array_expressions.py`): the
//! `_canonicalize` methods of ArrayTensorProduct, ArrayAdd, PermuteDims,
//! ArrayDiagonal and ArrayContraction, their flattening and denesting rules,
//! the contraction/diagonal index bookkeeping (`_push_indices_up/down`,
//! `_get_index_shifts`, `_lower_contraction_to_addends`,
//! `_sort_fully_contracted_args`) and `PermuteDims._nest_permutation`.
//!
//! The tree is index based: leaves (array symbols, matrices, scalars, ...)
//! are opaque entries of a leaf table with their rank and shape, and the
//! dimensions are opaque ids (equal ids are equal dimensions). The only
//! external decision is SymPy's `default_sort_key` order of fully
//! contracted tensor-product arguments, which the caller supplies.

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Node {
    Leaf(usize),
    Zero(Vec<usize>),
    Tensor(Vec<Node>),
    Add(Vec<Node>),
    Permute(Box<Node>, Vec<usize>),
    Contract(Box<Node>, Vec<Vec<usize>>),
    Diagonal(Box<Node>, Vec<Vec<usize>>),
    /// `Reshape(expr, shape)` (opaque to the canonical forms).
    Reshape(Box<Node>, Vec<usize>),
    /// `ArrayElementwiseApplyFunc(Lambda(d, d**e), expr)`, `e` the leaf index
    /// of the exponent (opaque to the canonical forms).
    PowApply(Box<Node>, usize),
    /// `ArrayElementwiseApplyFunc(function, expr)`, the function an opaque
    /// leaf-table entry (used by the array-to-matrix recognition).
    Apply(Box<Node>, usize),
}

#[derive(Clone, Debug)]
pub struct LeafInfo {
    pub rank: usize,
    pub shape: Option<Vec<usize>>,
    /// ZeroMatrix (treated as SymPy treats ZeroArray/ZeroMatrix).
    pub zero: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AErr {
    Value(String),
    NotImplemented(String),
    Type(String),
    Index(String),
    Attribute(String),
    /// An allocation the expression needs cannot be made (Python's MemoryError).
    Memory(String),
    /// The sort-key callback failed (the caller holds the error).
    Callback,
    Key(String),
    Assertion(String),
    /// A matrix-expression constructor failed.
    Matrix(crate::matrix::MatrixError),
    /// An input the native port does not model: the caller runs the stock
    /// SymPy algorithm instead (never a mathematical answer).
    Unsupported(String),
}

impl From<crate::matrix::MatrixError> for AErr {
    fn from(e: crate::matrix::MatrixError) -> Self { AErr::Matrix(e) }
}

pub type R<T> = Result<T, AErr>;

pub struct Ctx<'a> {
    pub leaves: &'a [LeafInfo],
    /// `ones[d]`: dimension id `d` is the integer 1.
    pub ones: &'a [bool],
    /// Stable order of the given nodes by SymPy's `default_sort_key`.
    pub order: &'a mut dyn FnMut(&[&Node]) -> R<Vec<usize>>,
}

pub fn invert(form: &[usize]) -> Vec<usize> {
    let mut out = vec![0; form.len()];
    for (i, &j) in form.iter().enumerate() { out[j] = i; }
    out
}

fn sort_groups(groups: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = groups.into_iter().map(|mut g| { g.sort_unstable(); g }).collect();
    groups.sort_by_key(|g| g.iter().copied().min());
    groups
}

/// `_build_push_indices_up_func_transformation`.
fn push_up_transform(flat: &[usize]) -> impl Fn(usize) -> Option<usize> + '_ {
    let mut shifts: Vec<(usize, usize)> = vec![(0, 0)];
    let (mut i, mut cumulative) = (0, 0);
    while i < flat.len() {
        let mut j = 1;
        while i + j < flat.len() && flat[i] + j == flat[i + j] { j += 1; }
        cumulative += j;
        match shifts.iter_mut().find(|(k, _)| *k == flat[i]) {
            Some(entry) => entry.1 = cumulative,
            None => shifts.push((flat[i], cumulative)),
        }
        i += j;
    }
    shifts.sort_unstable();
    move |idx| {
        if flat.contains(&idx) { return None; }
        let pos = shifts.partition_point(|(k, _)| *k <= idx);
        Some(idx - shifts[pos - 1].1)
    }
}

fn flat_sorted(groups: &[Vec<usize>]) -> Vec<usize> {
    let mut flat: Vec<usize> = groups.iter().flatten().copied().collect();
    flat.sort_unstable();
    flat
}

/// `ArrayContraction._push_indices_up` on a list of index lists.
pub fn contraction_push_up(groups: &[Vec<usize>], indices: &[Vec<usize>]) -> Vec<Vec<Option<usize>>> {
    let flat = flat_sorted(groups);
    let t = push_up_transform(&flat);
    indices.iter().map(|b| b.iter().map(|&j| t(j)).collect()).collect()
}

/// `ArrayContraction._push_indices_down` on a list of index lists.
pub fn contraction_push_down(groups: &[Vec<usize>], indices: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let flat = flat_sorted(groups);
    let Some(&last) = flat.last() else { return indices.to_vec() };
    let shifts: Vec<usize> = (0..last + 2).filter(|i| !flat.contains(i)).collect();
    let t = |j: usize| if j < shifts.len() { shifts[j] } else { j + shifts[shifts.len() - 1] + 1 - shifts.len() };
    indices.iter().map(|b| b.iter().map(|&j| t(j)).collect()).collect()
}

/// `ArrayDiagonal._get_positions_shape(range(rank), groups)` positions:
/// the free axes, then the diagonal groups.
#[derive(Clone, Debug, PartialEq)]
pub enum Position { Axis(usize), Group(Vec<usize>) }

pub fn diagonal_positions(rank: usize, groups: &[Vec<usize>]) -> Vec<Position> {
    let mut out: Vec<Position> = (0..rank).filter(|i| !groups.iter().any(|g| g.contains(i))).map(Position::Axis).collect();
    out.extend(groups.iter().cloned().map(Position::Group));
    out
}

/// `ArrayDiagonal._push_indices_down`.
pub fn diagonal_push_down(groups: &[Vec<usize>], indices: &[usize], rank: usize) -> Vec<Option<Position>> {
    let positions = diagonal_positions(rank, groups);
    indices.iter().map(|&x| positions.get(x).cloned()).collect()
}

/// `ArrayDiagonal._push_indices_up`.
pub fn diagonal_push_up(groups: &[Vec<usize>], indices: &[usize], rank: usize) -> Vec<Option<usize>> {
    let positions = diagonal_positions(rank, groups);
    indices.iter().map(|&x| positions.iter().position(|p| match p { Position::Axis(a) => *a == x, Position::Group(g) => g.contains(&x) })).collect()
}

/// SymPy `Permutation.cyclic_form` of an array form.
pub fn cyclic_form(form: &[usize]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; form.len()];
    let mut out = vec![];
    for start in 0..form.len() {
        if seen[start] { continue; }
        let mut cycle = vec![start];
        seen[start] = true;
        let mut j = form[start];
        while j != start { seen[j] = true; cycle.push(j); j = form[j]; }
        if cycle.len() > 1 { out.push(cycle); }
    }
    out
}

/// Array form of `Permutation(cycles, size=size)`.
fn from_cycles(cycles: &[Vec<usize>], size: usize) -> Vec<usize> {
    let size = cycles.iter().flatten().map(|&j| j + 1).max().unwrap_or(0).max(size);
    let mut form: Vec<usize> = (0..size).collect();
    for c in cycles {
        for (k, &a) in c.iter().enumerate() { form[a] = c[(k + 1) % c.len()]; }
    }
    form
}

fn cumulative(ranks: &[usize]) -> Vec<usize> {
    let mut out = vec![0];
    for r in ranks { out.push(out[out.len() - 1] + r); }
    out
}

impl Ctx<'_> {
    pub fn rank(&self, n: &Node) -> usize {
        match n {
            Node::Leaf(i) => self.leaves[*i].rank,
            Node::Zero(s) => s.len(),
            Node::Tensor(a) => a.iter().map(|x| self.rank(x)).sum(),
            Node::Add(a) => a.first().map_or(0, |x| self.rank(x)),
            Node::Permute(_, p) => p.len(),
            Node::Contract(e, g) => self.rank(e).saturating_sub(g.iter().map(Vec::len).sum::<usize>()),
            Node::Diagonal(e, g) => (self.rank(e) + g.len()).saturating_sub(g.iter().map(Vec::len).sum::<usize>()),
            Node::Reshape(_, s) => s.len(),
            Node::PowApply(e, _) | Node::Apply(e, _) => self.rank(e),
        }
    }

    pub fn shape(&self, n: &Node) -> Option<Vec<usize>> {
        match n {
            Node::Leaf(i) => self.leaves[*i].shape.clone(),
            Node::Zero(s) => Some(s.clone()),
            Node::Tensor(a) => {
                let mut out = vec![];
                for x in a { out.extend(self.shape(x)?); }
                Some(out)
            },
            Node::Add(a) => {
                let shapes: Vec<_> = a.iter().map(|x| self.shape(x)).collect();
                if shapes.iter().any(Option::is_none) { None } else { shapes.into_iter().next().flatten().or(Some(vec![])) }
            },
            Node::Permute(e, p) => { let s = self.shape(e)?; Some(p.iter().map(|&i| s[i]).collect()) },
            Node::Contract(e, g) => {
                let s = self.shape(e)?;
                Some(s.iter().enumerate().filter(|(i, _)| !g.iter().any(|x| x.contains(i))).map(|(_, d)| *d).collect())
            },
            Node::Diagonal(e, g) => {
                let s = self.shape(e)?;
                let mut out: Vec<usize> = s.iter().enumerate().filter(|(i, _)| !g.iter().any(|x| x.contains(i))).map(|(_, d)| *d).collect();
                out.extend(g.iter().map(|x| s[x[0]]));
                Some(out)
            },
            Node::Reshape(_, s) => Some(s.clone()),
            Node::PowApply(e, _) | Node::Apply(e, _) => self.shape(e),
        }
    }

    /// `_get_subranks`.
    pub fn subranks(&self, n: &Node) -> Vec<usize> {
        match n {
            Node::Tensor(a) => a.iter().map(|x| self.rank(x)).collect(),
            Node::Contract(e, _) | Node::Diagonal(e, _) | Node::Apply(e, _) => self.subranks(e),
            Node::Permute(e, _) => vec![self.rank(e)],
            _ => vec![self.rank(n)],
        }
    }

    pub fn subrank(&self, n: &Node) -> usize { self.subranks(n).iter().sum() }

    fn is_zero(&self, n: &Node) -> bool {
        match n { Node::Zero(_) => true, Node::Leaf(i) => self.leaves[*i].zero, _ => false }
    }

    // SymPy tests `isinstance(expr, (ZeroArray, ZeroMatrix))`: a zero-matrix leaf counts (`is_zero`).

    fn zero_shape(&self, n: &Node) -> R<Vec<usize>> {
        self.shape(n).ok_or_else(|| AErr::NotImplemented("cannot handle addition of ZeroMatrix/ZeroArray and undefined shape object".into()))
    }

    // ---------------------------------------------------------------- ATP

    /// `_array_tensor_product(*args)`.
    pub fn tensor(&mut self, args: Vec<Node>) -> R<Node> {
        let mut args: Vec<Node> = args.into_iter().flat_map(|a| match a { Node::Tensor(v) => v, o => vec![o] }).collect();
        let ranks: Vec<usize> = args.iter().map(|a| self.rank(a)).collect();
        let total: usize = ranks.iter().sum();
        let mut perm: Vec<usize> = (0..total).collect();
        let mut cycles = false;
        let mut offset = 0;
        for (i, a) in args.iter_mut().enumerate() {
            if let Node::Permute(e, p) = a {
                for (k, &pk) in p.iter().enumerate() { perm[offset + k] = offset + pk; }
                cycles |= p.iter().enumerate().any(|(k, &pk)| k != pk);
                *a = (**e).clone();
            }
            offset += ranks[i];
        }
        if cycles {
            let t = self.tensor(args)?;
            return self.permute(t, perm);
        }
        if args.len() == 1 { return Ok(args.pop().unwrap()); }
        if args.iter().any(|a| self.is_zero(a)) {
            let mut shape = vec![];
            for a in &args { shape.extend(self.shape(a).ok_or_else(|| AErr::Type("unsupported operand for shape concatenation".into()))?); }
            return Ok(Node::Zero(shape));
        }
        if args.iter().any(|a| matches!(a, Node::Contract(..))) {
            let ranks2: Vec<usize> = args.iter().map(|a| if matches!(a, Node::Contract(..)) { self.subrank(a) } else { self.rank(a) }).collect();
            let cum = cumulative(&ranks2);
            let mut groups = vec![];
            let mut inner = vec![];
            for (i, a) in args.into_iter().enumerate() {
                match a {
                    Node::Contract(e, g) => {
                        groups.extend(g.iter().map(|x| x.iter().map(|k| cum[i] + k).collect::<Vec<_>>()));
                        inner.push(*e);
                    },
                    o => inner.push(o),
                }
            }
            let t = self.tensor(inner)?;
            return self.contract(t, groups);
        }
        if args.iter().any(|a| matches!(a, Node::Diagonal(..))) {
            let cum = cumulative(&ranks);
            let mut inverse = vec![];
            let mut last = vec![];
            for (i, a) in args.iter().enumerate() {
                if let Node::Diagonal(_, g) = a {
                    let i1 = self.rank(a) - g.len();
                    inverse.extend((0..i1).map(|j| cum[i] + j));
                    last.extend((i1..i1 + g.len()).map(|j| cum[i] + j));
                } else {
                    inverse.extend((0..ranks[i]).map(|j| cum[i] + j));
                }
            }
            inverse.extend(last);
            let ranks2: Vec<usize> = args.iter().map(|a| if matches!(a, Node::Diagonal(..)) { self.subrank(a) } else { self.rank(a) }).collect();
            let cum2 = cumulative(&ranks2);
            let mut groups = vec![];
            let mut inner = vec![];
            for (i, a) in args.into_iter().enumerate() {
                match a {
                    Node::Diagonal(e, g) => {
                        groups.extend(g.iter().map(|x| x.iter().map(|k| cum2[i] + k).collect::<Vec<_>>()));
                        inner.push(*e);
                    },
                    o => inner.push(o),
                }
            }
            let t = self.tensor(inner)?;
            let d = self.diagonal(t, groups, false)?;
            return self.permute(d, invert(&inverse));
        }
        Ok(Node::Tensor(args))
    }

    // ---------------------------------------------------------------- Add

    /// `_array_add(*args)` (the rank and shape checks are the caller's).
    pub fn add(&mut self, args: Vec<Node>) -> R<Node> {
        let args: Vec<Node> = args.into_iter().flat_map(|a| match a { Node::Add(v) => v, o => vec![o] }).collect();
        let first = args.first().cloned();
        let mut kept: Vec<Node> = args.into_iter().filter(|a| !self.is_zero(a)).collect();
        match kept.len() {
            0 => {
                let first = first.ok_or_else(|| AErr::Value("summing arrays of different ranks".into()))?;
                Ok(Node::Zero(self.zero_shape(&first)?))
            },
            1 => Ok(kept.pop().unwrap()),
            _ => Ok(Node::Add(kept)),
        }
    }

    // ---------------------------------------------------------------- PermuteDims

    /// `_permute_dims(expr, permutation)`.
    pub fn permute(&mut self, expr: Node, perm: Vec<usize>) -> R<Node> {
        if perm.len() != self.rank(&expr) {
            return Err(AErr::Value("Permutation size must be the length of the shape of expr".into()));
        }
        let (mut expr, mut perm) = (expr, perm);
        if let Node::Permute(sub, subp) = expr {
            perm = perm.iter().map(|&i| subp[i]).collect();
            expr = *sub;
        }
        if matches!(expr, Node::Contract(..)) {
            (expr, perm) = self.permute_denest_contraction(expr, perm)?;
        }
        if let Node::Tensor(args) = expr {
            (expr, perm) = self.permute_denest_tensor(args, perm)?;
        }
        if self.is_zero(&expr) {
            let s = self.zero_shape(&expr)?;
            return Ok(Node::Zero(perm.iter().map(|&i| s[i]).collect()));
        }
        if perm.iter().enumerate().all(|(i, &p)| i == p) { return Ok(expr); }
        Ok(Node::Permute(Box::new(expr), perm))
    }

    fn permute_denest_tensor(&mut self, args: Vec<Node>, perm: Vec<usize>) -> R<(Node, Vec<usize>)> {
        let image = invert(&perm);
        let ranks: Vec<usize> = args.iter().map(|a| self.rank(a)).collect();
        let cum = cumulative(&ranks);
        let comps: Vec<Vec<usize>> = (0..args.len()).map(|i| image[cum[i]..cum[i + 1]].to_vec()).collect();
        let mut ps: Vec<(usize, Vec<usize>)> = comps.iter().enumerate().map(|(i, c)| { let mut s = c.clone(); s.sort_unstable(); (i, s) }).collect();
        ps.sort_by(|a, b| a.1.cmp(&b.1));
        let flat: Vec<usize> = ps.iter().flat_map(|(i, _)| comps[*i].clone()).collect();
        let sorted: Vec<Node> = ps.iter().map(|(i, _)| args[*i].clone()).collect();
        Ok((self.tensor(sorted)?, invert(&flat)))
    }

    fn permute_denest_contraction(&mut self, expr: Node, perm: Vec<usize>) -> R<(Node, Vec<usize>)> {
        let Node::Contract(inner, groups) = &expr else { return Ok((expr, perm)) };
        let Node::Tensor(args) = &**inner else { return Ok((expr, perm)) };
        let subranks: Vec<usize> = args.iter().map(|a| self.rank(a)).collect();
        let flat: Vec<usize> = groups.iter().flatten().copied().collect();
        let cum = cumulative(&subranks);
        let image = invert(&perm);
        let mut blocks_up = vec![];
        let mut counter = 0;
        for i in 0..subranks.len() {
            let mut current = vec![];
            for j in cum[i]..cum[i + 1] {
                if flat.contains(&j) { continue; }
                current.push(image[counter]);
                counter += 1;
            }
            blocks_up.push(current);
        }
        let index_blocks: Vec<Vec<usize>> = (0..subranks.len()).map(|i| (cum[i]..cum[i + 1]).collect()).collect();
        let up = contraction_push_up(groups, &index_blocks);
        let keys: Vec<Vec<usize>> = up.iter().map(|b| b.iter().flatten().map(|&j| image[j]).collect()).collect();
        let mut order: Vec<usize> = (0..args.len()).collect();
        order.sort_by(|&a, &b| keys[a].cmp(&keys[b]));
        let new_form = invert(&order.iter().flat_map(|&i| index_blocks[i].clone()).collect::<Vec<_>>());
        let new_groups: Vec<Vec<usize>> = groups.iter().map(|g| g.iter().map(|&j| new_form[j]).collect()).collect();
        let new_args: Vec<Node> = order.iter().map(|&i| args[i].clone()).collect();
        let t = self.tensor(new_args)?;
        let new_expr = self.contract(t, new_groups)?;
        let new_perm = invert(&order.iter().flat_map(|&k| blocks_up[k].clone()).collect::<Vec<_>>());
        Ok((new_expr, new_perm))
    }

    // ---------------------------------------------------------------- ArrayDiagonal

    fn validate_diagonal(&self, expr: &Node, groups: &[Vec<usize>], allow_trivial: bool) -> R<()> {
        let Some(shape) = self.shape(expr) else { return Ok(()) };
        for g in groups {
            if g.iter().any(|&j| j >= shape.len()) { return Err(AErr::Value("index is larger than expression shape".into())); }
            let mut dims: Vec<usize> = g.iter().map(|&j| shape[j]).collect();
            dims.dedup();
            dims.sort_unstable();
            dims.dedup();
            if dims.len() != 1 { return Err(AErr::Value("diagonalizing indices of different dimensions".into())); }
            if !allow_trivial && g.len() <= 1 { return Err(AErr::Value("need at least two axes to diagonalize".into())); }
            let mut s = g.clone();
            s.sort_unstable();
            s.dedup();
            if s.len() != g.len() { return Err(AErr::Value("axis index cannot be repeated".into())); }
        }
        Ok(())
    }

    /// `_array_diagonal(expr, *diagonal_indices)`.
    pub fn diagonal(&mut self, expr: Node, groups: Vec<Vec<usize>>, allow_trivial: bool) -> R<Node> {
        let groups: Vec<Vec<usize>> = groups.into_iter().map(|mut g| { g.sort_unstable(); g }).collect();
        self.validate_diagonal(&expr, &groups, allow_trivial)?;
        if groups.is_empty() { return Ok(expr); }
        if groups.iter().any(|g| g.len() == 1) {
            let short: Vec<Vec<usize>> = groups.iter().filter(|g| g.len() > 1).cloned().collect();
            let rank1 = self.rank(&Node::Diagonal(Box::new(expr.clone()), groups.clone()));
            let rank3 = rank1 - groups.len();
            let down = diagonal_push_down(&short, &(0..rank1).collect::<Vec<_>>(), self.rank(&expr));
            let mut inverse = vec![];
            let mut counter = 0;
            for p in down {
                match p {
                    Some(Position::Axis(a)) => match groups.iter().position(|g| g.len() == 1 && g[0] == a) {
                        Some(t) => inverse.push(rank3 + t),
                        None => { inverse.push(counter); counter += 1; },
                    },
                    Some(Position::Group(g)) => match groups.iter().position(|x| *x == g) {
                        Some(t) => inverse.push(rank3 + t),
                        None => return Err(AErr::Value("invalid diagonal indices".into())),
                    },
                    None => return Err(AErr::Value("invalid diagonal indices".into())),
                }
            }
            let perm = invert(&inverse);
            let inner = if short.is_empty() { expr } else { self.diagonal(expr, short, false)? };
            return self.permute(inner, perm);
        }
        match expr {
            Node::Add(args) => {
                let mut out = vec![];
                for a in args { out.push(self.diagonal(a, groups.clone(), false)?); }
                self.add(out)
            },
            Node::Diagonal(inner, inner_groups) => {
                let all_inner = flat_sorted(&inner_groups);
                let total = self.subrank(&Node::Diagonal(inner.clone(), inner_groups.clone()));
                let shifts = index_shifts(total, &all_inner);
                let mut all = inner_groups.clone();
                for g in &groups { all.push(g.iter().map(|&j| shifted(&shifts, j)).collect::<R<Vec<_>>>()?); }
                self.diagonal(*inner, all, false)
            },
            Node::Permute(inner, p) => {
                let rank = p.len();
                let back: Vec<Vec<usize>> = groups.iter().map(|g| g.iter().map(|&j| p[j]).collect()).collect();
                let nondiag: Vec<usize> = (0..rank).filter(|i| !groups.iter().any(|g| g.contains(i))).collect();
                let back_nondiag: Vec<usize> = nondiag.iter().map(|&i| p[i]).collect();
                let mut sorted = back_nondiag.clone();
                sorted.sort_unstable();
                let mut perm: Vec<usize> = back_nondiag.iter().map(|e| sorted.iter().position(|x| x == e).unwrap()).collect();
                let shift = perm.len();
                perm.extend((0..back.len()).map(|i| i + shift));
                let d = self.diagonal(*inner, back, false)?;
                self.permute(d, perm)
            },
            e if self.is_zero(&e) => {
                let d = Node::Diagonal(Box::new(e), groups);
                Ok(Node::Zero(self.zero_shape(&d)?))
            },
            e => Ok(Node::Diagonal(Box::new(e), groups)),
        }
    }

    // ---------------------------------------------------------------- ArrayContraction

    fn validate_contraction(&self, expr: &Node, groups: &[Vec<usize>]) -> R<()> {
        let Some(shape) = self.shape(expr) else { return Ok(()) };
        for g in groups {
            if g.iter().any(|&j| j >= shape.len()) { return Err(AErr::Index("tuple index out of range".into())); }
            let mut dims: Vec<usize> = g.iter().map(|&j| shape[j]).collect();
            dims.sort_unstable();
            dims.dedup();
            if dims.len() != 1 { return Err(AErr::Value("contracting indices of different dimensions".into())); }
        }
        // SymPy builds an inconsistent object for a repeated axis; it is refused here
        let mut flat: Vec<usize> = groups.iter().flatten().copied().collect();
        let n = flat.len();
        flat.sort_unstable();
        flat.dedup();
        if flat.len() != n { return Err(AErr::Value("contraction indices must be distinct".into())); }
        Ok(())
    }

    /// `_array_contraction(expr, *contraction_indices)`.
    pub fn contract(&mut self, expr: Node, groups: Vec<Vec<usize>>) -> R<Node> {
        let groups = sort_groups(groups);
        self.validate_contraction(&expr, &groups)?;
        if groups.is_empty() { return Ok(expr); }
        let mut expr = expr;
        let mut groups = groups;
        match expr {
            Node::Contract(inner, inner_groups) => {
                let all_inner = flat_sorted(&inner_groups);
                let total = self.subrank(&Node::Contract(inner.clone(), inner_groups.clone()));
                let shifts = index_shifts(total, &all_inner);
                let mut all = inner_groups.clone();
                for g in &groups { all.push(g.iter().map(|&j| shifted(&shifts, j)).collect::<R<Vec<_>>>()?); }
                return self.contract(*inner, all);
            },
            ref e if self.is_zero(e) => {
                let s = self.zero_shape(e)?;
                let flat: Vec<usize> = groups.iter().flatten().copied().collect();
                return Ok(Node::Zero(s.iter().enumerate().filter(|(i, _)| !flat.contains(i)).map(|(_, d)| *d).collect()));
            },
            Node::Permute(inner, p) => {
                let new_groups: Vec<Vec<usize>> = groups.iter().map(|g| g.iter().map(|&j| p[j]).collect()).collect();
                let plist: Vec<usize> = p.iter().copied().filter(|i| !new_groups.iter().any(|g| g.contains(i))).collect();
                let up: Vec<usize> = contraction_push_up(&new_groups, &[plist])[0].iter().map(|x| x.unwrap()).collect();
                let c = self.contract(*inner, new_groups)?;
                return self.permute(c, up);
            },
            Node::Tensor(args) => {
                let (e, g) = self.sort_fully_contracted(args, groups)?;
                let (e, g) = self.lower_to_addends(e, g)?;
                if g.is_empty() { return Ok(e); }
                expr = e;
                groups = g;
            },
            _ => {},
        }
        match expr {
            Node::Diagonal(inner, dgroups) => {
                let rank = self.rank(&inner);
                let positions = diagonal_positions(rank, &dgroups);
                let mut diags: Vec<Option<Vec<usize>>> = dgroups.into_iter().map(Some).collect();
                let mut new_groups = vec![];
                for g in &groups {
                    let mut down = vec![];
                    for &x in g {
                        match positions.get(x) {
                            Some(Position::Axis(a)) => down.push(*a),
                            Some(Position::Group(v)) => down.extend(v.iter().copied()),
                            None => return Err(AErr::Value("invalid contraction indices".into())),
                        }
                    }
                    let mut ind = down.clone();
                    for d in diags.iter_mut() {
                        if let Some(v) = d {
                            if down.iter().any(|i| v.contains(i)) { ind.extend(v.iter().copied()); *d = None; }
                        }
                    }
                    ind.sort_unstable();
                    ind.dedup();
                    new_groups.push(ind);
                }
                let remaining: Vec<Vec<usize>> = diags.into_iter().flatten().collect();
                let new_diags: Vec<Vec<usize>> = contraction_push_up(&new_groups, &remaining).into_iter()
                    .map(|v| v.into_iter().map(|x| x.ok_or_else(|| AErr::Value("invalid diagonal indices".into()))).collect::<R<Vec<_>>>())
                    .collect::<R<_>>()?;
                let c = self.contract(*inner, new_groups)?;
                self.diagonal(c, new_diags, false)
            },
            Node::Add(args) => {
                let mut out = vec![];
                for a in args { out.push(self.contract(a, groups.clone())?); }
                self.add(out)
            },
            e => {
                let shape = self.shape(&e);
                let groups: Vec<Vec<usize>> = groups.into_iter().filter(|g| g.len() > 1 || shape.as_ref().is_none_or(|s| !self.ones[s[g[0]]])).collect();
                if groups.is_empty() { Ok(e) } else { Ok(Node::Contract(Box::new(e), groups)) }
            },
        }
    }

    fn sort_fully_contracted(&mut self, args: Vec<Node>, groups: Vec<Vec<usize>>) -> R<(Node, Vec<Vec<usize>>)> {
        let expr = Node::Tensor(args);
        if self.shape(&expr).is_none() { return Ok((expr, groups)); }
        let Node::Tensor(args) = expr else { unreachable!() };
        let ranks: Vec<usize> = args.iter().map(|a| self.rank(a)).collect();
        let cum = cumulative(&ranks);
        let flat: Vec<usize> = groups.iter().flatten().copied().collect();
        let full: Vec<bool> = (0..args.len()).map(|i| (cum[i]..cum[i + 1]).all(|j| flat.contains(&j))).collect();
        let full_idx: Vec<usize> = (0..args.len()).filter(|&i| full[i]).collect();
        let mut new_pos: Vec<usize> = if full_idx.len() > 1 {
            let refs: Vec<&Node> = full_idx.iter().map(|&i| &args[i]).collect();
            let order = (self.order)(&refs)?;
            order.into_iter().map(|k| full_idx[k]).collect()
        } else { full_idx.clone() };
        new_pos.extend((0..args.len()).filter(|&i| !full[i]));
        let blocks: Vec<usize> = new_pos.iter().flat_map(|&i| cum[i]..cum[i + 1]).collect();
        let form = invert(&blocks);
        let new_groups = sort_groups(groups.iter().map(|g| g.iter().map(|&j| form[j]).collect()).collect());
        let new_args: Vec<Node> = new_pos.iter().map(|&i| args[i].clone()).collect();
        Ok((self.tensor(new_args)?, new_groups))
    }

    fn lower_to_addends(&mut self, expr: Node, groups: Vec<Vec<usize>>) -> R<(Node, Vec<Vec<usize>>)> {
        if matches!(expr, Node::Add(_)) { return Err(AErr::NotImplemented(String::new())); }
        let Node::Tensor(args) = expr else { return Ok((expr, groups)) };
        let ranks: Vec<usize> = args.iter().map(|a| self.rank(a)).collect();
        let cum = cumulative(&ranks);
        let mut remaining = vec![];
        let mut per_arg: Vec<Vec<Vec<usize>>> = vec![vec![]; args.len()];
        let mut backshift = vec![];
        for g in &groups {
            let target = (0..args.len()).find(|&j| matches!(args[j], Node::Add(_)) && g.iter().all(|&k| cum[j] <= k && k < cum[j + 1]));
            match target {
                Some(j) => { per_arg[j].push(g.iter().map(|k| k - cum[j]).collect()); backshift.extend(g.iter().copied()); },
                None => remaining.push(g.clone()),
            }
        }
        if remaining.len() == groups.len() { return Ok((Node::Tensor(args), groups)); }
        let total: usize = ranks.iter().sum();
        let mut shifts = vec![0; total];
        let mut acc = 0;
        for (i, s) in shifts.iter_mut().enumerate() { if backshift.contains(&i) { acc += 1; } *s = acc; }
        let remaining: Vec<Vec<usize>> = remaining.into_iter().map(|g| g.into_iter().map(|j| j - shifts[j]).collect()).collect();
        let mut new_args = vec![];
        for (a, c) in args.into_iter().zip(per_arg) { new_args.push(self.contract(a, c)?); }
        Ok((self.tensor(new_args)?, remaining))
    }

    // ---------------------------------------------------------------- canonical forms

    /// `_canonicalize` of the top node (its arguments are kept as given).
    pub fn canonicalize(&mut self, n: Node, allow_trivial: bool) -> R<Node> {
        match n {
            Node::Tensor(a) => self.tensor(a),
            Node::Add(a) => self.add(a),
            Node::Permute(e, p) => self.permute(*e, p),
            Node::Contract(e, g) => self.contract(*e, g),
            Node::Diagonal(e, g) => self.diagonal(*e, g, allow_trivial),
            o => Ok(o),
        }
    }

    /// `doit()`: the arguments first, then the node.
    pub fn doit(&mut self, n: Node, deep: bool) -> R<Node> {
        if !deep { return self.canonicalize(n, false); }
        let n = match n {
            Node::Tensor(a) => Node::Tensor(a.into_iter().map(|x| self.doit(x, true)).collect::<R<_>>()?),
            Node::Add(a) => Node::Add(a.into_iter().map(|x| self.doit(x, true)).collect::<R<_>>()?),
            Node::Permute(e, p) => Node::Permute(Box::new(self.doit(*e, true)?), p),
            Node::Contract(e, g) => Node::Contract(Box::new(self.doit(*e, true)?), g),
            Node::Diagonal(e, g) => Node::Diagonal(Box::new(self.doit(*e, true)?), g),
            Node::Reshape(e, s) => Node::Reshape(Box::new(self.doit(*e, true)?), s),
            Node::PowApply(e, x) => Node::PowApply(Box::new(self.doit(*e, true)?), x),
            Node::Apply(e, x) => Node::Apply(Box::new(self.doit(*e, true)?), x),
            o => o,
        };
        self.canonicalize(n, false)
    }

    /// `PermuteDims.nest_permutation()` (None: unchanged).
    pub fn nest_permutation(&mut self, n: Node) -> R<Node> {
        let Node::Permute(expr, perm) = n else { return Ok(n) };
        match *expr {
            Node::Tensor(args) => {
                let ranks: Vec<usize> = args.iter().map(|a| self.rank(a)).collect();
                let cum = cumulative(&ranks);
                let mut args = args;
                let mut keep = vec![];
                for cycle in cyclic_form(&perm) {
                    let (lo, hi) = (*cycle.iter().min().unwrap(), *cycle.iter().max().unwrap());
                    match (0..ranks.len()).find(|&j| lo >= cum[j] && hi < cum[j + 1]) {
                        Some(j) => {
                            let local: Vec<Vec<usize>> = vec![cycle.iter().map(|k| k - cum[j]).collect()];
                            let a = std::mem::replace(&mut args[j], Node::Tensor(vec![]));
                            // SymPy's `Permutation([[...]])` has the size of its largest point
                            args[j] = self.permute(a, from_cycles(&local, 0))?;
                        },
                        None => keep.push(cycle),
                    }
                }
                let t = self.tensor(args)?;
                self.permute(t, from_cycles(&keep, perm.len()))
            },
            Node::Contract(inner, groups) => {
                let total = self.subrank(&Node::Contract(inner.clone(), groups.clone()));
                let shifts = index_shifts(total, &flat_sorted(&groups));
                let cycles: Vec<Vec<usize>> = cyclic_form(&perm).into_iter().map(|c| c.into_iter().map(|j| shifted(&shifts, j)).collect::<R<Vec<_>>>()).collect::<R<_>>()?;
                let newperm = from_cycles(&cycles, 0);
                if newperm.len() != self.rank(&inner) {
                    return Err(AErr::Value("Permutation size must be the length of the shape of expr".into()));
                }
                let new_groups: Vec<Vec<usize>> = groups.iter().map(|g| g.iter().map(|&j| newperm[j]).collect()).collect();
                self.contract(Node::Permute(inner, newperm), new_groups)
            },
            Node::Add(args) => {
                let args = args.into_iter().map(|a| Node::Permute(Box::new(a), perm.clone())).collect();
                self.add(args)
            },
            e => Ok(Node::Permute(Box::new(e), perm)),
        }
    }
}

fn shifted(shifts: &[usize], j: usize) -> R<usize> {
    shifts.get(j).map(|s| s + j).ok_or_else(|| AErr::Index("list index out of range".into()))
}

/// `ArrayContraction._get_index_shifts` / the shifts of `ArrayDiagonal._flatten`.
pub fn index_shifts(total_rank: usize, all_inner: &[usize]) -> Vec<usize> {
    let inner_rank = all_inner.len();
    let outer_rank = total_rank.saturating_sub(inner_rank);
    let mut shifts = vec![0; outer_rank];
    let (mut counter, mut pointer) = (0, 0);
    for s in shifts.iter_mut() {
        while pointer < inner_rank && counter >= all_inner[pointer] { counter += 1; pointer += 1; }
        *s += pointer;
        counter += 1;
    }
    shifts
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx_run<T>(leaves: &[LeafInfo], f: impl FnOnce(&mut Ctx) -> T) -> T {
        let ones = vec![false, true];
        let mut order = |n: &[&Node]| Ok((0..n.len()).collect());
        let mut c = Ctx { leaves, ones: &ones, order: &mut order };
        f(&mut c)
    }

    fn mat() -> LeafInfo { LeafInfo { rank: 2, shape: Some(vec![0, 0]), zero: false } }

    #[test]
    fn tensor_lifts_permutations_and_contractions() {
        let leaves = vec![mat(), mat(), mat()];
        ctx_run(&leaves, |c| {
            // _array_tensor_product(PermuteDims(M, [1, 0]), N) -> PermuteDims(M x N, [1, 0, 2, 3])
            let r = c.tensor(vec![Node::Permute(Box::new(Node::Leaf(0)), vec![1, 0]), Node::Leaf(1)]).unwrap();
            assert_eq!(r, Node::Permute(Box::new(Node::Tensor(vec![Node::Leaf(0), Node::Leaf(1)])), vec![1, 0, 2, 3]));
            // contraction lifted out of a tensor product
            let r = c.tensor(vec![Node::Leaf(2), Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(0), Node::Leaf(1)])), vec![vec![1, 2]])]).unwrap();
            assert_eq!(r, Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(2), Node::Leaf(0), Node::Leaf(1)])), vec![vec![3, 4]]));
        });
    }

    #[test]
    fn contraction_flatten_and_trivial_axes() {
        let leaves = vec![mat(), mat(), LeafInfo { rank: 2, shape: Some(vec![0, 1]), zero: false }];
        ctx_run(&leaves, |c| {
            let inner = Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(0), Node::Leaf(1)])), vec![vec![1, 2]]);
            let r = c.contract(inner, vec![vec![0, 1]]).unwrap();
            assert_eq!(r, Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(0), Node::Leaf(1)])), vec![vec![0, 3], vec![1, 2]]));
            assert_eq!(c.contract(Node::Leaf(2), vec![vec![1]]).unwrap(), Node::Leaf(2));
            assert_eq!(c.contract(Node::Zero(vec![0, 0, 1]), vec![vec![0, 1]]).unwrap(), Node::Zero(vec![1]));
        });
    }

    #[test]
    fn permutation_bookkeeping() {
        assert_eq!(index_shifts(4, &[1, 2]), vec![0, 2]);
        assert_eq!(contraction_push_up(&[vec![1, 2]], &[vec![0, 1, 2, 3]]), vec![vec![Some(0), None, None, Some(1)]]);
        assert_eq!(contraction_push_down(&[vec![1, 2]], &[vec![0, 1]]), vec![vec![0, 3]]);
        assert_eq!(cyclic_form(&[1, 0, 3, 2]), vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(from_cycles(&[vec![2, 3]], 4), vec![0, 1, 3, 2]);
    }
}
