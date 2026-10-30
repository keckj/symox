//! SymPy 1.14 `convert_indexed_to_array` / `_convert_indexed_to_array`
//! (`sympy/tensor/array/expressions/from_indexed_to_array.py`) and the index
//! plan of `convert_array_to_indexed` (`from_array_to_indexed.py`).
//!
//! Indexed to array: the scalar expression is given as an `SNode` tree
//! (sums, products, sums of terms, powers, functions, Kronecker deltas,
//! matrix/array/indexed elements, other scalars). Indices are atom ids (the
//! caller interns them with Python equality, `int_atoms[a]` says whether
//! atom `a` is an `int`/`Integer`); a group of indices identified by
//! Kronecker deltas is an `It::Set`. The array results are canonical forms
//! of `arrayexpr.rs` over the leaf table; the `OneArray` objects of the
//! broadcasting rule are appended to that table as new leaves.
//!
//! An element-wise function `ArrayElementwiseApplyFunc(type(expr), sub)` is
//! the node `Node::PowApply(sub, f)` with `f` the leaf of the function (the
//! canonical forms only use its rank and shape; the caller decodes it).
//!
//! Array to indexed: `index_plan` computes the scalar expression as a
//! `Plan` (products, sums, `Sum`s with numbered dummies, leaf elements with
//! index expressions, applied functions); the caller builds the objects.
use crate::arrayexpr::{diagonal_push_down, invert, AErr, Ctx, LeafInfo, Node, Position};
use std::collections::BTreeSet;

/// An index: an atom or a set of atoms made equal by Kronecker deltas.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum It {
    A(usize),
    S(BTreeSet<usize>),
}

impl It {
    fn contains(&self, a: usize) -> bool { matches!(self, It::S(s) if s.contains(&a)) }
}

/// A limit `(ind, istart, iend)` of a Sum: the atom, whether `istart != 0`
/// and the dimension id of `iend + 1`; `None` when the limit is not a triple
/// (its length).
pub type Limit = Result<(usize, bool, usize), usize>;

#[derive(Clone, Debug)]
pub enum SNode {
    /// `Sum(function, *limits)`: the function, `expr.variables` and limits.
    Sum(Box<SNode>, Vec<usize>, Vec<Limit>),
    Mul(Vec<SNode>),
    Add(Vec<SNode>),
    /// `MatrixElement`, `ArrayElement`, `Indexed`: the leaf of the parent and
    /// the index atoms.
    Element(usize, Vec<usize>),
    /// `KroneckerDelta(i, j)`: its leaf and index atoms.
    Delta(usize, usize, usize),
    /// `Pow(base, exp)`: the base, and the integer exponent or the leaf of
    /// the whole power (non-integer exponent).
    Pow(Box<SNode>, Result<i64, usize>),
    /// A `Function`: the node of `args[0]` (`None`: no argument) and the leaf
    /// of `type(expr)`.
    Func(Option<Box<SNode>>, usize),
    /// An `IndexedBase` (the reference raises NotImplementedError).
    Base,
    /// Any other expression (returned with no indices).
    Leaf(usize),
}

/// Errors whose message the caller formats with the objects.
#[derive(Clone, Debug, PartialEq)]
pub enum IErr {
    A(AErr),
    /// `_get_argindex`: IndexError("%s not found in %s" % (ind, subindices)).
    NotFound(usize, Vec<It>),
    /// ValueError("summation index and array dimension mismatch: %s" % ind).
    Mismatch(usize),
    /// AttributeError: the leaf has no `shape` attribute.
    NoShape(usize),
}

impl From<AErr> for IErr {
    fn from(e: AErr) -> Self { IErr::A(e) }
}

pub type IR<T> = Result<T, IErr>;

fn type_err(m: &str) -> IErr { IErr::A(AErr::Type(m.into())) }
fn value_err(m: &str) -> IErr { IErr::A(AErr::Value(m.into())) }

pub struct IndexConv<'a> {
    /// The given leaves followed by the `OneArray` leaves made here.
    pub leaves: Vec<LeafInfo>,
    /// Leaves of undefined rank (IndexedBase without shape, `get_rank` -1).
    pub unknown: Vec<bool>,
    /// Leaves that have a `shape` attribute.
    pub has_shape: Vec<bool>,
    /// Number of given leaves (the following ones are `OneArray`s).
    pub given: usize,
    pub ones: &'a [bool],
    pub int_atoms: &'a [bool],
    /// `default_sort_key` order; receives the full leaf table.
    pub order: &'a mut dyn FnMut(&[&Node], &[LeafInfo]) -> Result<Vec<usize>, AErr>,
}

/// SymPy `Permutation(list)` validation (array form).
pub fn permutation(form: &[usize]) -> Result<Vec<usize>, AErr> {
    let mut seen = BTreeSet::new();
    for &i in form {
        if !seen.insert(i) { return Err(AErr::Value("there were repeated elements.".into())); }
    }
    if (0..seen.len()).any(|i| !seen.contains(&i)) {
        let max = seen.iter().next_back().copied().unwrap_or(0);
        return Err(AErr::Value(format!("Integers 0 through {max} must be present.")));
    }
    Ok(form.to_vec())
}

/// `_get_diagonal_indices(flattened_indices)`.
fn diagonal_indices(flat: &[It], int_atoms: &[bool]) -> (Vec<Vec<usize>>, Vec<It>) {
    let mut keys: Vec<It> = vec![];
    let mut axes: Vec<Vec<usize>> = vec![];
    for (i, ind) in flat.iter().enumerate() {
        if let It::A(a) = ind { if int_atoms[*a] { continue; } }
        match keys.iter().position(|k| k == ind) {
            Some(p) => axes[p].push(i),
            None => { keys.push(ind.clone()); axes.push(vec![i]); },
        }
    }
    let diag: Vec<(It, Vec<usize>)> = keys.into_iter().zip(axes).filter(|(_, v)| v.len() > 1).collect();
    let mut ret: Vec<It> = flat.iter().filter(|i| !diag.iter().any(|(k, _)| k == *i)).cloned().collect();
    // Sorted by first occurrence, which is the insertion order already.
    ret.extend(diag.iter().map(|(k, _)| k.clone()));
    (diag.into_iter().map(|(_, v)| v).collect(), ret)
}

/// `_get_argindex(subindices, ind)`.
fn argindex(sub: &[It], ind: usize) -> Option<usize> {
    sub.iter().position(|s| *s == It::A(ind) || s.contains(ind))
}

impl IndexConv<'_> {
    fn with_ctx<T>(&mut self, f: impl FnOnce(&mut Ctx) -> Result<T, AErr>) -> IR<T> {
        let leaves = &self.leaves;
        let order = &mut self.order;
        let mut o = |n: &[&Node]| order(n, leaves);
        let mut ctx = Ctx { leaves, ones: self.ones, order: &mut o };
        f(&mut ctx).map_err(IErr::A)
    }

    fn shape(&mut self, n: &Node) -> IR<Option<Vec<usize>>> { self.with_ctx(|c| Ok(c.shape(n))) }

    fn has_unknown(&self, n: &Node) -> bool {
        match n {
            Node::Leaf(i) => self.unknown[*i],
            Node::Zero(_) => false,
            Node::Tensor(a) | Node::Add(a) => a.iter().any(|x| self.has_unknown(x)),
            Node::Permute(e, _) | Node::Contract(e, _) | Node::Diagonal(e, _) | Node::Reshape(e, _) | Node::PowApply(e, _) => self.has_unknown(e),
        }
    }

    /// `get_rank`.
    fn get_rank(&mut self, n: &Node) -> IR<i64> {
        if let Node::Leaf(i) = n {
            return Ok(if self.unknown[*i] { -1 } else { self.leaves[*i].rank as i64 });
        }
        match self.shape(n)? {
            Some(s) => Ok(s.len() as i64),
            None => Err(type_err("object of type 'NoneType' has no len()")),
        }
    }

    fn tensor(&mut self, args: Vec<Node>) -> IR<Node> { self.with_ctx(|c| c.tensor(args)) }

    fn diagonal(&mut self, n: Node, groups: Vec<Vec<usize>>) -> IR<Node> {
        self.with_ctx(|c| c.diagonal(n, groups, false))
    }

    fn contract(&mut self, n: Node, groups: Vec<Vec<usize>>) -> IR<Node> {
        if groups.iter().any(|g| g.len() == 1) && self.shape(&n)?.is_none() {
            return Err(type_err("'NoneType' object is not subscriptable"));
        }
        self.with_ctx(|c| c.contract(n, groups))
    }

    /// `_permute_dims(expr, Permutation(form))` (the form is valid).
    fn permute(&mut self, n: Node, form: Vec<usize>) -> IR<Node> {
        if self.has_unknown(&n) {
            if matches!(n, Node::Leaf(_)) {
                return Err(value_err("Permutation size must be the length of the shape of expr"));
            }
            return Err(type_err("object of type 'NoneType' has no len()"));
        }
        self.with_ctx(|c| c.permute(n, form))
    }

    /// `_array_add(*args)` with the checks of `ArrayAdd.__new__`.
    fn add(&mut self, args: Vec<Node>) -> IR<Node> {
        let mut ranks = vec![];
        for a in &args {
            let r = self.get_rank(a)?;
            if !ranks.contains(&r) { ranks.push(r); }
        }
        if ranks.len() != 1 { return Err(value_err("summing arrays of different ranks")); }
        let mut shapes: Vec<Vec<usize>> = vec![];
        for a in &args {
            if let Node::Leaf(i) = a {
                if !self.has_shape[*i] { return Err(IErr::NoShape(*i)); }
            }
            if let Some(s) = self.shape(a)? {
                if !shapes.contains(&s) { shapes.push(s); }
            }
        }
        if shapes.len() > 1 { return Err(value_err("mismatching shapes in addition")); }
        self.with_ctx(|c| c.add(args))
    }

    fn one_array(&mut self, dims: Vec<usize>) -> Node {
        self.leaves.push(LeafInfo { rank: dims.len(), shape: Some(dims.clone()), zero: false });
        self.unknown.push(false);
        self.has_shape.push(!dims.is_empty());
        Node::Leaf(self.leaves.len() - 1)
    }

    fn map_kd(&self, i: &It, repl: &[(usize, It)]) -> It {
        match i {
            It::A(a) => repl.iter().find(|(k, _)| k == a).map_or_else(|| i.clone(), |(_, v)| v.clone()),
            o => o.clone(),
        }
    }

    /// `_convert_indexed_to_array(expr)`.
    pub fn convert(&mut self, e: &SNode) -> IR<(Node, Vec<It>)> {
        match e {
            SNode::Sum(function, variables, limits) => self.convert_sum(function, variables, limits),
            SNode::Mul(args) => {
                let mut conv = Vec::with_capacity(args.len());
                for a in args { conv.push(self.convert(a)?); }
                let mut repl: Vec<(usize, It)> = vec![];
                for a in args {
                    if let SNode::Delta(_, i, j) = a {
                        let mut k: BTreeSet<usize> = [*i, *j].into_iter().collect();
                        for x in [i, j] {
                            if let Some((_, It::S(s))) = repl.iter().find(|(key, _)| key == x) { k.extend(s.iter().copied()); }
                        }
                        let set = It::S(k.clone());
                        for idx in k {
                            match repl.iter_mut().find(|(key, _)| *key == idx) {
                                Some(entry) => entry.1 = set.clone(),
                                None => repl.push((idx, set.clone())),
                            }
                        }
                    }
                }
                let mut newargs = vec![];
                let mut flat = vec![];
                for (a, (n, ind)) in args.iter().zip(conv) {
                    if matches!(a, SNode::Delta(..)) { continue; }
                    newargs.push(n);
                    flat.extend(ind.iter().map(|j| self.map_kd(j, &repl)));
                }
                let (diag, ret) = diagonal_indices(&flat, self.int_atoms);
                let tp = self.tensor(newargs)?;
                if diag.is_empty() { Ok((tp, ret)) } else { Ok((self.diagonal(tp, diag)?, ret)) }
            },
            SNode::Element(parent, indices) => {
                let flat: Vec<It> = indices.iter().map(|&a| It::A(a)).collect();
                let (diag, ret) = diagonal_indices(&flat, self.int_atoms);
                if diag.is_empty() { Ok((Node::Leaf(*parent), ret)) } else { Ok((self.diagonal(Node::Leaf(*parent), diag)?, ret)) }
            },
            SNode::Base => Err(IErr::A(AErr::NotImplemented(String::new()))),
            SNode::Delta(leaf, i, j) => Ok((Node::Leaf(*leaf), vec![It::A(*i), It::A(*j)])),
            SNode::Add(args) => self.convert_add(args),
            SNode::Pow(base, exp) => {
                let (sub, subindices) = self.convert(base)?;
                match exp {
                    Ok(n) => {
                        let n = (*n).max(0) as usize;
                        let tp = self.tensor(vec![sub; n])?;
                        let groups: Vec<Vec<usize>> = if n == 0 { vec![] } else { (0..2).map(|r| (0..n).map(|i| 2 * i + r).collect()).collect() };
                        Ok((self.diagonal(tp, groups)?, subindices))
                    },
                    Err(leaf) => Ok((Node::Leaf(*leaf), vec![])),
                }
            },
            SNode::Func(arg, f) => {
                let arg = arg.as_ref().ok_or_else(|| IErr::A(AErr::Index("tuple index out of range".into())))?;
                let (sub, subindices) = self.convert(arg)?;
                Ok((Node::PowApply(Box::new(sub), *f), subindices))
            },
            SNode::Leaf(i) => Ok((Node::Leaf(*i), vec![])),
        }
    }

    fn convert_sum(&mut self, function: &SNode, variables: &[usize], limits: &[Limit]) -> IR<(Node, Vec<It>)> {
        let (mut sub, subindices) = self.convert(function)?;
        let sets: Vec<&BTreeSet<usize>> = subindices.iter().filter_map(|i| if let It::S(s) = i { Some(s) } else { None }).collect();
        let mut summation: Vec<It> = vec![];
        for &v in variables {
            // The last set holding `v` wins (dict comprehension).
            let it = sets.iter().rev().find(|s| s.contains(&v)).map_or(It::A(v), |s| It::S((*s).clone()));
            if !summation.contains(&it) { summation.push(it); }
        }
        let mut kronecker: BTreeSet<usize> = BTreeSet::new();
        if let SNode::Mul(args) = function {
            for a in args {
                if let SNode::Delta(_, i, j) = a {
                    if i != j { kronecker.insert(*i); kronecker.insert(*j); }
                }
            }
        }
        if let Some(shape) = self.shape(&sub)? {
            if !shape.is_empty() {
                for l in limits {
                    let (ind, start_nonzero, endp1) = l.map_err(|n| value_err(&format!("not enough values to unpack (expected 3, got {n})")))?;
                    let i = argindex(&subindices, ind).ok_or_else(|| IErr::NotFound(ind, subindices.clone()))?;
                    if start_nonzero { return Err(IErr::Mismatch(ind)); }
                    let d = *shape.get(i).ok_or_else(|| IErr::A(AErr::Index("tuple index out of range".into())))?;
                    if d != endp1 { return Err(IErr::Mismatch(ind)); }
                }
            }
        }
        let mut contraction: Vec<Vec<usize>> = vec![];
        let mut sub_ind: Vec<Option<It>> = subindices.into_iter().map(Some).collect();
        if let Node::Diagonal(inner, groups) = sub {
            let mut diag: Vec<Option<Vec<usize>>> = groups.into_iter().map(Some).collect();
            let cut = sub_ind.len().saturating_sub(diag.len());
            let dindices: Vec<Option<It>> = sub_ind[cut..].to_vec();
            sub_ind.truncate(cut);
            for index in &summation {
                if let Some(p) = dindices.iter().position(|d| d.as_ref() == Some(index)) {
                    contraction.push(diag[p].clone().expect("a diagonal group is taken once"));
                    diag[p] = None;
                }
            }
            let rest: Vec<Vec<usize>> = diag.into_iter().flatten().collect();
            sub = if rest.is_empty() { *inner } else { self.diagonal(*inner, rest)? };
        }
        let mut axes: Vec<(It, Vec<usize>)> = vec![];
        for i in 0..sub_ind.len() {
            let Some(ind) = sub_ind[i].clone() else { continue };
            let include = match &ind { It::S(s) => s.iter().all(|j| !kronecker.contains(j)), It::A(a) => !kronecker.contains(a) };
            if summation.contains(&ind) && include {
                match axes.iter_mut().find(|(k, _)| *k == ind) {
                    Some(e) => e.1.push(i),
                    None => axes.push((ind, vec![i])),
                }
                sub_ind[i] = None;
            }
        }
        for (k, v) in axes {
            let skip = match &k { It::S(s) => s.iter().any(|j| kronecker.contains(j)), It::A(a) => kronecker.contains(a) };
            if !skip { contraction.push(v); }
        }
        let free: Vec<It> = sub_ind.into_iter().flatten().collect();
        let mut ret: Vec<(usize, It)> = free.iter().map(|x| (free.iter().position(|y| y == x).unwrap(), x.clone())).collect();
        ret.sort_by_key(|(k, _)| *k);
        let out = self.contract(sub, contraction)?;
        Ok((out, ret.into_iter().map(|(_, x)| x).collect()))
    }

    fn convert_add(&mut self, terms: &[SNode]) -> IR<(Node, Vec<It>)> {
        let mut args = vec![];
        let mut indices = vec![];
        for t in terms {
            let (n, i) = self.convert(t)?;
            args.push(n);
            indices.push(i);
        }
        let mut index0: Vec<It> = vec![];
        let mut shape0: Vec<usize> = vec![];
        for (arg, ind) in args.iter().zip(&indices) {
            let missing: Vec<It> = ind.iter().filter(|i| !index0.contains(i)).cloned().collect();
            index0.extend(ind.iter().filter(|i| missing.contains(i)).cloned());
            let pos: Vec<usize> = ind.iter().enumerate().filter(|(_, e)| missing.contains(e)).map(|(i, _)| i).collect();
            if pos.is_empty() { continue; }
            let shape = self.get_shape(arg)?.ok_or_else(|| type_err("'NoneType' object is not subscriptable"))?;
            for i in pos {
                shape0.push(*shape.get(i).ok_or_else(|| IErr::A(AErr::Index("tuple index out of range".into())))?);
            }
        }
        let mut out = vec![];
        for (mut arg, ind) in args.into_iter().zip(indices) {
            let mut ind = ind;
            if ind.len() < index0.len() {
                let pos: Vec<usize> = index0.iter().enumerate().filter(|(_, e)| !ind.contains(e)).map(|(i, _)| i).collect();
                let shape: Vec<usize> = pos.iter().map(|&i| shape0[i]).collect();
                let mut new: Vec<It> = pos.iter().map(|&j| index0[j].clone()).collect();
                new.extend(ind);
                ind = new;
                let one = self.one_array(shape);
                arg = self.tensor(vec![one, arg])?;
            }
            let mut form = Vec::with_capacity(index0.len());
            for j in &index0 {
                form.push(ind.iter().position(|x| x == j).ok_or_else(|| value_err("tuple.index(x): x not in tuple"))?);
            }
            let form = permutation(&form)?;
            out.push(self.permute(arg, form)?);
        }
        Ok((self.add(out)?, index0))
    }

    /// `get_shape`: `expr.shape` if the attribute exists, else `()`.
    fn get_shape(&mut self, n: &Node) -> IR<Option<Vec<usize>>> {
        if let Node::Leaf(i) = n {
            if !self.has_shape[*i] { return Ok(Some(vec![])); }
        }
        self.shape(n)
    }

    /// `convert_indexed_to_array(expr, first_indices)` after the conversion,
    /// for a result without integer indices: `first` are the atoms of
    /// `first_indices` (empty: the result is returned as is).
    pub fn reorder(&mut self, result: Node, indices: &[It], first: &[usize]) -> IR<Node> {
        if first.is_empty() { return Ok(result); }
        let check_in = |elem: &It, list: &[It]| list.contains(elem) || matches!(elem, It::A(a) if list.iter().any(|i| i.contains(*a)));
        let repl = |a: usize| indices.iter().rev().find(|i| i.contains(a)).cloned().unwrap_or(It::A(a));
        let mut fi: Vec<It> = first.iter().map(|&a| repl(a)).collect();
        // `for i in first_indices: ... first_indices.remove(i)` (the
        // iterator advances over the shrinking list).
        let mut p = 0;
        while p < fi.len() {
            let x = fi[p].clone();
            if !check_in(&x, indices) {
                let q = fi.iter().position(|y| *y == x).expect("element of the list");
                fi.remove(q);
            }
            p += 1;
        }
        let extra: Vec<It> = indices.iter().filter(|i| !check_in(i, &fi)).cloned().collect();
        fi.extend(extra);
        let get_pos = |elem: &It| -> IR<usize> {
            if let Some(p) = fi.iter().position(|y| y == elem) { return Ok(p); }
            if let It::A(a) = elem {
                if let Some(p) = fi.iter().position(|e| e.contains(*a)) { return Ok(p); }
            }
            Err(value_err("not found"))
        };
        let mut pos = Vec::with_capacity(indices.len());
        for i in indices { pos.push(get_pos(i)?); }
        // `_af_invert(pos)`.
        let mut inv = vec![0usize; pos.len()];
        for (i, &a) in pos.iter().enumerate() {
            *inv.get_mut(a).ok_or_else(|| IErr::A(AErr::Index("list assignment index out of range".into())))? = i;
        }
        let form = permutation(&inv)?;
        if let Node::Add(args) = result {
            let mut out = vec![];
            for a in args { out.push(self.permute(a, form.clone())?); }
            self.add(out)
        } else {
            self.permute(result, form)
        }
    }
}

// ------------------------------------------------------------ array to indexed

/// An index expression of the array-to-indexed plan.
#[derive(Clone, Debug, PartialEq)]
pub enum IExpr {
    /// The `k`-th given index.
    Given(usize),
    /// The dummy `d{k}`.
    Dummy(usize),
    /// A dimension id.
    Dim(usize),
    Int(i64),
    Add(Vec<IExpr>),
    Mul(Vec<IExpr>),
    /// `a // b`.
    FloorDiv(Box<IExpr>, Box<IExpr>),
    /// `a % b`.
    Mod(Box<IExpr>, Box<IExpr>),
    /// No index (`None` in the reference's index lists).
    Missing,
}

/// The scalar expression of `convert_array_to_indexed`.
#[derive(Clone, Debug, PartialEq)]
pub enum Plan {
    /// `Mul.fromiter(...)`.
    Mul(Vec<Plan>),
    /// `Add.fromiter(...)`.
    Add(Vec<Plan>),
    /// `Sum(expr, *[(d{k}, 0, dim - 1)])`: dummy numbers and dimension ids.
    Sum(Box<Plan>, Vec<(usize, usize)>),
    /// `leaf[indices]` (`__getitem__` of an array symbol or constant array,
    /// `ArrayElement(leaf, indices)` otherwise).
    Get(usize, Vec<IExpr>),
    /// `ZeroArray(*dims)[indices]`.
    GetZero(Vec<usize>, Vec<IExpr>),
    /// `function(expr)` of an element-wise function node (its leaf).
    Apply(usize, Box<Plan>),
}

pub struct IndexPlan<'a> {
    pub leaves: &'a [LeafInfo],
    /// `ones[d]`: dimension id `d` is the integer 1.
    pub ones: &'a [bool],
    pub dummies: usize,
}

impl IndexPlan<'_> {
    fn ctx_rank(&self, n: &Node) -> usize {
        let mut o = |_: &[&Node]| -> Result<Vec<usize>, AErr> { Ok(vec![]) };
        Ctx { leaves: self.leaves, ones: self.ones, order: &mut o }.rank(n)
    }

    fn ctx_shape(&self, n: &Node) -> Option<Vec<usize>> {
        let mut o = |_: &[&Node]| -> Result<Vec<usize>, AErr> { Ok(vec![]) };
        Ctx { leaves: self.leaves, ones: self.ones, order: &mut o }.shape(n)
    }

    fn missing_shape() -> AErr { AErr::Type("'NoneType' object is not subscriptable".into()) }

    /// `_ConvertArrayToIndexed.do_convert(expr, indices)`.
    pub fn convert(&mut self, n: &Node, indices: &[IExpr]) -> Result<Plan, AErr> {
        match n {
            Node::Tensor(args) => {
                let mut out = vec![];
                let mut start = 0;
                for a in args {
                    let r = self.ctx_rank(a);
                    let lo = start.min(indices.len());
                    let hi = (start + r).min(indices.len());
                    out.push(self.convert(a, &indices[lo..hi])?);
                    start += r;
                }
                Ok(Plan::Mul(out))
            },
            Node::Contract(e, groups) => {
                let mut new = vec![IExpr::Missing; self.ctx_rank(e)];
                let bottom = self.ctx_shape(e);
                let mut limits = vec![];
                for g in groups {
                    let d = self.dummies;
                    self.dummies += 1;
                    let bottom = bottom.as_ref().ok_or_else(Self::missing_shape)?;
                    let dim = *bottom.get(g[0]).ok_or_else(|| AErr::Index("tuple index out of range".into()))?;
                    limits.push((d, dim));
                    for &i in g {
                        *new.get_mut(i).ok_or_else(|| AErr::Index("list assignment index out of range".into()))? = IExpr::Dummy(d);
                    }
                }
                let mut j = 0;
                for slot in new.iter_mut() {
                    if *slot == IExpr::Missing {
                        *slot = indices.get(j).cloned().ok_or_else(|| AErr::Index("list index out of range".into()))?;
                        j += 1;
                    }
                }
                let inner = self.convert(e, &new)?;
                Ok(Plan::Sum(Box::new(inner), limits))
            },
            Node::Diagonal(e, groups) => {
                let mut new = vec![IExpr::Missing; self.ctx_rank(e)];
                let rank = self.ctx_rank(n);
                let pos = diagonal_push_down(groups, &(0..indices.len()).collect::<Vec<_>>(), rank);
                for (p, index) in pos.into_iter().zip(indices) {
                    match p {
                        Some(Position::Group(g)) => for j in g {
                            *new.get_mut(j).ok_or_else(|| AErr::Index("list assignment index out of range".into()))? = index.clone();
                        },
                        Some(Position::Axis(j)) => {
                            *new.get_mut(j).ok_or_else(|| AErr::Index("list assignment index out of range".into()))? = index.clone();
                        },
                        None => return Err(AErr::Type("list indices must be integers or slices, not NoneType".into())),
                    }
                }
                self.convert(e, &new)
            },
            Node::Permute(e, perm) => {
                // `_apply_permutation_to_list`: new_list[perm(i)] = e.
                let mut new = vec![IExpr::Missing; perm.len()];
                for (i, x) in indices.iter().enumerate() {
                    let j = *perm.get(i).ok_or_else(|| AErr::Index("list index out of range".into()))?;
                    *new.get_mut(j).ok_or_else(|| AErr::Index("list assignment index out of range".into()))? = x.clone();
                }
                self.convert(e, &new)
            },
            Node::Add(args) => {
                let mut out = vec![];
                for a in args { out.push(self.convert(a, indices)?); }
                Ok(Plan::Add(out))
            },
            Node::PowApply(e, f) => Ok(Plan::Apply(*f, Box::new(self.convert(e, indices)?))),
            Node::Reshape(e, up) => {
                let down = self.ctx_shape(e).ok_or_else(Self::missing_shape)?;
                // cumul = accumulate([1] + reversed(shape_up), mul)
                let mut cumul: Vec<IExpr> = vec![IExpr::Int(1)];
                for &d in up.iter().rev() {
                    let last = cumul[cumul.len() - 1].clone();
                    cumul.push(if last == IExpr::Int(1) { IExpr::Dim(d) } else { IExpr::Mul(vec![last, IExpr::Dim(d)]) });
                }
                let one = IExpr::Add(indices.iter().rev().zip(&cumul).map(|(i, s)| IExpr::Mul(vec![i.clone(), s.clone()])).collect());
                let k = down.len();
                let mut dest = vec![IExpr::Missing; k];
                // `c` is the product of the dimensions passed; `c == 1` while
                // they are all 1.
                let mut c: Vec<usize> = vec![];
                for (i, &d) in down.iter().rev().enumerate() {
                    let e_ = IExpr::Dim(d);
                    let c_one = c.iter().all(|&x| self.ones.get(x).copied().unwrap_or(false));
                    let cc = IExpr::Mul(c.iter().map(|&x| IExpr::Dim(x)).collect());
                    dest[i] = if c_one {
                        if i == k - 1 { one.clone() } else { IExpr::Mod(Box::new(one.clone()), Box::new(e_)) }
                    } else if i == k - 1 {
                        IExpr::FloorDiv(Box::new(one.clone()), Box::new(cc))
                    } else {
                        IExpr::Mod(Box::new(IExpr::FloorDiv(Box::new(one.clone()), Box::new(cc))), Box::new(e_))
                    };
                    c.push(d);
                }
                dest.reverse();
                self.convert(e, &dest)
            },
            Node::Leaf(i) => Ok(Plan::Get(*i, indices.to_vec())),
            Node::Zero(d) => Ok(Plan::GetZero(d.clone(), indices.to_vec())),
        }
    }
}

/// The inverse permutation of a valid array form.
pub fn af_invert(form: &[usize]) -> Vec<usize> { invert(form) }

#[cfg(test)]
mod tests {
    use super::*;

    fn conv(leaves: Vec<LeafInfo>, ints: &[bool], e: &SNode) -> IR<(Node, Vec<It>)> {
        let n = leaves.len();
        let mut order = |nodes: &[&Node], _: &[LeafInfo]| Ok((0..nodes.len()).collect());
        let ones = vec![false; 4];
        let mut c = IndexConv { unknown: vec![false; n], has_shape: vec![true; n], given: n, leaves, ones: &ones, int_atoms: ints, order: &mut order };
        c.convert(e)
    }

    fn mat() -> LeafInfo { LeafInfo { rank: 2, shape: Some(vec![0, 0]), zero: false } }

    #[test]
    fn matrix_product_sum() {
        // Sum(M[i, j]*N[j, l], (j, 0, k-1)): atoms i=0, j=1, l=2; dim k-1+1 = 0.
        let e = SNode::Sum(Box::new(SNode::Mul(vec![SNode::Element(0, vec![0, 1]), SNode::Element(1, vec![1, 2])])), vec![1], vec![Ok((1, false, 0))]);
        let (n, ind) = conv(vec![mat(), mat()], &[false; 3], &e).unwrap();
        assert_eq!(n, Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(0), Node::Leaf(1)])), vec![vec![1, 2]]));
        assert_eq!(ind, vec![It::A(0), It::A(2)]);
        // A dimension mismatch.
        let e = SNode::Sum(Box::new(SNode::Element(0, vec![0, 0])), vec![0], vec![Ok((0, false, 3))]);
        assert_eq!(conv(vec![mat()], &[false; 3], &e).unwrap_err(), IErr::Mismatch(0));
    }

    #[test]
    fn diagonal_and_delta() {
        let e = SNode::Element(0, vec![0, 0]);
        let (n, ind) = conv(vec![mat()], &[false; 3], &e).unwrap();
        assert_eq!(n, Node::Diagonal(Box::new(Node::Leaf(0)), vec![vec![0, 1]]));
        assert_eq!(ind, vec![It::A(0)]);
        // KroneckerDelta(i, j)*M[i, k]
        let e = SNode::Mul(vec![SNode::Delta(1, 0, 1), SNode::Element(0, vec![0, 2])]);
        let (n, ind) = conv(vec![mat(), LeafInfo { rank: 0, shape: Some(vec![]), zero: false }], &[false; 3], &e).unwrap();
        assert_eq!(n, Node::Leaf(0));
        assert_eq!(ind, vec![It::S([0, 1].into_iter().collect()), It::A(2)]);
    }

    #[test]
    fn plan_contraction() {
        let leaves = vec![mat(), mat()];
        let mut p = IndexPlan { leaves: &leaves, ones: &[false], dummies: 0 };
        let n = Node::Contract(Box::new(Node::Tensor(vec![Node::Leaf(0), Node::Leaf(1)])), vec![vec![1, 2]]);
        let plan = p.convert(&n, &[IExpr::Given(0), IExpr::Given(1)]).unwrap();
        assert_eq!(plan, Plan::Sum(Box::new(Plan::Mul(vec![
            Plan::Get(0, vec![IExpr::Given(0), IExpr::Dummy(0)]),
            Plan::Get(1, vec![IExpr::Dummy(0), IExpr::Given(1)]),
        ])), vec![(0, 0)]));
    }
}
