//! Array-expression canonicalization (`symox.tensor.array.expressions`):
//! the SymPy 1.14 `_canonicalize`/`doit`/`nest_permutation` rules of
//! `crates/symox-core/src/arrayexpr.rs`.
//!
//! Wire format (the Python wrapper only encodes and decodes): a node is
//! `("L", i)`, `("Z", [dims])`, `("T", [nodes])`, `("A", [nodes])`,
//! `("P", node, [perm])`, `("C", node, [[group]])` or `("D", node, [[group]])`;
//! `leaves[i] = (rank, shape | None, is_zero_matrix)`, dimensions are ids and
//! `ones[d]` says whether dimension `d` is 1. `order(nodes)` returns the
//! stable `default_sort_key` order of the given encoded nodes.
use pyo3::exceptions::{PyAttributeError, PyIndexError, PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use symox::arrayexpr::{AErr, Ctx, LeafInfo, Node};

pub(crate) fn decode(o: &Bound<'_, PyAny>) -> PyResult<Node> {
    let t = o.downcast::<PyTuple>()?;
    let tag: String = t.get_item(0)?.extract()?;
    let b = t.get_item(1)?;
    let groups = |i: usize| -> PyResult<Vec<Vec<usize>>> { t.get_item(i)?.extract() };
    Ok(match tag.as_str() {
        "L" => Node::Leaf(b.extract()?),
        "Z" => Node::Zero(b.extract()?),
        "T" => Node::Tensor(b.try_iter()?.map(|x| decode(&x?)).collect::<PyResult<_>>()?),
        "A" => Node::Add(b.try_iter()?.map(|x| decode(&x?)).collect::<PyResult<_>>()?),
        "P" => Node::Permute(Box::new(decode(&b)?), t.get_item(2)?.extract()?),
        "C" => Node::Contract(Box::new(decode(&b)?), groups(2)?),
        "D" => Node::Diagonal(Box::new(decode(&b)?), groups(2)?),
        "R" => Node::Reshape(Box::new(decode(&b)?), t.get_item(2)?.extract()?),
        "E" => Node::PowApply(Box::new(decode(&b)?), t.get_item(2)?.extract()?),
        _ => return Err(PyValueError::new_err("unknown array-expression node")),
    })
}

pub(crate) fn encode<'py>(py: Python<'py>, n: &Node) -> PyResult<Bound<'py, PyAny>> {
    let list = |v: &[Node]| -> PyResult<Bound<'py, PyAny>> {
        Ok(PyList::new(py, v.iter().map(|x| encode(py, x)).collect::<PyResult<Vec<_>>>()?)?.into_any())
    };
    let groups = |g: &[Vec<usize>]| -> PyResult<Bound<'py, PyAny>> {
        Ok(PyList::new(py, g.iter().map(|x| PyList::new(py, x)).collect::<PyResult<Vec<_>>>()?)?.into_any())
    };
    let t = match n {
        Node::Leaf(i) => PyTuple::new(py, ["L".into_pyobject(py)?.into_any(), i.into_pyobject(py)?.into_any()])?,
        Node::Zero(s) => PyTuple::new(py, ["Z".into_pyobject(py)?.into_any(), PyList::new(py, s)?.into_any()])?,
        Node::Tensor(a) => PyTuple::new(py, ["T".into_pyobject(py)?.into_any(), list(a)?])?,
        Node::Add(a) => PyTuple::new(py, ["A".into_pyobject(py)?.into_any(), list(a)?])?,
        Node::Permute(e, p) => PyTuple::new(py, ["P".into_pyobject(py)?.into_any(), encode(py, e)?, PyList::new(py, p)?.into_any()])?,
        Node::Contract(e, g) => PyTuple::new(py, ["C".into_pyobject(py)?.into_any(), encode(py, e)?, groups(g)?])?,
        Node::Diagonal(e, g) => PyTuple::new(py, ["D".into_pyobject(py)?.into_any(), encode(py, e)?, groups(g)?])?,
        Node::Reshape(e, s) => PyTuple::new(py, ["R".into_pyobject(py)?.into_any(), encode(py, e)?, PyList::new(py, s)?.into_any()])?,
        Node::PowApply(e, x) => PyTuple::new(py, ["E".into_pyobject(py)?.into_any(), encode(py, e)?, x.into_pyobject(py)?.into_any()])?,
    };
    Ok(t.into_any())
}

pub(crate) fn to_pyerr(e: AErr, failure: Option<PyErr>) -> PyErr {
    match e {
        AErr::Value(m) => PyValueError::new_err(m),
        AErr::NotImplemented(m) => PyNotImplementedError::new_err(m),
        AErr::Type(m) => PyTypeError::new_err(m),
        AErr::Index(m) => PyIndexError::new_err(m),
        AErr::Attribute(m) => PyAttributeError::new_err(m),
        AErr::Callback => failure.unwrap_or_else(|| PyValueError::new_err("sort key failed")),
    }
}

/// Runs `f` with the `default_sort_key` oracle `order(nodes)` (a Python callable).
pub(crate) fn with_order<T>(py: Python<'_>, order: &Bound<'_, PyAny>, extra: Option<(&Bound<'_, PyAny>, &Bound<'_, PyAny>)>, f: impl FnOnce(&mut dyn FnMut(&[&Node]) -> Result<Vec<usize>, AErr>) -> Result<T, AErr>) -> PyResult<T> {
    let mut failure: Option<PyErr> = None;
    let result = {
        let mut sorter = |nodes: &[&Node]| -> Result<Vec<usize>, AErr> {
            let run = || -> PyResult<Vec<usize>> {
                let enc = PyList::new(py, nodes.iter().map(|n| encode(py, n)).collect::<PyResult<Vec<_>>>()?)?;
                match extra {
                    Some((leaves, dims)) => order.call1((enc, leaves, dims))?.extract(),
                    None => order.call1((enc,))?.extract(),
                }
            };
            run().map_err(|e| { failure = Some(e); AErr::Callback })
        };
        f(&mut sorter)
    };
    result.map_err(|e| to_pyerr(e, failure))
}

/// `kind`: "canon" (`_canonicalize` of the top node), "doit" (deep) or
/// "nest" (`nest_permutation`).
#[pyfunction]
#[pyo3(signature = (kind, tree, leaves, ones, order, allow_trivial=false))]
fn _arrayexpr_canon(
    py: Python<'_>, kind: &str, tree: &Bound<'_, PyAny>, leaves: Vec<(usize, Option<Vec<usize>>, bool)>,
    ones: Vec<bool>, order: &Bound<'_, PyAny>, allow_trivial: bool,
) -> PyResult<PyObject> {
    let node = decode(tree)?;
    if !matches!(kind, "canon" | "doit" | "nest") { return Err(PyValueError::new_err("unknown request")); }
    let leaves: Vec<LeafInfo> = leaves.into_iter().map(|(rank, shape, zero)| LeafInfo { rank, shape, zero }).collect();
    let out = with_order(py, order, None, |sorter| {
        let mut ctx = Ctx { leaves: &leaves, ones: &ones, order: sorter };
        match kind {
            "canon" => ctx.canonicalize(node, allow_trivial),
            "doit" => ctx.doit(node, true),
            _ => ctx.nest_permutation(node),
        }
    })?;
    Ok(encode(py, &out)?.unbind())
}

/// `convert_matrix_to_array(expr)`: `(tree, leaves, dims)` with the leaf
/// objects and the dimension of each dimension id (`order(nodes, leaves,
/// dims)` sorts), or None when `expr` is
/// neither a native matrix nor a native scalar (SymPy returns it unchanged).
#[pyfunction]
fn _arrayexpr_from_matrix(py: Python<'_>, expr: &Bound<'_, PyAny>, order: &Bound<'_, PyAny>) -> PyResult<Option<PyObject>> {
    use symox::arrayexpr_matrix::{Conversion, Source};
    let mut c = Conversion::default();
    let tree = if let Ok(m) = expr.downcast::<crate::matrixexpr::PyMatrixExpr>() {
        c.matrix(&m.get().m)
    } else if let Ok(e) = expr.downcast::<crate::PyExprNode>() {
        c.scalar(&e.get().e)
    } else {
        return Ok(None);
    };
    let tree = tree.map_err(|e| to_pyerr(e, None))?;
    let leaves = PyList::empty(py);
    for s in &c.sources {
        match s {
            Source::Matrix(m) => leaves.append(crate::matrixexpr::matrix_object(py, m)?)?,
            Source::Scalar(e) => leaves.append(crate::expr_object(py, e.clone())?)?,
        }
    }
    let dims = PyList::new(py, c.dims.iter().map(|d| crate::expr_object(py, d.clone())).collect::<PyResult<Vec<_>>>()?)?;
    let out = with_order(py, order, Some((leaves.as_any(), dims.as_any())), |sorter| c.evaluate(tree, sorter))?;
    Ok(Some(PyTuple::new(py, [encode(py, &out)?, leaves.into_any(), dims.into_any()])?.into_any().unbind()))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_arrayexpr_canon, m)?)?;
    m.add_function(wrap_pyfunction!(_arrayexpr_from_matrix, m)?)?;
    Ok(())
}
