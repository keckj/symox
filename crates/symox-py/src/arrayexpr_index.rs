//! Bindings of `array_derive` (`crates/symox-core/src/arrayexpr_derive.rs`),
//! `_convert_indexed_to_array`/`convert_indexed_to_array` and the plan of
//! `convert_array_to_indexed` (`crates/symox-core/src/arrayexpr_index.rs`).
//!
//! The Python wrappers (`symox/tensor/array/expressions/_native.py`) encode
//! the expressions as tuples and rebuild the objects; array nodes use the
//! wire format of `arrayexpr_canon.rs`.
use crate::arrayexpr_canon::{decode, encode, to_pyerr, with_order};
use pyo3::exceptions::{PyAttributeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use symox::arrayexpr::{AErr, LeafInfo, Node};
use symox::arrayexpr_derive::{array_derive, DErr, DNode, XInfo};
use symox::arrayexpr_index::{IErr, IExpr, IndexConv, IndexPlan, It, Limit, Plan, SNode};

type Leaves = Vec<(usize, Option<Vec<usize>>, bool)>;

fn leaf_infos(leaves: Leaves) -> Vec<LeafInfo> {
    leaves.into_iter().map(|(rank, shape, zero)| LeafInfo { rank, shape, zero }).collect()
}

/// `(True, value)` or `(False, error index)`.
fn result<'py, T: FromPyObject<'py>>(o: &Bound<'py, PyAny>) -> PyResult<Result<T, usize>> {
    let (ok, v): (bool, Bound<'py, PyAny>) = o.extract()?;
    Ok(if ok { Ok(v.extract()?) } else { Err(v.extract()?) })
}

fn dnode(o: &Bound<'_, PyAny>) -> PyResult<DNode> {
    let t = o.downcast::<PyTuple>()?;
    let tag: String = t.get_item(0)?.extract()?;
    let item = |i: usize| t.get_item(i);
    let sub = |i: usize| -> PyResult<Box<DNode>> { Ok(Box::new(dnode(&t.get_item(i)?)?)) };
    Ok(match tag.as_str() {
        "T" => DNode::Tensor(
            item(1)?.try_iter()?.map(|x| dnode(&x?)).collect::<PyResult<_>>()?,
            item(2)?.try_iter()?.map(|x| decode(&x?)).collect::<PyResult<_>>()?,
        ),
        "A" => DNode::Add(item(1)?.try_iter()?.map(|x| dnode(&x?)).collect::<PyResult<_>>()?),
        "P" => DNode::Permute(sub(1)?, item(2)?.extract()?),
        "C" => DNode::Contract(sub(1)?, item(2)?.extract()?),
        "D" => DNode::Diagonal(sub(1)?, item(2)?.extract()?),
        "R" => DNode::Reshape(sub(1)?, item(2)?.extract()?),
        "S" => DNode::Scalar,
        "ZC" => DNode::ZeroCat(item(1)?.extract()?),
        "SA" => DNode::SelfArray(item(1)?.extract()?),
        "SM" => DNode::SelfMatrix(item(1)?.extract()?, item(2)?.extract()?),
        "TR" => DNode::Transpose(sub(1)?),
        "INV" => DNode::Inverse(sub(1)?, item(2)?.extract()?, item(3)?.extract()?),
        "EWM" => DNode::EwMatrix { assert_err: item(1)?.extract()?, applied: result(&item(2)?)?, sub: sub(3)? },
        "EWA" => DNode::EwArray { applied: result(&item(1)?)?, sub: sub(2)?, rank: result(&item(3)?)? },
        "CV" => DNode::Converted(sub(1)?),
        "F" => DNode::Fail(item(1)?.extract()?),
        _ => return Err(PyValueError::new_err("unknown derivative node")),
    })
}

fn no_shape(objs: &Bound<'_, PyList>, i: usize) -> PyErr {
    let name = objs.get_item(i).and_then(|o| Ok(o.get_type().name()?.to_string())).unwrap_or_else(|_| "object".into());
    PyAttributeError::new_err(format!("'{name}' object has no attribute 'shape'"))
}

/// `array_derive(expr, x)`: `tree` the derivation nodes, `objs`/`leaves` the
/// leaf objects and their `(rank, shape, is_zero_matrix)`, `has_shape` the
/// leaves with a `shape` attribute, `szero` the leaf `S.Zero`, `x` the
/// results `(x.shape, get_shape(x), get_rank(x))`, `errors` the captured
/// exceptions. Returns the encoded array node.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn _arrayexpr_derive(
    py: Python<'_>, tree: &Bound<'_, PyAny>, objs: &Bound<'_, PyList>, leaves: Leaves, has_shape: Vec<bool>,
    ones: Vec<bool>, szero: usize, x: &Bound<'_, PyTuple>, errors: &Bound<'_, PyList>, order: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let tree = dnode(tree)?;
    let leaves = leaf_infos(leaves);
    let xi = XInfo { shape: result(&x.get_item(0)?)?, get_shape: result(&x.get_item(1)?)?, rank: result(&x.get_item(2)?)? };
    let mut derr: Option<DErr> = None;
    let out = with_order(py, order, None, |sorter| {
        array_derive(&tree, &leaves, &has_shape, &ones, szero, &xi, sorter).map_err(|e| match e {
            DErr::A(a) => a,
            other => { derr = Some(other); AErr::Value(String::new()) },
        })
    });
    match (out, derr) {
        (Ok(n), _) => Ok(encode(py, &n)?.unbind()),
        (Err(_), Some(DErr::Py(i))) => Err(PyErr::from_value(errors.get_item(i)?)),
        (Err(_), Some(DErr::NoShape(i))) => Err(no_shape(objs, i)),
        (Err(e), _) => Err(e),
    }
}

fn snode(o: &Bound<'_, PyAny>) -> PyResult<SNode> {
    let t = o.downcast::<PyTuple>()?;
    let tag: String = t.get_item(0)?.extract()?;
    let item = |i: usize| t.get_item(i);
    let list = |i: usize| -> PyResult<Vec<SNode>> { item(i)?.try_iter()?.map(|x| snode(&x?)).collect() };
    Ok(match tag.as_str() {
        "S" => {
            let limits: Vec<Limit> = item(3)?.try_iter()?.map(|l| result::<(usize, bool, usize)>(&l?)).collect::<PyResult<_>>()?;
            SNode::Sum(Box::new(snode(&item(1)?)?), item(2)?.extract()?, limits)
        },
        "M" => SNode::Mul(list(1)?),
        "A" => SNode::Add(list(1)?),
        "E" => SNode::Element(item(1)?.extract()?, item(2)?.extract()?),
        "K" => SNode::Delta(item(1)?.extract()?, item(2)?.extract()?, item(3)?.extract()?),
        "P" => SNode::Pow(Box::new(snode(&item(1)?)?), result(&item(2)?)?),
        "F" => {
            let a = item(1)?;
            SNode::Func(if a.is_none() { None } else { Some(Box::new(snode(&a)?)) }, item(2)?.extract()?)
        },
        "B" => SNode::Base,
        "L" => SNode::Leaf(item(1)?.extract()?),
        _ => return Err(PyValueError::new_err("unknown indexed node")),
    })
}

fn encode_it<'py>(py: Python<'py>, i: &It) -> PyResult<Bound<'py, PyAny>> {
    Ok(match i {
        It::A(a) => a.into_pyobject(py)?.into_any(),
        It::S(s) => PyList::new(py, s.iter())?.into_any(),
    })
}

/// `_convert_indexed_to_array(expr)` (`first` None) or
/// `convert_indexed_to_array(expr, first)`: returns `("ok", node, indices,
/// one_arrays)`, `("element", node, indices, one_arrays)` when the indices
/// have integers (the caller builds `ArrayElement(result, indices)`), or
/// `("error", kind, payload)` for the errors whose message needs the
/// objects. `one_arrays` are the dimension ids of the `OneArray` leaves
/// appended after the given ones; `order(nodes, one_arrays)` sorts.
#[pyfunction]
#[pyo3(signature = (tree, leaves, unknown, has_shape, ones, int_atoms, first, order))]
#[allow(clippy::too_many_arguments)]
fn _arrayexpr_from_indexed(
    py: Python<'_>, tree: &Bound<'_, PyAny>, leaves: Leaves, unknown: Vec<bool>, has_shape: Vec<bool>,
    ones: Vec<bool>, int_atoms: Vec<bool>, first: Option<Vec<usize>>, order: &Bound<'_, PyAny>,
) -> PyResult<PyObject> {
    let tree = snode(tree)?;
    let leaves = leaf_infos(leaves);
    let given = leaves.len();
    let mut failure: Option<PyErr> = None;
    let mut sorter = |nodes: &[&Node], table: &[LeafInfo]| -> Result<Vec<usize>, AErr> {
        let run = || -> PyResult<Vec<usize>> {
            let enc = PyList::new(py, nodes.iter().map(|n| encode(py, n)).collect::<PyResult<Vec<_>>>()?)?;
            let extra = PyList::new(py, table[given..].iter().map(|l| l.shape.clone().unwrap_or_default()))?;
            order.call1((enc, extra))?.extract()
        };
        run().map_err(|e| { failure = Some(e); AErr::Callback })
    };
    let mut conv = IndexConv { leaves, unknown, has_shape, given, ones: &ones, int_atoms: &int_atoms, order: &mut sorter };
    let out = (|| -> Result<(&'static str, Node, Vec<It>), IErr> {
        let (node, indices) = conv.convert(&tree)?;
        if indices.iter().any(|i| matches!(i, It::A(a) if int_atoms[*a])) {
            return Ok(("element", node, indices));
        }
        match &first {
            Some(f) => Ok(("ok", conv.reorder(node, &indices, f)?, indices)),
            None => Ok(("ok", node, indices)),
        }
    })();
    let extra: Vec<Vec<usize>> = conv.leaves[given..].iter().map(|l| l.shape.clone().unwrap_or_default()).collect();
    drop(conv);
    match out {
        Ok((tag, node, indices)) => {
            let ind = PyList::new(py, indices.iter().map(|i| encode_it(py, i)).collect::<PyResult<Vec<_>>>()?)?;
            Ok(PyTuple::new(py, [tag.into_pyobject(py)?.into_any(), encode(py, &node)?, ind.into_any(), PyList::new(py, extra)?.into_any()])?.into_any().unbind())
        },
        Err(IErr::A(e)) => Err(to_pyerr(e, failure)),
        Err(IErr::NotFound(a, sub)) => {
            let sub = PyList::new(py, sub.iter().map(|i| encode_it(py, i)).collect::<PyResult<Vec<_>>>()?)?;
            let payload = PyTuple::new(py, [a.into_pyobject(py)?.into_any(), sub.into_any()])?;
            Ok(PyTuple::new(py, ["error".into_pyobject(py)?.into_any(), "notfound".into_pyobject(py)?.into_any(), payload.into_any()])?.into_any().unbind())
        },
        Err(IErr::Mismatch(a)) => Ok(PyTuple::new(py, ["error".into_pyobject(py)?.into_any(), "mismatch".into_pyobject(py)?.into_any(), a.into_pyobject(py)?.into_any()])?.into_any().unbind()),
        Err(IErr::NoShape(i)) => Ok(PyTuple::new(py, ["error".into_pyobject(py)?.into_any(), "noshape".into_pyobject(py)?.into_any(), i.into_pyobject(py)?.into_any()])?.into_any().unbind()),
    }
}

fn encode_iexpr<'py>(py: Python<'py>, e: &IExpr) -> PyResult<Bound<'py, PyAny>> {
    let pair = |tag: &str, v: Bound<'py, PyAny>| -> PyResult<Bound<'py, PyAny>> { Ok(PyTuple::new(py, [tag.into_pyobject(py)?.into_any(), v])?.into_any()) };
    let list = |v: &[IExpr]| -> PyResult<Bound<'py, PyAny>> { Ok(PyList::new(py, v.iter().map(|x| encode_iexpr(py, x)).collect::<PyResult<Vec<_>>>()?)?.into_any()) };
    match e {
        IExpr::Given(k) => pair("g", k.into_pyobject(py)?.into_any()),
        IExpr::Dummy(k) => pair("d", k.into_pyobject(py)?.into_any()),
        IExpr::Dim(k) => pair("D", k.into_pyobject(py)?.into_any()),
        IExpr::Int(k) => pair("i", k.into_pyobject(py)?.into_any()),
        IExpr::Add(v) => pair("+", list(v)?),
        IExpr::Mul(v) => pair("*", list(v)?),
        IExpr::FloorDiv(a, b) => Ok(PyTuple::new(py, ["//".into_pyobject(py)?.into_any(), encode_iexpr(py, a)?, encode_iexpr(py, b)?])?.into_any()),
        IExpr::Mod(a, b) => Ok(PyTuple::new(py, ["%".into_pyobject(py)?.into_any(), encode_iexpr(py, a)?, encode_iexpr(py, b)?])?.into_any()),
        IExpr::Missing => Ok(py.None().into_bound(py)),
    }
}

fn encode_plan<'py>(py: Python<'py>, p: &Plan) -> PyResult<Bound<'py, PyAny>> {
    let s = |x: &str| -> PyResult<Bound<'py, PyAny>> { Ok(x.into_pyobject(py)?.into_any()) };
    let plans = |v: &[Plan]| -> PyResult<Bound<'py, PyAny>> { Ok(PyList::new(py, v.iter().map(|x| encode_plan(py, x)).collect::<PyResult<Vec<_>>>()?)?.into_any()) };
    let idx = |v: &[IExpr]| -> PyResult<Bound<'py, PyAny>> { Ok(PyList::new(py, v.iter().map(|x| encode_iexpr(py, x)).collect::<PyResult<Vec<_>>>()?)?.into_any()) };
    Ok(match p {
        Plan::Mul(v) => PyTuple::new(py, [s("M")?, plans(v)?])?.into_any(),
        Plan::Add(v) => PyTuple::new(py, [s("A")?, plans(v)?])?.into_any(),
        Plan::Sum(e, limits) => PyTuple::new(py, [s("S")?, encode_plan(py, e)?, PyList::new(py, limits)?.into_any()])?.into_any(),
        Plan::Get(l, v) => PyTuple::new(py, [s("G")?, l.into_pyobject(py)?.into_any(), idx(v)?])?.into_any(),
        Plan::GetZero(d, v) => PyTuple::new(py, [s("Z")?, PyList::new(py, d)?.into_any(), idx(v)?])?.into_any(),
        Plan::Apply(f, e) => PyTuple::new(py, [s("F")?, f.into_pyobject(py)?.into_any(), encode_plan(py, e)?])?.into_any(),
    })
}

/// `convert_array_to_indexed(expr, indices)`: the plan of the scalar
/// expression for `n_indices` given indices.
#[pyfunction]
fn _arrayexpr_to_indexed(py: Python<'_>, tree: &Bound<'_, PyAny>, leaves: Leaves, ones: Vec<bool>, n_indices: usize) -> PyResult<PyObject> {
    let tree = decode(tree)?;
    let leaves = leaf_infos(leaves);
    let mut p = IndexPlan { leaves: &leaves, ones: &ones, dummies: 0 };
    let given: Vec<IExpr> = (0..n_indices).map(IExpr::Given).collect();
    let plan = p.convert(&tree, &given).map_err(|e| to_pyerr(e, None))?;
    Ok(encode_plan(py, &plan)?.unbind())
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_arrayexpr_derive, m)?)?;
    m.add_function(wrap_pyfunction!(_arrayexpr_from_indexed, m)?)?;
    m.add_function(wrap_pyfunction!(_arrayexpr_to_indexed, m)?)?;
    Ok(())
}
