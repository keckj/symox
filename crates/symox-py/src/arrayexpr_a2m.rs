//! `convert_array_to_matrix`, its passes and the `_EditArrayContraction`
//! editor (`crates/symox-core/src/arrayexpr_a2m.rs`).
//!
//! Wire format: the trees of `arrayexpr_canon` (plus `("F", node, fn_leaf)`
//! for an element-wise application) over a leaf table of descriptors:
//! `("m", matrix)`, `("s", scalar)`, `("one", [dim ids])`,
//! `("a", object, rank, [dim ids] | None)` (an opaque array leaf) and
//! `("fn", lambda, original)`. `dims` holds the dimension of each id.
//! Results come back with the extended leaf table and dimension list; the
//! leaves the native port creates are `("m", ...)`, `("s", ...)` or
//! `("one", ...)`. `order(nodes, leaves, dims)` is the `default_sort_key`
//! oracle. `ArrayExprUnsupported` means the input is outside the native
//! port: the caller runs the stock algorithm.
use std::cell::RefCell;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use symox::arrayexpr::{AErr, LeafInfo, Node};
use symox::arrayexpr_a2m::{combine_removed, trace_normalized_arg, A2m, ArgE, Editor, Src};
use symox::expr::Expr;
use symox::matrixexpr::MatrixKind;

use crate::arrayexpr_canon::{decode, encode, to_pyerr};

type Parsed = (Vec<Src>, Vec<LeafInfo>, Vec<Expr>);

fn parse_leaves(descs: &Bound<'_, PyList>, dims_in: &Bound<'_, PyList>) -> PyResult<Parsed> {
    let mut dims: Vec<Expr> = dims_in.iter().map(|d| crate::to_expr(&d)).collect::<PyResult<_>>()?;
    let mut dim = |e: &Expr| -> usize {
        match dims.iter().position(|d| d == e) {
            Some(i) => i,
            None => { dims.push(e.clone()); dims.len() - 1 },
        }
    };
    let (mut srcs, mut info) = (vec![], vec![]);
    for d in descs.iter() {
        let t = d.downcast::<PyTuple>()?;
        let tag: String = t.get_item(0)?.extract()?;
        match tag.as_str() {
            "m" => {
                let m = crate::matrixexpr::to_matrixexpr(&t.get_item(1)?)?;
                let shape = vec![dim(m.rows()), dim(m.cols())];
                let zero = matches!(m.kind(), MatrixKind::Zero);
                srcs.push(Src::Matrix(m));
                info.push(LeafInfo { rank: 2, shape: Some(shape), zero });
            },
            "s" => {
                srcs.push(Src::Scalar(crate::to_expr(&t.get_item(1)?)?));
                info.push(LeafInfo { rank: 0, shape: Some(vec![]), zero: false });
            },
            "one" => {
                let ids: Vec<usize> = t.get_item(1)?.extract()?;
                info.push(LeafInfo { rank: ids.len(), shape: Some(ids.clone()), zero: false });
                srcs.push(Src::One(ids));
            },
            "a" => {
                let rank: usize = t.get_item(2)?.extract()?;
                let shape: Option<Vec<usize>> = t.get_item(3)?.extract()?;
                srcs.push(Src::Opaque);
                info.push(LeafInfo { rank, shape, zero: false });
            },
            "fn" => {
                srcs.push(Src::Fn(crate::to_expr(&t.get_item(1)?)?));
                info.push(LeafInfo { rank: 0, shape: Some(vec![]), zero: false });
            },
            _ => return Err(PyValueError::new_err("unknown array-expression leaf")),
        }
    }
    Ok((srcs, info, dims))
}

/// Appends to `out` the descriptors of `srcs[out.len()..]`.
fn sync_leaves(py: Python<'_>, descs: &Bound<'_, PyList>, srcs: &[Src], out: &Bound<'_, PyList>) -> PyResult<()> {
    for i in out.len()..srcs.len() {
        if i < descs.len() {
            out.append(descs.get_item(i)?)?;
            continue;
        }
        let d = match &srcs[i] {
            Src::Matrix(m) => PyTuple::new(py, ["m".into_pyobject(py)?.into_any(), crate::matrixexpr::matrix_object(py, m)?])?,
            Src::Scalar(e) => PyTuple::new(py, ["s".into_pyobject(py)?.into_any(), crate::expr_object(py, e.clone())?.into_bound(py)])?,
            Src::One(ids) => PyTuple::new(py, ["one".into_pyobject(py)?.into_any(), PyList::new(py, ids)?.into_any()])?,
            _ => return Err(PyValueError::new_err("array-expression leaf without a descriptor")),
        };
        out.append(d)?;
    }
    Ok(())
}

fn sync_dims(py: Python<'_>, dims: &[Expr], out: &Bound<'_, PyList>) -> PyResult<()> {
    for d in &dims[out.len()..] { out.append(crate::expr_object(py, d.clone())?)?; }
    Ok(())
}

enum Out {
    Tree(Node),
    Removed(Node, Vec<usize>),
    Editor(Editor, Option<Py<PyAny>>),
    EditorTree(Editor, Node),
}

type EdState = (Vec<Vec<Option<i64>>>, i64, Option<Vec<Vec<i64>>>);

fn editor_from(trees: &[Option<Node>], state: EdState) -> Editor {
    let (indices, n_contr, track) = state;
    let args: Vec<ArgE> = indices.into_iter().enumerate()
        .map(|(id, ind)| ArgE { id, element: trees.get(id).cloned().flatten(), indices: ind }).collect();
    let next_id = args.len();
    Editor { args, n_contr, track, next_id }
}

fn py_pairs(py: Python<'_>, v: &[(usize, usize)]) -> PyResult<Py<PyAny>> {
    Ok(PyList::new(py, v.iter().map(|&(a, b)| PyTuple::new(py, [a, b])).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind())
}

fn run_op(py: Python<'_>, a: &mut A2m<'_>, op: &str, trees: Vec<Option<Node>>, extra: Option<&Bound<'_, PyAny>>) -> PyResult<Result<Out, AErr>> {
    let first = || trees.first().cloned().flatten().ok_or_else(|| PyValueError::new_err("missing array expression"));
    if let Some(method) = op.strip_prefix("ed_") {
        let Some(extra) = extra else { return Err(PyValueError::new_err("missing editor state")) };
        let t = extra.downcast::<PyTuple>()?;
        if method == "new" {
            let base = first()?;
            return Ok(a.editor(&base).map(|e| Out::Editor(e, None)));
        }
        let state: EdState = (t.get_item(0)?.extract()?, t.get_item(1)?.extract()?, t.get_item(2)?.extract()?);
        let param = t.get_item(3)?;
        let mut ed = editor_from(&trees, state);
        let id = |p: &Bound<'_, PyAny>| -> PyResult<usize> { Ok(p.extract::<Option<usize>>()?.unwrap_or(usize::MAX)) };
        let res: Result<Option<Py<PyAny>>, AErr> = match method {
            "to_array" => return Ok(a.to_array(&mut ed).map(|n| Out::EditorTree(ed, n))),
            "merge" => a.merge_scalars(&mut ed).map(|_| None),
            "refresh" => { ed.refresh_indices(); Ok(None) },
            "new_index" => Ok(Some(ed.new_index().into_pyobject(py)?.into_any().unbind())),
            "contraction_indices" => match ed.contraction_indices() {
                Ok(v) => Ok(Some(PyList::new(py, v.iter().map(|g| PyList::new(py, g)).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind())),
                Err(e) => Err(e),
            },
            "mapping" => match ed.mapping_for_index(param.extract()?) { Ok(v) => Ok(Some(py_pairs(py, &v)?)), Err(e) => Err(e) },
            "ind_rel_pos" => match ed.ind_rel_pos() {
                Ok(v) => Ok(Some(PyList::new(py, v.iter().map(|g| py_pairs(py, g)).collect::<PyResult<Vec<_>>>()?)?.into_any().unbind())),
                Err(e) => Err(e),
            },
            "count" => Ok(Some(ed.count_with(param.extract()?).into_pyobject(py)?.into_any().unbind())),
            "args_with" => Ok(Some(PyList::new(py, ed.args_with(param.extract()?))?.into_any().unbind())),
            "n_diag" => Ok(Some(ed.n_diag().into_pyobject(py)?.into_any().unbind())),
            "track_start" => { ed.track_start(); Ok(None) },
            "track_merge" => {
                let (d, f): (usize, usize) = param.extract()?;
                ed.track_merge(d, f).map(|_| None)
            },
            "free_range" => match ed.abs_free_range(id(&param)?) { Ok(r) => Ok(Some(PyTuple::new(py, [r.0, r.1])?.into_any().unbind())), Err(e) => Err(e) },
            "range" => match ed.abs_range(id(&param)?) { Ok(r) => Ok(Some(PyTuple::new(py, [r.0, r.1])?.into_any().unbind())), Err(e) => Err(e) },
            _ => return Err(PyValueError::new_err("unknown editor request")),
        };
        return Ok(res.map(|r| Out::Editor(ed, r)));
    }
    Ok(match op {
        "convert" => a.convert(first()?).map(Out::Tree),
        "a2m" => a.a2m(first()?).map(Out::Tree),
        "rtd" => a.rtd(first()?).map(|(n, r)| Out::Removed(n, r)),
        "hadamard" => a.identify_hadamard_products(&first()?).map(Out::Tree),
        "removable" => a.identify_removable_identity_matrices(&first()?).map(Out::Tree),
        "remove_identity" => a.remove_identity_matrices(&first()?).map(|(n, r)| Out::Removed(n, r)),
        "multi_identity" => a.contraction_to_diagonal_multiple_identity(&first()?).map(|(n, r)| Out::Removed(n, r)),
        "diag2contr" => a.diag2contr_entry(first()?).map(Out::Tree),
        "split" => a.split_multiple_contractions(first()?).map(Out::Tree),
        "flatten" => a.flatten_contraction_of_diagonal(first()?).map(Out::Tree),
        "tp1" => {
            let ci: Vec<Vec<usize>> = extra.ok_or_else(|| PyValueError::new_err("missing contraction indices"))?.extract()?;
            let args: Option<Vec<Node>> = trees.into_iter().collect();
            let args = args.ok_or_else(|| PyValueError::new_err("missing array expression"))?;
            a.tp1(ci, args).map(Out::Tree)
        },
        _ => return Err(PyValueError::new_err("unknown request")),
    })
}

fn encode_opt<'py>(py: Python<'py>, n: &Option<Node>) -> PyResult<Bound<'py, PyAny>> {
    match n { Some(n) => encode(py, n), None => Ok(py.None().into_bound(py)) }
}

fn editor_out<'py>(py: Python<'py>, ed: &Editor, n_in: usize) -> PyResult<Bound<'py, PyAny>> {
    let args = PyList::empty(py);
    for a in &ed.args {
        let orig = if a.id < n_in { a.id.into_pyobject(py)?.into_any() } else { py.None().into_bound(py) };
        let ind = PyList::new(py, &a.indices)?;
        args.append(PyTuple::new(py, [orig, encode_opt(py, &a.element)?, ind.into_any()])?)?;
    }
    Ok(PyTuple::new(py, [args.into_any(), ed.n_contr.into_pyobject(py)?.into_any(), ed.track.clone().into_pyobject(py)?.into_any()])?.into_any())
}

/// `op` is a pass name (`convert`, `a2m`, `rtd`, `hadamard`, `removable`,
/// `remove_identity`, `multi_identity`, `diag2contr`, `split`, `flatten`,
/// `tp1`, `combine_removed`) or an editor method (`ed_*`, with
/// `extra = (indices, n_contr, track, param)` and the element trees in
/// `trees`). Returns `(result, leaves, dims)`; an editor method's result is
/// `((args, n_contr, track), value)` where each arg is
/// `(input position | None, tree | None, indices)`.
#[pyfunction]
#[pyo3(signature = (op, trees, leaves, dims, order, extra=None))]
fn _arrayexpr_a2m(
    py: Python<'_>, op: &str, trees: &Bound<'_, PyList>, leaves: &Bound<'_, PyList>, dims: &Bound<'_, PyList>,
    order: &Bound<'_, PyAny>, extra: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyObject> {
    if op == "combine_removed" {
        let (r1, r2): (Vec<usize>, Vec<usize>) = extra.ok_or_else(|| PyValueError::new_err("missing arguments"))?.extract()?;
        return Ok(PyList::new(py, combine_removed(&r1, &r2))?.into_any().unbind());
    }
    let nodes: Vec<Option<Node>> = trees.iter().map(|t| if t.is_none() { Ok(None) } else { decode(&t).map(Some) }).collect::<PyResult<_>>()?;
    let n_orig = if op == "ed_new" { 0 } else { nodes.len() };
    let (srcs, info, dimv) = parse_leaves(leaves, dims)?;
    let pleaves = PyList::empty(py);
    let pdims = PyList::empty(py);
    let failure: RefCell<Option<PyErr>> = RefCell::new(None);
    let mut order_fn = |nodes: &[&Node], srcs: &[Src], dims: &[Expr]| -> Result<Vec<usize>, AErr> {
        let run = || -> PyResult<Vec<usize>> {
            sync_leaves(py, leaves, srcs, &pleaves)?;
            sync_dims(py, dims, &pdims)?;
            let enc = PyList::new(py, nodes.iter().map(|n| encode(py, n)).collect::<PyResult<Vec<_>>>()?)?;
            order.call1((enc, &pleaves, &pdims))?.extract()
        };
        run().map_err(|e| { *failure.borrow_mut() = Some(e); AErr::Callback })
    };
    let (res, srcs, dimv) = {
        let mut a = A2m::new(srcs, info, dimv, &mut order_fn);
        let nodes = nodes.into_iter().map(|n| n.map(|n| a.unify(n))).collect();
        let r = run_op(py, &mut a, op, nodes, extra)?;
        (r, a.srcs, a.dims)
    };
    let out = res.map_err(|e| to_pyerr(e, failure.borrow_mut().take()))?;
    sync_leaves(py, leaves, &srcs, &pleaves)?;
    sync_dims(py, &dimv, &pdims)?;
    let value = match out {
        Out::Tree(n) => encode(py, &n)?,
        Out::Removed(n, r) => PyTuple::new(py, [encode(py, &n)?, PyList::new(py, r)?.into_any()])?.into_any(),
        Out::Editor(ed, v) => {
            let v = v.map_or_else(|| py.None().into_bound(py), |v| v.into_bound(py));
            PyTuple::new(py, [editor_out(py, &ed, n_orig)?, v])?.into_any()
        },
        Out::EditorTree(ed, n) => PyTuple::new(py, [editor_out(py, &ed, n_orig)?, encode(py, &n)?])?.into_any(),
    };
    Ok(PyTuple::new(py, [value, pleaves.into_any(), pdims.into_any()])?.into_any().unbind())
}

/// `Trace(m)._normalize()` for a native matrix expression.
#[pyfunction]
fn _matrix_trace_normalize(py: Python<'_>, m: &Bound<'_, PyAny>) -> PyResult<PyObject> {
    let m = crate::matrixexpr::to_matrixexpr(m)?;
    let t = trace_normalized_arg(&m).and_then(|a| a.trace_unevaluated().map_err(AErr::from)).map_err(|e| to_pyerr(e, None))?;
    crate::expr_object(py, t)
}

/// `ArrayContraction._push_indices_up/down` (`kind` "c_up", "c_down") and
/// `ArrayDiagonal._push_indices_up/down` (`kind` "d_up", "d_down", with
/// `rank`) on a flat list of axes. A diagonal group comes back as a tuple.
#[pyfunction]
#[pyo3(signature = (kind, groups, indices, rank=0))]
fn _arrayexpr_push_indices(py: Python<'_>, kind: &str, groups: Vec<Vec<usize>>, indices: Vec<usize>, rank: usize) -> PyResult<PyObject> {
    use pyo3::IntoPyObjectExt;
    use symox::arrayexpr::{contraction_push_down, contraction_push_up, diagonal_push_down, diagonal_push_up, Position};
    let out: Vec<PyObject> = match kind {
        "c_up" => contraction_push_up(&groups, &[indices]).remove(0).into_iter()
            .map(|v| v.into_py_any(py)).collect::<PyResult<_>>()?,
        "c_down" => contraction_push_down(&groups, &[indices]).remove(0).into_iter()
            .map(|v| v.into_py_any(py)).collect::<PyResult<_>>()?,
        "d_up" => diagonal_push_up(&groups, &indices, rank).into_iter()
            .map(|v| v.into_py_any(py)).collect::<PyResult<_>>()?,
        "d_down" => diagonal_push_down(&groups, &indices, rank).into_iter().map(|p| match p {
            None => Ok(py.None()),
            Some(Position::Axis(a)) => a.into_py_any(py),
            Some(Position::Group(g)) => PyTuple::new(py, g)?.into_py_any(py),
        }).collect::<PyResult<_>>()?,
        _ => return Err(PyValueError::new_err("unknown request")),
    };
    Ok(PyTuple::new(py, out)?.into_any().unbind())
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_arrayexpr_a2m, m)?)?;
    m.add_function(wrap_pyfunction!(_arrayexpr_push_indices, m)?)?;
    m.add_function(wrap_pyfunction!(_matrix_trace_normalize, m)?)?;
    Ok(())
}
