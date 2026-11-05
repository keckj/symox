"""SymPy 1.14 ``array_derive``, ``convert_indexed_to_array`` and
``convert_array_to_indexed``.

The three algorithms run in Rust (``crates/symox-core/src/arrayexpr_derive.rs``
and ``arrayexpr_index.rs``, bindings ``_arrayexpr_derive``,
``_arrayexpr_from_indexed`` and ``_arrayexpr_to_indexed``). This module only
encodes the expressions as tuples (the dispatch of the reference picks the
node kind of each object, the objects the reference builds from Python
methods are made here as leaves) and rebuilds the objects from the result.
"""
import symox as S
from . import (ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayElement, ArrayElementwiseApplyFunc, ArraySymbol,
               ArrayTensorProduct, OneArray, PermuteDims, Reshape, ZeroArray, _ArrayExpr, _ConstantArray,
               _sort_key)
from ._canon import _Codec
from ...._kind import matrix_tag
from .... import _core

_NO = object()


def _shape_attr(e):
    """``e.shape``, ``_NO`` without the attribute (``hasattr`` semantics)."""
    try:
        shape = getattr(e, 'shape', _NO)
    except Exception:
        return _NO
    if shape is _NO or shape is None:
        return shape
    try:
        return tuple(S.sympify(n) for n in shape)
    except TypeError:
        return _NO


def _no_shape_error(e):
    return AttributeError("'%s' object has no attribute 'shape'" % type(e).__name__)


class _NCodec(_Codec):
    """The array codec of ``_canon`` with element-wise functions (tag 'E'),
    reshapes (tag 'R') and leaves of unknown shape."""

    def __init__(self):
        super().__init__()
        self.has_shape = []
        self.unknown = []

    def leaf(self, e, info=None):
        if isinstance(e, (int, float, complex)):
            e = S.sympify(e)
        if info is None:
            shape = _shape_attr(e)
            if shape is _NO:
                info, has = (0, [], False), False
            elif shape is None:
                info, has = (0, None, False), True
            else:
                info, has = (len(shape), [self.dim(n) for n in shape], matrix_tag(e) == 'ZeroMatrix'), True
        else:
            has = False
        self.leaves.append(e)
        self.info.append(info)
        self.has_shape.append(has)
        self.unknown.append(False)
        return len(self.leaves) - 1

    def enc(self, e):
        t = type(e)
        if t in (ArrayTensorProduct, ArrayAdd, PermuteDims, ArrayContraction, ArrayDiagonal, ZeroArray):
            return super().enc(e)
        if t is Reshape:
            return ('R', self.enc(e.expr), [self.dim(n) for n in e.shape])
        if t is ArrayElementwiseApplyFunc:
            return ('E', self.enc(e.expr), self.leaf(e.function, (0, [], False)))
        return ('L', self.leaf(e))

    def dec(self, n):
        if n[0] == 'E':
            return ArrayElementwiseApplyFunc(self.leaves[n[2]], self.dec(n[1]))
        return super().dec(n)


def _ones(c):
    return [d == 1 for d in c.dimlist]


# ------------------------------------------------------------ array_derive

def _is_explicit_immutable(e):
    return isinstance(e, S.ImmutableMatrix) and not isinstance(e, S.MatrixExpr)


def _get_rank(e):
    """``get_rank`` (a leaf property computed for the encoding)."""
    if isinstance(e, S.MatrixExpr) or matrix_tag(e) is not None or isinstance(e, S.MatrixElement):
        return 2
    shape = _shape_attr(e)
    if shape is _NO:
        from ....codegen.ast import Indexed
        if isinstance(e, Indexed):
            return len(e.indices)
        return 0
    if shape is None:
        return -1
    return len(shape)


class _Deriver:
    def __init__(self, x):
        self.c = _NCodec()
        self.x = x
        self.errors = []

    def fail(self, exc):
        self.errors.append(exc)
        return ('F', len(self.errors) - 1)

    def result(self, fn):
        try:
            return (True, fn())
        except Exception as exc:
            self.errors.append(exc)
            return (False, len(self.errors) - 1)

    def child(self, e):
        # The structural constructors of the reference sympify their arguments.
        if isinstance(e, (int, float, complex)):
            e = S.sympify(e)
        return self.node(e)

    def node(self, e):
        c, t = self.c, type(e)
        if t is ArrayTensorProduct:
            args = [S.sympify(a) if isinstance(a, (int, float, complex)) else a for a in e.args]
            return ('T', [self.node(a) for a in args], [c.enc(a) for a in args])
        if t is ArrayAdd:
            return ('A', [self.child(a) for a in e.args])
        if t is PermuteDims:
            return ('P', self.child(e.expr), list(e.permutation.array_form))
        if t is ArrayContraction:
            return ('C', self.child(e.expr), [list(g) for g in e.contraction_indices])
        if t is ArrayDiagonal:
            return ('D', self.child(e.expr), [list(g) for g in e.diagonal_indices])
        if t is Reshape:
            return ('R', self.child(e.expr), [c.dim(n) for n in e.shape])
        if t is ArrayElementwiseApplyFunc:
            sub = e.expr
            fdiff = self.result(e._get_function_fdiff)
            applied = (True, c.leaf(ArrayElementwiseApplyFunc(fdiff[1], sub))) if fdiff[0] else fdiff
            shape = _shape_attr(sub)
            if shape is _NO or shape is None:
                self.errors.append(_no_shape_error(sub) if shape is _NO else TypeError("object of type 'NoneType' has no len()"))
                rank = (False, len(self.errors) - 1)
            else:
                rank = (True, len(shape))
            return ('EWA', applied, self.child(sub), rank)
        if t is ArraySymbol:
            if e == self.x:
                return ('SA', [c.leaf(S.Identity(n)) for n in e.shape])
            return ('ZC', [c.dim(n) for n in e.shape])
        tag = matrix_tag(e)
        if tag == 'MatrixSymbol':
            m, n = e.shape
            if e == self.x:
                return ('SM', c.leaf(S.Identity(m)), c.leaf(S.Identity(n)))
            return ('ZC', [c.dim(m), c.dim(n)])
        if tag in ('Identity', 'OneMatrix'):
            return ('ZC', [c.dim(n) for n in e.shape])
        if tag == 'Transpose':
            return ('TR', self.node(e.arg))
        if tag == 'Inverse':
            mat = e.I
            # MatrixExpr.__neg__ of the reference: MatMul(S.NegativeOne, expr).doit()
            return ('INV', self.node(mat), c.leaf(S.MatMul(S.Integer(-1), e).doit()), c.leaf(e))
        if tag == 'ElementwiseApplyFunction':
            assert_err = None
            if _get_rank(self.x) != 2:
                self.errors.append(AssertionError())
                assert_err = len(self.errors) - 1
            fdiff = self.result(e._get_function_fdiff)
            applied = (True, c.leaf(type(e)(fdiff[1], e.expr))) if fdiff[0] else fdiff
            return ('EWM', assert_err, applied, self.node(e.expr))
        if tag == 'HadamardProduct':
            return self.fail(NotImplementedError())
        if tag is not None or isinstance(e, S.MatrixExpr) or _is_explicit_immutable(e):
            from ._canon import convert_matrix_to_array
            cg = convert_matrix_to_array(e)
            if cg is e or cg == e:
                return self.fail(RecursionError('maximum recursion depth exceeded in comparison'))
            return ('CV', self.node(cg))
        if isinstance(e, (_ArrayExpr, ArrayElement, S.Expr)) and not isinstance(e, (int, float, complex)):
            return ('S',)
        return self.fail(NotImplementedError(f"not implemented for type {type(e)}"))

    def x_info(self):
        x, c = self.x, self.c
        shape = _shape_attr(x)
        if shape is _NO:
            self.errors.append(_no_shape_error(x))
            xs = gs = (False, len(self.errors) - 1)
        elif shape is None:
            self.errors.append(TypeError("ZeroArray() argument after * must be an iterable, not NoneType"))
            xs = (False, len(self.errors) - 1)
            self.errors.append(TypeError("object of type 'NoneType' has no len()"))
            gs = (False, len(self.errors) - 1)
        else:
            xs = gs = (True, [c.dim(n) for n in shape])
        if shape is _NO:
            gs = (True, [])
        return (xs, gs, (True, _get_rank(x)))


def array_derive(expr, x):
    """SymPy's ``array_derive(expr, x)``, computed in Rust
    (``crates/symox-core/src/arrayexpr_derive.rs``)."""
    d = _Deriver(x)
    tree = d.node(expr)
    xinfo = d.x_info()
    c = d.c
    szero = c.leaf(S.Integer(0))
    out = _core._arrayexpr_derive(tree, c.leaves, c.info, c.has_shape, _ones(c), szero, xinfo, d.errors, c.order)
    return c.dec(out)


# ------------------------------------------------- convert_indexed_to_array

class _IndexedEncoder:
    def __init__(self):
        self.c = _NCodec()
        self.atoms = []
        self.atom_ids = {}

    def atom(self, a):
        a = S.sympify(a)
        k = self.atom_ids.get(a)
        if k is None:
            k = self.atom_ids[a] = len(self.atoms)
            self.atoms.append(a)
        return k

    def parent(self, p, n_indices):
        c = self.c
        k = c.leaf(p)
        if c.info[k][1] is None:
            c.info[k] = (n_indices, None, False)
            c.unknown[k] = True
        return k

    def element(self, parent, indices):
        indices = list(indices)
        return ('E', self.parent(parent, len(indices)), [self.atom(i) for i in indices])

    def node(self, e):
        from ....codegen.ast import Indexed, IndexedBase, _indexed_reference
        c = self.c
        if isinstance(e, S.Sum):
            limits = []
            for lim in e.limits:
                if len(lim) != 3:
                    limits.append((False, len(lim)))
                    continue
                ind, start, end = lim
                limits.append((True, (self.atom(ind), bool(start != 0), c.dim(end + 1))))
            return ('S', self.node(e.function), [self.atom(v) for v in e.variables], limits)
        if isinstance(e, S.Mul):
            return ('M', [self.node(a) for a in e.args])
        if isinstance(e, S.MatrixElement):
            return self.element(e.args[0], e.args[1:])
        if isinstance(e, ArrayElement):
            return self.element(e.name, e.indices)
        if '_symox_element_parent' in type(e).__dict__:
            # Native scalar form of an ArrayElement: its class records the parent.
            return self.element(type(e)._symox_element_parent, e.args[1:])
        ref = _indexed_reference(e) if isinstance(e, S.Symbol) else None
        if ref is not None:
            e = ref
        if isinstance(e, Indexed):
            return self.element(e.base, e.indices)
        if isinstance(e, IndexedBase):
            return ('B',)
        if isinstance(e, S.KroneckerDelta):
            i, j = e.args
            return ('K', c.leaf(e), self.atom(i), self.atom(j))
        if isinstance(e, S.Add):
            return ('A', [self.node(a) for a in e.args])
        if isinstance(e, S.Pow):
            base = self.node(e.base)
            if isinstance(e.exp, S.Integer):
                return ('P', base, (True, int(e.exp)))
            return ('P', base, (False, c.leaf(e)))
        if isinstance(e, S.Function):
            arg = self.node(e.args[0]) if e.args else None
            return ('F', arg, c.leaf(type(e), (0, [], False)))
        return ('L', c.leaf(e))

    def index(self, i):
        return frozenset(self.atoms[a] for a in i) if isinstance(i, list) else self.atoms[i]

    def run(self, expr, first):
        c = self.c
        tree = self.node(expr)
        given = list(c.leaves)

        def extend(extra):
            c.leaves = given + [OneArray(*[c.dimlist[d] for d in dims]) if dims else S.Integer(1) for dims in extra]

        def order(nodes, extra):
            extend(extra)
            return c.order(nodes)

        first_ids = None if first is None else [self.atom(i) for i in first]
        int_atoms = [isinstance(a, S.Integer) for a in self.atoms]
        out = _core._arrayexpr_from_indexed(tree, c.info, c.unknown, c.has_shape, _ones(c), int_atoms, first_ids, order)
        if out[0] == 'error':
            kind, payload = out[1], out[2]
            if kind == 'notfound':
                a, sub = payload
                raise IndexError("%s not found in %s" % (self.atoms[a], tuple(self.index(i) for i in sub)))
            if kind == 'mismatch':
                raise ValueError("summation index and array dimension mismatch: %s" % self.atoms[payload])
            raise _no_shape_error(c.leaves[payload])
        tag, node, indices, extra = out
        extend(extra)
        return tag, c.dec(node), tuple(self.index(i) for i in indices)


def _convert_indexed_to_array(expr):
    """SymPy's ``_convert_indexed_to_array``: the array form and its free
    indices (Rust, ``crates/symox-core/src/arrayexpr_index.rs``)."""
    _, result, indices = _IndexedEncoder().run(expr, None)
    return result, indices


def convert_indexed_to_array(expr, first_indices=None):
    """SymPy's ``convert_indexed_to_array``, computed in Rust
    (``crates/symox-core/src/arrayexpr_index.rs``)."""
    first = list(first_indices) if first_indices else None
    tag, result, indices = _IndexedEncoder().run(expr, first)
    if tag == 'element':
        # ArrayElement(result, indices); the permutation that follows for
        # first_indices is the empty identity.
        return ArrayElement(result, indices)
    return result


# ------------------------------------------------- convert_array_to_indexed

def _get(obj, indices):
    """``_get_array_element_or_slice`` and ``_ArrayExpr.__getitem__``."""
    if isinstance(obj, (ArraySymbol, _ConstantArray)):
        ArrayElement._check_shape(obj, indices)
        if isinstance(obj, ArraySymbol):
            return ArrayElement(obj, indices)
        return S.Integer(obj.value)
    return ArrayElement(obj, indices)


def convert_array_to_indexed(expr, indices):
    """SymPy's ``convert_array_to_indexed``: the index plan is computed in
    Rust (``crates/symox-core/src/arrayexpr_index.rs``)."""
    c = _NCodec()
    tree = c.enc(expr)
    indices = list(indices)
    plan = _core._arrayexpr_to_indexed(tree, c.info, _ones(c), len(indices))
    dummies = {}

    def ie(x):
        if x is None:
            return None
        tag = x[0]
        if tag == 'g':
            return indices[x[1]]
        if tag == 'd':
            if x[1] not in dummies:
                dummies[x[1]] = S.Dummy(f"d{x[1]}")
            return dummies[x[1]]
        if tag == 'D':
            return c.dimlist[x[1]]
        if tag == 'i':
            return S.Integer(x[1])
        if tag == '+':
            return S.Add(*[ie(y) for y in x[1]])
        if tag == '*':
            return S.Mul(*[ie(y) for y in x[1]])
        if tag == '//':
            return ie(x[1]) // ie(x[2])
        return ie(x[1]) % ie(x[2])

    def build(p):
        tag = p[0]
        if tag == 'M':
            return S.Mul(*[build(q) for q in p[1]])
        if tag == 'A':
            return S.Add(*[build(q) for q in p[1]])
        if tag == 'S':
            body = build(p[1])
            return S.Sum(body, *[(ie(('d', k)), 0, c.dimlist[dim] - 1) for k, dim in p[2]])
        if tag == 'G':
            return _get(c.leaves[p[1]], tuple(ie(x) for x in p[2]))
        if tag == 'Z':
            return _get(ZeroArray(*[c.dimlist[d] for d in p[1]]), tuple(ie(x) for x in p[2]))
        return c.leaves[p[1]](build(p[2]))

    return build(plan)
