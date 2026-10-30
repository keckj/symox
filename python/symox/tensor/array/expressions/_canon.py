"""SymPy 1.14 array-expression canonical forms.

The ``_canonicalize``/``doit``/``nest_permutation`` rules of SymPy's
``array_expressions`` module run in Rust (``crates/symox-core/src/arrayexpr.rs``,
binding ``_arrayexpr_canon``). This module only encodes the expression tree
(tensor products, sums, permutations, contractions, diagonals and zero arrays
become nodes, every other object is an opaque leaf with its shape) and
rebuilds the objects from the result.
"""
import symox as S
from . import (ArrayAdd, ArrayContraction, ArrayElementwiseApplyFunc, Reshape, ArrayDiagonal, ArrayTensorProduct, PermuteDims,
               ZeroArray, _shape, _sort_key)
from ...._kind import matrix_tag
from .... import _core

_STRUCTURAL = (ArrayTensorProduct, ArrayAdd, PermuteDims, ArrayContraction, ArrayDiagonal, ZeroArray)


def _perm_form(permutation, rank=None):
    if hasattr(permutation, 'array_form'):
        return list(permutation.array_form)
    form = list(permutation)
    if form and isinstance(form[0], (list, tuple)):
        from ....combinatorics import Permutation
        return list(Permutation(form).array_form)
    return [int(i) for i in form]


class _Codec:
    def __init__(self, deep=False):
        self.deep = deep
        self.leaves = []
        self.info = []
        self.dims = {}
        self.dimlist = []

    def dim(self, d):
        if not isinstance(d, _core.Expr):
            d = S.sympify(d)
        k = self.dims.get(d)
        if k is None:
            k = self.dims[d] = len(self.dimlist)
            self.dimlist.append(d)
        return k

    def enc(self, e):
        t = type(e)
        if t is ArrayTensorProduct:
            return ('T', [self.enc(a) for a in e.args])
        if t is ArrayAdd:
            return ('A', [self.enc(a) for a in e.args])
        if t is PermuteDims:
            return ('P', self.enc(e.expr), list(e.permutation.array_form))
        if t is ArrayContraction:
            return ('C', self.enc(e.expr), [list(g) for g in e.groups])
        if t is ArrayDiagonal:
            return ('D', self.enc(e.expr), [list(g) for g in e.groups])
        if t is ZeroArray:
            return ('Z', [self.dim(n) for n in e.shape])
        if self.deep and hasattr(e, 'doit'):
            e = e.doit()
            if isinstance(e, _STRUCTURAL):
                return self.enc(e)
        if isinstance(e, (int, float, complex)):
            e = S.sympify(e)
        shape = _shape(e)
        self.leaves.append(e)
        self.info.append((len(shape), [self.dim(n) for n in shape], matrix_tag(e) == 'ZeroMatrix'))
        return ('L', len(self.leaves) - 1)

    def dec(self, n):
        tag = n[0]
        if tag == 'L':
            return self.leaves[n[1]]
        if tag == 'Z':
            return ZeroArray(*[self.dimlist[k] for k in n[1]])
        if tag == 'T':
            return ArrayTensorProduct(*[self.dec(a) for a in n[1]])
        if tag == 'A':
            return ArrayAdd(*[self.dec(a) for a in n[1]])
        if tag == 'P':
            return PermuteDims(self.dec(n[1]), n[2])
        if tag == 'C':
            return ArrayContraction(self.dec(n[1]), *[tuple(g) for g in n[2]])
        if tag == 'R':
            return Reshape(self.dec(n[1]), [self.dimlist[k] for k in n[2]])
        if tag == 'E':
            d = S.Dummy('d')
            return ArrayElementwiseApplyFunc(S.Lambda(d, d**self.leaves[n[2]]), self.dec(n[1]))
        return ArrayDiagonal(self.dec(n[1]), *[tuple(g) for g in n[2]])

    def order(self, nodes):
        keys = [_sort_key(self.dec(n)) for n in nodes]
        return sorted(range(len(nodes)), key=keys.__getitem__)

    def run(self, kind, tree, allow_trivial=False):
        ones = [d == 1 for d in self.dimlist]
        return self.dec(_core._arrayexpr_canon(kind, tree, self.info, ones, self.order, allow_trivial))


def _array_tensor_product(*args, **kwargs):
    c = _Codec()
    return c.run('canon', ('T', [c.enc(a) for a in args]))


def _array_add(*args, **kwargs):
    raw = ArrayAdd(*args)
    c = _Codec()
    return c.run('canon', c.enc(raw))


def _permute_dims(expr, permutation=None, index_order_old=None, index_order_new=None, **kwargs):
    if permutation is None:
        if index_order_old is None or index_order_new is None:
            raise ValueError("Permutation not defined")
        permutation = [list(index_order_old).index(i) for i in index_order_new]
    elif index_order_new is not None:
        raise ValueError("index_order_new cannot be defined with permutation")
    elif index_order_old is not None:
        raise ValueError("index_order_old cannot be defined with permutation")
    c = _Codec()
    return c.run('canon', ('P', c.enc(expr), _perm_form(permutation)))


def _array_contraction(expr, *contraction_indices, **kwargs):
    c = _Codec()
    return c.run('canon', ('C', c.enc(expr), [[int(j) for j in g] for g in contraction_indices]))


def _array_diagonal(expr, *diagonal_indices, **kwargs):
    c = _Codec()
    return c.run('canon', ('D', c.enc(expr), [[int(j) for j in g] for g in diagonal_indices]),
                 bool(kwargs.get('allow_trivial_diags', False)))


def _doit(expr, deep=True):
    c = _Codec(deep)
    return c.run('doit' if deep else 'canon', c.enc(expr))


def nest_permutation(expr):
    if not isinstance(expr, PermuteDims):
        return expr
    c = _Codec()
    return c.run('nest', c.enc(expr))


def _order_with(nodes, leaves, dims):
    c = _Codec()
    c.leaves, c.dimlist = list(leaves), list(dims)
    return c.order(nodes)


def convert_matrix_to_array(expr):
    """SymPy's ``convert_matrix_to_array``, computed in Rust
    (``crates/symox-core/src/arrayexpr_matrix.rs``)."""
    out = _core._arrayexpr_from_matrix(expr, _order_with)
    if out is None:
        return expr
    tree, leaves, dims = out
    c = _Codec()
    c.leaves, c.dimlist = list(leaves), list(dims)
    return c.dec(tree)
