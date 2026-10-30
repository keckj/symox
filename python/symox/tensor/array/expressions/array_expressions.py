"""``sympy.tensor.array.expressions.array_expressions`` import path.

Re-exports the array-expression classes of ``symox.tensor.array.expressions``
and SymPy's canonicalizing constructors, whose ``_canonicalize`` rules run in
Rust (``_canon``), and SymPy's contraction editor ``_EditArrayContraction``:
the editor's state (the ``_ArgE`` list, the contraction count and the
tracked permutation) lives here as in SymPy, so that callers can edit it,
and every method computes in Rust (``crates/symox-core/src/arrayexpr_a2m.rs``).
An element the native port does not model runs the stock SymPy method.
"""
from . import (ArrayAdd, ArrayContraction, ArrayDiagonal, ArrayElement,
               ArrayElementwiseApplyFunc, ArraySymbol, ArrayTensorProduct,
               OneArray, PermuteDims, Reshape, ZeroArray, _ArrayExpr, _shape)
from ._canon import (_array_add, _array_contraction, _array_diagonal,
                     _array_tensor_product, _permute_dims, nest_permutation)
from . import _a2m


def get_rank(expr):
    # SymPy's get_rank: a matrix element counts as rank 2.
    if getattr(getattr(expr, 'func', None), '__name__', None) == 'MatrixElement':
        return 2
    return len(_shape(expr))


class _ArgE:
    """
    The ``_ArgE`` object contains references to the array expression
    (``.element``) and a list containing the information about index
    contractions (``.indices``).

    Index contractions are numbered and contracted indices show the number of
    the contraction. Uncontracted indices have ``None`` value.
    """

    def __init__(self, element, indices=None):
        self.element = element
        if indices is None:
            self.indices = [None for i in range(get_rank(element))]
        else:
            self.indices = indices

    def __str__(self):
        return "_ArgE(%s, %s)" % (self.element, self.indices)

    __repr__ = __str__


class _IndPos:
    """
    Index position, requiring two integers in the constructor:

    - arg: the position of the argument in the tensor product,
    - rel: the relative position of the index inside the argument.
    """
    def __init__(self, arg, rel):
        self.arg = arg
        self.rel = rel

    def __str__(self):
        return "_IndPos(%i, %i)" % (self.arg, self.rel)

    __repr__ = __str__

    def __iter__(self):
        yield from [self.arg, self.rel]


def _stock_editor_class():
    from ...._arrayexpr import array_expressions as V
    return V


class _EditArrayContraction:
    """
    Utility class to help manipulate array contraction objects.

    This class takes as input an ``ArrayContraction`` object and turns it into
    an editable object.

    The field ``args_with_ind`` of this class is a list of ``_ArgE`` objects
    which can be used to easily edit the contraction structure of the
    expression.

    Once editing is finished, the ``ArrayContraction`` object may be recreated
    by calling the ``.to_array_contraction()`` method.
    """

    def __init__(self, base_array):
        if not isinstance(base_array, (ArrayContraction, ArrayDiagonal, ArrayTensorProduct)):
            raise NotImplementedError()
        self._track_permutation = None
        self.args_with_ind = []
        self.number_of_contraction_indices = 0
        try:
            c = _a2m._A2mCodec()
            tree = c.enc(base_array)
            (args, n, _), _ = c.call('ed_new', [tree], ((), 0, None, None))
        except (_a2m._Decline, _a2m.Unsupported):
            from ._stock import to_stock
            stock = _stock_editor_class()._EditArrayContraction(to_stock(base_array))
            self._from_stock(stock, {})
            return
        self.args_with_ind = [_ArgE(c.dec(t), list(ind)) for _, t, ind in args]
        self.number_of_contraction_indices = n

    # -- state transfer -------------------------------------------------

    def _state(self, param=None):
        return ([list(a.indices) for a in self.args_with_ind], self.number_of_contraction_indices,
                None if self._track_permutation is None else [list(t) for t in self._track_permutation], param)

    def _apply(self, c, state):
        args, n, track = state
        old = list(self.args_with_ind)
        new = []
        for orig, tree, ind in args:
            a = old[orig] if orig is not None else _ArgE.__new__(_ArgE)
            if tree is not None:
                a.element = c.dec(tree)
            a.indices = list(ind)
            new.append(a)
        self.args_with_ind[:] = new
        self.number_of_contraction_indices = n
        self._track_permutation = None if track is None else [list(t) for t in track]

    def _to_stock(self):
        from ._stock import to_stock
        V = _stock_editor_class()
        ed = V._EditArrayContraction.__new__(V._EditArrayContraction)
        back = {}
        ed.args_with_ind = []
        for a in self.args_with_ind:
            s = V._ArgE(to_stock(a.element), list(a.indices))
            back[id(s)] = a
            ed.args_with_ind.append(s)
        ed.number_of_contraction_indices = self.number_of_contraction_indices
        ed._track_permutation = None if self._track_permutation is None else [list(t) for t in self._track_permutation]
        return ed, back

    def _from_stock(self, ed, back):
        from ._stock import from_stock
        new = []
        for s in ed.args_with_ind:
            a = back.get(id(s))
            if a is None:
                a = _ArgE.__new__(_ArgE)
            a.element = from_stock(s.element)
            a.indices = list(s.indices)
            new.append(a)
        self.args_with_ind[:] = new
        self.number_of_contraction_indices = ed.number_of_contraction_indices
        self._track_permutation = ed._track_permutation

    def _run(self, method, param=None, elements=False):
        c = _a2m._A2mCodec()
        try:
            trees = [c.enc(a.element) if elements else None for a in self.args_with_ind]
            state, value = c.call('ed_' + method, trees, self._state(param))
        except (_a2m._Decline, _a2m.Unsupported):
            if not elements:
                raise
            ed, back = self._to_stock()
            name = {'merge': 'merge_scalars', 'to_array': 'to_array_contraction'}[method]
            out = getattr(ed, name)()
            self._from_stock(ed, back)
            from ._stock import from_stock
            return ('stock', from_stock(out))
        self._apply(c, state)
        if method == 'to_array':
            return ('tree', c.dec(value))
        return value

    def _position(self, arg):
        for i, a in enumerate(self.args_with_ind):
            if a == arg:
                return i
        return None

    # -- SymPy's methods ------------------------------------------------

    def insert_after(self, arg, new_arg):
        pos = self.args_with_ind.index(arg)
        self.args_with_ind.insert(pos + 1, new_arg)

    def get_new_contraction_index(self):
        return self._run('new_index')

    def refresh_indices(self):
        self._run('refresh')

    def merge_scalars(self):
        self._run('merge', elements=True)

    def to_array_contraction(self):
        return self._run('to_array', elements=True)[1]

    def get_contraction_indices(self):
        return [list(g) for g in self._run('contraction_indices')]

    def get_mapping_for_index(self, ind):
        return [_IndPos(i, j) for i, j in self._run('mapping', ind)]

    def get_contraction_indices_to_ind_rel_pos(self):
        return [[_IndPos(i, j) for i, j in g] for g in self._run('ind_rel_pos')]

    def count_args_with_index(self, index):
        """
        Count the number of arguments that have the given index.
        """
        return self._run('count', index)

    def get_args_with_index(self, index):
        """
        Get a list of arguments having the given index.
        """
        return [self.args_with_ind[i] for i in self._run('args_with', index)]

    @property
    def number_of_diagonal_indices(self):
        return self._run('n_diag')

    def track_permutation_start(self):
        self._run('track_start')

    def track_permutation_merge(self, destination, from_element):
        index_destination = self.args_with_ind.index(destination)
        index_element = self.args_with_ind.index(from_element)
        self._run('track_merge', (index_destination, index_element))

    def get_absolute_free_range(self, arg):
        """
        Return the range of the free indices of the arg as absolute positions
        among all free indices.
        """
        return self._run('free_range', self._position(arg))

    def get_absolute_range(self, arg):
        """
        Return the absolute range of indices for arg, disregarding dummy
        indices.
        """
        return self._run('range', self._position(arg))


__all__ = ['ArrayAdd', 'ArrayContraction', 'ArrayDiagonal', 'ArrayElement',
           'ArrayElementwiseApplyFunc', 'ArraySymbol', 'ArrayTensorProduct',
           'OneArray', 'PermuteDims', 'Reshape', 'ZeroArray', 'nest_permutation']
