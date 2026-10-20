"""Printer dispatch classes of SymPy 1.14 objects that are not native nodes.

Native nodes are instances of nominal classes whose MRO names SymPy's
classes. Matrices and matrix expressions are not: their SymPy MRO names,
taken from the reference SymPy, are listed here for ``Printer._print`` and
``precedence``.
"""

MATRIX_CLASSES = {
 'Adjoint': ('Adjoint', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'BlockDiagMatrix': ('BlockDiagMatrix', 'BlockMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'BlockMatrix': ('BlockMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'CompanionMatrix': ('CompanionMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'DFT': ('DFT', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'DenseMatrix': ('DenseMatrix', 'RepMatrix', 'MatrixBase', 'Printable'),
 'Determinant': ('Determinant', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'DiagMatrix': ('DiagMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'DiagonalMatrix': ('DiagonalMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'DiagonalOf': ('DiagonalOf', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'DotProduct': ('DotProduct', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'EigenValues': ('EigenValues', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'EigenVectors': ('EigenVectors', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'ElementwiseApplyFunction': ('ElementwiseApplyFunction', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'Factorization': ('Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'FunctionMatrix': ('FunctionMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'GenericIdentity': ('GenericIdentity', 'Identity', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'GenericZeroMatrix': ('GenericZeroMatrix', 'ZeroMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'HadamardPower': ('HadamardPower', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'HadamardProduct': ('HadamardProduct', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'IDFT': ('IDFT', 'DFT', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'Identity': ('Identity', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'ImmutableDenseMatrix': ('ImmutableDenseMatrix',
                          'DenseMatrix',
                          'ImmutableRepMatrix',
                          'RepMatrix',
                          'MatrixBase',
                          'MatrixExpr',
                          'Expr',
                          'Basic',
                          'Printable',
                          'EvalfMixin'),
 'ImmutableRepMatrix': ('ImmutableRepMatrix',
                        'RepMatrix',
                        'MatrixBase',
                        'MatrixExpr',
                        'Expr',
                        'Basic',
                        'Printable',
                        'EvalfMixin'),
 'ImmutableSparseMatrix': ('ImmutableSparseMatrix',
                           'SparseRepMatrix',
                           'ImmutableRepMatrix',
                           'RepMatrix',
                           'MatrixBase',
                           'MatrixExpr',
                           'Expr',
                           'Basic',
                           'Printable',
                           'EvalfMixin'),
 'Inverse': ('Inverse', 'MatPow', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'KroneckerProduct': ('KroneckerProduct', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'LofCholesky': ('LofCholesky', 'LofLU', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'LofLU': ('LofLU', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'MatAdd': ('MatAdd', 'MatrixExpr', 'Add', 'Expr', 'AssocOp', 'Basic', 'Printable', 'EvalfMixin'),
 'MatMul': ('MatMul', 'MatrixExpr', 'Mul', 'Expr', 'AssocOp', 'Basic', 'Printable', 'EvalfMixin'),
 'MatPow': ('MatPow', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'MatrixElement': ('MatrixElement', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'MatrixPermute': ('MatrixPermute', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'MatrixSlice': ('MatrixSlice', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'MatrixSymbol': ('MatrixSymbol', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'MutableDenseMatrix': ('MutableDenseMatrix',
                        'DenseMatrix',
                        'MutableRepMatrix',
                        'RepMatrix',
                        'MatrixBase',
                        'Printable'),
 'MutableRepMatrix': ('MutableRepMatrix', 'RepMatrix', 'MatrixBase', 'Printable'),
 'MutableSparseMatrix': ('MutableSparseMatrix',
                         'SparseRepMatrix',
                         'MutableRepMatrix',
                         'RepMatrix',
                         'MatrixBase',
                         'Printable'),
 'OneMatrix': ('OneMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'Permanent': ('Permanent', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'PermutationMatrix': ('PermutationMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'QofQR': ('QofQR', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'RepMatrix': ('RepMatrix', 'MatrixBase', 'Printable'),
 'RofQR': ('RofQR', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'SofSVD': ('SofSVD', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'SparseRepMatrix': ('SparseRepMatrix', 'RepMatrix', 'MatrixBase', 'Printable'),
 'Trace': ('Trace', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'Transpose': ('Transpose', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'UofCholesky': ('UofCholesky', 'UofLU', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'UofLU': ('UofLU', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'UofSVD': ('UofSVD', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'VofSVD': ('VofSVD', 'Factorization', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'),
 'ZeroMatrix': ('ZeroMatrix', 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin')}


# Symox implementation classes standing for a SymPy class of another name.
ALIASES = {'_CallableExpr': 'Lambda'}

# Symox implementation classes that have no SymPy counterpart.
INTERNAL = frozenset(('HostTree', '_HostAtom', '_HostFunction', '_IndexedArithmetic', 'object'))

import weakref

# Keyed by identity: nominal classes compare equal to their facades.
_names = {}


def _splice(names, more):
    for name in more:
        if name not in names:
            names.append(name)


def _type_names(cls):
    entry = _names.get(id(cls))
    if entry is not None and entry[0]() is cls:
        return entry[1]
    from .._node_types import NodeType, is_symox_class, nominal_anchor, node_class
    names = []
    for base in cls.__mro__:
        name = base.__name__
        if base is object:
            continue
        if isinstance(base, NodeType) or not is_symox_class(base):
            _splice(names, (name,))
            continue
        alias = ALIASES.get(name)
        if alias is None and is_symox_class(base):
            # A symox base standing for a SymPy class of another name
            # (``_HostFunction`` for ``Function``).
            alias = base.__dict__.get('_symox_nominal')
        anchor = node_class(alias) if alias is not None else nominal_anchor(base)
        if anchor is not None:
            _splice(names, [b.__name__ for b in anchor.__mro__ if isinstance(b, NodeType)])
        elif name in MATRIX_CLASSES:
            _splice(names, MATRIX_CLASSES[name])
        elif name not in INTERNAL and not name.startswith('_'):
            _splice(names, (name,))
    names = tuple(names)
    try:
        _names[id(cls)] = (weakref.ref(cls), names)
    except TypeError:
        pass
    return names


_UNDEFINED = ('UndefinedFunction', 'FunctionClass', 'type')


def _class_type_names(cls):
    """SymPy's metaclass MRO names for a class: ``FunctionClass`` for the
    applied-function classes, ``type`` otherwise."""
    from .._function_heads import Application
    try:
        if issubclass(cls, Application):
            return ('FunctionClass', 'type')
    except TypeError:
        pass
    return tuple(base.__name__ for base in type(cls).__mro__
                 if base is not object and not base.__name__.startswith('_') and base.__name__ != 'ABCMeta'
                 and base.__name__ != 'BuiltinHeadType')


def class_names(expr):
    """Names of the SymPy classes of ``expr``, in SymPy's MRO order."""
    from .. import _core
    from .._kind import matrix_tag
    from ..core.function import UndefinedFunction
    cls = type(expr)
    if isinstance(expr, UndefinedFunction):
        return _UNDEFINED
    if isinstance(expr, type):
        return _class_type_names(expr)
    if cls is _core._MatrixExpr:
        kind = matrix_tag(expr)
        return MATRIX_CLASSES.get(kind, (kind, 'MatrixExpr', 'Expr', 'Basic', 'Printable', 'EvalfMixin'))
    if cls is _core.Matrix:
        return MATRIX_CLASSES['ImmutableDenseMatrix' if expr._immutable else 'MutableDenseMatrix']
    if cls is _core.Expr:
        # A raw native node whose SymPy class is named by ``class_name``
        # (``MatrixElement``).
        name = getattr(expr, 'class_name', None)
        if name in MATRIX_CLASSES:
            return MATRIX_CLASSES[name]
    return _type_names(cls)
