"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Atom Basic')


def as_Basic(expr):
    from symox import sympify, SympifyError
    try:
        return sympify(expr, strict=True)
    except SympifyError:
        raise TypeError('Argument must be a Basic object, not `%s`' % type(expr).__name__) from None
__all__.append('as_Basic')


from functools import cmp_to_key as _cmp_to_key
from symox import _core
from symox._basic import Basic as _Basic
from symox._compat import _ExprType

# SymPy's class precedence for Basic.compare (sympy/core/basic.py, 1.14). The
# native argument order uses the same default table; native expressions keep
# that order even if downstream code inserts names for its own classes here.
ordering_of_classes = [
    # singleton numbers
    'Zero', 'One', 'Half', 'Infinity', 'NaN', 'NegativeOne', 'NegativeInfinity',
    # numbers
    'Integer', 'Rational', 'Float',
    # singleton symbols
    'Exp1', 'Pi', 'ImaginaryUnit',
    # symbols
    'Symbol', 'Wild',
    # arithmetic operations
    'Pow', 'Mul', 'Add',
    # function values
    'Derivative', 'Integral',
    # defined singleton functions
    'Abs', 'Sign', 'Sqrt',
    'Floor', 'Ceiling',
    'Re', 'Im', 'Arg',
    'Conjugate',
    'Exp', 'Log',
    'Sin', 'Cos', 'Tan', 'Cot', 'ASin', 'ACos', 'ATan', 'ACot',
    'Sinh', 'Cosh', 'Tanh', 'Coth', 'ASinh', 'ACosh', 'ATanh', 'ACoth',
    'RisingFactorial', 'FallingFactorial',
    'factorial', 'binomial',
    'Gamma', 'LowerGamma', 'UpperGamma', 'PolyGamma',
    'Erf',
    # special polynomials
    'Chebyshev', 'Chebyshev2',
    # undefined functions
    'Function', 'WildFunction',
    # anonymous functions
    'Lambda',
    # Landau O symbol
    'Order',
    # relational operations
    'Equality', 'Unequality', 'StrictGreaterThan', 'StrictLessThan',
    'GreaterThan', 'LessThan',
]


def _cmp_name(x, y):
    """Return -1, 0, 1 if the name of class x is before that of class y.

    A string comparison is done if either name does not appear in
    ``ordering_of_classes``. This is the helper for ``Basic.compare``.
    """
    # Public classes of native nodes (Integer, Add, ...) are Basic subclasses
    # in the reference although they emulate the hierarchy by name.
    y_is_basic = isinstance(y, _ExprType) or issubclass(y, _Basic)
    return _cmp_class_names(x.__name__, y.__name__, y_is_basic)


def _cmp_class_names(n1, n2, y_is_basic):
    if n1 == n2:
        return 0

    # If the other object is not a Basic subclass, then we are not equal to it.
    if not y_is_basic:
        return -1

    UNKNOWN = len(ordering_of_classes) + 1
    try:
        i1 = ordering_of_classes.index(n1)
    except ValueError:
        i1 = UNKNOWN
    try:
        i2 = ordering_of_classes.index(n2)
    except ValueError:
        i2 = UNKNOWN
    if i1 == UNKNOWN and i2 == UNKNOWN:
        return (n1 > n2) - (n1 < n2)
    return (i1 > i2) - (i1 < i2)


def _class_name(obj):
    # Native nodes share one Python type; their reference class name is
    # carried by the node.
    if isinstance(obj, _core.Expr):
        return obj.class_name
    return type(obj).__name__


def _basic_compare(self, other):
    """Basic.compare: -1, 0 or 1 in SymPy's canonical argument order."""
    if self is other:
        return 0
    if isinstance(self, _core.Expr) and isinstance(other, _core.Expr):
        return _core._st_compare(self, other)
    c = _cmp_class_names(_class_name(self), _class_name(other),
                         isinstance(other, _Basic))
    if c:
        return c
    st = self._hashable_content() if hasattr(self, '_hashable_content') else self.args
    ot = other._hashable_content() if hasattr(other, '_hashable_content') else other.args
    c = (len(st) > len(ot)) - (len(st) < len(ot))
    if c:
        return c
    for l, r in zip(st, ot):
        if isinstance(l, _Basic):
            c = _basic_compare(l, r)
        elif isinstance(l, frozenset):
            l = _Basic(*l)
            r = _Basic(*r) if isinstance(r, frozenset) else r
            c = _basic_compare(l, r)
        else:
            c = (l > r) - (l < r)
        if c:
            return c
    return 0


if not hasattr(_Basic, 'compare'):
    _Basic.compare = _basic_compare
_args_sortkey = _cmp_to_key(_Basic.compare)
