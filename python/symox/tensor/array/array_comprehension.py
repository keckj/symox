"""``ArrayComprehension`` and ``ArrayComprehensionMap`` (SymPy 1.14).

The classes hold their arguments; limit validation, shape, length, free
symbols, substitution and the expansion of the entries are native
(``symox._core._array_comprehension_*``).
"""
import inspect

from ... import _core
from ._array import ImmutableDenseNDimArray


def _limits(symbols):
    """``_check_limits_validity``: the sympified triples, checked natively."""
    from ... import Tuple, sympify
    if any(len(limit) != 3 or None for limit in symbols):
        raise ValueError('ArrayComprehension requires values lower and upper bound'
                         ' for the expression')
    triples = []
    for var, inf, sup in symbols:
        var = sympify(var, strict=True)
        inf = sympify(inf, strict=True)
        sup = Tuple(*sup) if isinstance(sup, list) else sympify(sup, strict=True)
        triples.append((var, inf, sup))
    return _core._array_comprehension_limits(tuple(triples))


class ArrayComprehension:
    """A comprehension over integer ranges, expanded by ``doit``."""
    __slots__ = ('_args',)

    def __new__(cls, function, *symbols, **assumptions):
        from ... import sympify
        if any(len(limit) != 3 or None for limit in symbols):
            raise ValueError('ArrayComprehension requires values lower and upper bound'
                             ' for the expression')
        function = sympify(function)
        limits = _limits(symbols)
        if assumptions:
            raise TypeError('Basic.__new__() got an unexpected keyword argument %r' % next(iter(assumptions)))
        result = object.__new__(cls)
        result._args = (function, *limits)
        return result

    _mapped = False

    @property
    def args(self): return self._args
    @property
    def function(self): return self._args[0]
    @property
    def limits(self): return self._args[1:]
    @property
    def func(self): return type(self)
    @property
    def variables(self): return [limit[0] for limit in self.limits]
    @property
    def bound_symbols(self): return [limit[0] for limit in self.limits if len(limit) != 1]
    @property
    def shape(self): return _core._array_comprehension_shape(self.limits)
    @property
    def is_shape_numeric(self): return _core._array_comprehension_is_shape_numeric(self.limits)
    @property
    def free_symbols(self): return _core._array_comprehension_free_symbols(self.function, self.limits)
    @property
    def kind(self):
        from ..._kind import UndefinedKind
        return UndefinedKind
    def rank(self): return len(self.limits)
    def __len__(self): return _core._array_comprehension_len(self.limits)

    def _callable(self): return self.function
    def _expand_array(self):
        entries = _core._array_comprehension_entries(self._callable(), self.limits, self._mapped)
        return ImmutableDenseNDimArray(entries, self.shape)
    def doit(self, **hints):
        return self._expand_array() if self.is_shape_numeric else self
    def tolist(self):
        if self.is_shape_numeric:
            return self._expand_array().tolist()
        raise ValueError("A symbolic array cannot be expanded to a list")
    def tomatrix(self):
        from ...matrices import Matrix
        if not self.is_shape_numeric:
            raise ValueError("A symbolic array cannot be expanded to a matrix")
        if self.rank() != 2:
            raise ValueError('Dimensions must be of size of 2')
        return Matrix(self._expand_array().tomatrix())

    def _rebuild(self, function, limits): return type(self)(function, *limits)
    def xreplace(self, rule):
        from ... import Tuple
        if self in rule:
            return rule[self]
        if self._mapped:
            return self._rebuild(None, Tuple(*self.limits).xreplace(rule).args)
        replaced = Tuple(*self.args).xreplace(rule).args
        return self._rebuild(replaced[0], replaced[1:])
    def subs(self, *args, **options):
        function, limits = _core._array_comprehension_subs(
            None if self._mapped else self.function, self.limits, args, options or None, self._mapped)
        return self._rebuild(function, limits)

    def __eq__(self, other): return type(self) is type(other) and self.args == other.args
    def __ne__(self, other): return not self == other
    def __hash__(self): return hash((type(self).__name__, self.args))
    def __repr__(self): return '%s(%s)' % (type(self).__name__, ', '.join(map(str, self.args)))
    __str__ = __repr__
    def __reduce__(self): return type(self), self.args
    def to_sympy(self):
        from sympy.tensor.array import ArrayComprehension
        return ArrayComprehension(self.function.to_sympy(),
            *(tuple(value.to_sympy() for value in limit) for limit in self.limits))


def _is_lambda(value):
    """SymPy's ``isLambda``."""
    return inspect.isfunction(value) and value.__name__ == '<lambda>'


class ArrayComprehensionMap(ArrayComprehension):
    """A comprehension whose entries are a lambda of the product of the indices."""
    __slots__ = ('_lambda',)
    _mapped = True

    def __new__(cls, function, *symbols, **assumptions):
        if any(len(limit) != 3 or None for limit in symbols):
            raise ValueError('ArrayComprehension requires values lower and upper bound'
                             ' for the expression')
        if not _is_lambda(function):
            raise ValueError('Data type not supported')
        limits = _limits(symbols)
        if assumptions:
            raise TypeError('Basic.__new__() got an unexpected keyword argument %r' % next(iter(assumptions)))
        result = object.__new__(cls)
        result._args = limits
        result._lambda = function
        return result

    # SymPy stores only the limits as args: ``function`` is ``args[0]``.
    @property
    def limits(self): return self._args
    @property
    def func(self):
        function = self._lambda
        return lambda *args, **kwargs: ArrayComprehensionMap(function, *args, **kwargs)
    def _callable(self): return self._lambda
    def _rebuild(self, function, limits): return ArrayComprehensionMap(self._lambda, *limits)
    def __reduce__(self): return ArrayComprehensionMap, (self._lambda, *self.limits)
    def to_sympy(self):
        from sympy.tensor.array import ArrayComprehensionMap
        return ArrayComprehensionMap(self._lambda,
            *(tuple(value.to_sympy() for value in limit) for limit in self.limits))


__all__ = ['ArrayComprehension', 'ArrayComprehensionMap']
