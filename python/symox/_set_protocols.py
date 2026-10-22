"""SymPy 1.14 set protocols missing on the native set wrappers.

``Set.kind`` (``SetKind``), ``Set.powerset`` and the iteration protocol of
``ProductSet`` follow ``sympy/sets/sets.py`` and ``fancysets.py``
(BSD-3-Clause, see LICENSE-SYMPY).
"""
from functools import reduce


def _set_kind(self):
    from .sets.sets import SetKind
    from ._kind import NumberKind, UndefinedKind
    from .core.containers import TupleKind
    import symox as s
    name = self.class_name
    if name in ('Interval', 'Reals', 'Integers', 'Naturals', 'Naturals0',
                'Rationals', 'Complexes', 'Range'):
        return SetKind(NumberKind)
    if name == 'EmptySet':
        return SetKind()
    if name == 'UniversalSet':
        return SetKind(UndefinedKind)
    args = self.args
    if name == 'FiniteSet':
        if not args:
            return SetKind()
        kinds = [a.kind for a in args]
        if all(k == kinds[0] for k in kinds):
            return SetKind(kinds[0])
        return SetKind(UndefinedKind)
    if name == 'ProductSet':
        return SetKind(TupleKind(*(a.kind.element_kind for a in args)))
    if name == 'Union':
        kinds = tuple(a.kind for a in args if a != s.EmptySet)
        if not kinds:
            return SetKind()
        if all(k == kinds[0] for k in kinds):
            return kinds[0]
        return SetKind(UndefinedKind)
    if name == 'Intersection':
        kinds = tuple(a.kind for a in args if a != s.UniversalSet)
        if not kinds:
            return SetKind(UndefinedKind)
        if all(k == kinds[0] for k in kinds):
            return kinds[0]
        return SetKind()
    if name == 'Complement':
        return args[0].kind
    if name == 'ImageSet':
        return SetKind(self.lamda.expr.kind)
    if name == 'ConditionSet':
        return SetKind(self.sym.kind)
    return SetKind(UndefinedKind)


def _eval_powerset(self):
    import symox as s
    if self.class_name == 'EmptySet':
        return s.FiniteSet(self)
    if self.class_name == 'FiniteSet':
        from .utilities.iterables import subsets
        return s.FiniteSet(*[s.FiniteSet(*sub) for sub in subsets(self.args)])
    from ._set_extensions import PowerSet
    return PowerSet(self)


def _powerset(self):
    """Find the power set of ``self`` (SymPy's ``Set.powerset``)."""
    return self._eval_powerset()


def _product_is_iterable(self):
    return all(a.is_iterable for a in self.sets)


def _product_iter(self):
    from .utilities.iterables import iproduct
    return iproduct(*self.sets)


def _product_len(self):
    return reduce(lambda a, b: a*b, (len(a) for a in self.args))


def install(namespace):
    Set = namespace['Set']
    if 'kind' not in Set.__dict__:
        Set.kind = property(_set_kind)
    Set._eval_powerset = _eval_powerset
    Set.powerset = _powerset
    product = namespace['ProductSet']
    product.is_iterable = property(_product_is_iterable)
    product.__iter__ = _product_iter
    product.__len__ = _product_len
