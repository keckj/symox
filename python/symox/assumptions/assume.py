"""Native predicates and mutable assumption contexts."""
from contextlib import contextmanager
from .. import _core


class _PredicateType(type):
    def __instancecheck__(cls, value):
        if cls in type(value).__mro__:return True
        return isinstance(value, _core.Expr) and (
            value.class_name.endswith("Predicate")
            and value.class_name != "AppliedPredicate")


class Predicate(metaclass=_PredicateType):
    def __new__(cls, name):
        if not isinstance(name, str):
            raise TypeError("Predicate expects a string name")
        return _core._predicate(name, undefined=True)


class _AppliedType(type):
    def __instancecheck__(cls, value):
        if cls in type(value).__mro__:return True
        return isinstance(value, _core.Expr) and value.class_name in (
            "AppliedPredicate", "AppliedBinaryRelation")


class AppliedPredicate(metaclass=_AppliedType):
    def __new__(cls, predicate, *args):
        return _core._applied_predicate((predicate, *args))


class AssumptionsContext(set):
    """Mutable set of default assumptions for ``ask`` and ``refine``."""
    def add(self, *assumptions):
        for assumption in assumptions:
            super().add(assumption)

    def __repr__(self):
        if not self:
            return "AssumptionsContext()"
        return "AssumptionsContext({" + ", ".join(sorted(map(str, self))) + "})"


global_assumptions = AssumptionsContext()


@contextmanager
def assuming(*assumptions):
    """Temporarily extend global assumptions, restoring them on exit."""
    previous = global_assumptions.copy()
    try:
        global_assumptions.add(*assumptions)
        yield
    finally:
        global_assumptions.clear()
        global_assumptions.update(previous)

PredicateMeta=type(Predicate)

class _UndefinedPredicateType(_PredicateType):
    def __instancecheck__(cls,value):return isinstance(value,_core.Expr) and value.class_name=='UndefinedPredicate'
class UndefinedPredicate(Predicate,metaclass=_UndefinedPredicateType):pass
