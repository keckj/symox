"""Native predicates and mutable assumption contexts."""
from contextlib import contextmanager
from .. import _core


class PredicateMeta(type):
    def __instancecheck__(cls, value):
        if cls in type(value).__mro__:return True
        return isinstance(value, _core.Expr) and (
            value.class_name.endswith("Predicate")
            and value.class_name != "AppliedPredicate")


# SymPy's metaclass name, shown by ``register`` on an undefined predicate.
_PredicateType = PredicateMeta


class Predicate(metaclass=_PredicateType):
    """``Predicate("P")`` builds an undefined predicate; a subclass with a
    ``name`` attribute builds a predicate whose handlers are registered with
    ``register`` (SymPy's dispatch-based predicates)."""

    def __new__(cls, *args, **kwargs):
        if cls is Predicate or cls.__dict__.get('_symox_undefined'):
            if len(args) != 1 or kwargs:
                raise TypeError("Predicate expects a string name")
            name, = args
            if not isinstance(name, str):
                raise TypeError("Predicate expects a string name")
            return _core._predicate(name, undefined=True)
        name = getattr(cls, 'name', None)
        if not isinstance(name, str):
            raise TypeError("Predicate subclasses define a string name")
        obj = object.__new__(cls)
        obj._symox_predicate = _core._predicate(name, undefined=True)
        obj._symox_args = args
        return obj

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if 'handler' not in cls.__dict__:
            cls.handler = _PredicateDispatcher(cls.__dict__.get('name', cls.__name__))

    # The applied form is the native applied predicate of the same name, so
    # ``ask`` reaches the registered handlers through ``register_handler``.
    def __call__(self, *args):
        return _core._applied_predicate((self._symox_predicate, *args))

    def _symox_(self):
        return self._symox_predicate

    @property
    def args(self):
        return self._symox_args

    @property
    def func(self):
        return type(self)

    def __eq__(self, other):
        if isinstance(other, Predicate) and hasattr(other, '_symox_predicate'):
            other = other._symox_predicate
        return self._symox_predicate == other

    def __hash__(self):
        return hash(self._symox_predicate)

    def __str__(self):
        return str(self._symox_predicate)

    __repr__ = __str__

    @classmethod
    def register(cls, *types, **kwargs):
        """Register a handler for the given argument types."""
        return cls.handler.register(*types, **kwargs)

    @classmethod
    def register_many(cls, *types, **kwargs):
        """Register the same handler for several type signatures."""
        def _(func):
            for t in types:
                if not isinstance(t, tuple):
                    t = (t,)
                cls.register(*t, **kwargs)(func)
            return func
        return _

    def eval(self, args, assumptions=True):
        return self.handler(*args, assumptions=assumptions)


class _PredicateDispatcher:
    """Type dispatch of predicate handlers, most specific signature first."""

    def __init__(self, name):
        self.name = name
        self.funcs = []
        self._installed = False

    def register(self, *types, **kwargs):
        def _(func):
            self.funcs.append((types, func))
            if not self._installed:
                from .._public_tools import register_handler
                register_handler(self.name, self)
                self._installed = True
            return func
        return _

    def dispatch_iter(self, *args):
        matches = [(types, func) for types, func in self.funcs
                   if len(types) == len(args)
                   and all(isinstance(a, t) for a, t in zip(args, types))]
        def specific(entry):
            return -sum(1 for other, _ in matches if other != entry[0]
                        and all(issubclass(a, b) for a, b in zip(entry[0], other)))
        for _types, func in sorted(matches, key=specific):
            yield func

    def __call__(self, *args, assumptions=True):
        from ..multipledispatch import MDNotImplementedError
        result = None
        for func in self.dispatch_iter(*args):
            try:
                result = func(*args, assumptions)
            except MDNotImplementedError:
                continue
            if result is not None:
                return result
        return result


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


class _UndefinedPredicateType(_PredicateType):
    def __instancecheck__(cls,value):return isinstance(value,_core.Expr) and value.class_name=='UndefinedPredicate'
class UndefinedPredicate(Predicate,metaclass=_UndefinedPredicateType):_symox_undefined=True
