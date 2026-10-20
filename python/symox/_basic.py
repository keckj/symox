"""Structural Basic objects and the native symbolic object category.

The category includes scalar nodes, symbolic matrices and immutable explicit
matrices. Virtual subclass registration does not change their concrete Python
MRO; nominal hierarchy compatibility remains a separate obligation.
"""
from abc import ABCMeta
from . import _core


class _BasicMeta(ABCMeta):
    def __instancecheck__(cls, value):
        if cls is Basic and isinstance(value, _core.Matrix):
            return value._immutable
        return super().__instancecheck__(value)

    # SymPy's Basic is the base of every native node's nominal class, so
    # attributes patched onto it after import (SageMath's ``_sage_``) must
    # reach native nodes as well as structural Basic subclasses.
    def __setattr__(cls, name, value):
        super().__setattr__(name, value)
        if globals().get('Basic') is cls:
            from ._node_types import forward_setattr
            forward_setattr(cls, 'Basic', name, value)

    def __delattr__(cls, name):
        super().__delattr__(name)
        if globals().get('Basic') is cls:
            from ._node_types import forward_delattr
            forward_delattr(cls, 'Basic', name)


class Basic(metaclass=_BasicMeta):
    """Immutable, reconstructible structural container for symbolic arguments."""
    # No instance slots: classes deriving from Basic and from the native
    # Expr type (Devito's Pointer) need a layout compatible with both.
    is_Atom = is_Symbol = is_Number = is_Add = is_Mul = is_Pow = False
    is_Function = is_Relational = is_Matrix = is_Boolean = False
    is_commutative = None

    def __new__(cls, *args):
        if issubclass(cls, _core.Expr):
            # An expression class: SymPy's Basic.__new__ builds the native
            # unevaluated node of this class (see _host.expression_class).
            return _core.Expr.__new__(cls, *args)
        result = object.__new__(cls)
        result._args = args
        return result

    @property
    def args(self):
        return self._args

    @property
    def func(self):
        return type(self)

    @property
    def free_symbols(self):
        return set().union(*(arg.free_symbols for arg in self.args))

    @property
    def assumptions0(self):
        return {}

    @classmethod
    def fromiter(cls, args, **assumptions):
        return cls(*tuple(args), **assumptions)

    def __repr__(self):
        return type(self).__name__ + '(' + ', '.join(map(str, self.args)) + ')'

    __str__ = __repr__

    def _hashable_content(self):
        return self._args

    def __eq__(self, other):
        return type(self) is type(other) and self._hashable_content() == other._hashable_content()

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash((type(self).__name__,) + tuple(self._hashable_content()))

    def __getnewargs__(self):
        return self.args

    def srepr(self):
        from . import srepr
        return type(self).__name__ + '(' + ', '.join(srepr(arg) for arg in self.args) + ')'

    def xreplace(self, rule):
        if self in rule:
            return rule[self]
        args = tuple(arg.xreplace(rule) if hasattr(arg, 'xreplace') else arg
                     for arg in self.args)
        return self if args == self.args else self.func(*args)

    def subs(self, arg1, arg2=None, **kwargs):
        from . import sympify, Dummy, default_sort_key
        simultaneous = kwargs.pop('simultaneous', False)
        if kwargs:
            raise NotImplementedError('Unsupported formal Basic substitution options: ' + ', '.join(kwargs))
        mapping = arg2 is None and hasattr(arg1, 'items')
        source = arg1.items() if mapping else arg1 if arg2 is None else [(arg1, arg2)]
        replacements = [(sympify(old), sympify(new)) for old, new in source]
        if mapping:
            from ._matrixexpr_methods import _walk
            replacements.sort(key=lambda pair: (-sum(1 for _ in _walk(pair[0])), default_sort_key(pair[0])))
        restore = {}
        if simultaneous:
            masked = []
            for old, new in replacements:
                temporary = Dummy('subs_basic')
                masked.append((old, temporary))
                restore[temporary] = new
            replacements = masked
        result = self
        for old, new in replacements:
            if result == old:
                result = new
            elif type(result) is not type(self) and Basic not in type(result).__mro__:
                result = result.subs(old, new)
            else:
                children = tuple(arg.subs(old, new) if hasattr(arg, 'subs') else arg
                                 for arg in result.args)
                result = result.func(*children)
        return result.xreplace(restore) if simultaneous else result

    def has(self, *patterns):
        from ._structural_queries import matcher
        tests = [matcher(pattern, matching=False) for pattern in patterns]
        from ._matrixexpr_methods import _walk
        return any(test(value) for value in _walk(self) for test in tests)

    def atoms(self, *types):
        if not self.args:
            return {self} if not types or isinstance(self, types) else set()
        return set().union(*(arg.atoms(*types) for arg in self.args))

    def doit(self, **hints):
        if not hints.get('deep', True):
            return self
        return self.func(*(arg.doit(**hints) if hasattr(arg, 'doit') else arg for arg in self.args))

    def _symox_(self):
        raise NotImplementedError('Formal Basic containers are not native scalar expressions')


Basic.register(_core.Expr)
Basic.register(_core._MatrixExpr)
Basic.register(_core.AlgebraicNumber)
