"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'true false to_cnf to_dnf to_nnf POSform SOPform simplify_logic bool_map')


def gateinputcount(expr):
    """Count wires entering Boolean gates, including nested gate occurrences."""
    from symox import sympify
    expr=sympify(expr)
    gates={'And','Or','Not','Xor','Nand','Nor','ITE','Implies','Equivalent'}
    if expr.class_name not in gates:return 0
    return len(expr.args)+sum(gateinputcount(a) for a in expr.args)
__all__ += ['gateinputcount']


from symox._basic import Basic, _BasicMeta


class _BooleanType(_BasicMeta):
    def __instancecheck__(cls, value):
        if type.__instancecheck__(cls, value):
            return True
        from symox import _core
        if not isinstance(value, _core.Expr):
            return False
        native_names = cls.__dict__.get('_native_names')
        if native_names is not None:
            return value.class_name in native_names
        if cls not in (Boolean, BooleanFunction):
            return False
        if cls is BooleanFunction:
            return value.is_Boolean and value.is_Function
        return value._is_boolean_category

    # The Boolean atoms are the nominal classes of the native singletons:
    # class dispatch on them accepts those node types, and attributes patched
    # onto them after import (SageMath's ``_sage_``) reach ``true``/``false``.
    def _native_counterpart(cls):
        if '_native_names' not in cls.__dict__ or globals().get(cls.__name__) is not cls:
            return None
        from symox._node_types import node_class
        return node_class(cls.__name__)

    def __subclasscheck__(cls, subclass):
        if super().__subclasscheck__(subclass):
            return True
        from symox._node_types import NodeType
        nominal = cls._native_counterpart() if isinstance(subclass, NodeType) else None
        return nominal is not None and issubclass(subclass, nominal)

    def __repr__(cls):
        # SymPy's Boolean functions are applications, whose FunctionClass
        # metaclass prints them by bare name.
        if 'BooleanFunction' in globals() and type.__subclasscheck__(BooleanFunction, cls):
            return cls.__name__
        return super().__repr__()

    def __setattr__(cls, name, value):
        super().__setattr__(name, value)
        if cls._native_counterpart() is not None:
            from symox._node_types import forward_setattr
            forward_setattr(cls, cls.__name__, name, value)

    def __delattr__(cls, name):
        super().__delattr__(name)
        if cls._native_counterpart() is not None:
            from symox._node_types import forward_delattr
            forward_delattr(cls, cls.__name__, name)


class Boolean(Basic, metaclass=_BooleanType):
    """Boolean category and formal Basic container.

    Native gate, relation and Boolean-atom operations remain in the native
    kernel. Algebra involving a directly constructed formal base container
    is not yet implemented by that kernel.
    """
    __slots__ = ()
    is_Boolean = False
    is_Function = False
    is_Atom = False
    is_Symbol = False
    is_Number = False
    is_Relational = False



class BooleanFunction(Boolean):
    """Category for native gates and formal Boolean function applications."""
    is_Boolean = True
    is_Function = True

    def __new__(cls, *args):
        from symox import sympify
        return Boolean.__new__(cls, *(sympify(arg) for arg in args))


class BooleanAtom(Boolean):
    """Base class of :py:class:`~.BooleanTrue` and :py:class:`~.BooleanFalse`.

    The two atoms are the native singletons ``true`` and ``false``.
    """
    _native_names = ('BooleanTrue', 'BooleanFalse')
    is_Boolean = True
    is_Atom = True
    _op_priority = 11  # higher than Expr


class BooleanTrue(BooleanAtom):
    """SymPy version of ``True``, a singleton that can be accessed via ``S.true``."""
    _native_names = ('BooleanTrue',)

    def __new__(cls):
        from symox import true
        return true


class BooleanFalse(BooleanAtom):
    """SymPy version of ``False``, a singleton that can be accessed via ``S.false``."""
    _native_names = ('BooleanFalse',)

    def __new__(cls):
        from symox import false
        return false


__all__ += ['Boolean', 'BooleanFunction', 'BooleanAtom', 'BooleanTrue', 'BooleanFalse']


# The gates are classes, as in SymPy: they key class dispatch
# (``singledispatch.register(And)``), ``make_args`` and ``expr.func``
# comparisons. Their applications are native nodes (instances of the
# nominal classes of the same name, which these classes accept in instance
# and subclass checks); construction evaluates in the native kernel.

from symox import _core  # noqa: E402
from symox.core.operations import LatticeOp, _OperationCategory  # noqa: E402


class _BooleanLatticeType(_BooleanType, _OperationCategory):
    """Metaclass of And and Or, which are both lattice operations and
    Boolean functions."""


def _gate_constructor(name):
    native = getattr(_core, name)

    def __new__(cls, *args, **options):
        if cls.__dict__.get('_symox_head_name') != name:
            # A user subclass: a formal application of that class.
            return BooleanFunction.__new__(cls, *args)
        if name in ('And', 'Or'):
            # LatticeOp.__new__ always evaluates and ignores its options.
            return native(*args)
        evaluate = options.pop('evaluate', None)
        if options:
            if name == 'Equivalent':
                raise TypeError("Equivalent.__new__() got an unexpected keyword argument '%s'"
                                % next(iter(options)))
            raise ValueError("Unknown options: %s" % options)
        if evaluate is False:
            from symox import sympify
            return _core._apply(name, *(sympify(arg) for arg in args), evaluate=False)
        return native(*args)

    __new__.__doc__ = native.__doc__
    return __new__


class And(LatticeOp, BooleanFunction, metaclass=_BooleanLatticeType):
    """Logical AND function.

    It evaluates its arguments in order, returning false immediately
    when an argument is false and true if they are all true.
    """
    _native_names = ('And',)
    _symox_head_name = 'And'
    __new__ = _gate_constructor('And')


class Or(LatticeOp, BooleanFunction, metaclass=_BooleanLatticeType):
    """Logical OR function.

    It evaluates its arguments in order, returning true immediately
    when an argument is true, and false if they are all false.
    """
    _native_names = ('Or',)
    _symox_head_name = 'Or'
    __new__ = _gate_constructor('Or')


class Not(BooleanFunction):
    """Logical Not function (negation)."""
    _native_names = ('Not',)
    _symox_head_name = 'Not'
    __new__ = _gate_constructor('Not')


class Xor(BooleanFunction):
    """Logical XOR (exclusive OR) function."""
    _native_names = ('Xor',)
    _symox_head_name = 'Xor'
    __new__ = _gate_constructor('Xor')


class Nand(BooleanFunction):
    """Logical NAND function."""
    _native_names = ('Nand',)
    _symox_head_name = 'Nand'
    __new__ = _gate_constructor('Nand')


class Nor(BooleanFunction):
    """Logical NOR function."""
    _native_names = ('Nor',)
    _symox_head_name = 'Nor'
    __new__ = _gate_constructor('Nor')


class Xnor(BooleanFunction):
    """Logical XNOR function."""
    _native_names = ('Xnor',)
    _symox_head_name = 'Xnor'
    __new__ = _gate_constructor('Xnor')


class Implies(BooleanFunction):
    """Logical implication: A implies B is equivalent to (not A) or B."""
    _native_names = ('Implies',)
    _symox_head_name = 'Implies'
    __new__ = _gate_constructor('Implies')


class Equivalent(BooleanFunction):
    """Equivalence relation: true when all the arguments are logically
    equivalent."""
    _native_names = ('Equivalent',)
    _symox_head_name = 'Equivalent'
    __new__ = _gate_constructor('Equivalent')


class ITE(BooleanFunction):
    """If-then-else clause: ``ITE(A, B, C)`` is B if A is true, else C."""
    _native_names = ('ITE',)
    _symox_head_name = 'ITE'
    __new__ = _gate_constructor('ITE')


def _install_lattice_elements():
    from symox import true, false
    And.zero, And.identity = false, true
    Or.zero, Or.identity = true, false


_install_lattice_elements()

__all__ += ['And', 'Or', 'Not', 'Xor', 'Nand', 'Nor', 'Xnor', 'Implies', 'Equivalent', 'ITE']


def simplify_univariate(expr):
    """return a simplified version of univariate boolean expression, else ``expr``"""
    from symox import S, Piecewise, Eq, Ne, And, Or, true, false
    from symox import _core
    if not isinstance(expr, BooleanFunction):
        return expr
    if expr.atoms(Eq, Ne):
        return expr
    c = expr
    free = c.free_symbols
    if len(free) != 1:
        return c
    x = free.pop()
    ok, i = _core._piecewise_intervals(Piecewise((0, c), evaluate=False), x, True)
    if not ok:
        return c
    if not i:
        return false
    args = []
    # the reference rebinds ``c`` inside this loop and so tests the
    # endpoints of later intervals on an earlier piece, which drops x = 0
    # from Or(x <= 0, x >= 1) when that piece comes first; the endpoints
    # are tested on the whole condition here
    for a, b, _, _ in i:
        if a is S.NegativeInfinity or a == S.NegativeInfinity:
            if b == S.Infinity:
                c = true
            else:
                if expr.subs(x, b) == True:  # noqa: E712
                    c = (x <= b)
                else:
                    c = (x < b)
        else:
            incl_a = (expr.subs(x, a) == True)  # noqa: E712
            incl_b = (expr.subs(x, b) == True)  # noqa: E712
            if incl_a and incl_b:
                if b.is_infinite:
                    c = (x >= a)
                else:
                    c = And(a <= x, x <= b)
            elif incl_a:
                c = And(a <= x, x < b)
            elif incl_b:
                if b.is_infinite:
                    c = (x > a)
                else:
                    c = And(a < x, x <= b)
            else:
                c = And(a < x, x < b)
        args.append(c)
    return Or(*args)
__all__ += ['simplify_univariate']
