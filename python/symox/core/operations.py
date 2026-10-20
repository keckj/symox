"""Associative and lattice operation categories.

These follow sympy/core/operations.py (SymPy 1.14). Native Add, Mul, And, Or,
Min and Max nodes, and the Python matrix, tensor and set operations that play
the same role, are instances of the corresponding category, so downstream
``isinstance`` and ``make_args`` checks see the reference hierarchy.

The categories carry no instance layout of their own. A Python class may then
combine them with the native expression type (``class M(Expr, LatticeOp)``) as
in the reference; constructing instances of such a class relies on the native
expression storage of that class.
"""
from symox import _core
from symox._basic import Basic as _Basic
from symox._compat import _ExprType
from symox.utilities.exceptions import sympy_deprecation_warning
from symox.utilities.iterables import sift

from .sympify import _sympify as _sympify_, sympify
from .cache import cacheit
from .sorting import ordered
from .logic import fuzzy_and
from .parameters import global_parameters


# Native class names of the reference AssocOp and LatticeOp subclasses.
_NATIVE_ASSOC_OPS = frozenset(('Add', 'Mul', 'And', 'Or', 'Min', 'Max'))
_NATIVE_LATTICE_OPS = frozenset(('And', 'Or', 'Min', 'Max'))

# Python classes that implement reference AssocOp subclasses (resolved lazily
# because they are defined after this module can be imported).
_PYTHON_ASSOC_OPS = ('symox:Add', 'symox:Mul', 'symox:MatAdd', 'symox:MatMul',
                     'symox:Union', 'symox:Intersection',
                     'symox.tensor.tensor:TensAdd', 'symox.tensor.tensor:TensMul')
_PYTHON_LATTICE_OPS = ('symox:Union', 'symox:Intersection')


def _resolve(paths):
    from importlib import import_module
    classes = []
    for path in paths:
        module, name = path.split(':')
        classes.append(getattr(import_module(module), name))
    return tuple(classes)


class _OperationCategory(_ExprType):
    """Metaclass giving a category reference instance and subclass checks.

    Real subclasses use the ordinary checks; only the category classes defined
    here recognize native nodes and the Python classes listed above.
    """
    def __instancecheck__(cls, value):
        if type.__instancecheck__(cls, value):
            return True
        category = cls.__dict__.get('_category')
        if category is None:
            return False
        native_names, python_paths = category
        if isinstance(value, _core.Expr):
            return value.class_name in native_names
        return isinstance(value, _resolve(python_paths))

    def __subclasscheck__(cls, sub):
        if type.__subclasscheck__(cls, sub):
            return True
        category = cls.__dict__.get('_category')
        if category is None or not isinstance(sub, type):
            return False
        from symox._node_types import NodeType
        if isinstance(sub, NodeType):
            # The nominal classes of native nodes (class dispatch keyed by
            # LatticeOp or Max accepts ``type(Max(x, y))``).
            return _ExprType.__subclasscheck__(cls, sub)
        return any(sub is c or type.__subclasscheck__(c, sub)
                   for c in _resolve(category[1]))


def _basic_new(cls, *args):
    # Reference: Basic.__new__(cls, *args) for a concrete subclass. The basic
    # storage of a Python operation class is the native expression type.
    if issubclass(cls, _core.Expr):
        return _core.Expr.__new__(cls, *args)
    raise TypeError(
        "%s defines no expression storage: derive it from Expr (or another "
        "native expression class) as well as from %s" % (
            cls.__name__, 'LatticeOp' if issubclass(cls, LatticeOp) else 'AssocOp'))


class AssocOp(metaclass=_OperationCategory):
    """ Associative operations, can separate noncommutative and
    commutative parts.

    (a op b) op c == a op (b op c) == a op b op c.

    Base class for Add and Mul.

    This is an abstract base class, concrete derived classes must define
    the attribute `identity`.

    Parameters
    ==========

    *args :
        Arguments which are operated

    evaluate : bool, optional
        Evaluate the operation. If not passed, refer to ``global_parameters.evaluate``.
    """
    _category = (_NATIVE_ASSOC_OPS, _PYTHON_ASSOC_OPS)

    _args_type = None

    @cacheit
    def __new__(cls, *args, evaluate=None, _sympify=True):
        # Allow faster processing by passing ``_sympify=False``, if all arguments
        # are already sympified.
        if _sympify:
            args = list(map(_sympify_, args))

        # Disallow non-Expr args in Add/Mul
        typ = cls._args_type
        if typ is not None:
            from .relational import Relational
            if any(isinstance(arg, Relational) for arg in args):
                raise TypeError("Relational cannot be used in %s" % cls.__name__)

            # This should raise TypeError once deprecation period is over:
            for arg in args:
                if not isinstance(arg, typ):
                    sympy_deprecation_warning(
                        f"""

Using non-Expr arguments in {cls.__name__} is deprecated (in this case, one of
the arguments has type {type(arg).__name__!r}).

If you really did intend to use a multiplication or addition operation with
this object, use the * or + operator instead.

                        """,
                        deprecated_since_version="1.7",
                        active_deprecations_target="non-expr-args-deprecated",
                        stacklevel=4,
                    )

        if evaluate is None:
            evaluate = global_parameters.evaluate
        if not evaluate:
            obj = cls._from_args(args)
            obj = cls._exec_constructor_postprocessors(obj)
            return obj

        args = [a for a in args if a is not cls.identity]

        if len(args) == 0:
            return cls.identity
        if len(args) == 1:
            return args[0]

        c_part, nc_part, order_symbols = cls.flatten(args)
        is_commutative = not nc_part
        obj = cls._from_args(c_part + nc_part, is_commutative)
        obj = cls._exec_constructor_postprocessors(obj)

        if order_symbols is not None:
            from symox import Order
            return Order(obj, *order_symbols)
        return obj

    @classmethod
    def _exec_constructor_postprocessors(cls, obj):
        # Native expressions carry no constructor postprocessors.
        return obj

    @classmethod
    def _from_args(cls, args, is_commutative=None):
        """Create new instance with already-processed args.
        If the args are not in canonical order, then a non-canonical
        result will be returned, so use with caution. The order of
        args may change if the sign of the args is changed."""
        if len(args) == 0:
            return cls.identity
        elif len(args) == 1:
            return args[0]

        obj = _basic_new(cls, *args)
        if is_commutative is None:
            is_commutative = fuzzy_and(a.is_commutative for a in args)
        obj.is_commutative = is_commutative
        return obj

    def _new_rawargs(self, *args, reeval=True, **kwargs):
        """Create new instance of own class with args exactly as provided by
        caller but returning the self class identity if args is empty.

        Commutativity is recomputed if self is non-commutative and
        ``reeval=False`` has not been passed.
        """
        if reeval and self.is_commutative is False:
            is_commutative = None
        else:
            is_commutative = self.is_commutative
        return self._from_args(args, is_commutative)

    @classmethod
    def flatten(cls, seq):
        """Return seq so that none of the elements are of type `cls`. This is
        the vanilla routine that will be used if a class derived from AssocOp
        does not define its own flatten routine."""
        # apply associativity, no commutativity property is used
        new_seq = []
        while seq:
            o = seq.pop()
            if o.__class__ is cls:  # classes must match exactly
                seq.extend(o.args)
            else:
                new_seq.append(o)
        new_seq.reverse()

        # c_part, nc_part, order_symbols
        return [], new_seq, None

    def _has_matcher(self):
        """Helper for .has() that checks for containment of
        subexpressions within an expr by using sets of args
        of similar nodes, e.g. x + 1 in x + y + 1 checks
        to see that {x, 1} & {x, y, 1} == {x, 1}
        """
        def _ncsplit(expr):
            # this is not the same as args_cnc because here
            # we don't assume expr is a Mul -- hence deal with args --
            # and always return a set.
            cpart, ncpart = sift(expr.args,
                lambda arg: arg.is_commutative is True, binary=True)
            return set(cpart), ncpart

        c, nc = _ncsplit(self)
        cls = self.__class__

        def is_in(expr):
            if isinstance(expr, cls):
                if expr == self:
                    return True
                _c, _nc = _ncsplit(expr)
                if (c & _c) == c:
                    if not nc:
                        return True
                    elif len(nc) <= len(_nc):
                        for i in range(len(_nc) - len(nc) + 1):
                            if _nc[i:i + len(nc)] == nc:
                                return True
            return False
        return is_in

    @classmethod
    def make_args(cls, expr):
        """
        Return a sequence of elements `args` such that cls(*args) == expr

        Examples
        ========

        >>> from sympy import Symbol, Mul, Add
        >>> x, y = map(Symbol, 'xy')

        >>> Mul.make_args(x*y)
        (x, y)
        >>> Add.make_args(x*y)
        (x*y,)
        >>> set(Add.make_args(x*y + y)) == set([y, x*y])
        True

        """
        if isinstance(expr, cls):
            return expr.args
        else:
            return (sympify(expr),)

    def doit(self, **hints):
        if hints.get('deep', True):
            terms = [term.doit(**hints) for term in self.args]
        else:
            terms = self.args
        return self.func(*terms, evaluate=True)


class ShortCircuit(Exception):
    pass


class LatticeOp(AssocOp):
    """
    Join/meet operations of an algebraic lattice[1].

    Explanation
    ===========

    These binary operations are associative (op(op(a, b), c) = op(a, op(b, c))),
    commutative (op(a, b) = op(b, a)) and idempotent (op(a, a) = op(a) = a).
    Common examples are AND, OR, Union, Intersection, max or min. They have an
    identity element (op(identity, a) = a) and an absorbing element
    conventionally called zero (op(zero, a) = zero).

    This is an abstract base class, concrete derived classes must declare
    attributes zero and identity. All defining properties are then respected.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Lattice_%28order%29
    """
    _category = (_NATIVE_LATTICE_OPS, _PYTHON_LATTICE_OPS)

    is_commutative = True

    def __new__(cls, *args, **options):
        args = (_sympify_(arg) for arg in args)

        try:
            # /!\ args is a generator and _new_args_filter
            # must be careful to handle as such; this
            # is done so short-circuiting can be done
            # without having to sympify all values
            _args = frozenset(cls._new_args_filter(args))
        except ShortCircuit:
            return sympify(cls.zero)
        if not _args:
            return sympify(cls.identity)
        elif len(_args) == 1:
            return set(_args).pop()
        else:
            # XXX in almost every other case for __new__, *_args is
            # passed along, but the expectation here is for _args
            obj = _basic_new(cls, *ordered(_args))
            obj._argset = _args
            return obj

    @classmethod
    def _new_args_filter(cls, arg_sequence, call_cls=None):
        """Generator filtering args"""
        ncls = call_cls or cls
        for arg in arg_sequence:
            if arg == ncls.zero:
                raise ShortCircuit(arg)
            elif arg == ncls.identity:
                continue
            elif arg.func == ncls:
                yield from arg.args
            else:
                yield arg

    @classmethod
    def make_args(cls, expr):
        """
        Return a set of args such that cls(*arg_set) == expr.
        """
        if isinstance(expr, cls):
            argset = getattr(expr, '_argset', None)
            # Native lattice nodes keep their (duplicate free) args only.
            return argset if argset is not None else frozenset(expr.args)
        else:
            return frozenset([sympify(expr)])


_Basic.register(AssocOp)

__all__ = ['AssocOp', 'ShortCircuit', 'LatticeOp']
