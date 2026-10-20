"""Identity, root and lattice functions (sympy/functions/elementary/miscellaneous.py)."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Id sqrt cbrt root real_root Rem')

from symox import _core  # noqa: E402
from symox._expr_facade import Expr  # noqa: E402
from symox.core.function import Application  # noqa: E402
from symox.core.operations import LatticeOp, ShortCircuit, _OperationCategory, _basic_new  # noqa: E402

# The native constructors: they evaluate native arguments in the kernel.
import symox as _symox  # noqa: E402
_NATIVE = {'Max': _symox.Max, 'Min': _symox.Min}


class _MinMaxType(_OperationCategory):
    """Metaclass of the lattice functions: SymPy's FunctionClass prints
    function classes by their bare name."""

    def __repr__(cls):
        return cls.__name__

    __str__ = __repr__


def _evaluation():
    from symox._public_core import _evaluation as evaluation
    return evaluation.get()


def _native_argument(arg):
    # The native kernel evaluates Max and Min over native arguments other
    # than nested Max and Min nodes, whose collapsing SymPy's algorithm
    # below performs; it rejects the arguments SymPy's filter rejects.
    if not isinstance(arg, _core.Expr) or arg.class_name in ('Max', 'Min') or _core._host_contains(arg)[0]:
        return False
    if arg.is_Number:
        return arg.class_name != 'NaN'
    return arg.is_extended_real is not False and not (arg.is_number and not arg.is_comparable)


class MinMaxBase(Expr, LatticeOp, metaclass=_MinMaxType):
    """Common base of :class:`Max` and :class:`Min`.

    Native nodes of either class are instances. Arguments holding instances
    of user classes (Torch's ``int_oo``), whose order relations are decided
    in Python, and user subclasses, are evaluated by SymPy's algorithm
    below. The native kernel evaluates the other arguments.
    """
    _category = (frozenset(('Max', 'Min')), ())
    _symox_expression_kind = 'expr'

    def __new__(cls, *args, **assumptions):
        evaluate = assumptions.pop('evaluate', None)
        if assumptions:
            raise TypeError(cls.__name__ + ' does not accept these keyword arguments')
        if evaluate is None:
            evaluate = _evaluation()
        native = _NATIVE.get(cls.__name__) if cls.__dict__.get('_symox_head_name') == cls.__name__ else None
        from symox import sympify
        from symox.core.sorting import ordered
        args = [sympify(arg) for arg in args]
        if native is not None and evaluate and all(_native_argument(arg) for arg in args):
            return native(*args)
        if evaluate:
            try:
                args = frozenset(cls._new_args_filter(args))
            except ShortCircuit:
                return cls.zero
            # remove redundant args that are easily identified
            args = cls._collapse_arguments(args, **assumptions)
            # find local zeros
            args = cls._find_localzeros(args, **assumptions)
        args = frozenset(args)

        if not args:
            return cls.identity

        if len(args) == 1:
            return list(args).pop()

        if native is not None:
            return _core._apply(cls.__name__, *ordered(args), evaluate=False)
        obj = _basic_new(cls, *ordered(args))
        try:
            obj._argset = args
        except AttributeError:
            pass
        return obj

    @classmethod
    def _collapse_arguments(cls, args, **assumptions):
        """Remove redundant args.

        Examples
        ========

        >>> from sympy import Min, Max
        >>> from sympy.abc import a, b, c, d, e

        Any arg in parent that appears in any
        parent-like function in any of the flat args
        of parent can be removed from that sub-arg:

        >>> Min(a, Max(b, Min(a, c, d)))
        Min(a, Max(b, Min(c, d)))

        If the arg of parent appears in an opposite-than parent
        function in any of the flat args of parent that function
        can be replaced with the arg:

        >>> Min(a, Max(b, Min(c, d, Max(a, e))))
        Min(a, Max(b, Min(a, c, d)))
        """
        from symox.core.sorting import ordered
        from symox.core.traversal import walk
        from symox.utilities.iterables import sift
        if not args:
            return args
        args = list(ordered(args))
        if cls == Min:
            other = Max
        else:
            other = Min

        # find global comparable max of Max and min of Min if a new
        # value is being introduced in these args at position 0 of
        # the ordered args
        if args[0].is_number:
            sifted = mins, maxs = [], []
            for i in args:
                for v in walk(i, Min, Max):
                    if v.args[0].is_comparable:
                        sifted[isinstance(v, Max)].append(v)
            small = Min.identity
            for i in mins:
                v = i.args[0]
                if v.is_number and (v < small) == True:  # noqa: E712
                    small = v
            big = Max.identity
            for i in maxs:
                v = i.args[0]
                if v.is_number and (v > big) == True:  # noqa: E712
                    big = v
            # at the point when this function is called from __new__,
            # there may be more than one numeric arg present since
            # local zeros have not been handled yet, so look through
            # more than the first arg
            if cls == Min:
                for arg in args:
                    if not arg.is_number:
                        break
                    if (arg < small) == True:  # noqa: E712
                        small = arg
            elif cls == Max:
                for arg in args:
                    if not arg.is_number:
                        break
                    if (arg > big) == True:  # noqa: E712
                        big = arg
            T = None
            if cls == Min:
                if small != Min.identity:
                    other = Max
                    T = small
            elif big != Max.identity:
                other = Min
                T = big
            if T is not None:
                # remove numerical redundancy
                for i in range(len(args)):
                    a = args[i]
                    if isinstance(a, other):
                        a0 = a.args[0]
                        if ((a0 > T) if other == Max else (a0 < T)) == True:  # noqa: E712
                            args[i] = cls.identity

        # remove redundant symbolic args
        def do(ai, a):
            if not isinstance(ai, (Min, Max)):
                return ai
            cond = a in ai.args
            if not cond:
                return ai.func(*[do(i, a) for i in ai.args],
                    evaluate=False)
            if isinstance(ai, cls):
                return ai.func(*[do(i, a) for i in ai.args if i != a],
                    evaluate=False)
            return a
        for i, a in enumerate(args):
            args[i + 1:] = [do(ai, a) for ai in args[i + 1:]]

        # factor out common elements as for
        # Min(Max(x, y), Max(x, z)) -> Max(x, Min(y, z))
        # and vice versa when swapping Min/Max -- do this only for the
        # easy case where all functions contain something in common;
        # trying to find some optimal subset of args to modify takes
        # too long

        def factor_minmax(args):
            is_other = lambda arg: isinstance(arg, other)  # noqa: E731
            other_args, remaining_args = sift(args, is_other, binary=True)
            if not other_args:
                return args

            # Min(Max(x, y, z), Max(x, y, u, v)) -> {x,y}, ({z}, {u,v})
            arg_sets = [set(arg.args) for arg in other_args]
            common = set.intersection(*arg_sets)
            if not common:
                return args

            new_other_args = list(common)
            arg_sets_diff = [arg_set - common for arg_set in arg_sets]

            # If any set is empty after removing common then all can be
            # discarded e.g. Min(Max(a, b, c), Max(a, b)) -> Max(a, b)
            if all(arg_sets_diff):
                other_args_diff = [other(*s, evaluate=False) for s in arg_sets_diff]
                new_other_args.append(cls(*other_args_diff, evaluate=False))

            other_args_factored = other(*new_other_args, evaluate=False)
            return remaining_args + [other_args_factored]

        if len(args) > 1:
            args = factor_minmax(args)

        return args

    @classmethod
    def _new_args_filter(cls, arg_sequence):
        """
        Generator filtering args.

        first standard filter, for cls.zero and cls.identity.
        Also reshape ``Max(a, Max(b, c))`` to ``Max(a, b, c)``,
        and check arguments for comparability
        """
        for arg in arg_sequence:
            # pre-filter, checking comparability of arguments (immutable
            # matrices, matrix expressions in SymPy, are not extended real)
            if not isinstance(arg, Expr) or getattr(arg, 'is_Matrix', False) or arg.is_extended_real is False or (
                    arg.is_number and
                    not arg.is_comparable):
                raise ValueError("The argument '%s' is not comparable." % arg)

            if arg == cls.zero:
                raise ShortCircuit(arg)
            elif arg == cls.identity:
                continue
            elif arg.func == cls:
                yield from arg.args
            else:
                yield arg

    @classmethod
    def _find_localzeros(cls, values, **options):
        """
        Sequentially allocate values to localzeros.

        When a value is identified as being more extreme than another member it
        replaces that member; if this is never true, then the value is simply
        appended to the localzeros.
        """
        localzeros = set()
        for v in values:
            is_newzero = True
            localzeros_ = list(localzeros)
            for z in localzeros_:
                if id(v) == id(z):
                    is_newzero = False
                else:
                    con = cls._is_connected(v, z)
                    if con:
                        is_newzero = False
                        if con is True or con == cls:
                            localzeros.remove(z)
                            localzeros.update([v])
            if is_newzero:
                localzeros.update([v])
        return localzeros

    @classmethod
    def _is_connected(cls, x, y):
        """
        Check if x and y are connected somehow.
        """
        from symox import S, factor_terms
        for i in range(2):
            if x == y:
                return True
            t, f = Max, Min
            for op in "><":
                for j in range(2):
                    try:
                        if op == ">":
                            v = x >= y
                        else:
                            v = x <= y
                    except TypeError:
                        return False  # non-real arg
                    if not v.is_Relational:
                        return t if v else f
                    t, f = f, t
                    x, y = y, x
                x, y = y, x  # run next pass with reversed order relative to start
            # simplification can be expensive, so be conservative
            # in what is attempted
            x = factor_terms(x - y)
            y = S.Zero

        return False


class Max(MinMaxBase, Application):
    """
    Return, if possible, the maximum value of the list.

    When number of arguments is equal one, then
    return this argument.

    When number of arguments is equal two, then
    return, if possible, the value from (a, b) that is $\\ge$ the other.

    In common case, when the length of list greater than 2, the task
    is more complicated. Return only the arguments, which are greater
    than others, if it is possible to determine directional relation.

    If is not possible to determine such a relation, return a partially
    evaluated result.

    Assumptions are used to make the decision too.

    Also, only comparable arguments are permitted.

    It is named ``Max`` and not ``max`` to avoid conflicts
    with the built-in function ``max``.

    See Also
    ========

    Min : find minimum values
    """
    _category = (frozenset(('Max',)), ())
    _symox_head_name = 'Max'


class Min(MinMaxBase, Application):
    """
    Return, if possible, the minimum value of the list.

    It is named ``Min`` and not ``min`` to avoid conflicts
    with the built-in function ``min``.

    See Also
    ========

    Max : find maximum values
    """
    _category = (frozenset(('Min',)), ())
    _symox_head_name = 'Min'


def _install_lattice_elements():
    from symox import S
    for cls, zero, identity in ((Max, S.Infinity, S.NegativeInfinity),
                                (Min, S.NegativeInfinity, S.Infinity)):
        cls.zero = zero
        cls.identity = identity
        cls.__signature__ = _NATIVE[cls.__name__].__signature__


_install_lattice_elements()

__all__ += ['MinMaxBase', 'Max', 'Min']
