"""SymPy expression constructors backed by the native Rust expression tree.

The public classes describe expression kinds. Their constructors return native
``Expr`` objects and their instance checks inspect those objects' native kind,
so arithmetic never crosses into a second Python expression implementation.
"""

from fractions import Fraction

from . import _core


def _matrix_operand(value):
    return (getattr(value, '_symox_matrixexpr', False)
            or isinstance(value, _core.Matrix)
            or getattr(value, '_symox_sparse', False))


class _ExprType(type):
    def __instancecheck__(cls, value):
        return isinstance(value, _core.Expr) and bool(cls._matches(value))


class Number(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Number)


class Rational(Number):
    _matches = staticmethod(lambda e: e.is_Rational)

    def __new__(cls, p, q=None):
        if q is None and isinstance(p, (str, float, Fraction)):
            value = Fraction(p)
            return _core.Rational(value.numerator, value.denominator)
        return _core.Rational(p, q)


class Integer(Rational):
    _matches = staticmethod(lambda e: e.is_Integer)

    def __new__(cls, value=0):
        if isinstance(value, (float, str, Fraction)):
            value = int(value)
        return _core.Integer(value)


class Float(Number):
    _matches = staticmethod(lambda e: e.is_Float)

    def __new__(cls, value=0, dps=None, precision=None):
        return _core.Float(value, dps, precision=precision)


class Symbol(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Symbol)

    def __new__(cls, name, **assumptions):
        return _core.Symbol(name, **assumptions)


class Dummy(Symbol):
    _matches = staticmethod(lambda e: e.is_Dummy)

    def __new__(cls, name=None, dummy_index=None, **assumptions):
        return _core.Dummy(name, dummy_index, **assumptions)


class Add(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Add)

    def __new__(cls, *args, evaluate=True):
        if any(_matrix_operand(arg) for arg in args):
            from ._matrixexpr import MatAdd
            return MatAdd(*args,evaluate=evaluate)
        return _core._node("Add", args, evaluate)

    @staticmethod
    def make_args(expr):
        expr = _core.sympify(expr)
        return expr.args if expr.is_Add else (expr,)


class Mul(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Mul)

    def __new__(cls, *args, evaluate=True):
        if any(_matrix_operand(arg) for arg in args):
            from ._matrixexpr import MatMul
            return MatMul(*args,evaluate=evaluate)
        return _core._node("Mul", args, evaluate)

    @staticmethod
    def make_args(expr):
        expr = _core.sympify(expr)
        return expr.args if expr.is_Mul else (expr,)


class Pow(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Pow)

    def __new__(cls, base, exp, evaluate=True):
        if _matrix_operand(base):
            from ._matrixexpr import MatPow
            return MatPow(base,exp,evaluate=evaluate)
        return _core._node("Pow", (base, exp), evaluate)


class Tuple(metaclass=_ExprType):
    """An immutable tuple of native symbolic objects."""
    _matches = staticmethod(lambda e: e.is_Tuple)

    def __new__(cls, *items):
        return _core._apply("Tuple", *(Tuple(*item) if isinstance(item, (tuple, list)) else item for item in items))


class UnevaluatedExpr(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.class_name == "UnevaluatedExpr")

    def __new__(cls, expr):
        return _core._node("UnevaluatedExpr", (expr,))


class _SingletonRegistry:
    Half = Rational(1, 2)
    Zero = Integer(0)
    One = Integer(1)
    NegativeOne = Integer(-1)

    def __call__(self, value):
        return _core.sympify(value)

    def __getattr__(self, name):
        if name == "IdentityFunction":
            from . import IdentityFunction
            return IdentityFunction()
        if name in ("EmptySet", "UniversalSet", "Reals", "Integers", "Naturals",
                    "Naturals0", "Complexes", "Rationals"):
            from . import __dict__ as namespace
            return namespace[name]
        names = {
            "Pi": "pi", "Exp1": "E", "ImaginaryUnit": "I",
            "Infinity": "oo", "ComplexInfinity": "zoo", "NaN": "nan",
            "EulerGamma": "EulerGamma", "true": "true", "false": "false",
        }
        if name == "NegativeInfinity":
            return -_core.constant("oo")
        try:
            return _core.constant(names[name])
        except KeyError:
            raise AttributeError(name) from None


S = _SingletonRegistry()


def count_ops(expr, visual=False):
    if visual:
        raise NotImplementedError("visual operation counts are not supported")
    return _core.sympify(expr).count_ops()


def default_sort_key(expr, order=None):
    return _core._default_sort_key(expr)


def _restore_atom(name, value):
    if name in ("Symbol", "Dummy"):
        return Symbol(value[0], **value[1])
    return _core.sympify(value)


def _restore_expr(name, args, evaluate=True):
    if name in ("Integral", "Sum"):
        from . import Integral, Sum
        return {"Integral":Integral,"Sum":Sum}[name](*args)
    if name in ("Add", "Mul", "Pow", "UnevaluatedExpr"):
        return _core._node(name, args, evaluate)
    relations = {
        "Equality": "==", "Unequality": "!=", "StrictLessThan": "<",
        "LessThan": "<=", "StrictGreaterThan": ">", "GreaterThan": ">=",
    }
    if name in relations:
        return _core.Rel(relations[name], *args)
    if name == "Derivative":
        return _core.Derivative(*args)
    if name == "Piecewise":
        return _core.Piecewise(*zip(args[::2], args[1::2]))
    return _core._apply(name, *args)
