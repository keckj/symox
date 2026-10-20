"""SymPy expression constructors backed by the native Rust expression tree.

The public classes describe expression kinds. Their constructors return native
``Expr`` objects and their instance checks inspect those objects' native kind,
so arithmetic never crosses into a second Python expression implementation.
Structured tensor operands dispatch to the indexed tensor layer; their scalar
coefficients still use these native expression nodes.
"""

from fractions import Fraction
from functools import lru_cache

from . import _core, _node_types


def _matrix_operand(value):
    return (getattr(value, '_symox_matrixexpr', False)
            or isinstance(value, _core.Matrix)
            or getattr(value, '_symox_sparse', False))


class _ExprType(type):
    def __getattr__(cls, name):
        # Attributes patched onto the class, and SymPy's own class-level
        # placeholders, live on the nominal node class of the same name.
        try:
            return _node_types.class_attribute(cls.__name__, name)
        except AttributeError:
            pass
        # Constructors return native Expr instances, so inherited unbound
        # methods must be accessible through the public constructor as well.
        value = getattr(_core.Expr, name)
        if callable(value):
            return value
        raise AttributeError(name)

    def fromiter(cls, args, **assumptions):
        return cls(*tuple(args), **assumptions)

    def __instancecheck__(cls, value):
        return isinstance(value, _core.Expr) and bool(cls._matches(value))

    def __subclasscheck__(cls, subclass):
        # Native nodes are instances of SymPy's nominal classes, so a class
        # dispatch keyed by this public class must accept their types.
        if type.__subclasscheck__(cls, subclass):
            return True
        nominal = _node_types.node_class(cls.__name__) if isinstance(subclass, _node_types.NodeType) else None
        return nominal is not None and issubclass(subclass, nominal)

    def __setattr__(cls, name, value):
        type.__setattr__(cls, name, value)
        _node_types.forward_setattr(cls, cls.__name__, name, value)

    def __delattr__(cls, name):
        type.__delattr__(cls, name)
        _node_types.forward_delattr(cls, cls.__name__, name)

    def __repr__(cls):
        if _node_types._facade_active():
            return _node_types.sympy_class_repr(cls) or type.__repr__(cls)
        return type.__repr__(cls)

    def class_key(cls):
        name = cls.__name__
        if issubclass(cls, Number):
            return (1, 0, 'Number')
        if name == 'NumberSymbol':
            return (2, 0, name)
        if issubclass(cls, Symbol):
            return (2, 0, name)
        if name in ('Mul', 'Add', 'Pow'):
            return (3, ('Mul', 'Add', 'Pow').index(name), name)
        return (5, 0, name)


class Number(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Number)


class NumberSymbol(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.class_name in (
        'Pi', 'Exp1', 'EulerGamma', 'Catalan', 'GoldenRatio', 'TribonacciConstant'))


class Rational(Number):
    _matches = staticmethod(lambda e: e.is_Rational)

    def __new__(cls, p, q=None):
        if q is None and isinstance(p, (str, float, Fraction)):
            value = Fraction(p)
            return _core.Rational(value.numerator, value.denominator)
        return _core.Rational(p, q)


class Integer(Rational):
    _matches = staticmethod(lambda e: e.is_Integer)

    def __new__(cls, i=0):
        return _core.Integer(int(i))


# Native nodes expose these logical numeric classes through __class__.
# Only genuine rational/integer values satisfy the corresponding Python ABC.
from numbers import Integral as _IntegralABC, Rational as _RationalABC
_RationalABC.register(Rational)
_IntegralABC.register(Integer)


class Float(Number):
    _matches = staticmethod(lambda e: e.is_Float)

    def __new__(cls, num=0, dps=None, precision=None):
        if dps is not None and precision is not None:
            raise ValueError('Supply decimal digits or binary precision, not both')
        infer = (dps is None and precision is None) or dps == '' or precision == ''
        exact_digits = dps == '' or precision == ''
        if (isinstance(num,int) and not isinstance(num,bool)) or getattr(num,'is_Integer',False):
            if infer:
                from decimal import Decimal
                digits=len(Decimal(int(num)).as_tuple().digits)
                dps=digits if exact_digits else max(15,digits)
                precision=None
            return _core.Float(num,dps,precision=precision)
        if isinstance(num, str):
            from decimal import Decimal, InvalidOperation
            source = num.strip().replace(' ', '_').lower()
            # Validate Python-style separators before removing them for the
            # native decimal parser. This conversion never supplies the value.
            try:
                float(source)
            except ValueError:
                raise ValueError('string-float not recognized: ' + num) from None
            try:
                decimal=Decimal(source)
            except InvalidOperation:
                if exact_digits:raise ValueError('string-float not recognized by Decimal: '+num) from None
                # Decimal's exponent range is narrower than the native binary
                # representation. The reference uses its default precision
                # when decimal-based inference is unavailable.
                return _core.Float(source.replace('_',''),15 if infer else dps,precision=precision)
            num = source.replace('_', '')
            if not decimal.is_finite():
                return _core.Float('nan' if decimal.is_nan() else '-inf' if decimal.is_signed() else 'inf')
            if infer and decimal.is_finite():
                digits = len(decimal.as_tuple().digits)
                if '.' not in source:
                    digits = max(digits, decimal.adjusted() + 1)
                dps = digits if exact_digits else max(15, digits)
                precision = None
        elif exact_digits:
            raise ValueError('Automatic exact precision requires a string or integer')
        elif (not isinstance(num, (float, tuple, _core.Expr)) and not hasattr(num, '_mpf_')
              and hasattr(num, '_mpmath_')):
            # SymPy converts foreign numbers through mpmath, never through
            # their _sympy_ hook (which may itself call Float).
            import mpmath
            from mpmath.libmp import dps_to_prec
            bits = int(precision) if precision is not None else dps_to_prec(15 if dps is None else int(dps))
            return _core.Float(mpmath.mpf(num, prec=bits), precision=bits)
        return _core.Float(num, dps, precision=precision)


@lru_cache(maxsize=1000, typed=True)
def _cached_symbol(cls, name, assumptions):
    return _core.Symbol(name, **dict(assumptions))


class Symbol(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Symbol)

    def __new__(cls, name, **assumptions):
        commutative = assumptions.get('commutative', True)
        if commutative is None or commutative not in (True, False):
            raise ValueError('Symbol commutativity must be True or False')
        given = tuple((key, bool(value)) for key,value in assumptions.items() if value is not None)
        try:
            hash((cls, name, given))
        except TypeError:
            return _core.Symbol(name, **dict(given))
        return _cached_symbol(cls, name, given)


class Dummy(Symbol):
    _matches = staticmethod(lambda e: e.is_Dummy)

    def __new__(cls, name=None, dummy_index=None, **assumptions):
        return _core.Dummy(name, dummy_index, **assumptions)


class Wild(Symbol):
    _matches = staticmethod(lambda e: e.is_Wild)

    def __new__(cls, name, exclude=(), properties=(), **assumptions):
        return _core.Wild(name,exclude,properties,**assumptions)


def _evaluation_state(value=None):
    from ._public_core import _evaluation
    return _evaluation.get() if value is None else value


class Add(metaclass=_ExprType):
    """Construct a symbolic sum, optionally preserving its unevaluated terms."""
    _matches = staticmethod(lambda e: e.is_Add)

    def __new__(cls, *args, evaluate=None):
        evaluate = _evaluation_state(evaluate)
        if any(getattr(arg,'_symox_tensexpr',False) for arg in args):
            from .tensor.tensor import TensAdd
            return TensAdd(*args,evaluate=evaluate)
        if any(_matrix_operand(arg) for arg in args):
            from ._matrixexpr import MatAdd
            return MatAdd(*args,evaluate=evaluate)
        return _core._node("Add", args, evaluate)

    @staticmethod
    def make_args(expr):
        if getattr(expr,'_symox_tensexpr',False):return (expr,)
        expr = _core.sympify(expr)
        return expr.args if expr.is_Add else (expr,)


class Mul(metaclass=_ExprType):
    """Construct a symbolic product with exact native scalar arithmetic."""
    _matches = staticmethod(lambda e: e.is_Mul)

    def __new__(cls, *args, evaluate=None):
        evaluate = _evaluation_state(evaluate)
        if any(getattr(arg,'_symox_tensexpr',False) for arg in args):
            from .tensor.tensor import TensMul
            return TensMul(*args,evaluate=evaluate)
        if any(_matrix_operand(arg) for arg in args):
            from ._matrixexpr import MatMul
            return MatMul(*args,evaluate=evaluate)
        return _core._node("Mul", args, evaluate)

    @staticmethod
    def make_args(expr):
        if getattr(expr,'_symox_tensexpr',False):return (expr,)
        expr = _core.sympify(expr)
        return expr.args if expr.is_Mul else (expr,)


class Pow(metaclass=_ExprType):
    """Construct a symbolic power, preserving the requested evaluation mode."""
    _matches = staticmethod(lambda e: e.is_Pow)

    def __new__(cls, b, e, evaluate=None):
        evaluate = _evaluation_state(evaluate)
        if _matrix_operand(b):
            from ._matrixexpr import MatPow
            return MatPow(b,e,evaluate=evaluate)
        return _core._node("Pow", (b, e), evaluate)


class Tuple(_core.Expr, metaclass=_ExprType):
    """An immutable tuple of native symbolic objects."""
    _matches = staticmethod(lambda e: e.is_Tuple)

    @classmethod
    def fromiter(cls, args, **assumptions):
        return cls(*tuple(args), **assumptions)

    @classmethod
    def class_key(cls):
        return (5, 0, cls.__name__)

    def __new__(cls, *args, **kwargs):
        if not kwargs.get('sympify', True) and any(not isinstance(item, _core.Expr) for item in args):
            raise NotImplementedError('Tuple(sympify=False) with non-symbolic items is not implemented')
        native = _core._apply("Tuple", *(Tuple(*item) if isinstance(item, (tuple, list)) else item for item in args))
        return cls._from_native(native)

    @classmethod
    def _from_native(cls, native):
        if type(native) is cls:
            return native
        return _core.Expr.__new__(cls, native)

    def subs(self, *args, **kwargs):
        return _tuple_result(super().subs(*args, **kwargs))

    def xreplace(self, rule):
        return _tuple_result(super().xreplace(rule))


def _tuple_result(value):
    return Tuple._from_native(value) if isinstance(value,_core.Expr) and value.is_Tuple else value


class UnevaluatedExpr(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.class_name == "UnevaluatedExpr")

    def __new__(cls, arg, **kwargs):
        from . import sympify
        arg=sympify(arg,strict=True)
        if kwargs:raise TypeError('UnevaluatedExpr does not accept expression options')
        return _core._node("UnevaluatedExpr", (arg,))


class Atom(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Atom)


class AtomicExpr(Atom):
    _matches = staticmethod(lambda e: e.is_Atom and e.class_name not in ('BooleanTrue', 'BooleanFalse'))


_REL_OPS = {
    "Equality": "==", "Unequality": "!=", "StrictLessThan": "<",
    "LessThan": "<=", "StrictGreaterThan": ">", "GreaterThan": ">=",
}
_REL_ALIASES = {
    "==": "==", "eq": "==", "!=": "!=", "<>": "!=", "ne": "!=",
    ">=": ">=", "ge": ">=", "<=": "<=", "le": "<=",
    ">": ">", "gt": ">", "<": "<", "lt": "<",
}


class Relational(metaclass=_ExprType):
    """``Relational(lhs, rhs, rop)``; instance checks match any relation."""
    _matches = staticmethod(lambda e: e.class_name in _REL_OPS)

    def __new__(cls, lhs, rhs=0, rop=None, **kwargs):
        op = _REL_OPS.get(cls.__name__)
        if op is None:
            try:
                op = _REL_ALIASES[rop if rop is not None else "=="]
            except (KeyError, TypeError):
                raise ValueError("Invalid relational operator symbol: %r" % (rop,)) from None
        return _core.Rel(op, lhs, rhs) if _evaluation_state(kwargs.get("evaluate")) else _core.Rel(op, lhs, rhs, evaluate=False)


Rel = Relational


def _relation_class(name, doc):
    def __new__(cls, lhs, rhs=0, **kwargs):
        return _core.Rel(_REL_OPS[name], lhs, rhs) if _evaluation_state(kwargs.get("evaluate")) else _core.Rel(_REL_OPS[name], lhs, rhs, evaluate=False)
    return _ExprType(name, (Relational,), {
        "__new__": __new__, "__doc__": doc, "__module__": __name__,
        "_matches": staticmethod(lambda e: e.class_name == name),
    })


Equality = _relation_class("Equality", "``lhs == rhs`` (``Eq``).")
Unequality = _relation_class("Unequality", "``lhs != rhs`` (``Ne``).")
GreaterThan = _relation_class("GreaterThan", "``lhs >= rhs`` (``Ge``).")
LessThan = _relation_class("LessThan", "``lhs <= rhs`` (``Le``).")
StrictGreaterThan = _relation_class("StrictGreaterThan", "``lhs > rhs`` (``Gt``).")
StrictLessThan = _relation_class("StrictLessThan", "``lhs < rhs`` (``Lt``).")


def _install_relation_classes(namespace):
    """Share constructor and instance-check protocols between relation aliases."""
    import inspect
    def constructor(factory):
        def __new__(cls, *args, **options):
            return factory(*args, **options)
        return staticmethod(__new__)
    for short, long in [('Eq', 'Equality'), ('Ne', 'Unequality'),
                        ('Ge', 'GreaterThan'), ('Le', 'LessThan'),
                        ('Gt', 'StrictGreaterThan'), ('Lt', 'StrictLessThan')]:
        factory, cls = namespace[short], namespace[long]
        cls.__new__ = constructor(factory)
        cls.__signature__ = inspect.signature(factory)
        cls.__doc__ = factory.__doc__
        namespace[short] = cls


class _SingletonRegistry:
    Half = Rational(1, 2)
    Zero = Integer(0)
    One = Integer(1)
    NegativeOne = Integer(-1)

    def __call__(self,a,locals=None,convert_xor=True,strict=False,rational=False,evaluate=None):
        from . import sympify
        if locals is None and convert_xor and not strict and not rational and evaluate is None:
            return sympify(a)
        return sympify(a,locals=locals,convert_xor=convert_xor,strict=strict,
                       rational=rational,evaluate=evaluate)

    def __getattr__(self, name):
        if name == 'EmptySequence':
            from . import EmptySequence
            return EmptySequence
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
            "Catalan": "Catalan", "GoldenRatio": "GoldenRatio",
            "TribonacciConstant": "TribonacciConstant",
        }
        if name == "NegativeInfinity":
            value = -_core.constant("oo")
            setattr(self, name, value)
            return value
        try:
            from . import __dict__ as namespace
            public = names[name]
            value = namespace.get(public)
            if value is None:
                value = _core.constant(public)
            setattr(self, name, value)
            return value
        except KeyError:
            raise AttributeError(name) from None


S = _SingletonRegistry()


def count_ops(expr, visual=False):
    from ._operation_counts import count_ops as implementation
    return implementation(expr,visual=visual)


def default_sort_key(item, order=None):
    from . import Set
    if isinstance(item,Set):return item._e.sort_key(order=order)
    if isinstance(item,_core.AlgebraicNumber):
        return item.sort_key(order=order)
    if isinstance(item,(_core.Expr,_core._MatrixExpr)):
        return item.sort_key(order=order)
    if callable(getattr(item,'sort_key',None)):
        return item.sort_key(order=order)
    if not isinstance(item,str):
        try:
            children=item.items() if isinstance(item,dict) else iter(item)
        except TypeError:
            try:
                return _core.sympify(item).sort_key(order=order)
            except (TypeError,ValueError):
                pass
        else:
            args=[default_sort_key(child,order) for child in children]
            if isinstance(item,(dict,set,frozenset)):
                args.sort()
            return ((10,0,type(item).__name__),(len(args),tuple(args)),S.One.sort_key(),S.One)
    return ((0,0,type(item).__name__),(1,(str(item),)),S.One.sort_key(),S.One)


def _restore_atom(name, value):
    if name in ("Symbol", "Dummy"):
        return Symbol(value[0], **value[1])
    return _core.sympify(value)


def _restore_expr(name, args, evaluate=True):
    if name in ("Integral", "Sum", "Product"):
        from . import Integral, Sum, Product
        return {"Integral":Integral,"Sum":Sum,"Product":Product}[name](*args)
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
        if args and all(isinstance(arg,(tuple,list)) or getattr(arg,'is_Tuple',False) for arg in args):
            return _core.Piecewise(*args,evaluate=evaluate)
        return _core.Piecewise(*zip(args[::2], args[1::2]))
    return _core._apply(name, *args)
