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
from ._basic import _BasicMeta


def _matrix_operand(value):
    return (getattr(value, '_symox_matrixexpr', False)
            or isinstance(value, _core.Matrix)
            or getattr(value, '_symox_sparse', False))


def _is_facade(cls):
    # The public Expr class (see _expr_facade).
    return '_symox_expr_facade' in type.__dict__['__dict__'].__get__(cls)


def _storage_bases(metaclass, name, bases, namespace, kwargs):
    """Replace the public Expr class by the native storage type in the bases
    of a user class; return the class when another metaclass now applies."""
    if not any(isinstance(base, _ExprType) and _is_facade(base) for base in bases):
        return bases, None
    bases = tuple(_core.Expr if isinstance(base, _ExprType) and _is_facade(base) else base
                  for base in bases)
    if metaclass is not _ExprType:
        return bases, None
    winner = type
    for base in bases:
        candidate = type(base)
        if issubclass(candidate, winner):
            winner = candidate
        elif not issubclass(winner, candidate):
            raise TypeError('metaclass conflict: the metaclass of a derived class must be '
                            'a (non-strict) subclass of the metaclasses of all its bases')
    if issubclass(winner, _ExprType):
        return bases, None
    return bases, winner(name, bases, namespace, **kwargs)


_FACADE_OWN = frozenset(('__doc__', '__module__', '__qualname__', '__name__', '__slots__',
                         '__annotations__', '__abstractmethods__', '__init_subclass__',
                         '__class_getitem__', '__dict__', '__weakref__'))


def _facade_target(name):
    # SymPy patches of Expr reach every expression class, so they are
    # stored on the nominal Expr class, ahead of the native methods.
    if name in _FACADE_OWN or name.startswith('_symox') or name.startswith('_abc_'):
        return None
    return _node_types.node_class('Expr')


def _facade_setattr(name, value):
    target = _facade_target(name)
    if target is not None:
        type.__setattr__(target, name, value)


def _facade_delattr(name):
    target = _facade_target(name)
    if target is not None and name in target.__dict__:
        type.__delattr__(target, name)


def _nominal_subclass(cls, subclass):
    # Whether SymPy's class standing for ``subclass`` (its nominal class, or
    # that of the symox class it stands for or derives from) derives from
    # the SymPy class this symox class stands for.
    if not _node_types.is_symox_class(cls):
        return False
    if isinstance(subclass, _node_types.NodeType):
        sub = subclass
    else:
        sub = _node_types.nominal_anchor(subclass)
        if sub is None:
            return False
    nominal = _node_types.node_class(cls.__name__)
    return nominal is not None and issubclass(sub, nominal)


class _ExprType(_BasicMeta):
    def __getattr__(cls, name):
        # Attributes patched onto the class, and SymPy's own class-level
        # placeholders, live on the nominal node class of the same name.
        try:
            return _node_types.anchored_class_attribute(cls, name)
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

    def __new__(metaclass, name, bases, namespace, **kwargs):
        # User subclasses of Symbol and Dummy are real native subclasses
        # whose instances keep their class and attributes (see _host).
        # User subclasses of Atom, AtomicExpr, Number and NumberSymbol are
        # native expression classes.
        module = namespace.get('__module__', '')
        if module != 'symox' and not module.startswith('symox.'):
            bases, created = _storage_bases(metaclass, name, bases, namespace, kwargs)
            if created is not None:
                return created
        if module != __name__ and _host_atom_bases(bases):
            from ._host import atom_class
            return atom_class(metaclass, name, bases, namespace, **kwargs)
        if module != 'symox' and not module.startswith('symox.'):
            from ._host import expression_bases
            bases = expression_bases(bases)
            return super().__new__(metaclass, name, bases, namespace, **kwargs)
        cls = super().__new__(metaclass, name, bases, namespace, **kwargs)
        # The public expression classes are Basic subclasses, as in SymPy.
        from ._basic import Basic
        Basic.register(cls)
        return cls

    def __instancecheck__(cls, value):
        namespace = cls.__dict__
        if '_symox_expr_facade' in namespace:
            from ._expr_facade import is_expr
            return is_expr(value)
        if '_symox_host' in namespace or '_symox_head' in namespace:
            # User classes are nominal: native nodes are never instances.
            return type.__instancecheck__(cls, value)
        if not isinstance(value, _core.Expr):
            return False
        kind = type(value)
        if kind is not _core.Expr and type.__subclasscheck__(cls, kind):
            # Instances of user expression classes (see _host).
            return True
        if kind is not _core.Expr and not isinstance(kind, _node_types.NodeType) \
                and _nominal_subclass(cls, kind):
            # Instances of user classes derived from a symox class that
            # SymPy derives from this one (a user subclass of ``sin`` is a
            # TrigonometricFunction).
            return True
        try:
            matches = cls._matches
        except AttributeError:
            matches = None
        if matches is None:
            # Classes without a structural test (SymPy's abstract classes
            # such as TrigonometricFunction) follow the nominal classes.
            return isinstance(kind, _node_types.NodeType) and _nominal_subclass(cls, kind)
        return bool(matches(value))

    def __subclasscheck__(cls, subclass):
        # Native nodes are instances of SymPy's nominal classes, so a class
        # dispatch keyed by this public class must accept their types.
        if type.__subclasscheck__(cls, subclass):
            return True
        if '_symox_expr_facade' in cls.__dict__:
            from ._expr_facade import is_expr_class
            return isinstance(subclass, type) and is_expr_class(subclass)
        return isinstance(subclass, type) and _nominal_subclass(cls, subclass)

    def __setattr__(cls, name, value):
        type.__setattr__(cls, name, value)
        if '_symox_expr_facade' in cls.__dict__:
            _facade_setattr(name, value)
        else:
            _node_types.forward_setattr(cls, cls.__name__, name, value)

    def __delattr__(cls, name):
        type.__delattr__(cls, name)
        if '_symox_expr_facade' in cls.__dict__:
            _facade_delattr(name)
        else:
            _node_types.forward_delattr(cls, cls.__name__, name)

    def __repr__(cls):
        if _node_types._facade_active() and _node_types.is_symox_class(cls):
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

    def __new__(cls, *obj):
        # SymPy's Number constructor dispatches on the Python value type.
        import decimal
        if len(obj) == 1:
            obj = obj[0]
        if isinstance(obj, Number):
            return obj
        if isinstance(obj, int) or type(obj).__name__ in ('mpz', 'MPZ'):
            return Integer(obj)
        if isinstance(obj, tuple) and len(obj) == 2:
            return Rational(*obj)
        if isinstance(obj, (float, decimal.Decimal)) or type(obj).__name__ == 'mpf':
            return Float(obj)
        if isinstance(obj, str):
            from . import S, sympify
            lowered = obj.lower()
            if lowered == 'nan':
                return S.NaN
            if lowered in ('inf', '+inf'):
                return S.Infinity
            if lowered == '-inf':
                return S.NegativeInfinity
            value = sympify(obj)
            if isinstance(value, Number):
                return value
            raise ValueError('String "%s" does not denote a Number' % obj)
        raise TypeError('expected str|int|long|float|Decimal|Number object but got %r' % obj.__class__.__name__)


class NumberSymbol(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.class_name in (
        'Pi', 'Exp1', 'EulerGamma', 'Catalan', 'GoldenRatio', 'TribonacciConstant'))


def _float_exact_ratio(value):
    if not getattr(value, 'is_Float', False):
        return value
    sign, man, exp, _ = value._mpf_
    if not man and exp:
        return value
    return Fraction(-man if sign else man) * (Fraction(2) ** exp)


class Rational(Number):
    _matches = staticmethod(lambda e: e.is_Rational)

    def __new__(cls, p, q=None):
        # The reference converts a Float to its exact binary ratio, in both
        # the numerator and the denominator.
        p, q = _float_exact_ratio(p), _float_exact_ratio(q)
        if (isinstance(p, Fraction) or isinstance(q, Fraction)) and isinstance(p, (int, Fraction)) \
                and isinstance(q, (int, Fraction)) and q != 0:
            value = Fraction(p) / Fraction(q)
            return _core.Rational(value.numerator, value.denominator)
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
        try:
            return _core.Float(num, dps, precision=precision)
        except NotImplementedError:
            # An object whose conversion hook declines (a geometric locus,
            # for instance) is not a number: the reference fails in mpmath
            # with a TypeError, which callers such as optics.deviation catch.
            if isinstance(num, (str, int, float, tuple, _core.Expr)):
                raise
            raise TypeError('cannot create mpf from %r' % (num,)) from None


def _host_atom_bases(bases):
    symbol = globals().get('Symbol')
    return symbol is not None and any(
        isinstance(base, type) and issubclass(base, symbol) for base in bases)


def _symbol_instance(cls, symbol):
    if '_symox_host' not in cls.__dict__:
        return symbol
    from ._host import host_atom
    return host_atom(cls, symbol)


@lru_cache(maxsize=1000, typed=True)
def _cached_symbol(cls, name, assumptions):
    return _symbol_instance(cls, _core.Symbol(name, **dict(assumptions)))


def _symbol_assumptions(assumptions):
    commutative = assumptions.get('commutative', True)
    if commutative is None or commutative not in (True, False):
        raise ValueError('Symbol commutativity must be True or False')
    return tuple((key, bool(value)) for key,value in assumptions.items() if value is not None)


class Symbol(metaclass=_ExprType):
    _matches = staticmethod(lambda e: e.is_Symbol)

    def __new__(cls, name, **assumptions):
        given = _symbol_assumptions(assumptions)
        try:
            hash((cls, name, given))
        except TypeError:
            return _symbol_instance(cls, _core.Symbol(name, **dict(given)))
        return _cached_symbol(cls, name, given)

    @staticmethod
    def __xnew__(cls, name, **assumptions):
        """SymPy's uncached construction of a symbol of class ``cls``."""
        if not isinstance(name, str):
            raise TypeError('name should be a string, not %s' % repr(type(name)))
        return _symbol_instance(cls, _core.Symbol(name, **dict(_symbol_assumptions(assumptions))))

    @staticmethod
    def __xnew_cached_(cls, name, **assumptions):
        """SymPy's cached construction of a symbol of class ``cls``."""
        return Symbol.__new__(cls, name, **assumptions)


class Dummy(Symbol):
    _matches = staticmethod(lambda e: e.is_Dummy)

    def __new__(cls, name=None, dummy_index=None, **assumptions):
        return _symbol_instance(cls, _core.Dummy(name, dummy_index, **assumptions))


class Wild(Symbol):
    _matches = staticmethod(lambda e: e.is_Wild)

    def __new__(cls, name, exclude=(), properties=(), **assumptions):
        # Function classes in exclude (Wild('a', exclude=[sin])) reject any
        # expression containing them, as expr.has(sin) does in the reference.
        expressions, heads = [], []
        for item in exclude:
            (heads if _is_exclusion_head(item) else expressions).append(item)
        if heads:
            properties = tuple(properties) + (_head_exclusion(tuple(heads)),)
        return _symbol_instance(cls, _core.Wild(name,expressions,properties,**assumptions))


def _is_exclusion_head(item):
    if isinstance(item, _core.Expr):
        return False
    try:
        _core.sympify(item)
    except Exception:
        return callable(item)
    return False


@lru_cache(maxsize=None)
def _head_exclusion(heads):
    def excluded(expr):
        return not expr.has(*heads)
    excluded.__name__ = 'excludes_' + '_'.join(getattr(h, '__name__', 'head') for h in heads)
    return excluded


# Set once a Number class with its own arithmetic exists (see _host).
_arithmetic_numbers = [False]


def _evaluation_state(value=None):
    from ._public_core import _evaluation
    return _evaluation.get() if value is None else value


class Add(metaclass=_ExprType):
    """Construct a symbolic sum, optionally preserving its unevaluated terms."""
    _matches = staticmethod(lambda e: e.is_Add)
    identity = _core.Integer(0)
    _op_priority = 10.0

    def __new__(cls, *args, evaluate=None):
        evaluate = _evaluation_state(evaluate)
        if cls is Add and evaluate is True:
            from ._add_memo import add
            cached=add(args)
            if cached is not None:return cached
        if any(getattr(arg,'_symox_tensexpr',False) for arg in args):
            from .tensor.tensor import TensAdd
            return TensAdd(*args,evaluate=evaluate)
        if any(_matrix_operand(arg) for arg in args):
            from ._matrixexpr import MatAdd
            return MatAdd(*args,evaluate=evaluate)
        if _arithmetic_numbers[0] and evaluate and any(getattr(arg,'_symox_arith',False) for arg in args):
            # Number classes with their own arithmetic fold the coefficient.
            from ._host import fold_numbers
            coeff, terms = fold_numbers("Add", args)
            if not terms:
                return coeff
            args = (coeff, *terms)
        return _core._node("Add", args, evaluate)

    @staticmethod
    def make_args(expr):
        if getattr(expr,'_symox_tensexpr',False):return (expr,)
        expr = _core.sympify(expr)
        return expr.args if expr.is_Add else (expr,)

    @classmethod
    def _from_args(cls, args, is_commutative=None):
        # Reference AssocOp._from_args: keep already processed args exactly.
        # Native nodes derive commutativity from their args.
        if len(args) == 0:
            return cls.identity
        if len(args) == 1:
            return args[0]
        return cls(*args, evaluate=False)


class Mul(metaclass=_ExprType):
    """Construct a symbolic product with exact native scalar arithmetic."""
    _matches = staticmethod(lambda e: e.is_Mul)
    identity = _core.Integer(1)
    _op_priority = 10.0

    def __new__(cls, *args, evaluate=None):
        evaluate = _evaluation_state(evaluate)
        if any(getattr(arg,'_symox_tensexpr',False) for arg in args):
            from .tensor.tensor import TensMul
            return TensMul(*args,evaluate=evaluate)
        if any(_matrix_operand(arg) for arg in args):
            from ._matrixexpr import MatMul
            return MatMul(*args,evaluate=evaluate)
        if _arithmetic_numbers[0] and evaluate and any(getattr(arg,'_symox_arith',False) for arg in args):
            # Number classes with their own arithmetic fold the coefficient.
            from ._host import fold_numbers
            coeff, terms = fold_numbers("Mul", args)
            if not terms:
                return coeff
            args = (coeff, *terms)
        return _core._node("Mul", args, evaluate)

    @staticmethod
    def make_args(expr):
        if getattr(expr,'_symox_tensexpr',False):return (expr,)
        expr = _core.sympify(expr)
        return expr.args if expr.is_Mul else (expr,)

    @classmethod
    def _from_args(cls, args, is_commutative=None):
        # Reference AssocOp._from_args: keep already processed args exactly.
        # Native nodes derive commutativity from their args.
        if len(args) == 0:
            return cls.identity
        if len(args) == 1:
            return args[0]
        return cls(*args, evaluate=False)


class Pow(metaclass=_ExprType):
    """Construct a symbolic power, preserving the requested evaluation mode."""
    _matches = staticmethod(lambda e: e.is_Pow)

    def __new__(cls, b, e, evaluate=None):
        evaluate = _evaluation_state(evaluate)
        if _matrix_operand(b):
            from ._matrixexpr import MatPow
            return MatPow(b,e,evaluate=evaluate)
        native = _core._node("Pow", (b, e), evaluate)
        if evaluate and '_symox_head' in type(b).__dict__ and hasattr(b, '_eval_power'):
            from ._host import eval_power
            return eval_power(b, e, native)
        return native


class Tuple(_core.Expr, metaclass=_ExprType):
    """An immutable tuple of native symbolic objects."""
    _matches = staticmethod(lambda e: e.is_Tuple)
    # The class flags of the reference Basic that the native classes take
    # from _relational_properties.
    is_AlgebraicNumber = is_Indexed = is_MatAdd = is_MatMul = is_Not = False
    is_NumberSymbol = is_Point = is_Poly = is_Vector = is_symbol = False

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

    def __new__(cls, *args):
        # SymPy's Basic.__new__ for user expression classes (see _host).
        if issubclass(cls, _core.Expr):
            return _core.Expr.__new__(cls, *args)
        return object.__new__(cls, *args)


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
            "BooleanTrue": "true", "BooleanFalse": "false",
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
    from ._function_subs import function_class_name
    name=function_class_name(item)
    if name is not None:
        # A function class sorts as a non-Basic object named by its
        # metaclass, as in the reference.
        return ((0,0,name),(1,(str(item),)),S.One.sort_key(),S.One)
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
