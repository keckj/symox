# Immutable code contract compiled for comparison only. Never executed.

class _BasicMeta:

    def __instancecheck__(cls, value):
        if cls is Basic and isinstance(value, _core.Matrix):
            return value._immutable
        if not super().__instancecheck__(value):
            return False
        # An undefined head ``Function('f')`` is a class in the reference,
        # not a Basic instance, although ``Function`` is a Basic subclass.
        return type(value) is not _undefined_head_type()


class ABCMeta:

    def __instancecheck__(cls, instance):
        """Override for isinstance(instance, cls)."""
        return _abc_instancecheck(cls, instance)


    def __subclasscheck__(cls, subclass):
        """Override for issubclass(subclass, cls)."""
        return _abc_subclasscheck(cls, subclass)


class _ExprType:

    def __subclasscheck__(cls, subclass):
        # Native nodes are instances of SymPy's nominal classes, so a class
        # dispatch keyed by this public class must accept their types.
        if type.__subclasscheck__(cls, subclass):
            return True
        if '_symox_expr_facade' in cls.__dict__:
            from ._expr_facade import is_expr_class
            return isinstance(subclass, type) and is_expr_class(subclass)
        return isinstance(subclass, type) and _nominal_subclass(cls, subclass)


class _OperationCategory:

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


class _BooleanType:

    def __subclasscheck__(cls, subclass):
        if super().__subclasscheck__(subclass):
            return True
        from symox._node_types import NodeType
        nominal = cls._native_counterpart() if isinstance(subclass, NodeType) else None
        return nominal is not None and issubclass(subclass, nominal)


class WeakKeyDictionary:

    def __getitem__(self, key):
        return self.data[ref(key)]


    def __setitem__(self, key, value):
        self.data[ref(key, self._remove)] = value


class DDM:

    def __init__(self, rowslist, shape, domain):
        if not (isinstance(rowslist, list) and all(type(row) is list for row in rowslist)):
            raise DMBadInputError("rowslist must be a list of lists")
        m, n = shape
        if len(rowslist) != m or any(len(row) != n for row in rowslist):
            raise DMBadInputError("Inconsistent row-list/shape")

        rowslist = _native_elements(rowslist, domain)
        super().__init__([i.copy() for i in rowslist])
        self.shape = (m, n)
        self.rows = m
        self.cols = n
        self.domain = domain


    def copy(self):
        copyrows = [row[:] for row in self]
        return DDM(copyrows, self.shape, self.domain)


    @classmethod
    def zeros(cls, shape, domain):
        z = domain.zero
        m, n = shape
        rowslist = [[z] * n for _ in range(m)]
        return DDM(rowslist, shape, domain)


    @classmethod
    def eye(cls, size, domain):
        if isinstance(size, tuple):
            m, n = size
        elif isinstance(size, int):
            m = n = size
        one = domain.one
        ddm = cls.zeros((m, n), domain)
        for i in range(min(m, n)):
            ddm[i][i] = one
        return ddm


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


def native_module(cls):
    """The module recorded in the class dictionary, ignoring the facade."""
    return type.__dict__['__module__'].__get__(cls)


def is_symox_class(cls):
    """Whether ``cls`` is defined by symox itself.

    Only symox's own public classes stand for the nominal class of their
    name: a user class reusing a SymPy name (Torch's ``Max``, Devito's
    ``Derivative``) must not see that node class's attributes, accept its
    node types or report its SymPy path.
    """
    module = native_module(cls) if isinstance(cls, type) else getattr(cls, '__module__', None)
    return isinstance(module, str) and (module == 'symox' or module.startswith('symox.'))


def nominal_anchor(cls):
    """The nominal class standing for ``cls`` in SymPy's MRO, or None.

    A symox class stands for the nominal class of its own name. A class
    defined outside symox (a user subclass of ``Symbol``, ``Function`` or
    ``Expr``) inherits, as in SymPy, the attributes of its nearest symox
    ancestor that has a SymPy counterpart: patches of ``Basic`` or
    ``Symbol`` reach it, those of an unrelated class of the same name do not.
    """
    try:
        return _anchors[cls]
    except (KeyError, TypeError):
        pass
    anchor = None
    if is_symox_class(cls):
        anchor = node_class(cls.__name__)
    else:
        for base in cls.__mro__:
            if base is object or isinstance(base, NodeType) or not is_symox_class(base):
                continue
            name = base.__dict__.get('_symox_nominal', base.__name__)
            if name in CLASSES:
                anchor = node_class(name)
                break
    try:
        _anchors[cls] = anchor
    except TypeError:
        pass
    return anchor


def node_class(name):
    """Return the nominal class of a SymPy class name, or None."""
    cls = _classes.get(name)
    if cls is not None or name not in CLASSES:
        return cls
    module, base_names = CLASSES[name]
    # HostTree is the root: a tree containing user subclass instances keeps
    # them alive in its slot while reporting its SymPy class.
    bases = tuple(node_class(base) for base in base_names) or (HostTree,)
    namespace = {'__module__': 'symox'+module[len('sympy'):], '__qualname__': name,
                 '__slots__': ()}
    if name == 'Basic':
        namespace['_sage_'] = _sage_placeholder
    cls = NodeType(name, bases, namespace)
    _classes[name] = cls
    return cls


def _resolve(paths):
    from importlib import import_module
    classes = []
    for path in paths:
        module, name = path.split(':')
        classes.append(getattr(import_module(module), name))
    return tuple(classes)


def import_module(name, package=None):
    """Import a module.

    The 'package' argument is required when performing a relative import. It
    specifies the package to use as the anchor point from which to resolve the
    relative import to an absolute import.

    """
    level = 0
    if name.startswith('.'):
        if not package:
            msg = ("the 'package' argument is required to perform a relative "
                   "import for {!r}")
            raise TypeError(msg.format(name))
        for character in name:
            if character != '.':
                break
            level += 1
    return _bootstrap._gcd_import(name[level:], package, level)


def _call(self, *args):
    """``K(a)`` builds ``K.dtype(a)``: ``ZZ(2.5)`` truncates like ``int``
    and ``QQ(2.5)`` is the exact binary fraction, as ``PythonMPQ``."""
    if len(args) == 1 and _inexact_value(args[0]) is not None:
        a = args[0]
        if getattr(self, 'is_ZZ', False):
            return _native_call(self, int(a))
        if getattr(self, 'is_QQ', False):
            if not isinstance(a, float):
                raise TypeError('PythonMPQ() requires numeric or string argument')
            if isinstance(a, float):
                from fractions import Fraction
                f = Fraction(a)
                return _native_call(self, f.numerator, f.denominator)
        if getattr(self, 'is_RealField', False) and (isinstance(a, complex) or type(a).__name__ == 'mpc'):
            raise TypeError('cannot create mpf from %r' % (a,))
        if getattr(self, 'is_FractionField', False) and not self.dom.is_Field:
            # FracField.ground_new retries through the field of the ground:
            # ZZ(y)(1.5) is 3/2 and ZZ_I(y)(1.5) reports QQ_I.
            from ... import CoercionFailed
            try:
                return self.convert(a)
            except CoercionFailed:
                field = self.dom.get_field()
                return self.convert(field.to_sympy(_convert_inexact(field, _inexact_value(a))))
        return self.convert(a)
    if len(args) == 1 and isinstance(args[0], bool):
        return _native_call(self, int(args[0]))
    if len(args) == 1 and isinstance(args[0], list) and getattr(self, 'is_AlgebraicField', False):
        # The reference ``ANP(rep, mod, dom)``: dense coefficients of the
        # primitive element, highest power first.
        from ._algebraic_methods import from_dense_list
        return from_dense_list(self, args[0])
    return _native_call(self, *args)


def convert_from(self, element, base):
    if getattr(base, '_symox_monogenic_extension', False):
        from ... import CoercionFailed
        result = from_MonogenicFiniteExtension(self, element, base)
        if result is None:
            raise CoercionFailed('there is no canonical embedding of this extension element')
        return result
    return self.convert(element, base)


def convert(self, element, base=None):
    from ... import CoercionFailed, Expr, Basic, Set
    # EXRAW accepts every Expr before the coefficient check (``of_type``).
    raw = getattr(self, 'is_EXRAW', False) and base is None and isinstance(element, Expr)
    if not raw and _not_a_coeff(element):
        raise CoercionFailed('%s is not in any domain' % (element,))
    if base is None and isinstance(element, (set, frozenset, dict)) and (
            getattr(self, 'is_EX', False) or getattr(self, 'is_EXRAW', False)):
        # Not a sequence: strict sympify gives a FiniteSet or a Dict.
        from ... import sympify
        element = sympify(element, strict=True)
        if getattr(self, 'is_EXRAW', False):
            raise CoercionFailed('Expecting an Expr instance but found: %s' % type(element).__name__)
        return _native_convert(self, element, None)
    if base is None and isinstance(element, (set, frozenset, dict)) and _sympy_failure(self, None) is not None:
        from ... import sympify
        a = sympify(element, strict=True)
        raise _sympy_failure(self, a)
    if base is None and (element is None or isinstance(
            element, (str, bytes, bytearray, list, tuple, set, frozenset, dict))):
        # The reference reads only strict sympify input: no strings, and
        # no sequences.
        raise CoercionFailed("Cannot convert %s of type %s to %s" % (element, type(element), self))
    if getattr(base, '_symox_monogenic_extension', False):
        return convert_from(self, element, base)
    if base is None:
        if isinstance(element, bool):
            # ``isinstance(True, int)``: booleans convert as integers.
            element = int(element)
        else:
            value = _inexact_value(element)
            if value is not None:
                if value[1] in _INFINITIES:
                    raise CoercionFailed('%s is not in any domain' % (element,))
                result = _convert_inexact(self, value)
                if result is not None:
                    return result
        if getattr(self, 'is_EXRAW', False) and isinstance(element, (Basic, Set)) and not isinstance(element, Expr):
            raise CoercionFailed('Expecting an Expr instance but found: %s' % type(element).__name__)
        result = _float_from_sympy(self, element)
        if result is not None:
            return result
    try:
        return _native_convert(self, element, base)
    except CoercionFailed:
        from fractions import Fraction
        if base is None and isinstance(element, (Basic, Fraction)):
            if isinstance(element, Fraction):
                from ... import Rational
                element = Rational(element.numerator, element.denominator)
            error = _sympy_failure(self, element)
            if error is not None:
                raise error from None
        raise
    except TypeError:
        if base is not None or isinstance(element, Basic):
            raise
        raise CoercionFailed("Cannot convert %s of type %s to %s" % (element, type(element), self)) from None


def _not_a_coeff(element):
    # The reference ``_not_a_coeff``: NaN and the infinities.
    import math
    if isinstance(element, float):
        return math.isnan(element) or math.isinf(element)
    from ... import Basic, S
    return isinstance(element, Basic) and any(
        element is value or element == value
        for value in (S.NaN, S.Infinity, S.NegativeInfinity, S.ComplexInfinity))


def _native_elements(rows, domain):
    """Native symox domain elements of ``rows``.

    The reference stores the given objects as they are: its ground types
    accept Python ints and SymPy numbers in place of ``ZZ`` and ``QQ``
    elements. Native symox elements are a distinct type, so entries that are
    not elements of ``domain`` yet are converted here, exactly.
    """
    of_type, convert = domain.of_type, domain.convert
    try:
        if isinstance(rows, list):
            if all(of_type(e) for row in rows for e in row):
                return rows
            return [[e if of_type(e) else convert(e) for e in row] for row in rows]
        if isinstance(rows, dict):
            if all(of_type(e) for row in rows.values() for e in row.values()):
                return rows
            return {i: {j: e if of_type(e) else convert(e) for j, e in row.items()}
                    for i, row in rows.items()}
    except (TypeError, AttributeError):
        return rows
    return rows

