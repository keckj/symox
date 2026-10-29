"""User subclasses of Symbol, Dummy, Function, Expr, Atom and Basic.

Libraries built on SymPy subclass its atoms and functions to attach state and
behaviour: instance attributes, renamable names, ``eval`` and ``fdiff``
classmethods, printing hooks. Such a class is made a real subclass of the
native ``Expr`` type. Native nodes carry an opaque handle to the class, so
``args``, ``subs``, ``xreplace``, matrices and every other native rebuild
return instances of the subclass, and ``F.eval`` runs again for each rebuilt
application as SymPy's ``Function.__new__`` does.

Subclasses of ``Expr``, ``Atom``, ``AtomicExpr``, ``Number`` or ``Basic``
that are not functions are expression classes: their instances are native
applications of a head owned by the class, built by ``Expr.__new__(cls,
*args)`` without evaluation as SymPy's ``Basic.__new__``.

Equality and hashing follow SymPy: the class is part of the identity, next to
the native content (name and assumptions, or head and arguments) and, for
classes defining ``_hashable_content``, the content it returns.
"""
import copyreg
import inspect
import operator

from . import _core

# Hook flags understood by the native host handles.
_STR, _LATEX, _RENAMES, _EVAL, _FDIFF = 1, 2, 4, 8, 16
_EXPR, _ATOM, _DERIV, _HASHABLE = 32, 64, 128, 256
_SUBS, _ASSUME, _EVALF, _FUNC = 512, 1024, 2048, 4096
_NONCOMMUTATIVE, _UNCOMMUTED, _COMMUTATIVE_HOOK = 8192, 16384, 32768
_PRETTY = 65536
_NUMBER = 131072
_REPR = 262144

# Hooks dispatched by the Python facades, enabled once a class defines one.
_hooks = {'expand': False, 'free_symbols': False}


class HostTree(_core.Expr):
    """Expression object whose node contains host nodes. The nodes own their
    Python instances; the native type reports them to the collector."""
    __slots__ = ()

    def _hashable_content(self):
        return _native_hashable_content(self)


def _float_mpf(value):
    import mpmath
    from mpmath.libmp import prec_to_dps
    context = mpmath.mp.clone()
    context.prec = value._prec
    return context.mpf(str(value.evalf(prec_to_dps(value._prec) + 12)))._mpf_


def _native_hashable_content(self):
    """``_hashable_content`` of the SymPy class of a native node, which the
    overrides of user subclasses extend through ``super()``."""
    if self.is_Symbol:
        content = (self.name,) + tuple(sorted(self.assumptions0.items()))
        if self.is_Dummy:
            return content + (self.dummy_index,)
        if self.is_Wild:
            return content + (self.exclude, self.properties)
        return content
    if self.is_Float:
        return (_float_mpf(self), self._prec)
    if self.is_Rational:
        return (self.p, self.q)
    if self.class_name == 'Subs':
        from .core.sorting import ordered
        return (self.expr.xreplace(self.canonical_variables),) + tuple(ordered(
            [(v, p) for v, p in zip(self.variables, self.point) if not self.expr.has(v)]))
    return tuple(self.args)


_core._host_register(HostTree)

_PRINTER_KINDS = ('Symbol', 'Dummy', 'Wild')


def _is_user(cls):
    # The module recorded in the class dictionary: under the routing facade
    # the nominal node and set classes report SymPy's module, and the native
    # application classes of undefined heads record none.
    module = type.__dict__['__module__'].__get__(cls)
    if '_symox_user_hooks' in cls.__dict__:
        # symox classes written as host classes (noncommutative adjoint)
        return True
    if not isinstance(module, str):
        return False
    return not module.startswith('symox.') and module != 'symox'


def _user_classes(cls):
    """Classes of ``cls.__mro__`` that do not belong to symox itself."""
    return [base for base in cls.__mro__ if base is not object and _is_user(base)]


def _defines(cls, name):
    # Defaults installed by expression_class are not user definitions.
    return any(name in base.__dict__ and name not in base.__dict__.get('_symox_defaults', ())
               for base in _user_classes(cls))


def _printer_kind(cls):
    # SymPy printers dispatch on class names along the MRO, so a Dummy
    # subclass named Symbol (optlang) prints like a Symbol.
    for base in cls.__mro__:
        if base.__name__ in _PRINTER_KINDS:
            return base.__name__
    return 'Symbol'


def _print(obj, format):
    """Host rendering of ``obj`` for the native printers, or None."""
    if format == 'str':
        method = getattr(obj, '_sympystr', None)
        if method is not None:
            from .printing.str import StrPrinter
            return str(method(StrPrinter()))
        if obj.is_Symbol:
            kind = _printer_kind(type(obj))
            name = obj.name
            return '_' + name if kind == 'Dummy' else name + '_' if kind == 'Wild' else name
        return None
    if format == 'repr':
        method = getattr(obj, '_sympyrepr', None)
        if method is not None:
            from .printing.repr import ReprPrinter
            return str(method(ReprPrinter()))
        return None
    if format == 'latex':
        method = getattr(obj, '_latex', None)
        if method is not None:
            from .printing.latex import LatexPrinter
            return str(method(LatexPrinter()))
    if format in ('pretty', 'pretty_ascii'):
        # A one-line _pretty form of an atom or nullary application; any
        # other form keeps the native layout.
        method = getattr(obj, '_pretty', None)
        if method is not None:
            from .printing.pretty.pretty import PrettyPrinter
            printer = PrettyPrinter({'use_unicode': format == 'pretty'})
            form = method(printer)
            picture = getattr(form, 'picture', None)
            if picture is not None and len(picture) == 1:
                return picture[0]
    return None


def _reduce(self, protocol):
    # Basic.__reduce_ex__: an overridden __reduce__ wins, otherwise the
    # constructor arguments are pickled and the constructor runs again.
    if protocol < 2:
        raise NotImplementedError("Only pickle protocol 2 or higher is supported by SymPy")
    cls = type(self)
    if cls.__reduce__ is not _default_reduce:
        return self.__reduce__()
    args, kwargs = self.__getnewargs_ex__()
    state = self.__getstate__()
    if kwargs:
        return (copyreg.__newobj_ex__, (cls, tuple(args), dict(kwargs)), state, None, None)
    return (copyreg.__newobj__, (cls,) + tuple(args), state, None, None)


def _default_reduce(self):
    return self.__reduce_ex__(2)


def _setstate(self, state):
    for name, value in state.items():
        setattr(self, name, value)


def _anchored_getattr(self, name):
    # User classes are not subclasses of the nominal node classes, so the
    # attributes patched onto SymPy's classes after import (``Basic``,
    # ``Symbol``, ``Function``...) are looked up on their nominal anchor.
    from ._node_types import anchored_instance_attribute
    return anchored_instance_attribute(self, name)


class _HostAtom:
    """Pickling protocol of SymPy's Symbol for native host atoms."""
    __slots__ = ()

    __getattr__ = _anchored_getattr

    __reduce_ex__ = _reduce
    __reduce__ = _default_reduce
    __setstate__ = _setstate

    def __getnewargs_ex__(self):
        args = (self.name, self.dummy_index) if self.is_Dummy else (self.name,)
        return args, dict(self._assumptions_orig)

    def __getnewargs__(self):
        # Basic.__getnewargs__ returns the (empty) args; Wild its options.
        return (self.name, self.exclude, self.properties) if self.is_Wild else ()

    def __getstate__(self):
        return None


_defined = [False]


def host_classes_defined():
    """Whether any host class exists; plain expressions skip host paths."""
    return _defined[0]


def atom_class(metaclass, name, bases, namespace, **kwargs):
    """Create a user subclass of the Symbol or Dummy facade."""
    if not any(issubclass(base, _core.Expr) for base in bases if isinstance(base, type)):
        bases = tuple(bases) + (_HostAtom, HostTree)
    cls = type.__new__(metaclass, name, bases, namespace, **kwargs)
    flags = _hook_flags(cls)
    flags |= _commutative_flags(cls, True)
    if _defines(cls, '_sympystr') or _printer_kind(cls) != ('Dummy' if _is_dummy_class(cls) else 'Symbol'):
        flags |= _STR
    if _defines(cls, '_latex'):
        flags |= _LATEX
    if _defines(cls, '_pretty'):
        flags |= _PRETTY
    if _defines(cls, '_sympyrepr'):
        flags |= _REPR
    if _defines(cls, 'name'):
        flags |= _RENAMES
    _install_content_hash(cls)
    cls._symox_host = _core._host_class(cls, flags)
    _defined[0] = True
    return cls


def _install_content_hash(cls):
    # Basic.__hash__ hashes the class name and _hashable_content(); a class
    # that also defines __eq__ without __hash__ stays unhashable.
    if _defines(cls, '_hashable_content') and not _defines(cls, '__hash__'):
        cls.__hash__ = _content_hash


def _is_dummy_class(cls):
    from ._compat import Dummy
    return issubclass(cls, Dummy)


def _hashable_flag(cls):
    return _HASHABLE if _defines(cls, '_hashable_content') else 0


def _assumption_handlers(cls):
    """The ``_eval_is_<fact>`` handlers that user classes define."""
    from .core.assumptions import _assume_defined
    return {fact: getattr(cls, '_eval_is_' + fact) for fact in sorted(_assume_defined)
            if _defines(cls, '_eval_is_' + fact)}


def _hook_flags(cls):
    """Flags of the hooks native algorithms consult, as SymPy's own code
    dispatches them: ``_eval_subs`` in ``Basic._subs``, ``_eval_is_<fact>``
    in ``_ask``, ``_eval_evalf`` in ``evalf`` and an overridden ``func`` in
    the rebuilds of ``xreplace`` and ``subs``."""
    flags = _hashable_flag(cls)
    if _defines(cls, '_eval_subs'):
        flags |= _SUBS
    cls._symox_handlers = _assumption_handlers(cls)
    if cls._symox_handlers:
        flags |= _ASSUME
    if _defines(cls, '_eval_evalf'):
        flags |= _EVALF
    if _defines(cls, 'func'):
        flags |= _FUNC
    if any(name.startswith('_eval_expand_') and _defines(cls, name) for name in dir(cls)):
        _hooks['expand'] = True
    if _defines(cls, 'free_symbols'):
        _install_free_symbols()
    return flags


def _commutative_flags(cls, default):
    """How instances of ``cls`` decide SymPy's ``is_commutative``: an
    ``is_commutative`` class attribute (False or None), an
    ``_eval_is_commutative`` handler, or the default of the SymPy base
    class: True for symbols and numbers, None for ``Expr``, ``Atom`` and
    ``Basic`` subclasses (``default`` None), and the arguments of a
    ``Function`` (``default`` 'arguments')."""
    for base in _user_classes(cls):
        if 'is_commutative' in base.__dict__:
            value = base.__dict__['is_commutative']
            if value is None:
                return _UNCOMMUTED
            if isinstance(value, (bool, int)):
                return 0 if value else _NONCOMMUTATIVE
            break
    # Other class facts imply it as in SymPy's default assumptions:
    # ``is_extended_positive = True`` makes the class commutative.
    from .core.assumptions import StdFactKB, _assume_defined
    declared = {}
    for base in reversed(_user_classes(cls)):
        defaults = base.__dict__.get('_symox_defaults', ())
        for fact in _assume_defined:
            value = base.__dict__.get('is_' + fact)
            if isinstance(value, bool) and 'is_' + fact not in defaults:
                declared[fact] = value
    if declared:
        try:
            implied = StdFactKB(declared).get('commutative')
        except Exception:
            implied = None
        if implied is not None:
            return 0 if implied else _NONCOMMUTATIVE
    if 'commutative' in cls._symox_handlers:
        return _COMMUTATIVE_HOOK
    return _UNCOMMUTED if default is None else 0


def _seed_facts(obj):
    facts = dict(_class_facts(type(obj)))
    if getattr(obj, 'is_Symbol', False):
        facts.update(getattr(obj, '_assumptions_orig', None) or {})
    return facts


def ask(obj, fact):
    """SymPy's ``_ask(fact, obj)`` over the ``_eval_is_<fact>`` handlers of
    a host instance. Every value found is deduced into the facts of the
    instance; None leaves the query to the native facts of its node."""
    from .core.assumptions import StdFactKB, _assume_rules
    handlers = type(obj)._symox_handlers
    namespace = getattr(obj, '__dict__', None)
    facts = None if namespace is None else namespace.get('_symox_facts')
    if facts is None:
        facts = StdFactKB(_seed_facts(obj))
        if namespace is not None:
            namespace['_symox_facts'] = facts
    if fact in facts:
        return facts[fact]
    queue, queued = [fact], {fact}
    for current in queue:
        if current in facts:
            continue
        handler = handlers.get(current)
        if handler is not None:
            value = handler(obj)
            if value is not None:
                facts.deduce_all_facts(((current, bool(value)),))
        value = facts.get(fact)
        if value is not None:
            return value
        new = sorted(_assume_rules.prereq[current] - queued)
        queue.extend(new)
        queued.update(new)
    if fact in facts:
        return facts[fact]
    facts._tell(fact, None)
    return None


def host_atom(cls, atom):
    """The instance of host class ``cls`` for the native atom ``atom``."""
    host = cls.__dict__.get('_symox_host')
    return atom if host is None else _core._host_atom(host, atom)


# --- functions ---------------------------------------------------------------

def _arity(cls):
    """SymPy's FunctionClass arity: argument counts accepted by ``eval``."""
    evaluator = getattr(cls, 'eval', None)
    if evaluator is None:
        return None
    try:
        parameters = inspect.signature(evaluator).parameters.values()
    except (TypeError, ValueError):
        return None
    if any(p.kind == p.VAR_POSITIONAL for p in parameters):
        return None
    positional = [p for p in parameters if p.kind == p.POSITIONAL_OR_KEYWORD]
    required = sum(1 for p in positional if p.default is p.empty)
    optional = len(positional) - required
    if not optional:
        return (required,)
    return tuple(range(required, required + optional + 1))


def _canonical_nargs(nargs):
    import operator
    if nargs is None:
        return None
    try:
        counts = (operator.index(nargs),)
    except TypeError:
        counts = tuple(operator.index(n) for n in nargs)
    if not counts:
        raise ValueError('nargs must not be empty')
    if any(n < 0 for n in counts):
        raise ValueError('nargs must contain nonnegative integers')
    return tuple(sorted(set(counts)))


def _class_facts(cls):
    """Head assumptions from ``is_<fact> = True/False`` class attributes."""
    given = {}
    for base in reversed(_user_classes(cls)):
        for key, value in base.__dict__.items():
            if key.startswith('is_') and key[3:].islower() and isinstance(value, bool):
                given[key[3:]] = value
    facts = {}
    for key, value in given.items():
        try:
            _core._function_facts({key: value})
        except ValueError:
            continue
        facts[key] = value
    return facts


class _HostFunctionType(type):
    """Metaclass of user Function subclasses (SymPy's FunctionClass)."""

    def __new__(metaclass, name, bases, namespace, **kwargs):
        cls = super().__new__(metaclass, name, bases, namespace, **kwargs)
        if namespace.get('_symox_host_base'):
            return cls
        nargs = namespace.get('nargs', None)
        if nargs is not None and not isinstance(nargs, property):
            nargs = _canonical_nargs(nargs)
        elif any(base.__dict__.get('eval') is not None for base in _user_classes(cls)):
            nargs = _arity(cls)
        else:
            nargs = next((base.__dict__['_nargs'] for base in cls.__mro__[1:]
                          if '_nargs' in base.__dict__), None)
        cls._nargs = nargs
        flags = _hook_flags(cls)
        flags |= _commutative_flags(cls, 'arguments')
        if _defines(cls, 'eval') or _defines(cls, '__new__'):
            flags |= _EVAL
        if _defines(cls, 'fdiff'):
            flags |= _FDIFF
        if _defines(cls, '_eval_derivative'):
            flags |= _DERIV
        if _defines(cls, '_sympystr'):
            flags |= _STR
        if _defines(cls, '_latex'):
            flags |= _LATEX
        if _defines(cls, '_pretty'):
            flags |= _PRETTY
        if _defines(cls, '_sympyrepr'):
            flags |= _REPR
        _install_content_hash(cls)
        host = _core._host_class(cls, flags)
        cls._symox_head = _core._host_head(host, name, _class_facts(cls), nargs)
        _defined[0] = True
        return cls

    @property
    def nargs(cls):
        from . import FiniteSet, Naturals0
        return Naturals0 if cls._nargs is None else FiniteSet(*cls._nargs)

    def __repr__(cls):
        return cls.__name__


def _function_class_key(cls):
    return (4, 0 if cls._nargs is None else 10000, cls.__name__)


class _HostFunction(HostTree, metaclass=_HostFunctionType):
    """Base of user Function subclasses, mirroring SymPy's Function."""
    _symox_host_base = True
    _symox_nominal = 'Function'
    __getattr__ = _anchored_getattr
    _nargs = None

    def __new__(cls, *args, **options):
        evaluate = options.pop('evaluate', None)
        if options:
            raise TypeError('%s() got unexpected keyword arguments %s' % (cls.__name__, sorted(options)))
        n = len(args)
        if cls._nargs is not None and n not in cls._nargs:
            counts = cls._nargs
            if min(counts) == max(counts):
                qualifier, count = 'exactly', min(counts)
            elif n < min(counts):
                qualifier, count = 'at least', min(counts)
            else:
                qualifier, count = 'at most', max(counts)
            raise TypeError('%s takes %s %s argument%s (%s given)' % (
                cls.__name__, qualifier, count, '' if count == 1 else 's', n))
        from . import sympify
        from ._public_core import _evaluation
        values = tuple(sympify(arg) for arg in args)
        if evaluate is None:
            evaluate = _evaluation.get()
        if evaluate:
            value = cls.eval(*values)
            if value is not None:
                return sympify(value)
        return _core._host_apply(cls._symox_head, values)

    @classmethod
    def eval(cls, *args):
        return None

    @classmethod
    def class_key(cls):
        return _function_class_key(cls)

    @property
    def nargs(self):
        return type(self).nargs

    def fdiff(self, argindex=1):
        """SymPy's default partial derivative of an application."""
        from . import Derivative, Dummy, Subs
        from .core.function import ArgumentIndexError
        if not 1 <= argindex <= len(self.args):
            raise ArgumentIndexError(self, argindex)
        index = argindex - 1
        argument = self.args[index]
        if getattr(argument, '_diff_wrt', False):
            if len(self.args) == 1 or not argument.is_Symbol or not any(
                    argument in other.free_symbols for i, other in enumerate(self.args) if i != index):
                return Derivative(self, argument)
        dummy = Dummy('xi_%i' % argindex, dummy_index=hash(argument))
        arguments = list(self.args)
        arguments[index] = dummy
        return Subs(Derivative(self.func(*arguments), dummy), dummy, argument)

    __reduce_ex__ = _reduce
    __reduce__ = _default_reduce
    __setstate__ = _setstate

    def __getnewargs_ex__(self):
        return tuple(self.args), {}

    def __getnewargs__(self):
        return tuple(self.args)

    def __getstate__(self):
        return None


def rewrite(expr, pattern, rule, method, native, deep, hints):
    """SymPy's ``Basic._rewrite`` for trees containing host nodes.

    Host nodes dispatch to their ``_eval_rewrite_as_<rule>`` methods after
    their arguments were rewritten; other nodes keep the native rewrite.
    """
    contains, is_host = _core._host_contains(expr)
    if not contains:
        return native(expr, deep)
    args = tuple(rewrite(a, pattern, rule, method, native, deep, hints) for a in expr.args) if deep else expr.args
    if not is_host:
        rebuilt = expr if args == expr.args else expr.func(*args)
        return native(rebuilt, False)
    if not pattern or any(isinstance(expr, p) for p in pattern if isinstance(p, type)):
        evaluator = getattr(expr, method, None)
        if evaluator is not None:
            value = evaluator(*args, **hints)
        else:
            evaluator = getattr(expr, '_eval_rewrite', None)
            value = None if evaluator is None else evaluator(rule, args, **hints)
        if value is not None:
            return value
    if not args:
        return expr
    return expr.func(*args)


def function_class(name, bases, namespace, **kwargs):
    """Create ``class name(Function)`` as a native host function class."""
    from . import Function
    bases = tuple(_HostFunction if base is Function else base for base in bases)
    return _HostFunctionType(name, bases, namespace, **kwargs)


def is_host_function_class(value):
    return isinstance(value, _HostFunctionType) and not value.__dict__.get('_symox_host_base')


# --- expression classes -------------------------------------------------------

_abstract = {}


def _abstract_kinds():
    """Kinds of the symox classes an expression class may derive from."""
    if not _abstract:
        from ._basic import Basic
        from ._compat import Atom, AtomicExpr, Number, NumberSymbol
        _abstract.update({_core.Expr: 'expr', HostTree: 'expr', Basic: 'expr',
                          Atom: 'atom', AtomicExpr: 'atom', NumberSymbol: 'atom',
                          Number: 'number'})
    return _abstract


def expression_bases(bases):
    """Bases of a user class deriving from the Atom, AtomicExpr, Number or
    NumberSymbol facades, completed with the native expression storage."""
    if any(isinstance(base, type) and issubclass(base, _core.Expr) for base in bases):
        return bases
    kinds = _abstract_kinds()
    for base in bases:
        if isinstance(base, type) and any(kinds.get(b) in ('atom', 'number') for b in base.__mro__):
            return tuple(bases) + (HostTree,)
    return bases


def _expression_kind(cls):
    kinds = _abstract_kinds()
    for base in cls.__mro__:
        if base is not object and not _is_user(base):
            # Symox expression classes evaluated in Python (Max, Min) are
            # expression storage for their user subclasses.
            return kinds.get(base, getattr(base, '_symox_expression_kind', None))
    return None


def _atom_class_key(cls):
    return (2, 0, cls.__name__)


def _number_class_key(cls):
    return (1, 0, 'Number')


def _expr_class_key(cls):
    return (5, 0, cls.__name__)


def _expr_is_number(self):
    return all(arg.is_number for arg in self.args)


def _getnewargs(self):
    return tuple(self.args)


def _getnewargs_ex(self):
    return tuple(self.__getnewargs__()), {}


def _getstate(self):
    return None


def _hashable_args(self):
    return tuple(self.args)


def _content_hash(self):
    return hash((type(self).__name__,) + tuple(self._hashable_content()))


class _InstanceAttribute:
    """A read-only native attribute that instances of an expression class
    set, as SymPy objects do (Devito's ``self.name = name``)."""

    def __init__(self, name, native):
        self.name, self.native = name, native

    def __get__(self, obj, owner=None):
        if obj is None:
            return self
        try:
            return obj.__dict__[self.name]
        except (AttributeError, KeyError):
            return self.native.__get__(obj, owner)

    def __set__(self, obj, value):
        obj.__dict__[self.name] = value

    def __delete__(self, obj):
        try:
            del obj.__dict__[self.name]
        except (AttributeError, KeyError):
            raise AttributeError(self.name) from None


def _expression_setattr(self, name, value):
    native = _core.Expr.__dict__.get(name)
    if native is not None and hasattr(type(native), '__set__'):
        cls = type(self)
        for base in cls.__mro__:
            if name in base.__dict__:
                if base.__dict__[name] is native:
                    type.__setattr__(cls, name, _InstanceAttribute(name, native))
                break
    object.__setattr__(self, name, value)


def eval_power(base, exp, native):
    """Return the power of an expression class instance as the class's
    ``_eval_power`` evaluates it, as the reference ``Pow.__new__`` does once
    the generic exponent rules leave ``native`` unevaluated."""
    hook = getattr(type(base), '_eval_power', None)
    if hook is None or not getattr(native, 'is_Pow', False):
        return native
    native_base, native_exp = native.args
    if native_base is not base and native_base != base:
        return native
    result = hook(base, native_exp)
    return native if result is None else result


def _number_float(self):
    from mpmath.libmp import to_float
    return to_float(self._as_mpf_val(53))


def _number_eval_evalf(self, prec):
    # Reference Number._eval_evalf: Float._new(self._as_mpf_val(prec), prec),
    # which maps the special values to the singletons.
    import mpmath
    from mpmath.libmp import finf, fnan, fninf, fzero
    from . import Float, S
    value = self._as_mpf_val(prec)
    special = {fzero: S.Zero, fnan: S.NaN, finf: S.Infinity, fninf: S.NegativeInfinity}
    if value in special:
        return special[value]
    with mpmath.workprec(prec):
        return Float(mpmath.mpf(value), precision=prec)


def fold_numbers(op, args):
    """Fold the Number terms of an evaluated Add or Mul in the order of the
    reference ``AssocOp.flatten``, so that the arithmetic of a Number class
    combines them, and return the coefficient with the remaining terms."""
    combine = operator.add if op == 'Add' else operator.mul
    coeff = _core.Integer(0 if op == 'Add' else 1)
    nested = 'is_' + op
    seq = [_core.sympify(arg) for arg in args]
    terms = []
    # The reference appends the terms of nested operations to the sequence
    # it is iterating.
    for term in seq:
        if getattr(term, nested, False):
            seq.extend(term.args)
        elif term.is_Number:
            coeff = combine(coeff, term)
        else:
            terms.append(term)
    return coeff, terms


def _expression_pow(self, other, mod=None):
    if mod is not None:
        return _core.Expr.__pow__(self, other, mod)
    native = _core.Expr.__pow__(self, other)
    if native is NotImplemented:
        return native
    return eval_power(self, other, native)


_NUMBER_ARITHMETIC = ('__add__', '__radd__', '__mul__', '__rmul__')


# Methods of the formal Basic container that the native implementations
# replace for expression classes: their instances are native nodes.
_NATIVE_BASIC = ('args', 'func', 'free_symbols', '__repr__', '__str__', '__eq__', '__ne__',
                 '__hash__', 'srepr', 'xreplace', 'subs', 'has', 'atoms', 'doit',
                 'is_Atom', 'is_Function')


def expression_class(cls):
    """Register a user subclass of the native Expr type that is neither a
    Symbol nor a Function subclass as an expression class."""
    if not _is_user(cls) or isinstance(cls, _HostFunctionType) or issubclass(cls, _HostAtom):
        return
    kind = _expression_kind(cls)
    if kind is None:
        return
    from ._basic import Basic
    from ._compat import NumberSymbol
    atomic = kind != 'expr'
    number = kind == 'number' or issubclass(cls, NumberSymbol)
    defaults = {
        'is_Function': False,
        'is_Atom': atomic,
        'class_key': classmethod(_number_class_key if kind == 'number' else
                                 _atom_class_key if atomic else _expr_class_key),
        # Number and NumberSymbol subclasses are numbers (SymPy is_number).
        'is_number': True if number else False if atomic else property(_expr_is_number),
        '__reduce_ex__': _reduce,
        '__reduce__': _default_reduce,
        '__getnewargs__': _getnewargs,
        '__getnewargs_ex__': _getnewargs_ex,
        '__getstate__': _getstate,
        '__setstate__': _setstate,
        '_hashable_content': _hashable_args,
        '__setattr__': _expression_setattr,
        '__getattr__': _anchored_getattr,
    }
    if kind == 'number':
        defaults['is_Number'] = True
        if _defines(cls, '_as_mpf_val'):
            # Reference Number.__float__ and Number._eval_evalf.
            defaults['__float__'] = _number_float
            defaults['_eval_evalf'] = _number_eval_evalf
    if issubclass(cls, Basic):
        for name in _NATIVE_BASIC:
            if name not in defaults and name in _core.Expr.__dict__:
                defaults[name] = _core.Expr.__dict__[name]
    if _defines(cls, '_eval_power'):
        defaults['__pow__'] = _expression_pow
    if kind == 'number' and any(_defines(cls, name) for name in _NUMBER_ARITHMETIC):
        cls._symox_arith = True
        from ._compat import _arithmetic_numbers
        _arithmetic_numbers[0] = True
    hashable = _defines(cls, '_hashable_content')
    if hashable:
        defaults['__hash__'] = _content_hash
    installed = [name for name in defaults if not _defines(cls, name)]
    for name in installed:
        setattr(cls, name, defaults[name])
    cls._symox_defaults = frozenset(installed)
    flags = _EXPR | (_ATOM if atomic else 0) | (_NUMBER if number else 0) | _hook_flags(cls)
    if '_eval_evalf' in installed:
        flags |= _EVALF
    from ._compat import NumberSymbol
    flags |= _commutative_flags(cls, True if kind == 'number' or issubclass(cls, NumberSymbol) else None)
    if _defines(cls, '__new__') or _defines(cls, '__init__'):
        flags |= _EVAL
    if _defines(cls, '_eval_derivative'):
        flags |= _DERIV
    if _defines(cls, '_sympystr'):
        flags |= _STR
    if _defines(cls, '_latex'):
        flags |= _LATEX
    if _defines(cls, '_pretty'):
        flags |= _PRETTY
    if _defines(cls, '_sympyrepr'):
        flags |= _REPR
    key = cls.class_key()
    host = _core._host_class(cls, flags, (key[0], key[1]))
    cls._symox_head = _core._host_head(host, cls.__name__, _class_facts(cls), None)
    _defined[0] = True


def _init_subclass(cls, **kwargs):
    super(_core.Expr, cls).__init_subclass__(**kwargs)
    expression_class(cls)


_core.Expr.__init_subclass__ = classmethod(_init_subclass)


# --- hooks dispatched by the Python facades ---------------------------------

def _tree_free_symbols(self):
    # Basic.free_symbols: the union over the arguments, so that a host class
    # overriding free_symbols decides for its own nodes. Binding nodes keep
    # their native rules (bound variables).
    native = _core.Expr.free_symbols.__get__(self)
    if not _core._host_contains(self)[0] or self.is_Atom or not (
            self.is_Add or self.is_Mul or self.is_Pow or self.is_Function or self.is_Relational):
        return native
    return set().union(*[arg.free_symbols for arg in self.args])


def _install_free_symbols():
    if not _hooks['free_symbols']:
        _hooks['free_symbols'] = True
        HostTree.free_symbols = property(_tree_free_symbols)


def simplify_hook(expr):
    """The ``_eval_simplify`` method of a host instance whose class defines
    one, which SymPy's ``simplify`` returns at once, or None."""
    if not _defined[0] or not isinstance(expr, HostTree) or not _core._host_contains(expr)[1]:
        return None
    if not _defines(type(expr), '_eval_simplify'):
        return None
    return expr._eval_simplify


# Hints the native expansion implements, and SymPy's default hints.
_NATIVE_HINTS = frozenset(('basic', 'mul', 'multinomial', 'power_base', 'power_exp', 'log',
                           'trig', 'func', 'complex'))
_DEFAULT_HINTS = ('power_base', 'power_exp', 'mul', 'log', 'multinomial', 'basic')
_OPTIONS = frozenset(('force', 'max_terms', 'factor'))


def _expand_node(expr, hint, hints, native):
    # SymPy's ``getattr(expr, '_eval_expand_' + hint)(**hints)`` for one node.
    if _core._host_contains(expr)[1]:
        method = getattr(expr, '_eval_expand_' + hint, None)
        return None if method is None else method(**hints)
    if hint not in _NATIVE_HINTS:
        return None
    single = {name: False for name in _DEFAULT_HINTS}
    single[hint] = True
    for name in ('force', 'max_terms'):
        if name in hints:
            single[name] = hints[name]
    return native(expr, False, single)


def _expand_hint(expr, hint, deep, hints, native):
    # SymPy's Expr._expand_hint. Subtrees without host nodes expand natively.
    if not _core._host_contains(expr)[0]:
        if hint not in _NATIVE_HINTS:
            return expr, False
        single = {name: False for name in _DEFAULT_HINTS}
        single[hint] = True
        for name in ('force', 'max_terms'):
            if name in hints:
                single[name] = hints[name]
        result = native(expr, deep, single)
        return result, result != expr
    hit = False
    if deep and expr.args and not expr.is_Atom:
        args = []
        for arg in expr.args:
            arg, arg_hit = _expand_hint(arg, hint, True, hints, native)
            hit |= arg_hit
            args.append(arg)
        if hit:
            expr = expr.func(*args)
    result = _expand_node(expr, hint, hints, native)
    if result is not None and result != expr:
        return result, True
    return expr, hit


def expand(expr, deep, hints, native):
    """SymPy's ``Expr.expand`` for a tree containing host nodes, when a host
    class defines ``_eval_expand_<hint>`` methods; ``native(expr, deep,
    hints)`` expands nodes without host nodes. Returns None to leave the
    expansion to ``native``."""
    if not _hooks['expand'] or not isinstance(expr, _core.Expr) or not _core._host_contains(expr)[0]:
        return None
    hints = dict(hints)
    if any(hints.get(name) for name in ('frac', 'numer', 'denom')) or hints.get('modulus') is not None:
        return None
    hints.pop('modulus', None)
    for name in ('frac', 'numer', 'denom'):
        hints.pop(name, None)
    for name in _DEFAULT_HINTS:
        hints.setdefault(name, True)
    # Hints apply once each, multinomial before mul, then mul, multinomial
    # and log repeat until nothing changes.
    for hint in sorted(hints, key=lambda name: 'mulz' if name == 'mul' else name):
        if hints[hint] and hint not in _OPTIONS:
            expr, _ = _expand_hint(expr, hint, deep, hints, native)
    while True:
        was = expr
        for hint in ('multinomial', 'mul', 'log'):
            if hints.get(hint):
                expr, _ = _expand_hint(expr, hint, deep, hints, native)
        if expr == was:
            return expr


# SymPy's expand_<hint> functions: the hints each one enables.
_EXPAND_FUNCTIONS = {
    'expand_func': {'func': True},
    'expand_mul': {'mul': True},
    'expand_multinomial': {'multinomial': True},
    'expand_log': {'log': True},
    'expand_trig': {'trig': True},
    'expand_power_base': {'power_base': True},
    'expand_power_exp': {'power_exp': True},
}


def install_expand(namespace):
    """Route ``Expr.expand``, ``expand`` and the ``expand_<hint>`` functions
    of ``namespace`` through ``expand`` for trees with expansion hooks."""
    import functools
    import sys
    method = _core.Expr.expand

    def native(expr, deep, hints):
        return method(expr, deep=deep, **hints)

    @functools.wraps(method)
    def expand_method(self, deep=True, **hints):
        if _hooks['expand']:
            result = expand(self, deep, hints, native)
            if result is not None:
                return result
        return method(self, deep=deep, **hints)

    # The native method, for callers that check ``_hooks['expand']`` first.
    expand_method._symox_native_expand = method
    _core.Expr.expand = expand_method
    replaced = {}
    root = namespace['expand']

    @functools.wraps(root)
    def expand_function(e, deep=True, **hints):
        if _hooks['expand']:
            from . import sympify
            value = sympify(e)
            if isinstance(value, _core.Expr):
                result = expand(value, deep, hints, native)
                if result is not None:
                    return result
        return root(e, deep=deep, **hints)

    replaced[id(root)] = ('expand', root, expand_function)
    for name, enabled in _EXPAND_FUNCTIONS.items():
        original = namespace.get(name)
        if original is None:
            continue

        def make(original, enabled):
            @functools.wraps(original)
            def function(expr, deep=True, **options):
                if _hooks['expand']:
                    from . import sympify
                    value = sympify(expr)
                    if isinstance(value, _core.Expr):
                        hints = {hint: False for hint in _DEFAULT_HINTS}
                        hints.update(enabled)
                        hints.update(options)
                        result = expand(value, deep, hints, native)
                        if result is not None:
                            return result
                return original(expr, deep=deep, **options)
            return function
        replaced[id(original)] = (name, original, make(original, enabled))
    for name, original, wrapper in replaced.values():
        namespace[name] = wrapper
    for module_name, module in tuple(sys.modules.items()):
        if not module_name.startswith('symox.') or module_name == 'symox._core' or module is None:
            continue
        for name, original, wrapper in replaced.values():
            if module.__dict__.get(name) is original:
                setattr(module, name, wrapper)
