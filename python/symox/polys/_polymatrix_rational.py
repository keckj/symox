"""Admission for exact rational entries of ``PolyMatrix`` (wave 30, prdefrac).

``PolyMatrix`` keeps its entries as expressions and builds a ``Poly`` in
``self.gens`` over ``self.domain`` each time an entry is read. For the
constant matrices of the parametric Risch DE (exact Integer/Rational entries
over QQ, or over a fraction field ``ZZ(a, ...)``/``QQ(a, ...)`` of symbols)
that public ``Poly(e, *gens, domain=K)`` call only parses options before
reaching the native constructor ``Poly._from_normalized``. These helpers call
the native constructor directly, and prove the ground check and the ``unify``
chain of ``from_list`` without building the entries, while the stock
protocol captured at import is unchanged:

- the ``Poly`` and ``QQ`` names of ``symox.polys.polymatrix`` and the
  ``Poly`` class entries that path reads (``__new__``, ``_from_normalized``,
  ``unify``, ``is_ground``, ``domain``, ``gens``, ``as_expr``, ``__add__``);
- the module ``symox._poly_constructor`` that ``Poly.__new__`` imports, its
  ``construct``/``_options``/``_domain_alias``/``Poly`` names, and the
  functions ``construct`` reaches for a scalar entry (``sympify``,
  ``ring_domain``, ``iterable``), with their code, defaults and closures
  (read by the native ``_function_states_unchanged`` through the public
  CPython getters, else by the ``_polymatrix_native`` helpers);
- the ``PolyMatrix`` class entries ``_poly``, ``__getitem__``, ``__iter__``
  (``_ground_check`` iterates ``self``).

Every helper returns ``None``/``False`` to request the unchanged stock code.
"""
from sys import modules as _MODULES

_TYPE = type
_LEN = len
_ROOT = None
_MISSING = object()
_FUNCTION = (lambda: None).__class__
_BUILTIN = [].append.__class__


def _states_unchanged(records):
    # The native check (bound at import) reads the live function fields
    # through the public CPython getters (no ``object.__getattr__`` audit
    # event). Without those getters it returns None, and the shared
    # ``_polymatrix_native`` helpers run instead, with the same reader.
    native = _ROOT[14]
    if native is not None:
        result = native(records)
        if result is not None:
            return result
    ground = _ROOT[13]
    reader = ground._metadata_reader()
    for record in records:
        if not ground._unchanged(record, reader):
            return False
    return True


def register(cls, namespace):
    """Capture the stock protocol once, after ``polymatrix`` is executed."""
    global _ROOT
    if _ROOT is not None:
        return
    import symox
    from symox import _core
    from symox import _poly_constructor as constructor
    from symox import _modular_poly
    from symox.utilities import iterables
    poly = namespace.get('Poly')
    qq, zz = namespace.get('QQ'), getattr(symox, 'ZZ', None)
    if poly is not _core.Poly or qq is None or zz is None or _TYPE(qq) is not _core.Domain:
        return
    names = ('__new__', '_from_normalized', 'unify', 'is_ground', 'domain', 'gens', 'as_expr', '__add__')
    entries = tuple((name, poly.__dict__.get(name, _MISSING)) for name in names)
    if any(v is _MISSING for _, v in entries):
        return
    cnames = ('construct', '_options', '_domain_alias', 'Poly')
    cvalues = tuple((name, constructor.__dict__.get(name, _MISSING)) for name in cnames)
    if any(v is _MISSING for _, v in cvalues) or constructor.__dict__['Poly'] is not poly:
        return
    from . import _polymatrix_native as ground
    functions = (constructor.construct, constructor._options, constructor._domain_alias,
                 symox.sympify, _modular_poly.ring_domain, iterables.iterable,
                 poly.__dict__['unify'], poly.__dict__['as_expr'])
    states = tuple(ground._function(f) for f in functions if _TYPE(f) is _FUNCTION)
    native = _core.__dict__.get('_function_states_unchanged')
    if (_TYPE(native) is not _BUILTIN or native.__self__ is not _core
            or native.__name__ != '_function_states_unchanged'):
        native = None
    links = ((symox.__dict__, 'sympify', symox.sympify),
             (_modular_poly.__dict__, 'ring_domain', _modular_poly.ring_domain),
             (iterables.__dict__, 'iterable', iterables.iterable),
             (namespace, 'Poly', poly), (namespace, 'QQ', qq))
    methods = tuple((name, cls.__dict__.get(name, _MISSING)) for name in ('_poly', '__getitem__', '__iter__'))
    S = symox.S
    rationals = frozenset(_TYPE(v) for v in (S.Zero, S.One, S.NegativeOne, S.Half,
                                              symox.Integer(2), symox.Rational(1, 3)))
    symbols = frozenset((symox.Symbol, symox.Dummy))
    _ROOT = (poly, qq, zz, entries, constructor, cvalues, states,
             links, cls, methods, rationals, symbols, _TYPE(qq), ground, native)


def _stock():
    root = _ROOT
    if root is None:
        return None
    poly, _, _, entries, constructor, cvalues, functions, links, cls, methods, _, _, _, _, _ = root
    d = poly.__dict__
    for name, value in entries:
        if d.get(name, _MISSING) is not value:
            return None
    if _MODULES.get('symox._poly_constructor') is not constructor:
        return None
    d = constructor.__dict__
    for name, value in cvalues:
        if d.get(name, _MISSING) is not value:
            return None
    for ns, name, value in links:
        if ns.get(name, _MISSING) is not value:
            return None
    if not _states_unchanged(functions):
        return None
    d = cls.__dict__
    for name, value in methods:
        if d.get(name, _MISSING) is not value:
            return None
    return root


def _gens(gens, symbols):
    # One commutative stock symbol: ``_options`` keeps it as is and cannot raise.
    if _TYPE(gens) is not tuple or _LEN(gens) != 1:
        return False
    g = gens[0]
    return _TYPE(g) in symbols and g.is_commutative is True


def _field(domain, root, gen):
    """True for QQ, and for a fraction field ``ZZ(a, ...)``/``QQ(a, ...)`` of
    stock symbols disjoint from ``gen``. For those, stock ``construct`` with
    an exact rational entry is ``Poly._from_normalized(e, gens, domain, None,
    None)`` (not a ``str``, not a ``ModularRing``, not ``RR``, no option), the
    constant ``Poly`` cannot raise and is ground, and the domain is a field."""
    if domain is root[1]:
        return True
    if _TYPE(domain) is not root[12] or domain.is_FractionField is not True or domain.is_Field is not True:
        return False
    ground = domain.dom
    if _TYPE(ground) is not root[12] or (ground.is_ZZ is not True and ground.is_QQ is not True):
        return False
    symbols = root[11]
    names = domain.symbols
    if _TYPE(names) is not tuple or not names:
        return False
    for s in names:
        if _TYPE(s) not in symbols or (s == gen) is not False:
            return False
    return True


def ground_poly(matrix, e):
    """``Poly(e, *matrix.gens, domain=matrix.domain)`` for an exact rational
    entry over an admitted field (see ``_field``), or None."""
    root = _ROOT
    if root is None or _TYPE(e) not in root[10] or _TYPE(matrix) is not root[8]:
        return None
    gens, domain = matrix.gens, matrix.domain
    if not _gens(gens, root[11]) or not _field(domain, root, gens[0]) or _stock() is None:
        return None
    return root[0]._from_normalized(e, gens, domain, None, None)


def ground_items(matrix):
    """True when every entry is an exact rational over an admitted field (the
    stock ground check then passes: the domain is a field and each ``Poly`` is
    a constant)."""
    root = _ROOT
    if root is None or _TYPE(matrix) is not root[8]:
        return False
    rationals = root[10]
    for e in matrix._items:
        if _TYPE(e) not in rationals:
            return False
    gens = matrix.gens
    return (_gens(gens, root[11]) and _field(matrix.domain, root, gens[0])
            and _stock() is not None)


def common_ring(items):
    """``(domain, gens)`` of the stock ``unify`` chain of ``from_list`` when all
    items are native ``Poly`` objects over the same ZZ or QQ domain object and
    the same generators (each ``unify`` then keeps that domain), else None."""
    root = _ROOT
    if root is None or not items:
        return None
    poly = root[0]
    first = items[0]
    if _TYPE(first) is not poly:
        return None
    domain, gens = first.domain, first.gens
    if domain is not root[1] and domain is not root[2]:
        return None
    if not _gens(gens, root[11]):
        return None
    g = gens[0]
    for p in items:
        if _TYPE(p) is not poly or p.domain is not domain:
            return None
        other = p.gens
        if _LEN(other) != 1 or _TYPE(other[0]) is not _TYPE(g) or (other[0] == g) is not True:
            return None
    if _stock() is None:
        return None
    return domain, tuple(gens)
