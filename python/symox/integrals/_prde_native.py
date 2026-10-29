"""Admission for the native base-case parametric Risch DE kernels.

The kernels in ``symox._core._prde_kernel`` reproduce the stock bodies of
``symox.integrals.prde`` for ``DE.case == 'base'`` (``D = [Poly(1, x)]``)
with plain ZZ/QQ polynomials. They are used only while every Python
function those bodies reach (prde, rde, risch, PolyMatrix), every global
or builtin name they read, and the Poly, PolyMatrix and
DifferentialExtension class dictionaries are the ones captured at import.
Each entry point returns ``None`` to request the unchanged Python body.
"""
import builtins as _builtins
from types import FunctionType, BuiltinFunctionType
from .. import _core

_TYPE = type
_LEN = len
_STR = str
_INT = int
_LIST = list
_MISSING = object()
_ROOT = None


def _state(f):
    return (f, f.__code__, f.__defaults__, f.__kwdefaults__, f.__closure__,
            tuple(f.__dict__.items()))


def register_root(namespace):
    """Capture the stock prde protocol after package import."""
    global _ROOT
    if _ROOT is not None:
        return
    kernel = _core.__dict__.get('_prde_kernel')
    if _TYPE(kernel) is not BuiltinFunctionType or kernel.__self__ is not _core:
        return
    from importlib import import_module
    prde = import_module('symox.integrals.prde')
    rde = import_module('symox.integrals.rde')
    risch = import_module('symox.integrals.risch')
    polymatrix = import_module('symox.polys.polymatrix')
    cls = risch.DifferentialExtension
    matrix = polymatrix.PolyMatrix
    try:
        slots = tuple((name, cls.__dict__[name]) for name in ('case', 't', 'd', 'T', 'D', 'level'))
    except KeyError:
        return
    if any(_TYPE(s).__name__ != 'member_descriptor' for _, s in slots):
        return
    new = matrix.__dict__.get('_new')
    qq = polymatrix.__dict__.get('QQ')
    if _TYPE(new) is not classmethod or qq is None:
        return
    modules = (prde.__dict__, rde.__dict__, risch.__dict__, polymatrix.__dict__)
    functions, lookups, seen = [], [], set()

    def visit(f):
        if _TYPE(f) is not FunctionType or id(f) in seen:
            return
        if not any(f.__globals__ is m for m in modules):
            return
        seen.add(id(f))
        functions.append(_state(f))
        for cell in f.__closure__ or ():
            visit(cell.cell_contents)
        ns = f.__globals__
        for name in f.__code__.co_names:
            value = ns.get(name, _MISSING)
            lookups.append((ns, name, value))
            if value is _MISSING:
                lookups.append((_builtins.__dict__, name, _builtins.__dict__.get(name, _MISSING)))
            visit(value)
        for const in f.__code__.co_consts:
            if _TYPE(const) is _TYPE(f.__code__):
                for name in const.co_names:
                    value = ns.get(name, _MISSING)
                    lookups.append((ns, name, value))
                    if value is _MISSING:
                        lookups.append((_builtins.__dict__, name, _builtins.__dict__.get(name, _MISSING)))
                    visit(value)

    for f in list(prde.__dict__.values()):
        visit(f)
    for f in (rde.weak_normalizer, rde.bound_degree, rde.order_at_oo,
              risch.splitfactor, risch.derivation, risch.gcdex_diophantine):
        visit(f)
    for value in list(matrix.__dict__.values()):
        visit(value.__func__ if _TYPE(value) in (classmethod, staticmethod) else value)
    classes = tuple((c, tuple(c.__dict__.items())) for c in (_core.Poly, matrix, cls))
    _ROOT = (kernel, cls, slots, tuple(functions), tuple(lookups), classes, matrix, new.__func__, qq)


def _verified():
    root = _ROOT
    if root is None or _core.__dict__.get('_prde_kernel') is not root[0]:
        return None
    _, _, _, functions, lookups, classes, _, _, _ = root
    for f, code, defaults, kwdefaults, closure, attrs in functions:
        if (f.__code__ is not code or f.__defaults__ is not defaults
                or f.__kwdefaults__ is not kwdefaults or f.__closure__ is not closure):
            return None
        d = f.__dict__
        if _LEN(d) != _LEN(attrs) or any(d.get(k, _MISSING) is not v for k, v in attrs):
            return None
    for ns, name, value in lookups:
        if ns.get(name, _MISSING) is not value:
            return None
    for c, items in classes:
        d = c.__dict__
        if _LEN(d) != _LEN(items) or any(d.get(k, _MISSING) is not v for k, v in items):
            return None
    return root


def _fields(de):
    root = _ROOT
    if root is None or _TYPE(de) is not root[1]:
        return None
    cls, slots = root[1], root[2]
    try:
        values = [slot.__get__(de, cls) for _, slot in slots]
    except AttributeError:
        return None
    case = values[0]
    if _TYPE(case) is not _STR or case != 'base' or _verified() is None:
        return None
    return values[1:]


def _matrix(spec, gens):
    root = _ROOT
    rows, cols, items = spec
    return root[7](root[6], rows, cols, items, gens, root[8])


def _call(kind, args, fields):
    t, d, T, D, level = fields
    return _ROOT[0](kind, args, t, d, T, D, level)


def param_rischDE(fa, fd, G, DE):
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call('param_rischDE', (fa, fd, G), fields)
    if r is None:
        return None
    h, spec = r
    return h, _matrix(spec, (fields[0],))


def param_poly_rischDE(a, b, q, n, DE):
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call('param_poly_rischDE', (a, b, q, n), fields)
    if r is None:
        return None
    h, spec = r
    return h, _matrix(spec, (fields[0],))


def _no_cancel(kind, b, Q, n, DE):
    if _TYPE(Q) is not _LIST or _TYPE(n) is not _INT:
        return None
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call(kind, (b, Q, n), fields)
    if r is None:
        return None
    h, spec, final = r
    if kind == 'prde_no_cancel_b_large':
        ran = n >= 0
    else:
        ran = n >= 1 or (not b.is_zero and b.degree() > 0)
    if ran:
        Q[:] = final
    return h, _matrix(spec, (fields[0],))


def prde_no_cancel_b_large(b, Q, n, DE):
    return _no_cancel('prde_no_cancel_b_large', b, Q, n, DE)


def prde_no_cancel_b_small(b, Q, n, DE):
    return _no_cancel('prde_no_cancel_b_small', b, Q, n, DE)


def prde_spde(a, b, Q, n, DE):
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call('prde_spde', (a, b, Q, n), fields)
    if r is None:
        return None
    B, Qq, R, n1 = r
    return (a, B, Qq, R, n1)


def prde_normal_denom(fa, fd, G, DE):
    fields = _fields(DE)
    if fields is None:
        return None
    return _call('prde_normal_denom', (fa, fd, G), fields)


def prde_special_denom(a, ba, bd, G, DE, case):
    if _TYPE(case) is not _STR or case not in ('auto', 'base', 'primitive'):
        return None
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call('prde_special_denom', (ba, bd), fields)
    if r is None:
        return None
    B, one = r
    return (a, B, G, one)


def prde_linear_constraints(a, b, G, DE):
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call('prde_linear_constraints', (G,), fields)
    if r is None:
        return None
    qs, spec = r
    return qs, _matrix(spec, (fields[0],))


def poly_linear_constraints(p, d):
    if _TYPE(d) is not _core.Poly or _verified() is None:
        return None
    gens = d.gens
    if _TYPE(gens) is not tuple or _LEN(gens) != 1:
        return None
    r = _ROOT[0]('poly_linear_constraints', (p,), gens[0], d, None, None, None)
    if r is None:
        return None
    q, spec = r
    return q, _matrix(spec, gens)


def constant_system(A, u, DE):
    root = _ROOT
    if root is None or _TYPE(A) is not root[6] or _TYPE(u) is not root[6]:
        return None
    if not (A._rows * A._cols) or A.domain is not root[8] or u.domain is not root[8]:
        return None
    gens = A.gens
    fields = _fields(DE)
    if fields is None or gens != (fields[0],):
        return None
    r = _call('constant_system', (A._rows, A._cols, _LIST(A._items), u._rows, u._cols, _LIST(u._items)), fields)
    if r is None:
        return None
    a, v = r
    return _matrix(a, gens), _matrix(v, gens)
