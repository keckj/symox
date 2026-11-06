"""Admission for the native base-case parametric Risch DE kernels.

The kernels in ``symox._core._prde_kernel`` reproduce the stock bodies of
``symox.integrals.prde`` for ``DE.case == 'base'`` (``D = [Poly(1, x)]``)
with plain ZZ/QQ polynomials, and ``constant_system`` for rational QQ
matrices in any DE state. They are used only while every Python
function those bodies reach (prde, rde, risch, PolyMatrix), every global
or builtin name they read, and the Poly, PolyMatrix and
DifferentialExtension class dictionaries are the ones captured at import.
Each entry point returns ``None`` to request the unchanged Python body.
"""
import builtins as _builtins
from types import FunctionType, BuiltinFunctionType, ModuleType
from .. import _core

_TYPE = type
_LEN = len
_ZIP = zip
_STR = str
_INT = int
_LIST = list
_MISSING = object()
_ROOT = None
# Number of calls answered by the native kernels (per entry point).
STATS = {'constant_system_field': 0, 'derivation': 0, 'is_deriv_k': 0,
         'is_log_deriv_k_t_radical_in_field': 0, 'parametric_log_deriv_heu': 0,
         'limited_integrate': 0, 'is_log_deriv_k_t_radical_in_field_tower': 0,
         'parametric_log_deriv_heu_tower': 0, 'param_rischDE_tower': 0,
         'param_poly_rischDE_tower': 0, 'prde_cancel_liouvillian_tower': 0,
         'limited_integrate_tower': 0, 'prde_no_cancel_b_large_tower': 0,
         'prde_no_cancel_b_small_tower': 0}
_POLY = _core.Poly
_EXPR = _core.Expr
_AS_POLY = None
_DERIV_K = None
# The ``DifferentialExtension.cases`` slot descriptor (tower kernels).
_CASES = None
_NOT_IMPLEMENTED = NotImplementedError
_VALUE_ERROR = ValueError
_POLYNOMIAL_ERROR = None


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
    global _AS_POLY
    reader = _core.__dict__.get('_ddm_function_metadata')
    if _TYPE(reader) is not BuiltinFunctionType or reader.__self__ is not _core:
        return
    from .. import _poly_constructor
    _AS_POLY = (_core.Poly.__dict__.get('as_poly'), _core.Poly.__dict__.get('diff'),
                _core.Poly.__dict__.get('as_expr'), _core.Poly.__dict__.get('__add__'),
                _poly_constructor.__dict__.get('construct'), _core.Expr.__dict__.get('as_poly'),
                risch.__dict__, risch.derivation, risch.derivation.__code__, reader,
                _poly_constructor.__dict__)
    global _DERIV_K
    from .. import _compat, _public_core
    import symox
    try:
        ext_slots = (cls.__dict__['exts'], cls.__dict__['extargs'])
    except KeyError:
        return
    if any(_TYPE(s).__name__ != 'member_descriptor' for s in ext_slots):
        return
    # Once an AccumBounds exists, Add and Mul consult the stock coefficient
    # fold, which returns None unless an argument is (or directly holds) an
    # AccumulationBounds. The kernel only builds host-free expressions, so the
    # fold is admitted when it is the stock function with its stock globals.
    from importlib import import_module
    bounds = import_module('symox.calculus.accumulationbounds')
    fold = bounds.__dict__.get('_fold_coefficient')
    if _TYPE(fold) is not FunctionType or fold.__globals__ is not bounds.__dict__:
        return
    folding = (fold, reader(fold), bounds.__dict__, bounds.__dict__.get('AccumulationBounds'),
               tuple((name, _builtins.__dict__[name]) for name in ('type', 'getattr', 'any')))
    if folding[1] is None or folding[3] is None:
        return
    _DERIV_K = (_public_core.__dict__.get('_evaluation'), _compat.__dict__.get('_accum_fold'),
                prde.__dict__.get('Dummy'), ext_slots, symox.__dict__, symox.__dict__.get('cancel'),
                prde.__dict__.get('cancel'), folding, _compat.__dict__)
    global _CASES, _POLYNOMIAL_ERROR
    # The stock residue_reduce raises the PolynomialError of the polys
    # package (kkinv gcd, status 3 of the tower kernels).
    _POLYNOMIAL_ERROR = import_module('symox.polys.polyerrors').PolynomialError
    cases = cls.__dict__.get('cases')
    if _TYPE(cases).__name__ == 'member_descriptor':
        _CASES = cases


def _stock_fold(fold, folding):
    """``fold`` (the installed ``_compat._accum_fold[0]``) is the stock
    ``_fold_coefficient`` with its captured code, defaults, closure, globals
    and the builtins it reads."""
    function, record, ns, bounds, names = folding
    if fold is not function:
        return False
    values = _AS_POLY[9](function)
    if (values is None or _LEN(values) != _LEN(record)
            or any(a is not b for a, b in _ZIP(values, record))
            or ns.get('AccumulationBounds') is not bounds):
        return False
    b = ns.get('__builtins__')
    if _TYPE(b) is ModuleType:
        b = b.__dict__
    if _TYPE(b) is not dict:
        return False
    return all(name not in ns and b.get(name) is value for name, value in names)


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


def _derivation_protocol():
    """The Poly methods, the constructor and ``risch.derivation`` (code,
    globals, builtins) are the ones captured at import."""
    c = _POLY.__dict__
    captured = _AS_POLY
    if (captured is None or c.get('as_poly') is not captured[0] or c.get('diff') is not captured[1]
            or c.get('as_expr') is not captured[2] or c.get('__add__') is not captured[3]
            or captured[10].get('construct') is not captured[4]
            or _EXPR.__dict__.get('as_poly') is not captured[5]):
        return False
    # The stock body reads Poly, len and zip: a module global or a builtin
    # replaced for risch.derivation keeps the stock body (and its calls).
    ns, f = captured[6], captured[6].get('derivation')
    if f is not captured[7] or ns.get('Poly') is not _POLY or 'len' in ns or 'zip' in ns:
        return False
    # The code object is read by the native metadata reader (no audit event).
    values = captured[9](f)
    if values is None or values[0] is not captured[8]:
        return False
    b = f.__builtins__
    return b.get('len') is _LEN and b.get('zip') is _ZIP


def _expr_state():
    """The kernels that build expressions (``Add``, ``Mul``, ``Pow``,
    ``cancel``) need the stock ``derivation`` protocol, global evaluation
    switched on, the stock ``cancel`` and the stock AccumBounds fold."""
    captured = _DERIV_K
    if captured is None or not _derivation_protocol():
        return False
    evaluation, accum, _, _, public, cancel, prde_cancel, folding, compat = captured
    if (evaluation.get() is not True or public.get('cancel') is not cancel
            or cancel is not prde_cancel or compat.get('_accum_fold') is not accum
            or _TYPE(accum) is not _LIST or _LEN(accum) != 1):
        return False
    fold = accum[0]
    return fold is None or _stock_fold(fold, folding)


def is_deriv_k(fa, fd, DE):
    """``is_deriv_k(fa, fd, DE)`` for stock Polys over ZZ, QQ, ZZ[...],
    QQ[...], ZZ(...) or QQ(...) whose tower is elementary (``len(DE.exts)
    == len(DE.D)``; otherwise the stock body raises). The native kernel
    reproduces the whole body (``Df/f``, the PolyMatrix row, the
    fraction-field ``constant_system``, the constant and rational checks and
    the answer). It is used only while the stock prde protocol, the
    ``derivation`` protocol and the global evaluation state are the ones
    captured at import, and the Dummy generator is created here exactly as
    the stock body does.

    Returns ``None`` (stock body, which creates its own Dummy),
    ``(True, result)``, or ``(False, dum)`` (stock body with the Dummy
    already created).
    """
    root = _ROOT
    if root is None or _TYPE(fa) is not _POLY or _TYPE(fd) is not _POLY or _TYPE(DE) is not root[1]:
        return None
    if not _expr_state() or _verified() is None:
        return None
    dummy, (exts_slot, extargs_slot) = _DERIV_K[2], _DERIV_K[3]
    slots, cls = root[2], root[1]
    try:
        t = slots[1][1].__get__(DE, cls)
        T = slots[3][1].__get__(DE, cls)
        D = slots[4][1].__get__(DE, cls)
        level = slots[5][1].__get__(DE, cls)
        exts = exts_slot.__get__(DE, cls)
        extargs = extargs_slot.__get__(DE, cls)
    except AttributeError:
        return None
    if (_TYPE(T) is not _LIST or _TYPE(D) is not _LIST or _TYPE(level) is not _INT
            or _TYPE(exts) is not _LIST or _TYPE(extargs) is not _LIST or _LEN(exts) != _LEN(D)):
        return None
    dum = dummy()
    r = root[0]('is_deriv_k', (fa, fd, dum, exts, extargs), t, None, T, D, level)
    if r is None:
        return (False, dum)
    STATS['is_deriv_k'] += 1
    return (True, r[0])


def derivation(p, DE, coefficientD):
    """``derivation(p, DE, coefficientD, basic=False)`` for a stock Poly
    ``p`` over ZZ, QQ, ZZ[...], QQ[...], ZZ(...) or QQ(...) (exp, primitive,
    base or any other case: the stock body does not depend on it). The
    native kernel runs the stock loop with the native Poly constructor; it
    is used only while the Poly methods and the constructor the stock body
    reaches are the ones captured at import. ``None`` requests the Python
    body.
    """
    root = _ROOT
    if root is None or _TYPE(p) is not _POLY or _TYPE(DE) is not root[1]:
        return None
    if not _derivation_protocol():
        return None
    slots, cls = root[2], root[1]
    try:
        t = slots[1][1].__get__(DE, cls)
        T = slots[3][1].__get__(DE, cls)
        D = slots[4][1].__get__(DE, cls)
        level = slots[5][1].__get__(DE, cls)
    except AttributeError:
        return None
    if _TYPE(T) is not _LIST or _TYPE(D) is not _LIST or _TYPE(level) is not _INT:
        return None
    r = root[0]('derivation', (p, coefficientD), t, None, T, D, level)
    if r is not None:
        STATS['derivation'] += 1
    return r


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


def _above_base(DE):
    """True when the active level of ``DE`` is an extension (not 'base')."""
    root = _ROOT
    if root is None or _TYPE(DE) is not root[1]:
        return False
    try:
        current = root[2][0][1].__get__(DE, root[1])
    except AttributeError:
        return False
    return _TYPE(current) is _STR and current != 'base'


def _tower_matrix(kind, args, DE):
    """Tower kernels answering ``(list, PolyMatrix)``: ``None`` (stock body)
    or the stock result, the matrix rebuilt from ``(rows, cols, items, gens,
    domain)``."""
    r = _tower(kind, args, DE, kind)
    if r is None:
        return None
    h, (rows, cols, items, gens, domain) = r[0]
    root = _ROOT
    return h, root[7](root[6], rows, cols, items, gens, domain)


def param_rischDE(fa, fd, G, DE):
    if _above_base(DE):
        if _TYPE(G) is not _LIST or _TYPE(fa) is not _POLY or _TYPE(fd) is not _POLY:
            return None
        return _tower_matrix('param_rischDE_tower', (fa, fd, G), DE)
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call('param_rischDE', (fa, fd, G), fields)
    if r is None:
        return None
    # The kernel answers when the resultant of weak_normalizer has no z:
    # the stock body creates its Dummy('z') and no other Dummy.
    _create(('z',))
    h, spec = r
    return h, _matrix(spec, (fields[0],))


def param_poly_rischDE(a, b, q, n, DE):
    if _above_base(DE):
        if (_TYPE(q) is not _LIST or _TYPE(n) is not _INT or _TYPE(a) is not _POLY
                or _TYPE(b) is not _POLY):
            return None
        return _tower_matrix('param_poly_rischDE_tower', (a, b, q, n), DE)
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
    if _above_base(DE):
        if _TYPE(b) is not _POLY:
            return None
        r = _tower(kind + '_tower', (b, Q, n), DE, kind + '_tower')
        if r is None:
            return None
        h, (rows, cols, items, gens, domain), final = r[0]
        root = _ROOT
        matrix = root[7](root[6], rows, cols, items, gens, domain)
    else:
        fields = _fields(DE)
        if fields is None:
            return None
        r = _call(kind, (b, Q, n), fields)
        if r is None:
            return None
        h, spec, final = r
        matrix = _matrix(spec, (fields[0],))
    if kind == 'prde_no_cancel_b_large':
        ran = n >= 0
    else:
        ran = n >= 1 or (not b.is_zero and b.degree() > 0)
    if ran:
        Q[:] = final
    return h, matrix


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
    """Rational QQ systems in any DE state (base, exp, primitive, or a tower
    still being built): the derivation loop of the stock body only runs for
    entries that involve ``DE.T``, which rational entries never do. The kernel
    also requires the matrix generator and every ``DE.T`` element to be plain
    stock symbols, so ``rref``'s ground check and ``has(*DE.T)`` cannot raise.
    """
    root = _ROOT
    if root is None or _TYPE(A) is not root[6] or _TYPE(u) is not root[6]:
        return None
    if not (A._rows * A._cols):
        return None
    gens = A.gens
    if _TYPE(gens) is not tuple or _LEN(gens) != 1 or _TYPE(DE) is not root[1]:
        return None
    if A.domain is root[8] and u.domain is root[8]:
        try:
            T = root[2][3][1].__get__(DE, root[1])
        except AttributeError:
            return None
        if _TYPE(T) is not _LIST or _verified() is None:
            return None
        r = _call('constant_system', (A._rows, A._cols, _LIST(A._items), u._rows, u._cols, _LIST(u._items)),
                  (gens[0], None, T, None, None))
        if r is None:
            return None
        a, v = r
        return _matrix(a, gens), _matrix(v, gens)
    return _constant_system_field(A, u, DE, gens)


def _constant_system_field(A, u, DE, gens):
    """Fraction-field systems (ZZ(x), QQ(x, t0), ...) in any DE state: the
    rref, the derivation loop and the PolyMatrix domains are computed by the
    native kernel, which calls the Poly constructor for the domain inference
    of ``PolyMatrix(items, gens)``. As in the stock body, the result uses the
    generators of ``A`` (``A.row_join(u)``).
    """
    root = _ROOT
    domain_type = _TYPE(root[8])
    if _TYPE(A.domain) is not domain_type or _TYPE(u.domain) is not domain_type:
        return None
    if _TYPE(A._items) is not _LIST or _TYPE(u._items) is not _LIST:
        return None
    cls = root[1]
    try:
        T = root[2][3][1].__get__(DE, cls)
        D = root[2][4][1].__get__(DE, cls)
        level = root[2][5][1].__get__(DE, cls)
    except AttributeError:
        return None
    if _TYPE(T) is not _LIST or _TYPE(D) is not _LIST or _TYPE(level) is not _INT or _verified() is None:
        return None
    r = root[0]('constant_system_field',
                (A._rows, A._cols, _LIST(A._items), A.domain, u._rows, u._cols, _LIST(u._items), u.domain),
                gens[0], None, T, D, level)
    if r is None:
        return None
    (ar, ac, ai, ad), (ur, uc, ui, ud) = r
    if ad == root[8]:
        ad = root[8]
    if ud == root[8]:
        ud = root[8]
    new, matrix = root[7], root[6]
    STATS['constant_system_field'] += 1
    return new(matrix, ar, ac, ai, gens, ad), new(matrix, ur, uc, ui, gens, ud)


def limited_integrate_reduce(fa, fd, G, DE):
    if _TYPE(G) is not _LIST:
        return None
    fields = _fields(DE)
    if fields is None:
        return None
    r = _call('limited_integrate_reduce', (fa, fd, G), fields)
    if r is None:
        return None
    a, b, N, g, V = r
    return (a, b, a, N, g, V)


def _log_fields(fa, fd, DE):
    """Base-case fields for the expression-building kernels, or ``None``."""
    if _TYPE(fa) is not _POLY or _TYPE(fd) is not _POLY:
        return None
    fields = _fields(DE)
    if fields is None or not _expr_state():
        return None
    return fields


def _create(trace):
    """Create the Dummies of the stock body, in its order."""
    dummy = _DERIV_K[2]
    for name in trace:
        if name:
            dummy(name)
        else:
            dummy()


def _tower(kind, args, DE, stat):
    """Tower kernels (levels above the base): ``None`` (stock body) or
    ``(answer,)``. The stock NotImplementedError, ValueError and
    PolynomialError (status 1, 2, 3) are raised here after the Dummies
    created before them;
    ``STATS[stat]`` counts every native answer, raised ones included."""
    root = _ROOT
    if root is None or _CASES is None or _TYPE(DE) is not root[1]:
        return None
    cls, slots = root[1], root[2]
    try:
        case, t, d, T, D, level = [slot.__get__(DE, cls) for _, slot in slots]
        cases = _CASES.__get__(DE, cls)
    except AttributeError:
        return None
    if (_TYPE(T) is not _LIST or _TYPE(D) is not _LIST or _TYPE(level) is not _INT
            or _TYPE(cases) is not _LIST or _TYPE(case) is not _STR
            or _LEN(cases) != _LEN(T) or _LEN(D) != _LEN(T)):
        return None
    active = _LEN(T) + level + 1
    if active < 2 or active > _LEN(T):
        return None
    top = cases[active - 1]
    if _TYPE(top) is not _STR or top != case:
        return None
    if not _expr_state() or _verified() is None:
        return None
    r = root[0](kind, args + (_LIST(cases),), t, d, T, D, level)
    if r is None:
        return None
    status, payload, trace = r
    _create(trace)
    STATS[stat] += 1
    if status == 1:
        raise _NOT_IMPLEMENTED(payload)
    if status == 2:
        raise _VALUE_ERROR(payload)
    if status == 3:
        raise _POLYNOMIAL_ERROR(payload)
    return (payload,)


def is_log_deriv_k_t_radical_in_field(fa, fd, DE, case, z):
    """``is_log_deriv_k_t_radical_in_field(fa, fd, DE, case, z)`` for
    ``z=None``. Base case (``case`` 'auto' or 'base'): the kernel reproduces
    ``cancel``, ``residue_reduce`` (subresultant PRS, square-free resultant,
    the ``invert`` normalization), the rational root check, ``p in k`` and
    the base branch. Exp/primitive towers (any ``case`` string): the same
    steps over the rational function field of the tower, the exp branch
    (``parametric_log_deriv`` one level down), the primitive branch
    (recursive call), the base branch, the stock exceptions and the final
    ``cancel``. The Dummies (``z``, one per ``real_roots`` call, ``i``,
    ``c1`` of nested calls) are created as the stock body does.

    Returns ``None`` (stock body) or ``(answer,)`` with ``answer`` the stock
    result (``None`` or ``(n, u)``).
    """
    if z is not None or _TYPE(case) is not _STR or _TYPE(fa) is not _POLY or _TYPE(fd) is not _POLY:
        return None
    root = _ROOT
    if root is None or _TYPE(DE) is not root[1]:
        return None
    try:
        current = root[2][0][1].__get__(DE, root[1])
    except AttributeError:
        return None
    if _TYPE(current) is _STR and current != 'base':
        return _tower('is_log_deriv_tower', (fa, fd, case), DE,
                      'is_log_deriv_k_t_radical_in_field_tower')
    if case not in ('auto', 'base'):
        return None
    fields = _log_fields(fa, fd, DE)
    if fields is None:
        return None
    r = _call('is_log_deriv_in_field', (fa, fd), fields)
    if r is None:
        return None
    answer, trace = r
    _create(trace)
    STATS['is_log_deriv_k_t_radical_in_field'] += 1
    return (answer,)


def parametric_log_deriv_heu(fa, fd, wa, wd, DE, c1):
    """``parametric_log_deriv_heu(fa, fd, wa, wd, DE, c1=None)``. Base case
    (``B = 0``): the linear equations in ``c1``, ``solve`` for one unknown
    and the base-case ``is_log_deriv_k_t_radical_in_field``. Exp/primitive
    towers: the same steps over the rational function field of the tower
    (``splitfactor`` included). The Dummies (``c1``, then those of
    ``is_log_deriv_k_t_radical_in_field``) are created as the stock body
    does. Returns ``None`` (stock body) or ``(answer,)``.
    """
    if (c1 is not None or _TYPE(fa) is not _POLY or _TYPE(fd) is not _POLY
            or _TYPE(wa) is not _POLY or _TYPE(wd) is not _POLY):
        return None
    root = _ROOT
    if root is None or _TYPE(DE) is not root[1]:
        return None
    try:
        current = root[2][0][1].__get__(DE, root[1])
    except AttributeError:
        return None
    if _TYPE(current) is _STR and current != 'base':
        return _tower('heu_tower', (fa, fd, wa, wd), DE,
                      'parametric_log_deriv_heu_tower')
    fields = _log_fields(fa, fd, DE)
    if fields is None:
        return None
    r = _call('parametric_log_deriv_heu', (fa, fd, wa, wd), fields)
    if r is None:
        return None
    answer, trace = r
    _create(trace)
    STATS['parametric_log_deriv_heu'] += 1
    return (answer,)


def limited_integrate(fa, fd, G, DE):
    """Base-case ``limited_integrate(fa, fd, G, DE)``: the native
    ``param_rischDE``, the nullspace, ``y`` built with ``Add``/``Mul`` and
    ``as_numer_denom`` and the normalized ``(Ya, Yd)``; ``C`` holds
    ``Poly(c, x, domain=QQ)``. Returns ``None`` (stock body) or
    ``(answer,)``.
    """
    if _TYPE(G) is not _LIST:
        return None
    if _above_base(DE):
        if _TYPE(fa) is not _POLY or _TYPE(fd) is not _POLY:
            return None
        return _tower('limited_integrate_tower', (fa, fd, G), DE, 'limited_integrate_tower')
    fields = _log_fields(fa, fd, DE)
    if fields is None:
        return None
    r = _call('limited_integrate', (fa, fd, G), fields)
    if r is None:
        return None
    # Dummy('z') of weak_normalizer in the inner param_rischDE.
    _create(('z',))
    STATS['limited_integrate'] += 1
    return r


def prde_cancel_liouvillian(b, Q, n, DE):
    """``prde_cancel_liouvillian(b, Q, n, DE)`` on exp/primitive towers: the
    loop over ``i``, the lower-level ``param_rischDE`` (native base case or
    tower), ``set_gens`` and the block matrix. Returns ``None`` (stock body)
    or the stock ``(H, M)``."""
    if (_TYPE(Q) is not _LIST or _TYPE(n) is not _INT or _TYPE(b) is not _POLY
            or not _above_base(DE)):
        return None
    return _tower_matrix('prde_cancel_liouvillian_tower', (b, Q, n), DE)
