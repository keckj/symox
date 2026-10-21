"""Ported ``sympy.simplify`` helpers and the matching ``Expr`` methods.

The algorithms run natively (``crates/symox-core/src/simptools``); this module
only adapts the SymPy call signatures and installs the ``Expr`` methods.
"""
import math as _math

from . import _core
from . import _noncommutative as _nc
from ._core import Expr

__all__ = [
    "fraction", "numer", "denom", "collect_const", "rcollect", "radsimp",
    "signsimp", "posify", "separatevars", "hypersimp", "hypersimilar",
    "powsimp", "powdenest", "factor_terms", "gcd_terms", "factor_nc",
    "sqrtdenest", "combsimp", "gammasimp", "ratsimp", "ratsimpmodprime",
    "fu", "exptrigsimp", "trigsimp",
]


def _sympify(e):
    from . import sympify
    return sympify(e)


def _is_container(e):
    return isinstance(e, (list, tuple, set, frozenset, dict))


# --- radsimp ---------------------------------------------------------------

def fraction(expr, exact=False):
    """``(numerator, denominator)`` of ``expr`` without any simplification."""
    return _core._st_fraction(_sympify(expr), bool(exact))


def numer(expr, exact=False):
    return fraction(expr, exact)[0]


def denom(expr, exact=False):
    return fraction(expr, exact)[1]


def collect_const(expr, *vars, Numbers=True):
    """Collect terms of an Add sharing a constant (or ``vars``) coefficient."""
    return _core._st_collect_const(_sympify(expr), [_sympify(v) for v in vars], bool(Numbers))


def rcollect(expr, *vars):
    """Recursively collect sums in ``expr`` with respect to ``vars``."""
    return _core._st_rcollect(_sympify(expr), [_sympify(v) for v in vars])


def radsimp(expr, symbolic=True, max_terms=4):
    """Rationalize denominators by removing square roots (SymPy ``radsimp``)."""
    expr = _sympify(expr)
    if _nc.is_nc(expr):
        return _nc.masked(lambda e: _core._st_radsimp(e, bool(symbolic), int(max_terms)), expr)
    return _core._st_radsimp(expr, bool(symbolic), int(max_terms))


# --- simplify.simplify -----------------------------------------------------

def signsimp(expr, evaluate=None):
    """Make the sign of Add terms canonical (SymPy ``signsimp``)."""
    if _is_container(expr):
        return expr
    return _core._st_signsimp(_sympify(expr), True if evaluate is None else bool(evaluate))


def posify(eq):
    """Replace symbols without a sign assumption by positive dummies."""
    if _is_container(eq) or (hasattr(eq, "__iter__") and not isinstance(eq, (Expr, str))):
        kind = type(eq)
        items = [_sympify(e) for e in eq]
        out, reps = _core._st_posify_seq(items)
        try:
            out = kind(out)
        except TypeError:
            pass
        return out, dict(reps)
    out, reps = _core._st_posify(_sympify(eq))
    return out, dict(reps)


def separatevars(expr, symbols=[], dict=False, force=False):  # noqa: B006
    """Separate the variables of ``expr`` into a product, if possible."""
    expr = _sympify(expr)
    if not dict:
        return _core._st_separatevars(expr, bool(force))
    syms = None if symbols is None else [_sympify(s) for s in symbols]
    out = _core._st_separatevars_dict(expr, syms, bool(force))
    if out is None:
        return None
    return {("coeff" if k is None else k): v for k, v in out}


def hypersimp(f, k):
    """Simplified ``f(k+1)/f(k)`` if ``f`` is hypergeometric, else ``None``."""
    return _core._st_hypersimp(_sympify(f), _sympify(k))


def hypersimilar(f, g, k):
    """True if ``f/g`` is a rational function of ``k``."""
    return _core._st_hypersimilar(_sympify(f), _sympify(g), _sympify(k))


# --- powsimp ---------------------------------------------------------------

def powsimp(expr, deep=False, combine="all", force=False, measure=None):
    """Combine powers with similar bases and exponents (SymPy ``powsimp``)."""
    if combine not in ("all", "exp", "base"):
        raise ValueError("combine must be one of ('all', 'exp', 'base').")
    if _is_container(expr):
        return expr
    return _core._st_powsimp(_sympify(expr), bool(deep), combine, bool(force))


def powdenest(eq, force=False, polar=False):
    """Collect exponents of nested powers (SymPy ``powdenest``)."""
    if polar:
        raise NotImplementedError("powdenest(polar=True) is not implemented")
    return _core._st_powdenest(_sympify(eq), bool(force))


# --- exprtools -------------------------------------------------------------

def factor_terms(expr, radical=False, clear=False, fraction=False, sign=True):
    """Remove common factors from terms without expanding (SymPy ``factor_terms``)."""
    if isinstance(expr, (list, tuple, set)):
        return type(expr)([factor_terms(e, radical, clear, fraction, sign) for e in expr])
    if isinstance(expr, dict):
        return expr
    expr = _sympify(expr)
    if _nc.is_nc(expr):
        return _nc.factor_terms(expr, radical, clear, fraction, sign)
    return _core._st_factor_terms(expr, bool(radical), bool(clear), bool(fraction), bool(sign))


def gcd_terms(terms, isprimitive=False, clear=True, fraction=True):
    """Factor the gcd out of an Add or a sequence of terms (SymPy ``gcd_terms``)."""
    if isinstance(terms, (list, tuple, set)):
        terms = [_sympify(t) for t in terms]
        if _nc._any_nc(terms):
            return _nc.gcd_terms(terms, isprimitive, clear, fraction)
        return _core._st_gcd_terms_list(terms, bool(isprimitive), bool(clear), bool(fraction))
    if isinstance(terms, dict):
        return {k: gcd_terms(v, isprimitive, clear, fraction) for k, v in terms.items()}
    terms = _sympify(terms)
    if _nc.is_nc(terms):
        return _nc.gcd_terms(terms, isprimitive, clear, fraction)
    return _core._st_gcd_terms(terms, bool(isprimitive), bool(clear), bool(fraction))


def factor_nc(expr):
    """Factor a noncommutative expression (SymPy ``factor_nc``)."""
    expr = _sympify(expr)
    if _nc.is_nc(expr):
        return _nc.factor_nc(expr)
    return _core._st_factor_nc(expr)


# --- sqrtdenest, gammasimp, ratsimp, fu ------------------------------------

def sqrtdenest(expr, max_iter=3):
    """Denest square roots in ``expr`` (SymPy ``sqrtdenest``)."""
    return _core._st_sqrtdenest(_sympify(expr), int(max_iter))


def gammasimp(expr):
    """Simplify gamma functions and combinatorial functions as gammas."""
    expr = _sympify(expr)
    if _nc.is_nc(expr):
        return _nc.combsimp(expr, _core._st_gammasimp)
    return _core._st_gammasimp(expr)


def combsimp(expr):
    """Simplify combinatorial expressions (SymPy ``combsimp``)."""
    expr = _sympify(expr)
    if _nc.is_nc(expr):
        return _nc.combsimp(expr, _core._st_combsimp)
    return _core._st_combsimp(expr)


def ratsimp(expr):
    """Put an expression over a common denominator, then cancel and reduce."""
    expr = _sympify(expr)
    if _nc.is_nc(expr):
        return _nc.ratsimp(expr)
    return _core._st_ratsimp(expr)


def ratsimpmodprime(expr, G, *gens, quick=True, polynomial=False, **args):
    """Simplify a rational expression modulo the prime ideal generated by ``G``."""
    for key in args:
        if key not in ("order", "domain", "modulus", "field"):
            raise NotImplementedError("ratsimpmodprime option %r is not supported" % key)
    order = args.get("order", "grevlex")
    order = str(getattr(order, "alias", order))
    modulus, field = int(args.get("modulus") or 0), bool(args.get("field", False))
    domain = args.get("domain")
    if domain is not None:
        if isinstance(domain, str):
            import re
            m = re.fullmatch(r"\s*(?:GF|FF)\((\d+)\)\s*", domain)
            domain = _core.GF(int(m.group(1))) if m else {"QQ": _core.QQ, "ZZ": _core.ZZ}.get(domain.strip(), domain)
        if getattr(domain, "is_FiniteField", False):
            modulus = int(domain.characteristic())
        elif getattr(domain, "is_QQ", False):
            field = True
        elif not getattr(domain, "is_ZZ", False):
            raise NotImplementedError("ratsimpmodprime over %s is not supported" % (domain,))
    if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
        gens = tuple(gens[0])
    try:
        return _core._st_ratsimpmodprime(_sympify(expr), [_sympify(g) for g in G],
                                         [_sympify(g) for g in gens], bool(quick), bool(polynomial), order,
                                         modulus, field)
    except ValueError as err:
        if "generator" in str(err):
            from . import PolynomialError
            raise PolynomialError(str(err)) from None
        raise


def fu(rv, measure=None):
    """Fu-Harrington-Wang trigonometric simplification (SymPy ``fu``)."""
    return _core._st_fu(_sympify(rv), measure)


def exptrigsimp(expr):
    """Simplify exponential/trigonometric/hyperbolic mixtures."""
    return _core._st_exptrigsimp(_sympify(expr))


def futrig(e, *, hyper=True, **kwargs):
    """Fu-based trigonometric simplification (``trigsimp`` method 'matching')."""
    e = _sympify(e)
    if not isinstance(e, Expr):
        return e
    return _core._st_futrig(e, bool(hyper))


def trigsimp(expr, inverse=False, **opts):
    """Simplify trigonometric expressions (SymPy ``trigsimp``).

    Methods include native matching/Fu and polynomial Groebner reduction.
    """
    if hasattr(expr, "applyfunc") and not isinstance(expr, Expr):
        return expr.applyfunc(lambda x: trigsimp(x, inverse=inverse, **opts))
    expr = _sympify(expr)
    if getattr(expr, "is_Relational", False):
        return expr.func(trigsimp(expr.lhs, inverse=inverse, **opts),
                         trigsimp(expr.rhs, inverse=inverse, **opts))
    legacy = opts.pop("old", False)
    opts.pop("deep", None)
    opts.pop("recursive", None)
    method = opts.pop("method", "matching")
    if legacy:method='old'
    if method == "fu":
        result = fu(expr, **opts)
    elif method == "matching":
        result = futrig(expr)
        if result != 0:
            from .simplify.trigsimp import _common_angles, trigsimp_groebner
            if _common_angles(result) is not None:
                try:result = trigsimp_groebner(result)
                except (ValueError,NotImplementedError):pass
    elif method in ("groebner", "combined"):
        from .simplify.trigsimp import trigsimp_groebner
        result=trigsimp_groebner(expr,**opts)
        if method=='combined':result=futrig(result)
    elif method=='old':
        from .simplify.trigsimp import trigsimp_old
        result=trigsimp_old(expr,**opts)
    else:
        raise KeyError(method)
    if inverse:
        result = _core._st_trigsimp_inverse(result)
    return result


def _rule(name, *params):
    def tr(rv, *args, **kw):
        values = dict(zip(params, args))
        values.update(kw)
        opts = {}
        for k, v in values.items():
            if k not in params:
                raise TypeError("%s() got an unexpected keyword argument %r" % (name, k))
            if k in ("max", "base"):
                opts[k] = _sympify(v)
            else:
                opts[k] = bool(v)
        return _core._st_tr(name, _sympify(rv), **opts)
    tr.__name__ = name
    tr.__doc__ = "Fu transformation rule ``%s``." % name
    return tr


_TR_PARAMS = {
    "TR0": (), "TR1": (), "TR2": (), "TR2i": ("half",), "TR3": (), "TR4": (),
    "TR5": ("max", "pow"), "TR6": ("max", "pow"), "TR7": (), "TR8": ("first",), "TR9": (),
    "TR10": ("first",), "TR10i": (), "TR11": ("base",), "_TR11": (), "TR12": ("first",),
    "TR12i": (), "TR13": (), "TR14": ("first",), "TR15": ("max", "pow"),
    "TR16": ("max", "pow"), "TR111": (), "TR22": ("max", "pow"), "TRmorrie": (),
    "TRpower": (), "sincos_to_sum": (),
}


def L(rv):
    """Number of trigonometric functions in ``rv``."""
    from . import Integer
    return Integer(_core._st_fu_l(_sympify(rv)))


def _TR56(rv, f, g, h, max, pow):
    """Shared helper of TR5 and TR6: ``f(x)**n -> h(g(x)**2)``."""
    name = lambda fn: getattr(fn, "__name__", str(fn))
    return _core._st_tr56(_sympify(rv), name(f), name(g), lambda x: _sympify(h(x)), _sympify(max), bool(pow))


def process_common_addends(rv, do, key2=None, key1=True):
    """Apply ``do`` to the addends of ``rv`` sharing a common coefficient."""
    from . import Add, S
    absc = {}
    for a in rv.args:
        if key1:
            c, a = a.as_coeff_Mul()
            if c < 0:
                c, a = -c, -a
        else:
            c = S.One
        absc.setdefault((c, key2(a) if key2 else 1), []).append(a)
    args, hit = [], False
    for (c, _), v in absc.items():
        if len(v) > 1:
            e = Add(*v, evaluate=False)
            new = do(e)
            v = new if new != e else Add(*v)
            hit = hit or new != e
        else:
            v = v[0]
        args.append(c * v)
    return Add(*args) if hit else rv


def trig_split(a, b, two=False):
    """``(gcd, n1, n2, a, b, flag)`` when ``a + b`` is a sin/cos sum pattern."""
    r = _core._st_trig_split(_sympify(a), _sympify(b), bool(two))
    return None if r is None else r


def as_f_sign_1(e):
    """``(g, f, s)`` with ``e = g*(f + s)`` and ``s = +/-1``, else None."""
    return _core._st_as_f_sign_1(_sympify(e))


def _osborne(e, d):
    return _core._st_osborne(_sympify(e), _sympify(d), False)


def _osbornei(e, d):
    return _core._st_osborne(_sympify(e), _sympify(d), True)


def hyper_as_trig(rv):
    """``(trig_expr, inverse)``: hyperbolics rewritten with Osborne's rule."""
    t, d, back = _core._st_hyper_as_trig(_sympify(rv))
    return t, lambda x: _core._st_hyper_as_trig_back(_sympify(x), d, back)


_TR_NAMES = list(_TR_PARAMS) + ["L", "_TR56", "process_common_addends", "trig_split",
                               "as_f_sign_1", "hyper_as_trig", "_osborne", "_osbornei"]


class _FuNamespace:
    """``sympy.simplify.fu`` lookalike holding the TR rules."""

    def __init__(self):
        for n, params in _TR_PARAMS.items():
            setattr(self, n, _rule(n, *params))
        for n in _TR_NAMES:
            if n not in _TR_PARAMS:
                setattr(self, n, globals()[n])
        self.fu = fu


fu_rules = _FuNamespace()


# --- Expr methods ----------------------------------------------------------

def _decompose_power(expr):
    from . import Pow, Rational, Mul, Add
    base, exp = expr.as_base_exp()
    if exp.is_Number:
        if exp.is_Rational:
            if not exp.is_Integer:
                base = Pow(base, Rational(1, exp.q))
            return base, int(exp.p)
        return expr, 1
    c, tail = exp.as_coeff_Mul(rational=True)
    if c == -1:
        return Pow(base, tail), -1
    if c != 1:
        q = Rational(1, c.q)
        tail = Mul(q, tail, evaluate=False) if isinstance(tail, Add) else q * tail
        return Pow(base, tail), int(c.p)
    return expr, 1


def _as_terms(self):
    from . import Add, Mul, default_sort_key
    gens, terms = set(), []
    for term in (self.args if isinstance(self, Add) else (self,)):
        coeff, _term = term.as_coeff_Mul()
        coeff = complex(coeff)
        cpart, ncpart = {}, []
        if _term != 1:
            for factor in (_term.args if isinstance(_term, Mul) else (_term,)):
                if factor.is_number:
                    try:
                        coeff *= complex(factor)
                    except (TypeError, ValueError):
                        pass
                    else:
                        continue
                if factor.is_commutative and not factor.is_Order:
                    base, exp = _decompose_power(factor)
                    cpart[base] = exp
                    gens.add(base)
                else:
                    ncpart.append(factor)
        terms.append((term, ((coeff.real, coeff.imag), cpart, tuple(ncpart))))
    gens = sorted(gens, key=default_sort_key)
    indices = {g: i for i, g in enumerate(gens)}
    result = []
    for term, (coeff, cpart, ncpart) in terms:
        monom = [0] * len(gens)
        for base, exp in cpart.items():
            monom[indices[base]] = exp
        result.append((term, (coeff, tuple(monom), ncpart)))
    return result, gens


def _equals(self, other, failing_expression=False):
    other = _sympify(other)
    if _nc.is_nc(self) or _nc.is_nc(other):
        # Expr.equals: a zero difference proves equality; the numerical
        # constancy tests cannot decide a noncommutative difference.
        if self == other:
            return True
        from . import simplify
        diff = factor_terms(simplify(self - other), radical=True)
        if not diff:
            return True
        return diff if failing_expression else None
    rv = _core._st_equals(self, other)
    if rv is None and failing_expression:
        from . import simplify
        return factor_terms(simplify(self - other), radical=True)
    return rv


def _is_constant(self, *wrt, **flags):
    for key in flags:
        if key != "simplify":
            raise NotImplementedError("is_constant flag %r is not supported" % key)
    return _core._st_is_constant(self, [_sympify(w) for w in wrt], bool(flags.get("simplify", True)))


def _invert(self, g, *gens, **args):
    from . import mod_inverse, invert
    if self.is_number and getattr(g, "is_number", True):
        return mod_inverse(self, g)
    return invert(self, g, *gens, **args)


def _mag(x):
    from . import Integer
    xpos = abs(x.n())
    if not xpos:
        return 0
    try:
        mag = int(_math.ceil(_math.log10(float(xpos))))
    except (ValueError, OverflowError):
        from . import log
        mag = int(_math.ceil(float((log(xpos) / log(10)).n(30))))
    if (xpos / Integer(10) ** mag) >= 1:
        mag += 1
    return mag


def _prec_to_dps(prec):
    return max(1, int(round(int(prec) / 3.3219280948873626) - 1))


def _round(self, n=None):
    from . import Float, Integer, Rational, Pow, re, im, I, nan, oo, zoo
    x = self
    if not x.is_number:
        raise TypeError("Cannot round symbolic expression")
    if not x.is_Atom:
        h, t = x.n(2).as_coeff_Add()
        ok = (not t) or (t.as_coeff_Mul()[1] == I)
        if not ok:
            raise TypeError("Expected a number but got %s:" % type(x).__name__)
    elif x in (nan, oo, -oo, zoo):
        return x
    xr = x.is_extended_real
    if not xr:
        r, i = re(x), im(x)
        if xr is None and i.is_number and i.equals(0) is False:
            xr = False
        if xr is False:
            return _round(r, n) + I * _round(i, n)
        if i.equals(0):
            return _round(r, n)
    if x == 0:
        return Integer(0) if n is None else x
    p = int(n or 0)
    if x.is_Integer:
        return Integer(round(int(x), p))
    digits_to_decimal = _mag(x)
    allow = digits_to_decimal + p
    precs = [f._prec for f in x.atoms(Float)]
    dps = _prec_to_dps(max(precs)) if precs else None
    if dps is None:
        dps = max(15, allow)
    else:
        allow = min(allow, dps)
    shift = -digits_to_decimal + dps
    extra = 1
    xf = x.n(dps + extra) * Pow(10, shift)
    xi = Integer(int(xf))
    sign = 1 if x > 0 else -1
    dif2 = sign * (xf - xi).n(extra)
    if dif2 < 0:
        raise NotImplementedError("not expecting int(x) to round away from 0")
    if dif2 > .5:
        xi += sign
    elif dif2 == .5:
        xi += sign if int(xi) % 2 else -sign
    ip = p - shift
    xr = round(int(xi), ip)
    rv = Rational(xr, 10 ** shift)
    if rv.is_Integer:
        if n is None:
            return rv
        return Float(str(rv), dps)
    if not allow and rv > self:
        allow += 1
    return Float(rv, allow)


def _preorder(e):
    yield e
    for a in e.args:
        yield from _preorder(a)


def _head_name(query):
    """Name of a symox function constructor such as ``sin``, else None."""
    if isinstance(query, (type, Expr)) or not callable(query):
        return None
    if not str(getattr(query, "__module__", "")).startswith("symox"):
        return None
    return getattr(query, "__name__", None)


def _head_test(name):
    def test(expr):
        return not expr.is_Atom and str(expr.func) == name
    return test


def _count(self, query):
    if not isinstance(query, type) and not callable(query) or isinstance(query, Expr):
        try:
            query = _sympify(query)
        except Exception:
            pass
    if _head_name(query):
        test = _head_test(_head_name(query))
    elif isinstance(query, type):
        def test(expr):
            return isinstance(expr, query)
    elif isinstance(query, Expr):
        def test(expr):
            return expr.match(query) is not None
    else:
        test = query
    return sum(bool(test(sub)) for sub in _preorder(self))


def _separate(self, deep=False, force=False):
    from . import expand_power_base
    return expand_power_base(self, deep=deep, force=force)


def _compare(self, other):
    if not isinstance(other, Expr):
        return -1
    return _core._st_compare(self, other)


def _dummy_eq(self, other, symbol=None):
    from . import Dummy
    o = _sympify(other)
    dummies = [i for i in self.free_symbols if i.is_Dummy]
    if len(dummies) != 1:
        return self == o
    dummy = dummies[0]
    if symbol is None:
        symbols = o.free_symbols
        if len(symbols) != 1:
            return self == o
        symbol = next(iter(symbols))
    tmp = Dummy()
    return self.xreplace({dummy: tmp}) == o.xreplace({symbol: tmp})


def _iterfreeargs(e):
    def children(value):
        return value if isinstance(value, (tuple, list)) else getattr(value, 'args', ())

    queue = [e]
    for node in queue:
        yield node
        bound = getattr(node, 'bound_symbols', ())
        if bound:
            # Canonicalize the binder so occurrences of a bound variable
            # cannot match a query for the original free symbol. Continue
            # through dependent ancestors to retain their free descendants.
            canonical = node.as_dummy()
            hidden = set(node.canonical_variables.values())
            pending = list(children(canonical))
            for value in pending:
                pending.extend(children(value))
                if isinstance(value, (tuple, list)):
                    continue
                if value not in hidden and not any(
                        getattr(value, 'has', lambda *a: False)(h) for h in hidden):
                    yield value
        else:
            queue.extend(children(node))


def _has_xfree(self, s):
    if type(s) is not set:
        raise TypeError("expecting set argument")
    return any(a in s for a in _iterfreeargs(self))


def _has_matcher(pattern):
    cls = pattern.func
    c = set(pattern.args)

    def is_in(expr):
        if expr.func == cls:
            return expr == pattern or (c & set(expr.args)) == c
        return False
    return is_in


def _has_free(self, *patterns):
    from . import Add, Mul
    if not patterns:
        return False
    p0 = patterns[0]
    if len(patterns) == 1 and hasattr(p0, "__iter__") and not isinstance(p0, (Expr, str)):
        raise TypeError("Expecting 1 or more Basic args, not a single non-Basic iterable. "
                        "Don't forget to unpack iterables: `eq.has_free(*patterns)`")
    types, values = [], set()
    for p in patterns:
        if _head_name(p):
            if any(_head_test(_head_name(p))(i) for i in _iterfreeargs(self)):
                return True
            continue
        if isinstance(p, type):
            types.append(p)
            continue
        try:
            values.add(_sympify(p))
        except Exception:
            continue
    if _has_xfree(self, values):
        return True
    types = tuple(types)
    for i in _iterfreeargs(self):
        if i in values or (types and isinstance(i, types)):
            return True
    for p in values:
        if isinstance(p, (Add, Mul)):
            match = _has_matcher(p)
            if any(match(a) for a in _iterfreeargs(self)):
                return True
    return False


def _postorder(e):
    for a in e.args:
        yield from _postorder(a)
    yield e


def _is_same(a, b, approx=None):
    from itertools import zip_longest
    from . import Number
    for x, y in zip_longest(_postorder(a), _postorder(_sympify(b))):
        if x is None or y is None:
            return False
        if isinstance(x, Number):
            if not isinstance(y, Number):
                return False
            if approx:
                return approx(x, y)
        if not (x == y and x.func == y.func):
            return False
    return True


def _rcall(self, *args):
    from . import Lambda

    def rec(e):
        if isinstance(e, Lambda):
            return e(*args)
        if e.args:
            return e.func(*[rec(a) for a in e.args])
        return e
    return rec(self)


_PERIODIC = ("sin", "cos", "tan", "cot", "sec", "csc")


def _as_set(self):
    from . import S, Intersection, Union, solve_univariate_inequality, ConditionSet
    if self == S.true:
        return S.UniversalSet
    if self == S.false:
        return S.EmptySet
    free = self.free_symbols
    if len(free) != 1:
        raise NotImplementedError("Sorry, as_set has not yet been implemented for multivariate expressions")
    x = next(iter(free))

    def rec(e):
        head = "" if e.is_Atom else str(e.func)
        if head == "And":
            return Intersection(*[rec(a) for a in e.args])
        if head == "Or":
            return Union(*[rec(a) for a in e.args])
        if head == "Not":
            return rec(e.args[0]).complement(S.Reals)
        if e == S.true:
            return S.UniversalSet
        if e == S.false:
            return S.EmptySet
        if not e.is_Relational:
            raise NotImplementedError("as_set is not implemented for %s" % e)
        if any(str(a.func) in _PERIODIC for a in _preorder(e)):
            raise NotImplementedError("as_set is not implemented for relationals with periodic solutions")
        try:
            return solve_univariate_inequality(e, x, relational=False)
        except NotImplementedError:
            return ConditionSet(x, e, S.Reals)
    return rec(self)


def _install():
    methods = {
        "ratsimp": lambda self: ratsimp(self),
        "trigsimp": lambda self, **kw: trigsimp(self, **kw),
        "combsimp": lambda self: combsimp(self),
        "gammasimp": lambda self: gammasimp(self),
        "powsimp": lambda self, *args, **kw: powsimp(self, *args, **kw),
        "radsimp": lambda self, **kw: radsimp(self, **kw),
        "normal": lambda self: _core._st_normal(self),
        "primitive": lambda self: _core._st_primitive(self),
        "extract_additively": lambda self, c: _core._st_extract_additively(self, _sympify(c)),
        "extract_multiplicatively": lambda self, c: _core._st_extract_multiplicatively(self, _sympify(c)),
        "as_coefficient": lambda self, expr: _core._st_as_coefficient(self, _sympify(expr)),
        "as_terms": _as_terms,
        "equals": _equals,
        "is_constant": _is_constant,
        "invert": _invert,
        "round": _round,
        "__round__": _round,
        "count": _count,
        "separate": _separate,
        "compare": _compare,
        "dummy_eq": _dummy_eq,
        "has_free": _has_free,
        "has_xfree": _has_xfree,
        "is_same": _is_same,
        "rcall": _rcall,
        "as_set": _as_set,
    }
    for name, fn in methods.items():
        try:
            setattr(Expr, name, fn)
        except (TypeError, AttributeError):  # pragma: no cover
            pass


_install()


# Public Fu rule registry, backed by the existing native transformations.
FU = {n: getattr(fu_rules, n) for n in (
    'TR0 TR1 TR2 TR3 TR4 TR5 TR6 TR7 TR8 TR9 TR10 TR10i TR11 '
    'TR12 TR13 L TR2i TRmorrie TR12i TR14 TR15 TR16 TR111 TR22').split()}


def kroneckersimp(expr):
    """Remove products with incompatible Kronecker equality constraints."""
    from itertools import combinations
    from . import Eq, Integer, sympify
    expr = sympify(expr)

    def visit(e):
        if not e.args:
            return e
        args = tuple(visit(a) for a in e.args)
        value = e if args == e.args else e.func(*args)
        if value.is_Mul:
            deltas = [a for a in value.args if a.class_name == 'KroneckerDelta']
            for a, b in combinations(deltas, 2):
                for i in range(2):
                    for j in range(2):
                        if Eq(a.args[i], b.args[j]) == True and Eq(a.args[1-i], b.args[1-j]) == False:
                            return Integer(0)*value
        return value
    return visit(expr)
