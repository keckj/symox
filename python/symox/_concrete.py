"""SymPy ``ExprWithIntLimits`` protocol for held ``Sum``/``Product`` nodes.

Convergence tests and limit manipulations are native; this module only
resolves symbolic indices, attaches the methods to ``Expr`` and raises
SymPy's ``ReorderError``.
"""
from . import _core
from ._core import Expr


class ReorderError(NotImplementedError):
    """Exception raised when trying to reorder dependent limits."""
    def __init__(self, expr, msg):
        super().__init__("%s could not be reordered: %s." % (expr, msg))


def _require(self, name, *heads):
    if not any(getattr(self, 'is_' + h) for h in heads):
        raise AttributeError("'%s' object has no attribute '%s'" % (type(self).__name__, name))


def _is_convergent(self):
    """Convergence of a one-limit Sum or Product."""
    return _core._is_convergent(self)


def _is_absolutely_convergent(self):
    """``Sum(Abs(f)).is_convergent()`` for a one-limit Sum."""
    return _core._is_absolutely_convergent(self)


def _index(self, x):
    _require(self, 'index', 'Sum', 'Product')
    return _core._limits_index(self, x)


def _resolve(self, indices):
    return [i if isinstance(i, int) else self.index(i) for i in indices]


def _reverse_order(self, *indices):
    _require(self, 'reverse_order', 'Sum', 'Product')
    return _core._reverse_order(self, _resolve(self, indices))


def _change_index(self, var, trafo, newvar=None):
    _require(self, 'change_index', 'Sum', 'Product')
    return _core._change_index(self, var, trafo, newvar)


def _reorder_limit(self, x, y):
    _require(self, 'reorder_limit', 'Sum', 'Product')
    try:
        return _core._reorder_limit(self, x, y)
    except NotImplementedError as error:
        raise ReorderError(self, str(error)) from None


def _reorder(self, *arg):
    _require(self, 'reorder', 'Sum', 'Product')
    new_expr = self
    for r in arg:
        if len(r) != 2:
            raise ValueError(r, "Invalid number of arguments")
        index1, index2 = _resolve(self, r)
        new_expr = new_expr.reorder_limit(index1, index2)
    return new_expr


def _limits_property(name, *heads):
    def getter(self):
        _require(self, name, *heads)
        return _core._limits_property(self, name)
    getter.__name__ = name
    return property(getter)


def _term(self):
    _require(self, 'term', 'Product')
    return self.args[0]


Expr.is_convergent = _is_convergent
Expr.is_absolutely_convergent = _is_absolutely_convergent
Expr.index = _index
Expr.reverse_order = _reverse_order
Expr.change_index = _change_index
Expr.reorder_limit = _reorder_limit
Expr.reorder = _reorder
Expr.has_finite_limits = _limits_property('has_finite_limits', 'Sum', 'Product', 'Integral')
Expr.has_reversed_limits = _limits_property('has_reversed_limits', 'Sum', 'Product', 'Integral')
Expr.has_empty_sequence = _limits_property('has_empty_sequence', 'Sum', 'Product')
Expr.term = property(_term)


# Indexed objects enter native arithmetic as symbols named like ``f[n]``
# (codegen.ast.Indexed._symox_), which hides their dependence on the index.
# The symbols are recorded by _public_matrix, and the evaluation of held sums and products
# replaces those whose indices are not constant by applications of an
# undefined function of the indices, so the summation engine substitutes the
# summation variable into them as SymPy substitutes into Indexed.
_INDEXED_HEADS = {}
_REGISTRY = []


def _indexed_mapping(expr):
    if not _REGISTRY:
        from ._public_matrix import _indexed
        _REGISTRY.append(_indexed)
    _INDEXED_SYMBOLS = _REGISTRY[0]
    if not _INDEXED_SYMBOLS:
        return None
    if type(expr).__name__ == 'Indexed' and hasattr(expr, '_symox_'):
        expr = expr._symox_()
    if not isinstance(expr, Expr):
        return None
    try:
        free = expr.free_symbols
    except Exception:
        return None
    flattened = [s for s in free if s in _INDEXED_SYMBOLS]
    if not flattened:
        return None
    from . import Function
    forward, heads = {}, {}
    for symbol in flattened:
        indexed = _INDEXED_SYMBOLS[symbol]
        if not any(getattr(i, 'free_symbols', None) for i in indexed.indices):
            continue
        key = (indexed.base, len(indexed.indices))
        head = _INDEXED_HEADS.get(key)
        if head is None:
            head = _INDEXED_HEADS[key] = Function('_symox_indexed_%d' % len(_INDEXED_HEADS))
        heads[head] = key[0]
        forward[symbol] = head(*indexed.indices)
    return (forward, heads, expr) if forward else None


def _indexed_restore(value, heads):
    from .codegen.ast import Indexed
    if not isinstance(value, Expr):
        return value
    for head, base in heads.items():
        value = value.replace(head, lambda *indices, _base=base: Indexed(_base, *indices)._symox_())
    return value


def _through_indexed(native):
    def call(expr, *args, **kwargs):
        mapping = _indexed_mapping(expr)
        if mapping is None:
            return native(expr, *args, **kwargs)
        forward, heads, expr = mapping
        return _indexed_restore(native(expr.xreplace(forward), *args, **kwargs), heads)
    call.__name__ = getattr(native, '__name__', 'call')
    call.__doc__ = getattr(native, '__doc__', None)
    call._symox_native = native
    return call


Expr.doit = _through_indexed(Expr.doit)
summation = _through_indexed(_core.summation)
product = _through_indexed(_core.product)


# Numerical evaluation of infinite sums (SymPy's evalf_sum): hypergeometric
# terms use exact fixed-point summation or Richardson extrapolation (hypsum),
# other terms use Euler-Maclaurin summation.

def _check_convergence(numer, denom, n):
    from . import Poly
    npol = Poly(numer, n)
    dpol = Poly(denom, n)
    p = npol.degree()
    q = dpol.degree()
    rate = q - p
    if rate:
        return rate, None, None
    constant = dpol.LC() / npol.LC()
    if abs(constant) != 1:
        return rate, constant, None
    if npol.degree() == dpol.degree() == 0:
        return rate, constant, 0
    pc = npol.all_coeffs()[1]
    qc = dpol.all_coeffs()[1]
    return rate, constant, (qc - pc)/dpol.LC()


def _integer_polynomial(num, den, n):
    # The term ratio num/den with integer coefficients, as Horner tables.
    from . import Poly, ilcm
    pn, pd = Poly(num, n), Poly(den, n)
    coeffs = pn.all_coeffs() + pd.all_coeffs()
    if not all(c.is_Rational for c in coeffs):
        raise NotImplementedError("a rational term ratio is required")
    m = 1
    for c in coeffs:
        m = ilcm(m, c.q)
    to_int = lambda p: [int(c*m) for c in p.all_coeffs()]
    return to_int(pn), to_int(pd)


def _horner(coeffs, k):
    v = 0
    for c in coeffs:
        v = v*k + c
    return v


def _hypsum(expr, n, start, prec):
    from . import Float, hypersimp
    from mpmath import nsum, inf as mpmath_inf, workprec, make_mpf
    from mpmath.libmp import from_man_exp, prec_to_dps
    if prec == float('inf'):
        raise NotImplementedError('does not support inf prec')
    if start:
        expr = expr.subs(n, n + start)
    hs = hypersimp(expr, n)
    if hs is None:
        raise NotImplementedError("a hypergeometric series is required")
    num, den = hs.as_numer_denom()
    h, g, p = _check_convergence(num, den, n)
    if h < 0:
        raise ValueError("Sum diverges like (n!)^%i" % (-h))
    term = expr.subs(n, 0)
    if not term.is_Rational:
        raise NotImplementedError("Non rational term functionality is not implemented.")
    c1, c2 = _integer_polynomial(num, den, n)
    if h > 0 or (h == 0 and abs(g) > 1):
        term = (int(term.p) << prec) // int(term.q)
        s = term
        k = 1
        while abs(term) > 5:
            term *= _horner(c1, k - 1)
            term //= _horner(c2, k - 1)
            s += term
            k += 1
        return from_man_exp(s, -prec)
    alt = g < 0
    if abs(g) < 1:
        raise ValueError("Sum diverges like (%i)^n" % abs(1/g))
    if p < 1 or (p == 1 and not alt):
        raise ValueError("Sum diverges like n^%i" % (-p))
    vold = None
    ndig = prec_to_dps(prec)
    while True:
        prec2 = 4*prec
        term0 = (int(term.p) << prec2) // int(term.q)

        def summand(k, _term=[term0]):
            if k:
                k = int(k)
                _term[0] *= _horner(c1, k - 1)
                _term[0] //= _horner(c2, k - 1)
            return make_mpf(from_man_exp(_term[0], -prec2))

        with workprec(prec):
            v = nsum(summand, [0, mpmath_inf], method='richardson')
        vf = Float(v, ndig)
        if vold is not None and vold == vf:
            break
        prec += prec
        vold = vf
    return v._mpf_


def _euler_maclaurin(self, m=0, n=0, eps=0, eval_integral=True):
    """Euler-Maclaurin approximation ``(s, e)`` of a one-limit Sum."""
    from . import S, Add, Dummy, Integral, bernoulli, factorial
    m = int(m)
    n = int(n)
    f = self.function
    if len(self.limits) != 1:
        raise ValueError("More than 1 limit")
    i, a, b = self.limits[0]
    if (a > b) == True:
        if a - b == 1:
            return S.Zero, S.Zero
        a, b = b + 1, a - 1
        f = -f
    s = S.Zero
    if m:
        if b.is_Integer and a.is_Integer:
            m = min(m, b - a + 1)
        if not eps or f.is_polynomial(i):
            s = Add(*[f.subs(i, a + k) for k in range(m)])
        else:
            term = f.subs(i, a)
            if term:
                test = abs(term.evalf(3)) < eps
                if test == True:
                    return s, abs(term)
                elif not (test == False):
                    return term, S.Zero
            s = term
            for k in range(1, m):
                term = f.subs(i, a + k)
                if abs(term.evalf(3)) < eps and term != 0:
                    return s, abs(term)
                s += term
        if b - a + 1 == m:
            return s, S.Zero
        a += m
    x = Dummy('x')
    integral = Integral(f.subs(i, x), (x, a, b))
    if eval_integral:
        integral = integral.doit()
    s += integral

    def fpoint(expr):
        if b is S.Infinity:
            return expr.subs(i, a), 0
        return expr.subs(i, a), expr.subs(i, b)
    fa, fb = fpoint(f)
    iterm = (fa + fb)/2
    g = f.diff(i)
    for k in range(1, n + 2):
        ga, gb = fpoint(g)
        term = bernoulli(2*k)/factorial(2*k)*(gb - ga)
        if k > n:
            break
        if eps and term:
            term_evalf = term.evalf(3)
            if term_evalf is S.NaN:
                return S.NaN, S.NaN
            if abs(term_evalf) < eps:
                break
        s += term
        g = g.diff(i, 2)
    return s + iterm, abs(term)


def evalf_infinite_sum(expr, dps):
    """The value of a closed Sum over ``[a, oo)`` at ``dps`` digits, or None."""
    from . import S, Float
    from mpmath.libmp import dps_to_prec
    if len(expr.limits) != 1 or len(expr.limits[0]) != 3:
        return None
    func = expr.function
    if func.is_zero:
        return S.Zero
    n, a, b = expr.limits[0]
    if b is not S.Infinity or a is S.NegativeInfinity or not a.is_Integer:
        return None
    prec = dps_to_prec(dps)
    prec2 = prec + 10
    try:
        v = _hypsum(func, n, int(a), prec2)
        if _fastlog(v) < -10:
            v = _hypsum(func, n, int(a), prec - _fastlog(v))
        return Float(_mpf(v), dps)
    except NotImplementedError:
        pass
    eps = Float(2.0)**(-prec)
    for i in range(1, 5):
        m = nn = 2**i * prec
        s, err = _euler_maclaurin(expr, m=m, n=nn, eps=eps, eval_integral=False)
        err = err.evalf()
        if err is S.NaN:
            return None
        if err <= eps:
            break
    value = s.evalf(dps + 3)
    if not value.is_number or value.has(S.NaN) or value.free_symbols:
        return None
    return value.evalf(dps)


def _mpf(value):
    from mpmath import make_mpf
    return make_mpf(value)


def _fastlog(x):
    # mpmath's fastlog: the binary magnitude of a raw mpf.
    if not x or x[1] == 0 and x[3] == 0:
        return float('-inf')
    return x[2] + x[3]


def _install_euler_maclaurin():
    from . import _node_types
    _node_types.node_class('Sum').euler_maclaurin = _euler_maclaurin


_install_euler_maclaurin()
