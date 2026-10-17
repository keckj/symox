"""High-level polynomial manipulation (``sympy.polys.polyfuncs``).

``symmetrize`` follows SymPy's reduction by elementary symmetric polynomials
on native ``Poly`` arithmetic; ``horner``, ``interpolate``,
``rational_interpolate`` and ``viete`` port the corresponding SymPy routines
with exact symox arithmetic.
"""
from itertools import islice

from .._core import Poly
from .polyerrors import (
    ComputationFailed, MultivariatePolynomialError, OptionError,
    PolificationFailed,
)
from .polytools import allowed_flags, poly_from_expr, parallel_poly_from_expr
from .specialpolys import interpolating_poly, symmetric_poly

__all__ = ['symmetrize', 'horner', 'interpolate', 'rational_interpolate', 'viete']


def _is_iterable(value):
    return isinstance(value, (list, tuple, set)) or (
        hasattr(value, '__iter__') and not isinstance(value, str)
        and not hasattr(value, 'is_Number') and not isinstance(value, Poly))


def _terms_expr(poly, gens):
    from .. import Add, Mul
    return Add(*[c*Mul(*[g**e for g, e in zip(gens, m)]) for m, c in poly.terms()])


def _symmetrize_poly(f, gens):
    """``(symmetric terms, remainder, mapping)`` for a Poly ``f`` over ``gens``."""
    n = len(gens)
    kw = {'domain': f.domain} if f.modulus is None else {'modulus': f.modulus}
    polys = [Poly(symmetric_poly(i + 1, *gens), *gens, **kw) for i in range(n)]
    powers = {}

    def power(i, k):
        if (i, k) not in powers:
            powers[(i, k)] = polys[i]**k
        return powers[(i, k)]

    indices = list(range(n - 1))
    weights = list(range(n, 0, -1))
    symmetric = []
    while not f.is_zero:
        best_height, best = -1, None
        for monom, coeff in f.terms():
            if all(monom[i] >= monom[i + 1] for i in indices):
                height = max(w*m for w, m in zip(weights, monom))
                if height > best_height:
                    best_height, best = height, (monom, coeff)
        if best is None:
            break
        monom, coeff = best
        exponents = [m1 - m2 for m1, m2 in zip(monom, monom[1:] + (0,))]
        symmetric.append((coeff, exponents))
        product = Poly(coeff, *gens, **kw)
        for i, k in enumerate(exponents):
            if k:
                product = product*power(i, k)
        f = f - product
    return symmetric, f, polys


def symmetrize(F, *gens, **args):
    """Rewrite a polynomial in terms of elementary symmetric polynomials.

    Returns ``(symmetric part, remainder)``; with ``formal=True`` the
    symmetric part is expressed in new symbols and the mapping to the
    elementary symmetric polynomials is appended.
    """
    from .. import Add, Mul, numbered_symbols, sympify
    allowed_flags(args, ['formal', 'symbols'])
    formal = args.pop('formal', False)
    symbols = args.pop('symbols', None)
    iterable = True
    if not _is_iterable(F):
        iterable = False
        F = [F]
    F = [f.as_expr() if isinstance(f, Poly) else sympify(f) for f in F]
    if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
        gens = tuple(gens[0])
    try:
        polys, opt = parallel_poly_from_expr(F, *gens, **args)
        gens = tuple(opt.gens)
    except PolificationFailed:
        polys, gens = None, ()
    symbols = iter(symbols) if symbols is not None else numbered_symbols('s', start=1)
    symbols = [next(symbols) for _ in range(len(gens))]
    result, mapping = [], []
    if polys is None:
        result = [(f, sympify(0)) for f in F]
    else:
        for f in polys:
            sym, rem, elementary = _symmetrize_poly(f, gens)
            sym_expr = Add(*[c*Mul(*[s**k for s, k in zip(symbols, ex)]) for c, ex in sym])
            result.append((sym_expr, _terms_expr(rem, gens)))
            mapping = elementary
    polys_map = [(s, p.as_expr()) for s, p in zip(symbols, mapping)]
    if not formal:
        result = [(sym.subs(polys_map), non_sym) for sym, non_sym in result]
    if not iterable:
        result, = result
    if not formal:
        return result
    if iterable:
        return result, polys_map
    return result + (polys_map,)


def horner(f, *gens, **args):
    """Rewrite a polynomial in Horner form."""
    from .. import S
    allowed_flags(args, [])
    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        return exc.expr
    form, gen = S.Zero, F.gen
    if F.is_univariate:
        for coeff in F.all_coeffs():
            form = form*gen + coeff
    else:
        rest = gens[1:]
        F = Poly(F.as_expr(), gen)
        for coeff in F.all_coeffs():
            form = form*gen + horner(coeff, *rest, **args)
    return form


def interpolate(data, x):
    """Lagrange interpolating polynomial through ``data`` evaluated at ``x``."""
    from .. import Dummy, S, sympify
    n = len(data)
    x = sympify(x)
    if isinstance(data, dict):
        if x in data:
            return S(data[x])
        X, Y = list(zip(*data.items()))
    else:
        if isinstance(data[0], tuple):
            X, Y = list(zip(*data))
            if x in X:
                return S(Y[X.index(x)])
        else:
            if x.is_Integer and 1 <= int(x) <= n:
                return S(data[int(x) - 1])
            Y = list(data)
            X = list(range(1, n + 1))
    X = [sympify(v) for v in X]
    Y = [sympify(v) for v in Y]
    try:
        return interpolating_poly(n, x, X, Y).expand()
    except ValueError:
        d = Dummy()
        return interpolating_poly(n, d, X, Y).expand().subs(d, x)


def rational_interpolate(data, degnum, X=None):
    """Rational interpolation with numerator degree ``degnum``."""
    from .. import Symbol, sympify, Matrix
    if X is None:
        X = Symbol('x')
    xdata, ydata = list(zip(*data))
    xdata = [sympify(v) for v in xdata]
    ydata = [sympify(v) for v in ydata]
    k = len(xdata) - degnum - 1
    if k < 0:
        raise OptionError("Too few values for the required degree.")
    rows, cols = degnum + k + 1, degnum + k + 2
    c = [[sympify(1)]*cols for _ in range(rows)]
    for j in range(max(degnum, k)):
        for i in range(rows):
            c[i][j + 1] = c[i][j]*xdata[i]
    for j in range(k + 1):
        for i in range(rows):
            c[i][degnum + k + 1 - j] = -c[i][k - j]*ydata[i]
    r = Matrix(c).nullspace()[0]
    return (sum(r[i]*X**i for i in range(degnum + 1))
            / sum(r[i + degnum + 1]*X**i for i in range(k + 1)))


def viete(f, roots=None, *gens, **args):
    """Vieta's formulas ``[(elementary symmetric polynomial, coefficient), ...]``."""
    from .. import Basic, numbered_symbols
    allowed_flags(args, [])
    if isinstance(roots, Basic):
        gens, roots = (roots,) + gens, None
    try:
        f, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('viete', 1, exc) from None
    if f.is_multivariate:
        raise MultivariatePolynomialError("multivariate polynomials are not allowed")
    n = f.degree()
    if n < 1:
        raise ValueError("Cannot derive Viete's formulas for a constant polynomial")
    if roots is None:
        roots = numbered_symbols('r', start=1)
    roots = list(islice(roots, n))
    if n != len(roots):
        raise ValueError("required %s roots, got %s" % (n, len(roots)))
    lc, coeffs = f.LC(), f.all_coeffs()
    result, sign = [], -1
    for i, coeff in enumerate(coeffs[1:]):
        poly = symmetric_poly(i + 1, roots)
        result.append((poly, sign*(coeff/lc)))
        sign = -sign
    return result
