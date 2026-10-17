"""Top-level polynomial helpers (``sympy.polys.polytools``).

The functions follow SymPy's calling conventions: an expression is first
converted with :func:`poly_from_expr`, a constant without generators raises
``ComputationFailed``, and the result is an expression unless the input was
a ``Poly`` or ``polys=True`` was passed.
"""
from .._core import (
    Poly, GroebnerBasis, groebner, reduced, factor, factor_list, sqf_list,
    gcd, lcm, resultant, subresultants, discriminant, sturm, count_roots, intervals, refine_root, ground_roots, gff_list,
)
from .polyerrors import (
    ComputationFailed, FlagError, MultivariatePolynomialError, OptionError,
    PolificationFailed, PolynomialError,
)

__all__ = [
    'Poly', 'GroebnerBasis', 'groebner', 'reduced', 'factor', 'factor_list',
    'sqf_list', 'gcd', 'lcm', 'resultant', 'subresultants', 'discriminant',
    'sturm', 'count_roots', 'intervals', 'refine_root', 'ground_roots', 'gff_list',
    'poly_from_expr', 'parallel_poly_from_expr', 'degree_list', 'total_degree',
    'trunc', 'monic', 'content', 'gff', 'sqf_norm', 'all_roots', 'real_roots',
    'nth_power_roots_poly', 'LC', 'LM', 'LT',
]

_FLAGS = ('auto', 'frac', 'formal', 'polys', 'include', 'all', 'gen', 'series', 'symbols', 'method')
_DEFAULTS = {
    'gens': (), 'domain': None, 'polys': None, 'auto': True, 'expand': True,
    'frac': False, 'formal': False, 'include': False, 'all': False,
    'greedy': None, 'modulus': None, 'extension': None, 'field': None,
    'strict': True, 'split': None, 'gaussian': None, 'wrt': None,
    'sort': None, 'series': False, 'gen': 0, 'method': None, 'symbols': None,
}


class Options(dict):
    """Polynomial options as a dictionary with attribute access.

    Only the options that were given or computed are stored; reading any
    other known option returns SymPy's default.
    """

    def __getattr__(self, name):
        if name in self:
            return self[name]
        if name == 'order':
            from .orderings import lex
            return lex
        if name in _DEFAULTS:
            return _DEFAULTS[name]
        raise AttributeError(name)

    def __setattr__(self, name, value):
        self[name] = value

    def clone(self, updates={}):
        out = Options(self)
        out.update(updates)
        return out


def allowed_flags(args, flags):
    """Raise ``FlagError`` for flags not allowed in this context."""
    for arg in args:
        if arg in _FLAGS and arg not in flags:
            raise FlagError("'%s' flag is not allowed in this context" % arg)


def _build_options(gens, args):
    from .. import sympify
    if len(gens) == 1 and isinstance(gens[0], (list, tuple)):
        gens = tuple(gens[0])
    if gens and args.get('gens'):
        raise OptionError("both '*gens' and keyword argument 'gens' supplied")
    opt = Options()
    for key, value in args.items():
        if key == 'order':
            from .orderings import monomial_key
            value = monomial_key(value)
        elif key == 'gens':
            value = tuple(sympify(g) for g in value)
        opt[key] = value
        if key == 'domain' and 'auto' not in args:
            opt['auto'] = False
    ext = args.get('extension')
    if ext is not None and ext is not True and ext is not False:
        opt['extension'] = set(ext) if isinstance(ext, (list, tuple, set)) else {sympify(ext)}
    if ('modulus' in args or (ext is not None and ext is not True)) and 'domain' not in args:
        opt['domain'] = None  # filled in with the computed domain, as in SymPy
        if ext is not None and 'auto' not in args:
            opt['auto'] = False
    if gens:
        opt['gens'] = tuple(sympify(g) for g in gens)
    return opt


def _poly_kwargs(opt):
    out = {k: opt[k] for k in ('domain', 'modulus', 'extension') if k in opt and opt[k] is not None}
    if isinstance(out.get('extension'), set):
        ext = sorted(out['extension'], key=lambda e: e.sort_key())
        out['extension'] = ext[0] if len(ext) == 1 else ext
    return out


def _without_generators(exc):
    return "without generators" in str(exc) or "without generator" in str(exc)


def _poly_from_expr(expr, opt):
    from .. import sympify
    orig = expr
    if isinstance(expr, Poly):
        poly = Poly(expr, *opt.gens, **_poly_kwargs(opt)) if (opt.gens or _poly_kwargs(opt)) else expr
        opt['gens'] = poly.gens
        opt['domain'] = poly.domain
        if 'polys' not in opt:
            opt['polys'] = True
        return poly, opt
    if isinstance(expr, (list, tuple, set, dict)):
        raise PolificationFailed(opt, orig, expr)
    try:
        expr = sympify(expr)
    except (TypeError, ValueError):
        raise PolificationFailed(opt, orig, expr) from None
    try:
        poly = Poly(expr, *opt.gens, **_poly_kwargs(opt))
    except ValueError as exc:
        if not opt.gens and _without_generators(exc):
            raise PolificationFailed(opt, orig, expr) from None
        raise
    opt['gens'] = poly.gens
    if 'domain' not in opt or opt['domain'] is None:
        try:
            opt['domain'] = poly.domain
        except NotImplementedError:
            opt['domain'] = poly.domain_descriptor
    if 'polys' not in opt:
        opt['polys'] = False
    return poly, opt


def poly_from_expr(expr, *gens, **args):
    """Construct a polynomial from an expression, returning ``(Poly, options)``."""
    return _poly_from_expr(expr, _build_options(gens, args))


def _gens_of(expr):
    try:
        return Poly(expr).gens
    except ValueError as exc:
        if _without_generators(exc):
            return ()
        raise


def _ordered_gens(gens):
    from .. import Add
    if len(gens) <= 1:
        return tuple(gens)
    try:
        ordered = Poly(Add(*gens)).gens
    except Exception:
        ordered = ()
    if len(ordered) == len(gens) and set(ordered) == set(gens):
        return tuple(ordered)
    return tuple(sorted(gens, key=lambda g: g.sort_key()))


def parallel_poly_from_expr(exprs, *gens, **args):
    """Construct polynomials over common generators and a common domain."""
    from .. import sympify
    opt = _build_options(gens, args)
    origs = list(exprs)
    converted, have_polys, failed = [], False, False
    for expr in origs:
        if isinstance(expr, Poly):
            have_polys = True
            converted.append(expr.as_expr())
            continue
        try:
            converted.append(sympify(expr))
        except (TypeError, ValueError):
            failed = True
            converted.append(expr)
    if failed:
        raise PolificationFailed(opt, origs, converted, True)
    kwargs = _poly_kwargs(opt)
    if opt.gens:
        all_gens = tuple(opt.gens)
    else:
        found = []
        for expr in converted:
            for g in _gens_of(expr):
                if g not in found:
                    found.append(g)
        if not found:
            raise PolificationFailed(opt, origs, converted, True)
        all_gens = _ordered_gens(found)
    polys = [Poly(e, *all_gens, **kwargs) for e in converted]
    if not kwargs:
        domain = None
        try:
            for p in polys:
                domain = p.domain if domain is None else domain.unify(p.domain)
        except NotImplementedError:
            domain = None
        if domain is not None:
            polys = [Poly(e, *all_gens, domain=domain) for e in converted]
    opt['gens'] = all_gens
    try:
        opt['domain'] = polys[0].domain
    except NotImplementedError:
        opt['domain'] = polys[0].domain_descriptor
    if 'polys' not in opt:
        opt['polys'] = have_polys
    return polys, opt


def _poly_or_fail(name, f, gens, args):
    try:
        return poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed(name, 1, exc) from None


def degree_list(f, *gens, **args):
    """Degrees of ``f`` in all of its generators."""
    allowed_flags(args, ['polys'])
    F, opt = _poly_or_fail('degree_list', f, gens, args)
    return F.degree_list()


def total_degree(f, *gens):
    """Total degree of ``f`` in ``gens`` (all of its generators by default)."""
    from .. import Integer, sympify
    p = sympify(f) if not isinstance(f, Poly) else f
    if isinstance(p, Poly):
        gens = gens or p.gens
        p = p.as_expr()
    if p.is_Number:
        return Integer(0)
    return Integer(Poly(p, *gens).total_degree())


def trunc(f, p, *gens, **args):
    """Reduce the integer coefficients of ``f`` to symmetric residues modulo ``p``."""
    from .. import sympify
    allowed_flags(args, ['auto', 'polys'])
    F, opt = _poly_or_fail('trunc', f, gens, args)
    result = F.trunc(sympify(p))
    return result if opt.polys else result.as_expr()


def monic(f, *gens, **args):
    """Divide all coefficients of ``f`` by its leading coefficient."""
    allowed_flags(args, ['auto', 'polys'])
    F, opt = _poly_or_fail('monic', f, gens, args)
    if not opt.auto:
        K = F.domain
        if not K.is_Field and F.LC() not in (1, -1) and any(not (c/F.LC()).is_integer for c in F.coeffs()):
            from .polyerrors import ExactQuotientFailed
            raise ExactQuotientFailed("%s does not divide %s in %s" % (F.LC(), F.coeffs()[-1], K))
    result = F.monic()
    return result if opt.polys else result.as_expr()


def content(f, *gens, **args):
    """GCD of the coefficients of ``f``."""
    allowed_flags(args, ['polys'])
    F, opt = _poly_or_fail('content', f, gens, args)
    return F.content()


def gff(f, *gens, **args):
    """Greatest factorial factorization (not implemented, as in SymPy)."""
    raise NotImplementedError('symbolic falling factorial')


def sqf_norm(f, *gens, **args):
    """``(s, g, r)`` with ``r = norm(g)`` square-free and ``g = f(x - s*a)``."""
    from .. import Integer
    allowed_flags(args, ['polys'])
    F, opt = _poly_or_fail('sqf_norm', f, gens, args)
    s, g, r = F.sqf_norm()
    s = [Integer(k) for k in s]
    if not opt.polys:
        return s, g.as_expr(), r.as_expr()
    return s, g, r


def _univariate_for_roots(f, extension):
    from .. import sympify
    kwargs = {'extension': extension} if extension else {}
    try:
        F = Poly(f, **kwargs) if isinstance(f, Poly) else Poly(sympify(f), **kwargs)
    except ValueError as exc:
        if _without_generators(exc):
            raise PolynomialError("Cannot compute real roots of %s, not a polynomial" % f) from None
        raise
    if not F.is_univariate:
        raise PolynomialError("only univariate polynomials are allowed")
    return F


def all_roots(f, multiple=True, radicals=True, extension=False):
    """All complex roots of a univariate polynomial, real roots first."""
    F = _univariate_for_roots(f, extension)
    return F.all_roots(multiple=multiple, radicals=radicals)


def real_roots(f, multiple=True, radicals=True, extension=False):
    """Real roots of a univariate polynomial in increasing order."""
    F = _univariate_for_roots(f, extension)
    return F.real_roots(multiple=multiple, radicals=radicals)


def nth_power_roots_poly(f, n, *gens, **args):
    """Polynomial whose roots are the ``n``-th powers of the roots of ``f``."""
    allowed_flags(args, [])
    F, opt = _poly_or_fail('nth_power_roots_poly', f, gens, args)
    result = F.nth_power_roots_poly(n)
    return result if opt.polys else result.as_expr()


def LC(f, *gens, **args):
    """Leading coefficient of ``f`` with respect to ``order``."""
    allowed_flags(args, ['polys'])
    F, opt = _poly_or_fail('LC', f, gens, args)
    return F.LC(order=opt.order)


def LM(f, *gens, **args):
    """Leading monomial of ``f`` with respect to ``order``."""
    allowed_flags(args, ['polys'])
    F, opt = _poly_or_fail('LM', f, gens, args)
    return F.LM(order=opt.order).as_expr()


def LT(f, *gens, **args):
    """Leading term of ``f`` with respect to ``order``."""
    allowed_flags(args, ['polys'])
    F, opt = _poly_or_fail('LT', f, gens, args)
    monom, coeff = F.LT(order=opt.order)
    return coeff*monom.as_expr()
