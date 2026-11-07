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
    'Poly', 'PurePoly', 'GroebnerBasis', 'groebner', 'reduced', 'factor', 'factor_list',
    'sqf_list', 'gcd', 'lcm', 'resultant', 'subresultants', 'discriminant',
    'sturm', 'count_roots', 'intervals', 'refine_root', 'ground_roots', 'gff_list',
    'poly_from_expr', 'parallel_poly_from_expr', 'degree_list', 'total_degree',
    'trunc', 'monic', 'content', 'gff', 'sqf_norm', 'all_roots', 'real_roots',
    'nth_power_roots_poly', 'LC', 'LM', 'LT', 'decompose', 'terms_gcd',
]


# Remaining reference polytools functions, implemented by the public layer.
_PUBLIC_FUNCTIONS = frozenset((
    'degree', 'pdiv', 'prem', 'pquo', 'pexquo', 'div', 'rem', 'quo', 'exquo',
    'half_gcdex', 'gcdex', 'invert', 'cofactors', 'gcd_list', 'lcm_list',
    'primitive', 'compose', 'sqf_part', 'sqf', 'nroots', 'cancel',
    'is_zero_dimensional', 'poly',
))
__all__ += sorted(_PUBLIC_FUNCTIONS)


def __getattr__(name):
    if name=='PurePoly' or name in _PUBLIC_FUNCTIONS:
        import symox
        return getattr(symox, name)
    raise AttributeError(name)

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
    # The reference stores the given options (positional gens appended to
    # them) before the postprocess steps add the derived ones like auto.
    for key, value in args.items():
        if key == 'order':
            from .orderings import monomial_key
            value = monomial_key(value)
        elif key == 'gens':
            value = tuple(sympify(g) for g in value)
        opt[key] = value
    if gens:
        opt['gens'] = tuple(sympify(g) for g in gens)
    if ('domain' in args or 'field' in args) and 'auto' not in args:
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
    if 'expand' in opt and opt['expand'] is False:
        out['expand'] = False
    # The generator ordering options select the main generator (horner wrt).
    for key in ('wrt', 'sort'):
        if opt.get(key):
            out[key] = opt[key]
    if opt.get('field') and 'domain' not in out:
        out['field'] = True
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
            try:
                opt['domain'] = poly.domain_descriptor()
            except NotImplementedError:
                # No typed descriptor (an algebraic extension such as
                # QQ<I>); callers must use the polynomial itself.
                opt['domain'] = None
    if 'polys' not in opt:
        opt['polys'] = False
    return poly, opt


def poly_from_expr(expr, *gens, **args):
    """Construct a polynomial from an expression, returning ``(Poly, options)``."""
    return _poly_from_expr(expr, _build_options(gens, args))


def _gens_of(expr, expand=True):
    try:
        # With expand=False the reference keeps unexpanded factors as
        # generators (``_parallel_dict_from_expr_no_gens``).
        return Poly(expr, expand=False).gens if expand is False else Poly(expr).gens
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


def _joint_domain(exprs, gens, kwargs):
    """The domain inferred from the coefficients of all ``exprs`` together.

    A fresh generator separates the expressions in one carrier polynomial,
    whose coefficients in ``gens`` are exactly their union.
    """
    from .. import Add, Dummy
    t = Dummy('t')
    try:
        return Poly(Add(*[e*t**i for i, e in enumerate(exprs)]), t, *gens, **kwargs).domain
    except (NotImplementedError, ValueError, TypeError):
        return None


def _separate_domains_joint(polys, exprs, gens):
    """The polynomials over one domain built from all coefficients.

    Only needed when the separately inferred domains differ: the reference
    calls ``construct_domain`` once on the union of the coefficients, so
    ``EX`` with ``ZZ[y]`` (shared symbols) gives ``EX`` and ``RR`` with
    ``ZZ[y]`` gives ``RR[y]``.
    """
    from .constructor import construct_domain
    try:
        domains = {str(p.domain) for p in polys}
    except NotImplementedError:
        domains = None
    if domains is not None and len(domains) == 1:
        return None
    try:
        coefficients = [c for e in exprs for c in Poly(e, *gens, domain='EX').coeffs()]
        domain, _ = construct_domain(coefficients)
        return [Poly(e, *gens, domain=domain) for e in exprs]
    except Exception:  # any failure keeps the separately inferred domains
        return None


def _domains_differ(polys):
    try:
        return len({str(p.domain) for p in polys}) > 1
    except NotImplementedError:
        return False


def _field_joint(polys, exprs, gens, kwargs):
    from .constructor import construct_domain
    rest = {k: v for k, v in kwargs.items() if k != 'field'}
    try:
        coefficients = [c for e in exprs for c in Poly(e, *gens, domain='EX', **rest).coeffs()]
        domain, _ = construct_domain(coefficients, field=True)
        return [Poly(e, *gens, domain=domain, **rest) for e in exprs]
    except Exception:  # any failure keeps the separately inferred domains
        return polys


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
            for g in _gens_of(expr, kwargs.get('expand', True)):
                if g not in found:
                    found.append(g)
        if not found:
            raise PolificationFailed(opt, origs, converted, True)
        all_gens = _ordered_gens(found)
    polys = [Poly(e, *all_gens, **kwargs) for e in converted]
    if len(polys) > 1 and kwargs.get('extension') is not True and kwargs.get('field') is True \
            and _domains_differ(polys):
        # With ``field=True`` the reference also builds one domain from the
        # union of the coefficients: ``[0, x]`` in ``t`` is over ``ZZ(x)``,
        # ``[1/2, x]`` too, where separate domains give ``QQ`` and ``QQ(x)``.
        polys = _field_joint(polys, converted, all_gens, kwargs)
    if kwargs.get('extension') is True and len(polys) > 1:
        # The reference infers one domain from the union of all
        # coefficients; separately inferred domains may not unify to it
        # (QQ<sqrt(2)> and ZZ[y] jointly give EX in x).
        joint = _joint_domain(converted, all_gens, kwargs)
        if joint is not None:
            rest = {k: v for k, v in kwargs.items() if k != 'extension'}
            try:
                polys = [Poly(e, *all_gens, domain=joint, **rest) for e in converted]
            except (NotImplementedError, ValueError, TypeError):
                pass
    joint = None
    if not kwargs and len(polys) > 1:
        joint = _separate_domains_joint(polys, converted, all_gens)
        if joint is not None:
            polys = joint
    if not kwargs and joint is None:
        domain = None
        try:
            for p in polys:
                domain = p.domain if domain is None else domain.unify(p.domain)
        except NotImplementedError:
            domain = None
        if domain is not None:
            polys = [Poly(e, *all_gens, domain=domain) for e in converted]
    opt['gens'] = all_gens
    if not polys:
        # The reference reports the default ground domain of no polynomials.
        from .domains import QQ, ZZ
        opt['domain'] = QQ if opt.get('field') else ZZ
    else:
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


def decompose(f, *gens, **args):
    """Exact univariate functional decomposition in the selected coefficient domain."""
    allowed_flags(args, ['polys'])
    F, opt = _poly_or_fail('decompose', f, gens, args)
    result = F.decompose()
    return result if opt.polys else [part.as_expr() for part in result]


def total_degree(f, *gens):
    """Total degree of ``f`` in ``gens`` (all of its generators by default)."""
    from .. import Integer, sympify
    p = sympify(f) if not isinstance(f, Poly) else f
    if isinstance(p, Poly):
        if not gens or tuple(gens) == p.gens:
            return Integer(p.total_degree())
        p = p.as_expr()
    if p.is_Number:
        return Integer(0)
    return Integer(Poly(p, *gens).total_degree())


def terms_gcd(f, *gens, **args):
    """Remove the GCD of the terms of ``f``, SymPy's ``terms_gcd``.

    ``deep`` applies it to every subexpression first and ``clear=False``
    keeps a fractional coefficient inside a sum (``x*(1 - x**3/12)``).
    Adapted from SymPy 1.14's sympy/polys/polytools.py.
    """
    from .. import sympify, Equality, Relational, Expr, Mul, S
    from ..core.mul import _keep_coeff
    from ..core.numbers import equal_valued
    orig = sympify(f)
    if isinstance(f, Equality):
        return Equality(*(terms_gcd(s, *gens, **args) for s in [f.lhs, f.rhs]))
    elif isinstance(f, Relational):
        raise TypeError("Inequalities cannot be used with terms_gcd. Found: %s" % (f,))
    f = orig
    if not isinstance(f, Expr) or f.is_Atom:
        return orig
    if args.get('deep', False):
        new = f.func(*[terms_gcd(a, *gens, **args) for a in f.args])
        args.pop('deep')
        args['expand'] = False
        return terms_gcd(new, *gens, **args)
    clear = args.pop('clear', True)
    allowed_flags(args, ['polys'])
    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        return exc.expr
    J, f = F.terms_gcd()
    gaussian = _gaussian_primitive(f) if opt.domain is None else None
    if gaussian is not None:
        coeff, expr = gaussian
    elif opt.domain is not None and opt.domain.is_Exact is False:
        # RR and CC: their gcd and denominators are one.
        coeff = S.One
        expr = f.as_expr()
    elif opt.domain is None or opt.domain.is_Ring:
        if opt.domain is not None and opt.domain.is_Field:
            denom, f = f.clear_denoms(convert=True)
        coeff, f = f.primitive()
        if opt.domain is not None and opt.domain.is_Field:
            coeff /= denom
        expr = f.as_expr()
    else:
        coeff = S.One
        expr = f.as_expr()
    term = Mul(*[x**j for x, j in zip(f.gens, J)])
    if equal_valued(coeff, 1):
        coeff = S.One
        if term == 1:
            return orig
    if clear:
        return _keep_coeff(coeff, term*expr)
    # base the clearing on the form of the original expression, not
    # the (perhaps) Mul that we have now
    coeff, f = _keep_coeff(coeff, expr, clear=False).as_coeff_Mul()
    return _keep_coeff(coeff, term*f, clear=False)


def _gaussian_primitive(f):
    """SymPy's QQ_I ``clear_denoms(convert=True)`` and ZZ_I ``primitive``.

    Returns ``(content, primitive expression)``, or ``None`` when a
    coefficient is not a Gaussian rational.
    """
    from .. import Add, Mul, ZZ_I, ilcm
    terms = f.terms()
    parts = []
    for _, c in terms:
        re, im = c.as_real_imag()
        if not (re.is_Rational and im.is_Rational):
            return None
        parts.append((re, im))
    denom = 1
    for re, im in parts:
        denom = ilcm(denom, re.q, im.q)
    elements = [ZZ_I(int(re*denom), int(im*denom)) for re, im in parts]
    content = ZZ_I.zero
    for c in elements:
        content = ZZ_I.gcd(content, c)
        if content == ZZ_I.one:
            break
    if content == ZZ_I.zero:
        return None
    expr = Add(*[ZZ_I.to_sympy(ZZ_I.exquo(c, content))*Mul(*[g**e for g, e in zip(f.gens, monom)])
                 for (monom, _), c in zip(terms, elements)])
    return ZZ_I.to_sympy(content)/denom, expr


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
    result = F.monic(auto=opt.auto)
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
    # SymPy builds Poly(f, greedy=False): algebraic or transcendental
    # coefficients stay in the ground domain (EX, ZZ[pi]) instead of becoming
    # generators, so sorting the roots fails with NotImplementedError.
    # With ``extension`` set the reference rebuilds a Poly over a
    # non-algebraic domain from its expression with ``extension=True``,
    # whatever extension was given: the natural generators of
    # ``Poly(x - sqrt(2))`` are ``(x, sqrt(2))``.
    kwargs = {'extension': True} if extension else {'greedy': False}
    try:
        if isinstance(f, Poly):
            F = Poly(f.as_expr(), **kwargs) if extension and not f.domain.is_AlgebraicField else f
        else:
            F = Poly(sympify(f), **kwargs)
    except ValueError as exc:
        if _without_generators(exc):
            raise PolynomialError("Cannot compute real roots of %s, not a polynomial" % f) from None
        raise
    if not isinstance(f, Poly) and F.is_univariate and not F.gen.is_Symbol:
        raise PolynomialError("generator must be a Symbol")
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
