"""Expression-level polynomial functions over the ring Z/nZ.

When the options select a ``GF(n)`` outside the native prime field path
(composite or large modulus, ``symmetric=False``) or an argument is a
``ModularPoly``, the reference ``polytools`` functions build the
polynomials with ``parallel_poly_from_expr`` and call the ``Poly`` method;
this module does the same with ``_modular_poly.ModularPoly``.
"""
from .._modular_poly import ModularPoly, ring_domain

_BINARY = {'gcd', 'lcm', 'cofactors', 'div', 'rem', 'quo', 'exquo', 'prem', 'pquo', 'pexquo',
           'pdiv', 'half_gcdex', 'gcdex', 'invert', 'subresultants', 'resultant', 'compose'}
_UNARY = {'factor_list', 'sqf_list', 'factor', 'sqf', 'sqf_part', 'discriminant', 'decompose'}
NAMES = tuple(sorted(_BINARY | _UNARY | {'groebner', 'reduced', 'roots', 'real_roots', 'nroots',
                                         'gcd_list', 'lcm_list'}))


def selected(fargs, args):
    """The ``ModularRing`` the call runs over, or None for the native path."""
    for arg in fargs:
        if isinstance(arg, ModularPoly):
            return arg.domain
        if isinstance(arg, (list, tuple)):
            for item in arg:
                if isinstance(item, ModularPoly):
                    return item.domain
    return ring_domain(args.get('domain'), args.get('modulus'), args.get('symmetric', True))


def _polys(exprs, gens, dom):
    from .. import parallel_poly_from_expr
    plain = [e.as_expr() if isinstance(e, ModularPoly) else e for e in exprs]
    if not gens:
        for e in exprs:
            if isinstance(e, ModularPoly):
                gens = e.gens
                break
    polys, opt = parallel_poly_from_expr(plain, *gens)
    return [ModularPoly.from_expr(p, opt['gens'], dom) for p in polys]


def _out(value, polys):
    if isinstance(value, ModularPoly):
        return value if polys else value.as_expr()
    if isinstance(value, tuple):
        return tuple(_out(v, polys) for v in value)
    if isinstance(value, list):
        return [_out(v, polys) for v in value]
    return value


def _seeded(name):
    # The reference equal-degree factorization draws random polynomials, so
    # over a composite modulus it may or may not meet a zero divisor. Symox
    # uses a fixed seed so the outcome is reproducible.
    from ._dense import _support
    _support._rng.seed(0)


def call(name, dom, fargs, args):
    from .. import Mul, sympify
    from ..core.mul import _keep_coeff
    args = dict(args)
    polys = args.pop('polys', None)
    for key in ('modulus', 'domain', 'symmetric', 'auto', 'all', 'include', 'frac'):
        args.pop(key, None)
    if name in ('roots', 'real_roots', 'nroots'):
        raise NotImplementedError('%s of polynomials over %s is not implemented' % (name, dom))
    if name in ('groebner', 'reduced'):
        if name == 'reduced' and not dom.is_Field:
            from .polyerrors import DomainError
            raise DomainError('there is no ring associated with %s' % dom)
        # The reference runs Buchberger or F5B over the ring with ``//`` on
        # the coefficients: it returns a basis when every division is by one
        # and raises TypeError otherwise. Symox has no such ring algorithm;
        # it only reproduces the case decided by the first ``monic`` call.
        if not dom.is_Field and name == 'groebner':
            seq = fargs[0] if fargs and isinstance(fargs[0], (list, tuple)) else ()
            ps = [p for p in _polys(list(seq), fargs[1:], dom) if not p.is_zero] if seq else []
            lex = str(args.get('order', 'lex')) == 'lex' or (ps and len(ps[0].gens) == 1)
            if len(ps) == 1 and lex and not ps[0].is_monic:
                t = dom.dtype.__name__
                raise TypeError("unsupported operand type(s) for //: '%s' and '%s'" % (t, t))
        raise NotImplementedError('%s over %s is not implemented' % (name, dom))
    if name in ('gcd_list', 'lcm_list'):
        seq = list(fargs[0])
        gens = fargs[1:]
        if polys is None:
            polys = bool(seq) and all(isinstance(p, ModularPoly) for p in seq)
        ps = _polys(seq, gens, dom)
        result = ps[0]
        for p in ps[1:]:
            result = result.gcd(p) if name == 'gcd_list' else result.lcm(p)
        return _out(result, polys)
    count = 2 if name in _BINARY else 1
    exprs, gens = fargs[:count], fargs[count:]
    if polys is None:
        polys = all(isinstance(e, ModularPoly) for e in exprs)
    if name in ('factor', 'sqf', 'factor_list', 'sqf_list'):
        _seeded(name)
        method = 'factor_list' if name.startswith('factor') else 'sqf_list'
        if isinstance(exprs[0], ModularPoly):
            coeff, factors = getattr(exprs[0], method)()
        else:
            coeff, factors = _symbolic_list(sympify(exprs[0]), gens, dom, method)
        if name.endswith('_list'):
            return coeff, [(_out(f, polys), k) for f, k in factors]
        return _keep_coeff(sympify(coeff), Mul(*[f.as_expr()**k for f, k in factors]))
    ps = _polys(exprs, gens, dom)
    method = getattr(ps[0], name)
    result = method(*ps[1:])
    return _out(result, polys)


def _symbolic_list(expr, gens, dom, method):
    """The reference ``_symbolic_factor_list`` for one polynomial expression."""
    from .. import Mul, S, together
    coeff, factors = S.One, []
    for arg in Mul.make_args(together(expr)):
        if arg.is_Number:
            coeff *= arg
            continue
        base, exp = arg, S.One
        if arg.is_Pow and arg.exp.is_Integer and arg.exp > 0:
            base, exp = arg.base, arg.exp
        poly = _polys([base], gens, dom)[0]
        c, fs = getattr(poly, method)()
        if c is not S.One and c != 1:
            coeff *= c**exp
        factors.extend((f, k*exp) for f, k in fs)
    return coeff, factors
