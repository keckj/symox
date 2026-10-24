"""Expression-level polynomial functions over ``RR``, ``CC`` and their composites.

The reference ``polytools`` functions build the polynomials with
``parallel_poly_from_expr`` (or ``_poly_from_expr`` for each factor in
``_symbolic_factor_list``) and call the ``Poly`` method. The native
expression functions only cover exact coefficient domains and report an
explicit inexact ``domain=`` as unsupported, so the call is redone here the
same way with the ``Poly`` methods, which do support these domains.
"""
import re

NAMES = ('factor', 'factor_list', 'sqf', 'sqf_list', 'gcd', 'lcm', 'div',
         'quo', 'rem', 'gcdex', 'invert')
EUCLID = ('gcdex', 'invert')
_BINARY = {'gcd', 'lcm', 'div', 'quo', 'rem', 'gcdex', 'invert'}
_STRING = re.compile(r'^\s*(?:RR|CC|RR_\d+|CC_\d+)\s*(?:[\[(].*)?$')
_ALLOWED = {'domain', 'polys', 'auto', 'frac', 'fraction', 'expand'}


def _inexact(domain):
    if getattr(domain, 'is_RR', False) or getattr(domain, 'is_CC', False):
        return True
    try:
        ground = getattr(domain, 'domain', None)
    except Exception:
        return False
    return ground is not None and ground is not domain and _inexact(ground)


def selected(fargs, args):
    """Whether the options select an inexact ground domain."""
    if args.get('modulus') is not None or set(args) - _ALLOWED:
        return False
    domain = args.get('domain')
    if domain is None:
        return False
    if isinstance(domain, str):
        return bool(_STRING.match(domain))
    return _inexact(domain)


def _factor_list(expr, gens, options, method):
    """The reference ``_symbolic_factor_list`` with explicit generators."""
    from functools import reduce
    from operator import mul
    from collections import defaultdict
    from .. import Mul, Poly, S, Expr
    from ..core.evalf import pure_complex
    from .polyerrors import PolificationFailed, GeneratorsNeeded
    coeff, factors = S.One, []
    args = [i._eval_factor() if hasattr(i, '_eval_factor') else i
            for i in Mul.make_args(expr)]
    for arg in args:
        if arg.is_Number or (isinstance(arg, Expr) and pure_complex(arg)):
            coeff *= arg
            continue
        elif arg.is_Pow and arg.base != S.Exp1:
            base, exp = arg.args
            if base.is_Number and exp.is_Number:
                coeff *= arg
                continue
            if base.is_Number:
                factors.append((base, exp))
                continue
        else:
            base, exp = arg, S.One
        try:
            poly = Poly(base, *gens, **options)
        except (PolificationFailed, GeneratorsNeeded):
            factors.append((base, exp))
            continue
        _coeff, _factors = getattr(poly, method + '_list')()
        if _coeff is not S.One:
            if exp.is_Integer:
                coeff *= _coeff**exp
            elif _coeff.is_positive:
                factors.append((_coeff, exp))
            else:
                _factors.append((_coeff, S.One))
        if exp is S.One:
            factors.extend(_factors)
        elif exp.is_integer:
            factors.extend([(f, k*exp) for f, k in _factors])
        else:
            other = []
            for f, k in _factors:
                if f.as_expr().is_positive:
                    factors.append((f, k*exp))
                else:
                    other.append((f, k))
            factors.append((Mul(*[f.as_expr()**k for f, k in other]), exp))
    if method == 'sqf':
        factors = [(reduce(mul, (f for f, _ in factors if _ == k)), k)
                   for k in {i for _, i in factors}]
    rv = defaultdict(int)
    for k, v in factors:
        rv[k] += v
    return coeff, list(rv.items())


def call(name, fargs, args):
    from .. import Mul, Poly, sympify, together
    from ..core.mul import _keep_coeff
    from .polyerrors import PolynomialError
    from ._wrt_option import _sorted_factors
    args = dict(args)
    polys = args.pop('polys', None)
    auto = args.pop('auto', True)
    frac = args.pop('frac', False)
    fraction = args.pop('fraction', True)
    options = {'domain': args['domain']} if 'domain' in args else {}
    if name in ('factor', 'sqf', 'factor_list', 'sqf_list'):
        expr, gens = sympify(fargs[0]), fargs[1:]
        method = 'sqf' if name.startswith('sqf') else 'factor'
        if name in ('factor', 'sqf'):
            coeff, factors = _factor_list(together(expr, fraction=fraction), gens, options, method)
            return _keep_coeff(coeff, Mul(*[f.as_expr()**k for f, k in factors]))
        numer, denom = together(expr).as_numer_denom()
        cp, fp = _factor_list(numer, gens, options, method)
        cq, fq = _factor_list(denom, gens, options, method)
        if fq and not frac:
            raise PolynomialError("a polynomial expected, got %s" % expr)
        for fs in (fp, fq):
            for i, (f, k) in enumerate(fs):
                if not isinstance(f, Poly):
                    fs[i] = (Poly(f, *gens, expand=True, **options), k)
        fp, fq = _sorted_factors(fp, method), _sorted_factors(fq, method)
        if not polys:
            fp = [(f.as_expr(), k) for f, k in fp]
            fq = [(f.as_expr(), k) for f, k in fq]
        if not frac:
            return cp/cq, fp
        return cp/cq, fp, fq
    from .polytools import parallel_poly_from_expr
    exprs, gens = fargs[:2], fargs[2:]
    if not gens and 'domain' not in options:
        # Only symbols are generators: `I` is a coefficient of CC.
        from .. import default_sort_key
        free = set().union(*(sympify(e).free_symbols for e in exprs))
        gens = tuple(sorted(free, key=default_sort_key))
    if polys is None:
        polys = all(isinstance(e, Poly) for e in exprs)
    (F, G), _ = parallel_poly_from_expr(list(exprs), *gens, **options)
    if name in ('gcdex', 'invert'):
        result = getattr(F, name)(G, auto=auto)
    else:
        result = getattr(F, name)(G)
    if polys:
        return result
    if isinstance(result, tuple):
        return tuple(r.as_expr() for r in result)
    return result.as_expr()
