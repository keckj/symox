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
         'quo', 'rem', 'gcdex', 'invert', 'cancel')
EUCLID = ('gcdex', 'invert')
_BINARY = {'gcd', 'lcm', 'div', 'quo', 'rem', 'gcdex', 'invert'}
_STRING = re.compile(r'^\s*(?:RR|CC|RR_\d+|CC_\d+)\s*(?:[\[(].*)?$')
_CANCEL_STRING = re.compile(r'^\s*(?:RR|CC|R|C)(?:_\d+)?\s*(?:[\[(].*)?$')
_ALLOWED = {'domain', 'polys', 'auto', 'frac', 'fraction', 'expand'}
_CANCEL_ALLOWED = _ALLOWED | {'_signsimp', 'gens', 'field', 'greedy'}


def _inexact(domain):
    if getattr(domain, 'is_RR', False) or getattr(domain, 'is_CC', False):
        return True
    try:
        ground = getattr(domain, 'domain', None)
    except Exception:
        return False
    return ground is not None and ground is not domain and _inexact(ground)


def selected(fargs, args, name=None):
    """Whether the options select an inexact ground domain."""
    allowed = _CANCEL_ALLOWED if name == 'cancel' else _ALLOWED
    if args.get('modulus') is not None or set(args) - allowed:
        return False
    domain = args.get('domain')
    if domain is None:
        return False
    if isinstance(domain, str):
        return bool((_CANCEL_STRING if name == 'cancel' else _STRING).match(domain))
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


def _cancel_nonpolynomial(f):
    """The reference cancel fallback for Piecewise or noncommutative input."""
    from .. import Expr, Piecewise, cancel, preorder_traversal
    from ..logic.boolalg import BooleanAtom
    from ..utilities.iterables import sift
    if f.is_Add or f.is_Mul:
        c, nc = sift(f.args, lambda x: x.is_commutative is True and not x.has(Piecewise), binary=True)
        return f.func(cancel(f.func(*c)), *[cancel(i) for i in nc])
    replacements = []
    traversal = preorder_traversal(f)
    next(traversal)
    for e in traversal:
        if isinstance(e, BooleanAtom) or not isinstance(e, Expr):
            continue
        try:
            replacements.append((e, cancel(e)))
            traversal.skip()
        except NotImplementedError:
            pass
    return f.xreplace(dict(replacements))


def _cancel(fargs, args):
    """Reference cancel protocol over the selected inexact polynomial ring.

    The Poly cancellation helper follows ``dmp_cancel(include=True)`` at
    the selected precision, including composite grounds. Tuple inputs keep
    the reference (1, numerator, denominator) return convention.
    """
    from .. import Expr, Piecewise, Poly, S, Tuple, factor_terms, signsimp, sympify
    from .._noncommutative import _poly_cancel
    from .polytools import allowed_flags, parallel_poly_from_expr
    from .polyerrors import GeneratorsNeeded, OptionError
    from .._poly_constructor import _options, _real_coefficients, _sort, _unexpanded

    args = dict(args)
    do_signsimp = args.pop('_signsimp', True)
    allowed_flags(args, ['polys'])
    original, positional_gens = fargs[0], fargs[1:]
    # sympify embeds native Poly objects as expressions, so retain this
    # public return-type information before normalizing the input tuple.
    input_polys = (isinstance(original, (tuple, Tuple)) and len(original) == 2
                   and all(isinstance(p, Poly) for p in original))
    output_gens = original[0].gens if input_polys else None
    output_domain = original[0].domain if input_polys else None
    if isinstance(original, Poly):
        return original
    f = sympify(original)
    if do_signsimp:
        f = signsimp(f)
    pair = isinstance(f, (tuple, Tuple))
    if pair:
        if len(f) != 2:
            raise ValueError('unexpected argument: %s' % f)
        p, q = f
        polys = args.get('polys', input_polys)
        p, q = p.as_expr(), q.as_expr()
        if (f.has(Piecewise) or p.is_commutative is False or q.is_commutative is False):
            return _cancel_nonpolynomial(f)
    else:
        if f.is_Number or f.is_Relational or not isinstance(f, Expr):
            return f
        f = factor_terms(f, radical=True)
        if f.is_commutative is False or f.has(Piecewise):
            return _cancel_nonpolynomial(f)
        p, q = f.as_numer_denom()
        polys = False
    gens, args = _options(positional_gens, None, None, None, args)
    domain = args['domain']
    if isinstance(domain, str):
        from .. import ComplexField, Dummy, RealField
        numeric = re.fullmatch(r'\s*(R|RR|C|CC)_(\d+)\s*', domain)
        if numeric:
            domain = (ComplexField if numeric[1].startswith('C') else RealField)(int(numeric[2]))
        else:
            domain = Poly(S.One, Dummy('g'), domain=domain).domain
        args['domain'] = domain
    options = {key: value for key, value in args.items() if key not in {'polys', 'gens'}}
    if not gens:
        # Generator inference is relative to the *selected* domain. Numeric
        # constants such as sqrt(2), sin(1), or I over CC are coefficients,
        # rather than generators inferred over the default exact ground.
        from .. import expand
        found, parts = set(), []
        for part in (p, q):
            if args.get('expand', True):
                part = expand(part)
            if domain.is_RR:
                part = _real_coefficients(part, (), args)
            parts.append(part)
            try:
                found.update(_unexpanded(part, (), domain, None, None, args).gens)
            except GeneratorsNeeded:
                pass
        gens = _sort(found, args)
        if not gens:
            return (S.One, p, q) if pair else f.expand()
        p, q = parts
    (F, G), _ = parallel_poly_from_expr((p, q), *gens, **options)
    P, Q = _poly_cancel(F, G)
    if not pair:
        return P/Q
    if polys:
        if positional_gens:
            raise OptionError("both '*gens' and keyword argument 'gens' supplied")
        output_options = {} if output_domain is None else {'domain': output_domain}
        return (S.One, Poly(P, *(output_gens or F.gens), **output_options),
                Poly(Q, *(output_gens or G.gens), **output_options))
    return S.One, P, Q


def call(name, fargs, args):
    if name == 'cancel':
        return _cancel(fargs, args)
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
