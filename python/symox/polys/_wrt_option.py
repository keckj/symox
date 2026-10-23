"""The ``wrt`` and ``sort`` generator-ordering options of the polys functions.

The reference ``Options`` only use them to order automatically detected
generators (``_sort_gens``): explicit generators ignore them. The native
functions take explicit generators, so the wrapper detects the generators
of all inputs jointly, orders them like the reference and passes them on.
"""
import functools

# name -> how the polynomial inputs are taken from the positional arguments:
# an int n is the n leading arguments, 'reduced' is (f, G), 'list' is (F,)
_ARITY = {
    'factor_list': 1, 'sqf_list': 1, 'sqf_part': 1, 'discriminant': 1,
    'decompose': 1, 'sturm': 1, 'gff_list': 1, 'ground_roots': 1,
    'div': 2, 'rem': 2, 'quo': 2, 'pdiv': 2, 'prem': 2, 'pquo': 2,
    'exquo': 2, 'pexquo': 2, 'resultant': 2, 'subresultants': 2,
    'compose': 2, 'reduced': 'reduced', 'groebner': 'list',
}

_DOMAIN_KEYS = ('domain', 'extension', 'modulus', 'gaussian', 'greedy',
                'field', 'composite', 'expand')


def _split(name, fargs):
    arity = _ARITY[name]
    if arity == 'reduced':
        if len(fargs) < 2:
            return None
        return [fargs[0]] + list(fargs[1]), fargs[2:]
    if arity == 'list':
        if not fargs:
            return None
        return list(fargs[0]), fargs[1:]
    if len(fargs) < arity:
        return None
    return list(fargs[:arity]), fargs[arity:]


def _is_poly(value):
    from .. import Poly
    return isinstance(value, Poly)


def _ordered_gens(exprs, args, order):
    from .. import sympify, Poly
    from .._poly_constructor import _sort
    from .polytools import parallel_poly_from_expr
    parts = []
    for e in exprs:
        if isinstance(e, Poly):
            return None
        e = sympify(e)
        if getattr(e, 'is_Relational', False):
            e = e.lhs - e.rhs
        n, d = e.as_numer_denom()
        parts.extend((n, d))
    options = {k: args[k] for k in _DOMAIN_KEYS if k in args}
    try:
        _, opt = parallel_poly_from_expr(parts, **options)
    except Exception:
        return None
    gens = tuple(opt.get('gens', ()) if isinstance(opt, dict) else opt.gens)
    if not gens:
        return None
    return _sort(gens, order)


def _sorted_factors(factors, method):
    if method == 'sqf':
        def key(obj):
            poly, exp = obj
            rep = poly.as_list()
            return (exp, len(rep), len(poly.gens), str(poly.domain), rep)
    else:
        def key(obj):
            poly, exp = obj
            rep = poly.as_list()
            return (len(rep), len(poly.gens), exp, str(poly.domain), rep)
    return sorted(factors, key=key)


def _symbolic_factor_list(expr, options, method):
    """Port of the reference ``_symbolic_factor_list`` with automatic,
    ``wrt``/``sort`` ordered generators for each factor."""
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
            poly = Poly(base, **options)
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


def _generic_factor_list(expr, args, method):
    """Port of the reference ``_generic_factor_list`` (automatic gens)."""
    from .. import Poly, sympify, together, Expr
    from .polyerrors import PolynomialError
    frac = args.pop('frac', False)
    polys = args.pop('polys', False)
    expr = sympify(expr)
    if not isinstance(expr, Expr):
        raise PolynomialError("a polynomial expected, got %s" % expr)
    numer, denom = together(expr).as_numer_denom()
    cp, fp = _symbolic_factor_list(numer, args, method)
    cq, fq = _symbolic_factor_list(denom, args, method)
    if fq and not frac:
        raise PolynomialError("a polynomial expected, got %s" % expr)
    _args = dict(args, expand=True)
    for factors in (fp, fq):
        for i, (f, k) in enumerate(factors):
            if not isinstance(f, Poly):
                factors[i] = (Poly(f, **_args), k)
    fp = _sorted_factors(fp, method)
    fq = _sorted_factors(fq, method)
    if not polys:
        fp = [(f.as_expr(), k) for f, k in fp]
        fq = [(f.as_expr(), k) for f, k in fq]
    coeff = cp/cq
    if not frac:
        return coeff, fp
    return coeff, fp, fq


def _wrap(native, name):
    @functools.wraps(native)
    def wrapper(*fargs, **args):
        if args.get('wrt') is None and args.get('sort') is None:
            args.pop('wrt', None)
            args.pop('sort', None)
            return native(*fargs, **args)
        order = {k: args.pop(k) for k in ('wrt', 'sort') if args.get(k) is not None}
        args.pop('wrt', None)
        args.pop('sort', None)
        from .._poly_constructor import _sort
        _sort((), order)  # validates the option values as the reference
        split = _split(name, fargs)
        if (name in ('factor_list', 'sqf_list') and split is not None and not split[1]
                and args.get('gens') is None and not _is_poly(split[0][0])):
            unknown = set(args) - set(_DOMAIN_KEYS) - {'frac', 'polys'}
            if not unknown:
                return _generic_factor_list(split[0][0], dict(args, **order),
                                            name[:-5])
        if split is None or split[1] or args.get('gens') is not None:
            return native(*fargs, **args)
        gens = _ordered_gens(split[0], args, order)
        if gens is None:
            return native(*fargs, **args)
        return native(*fargs, *gens, **args)
    wrapper._symox_wrt_option = True
    return wrapper


def install(*namespaces):
    done = {}
    for namespace in namespaces:
        for name in _ARITY:
            native = namespace.get(name)
            if native is None or getattr(native, '_symox_wrt_option', False):
                continue
            key = (name, id(native))
            if key not in done:
                done[key] = _wrap(native, name)
            namespace[name] = done[key]


# The reference errors of the univariate-only operations on a multivariate
# polynomial (DMP.decompose / DMP.gff_list and Poly.ground_roots).
_UNIVARIATE = {
    'decompose': ('univariate polynomial expected',),
    'gff_list': ('greatest factorial factorization requires one polynomial generator',),
    'ground_roots': ('ground_roots requires one polynomial generator',),
}


def _univariate_error(name, poly):
    from .polyerrors import MultivariatePolynomialError
    if name == 'ground_roots':
        return MultivariatePolynomialError("Cannot compute ground roots of %s" % poly)
    return ValueError('univariate polynomial expected')


def _wrap_univariate(native, name, method):
    messages = _UNIVARIATE[name]

    @functools.wraps(native)
    def wrapper(*fargs, **args):
        try:
            return native(*fargs, **args)
        except (ValueError, NotImplementedError) as error:
            if type(error) not in (ValueError, NotImplementedError) or str(error) not in messages:
                raise
            from .. import Poly
            if method:
                poly = fargs[0]
            else:
                gens = fargs[1:]
                options = {k: v for k, v in args.items() if k not in ('polys', 'multiple')}
                try:
                    poly = fargs[0] if isinstance(fargs[0], Poly) else Poly(fargs[0], *gens, **options)
                except Exception:
                    raise error from None
            if len(poly.gens) <= 1:
                raise
            raise _univariate_error(name, poly) from None
    wrapper._symox_univariate_error = True
    return wrapper


def install_univariate_errors(*namespaces):
    from .. import Poly
    for name in _UNIVARIATE:
        native = Poly.__dict__.get(name)
        if native is not None and not getattr(native, '_symox_univariate_error', False):
            setattr(Poly, name, _wrap_univariate(native, name, True))
    done = {}
    for namespace in namespaces:
        for name in _UNIVARIATE:
            native = namespace.get(name)
            if native is None or getattr(native, '_symox_univariate_error', False):
                continue
            key = (name, id(native))
            if key not in done:
                done[key] = _wrap_univariate(native, name, False)
            namespace[name] = done[key]
