"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'collect collect_const denom fraction numer radsimp rcollect')

import symox as S


def fraction_expand(expr,**hints):
    numerator,denominator=S.fraction(expr)
    return S.expand(numerator,**hints)/S.expand(denominator,**hints)
def numer_expand(expr,**hints):
    numerator,denominator=S.fraction(expr)
    return S.expand(numerator,**hints)/denominator
def denom_expand(expr,**hints):
    numerator,denominator=S.fraction(expr)
    return numerator/S.expand(denominator,**hints)
expand_fraction=fraction_expand
expand_numer=numer_expand
expand_denom=denom_expand


def collect_abs(expr):
    # Port of the reference collect_abs.
    def _abs(mul):
        c, nc = mul.args_cnc()
        a = []
        o = []
        for i in c:
            if i.class_name == 'Abs':
                a.append(i.args[0])
            elif i.is_Pow and i.base.class_name == 'Abs' and i.exp.is_real:
                a.append(i.base.args[0]**i.exp)
            else:
                o.append(i)
        if len(a) < 2 and not any(i.exp.is_negative for i in a if i.is_Pow):
            return mul
        absarg = S.Mul(*a)
        A = S.Abs(absarg)
        args = [A]
        args.extend(o)
        if not A.has(S.Abs):
            args.extend(nc)
            return S.Mul(*args)
        if A.class_name != 'Abs':
            A = S.Abs(absarg, evaluate=False)
        args[0] = A
        args.extend(nc)
        return S.Mul(*args, evaluate=False)

    return S.sympify(expr).replace(
        lambda x: x.is_Mul,
        lambda x: _abs(x)).replace(
            lambda x: x.is_Pow,
            lambda x: _abs(x))


def collect_sqrt(expr, evaluate=None):
    # Port of the reference collect_sqrt; the radicals are passed to
    # collect_const in sorted order instead of hash-set order.
    if evaluate is None:
        from .._public_core import _evaluation
        evaluate = _evaluation.get()
    expr = S.sympify(expr)
    coeff, expr = expr.as_content_primitive()
    vars = set()
    for a in S.Add.make_args(expr):
        for m in a.args_cnc()[0]:
            if m.is_number and (
                    m.is_Pow and m.exp.is_Rational and m.exp.q == 2 or
                    m == S.I):
                vars.add(m)
    d = S.collect_const(expr, *sorted(vars, key=S.default_sort_key), Numbers=False)
    hit = expr != d
    if not evaluate:
        nrad = 0
        args = list(S.ordered(S.Add.make_args(d)))
        for i, m in enumerate(args):
            c, nc = m.args_cnc()
            for ci in c:
                if ci.is_Pow and ci.exp.is_Rational and ci.exp.q == 2 or \
                        ci == S.I:
                    nrad += 1
                    break
            args[i] *= coeff
        if not (hit or nrad):
            args = [S.Add(*args)]
        return tuple(args), nrad
    return coeff*d


def rad_rationalize(num, den):
    """SymPy's ``rad_rationalize``: multiply by the conjugate of the
    ``split_surds`` decomposition until the denominator is not an Add."""
    from .sqrtdenest import _mexpand
    num, den = S.sympify(num), S.sympify(den)
    if not den.is_Add:
        return num, den
    g, a, b = split_surds(den)
    a = a*S.sqrt(g)
    num = _mexpand((a - b)*num)
    den = _mexpand(a**2 - b**2)
    return rad_rationalize(num, den)


def _split_gcd(*a):
    """SymPy's ``_split_gcd``."""
    g = a[0]
    b1 = [g]
    b2 = []
    for x in a[1:]:
        g1 = S.gcd(g, x)
        if g1 == 1:
            b2.append(x)
        else:
            g = g1
            b1.append(x)
    return g, b1, b2


def split_surds(expr):
    """SymPy's ``split_surds``."""
    from ..core.sorting import default_sort_key
    expr = S.sympify(expr)
    args = sorted(expr.args, key=default_sort_key)
    coeff_muls = [x.as_coeff_Mul() for x in args]
    surds = [x[1]**2 for x in coeff_muls if x[1].is_Pow]
    if not surds:
        raise ValueError('surds must include a square root of a positive rational')
    surds.sort(key=default_sort_key)
    g, b1, b2 = _split_gcd(*surds)
    g2 = g
    if not b2 and len(b1) >= 2:
        b1n = [x/g for x in b1]
        b1n = [x for x in b1n if x != 1]
        # only a common factor has been factored; split again
        g1, b1n, b2 = _split_gcd(*b1n)
        g2 = g*g1
    a1v, a2v = [], []
    for c, s in coeff_muls:
        if s.is_Pow and s.exp == S.S.Half:
            s1 = s.base
            if s1 in b1:
                a1v.append(c*S.sqrt(s1/g2))
            else:
                a2v.append(c*s)
        else:
            a2v.append(c*s)
    return g2, S.Add(*a1v), S.Add(*a2v)

__all__ += ['collect_abs','collect_sqrt','denom_expand','expand_denom','expand_fraction','expand_numer','fraction_expand',
 'numer_expand','rad_rationalize','split_surds']
