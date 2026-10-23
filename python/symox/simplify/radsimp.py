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


def rad_rationalize(num,den):return S.fraction(S.radsimp(S.sympify(num)/S.sympify(den)))


def split_surds(expr):
    expr=S.sympify(expr);terms=S.Add.make_args(expr);radicands=[]
    for term in terms:
        coefficient,radical=term.as_coeff_Mul()
        if radical.is_Pow and radical.exp==S.Rational(1,2) and radical.base.is_Rational and radical.base.is_positive:radicands.append(radical.base)
    if not radicands:raise ValueError('surds must include a square root of a positive rational')
    g=min(radicands);selected=[];rest=[]
    for term in terms:
        coefficient,radical=term.as_coeff_Mul()
        if radical.is_Pow and radical.exp==S.Rational(1,2) and radical.base.is_Rational and (radical.base/g).is_Integer:
            selected.append(coefficient*S.sqrt(radical.base/g))
        else:rest.append(term)
    return g,S.Add(*selected),S.Add(*rest)

__all__ += ['collect_abs','collect_sqrt','denom_expand','expand_denom','expand_fraction','expand_numer','fraction_expand',
 'numer_expand','rad_rationalize','split_surds']
