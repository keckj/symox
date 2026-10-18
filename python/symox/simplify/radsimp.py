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
    def visit(node):
        if not node.is_Mul:return node
        inside=[];outside=[]
        for factor in node.args:
            if factor.class_name=='Abs':inside.append(factor.args[0])
            elif factor.is_Pow and factor.base.class_name=='Abs' and factor.exp.is_Integer:inside.append(factor.base.args[0]**factor.exp)
            else:outside.append(factor)
        return S.Mul(*outside)*S.Abs(S.Mul(*inside)) if len(inside)>1 else node
    return S.bottom_up(expr,visit)


def collect_sqrt(expr,evaluate=None):
    expr=S.sympify(expr);radicals=set()
    for term in S.Add.make_args(expr):
        for factor in S.Mul.make_args(term):
            if factor==S.I or factor.is_Pow and factor.exp.is_Rational and factor.exp.q==2 and factor.is_number:radicals.add(factor)
    result=S.collect_const(expr,*sorted(radicals,key=S.default_sort_key),Numbers=False) if radicals else expr
    if evaluate is not False:return result
    terms=S.Add.make_args(result)
    count=sum(any(factor==S.I or factor.is_Pow and factor.exp.is_Rational and factor.exp.q==2 for factor in S.Mul.make_args(term)) for term in terms)
    return (tuple(sorted(terms,key=S.default_sort_key)) if count else (result,)),count


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
