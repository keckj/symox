"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'exptrigsimp trigsimp')
_export(globals(), 'symox._simptools', 'futrig')

import symox as S


def trigsimp_old(expr,*,first=True,**opts):
    from .._simptools import futrig
    return futrig(S.sympify(expr))


def _common_angles(expr):
    """Expand commensurate angles in one auxiliary variable per exact base."""
    from math import gcd, lcm
    groups={}
    for node in S.preorder_traversal(expr):
        if node.class_name not in ('sin','cos','tan','cot','sec','csc','sinh','cosh'):
            continue
        coefficient,base=node.args[0].as_coeff_Mul(rational=True)
        if coefficient.is_Rational and base!=1:
            groups.setdefault(base,set()).add(coefficient)
    replacements={};back={}
    for base,coefficients in groups.items():
        if len(coefficients)<2:continue
        denominator=lcm(*(int(c.q) for c in coefficients))
        numerator=0
        for c in coefficients:numerator=gcd(numerator,abs(int(c.p)*(denominator//int(c.q))))
        unit=S.Rational(numerator,denominator)
        if any(abs(int(c/unit))>16 for c in coefficients):continue
        local=S.Dummy('common_angle');back[local]=unit*base
        for node in S.preorder_traversal(expr):
            if node.class_name in ('sin','cos','tan','cot','sec','csc','sinh','cosh'):
                c,b=node.args[0].as_coeff_Mul(rational=True)
                if b==base and c in coefficients:
                    replacements[node]=node.func((c/unit)*local)
    if not replacements:return None
    return S.expand_trig(expr.xreplace(replacements)).xreplace(back)


def trigsimp_groebner(expr,hints=None,quick=False,order='grlex',polynomial=False):
    """Reduce rational trig expressions modulo their quadratic identities."""
    expr=S.sympify(expr);common=_common_angles(expr)
    expanded=S.expand_trig(expr) if common is None else common
    def ratios(node):
        name=getattr(node,'class_name','')
        if name in ('tan','cot','sec','csc'):
            arg=node.args[0];s,c=S.sin(arg),S.cos(arg)
            return {'tan':s/c,'cot':c/s,'sec':1/c,'csc':1/s}[name]
        return node
    expanded=S.bottom_up(expanded,ratios)
    arguments={node.args[0] for node in S.preorder_traversal(expanded) if node.class_name in ('sin','cos','sinh','cosh')}
    if not arguments:return expr
    if len(arguments)>8:raise ValueError('trigonometric Groebner reduction is limited to eight arguments')
    replacements={};relations=[];pairs=[]
    for arg in sorted(arguments,key=S.default_sort_key):
        for sine,cosine,hyperbolic in [(S.sin,S.cos,False),(S.sinh,S.cosh,True)]:
            a,b=sine(arg),cosine(arg)
            if not expanded.has(a,b):continue
            s,c=S.Dummy('s'),S.Dummy('c');replacements[a]=s;replacements[b]=c;pairs.append((s,c))
            relations.append(c*c-s*s-1 if hyperbolic else s*s+c*c-1)
    native=S.cancel(expanded.xreplace(replacements));numerator,denominator=S.fraction(native);back={new:old for old,new in replacements.items()};candidates=[expr]
    for reverse in ([False] if quick else [False,True]):
        generators=[value for pair in pairs for value in (reversed(pair) if reverse else pair)]
        try:
            basis=S.groebner(relations,*generators,order=order)
            n=basis.reduce(numerator)[1];d=basis.reduce(denominator)[1]
            if d==0:continue
            candidates.append(S.cancel(n/d).xreplace(back))
        except (ValueError,TypeError,NotImplementedError):continue
    return min(candidates,key=lambda value:(S.count_ops(value),len(str(value))))

__all__ += ['trigsimp_old','trigsimp_groebner']
