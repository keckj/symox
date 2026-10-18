"""Integration adapters using the native engine's resolved exact results.

These functions preserve the reference result/condition conventions. They do
not implement a second independent Meijer-G formula catalog.
"""
import symox as S


def _resolved(expr):return not expr.has(S.Integral)

def meijerint_indefinite(f,x):
    result=S.integrate(f,x)
    return result if _resolved(result) else None

def meijerint_definite(f,x,a,b):
    integral=S.Integral(f,(x,a,b))
    try:converges=integral.is_convergent()
    except (NotImplementedError,ValueError,TypeError,AttributeError):converges=None
    if converges is False or converges==S.false:return None
    expression=S.cancel(f)
    if expression.is_rational_function(x):
        denominator=S.denom(expression)
        zeros=S.solveset(denominator,x,domain=S.Interval(a,b))
        if getattr(zeros,'is_finite_set',False) and zeros!=S.EmptySet:return None
    result=S.integrate(f,(x,a,b))
    if not _resolved(result):return None
    if getattr(result,'class_name','')=='Piecewise':
        from ..printing.mathml import _pairs
        for value,condition in _pairs(result):
            if _resolved(value) and not value.has(S.nan,S.oo,-S.oo):return value,condition
        return None
    if result.has(S.nan,S.oo,-S.oo):return None
    if (S.sympify(f).free_symbols|S.sympify(a).free_symbols|S.sympify(b).free_symbols)-{x} and converges not in (True,S.true):return None
    return result,S.true

def meijerint_inversion(f,x,t):
    result=S.inverse_laplace_transform(f,x,t)
    return None if result.has(S.InverseLaplaceTransform) else result

__all__=['meijerint_indefinite','meijerint_definite','meijerint_inversion']
