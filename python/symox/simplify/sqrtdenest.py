"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'sqrtdenest')

import symox as S

class SqrtdenestStopIteration(Exception):pass

def is_sqrt(expr):return getattr(expr,'is_Pow',False) and expr.exp==S.Rational(1,2)
def is_algebraic(p):return S.sympify(p).is_algebraic is True

def sqrt_depth(p):
    p=S.sympify(p)
    return int(is_sqrt(p))+max((sqrt_depth(arg) for arg in p.args),default=0)
def sqrt_biquadratic_denest(expr,a,b,r,d2):
    expr=S.sympify(expr);candidate=S.sqrtdenest(expr)
    if sqrt_depth(candidate)<sqrt_depth(expr) and S.simplify(candidate**2-expr**2)==0:return candidate
    return None

__all__ += ['SqrtdenestStopIteration','is_algebraic','is_sqrt','sqrt_biquadratic_denest','sqrt_depth']
