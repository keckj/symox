from .. import Number,NumberSymbol,Integer,Rational,Float,RealNumber,comp,seterr
__all__ = ['Number','NumberSymbol','Integer','Rational','Float','RealNumber','comp','seterr']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'AlgebraicNumber Catalan EulerGamma GoldenRatio TribonacciConstant mod_inverse')

from .. import S as _S, sympify as _sympify

class IntegerConstant(Integer):
    _matches=staticmethod(lambda e:e.is_Integer and int(e) in (-1,0,1))
class RationalConstant(Rational):
    _matches=staticmethod(lambda e:e==_S.Half)


def _constant_class(name,base,value):
    return type(base)(name,(base,),{'__module__':__name__,
        '_matches':staticmethod(lambda e:e==value),'__new__':staticmethod(lambda cls:value)})

Zero=_constant_class('Zero',IntegerConstant,_S.Zero)
One=_constant_class('One',IntegerConstant,_S.One)
NegativeOne=_constant_class('NegativeOne',IntegerConstant,_S.NegativeOne)
Half=_constant_class('Half',RationalConstant,_S.Half)
Infinity=_constant_class('Infinity',Number,_S.Infinity)
NegativeInfinity=_constant_class('NegativeInfinity',Number,_S.NegativeInfinity)
ComplexInfinity=_constant_class('ComplexInfinity',Number,_S.ComplexInfinity)
NaN=_constant_class('NaN',Number,_S.NaN)
ImaginaryUnit=_constant_class('ImaginaryUnit',Number,_S.ImaginaryUnit)
Pi=_constant_class('Pi',NumberSymbol,_S.Pi)
Exp1=_constant_class('Exp1',NumberSymbol,_S.Exp1)


def int_valued(x):
    import math
    if isinstance(x,int):return True
    if isinstance(x,float):return math.isfinite(x) and x.is_integer()
    if not getattr(x,'is_Number',False):return False
    if x.is_Integer:return True
    if x.is_Float:
        try:return _float_fraction(x).denominator==1
        except (ValueError,OverflowError):return False
    return False


def equal_valued(x,y):
    from fractions import Fraction
    x,y=_sympify(x),_sympify(y)
    def exact(value):
        if value.is_Float:
            # Preserve the exact binary value, rather than a decimal tolerance.
            p=getattr(value,'_mpf_',None)
            if p is not None:
                sign,mantissa,exponent,bits=p
                return Fraction((-1)**sign*mantissa)*(Fraction(2)**exponent)
            return _float_fraction(value)
        if value.is_Rational:return Fraction(int(value.p),int(value.q))
        return value
    return exact(x)==exact(y)

def _float_fraction(value):
    from fractions import Fraction
    import mpmath
    from mpmath.libmp import prec_to_dps
    context=mpmath.mp.clone();context.prec=value._prec
    sign,mantissa,exponent,bits=context.mpf(str(value.evalf(prec_to_dps(value._prec)+12)))._mpf_
    return Fraction((-1)**sign*int(mantissa))*(Fraction(2)**exponent)


def all_close(expr1,expr2,rtol=1e-5,atol=1e-8):
    import math
    if isinstance(expr1,(list,tuple)) or isinstance(expr2,(list,tuple)):
        return isinstance(expr1,(list,tuple)) and isinstance(expr2,(list,tuple)) and len(expr1)==len(expr2) and all(all_close(a,b,rtol,atol) for a,b in zip(expr1,expr2))
    if hasattr(expr1,'shape') or hasattr(expr2,'shape'):
        return getattr(expr1,'shape',None)==getattr(expr2,'shape',None) and all(all_close(a,b,rtol,atol) for a,b in zip(expr1,expr2))
    a,b=_sympify(expr1),_sympify(expr2)
    if a==b:return True
    if a.is_number and b.is_number:
        try:return abs(complex((a-b).evalf(30)))<=atol+rtol*abs(complex(b.evalf(30)))
        except (TypeError,ValueError,OverflowError):return False
    return a.func==b.func and len(a.args)==len(b.args) and bool(a.args) and all(all_close(x,y,rtol,atol) for x,y in zip(a.args,b.args))


def mpf_norm(mpf,prec):
    from mpmath.libmp import mpf_pos
    return mpf_pos(mpf,prec)
def sympify_fractions(value):return Rational(value.numerator,value.denominator)
def sympify_complex(value):return _sympify(value.real)+_S.ImaginaryUnit*_sympify(value.imag)
def sympify_mpmath(value):
    from mpmath.libmp import prec_to_dps,to_str
    prec=value.context.prec
    if hasattr(value,'_mpc_'):
        a,b=value._mpc_;return Float(to_str(a,prec_to_dps(prec)+3),precision=prec)+_S.ImaginaryUnit*Float(to_str(b,prec_to_dps(prec)+3),precision=prec)
    return Float(to_str(value._mpf_,prec_to_dps(prec)+3),precision=prec)

__all__ += ['IntegerConstant','RationalConstant','Zero','One','NegativeOne','Half','Infinity','NegativeInfinity',
 'ComplexInfinity','NaN','ImaginaryUnit','Pi','Exp1','int_valued','equal_valued','all_close','mpf_norm',
 'sympify_fractions','sympify_complex','sympify_mpmath']
