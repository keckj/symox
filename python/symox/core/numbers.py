from .. import Number,NumberSymbol,Integer,Rational,Float,RealNumber,comp,seterr
__all__ = ['Number','NumberSymbol','Integer','Rational','Float','RealNumber','comp','seterr']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'AlgebraicNumber mod_inverse oo zoo nan E pi I')

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
GoldenRatio=_constant_class('GoldenRatio',NumberSymbol,_S.GoldenRatio)
EulerGamma=_constant_class('EulerGamma',NumberSymbol,_S.EulerGamma)
Catalan=_constant_class('Catalan',NumberSymbol,_S.Catalan)
TribonacciConstant=_constant_class('TribonacciConstant',NumberSymbol,_S.TribonacciConstant)


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


def all_close(expr1, expr2, rtol=1e-5, atol=1e-8):
    # sympy/core/numbers.py all_close (SymPy 1.14): Add and Mul compare their
    # coefficients with the tolerances and match the remaining factors.
    from symox import Rational, Float
    from symox.core.sorting import ordered
    NUM_TYPES = (Rational, Float)

    def _all_close(obj1, obj2):
        if type(obj1) == type(obj2) and isinstance(obj1, (list, tuple)):
            if len(obj1) != len(obj2):
                return False
            return all(_all_close(e1, e2) for e1, e2 in zip(obj1, obj2))
        else:
            return _all_close_expr(_sympify(obj1), _sympify(obj2))

    def _all_close_expr(expr1, expr2):
        num1 = isinstance(expr1, NUM_TYPES)
        num2 = isinstance(expr2, NUM_TYPES)
        if num1 != num2:
            return False
        elif num1:
            return _close_num(expr1, expr2)
        if expr1.is_Add or expr1.is_Mul or expr2.is_Add or expr2.is_Mul:
            return _all_close_ac(expr1, expr2)
        if expr1.func != expr2.func or len(expr1.args) != len(expr2.args):
            return False
        args = zip(expr1.args, expr2.args)
        return all(_all_close_expr(a1, a2) for a1, a2 in args)

    def _close_num(num1, num2):
        return bool(abs(num1 - num2) <= atol + rtol*abs(num2))

    def _all_close_ac(expr1, expr2):
        # compare expressions with associative commutative operators for
        # approximate equality by seeing that all terms have equivalent
        # coefficients (which are always Rational or Float)
        if expr1.is_Mul or expr2.is_Mul:
            # as_coeff_mul automatically will supply coeff of 1
            c1, e1 = expr1.as_coeff_mul(rational=False)
            c2, e2 = expr2.as_coeff_mul(rational=False)
            if not _close_num(c1, c2):
                return False
            s1 = set(e1)
            s2 = set(e2)
            common = s1 & s2
            s1 -= common
            s2 -= common
            if not s1:
                return True
            if not any(i.has(Float) for j in (s1, s2) for i in j):
                return False
            # factors might not be matching, e.g.
            # x != x**1.0, exp(x) != exp(1.0*x), etc...
            s1 = [i.as_base_exp() for i in ordered(s1)]
            s2 = [i.as_base_exp() for i in ordered(s2)]
            unmatched = list(range(len(s1)))
            for be1 in s1:
                for i in unmatched:
                    be2 = s2[i]
                    if _all_close(be1, be2):
                        unmatched.remove(i)
                        break
                else:
                    return False
            return not(unmatched)
        assert expr1.is_Add or expr2.is_Add
        cd1 = expr1.as_coefficients_dict()
        cd2 = expr2.as_coefficients_dict()
        # this test will assure that the key of 1 is in
        # each dict and that they have equal values
        if not _close_num(cd1[1], cd2[1]):
            return False
        if len(cd1) != len(cd2):
            return False
        for k in list(cd1):
            if k in cd2:
                if not _close_num(cd1.pop(k), cd2.pop(k)):
                    return False
            # k (or a close version in cd2) might have
            # Floats in a factor of the term which will
            # be handled below
        else:
            if not cd1:
                return True
        for k1 in cd1:
            for k2 in cd2:
                if _all_close_expr(k1, k2):
                    # found a matching key
                    # XXX there could be a corner case where
                    # more than 1 might match and the numbers are
                    # such that one is better than the other
                    # that is not being considered here
                    if not _close_num(cd1[k1], cd2[k2]):
                        return False
                    break
            else:
                # no key matched
                return False
        return True

    return _all_close(expr1, expr2)


def mpf_norm(mpf,prec):
    from mpmath.libmp import mpf_pos
    return mpf_pos(mpf,prec)
def sympify_fractions(f):return Rational(f.numerator,f.denominator)
def sympify_complex(a):return _sympify(a.real)+_S.ImaginaryUnit*_sympify(a.imag)
def sympify_mpmath(x):
    from mpmath.libmp import prec_to_dps,to_str
    prec=x.context.prec
    if hasattr(x,'_mpc_'):
        a,b=x._mpc_;return Float(to_str(a,prec_to_dps(prec)+3),precision=prec)+_S.ImaginaryUnit*Float(to_str(b,prec_to_dps(prec)+3),precision=prec)
    return Float(to_str(x._mpf_,prec_to_dps(prec)+3),precision=prec)

__all__ += ['IntegerConstant','RationalConstant','Zero','One','NegativeOne','Half','Infinity','NegativeInfinity',
 'ComplexInfinity','NaN','ImaginaryUnit','Pi','Exp1','GoldenRatio','EulerGamma','Catalan','TribonacciConstant','int_valued','equal_valued','all_close','mpf_norm',
 'sympify_fractions','sympify_complex','sympify_mpmath']
