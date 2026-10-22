"""Symbolic C99 intrinsic calls with exact evaluation of constant arguments."""
import symox as S
from .._compat import _ExprType


class _Intrinsic(metaclass=_ExprType):
    _arity=(1,)
    @classmethod
    def _matches(cls,value):return isinstance(value,S.Expr) and value.class_name==cls.__name__
    def __new__(cls,*args,evaluate=True):
        if len(args) not in cls._arity:raise TypeError(cls.__name__+' received the wrong number of arguments')
        args=tuple(map(S.sympify,args));name=cls.__name__
        if evaluate and all(a.is_number for a in args):
            a=args[0]
            if name=='isnan':return S.true if a==S.nan else S.false
            if name=='isinf':return S.true if a in (S.oo,-S.oo) else S.false
            formulas={'Sqrt':lambda:S.sqrt(a),'Cbrt':lambda:S.real_root(a,3),
                'exp2':lambda:2**a,'expm1':lambda:S.exp(a)-1,'log2':lambda:S.log(a,2),
                'log10':lambda:S.log(a,10),'log1p':lambda:S.log(1+a),
                'fma':lambda:args[0]*args[1]+args[2],
                'hypot':lambda:S.sqrt(sum(v*v for v in args))}
            return formulas[name]()
        return S._core._apply(name,*args,evaluate=False)

class Sqrt(_Intrinsic):pass
class Cbrt(_Intrinsic):pass
class exp2(_Intrinsic):pass
def _has_float(value):
    return bool(value.atoms(S.Float))


def _mpmath_number(value,prec):
    # A numeric argument at ``prec`` bits, with the guard digits of evalf.
    import mpmath
    from mpmath.libmp import prec_to_dps
    dps=prec_to_dps(prec)+5
    real,imag=S.N(value,dps).as_real_imag()
    if not (real.is_number and imag.is_number):return None
    try:
        real,imag=mpmath.mpf(str(S.N(real,dps))),mpmath.mpf(str(S.N(imag,dps)))
    except (TypeError,ValueError):return None
    return real if imag==0 else mpmath.mpc(real,imag)


def _float_result(value,prec):
    import mpmath
    if isinstance(value,mpmath.mpc):
        return S.Float(value.real,precision=prec)+S.I*S.Float(value.imag,precision=prec)
    return S.Float(value,precision=prec)


def _evaluated(function,arg,prec):
    import mpmath
    with mpmath.workprec(prec+20):
        value=_mpmath_number(arg,prec+20)
        if value is None:return None
        try:value=function(value)
        except (ValueError,ZeroDivisionError):return None
        return _float_result(value,prec)


def _numeric(name):
    # Numerical implementation used by lambdify: SymPy's NumPy and math
    # printers emit numpy.expm1/numpy.log1p or math.expm1/math.log1p.
    def implementation(value):
        try:
            import numpy
        except ImportError:
            import math
            return getattr(math,name)(value)
        return getattr(numpy,name)(value)
    implementation.__name__=name
    return staticmethod(implementation)


class expm1(S.Function):
    """``exp(x) - 1``, evaluated without cancellation near ``x = 0``."""
    nargs=1
    _symox_user_hooks=True
    _imp_=_numeric('expm1')
    def fdiff(self,argindex=1):
        if argindex==1:return S.exp(*self.args)
        from ..core.function import ArgumentIndexError;raise ArgumentIndexError(self,argindex)
    def _eval_expand_func(self,**hints):return S.exp(self.args[0])-1
    def _eval_rewrite_as_exp(self,arg,**kwargs):return S.exp(arg)-1
    _eval_rewrite_as_tractable=_eval_rewrite_as_exp
    @classmethod
    def eval(cls,arg):
        if arg.is_number and _has_float(arg):
            # SymPy's Function.__new__ evaluates Float arguments at their
            # precision with mpmath.expm1.
            import mpmath
            prec=max(f._prec for f in arg.atoms(S.Float))
            return _evaluated(mpmath.expm1,arg,prec)
        value=S.exp(arg)
        if not (getattr(value,'class_name','')=='exp' and len(value.args)==1):return value-1
    def _eval_evalf(self,prec):
        import mpmath
        return _evaluated(mpmath.expm1,self.args[0],prec)
    def _eval_is_real(self):return self.args[0].is_real
    def _eval_is_finite(self):return self.args[0].is_finite
class log2(_Intrinsic):pass
class log10(_Intrinsic):pass
class log1p(S.Function):
    """``log(x + 1)``, evaluated without cancellation near ``x = 0``."""
    nargs=1
    _symox_user_hooks=True
    _imp_=_numeric('log1p')
    def fdiff(self,argindex=1):
        if argindex==1:return 1/(self.args[0]+1)
        from ..core.function import ArgumentIndexError;raise ArgumentIndexError(self,argindex)
    def _eval_expand_func(self,**hints):return S.log(self.args[0]+1)
    def _eval_rewrite_as_log(self,arg,**kwargs):return S.log(arg+1)
    _eval_rewrite_as_tractable=_eval_rewrite_as_log
    @classmethod
    def eval(cls,arg):
        if arg.is_Rational:return S.log(arg+1)
        if arg.is_Float:
            # Adding 1 to a Float loses its digits: SymPy converts it exactly.
            return S.log(S.Rational(arg)+1)
        if arg.is_number and _has_float(arg):
            import mpmath
            prec=max(f._prec for f in arg.atoms(S.Float))
            return _evaluated(mpmath.log1p,arg,prec)
        value=S.log(arg+1)
        if not (getattr(value,'class_name','')=='log' and len(value.args)==1 and value.args[0]==arg+1):return value
    def _eval_evalf(self,prec):
        import mpmath
        return _evaluated(mpmath.log1p,self.args[0],prec)
    def _eval_is_real(self):return (self.args[0]+1).is_nonnegative
    def _eval_is_finite(self):
        if (self.args[0]+1).is_zero:return False
        return self.args[0].is_finite
    def _eval_is_positive(self):return self.args[0].is_positive
    def _eval_is_zero(self):return self.args[0].is_zero
    def _eval_is_nonnegative(self):return self.args[0].is_nonnegative
class isnan(_Intrinsic):pass
class isinf(_Intrinsic):pass
class fma(_Intrinsic):_arity=(3,)
class hypot(_Intrinsic):_arity=(2,)
__all__=['Sqrt','Cbrt','exp2','expm1','log2','log10','log1p','isnan','isinf','fma','hypot']
