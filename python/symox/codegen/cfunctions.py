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
        if evaluate:
            a=args[0]
            # The ``eval`` of SymPy's isnan and isinf BooleanFunctions.
            if name=='isnan':
                if a is S.nan or a==S.nan:return S.true
                if a.is_number:return S.false
            if name=='isinf':
                if a.is_infinite:return S.true
                if a.is_finite:return S.false
        return S._core._apply(name,*args,evaluate=False)

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
class isnan(S.Function):
    """The C99 ``isnan`` classification (a boolean function in SymPy)."""
    nargs=1
    def _latex(self,printer):
        # SymPy prints boolean functions without the function-argument braces.
        return r'\operatorname{isnan}\left(%s\right)'%printer._print(self.args[0])
    _symox_user_hooks=True
    @classmethod
    def eval(cls,arg):
        if arg is S.nan or arg==S.nan:return S.true
        if arg.is_number:return S.false
class isinf(S.Function):
    """The C99 ``isinf`` classification (a boolean function in SymPy)."""
    nargs=1
    def _latex(self,printer):
        # SymPy prints boolean functions without the function-argument braces.
        return r'\operatorname{isinf}\left(%s\right)'%printer._print(self.args[0])
    _symox_user_hooks=True
    @classmethod
    def eval(cls,arg):
        if arg.is_infinite:return S.true
        if arg.is_finite:return S.false
_Two=S.Integer(2);_Ten=S.Integer(10)


def _argindex_error(self,argindex):
    from ..core.function import ArgumentIndexError
    raise ArgumentIndexError(self,argindex)


def _log_base_atom(arg,base):
    # SymPy's ``log.eval(arg, base)`` result, kept only when it is an atom.
    result=S.log(arg,base)
    return result if result.is_Atom else None


class exp2(S.Function):
    """``2**x``, as the C99 ``exp2``."""
    nargs=1
    _symox_user_hooks=True
    def fdiff(self,argindex=1):
        if argindex==1:return self*S.log(_Two)
        _argindex_error(self,argindex)
    def _eval_rewrite_as_Pow(self,arg,**kwargs):return S.Pow(_Two,arg)
    _eval_rewrite_as_tractable=_eval_rewrite_as_Pow
    def _eval_expand_func(self,**hints):return S.Pow(_Two,self.args[0])
    @classmethod
    def eval(cls,arg):
        if arg.is_number:return S.Pow(_Two,arg)


class log2(S.Function):
    """``log(x)/log(2)``, as the C99 ``log2``."""
    nargs=1
    _symox_user_hooks=True
    def fdiff(self,argindex=1):
        if argindex==1:return S.Integer(1)/(S.log(_Two)*self.args[0])
        _argindex_error(self,argindex)
    @classmethod
    def eval(cls,arg):
        if arg.is_number:return _log_base_atom(arg,_Two)
        elif arg.is_Pow and arg.base==_Two:return arg.exp
    def _eval_evalf(self,prec):
        from ..core.evalf import prec_to_dps
        return (S.log(self.args[0])/S.log(_Two)).evalf(prec_to_dps(prec))
    def _eval_expand_func(self,**hints):return S.log(self.args[0])/S.log(_Two)
    def _eval_rewrite_as_log(self,arg,**kwargs):return S.log(arg)/S.log(_Two)
    _eval_rewrite_as_tractable=_eval_rewrite_as_log


class fma(S.Function):
    """``x*y + z``, as the C99 fused multiply-add."""
    nargs=3
    _symox_user_hooks=True
    def fdiff(self,argindex=1):
        if argindex in (1,2):return self.args[2-argindex]
        elif argindex==3:return S.Integer(1)
        _argindex_error(self,argindex)
    def _eval_expand_func(self,**hints):
        x,y,z=self.args;return x*y+z
    def _eval_rewrite_as_tractable(self,arg,limitvar=None,**kwargs):
        # The reference passes a single argument to the three-argument
        # helper, so this rewrite raises TypeError there as well.
        raise TypeError('_fma() missing 2 required positional arguments')


class log10(S.Function):
    """``log(x)/log(10)``, as the C99 ``log10``."""
    nargs=1
    _symox_user_hooks=True
    _imp_=_numeric('log10')
    def fdiff(self,argindex=1):
        if argindex==1:return S.Integer(1)/(S.log(_Ten)*self.args[0])
        _argindex_error(self,argindex)
    @classmethod
    def eval(cls,arg):
        if arg.is_number:return _log_base_atom(arg,_Ten)
        elif arg.is_Pow and arg.base==_Ten:return arg.exp
    def _eval_evalf(self,prec):
        import mpmath
        return _evaluated(mpmath.log10,self.args[0],prec)
    def _eval_expand_func(self,**hints):return S.log(self.args[0])/S.log(_Ten)
    def _eval_rewrite_as_log(self,arg,**kwargs):return S.log(arg)/S.log(_Ten)
    _eval_rewrite_as_tractable=_eval_rewrite_as_log


class Sqrt(S.Function):
    """The square root, as the C99 ``sqrt``; it is never evaluated."""
    nargs=1
    _symox_user_hooks=True
    def fdiff(self,argindex=1):
        if argindex==1:return S.Pow(self.args[0],S.Rational(-1,2))/_Two
        _argindex_error(self,argindex)
    def _eval_expand_func(self,**hints):return S.Pow(self.args[0],S.Rational(1,2))
    def _eval_rewrite_as_Pow(self,arg,**kwargs):return S.Pow(arg,S.Rational(1,2))
    _eval_rewrite_as_tractable=_eval_rewrite_as_Pow
    def _eval_evalf(self,prec):return self


class Cbrt(S.Function):
    """The cube root, as the C99 ``cbrt``; it is never evaluated."""
    nargs=1
    _symox_user_hooks=True
    def fdiff(self,argindex=1):
        if argindex==1:return S.Pow(self.args[0],S.Rational(-2,3))/3
        _argindex_error(self,argindex)
    def _eval_expand_func(self,**hints):return S.Pow(self.args[0],S.Rational(1,3))
    def _eval_rewrite_as_Pow(self,arg,**kwargs):return S.Pow(arg,S.Rational(1,3))
    _eval_rewrite_as_tractable=_eval_rewrite_as_Pow
    def _eval_evalf(self,prec):return self


class hypot(S.Function):
    """``sqrt(x**2 + y**2)``, as the C99 ``hypot``."""
    nargs=2
    _symox_user_hooks=True
    def fdiff(self,argindex=1):
        if argindex in (1,2):return 2*self.args[argindex-1]/(_Two*self.func(*self.args))
        _argindex_error(self,argindex)
    def _eval_expand_func(self,**hints):
        x,y=self.args;return S.sqrt(S.Pow(x,2)+S.Pow(y,2))
    def _eval_rewrite_as_Pow(self,arg,**kwargs):
        # The reference helper takes two arguments, so the rewrite raises.
        raise TypeError('hypot._eval_rewrite_as_Pow() takes 2 positional arguments but 3 were given')
    _eval_rewrite_as_tractable=_eval_rewrite_as_Pow
    def _eval_evalf(self,prec):
        import mpmath
        x,y=self.args
        with mpmath.workprec(prec+20):
            a,b=_mpmath_number(x,prec+20),_mpmath_number(y,prec+20)
            if a is None or b is None:return None
            return _float_result(mpmath.hypot(a,b),prec)


__all__=['Sqrt','Cbrt','exp2','expm1','log2','log10','log1p','isnan','isinf','fma','hypot']
