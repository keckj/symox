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
class expm1(_Intrinsic):pass
class log2(_Intrinsic):pass
class log10(_Intrinsic):pass
class log1p(_Intrinsic):pass
class isnan(_Intrinsic):pass
class isinf(_Intrinsic):pass
class fma(_Intrinsic):_arity=(3,)
class hypot(_Intrinsic):_arity=(2,)
__all__=['Sqrt','Cbrt','exp2','expm1','log2','log10','log1p','isnan','isinf','fma','hypot']
