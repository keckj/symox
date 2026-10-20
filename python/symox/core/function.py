from .. import Function,FunctionClass,WildFunction,Derivative,Lambda,Subs,arity,expand_multinomial
__all__ = ['Function','FunctionClass','WildFunction','Derivative','Lambda','Subs','arity','expand_multinomial']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'BadArgumentsError BadSignatureError PoleError count_ops diff expand expand_complex expand_func expand_log expand_mul expand_power_base expand_power_exp expand_trig nfloat')

from .._compat import _ExprType

class ArgumentIndexError(ValueError):
    def __init__(self,function,index):self.function,self.index=function,index;super().__init__('invalid argument index %s for %s'%(index,function))

class Application(metaclass=_ExprType):
    _matches=staticmethod(lambda e:e.is_Function)
    def __new__(cls,*args):
        if cls is Application:raise TypeError('instantiate a concrete function application')
        return Function(cls.__name__)(*args)
class DefinedFunction(Application):
    _matches=staticmethod(lambda e:e.is_Function and not e.is_UndefinedFunction)
    def __new__(cls,*args):
        from symox import sympify
        values=tuple(sympify(a) for a in args)
        evaluator=getattr(cls,'eval',None)
        if evaluator:
            result=evaluator(*values)
            if result is not None:return sympify(result)
        return Function(cls.__name__)(*values)
class AppliedUndef(Application):_matches=staticmethod(lambda e:e.is_UndefinedFunction)

class _UndefinedFunctionType(type):
    def __instancecheck__(cls,value):return type.__instancecheck__(cls,value) or isinstance(value,Function)
class UndefinedFunction(type,metaclass=_UndefinedFunctionType):
    def __new__(mcls,name,bases=(),namespace=None,**kwargs):
        return Function(name,**kwargs)

class UndefSageHelper:
    """Helper to facilitate Sage conversion of undefined functions."""
    def __get__(self,ins,typ):
        import sage.all as sage
        if ins is None:
            return lambda:sage.function(typ.__name__)
        else:
            args=[arg._sage_() for arg in ins.args]
            return lambda:sage.function(ins.__class__.__name__)(*args)

# Applied undefined-function classes read this at creation, as in SymPy.
_undef_sage_helper=UndefSageHelper()

__all__ += ['Application','DefinedFunction','AppliedUndef','UndefinedFunction','ArgumentIndexError','UndefSageHelper']
