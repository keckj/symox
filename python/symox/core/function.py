from .. import Function,FunctionClass,WildFunction,Derivative,Lambda,Subs,arity,expand_multinomial
__all__ = ['Function','FunctionClass','WildFunction','Derivative','Lambda','Subs','arity','expand_multinomial']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'BadArgumentsError BadSignatureError PoleError count_ops diff expand expand_complex expand_func expand_log expand_mul expand_power_base expand_power_exp expand_trig nfloat')

from .._compat import _ExprType

class ArgumentIndexError(ValueError):
    def __init__(self,function,index):self.function,self.index=function,index;super().__init__('Invalid operation with argument number %s for Function %s'%(index,function))

from .._function_heads import Application, DefinedFunction  # noqa: E402
class AppliedUndef(Application):_matches=staticmethod(lambda e:e.is_UndefinedFunction)

class _UndefinedFunctionType(type):
    def __instancecheck__(cls,value):
        from .. import _is_function_head
        return type.__instancecheck__(cls,value) or _is_function_head(value)
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
