"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'ccode cxxcode fcode print_ccode print_fcode rust_code')

from functools import wraps as _wraps
from .cxx import CXX17CodePrinter
from ._nodes import NodePrinter

class AssignmentError(ValueError):pass
class PrintMethodNotImplementedError(NotImplementedError):pass

def requires(**requirements):
    def decorate(function):
        @_wraps(function)
        def wrapped(self,*args,**kwargs):
            for attribute,values in requirements.items():
                current=getattr(self,attribute,None)
                if current is None:current=set();setattr(self,attribute,current)
                current.update(values)
            return function(self,*args,**kwargs)
        return wrapped
    return decorate

class CodePrinter(NodePrinter,CXX17CodePrinter):
    language='C';printmethod='_ccode'
    def __init__(self,settings=None,**kwargs):
        super().__init__(settings,**kwargs);self.headers=set();self.libraries=set();self.macros=set()
    def _print_Pow(self,expr):
        from symox import Rational
        base,exponent=expr.args
        if exponent==Rational(1,2):return 'sqrt('+self._print(base)+')'
        if exponent==Rational(1,3):return 'cbrt('+self._print(base)+')'
        return 'pow('+self._print(base)+', '+self._print(exponent)+')'
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Function',False):
            name=expr.class_name;function=self.known_functions.get(name)
            if function is None:
                if getattr(expr,'is_UndefinedFunction',False):raise PrintMethodNotImplementedError('unknown function '+name+' requires user_functions')
                function={'Abs':'fabs','ceiling':'ceil','Max':'fmax','Min':'fmin'}.get(name,name)
            if not isinstance(function,str):function=next((value for condition,value in function if condition(*expr.args)),name)
            return function+'('+', '.join(self._print(a) for a in expr.args)+')'
        from symox import _core
        return _core.ccode(expr)

def print_rust_code(expr,**settings):
    from .rust import rust_code
    print(rust_code(expr,**settings))
__all__ += ['CodePrinter','AssignmentError','PrintMethodNotImplementedError','requires','print_rust_code']


def _c_relation(self,expr,op):return '('+self._print(expr.lhs)+' '+op+' '+self._print(expr.rhs)+')'
CodePrinter._print_And=lambda self,e:'('+' && '.join(self._print(a) for a in e.args)+')'
CodePrinter._print_Or=lambda self,e:'('+' || '.join(self._print(a) for a in e.args)+')'
CodePrinter._print_Not=lambda self,e:'!('+self._print(e.args[0])+')'
CodePrinter._print_Xor=lambda self,e:'('+' ^ '.join(self._print(a) for a in e.args)+')'
CodePrinter._print_ITE=lambda self,e:'('+self._print(e.args[0])+' ? '+self._print(e.args[1])+' : '+self._print(e.args[2])+')'
for _name,_op in [('Equality','=='),('Unequality','!='),('LessThan','<='),('StrictLessThan','<'),('GreaterThan','>='),('StrictGreaterThan','>')]:
    setattr(CodePrinter,'_print_'+_name,lambda self,e,op=_op:_c_relation(self,e,op))
def _c_minmax(self,expr):
    function='fmin' if expr.class_name=='Min' else 'fmax';value=self._print(expr.args[0])
    for arg in expr.args[1:]:value=function+'('+value+', '+self._print(arg)+')'
    return value
CodePrinter._print_Min=CodePrinter._print_Max=_c_minmax
