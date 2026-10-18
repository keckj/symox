"""Native rust code generation entry point."""
from symox._core import rust_code
__all__ = ['rust_code']

from .codeprinter import CodePrinter
class RustCodePrinter(CodePrinter):
    language='Rust';printmethod='_rust_code';_style='rust'
    def _print_Integer(self,expr):return str(int(expr))+'.0'
    def _print_Pow(self,expr):
        b,e=expr.args
        return '('+self._print(b)+').'+('powi' if e.is_Integer else 'powf')+'('+(str(int(e)) if e.is_Integer else self._print(e))+')'
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Function',False):
            name={'Abs':'abs','ceiling':'ceil','log':'ln'}.get(expr.class_name,expr.class_name)
            if len(expr.args)==1 and name in known_functions:
                return '('+self._print(expr.args[0])+').'+name+'()'
        from symox import _core
        return _core.rust_code(expr)
class TypeCast:
    def __init__(self,value,type):self.value,self.type=value,type
    @property
    def args(self):return self.value,self.type
    def _rust_code(self,printer):return '('+printer._print(self.value)+' as '+str(self.type)+')'
def float_floor(value):
    from symox import floor
    return floor(value)
def float_ceiling(value):
    from symox import ceiling
    return ceiling(value)
known_functions={name:name for name in 'sin cos tan asin acos atan sinh cosh tanh exp ln log10 log2 sqrt cbrt abs floor ceil'.split()}
__all__ += ['RustCodePrinter','TypeCast','float_floor','float_ceiling','known_functions']
