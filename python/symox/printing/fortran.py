"""Native fortran code generation entry point."""
from symox._core import fcode
__all__ = ['fcode']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_fcode')

from ._legacy import LegacyCodePrinter as CodePrinter
class FCodePrinter(CodePrinter):
    language='Fortran';printmethod='_fcode';_style='fortran'
    def _print_Pow(self,expr):return '('+self._print(expr.base)+'**'+self._print(expr.exp)+')'
    def _print_Rational(self,expr):return str(expr.p)+'.0d0/'+str(expr.q)+'.0d0'
    def _print_Assignment(self,expr):return self._print(expr.lhs)+' = '+self._print(expr.rhs)
    def _print_CodeBlock(self,expr):return '\n'.join(self._print(v) for v in expr.args)
    def _print_Comment(self,expr):return '! '+expr.text
    def emptyPrinter(self,expr):
        from symox import _core
        if getattr(expr,'is_Function',False):return {'Abs':'abs','ceiling':'ceiling','Min':'min','Max':'max'}.get(expr.class_name,expr.class_name)+'('+', '.join(self._print(a) for a in expr.args)+')'
        return _core.fcode(expr).strip()
__all__ += ['FCodePrinter']
FCodePrinter._print_And=lambda self,e:'('+' .and. '.join(self._print(a) for a in e.args)+')'
FCodePrinter._print_Or=lambda self,e:'('+' .or. '.join(self._print(a) for a in e.args)+')'
FCodePrinter._print_Not=lambda self,e:'.not. ('+self._print(e.args[0])+')'
FCodePrinter._print_Xor=lambda self,e:'('+' .neqv. '.join(self._print(a) for a in e.args)+')'
FCodePrinter._print_Min=lambda self,e:'min('+', '.join(self._print(a) for a in e.args)+')'
FCodePrinter._print_Max=lambda self,e:'max('+', '.join(self._print(a) for a in e.args)+')'
