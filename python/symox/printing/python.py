"""Reconstruct native symbolic expressions as executable Python programs."""
from .repr import srepr

def python(expr,**settings):
    return 'from symox import *\ne = '+srepr(expr,**settings)

def print_python(expr,**settings):
    print(python(expr,**settings))

from .str import StrPrinter

class PythonPrinter(StrPrinter):
    printmethod='_pythoncode'
    def __init__(self,settings=None):super().__init__(settings);self.symbols=[];self.functions=[]
    def _print_Symbol(self,expr):
        if expr not in self.symbols:self.symbols.append(expr)
        return expr.name
    def emptyPrinter(self,expr):
        self.symbols.extend(s for s in getattr(expr,'free_symbols',()) if s not in self.symbols)
        return super().emptyPrinter(expr)
