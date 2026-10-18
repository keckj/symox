"""Render Fortran procedures in a module."""
from ..printing.fortran import FCodePrinter
def render_as_module(content,name,Printer=FCodePrinter,settings=None):
    from .fnodes import Module
    nodes=list(content) if isinstance(content,(list,tuple)) else [content]
    return Printer(settings or {}).doprint(Module(name,(),nodes))+'\n'
__all__=['render_as_module']
