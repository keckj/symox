"""Render executable Python modules from statement nodes."""
from ..printing.pycode import PythonCodePrinter
def render_as_module(content,Printer=PythonCodePrinter,settings=None):
    from .ast import CodeBlock
    if isinstance(content,(list,tuple)):content=CodeBlock(*content)
    return 'import math\n\n'+Printer(settings or {}).doprint(content)+'\n'
__all__=['render_as_module']
