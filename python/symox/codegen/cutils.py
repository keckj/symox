"""Render complete C source from program nodes."""
from ..printing.c import C99CodePrinter
def render_as_source_file(content,Printer=C99CodePrinter,settings=None):
    from .ast import CodeBlock
    if isinstance(content,(list,tuple)):content=CodeBlock(*content)
    printer=Printer(settings or {})
    return '#include <math.h>\n#include <stdbool.h>\n#include <stdint.h>\n#include <stdio.h>\n#include <stdlib.h>\n\n'+printer.doprint(content)+'\n'
__all__=['render_as_source_file']
