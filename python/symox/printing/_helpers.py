"""Compatibility entry points for native printers and output streams."""
from .latex import latex,print_latex,multiline_latex
from .repr import srepr
from .python import python,print_python
from .tree import print_tree
from .tableform import TableForm
from .dot import dotprint
from .mathml import mathml,print_mathml
from .preview import preview
from .gtk import print_gtk
from .smtlib import smtlib_code
from .jscode import jscode
from .julia import julia_code
from .octave import octave_code
from .rcode import rcode
from .maple import maple_code
from .mathematica import mathematica_code
from .glsl import glsl_code
from .str import sstr
from .pretty import pretty,pprint,pretty_print
from .pretty.pretty import pprint_use_unicode,pprint_try_use_unicode
from symox._core import cxxcode
from .c import ccode
from .fortran import fcode
from .pycode import pycode
from .rust import rust_code


def sstrrepr(expr,**settings):
    return repr(expr) if isinstance(expr,str) else sstr(expr,**settings)


def pager_print(expr,**settings):
    import pydoc
    pydoc.pager(pretty(expr,**settings))


def print_fcode(expr,**settings):print(fcode(expr,**settings))
def print_glsl(expr,**settings):print(glsl_code(expr,**settings))
def print_jscode(expr,**settings):print(jscode(expr,**settings))
def print_maple_code(expr,**settings):print(maple_code(expr,**settings))
def print_rcode(expr,**settings):print(rcode(expr,**settings))
def print_ccode(expr,**settings):print(ccode(expr,**settings))
