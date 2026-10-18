"""Python printer classes compatible with SymPy subclassing."""

from .cxx import CXX17CodePrinter, cxxcode
from .pretty import PrettyPrinter, pretty, pprint, pretty_print
from .str import StrPrinter, sstr
from .pytorch import PyTorchPrinter, TorchPrinter, torchcode
from .jax import JaxPrinter, JAXPrinter, jaxcode
from .slang import SlangCodePrinter, SlangPrinter, slangcode
from .sycl import SYCLCodePrinter, AdaptiveCppCodePrinter, syclcode, adaptivecppcode

__all__ = ["CXX17CodePrinter", "PrettyPrinter", "StrPrinter", "sstr", "cxxcode", "pretty", "pprint", "pretty_print", "PyTorchPrinter", "TorchPrinter", "torchcode", "JaxPrinter", "JAXPrinter", "jaxcode", "SlangCodePrinter", "SlangPrinter", "slangcode", "SYCLCodePrinter", "AdaptiveCppCodePrinter", "syclcode", "adaptivecppcode"]

from ._helpers import (
    TableForm, ccode, dotprint, fcode, glsl_code, jscode, julia_code, latex,
    maple_code, mathematica_code, mathml, multiline_latex, octave_code,
    pager_print, pprint_try_use_unicode, pprint_use_unicode, preview,
    print_ccode, print_fcode, print_glsl, print_gtk, print_jscode, print_latex,
    print_maple_code, print_mathml, print_python, print_rcode, print_tree,
    pycode, python, rcode, rust_code, smtlib_code, srepr, sstrrepr,
)
__all__ += [
    'TableForm','ccode','dotprint','fcode','glsl_code','jscode','julia_code','latex',
    'maple_code','mathematica_code','mathml','multiline_latex','octave_code',
    'pager_print','pprint_try_use_unicode','pprint_use_unicode','preview',
    'print_ccode','print_fcode','print_glsl','print_gtk','print_jscode','print_latex',
    'print_maple_code','print_mathml','print_python','print_rcode','print_tree',
    'pycode','python','rcode','rust_code','smtlib_code','srepr','sstrrepr',
]
