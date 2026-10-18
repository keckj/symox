"""Native c code generation entry point."""
from symox._core import ccode
__all__ = ['ccode']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.printing', 'print_ccode')

from .codeprinter import CodePrinter,requires
class C89CodePrinter(CodePrinter):standard='C89'
class C99CodePrinter(CodePrinter):
    standard='C99'
    def __init__(self,settings=None,**kwargs):
        super().__init__(settings,**kwargs);self.known_functions={**known_functions_C99,**self.known_functions}
class C11CodePrinter(C99CodePrinter):standard='C11'
known_functions_C89={n:n for n in 'sin cos tan asin acos atan atan2 sinh cosh tanh exp log log10 sqrt pow floor ceil fabs'.split()}
known_functions_C99={**known_functions_C89,**{n:n for n in 'cbrt exp2 expm1 log2 log1p erf erfc tgamma lgamma fma hypot'.split()}}
def get_math_macros():return {'E':'M_E','Pi':'M_PI','log2(E)':'M_LOG2E','log(2)':'M_LN2','sqrt(2)':'M_SQRT2'}
__all__ += ['C89CodePrinter','C99CodePrinter','C11CodePrinter','CodePrinter','known_functions_C89','known_functions_C99','get_math_macros','requires']

for _name,_target in {'Sqrt':'sqrt','Cbrt':'cbrt','gamma':'tgamma','loggamma':'lgamma'}.items():known_functions_C99[_name]=_target
