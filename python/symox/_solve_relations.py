"""Route relational solve inputs through the native inequality engine."""
from functools import wraps
import sys


def install(namespace):
    native=namespace['solve']

    @wraps(native)
    def solve(f, *symbols, **flags):
        from . import Tuple, reduce_inequalities, sympify
        sequences=(list,tuple,set,frozenset,Tuple)
        if hasattr(f,'__next__'):
            f=list(f)
        equations=f if isinstance(f,sequences) else (f,)
        prepared=[];relational=False
        for equation in equations:
            value=sympify(equation) if isinstance(equation,str) else equation
            prepared.append(value)
            if getattr(value,'is_Relational',False) and value.rel_op!='==':
                relational=True
        if relational:
            targets=symbols[0] if len(symbols)==1 else symbols
            return reduce_inequalities(prepared,targets)
        return native(f,*symbols,**flags)

    namespace['solve']=solve
    # Keep public import paths coherent without replacing the private binding.
    for name,module in tuple(sys.modules.items()):
        if (name.startswith('symox.') and name!='symox._core' and module is not None
                and module.__dict__.get('solve') is native):
            module.solve=solve
