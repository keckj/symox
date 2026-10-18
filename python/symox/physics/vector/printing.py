"""Text and LaTeX printing for physics vector expressions."""
import symox as S
from ._core import Vector,Dyadic


def vsprint(expr,**settings):return str(expr)
def vsstrrepr(expr,**settings):return str(expr)
def vprint(expr,**settings):print(vsprint(expr,**settings))
def vpprint(expr,**settings):print(vsprint(expr,**settings))


def vlatex(expr,**settings):
    if isinstance(expr,Vector):
        terms = []
        for frame,components in expr._components.items():
            for i,value in enumerate(components):
                if value!=0:terms.append(S.latex(value,**settings)+r'\mathbf{\hat{'+frame.name.lower()+'}_{'+frame._indices[i]+'}}')
        return ' + '.join(terms) or '0'
    if isinstance(expr,Dyadic):return ' + '.join(S.latex(c,**settings)+vlatex(a)+r'\otimes '+vlatex(b) for c,a,b in expr._terms()) or '0'
    return S.latex(expr,**settings)


def init_vprinting(**kwargs):
    try:from IPython import get_ipython
    except ImportError:return None
    shell = get_ipython()
    if shell is not None:shell.display_formatter.formatters['text/latex'].for_type(Vector,lambda value:'$'+vlatex(value)+'$')
