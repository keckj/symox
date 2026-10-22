"""Held and exact Clebsch-Gordan and Wigner coupling expressions."""
import symox as S
from .. import wigner as W


class _Coupling:
    nargs=6
    def __init__(self,*args):
        if len(args)!=self.nargs:raise TypeError('incorrect number of coupling indices')
        self.args=tuple(map(S.sympify,args));self._expr=S._core._apply(type(self).__name__,*self.args,evaluate=False)
    def _symox_(self):return self._expr
    def __repr__(self):return str(self._expr)
    @property
    def free_symbols(self):return self._expr.free_symbols
    def subs(self,*args,**kwargs):return type(self)(*(a.subs(*args,**kwargs) for a in self.args))
    def __mul__(self,other):return self._expr*other
    __rmul__=__mul__
    def __add__(self,other):return self._expr+other
    __radd__=__add__
    def __pow__(self,power):return self._expr**power
class CG(_Coupling):
    def doit(self,**hints):
        if self.free_symbols:return self
        a,p,b,q,c,r=self.args;return W.clebsch_gordan(a,b,c,p,q,r)
class Wigner3j(_Coupling):
    # Arguments are (j1, m1, j2, m2, j3, m3), unlike wigner_3j(j1, j2, j3, m1, m2, m3).
    def doit(self,**hints):
        if self.free_symbols:return self
        a,p,b,q,c,r=self.args;return W.wigner_3j(a,b,c,p,q,r)
class Wigner6j(_Coupling):
    def doit(self,**hints):return self if self.free_symbols else W.wigner_6j(*self.args)
class Wigner9j(_Coupling):
    nargs=9
    def doit(self,**hints):return self if self.free_symbols else W.wigner_9j(*self.args)


def cg_simp(e):
    if isinstance(e,_Coupling):return e.doit()
    expression=S.sympify(e);replacements={};stack=[expression]
    while stack:
        node=stack.pop()
        cls={'CG':CG,'Wigner3j':Wigner3j,'Wigner6j':Wigner6j,'Wigner9j':Wigner9j}.get(node.class_name)
        if cls is not None and not node.free_symbols:replacements[node]=cls(*node.args).doit()
        else:stack.extend(a for a in node.args if isinstance(a,S.Expr))
    return expression.xreplace(replacements)


__all__=['CG','Wigner3j','Wigner6j','Wigner9j','cg_simp']
