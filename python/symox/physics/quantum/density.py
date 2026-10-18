"""Finite mixed states as weighted outer products."""
import symox as S
from ._core import Dagger,represent,qapply


class Density:
    def __init__(self,*args):self.args=tuple((state,S.sympify(weight)) for state,weight in args)
    def doit(self,**hints):return sum((weight*state*Dagger(state) for state,weight in self.args),S.Integer(0))
    def _symox_(self):return S.sympify(self.doit())
    def entropy(self):
        matrix=self.to_matrix();values=matrix.eigenvals();return -S.Add(*(multiplicity*value*S.log(value) for value,multiplicity in values.items() if value!=0))
    def to_matrix(self,**options):return sum((weight*represent(state,**options)*represent(Dagger(state),**options) for state,weight in self.args),S.zeros(*represent(self.args[0][0],**options).shape[:1]*2))
    def apply_op(self,operator):return Density(*[(qapply(operator*state),weight) for state,weight in self.args])
