"""Immutable propositional literals for the known-predicate CNF inventory."""
from dataclasses import dataclass

@dataclass(frozen=True)
class Literal:
    lit:object
    is_Not:bool=False
    def __invert__(self):return type(self)(self.lit,not self.is_Not)
    def __repr__(self):return 'Literal(%s, %s)'%(self.lit,self.is_Not)
    def rcall(self,*args):return type(self)(self.lit(*args),self.is_Not)
