"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'satisfiable')

from abc import ABC, abstractmethod
import symox as S


def literal_symbol(literal):
    literal=S.sympify(literal)
    return literal.args[0] if literal.class_name=='Not' else literal


def pl_true(expr,model=None,deep=False):
    value=S.sympify(expr).subs({symbol:S.true if truth else S.false for symbol,truth in (model or {}).items()})
    if value==S.true:return True
    if value==S.false:return False
    if deep:
        if S.satisfiable(value) is False:return False
        if S.satisfiable(S.Not(value)) is False:return True
    return None

def valid(expr):return S.satisfiable(S.Not(expr)) is False

def entails(expr,formula_set=None):
    formula_set=list(formula_set) if formula_set else []
    formula_set.append(S.Not(expr))
    return not S.satisfiable(S.And(*formula_set))


class KB(ABC):
    def __init__(self,sentence=None):
        self.clauses_=set()
        if sentence:self.tell(sentence)
    @abstractmethod
    def tell(self,sentence):raise NotImplementedError
    @abstractmethod
    def ask(self,query):raise NotImplementedError
    @abstractmethod
    def retract(self,sentence):raise NotImplementedError
    @property
    def clauses(self):
        return list(S.ordered(self.clauses_))


class PropKB(KB):
    def tell(self,sentence):
        for c in S.And.make_args(S.to_cnf(sentence)):self.clauses_.add(c)
    def ask(self,query):return entails(query,self.clauses)
    def retract(self,sentence):
        for c in S.And.make_args(S.to_cnf(sentence)):self.clauses_.discard(c)
__all__ += ['KB','PropKB','entails','literal_symbol','pl_true','valid']
