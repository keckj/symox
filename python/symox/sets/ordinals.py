"""Exact ordinal arithmetic in finite Cantor normal form, below epsilon-zero."""
from functools import total_ordering
import operator


@total_ordering
class OmegaPower:
    def __init__(self,a,b):
        self.exp=Ordinal.convert(a)
        self.mult=operator.index(b)
        if self.mult<=0:
            raise TypeError('multiplicity must be a positive integer')

    @property
    def args(self):
        return self.exp,self.mult

    def __eq__(self,other):
        if not isinstance(other,OmegaPower):
            try:
                other=OmegaPower(0,other)
            except (TypeError,ValueError):
                return NotImplemented
        return self.args==other.args

    def __lt__(self,other):
        if not isinstance(other,OmegaPower):
            other=OmegaPower(0,other)
        return self.args<other.args

    def __hash__(self):
        return hash(self.mult) if self.exp==0 else hash(self.args)

    def __repr__(self):
        return 'OmegaPower(%s, %s)' % self.args

    def to_sympy(self):
        from sympy import OmegaPower
        return OmegaPower(self.exp.to_sympy(),self.mult)


@total_ordering
class Ordinal:
    def __init__(self,*terms):
        if any(not isinstance(t,OmegaPower) for t in terms):
            raise TypeError('Cantor terms must be OmegaPower values')
        if any(a.exp<=b.exp for a,b in zip(terms,terms[1:])):
            raise ValueError('Cantor powers must be strictly decreasing')
        self.terms=tuple(terms)

    @property
    def args(self):
        return self.terms

    @classmethod
    def convert(cls,value):
        if isinstance(value,Ordinal):
            return value
        value=operator.index(value)
        if value<0:
            raise TypeError('ordinals are nonnegative')
        return Ordinal() if value==0 else Ordinal(OmegaPower(Ordinal(),value))

    def __eq__(self,other):
        try:
            return self.terms==Ordinal.convert(other).terms
        except (TypeError,ValueError):
            return NotImplemented

    def __lt__(self,other):
        return self.terms<Ordinal.convert(other).terms

    def __hash__(self):
        if not self.terms:
            return hash(0)
        if len(self.terms)==1 and self.terms[0].exp==0:
            return hash(self.terms[0].mult)
        return hash(self.terms)

    @property
    def leading_term(self):
        if not self.terms:
            raise ValueError('zero has no leading term')
        return self.terms[0]

    @property
    def trailing_term(self):
        if not self.terms:
            raise ValueError('zero has no trailing term')
        return self.terms[-1]

    @property
    def degree(self):
        return self.leading_term.exp

    @property
    def is_successor_ordinal(self):
        return bool(self.terms) and self.trailing_term.exp==0

    @property
    def is_limit_ordinal(self):
        return bool(self.terms) and self.trailing_term.exp!=0

    def __add__(self,other):
        other=Ordinal.convert(other)
        if not other.terms:
            return self
        result=list(self.terms)
        head=other.leading_term
        while result and result[-1].exp<head.exp:
            result.pop()
        if result and result[-1].exp==head.exp:
            previous=result.pop()
            return Ordinal(*result,OmegaPower(head.exp,previous.mult+head.mult),*other.terms[1:])
        return Ordinal(*result,*other.terms)

    def __radd__(self,other):
        return Ordinal.convert(other)+self

    def __mul__(self,other):
        other=Ordinal.convert(other)
        if not self.terms or not other.terms:
            return ord0
        result=ord0
        for term in other.terms:
            if term.exp==0:
                finite=Ordinal(OmegaPower(self.degree,self.leading_term.mult*term.mult),*self.terms[1:])
                result+=finite
            else:
                result+=Ordinal(OmegaPower(self.degree+term.exp,term.mult))
        return result

    def __rmul__(self,other):
        return Ordinal.convert(other)*self

    def __pow__(self,other):
        if self==omega:
            return Ordinal(OmegaPower(other,1))
        n=operator.index(other)
        if n<0:
            raise ValueError('ordinal exponents must be nonnegative')
        result=Ordinal.convert(1)
        value=self
        while n:
            if n&1:
                result=result*value
            value=value*value
            n//=2
        return result

    def __repr__(self):
        if not self.terms:
            return 'ord0'
        parts=[]
        for t in self.terms:
            if t.exp==0:
                parts.append(str(t.mult))
            else:
                text='w' if t.exp==1 else 'w**(%s)'%t.exp
                parts.append(text if t.mult==1 else text+'*%d'%t.mult)
        return ' + '.join(parts)
    __str__=__repr__

    def to_sympy(self):
        from sympy import Ordinal
        return Ordinal(*(t.to_sympy() for t in self.terms))


class OrdinalZero(Ordinal):
    pass


ord0=OrdinalZero()
omega=Ordinal(OmegaPower(1,1))
