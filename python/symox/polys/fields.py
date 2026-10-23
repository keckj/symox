"""Rational function fields with native cancellation and coefficient domains."""
import inspect
import operator
from .rings import PolyRing, PolyElement, _inferred
from .orderings import lex

__all__ = ['field','xfield','vfield','sfield']


class FracField:
    def __init__(self, symbols, domain, order=lex):
        self.ring = PolyRing(symbols,domain,order)
        self.symbols,self.ngens,self.domain,self.order = self.ring.symbols,self.ring.ngens,self.ring.domain,self.ring.order
        self._native = self.domain.frac_field(*self.symbols,order=str(self.order)) if self.ngens else self.domain.get_field()
        self.gens = tuple(self(s) for s in self.symbols)
        self.zero,self.one = self(0),self(1)
        self.dtype = self.__call__
        for symbol,gen in zip(self.symbols,self.gens):
            name=getattr(symbol,'name',None)
            if name and not hasattr(self,name):
                setattr(self,name,gen)

    def __repr__(self):
        return 'Rational function field in %s over %s with %s order' % (', '.join(map(str,self.symbols)),self.domain,self.order)

    def __eq__(self, other):
        return isinstance(other,FracField) and self.ring==other.ring

    def __hash__(self):
        return hash(('FracField',self.ring))

    def __call__(self, numer=0, denom=None):
        if denom is None and isinstance(numer,tuple) and len(numer)==2:
            numer,denom = numer
        if isinstance(numer,FracElement):
            if numer.field==self and denom is None:
                return numer
            numer=numer.as_expr()
        elif isinstance(numer,(PolyElement,dict)):
            numer=self.ring(numer).as_expr()
        value=self._native.convert(numer)
        if denom is not None:
            if isinstance(denom,(PolyElement,dict)):
                denom=self.ring(denom).as_expr()
            value=value/self._native.convert(denom)
        return FracElement(self,value)

    from_expr = __call__
    new = __call__

    def to_domain(self):
        return self._native

    def to_ring(self):
        return self.ring

    def is_element(self,value):
        return isinstance(value,FracElement) and value.field==self


class FracElement:
    # Expressions on the left defer to the element, as the reference
    # ``Expr`` binary methods return ``NotImplemented`` for it.
    _op_priority = 10.1
    _symox_ring_element = True

    def __init__(self,field,value):
        self.field,self._value=field,value

    def parent(self):
        return self.field

    def as_expr(self):
        return self._value.as_expr()
    _symox_=as_expr

    def to_sympy(self):
        return self.as_expr().to_sympy()

    def __repr__(self):
        # The reference ``StrPrinter._print_FracElement``.
        from ..printing.precedence import PRECEDENCE
        from .rings import _parenthesize
        if self.denom == 1:
            return repr(self.numer)
        return (_parenthesize(self.numer, PRECEDENCE["Mul"], True) + "/" +
                _parenthesize(self.denom, PRECEDENCE["Atom"], True))

    def __hash__(self):
        return hash((self.field,self._value))

    def __eq__(self,other):
        try:
            return self._value==self.field(other)._value
        except (TypeError,ValueError):
            return False

    def __bool__(self):
        return self._value!=self.field._native.zero

    def __neg__(self):
        return self.field(-self._value)

    def __add__(self,other):
        return self.field(self._value+self.field(other)._value)
    __radd__=__add__

    def __sub__(self,other):
        return self+(-self.field(other))

    def __rsub__(self,other):
        return self.field(other)-self

    def __mul__(self,other):
        return self.field(self._value*self.field(other)._value)
    __rmul__=__mul__

    def __truediv__(self,other):
        return self.field(self._value/self.field(other)._value)

    def __rtruediv__(self,other):
        return self.field(other)/self

    def __pow__(self,n):
        return self.field(self._value**operator.index(n))

    @property
    def numer(self):
        return self.field.ring(self._value.numer.as_expr())

    @property
    def denom(self):
        return self.field.ring(self._value.denom.as_expr())

    def to_poly(self):
        if self.denom!=1:
            raise ValueError('rational function is not a polynomial')
        return self.numer

    def diff(self,x):
        return self.field(self.as_expr().diff(self.field.symbols[self.field.ring.index(x)]))

    def evaluate(self,x,a=None):
        assignments=x if a is None else [(x,a)]
        substitutions={self.field.symbols[self.field.ring.index(g)]:v.as_expr() if isinstance(v,(PolyElement,FracElement)) else v for g,v in assignments}
        remaining=tuple(g for g in self.field.symbols if g not in substitutions)
        value=self.as_expr().subs(substitutions,simultaneous=True)
        return FracField(remaining,self.field.domain,self.field.order)(value) if remaining else self.field.domain.get_field().convert(value)


def field(symbols, domain, order=lex):
    F=FracField(symbols,domain,order)
    return (F,)+F.gens


def xfield(symbols, domain, order=lex):
    F=FracField(symbols,domain,order)
    return F,F.gens


def vfield(symbols, domain, order=lex):
    F=FracField(symbols,domain,order)
    inspect.currentframe().f_back.f_globals.update({str(s):g for s,g in zip(F.symbols,F.gens)})
    return F


def sfield(exprs,*symbols,**options):
    values,symbols,domain,order,single=_inferred(exprs,symbols,options,True)
    F=FracField(symbols,domain,order)
    result=[F(v) for v in values]
    return F,result[0] if single else result
