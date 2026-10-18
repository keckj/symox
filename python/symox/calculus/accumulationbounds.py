"""Real accumulation intervals used for oscillatory sequence limits."""


class AccumulationBounds:
    __slots__=('_args',)
    def __new__(cls,min,max):
        from .. import sympify,S
        lower,upper=map(sympify,(min,max))
        if lower.is_extended_real is not True or upper.is_extended_real is not True or lower in (S.NaN,S.ComplexInfinity) or upper in (S.NaN,S.ComplexInfinity):
            raise ValueError('Accumulation bounds must be extended real numbers')
        if (upper-lower).is_negative is True:
            raise ValueError('The lower accumulation bound must not exceed the upper')
        if lower==upper:return lower
        result=object.__new__(cls)
        result._args=(lower,upper)
        return result
    @property
    def args(self):return self._args
    @property
    def func(self):return type(self)
    @property
    def min(self):return self.args[0]
    @property
    def max(self):return self.args[1]
    @property
    def delta(self):return self.max-self.min
    @property
    def mid(self):return (self.max+self.min)/2
    @property
    def free_symbols(self):return self.min.free_symbols|self.max.free_symbols
    def __repr__(self):return 'AccumBounds(%s, %s)' % self.args
    def __hash__(self):return hash((type(self),self.args))
    def __eq__(self,other):return isinstance(other,AccumulationBounds) and self.args==other.args
    def __reduce__(self):return type(self),self.args
    def __neg__(self):return type(self)(-self.max,-self.min)
    def __add__(self,other):
        from .. import sympify
        if isinstance(other,AccumulationBounds):return type(self)(self.min+other.min,self.max+other.max)
        other=sympify(other)
        return type(self)(self.min+other,self.max+other)
    __radd__=__add__
    def __sub__(self,other):return self+-other
    def __rsub__(self,other):return -self+other
    def __mul__(self,other):
        from .. import Min,Max,sympify
        bounds=other.args if isinstance(other,AccumulationBounds) else (sympify(other),)*2
        products=[0 if a==0 or b==0 else a*b for a in self.args for b in bounds]
        return type(self)(Min(*products),Max(*products))
    __rmul__=__mul__
    def __truediv__(self,other):
        if isinstance(other,AccumulationBounds):return self*other**-1
        from .. import sympify
        return self*(1/sympify(other))
    def __rtruediv__(self,other):return other*self**-1
    def __pow__(self,exponent):
        from .. import sympify,Min,Max,oo,S
        exponent=sympify(exponent)
        if not exponent.is_Integer:
            raise NotImplementedError('Accumulation interval powers require an integer exponent')
        power=int(exponent)
        if not power:return S.One
        if power<0:
            if self.min.is_positive is True or self.max.is_negative is True:
                inverse=type(self)(1/self.max,1/self.min)
            elif self.min==0:inverse=type(self)(1/self.max,oo)
            elif self.max==0:inverse=type(self)(-oo,1/self.min)
            else:inverse=type(self)(-oo,oo)
            return inverse if power==-1 else inverse**(-power)
        lower,upper=self.min**power,self.max**power
        if power%2:return type(self)(lower,upper)
        minimum=Min(lower,upper)
        if self.min.is_extended_nonpositive is True and self.max.is_extended_nonnegative is True:
            minimum=S.Zero
        return type(self)(minimum,Max(lower,upper))
    def subs(self,*args,**kwargs):return type(self)(*(value.subs(*args,**kwargs) for value in self.args))
    def xreplace(self,mapping):
        return mapping[self] if self in mapping else type(self)(*(value.xreplace(mapping) for value in self.args))
    def to_sympy(self):
        from sympy import AccumBounds
        return AccumBounds(self.min.to_sympy(),self.max.to_sympy())


AccumBounds=AccumulationBounds
__all__=['AccumBounds','AccumulationBounds']
