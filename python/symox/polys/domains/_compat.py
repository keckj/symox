"""Legacy backend names constructed on symox's native coefficient domains."""
import operator
from fractions import Fraction

from ._types import IntegerRing,RationalField,FiniteField


class PythonIntegerRing(IntegerRing):
    """Python-domain constructor syntax returning the native integer ring."""


class GMPYIntegerRing(IntegerRing):
    """GMPY-domain constructor syntax returning the native integer ring."""


class PythonRationalField(RationalField):
    """Python-domain constructor syntax returning the native rational field."""


class GMPYRationalField(RationalField):
    """GMPY-domain constructor syntax returning the native rational field."""


class PythonFiniteField(FiniteField):
    def __new__(cls,mod,symmetric=True):
        return super().__new__(cls,mod,symmetric=symmetric)


class GMPYFiniteField(FiniteField):
    def __new__(cls,mod,symmetric=True):
        return super().__new__(cls,mod,symmetric=symmetric)


class PythonRational(Fraction):
    """Exact Python integer rational interoperating with native coefficients."""
    __slots__=()
    def __new__(cls,numerator,denominator=None):
        if denominator is not None:
            if not isinstance(numerator,int) or not isinstance(denominator,int):
                raise TypeError('PythonRational numerator and denominator must be Python integers')
            return super().__new__(cls,numerator,denominator)
        if hasattr(numerator,'_symox_') and not isinstance(numerator,PythonRational):
            numerator=numerator._symox_()
        if getattr(numerator,'is_Rational',False):
            return super().__new__(cls,int(numerator.p),int(numerator.q))
        return super().__new__(cls,numerator)
    def __repr__(self):return 'MPQ(%s,%s)' % (self.numerator,self.denominator)
    def __pos__(self):return self
    def __neg__(self):return type(self)(-self.numerator,self.denominator)
    def __abs__(self):return type(self)(abs(self.numerator),self.denominator)
    def _binary(self,other,operation):
        if not isinstance(other,(int,Fraction)):return NotImplemented
        if isinstance(other,Fraction):other=Fraction(other.numerator,other.denominator)
        result=operation(Fraction(self.numerator,self.denominator),other)
        return type(self)(result.numerator,result.denominator)
    def __add__(self,other):return self._binary(other,operator.add)
    __radd__=__add__
    def __sub__(self,other):return self._binary(other,operator.sub)
    def __rsub__(self,other):return -self+other
    def __mul__(self,other):return self._binary(other,operator.mul)
    __rmul__=__mul__
    def __truediv__(self,other):return self._binary(other,operator.truediv)
    def __rtruediv__(self,other):
        if not isinstance(other,(int,Fraction)):return NotImplemented
        return type(self)(other)/self
    def __pow__(self,power):
        power=operator.index(power)
        result=Fraction(self.numerator,self.denominator)**power
        return type(self)(result.numerator,result.denominator)
    def __floordiv__(self,other):raise TypeError('Rational domain elements do not support floor division')
    __rfloordiv__=__floordiv__
    def __mod__(self,other):raise TypeError('Rational domain elements do not support remainder')
    __rmod__=__mod__
    def __reduce__(self):return type(self),(self.numerator,self.denominator)
    def _symox_(self):
        from ... import Rational
        return Rational(self.numerator,self.denominator)
    def to_sympy(self):
        from sympy import Rational
        return Rational(self.numerator,self.denominator)
    _sympy_=to_sympy


ZZ_python=PythonIntegerRing
ZZ_gmpy=GMPYIntegerRing
QQ_python=PythonRationalField
QQ_gmpy=GMPYRationalField
FF_python=PythonFiniteField
FF_gmpy=GMPYFiniteField

__all__=['PythonIntegerRing','GMPYIntegerRing','PythonRationalField','GMPYRationalField',
         'PythonFiniteField','GMPYFiniteField','PythonRational','ZZ_python','ZZ_gmpy',
         'QQ_python','QQ_gmpy','FF_python','FF_gmpy']
