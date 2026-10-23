"""De Groote et al. (2016) characteristic curves, evaluated by the native engine."""
from abc import ABC, abstractmethod
from dataclasses import dataclass
import symox as S


class CharacteristicCurveFunction(ABC):
    defaults = ()
    def __init__(self, value, *constants):
        if len(constants) != len(self.defaults):
            raise TypeError(f'{type(self).__name__} requires {len(self.defaults)} curve constants')
        self.args = tuple(map(S.sympify, (value, *constants)))
    @classmethod
    def with_defaults(cls, value):return cls(value, *cls.defaults)
    @abstractmethod
    def _formula(self, x, c):raise NotImplementedError
    def doit(self, **hints):return self._formula(self.args[0], self.args[1:])
    def _symox_(self):return self.doit()
    @property
    def free_symbols(self):return self.doit().free_symbols
    def fdiff(self, argindex=1):
        if not 1 <= argindex <= len(self.args):raise ValueError('invalid argument index')
        dummies = tuple(S.Dummy('curve') for _ in self.args)
        formula = self._formula(dummies[0], dummies[1:])
        return S.diff(formula, dummies[argindex-1]).xreplace(dict(zip(dummies, self.args)))
    def diff(self, *args, **kwargs):return self.doit().diff(*args, **kwargs)
    def subs(self, *args, **kwargs):return type(self)(*(x.subs(*args, **kwargs) for x in self.args))
    def _sympystr(self, printer):
        return '%s(%s)' % (type(self).__name__, ', '.join(printer._print(a) for a in self.args))
    def __repr__(self):return S.sstr(self)
    __str__ = __repr__
    def __mul__(self, other):return self.doit()*S.sympify(other)
    __rmul__ = __mul__
    def __add__(self, other):return self.doit()+S.sympify(other)
    __radd__ = __add__
    def __sub__(self, other):return self.doit()-S.sympify(other)
    def __rsub__(self, other):return S.sympify(other)-self.doit()
    def __truediv__(self, other):return self.doit()/S.sympify(other)
    def __pow__(self, power):return self.doit()**power
    def __neg__(self):return -self.doit()


class FiberForceLengthActiveDeGroote2016(CharacteristicCurveFunction):
    defaults = (.814, 1.06, .162, .0633, .433, .717, -.0299, .2, .1, 1., .354, 0.)
    def _formula(self, x, c):return sum(c[i]*S.exp(-S.Rational(1, 2)*((x-c[i+1])/(c[i+2]+c[i+3]*x))**2) for i in (0, 4, 8))
class FiberForceLengthPassiveDeGroote2016(CharacteristicCurveFunction):
    defaults = (.6, 4.)
    def _formula(self, x, c):return (S.exp(c[1]*(x-1)/c[0])-1)/(S.exp(c[1])-1)
class FiberForceLengthPassiveInverseDeGroote2016(CharacteristicCurveFunction):
    defaults = (.6, 4.)
    def _formula(self, x, c):return 1+c[0]/c[1]*S.log(1+x*(S.exp(c[1])-1))
class FiberForceVelocityDeGroote2016(CharacteristicCurveFunction):
    defaults = (-.318, -8.149, -.374, .886)
    def _formula(self, x, c):return c[0]*S.log(c[1]*x+c[2]+S.sqrt((c[1]*x+c[2])**2+1))+c[3]
class FiberForceVelocityInverseDeGroote2016(CharacteristicCurveFunction):
    defaults = (-.318, -8.149, -.374, .886)
    def _formula(self, x, c):return (S.sinh((x-c[3])/c[0])-c[2])/c[1]
class TendonForceLengthDeGroote2016(CharacteristicCurveFunction):
    defaults = (.2, .995, .25, S.Float('33.93669377311689'))
    def _formula(self, x, c):return c[0]*S.exp(c[3]*(x-c[1]))-c[2]
class TendonForceLengthInverseDeGroote2016(CharacteristicCurveFunction):
    defaults = (.2, .995, .25, S.Float('33.93669377311689'))
    def _formula(self, x, c):return c[1]+S.log((x+c[2])/c[0])/c[3]


@dataclass(frozen=True)
class CharacteristicCurveCollection:
    tendon_force_length: type
    tendon_force_length_inverse: type
    fiber_force_length_passive: type
    fiber_force_length_passive_inverse: type
    fiber_force_length_active: type
    fiber_force_velocity: type
    fiber_force_velocity_inverse: type
