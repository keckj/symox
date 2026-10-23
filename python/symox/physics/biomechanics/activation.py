"""Activation dynamics for musculotendon models, ported from SymPy 1.14."""
from abc import ABC, abstractmethod
from functools import cached_property
from symox.core.symbol import Symbol
from symox.core.numbers import Float, Integer, Rational
from symox.functions.elementary.hyperbolic import tanh
from symox.matrices.dense import MutableDenseMatrix as Matrix, zeros
from symox.physics.biomechanics._mixin import _NamedMixin
from symox.physics.mechanics import dynamicsymbols
__all__ = ['ActivationBase', 'FirstOrderActivationDeGroote2016', 'ZerothOrderActivation']

class ActivationBase(ABC, _NamedMixin):

    def __init__(self, name):
        self.name = str(name)
        self._e = dynamicsymbols(f'e_{name}')
        self._a = dynamicsymbols(f'a_{name}')

    @classmethod
    @abstractmethod
    def with_defaults(cls, name):
        pass

    @property
    def excitation(self):
        return self._e

    @property
    def e(self):
        return self._e

    @property
    def activation(self):
        return self._a

    @property
    def a(self):
        return self._a

    @property
    @abstractmethod
    def order(self):
        pass

    @property
    @abstractmethod
    def state_vars(self):
        pass

    @property
    @abstractmethod
    def x(self):
        pass

    @property
    @abstractmethod
    def input_vars(self):
        pass

    @property
    @abstractmethod
    def r(self):
        pass

    @property
    @abstractmethod
    def constants(self):
        pass

    @property
    @abstractmethod
    def p(self):
        pass

    @property
    @abstractmethod
    def M(self):
        pass

    @property
    @abstractmethod
    def F(self):
        pass

    @abstractmethod
    def rhs(self):
        pass

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        if self.name != other.name:
            return False
        return True

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name!r})'

class ZerothOrderActivation(ActivationBase):

    def __init__(self, name):
        super().__init__(name)
        self._a = self._e

    @classmethod
    def with_defaults(cls, name):
        return cls(name)

    @property
    def order(self):
        return 0

    @property
    def state_vars(self):
        return zeros(0, 1)

    @property
    def x(self):
        return zeros(0, 1)

    @property
    def input_vars(self):
        return Matrix([self._e])

    @property
    def r(self):
        return Matrix([self._e])

    @property
    def constants(self):
        return zeros(0, 1)

    @property
    def p(self):
        return zeros(0, 1)

    @property
    def M(self):
        return Matrix([])

    @property
    def F(self):
        return zeros(0, 1)

    def rhs(self):
        return zeros(0, 1)

class FirstOrderActivationDeGroote2016(ActivationBase):

    def __init__(self, name, activation_time_constant=None, deactivation_time_constant=None, smoothing_rate=None):
        super().__init__(name)
        self.activation_time_constant = activation_time_constant
        self.deactivation_time_constant = deactivation_time_constant
        self.smoothing_rate = smoothing_rate

    @classmethod
    def with_defaults(cls, name):
        tau_a = Float('0.015')
        tau_d = Float('0.060')
        b = Float('10.0')
        return cls(name, tau_a, tau_d, b)

    @property
    def activation_time_constant(self):
        return self._tau_a

    @activation_time_constant.setter
    def activation_time_constant(self, tau_a):
        if hasattr(self, '_tau_a'):
            msg = f"Can't set attribute `activation_time_constant` to {repr(tau_a)} as it is immutable and already has value {self._tau_a}."
            raise AttributeError(msg)
        self._tau_a = Symbol(f'tau_a_{self.name}') if tau_a is None else tau_a

    @property
    def tau_a(self):
        return self._tau_a

    @property
    def deactivation_time_constant(self):
        return self._tau_d

    @deactivation_time_constant.setter
    def deactivation_time_constant(self, tau_d):
        if hasattr(self, '_tau_d'):
            msg = f"Can't set attribute `deactivation_time_constant` to {repr(tau_d)} as it is immutable and already has value {self._tau_d}."
            raise AttributeError(msg)
        self._tau_d = Symbol(f'tau_d_{self.name}') if tau_d is None else tau_d

    @property
    def tau_d(self):
        return self._tau_d

    @property
    def smoothing_rate(self):
        return self._b

    @smoothing_rate.setter
    def smoothing_rate(self, b):
        if hasattr(self, '_b'):
            msg = f"Can't set attribute `smoothing_rate` to {b!r} as it is immutable and already has value {self._b!r}."
            raise AttributeError(msg)
        self._b = Symbol(f'b_{self.name}') if b is None else b

    @property
    def b(self):
        return self._b

    @property
    def order(self):
        return 1

    @property
    def state_vars(self):
        return Matrix([self._a])

    @property
    def x(self):
        return Matrix([self._a])

    @property
    def input_vars(self):
        return Matrix([self._e])

    @property
    def r(self):
        return Matrix([self._e])

    @property
    def constants(self):
        constants = [self._tau_a, self._tau_d, self._b]
        symbolic_constants = [c for c in constants if not c.is_number]
        return Matrix(symbolic_constants) if symbolic_constants else zeros(0, 1)

    @property
    def p(self):
        constants = [self._tau_a, self._tau_d, self._b]
        symbolic_constants = [c for c in constants if not c.is_number]
        return Matrix(symbolic_constants) if symbolic_constants else zeros(0, 1)

    @property
    def M(self):
        return Matrix([Integer(1)])

    @property
    def F(self):
        return Matrix([self._da_eqn])

    def rhs(self):
        return Matrix([self._da_eqn])

    @cached_property
    def _da_eqn(self):
        HALF = Rational(1, 2)
        a0 = HALF * tanh(self._b * (self._e - self._a))
        a1 = HALF + Rational(3, 2) * self._a
        a2 = (HALF + a0) / (self._tau_a * a1)
        a3 = a1 * (HALF - a0) / self._tau_d
        activation_dynamics_equation = (a2 + a3) * (self._e - self._a)
        return activation_dynamics_equation

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        self_attrs = (self.name, self.tau_a, self.tau_d, self.b)
        other_attrs = (other.name, other.tau_a, other.tau_d, other.b)
        if self_attrs == other_attrs:
            return True
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name!r}, activation_time_constant={self.tau_a!r}, deactivation_time_constant={self.tau_d!r}, smoothing_rate={self.b!r})'
