"""Excitation-to-activation dynamics for musculotendon systems."""
from abc import ABC, abstractmethod
import symox as S
from ..vector import dynamicsymbols


def _column(values):return S.Matrix(values) if values else S.zeros(0, 1)
class ActivationBase(ABC):
    def __init__(self, name):
        if not isinstance(name, str) or not name:raise ValueError('activation name must be a nonempty string')
        self.name = name
        self._e = dynamicsymbols('e_'+name)
        self._a = dynamicsymbols('a_'+name)
    @property
    def e(self):return self._e
    @property
    def a(self):return self._a
    @property
    def r(self):return _column([self.e])
    @property
    def input_vars(self):return self.r
    @property
    def state_vars(self):return self.x
    @property
    def constants(self):return self.p
    @property
    def M(self):return S.eye(self.x.rows)
    @property
    def F(self):return self.rhs()
    @abstractmethod
    def rhs(self):raise NotImplementedError
class ZerothOrderActivation(ActivationBase):
    def __init__(self, name):super().__init__(name);self._a = self._e
    @property
    def x(self):return S.zeros(0, 1)
    @property
    def p(self):return S.zeros(0, 1)
    def rhs(self):return S.zeros(0, 1)
class FirstOrderActivationDeGroote2016(ActivationBase):
    def __init__(self, name, activation_time_constant=None, deactivation_time_constant=None, smoothing_rate=None):
        super().__init__(name)
        # SymPy names the default constants without assumptions.
        self.tau_a = S.Symbol('tau_a_'+name) if activation_time_constant is None else S.sympify(activation_time_constant)
        self.tau_d = S.Symbol('tau_d_'+name) if deactivation_time_constant is None else S.sympify(deactivation_time_constant)
        self.b = S.Symbol('b_'+name) if smoothing_rate is None else S.sympify(smoothing_rate)
    @classmethod
    def with_defaults(cls, name):return cls(name, .015, .060, 10.)
    @property
    def activation_time_constant(self):return self.tau_a
    @property
    def deactivation_time_constant(self):return self.tau_d
    @property
    def smoothing_rate(self):return self.b
    @property
    def x(self):return _column([self.a])
    @property
    def p(self):return _column([c for c in (self.tau_a, self.tau_d, self.b) if c.free_symbols])
    def rhs(self):
        switch = (1+S.tanh(self.b*(self.e-self.a)))/2
        factor = (1+3*self.a)/2
        return _column([(self.e-self.a)*(switch/(self.tau_a*factor)+(1-switch)*factor/self.tau_d)])
