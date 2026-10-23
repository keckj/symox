"""Monochromatic transverse waves and coherent superposition, ported from SymPy 1.14."""
__all__ = ['TWave']
from symox.core.basic import Basic
from symox._expr_facade import Expr
from symox.core.function import Derivative, Function
from symox.core.numbers import Number, pi, I
from symox.core.singleton import S
from symox.core.symbol import Symbol, symbols
from symox.core.sympify import _sympify, sympify
from symox.functions.elementary.exponential import exp
from symox.functions.elementary.miscellaneous import sqrt
from symox.functions.elementary.trigonometric import atan2, cos, sin
from symox.physics.units import speed_of_light, meter, second
c = speed_of_light.convert_to(meter / second)

class TWave(Expr):
    _symox_user_hooks = True
    __module__ = 'sympy.physics.optics.waves'

    def __new__(cls, amplitude, frequency=None, phase=S.Zero, time_period=None, n=Symbol('n')):
        if time_period is not None:
            time_period = _sympify(time_period)
            _frequency = S.One / time_period
        if frequency is not None:
            frequency = _sympify(frequency)
            _time_period = S.One / frequency
            if time_period is not None:
                if frequency != S.One / time_period:
                    raise ValueError('frequency and time_period should be consistent.')
        if frequency is None and time_period is None:
            raise ValueError('Either frequency or time period is needed.')
        if frequency is None:
            frequency = _frequency
        if time_period is None:
            time_period = _time_period
        amplitude = _sympify(amplitude)
        phase = _sympify(phase)
        n = sympify(n)
        obj = Expr.__new__(cls, amplitude, frequency, phase, time_period, n)
        return obj

    @property
    def amplitude(self):
        return self.args[0]

    @property
    def frequency(self):
        return self.args[1]

    @property
    def phase(self):
        return self.args[2]

    @property
    def time_period(self):
        return self.args[3]

    @property
    def n(self):
        return self.args[4]

    @property
    def wavelength(self):
        return c / (self.frequency * self.n)

    @property
    def speed(self):
        return self.wavelength * self.frequency

    @property
    def angular_velocity(self):
        return 2 * pi * self.frequency

    @property
    def wavenumber(self):
        return 2 * pi / self.wavelength

    def __str__(self):
        from symox.printing import sstr
        return type(self).__name__ + sstr(self.args)
    __repr__ = __str__

    def __add__(self, other):
        if isinstance(other, TWave):
            if self.frequency == other.frequency and self.wavelength == other.wavelength:
                return TWave(sqrt(self.amplitude ** 2 + other.amplitude ** 2 + 2 * self.amplitude * other.amplitude * cos(self.phase - other.phase)), self.frequency, atan2(self.amplitude * sin(self.phase) + other.amplitude * sin(other.phase), self.amplitude * cos(self.phase) + other.amplitude * cos(other.phase)))
            else:
                raise NotImplementedError('Interference of waves with different frequencies has not been implemented.')
        else:
            raise TypeError(type(other).__name__ + ' and TWave objects cannot be added.')

    def __mul__(self, other):
        other = sympify(other)
        if isinstance(other, Number):
            return TWave(self.amplitude * other, *self.args[1:])
        else:
            raise TypeError(type(other).__name__ + ' and TWave objects cannot be multiplied.')

    def __sub__(self, other):
        return self.__add__(-1 * other)

    def __neg__(self):
        return self.__mul__(-1)

    def __radd__(self, other):
        return self.__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rsub__(self, other):
        return (-self).__radd__(other)

    def _eval_rewrite_as_sin(self, *args, **kwargs):
        return self.amplitude * sin(self.wavenumber * Symbol('x') - self.angular_velocity * Symbol('t') + self.phase + pi / 2, evaluate=False)

    def _eval_rewrite_as_cos(self, *args, **kwargs):
        return self.amplitude * cos(self.wavenumber * Symbol('x') - self.angular_velocity * Symbol('t') + self.phase)

    def _eval_rewrite_as_pde(self, *args, **kwargs):
        (mu, epsilon, x, t) = symbols('mu, epsilon, x, t')
        E = Function('E')
        return Derivative(E(x, t), x, 2) + mu * epsilon * Derivative(E(x, t), t, 2)

    def _eval_rewrite_as_exp(self, *args, **kwargs):
        return self.amplitude * exp(I * (self.wavenumber * Symbol('x') - self.angular_velocity * Symbol('t') + self.phase))
