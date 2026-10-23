"""Optical media, ported from SymPy 1.14."""
from symox.physics.units import second, meter, kilogram, ampere
__all__ = ['Medium']
from symox.core.basic import Basic
from symox.core.symbol import Str
from symox.core.sympify import _sympify
from symox.functions.elementary.miscellaneous import sqrt
from symox.physics.units import speed_of_light, u0, e0
c = speed_of_light.convert_to(meter / second)
_e0mksa = e0.convert_to(ampere ** 2 * second ** 4 / (kilogram * meter ** 3))
_u0mksa = u0.convert_to(meter * kilogram / (ampere ** 2 * second ** 2))

class Medium(Basic):

    def __new__(cls, name, permittivity=None, permeability=None, n=None):
        if not isinstance(name, Str):
            name = Str(name)
        permittivity = _sympify(permittivity) if permittivity is not None else permittivity
        permeability = _sympify(permeability) if permeability is not None else permeability
        n = _sympify(n) if n is not None else n
        if n is not None:
            if permittivity is not None and permeability is None:
                permeability = n ** 2 / (c ** 2 * permittivity)
                return MediumPP(name, permittivity, permeability)
            elif permeability is not None and permittivity is None:
                permittivity = n ** 2 / (c ** 2 * permeability)
                return MediumPP(name, permittivity, permeability)
            elif permittivity is not None and permittivity is not None:
                raise ValueError('Specifying all of permittivity, permeability, and n is not allowed')
            else:
                return MediumN(name, n)
        elif permittivity is not None and permeability is not None:
            return MediumPP(name, permittivity, permeability)
        elif permittivity is None and permeability is None:
            return MediumPP(name, _e0mksa, _u0mksa)
        else:
            raise ValueError('Arguments are underspecified. Either specify n or any two of permittivity, permeability, and n')

    @property
    def name(self):
        return self.args[0]

    @property
    def speed(self):
        return c / self.n

    @property
    def refractive_index(self):
        return c / self.speed

class MediumN(Medium):

    def __new__(cls, name, n):
        obj = super(Medium, cls).__new__(cls, name, n)
        return obj

    @property
    def n(self):
        return self.args[1]

class MediumPP(Medium):

    def __new__(cls, name, permittivity, permeability):
        obj = super(Medium, cls).__new__(cls, name, permittivity, permeability)
        return obj

    @property
    def intrinsic_impedance(self):
        return sqrt(self.permeability / self.permittivity)

    @property
    def permittivity(self):
        return self.args[1]

    @property
    def permeability(self):
        return self.args[2]

    @property
    def n(self):
        return c * sqrt(self.permittivity * self.permeability)
