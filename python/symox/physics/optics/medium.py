"""Homogeneous linear optical media."""
import symox as S
from ..units import vacuum_permittivity as epsilon0,vacuum_permeability as mu0,speed_of_light


class Medium:
    def __init__(self,name,permittivity=None,permeability=None,n=None):
        self.name=str(name);self._epsilon=S.sympify(permittivity) if permittivity is not None else None;self._mu=S.sympify(permeability) if permeability is not None else S.sympify(mu0)
        self._n=S.sympify(n) if n is not None else None
        if self._n is None and self._epsilon is None:self._n=S.Integer(1)
        if self._epsilon is None:self._epsilon=self._n**2*S.sympify(epsilon0)*S.sympify(mu0)/self._mu
        if self._n is None:self._n=S.sqrt(self._epsilon*self._mu/(S.sympify(epsilon0)*S.sympify(mu0)))
        if self._n.is_nonpositive is True:raise ValueError('refractive index must be positive')
    def __repr__(self):return self.name
    @property
    def refractive_index(self):return self._n
    @property
    def permittivity(self):return self._epsilon
    @property
    def permeability(self):return self._mu
    @property
    def speed(self):return S.sympify(speed_of_light)/self._n
    @property
    def intrinsic_impedance(self):return S.sqrt(self._mu/self._epsilon)


__all__=['Medium']
