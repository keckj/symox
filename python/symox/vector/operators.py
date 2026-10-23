"""Differential operators in orthogonal three-dimensional coordinates."""
from .. import Integer, sympify, simplify
from ._basis import Vector, Dyadic, _Held
from .coordsysrect import _scalar_systems


def _get_coord_systems(expression):
    result = set()
    if isinstance(expression,Vector):
        result.update(b.system for b in expression.components)
    elif isinstance(expression,Dyadic):
        result.update(b.vector1.system for b in expression.components)
        result.update(b.vector2.system for b in expression.components)
    for symbol in getattr(expression,'free_symbols',()):
        if symbol in _scalar_systems:
            result.add(_scalar_systems[symbol][0])
    return result


def _system(expression):
    systems = _get_coord_systems(expression)
    if not systems:
        return None
    if len(systems) == 1:
        return next(iter(systems))
    roots = {c._root for c in systems}
    if len(roots) != 1:
        raise ValueError('field contains disconnected coordinate systems')
    return next(iter(roots))


def gradient(scalar_field,doit=True):
    if isinstance(scalar_field,(Vector,Dyadic)):
        raise TypeError('gradient requires a scalar field')
    scalar_field = sympify(scalar_field)
    # Reference gradient: a field free of base scalars has the zero
    # gradient, whether or not doit is requested.
    if not _get_coord_systems(scalar_field):
        return Vector.zero
    if not doit:
        return Gradient(scalar_field)
    system = _system(scalar_field)
    from .functions import express
    scalar_field = express(scalar_field,system,variables=True)
    return sum((simplify(scalar_field.diff(x)/h)*b for x,h,b in zip(system.base_scalars(),system.lame_coefficients(),system.base_vectors())),Vector.zero)


def divergence(vect,doit=True):
    if not isinstance(vect,Vector):
        raise TypeError('divergence requires a vector field')
    if not doit:
        return Divergence(vect)
    system = _system(vect)
    if system is None:
        return Integer(0)
    from .functions import express
    values = express(vect,system,variables=True).to_matrix(system)
    h = system.lame_coefficients()
    product = h[0]*h[1]*h[2]
    return simplify(sum((product*values[i]/h[i]).diff(x) for i,x in enumerate(system.base_scalars()))/product)


def curl(vect,doit=True):
    if not isinstance(vect,Vector):
        raise TypeError('curl requires a vector field')
    if not doit:
        return Curl(vect)
    system = _system(vect)
    if system is None:
        return Vector.zero
    from .functions import express
    a = express(vect,system,variables=True).to_matrix(system)
    x,h,b = system.base_scalars(),system.lame_coefficients(),system.base_vectors()
    return sum((simplify(((h[(i+2)%3]*a[(i+2)%3]).diff(x[(i+1)%3])-(h[(i+1)%3]*a[(i+1)%3]).diff(x[(i+2)%3]))/(h[(i+1)%3]*h[(i+2)%3]))*b[i] for i in range(3)),Vector.zero)


class Gradient(_Held):
    def __init__(self,expr):
        if isinstance(expr,(Vector,Dyadic)):
            raise TypeError('Gradient requires a scalar field')
        super().__init__(sympify(expr))

    def doit(self,**kwargs):
        return gradient(self.args[0])


class Divergence(_Held):
    def __init__(self,expr):
        if not isinstance(expr,Vector):
            raise TypeError('Divergence requires a vector')
        super().__init__(expr)

    def doit(self,**kwargs):
        return divergence(self.args[0])


class Curl(Divergence):
    def doit(self,**kwargs):
        return curl(self.args[0])


class Laplacian(_Held):
    def doit(self,**kwargs):
        from .functions import laplacian
        return laplacian(self.args[0])
