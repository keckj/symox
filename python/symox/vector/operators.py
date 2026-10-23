"""Differential operators in orthogonal three-dimensional coordinates."""
from .. import Integer, sympify, simplify, Derivative
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
    systems = _get_coord_systems(scalar_field)
    if len(systems) == 1:
        # SymPy 1.14: unevaluated derivatives over the Lame coefficients.
        system = next(iter(systems))
        h1, h2, h3 = system.lame_coefficients()
        i, j, k = system.base_vectors()
        x, y, z = system.base_scalars()
        result = (Derivative(scalar_field, x)/h1*i + Derivative(scalar_field, y)/h2*j
                  + Derivative(scalar_field, z)/h3*k)
        return result.doit() if doit else result
    split = _split_gradient(scalar_field)
    if split is not None:
        return split
    if not doit:
        return Gradient(scalar_field)
    system = _system(scalar_field)
    from .functions import express
    scalar_field = express(scalar_field,system,variables=True)
    return sum((simplify(scalar_field.diff(x)/h)*b for x,h,b in zip(system.base_scalars(),system.lame_coefficients(),system.base_vectors())),Vector.zero)


def _split_gradient(scalar_field):
    """SymPy's multi-system gradient: split sums and products by system.

    None when a piece still mixes systems (SymPy holds Gradient there; the
    native fallback below evaluates it through a common root instead)."""
    from .. import Add, Mul, S
    if scalar_field.is_Add:
        pieces = [gradient(term) for term in scalar_field.args]
        if all(isinstance(piece,Vector) for piece in pieces):
            return sum(pieces,Vector.zero)
        return None
    if scalar_field.is_Mul:
        groups = {}
        for factor in scalar_field.args:
            key = frozenset(_get_coord_systems(factor))
            groups[key] = groups.get(key,S.One)*factor
        if any(len(key) > 1 for key in groups):
            return None
        return sum((scalar_field/factor*gradient(factor) for factor in groups.values()),Vector.zero)
    return None


def _split_divergence(vect):
    """SymPy's evaluated multi-system divergence, term by term."""
    from .. import Add
    if len(vect.components) > 1:
        pieces = [_split_divergence(c*b) if len(_get_coord_systems(c*b)) > 1 else divergence(c*b)
                  for b,c in vect.components.items()]
        if any(piece is None for piece in pieces):
            return None
        return Add(*pieces)
    (base,scalar), = vect.components.items()
    grad = gradient(scalar)
    if not isinstance(grad,Vector):
        return None
    try:
        projection = base.dot(grad)
    except ValueError:
        # Disconnected systems: SymPy holds Dot(base, grad) here.
        return None
    return projection + scalar*divergence(base)


def divergence(vect,doit=True):
    if not isinstance(vect,Vector):
        raise TypeError('divergence requires a vector field')
    systems = _get_coord_systems(vect)
    if not systems:
        return Integer(0)
    if len(systems) == 1:
        system = next(iter(systems))
        i, j, k = system.base_vectors()
        x, y, z = system.base_scalars()
        h1, h2, h3 = system.lame_coefficients()
        def conditional(component, scalar, coeff_1, coeff_2):
            arg = coeff_1*coeff_2*component
            return Derivative(arg, scalar) if arg else Integer(0)
        result = (conditional(vect.dot(i), x, h2, h3)/(h1*h2*h3)
                  + conditional(vect.dot(j), y, h3, h1)/(h1*h2*h3)
                  + conditional(vect.dot(k), z, h1, h2)/(h1*h2*h3))
        return result.doit() if doit else result
    if doit:
        split = _split_divergence(vect)
        if split is not None:
            return split
    if not doit:
        return Divergence(vect)
    system = _system(vect)
    from .functions import express
    values = express(vect,system,variables=True).to_matrix(system)
    h = system.lame_coefficients()
    product = h[0]*h[1]*h[2]
    return simplify(sum((product*values[i]/h[i]).diff(x) for i,x in enumerate(system.base_scalars()))/product)


def curl(vect,doit=True):
    if not isinstance(vect,Vector):
        raise TypeError('curl requires a vector field')
    systems = _get_coord_systems(vect)
    if not systems:
        return Vector.zero
    if len(systems) == 1:
        system = next(iter(systems))
        i, j, k = system.base_vectors()
        x, y, z = system.base_scalars()
        h1, h2, h3 = system.lame_coefficients()
        vx, vy, vz = vect.dot(i), vect.dot(j), vect.dot(k)
        result = Vector.zero
        result += (Derivative(vz*h3, y) - Derivative(vy*h2, z))*i/(h2*h3)
        result += (Derivative(vx*h1, z) - Derivative(vz*h3, x))*j/(h1*h3)
        result += (Derivative(vy*h2, x) - Derivative(vx*h1, y))*k/(h2*h1)
        return result.doit() if doit else result
    if not doit:
        return Curl(vect)
    system = _system(vect)
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
