"""Coordinate conversion, vector calculus identities and scalar potentials."""
from .. import Integer, Matrix, sympify, simplify, integrate
from ._basis import Vector, Dyadic
from .operators import _get_coord_systems, _system, gradient, divergence, curl


def express(expr,system,system2=None,variables=False):
    from .coordsysrect import CoordSys3D
    if not isinstance(system,CoordSys3D) or (system2 is not None and not isinstance(system2,CoordSys3D)):
        raise TypeError('expected coordinate systems')
    if isinstance(expr,Vector):
        if system2 is not None:
            raise ValueError('system2 is only meaningful for dyadics')
        result = matrix_to_vector(expr.to_matrix(system),system)
    elif isinstance(expr,Dyadic):
        second = system if system2 is None else system2
        matrix = expr.to_matrix(system,second)
        result = sum((matrix[i,j]*system.base_vectors()[i].outer(second.base_vectors()[j]) for i in range(3) for j in range(3)),Dyadic.zero)
    else:
        result = sympify(expr)
    if variables:
        replacements = {}
        for source in _get_coord_systems(result):
            if source != system:
                replacements.update(source.scalar_map(system))
        result = result.subs(replacements,simultaneous=True)
    return result


def matrix_to_vector(matrix,system):
    matrix = Matrix(matrix)
    if matrix.shape not in ((3,1),(1,3)):
        raise ValueError('expected a three-component row or column matrix')
    return sum((v*b for v,b in zip(tuple(matrix),system.base_vectors())),Vector.zero)


def directional_derivative(field,direction_vector):
    if not isinstance(direction_vector,Vector):
        raise TypeError('direction must be a vector')
    if not isinstance(field,Vector):
        return direction_vector.dot(gradient(field))
    system = _system(field)
    if system is None:
        return Vector.zero
    # Differentiate physical Cartesian components, then project back. This
    # includes derivatives of the moving cylindrical/spherical basis vectors.
    frame = system._frame
    physical = frame*express(field,system,variables=True).to_matrix(system)
    direction = express(direction_vector,system,variables=True).to_matrix(system)
    variables,h = system.base_scalars(),system.lame_coefficients()
    derivative = Matrix([sum(direction[j]/h[j]*physical[i].diff(variables[j]) for j in range(3)) for i in range(3)])
    return matrix_to_vector((frame.T*derivative).applyfunc(simplify),system)


def laplacian(expr):
    if isinstance(expr,Vector):
        return (gradient(divergence(expr))-curl(curl(expr))).doit()
    if len(_get_coord_systems(sympify(expr))) == 1:
        # SymPy: Del().dot(Del()(expr)).doit(), over held derivatives.
        return divergence(gradient(expr,doit=False),doit=False).doit()
    return divergence(gradient(expr))


def is_conservative(field):
    if not isinstance(field,Vector):
        raise TypeError('expected a vector field')
    # This is the local curl test; topology is not inferred.
    return not curl(field).simplify()


def is_solenoidal(field):
    if not isinstance(field,Vector):
        raise TypeError('expected a vector field')
    return simplify(divergence(field)) == 0


def scalar_potential(field,coord_sys):
    if not isinstance(field,Vector):
        raise TypeError('expected a vector field')
    converted = express(field,coord_sys,variables=True)
    if not is_conservative(converted):
        raise ValueError('field is not conservative')
    values = converted.to_matrix(coord_sys)
    potential = Integer(0)
    for x,h,value in zip(coord_sys.base_scalars(),coord_sys.lame_coefficients(),tuple(values)):
        residual = simplify(h*value-potential.diff(x))
        potential += integrate(residual,x)
    if any(simplify(v) != 0 for v in (gradient(potential)-converted).components.values()):
        raise NotImplementedError('unable to construct a scalar potential with the native integrator')
    return simplify(potential)


def scalar_potential_difference(field,coord_sys,point1,point2):
    potential = scalar_potential(field,coord_sys) if isinstance(field,Vector) else sympify(field)
    first = point1.express_coordinates(coord_sys)
    second = point2.express_coordinates(coord_sys)
    return simplify(potential.subs(dict(zip(coord_sys.base_scalars(),second)),simultaneous=True)-potential.subs(dict(zip(coord_sys.base_scalars(),first)),simultaneous=True))
