"""Line, flux, surface and volume integrals over explicit parametrizations."""
from .. import Integer, Matrix, sqrt, sympify, simplify, trigsimp, integrate, pi
from ._basis import Vector
from .coordsysrect import CoordSys3D
from .operators import _system
from .parametricregion import ParametricRegion,parametric_region_list


class ParametricIntegral:
    def __new__(cls,field,parametricregion):
        if not isinstance(parametricregion,ParametricRegion):
            raise TypeError('expected a ParametricRegion')
        region = parametricregion
        if region.dimensions == 0:
            return Integer(0)
        if not isinstance(field,Vector):
            field = sympify(field)
        system = _system(field) or CoordSys3D('_integration_frame')
        from .functions import express
        field = express(field,system,variables=True)
        definitions = region.definition+(Integer(0),)*(3-len(region.definition))
        substitutions = dict(zip(system.base_scalars(),definitions))
        # The physical parametrization includes curvilinear metric factors.
        position = system._position.subs(substitutions,simultaneous=True)
        parameters = cls._bounds_case(region.parameters,region.limits)
        tangents = [position.diff(p) for p in parameters]
        if isinstance(field,Vector):
            physical = (system._frame*field.to_matrix(system)).subs(substitutions,simultaneous=True)
        else:
            physical = field.subs(substitutions,simultaneous=True)
        if len(parameters) == 1:
            tangent = tangents[0]
            integrand = sum(physical[i]*tangent[i] for i in range(3)) if isinstance(field,Vector) else physical*sqrt(trigsimp(sum(v*v for v in tangent)))
        elif len(parameters) == 2:
            normal = tangents[0].cross(tangents[1]).applyfunc(simplify)
            integrand = sum(physical[i]*normal[i] for i in range(3)) if isinstance(field,Vector) else physical*sqrt(trigsimp(sum(v*v for v in normal)))
        else:
            if isinstance(field,Vector):
                raise TypeError('volume integration requires a scalar field')
            determinant = Matrix(3,3,lambda i,j:tangents[j][i]).det()
            integrand = physical*determinant
        result = simplify(integrand)
        for p in parameters:
            result = _bounded_metric(result,p,*region.limits[p])
            result = integrate(result,(p,*region.limits[p]))
        return result

    @staticmethod
    def _bounds_case(parameters,limits):
        pending = [p for p in parameters if p in limits]
        ordered = []
        # Integrate dependent bounds first: an inner upper bound may contain an
        # outer variable, never a variable that has already been integrated.
        while pending:
            candidates = [p for p in pending if not any(p in (limits[q][0].free_symbols|limits[q][1].free_symbols) for q in pending if p != q)]
            if not candidates:
                raise ValueError('cyclic dependencies in region bounds')
            p = candidates[0]
            ordered.append(p)
            pending.remove(p)
        return tuple(ordered)


def _bounded_metric(expression,parameter,lower,upper):
    """Resolve the elementary chart metric signs on standard angular ranges."""
    from .. import Abs,sin,cos
    replacements = {}
    if (lower,upper) in ((0,pi),(pi,0)):
        replacements[Abs(sin(parameter))] = sin(parameter)
    if (lower,upper) in ((-pi/2,pi/2),(pi/2,-pi/2)):
        replacements[Abs(cos(parameter))] = cos(parameter)
    return expression.xreplace(replacements)


def vector_integrate(field,*region):
    if len(region) == 1 and isinstance(region[0],ParametricRegion):
        return ParametricIntegral(field,region[0])
    if len(region) == 1:
        from .implicitregion import ImplicitRegion
        from ..geometry.entity import GeometryEntity
        if isinstance(region[0],(ImplicitRegion,GeometryEntity)):
            return sum((ParametricIntegral(field,r) for r in parametric_region_list(region[0])),Integer(0))
    if isinstance(field,Vector):
        return field.integrate(*region)
    return integrate(field,*region)
