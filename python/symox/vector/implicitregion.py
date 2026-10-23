"""Polynomial implicit regions and rational parametrizations of quadrics."""
from itertools import product
from .. import (sympify, Symbol, Poly, solve, simplify, diff, FiniteSet,
                EmptySet, Tuple, sqrt, sin, cos, pi, Equality)


from ._implicitregion_ref import _ReferenceImplicitRegion, conic_coeff  # noqa: F401


class ImplicitRegion(_ReferenceImplicitRegion):
    # regular_point, singular_points, multiplicity and
    # rational_parametrization follow the reference algorithms.
    def __init__(self,variables,equation):
        self.variables = tuple(sympify(x) for x in variables)
        if not self.variables or len(self.variables) > 3 or any(not x.is_Symbol for x in self.variables) or len(set(self.variables)) != len(self.variables):
            raise ValueError('provide one to three distinct coordinate symbols')
        equation = sympify(equation)
        self.equation = equation.lhs-equation.rhs if isinstance(equation,Equality) else equation
        self.args = (self.variables,self.equation)

    def __repr__(self):
        return 'ImplicitRegion(%s, %s)' % (self.variables,self.equation)

    def rational_parametrization(self,parameters=('t','s'),reg_point=None):
        # The reference raises NotImplementedError for a plane in three
        # variables; solving the linear equation for its last variable is
        # an exact extension kept from the previous implementation.
        if self.degree == 1 and len(self.variables) == 3:
            if isinstance(parameters,str) or getattr(parameters,'is_Symbol',False):
                parameters = (parameters,)
            if len(parameters) < 2:
                raise ValueError('insufficient parameters')
            params = tuple(Symbol(p,real=True) if isinstance(p,str) else sympify(p) for p in parameters[:2])
            if set(params) & set(self.variables) or len(set(params)) != 2:
                raise ValueError('parameters must be distinct from region coordinates')
            variable = next(x for x in reversed(self.variables) if self.equation.diff(x) != 0)
            others = [x for x in self.variables if x != variable]
            mapping = dict(zip(others,params))
            mapping[variable] = solve(self.equation,variable)[0].subs(mapping,simultaneous=True)
            return tuple(mapping[x] for x in self.variables)
        return super().rational_parametrization(parameters,reg_point)

    @property
    def degree(self):
        from ..polys import total_degree
        return total_degree(self.equation)

    def _parametric_regions(self,parameters=None):
        from .parametricregion import ParametricRegion
        params = ('t','s') if parameters is None else (parameters,) if isinstance(parameters,str) else parameters
        params = tuple(Symbol(p,real=True) if isinstance(p,str) else sympify(p) for p in params)
        if len(self.variables) == 1:
            return [ParametricRegion((root,)) for root in solve(self.equation,self.variables[0])]
        if self.degree == 2:
            polynomial = Poly(self.equation,*self.variables)
            n = len(self.variables)
            weights = [polynomial.coeff_monomial(x**2) for x in self.variables]
            mixed = [polynomial.coeff_monomial(self.variables[i]*self.variables[j]) for i in range(n) for j in range(i)]
            if all(w != 0 for w in weights) and all(c == 0 for c in mixed):
                center = tuple(simplify(-polynomial.coeff_monomial(x)/(2*w)) for x,w in zip(self.variables,weights))
                scale = simplify(-self.equation.subs(dict(zip(self.variables,center)),simultaneous=True))
                radii = [sqrt(simplify(scale/w)) for w in weights]
                if all(r.is_real is not False for r in radii):
                    t = params[0]
                    if n == 2:
                        return [ParametricRegion((center[0]+radii[0]*cos(t),center[1]+radii[1]*sin(t)),(t,0,2*pi))]
                    if len(params) < 2:
                        raise ValueError('a surface parametrization needs two parameters')
                    s = params[1]
                    return [ParametricRegion((center[0]+radii[0]*sin(t)*cos(s),center[1]+radii[1]*sin(t)*sin(s),center[2]+radii[2]*cos(t)),(t,0,pi),(s,0,2*pi))]
        raise NotImplementedError('automatic bounded integration is implemented for axis-aligned ellipses and ellipsoids; use an explicit ParametricRegion')
