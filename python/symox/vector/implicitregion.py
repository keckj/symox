"""Polynomial implicit regions and rational parametrizations of quadrics."""
from itertools import product
from .. import (sympify, Symbol, Poly, solve, simplify, diff, FiniteSet,
                EmptySet, Tuple, sqrt, sin, cos, pi, Equality)


class ImplicitRegion:
    def __init__(self,variables,equation):
        self.variables = tuple(sympify(x) for x in variables)
        if not self.variables or len(self.variables) > 3 or any(not x.is_Symbol for x in self.variables) or len(set(self.variables)) != len(self.variables):
            raise ValueError('provide one to three distinct coordinate symbols')
        equation = sympify(equation)
        self.equation = equation.lhs-equation.rhs if isinstance(equation,Equality) else equation
        self.args = (self.variables,self.equation)

    def __repr__(self):
        return 'ImplicitRegion(%s, %s)' % (self.variables,self.equation)

    @property
    def degree(self):
        return Poly(self.equation,*self.variables).total_degree()

    def singular_points(self):
        equations = [self.equation]+[self.equation.diff(x) for x in self.variables]
        if all(simplify(e) == 0 for e in equations):
            raise NotImplementedError('positive-dimensional singular loci are not represented')
        if any(e.is_number and e != 0 for e in equations):
            return EmptySet
        solutions = solve(equations,self.variables,dict=True)
        if any(any(x not in solution for x in self.variables) for solution in solutions):
            raise NotImplementedError('positive-dimensional singular loci are not represented')
        return FiniteSet(*(Tuple(*(solution[x] for x in self.variables)) for solution in solutions))

    def regular_point(self):
        # Deterministic finite search followed by an exact univariate solve.
        # Failure is explicit rather than claiming the real locus is empty.
        last = self.variables[-1]
        for values in product((0,1,-1,2,-2,3,-3),repeat=len(self.variables)-1):
            substitutions = dict(zip(self.variables[:-1],values))
            reduced = self.equation.subs(substitutions,simultaneous=True)
            roots = solve(reduced,last) if reduced != 0 else [sympify(0),sympify(1)]
            for root in roots:
                if root.is_real is False or root.free_symbols & set(self.variables):
                    continue
                candidate = values+(root,)
                mapping = dict(zip(self.variables,candidate))
                if simplify(self.equation.subs(mapping,simultaneous=True)) == 0 and any(simplify(self.equation.diff(x).subs(mapping,simultaneous=True)) != 0 for x in self.variables):
                    return tuple(sympify(v) for v in candidate)
        raise NotImplementedError('no regular point found by the finite exact search')

    def rational_parametrization(self,parameters=('t','s')):
        if self.degree > 2:
            raise NotImplementedError('rational parametrization is implemented for linear equations and quadrics')
        if isinstance(parameters,(str,)) or getattr(parameters,'is_Symbol',False):
            parameters = (parameters,)
        count = len(self.variables)-1
        if len(parameters) < count:
            raise ValueError('insufficient parameters')
        params = tuple(Symbol(p,real=True) if isinstance(p,str) else sympify(p) for p in parameters[:count])
        if set(params) & set(self.variables) or len(set(params)) != len(params):
            raise ValueError('parameters must be distinct from region coordinates')
        if self.degree == 1:
            variable = next(x for x in reversed(self.variables) if self.equation.diff(x) != 0)
            others = [x for x in self.variables if x != variable]
            mapping = dict(zip(others,params))
            values = solve(self.equation,variable)
            mapping[variable] = values[0].subs(mapping,simultaneous=True)
            return tuple(mapping[x] for x in self.variables)
        point = self.regular_point()
        u = Symbol('_line_parameter')
        # Select a direction chart where the second intersection is not
        # identically degenerate (e.g. at a vertex of a parabola).
        for fixed in range(len(self.variables)-1,-1,-1):
            direction = list(params)
            direction.insert(fixed,sympify(1))
            line = [p+u*d for p,d in zip(point,direction)]
            polynomial = Poly(self.equation.subs(dict(zip(self.variables,line)),simultaneous=True).expand(),u)
            quadratic = polynomial.nth(2)
            if quadratic == 0:
                continue
            distance = simplify(-polynomial.nth(1)/quadratic)
            return tuple(simplify(p+distance*d) for p,d in zip(point,direction))
        raise NotImplementedError('the chosen quadric has no nondegenerate line parametrization')

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
