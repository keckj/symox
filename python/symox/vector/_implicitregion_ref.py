"""Reference ImplicitRegion algorithms (SymPy 1.14 sympy/vector/implicitregion.py).

Ported verbatim: regular points of conics through the Diophantine
classification, singular points, multiplicities and the stereographic
rational parametrization.
"""
from ..core.numbers import Rational
from ..core.singleton import S
from ..core.symbol import symbols
from ..functions.elementary.complexes import sign
from ..functions.elementary.miscellaneous import sqrt
from ..polys.polytools import gcd
from ..sets.sets import Complement
from ..core import Tuple, diff, expand, Eq, Integer
from ..core.sorting import ordered
from ..core.symbol import Symbol


def _symbol(s, matching_symbol=None, **assumptions):
    # sympy.core.symbol._symbol
    if isinstance(s, str):
        if matching_symbol and matching_symbol.name == s:
            return matching_symbol
        return Symbol(s, **assumptions)
    elif isinstance(s, Symbol):
        return s
    else:
        raise ValueError('symbol must be string for symbol name or Symbol')
from ..solvers import solveset, nonlinsolve, diophantine
from ..polys import total_degree
from ..geometry import Point
from ..ntheory.factor_ import core


def _first_element(s):
    """``next(iter(s))`` with the reference ``Complement.__iter__`` order.

    The reference iterates ``Complement(A, B)`` as ``A`` with the members of
    ``B`` skipped; here ``A`` is ``Integers`` (0, 1, -1, 2, -2, ...).
    """
    if isinstance(s, Complement) and s.args[0] == S.Integers:
        excluded = s.args[1]
        if excluded.is_FiniteSet or excluded == S.EmptySet:
            for i in S.Integers:
                if i not in excluded:
                    return i
    return next(iter(s))


class _ReferenceImplicitRegion:
    def regular_point(self):
        equation = self.equation

        if len(self.variables) == 1:
            return (list(solveset(equation, self.variables[0], domain=S.Reals))[0],)
        elif len(self.variables) == 2:

            if self.degree == 2:
                coeffs = a, b, c, d, e, f = conic_coeff(self.variables, equation)

                if b**2 == 4*a*c:
                    x_reg, y_reg = self._regular_point_parabola(*coeffs)
                else:
                    x_reg, y_reg = self._regular_point_ellipse(*coeffs)
                return x_reg, y_reg

        if len(self.variables) == 3:
            x, y, z = self.variables

            for x_reg in range(-10, 10):
                for y_reg in range(-10, 10):
                    if not solveset(equation.subs({x: x_reg, y: y_reg}), self.variables[2], domain=S.Reals).is_empty:
                        point = (x_reg, y_reg, list(solveset(equation.subs({x: x_reg, y: y_reg})))[0])
                        # Divergence: the reference returns the first point
                        # on the surface even where the gradient vanishes;
                        # a singular point is skipped here.
                        if all(diff(equation, v).subs(dict(zip(self.variables, point))) == 0
                               for v in self.variables):
                            continue
                        return point

        # The reference ends with ``return list[self.singular_points()][0]``,
        # which raises TypeError (``list[...]`` is a generic alias) and would
        # give a singular point anyway; no regular point was found.
        raise NotImplementedError()

    def _regular_point_parabola(self, a, b, c, d, e, f):
            ok = (a, d) != (0, 0) and (c, e) != (0, 0) and b**2 == 4*a*c and (a, c) != (0, 0)

            if not ok:
                raise ValueError("Rational Point on the conic does not exist")

            if a != 0:
                d_dash, f_dash = (4*a*e - 2*b*d, 4*a*f - d**2)
                if d_dash != 0:
                    y_reg = -f_dash/d_dash
                    x_reg = -(d + b*y_reg)/(2*a)
                else:
                    ok = False
            elif c != 0:
                d_dash, f_dash = (4*c*d - 2*b*e, 4*c*f - e**2)
                if d_dash != 0:
                    x_reg = -f_dash/d_dash
                    y_reg = -(e + b*x_reg)/(2*c)
                else:
                    ok = False

            if ok:
                return x_reg, y_reg
            else:
                raise ValueError("Rational Point on the conic does not exist")

    def _regular_point_ellipse(self, a, b, c, d, e, f):
            D = 4*a*c - b**2
            ok = D

            if not ok:
                raise ValueError("Rational Point on the conic does not exist")

            if a == 0 and c == 0:
                K = -1
                L = 4*(d*e - b*f)
            elif c != 0:
                K = D
                L = 4*c**2*d**2 - 4*b*c*d*e + 4*a*c*e**2 + 4*b**2*c*f - 16*a*c**2*f
            else:
                K = D
                L = 4*a**2*e**2 - 4*b*a*d*e + 4*b**2*a*f

            ok = L != 0 and not(K > 0 and L < 0)
            if not ok:
                raise ValueError("Rational Point on the conic does not exist")

            K = Rational(K).limit_denominator(10**12)
            L = Rational(L).limit_denominator(10**12)

            k1, k2 = K.p, K.q
            l1, l2 = L.p, L.q
            g = gcd(k2, l2)

            a1 = (l2*k2)/g
            b1 = (k1*l2)/g
            c1 = -(l1*k2)/g
            a2 = sign(a1)*core(abs(a1), 2)
            r1 = sqrt(a1/a2)
            b2 = sign(b1)*core(abs(b1), 2)
            r2 = sqrt(b1/b2)
            c2 = sign(c1)*core(abs(c1), 2)
            r3 = sqrt(c1/c2)

            g = gcd(gcd(a2, b2), c2)
            a2 = a2/g
            b2 = b2/g
            c2 = c2/g

            g1 = gcd(a2, b2)
            a2 = a2/g1
            b2 = b2/g1
            c2 = c2*g1

            g2 = gcd(a2,c2)
            a2 = a2/g2
            b2 = b2*g2
            c2 = c2/g2

            g3 = gcd(b2, c2)
            a2 = a2*g3
            b2 = b2/g3
            c2 = c2/g3

            x, y, z = symbols("x y z")
            eq = a2*x**2 + b2*y**2 + c2*z**2

            # The reference parametric answer of diophantine(eq) (the
            # public diophantine refuses this incomplete family).
            from ._ternary_reference import ternary_normal_diophantine
            solutions = ternary_normal_diophantine(a2, b2, c2)

            if len(solutions) == 0:
                raise ValueError("Rational Point on the conic does not exist")

            flag = False
            for sol in solutions:
                syms = Tuple(*sol).free_symbols
                rep = dict.fromkeys(syms, 3)
                sol_z = sol[2]

                if sol_z == 0:
                    flag = True
                    continue

                if not isinstance(sol_z, (int, Integer)):
                    syms_z = sol_z.free_symbols

                    if len(syms_z) == 1:
                        p = next(iter(syms_z))
                        p_values = Complement(S.Integers, solveset(Eq(sol_z, 0), p, S.Integers))
                        rep[p] = _first_element(p_values)

                    if len(syms_z) == 2:
                        p, q = list(ordered(syms_z))

                        for i in S.Integers:
                            subs_sol_z = sol_z.subs(p, i)
                            q_values = Complement(S.Integers, solveset(Eq(subs_sol_z, 0), q, S.Integers))

                            if not q_values.is_empty:
                                rep[p] = i
                                rep[q] = _first_element(q_values)
                                break

                    if len(syms) != 0:
                        x, y, z = tuple(s.subs(rep) for s in sol)
                    else:
                        x, y, z =   sol
                    flag = False
                    break

            if flag:
                raise ValueError("Rational Point on the conic does not exist")

            x = (x*g3)/r1
            y = (y*g2)/r2
            z = (z*g1)/r3
            x = x/z
            y = y/z

            if a == 0 and c == 0:
                x_reg = (x + y - 2*e)/(2*b)
                y_reg = (x - y - 2*d)/(2*b)
            elif c != 0:
                x_reg = (x - 2*d*c + b*e)/K
                y_reg = (y - b*x_reg - e)/(2*c)
            else:
                y_reg = (x - 2*e*a + b*d)/K
                x_reg = (y - b*y_reg - d)/(2*a)

            return x_reg, y_reg

    def singular_points(self):
        eq_list = [self.equation]
        for var in self.variables:
            eq_list += [diff(self.equation, var)]

        return nonlinsolve(eq_list, list(self.variables))

    def multiplicity(self, point):
        if isinstance(point, Point):
            point = point.args

        modified_eq = self.equation

        for i, var in enumerate(self.variables):
            modified_eq = modified_eq.subs(var, var + point[i])
        modified_eq = expand(modified_eq)

        if len(modified_eq.args) != 0:
            terms = modified_eq.args
            m = min(total_degree(term) for term in terms)
        else:
            terms = modified_eq
            m = total_degree(terms)

        return m

    def rational_parametrization(self, parameters=('t', 's'), reg_point=None):
        equation = self.equation
        degree = self.degree

        if degree == 1:
            if len(self.variables) == 1:
                return (equation,)
            elif len(self.variables) == 2:
                x, y = self.variables
                y_par = list(solveset(equation, y))[0]
                return x, y_par
            else:
                raise NotImplementedError()

        point = ()

        # Finding the (n - 1) fold point of the monoid of degree
        if degree == 2:
            # For degree 2 curves, either a regular point or a singular point can be used.
            if reg_point is not None:
                # Using point provided by the user as regular point
                point = reg_point
            else:
                if len(self.singular_points()) != 0:
                    point = list(self.singular_points())[0]
                else:
                    point = self.regular_point()

        if len(self.singular_points()) != 0:
            singular_points = self.singular_points()
            for spoint in singular_points:
                syms = Tuple(*spoint).free_symbols
                rep = dict.fromkeys(syms, 2)

                if len(syms) != 0:
                    spoint = tuple(s.subs(rep) for s in spoint)

                if self.multiplicity(spoint) == degree - 1:
                    point = spoint
                    break

        if len(point) == 0:
            # The region in not a monoid
            raise NotImplementedError()

        modified_eq = equation

        # Shifting the region such that fold point moves to origin
        for i, var in enumerate(self.variables):
            modified_eq = modified_eq.subs(var, var + point[i])
        modified_eq = expand(modified_eq)

        hn = hn_1 = 0
        for term in modified_eq.args:
            if total_degree(term) == degree:
                hn += term
            else:
                hn_1 += term

        hn_1 = -1*hn_1

        if not isinstance(parameters, tuple):
            parameters = (parameters,)

        if len(self.variables) == 2:

            parameter1 = parameters[0]
            if parameter1 == 's':
                # To avoid name conflict between parameters
                s = _symbol('s_', real=True)
            else:
                s = _symbol('s', real=True)
            t = _symbol(parameter1, real=True)

            hn = hn.subs({self.variables[0]: s, self.variables[1]: t})
            hn_1 = hn_1.subs({self.variables[0]: s, self.variables[1]: t})

            x_par = (s*(hn_1/hn)).subs(s, 1) + point[0]
            y_par = (t*(hn_1/hn)).subs(s, 1) + point[1]

            return x_par, y_par

        elif len(self.variables) == 3:

            parameter1, parameter2 = parameters
            if 'r' in parameters:
                # To avoid name conflict between parameters
                r = _symbol('r_', real=True)
            else:
                r = _symbol('r', real=True)
            s = _symbol(parameter2, real=True)
            t = _symbol(parameter1, real=True)

            hn = hn.subs({self.variables[0]: r, self.variables[1]: s, self.variables[2]: t})
            hn_1 = hn_1.subs({self.variables[0]: r, self.variables[1]: s, self.variables[2]: t})

            x_par = (r*(hn_1/hn)).subs(r, 1) + point[0]
            y_par = (s*(hn_1/hn)).subs(r, 1) + point[1]
            z_par = (t*(hn_1/hn)).subs(r, 1) + point[2]

            return x_par, y_par, z_par

        raise NotImplementedError()


def conic_coeff(variables, equation):
    if total_degree(equation) != 2:
        raise ValueError()
    x = variables[0]
    y = variables[1]

    equation = expand(equation)
    a = equation.coeff(x**2)
    b = equation.coeff(x*y)
    c = equation.coeff(y**2)
    d = equation.coeff(x, 1).coeff(y, 0)
    e = equation.coeff(y, 1).coeff(x, 0)
    f = equation.coeff(x, 0).coeff(y, 0)
    return a, b, c, d, e, f
