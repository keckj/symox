"""Euler–Bernoulli beams using singularity functions and equilibrium equations."""
import symox as S
from math import factorial


def _sf(x, start, order):return S.SingularityFunction(x, start, order)
def _solve(equations, unknowns):
    solutions = S.solve(equations, unknowns, dict=True)
    if not solutions:raise ValueError('boundary conditions or equilibrium have no solution')
    solution = solutions[0]
    if any(v not in solution or solution[v].has(*unknowns) for v in unknowns):raise ValueError('insufficient independent boundary conditions')
    return solution
class Beam:
    def __init__(self, length, elastic_modulus, second_moment, area=S.Symbol('A'), variable=S.Symbol('x'), base_char='C'):
        self.length, self.elastic_modulus, self.second_moment, self.area, self.variable = map(S.sympify, (length, elastic_modulus, second_moment, area, variable))
        if self.length.is_nonpositive is True:raise ValueError('beam length must be positive')
        self.base_char = base_char;self._loads = [];self._reaction_loads = {};self._solution = {}
        self._bc = {'deflection':[], 'slope':[], 'bending_moment':[], 'shear_force':[]}
    @property
    def bc_deflection(self):return self._bc['deflection']
    @bc_deflection.setter
    def bc_deflection(self, value):self._bc['deflection'] = list(value)
    @property
    def bc_slope(self):return self._bc['slope']
    @bc_slope.setter
    def bc_slope(self, value):self._bc['slope'] = list(value)
    @property
    def bc_bending_moment(self):return self._bc['bending_moment']
    @bc_bending_moment.setter
    def bc_bending_moment(self, value):self._bc['bending_moment'] = list(value)
    @property
    def bc_shear_force(self):return self._bc['shear_force']
    @bc_shear_force.setter
    def bc_shear_force(self, value):self._bc['shear_force'] = list(value)
    @property
    def applied_loads(self):return list(self._loads)
    @property
    def reaction_loads(self):return dict(self._reaction_loads)
    def apply_load(self, value, start, order, end=None):
        if int(order)!=order or order < -2:raise ValueError('load order must be an integer at least -2')
        if end is not None and order < 0:raise ValueError('end applies only to distributed loads')
        self._loads.append((S.sympify(value), S.sympify(start), int(order), None if end is None else S.sympify(end)))
        self._solution = {}
    def remove_load(self, value, start, order, end=None):
        self._loads.remove((S.sympify(value), S.sympify(start), int(order), None if end is None else S.sympify(end)));self._solution = {}
    def apply_support(self, loc, type='fixed'):
        loc = S.sympify(loc)
        if type not in ('fixed', 'pin', 'roller'):raise ValueError('support type must be fixed, pin or roller')
        reaction = S.Symbol('R_'+str(loc));self.apply_load(reaction, loc, -1);self.bc_deflection.append((loc, 0))
        if type != 'fixed':return reaction
        moment = S.Symbol('M_'+str(loc));self.apply_load(moment, loc, -2);self.bc_slope.append((loc, 0))
        return reaction, moment
    def _pieces(self):
        for value, start, order, end in self._loads:
            yield value, start, order
            if end is not None:
                for k in range(order+1):yield -value*S.binomial(order, k)*(end-start)**(order-k), end, k
    def _integrated_load(self, times):
        x = self.variable;result = S.Integer(0)
        for value, start, order in self._pieces():
            factor = S.Integer(1);power = order
            for _ in range(times):
                power += 1
                if power > 0:factor /= power
            result += value*factor*_sf(x, start, power)
        return result
    @property
    def load(self):return self._integrated_load(0).subs(self._reaction_loads)
    def shear_force(self):return -self._integrated_load(1).subs(self._reaction_loads)
    def bending_moment(self):return -self._integrated_load(2).subs(self._reaction_loads)
    def _kinematics(self):
        c1, c2 = S.symbols(self.base_char+'3 '+self.base_char+'4');x = self.variable
        rigidity = self.elastic_modulus*self.second_moment
        if x in rigidity.free_symbols:
            slope = S.integrate(-self.bending_moment()/rigidity, x)+c1
            deflection = S.integrate(slope, x)+c2
        else:
            slope = self._integrated_load(3).subs(self._reaction_loads)/rigidity+c1
            deflection = self._integrated_load(4).subs(self._reaction_loads)/rigidity+c1*x+c2
        return slope, deflection, (c1, c2)
    def _at(self, expr, position, outside=False):
        x = self.variable;position = S.sympify(position);mapping = {}
        for node in S.preorder_traversal(expr):
            if getattr(node, 'class_name', None) != 'SingularityFunction':continue
            _, start, power = node.args;delta = S.simplify(position-start)
            if power < 0:val = 0
            elif delta.is_negative is True:val = 0
            elif delta.is_nonnegative is True or outside:val = delta**power
            else:continue
            mapping[node] = val
        return expr.xreplace(mapping).subs(x, position)
    def solve_for_reaction_loads(self, *reactions):
        self._reaction_loads = {};self._solution = {}
        slope, deflection, constants = self._kinematics()
        equations = [self._at(self.shear_force(), self.length, True), self._at(self.bending_moment(), self.length, True)]
        curves = {'slope':slope, 'deflection':deflection, 'shear_force':self.shear_force(), 'bending_moment':self.bending_moment()}
        equations += [self._at(curves[kind], position)-S.sympify(value) for kind, conditions in self._bc.items() for position, value in conditions]
        # Integration constants need boundary conditions only when solving kinematics.
        unknowns = list(reactions)+[c for c in constants if any(e.has(c) for e in equations)]
        self._solution = _solve(equations, unknowns)
        self._reaction_loads = {r:self._solution[r] for r in reactions}
    def _curve(self, which):
        slope, deflection, constants = self._kinematics();curve = slope if which == 'slope' else deflection
        mapping = dict(self._solution)
        if not all(c in mapping for c in constants):
            equations = [self._at(slope, x)-S.sympify(v) for x, v in self.bc_slope]+[self._at(deflection, x)-S.sympify(v) for x, v in self.bc_deflection]
            if equations:
                solutions = S.solve(equations, constants, dict=True)
                if solutions:mapping.update(solutions[0])
        return curve.subs(mapping)
    def slope(self):return self._curve('slope')
    def deflection(self):return self._curve('deflection')
    def point_cflexure(self):
        solutions = S.solve(self.bending_moment(), self.variable)
        return [v for v in solutions if (v > 0) is True and (v < self.length) is True]
    def _plot(self, expr, *args, **kwargs):return S.plot(expr, (self.variable, 0, self.length), *args, **kwargs)
    def plot_loading_results(self, *args, **kwargs):
        return [self._plot(expr, *args, **kwargs) for expr in (self.load, self.shear_force(), self.bending_moment(), self.slope(), self.deflection())]
    def plot_shear_force(self, *args, **kwargs):return self._plot(self.shear_force(), *args, **kwargs)
    def plot_bending_moment(self, *args, **kwargs):return self._plot(self.bending_moment(), *args, **kwargs)
    def plot_slope(self, *args, **kwargs):return self._plot(self.slope(), *args, **kwargs)
    def plot_deflection(self, *args, **kwargs):return self._plot(self.deflection(), *args, **kwargs)


class Beam3D:
    """Timoshenko beam in three dimensions, fixed at both ends (reference formulation).

    Loads enter the governing equations through their magnitudes (load_vector and
    moment_load_vector); slopes and deflections solve the coupled shear/bending ODEs
    with zero slope and deflection at x = 0 and x = length.
    """
    def __init__(self, length, elastic_modulus, shear_modulus, second_moment, area, variable=S.Symbol('x')):
        self.length, self.elastic_modulus, self.variable = map(S.sympify, (length, elastic_modulus, variable))
        self.shear_modulus = shear_modulus;self.second_moment = second_moment;self.area = area
        self._boundary_conditions = {'deflection':[], 'slope':[]};self._reaction_loads = {}
        self._load_vector = [0, 0, 0];self._moment_load_vector = [0, 0, 0];self._torsion_moment = {}
        self._load_Singularity = [0, 0, 0];self._slope = [0, 0, 0];self._deflection = [0, 0, 0];self._angular_deflection = 0
    shear_modulus = property(lambda self:self._shear_modulus, lambda self,e:setattr(self, '_shear_modulus', S.sympify(e)))
    area = property(lambda self:self._area, lambda self,a:setattr(self, '_area', S.sympify(a)))
    @property
    def second_moment(self):return self._second_moment
    @second_moment.setter
    def second_moment(self, i):self._second_moment = [S.sympify(v) for v in i] if isinstance(i, list) else S.sympify(i)
    load_vector = property(lambda self:self._load_vector)
    moment_load_vector = property(lambda self:self._moment_load_vector)
    boundary_conditions = property(lambda self:self._boundary_conditions)
    reaction_loads = property(lambda self:self._reaction_loads)
    bc_deflection = property(lambda self:self._boundary_conditions['deflection'], lambda self,v:self._boundary_conditions.__setitem__('deflection', list(v)))
    bc_slope = property(lambda self:self._boundary_conditions['slope'], lambda self,v:self._boundary_conditions.__setitem__('slope', list(v)))
    def polar_moment(self):
        if not isinstance(self.second_moment, (list, tuple)):return 2*self.second_moment
        return S.Add(*self.second_moment)
    def apply_load(self, value, start, order, dir="y"):
        x = self.variable;value, start, order = map(S.sympify, (value, start, order))
        k = {'x':0, 'y':1}.get(dir, 2)
        if not order == -1:self._load_vector[k] += value
        self._load_Singularity[k] += value*_sf(x, start, order)
    def apply_moment_load(self, value, start, order, dir="y"):
        x = self.variable;value, start, order = map(S.sympify, (value, start, order))
        k = {'x':0, 'y':1}.get(dir, 2)
        if not order == -2:self._moment_load_vector[k] += value
        elif k == 0:self._torsion_moment[start] = self._torsion_moment.get(start, 0)+value
        self._load_Singularity[0] += value*_sf(x, start, order)
    def apply_support(self, loc, type="fixed"):
        reaction_load = S.Symbol('R_'+str(loc))
        if type in ("pin", "roller"):
            self._reaction_loads[reaction_load] = reaction_load;self.bc_deflection.append((loc, [0, 0, 0]))
        else:
            self._reaction_loads[reaction_load] = [reaction_load, S.Symbol('M_'+str(loc))]
            self.bc_deflection.append((loc, [0, 0, 0]));self.bc_slope.append((loc, [0, 0, 0]))
    def solve_for_reaction_loads(self, *reaction):
        x, l = self.variable, self.length
        shear_curves = [S.integrate(load, x) for load in self._load_Singularity]
        moment_curves = [S.integrate(shear, x) for shear in shear_curves]
        for i in range(3):
            react = [r for r in reaction if (shear_curves[i].has(r) or moment_curves[i].has(r))]
            if len(react) == 0:continue
            shear_curve = S.limit(shear_curves[i], x, l);moment_curve = S.limit(moment_curves[i], x, l)
            sol_dict = dict(zip(react, list(S.linsolve([shear_curve, moment_curve], react).args[0])))
            for key in sol_dict:
                if key in self._reaction_loads and sol_dict[key] != self._reaction_loads[key]:
                    raise ValueError("Ambiguous solution for %s in different directions." % key)
            self._reaction_loads.update(sol_dict)
    def shear_force(self):
        x, q = self.variable, self._load_vector
        return [S.integrate(-q[0], x), S.integrate(-q[1], x), S.integrate(-q[2], x)]
    def axial_force(self):return self.shear_force()[0]
    def shear_stress(self):return [v/self._area for v in self.shear_force()]
    def axial_stress(self):return self.axial_force()/self._area
    def bending_moment(self):
        x, m = self.variable, self._moment_load_vector;shear = self.shear_force()
        return [S.integrate(-m[0], x), S.integrate(-m[1]+shear[2], x), S.integrate(-m[2]-shear[1], x)]
    def torsional_moment(self):return self.bending_moment()[0]
    def solve_for_torsion(self):
        x = self.variable;sum_moments = S.Add(*self._torsion_moment.values());points = list(self._torsion_moment)
        torque_diagram = S.Piecewise((sum_moments, x <= points[0]), (0, x >= points[0]))
        for i in range(1, len(points)):
            sum_moments -= self._torsion_moment[points[i-1]]
            torque_diagram += S.Piecewise((0, x <= points[i-1]), (sum_moments, x <= points[i]), (0, x >= points[i]))
        self._angular_deflection = S.integrate(torque_diagram, x)/(self.shear_modulus*self.polar_moment())
    def solve_slope_deflection(self):
        x, l, E, G, A = self.variable, self.length, self.elastic_modulus, self.shear_modulus, self._area
        I = self.second_moment
        I_y, I_z = (I[0], I[1]) if isinstance(I, list) else (I, I)
        load, moment = self._load_vector, self._moment_load_vector
        defl, theta = S.Function('defl'), S.Function('theta')
        C1, C2, C_i = S.Symbol('C1'), S.Symbol('C2'), S.Symbol('C_i')
        def _consts(expr, a, b):return list(S.linsolve([expr.subs(x, 0), expr.subs(x, l)], a, b).args[0])
        # Axial: (E*A*u')' + load_x = 0
        def_x = S.dsolve(S.Eq(S.Derivative(E*A*S.Derivative(defl(x), x), x)+load[0], 0), defl(x)).args[1]
        c = _consts(def_x, C1, C2);def_x = def_x.subs({C1:c[0], C2:c[1]})
        self._deflection[0] = def_x;self._slope[0] = def_x.diff(x)
        # y deflection and z slope
        eq1 = S.Derivative(E*I_z*S.Derivative(theta(x), x), x)+(S.integrate(-load[1], x)+C_i)+moment[2]
        slope_z = S.dsolve(S.Eq(eq1, 0), theta(x)).args[1]
        c = _consts(slope_z, C1, C2);slope_z = slope_z.subs({C1:c[0], C2:c[1]})
        eq2 = G*A*(S.Derivative(defl(x), x))+load[1]*x-C_i-G*A*slope_z
        def_y = S.dsolve(S.Eq(eq2, 0), defl(x)).args[1]
        c = _consts(def_y, C1, C_i)
        self._deflection[1] = def_y.subs({C1:c[0], C_i:c[1]});self._slope[2] = slope_z.subs(C_i, c[1])
        # z deflection and y slope
        eq1 = S.Derivative(E*I_y*S.Derivative(theta(x), x), x)+(S.integrate(load[2], x)-C_i)+moment[1]
        slope_y = S.dsolve(S.Eq(eq1, 0), theta(x)).args[1]
        c = _consts(slope_y, C1, C2);slope_y = slope_y.subs({C1:c[0], C2:c[1]})
        eq2 = G*A*(S.Derivative(defl(x), x))+load[2]*x-C_i+G*A*slope_y
        def_z = S.dsolve(S.Eq(eq2, 0), defl(x)).args[1]
        c = _consts(def_z, C1, C_i)
        self._deflection[2] = def_z.subs({C1:c[0], C_i:c[1]});self._slope[1] = slope_y.subs(C_i, c[1])
    def slope(self):return self._slope
    def deflection(self):return self._deflection
    def angular_deflection(self):return self._angular_deflection

__all__ = ['Beam', 'Beam3D']
