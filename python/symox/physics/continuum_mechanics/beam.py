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
    """Independent Euler–Bernoulli bending planes, axial extension and Saint-Venant torsion."""
    def __init__(self, length, elastic_modulus, shear_modulus, second_moment, area=S.Symbol('A'), variable=S.Symbol('x')):
        self.length, self.elastic_modulus, self.shear_modulus, self.area, self.variable = map(S.sympify, (length, elastic_modulus, shear_modulus, area, variable))
        moments = list(second_moment) if isinstance(second_moment, (list, tuple)) else [second_moment]*3
        if len(moments)!=3:raise ValueError('second_moment must be a scalar or three moments')
        self.second_moment = tuple(map(S.sympify, moments));self._beams = [Beam(length, elastic_modulus, i, area, variable, 'C'+str(k)) for k, i in enumerate(moments)]
        self._torsion = [];self._bc_slope = [];self._bc_deflection = []
    def apply_load(self, value, start, order, dir='y'):
        if dir not in 'xyz' or len(dir)!=1:raise ValueError('direction must be x, y or z')
        self._beams['xyz'.index(dir)].apply_load(value, start, order)
    def apply_moment_load(self, value, start, order, dir='y'):
        if dir == 'x':self._torsion.append((S.sympify(value), S.sympify(start), order));return
        # M_y bends z; M_z bends y.
        index = {'y':2, 'z':1}.get(dir)
        if index is None:raise ValueError('direction must be x, y or z')
        self._beams[index].apply_load(value if dir == 'z' else -S.sympify(value), start, order-1)
    bc_slope = property(lambda self:self._bc_slope, lambda self,v:setattr(self, '_bc_slope', list(v)))
    bc_deflection = property(lambda self:self._bc_deflection, lambda self,v:setattr(self, '_bc_deflection', list(v)))
    def shear_force(self):return [b.shear_force() for b in self._beams]
    def bending_moment(self):return [S.Add(*(v*_sf(self.variable, x, n+1) for v, x, n in self._torsion)), -self._beams[2].bending_moment(), self._beams[1].bending_moment()]
    def solve_for_reaction_loads(self, *reactions):
        equations = [b._at(b.shear_force(), self.length, True) for b in self._beams]
        equations += [b._at(b.bending_moment(), self.length, True) for b in self._beams[1:]]
        sol = _solve(equations, list(reactions))
        for b in self._beams:b._reaction_loads = sol
    def _kinematics(self):
        results = []
        for i, b in enumerate(self._beams):
            b.bc_deflection = [(p, values[i]) for p, values in self.bc_deflection]
            b.bc_slope = [(p, values[i]) for p, values in self.bc_slope]
            if i == 0:
                c = S.Symbol('C_axial');expr = b._integrated_load(2).subs(b._reaction_loads)/(self.elastic_modulus*self.area)+c
                equations = [b._at(expr, p)-values[0] for p, values in self.bc_deflection]
                sol = S.solve(equations, [c], dict=True) if equations else []
                results.append((S.diff(expr, self.variable), expr.subs(sol[0]) if sol else expr))
            else:results.append((b.slope(), b.deflection()))
        return results
    def slope(self):return [v[0] for v in self._kinematics()]
    def deflection(self):return [v[1] for v in self._kinematics()]
    def angular_deflection(self):
        x = self.variable;return S.integrate(self.bending_moment()[0]/(self.shear_modulus*self.second_moment[0]), x)
    def solve_slope_deflection(self):return self.slope(), self.deflection()

__all__ = ['Beam', 'Beam3D']
