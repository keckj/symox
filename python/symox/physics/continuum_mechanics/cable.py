"""Equilibrium of polygonal loaded cables and uniformly loaded parabolic cables."""
import symox as S
from .beam import _solve


class Cable:
    def __init__(self, support_1, support_2):
        if support_1[0] == support_2[0] or support_1[1] == support_2[1]:raise ValueError('supports require distinct labels and horizontal positions')
        supports = sorted((support_1, support_2), key=lambda p:p[1])
        self._labels = [p[0] for p in supports];self.supports = {p[0]:tuple(map(S.sympify, p[1:])) for p in supports}
        self.loads = {'distributed':{}, 'point_load':{}};self.loads_position = {};self.reaction_loads = {};self.tension = {};self.length = S.Integer(0)
        self._tension_func = None;self._cable_eqn = None;self._given_length = None
    def apply_load(self, order, load):
        label = load[0]
        if label in self.loads['distributed'] or label in self.loads['point_load']:raise ValueError('duplicate load label')
        if order == -1:
            if self.loads['distributed']:raise ValueError('point and distributed loads cannot be combined')
            x, y, mag, direction = map(S.sympify, load[1:]);a, b = [self.supports[n][0] for n in self._labels]
            if (x <= a) is S.true or (x >= b) is S.true:raise ValueError('point load must lie strictly between the supports')
            if mag.is_negative is True:raise ValueError('load magnitude must be nonnegative')
            self.loads_position[label] = (x, y);self.loads['point_load'][label] = (mag, direction)
        elif order == 0:
            if self.loads_position:raise ValueError('point and distributed loads cannot be combined')
            mag = S.sympify(load[1])
            if mag.is_nonpositive is True:raise ValueError('distributed downward load must be positive')
            self.loads['distributed'][label] = mag
        else:raise ValueError('load order must be -1 or 0')
    def remove_loads(self, *labels):
        for label in labels:
            self.loads_position.pop(label, None)
            for kind in self.loads.values():kind.pop(label, None)
    def apply_length(self, length):self._given_length = S.sympify(length)
    def solve(self, *args):
        a, b = self._labels;left, right = self.supports[a], self.supports[b];x = S.Symbol('x')
        if self.loads_position:
            names = [a]+sorted(self.loads_position, key=lambda n:self.loads_position[n][0])+[b]
            positions = [left]+[self.loads_position[n] for n in names[1:-1]]+[right]
            tensions = [S.Dummy('T') for _ in range(len(names)-1)];directions = [];length = S.Integer(0)
            for p, q in zip(positions, positions[1:]):
                delta = S.Matrix([q[i]-p[i] for i in range(2)]);size = S.sqrt(delta.dot(delta));directions.append(delta/size);length += size
            equations = []
            for i, name in enumerate(names[1:-1], 1):
                mag, angle = self.loads['point_load'][name];force = S.Matrix([mag*S.cos(angle*S.pi/180), mag*S.sin(angle*S.pi/180)])
                equations.extend(-tensions[i-1]*directions[i-1]+tensions[i]*directions[i]+force)
            solution = _solve(equations, tensions)
            if any(v.is_negative is True for v in solution.values()):raise ValueError('specified geometry requires compressive cable forces')
            self.tension = {S.Symbol(str(p)+'_'+str(q)):solution[t] for p, q, t in zip(names, names[1:], tensions)}
            self.length = length
            self.reaction_loads = {S.Symbol('R_'+str(a)+'_'+axis):-solution[tensions[0]]*directions[0][i] for i, axis in enumerate('xy')}
            self.reaction_loads.update({S.Symbol('R_'+str(b)+'_'+axis):solution[tensions[-1]]*directions[-1][i] for i, axis in enumerate('xy')})
            self._tension_func = S.Piecewise(*[(solution[t], x <= positions[i+1][0]) for i, t in enumerate(tensions)])
            self._cable_eqn = S.Piecewise(*[(p[1]+(q[1]-p[1])/(q[0]-p[0])*(x-p[0]), x <= q[0]) for p, q in zip(positions, positions[1:])])
        elif self.loads['distributed']:
            if not 1 <= len(args) <= 2:raise ValueError('uniform loading requires the lowest-point x coordinate, and optionally its y coordinate')
            x0 = S.sympify(args[0]);denom = (right[0]-x0)**2-(left[0]-x0)**2
            if len(args)==2:
                y0 = S.sympify(args[1]);coefficient = (left[1]-y0)/(left[0]-x0)**2
                if S.simplify(y0+coefficient*(right[0]-x0)**2-right[1])!=0:raise ValueError('lowest point is inconsistent with the supports')
            elif denom != 0:
                coefficient = (right[1]-left[1])/denom;y0 = left[1]-coefficient*(left[0]-x0)**2
            else:raise ValueError('equal-height supports require the lowest-point y coordinate to determine sag')
            if coefficient.is_nonpositive is True:raise ValueError('cable must sag below its supports')
            w = sum(self.loads['distributed'].values());horizontal = w/(2*coefficient)
            slope = 2*coefficient*(x-x0);self._cable_eqn = y0+coefficient*(x-x0)**2
            self._tension_func = horizontal*S.sqrt(1+slope**2);self.tension = {'distributed':self._tension_func}
            self.reaction_loads = {S.Symbol('R_'+str(a)+'_x'):-horizontal,S.Symbol('R_'+str(a)+'_y'):-horizontal*slope.subs(x,left[0]),S.Symbol('R_'+str(b)+'_x'):horizontal,S.Symbol('R_'+str(b)+'_y'):horizontal*slope.subs(x,right[0])}
            u, v = 2*coefficient*(left[0]-x0),2*coefficient*(right[0]-x0)
            primitive = lambda z:(z*S.sqrt(1+z*z)+S.asinh(z))/(4*coefficient)
            self.length = primitive(v)-primitive(u)
        else:raise ValueError('apply a load before solving the cable')
        if self._given_length is not None and S.simplify(self.length-self._given_length).is_zero is False:raise ValueError('specified cable length is inconsistent with the loaded geometry')
    def tension_at(self, x):
        if self._tension_func is None:raise ValueError('solve the cable first')
        return self._tension_func.subs(S.Symbol('x'),x)
    def plot_tension(self, **kwargs):
        if self._tension_func is None:raise ValueError('solve the cable first')
        return S.plot(self._tension_func,(S.Symbol('x'),self.supports[self._labels[0]][0],self.supports[self._labels[1]][0]),**kwargs)

__all__ = ['Cable']
