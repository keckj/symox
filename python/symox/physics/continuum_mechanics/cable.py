"""Equilibrium of polygonal loaded cables and uniformly loaded parabolic cables."""
import symox as S
from ._statics import _solve


class Cable:
    def __init__(self, support_1, support_2):
        if support_1[0] == support_2[0] or support_1[1] == support_2[1]:raise ValueError('supports require distinct labels and horizontal positions')
        supports = sorted((support_1, support_2), key=lambda p:p[1])
        self._labels = [p[0] for p in supports];self.supports = {p[0]:list(map(S.sympify, p[1:])) for p in (support_1, support_2)}
        self.loads = {'distributed':{}, 'point_load':{}};self.loads_position = {};self.reaction_loads = {S.Symbol('R_'+str(n)+'_'+a):S.Integer(0) for n in self._labels for a in 'xy'};self.tension = {};self.length = S.Integer(0)
        self._tension_func = None;self._cable_eqn = None;self._given_length = None
        self._lowest_x_global = S.Integer(0);self._lowest_y_global = S.Integer(0)
    def apply_load(self, order, load):
        label = load[0]
        if label in self.loads['distributed'] or label in self.loads['point_load']:raise ValueError('duplicate load label')
        if order == -1:
            if self.loads['distributed']:raise ValueError('point and distributed loads cannot be combined')
            x, y, mag, direction = map(S.sympify, load[1:]);a, b = [self.supports[n][0] for n in self._labels]
            if (x <= a) is S.true or (x >= b) is S.true:raise ValueError('point load must lie strictly between the supports')
            if mag.is_negative is True:raise ValueError('load magnitude must be nonnegative')
            self.loads_position[label] = [x, y];self.loads['point_load'][label] = [mag, direction]
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
    def change_support(self, label, new_support):
        """Replace the support ``label`` by ``new_support = (new_label, x, y)``.

        Point loads must stay strictly between the supports. (The reference
        compares the first character of each load label instead of the load
        position, a TypeError as soon as a point load exists.)
        """
        if label not in self.supports:raise ValueError("No support exists with the given label")
        kept = self._labels[(self._labels.index(label)+1) % 2];x1, y1 = self.supports[kept]
        x, y = S.sympify(new_support[1]), S.sympify(new_support[2])
        for px, _ in self.loads_position.values():
            if (px >= S.Max(x, x1)) is S.true or (px <= S.Min(x, x1)) is S.true:
                raise ValueError("The change in support will throw an existing load out of range")
        self.supports.pop(label);self.supports[new_support[0]] = [x, y]
        self._labels = [kept, new_support[0]] if (x1 < x) is S.true else [new_support[0], kept]
        self.reaction_loads = {S.Symbol('R_'+str(n)+'_'+a):S.Integer(0) for n in self._labels for a in 'xy'}
    def apply_length(self, length):
        # The reference compares against sqrt(dx**2 - dy**2); the support distance is sqrt(dx**2 + dy**2).
        (x1, y1), (x2, y2) = (self.supports[n] for n in self._labels)
        if length < S.sqrt((x1-x2)**2+(y1-y2)**2):raise ValueError("length should not be less than the distance between the supports")
        self._given_length = S.sympify(length);self.length = self._given_length
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
            self._lowest_x_global, self._lowest_y_global = x0, y0
            w = sum(self.loads['distributed'].values());horizontal = w/(2*coefficient)
            slope = 2*coefficient*(x-x0);self._cable_eqn = y0+coefficient*(x-x0)**2
            self._tension_func = horizontal*S.sqrt(1+slope**2);self.tension = {'distributed':self._tension_func}
            self.reaction_loads = {S.Symbol('R_'+str(a)+'_x'):-horizontal,S.Symbol('R_'+str(a)+'_y'):-horizontal*slope.subs(x,left[0]),S.Symbol('R_'+str(b)+'_x'):horizontal,S.Symbol('R_'+str(b)+'_y'):horizontal*slope.subs(x,right[0])}
            u, v = 2*coefficient*(left[0]-x0),2*coefficient*(right[0]-x0)
            primitive = lambda z:(z*S.sqrt(1+z*z)+S.asinh(z))/(4*coefficient)
            # The reference leaves an applied length untouched for distributed loads.
            if self._given_length is None:self.length = primitive(v)-primitive(u)
        else:raise ValueError('apply a load before solving the cable')
    def tension_at(self, x):
        if self._tension_func is None:raise ValueError('solve the cable first')
        return self._tension_func.subs(S.Symbol('x'),x)
    def plot_tension(self, **kwargs):
        if self._tension_func is None:raise ValueError('solve the cable first')
        return S.plot(self._tension_func,(S.Symbol('x'),self.supports[self._labels[0]][0],self.supports[self._labels[1]][0]),**kwargs)
    def draw(self):
        """Plot of the cable with its supports, shape and loads, built as the
        reference draw (support rectangles, load annotations, no axis)."""
        x = S.Symbol('x');(lx, ly), (rx, ry) = (self.supports[n] for n in self._labels)
        xy_min = S.Min(lx, self._lowest_y_global);xy_max = S.Max(rx, S.Max(ry, ly));max_diff = xy_max-xy_min
        width = 0.075*max_diff
        rectangles = [{'xy':(lx-width, ly), 'width':width, 'height':width, 'color':'brown', 'fill':False},
                      {'xy':(rx, ry), 'width':width, 'height':width, 'color':'brown', 'fill':False}]
        annotations = []
        if self.loads_position:
            ordered = sorted(self.loads_position.values(), key=lambda p:p[0]);points = [(lx, ly)]+ordered
            pieces = [((q[1]-p[1])*(x-p[0])/(q[0]-p[0])+p[1], x <= q[0]) for p, q in zip(points, points[1:])]
            # The reference closes the polygon with the last load point's abscissa.
            last = ordered[-1];pieces.append(((last[1]-ry)*(x-rx)/(last[0]-rx)+ry, x <= rx))
            curves = [S.Piecewise(*pieces)];length = max_diff*0.1
            for key, (mag, angle) in self.loads['point_load'].items():
                px, py = self.loads_position[key];c, s = S.cos(S.rad(angle)), S.sin(S.rad(angle))
                annotations.append({'text':'', 'xy':(px+length*c, py+length*s), 'xytext':(px, py),
                                    'arrowprops':{'width':1, 'headlength':3, 'headwidth':3, 'facecolor':'black'}})
                annotations.append({'text':f'{mag}N', 'xy':(px+length*1.6*c, py+length*1.6*s)})
        elif self.loads['distributed']:
            x0 = self._lowest_x_global;a, c = S.Dummy('a'), S.Dummy('c')
            solution = S.solve([a*(px-x0)**2+c-py for px, py in ((lx, ly), (rx, ry))], (a, c))
            curves = [solution[a]*(x-x0)**2+solution[c], self._lowest_y_global-max_diff*0.075]
            for i in range(1, 10):
                xi = lx+((rx-lx)/10)*i
                annotations.append({'text':'', 'xytext':(xi, curves[0].subs(x, xi)), 'xy':(xi, curves[1].subs(x, xi)),
                                    'arrowprops':{'width':1, 'headlength':3.5, 'headwidth':3.5, 'facecolor':'black'}})
            annotations.append({'text':f"{sum(self.loads['distributed'].values())} N/m",
                                'xy':((lx+rx)/2, self._lowest_y_global-max_diff*0.15)})
        else:raise ValueError("solve method not called and/or values provided for loads and supports not adequate")
        lim = (xy_min-0.5*max_diff, xy_max+0.5*max_diff)
        return S.plot(*curves, (x, lx, rx), xlim=lim, ylim=lim, rectangles=rectangles, show=False,
                      annotations=annotations, axis=False)

__all__ = ['Cable']
