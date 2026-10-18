"""Planar pin-jointed trusses solved by joint equilibrium."""
import symox as S


class Truss:
    def __init__(self):
        self._nodes = {};self._members = {};self._loads = {};self._supports = {};self._reaction_loads = {};self._internal_forces = {}
    nodes = property(lambda self:[(name, *xy) for name, xy in self._nodes.items()])
    members = property(lambda self:dict(self._members))
    loads = property(lambda self:{n:list(v) for n, v in self._loads.items()})
    supports = property(lambda self:dict(self._supports))
    reaction_loads = property(lambda self:dict(self._reaction_loads))
    internal_forces = property(lambda self:dict(self._internal_forces))
    def add_node(self, *args):
        for label, x, y in args:
            if label in self._nodes:raise ValueError('duplicate node label')
            xy = tuple(map(S.sympify, (x, y)))
            if xy in self._nodes.values():raise ValueError('nodes must have distinct positions')
            self._nodes[label] = xy
    def remove_node(self, *labels):
        for label in labels:
            if any(label in pair for pair in self._members.values()):raise ValueError('remove attached members first')
            del self._nodes[label];self._loads.pop(label, None);self._supports.pop(label, None)
    def add_member(self, *args):
        for label, start, end in args:
            if label in self._members or start == end or start not in self._nodes or end not in self._nodes:raise ValueError('invalid member or duplicate label')
            if any(set(pair)=={start, end} for pair in self._members.values()):raise ValueError('duplicate member')
            self._members[label] = (start, end)
    def remove_member(self, *labels):
        for label in labels:del self._members[label]
    def apply_load(self, *args):
        for node, magnitude, direction in args:
            if node not in self._nodes:raise ValueError('unknown load node')
            self._loads.setdefault(node, []).append(tuple(map(S.sympify, (magnitude, direction))))
    def remove_load(self, *args):
        for node, magnitude, direction in args:self._loads[node].remove(tuple(map(S.sympify, (magnitude, direction))))
    def apply_support(self, *args):
        for node, kind in args:
            if node not in self._nodes or kind not in ('pinned', 'roller'):raise ValueError('supports require an existing node and pinned or roller type')
            self._supports[node] = kind
    def remove_support(self, *nodes):
        for node in nodes:del self._supports[node]
    def solve(self):
        equations = {n:[S.Integer(0), S.Integer(0)] for n in self._nodes};unknowns = []
        for n, loads in self._loads.items():
            for mag, angle in loads:
                equations[n][0] += mag*S.cos(angle*S.pi/180);equations[n][1] += mag*S.sin(angle*S.pi/180)
        reactions = []
        for n, kind in self._supports.items():
            for axis in ('xy' if kind == 'pinned' else 'y'):
                symbol = S.Symbol(f'R_{n}_{axis}');unknowns.append(symbol);reactions.append(symbol);equations[n]['xy'.index(axis)] += symbol
        forces = {}
        for name, (a, b) in self._members.items():
            force = S.Symbol('_member_'+str(name));unknowns.append(force);forces[name] = force
            delta = [self._nodes[b][i]-self._nodes[a][i] for i in range(2)];length = S.sqrt(sum(d*d for d in delta))
            for i in range(2):equations[a][i] += force*delta[i]/length;equations[b][i] -= force*delta[i]/length
        solution = S.solve([v for pair in equations.values() for v in pair], unknowns, dict=True)
        if not solution or any(s not in solution[0] or solution[0][s].has(*unknowns) for s in unknowns):raise ValueError('truss is unstable, inconsistent or statically indeterminate')
        self._reaction_loads = {str(s):solution[0][s] for s in reactions}
        self._internal_forces = {n:solution[0][s] for n, s in forces.items()}

__all__ = ['Truss']
