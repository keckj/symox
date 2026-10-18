"""Manifolds, local charts and coordinate identities over native scalars."""
from collections import deque
import operator

from .. import Symbol,Dummy,Integer,Matrix,Tuple,sympify,solve,simplify

_charts = {}
_fields = {}
_coordinates = {}
_atoms = {}


def _native(value):
    return value._atom if isinstance(value,_ScalarAtom) else sympify(value)


def _sequence(value):
    return tuple(value.args) if getattr(value,'is_Tuple',False) else tuple(value)


class _ScalarAtom:
    def _symox_(self):
        return self._atom

    as_expr = _symox_

    def __repr__(self):
        return self._name

    __str__ = __repr__

    def __hash__(self):
        return hash(self._atom)

    def __eq__(self,other):
        try:
            return self._atom == _native(other)
        except (TypeError,ValueError):
            return False

    def __getattr__(self,name):
        return getattr(self._atom,name)

    def __add__(self,other):
        return self._atom+other

    def __radd__(self,other):
        return other+self._atom

    def __sub__(self,other):
        return self._atom-other

    def __rsub__(self,other):
        return other-self._atom

    def __mul__(self,other):
        return self._atom*other

    def __rmul__(self,other):
        return other*self._atom

    def __truediv__(self,other):
        return self._atom/other

    def __rtruediv__(self,other):
        return other/self._atom

    def __pow__(self,other):
        return self._atom**other

    def __rpow__(self,other):
        return other**self._atom

    def __neg__(self):
        return -self._atom

    def _latex(self,printer):
        return printer.doprint(Symbol(self._name))

    patch = property(lambda self:self.coord_sys.patch)
    manifold = property(lambda self:self.coord_sys.manifold)
    dim = property(lambda self:self.coord_sys.dim)


class Manifold:
    def __init__(self,name,dim,**kwargs):
        if kwargs:
            raise TypeError('unsupported manifold options')
        dim = operator.index(dim)
        if not isinstance(name,str) or not name or dim < 1:
            raise ValueError('a manifold requires a name and positive dimension')
        self.name,self.dim = name,dim
        self.args = (name,Integer(dim))
        self.patches = []

    def __hash__(self):
        return hash((Manifold,self.args))

    def __eq__(self,other):
        return isinstance(other,Manifold) and self.args == other.args

    def __repr__(self):
        return self.name


class Patch:
    def __init__(self,name,manifold,**kwargs):
        if kwargs or not isinstance(manifold,Manifold):
            raise TypeError('a patch requires a Manifold')
        if not isinstance(name,str) or not name:
            raise ValueError('a patch requires a nonempty name')
        self.name,self.manifold = name,manifold
        self.dim = manifold.dim
        self.args = (name,manifold)
        self.coord_systems = []
        manifold.patches.append(self)

    def __hash__(self):
        return hash((Patch,self.args))

    def __eq__(self,other):
        return isinstance(other,Patch) and self.args == other.args

    def __repr__(self):
        return self.name


class CoordinateSymbol(_ScalarAtom):
    def __new__(cls,coord_sys,index,**assumptions):
        index = operator.index(index)
        if not 0 <= index < coord_sys.dim:
            raise ValueError('coordinate index is out of range')
        key = (coord_sys,index)
        if key not in _coordinates:
            obj = super().__new__(cls)
            obj.coord_sys,obj.index = coord_sys,index
            source = coord_sys._input_symbols[index]
            obj._name = source.name
            obj._atom = Dummy(source.name,**(assumptions or source.assumptions0))
            _coordinates[key] = obj
            _atoms[obj._atom] = obj
        return _coordinates[key]

    def __init__(self,*args,**kwargs):
        pass


class BaseScalarField(_ScalarAtom):
    def __new__(cls,coord_sys,index,**kwargs):
        if kwargs:
            raise TypeError('unsupported scalar field options')
        index = operator.index(index)
        if not 0 <= index < coord_sys.dim:
            raise ValueError('scalar field index is out of range')
        key = (coord_sys,index)
        if key not in _fields:
            obj = super().__new__(cls)
            obj.coord_sys,obj.index = coord_sys,index
            source = coord_sys._input_symbols[index]
            obj._name = source.name
            obj._atom = Dummy(source.name,**source.assumptions0)
            _fields[key] = obj
            _atoms[obj._atom] = obj
        return _fields[key]

    def __init__(self,*args,**kwargs):
        pass

    def __call__(self,point):
        if not isinstance(point,Point):
            raise TypeError('scalar fields are evaluated at differential-geometric points')
        return point.coords(self.coord_sys)[self.index]

    rcall = __call__


class _CoordinateMap:
    def __init__(self,variables,expression):
        self.signature = self.variables = tuple(_native(v) for v in variables)
        self.expr = Matrix(expression)

    def __call__(self,*args):
        if len(args) != len(self.variables):
            raise ValueError('wrong coordinate dimension')
        return self.expr.subs(dict(zip(self.variables,map(_native,args))),simultaneous=True)


class CoordSystem:
    def __init__(self,name,patch,symbols=None,relations=None,**kwargs):
        if not isinstance(patch,Patch):
            raise TypeError('coordinate systems require a Patch')
        if symbols is None:
            symbols = kwargs.pop('names',None) or [Symbol('%s_%s' % (name,i),real=True) for i in range(patch.dim)]
        if kwargs:
            raise TypeError('unsupported coordinate-system options')
        self.name,self.patch,self.manifold,self.dim = str(name),patch,patch.manifold,patch.dim
        self._input_symbols = tuple(Symbol(s,real=True) if isinstance(s,str) else _native(s) for s in _sequence(symbols))
        if len(self._input_symbols) != self.dim or len(set(self._input_symbols)) != self.dim or any(not s.is_Symbol for s in self._input_symbols):
            raise ValueError('supply one distinct symbol per coordinate')
        self._relations = {}
        for (source,target),relation in (relations or {}).items():
            if hasattr(relation,'signature'):
                variables,expressions = relation.signature,relation.expr
            else:
                variables,expressions = relation
            if getattr(expressions,'is_Tuple',False):
                expressions = expressions.args
            variables,expressions = tuple(map(_native,_sequence(variables))),tuple(map(_native,_sequence(expressions)))
            if len(variables) != self.dim or len(expressions) != self.dim:
                raise ValueError('coordinate relation has the wrong dimension')
            self._relations[str(source),str(target)] = (variables,expressions)
        self._key = (self.name,patch,self._input_symbols,tuple(sorted(self._relations.items())))
        self.symbols = tuple(CoordinateSymbol(self,i) for i in range(self.dim))
        self.transforms = {}
        _charts[patch,self.name] = self
        patch.coord_systems.append(self)

    def __hash__(self):
        return hash(self._key)

    def __eq__(self,other):
        return isinstance(other,CoordSystem) and self._key == other._key

    def __repr__(self):
        return self.name

    relations = property(lambda self:dict(self._relations))

    def _path(self,target):
        if not isinstance(target,CoordSystem) or target.manifold != self.manifold:
            raise ValueError('coordinate systems must belong to the same manifold')
        if target == self:
            return tuple(_native(s) for s in self.symbols)
        edges = {}
        for chart in list(_charts.values()):
            if chart.manifold == self.manifold:
                for key,value in chart._relations.items():
                    if key in edges and edges[key] != value:
                        raise ValueError('conflicting coordinate transformation relations')
                    edges[key] = value
                edges.update(chart.transforms)
        def search():
            queue = deque([(self.name,tuple(_native(s) for s in self.symbols))])
            visited = {self.name}
            while queue:
                name,current = queue.popleft()
                for (source,destination),(variables,expressions) in edges.items():
                    if source != name or destination in visited:
                        continue
                    new = tuple(simplify(e.subs(dict(zip(variables,current)),simultaneous=True)) for e in expressions)
                    if destination == target.name:
                        return new
                    visited.add(destination)
                    queue.append((destination,new))
            return None
        result = search()
        if result is None:
            for (source,destination),(variables,expressions) in list(edges.items()):
                if (destination,source) in edges:
                    continue
                output = tuple(Dummy('inverse_%s' % i) for i in range(self.dim))
                try:
                    inverses = solve([e-y for e,y in zip(expressions,output)],variables,dict=True)
                except (NotImplementedError,ValueError):
                    continue
                if len(inverses) == 1 and all(v in inverses[0] for v in variables):
                    edges[destination,source] = (output,tuple(inverses[0][v] for v in variables))
            result = search()
        if result is None:
            raise NotImplementedError('no unambiguous coordinate transformation path; provide explicit inverse relations')
        return result

    def transformation(self,sys):
        return _CoordinateMap(self.symbols,self._path(sys))

    def transform(self,sys,coordinates=None):
        coordinates = self.symbols if coordinates is None else _sequence(coordinates)
        return self.transformation(sys)(*coordinates)

    coord_tuple_transform_to = transform

    def jacobian(self,sys,coordinates=None):
        result = Matrix(self._path(sys)).jacobian(tuple(_native(s) for s in self.symbols))
        if coordinates is not None:
            if len(coordinates) != self.dim:
                raise ValueError('wrong coordinate dimension')
            result = result.subs(dict(zip(map(_native,self.symbols),map(_native,coordinates))),simultaneous=True)
        return result

    def jacobian_determinant(self,sys,coordinates=None):
        return simplify(self.jacobian(sys,coordinates).det())

    def connect_to(self,to_sys,from_coords,to_exprs,inverse=True,fill_in_gaps=False):
        if fill_in_gaps:
            raise NotImplementedError('transformations are composed on demand; explicit gap filling is not supported')
        if to_sys.manifold != self.manifold or len(from_coords) != self.dim or len(to_exprs) != self.dim:
            raise ValueError('incompatible coordinate systems or dimensions')
        self.transforms[self.name,to_sys.name] = (tuple(map(_native,from_coords)),tuple(map(_native,to_exprs)))
        if inverse:
            to_sys._path(self)

    def point(self,coords):
        return Point(self,coords)

    def point_to_coords(self,point):
        return point.coords(self)

    def base_scalar(self,index):
        return BaseScalarField(self,index)

    def base_scalars(self):
        return tuple(self.base_scalar(i) for i in range(self.dim))

    coord_functions = base_scalars

    def base_vector(self,index):
        from ._tensors import BaseVectorField
        return BaseVectorField(self,index)

    def base_vectors(self):
        return tuple(self.base_vector(i) for i in range(self.dim))

    def base_oneform(self,index):
        from ._tensors import Differential
        return Differential(self.base_scalar(index))

    def base_oneforms(self):
        return tuple(self.base_oneform(i) for i in range(self.dim))


class Point:
    def __init__(self,coord_sys,coords,**kwargs):
        if kwargs or not isinstance(coord_sys,CoordSystem):
            raise TypeError('a point requires a coordinate system')
        self.coord_sys,self._coords = coord_sys,tuple(map(_native,_sequence(coords)))
        if len(self._coords) != coord_sys.dim:
            raise ValueError('wrong point dimension')
        self.patch,self.manifold,self.dim = coord_sys.patch,coord_sys.manifold,coord_sys.dim
        self.args = (coord_sys,self._coords)

    def coords(self,sys=None):
        return Matrix(self._coords) if sys is None else self.coord_sys.transform(sys,self._coords)

    free_symbols = property(lambda self:set().union(*(c.free_symbols for c in self._coords)))

    def __eq__(self,other):
        return isinstance(other,Point) and self.args == other.args

    def __hash__(self):
        return hash((Point,self.args))

    def __repr__(self):
        return 'Point(%s, %s)' % (self.coord_sys,self._coords)


def _scalar_at(expression,point):
    expression = _native(expression)
    replacements = {}
    for atom in expression.free_symbols:
        field = _atoms.get(atom)
        if isinstance(field,BaseScalarField):
            replacements[atom] = field(point)
    return expression.subs(replacements,simultaneous=True)


def _scalar_in_chart(expression,chart):
    expression = _native(expression)
    replacements = {}
    for atom in expression.free_symbols:
        field = _atoms.get(atom)
        if isinstance(field,BaseScalarField) and field.coord_sys != chart:
            values = chart.transform(field.coord_sys,chart.base_scalars())
            replacements[atom] = values[field.index]
    return expression.subs(replacements,simultaneous=True)


def _field_charts(expression):
    return {field.coord_sys for atom in _native(expression).free_symbols
            if isinstance((field := _atoms.get(atom)),BaseScalarField)}
