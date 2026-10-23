"""Finite vector and dyadic expressions with native scalar coefficients."""
from .. import Integer, Matrix, sqrt, sympify, NumberKind
from .._kind import Kind


class VectorKind(Kind):
    _cache = {}

    def __new__(cls, element_kind=NumberKind):
        if element_kind not in cls._cache:
            item = super().__new__(cls)
            item.element_kind = element_kind
            cls._cache[element_kind] = item
        return cls._cache[element_kind]

    def __repr__(self):
        return 'VectorKind(%s)' % self.element_kind


class _BasisDependent:
    def _new(self, components):
        cls = Vector if isinstance(self, Vector) else Dyadic
        result = object.__new__(cls)
        result._components = {b:sympify(c) for b,c in components.items() if c != 0}
        return result

    @property
    def components(self):
        return dict(self._components)

    @property
    def free_symbols(self):
        return set().union(*(c.free_symbols for c in self._components.values())) if self._components else set()

    @property
    def args(self):
        if isinstance(self,(BaseVector,BaseDyadic)):
            return self._args
        if hasattr(self,'_factors'):
            return self._factors
        mul = VectorMul if isinstance(self,Vector) else DyadicMul
        return tuple(b if c == 1 else mul(c,b) for b,c in self._components.items())

    def __bool__(self):
        return bool(self._components)

    def __eq__(self, other):
        if isinstance(other, type(0)) and other == 0:
            return not self._components
        if not isinstance(other, Vector if isinstance(self,Vector) else Dyadic):
            return False
        return self._components == other._components

    __hash__ = None

    def __add__(self, other):
        if isinstance(other, (int, Integer)) and other == 0:
            return self
        if not isinstance(other, Vector if isinstance(self,Vector) else Dyadic):
            raise TypeError('addition requires the same vector or dyadic type')
        components = self.components
        for b,c in other._components.items():
            components[b] = components.get(b,Integer(0))+c
        return self._new(components)

    __radd__ = __add__

    def __neg__(self):
        return self._new({b:-c for b,c in self._components.items()})

    def __sub__(self, other):
        return self+-other

    def __rsub__(self, other):
        return -self+other

    def __mul__(self, other):
        if isinstance(other,_BasisDependent):
            raise TypeError('use dot, cross or outer for products of vectors and dyadics')
        scalar = sympify(other)
        return self._new({b:c*scalar for b,c in self._components.items()})

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = sympify(scalar)
        if scalar == 0:
            raise ZeroDivisionError('division by zero')
        return self*(1/scalar)

    def __repr__(self):
        # Reference str printing: a unit coefficient is omitted and a
        # coefficient is parenthesized only when its text has '(', '-' or
        # '+'. Vector terms are grouped by system name, then base order;
        # both zero objects print as 0.
        def term(c,label):
            if getattr(c,'is_Integer',False) and c == 1:
                return label
            if getattr(c,'is_Float',False):
                from .. import sstr
                measure = sstr(c,full_prec=False)
            else:
                measure = str(c)
            if '(' in measure or '-' in measure or '+' in measure:
                measure = '('+measure+')'
            return measure+'*'+label
        if isinstance(self,Vector):
            key = lambda item:(str(item[0].system._name),item[0].index)
        else:
            key = lambda item:item[0]._label
        return ' + '.join(term(c,b._label) for b,c in sorted(self._components.items(),key=key)) or '0'

    __str__ = __repr__

    def _map(self, method, *args, **kwargs):
        return self._new({b:getattr(c,method)(*args,**kwargs) for b,c in self._components.items()})

    def subs(self, *args, **kwargs):
        return self._map('subs',*args,**kwargs)

    def xreplace(self, rule):
        return self._map('xreplace',rule)

    def diff(self, *args, **kwargs):
        return self._map('diff',*args,**kwargs)

    def integrate(self, *args, **kwargs):
        # SymPy's integrate() builds Integral(self), whose _eval_Integral holds
        # every component, then re-evaluates only the Integral arguments of
        # the result: a single-term field is integrated, a sum stays held.
        from .. import integrate, Integral
        flags = ('meijerg', 'conds', 'risch', 'heurisch', 'manual')
        if len(self._components) == 1:
            return self._new({b:integrate(c,*args,**kwargs) for b,c in self._components.items()})
        options = {k:v for k,v in kwargs.items() if k not in flags}
        return self._new({b:Integral(c,*args,**options) for b,c in self._components.items()})

    def _eval_Integral(self, *symbols, **assumptions):
        from .. import Integral
        return self._new({b:Integral(c,*symbols,**assumptions) for b,c in self._components.items()})

    def simplify(self, **kwargs):
        return self._map('simplify',**kwargs)

    def trigsimp(self, **kwargs):
        from .. import trigsimp
        return self._new({b:trigsimp(c,**kwargs) for b,c in self._components.items()})

    def expand(self, **kwargs):
        return self._map('expand',**kwargs)

    def evalf(self, *args, **kwargs):
        return self._map('evalf',*args,**kwargs)

    def doit(self, **kwargs):
        return self._map('doit',**kwargs)


class Vector(_BasisDependent):
    kind = VectorKind()
    is_Vector = True
    # Add(v1, v2) builds a VectorAdd, as the reference constructor
    # postprocessors do.
    _symox_vector_sum = True

    def __init__(self, *args):
        if not args:
            self._components = {}
        elif len(args) == 1 and isinstance(args[0],Vector):
            self._components = args[0].components
        else:
            raise TypeError('construct vectors from coordinate base vectors')

    def to_matrix(self, system):
        result = Matrix([0,0,0])
        for b,c in self._components.items():
            if b.system == system:
                result[b.index,0] += c
            else:
                rotation = system.rotation_matrix(b.system)
                for i in range(3):
                    result[i,0] += rotation[i,b.index]*c
        return result

    @property
    def _projections(self):
        if not self:
            return (Integer(0),)*3
        return tuple(self.to_matrix(next(iter(self._components)).system))

    def separate(self):
        result = {}
        for b,c in self._components.items():
            result[b.system] = result.get(b.system,Vector.zero)+c*b
        return result

    def dot(self, other):
        from .deloperator import Del
        from .functions import directional_derivative
        if isinstance(other,Del):
            return lambda field:directional_derivative(field,self)
        if isinstance(other,Dyadic):
            return sum((c*self.dot(b.vector1)*b.vector2 for b,c in other._components.items()),Vector.zero)
        if not isinstance(other,Vector):
            raise TypeError('dot requires a vector or dyadic')
        if not self or not other:
            return Integer(0)
        system = next(iter(self._components)).system
        a,b = self.to_matrix(system), other.to_matrix(system)
        return sum((a[i]*b[i] for i in range(3)),Integer(0))

    __and__ = dot

    def cross(self, other):
        if isinstance(other,Dyadic):
            return sum((c*self.cross(b.vector1).outer(b.vector2) for b,c in other._components.items()),Dyadic.zero)
        if not isinstance(other,Vector):
            raise TypeError('cross requires a vector or dyadic')
        if not self or not other:
            return Vector.zero
        system = next(iter(self._components)).system
        a,b = self.to_matrix(system), other.to_matrix(system)
        return sum(((a[(i+1)%3]*b[(i+2)%3]-a[(i+2)%3]*b[(i+1)%3])*system.base_vectors()[i] for i in range(3)),Vector.zero)

    __xor__ = cross

    def outer(self, other):
        if not isinstance(other,Vector):
            raise TypeError('outer product requires a vector')
        return sum((a*b*BaseDyadic(u,v) for u,a in self._components.items() for v,b in other._components.items()),Dyadic.zero)

    __or__ = outer

    def magnitude(self):
        return sqrt(self.dot(self))

    def normalize(self):
        if not self:
            raise ValueError('cannot normalize the zero vector')
        return self/self.magnitude()

    def projection(self, other, scalar=False):
        # Reference semantics: the scalar form is the coefficient
        # dot(other)/dot(self) of the vector projection, and projecting onto
        # the zero vector gives zero.
        if not isinstance(other,Vector):
            raise TypeError('projection requires a vector')
        if not self:
            return Integer(0) if scalar else Vector.zero
        coefficient = self.dot(other)/self.dot(self)
        return coefficient if scalar else coefficient*self

    def equals(self, other):
        # Reference: compare through the squared magnitude of the difference,
        # which expresses components from different (rotated) frames.
        if not isinstance(other,Vector):
            return False
        diff = self-other
        return diff.dot(diff).equals(0)


class BaseVector(Vector):
    def __init__(self, index, system, pretty_str=None, latex_str=None):
        if index not in (0,1,2):
            raise ValueError('base vector index must be 0, 1 or 2')
        self.index, self.system = index, system
        self._label = '%s.%s' % (system._name,system._vector_names[index])
        self._args = (Integer(index),system)
        self._components = {self:Integer(1)}

    def __hash__(self):
        return hash((BaseVector,self.index,self.system))

    def __eq__(self, other):
        if isinstance(other,BaseVector):
            return self.index == other.index and self.system == other.system
        return super().__eq__(other)

    def __repr__(self):
        return self._label

    __str__ = __repr__


class VectorAdd(Vector):
    def __init__(self, *args, **options):
        options.pop('evaluate',None)
        if options:
            raise TypeError('unknown vector addition options')
        self._components = sum(args,Vector.zero).components


class VectorMul(Vector):
    def __init__(self, *args, **options):
        options.pop('evaluate',None)
        if options:
            raise TypeError('unknown vector multiplication options')
        scalar, vector = Integer(1),None
        for arg in args:
            if isinstance(arg,Vector):
                if vector is not None:
                    raise TypeError('VectorMul accepts exactly one vector factor')
                vector = arg
            else:
                scalar *= sympify(arg)
        if vector is None:
            raise TypeError('VectorMul requires a vector factor')
        self._components = (vector*scalar).components
        self._factors = (scalar,vector)


class VectorZero(Vector):
    def __init__(self):
        super().__init__()


class Dyadic(_BasisDependent):
    is_Dyadic = True

    def __init__(self, *args):
        if not args:
            self._components = {}
        elif len(args) == 1 and isinstance(args[0],Dyadic):
            self._components = args[0].components
        else:
            raise TypeError('construct dyadics with vector.outer(vector)')

    def dot(self, other):
        if isinstance(other,Vector):
            return sum((c*b.vector2.dot(other)*b.vector1 for b,c in self._components.items()),Vector.zero)
        if isinstance(other,Dyadic):
            return sum((a*b*u.vector2.dot(v.vector1)*u.vector1.outer(v.vector2)
                        for u,a in self._components.items() for v,b in other._components.items()),Dyadic.zero)
        raise TypeError('dot requires a vector or dyadic')

    __and__ = dot

    def cross(self, other):
        if not isinstance(other,Vector):
            raise TypeError('dyadic cross product requires a vector')
        return sum((c*b.vector1.outer(b.vector2.cross(other)) for b,c in self._components.items()),Dyadic.zero)

    __xor__ = cross

    def to_matrix(self, system, second_system=None):
        second_system = system if second_system is None else second_system
        result = Matrix.zeros(3,3)
        for b,c in self._components.items():
            result += c*b.vector1.to_matrix(system)*b.vector2.to_matrix(second_system).T
        return result

    @property
    def T(self):
        return sum((c*b.vector2.outer(b.vector1) for b,c in self._components.items()),Dyadic.zero)

    def transpose(self):
        return self.T


class BaseDyadic(Dyadic):
    def __init__(self, vector1, vector2):
        if not isinstance(vector1,BaseVector) or not isinstance(vector2,BaseVector):
            raise TypeError('BaseDyadic requires two base vectors')
        self.vector1, self.vector2 = vector1,vector2
        self._label = '(%s|%s)' % (vector1,vector2)
        self._args = (vector1,vector2)
        self._components = {self:Integer(1)}

    def __hash__(self):
        return hash((BaseDyadic,self.vector1,self.vector2))

    def __eq__(self, other):
        if isinstance(other,BaseDyadic):
            return self.vector1 == other.vector1 and self.vector2 == other.vector2
        return super().__eq__(other)


class DyadicAdd(Dyadic):
    def __init__(self,*args,**options):
        if set(options)-{'evaluate'}:
            raise TypeError('unknown dyadic addition options')
        self._components = sum(args,Dyadic.zero).components


class DyadicMul(Dyadic):
    def __init__(self,*args,**options):
        if set(options)-{'evaluate'}:
            raise TypeError('unknown dyadic multiplication options')
        factors = [a for a in args if isinstance(a,Dyadic)]
        if len(factors) != 1:
            raise TypeError('DyadicMul requires exactly one dyadic')
        result = factors[0]
        for arg in args:
            if not isinstance(arg,Dyadic):
                result = result*arg
        self._components = result.components
        self._factors = tuple(args)


class DyadicZero(Dyadic):
    def __init__(self):
        super().__init__()


Vector.zero = VectorZero()
Dyadic.zero = DyadicZero()


def cross(vect1,vect2):
    return vect1.cross(vect2)


def dot(vect1,vect2):
    return vect1.dot(vect2)


class _Held:
    def __init__(self,*args):
        self.args = args

    def __repr__(self):
        return '%s(%s)' % (type(self).__name__,', '.join(map(str,self.args)))

    def __eq__(self,other):
        return type(self) is type(other) and self.args == other.args

    def subs(self,*args,**kwargs):
        return type(self)(*(a.subs(*args,**kwargs) for a in self.args))


class Cross(_Held):
    def __init__(self,expr1,expr2):
        if not isinstance(expr1,Vector) or not isinstance(expr2,Vector):
            raise TypeError('Cross requires vectors')
        super().__init__(expr1,expr2)

    def doit(self,**kwargs):
        return cross(*self.args)


class Dot(Cross):
    def doit(self,**kwargs):
        return dot(*self.args)
