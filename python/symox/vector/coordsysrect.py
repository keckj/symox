"""Orthogonal coordinate charts, connected frames and position vectors."""
from .. import (Symbol, Dummy, Integer, Matrix, sympify, sin, cos, sqrt,
                atan2, acos, simplify)
from ._basis import Vector, BaseVector

_scalar_systems = {}
_scalar_cache = {}


class _ScalarMeta(type):
    def __instancecheck__(cls,value):
        try:
            return value in _scalar_systems
        except TypeError:
            return False


class BaseScalar(metaclass=_ScalarMeta):
    """Factory for native scalar symbols carrying a coordinate registry entry."""
    def __new__(cls,index,system,pretty_str=None,latex_str=None):
        if index not in (0,1,2):
            raise ValueError('base scalar index must be 0, 1 or 2')
        key = (system,index)
        if key not in _scalar_cache:
            label = '%s.%s' % (system._name,system._variable_names[index])
            scalar = Symbol(label,real=True)
            if scalar in _scalar_systems and _scalar_systems[scalar] != key:
                scalar = Dummy(label,real=True)
            _scalar_cache[key] = scalar
            _scalar_systems[scalar] = key
        return _scalar_cache[key]


class CoordSys3D:
    def __init__(self,name,transformation=None,parent=None,location=None,rotation_matrix=None,
                 vector_names=None,variable_names=None):
        from .point import Point
        if not isinstance(name,str):
            raise TypeError('coordinate system name must be a string')
        if parent is not None and not isinstance(parent,CoordSys3D):
            raise TypeError('parent must be a CoordSys3D')
        if transformation is None:
            transformation = 'cartesian'
        if not isinstance(transformation,str) and not callable(transformation) and not isinstance(transformation,(tuple,list)):
            raise TypeError('transformation must be a chart name, callable or (variables, expressions)')
        if isinstance(transformation,str) and transformation not in ('cartesian','cylindrical','spherical'):
            raise ValueError('unknown coordinate chart: '+transformation)
        if parent is not None and parent._chart != 'cartesian':
            raise ValueError('the parent chart must be Cartesian')
        self._name,self._parent,self._chart = name,parent,transformation
        self._rotation = Matrix.eye(3) if rotation_matrix is None else Matrix(rotation_matrix)
        if self._rotation.shape != (3,3):
            raise ValueError('rotation matrix must be 3 by 3')
        if any(simplify((self._rotation*self._rotation.T)[i,j]-(1 if i == j else 0)) != 0 for i in range(3) for j in range(3)) or simplify(self._rotation.det()-1) != 0:
            raise ValueError('rotation matrix must be a proper orthogonal matrix')
        self._location = Vector.zero if location is None else location
        if not isinstance(self._location,Vector):
            raise TypeError('location must be a Vector')
        if any(isinstance(s,BaseScalar) for s in self._location.free_symbols):
            raise ValueError('location must be independent of coordinate scalars')
        if parent is None and self._location:
            raise ValueError('a translated system needs a parent')
        self._vector_names = tuple(vector_names or ('i','j','k'))
        default = ('r','theta','phi') if transformation == 'spherical' else ('r','theta','z') if transformation == 'cylindrical' else ('x','y','z')
        self._variable_names = tuple(variable_names or default)
        names = self._vector_names+self._variable_names
        if len(self._vector_names) != 3 or len(self._variable_names) != 3 or any(not isinstance(n,str) or not n.isidentifier() for n in names) or len(set(names)) != 6:
            raise ValueError('provide three distinct vector names and three distinct scalar names')
        for n in names:
            if hasattr(self,n):
                raise ValueError('coordinate name conflicts with a system attribute: '+n)
        self._key = (name,parent,str(transformation),tuple(self._rotation),frozenset(self._location.components.items()))
        self._base_vectors = tuple(BaseVector(i,self) for i in range(3))
        self._base_scalars = tuple(BaseScalar(i,self) for i in range(3))
        for n,v in zip(names,self._base_vectors+self._base_scalars):
            setattr(self,n,v)
        x,y,z = self._base_scalars
        if transformation == 'cartesian':
            self._local = Matrix([x,y,z])
            self._lame = (Integer(1),)*3
            self._chart_frame = Matrix.eye(3)
        elif transformation == 'cylindrical':
            self._local = Matrix([x*cos(y),x*sin(y),z])
            self._lame = (Integer(1),x,Integer(1))
            self._chart_frame = Matrix([[cos(y),-sin(y),0],[sin(y),cos(y),0],[0,0,1]])
        elif transformation == 'spherical':
            self._local = Matrix([x*sin(y)*cos(z),x*sin(y)*sin(z),x*cos(y)])
            self._lame = (Integer(1),x,x*sin(y))
            self._chart_frame = Matrix([[sin(y)*cos(z),cos(y)*cos(z),-sin(z)],
                                       [sin(y)*sin(z),cos(y)*sin(z),cos(z)],[cos(y),-sin(y),0]])
        else:
            if callable(transformation):
                values = transformation(x,y,z)
            else:
                variables,values = transformation
                substitutions = dict(zip(variables,self._base_scalars))
                values = [sympify(v).subs(substitutions,simultaneous=True) for v in values]
            values = tuple(values.args) if getattr(values,'is_Tuple',False) else tuple(values)
            if len(values) != 3:
                raise ValueError('coordinate transformations need three expressions')
            self._local = Matrix(values)
            jacobian = self._local.jacobian(self._base_scalars)
            gram = jacobian.T*jacobian
            if any(simplify(gram[i,j]) != 0 for i in range(3) for j in range(i)):
                raise ValueError('coordinate chart is not orthogonal')
            self._lame = tuple(sqrt(simplify(gram[i,i])) for i in range(3))
            if any(h == 0 for h in self._lame):
                raise ValueError('degenerate coordinate chart')
            self._chart_frame = Matrix(3,3,lambda i,j:simplify(jacobian[i,j]/self._lame[j]))
        if parent is None:
            self._root = self
            self._frame = self._rotation.T*self._chart_frame
            self._position = self._rotation.T*self._local
            self.origin = Point(name+'.origin')
        else:
            self._root = parent._root
            self._frame = parent._frame*self._rotation.T*self._chart_frame
            self._position = parent._position.subs(dict(zip(parent.base_scalars(),(0,0,0)))) + self._location.to_matrix(parent._root) + parent._frame*self._rotation.T*self._local
            self.origin = parent.origin.locate_new(name+'.origin',self._location)

    def __repr__(self):
        return self._name

    def __hash__(self):
        return hash(self._key)

    def __eq__(self,other):
        return isinstance(other,CoordSys3D) and self._key == other._key

    def base_vectors(self):
        return self._base_vectors

    def base_scalars(self):
        return self._base_scalars

    def lame_coefficients(self):
        return self._lame

    def rotation_matrix(self,other):
        if not isinstance(other,CoordSys3D):
            raise TypeError('expected a coordinate system')
        if self._root != other._root:
            raise ValueError('coordinate systems have no connecting path')
        return (self._frame.T*other._frame).applyfunc(simplify)

    def scalar_map(self,other):
        if self._root != other._root:
            raise ValueError('coordinate systems have no connecting path')
        if self == other:
            return dict(zip(self.base_scalars(),self.base_scalars()))
        # Convert the other chart to Cartesian coordinates in this system's
        # unrotated parent frame, then apply the inverse chart.
        origin = self._position.subs(dict(zip(self.base_scalars(),(0,0,0))))
        parent_frame = Matrix.eye(3) if self._parent is None else self._parent._frame
        local = self._rotation*parent_frame.T*(other._position-origin)
        a,b,c = tuple(local)
        if self._chart == 'cartesian':
            values = local
        elif self._chart == 'cylindrical':
            values = (sqrt(a*a+b*b),atan2(b,a),c)
        elif self._chart == 'spherical':
            radius = sqrt(a*a+b*b+c*c)
            values = (radius,acos(c/radius),atan2(b,a))
        else:
            raise NotImplementedError('inversion of custom coordinate charts is not implemented')
        return {x:simplify(v) for x,v in zip(self.base_scalars(),values)}

    def transformation_to_parent(self):
        values = self._rotation.T*self._local
        if self._parent is not None:
            values += self._location.to_matrix(self._parent)
        return tuple(values)

    def transformation_from_parent(self):
        if self._parent is None:
            raise ValueError('root coordinate system has no parent')
        mapping = self.scalar_map(self._parent)
        return tuple(mapping[x] for x in self.base_scalars())

    def position_wrt(self,other):
        return self.origin.position_wrt(other.origin if isinstance(other,CoordSys3D) else other)

    def locate_new(self,name,position,vector_names=None,variable_names=None):
        return CoordSys3D(name,parent=self,location=position,vector_names=vector_names,variable_names=variable_names)

    def create_new(self,name,transformation,variable_names=None,vector_names=None):
        return CoordSys3D(name,parent=self,transformation=transformation,variable_names=variable_names,vector_names=vector_names)

    def orient_new(self,name,orienters,location=None,vector_names=None,variable_names=None):
        from .orienters import AxisOrienter
        if not isinstance(orienters,(tuple,list)):
            orienters = [orienters]
        rotation = Matrix.eye(3)
        for orienter in orienters:
            rotation = orienter.rotation_matrix(self)*rotation if isinstance(orienter,AxisOrienter) else orienter.rotation_matrix()*rotation
        return CoordSys3D(name,parent=self,rotation_matrix=rotation,location=location,vector_names=vector_names,variable_names=variable_names)

    def orient_new_axis(self,name,angle,axis,location=None,vector_names=None,variable_names=None):
        from .orienters import AxisOrienter
        return self.orient_new(name,AxisOrienter(angle,axis),location,vector_names,variable_names)

    def orient_new_body(self,name,angle1,angle2,angle3,rotation_order,location=None,vector_names=None,variable_names=None):
        from .orienters import BodyOrienter
        return self.orient_new(name,BodyOrienter(angle1,angle2,angle3,rotation_order),location,vector_names,variable_names)

    def orient_new_space(self,name,angle1,angle2,angle3,rotation_order,location=None,vector_names=None,variable_names=None):
        from .orienters import SpaceOrienter
        return self.orient_new(name,SpaceOrienter(angle1,angle2,angle3,rotation_order),location,vector_names,variable_names)

    def orient_new_quaternion(self,name,q0,q1,q2,q3,location=None,vector_names=None,variable_names=None):
        from .orienters import QuaternionOrienter
        return self.orient_new(name,QuaternionOrienter(q0,q1,q2,q3),location,vector_names,variable_names)
