"""Moving orthonormal frames, vectors, dyadics and point kinematics."""
from collections import deque
import symox as S
from ...vector._basis import _BasisDependent

_TIME = S.Symbol('t')
_COORDINATES = {}


def _matrix(value,rows=3,cols=1):
    value = S.Matrix(value)
    if value.shape != (rows,cols):raise ValueError('incorrect component matrix shape')
    return value


def _path(start,end,attribute):
    if start is end:return []
    queue = deque([(start,[])]);visited = {start}
    while queue:
        node,path = queue.popleft()
        for neighbor,value in getattr(node,attribute).items():
            if neighbor in visited:continue
            next_path = path+[(node,neighbor,value)]
            if neighbor is end:return next_path
            visited.add(neighbor);queue.append((neighbor,next_path))
    raise ValueError('objects have no connecting path')


def dynamicsymbols(names,level=0,**assumptions):
    values = S.symbols(names,**assumptions)
    multiple = isinstance(values,(list,tuple))
    values = values if multiple else (values,)
    result = [S.diff(S.Function(str(v))(_TIME),_TIME,level) if level else S.Function(str(v))(_TIME) for v in values]
    return result if multiple else result[0]


dynamicsymbols._t = _TIME
dynamicsymbols._str = "'"


class _CoordinateMeta(type):
    def __instancecheck__(cls,value):
        try:return value in _COORDINATES
        except TypeError:return False


class CoordinateSym(metaclass=_CoordinateMeta):
    def __new__(cls,name,frame,index):
        if index not in (0,1,2):raise ValueError('coordinate index must be 0, 1 or 2')
        value = S.Dummy(str(name),real=True)
        _COORDINATES[value] = frame,index
        return value


class Vector(_BasisDependent):
    is_Vector = True
    def __init__(self,inlist=0):
        if isinstance(inlist,Vector):self._components = dict(inlist._components);return
        self._components = {}
        if isinstance(inlist,(int,S.Expr)) and inlist == 0:return
        data = inlist.items() if isinstance(inlist,dict) else ((frame,value) for value,frame in inlist)
        for frame,value in data:
            if not isinstance(frame,ReferenceFrame):raise TypeError('vector components need a ReferenceFrame')
            value = _matrix(value)
            if any(v != 0 for v in value):self._components[frame] = value
    def _new(self,components):return Vector(components)
    @property
    def components(self):return dict(self._components)
    @property
    def args(self):return [(matrix,frame) for frame,matrix in self._components.items()]
    @property
    def free_symbols(self):return set().union(*(v.free_symbols for m in self._components.values() for v in m)) if self._components else set()
    def __repr__(self):
        terms = []
        for frame,values in self._components.items():
            terms.extend('%s*%s.%s'%(value,frame.name,frame._indices[i]) for i,value in enumerate(values) if value != 0)
        return ' + '.join(terms) or '0'
    __str__ = __repr__
    def __bool__(self):return bool(self._components)
    def __eq__(self,other):
        if isinstance(other,(int,S.Expr)) and other == 0:other = Vector(0)
        if not isinstance(other,Vector):return False
        if not self._components and not other._components:return True
        frame = next(iter(self._components or other._components))
        try:return all(S.trigsimp(v)==0 for v in self.to_matrix(frame)-other.to_matrix(frame))
        except ValueError:return False
    __hash__ = None
    def __add__(self,other):
        if isinstance(other,(int,S.Expr)) and other == 0:return self
        if not isinstance(other,Vector):raise TypeError('vector addition requires another vector')
        components = dict(self._components)
        for f,m in other._components.items():components[f] = components.get(f,S.zeros(3,1))+m
        return Vector(components)
    __radd__ = __add__
    def __neg__(self):return Vector({f:-m for f,m in self._components.items()})
    def __sub__(self,other):return self+-other
    def __rsub__(self,other):return -self+other
    def __mul__(self,other):
        if isinstance(other,(Vector,Dyadic)):raise TypeError('use dot, cross or outer for vector products')
        return Vector({f:m*S.sympify(other) for f,m in self._components.items()})
    __rmul__ = __mul__
    def __truediv__(self,other):return self*(1/S.sympify(other))
    def to_matrix(self,frame):
        if not isinstance(frame,ReferenceFrame):raise TypeError('expected a ReferenceFrame')
        result = S.zeros(3,1)
        for f,m in self._components.items():result += frame.dcm(f)*m
        return result.applyfunc(S.trigsimp)
    def express(self,frame,variables=False):
        result = Vector({frame:self.to_matrix(frame)})
        if variables:
            mapping = {}
            for f in self._components:mapping.update(f.variable_map(frame))
            result = result.subs(mapping,simultaneous=True)
        return result
    def dot(self,other):
        if isinstance(other,Dyadic):return other._left_dot(self)
        if not isinstance(other,Vector):raise TypeError('dot product requires a vector or dyadic')
        if not self or not other:return S.Integer(0)
        frame = next(iter(self._components));return S.trigsimp((self.to_matrix(frame).T*other.to_matrix(frame))[0])
    __and__ = dot
    def cross(self,other):
        if isinstance(other,Dyadic):return sum((coefficient*self.cross(a).outer(b) for coefficient,a,b in other._terms()),Dyadic(0))
        if not isinstance(other,Vector):raise TypeError('cross product requires a vector or dyadic')
        if not self or not other:return Vector(0)
        frame = next(iter(self._components));return Vector({frame:self.to_matrix(frame).cross(other.to_matrix(frame))})
    __xor__ = cross
    def outer(self,other):
        if not isinstance(other,Vector):raise TypeError('outer product requires another vector')
        return Dyadic({(a,b):left*right.T for a,left in self._components.items() for b,right in other._components.items()})
    __or__ = outer
    def magnitude(self):return S.sqrt(self.dot(self))
    def normalize(self):
        magnitude = self.magnitude()
        if magnitude == 0:raise ValueError('cannot normalize a zero vector')
        return self/magnitude
    def angle_between(self,other):return S.acos(self.dot(other)/(self.magnitude()*other.magnitude()))
    def projection(self,other,scalar=False):
        coefficient = self.dot(other)/other.dot(other)
        return coefficient if scalar else coefficient*other
    def _map(self,method,*args,**kwargs):return Vector({f:m.applyfunc(lambda v:getattr(v,method)(*args,**kwargs)) for f,m in self._components.items()})
    def subs(self,*args,**kwargs):return self._map('subs',*args,**kwargs)
    def xreplace(self,rule):return self._map('xreplace',rule)
    def simplify(self,**kwargs):return self._map('simplify',**kwargs)
    def trigsimp(self,**kwargs):return Vector({f:m.applyfunc(lambda v:S.trigsimp(v,**kwargs)) for f,m in self._components.items()})
    def expand(self,**kwargs):return self._map('expand',**kwargs)
    def evalf(self,*args,**kwargs):return self._map('evalf',*args,**kwargs)
    def doit(self,**kwargs):return self._map('doit',**kwargs)
    def diff(self,var,frame=None,var_in_dcm=True):
        if frame is not None:return Vector({frame:self.to_matrix(frame).applyfunc(lambda v:S.diff(v,var))})
        return Vector({f:m.applyfunc(lambda v:S.diff(v,var)) for f,m in self._components.items()})
    def dt(self,frame):
        result = Vector(0)
        for f,m in self._components.items():
            component = Vector({f:m});result += Vector({f:m.applyfunc(lambda v:S.diff(v,_TIME))})+f.ang_vel_in(frame).cross(component)
        return result
    def integrate(self,*args,**kwargs):return Vector({f:m.applyfunc(lambda v:S.integrate(v,*args,**kwargs)) for f,m in self._components.items()})
    def separate(self):return {f:Vector({f:m}) for f,m in self._components.items()}
    def free_dynamicsymbols(self,reference_frame):return set().union(*(v.atoms(S.Function,S.Derivative) for v in self.to_matrix(reference_frame)))


class Dyadic(_BasisDependent):
    def __init__(self,inlist=0):
        self._components = {}
        if isinstance(inlist,Dyadic):self._components = dict(inlist._components)
        elif isinstance(inlist,dict):self._components = {key:_matrix(value,3,3) for key,value in inlist.items() if any(v!=0 for v in value)}
        elif not(isinstance(inlist,(int,S.Expr)) and inlist==0):
            for coefficient,left,right in inlist:
                contribution = left.outer(right)*coefficient
                self._components = (self+contribution)._components
    def __repr__(self):return ' + '.join('%s*(%s|%s)'%(c,a,b) for c,a,b in self._terms()) or '0'
    __str__ = __repr__
    @property
    def args(self):return list(self._terms())
    @property
    def free_symbols(self):return set().union(*(v.free_symbols for m in self._components.values() for v in m)) if self._components else set()
    def _terms(self):
        for (a,b),matrix in self._components.items():
            for i in range(3):
                for j in range(3):
                    if matrix[i,j] != 0:yield matrix[i,j],a._axes[i],b._axes[j]
    def __bool__(self):return bool(self._components)
    def __add__(self,other):
        if isinstance(other,(int,S.Expr)) and other==0:return self
        if not isinstance(other,Dyadic):raise TypeError('dyadic addition requires another dyadic')
        result = dict(self._components)
        for key,m in other._components.items():result[key] = result.get(key,S.zeros(3,3))+m
        return Dyadic(result)
    __radd__ = __add__
    def __neg__(self):return self*-1
    def __sub__(self,other):return self+-other
    def __rsub__(self,other):return -self+other
    def __mul__(self,scalar):return Dyadic({key:m*S.sympify(scalar) for key,m in self._components.items()})
    __rmul__ = __mul__
    def __truediv__(self,scalar):return self*(1/S.sympify(scalar))
    def to_matrix(self,frame,second_frame=None):
        second_frame = second_frame or frame;result = S.zeros(3,3)
        for (a,b),m in self._components.items():result += frame.dcm(a)*m*b.dcm(second_frame)
        return result.applyfunc(S.trigsimp)
    def express(self,frame,frame2=None):return Dyadic({(frame,frame2 or frame):self.to_matrix(frame,frame2)})
    def __eq__(self,other):
        if isinstance(other,(int,S.Expr)) and other==0:other = Dyadic(0)
        if not isinstance(other,Dyadic):return False
        if not self._components and not other._components:return True
        f,g = next(iter(self._components or other._components))
        try:return all(S.trigsimp(v)==0 for v in self.to_matrix(f,g)-other.to_matrix(f,g))
        except ValueError:return False
    __hash__ = None
    def dot(self,other):
        if isinstance(other,Vector):return sum((c*a*b.dot(other) for c,a,b in self._terms()),Vector(0))
        if isinstance(other,Dyadic):return sum((c*d*b.dot(x)*a.outer(y) for c,a,b in self._terms() for d,x,y in other._terms()),Dyadic(0))
        raise TypeError('dyadic dot requires a vector or dyadic')
    __and__ = dot
    def _left_dot(self,vector):return sum((c*vector.dot(a)*b for c,a,b in self._terms()),Vector(0))
    def cross(self,vector):return sum((c*a.outer(b.cross(vector)) for c,a,b in self._terms()),Dyadic(0))
    __xor__ = cross
    def _map(self,method,*args,**kwargs):return Dyadic({key:m.applyfunc(lambda v:getattr(v,method)(*args,**kwargs)) for key,m in self._components.items()})
    def subs(self,*args,**kwargs):return self._map('subs',*args,**kwargs)
    def xreplace(self,rule):return self._map('xreplace',rule)
    def simplify(self,**kwargs):return self._map('simplify',**kwargs)
    def expand(self,**kwargs):return self._map('expand',**kwargs)
    def trigsimp(self,**kwargs):return Dyadic({key:m.applyfunc(lambda v:S.trigsimp(v,**kwargs)) for key,m in self._components.items()})
    def doit(self,**kwargs):return self._map('doit',**kwargs)
    def evalf(self,*args,**kwargs):return self._map('evalf',*args,**kwargs)
    def diff(self,var,frame=None):
        if frame is not None:return Dyadic({(frame,frame):self.to_matrix(frame).applyfunc(lambda v:S.diff(v,var))})
        return self._map('diff',var)
    def dt(self,frame):return sum((S.diff(c,_TIME)*a.outer(b)+c*a.dt(frame).outer(b)+c*a.outer(b.dt(frame)) for c,a,b in self._terms()),Dyadic(0))
    def integrate(self,*args,**kwargs):return Dyadic({key:m.applyfunc(lambda v:S.integrate(v,*args,**kwargs)) for key,m in self._components.items()})


class ReferenceFrame:
    def __init__(self,name,indices=None,latexs=None,variables=None):
        self.name = str(name);self._indices = tuple(indices or ('x','y','z'))
        if len(self._indices)!=3 or len(set(self._indices))!=3:raise ValueError('provide three distinct frame indices')
        self._dcm_dict = {};self._ang_vel_dict = {};self._ang_acc_dict = {}
        self.varlist = tuple(CoordinateSym(n,self,i) for i,n in enumerate(variables or tuple(self.name+'_'+v for v in ('x','y','z'))))
        self._axes = tuple(Vector({self:S.Matrix([int(i==j) for j in range(3)])}) for i in range(3))
        self.x,self.y,self.z = self._axes
        self._latex_vecs = latexs
    def __repr__(self):return self.name
    def __iter__(self):return iter(self._axes)
    def __getitem__(self,index):return self._axes[self._indices.index(index)] if isinstance(index,str) else self.varlist[index]
    def dcm(self,other):
        if not isinstance(other,ReferenceFrame):raise TypeError('expected a ReferenceFrame')
        result = S.eye(3)
        for a,b,value in _path(self,other,'_dcm_dict'):result *= value
        return result.applyfunc(S.trigsimp)
    def variable_map(self,other):
        coordinates = self.dcm(other)*S.Matrix(other.varlist)
        return dict(zip(self.varlist,coordinates))
    def orient_dcm(self,parent,dcm):
        if parent is self:raise ValueError('cannot orient a frame relative to itself')
        dcm = _matrix(dcm,3,3)
        self._dcm_dict[parent] = dcm;parent._dcm_dict[self] = dcm.T
        skew = dcm.T.diff(_TIME)*dcm
        self.set_ang_vel(parent,sum((skew[i,j]*parent._axes[k] for i,j,k in ((2,1,0),(0,2,1),(1,0,2))),Vector(0)).simplify())
    def orient_explicit(self,parent,dcm):self.orient_dcm(parent,_matrix(dcm,3,3).T)
    def orient_axis(self,parent,axis,angle):
        if not isinstance(axis,Vector):raise TypeError('rotation axis must be a vector')
        if axis.dt(parent)!=0:raise ValueError('axis must be fixed in the parent frame')
        angle = S.sympify(angle);n = axis.normalize().to_matrix(parent);c,s = S.cos(angle),S.sin(angle)
        skew = S.Matrix([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
        self.orient_dcm(parent,c*S.eye(3)+(1-c)*n*n.T-s*skew)
        self.set_ang_vel(parent,S.diff(angle,_TIME)*axis.normalize())
    @staticmethod
    def _axis_matrix(index,angle):
        c,s = S.cos(angle),S.sin(angle)
        return {1:S.Matrix([[1,0,0],[0,c,s],[0,-s,c]]),2:S.Matrix([[c,0,-s],[0,1,0],[s,0,c]]),3:S.Matrix([[c,s,0],[-s,c,0],[0,0,1]])}[index]
    def _orient_angles(self,parent,angles,order,body):
        order = str(order).upper().translate(str.maketrans('XYZ','123'))
        if len(order)!=3 or any(x not in '123' for x in order) or order[0]==order[1] or order[1]==order[2]:raise ValueError('invalid Euler rotation order')
        if len(angles)!=3:raise ValueError('three rotation angles are required')
        matrices = [self._axis_matrix(int(i),S.sympify(q)) for i,q in zip(order,angles)]
        rotation = matrices[2]*matrices[1]*matrices[0] if body else matrices[0]*matrices[1]*matrices[2]
        self.orient_dcm(parent,rotation)
    def orient_body_fixed(self,parent,angles,rotation_order):self._orient_angles(parent,angles,rotation_order,True)
    def orient_space_fixed(self,parent,angles,rotation_order):self._orient_angles(parent,angles,rotation_order,False)
    def orient_quaternion(self,parent,numbers):
        if len(numbers)!=4:raise ValueError('four quaternion components are required')
        q0,q1,q2,q3 = map(S.sympify,numbers)
        self.orient_dcm(parent,S.Matrix([[q0*q0+q1*q1-q2*q2-q3*q3,2*(q1*q2+q0*q3),2*(q1*q3-q0*q2)],
         [2*(q1*q2-q0*q3),q0*q0-q1*q1+q2*q2-q3*q3,2*(q2*q3+q0*q1)],
         [2*(q1*q3+q0*q2),2*(q2*q3-q0*q1),q0*q0-q1*q1-q2*q2+q3*q3]]))
    def orient(self,parent,rot_type,amounts,rot_order=''):
        kind = str(rot_type).lower()
        if kind == 'axis':
            first,second = amounts;axis,angle = (first,second) if isinstance(first,Vector) else (second,first)
            return self.orient_axis(parent,axis,angle)
        if kind == 'body':return self.orient_body_fixed(parent,amounts,rot_order)
        if kind == 'space':return self.orient_space_fixed(parent,amounts,rot_order)
        if kind == 'quaternion':return self.orient_quaternion(parent,amounts)
        if kind == 'dcm':return self.orient_dcm(parent,amounts)
        raise ValueError('unknown rotation type')
    def orientnew(self,newname,rot_type,amounts,rot_order='',variables=None,indices=None,latexs=None):
        child = ReferenceFrame(newname,indices,latexs,variables);child.orient(self,rot_type,amounts,rot_order);return child
    def set_ang_vel(self,other,value):
        value = _vector(value);self._ang_vel_dict[other] = value;other._ang_vel_dict[self] = -value
    def set_ang_acc(self,other,value):
        value = _vector(value);self._ang_acc_dict[other] = value;other._ang_acc_dict[self] = -value
    def ang_vel_in(self,other):return sum((value for a,b,value in _path(self,other,'_ang_vel_dict')),Vector(0))
    def ang_acc_in(self,other):return self._ang_acc_dict[other] if other in self._ang_acc_dict else self.ang_vel_in(other).dt(other)
    def partial_velocity(self,frame,*speeds):
        result = tuple(self.ang_vel_in(frame).diff(u,frame) for u in speeds)
        return result[0] if len(result)==1 else result
    def scalar_map(self,other):return self.variable_map(other)


def _vector(value):
    if isinstance(value,Vector):return value
    if isinstance(value,(int,S.Expr)) and value==0:return Vector(0)
    raise TypeError('expected a physics vector')


class Point:
    def __init__(self,name):self.name = str(name);self._pos_dict = {};self._vel_dict = {};self._acc_dict = {}
    def __repr__(self):return self.name
    def set_pos(self,other,value):
        if not isinstance(other,Point):raise TypeError('position reference must be a Point')
        value = _vector(value);self._pos_dict[other] = value;other._pos_dict[self] = -value
    def pos_from(self,other):return sum((value for a,b,value in _path(self,other,'_pos_dict')),Vector(0))
    def locatenew(self,name,value):
        child = Point(name);child.set_pos(self,value);return child
    def set_vel(self,frame,value):self._vel_dict[frame] = _vector(value)
    def set_acc(self,frame,value):self._acc_dict[frame] = _vector(value)
    def vel(self,frame):
        if frame in self._vel_dict:return self._vel_dict[frame]
        queue = deque([self]);visited = {self}
        while queue:
            point = queue.popleft()
            if point is not self and frame in point._vel_dict:return point._vel_dict[frame]+self.pos_from(point).dt(frame)
            for neighbor in point._pos_dict:
                if neighbor not in visited:visited.add(neighbor);queue.append(neighbor)
        raise ValueError('velocity has not been defined in this frame')
    def acc(self,frame):return self._acc_dict[frame] if frame in self._acc_dict else self.vel(frame).dt(frame)
    def v2pt_theory(self,otherpoint,outframe,fixedframe):
        value = otherpoint.vel(outframe)+fixedframe.ang_vel_in(outframe).cross(self.pos_from(otherpoint));self.set_vel(outframe,value);return value
    def v1pt_theory(self,otherpoint,outframe,interframe):
        value = self.vel(interframe)+otherpoint.vel(outframe)+interframe.ang_vel_in(outframe).cross(self.pos_from(otherpoint));self.set_vel(outframe,value);return value
    def a2pt_theory(self,otherpoint,outframe,fixedframe):
        r = self.pos_from(otherpoint);w = fixedframe.ang_vel_in(outframe)
        value = otherpoint.acc(outframe)+fixedframe.ang_acc_in(outframe).cross(r)+w.cross(w.cross(r));self.set_acc(outframe,value);return value
    def a1pt_theory(self,otherpoint,outframe,interframe):
        r = self.pos_from(otherpoint);w = interframe.ang_vel_in(outframe)
        value = self.acc(interframe)+otherpoint.acc(outframe)+interframe.ang_acc_in(outframe).cross(r)+w.cross(w.cross(r))+2*w.cross(self.vel(interframe));self.set_acc(outframe,value);return value
    def partial_velocity(self,frame,*speeds):
        result = tuple(self.vel(frame).diff(u,frame) for u in speeds)
        return result[0] if len(result)==1 else result
