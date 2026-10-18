"""Hamilton quaternions over native symbolic scalars."""
import operator


class Quaternion:
    def __init__(self,a=0,b=0,c=0,d=0,real_field=True,norm=None):
        from .. import sympify,simplify
        self.a,self.b,self.c,self.d = map(sympify,(a,b,c,d))
        self.real_field = bool(real_field)
        self._norm = None if norm is None else sympify(norm)
        if self._norm is not None and simplify(self._norm**2-sum(v*v for v in self.args)) != 0:
            raise ValueError('the supplied norm does not match the quaternion')

    @property
    def args(self):return self.a,self.b,self.c,self.d
    @property
    def free_symbols(self):return set().union(*(v.free_symbols for v in self.args))
    def __repr__(self):return 'Quaternion(%s)' % ', '.join(map(str,self.args))
    def __hash__(self):return hash((Quaternion,self.args))
    def __eq__(self,other):return isinstance(other,Quaternion) and self.args == other.args
    def _coerce(self,value):return value if isinstance(value,Quaternion) else Quaternion(value,real_field=self.real_field)
    def __add__(self,other):
        other = self._coerce(other)
        return Quaternion(*(a+b for a,b in zip(self.args,other.args)),real_field=self.real_field)
    __radd__ = __add__
    def __neg__(self):return Quaternion(*(-v for v in self.args),real_field=self.real_field)
    def __sub__(self,other):return self+-self._coerce(other)
    def __rsub__(self,other):return self._coerce(other)+-self
    def __mul__(self,other):
        e,f,g,h = self._coerce(other).args
        a,b,c,d = self.args
        return Quaternion(a*e-b*f-c*g-d*h,a*f+b*e+c*h-d*g,a*g-b*h+c*e+d*f,a*h+b*g-c*f+d*e,real_field=self.real_field)
    def __rmul__(self,other):return self._coerce(other)*self
    def __truediv__(self,other):return self*self._coerce(other).inverse()
    def __rtruediv__(self,other):return self._coerce(other)*self.inverse()
    def __pow__(self,n):
        n = operator.index(n)
        if n < 0:return self.inverse()**(-n)
        result,base = Quaternion(1),self
        while n:
            if n&1:result *= base
            base *= base;n >>= 1
        return result
    def conjugate(self):return Quaternion(self.a,-self.b,-self.c,-self.d,real_field=self.real_field)
    def norm(self):
        from .. import sqrt,simplify
        return self._norm if self._norm is not None else sqrt(simplify(sum(v*v for v in self.args)))
    def inverse(self):
        from .. import simplify
        squared = simplify(sum(v*v for v in self.args))
        if squared == 0:raise ValueError('zero quaternion has no inverse')
        return Quaternion(*(v/squared for v in self.conjugate().args),real_field=self.real_field)
    def normalize(self):
        n = self.norm()
        if n == 0:raise ValueError('zero quaternion cannot be normalized')
        return Quaternion(*(v/n for v in self.args),real_field=self.real_field,norm=1)
    def subs(self,*args,**kwargs):return Quaternion(*(v.subs(*args,**kwargs) for v in self.args),real_field=self.real_field)
    def simplify(self,**kwargs):return Quaternion(*(v.simplify(**kwargs) for v in self.args),real_field=self.real_field)
    def expand(self,**kwargs):return Quaternion(*(v.expand(**kwargs) for v in self.args),real_field=self.real_field)
    def diff(self,*args,**kwargs):return Quaternion(*(v.diff(*args,**kwargs) for v in self.args),real_field=self.real_field)
    def integrate(self,*args,**kwargs):
        from .. import integrate
        return Quaternion(*(integrate(v,*args,**kwargs) for v in self.args),real_field=self.real_field)
    def evalf(self,*args,**kwargs):return Quaternion(*(v.evalf(*args,**kwargs) for v in self.args),real_field=self.real_field)
    def scalar_part(self):return self.a
    def vector_part(self):return Quaternion(0,self.b,self.c,self.d,real_field=self.real_field)
    def is_zero_quaternion(self):return all(v == 0 for v in self.args)
    def is_pure(self):return self.a == 0
    def exp(self):
        from .. import sqrt,exp,cos,sin
        v = sqrt(self.b**2+self.c**2+self.d**2)
        if v == 0:return Quaternion(exp(self.a))
        factor = exp(self.a)*sin(v)/v
        return Quaternion(exp(self.a)*cos(v),self.b*factor,self.c*factor,self.d*factor)
    def log(self):
        from .. import log,sqrt,acos
        v,n = sqrt(self.b**2+self.c**2+self.d**2),self.norm()
        if n == 0:raise ValueError('zero quaternion has no logarithm')
        if v == 0:
            if self.a.is_positive is True:return Quaternion(log(self.a))
            raise ValueError('a negative real quaternion requires a choice of logarithm axis')
        factor = acos(self.a/n)/v
        return Quaternion(log(n),self.b*factor,self.c*factor,self.d*factor)
    def pow(self,p):return (self.log()*p).exp()
    @classmethod
    def from_axis_angle(cls,pair):
        from .. import sympify,sqrt,sin,cos
        axis,angle = pair
        axis = tuple(map(sympify,axis));angle = sympify(angle)
        if len(axis) != 3:raise ValueError('rotation axis must have three components')
        norm = sqrt(sum(v*v for v in axis))
        if norm == 0:raise ValueError('rotation axis must be nonzero')
        return cls(cos(angle/2),*(sin(angle/2)*v/norm for v in axis))
    def to_axis_angle(self):
        from .. import acos,sqrt
        q = self.normalize();length = sqrt(q.b**2+q.c**2+q.d**2)
        return ((1,0,0),0) if length == 0 else ((q.b/length,q.c/length,q.d/length),2*acos(q.a))
    def to_rotation_matrix(self,v=None,homogeneous=True):
        from .. import Matrix,simplify,eye
        a,b,c,d = self.args;n = self.norm()**2
        if n == 0:raise ValueError('zero quaternion does not define a rotation')
        matrix = Matrix([[a*a+b*b-c*c-d*d,2*(b*c-a*d),2*(b*d+a*c)],
                         [2*(b*c+a*d),a*a-b*b+c*c-d*d,2*(c*d-a*b)],
                         [2*(b*d-a*c),2*(c*d+a*b),a*a-b*b-c*c+d*d]])/n
        matrix = matrix.applyfunc(simplify)
        if v is None:return matrix
        vector = Matrix(v)
        if vector.shape != (3,1):raise ValueError('rotation center must be a three-vector')
        result = eye(4)
        translation = vector-matrix*vector
        for i in range(3):
            for j in range(3):result[i,j] = matrix[i,j]
            result[i,3] = translation[i]
        return result
    @classmethod
    def from_rotation_matrix(cls,M):
        from .. import sqrt
        if M.shape != (3,3):raise ValueError('expected a 3 by 3 rotation matrix')
        trace = M.trace()
        # Choose a nonzero diagonal branch, including rotations through pi.
        candidates = [1+trace,1+2*M[0,0]-trace,1+2*M[1,1]-trace,1+2*M[2,2]-trace]
        k = next((i for i,c in enumerate(candidates) if c.is_positive is True),0)
        r = sqrt(candidates[k])/2
        if r == 0:raise ValueError('degenerate rotation matrix')
        if k == 0:return cls(r,(M[2,1]-M[1,2])/(4*r),(M[0,2]-M[2,0])/(4*r),(M[1,0]-M[0,1])/(4*r))
        if k == 1:return cls((M[2,1]-M[1,2])/(4*r),r,(M[0,1]+M[1,0])/(4*r),(M[0,2]+M[2,0])/(4*r))
        if k == 2:return cls((M[0,2]-M[2,0])/(4*r),(M[0,1]+M[1,0])/(4*r),r,(M[1,2]+M[2,1])/(4*r))
        return cls((M[1,0]-M[0,1])/(4*r),(M[0,2]+M[2,0])/(4*r),(M[1,2]+M[2,1])/(4*r),r)
    @staticmethod
    def rotate_point(pin,r):
        q = Quaternion.from_axis_angle(r) if isinstance(r,(tuple,list)) else r
        if not isinstance(q,Quaternion):raise TypeError('rotation must be a quaternion or axis-angle pair')
        point = tuple(pin)
        if len(point) != 3:raise ValueError('point must have three components')
        result = q*Quaternion(0,*point)*q.inverse()
        return tuple(v.simplify() for v in result.args[1:])
