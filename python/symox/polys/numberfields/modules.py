"""Exact lattices in a rational power basis, for integral orders and ideals."""
from functools import cached_property
from math import gcd,lcm
from ... import Matrix, Rational, Integer, eye, Poly, ZZ, QQ


def _egcd(a,b):
    x,y,u,v=1,0,0,1
    while b:
        q,a,b=a//b,b,a%b;x,u=u,x-q*u;y,v=v,y-q*v
    return (a,x,y) if a>=0 else (-a,-x,-y)


def _hnf(matrix):
    rows,cols=matrix.rows,matrix.cols
    a=[[int(matrix[i,j]) for j in range(cols)] for i in range(rows)]
    k=cols
    for i in range(rows-1,-1,-1):
        k-=1
        pivot=next((j for j in range(k,-1,-1) if a[i][j]),None)
        if pivot is None:raise ValueError('module lattice must have full rank')
        for row in a:row[k],row[pivot]=row[pivot],row[k]
        for j in range(k):
            if not a[i][j]:continue
            av,bv=a[i][k],a[i][j];g,s,t=_egcd(av,bv)
            for row in a:
                u,v=row[k],row[j]
                row[k],row[j]=s*u+t*v,-(bv//g)*u+(av//g)*v
        if a[i][k]<0:
            for row in a:row[k]=-row[k]
        for j in range(k+1,cols):
            q=a[i][j]//a[i][k]
            for row in a:row[j]-=q*row[k]
    return Matrix([row[k:] for row in a])


def _normalize(matrix):
    denom=1
    for value in matrix:denom=lcm(denom,int(value.q))
    numerator=_hnf(matrix*denom)
    common=denom
    for value in numerator:common=gcd(common,int(value))
    return numerator/common,denom//common


class _LatticeMatrix:
    """Small matrix view with the DomainMatrix conversion protocol."""
    def __init__(self,value,domain):self._value,self.domain=value,domain
    @property
    def shape(self):return self._value.shape
    def to_Matrix(self):return self._value.copy()
    def __getitem__(self,key):return self._value[key]
    def det(self):return self._value.det()
    def inv(self):return _LatticeMatrix(self._value.inv(),QQ)
    def __mul__(self,other):return _LatticeMatrix(self._value*other._value,QQ)
    def __repr__(self):return repr(self._value)


class PowerBasis:
    def __init__(self,T):
        self.T=T;self.n=int(T.degree());self._coeffs=list(reversed(T.all_coeffs()))
    def __eq__(self,other):return isinstance(other,PowerBasis) and self.T==other.T
    def __hash__(self):return hash(self.T)
    def __repr__(self):return 'PowerBasis(%s)'%self.T.as_expr()
    def __call__(self,coeffs,denom=1):
        from ... import sympify
        values=[sympify(v)/denom for v in coeffs]
        if len(values)>self.n:raise ValueError('too many power basis coefficients')
        return PowerBasisElement(self,values+[Integer(0)]*(self.n-len(values)))
    def element_from_poly(self,poly):
        poly=Poly(poly,self.T.gens[0],domain=QQ).rem(self.T)
        return self([poly.nth(i) for i in range(self.n)])
    def one(self):return self([1])
    def basis_elements(self):return [self([0]*i+[1]) for i in range(self.n)]
    def whole_submodule(self):return Submodule(self,eye(self.n))
    def _mul(self,a,b):
        n=self.n;coeff=[Integer(0)]*(2*n-1)
        for i,av in enumerate(a):
            for j,bv in enumerate(b):coeff[i+j]+=av*bv
        for i in range(2*n-2,n-1,-1):
            c=coeff[i]
            for j in range(n):coeff[i-n+j]-=c*self._coeffs[j]
        return coeff[:n]
    def mult_matrix(self,vector):
        return Matrix.hstack(*(Matrix(self._mul(vector,[Integer(int(j==i)) for j in range(self.n)])) for i in range(self.n)))


class PowerBasisElement:
    def __init__(self,module,vector):self.module,self._vector=module,tuple(vector)
    @property
    def denom(self):return lcm(*(int(v.q) for v in self._vector))
    @property
    def coeffs(self):return [int(v*self.denom) for v in self._vector]
    def as_expr(self):return sum((v*self.module.T.gens[0]**i for i,v in enumerate(self._vector)),Integer(0))
    def __repr__(self):return str(self.coeffs)+('/%d'%self.denom if self.denom!=1 else '')
    def __eq__(self,other):return isinstance(other,PowerBasisElement) and self.module==other.module and self._vector==other._vector
    def __hash__(self):return hash((self.module,self._vector))
    def __add__(self,other):
        if not isinstance(other,PowerBasisElement):other=self.module([other])
        if self.module!=other.module:raise ValueError('incompatible power bases')
        return self.module([a+b for a,b in zip(self._vector,other._vector)])
    __radd__=__add__
    def __neg__(self):return self.module([-v for v in self._vector])
    def __sub__(self,other):return self+(-other)
    def __mul__(self,other):
        if isinstance(other,Submodule):return other*self
        if not isinstance(other,PowerBasisElement):return self.module([v*other for v in self._vector])
        if self.module!=other.module:raise ValueError('incompatible power bases')
        return self.module(self.module._mul(self._vector,other._vector))
    __rmul__=__mul__
    def __truediv__(self,other):return self.module([v/other for v in self._vector])
    def __pow__(self,n):
        import operator
        n=operator.index(n)
        if n<0:raise NotImplementedError('negative powers of power basis elements')
        result=self.module.one();base=self
        while n:
            if n&1:result*=base
            base*=base;n//=2
        return result


class Submodule:
    def __init__(self,parent,matrix):
        self.parent=parent;self.n=parent.n
        numerator,self.denom=_normalize(matrix)
        self._matrix=numerator
        self.matrix=_LatticeMatrix(numerator,ZZ)
        self.QQ_matrix=_LatticeMatrix(numerator/self.denom,QQ)
    def __repr__(self):return 'Submodule%s/%s'%(self._matrix.T.tolist(),self.denom)
    def __eq__(self,other):return isinstance(other,Submodule) and self.parent==other.parent and self.denom==other.denom and self._matrix==other._matrix
    def basis_elements(self):return [self.parent(list(self._matrix[:,j]),self.denom) for j in range(self.n)]
    @cached_property
    def _inverse(self):return (self._matrix/self.denom).inv()
    def represent(self,element):
        if element.module!=self.parent:raise ValueError('element belongs to another power basis')
        vector=self._inverse*Matrix(element._vector)
        if any(not v.is_Integer for v in vector):raise ValueError('element is outside the module')
        return vector
    def contains(self,element):
        try:self.represent(element);return True
        except ValueError:return False
    def __contains__(self,element):return self.contains(element)
    def __add__(self,other):
        if self.parent!=other.parent:raise ValueError('incompatible power bases')
        return Submodule(self.parent,Matrix.hstack(self._matrix/self.denom,other._matrix/other.denom))
    def __mul__(self,other):
        if hasattr(other,'as_submodule'):other=other.as_submodule()
        if isinstance(other,Submodule):
            if self.parent!=other.parent:raise ValueError('incompatible power bases')
            columns=[Matrix((a*b)._vector) for a in self.basis_elements() for b in other.basis_elements()]
            return Submodule(self.parent,Matrix.hstack(*columns))
        if isinstance(other,PowerBasisElement):
            return Submodule(self.parent,Matrix.hstack(*(Matrix((a*other)._vector) for a in self.basis_elements())))
        if other==0:raise NotImplementedError('zero-rank ideals are not represented as full lattices')
        return Submodule(self.parent,(self._matrix/self.denom)*other)
    __rmul__=__mul__
    def __truediv__(self,other):return self*Rational(1,other)
    def __pow__(self,n):
        import operator
        n=operator.index(n)
        if n<1:raise ValueError('module powers require a positive integer')
        result=self
        for _ in range(n-1):result=result*self
        return result
