"""Exact lattices in a rational power basis, for integral orders and ideals."""
from functools import cached_property
from math import gcd,lcm
from ... import Matrix, Rational, Integer, eye, Poly, ZZ, QQ


def _egcd(a,b):
    x,y,u,v=1,0,0,1
    while b:
        q,a,b=a//b,b,a%b;x,u=u,x-q*u;y,v=v,y-q*v
    return (a,x,y) if a>=0 else (-a,-x,-y)


def _hnf(matrix,full_rank=True):
    rows,cols=matrix.rows,matrix.cols
    a=[[int(matrix[i,j]) for j in range(cols)] for i in range(rows)]
    k=cols
    for i in range(rows-1,-1,-1):
        if k==0:
            if full_rank:raise ValueError('module lattice must have full rank')
            break
        k-=1
        pivot=next((j for j in range(k,-1,-1) if a[i][j]),None)
        if pivot is None:
            if full_rank:raise ValueError('module lattice must have full rank')
            k+=1;continue
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


from abc import ABC, abstractmethod


class Module(ABC):
    @abstractmethod
    def basis_elements(self):raise NotImplementedError
    def __call__(self,col,denom=1):return make_mod_elt(self,col,denom)


class ModuleElement:
    def __init__(self,module,col,denom=1):
        if hasattr(col,'to_Matrix'):col=col.to_Matrix()
        self.module=module;self._vector=tuple(Rational(v)/denom for v in col)
        if len(self._vector)!=module.n:raise ValueError('coefficient column has the wrong size')
    @property
    def denom(self):return lcm(*(int(v.q) for v in self._vector))
    @property
    def coeffs(self):return [int(v*self.denom) for v in self._vector]
    @property
    def col(self):return to_col(self.coeffs)
    def to_parent(self):
        if isinstance(self.module,PowerBasis):return self
        basis=self.module.basis_elements();result=basis[0].module([0])
        for a,b in zip(self._vector,basis):result+=a*b
        return result
    def as_expr(self):return self.to_parent().as_expr()
    def __add__(self,other):
        if not isinstance(other,ModuleElement) or other.module!=self.module:raise ValueError('incompatible modules')
        return make_mod_elt(self.module,[a+b for a,b in zip(self._vector,other._vector)])
    def __neg__(self):return make_mod_elt(self.module,[-a for a in self._vector])
    def __sub__(self,other):return self+-other
    def __mul__(self,other):
        if isinstance(other,ModuleElement):
            product=self.to_parent()*other.to_parent();return make_mod_elt(self.module,list(self.module.represent(product)))
        return make_mod_elt(self.module,[a*other for a in self._vector])
    __rmul__=__mul__
    def __eq__(self,other):return isinstance(other,ModuleElement) and self.module==other.module and self._vector==other._vector


class PowerBasis(Module):
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


class PowerBasisElement(ModuleElement):
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


class Submodule(Module):
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


def to_col(coeffs):
    from ..matrices.domainmatrix import DomainMatrix
    values=list(coeffs)
    return DomainMatrix([[ZZ.convert(v)] for v in values],(len(values),1),ZZ)

def make_mod_elt(module,col,denom=1):
    if hasattr(col,'to_Matrix'):col=col.to_Matrix()
    return module(list(col),denom) if isinstance(module,PowerBasis) else ModuleElement(module,list(col),denom)

def is_sq_maxrank_HNF(dm):
    matrix=dm.to_Matrix() if hasattr(dm,'to_Matrix') else Matrix(dm)
    if matrix.rows!=matrix.cols:return False
    for i in range(matrix.rows):
        if matrix[i,i].is_positive is not True:return False
        for j in range(matrix.cols):
            if i>j and matrix[i,j]!=0:return False
            if j>i and not (bool(matrix[i,j]>=0) and bool(matrix[i,j]<matrix[i,i])):return False
    return True


class _SpanModule(Module):
    def __init__(self,parent,basis):
        self.parent=parent;basis=list(basis)
        if basis:
            columns=Matrix.hstack(*(Matrix(e._vector) if e.module==parent else parent.represent(e) for e in basis))
            denominator=lcm(*(int(v.q) for v in columns))
            columns=_hnf(columns*denominator,full_rank=False)/denominator
            self._basis=[make_mod_elt(parent,list(columns[:,j])) for j in range(columns.cols)]
        else:self._basis=[]
        self.n=len(self._basis)
    def basis_elements(self):return list(self._basis)
    def represent(self,element):
        if not self._basis:
            if any(element._vector):raise ValueError('element is outside the zero module')
            return Matrix.zeros(0,1)
        matrix=Matrix.hstack(*(Matrix(e._vector) for e in self._basis))
        solutions=matrix.gauss_jordan_solve(Matrix(element._vector))[0]
        if any(v.is_Integer is not True for v in solutions):raise ValueError('element is outside the integral span')
        return solutions


class ModuleHomomorphism:
    def __init__(self,domain,codomain,mapping):
        if not callable(mapping):raise TypeError('module mapping must be callable')
        self.domain,self.codomain,self.mapping=domain,codomain,mapping
    def __call__(self,element):
        if element.module==self.domain:return self.mapping(element)
        if hasattr(self.domain,'represent'):
            return self.mapping(make_mod_elt(self.domain,list(self.domain.represent(element))))
        raise ValueError('element is outside the homomorphism domain')
    def matrix(self):
        from ..matrices.domainmatrix import DomainMatrix
        columns=[]
        for e in self.domain.basis_elements():
            image=self.mapping(e)
            col=image._vector if image.module==self.codomain else self.codomain.represent(image)
            columns.append(Matrix(col))
        matrix=Matrix.hstack(*columns) if columns else Matrix.zeros(self.codomain.n,0)
        return DomainMatrix([[QQ.convert(v) for v in row] for row in matrix.tolist()],matrix.shape,QQ)
    def kernel(self,modulus=None):
        if modulus is not None:
            matrix=self.matrix().to_Matrix();m=int(modulus)
            if m<2:raise ValueError('kernel modulus must exceed one')
            if m**self.domain.n>1000000:raise ValueError('explicit residue kernel exceeds one million elements')
            # Enumerate a finite residue module exactly, then include m times
            # each basis vector to describe its full integral preimage.
            from itertools import product
            vectors=[]
            for values in product(range(m),repeat=self.domain.n):
                if all(int(v)%m==0 for v in matrix*Matrix(values)):vectors.append(make_mod_elt(self.domain,values))
            vectors += [m*e for e in self.domain.basis_elements()]
            return _SpanModule(self.domain,vectors)
        kernel=_integer_kernel(self.matrix().to_Matrix())
        return _SpanModule(self.domain,[make_mod_elt(self.domain,list(kernel[:,j])) for j in range(kernel.cols)])
    def image(self):
        images=[self.mapping(e) for e in self.domain.basis_elements()]
        return _SpanModule(self.codomain,images)
    def __add__(self,other):return ModuleHomomorphism(self.domain,self.codomain,lambda e:self(e)+other(e))
    def __mul__(self,other):
        if other.codomain!=self.domain:raise ValueError('homomorphism composition domains differ')
        return ModuleHomomorphism(other.domain,self.codomain,lambda e:self(other(e)))


class ModuleEndomorphism(ModuleHomomorphism):
    def __init__(self,domain,mapping):super().__init__(domain,domain,mapping)
class InnerEndomorphism(ModuleEndomorphism):
    def __init__(self,domain,multiplier):self.multiplier=multiplier;super().__init__(domain,lambda e:multiplier*e)
class EndomorphismRing:
    def __init__(self,domain):self.domain=domain
    def __call__(self,mapping):return mapping if isinstance(mapping,ModuleEndomorphism) and mapping.domain==self.domain else ModuleEndomorphism(self.domain,mapping)
    def one(self):return ModuleEndomorphism(self.domain,lambda e:e)
    def zero(self):return ModuleEndomorphism(self.domain,lambda e:0*e)


def find_min_poly(alpha,domain,x=None,powers=None):
    from ... import Dummy
    x=Dummy('x') if x is None else x
    basis=alpha.module if isinstance(alpha.module,PowerBasis) else alpha.module.parent
    value=alpha if isinstance(alpha.module,PowerBasis) else alpha.to_parent()
    terms=[basis.one()];current=terms[0]
    for degree in range(1,basis.n+1):
        current=current*value;terms.append(current)
        if powers is not None:powers[degree]=current
        columns=Matrix.hstack(*(Matrix(e._vector) for e in terms));null=columns.nullspace()
        relation=next((v for v in null if v[-1]!=0),None)
        if relation is not None:
            polynomial=sum((relation[i]/relation[-1])*x**i for i in range(degree+1))
            return Poly(polynomial,x,domain=domain)
    raise ValueError('element powers did not produce a finite algebraic relation')

__all__=['Module','ModuleElement','PowerBasis','PowerBasisElement','Submodule','ModuleHomomorphism',
 'ModuleEndomorphism','InnerEndomorphism','EndomorphismRing','find_min_poly','is_sq_maxrank_HNF','make_mod_elt','to_col']


def _integer_kernel(matrix):
    """A saturated integral kernel from unimodular column transformations."""
    basis=eye(matrix.cols)
    for row in range(matrix.rows):
        values=list((matrix[row,:]*basis))
        if not any(values):continue
        denominator=lcm(*(int(v.q) for v in values));values=[int(v*denominator) for v in values]
        transform=eye(len(values))
        for j in range(1,len(values)):
            a,b=values[0],values[j]
            if b==0:continue
            g,s,t=_egcd(a,b)
            first,other=transform[:,0],transform[:,j]
            transform[:,0]=s*first+t*other
            transform[:,j]=-(b//g)*first+(a//g)*other
            values[0],values[j]=g,0
        basis=basis*transform[:,1:]
        if not basis.cols:break
    return basis
