"""ABCD ray optics and Gaussian beam propagation."""
import symox as S


class RayTransferMatrix:
    _symox_matrixexpr = True
    def __init__(self,*args):
        # Reference constructor: the four elements A, B, C, D, or one 2x2 matrix.
        if len(args) == 4:
            args = ([[args[0], args[1]], [args[2], args[3]]],)
        self._matrix = S.Matrix(*args)
        if self._matrix.shape!=(2,2):raise ValueError('ray transfer matrix must be 2 by 2')
    def _symox_(self):return S.MatrixExpr(self._matrix)
    def __repr__(self):return repr(self._matrix)
    @property
    def shape(self):return (2,2)
    rows=cols=2
    @property
    def A(self):return self._matrix[0,0]
    @property
    def B(self):return self._matrix[0,1]
    @property
    def C(self):return self._matrix[1,0]
    @property
    def D(self):return self._matrix[1,1]
    def __getitem__(self,key):return self._matrix[key]
    def __iter__(self):return iter(self._matrix)
    def __eq__(self,other):return self._matrix==(other._matrix if isinstance(other,RayTransferMatrix) else other)
    def __mul__(self,other):
        if isinstance(other,BeamParameter):return BeamParameter._from_q(other.wavelen,(self.A*other.q+self.B)/(self.C*other.q+self.D),other.n)
        if isinstance(other,GeometricRay):return GeometricRay(self._matrix*other._matrix)
        if isinstance(other,RayTransferMatrix):return RayTransferMatrix(self._matrix*other._matrix)
        return self._matrix*other
    def __rmul__(self,other):return other*self._matrix
    def __pow__(self,power):return RayTransferMatrix(self._matrix**power)
    def det(self):return self._matrix.det()
    def inv(self,*args,**kwargs):return RayTransferMatrix(self._matrix.inv(*args,**kwargs))
    def as_explicit(self):return self._matrix.copy()
    def subs(self,*args,**kwargs):return RayTransferMatrix(self._matrix.subs(*args,**kwargs))
    def simplify(self):return RayTransferMatrix(self._matrix.applyfunc(S.simplify))


class FreeSpace(RayTransferMatrix):
    def __init__(self,d):super().__init__([[1,d],[0,1]])
class FlatRefraction(RayTransferMatrix):
    def __init__(self,n1,n2):super().__init__([[1,0],[0,S.sympify(n1)/n2]])
class CurvedRefraction(RayTransferMatrix):
    def __init__(self,R,n1,n2):
        R,n1,n2=map(S.sympify,(R,n1,n2));super().__init__([[1,0],[(n1-n2)/(R*n2),n1/n2]])
class FlatMirror(RayTransferMatrix):
    def __init__(self):super().__init__(S.eye(2))
class CurvedMirror(RayTransferMatrix):
    def __init__(self,R):super().__init__([[1,0],[-2/S.sympify(R),1]])
class ThinLens(RayTransferMatrix):
    def __init__(self,f):super().__init__([[1,0],[-1/S.sympify(f),1]])


class GeometricRay:
    _symox_matrixexpr = True
    def __init__(self,*args):
        self._matrix = S.Matrix(list(args) if len(args)==2 and not isinstance(args[0],(list,tuple,S.Matrix)) else args[0])
        if self._matrix.shape!=(2,1):raise ValueError('a geometric ray has height and angle components')
    def _symox_(self):return S.MatrixExpr(self._matrix)
    def __repr__(self):return repr(self._matrix)
    rows,cols=2,1
    @property
    def shape(self):return (2,1)
    @property
    def height(self):return self._matrix[0]
    @property
    def angle(self):return self._matrix[1]
    def __getitem__(self,key):return self._matrix[key]
    def as_explicit(self):return self._matrix.copy()
    def __eq__(self,other):return self._matrix==(other._matrix if isinstance(other,GeometricRay) else other)


def waist2rayleigh(w,wavelen,n=1):return S.pi*S.sympify(n)*S.sympify(w)**2/S.sympify(wavelen)
def rayleigh2waist(z_r,wavelen):return S.sqrt(S.sympify(z_r)*S.sympify(wavelen)/S.pi)


class BeamParameter:
    def __init__(self,wavelen,z,z_r=None,w=None,n=1):
        self.wavelen,self.z,self.n = map(S.sympify,(wavelen,z,n))
        if z_r is None and w is None:raise ValueError('a Rayleigh range or waist radius is required')
        self.z_r = S.sympify(z_r) if z_r is not None else waist2rayleigh(w,self.wavelen,self.n)
        if any(value.is_nonpositive is True for value in (self.wavelen,self.z_r,self.n)):raise ValueError('wavelength, Rayleigh range and refractive index must be positive')
    @classmethod
    def _from_q(cls,wavelen,q,n):
        q=S.expand_complex(S.cancel(q));return cls(wavelen,S.simplify(S.re(q)),S.simplify(S.im(q)),n=n)
    @property
    def q(self):return self.z+S.I*self.z_r
    @property
    def w_0(self):return S.sqrt(self.wavelen*self.z_r/(S.pi*self.n))
    @property
    def w(self):return self.w_0*S.sqrt(1+(self.z/self.z_r)**2)
    @property
    def radius(self):return S.oo if self.z==0 else self.z+self.z_r**2/self.z
    @property
    def gouy(self):return S.atan(self.z/self.z_r)
    @property
    def waist_approximation_limit(self):return 2*self.wavelen/S.pi
    @property
    def divergence(self):return self.wavelen/(S.pi*self.n*self.w_0)
    @property
    def args(self):return self.wavelen,self.z,self.z_r,self.n
    def __repr__(self):return 'BeamParameter(%s, %s, %s, n=%s)'%self.args


def geometric_conj_ab(a,b):
    a,b=map(S.sympify,(a,b))
    if a.is_infinite or b.is_infinite:return a if b.is_infinite else b
    return a*b/(a+b)
def geometric_conj_af(a,f):return S.sympify(a)*S.sympify(f)/(S.sympify(a)-S.sympify(f))
def geometric_conj_bf(a,f):return geometric_conj_af(a,f)


def gaussian_conj(s_in,z_r_in,f):
    s,z,f=map(S.sympify,(s_in,z_r_in,f));denominator=(s-f)**2+z*z
    magnification=S.Abs(f)/S.sqrt(denominator)
    # s_out is SymPy's expression; the Rayleigh range and magnification use
    # (1 - s/f)**2 where SymPy has 1 - (s/f)**2 (see docs/algorithmic-parity.md).
    return 1/(-1/(s+z**2/(s-f))+1/f),S.simplify(z*f*f/denominator),S.simplify(magnification)


def conjugate_gauss_beams(wavelen,waist_in,waist_out,**kwargs):
    wavelen,waist_in,waist_out=map(S.sympify,(wavelen,waist_in,waist_out))
    m=waist_out/waist_in;z=waist2rayleigh(waist_in,wavelen)
    if len(kwargs)!=1:raise ValueError('The function expects only one named argument')
    elif 'dist' in kwargs:raise NotImplementedError('Currently only focal length is supported as a parameter')
    elif 'f' in kwargs:
        f=S.sympify(kwargs['f']);s_in=f*(1-S.sqrt(1/m**2-z**2/f**2));s_out=gaussian_conj(s_in,z,f)[0]
    elif 's_in' in kwargs:
        # Symox extension (SymPy 1.14 raises NotImplementedError): the focal
        # length placing the output waist for a given input distance.
        s_in=S.sympify(kwargs['s_in']);f=(s_in**2+z**2)/(s_in+S.sqrt((s_in**2+z**2)/m**2-z**2));s_out=gaussian_conj(s_in,z,f)[0]
    else:raise ValueError('The function expects only one named argument')
    return (s_in,s_out,f)


__all__=['RayTransferMatrix','FreeSpace','FlatRefraction','CurvedRefraction','FlatMirror','CurvedMirror','ThinLens','GeometricRay','BeamParameter','waist2rayleigh','rayleigh2waist','geometric_conj_ab','geometric_conj_af','geometric_conj_bf','gaussian_conj','conjugate_gauss_beams']
