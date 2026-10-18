"""Exact angular-momentum coupling from finite factorial sums."""
from functools import lru_cache
from itertools import product
import symox as S

H = S.Rational(1,2)


def _half(value):
    value = S.sympify(value)
    if not (2*value).is_Integer:raise ValueError('angular momenta must be explicit integers or half-integers')
    return value


def _fact(value):
    value = S.sympify(value)
    if not value.is_Integer or value<0:return S.Integer(0)
    return S.factorial(value)


def _triangle(a,b,c):return a>=0 and b>=0 and c>=0 and abs(a-b)<=c<=a+b and (a+b+c).is_Integer


def _delta(a,b,c):
    if not _triangle(a,b,c):return S.Integer(0)
    return S.sqrt(_fact(a+b-c)*_fact(a-b+c)*_fact(-a+b+c)/_fact(a+b+c+1))


@lru_cache(maxsize=32768)
def wigner_3j(j_1,j_2,j_3,m_1,m_2,m_3):
    a,b,c,p,q,r = map(_half,(j_1,j_2,j_3,m_1,m_2,m_3))
    if not _triangle(a,b,c) or p+q+r!=0 or any(abs(m)>j or not(j-m).is_Integer for j,m in ((a,p),(b,q),(c,r))):return S.Integer(0)
    lower = max(0,int(b-c-p),int(a-c+q));upper = min(int(a+b-c),int(a-p),int(b+q))
    coefficient = S.Integer(-1)**int(a-b-r)*_delta(a,b,c)*S.sqrt(S.Mul(*(_fact(j+m)*_fact(j-m) for j,m in ((a,p),(b,q),(c,r)))))
    total = S.Add(*(S.Integer(-1)**z/(_fact(z)*_fact(a+b-c-z)*_fact(a-p-z)*_fact(b+q-z)*_fact(c-b+p+z)*_fact(c-a-q+z)) for z in range(lower,upper+1)))
    return S.simplify(coefficient*total)


def clebsch_gordan(j_1,j_2,j_3,m_1,m_2,m_3):
    a,b,c,p,q,r = map(_half,(j_1,j_2,j_3,m_1,m_2,m_3))
    if p+q!=r:return S.Integer(0)
    return S.simplify(S.Integer(-1)**int(a-b+r)*S.sqrt(2*c+1)*wigner_3j(a,b,c,p,q,-r))


@lru_cache(maxsize=32768)
def wigner_6j(j_1,j_2,j_3,j_4,j_5,j_6,prec=None):
    a,b,c,d,e,f = map(_half,(j_1,j_2,j_3,j_4,j_5,j_6))
    triangles = ((a,b,c),(a,e,f),(d,b,f),(d,e,c))
    if any(not(x+y+z).is_Integer for x,y,z in triangles):raise ValueError('triangle sums must be integers')
    coefficient = S.Mul(*(_delta(*triple) for triple in triangles))
    if coefficient==0:return coefficient
    first = (a+b+c,a+e+f,d+b+f,d+e+c);second = (a+b+d+e,a+c+d+f,b+c+e+f)
    total = S.Add(*(S.Integer(-1)**z*_fact(z+1)/(S.Mul(*(_fact(z-v) for v in first))*S.Mul(*(_fact(v-z) for v in second))) for z in range(int(max(first)),int(min(second))+1)))
    value = S.simplify(coefficient*total)
    return value.evalf(prec) if prec else value


def racah(aa,bb,cc,dd,ee,ff,prec=None):
    a,b,c,d,e,f = map(_half,(aa,bb,cc,dd,ee,ff))
    value = S.Integer(-1)**int(a+b+c+d)*wigner_6j(a,b,e,d,c,f)
    return value.evalf(prec) if prec else value


def wigner_9j(j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8,j_9,prec=None):
    a,b,c,d,e,f,g,h,i = map(_half,(j_1,j_2,j_3,j_4,j_5,j_6,j_7,j_8,j_9))
    upper = int(2*min(a+i,b+f,d+h));terms = []
    for twice in range(upper%2,upper+1,2):
        x = S.Rational(twice,2)
        try:terms.append((twice+1)*racah(a,b,i,f,c,x)*racah(d,f,h,b,e,x)*racah(a,d,i,h,g,x))
        except ValueError:continue
    value = S.simplify(S.Add(*terms))
    return value.evalf(prec) if prec else value


def gaunt(l_1,l_2,l_3,m_1,m_2,m_3,prec=None):
    values = tuple(map(S.sympify,(l_1,l_2,l_3,m_1,m_2,m_3)))
    if any(not v.is_Integer for v in values):raise ValueError('Gaunt indices must be integers')
    a,b,c,p,q,r = values
    value = S.sqrt((2*a+1)*(2*b+1)*(2*c+1)/(4*S.pi))*wigner_3j(a,b,c,0,0,0)*wigner_3j(a,b,c,p,q,r)
    value = S.simplify(value)
    return value.evalf(prec) if prec else value


def real_gaunt(l_1,l_2,l_3,m_1,m_2,m_3,prec=None):
    degrees = tuple(map(S.sympify,(l_1,l_2,l_3)));orders = tuple(map(S.sympify,(m_1,m_2,m_3)))
    if any(not v.is_Integer for v in degrees+orders):raise ValueError('real Gaunt indices must be integers')
    def basis(m):
        if m==0:return [(0,S.Integer(1))]
        if m>0:return [(m,1/S.sqrt(2)),(-m,S.Integer(-1)**int(m)/S.sqrt(2))]
        return [(m,1/(S.I*S.sqrt(2))),(-m,-S.Integer(-1)**int(m)/(S.I*S.sqrt(2)))]
    value = S.simplify(S.Add(*(S.Mul(*(v[1] for v in combination))*gaunt(*degrees,*(v[0] for v in combination)) for combination in product(*(basis(m) for m in orders)))))
    return value.evalf(prec) if prec else value


def wigner_d_small(J,beta):
    J,beta = _half(J),S.sympify(beta)
    if J<0:raise ValueError('spin must be nonnegative')
    states = [J-k for k in range(int(2*J)+1)]
    def element(row,col):
        p,m = states[row],states[col];terms = []
        coefficient = S.sqrt(_fact(J+m)*_fact(J-m)*_fact(J+p)*_fact(J-p))
        for s in range(max(0,int(m-p)),min(int(J+m),int(J-p))+1):
            terms.append(S.Integer(-1)**s*S.cos(beta/2)**(2*J+m-p-2*s)*S.sin(beta/2)**(p-m+2*s)/(_fact(J+m-s)*_fact(s)*_fact(p-m+s)*_fact(J-p-s)))
        return S.simplify(coefficient*S.Add(*terms))
    return S.Matrix(len(states),len(states),element)


def wigner_d(J,alpha,beta,gamma):
    J = _half(J);small = wigner_d_small(J,beta)
    return S.Matrix(small.rows,small.cols,lambda i,j:S.exp(S.I*(J-i)*alpha)*small[i,j]*S.exp(S.I*(J-j)*gamma))


class Wigner3j:
    def __init__(self,*args):
        if len(args)!=6:raise TypeError('Wigner3j requires six indices')
        self.args = tuple(map(S.sympify,args));self._expr = S._core._apply('Wigner3j',*self.args,evaluate=False)
    def _symox_(self):return self._expr
    def __repr__(self):return str(self._expr)
    @property
    def free_symbols(self):return self._expr.free_symbols
    def doit(self,**hints):return self if self.free_symbols else wigner_3j(*self.args)
    def subs(self,*args,**kwargs):return Wigner3j(*(v.subs(*args,**kwargs) for v in self.args))
    def __mul__(self,other):return self._expr*other
    def __rmul__(self,other):return other*self._expr


def dot_rot_grad_Ynm(j,p,l,m,theta,phi):
    j,p,l,m,theta,phi = map(S.sympify,(j,p,l,m,theta,phi))
    def term(k):
        factor = S.sqrt((2*l+1)*(2*j+1)*(2*k+1)/(4*S.pi))
        first,second = Wigner3j(j,l,k,0,0,0),Wigner3j(j,l,k,p,m,-m-p)
        a,b = first.doit(),second.doit()
        return (-1)**(m+p)*S.Ynm(k,m+p,theta,phi)*factor*a*b*(k*(k+1)-j*(j+1)-l*(l+1))/2
    if j.is_Integer and l.is_Integer:return S.simplify(S.Add(*(term(S.Integer(k)) for k in range(abs(int(l-j)),int(l+j)+1))))
    k = S.Dummy('k',integer=True)
    return S.Sum(term(k),(k,S.Abs(l-j),l+j))


__all__ = ['Wigner3j','clebsch_gordan','dot_rot_grad_Ynm','gaunt','racah','real_gaunt','wigner_3j','wigner_6j','wigner_9j','wigner_d','wigner_d_small']
