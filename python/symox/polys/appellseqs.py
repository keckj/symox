"""Exact Appell polynomial sequences built from native Bernoulli/Euler numbers."""
import operator
import symox as S


def _argument(n,x):
    n=operator.index(n)
    if n<0:raise ValueError('polynomial degree must be nonnegative')
    return n,S.Dummy('x') if x is None else S.sympify(x)

def _bernoulli(n,x):
    return S.expand(sum(S.binomial(n,k)*(-S.Rational(1,2) if k==1 else S.bernoulli(k))*x**(n-k) for k in range(n+1)))

def bernoulli_poly(n,x=None,polys=False):
    n,x=_argument(n,x);result=_bernoulli(n,x)
    return S.Poly(result,x) if polys else result

def bernoulli_c_poly(n,x=None,polys=False):
    n,x=_argument(n,x);result=S.expand(2**n*_bernoulli(n,(x+1)/2))
    return S.Poly(result,x) if polys else result

def euler_poly(n,x=None,polys=False):
    n,x=_argument(n,x);result=S.expand(2*S.Rational(1,n+1)*(_bernoulli(n+1,x)-2**(n+1)*_bernoulli(n+1,x/2)))
    return S.Poly(result,x) if polys else result

def genocchi_poly(n,x=None,polys=False):
    n,x=_argument(n,x);result=S.Integer(0) if n==0 else S.expand(-n*euler_poly(n-1,x))
    return S.Poly(result,x) if polys else result

def andre_poly(n,x=None,polys=False):
    n,x=_argument(n,x);result=S.expand(sum(S.binomial(n,k)*S.euler(k)*x**(n-k) for k in range(n+1)))
    return S.Poly(result,x) if polys else result

__all__=['bernoulli_poly','bernoulli_c_poly','euler_poly','genocchi_poly','andre_poly']
