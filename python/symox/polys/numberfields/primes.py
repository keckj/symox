"""Prime ideals for quadratic orders and Dedekind's good-index primes."""
import operator
from ... import Poly, ZZ, Integer, Matrix, eye, Symbol
from .basis import round_two, _polynomial
from .modules import PowerBasis, PowerBasisElement, Submodule


class PrimeIdeal:
    def __init__(self,ZK,p,alpha,f,e=None):
        from ... import isprime
        self.ZK,self.p,self.alpha,self.f=ZK,operator.index(p),alpha,operator.index(f)
        self.e=operator.index(e) if e is not None else None
        if not isprime(self.p) or self.f<1 or self.e is None or self.e<1:
            raise ValueError('prime ideal requires prime p and positive e and f')
        self._ideal=self.p*ZK if not any(alpha._vector) else self.p*ZK+alpha*ZK
        relative=abs(self._ideal.QQ_matrix.det()/ZK.QQ_matrix.det())
        if relative!=self.p**self.f:raise ArithmeticError('prime ideal norm certificate failed')

    def as_submodule(self):return self._ideal
    def is_inert(self):return self.f==self.ZK.n
    def __repr__(self):return '(%s, %s), e=%s, f=%s'%(self.p,self.alpha.as_expr(),self.e,self.f)
    def __eq__(self,other):return isinstance(other,PrimeIdeal) and self._ideal==other._ideal
    def __mul__(self,other):return self._ideal*other
    def __rmul__(self,other):return self._ideal*other
    def __pow__(self,n):
        n=operator.index(n)
        if n<0:raise NotImplementedError('negative prime ideal powers require inverse lattices')
        return self.ZK if n==0 else self._ideal**n
    def __add__(self,other):return self._ideal+(other._ideal if isinstance(other,PrimeIdeal) else other)
    def valuation(self,I):return prime_valuation(I,self)


def prime_decomp(p,T=None,ZK=None,dK=None,radical=None):
    """Factor p in the maximal order, with exact ideal-product certificates.

    Quadratic fields use an integral generator, including primes dividing
    the power-basis index. Higher degrees currently require a prime coprime
    to that index. Supplying a radical cache is not yet supported.
    """
    from ... import isprime
    p=operator.index(p)
    if not isprime(p):raise ValueError('p must be prime')
    if radical is not None:raise NotImplementedError('precomputed prime radicals are not supported')
    if T is None and ZK is None:raise ValueError('supply a defining polynomial or maximal order')
    if T is None:T=ZK.parent.T
    T=_polynomial(T)
    computed,disc=round_two(T)
    if ZK is None:ZK=computed
    elif ZK!=computed:raise ValueError('ZK is not the certified maximal order')
    if dK is not None and dK!=disc:raise ValueError('dK disagrees with the field discriminant')
    n=int(T.degree());basis=ZK.parent
    index_square=abs(int(T.discriminant()/disc))
    if n==2:
        theta=ZK.basis_elements()[1]
        characteristic=basis.mult_matrix(theta._vector).charpoly()
        relation=Poly(characteristic.as_expr(),characteristic.gen,domain=ZZ)
        # charpoly's private variable can be used directly by native Poly.
        generator=relation.gens[0]
    else:
        if index_square%(p*p)==0:
            raise NotImplementedError('higher-degree decomposition at primes dividing the power-basis index is not implemented')
        theta=basis([0,1]) if n>1 else basis([-T.nth(0)])
        relation=T;generator=T.gens[0]
    factors=relation.set_modulus(p).factor_list()[1]
    result=[]
    for factor,multiplicity in factors:
        f=int(factor.degree());alpha=basis([0])
        if f!=n:
            for coefficient in factor.all_coeffs():
                coefficient=int(coefficient)%p
                if 2*coefficient>p:coefficient-=p
                alpha=alpha*theta+coefficient
        result.append(PrimeIdeal(ZK,p,alpha,f,int(multiplicity)))
    if sum(P.e*P.f for P in result)!=n:raise ArithmeticError('prime degree certificate failed')
    reconstructed=ZK
    for P in result:reconstructed=reconstructed*(P**P.e)
    if reconstructed!=p*ZK:raise ArithmeticError('prime ideal factorization certificate failed')
    return result


def prime_valuation(I,P):
    """Exact exponent of P in a full fractional ideal lattice."""
    if not isinstance(P,PrimeIdeal):raise TypeError('P must be a prime ideal')
    if isinstance(I,PrimeIdeal):I=I.as_submodule()
    if isinstance(I,PowerBasisElement):I=I*P.ZK
    if isinstance(I,int) or getattr(I,'is_Integer',False):I=P.ZK*I
    if not isinstance(I,Submodule) or I.parent!=P.ZK.parent:
        raise ValueError('ideal and prime must have the same power basis')
    if I*P.ZK!=I:
        raise ValueError('lattice is not an ideal of the selected maximal order')
    # Clear rational lattice denominators using the exact integer denominator;
    # rational p-powers contribute e times their ordinary p-adic valuation.
    denominator=I.denom;shift=0;value=denominator
    while value%P.p==0:value//=P.p;shift+=P.e
    integral=I*denominator
    if any(not P.ZK.contains(b) for b in integral.basis_elements()):
        raise ValueError('lattice is not a fractional ideal of the selected order')
    count=0;power=P.as_submodule()
    for _ in range(1024):
        if not all(power.contains(b) for b in integral.basis_elements()):return count-shift
        count+=1;power=power*P.as_submodule()
    raise NotImplementedError('prime valuation exceeds 1024 lattice steps')
