"""Certified integral lattices by bounded prime saturation.

This exact enumeration complements the existing native number-field kernels.
It is intended for small degrees, not as a large-degree Round Two backend.
"""
from itertools import product
from ... import Poly, QQ, ZZ, Matrix, Integer, Rational, eye
from .modules import PowerBasis, Submodule


def _polynomial(T):
    T=T if isinstance(T,Poly) else Poly(T)
    if len(T.gens)!=1 or T.domain not in (ZZ,QQ) or not T.is_irreducible:
        raise ValueError('an irreducible polynomial over ZZ or QQ is required')
    T=T.monic()
    if not all(c.is_Integer for c in T.all_coeffs()):
        raise NotImplementedError('a monic integral defining polynomial is required')
    if not 1<=T.degree()<=8:raise NotImplementedError('integral bases are limited to degrees 1 through 8')
    return T.set_domain(ZZ)


def round_two(T,radicals=None):
    """Return an exact maximal order lattice and its field discriminant.

    The order is saturated by testing integral classes in (1/p)L/L. Every
    enlargement is certified by an integral characteristic polynomial.
    Projective residue enumeration is limited to 100000 classes per prime.
    ``radicals`` is accepted as a cache container and remains unmodified.
    """
    from ... import factorint
    T=_polynomial(T);n=int(T.degree());basis=PowerBasis(T)
    B=eye(n)
    monomials=basis.basis_elements()
    trace=Matrix([[basis.mult_matrix((a*b)._vector).trace() for b in monomials] for a in monomials])
    disc=int(T.discriminant())
    primes=[int(p) for p,e in factorint(abs(disc)).items() if e>=2]
    for p in primes:
        if (p**n-1)//(p-1)>100000:
            raise NotImplementedError('integral basis residue enumeration exceeds 100000 projective classes')
        while disc%(p*p)==0:
            enlarged=None
            # Projective normalization: the first nonzero residue is one.
            for pivot in range(n):
                for tail in product(range(p),repeat=n-pivot-1):
                    vector=Matrix([0]*pivot+[1]+list(tail))
                    beta=B*vector/p
                    cp=basis.mult_matrix(list(beta)).charpoly().all_coeffs()
                    if all(c.is_Integer for c in cp):
                        enlarged=(pivot,beta);break
                if enlarged is not None:break
            if enlarged is None:break
            pivot,beta=enlarged
            B=Matrix.hstack(*(beta if j==pivot else B[:,j] for j in range(n)))
            disc//=p*p
    order=Submodule(basis,B)
    certified=(order.QQ_matrix.to_Matrix().T*trace*order.QQ_matrix.to_Matrix()).det()
    if certified!=disc:raise ArithmeticError('integral basis discriminant certificate failed')
    return order,Integer(disc)
