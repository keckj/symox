"""Embeddings, maximal orders and prime factorizations need exact certificates."""
import pytest
import symox as s
import sympy as sp
from sympy.polys.numberfields import round_two as oracle_order, prime_decomp as oracle_primes
from symox.polys.numberfields._normalforms import hermite_normal_form as _dm_hnf
from symox.polys.matrices import DomainMatrix
from symox.polys.numberfields.modules import find_min_poly


def _hnf(M):
    return _dm_hnf(DomainMatrix.from_Matrix(M).convert_to(s.ZZ)).to_Matrix()

x=s.Symbol('x')
CASES=[x*x-5,x*x-7,x*x+1,x*x-12,x**3-x-1,
 x**3+x*x-2*x+8,x**3+3*x*x-4*x+4,x**3+3*x*x+3*x-3,
 x**3+5*x*x-x+3,x**3+5*x*x-5*x-5,x**3+3*x*x+5,
 x**3+6*x*x+3*x-1,x**3+6*x*x+4,x**3+7*x*x+7*x-7,
 x**3+7*x*x-x+5,x**3+7*x*x-5*x+5,x**3+4*x*x-3*x+7,
 x**3+8*x*x+5*x-1,x**3+8*x*x-2*x+6,x**3+6*x*x-3*x+8,
 x**3+9*x*x+6*x-8,x**3+15*x*x-9*x+13,
 x**4+x**3+x*x+x+1]


@pytest.mark.parametrize('expr',CASES)
def test_maximal_orders_match_oracle_lattices(expr):
    T=s.Poly(expr,x);Z,d=s.round_two(T)
    expected,disc=oracle_order(sp.Poly(expr.to_sympy(),x.to_sympy()))
    assert int(d)==disc
    B=Z.QQ_matrix.to_Matrix().to_sympy()
    comparison=expected.QQ_matrix.to_Matrix().inv()*B
    assert all(c.q==1 for c in comparison) and abs(comparison.det())==1
    # Module.represent raises ClosureFailure outside the order.
    assert Z.represent(Z.parent.one()).domain.is_ZZ
    for a in Z.basis_elements():
        for b in Z.basis_elements():
            assert Z.represent((a*b).to_parent()).domain.is_ZZ
        assert all(c.is_Integer for c in find_min_poly(a, s.QQ).all_coeffs())


@pytest.mark.parametrize('expr',[x*x-5,x*x+1,x*x-2,x*x-12,x**3-x-1,x**4+x**3+x*x+x+1])
@pytest.mark.parametrize('p',[2,3,5,7,11])
def test_prime_decomposition_and_fractional_valuations(expr,p):
    T=s.Poly(expr,x);Z,d=s.round_two(T)
    Ps=s.prime_decomp(p,ZK=Z,dK=d)
    expected=oracle_primes(p,T=sp.Poly(expr.to_sympy(),x.to_sympy()))
    assert sorted((P.e,P.f) for P in Ps)==sorted((P.e,P.f) for P in expected)
    product=Z
    for P in Ps:
        product=product*(P**P.e)
        assert P.as_submodule().QQ_matrix.det()/Z.QQ_matrix.det()==p**P.f
        assert s.prime_valuation(P**3,P)==3
        assert s.prime_valuation(p*Z,P)==P.e
        for Q in Ps:
            assert s.prime_valuation(Q.as_submodule(),P)==int(P==Q)
    assert product==p*Z


@pytest.mark.parametrize('a,b',[(s.sqrt(2),s.sqrt(2)+s.sqrt(3)),(s.I*s.sqrt(3),s.I*s.sqrt(3)/2),(3,s.sqrt(2)),(s.sqrt(2),s.sqrt(3)),(s.AlgebraicNumber(s.sqrt(2),[2,1]),s.sqrt(2))])
def test_field_embeddings(a,b):
    coeffs=s.field_isomorphism(a,b)
    expected=sp.polys.numberfields.field_isomorphism(s.sympify(a).to_sympy(),s.sympify(b).to_sympy())
    assert (coeffs is None)==(expected is None)
    if coeffs is not None:
        expr=s.Integer(0)
        for c in coeffs:expr=expr*b+c
        target=a.as_expr() if getattr(a,'is_AlgebraicNumber',False) else a
        assert s.simplify(expr-target)==0
    assert s.field_isomorphism(a,b,fast=False)==coeffs


@pytest.mark.parametrize('a',[s.sqrt(2),-s.sqrt(5),s.sqrt(2)-s.sqrt(3),s.real_root(2,3),s.Rational(2,3)])
def test_root_isolation_certificates(a):
    left,right=s.isolate(a,eps=s.Rational(1,10000))
    assert right-left<=s.Rational(1,10000)
    if left==right:assert a==left
    else:
        assert (a-left).is_positive is True and (right-a).is_positive is True
        poly=s.minpoly(a,polys=True)
        assert poly.count_roots(left,right)==1


def test_numberfield_errors_and_hnf_certificate():
    from sympy.matrices.normalforms import hermite_normal_form
    matrices=[[[2,0,1],[0,2,1]],[[3,5,7],[4,2,6]],[[2,0,0,1],[0,3,0,2],[0,0,5,3]]]
    for entries in matrices:
        assert _hnf(s.Matrix(entries)).to_sympy()==hermite_normal_form(sp.Matrix(entries))
    with pytest.raises(ValueError):s.round_two(s.Poly(x*x-1,x))
    with pytest.raises(ValueError):s.prime_decomp(4,T=s.Poly(x*x-2,x))
    # 2 is a common index divisor of this field: the reference splits it
    # with idempotents instead of Kummer-Dedekind.
    Ps=s.prime_decomp(2,T=s.Poly(x**3+x*x-2*x+8,x))
    assert [(P.e,P.f) for P in Ps]==[(1,1),(1,1),(1,1)]
    with pytest.raises(NotImplementedError):s.isolate(s.I)
