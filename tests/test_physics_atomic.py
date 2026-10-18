"""Exact angular momentum, operator ordering and elementary quantum systems."""
import itertools
import random
import mpmath as mp
import pytest
import symox as S
from symox.physics import matrices as P
from symox.physics import wigner as W
from symox.physics import hydrogen as H
from symox.physics import qho_1d as Q
from symox.physics import sho as SHO
from symox.physics.paulialgebra import Pauli,evaluate_pauli_product
import sympy as R
from sympy.physics import wigner as RW
from sympy.physics import hydrogen as RH
from sympy.physics import sho as RSHO


def same(a,b):
    a=S.sympify(a)
    b=S.parse_expr(b,local_dict={str(v):v for v in a.free_symbols}) if isinstance(b,str) else S.sympify(b)
    difference=S.trigsimp(S.simplify(a-b))
    if difference!=0:assert R.simplify(R.trigsimp(R.sympify(str(difference),locals={str(v):R.Symbol(str(v),real=True) for v in difference.free_symbols}),method='fu'))==0


def normalized(expression,variable,domain):
    f=R.lambdify(R.Symbol(str(variable)),R.sympify(str(expression)),'mpmath')
    assert float(mp.quad(f,domain))==pytest.approx(1,abs=1e-10)


@pytest.mark.parametrize('a,b,c',[(a,b,c) for a in range(4) for b in range(4) for c in range(abs(a-b),a+b+1)])
def test_wigner_3j_and_clebsch(a,b,c):
    for p,q in [(0,0),(a,-b),(-a,b)]:
        r=-p-q
        ours=W.wigner_3j(a,b,c,p,q,r);theirs=RW.wigner_3j(a,b,c,p,q,r)
        same(ours,str(theirs))
        same(W.clebsch_gordan(a,b,c,p,q,-r),str(RW.clebsch_gordan(a,b,c,p,q,-r)))


@pytest.mark.parametrize('a,b,c',[(S.Rational(1,2),S.Rational(1,2),0),(S.Rational(1,2),1,S.Rational(3,2)),(S.Rational(3,2),2,S.Rational(5,2))])
def test_half_integer_coupling(a,b,c):
    same(W.wigner_3j(a,b,c,a,-b,b-a),str(RW.wigner_3j(*map(lambda x:R.sympify(str(x)),[a,b,c,a,-b,b-a]))))


@pytest.mark.parametrize('indices',[(1,1,1,1,1,1),(2,2,2,2,2,2),(1,2,2,2,1,1),(S.Rational(1,2),S.Rational(1,2),1,S.Rational(1,2),S.Rational(1,2),1)])
def test_six_j(indices):same(W.wigner_6j(*indices),str(RW.wigner_6j(*(R.sympify(str(x)) for x in indices))))


@pytest.mark.parametrize('indices',[(1,1,1,1,1,1,1,1,1),(1,1,0,1,1,0,0,0,0),(2,2,2,2,2,2,2,2,2),(1,1,1,1,1,2,1,2,1)])
def test_nine_j(indices):same(W.wigner_9j(*indices),str(RW.wigner_9j(*indices)))


@pytest.mark.parametrize('degrees,orders',[((0,0,0),(0,0,0)),((1,1,0),(1,-1,0)),((2,2,2),(1,1,-2)),((2,2,2),(-1,-1,2)),((2,2,2),(0,0,0)),((1,1,2),(1,1,0))])
def test_gaunt_and_real_gaunt(degrees,orders):
    same(W.gaunt(*degrees,*orders),str(RW.gaunt(*degrees,*orders)))
    same(W.real_gaunt(*degrees,*orders),str(RW.real_gaunt(*degrees,*orders)))


@pytest.mark.parametrize('spin',[0,S.Rational(1,2),1,S.Rational(3,2),2])
def test_wigner_rotation(spin):
    beta=S.Symbol('beta',real=True);rbeta=R.Symbol('beta',real=True)
    ours=W.wigner_d_small(spin,beta);theirs=RW.wigner_d_small(R.sympify(str(spin)),rbeta)
    for a,b in zip(ours,theirs):same(S.trigsimp(a),str(R.trigsimp(b)))
    for v in ours*ours.T-S.eye(ours.rows):same(S.trigsimp(v),0)


def test_held_wigner_and_gradient():
    same(W.Wigner3j(1,1,0,0,0,0).doit(),-1/S.sqrt(3))
    x=S.Symbol('x');assert W.Wigner3j(1,1,x,0,0,0).doit().free_symbols=={x}
    theta,phi=S.symbols('theta phi')
    value=W.dot_rot_grad_Ynm(3,2,2,0,theta,phi)
    same(value,3*S.sqrt(55)*S.Ynm(5,2,theta,phi)/(11*S.sqrt(S.pi)))
    with pytest.raises(ValueError):W.wigner_3j(S.Rational(1,3),1,1,0,0,0)


@pytest.mark.parametrize('i,j',list(itertools.product(range(1,4),repeat=2)))
def test_pauli_and_matrix_algebra(i,j):
    actual=evaluate_pauli_product(Pauli(i)*Pauli(j))
    matrix=P.msigma(i)*P.msigma(j)
    if i==j:assert actual==1 and matrix==S.eye(2)
    else:
        k=6-i-j;phase=S.I if (i,j) in [(1,2),(2,3),(3,1)] else -S.I
        assert actual==phase*Pauli(k)
        assert matrix==phase*P.msigma(k)
    assert evaluate_pauli_product(S.I*Pauli(1)*Pauli(2)*Pauli(3))==-1
    assert (Pauli(1)+Pauli(2))**2==2


def test_dirac_clifford_fourier_and_parallel_axis():
    metric=[1,-1,-1,-1]
    for i,j in itertools.product(range(4),repeat=2):assert P.mgamma(i)*P.mgamma(j)+P.mgamma(j)*P.mgamma(i)==(2*metric[i]*S.eye(4) if i==j else S.zeros(4))
    assert P.mgamma(5)*P.mgamma(5)==S.eye(4)
    for i in range(4):assert P.mgamma(5)*P.mgamma(i)+P.mgamma(i)*P.mgamma(5)==S.zeros(4)
    assert P.mdft(2)==S.Matrix([[1,1],[1,-1]])/S.sqrt(2)
    assert P.pat_matrix(2,1,0,0)==S.diag(0,2,2)


@pytest.mark.parametrize('n,l',[(n,l) for n in range(1,5) for l in range(n)])
def test_hydrogen_radial_and_energy(n,l):
    r=S.Symbol('r',positive=True);rr=R.Symbol('r',positive=True)
    same(H.R_nl(n,l,r),str(RH.R_nl(n,l,rr)))
    same(H.E_nl(n),str(RH.E_nl(n)))
    assert float(H.E_nl_dirac(n,l).evalf())==pytest.approx(float(RH.E_nl_dirac(n,l)),abs=1e-10)
    normalized(H.R_nl(n,l,r)**2*r**2,r,[0,1,mp.inf])


@pytest.mark.parametrize('n,l',[(0,0),(1,0),(1,1),(2,2)])
def test_spherical_oscillator(n,l):
    r=S.Symbol('r',positive=True);rr=R.Symbol('r',positive=True)
    same(SHO.R_nl(n,l,1,r),str(RSHO.R_nl(n,l,1,rr)))
    normalized(SHO.R_nl(n,l,1,r)**2*r*r,r,[0,1,mp.inf])


def test_one_dimensional_oscillator():
    x=S.Symbol('x',real=True)
    same(Q.E_n(0,2),Q.hbar)
    same(Q.coherent_state(0,0),1)
    for n in range(3):
        state=Q.psi_n(n,x,1,Q.hbar)
        normalized(state**2,x,[-mp.inf,0,mp.inf])
