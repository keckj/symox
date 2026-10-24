"""simp-rest6 parity regressions; expected values from SymPy 1.14."""
from symox import *
del test
from symox.simplify.simplify import nc_simplify
from symox.simplify.gammasimp import gammasimp as gs2
from symox.polys.domains.realfield import RR as RR2

x,s,y,z,a=symbols('x s y z a')
n=Symbol('n', odd=True); m=Symbol('m', even=True); k=Symbol('k', integer=True)
A,B,C=[MatrixSymbol(q,3,3) for q in 'ABC']
an,bn,cn=symbols('a b c', commutative=False)
X=Symbol('X', commutative=False)
M=MatrixSymbol('M',10,10)
p=Symbol('p', positive=True)


def test_pd1():
    assert str(powdenest(2*4**(s/2-1), polar=True)) == '2**(s - 1)'


def test_pd2():
    assert str(powdenest((x**2)**y, polar=True)) == 'x**(2*y)'


def test_pd3():
    assert str(powdenest(sqrt(x**2), polar=True)) == 'x'


def test_pd4():
    assert str(powdenest(4**(s/2-1)*x, polar=True)) == '2**(s - 2)*x'


def test_pd5():
    assert str(powdenest((2**x)**y, force=True, polar=True)) == '2**(x*y)'


def test_mellin():
    assert str(mellin_transform(besselj(0,x),x,s)) == '(2**(s - 1)*gamma(s/2)/gamma(1 - s/2), (0, 3/2), True)'


def test_sim1():
    assert str(simplify(2**(s-1)*y)) == '2**(s - 1)*y'


def test_sim2():
    assert str(simplify(2*3**(x-1))) == '2*3**(x - 1)'


def test_sim3():
    assert str(simplify(sqrt(4*p*x+4*p*y))) == '2*sqrt(p)*sqrt(x + y)'


def test_pw1():
    assert str(powsimp(3**x*Rational(2, 3))) == '2*3**(x - 1)'


def test_pw2():
    assert str(powsimp(2*3**x/3)) == '2*3**(x - 1)'


def test_pw3():
    assert str(powsimp(2**x*Rational(3,4))) == '3*2**(x - 2)'


def test_pw4():
    assert str(powsimp(2**x/4)) == '2**(x - 2)'


def test_pw5():
    assert str(powsimp(6**x/4)) == '6**x/4'


def test_pwm1():
    assert str(powsimp(M, deep=True) == M) == 'True'


def test_pwm2():
    assert str(type(powsimp(M, deep=True)).__name__) == 'MatrixSymbol'


def test_pwm3():
    assert str(type(powsimp(A**2*A, deep=True)).__name__) == 'MatPow'


def test_nce1():
    assert str(expand(cn**-1*an**-1)) == 'c**(-1)*a**(-1)'


def test_nce2():
    assert str(expand(an**-1*bn*cn**-1)) == 'a**(-1)*b*c**(-1)'


def test_ncs1():
    assert str(nc_simplify(an**-1*bn*cn**-1)) == '(c*b**(-1)*a)**(-1)'


def test_ncr1():
    assert str((2.0*cos(pi/3)*X).simplify()) == 'X'


def test_ncr2():
    assert str((2*cos(pi/3)*X).simplify()) == 'X'


def test_ncm1():
    assert str(nc_simplify(MatMul(A**-1, B, C**-1))) == '(C*B**(-1)*A)**(-1)'


def test_ncm2():
    assert str(nc_simplify(MatMul(B**-1, A**-1, MatPow(MatMul(A, B), 2)))) == 'A*B'


def test_ncm3():
    assert str(MatPow(MatMul(C, B**-1, A), -1).doit(inv_expand=False)) == '(C*B**(-1)*A)**(-1)'


def test_ncm4():
    assert str(MatPow(MatMul(A, B), 2)**-1) == '(A*B)**(-2)'


def test_ncm5():
    assert str(MatPow(MatMul(C, B**-1, A), -1).doit()) == 'A**(-1)*B*C**(-1)'


def test_bs1():
    assert str(besselsimp(-besselj(1,x) - besselj(-1,x))) == '0'


def test_zoo1():
    assert str(zoo*sqrt(5-sqrt(5))) == 'zoo'


def test_zoo2():
    assert str(expand_mul(zoo*(sqrt(5-sqrt(5))*sqrt(40-8*sqrt(5))+1))) == 'zoo'


def test_neg1():
    assert str((-I)**n) == '-I**n'


def test_neg2():
    assert str((-I)**m) == 'I**m'


def test_neg3():
    assert str((-I)**k) == '(-I)**k'


def test_neg4():
    assert str((-2)**n) == '-2**n'


def test_neg5():
    assert str((-sqrt(2))**n) == '-2**(n/2)'


def test_neg6():
    assert str((-2*I)**m) == '(2*I)**m'


def test_neg7():
    assert str(powsimp((-1)**(n*Rational(3, 2)))) == '-I**n'


def test_neg8():
    assert str((-x)**n) == '(-x)**n'


def test_gs():
    assert str(gs2(gamma(x+1)/gamma(x))) == 'x'


def test_rr():
    assert str(RR2(1.5)) == '1.5'
