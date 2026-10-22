"""Wave-4 polys domains and number fields: values hard-coded from SymPy 1.14."""
import pytest

import symox as S
from symox import (QQ, ZZ, GF, FF, QQ_I, ZZ_I, Poly, sqrt, factor, cancel, gcd,
                   round_two, prime_decomp, prime_valuation, cyclotomic_poly)
from symox.abc import x
from symox.polys.matrices import DomainMatrix
from symox.polys.numberfields.modules import PowerBasis, find_min_poly


def col(*c):
    return DomainMatrix([[ZZ(v)] for v in c], (len(c), 1), ZZ)


# Finite fields

def test_finite_field_elements_print_and_compare_like_modular_integers():
    F = GF(5)
    a = F(3)
    assert repr(a) == 'SymmetricModularIntegerMod5(3)'
    assert str(a) == '3 mod 5'
    assert repr(F(-2)) == 'SymmetricModularIntegerMod5(3)'
    assert F.to_sympy(a) == -2
    assert a == 3 and a == 8 and a == F(-2) and not (a != 3) and a == S.Integer(3)
    assert (F.mod, F.dom, F.sym) == (5, ZZ, True)
    assert repr(FF(7)(10)) == 'SymmetricModularIntegerMod7(3)'


@pytest.mark.parametrize('domain', [GF(2), 'GF(2)', FF(2), 'FF(2)'])
def test_finite_field_domain_option_is_the_modulus(domain):
    assert factor(x**2 + 1, domain=domain) == (x + 1)**2
    assert cancel((x**2 + 1)/(x + 1), domain=domain) == x + 1
    assert gcd(x**2 + 1, x + 1, domain=domain) == x + 1


@pytest.mark.parametrize('expr,p,expected', [
    ((x**2 + 1)/(x + 1), 2, x + 1),
    ((x**2 - 1)/(x + 1), 3, x - 1),
    ((x**3 + 1)/(x**2 + 1), 2, (x**2 + x + 1)/(x + 1)),
    ((x**2 + 1)/(x**3 + x), 2, 1/x),
    ((2*x**2 + 1)/(x + 1), 3, 1 - x),
])
def test_cancel_modulo_p(expr, p, expected):
    assert cancel(expr, modulus=p) == expected


def test_prime_power_finite_fields_are_the_ring_modulo_n():
    # SymPy 1.14 builds the ring Z/9Z for GF(9), which is not a field.
    assert GF(9).is_Field is False
    assert repr(GF(9)(3)*GF(9)(3)) == 'SymmetricModularIntegerMod9(0)'


# Gaussian domains

def test_gaussian_integer_division_lands_in_gaussian_rationals():
    z1, z2 = ZZ_I(2, 3), ZZ_I(1, 3)
    assert repr(ZZ_I(5, 0)/ZZ_I(1, 3)) == 'QQ_I(1/2, -3/2)'
    assert repr(ZZ_I(2, 0)/ZZ_I(1, 1)) == 'QQ_I(1, -1)'
    assert repr(ZZ_I(5, 0)/2) == 'QQ_I(5/2, 0)'
    assert repr(2/ZZ_I(1, 3)) == 'QQ_I(1/5, -3/5)'
    assert repr(ZZ_I(5, 0)//ZZ_I(1, 3)) == 'ZZ_I(1, -1)'
    assert repr(ZZ_I(5, 0) % ZZ_I(1, 3)) == 'ZZ_I(1, -2)'
    assert repr(divmod(ZZ_I(5, 0), ZZ_I(1, 3))) == '(ZZ_I(1, -1), ZZ_I(1, -2))'
    assert repr(ZZ_I(5, 0)//ZZ_I(0, 2)) == 'ZZ_I(0, -2)'
    assert repr(7//ZZ_I(1, 3)) == 'ZZ_I(1, -2)'
    assert repr(7 % ZZ_I(1, 3)) == 'ZZ_I(0, -1)'
    assert (z1//z2)*z2 + z1 % z2 == z1
    with pytest.raises(S.ExactQuotientFailed, match=r'\(1 \+ 3\*I\) does not divide \(5 \+ 0\*I\) in ZZ_I'):
        ZZ_I.exquo(ZZ_I(5, 0), ZZ_I(1, 3))
    with pytest.raises(S.ExactQuotientFailed, match=r'\(1 \+ -3\*I\) does not divide'):
        ZZ_I.exquo(ZZ_I(5, 0), ZZ_I(1, -3))
    with pytest.raises(ZeroDivisionError):
        ZZ_I(5, 1)/ZZ_I(0, 0)


def test_gaussian_rational_floor_division_is_exact():
    z1, z2 = QQ_I(5, 1), QQ_I(2, 1)
    assert repr(z1//z2) == 'QQ_I(11/5, -3/5)'
    assert repr(z1 % z2) == 'QQ_I(0, 0)'
    assert repr(QQ_I(1, 2)//QQ_I(0, 1)) == 'QQ_I(2, -1)'
    assert (z1//z2)*z2 + z1 % z2 == z1


# Algebraic fields

def test_algebraic_field_dense_construction_and_orig_ext():
    K = QQ.algebraic_field(sqrt(2))
    assert K([1, 2]) == K.convert(sqrt(2) + 2)
    assert K.orig_ext == (sqrt(2),)
    assert QQ.algebraic_field(sqrt(2), sqrt(3)).orig_ext == (sqrt(2), sqrt(3))
    assert str(QQ.cyclotomic_field(5).orig_ext[0]) == 'zeta'
    assert K.to_alg_num(K([1, 2])).as_expr() == sqrt(2) + 2


def test_integral_bases_and_discriminants():
    K = QQ.algebraic_field(sqrt(5))
    assert K.integral_basis(fmt='sympy') == [1, S.Rational(1, 2) + sqrt(5)/2]
    assert K.discriminant() == 5
    assert [K.to_sympy(b) for b in K.integral_basis()] == [1, S.Rational(1, 2) + sqrt(5)/2]
    L = QQ.alg_field_from_poly(x**3 + x**2 - 2*x + 8)
    assert L.discriminant() == -503
    t = 'CRootOf(x**3 + x**2 - 2*x + 8, 2)'
    assert str(L.integral_basis(fmt='sympy')) == '[1, %s, %s/2 + %s**2/2]' % (t, t, t)
    Z = QQ.cyclotomic_field(5)
    assert Z.discriminant() == 125


def test_primes_above_common_index_divisor():
    K = QQ.alg_field_from_poly(x**3 + x**2 - 2*x + 8)
    assert str(K.primes_above(2)) == ('[[ (2, _x**2 + 1) e=1, f=1 ], [ (2, (_x**2 + 3*_x + 2)/2) e=1, f=1 ], '
                                      '[ (2, (3*_x**2 + 3*_x)/2) e=1, f=1 ]]')
    assert str(K.primes_above(3)) == '[[ (3, 0) e=1, f=3 ]]'


def test_galois_groups_through_the_field():
    K = QQ.algebraic_field(sqrt(2) + sqrt(3))
    assert K.galois_group(by_name=True)[0].name == 'V'
    # Degree 6 (e.g. QQ(zeta7)) needs resolvents not ported yet.
    with pytest.raises(NotImplementedError):
        QQ.cyclotomic_field(7).galois_group(by_name=True)


# Number field modules, Round Two, primes

def test_round_two_and_prime_decomposition():
    T = Poly(x**3 + x**2 - 2*x + 8)
    ZK, dK = round_two(T)
    assert dK == -503
    assert ZK.QQ_matrix.to_Matrix() == S.Matrix([[1, 0, 0], [0, 1, S.Rational(1, 2)], [0, 0, S.Rational(1, 2)]])
    Ps = prime_decomp(2, T)
    assert [(P.e, P.f) for P in Ps] == [(1, 1), (1, 1), (1, 1)]
    T = Poly(cyclotomic_poly(7, x))
    P = prime_decomp(7, T)
    assert [(Q.e, Q.f) for Q in P] == [(6, 1)]
    Ps = prime_decomp(2, Poly(x**2 + 1))
    assert [(Q.e, Q.f) for Q in Ps] == [(2, 1)]


def test_prime_valuation_reference_and_containment_paths():
    T = Poly(x**2 + 1, x)
    ZK, dK = round_two(T)
    P = prime_decomp(5, ZK=ZK, dK=dK)[0]
    assert prime_valuation(5*ZK, P) == 1
    assert prime_valuation(P**3, P) == 3
    # SymPy 1.14 raises CoercionFailed on these integral ideals (its loop
    # reconverts a non-integral matrix to ZZ); symox returns the valuation.
    ZK, dK = round_two(Poly(x**2 - 2, x))
    P = prime_decomp(2, ZK=ZK, dK=dK)[0]
    assert prime_valuation(2*ZK, P) == 2
    ZK, dK = round_two(Poly(x**3 - x - 1, x))
    Ps = prime_decomp(5, ZK=ZK, dK=dK)
    assert [prime_valuation(Q.as_submodule(), Ps[1]) for Q in Ps] == [0, 1]
    assert prime_valuation(Ps[1]**3, Ps[1]) == 3


def test_module_element_from_rational_and_represent():
    T = Poly(cyclotomic_poly(5, x))
    A = PowerBasis(T)
    B = A.submodule_from_matrix(2 * DomainMatrix.eye(4, ZZ))
    a = A.element_from_rational(QQ(2, 3))
    assert str(a) == '[2, 0, 0, 0]/3'
    b = B.element_from_rational(QQ(2, 3))
    assert str(b) == '[2, 0, 0, 0]/3'
    e = A(col(4, 6, 8, 10))
    assert B.represent(e).to_Matrix() == S.Matrix([2, 3, 4, 5])
    alpha = A(col(0, 1, 0, 0))
    assert find_min_poly(alpha, QQ, x) == Poly(x**4 + x**3 + x**2 + x + 1, x, domain='QQ')
    assert find_min_poly(alpha**2, FF(5), x) == Poly(x**4 + x**3 + x**2 + x + 1, x, modulus=5)
