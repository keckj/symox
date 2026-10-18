"""Exact arithmetic and domain protocols, compared with SymPy where defined."""
from fractions import Fraction
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('rational', [False, True])
@pytest.mark.parametrize('p,q', [(0, 1), (1, 1), (-1, 1), (4, 9), (9, 4),
                               (5, 7), (10**80, 49), (-4, 9)])
def test_ground_sign_and_square_roots(rational, p, q):
    if not rational:
        q = 1
    k, ref = (sx.QQ, sp.QQ) if rational else (sx.ZZ, sp.ZZ)
    a, b = k(p, q), ref(p, q) if rational else ref(p)
    for name in ['is_positive', 'is_negative', 'is_nonpositive', 'is_nonnegative',
                 'is_square', 'canonical_unit', 'abs', 'pos', 'real', 'imag']:
        actual, expected = getattr(k, name)(a), getattr(ref, name)(b)
        if isinstance(expected, bool):
            assert actual is expected
        else:
            assert k.to_sympy(actual).to_sympy() == ref.to_sympy(expected)
    actual, expected = k.exsqrt(a), ref.exsqrt(b)
    assert (actual is None) == (expected is None)
    if actual is not None:
        assert actual.parent() == k and actual*actual == a
        assert k.to_sympy(actual).to_sympy() == ref.to_sympy(expected)


@pytest.mark.parametrize('prime', [2, 3, 5, 7, 41, 101, 65537])
def test_finite_square_roots_and_sign_convention(prime):
    k, ref = sx.GF(prime), sp.GF(prime)
    rng = random.Random(prime)
    for n in [0, 1, prime-1] + [rng.randrange(prime) for _ in range(20)]:
        a, b = k(n), ref(n)
        root, expected = k.exsqrt(a), ref.exsqrt(b)
        assert k.is_square(a) == ref.is_square(b)
        assert (root is None) == (expected is None)
        if root is not None:
            assert root*root == a and int(root) == int(expected)
        for name in ['is_positive', 'is_negative', 'is_nonpositive', 'is_nonnegative']:
            assert getattr(k, name)(a) == getattr(ref, name)(b)


@pytest.mark.parametrize('rational', [False, True])
def test_bezout_cofactors_and_reversible_elements(rational):
    k = sx.QQ if rational else sx.ZZ
    rng = random.Random(193)
    for _ in range(60):
        a, b = k(rng.randrange(-1000, 1001)), k(rng.randrange(1, 1001))
        g, u, v = k.cofactors(a, b)
        assert g*u == a and g*v == b
        s, h = k.half_gcdex(a, b)
        assert h == k.gcd(a, b)
        assert k.rem(h-s*a, b) == k.zero
        if k.is_unit(a):
            assert k.revert(a)*a == k.one
        else:
            with pytest.raises(sx.NotReversible):
                k.revert(a)


@pytest.mark.parametrize('kind', ['ZZring', 'QQring', 'field'])
def test_polynomial_signs_units_and_squares(kind):
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    if kind == 'ZZring':
        k, ref = sx.ZZ.poly_ring(x, y), sp.ZZ.poly_ring(X, Y)
    elif kind == 'QQring':
        k, ref = sx.QQ.poly_ring(x, y), sp.QQ.poly_ring(X, Y)
    else:
        k, ref = sx.QQ.frac_field(x, y), sp.QQ.frac_field(X, Y)
    for expression in [0, -3, x-y, -y*x*x+3*x, -(x+y)**2]:
        a = k(expression)
        b = ref.from_sympy(sx.sympify(expression).to_sympy())
        for name in ['is_positive', 'is_negative', 'canonical_unit']:
            got, expected = getattr(k, name)(a), getattr(ref, name)(b)
            if isinstance(expected, bool):
                assert got == expected
            else:
                assert k.to_sympy(got).to_sympy() == ref.to_sympy(expected)
    square = k((x-2*y)**2 * (x+y)**4)
    root = k.exsqrt(square)
    assert root.parent() == k and root*root == square
    assert k.exsqrt(k(x+y)) is None
    if k.is_Field:
        square /= k((x-y)**6)
        root = k.exsqrt(square)
        assert root*root == square


def test_generator_operations_and_unification():
    x, y, z = sx.symbols('x y z')
    for ground in [sx.ZZ, sx.QQ]:
        assert ground.drop(x) == ground
        assert ground.inject(x) == ground.poly_ring(x)
        for factory in [ground.poly_ring, ground.frac_field]:
            k = factory(x, y, order='grevlex')
            assert k.drop(x) == factory(y, order='grevlex')
            assert k.drop(x, y) == ground
            assert k.inject(z) == factory(x, y, z, order='grevlex')
            assert k.unify_composite(sx.QQ) == k.unify(sx.QQ)
            assert k.unify_with_symbols(sx.QQ, [z]) == k.unify(sx.QQ)
            with pytest.raises(sx.UnificationFailed):
                k.unify_with_symbols(sx.QQ, [x])
            with pytest.raises(sx.GeneratorsError):
                k.inject(x)


def test_conversion_hooks_and_nested_sequences():
    x = sx.Symbol('x')
    k = sx.QQ.algebraic_field(sx.sqrt(2))
    assert sx.QQ.from_AlgebraicField(k(3), k) == sx.QQ(3)
    assert sx.QQ.from_AlgebraicField(k.unit, k) is None
    assert sx.QQ.from_FF(sx.GF(7)(3), sx.GF(7)) is None
    assert sx.ZZ.from_QQ(sx.QQ(1, 2), sx.QQ) is None
    assert sx.QQ.from_ZZ(sx.ZZ(3), sx.ZZ) == sx.QQ(3)
    assert sx.QQ.from_PolynomialRing(sx.QQ[x](3), sx.QQ[x]) == sx.QQ(3)
    assert sx.QQ.from_PolynomialRing(sx.QQ[x](x), sx.QQ[x]) is None
    with pytest.raises(sx.CoercionFailed):
        sx.ZZ.convert_from(sx.QQ(1, 2), sx.QQ)
    with pytest.raises(sx.CoercionFailed):
        sx.QQ.convert_from(sx.ZZ(3), sx.GF(5))
    assert sx.QQ.map([1, [Fraction(2, 3), -4]]) == [sx.QQ(1), [sx.QQ(2, 3), sx.QQ(-4)]]
    assert sx.QQ.sum(Fraction(n, 6) for n in range(7)) == sx.QQ(7, 2)
    assert sx.QQ.sum([]) == sx.QQ.zero
    for name in ['new', 'normal', 'tp']:
        assert getattr(sx.QQ, name)(2, 3) == sx.QQ(2, 3)
    assert sx.QQ.almosteq(sx.QQ(1, 3), sx.QQ(1, 3))
    assert not sx.QQ.almosteq(sx.QQ(1, 3), sx.QQ(1, 2), tolerance=1)
    assert sx.QQ.get_exact() is sx.QQ
    assert sx.QQ.evalf(sx.QQ(1, 3), 30) == sx.Rational(1, 3).evalf(30)


def test_integer_sqrt_and_exact_large_logarithm():
    for n in [0, 1, 2, 99, 10**200+1]:
        root = sx.ZZ.sqrt(sx.ZZ(n))
        assert int(root)**2 <= n < (int(root)+1)**2
    for base in [2, 3, 7, 10, 10**20]:
        for delta in [-1, 0, 1]:
            n = base**100+delta
            result = int(sx.ZZ.log(sx.ZZ(n), sx.ZZ(base)))
            assert base**result <= n < base**(result+1)
    with pytest.raises(ValueError):
        sx.ZZ.log(2, 1)
    with pytest.raises(ValueError):
        sx.ZZ.sqrt(-1)
    with pytest.raises(NotImplementedError):
        sx.QQ.sqrt(sx.QQ(2))


def test_domain_protocols_without_sympy():
    code = '''
import sys
sys.modules['sympy'] = None
import symox as s
x=s.Symbol('x')
assert s.QQ.exsqrt(s.QQ(49,81)) == s.QQ(7,9)
assert s.GF(41).exsqrt(s.GF(41)(5))**2 == s.GF(41)(5)
k=s.QQ.poly_ring(x)
assert k.exsqrt(k((x-2)**4))**2 == k((x-2)**4)
assert s.QQ.from_ZZ(s.ZZ(7),s.ZZ) == s.QQ(7)
assert s.ZZ.log(2**200,2) == s.ZZ(200)
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
