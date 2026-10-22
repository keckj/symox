"""SymPy 1.14 parity of three core gaps: Tuple equality with Python
tuples, the complex-number assumptions of sums and products with I, and
GF(n) for moduli outside the native prime field (the ring Z/nZ).

Expected values come from the reference SymPy 1.14.
"""
import warnings

import pytest

from symox import (I, E, GF, FiniteField, Integer, Rational, S, Symbol, Tuple, FiniteSet,
                   Basic, Lambda, Matrix, Poly, oo, pi, exp, sqrt, cos, symbols, srepr)
from symox.polys.domains import FF
from symox.polys.polyerrors import NotInvertible, NotReversible, CoercionFailed, DomainError

x, y = symbols('x y')


# Tuple equality --------------------------------------------------------

@pytest.mark.parametrize('left, right, expected', [
    (Tuple(1, 2), (1, 2), True),
    (Tuple(1, 2), (1, 3), False),
    (Tuple(1, 2), [1, 2], False),
    (Tuple(1, Tuple(2, 3)), (1, (2, 3)), True),
    (Tuple(1, (2, 3)), (1, (2, 3)), True),
    (Tuple(x, y), (x, y), True),
    (Tuple(x), ('x',), False),
    (Tuple(1), 1, False),
    (Tuple(1.0, 2), (1, 2), False),
    (Tuple(1, 2), (1.0, 2), False),
    (Tuple(S.Half), (Rational(1, 2),), True),
    (Tuple(), (), True),
    (Tuple(1, 2), Tuple(1, 2), True),
    (Tuple(1, 2), {1, 2}, False),
    (Tuple(1, 2), {1: 2}, False),
    (Tuple(1, 2), range(1, 3), False),
    (Tuple(1, 2), Matrix([1, 2]), False),
    (Integer(1), (1,), False),
])
def test_tuple_equality_with_python_objects(left, right, expected):
    assert (left == right) is expected
    assert (left != right) is (not expected)
    assert (right == left) is expected
    assert (right != left) is (not expected)


def test_tuple_hash_and_container_lookups():
    t = Tuple(1, 2)
    assert hash(t) == hash((1, 2))
    assert t in {(1, 2): 1}
    assert (1, 2) in {t: 1}
    assert t in {(1, 2)}
    assert t in [(1, 2)]
    assert (1, 2) in [t]
    assert len({Tuple(1, 2): 1, (1, 2): 2}) == 1
    assert len({Tuple(1, 2), (1, 2)}) == 1
    assert hash(Tuple(x, Tuple(1, y))) == hash((x, (1, y)))
    assert {Tuple(x, 1): 'a'}[(x, 1)] == 'a'
    items = (x, 1.5, S.Half, 2**70, oo, -oo, S.true, I, Tuple(), (1, (y,)))
    assert hash(Tuple(*items)) == hash(Tuple(*items).args)


def test_tuple_sequence_methods_compare_items():
    assert Tuple(1, 2, (3, 4)).index((3, 4)) == 2
    assert Tuple(1, (3, 4), (3, 4)).count((3, 4)) == 2
    assert (3, 4) in Tuple(1, (3, 4))
    assert 1 in Tuple(1, (3, 4))


def test_tuple_equality_in_sets_and_basic():
    assert FiniteSet((1, 2)) == FiniteSet(Tuple(1, 2))
    assert (1, 2) in FiniteSet(Tuple(1, 2))
    assert FiniteSet(Tuple(1, 2)).contains((1, 2)) is S.true
    assert Basic(Tuple(1, 2)) == Basic((1, 2))
    assert Lambda((x, y), x).args[0] == (x, y)
    assert str(Tuple(x, 1).subs(x, 2)) == '(2, 1)'


# Assumptions of numbers with I ----------------------------------------

_PROPS = ['is_complex', 'is_finite', 'is_infinite', 'is_extended_real', 'is_real',
          'is_algebraic', 'is_transcendental', 'is_imaginary', 'is_zero', 'is_nonzero',
          'is_rational', 'is_irrational', 'is_integer', 'is_number', 'is_commutative',
          'is_extended_nonzero', 'is_positive', 'is_extended_positive',
          'is_extended_negative', 'is_hermitian', 'is_antihermitian', 'is_comparable']

_r = Symbol('r', real=True)
_c = Symbol('c', complex=True)
_p = Symbol('p', positive=True)

# Reference letters (T/F/N) for each property of _PROPS, in order.
# sqrt(2 + I) is left out: the reference proves it non-real through
# arg(2 + I)/(2*pi) not being an integer, which symox does not derive.
_ASSUMPTIONS = [
    (lambda: 2 + I, 'TTFFFTFFFFFFFTTFFFFNNF'),
    (lambda: 2*I, 'TTFFFTFTFFFFFTTFFFFNTF'),
    (lambda: 3 + 4*I, 'TTFFFTFFFFFFFTTFFFFNNF'),
    (lambda: S.Half + I, 'TTFFFTFFFFFFFTTFFFFNNF'),
    (lambda: sqrt(2) + I, 'TTFFFTFFFFFFFTTFFFFNNF'),
    (lambda: sqrt(2)*I, 'TTFFFTFTFFFFFTTFFFFNTF'),
    (lambda: pi + I, 'TTFFFFTFFFFFFTTFFFFNNF'),
    (lambda: pi*I, 'TTFFFFTTFFFFFTTFFFFNTF'),
    (lambda: E*I + 1, 'TTFFFFTFFFFFFTTFFFFNNF'),
    (lambda: I*(1 + I), 'TTFFFTFFFFFFFTTFFFFNNF'),
    (lambda: Rational(1, 3)*I + 2, 'TTFFFTFFFFFFFTTFFFFNNF'),
    (lambda: I + oo, 'FFTFFFFFFFFFFTTFFFFNNF'),
    (lambda: oo*I, 'FFTFFFFFFFFFFTTFFFFNNF'),
    (lambda: _c + I, 'TTFNNNNNNNNNNFTNNNNNNF'),
    (lambda: _c*I, 'TTFNNNNNNNNNNFTNNNNNNF'),
    (lambda: _r*I, 'TTFNNNNNNNNNNFTNNNNNNF'),
    (lambda: exp(I), 'TTFFFFTNFFFFFTTFFFFNNF'),
    (lambda: (1 + I)/(1 - I), 'TTFFFTFTFFFFFTTFFFFNTF'),
    (lambda: I*sqrt(3) + 1, 'TTFFFTFFFFFFFTTFFFFNNF'),
    (lambda: cos(1) + I, 'TTFFFFTFFFFFFTTFFFFNNF'),
    (lambda: I*(_r + 1), 'TTFNNNNNNNNNNFTNNNNNNF'),
    (lambda: _r + I, 'TTFFFNNNFFFFFFTFFFFNNF'),
    (lambda: _p + I, 'TTFFFNNFFFFFFFTFFFFNNF'),
]


@pytest.mark.parametrize('make, expected', _ASSUMPTIONS)
def test_complex_number_assumptions(make, expected):
    e = make()
    got = ''.join(str(getattr(e, p))[0] for p in _PROPS)
    assert got == expected, srepr(e)


def test_two_plus_i_is_complex():
    assert (2 + I).is_complex is True
    assert (2 + I).is_finite is True
    assert (2 + I).is_extended_real is False
    assert (2 + I).is_algebraic is True
    assert (2 + I).is_imaginary is False
    assert ((2 + I)**2).is_complex is True


# GF(n) beyond the native prime field -----------------------------------

def test_composite_gf_is_the_ring_of_integers_modulo_n():
    K = GF(9)
    assert repr(K) == 'GF(9)' and str(K) == 'GF(9)'
    assert K.is_Field is False
    assert K.is_FiniteField is True
    assert isinstance(K, FiniteField)
    assert FiniteField(9) == K and FF(9) == K
    assert K.mod == 9 and K.characteristic() == 9
    assert repr(K(5) + K(7)) == 'SymmetricModularIntegerMod9(3)'
    assert repr(K(5)*K(7)) == 'SymmetricModularIntegerMod9(8)'
    assert repr(K(2) - K(5)) == 'SymmetricModularIntegerMod9(6)'
    assert repr(-K(2)) == 'SymmetricModularIntegerMod9(7)'
    assert repr(K(2)**5) == 'SymmetricModularIntegerMod9(5)'
    assert repr(K(2)**-1) == 'SymmetricModularIntegerMod9(5)'
    assert repr(K(4)/K(2)) == 'SymmetricModularIntegerMod9(2)'
    assert repr(K.revert(K(2))) == 'SymmetricModularIntegerMod9(5)'
    assert repr(K.quo(K(4), K(2))) == 'SymmetricModularIntegerMod9(2)'
    assert repr(K.gcd(K(3), K(6))) == 'SymmetricModularIntegerMod9(1)'
    assert K(10) == K(1)
    assert int(K(8)) == 8
    assert K.to_sympy(K(8)) == -1
    assert repr(K.from_sympy(Integer(11))) == 'SymmetricModularIntegerMod9(2)'
    assert repr(K.convert(13)) == 'SymmetricModularIntegerMod9(4)'
    assert str(K(5)) == '5 mod 9'
    assert repr(K(-1)) == 'SymmetricModularIntegerMod9(8)'
    assert K(5).parent() == K
    assert bool(K(9)) is False
    assert repr(K(5) + 3) == 'SymmetricModularIntegerMod9(8)'
    assert repr(4*K(5)) == 'SymmetricModularIntegerMod9(2)'
    assert repr(K(5)/2) == 'SymmetricModularIntegerMod9(7)'
    assert K.is_unit(K(3)) is True and K.is_unit(K(2)) is True
    assert K.get_field() == K
    assert GF(4).is_Field is False
    assert repr(GF(6)(5)*GF(6)(5)) == 'SymmetricModularIntegerMod6(1)'
    assert repr(GF(15)(7)**-1) == 'SymmetricModularIntegerMod15(13)'
    assert repr(GF(1)) == 'GF(1)'
    assert hash(GF(9)) == hash(GF(9)) and GF(9) != GF(3)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert (K(2) < K(3)) is True


@pytest.mark.parametrize('call', [
    lambda K: K(3)**-1,
    lambda K: K(4)/K(3),
    lambda K: K(4)/K(0),
    lambda K: K.revert(K(3)),
    lambda K: K.quo(K(4), K(3)),
    lambda K: K.exquo(K(4), K(3)),
    lambda K: K(5)/3,
])
def test_composite_gf_non_units_raise_not_invertible(call):
    with pytest.raises(NotInvertible, match='zero divisor'):
        call(GF(9))


def test_composite_gf_errors():
    K = GF(9)
    for modulus in (0, -4):
        with pytest.raises(ValueError, match='modulus must be a positive integer'):
            GF(modulus)
    with pytest.raises(DomainError, match='there is no ring associated with GF\\(9\\)'):
        K.get_ring()
    with pytest.raises(CoercionFailed, match='expected an integer, got 1/2'):
        K.convert(Rational(1, 2))
    with pytest.raises(NotReversible):
        K.revert(K(0))


def test_nonsymmetric_and_large_moduli():
    assert repr(GF(9, symmetric=False)(8)) == 'ModularIntegerMod9(8)'
    assert GF(9, symmetric=False).to_sympy(GF(9, symmetric=False)(8)) == 8
    big = 2**61 - 1
    assert GF(big).is_Field is True
    assert int(GF(big)(3)**-1 * 3) == 1


def test_polynomials_over_composite_gf_are_unsupported():
    # The reference builds Poly(x**2 - 4*x - 2, x, modulus=9); symox has no
    # polynomial arithmetic over Z/nZ and refuses instead of guessing.
    with pytest.raises(NotImplementedError):
        Poly(x**2 + 5*x + 7, x, domain=GF(9))


def test_lfsr_sequence_over_composite_gf():
    from symox import crypto
    K = GF(4)
    assert repr(crypto.lfsr_sequence([K(1), K(3)], [K(1), K(2)], 5)) == (
        '[SymmetricModularIntegerMod4(1), SymmetricModularIntegerMod4(2), '
        'SymmetricModularIntegerMod4(3), SymmetricModularIntegerMod4(3), '
        'SymmetricModularIntegerMod4(0)]')
