"""SymPy 1.14 documented "c" objects (wave 4 final): expected values from SymPy."""
from symox import (Dict, FiniteSet, KroneckerProduct, MatrixSymbol, Poly, QQ,
                   Set, symbols, nonlinsolve, EmptySet)
from symox.abc import theta, x as X
from symox.utilities.iterables import has_dups


def test_has_dups_set_and_dict_shortcut():
    assert all(has_dups(c) is False for c in (set(), Set(), dict(), Dict()))
    assert has_dups([1, 1]) is True
    assert has_dups(FiniteSet(1, 2)) is False
    assert str(Set()) == 'Set()'
    assert Set() == Set()
    assert Set().args == ()


def test_nonlinsolve_float_coefficients():
    x, y = symbols('x y')
    assert str(nonlinsolve([x**2 - 0.25, y - x], [x, y])) == '{(-0.5, -0.5), (0.5, 0.5)}'
    assert str(nonlinsolve([x**2 - 0.2, y - x*1.5], [x, y])) == \
        '{(-0.447213595499958, -0.670820393249937), (0.447213595499958, 0.670820393249937)}'
    assert nonlinsolve([x**2 + y**2 - 0.0625, 2*x - 1.0, 2*y], [x, y]) == EmptySet
    assert str(nonlinsolve([x*y - 0.5, x - y], [x, y])) == \
        '{(-0.707106781186548, -0.707106781186548), (0.707106781186548, 0.707106781186548)}'


def _gamma_setup():
    from symox.physics.hep.gamma_matrices import GammaMatrix as G, LorentzIndex
    from symox.tensor.tensor import tensor_indices, tensor_heads
    p, q = tensor_heads('p, q', [LorentzIndex])
    i0, i1 = tensor_indices('i0:2', LorentzIndex)
    return G, p, q, i0, i1


def test_simplify_gpgp_reference_examples():
    from symox.physics.hep.gamma_matrices import simplify_gpgp
    G, p, q, i0, i1 = _gamma_setup()
    ps = p(i0)*G(-i0)
    qs = q(i0)*G(-i0)
    assert str(simplify_gpgp(ps*qs*qs)) == 'GammaMatrix(-L_0)*p(L_0)*q(L_1)*q(-L_1)'
    assert str(simplify_gpgp(ps*qs)) == 'GammaMatrix(-L_0)*GammaMatrix(-L_1)*p(L_0)*q(L_1)'
    assert str(simplify_gpgp(ps*ps*qs*qs*ps)) == \
        'GammaMatrix(-L_0)*p(L_0)*p(L_1)*p(-L_1)*q(L_2)*q(-L_2)'


def test_simplify_gpgp_lone_pair_is_not_traced():
    # Intentional divergence: the reference returns 4*p(L_0)*p(-L_0), the
    # trace, although p_slash*p_slash is p.p times the identity.
    from symox.physics.hep.gamma_matrices import simplify_gpgp
    G, p, q, i0, i1 = _gamma_setup()
    ps = p(i0)*G(-i0)
    qs = q(i0)*G(-i0)
    assert str(simplify_gpgp(ps*ps)) == 'p(L_0)*p(-L_0)'
    assert str(simplify_gpgp(qs*ps*ps*qs)) == 'p(L_0)*p(-L_0)*q(L_1)*q(-L_1)'


def test_prime_decomp_symbolic_generator():
    T = Poly(X**3 + X**2 - 2*X + 8)
    K = QQ.algebraic_field((T, theta))
    assert str(K) == 'QQ<theta>'
    assert K.orig_ext == (theta,)
    assert str(K.primes_above(2)) == ('[[ (2, x**2 + 1) e=1, f=1 ], [ (2, (x**2 + 3*x + 2)/2) e=1, f=1 ], '
                                      '[ (2, (3*x**2 + 3*x)/2) e=1, f=1 ]]')
    assert str(K.primes_above(3)) == '[[ (3, 0) e=1, f=3 ]]'
    assert K.discriminant() == -503
    assert str(K.galois_group(by_name=True)[0]) == 'S3TransitiveSubgroups.S3'


def test_prime_decomp_symbolic_generator_reducible_fails_loudly():
    import pytest
    with pytest.raises(NotImplementedError):
        QQ.algebraic_field((Poly(X**2 - 1), theta))


def test_combine_kronecker_sum_is_not_factored():
    # The reference returns KroneckerProduct(A + B.T, B + A.T), which drops
    # the cross terms A x A.T + B.T x B (1x1 example: 2 != 4).
    from symox.matrices.expressions import combine_kronecker
    m, n = symbols('m, n', integer=True)
    A = MatrixSymbol('A', m, n)
    B = MatrixSymbol('B', n, m)
    assert str(combine_kronecker(KroneckerProduct(A, B)*KroneckerProduct(B, A))) == 'KroneckerProduct(A*B, B*A)'
    e = KroneckerProduct(A, B) + KroneckerProduct(B.T, A.T)
    assert combine_kronecker(e) == e
