"""mat-rest6 parity: DFT/IDFT and generic identity nodes, Matrix == RawMatrix,
held matrix derivatives. Expected values come from SymPy 1.14."""
import warnings

import pytest

from symox import *
del test
from symox.matrices.expressions.fourier import DFT, IDFT
from symox.matrices.expressions.special import GenericIdentity, GenericZeroMatrix

n = Symbol('n', integer=True, positive=True)
k = Symbol('k')
A = MatrixSymbol('A', n, n)
X = MatrixSymbol('X', k, k)
B = MatrixSymbol('B', 2, 3)


def test_generic_identities():
    G, Z = GenericIdentity(), GenericZeroMatrix()
    assert srepr(MatMul()) == 'GenericIdentity()'
    assert srepr(MatAdd()) == 'GenericZeroMatrix()'
    assert srepr(G**-1) == 'GenericIdentity()'
    assert srepr(G.T) == 'GenericIdentity()'
    assert MatMul(G, G) == 1
    assert str(G) == 'I' and str(Z) == '0'
    assert latex(G) == r'\mathbb{I}' and latex(Z) == '0'
    assert (G == Identity(n)) is False
    assert G.free_symbols == set()
    assert isinstance(G, Identity) and isinstance(Z, ZeroMatrix)
    assert hash(G) == hash(GenericIdentity())
    assert MatMul.identity == G and MatAdd.identity == Z
    assert G * A == A and A + Z == A
    assert 2 * G == 2 and G * 2 == 2 and G * G == 1
    assert (2 * G).is_Integer and isinstance(2 * G, Integer)
    x = Symbol('x')
    assert x * G == x and G * x == x and G @ A == A and A @ G == A
    with pytest.raises(TypeError, match='does not have a specified shape'):
        G.shape
    with pytest.raises(TypeError):
        G[0, 0]
    with pytest.raises(ValueError, match='not enough values to unpack'):
        MatAdd(Z, Z)


def test_dft():
    assert srepr(DFT(n) * IDFT(n)) == "Identity(Symbol('n', integer=True, positive=True))"
    assert srepr(DFT(n).I) == "IDFT(Symbol('n', integer=True, positive=True))"
    assert srepr(IDFT(n).I) == "DFT(Symbol('n', integer=True, positive=True))"
    assert srepr(DFT(n).conjugate()) == "Adjoint(Transpose(DFT(Symbol('n', integer=True, positive=True))))"
    assert srepr(DFT(n).T) == "Transpose(DFT(Symbol('n', integer=True, positive=True)))"
    assert srepr(DFT(n).adjoint()) == "Adjoint(DFT(Symbol('n', integer=True, positive=True)))"
    r = sqrt(3) / 3
    assert DFT(3).as_explicit() == Matrix([
        [r, r, r],
        [r, sqrt(3) * exp(-2 * I * pi / 3) / 3, sqrt(3) * exp(2 * I * pi / 3) / 3],
        [r, sqrt(3) * exp(2 * I * pi / 3) / 3, sqrt(3) * exp(-2 * I * pi / 3) / 3]])
    assert IDFT(2)[1, 1] == -sqrt(2) / 2
    assert latex(DFT(n)) == r'\text{DFT}_{n}'
    assert latex(IDFT(3)) == r'\text{IDFT}_{3}'
    assert ask(Q.unitary(DFT(n))) is True
    assert DFT(n).xreplace({n: 3}) == DFT(3)
    assert DFT(n).subs(n, 3) == DFT(3)
    assert isinstance(IDFT(n), DFT)
    assert DFT(n).func is DFT and IDFT(n).func is IDFT
    assert DFT(n).args == (n,)
    assert (DFT(n) == IDFT(n)) is False
    assert (DFT(n) * A).args == (DFT(n), A)


def test_matrix_eq_rawmatrix():
    # test_smith_normal_deprecated: Matrix.__eq__ returns NotImplemented.
    from symox.matrices.normalforms import smith_normal_form
    from symox.polys.domains import ZZ
    from symox.polys.solvers import RawMatrix
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        m = RawMatrix([[2, 4]])
        smf = RawMatrix([[2, 0]])
        raw = RawMatrix([[1, 2]])
    m.ring = ZZ
    assert smith_normal_form(m) == smf
    assert (Matrix([[1, 2]]) == raw) is True
    assert (raw == Matrix([[1, 2]])) is True
    assert (Matrix([[1, 2]]) == 3) is False
    assert (Matrix([[1, 2]]) == [[1, 2]]) is False


def test_held_matrix_derivative():
    assert str(X.diff(B)) == 'Derivative(X, B)'
    assert X.diff(B).shape == (2, 3, k, k)
    assert str(X.diff(B, B)) == 'Derivative(X, (B, 2))'
    assert str((X * X).diff(B)) == 'Derivative(X**2, B)'
    assert isinstance(X.diff(B), Derivative)
    assert latex(X.diff(B)) == r'\frac{d}{d B} X'
    assert X.diff(B).expr == X
    assert X.diff(B).variables == (B,)
    assert str(Derivative(X, B)) == 'Derivative(X, B)'
    assert Derivative(X, B).doit() == ZeroMatrix(k, k)
    assert Derivative(X, B, evaluate=True) == ZeroMatrix(k, k)
    assert Derivative(X, B).free_symbols == {X}
    assert (X.diff(B) == Derivative(X, B)) is False
