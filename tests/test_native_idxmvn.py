"""Matrix heads on abstract matrices, differential against SymPy 1.14.

Every expected value below was produced by SymPy 1.14.0 and is hard-coded.
"""
import pytest

from symox import *  # noqa: F401,F403
from symox import (Abs, Adjoint, DiracDelta, Heaviside, I, Identity, ImmutableMatrix, MatAdd,
                   MatrixSymbol, OneMatrix, Symbol, Transpose, ZeroMatrix, arg, conjugate,
                   cos, exp, im, log, re, sin, sqrt, srepr)

del test

A = MatrixSymbol('A', 2, 2)
B = MatrixSymbol('B', 2, 2)
x = Symbol('x')
y = Symbol('y')


def _shape(e):
    """Type names of the tree, a structural view independent of srepr's
    spelling of MatrixSymbol names."""
    if getattr(e, 'is_Atom', False) or type(e).__name__ in ('MatrixSymbol', 'Symbol'):
        return type(e).__name__
    args = getattr(e, 'args', ())
    return (type(e).__name__,) + tuple(_shape(a) for a in args)


# conjugate(M): MatrixExpr._eval_conjugate is Adjoint(Transpose(M)), with the
# kind-specific rules of Transpose, Adjoint, Inverse, MatPow, Kronecker and
# the self-conjugate constants.
CONJUGATE = [
    (lambda: A, 'Adjoint(A.T)', 'Adjoint'),
    (lambda: A.T, 'Adjoint(A)', 'Adjoint'),
    (lambda: Adjoint(A), 'A.T', 'Transpose'),
    (lambda: Adjoint(A.T), 'A', 'MatrixSymbol'),
    (lambda: A*B, 'Adjoint((A*B).T)', 'Adjoint'),
    (lambda: A + B, 'Adjoint((A + B).T)', 'Adjoint'),
    (lambda: 2*A, 'Adjoint((2*A).T)', 'Adjoint'),
    (lambda: x*A, 'Adjoint((x*A).T)', 'Adjoint'),
    (lambda: A.I, 'Adjoint(A.T)**(-1)', 'Inverse'),
    (lambda: A**2, 'Adjoint(A.T)**2', 'MatPow'),
    (lambda: A**x, 'Adjoint(A.T)**x', 'MatPow'),
    (lambda: KroneckerProduct(A, B), 'KroneckerProduct(Adjoint(A.T), Adjoint(B.T))',
     'KroneckerProduct'),
    (lambda: Identity(2), 'I', 'Identity'),
    (lambda: ZeroMatrix(2, 2), '0', 'ZeroMatrix'),
    (lambda: OneMatrix(2, 2), '1', 'OneMatrix'),
    (lambda: HadamardProduct(A, B), 'Adjoint((HadamardProduct(A, B)).T)', 'Adjoint'),
    (lambda: DiagonalOf(A), 'Adjoint(DiagonalOf(A).T)', 'Adjoint'),
]


@pytest.mark.parametrize('build,text,kind', CONJUGATE)
def test_conjugate_of_matrix_expressions(build, text, kind):
    M = build()
    value = conjugate(M)
    assert str(value) == text
    assert type(value).__name__ == kind
    assert str(M.conjugate()) == text


def test_conjugate_structure():
    assert _shape(conjugate(A)) == ('Adjoint', ('Transpose', 'MatrixSymbol'))
    assert conjugate(conjugate(A)) == A
    assert _shape(conjugate(A.I)) == ('Inverse', ('Adjoint', ('Transpose', 'MatrixSymbol')),
                                       'NegativeOne')


# re, im: MatrixExpr.as_real_imag, S.Half*(M + conj) and (M - conj)/(2*I),
# with the ``.doit()`` of MatrixExpr.__add__.
REAL_IMAG = [
    (lambda: A, '(1/2)*Adjoint(A.T) + (1/2)*A', '-(I/2)*(-Adjoint(A.T) + A)'),
    (lambda: x*A, '(x/2)*A + (conjugate(x)/2)*Adjoint(A.T)',
     '-(I/2)*(x*A + (-conjugate(x))*Adjoint(A.T))'),
    (lambda: I*A, '-(I/2)*Adjoint(A.T) + (I/2)*A', '-(I/2)*(I*Adjoint(A.T) + I*A)'),
    (lambda: A**2, '(1/2)*Adjoint(A.T)**2 + (1/2)*A**2', '-(I/2)*(-Adjoint(A.T)**2 + A**2)'),
    (lambda: A.I, '(1/2)*Adjoint(A.T)**(-1) + (1/2)*A**(-1)',
     '-(I/2)*(-Adjoint(A.T)**(-1) + A**(-1))'),
    (lambda: Adjoint(A), '(1/2)*Adjoint(A) + (1/2)*A.T', '-(I/2)*(Adjoint(A) - A.T)'),
    (lambda: Transpose(A*B), '(1/2)*Adjoint(B)*Adjoint(A) + (1/2)*B.T*A.T',
     '-(I/2)*(-Adjoint(B)*Adjoint(A) + B.T*A.T)'),
    (lambda: A*B, '(1/2)*Adjoint(A.T)*Adjoint(B.T) + (1/2)*A*B',
     '-(I/2)*(-Adjoint(A.T)*Adjoint(B.T) + A*B)'),
    (lambda: A + B, '(1/2)*Adjoint(A.T) + (1/2)*Adjoint(B.T) + (1/2)*A + (1/2)*B',
     '-(I/2)*(-Adjoint(A.T) - Adjoint(B.T) + A + B)'),
    (lambda: 2*A, 'Adjoint(A.T) + A', '-(I/2)*(-2*Adjoint(A.T) + 2*A)'),
    (lambda: -A, '-(1/2)*Adjoint(A.T) - (1/2)*A', '-(I/2)*(Adjoint(A.T) - A)'),
    (lambda: Identity(2), 'I', '0'),
    (lambda: ZeroMatrix(2, 2), '0', '0'),
    (lambda: OneMatrix(2, 2), '1', '0'),
]


@pytest.mark.parametrize('build,real,imag', REAL_IMAG)
def test_re_im_of_matrix_expressions(build, real, imag):
    M = build()
    assert str(re(M)) == real
    assert str(im(M)) == imag
    assert tuple(str(part) for part in M.as_real_imag()) == (real, imag)


def test_re_im_types():
    assert type(re(A)).__name__ == 'MatAdd'
    assert type(im(A)).__name__ == 'MatMul'
    assert type(re(Identity(2))).__name__ == 'Identity'
    assert type(im(Identity(2))).__name__ == 'ZeroMatrix'
    assert type(im(OneMatrix(2, 2))).__name__ == 'ZeroMatrix'
    assert str(re(A, evaluate=False)) == 're(A)'


# Abs: Abs.eval of a MatMul takes Abs of each argument; otherwise
# sqrt(expand_mul(M*conjugate(M))), unevaluated when M is its own conjugate.
ABS = [
    (lambda: A, '(A*Adjoint(A.T))**(1/2)', 'MatPow'),
    (lambda: A.T, '(A.T*Adjoint(A))**(1/2)', 'MatPow'),
    (lambda: Adjoint(A), '(Adjoint(A)*A.T)**(1/2)', 'MatPow'),
    (lambda: A**2, '(A**2*Adjoint(A.T)**2)**(1/2)', 'MatPow'),
    (lambda: A**3, '(A**3*Adjoint(A.T)**3)**(1/2)', 'MatPow'),
    (lambda: A.I, '(A**(-1)*Adjoint(A.T)**(-1))**(1/2)', 'MatPow'),
    (lambda: x*A, 'Abs(x)*(A*Adjoint(A.T))**(1/2)', 'MatMul'),
    (lambda: I*A, '(A*Adjoint(A.T))**(1/2)', 'MatPow'),
    (lambda: 2*A, '2*(A*Adjoint(A.T))**(1/2)', 'MatMul'),
    (lambda: -A, '(A*Adjoint(A.T))**(1/2)', 'MatPow'),
    (lambda: -I*x*A.T, 'Abs(x)*(A.T*Adjoint(A))**(1/2)', 'MatMul'),
    (lambda: A*B, '(A*Adjoint(A.T))**(1/2)*(B*Adjoint(B.T))**(1/2)', 'MatMul'),
    (lambda: A*B*x, 'Abs(x)*(A*Adjoint(A.T))**(1/2)*(B*Adjoint(B.T))**(1/2)', 'MatMul'),
    (lambda: A.I*B, '(A**(-1)*Adjoint(A.T)**(-1))**(1/2)*(B*Adjoint(B.T))**(1/2)', 'MatMul'),
    (lambda: Identity(2), 'Abs(I)', 'Abs'),
    (lambda: ZeroMatrix(2, 2), 'Abs(0)', 'Abs'),
    (lambda: OneMatrix(2, 2), 'Abs(1)', 'Abs'),
]


@pytest.mark.parametrize('build,text,kind', ABS)
def test_abs_of_matrix_expressions(build, text, kind):
    value = Abs(build())
    assert str(value) == text
    assert type(value).__name__ == kind


def test_abs_structure_and_explicit_factor():
    assert _shape(Abs(A)) == ('MatPow', ('MatMul', 'MatrixSymbol',
                                         ('Adjoint', ('Transpose', 'MatrixSymbol'))), 'Half')
    value = Abs(x*y*A)
    assert str(value) == 'Abs(x*y)*(A*Adjoint(A.T))**(1/2)'
    value = Abs(A*ImmutableMatrix([[1, 2], [3, 4]]))
    assert type(value).__name__ == 'MatMul'
    assert str(value.args[0]) == '(A*Adjoint(A.T))**(1/2)'
    assert value.args[1] == ImmutableMatrix([[1, 2], [3, 4]])
    assert str(Abs(A, evaluate=False)) == 'Abs(A)'


def test_heaviside_and_diracdelta_of_a_matrix():
    with pytest.raises(ValueError) as error:
        Heaviside(A)
    assert str(error.value) == ('Function defined only for Real Values. Complex part: '
                                '-(I/2)*(-Adjoint(A.T) + A)  found in A .')
    with pytest.raises(ValueError) as error:
        DiracDelta(A)
    assert str(error.value) == ('\nFunction defined only for Real Values. Complex part:\n'
                                '-(I/2)*(-Adjoint(A.T) + A)  found in A .')
    with pytest.raises(ValueError) as error:
        Heaviside(Identity(2))
    assert str(error.value) == 'Function defined only for Real Values. Complex part: 0  found in I .'
    with pytest.raises(ValueError) as error:
        Heaviside(A**2)
    assert str(error.value) == ('Function defined only for Real Values. Complex part: '
                                '-(I/2)*(-Adjoint(A.T)**2 + A**2)  found in A**2 .')
    assert str(Heaviside(A, evaluate=False)) == 'Heaviside(A)'


def test_re_im_of_functions_of_a_matrix():
    assert str(re(log(A))) == 'log((A*Adjoint(A.T))**(1/2))'
    assert str(im(log(A))) == 'arg(A)'
    assert str(re(log(A.T))) == 'log((A.T*Adjoint(A))**(1/2))'
    assert str(im(log(A.T))) == 'arg(A.T)'
    assert str(re(arg(A))) == 'arg(A)'
    assert im(arg(A)) == 0
    for head in (exp,):
        for part in (re, im):
            with pytest.raises(TypeError, match=r'MatAdd\.expand\(\) takes 1 positional argument'):
                part(head(A))
    # cosh is left out: its reference result depends on the SymPy cache state.
    for head in (sec, csc, tanh, coth, acosh, atanh):
        for part in (re, im):
            with pytest.raises(TypeError, match='Mix of Matrix and Scalar symbols'):
                part(head(A))


def test_scalar_matrix_mix_errors():
    for build in (lambda: x + A, lambda: A + x, lambda: sin(A) + A, lambda: A + sin(A),
                  lambda: A - sin(A), lambda: sin(A) - A, lambda: 2 + A, lambda: A + 2,
                  lambda: A - 2, lambda: 2 - A, lambda: MatAdd(sin(A), A), lambda: MatAdd(x, A),
                  lambda: x*A + sin(A), lambda: A**2 + sin(A)):
        with pytest.raises(TypeError) as error:
            build()
        assert str(error.value) == 'Mix of Matrix and Scalar symbols'
    assert str(sin(A) + x) == 'x + sin(A)'


def test_matadd_order_inside_functions():
    C = MatrixSymbol('C', 3, 3)
    assert str(sin(C - C.T)) == 'sin(C - C.T)'
    assert str(sin(2*A - B)) == 'sin(2*A - B)'
    assert str(cos(B + A)) == 'cos(A + B)'
    assert str(C - C.T) == 'C - C.T'
    assert _shape(C - C.T) == ('MatAdd', ('MatMul', 'NegativeOne', ('Transpose', 'MatrixSymbol')),
                               'MatrixSymbol')
    assert srepr(sin(C - C.T)).replace('"C"', "Str('C')") == (
        "sin(MatAdd(MatrixSymbol(Str('C'), Integer(3), Integer(3)), "
        "MatMul(Integer(-1), Transpose(MatrixSymbol(Str('C'), Integer(3), Integer(3))))))")
