"""Wave 15 pw-nc: noncommutative Piecewise, partial evalf of -1 and I,
adjoint printing and matrix algorithms on noncommutative entries, checked
against SymPy 1.14 (expected values hard-coded from SymPy)."""
import numpy as np
import pytest

import symox as sp
from symox import *  # noqa: F401,F403
from symox.polys.polyerrors import PolynomialError

del test

A, B, C = symbols('A B C', commutative=False)
x, y, z = symbols('x y z')
a, b = symbols('a b', real=True)
f = Function('f')
NS = dict(globals())


EVALF = [
    ('(-x).evalf()', '-x'),
    ('(-x*y).evalf()', '-x*y'),
    ('(x - y*z).evalf()', 'x - y*z'),
    ('(-2*x*y).evalf()', '-2.0*x*y'),
    ('(-pi*x).evalf()', '-3.14159265358979*x'),
    ('(-x*y).n(5)', '-x*y'),
    ('(-x**2).evalf()', '-x**2'),
    ('(-sin(x)*y).evalf()', '-y*sin(x)'),
    ('(x*y - y*z).evalf()', 'x*y - y*z'),
    ('(-A*B).evalf()', '-A*B'),
    ('(A*B - B*A).evalf()', 'A*B - B*A'),
    ('(-x*A).evalf()', '-x*A'),
    ('(-A).evalf()', '-A'),
    ('(x - A*B).evalf()', 'x - A*B'),
    ('(-x*y + 1).evalf()', '-x*y + 1.0'),
    ('(-x*y*z).evalf()', '-x*y*z'),
    ('(-sqrt(2)*x*y).evalf()', '-1.4142135623731*x*y'),
    ('(-I*x*y).evalf()', '-I*x*y'),
    ('(-x*y/3).evalf()', '-0.333333333333333*x*y'),
    ('Mul(-1, x, y, evaluate=False).evalf()', '-x*y'),
    ('N(-x*y)', '-x*y'),
    ('(I*x*y).evalf()', 'I*x*y'),
    ('(I*x).evalf()', 'I*x'),
    ('(2*I*x).evalf()', '2.0*I*x'),
    ('(-2*I*x).evalf()', '-2.0*I*x'),
    ('(-I*x).evalf()', '-I*x'),
    ('(-I*A).evalf()', '-I*A'),
    ('(-pi*x*A).evalf()', '-3.14159265358979*x*A'),
    ('(-x*sin(1)).evalf()', '-0.841470984807897*x'),
    ('(-x*y*pi).evalf()', '-3.14159265358979*x*y'),
    ('(-(x + 1)).evalf()', '-x - 1.0'),
    ('Mul(-1, x + 1, evaluate=False).evalf()', '-x - 1.0'),
    ('(-x*y).evalf(chop=True)', '-x*y'),
    ('(-sin(x)).evalf()', '-sin(x)'),
    ('(-sin(x)).evalf(subs={x: 1})', '-0.841470984807897'),
    ('(-x - y).evalf()', '-x - y'),
    ('(-oo*x).evalf()', '-oo*x'),
    ('(-zoo*x).evalf()', 'zoo*x'),
    ('(x + I).evalf()', 'x + I'),
    ('(x + 1 + I).evalf()', 'x + 1.0 + I'),
    ('(x + 2*I).evalf()', 'x + 2.0*I'),
    ('(x**I).evalf()', 'x**I'),
    ('sin(I*x).evalf()', 'I*sinh(x)'),
    ('exp(I*x).evalf()', 'exp(I*x)'),
    ('exp(-x).evalf()', 'exp(-x)'),
    ('sin(-x).evalf()', '-sin(x)'),
    ('(x**(-y)).evalf()', 'x**(-y)'),
    ('(2**(-x)).evalf()', '2.0**(-x)'),
    ('(x*(-y + 1)).evalf()', 'x*(1.0 - y)'),
    ('(x + pi*I).evalf()', 'x + 3.14159265358979*I'),
    ('(I*pi*x).evalf()', '3.14159265358979*I*x'),
    ('(-I*pi*x).evalf()', '-3.14159265358979*I*x'),
    ('(-x*exp(-y)).evalf()', '-x*exp(-y)'),
    ('(-x/y).evalf()', '-x/y'),
    ('(x - y/2).evalf()', 'x - 0.5*y'),
    ('(-x**2*y).evalf(3)', '-x**2*y'),
    ('Piecewise((-x, y > 0), (x, True)).evalf()', 'Piecewise((-x, y > 0), (x, True))'),
    ('Eq(-x, y).evalf()', 'Eq(-x, y)'),
    ('(-x).evalf(subs={y: 1})', '-x'),
    ('(-sqrt(x)).evalf()', '-x**0.5'),
    ('log(-x).evalf()', 'log(-x)'),
]

MISC = [
    ('Piecewise((A, y < 1), (B, True)).is_commutative', 'False'),
    ('Piecewise((A, y < 1), (x, True)).is_commutative', 'False'),
    ('(Piecewise((A, y < 1), (B, True))*C).args', '(Piecewise((A, y < 1), (B, True)), C)'),
    ('(C*Piecewise((A, y < 1), (B, True))*A).args', '(C, Piecewise((A, y < 1), (B, True)), A)'),
    ('(A*Piecewise((A, y < 1), (B, True))).args', '(A, Piecewise((A, y < 1), (B, True)))'),
    ('piecewise_fold(C*Piecewise((A, y < 1), (B, True))*A)', 'Piecewise((C*A**2, y < 1), (C*B*A, True))'),
    ('piecewise_fold(Piecewise((A, y < 1), (B, True))*C + Piecewise((A, y < 1), (B, True)))', 'Piecewise((A + A*C, y < 1), (B + B*C, True))'),
    ('piecewise_fold(Piecewise((A, y < 1), (B, True))*Piecewise((B, y < 1), (A, True)))', 'Piecewise((A*B, y < 1), (B*A, True))'),
    ('(C*Piecewise((A, y < 1), (B, True))*A).subs(y, 2)', 'C*B*A'),
    ('expand(C*Piecewise((A, y < 1), (B, True))*(A + B))', 'C*Piecewise((A, y < 1), (B, True))*A + C*Piecewise((A, y < 1), (B, True))*B'),
    ('Piecewise((A, y < 1), (B, True))**2', 'Piecewise((A**2, y < 1), (B**2, True))'),
    ('latex(adjoint(A)**2)', '\\left(A^{\\dagger}\\right)^{2}'),
    ('latex(adjoint(A)**x)', '\\left(A^{\\dagger}\\right)^{x}'),
    ('latex(adjoint(f(A)))', 'f{\\left(A \\right)}^{\\dagger}'),
    ('latex(adjoint(A) + B)', 'B + A^{\\dagger}'),
    ('latex(adjoint(A)*B)', 'A^{\\dagger} B'),
    ('latex(sin(adjoint(A)))', '\\sin{\\left(A^{\\dagger} \\right)}'),
    ('latex(adjoint(A*B))', 'B^{\\dagger} A^{\\dagger}'),
    ('pretty(adjoint(A))', ' †\nA '),
    ('pretty(adjoint(A)**2)', '    2\n⎛ †⎞ \n⎝A ⎠ '),
    ('pretty(adjoint(A*B))', ' †  †\nB ⋅A '),
    ('pretty(adjoint(A) + A)', '     †\nA + A '),
    ('pretty(adjoint(sin(A)))', '         ⎛   ⎛_⎞⎞\ntranspose⎝sin⎝A⎠⎠'),
    ('pretty(transpose(A))', 'transpose(A)'),
    ('pretty(adjoint(A), use_unicode=False)', ' +\nA '),
]



@pytest.mark.parametrize('source, expected', EVALF)
def test_partial_evalf_keeps_minus_one_and_I(source, expected):
    assert str(eval(source, NS)) == expected


@pytest.mark.parametrize('source, expected', MISC)
def test_noncommutative_piecewise_and_adjoint_printing(source, expected):
    assert str(eval(source, NS)) == expected


@pytest.mark.parametrize('source', [
    'Matrix([[A, B], [C, A]]).charpoly()',
    'Matrix([[A, B], [C, A]]).charpoly(x)',
    'Matrix([[A, 0], [0, A]]).eigenvals()',
    'Matrix([[A, 0], [0, A]]).eigenvects()',
    'Matrix([[A, 0], [0, B]]).is_diagonalizable()',
    'Matrix([[A, 0], [0, B]]).diagonalize()',
    'Matrix([[A, 0], [0, B]]).jordan_form()',
    'Matrix([[A, 0], [0, B]]).exp()',
    'Matrix([[A, B], [C, A]]).pinv()',
    'Matrix([[A, B], [C, A]]).singular_values()',
    'Matrix([[A, B], [C, A]]).is_nilpotent()',
    'ImmutableMatrix([[A, B], [C, A]]).charpoly()',
])
def test_polynomial_algorithms_reject_noncommutative_entries(source):
    with pytest.raises(PolynomialError, match='non-commutative expressions are not supported'):
        eval(source, NS)


@pytest.mark.parametrize('source', [
    'Matrix([[A, B], [C, A]]).echelon_form()',
    'Matrix([[A, B], [C, A]]).QRdecomposition()',
])
def test_unsupported_noncommutative_matrix_algorithms(source):
    with pytest.raises(NotImplementedError):
        eval(source, NS)


def _numeric(e, values):
    e = sympify(e)
    if e.is_Symbol:
        return values[e]
    if e.is_Number:
        return float(e)*np.eye(2)
    if e.is_Add:
        return sum(_numeric(t, values) for t in e.args)
    if e.is_Mul:
        r = np.eye(2)
        for t in e.args:
            r = r @ _numeric(t, values)
        return r
    if e.is_Pow:
        base, n = _numeric(e.base, values), int(e.exp)
        return np.linalg.matrix_power(base if n >= 0 else np.linalg.inv(base), abs(n))
    raise ValueError(e)


@pytest.mark.parametrize('method', ['GE', 'ADJ', 'LU'])
def test_noncommutative_inverse_is_the_block_inverse(method):
    # Intentional correction: SymPy's GE and ADJ inverses apply commutative
    # formulas (A*(A**2 - B*C)**-1); the block inverse is the LU one.
    M = Matrix([[A, B], [C, A]])
    inverse = M.inv(method=method)
    assert inverse == M.inv(method='LU')
    rng = np.random.default_rng(1)
    values = {s: rng.standard_normal((2, 2)) for s in (A, B, C)}
    block = np.block([[values[A], values[B]], [values[C], values[A]]])
    numeric = np.block([[_numeric(inverse[i, j], values) for j in range(2)] for i in range(2)])
    assert np.abs(numeric @ block - np.eye(4)).max() < 1e-9


def test_noncommutative_solve_uses_lu():
    M = Matrix([[A, B], [C, A]])
    rhs = Matrix([1, 2])
    assert M.solve(rhs) == M.LUsolve(rhs)
    rng = np.random.default_rng(2)
    values = {s: rng.standard_normal((2, 2)) for s in (A, B, C)}
    block = np.block([[values[A], values[B]], [values[C], values[A]]])
    solution = M.solve(rhs)
    numeric = np.vstack([_numeric(solution[i], values) for i in range(2)])
    target = np.vstack([np.eye(2), 2*np.eye(2)])
    assert np.abs(block @ numeric - target).max() < 1e-9


def test_commutative_matrix_algorithms_unchanged():
    M = Matrix([[1, 2], [3, 4]])
    assert str(M.charpoly()) == "PurePoly(lambda**2 - 5*lambda - 2, lambda, domain='ZZ')"
    assert M.inv(method='GE') == Matrix([[-2, 1], [Rational(3, 2), Rational(-1, 2)]])
    assert M.solve(Matrix([1, 2])) == Matrix([0, Rational(1, 2)])
    assert M.echelon_form() == Matrix([[1, 2], [0, -2]])


def test_polarify_subs_mapping():
    # SymPy's mapping follows the hash order of free_symbols, which changes
    # between interpreter runs; the expression and the mapping are compared.
    expr, reps = polarify(x**2*y + z, subs=True)
    assert {str(k): v for k, v in reps.items()} == {'_x': x, '_y': y, '_z': z}
    assert expr.subs(reps) == x**2*y + z
