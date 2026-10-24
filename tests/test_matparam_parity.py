"""Matrix domain inference, Float/symbolic parameter routes, identity
entries and matrix exp/log forms (wave 16 matparam).

Expected values are the printed results of SymPy 1.14 for the same
statements.
"""
import pytest

import symox as sx


def _namespace():
    names = {name: getattr(sx, name) for name in dir(sx)
             if not name.startswith('_')}
    a, x, r, p = sx.symbols('a x r p')
    names.update(a=a, x=x, r=r, p=p,
                 W=sx.sin(x)**2 + sx.cos(x)**2)
    return names


CASES = [
    ('Matrix([[sqrt(a), 1],[1, sqrt(a)]]).to_DM().domain',
     'ZZ[sqrt(a)]'),
    ('Matrix([[exp(a), 1],[1, 1/exp(a)]]).to_DM().domain',
     'ZZ(exp(a))'),
    ('Matrix([[I*a, 1],[1, a]]).to_DM().domain',
     'ZZ_I[a]'),
    ('Matrix([[1.0*a, 1],[1, a]]).to_DM().domain',
     'RR[a]'),
    ('Matrix([[sqrt(2), 0],[0, 1.0]]).to_DM().domain',
     'EXRAW'),
    ('Matrix([[pi, 1],[1, pi]]).to_DM().domain',
     'ZZ[pi]'),
    ('Matrix([[sin(x), 1],[0, 1]]).to_DM().domain',
     'ZZ[sin(x)]'),
    ('Matrix([[1.0*a, 1],[1, 1.0*a]]).inv()',
     'Matrix([[1.0*a/(1.0*a**2 - 1.0), -1.0/(1.0*a**2 - 1.0)], [-1.0/(1.0*a**2 - 1.0), 1.0*a/(1.0*a**2 - 1.0)]])'),
    ("Matrix([[1.0*a, 1],[1, 1.0*a]]).inv(method='DM')",
     'Matrix([[1.0*a/(1.0*a**2 - 1.0), -1.0/(1.0*a**2 - 1.0)], [-1.0/(1.0*a**2 - 1.0), 1.0*a/(1.0*a**2 - 1.0)]])'),
    ('Matrix([[0.5*a, 1],[2, 1.5*a]]).inv()',
     'Matrix([[1.5*a/(0.75*a**2 - 2.0), -1.0/(0.75*a**2 - 2.0)], [-2.0/(0.75*a**2 - 2.0), 0.5*a/(0.75*a**2 - 2.0)]])'),
    ('Matrix([[1.0*a, 1.0*a],[1, 1]]).inv()',
     'NonInvertibleMatrixError: Matrix det == 0; not invertible.'),
    ('Matrix([[1.0*a, 2],[1, a]]).charpoly(x).domain',
     'RR[a]'),
    ('Matrix([[1.0*a, 2],[1, a]]).charpoly(x).as_expr()',
     '1.0*a**2 - 2.0*a*x + 1.0*x**2 - 2.0'),
    ('Matrix([[1.0*a, 1],[1, 1.0*a]]).eigenvals()',
     '{1.0*a - 1.0: 1, 1.0*a + 1.0: 1}'),
    ('Matrix([[W, 0],[0, W]]).charpoly(x).as_expr()',
     '_x**2 - 2*_x + 1'),
    ('Matrix([[W, 1],[1, W]]).eigenvals()',
     '{2: 1, 0: 1}'),
    ('Matrix([[W, 0],[0, 1]]).singular_values()',
     '[1, 1]'),
    ('Matrix([[W, 1],[1, 1]]).singular_values()',
     '[2, 0]'),
    ('cancel(log(1.0*I))',
     'I*pi/2'),
    ('cancel(log(2.0*I))',
     '0.693147180559945 + 0.5*I*pi'),
    ('cancel(log(1.0*I)+x)',
     'x + I*pi/2'),
    ('Matrix([[0, 1.0],[-1.0, 0]]).log()',
     'Matrix([[0, 0.5*pi], [-0.5*pi, 0]])'),
    ('Matrix([[0, r],[-r, 0]]).exp()',
     'Matrix([[exp(I*r)/2 + exp(-I*r)/2, -I*exp(I*r)/2 + I*exp(-I*r)/2], [I*exp(I*r)/2 - I*exp(-I*r)/2, exp(I*r)/2 + exp(-I*r)/2]])'),
    ('Matrix([[0, x],[-x, 0]]).exp()',
     'Matrix([[exp(I*x)/2 + exp(-I*x)/2, -I*exp(I*x)/2 + I*exp(-I*x)/2], [I*exp(I*x)/2 - I*exp(-I*x)/2, exp(I*x)/2 + exp(-I*x)/2]])'),
    ('Matrix([[a, 1],[0, a]]).exp()',
     'Matrix([[exp(a), exp(a)], [0, exp(a)]])'),
    ('Matrix([[x, 1],[1, x]]).exp()',
     'Matrix([[exp(x - 1)/2 + exp(x + 1)/2, -exp(x - 1)/2 + exp(x + 1)/2], [-exp(x - 1)/2 + exp(x + 1)/2, exp(x - 1)/2 + exp(x + 1)/2]])'),
    ('Matrix([[0, 1],[-x**2, 0]]).exp()',
     'Matrix([[exp(I*x)/2 + exp(-I*x)/2, -I*exp(I*x)/(2*x) + I*exp(-I*x)/(2*x)], [I*x*exp(I*x)/2 - I*x*exp(-I*x)/2, exp(I*x)/2 + exp(-I*x)/2]])'),
]


@pytest.mark.parametrize('statement,expected', CASES)
def test_matches_reference(statement, expected):
    try:
        value = str(eval(statement, _namespace()))
    except Exception as error:  # the reference error class and message
        value = '%s: %s' % (type(error).__name__, error)
    assert value == expected


def test_inexact_composite_inverse_round_trips():
    a = sx.Symbol('a')
    M = sx.Matrix([[1.0*a, 1], [1, 1.0*a]])
    error = (M*M.inv() - sx.eye(2)).subs(a, 3)
    assert all(abs(complex(entry)) < 1e-12 for entry in error)


def test_singular_inexact_composite_raises():
    a = sx.Symbol('a')
    with pytest.raises(sx.NonInvertibleMatrixError):
        sx.Matrix([[1.0*a, 1.0*a], [1, 1]]).inv()


def test_charpoly_evaluate_context_independent():
    x = sx.Symbol('x')
    W = sx.sin(x)**2 + sx.cos(x)**2
    M = sx.Matrix([[W, 0], [0, W]])
    expected = M.charpoly().as_expr()
    with sx.evaluate(False):
        value = M.charpoly().as_expr()
    assert str(value) == str(expected)
