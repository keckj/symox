"""Public arithmetic structure must remain visible to exact replacement."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('expression,old,new', [
    ('sqrt(2)/2', 'sqrt(2)', 'sqrt(3)'),
    ('-sqrt(2)', 'sqrt(2)', 'sqrt(3)'),
    ('2*x+y', '2*x', 'z'), ('2*x', '2', '3'),
    ('x**2*y', 'x**2', 'z'), ('2*x**3', '3', '4'),
    ('2*x+3*y', '2', '5'), ('x/3', '1/3', '2/3'),
    ('sin(2*x)+cos(2*x)', '2', '3'),
])
def test_xreplace_public_coefficients_and_powers(expression, old, new):
    e, a, b = map(sx.sympify, (expression, old, new))
    result = e.xreplace({a: b})
    assert result.to_sympy() == sp.sympify(expression).xreplace({sp.sympify(old): sp.sympify(new)})


@pytest.mark.parametrize('expression,replacements', [
    ('sqrt(2)+2', {'sqrt(2)': 'sqrt(3)', '2': '7'}),
    ('x**2+x', {'x**2': 'y', 'x': '0'}),
    ('(x+y)/(x-y)', {'x+y': '0', 'x-y': '0'}),
    ('x**4+x**2', {'x**2': 'y', 'y': '0'}),
    ('sin(x)+x', {'sin(x)': 'y', 'x': '2'}),
])
def test_simultaneous_algebraic_patterns(expression, replacements):
    mapping = {sx.sympify(a): sx.sympify(b) for a, b in replacements.items()}
    expected_mapping = {sp.sympify(a): sp.sympify(b) for a, b in replacements.items()}
    actual = sx.sympify(expression).subs(mapping, simultaneous=True)
    expected = sp.sympify(expression).subs(expected_mapping, simultaneous=True)
    assert actual.to_sympy() == expected


@pytest.mark.parametrize('positive', [False, True])
def test_fractional_power_substitution_respects_real_branch(positive):
    x = sx.Symbol('x', positive=True) if positive else sx.Symbol('x')
    X = sp.Symbol('x', positive=True) if positive else sp.Symbol('x')
    y = sx.Symbol('y')
    assert (x*x).subs(x**3, y).to_sympy() == (X*X).subs(X**3, sp.Symbol('y'))
    assert (2**sx.Rational(1, 3)).subs(sx.sqrt(2), sx.sqrt(3)) == 3**sx.Rational(1, 3)


@pytest.mark.parametrize('expression', ['3+2*sqrt(2)', 'pi+sqrt(3)', 'exp(1/3)',
                                       '2**(1/3)+I', 'sin(1)', 'sqrt(2)/3'])
def test_nsimplify_does_not_approximate_exact_input(expression):
    original = sx.sympify(expression)
    assert sx.nsimplify(original) == original


def test_scalar_poly_rejects_matrix_input():
    x = sx.Symbol('x')
    a = sx.MatrixSymbol('A', 2, 2)
    with pytest.raises(TypeError, match='matrix'):
        sx.Poly(a, x)
    with pytest.raises(TypeError, match='[Mm]atrix'):
        sx.Poly(sx.Matrix([[1, 2], [3, 4]]), x)
    assert a.as_poly(x) is None
