"""Powers of powers combine like SymPy's Pow._eval_power."""
import pytest
import sympy as sp
import symox as sx

SOURCES = [
    'sqrt(sqrt(t))', 'sqrt(t)**y', '(t**(S(1)/3))**y', '(t**(S(2)/3))**(S(1)/2)',
    '(t**(S(1)/2))**(-S(1)/2)', '(t**(-S(1)/2))**(S(1)/2)', '(t**2)**(S(1)/2)',
    '(t**(S(3)/2))**(S(1)/2)', '(t**-1)**(S(1)/2)', '(t**y)**(S(1)/2)', '(t**0.5)**y',
    '(t**-0.5)**(S(1)/3)', 'sqrt(x + 1)**(S(1)/3)', '((-t)**(S(1)/2))**(S(2)/3)',
    'sqrt(1/sqrt(t))', '(sqrt(t))**(S(1)/3)', 'sqrt(sqrt(t) + 1)', '(t**(S(1)/2))**3',
]


def _evaluate(module, source):
    t, x, y = module.symbols('t x y')
    return eval(source, {'t': t, 'x': x, 'y': y, 'S': module.S, 'sqrt': module.sqrt})


@pytest.mark.parametrize('source', SOURCES)
def test_power_of_power_matches_sympy(source):
    actual, expected = _evaluate(sx, source), _evaluate(sp, source)
    # str rather than srepr: Float digits are spelled differently by srepr.
    assert str(actual) == str(expected)
    assert [type(a).__name__ for a in sx.preorder_traversal(actual)] == \
        [type(a).__name__ for a in sp.preorder_traversal(expected)]
