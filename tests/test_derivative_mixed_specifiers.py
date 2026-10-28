"""A scalar variable followed by a pair starts two differentiation specs."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('carrier', ['tuple', 'list', 'Tuple'])
@pytest.mark.parametrize('pattern', [0, 1, 2])
@pytest.mark.parametrize('symbolic', [False, True])
@pytest.mark.parametrize('evaluate', [False, True])
def test_mixed_variable_and_pair_specifiers(carrier, pattern, symbolic, evaluate):
    def construct(module):
        x, y, z, n = module.symbols('x y z n')
        count = n if symbolic else module.Integer(2)
        pair = {'tuple': tuple, 'list': list}.get(carrier)
        if pair is None:
            pair = lambda values: module.Tuple(*values)
        specs = [(x, pair((y, count))),
                 (x, pair((y, count)), z),
                 (pair((x, 2)), y, pair((z, count)))][pattern]
        body = module.Function('f')(x, y, z)
        return module.Derivative(body, *specs, evaluate=evaluate)

    assert construct(sx).to_sympy() == construct(sp)


@pytest.mark.parametrize('carrier', [tuple, list, sx.Tuple])
def test_diff_uses_the_same_mixed_specifier_parser(carrier):
    x, y = sx.symbols('x y')
    pair = carrier(y, 2) if carrier is sx.Tuple else carrier((y, 2))
    assert sx.diff(x*y**3, x, pair) == 6*y
