"""A zero product describes a union of entire solution components."""
import pytest
import symox as sx
import sympy as sp


def foreign(value):
    if isinstance(value, dict):
        return {foreign(k): foreign(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return type(value)(foreign(v) for v in value)
    return value.to_sympy() if hasattr(value, 'to_sympy') else value


@pytest.mark.parametrize('flags', [{}, {'dict': True}, {'set': True}, {'check': False}])
@pytest.mark.parametrize('source', [
    '(x-1)*(y-2)', '(x-a)*(y-b)', '(x-1)**2*(y-2)**3',
    '[(x-1)*(y-2),x+y-4]',
    '[(x-1)*(y-2),(x-3)*(y-4)]',
])
@pytest.mark.parametrize('inferred', [False, True])
@pytest.mark.parametrize('expanded', [False, True])
def test_every_polynomial_factor_component_is_retained(source, flags, inferred, expanded):
    x, y, a, b = sx.symbols('x y a b')
    equations = eval(source, dict(x=x, y=y, a=a, b=b))
    if expanded:
        equations = [sx.expand(e) for e in equations] if isinstance(equations,list) else sx.expand(equations)
    variables = () if inferred else (x, y)
    actual = sx.solve(equations, *variables, **flags)
    expected = sp.solve(foreign(equations), *foreign(variables), **flags)
    assert foreign(actual) == expected


def test_each_family_contains_an_unbounded_free_parameter():
    x, y = sx.symbols('x y')
    # Exact symbolic families, not a finite list of numerically checked points.
    assert sx.solve((x-1)*(y-2)) == [{x: 1}, {y: 2}]
    assert sx.solve((x-1)*(y-2), x, y) == [(1, y), (x, 2)]
    assert sx.solve((x-1)*(y-2), set=True) == ([x, y], {(1, y), (x, 2)})


def test_regular_pivot_covers_the_whole_graph_without_division_by_a_free_parameter():
    x, y, z = sx.symbols('x y z')
    assert sx.solve([x*y-z]) == [{z: x*y}]
    assert sx.solve([x*y-z], x, y, z) == [(x,y,x*y)]
    # Both zero and nonzero x/y are included by this polynomial parametrization.
    assert sx.solve([x*y-z], set=True) == ([z], {(x*y,)})
