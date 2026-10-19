"""Solve output depends on scalar/system input and symbol ordering."""
import pytest
import symox as sx
import sympy as sp


def foreign(value):
    if isinstance(value, dict):
        return {foreign(k): foreign(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return type(value)(foreign(v) for v in value)
    return value.to_sympy() if hasattr(value, 'to_sympy') else value


@pytest.mark.parametrize('source', [
    '[x-3,y**2-1]', '[x-3,y-1]', '[x+y-3]', 'x+y-3',
])
@pytest.mark.parametrize('arguments', ['inferred', 'positional', 'list', 'set', 'frozenset'])
@pytest.mark.parametrize('flags', [{}, {'dict': True}, {'set': True}])
def test_reference_output_shape_and_free_variables(source, arguments, flags):
    x, y = sx.symbols('x y')
    expressions = eval(source, {'x': x, 'y': y})
    variables = {
        'inferred': (), 'positional': (x, y), 'list': ([x, y],),
        'set': ({x, y},), 'frozenset': (frozenset([x, y]),),
    }[arguments]
    expected = sp.solve(foreign(expressions), *foreign(variables), **flags)
    actual = sx.solve(expressions, *variables, **flags)
    assert foreign(actual) == expected


@pytest.mark.parametrize('flags', [{}, {'dict': True}, {'set': True}])
def test_explicit_variable_order_is_preserved(flags):
    x, y = sx.symbols('x y')
    equations = [x-3, y*y-1]
    assert foreign(sx.solve(equations, y, x, **flags)) == sp.solve(
        foreign(equations), y.to_sympy(), x.to_sympy(), **flags)


def test_polynomial_graph_can_have_a_simpler_complete_parametrization():
    x, y = sx.symbols('x y')
    # Native linear elimination parametrizes this whole graph by x. The
    # reference chooses two square-root branches parametrized by y. Comparing
    # their raw list lengths would falsely classify an equivalent family.
    assert sx.solve([x*x-y]) == [{y: x*x}]
    rx, ry = foreign((x, y))
    branches = sp.solve([rx*rx-ry])
    assert len(branches) == 2
    t = sp.Symbol('t')
    assert sp.expand(sp.prod(t-row[rx] for row in branches)) == t*t-ry
    # Thus their roots are exactly all x with x²=y, for every complex y;
    # conversely every complex x occurs in native (x,x²), with no denominator.
    assert sx.solve([x*x-y], x, y) == [(x, x*x)]
    assert sx.solve([x*x-y], set=True) == ([y], {(x*x,)})


@pytest.mark.parametrize('flags', [{},{'dict':True},{'set':True}])
def test_native_tuples_are_equation_and_unknown_sequences(flags):
    x,y=sx.symbols('x y')
    equations=sx.Tuple(x-3,y*y-1)
    variables=sx.Tuple(x,y)
    assert foreign(sx.solve(equations,variables,**flags))==sp.solve(
        equations.to_sympy(),variables.to_sympy(),**flags)


def test_nsolve_native_tuple_inputs_are_sequences_too():
    x,y=sx.symbols('x y')
    vector=sx.nsolve(sx.Tuple(x+y-3,x-y-1),sx.Tuple(x,y),sx.Tuple(0,0))
    assert vector.shape==(2,1)
    assert abs(float(vector[0])-2)<1e-13 and abs(float(vector[1])-1)<1e-13
    root=sx.nsolve(x*x-2,x,sx.Tuple(1,2))
    assert abs(float(root)-2**0.5)<1e-13
