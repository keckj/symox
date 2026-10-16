"""Complete finite polynomial systems without a SymPy runtime dependency."""

import subprocess
import sys

import pytest
import sympy as sp
import symox as sm


def oracle(e):
    return sp.sympify(str(e))


def tuples_from_set(result):
    return [] if result == sm.EmptySet else [tuple(t.args) for t in result.args]


def compare_solutions(native, expected):
    assert len(native) == len(expected)
    remaining = [tuple(oracle(e) for e in row) for row in native]
    for row in expected:
        found = next((i for i, values in enumerate(remaining)
                      if all(sp.simplify(a-b) == 0 for a, b in zip(values, row))), None)
        if found is None:
            found = next((i for i, values in enumerate(remaining)
                          if all(not (a-b).free_symbols and abs(complex(sp.N(a-b))) < 1e-8
                                 for a, b in zip(values, row))), None)
        assert found is not None, (remaining, row)
        remaining.pop(found)


@pytest.mark.parametrize("system", [
    ["x**2-1", "y**2-1"],
    ["x**2+y**2-5", "x-y-1"],
    ["x*y-1", "y**2-x"],
    ["x**2-2", "y-x"],
    ["x**2+1", "y**2-x"],
    ["x**2", "y**3"],
    ["(x-1)*(x+2)*(x-3)", "y-x*x"],
    ["x**2-1", "y-x", "x+y-2"],
    ["x/2+y/3-1", "x-y"],
])
@pytest.mark.parametrize("reverse", [False, True])
def test_polynomial_systems_match_sympy_and_certify_substitution(system, reverse):
    variables = sm.symbols("x y")
    symbols = sp.symbols("x y")
    if reverse:
        variables, symbols = variables[::-1], symbols[::-1]
    equations = [sm.sympify(e) for e in system]
    expected = sp.solve_poly_system([sp.sympify(e) for e in system], *symbols)
    native = sm.solve_poly_system(equations, *variables)
    assert native is not None
    compare_solutions(native, expected)
    compare_solutions(tuples_from_set(sm.nonlinsolve(equations, variables)), expected)
    for solution in native:
        mapping = dict(zip(symbols, map(oracle, solution)))
        for equation in system:
            residual = sp.simplify(sp.sympify(equation).subs(mapping))
            assert residual == 0 or (not residual.free_symbols and abs(complex(sp.N(residual))) < 1e-8)


@pytest.mark.parametrize("system", [
    ["x**2-a", "y-x"],
    ["x+y-a", "y**2-b"],
    ["a*x+y-1", "x-y"],
    ["x**2-a", "y**2-b"],
])
def test_symbolic_parameters_in_triangular_systems(system):
    x, y = sm.symbols("x y")
    sx, sy = sp.symbols("x y")
    equations = [sm.sympify(e) for e in system]
    expected = sp.solve_poly_system([sp.sympify(e) for e in system], sx, sy)
    native = sm.solve_poly_system(equations, x, y)
    compare_solutions(native, expected)
    for row in native:
        for equation in equations:
            assert sp.simplify(oracle(equation).subs(dict(zip((sx, sy), map(oracle, row))))) == 0


def test_three_variable_system_inference_and_equations():
    x, y, z = sm.symbols("x y z")
    equations = [sm.Eq(x*x, 1), sm.Eq(y, x+1), sm.Eq(z, 2*y)]
    expected = [(-1, 0, 0), (1, 2, 4)]
    assert sm.solve_poly_system(equations, x, y, z) == expected
    assert sm.solve_poly_system(equations) == expected
    assert sm.solve_poly_system([x*x-1], x) == [(-1,), (1,)]


def test_inconsistency_is_distinct_from_unsupported_systems():
    x, y = sm.symbols("x y")
    for equations in ([1], [x, x-1], [x*x+1, x*x+2], [x, x*y-1]):
        assert sm.solve_poly_system(equations, x, y) is None
        assert sm.nonlinsolve(equations, x, y) == sm.EmptySet
    for equations in ([], [0], [x], [x*y-1]):
        with pytest.raises(NotImplementedError, match="zero-dimensional"):
            sm.solve_poly_system(equations, x, y)
    with pytest.raises(NotImplementedError, match="zero-dimensional"):
        sm.nonlinsolve([x*y-1], x, y)
    for strict in (False, True):
        with pytest.raises(NotImplementedError, match="complete solution"):
            sm.solve_poly_system([x**5-x+1, y-x], x, y, strict=strict)


def test_nonlinsolve_retains_parametric_linear_solutions():
    x, y = sm.symbols("x y")
    assert sm.nonlinsolve([x-y], x, y) == sm.linsolve([x-y], x, y)
    assert sm.nonlinsolve([x], x, y) == sm.linsolve([x], x, y)
    assert sm.nonlinsolve([], x, y) == sm.EmptySet
    assert sm.nonlinsolve([x-y, x-y-1], x, y) == sm.EmptySet


def test_invalid_polynomial_system_inputs_are_explicit():
    x, y, a = sm.symbols("x y a")
    with pytest.raises(ValueError, match="unique"):
        sm.solve_poly_system([x], x, x)
    with pytest.raises(ValueError, match="generators"):
        sm.solve_poly_system([1])
    with pytest.raises(ValueError, match="symbols"):
        sm.nonlinsolve([x*x-1])
    with pytest.raises(ValueError, match="generators"):
        sm.solve_poly_system([x], 1)
    with pytest.raises(ValueError):
        sm.solve_poly_system([sm.sin(x), y], x, y)
    with pytest.raises(NotImplementedError, match="exact coefficients"):
        sm.solve_poly_system([x*x-1.5, y-x], x, y)
    with pytest.raises(NotImplementedError, match="coefficient-field"):
        sm.solve_poly_system([x*x+y*y-a, x*y-1], x, y)
    with pytest.raises(NotImplementedError, match="options"):
        sm.solve_poly_system([x], x, extension=True)


def test_polynomial_solvers_do_not_import_sympy_at_runtime():
    subprocess.run([sys.executable, "-c", """
import sys
sys.modules['sympy'] = None
import symox as s
x,y=s.symbols('x y')
assert s.solve_poly_system([x*x-1,y-x],x,y) == [(-1,-1),(1,1)]
assert len(s.nonlinsolve([x*x-1,y-x],x,y).args) == 2
"""], check=True, timeout=30)
