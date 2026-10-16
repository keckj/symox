"""Exact linear coefficient extraction and parametric solution compatibility."""

import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize("equations,symbols", [
    (["x+y-1"], "x y"),
    (["2*x+y-3", "x-y"], "x y"),
    (["x+y-1", "2*x+2*y-2"], "x y z"),
    (["x-1", "x-2"], "x"),
    (["0", "0"], "x y"),
    (["a*x+y-1", "x+a*y-2"], "x y"),
    (["a*(x+2*y)-3"], "x y"),
    (["x/y-1"], "x"),
    (["sin(a)*x+cos(a)*y-1"], "x y"),
    ([], "x y"),
])
def test_coefficients_and_solutions_against_sympy(equations, symbols):
    variables = tuple(sx.sympify(name) for name in symbols.split())
    expressions = [sx.sympify(value) for value in equations]
    a, b = sx.linear_eq_to_matrix(expressions, variables)
    pa, pb = sp.linear_eq_to_matrix([value.to_sympy() for value in expressions], [v.to_sympy() for v in variables])
    assert a.to_sympy() == pa and b.to_sympy() == pb
    actual = sx.linsolve(expressions, variables)
    expected = sp.linsolve([value.to_sympy() for value in expressions], [v.to_sympy() for v in variables])
    assert actual.to_sympy() == expected
    for solution in actual:
        mapping = dict(zip(variables, solution.args))
        for expression in expressions:
            assert sx.simplify(expression.xreplace(mapping)) == 0


@pytest.mark.parametrize("expression", ["x**2", "x*y", "sin(x)", "1/x", "(x**2-1)/(x-1)", "y**2-y*(y-1)+x"])
def test_strict_nonlinearity_is_not_cancelled(expression):
    x, y = sx.symbols("x y")
    for solver in (sx.linear_eq_to_matrix, sx.linsolve):
        with pytest.raises(sx.NonlinearError):
            solver([sx.sympify(expression)], x, y)


def test_matrix_formats_generated_parameters_and_collision_avoidance():
    x, y, z, tau0 = sx.symbols("x y z tau0")
    a, b = sx.Matrix([[1, 2, 3], [2, 4, 6]]), sx.Matrix([1, 2])
    for system in ((a, b), a.row_join(b)):
        assert sx.linsolve(system, x, y, z).to_sympy() == sp.linsolve((a.to_sympy(), b.to_sympy()), sp.symbols("x y z"))
        assert sx.linsolve(system).to_sympy() == sp.linsolve((a.to_sympy(), b.to_sympy()))
        assert sx.linsolve(system, iter([x, y, z])) == sx.linsolve(system, x, y, z)
    collision = sx.linsolve((sx.Matrix([[1, 1]]), sx.Matrix([tau0])))
    solution = next(iter(collision)).args
    assert sx.expand(sum(solution)-tau0) == 0
    assert tau0 in collision.free_symbols and len(collision.free_symbols) == 2
    assert sx.linsolve((sx.Matrix(0, 2, []), sx.Matrix(0, 1, [])), x, y).to_sympy() == sp.FiniteSet(sp.Tuple(*sp.symbols("x y")))


def test_equalities_functions_composite_unknowns_and_assumptions():
    x, y = sx.symbols("x y")
    a, b = sx.Symbol("same", positive=True), sx.Symbol("same", negative=True)
    assert next(iter(sx.linsolve([a+b-1], a, b))).args == (1-b, b)
    f = sx.Function("f")
    for variable in [f(x), x*x]:
        result = sx.linsolve([variable+y-1], variable, y)
        assert next(iter(result)).args == (1-y, y)
    assert sx.linsolve([sx.Eq(2*x+y, 3), sx.Eq(x-y, 0)], x, y).to_sympy() == sp.FiniteSet(sp.Tuple(1, 1))


def test_solve_linear_system_dict_and_errors():
    x, y = sx.symbols("x y")
    for values in [[[1, 1, 1]], [[1, 1, 1], [1, -1, 3]], [[0, 0, 0]], [[1, 0, 1], [1, 0, 2]]]:
        matrix = sx.Matrix(values)
        actual = sx.solve_linear_system(matrix, x, y)
        expected = sp.solve_linear_system(matrix.to_sympy(), *sp.symbols("x y"))
        assert actual is None if expected is None else {k.to_sympy(): v.to_sympy() for k, v in actual.items()} == expected
    for solver in (sx.linear_eq_to_matrix, sx.linsolve):
        with pytest.raises(ValueError):
            solver([x+y], x, x)
        with pytest.raises(ValueError):
            solver([x+y])
        with pytest.raises(TypeError):
            solver([x+y], {x, y})
    with pytest.raises(ValueError):
        sx.linsolve((sx.Matrix([[1, 2]]), sx.Matrix([1])), x)
    with pytest.raises(ValueError):
        sx.linsolve((sx.eye(2), sx.ones(2, 2)), x, y)


def test_native_linear_functions_do_not_import_sympy():
    script = "import sys; import symox as s; x,y=s.symbols('x y'); assert str(s.linsolve([x+y-1],x,y))=='{(1 - y, y)}'; assert 'sympy' not in sys.modules"
    subprocess.run([sys.executable, "-c", script], check=True, capture_output=True, text=True)


def test_conjugate_integer_powers_preserves_complex_semantics():
    x, n = sx.Symbol("x"), sx.Symbol("n", integer=True)
    for exponent in [sx.Integer(2), sx.Integer(-3), n]:
        assert sx.conjugate(x**exponent) == sx.conjugate(x)**exponent
    assert sx.conjugate(sx.sqrt(x)) != sx.sqrt(sx.conjugate(x))
