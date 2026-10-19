"""Dense symbolic factorizations retain exact arithmetic and usable expressions."""
import numpy as np
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('method', ['cholesky_solve', 'LDLsolve'])
def test_dense_symbolic_solve_exact_specialization_and_complex_values(method):
    n = 6
    coefficients = sx.symbols('a:21')
    entries = iter(coefficients)
    matrix = sx.zeros(n)
    for i in range(n):
        for j in range(i, n):
            matrix[i, j] = matrix[j, i] = next(entries)
    values = sx.symbols('b:12')
    rhs = sx.Matrix(n, 2, values)
    solution = getattr(matrix, method)(rhs)
    rational = sx.Matrix(n, n, lambda i, j:
        10+i if i == j else sx.Rational(i+j+1, 10))
    concrete_rhs = sx.Matrix(n, 2, range(1, 13))
    substitutions = {matrix[i, j]: rational[i, j]
                     for i in range(n) for j in range(i, n)}
    substitutions.update(zip(values, concrete_rhs))
    expected = rational.to_sympy().inv() * concrete_rhs.to_sympy()
    actual = solution.xreplace(substitutions).to_sympy()
    assert actual.applyfunc(sp.simplify) == expected
    evaluate = sx.lambdify([*coefficients, *values], solution, 'numpy', cse=True)
    for factor in [1.0, 1.0+0.5j]:
        numeric = np.array(rational.tolist(), dtype=float)*factor
        input_values = [numeric[i, j] for i in range(n) for j in range(i, n)]
        calculated = evaluate(*input_values, *range(1, 13))
        np.testing.assert_allclose(calculated,
            np.linalg.solve(numeric, np.arange(1,13).reshape(n,2)), rtol=2e-13, atol=2e-13)


@pytest.mark.parametrize('method', ['cholesky_solve', 'LDLsolve'])
def test_symbolic_singular_system_retains_noninvertible_error(method):
    x = sx.Symbol('x', nonzero=True)
    matrix = sx.Matrix([[x, x], [x, x]])
    with pytest.raises(sx.NonInvertibleMatrixError):
        getattr(matrix, method)(sx.Matrix([1, 2]))


@pytest.mark.parametrize('method', ['cholesky_solve', 'LDLsolve'])
def test_polynomial_rank_deficiency_is_not_hidden_by_compact_pivots(method):
    x,y = sx.symbols('x y', real=True)
    rectangular = sx.Matrix([[x, 1], [y, x], [1, y]])
    matrix = rectangular*rectangular.T
    assert matrix.to_sympy().rank() == 2
    with pytest.raises(sx.NonInvertibleMatrixError):
        getattr(matrix, method)(sx.ones(3,1))


@pytest.mark.parametrize('method', ['cholesky_solve', 'LDLsolve'])
def test_rank_probes_do_not_expand_enormous_symbolic_powers(method):
    x = sx.Symbol('x', positive=True)
    power = x**(2**80)
    matrix = sx.diag(power, 1)
    solution = getattr(matrix, method)(sx.ones(2,1))
    assert solution[0]*power == 1
    assert solution[1] == 1
