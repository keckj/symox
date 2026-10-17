"""Public matrix coefficient views preserve matrix order and raw tails."""
import pickle

import pytest
import sympy as sp
import symox as sx


def converted(value):
    if isinstance(value, (tuple, list)):
        return type(value)(converted(item) for item in value)
    return value.to_sympy()


def examples(module):
    A, B, C = (module.MatrixSymbol(name, 2, 2) for name in "ABC")
    x = module.Symbol("x")
    return [A, A+B, A*B, 2*A, -A, x*A, 2*x*A, A**2, A**-1,
            A.T, A.adjoint(), module.Identity(2), module.ZeroMatrix(2, 2),
            module.MatMul(A), module.MatMul(A, B, C),
            module.MatMul(2, A, B, C), module.MatAdd(A),
            module.MatAdd(A, A, B), module.MatMul(A, module.Identity(2))]


@pytest.mark.parametrize("index", range(19))
@pytest.mark.parametrize("method", ["as_coeff_mmul", "as_coeff_Mul", "as_base_exp"])
def test_public_matrix_views_against_sympy(index, method):
    actual, reference = examples(sx)[index], examples(sp)[index]
    result, expected = getattr(actual, method)(), getattr(reference, method)()
    assert converted(result) == expected
    if method == "as_coeff_mmul":
        coefficient, tail = result
        assert isinstance(tail, sx.MatMul)
        assert tail.args == actual.as_coeff_mmul()[1].args
        assert coefficient*tail.doit() == actual.doit()
        assert pickle.loads(pickle.dumps(tail)) == tail
    if method == "as_coeff_Mul":
        assert actual.as_coeff_Mul(rational=True) == result


@pytest.mark.parametrize("coefficient", ["2", "-2", "2.5", "-2.5", "oo", "-oo", "zoo", "nan", "pi", "I"])
@pytest.mark.parametrize("rational", [True, False])
def test_numeric_matrix_coefficients(coefficient, rational):
    a, A = sx.MatrixSymbol("A", 2, 2), sp.MatrixSymbol("A", 2, 2)
    actual = sx.MatMul(sx.sympify(coefficient), a)
    reference = sp.MatMul(sp.sympify(coefficient), A)
    assert converted(actual.as_coeff_mul(rational=rational)) == reference.as_coeff_mul(rational=rational)
    scalar, matrices = actual.as_coeff_matrices()
    assert isinstance(matrices, list)
    assert converted((scalar, matrices)) == reference.as_coeff_matrices()


@pytest.mark.parametrize("index", range(19))
def test_matrix_factor_tuple_and_two_terms(index):
    actual, reference = examples(sx)[index], examples(sp)[index]
    assert converted(actual.as_coeff_mul()) == reference.as_coeff_mul()
    if index == 16:
        with pytest.raises(NotImplementedError, match="GenericZeroMatrix"):
            actual.as_two_terms()
    elif isinstance(actual, (sx.MatAdd, sx.MatMul)):
        assert converted(actual.as_two_terms()) == reference.as_two_terms()


def test_symbolic_coefficients_and_noncommutative_matrix_order():
    A, B = sx.MatrixSymbol("A", 2, 3), sx.MatrixSymbol("B", 3, 4)
    x = sx.Symbol("x")
    expression = sx.MatMul(2*x, A, B)
    assert expression.as_coeff_matrices() == (2*x, [A, B])
    coefficient, tail = expression.as_coeff_mmul()
    assert coefficient == 2*x
    assert tail.args == (A, B)
    assert tail.shape == (2, 4)
    assert expression.as_coeff_mul() == (1, (2*x, A, B))
    assert expression.as_two_terms()[1].args == (A, B)


def test_inverse_alias_and_failures_are_explicit():
    A = sx.MatrixSymbol("A", 2, 2)
    assert A.I == A.inverse()
    assert A.I.I == A
    with pytest.raises(sx.ShapeError):
        sx.MatrixSymbol("R", 2, 3).I
    with pytest.raises(AttributeError):
        A.as_coeff_matrices()
    with pytest.raises(AttributeError):
        A.as_two_terms()
    assert (2*A).as_coeff_mul(A) == (2, (A,))
    with pytest.raises(NotImplementedError, match="GenericIdentity"):
        sx.MatMul(A).as_coeff_mul(A)
    with pytest.raises(NotImplementedError):
        (2*A).as_coeff_mul(extra=True)


def test_coefficient_views_do_not_import_sympy():
    import subprocess
    import sys
    program = """
import builtins
original = builtins.__import__
def blocked(name, *args, **kwargs):
    if name == 'sympy' or name.startswith('sympy.'):
        raise AssertionError('unexpected SymPy runtime dependency')
    return original(name, *args, **kwargs)
builtins.__import__ = blocked
import symox as s
A = s.MatrixSymbol('A', 2, 2)
x = s.Symbol('x')
e = x*A
assert e.as_coeff_matrices() == (x, [A])
assert e.as_coeff_mmul()[1].args == (A,)
assert e.as_coeff_Mul() == (1, e)
assert e.as_coeff_mul() == (1, (x, A))
assert (A**2).as_base_exp() == (A**2, 1)
assert A.I == A**-1
"""
    subprocess.run([sys.executable, "-c", program], check=True)


@pytest.mark.parametrize("dependency", ["x", "y", "A", "B", "A*B", "x*y", "2"])
@pytest.mark.parametrize("index", range(8))
def test_dependency_partition_against_sympy(index, dependency):
    def build(module):
        A, B = (module.MatrixSymbol(name, 2, 2) for name in "AB")
        x, y = module.symbols("x y")
        dependency_value = {"x": x, "y": y, "A": A, "B": B,
                            "A*B": A*B, "x*y": x*y, "2": module.Integer(2)}[dependency]
        values = [A, 2*A, x*A, 2*x*A, x*y*A*B, A+B, A**2, module.Trace(A)*B]
        return values[index], dependency_value
    actual, dep = build(sx)
    reference, sp_dep = build(sp)
    expected = reference.as_coeff_mul(sp_dep)
    from sympy.matrices.expressions.special import GenericIdentity
    if isinstance(expected[0], GenericIdentity):
        with pytest.raises(NotImplementedError, match="GenericIdentity"):
            actual.as_coeff_mul(dep)
    else:
        assert converted(actual.as_coeff_mul(dep)) == expected


def test_structural_queries_cover_typed_parents_and_matrix_subproducts():
    def build(module):
        A, B, C = (module.MatrixSymbol(name, 2, 2) for name in "ABC")
        x, y = module.symbols("x y")
        expression = module.MatMul(2*x*y + module.Trace(B), A, B, C)
        patterns = [A, B, C, A*B, B*C, A*C, C*B, x, y, x*y, 2*x,
                    module.Integer(2), module.Integer(1), module.Integer(7)]
        return expression, patterns
    actual, patterns = build(sx)
    reference, expected_patterns = build(sp)
    assert [actual.has(p) for p in patterns] == [reference.has(p) for p in expected_patterns]
    assert not actual.has()
    assert actual.has(sx.Symbol("missing"), sx.Symbol("x"))
    A, B = sx.MatrixSymbol("A", 2, 1), sx.MatrixSymbol("B", 1, 2)
    assert (sx.DotProduct(A, B)*sx.Identity(2)).has(B)


def test_dimension_arguments_and_reconstruction_are_structural():
    n, m = sx.symbols("n m", integer=True, positive=True)
    A = sx.MatrixSymbol("A", n, m)
    assert A.has(n, m)
    assert not A.has(sx.Symbol("A"))
    square = sx.MatrixSymbol("S", n, n)
    assert not (square*square).has(sx.Integer(1))
    assert not sx.MatMul(square, square).has(sx.Integer(2))
    assert (square**2).has(sx.Integer(2))
    assert (square + square.T).has(n)


@pytest.mark.parametrize("shape", [(0, 0), (2, 2)])
def test_inverse_public_arguments_include_minus_one(shape):
    A = sx.MatrixSymbol("InverseParent", *shape)
    inverse = sx.Inverse(A)
    assert inverse.args == (A, sx.Integer(-1))
    assert inverse.func(*inverse.args) == inverse
    assert inverse.has(-1)
    assert not inverse.has(1)
    assert sx.sympify(sx.srepr(inverse)) == inverse
    assert pickle.loads(pickle.dumps(inverse)).args == inverse.args



def test_dependency_partition_rejects_invalid_rectangular_reconstruction():
    A = sx.MatrixSymbol("A", 2, 3)
    B = sx.MatrixSymbol("B", 3, 4)
    C = sx.MatrixSymbol("C", 4, 5)
    # Dropping B would create an invalid A*C product, even without evaluation.
    with pytest.raises(sx.ShapeError):
        (A*B*C).as_coeff_mul(B)
