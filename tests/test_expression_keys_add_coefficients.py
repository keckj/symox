"""Public expression class/sort keys and leading Add coefficient views."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def normalized(value):
    if isinstance(value, tuple):
        return tuple(normalized(item) for item in value)
    return value.to_sympy() if hasattr(value, "to_sympy") else value


def expressions(module):
    x, y = module.symbols("x y")
    A = module.MatrixSymbol("A", 2, 2)
    return [x, 2*x, x+y, x**3, module.Integer(2), module.Rational(1, 3),
            module.Float("2.5"), module.I, module.oo, -module.oo, module.zoo,
            module.nan, module.pi, module.E, module.true, module.false,
            module.exp(x), module.sin(x), module.Function("f")(x),
            module.Derivative(module.Function("f")(x), x), module.Integral(x, x),
            module.Sum(x, (x, 0, 3)), module.UnevaluatedExpr(2*x),
            module.Tuple(x, y), A, 2*A, A+A.T, A**2, A**-1, A.T,
            module.Trace(A), module.Determinant(A), A[0, 1],
            module.DotProduct(module.MatrixSymbol("a", 2, 1), module.MatrixSymbol("b", 1, 2))]


@pytest.mark.parametrize("index", range(34))
def test_instance_class_keys_distinguish_full_class_from_sort_base(index):
    actual, expected = expressions(sx)[index], expressions(sp)[index]
    assert actual.class_key() == expected.class_key()
    assert normalized(actual.sort_key()) == expected.sort_key()
    assert normalized(actual.sort_key(order="lex")) == expected.sort_key(order="lex")
    assert actual.sort_key() == sx.default_sort_key(actual)
    if index == 1:
        assert actual.class_key() == (3, 0, "Mul")
        assert actual.sort_key()[0] == (2, 0, "Symbol")


@pytest.mark.parametrize("name", ["Number", "Integer", "Rational", "Float", "Symbol", "Dummy",
                                  "Add", "Mul", "Pow", "Tuple", "UnevaluatedExpr",
                                  "MatrixExpr", "MatrixSymbol", "MatAdd", "MatMul", "MatPow",
                                  "Inverse", "Transpose", "Identity", "ZeroMatrix"])
def test_constructor_class_keys(name):
    assert getattr(sx, name).class_key() == getattr(sp, name).class_key()


@pytest.mark.parametrize("source", ["0", "2", "-2", "Rational(2,3)", "Float('2.5')",
                                    "Float('-2.5')", "oo", "-oo", "nan", "zoo", "pi",
                                    "x", "2*x", "2+x+y", "Float('2.5')+x+y",
                                    "x+y", "x+oo"])
@pytest.mark.parametrize("rational", [False, True])
def test_add_coefficient_head_matches_sympy(source, rational):
    actual, reference = sx.sympify(source), sp.sympify(source)
    result = actual.as_coeff_Add(rational=rational)
    assert normalized(result) == reference.as_coeff_Add(rational=rational)
    assert pickle.loads(pickle.dumps(result)) == result


@pytest.mark.parametrize("items", [("x", "2", "y"), ("2", "x", "y"),
                                  ("2.5", "x", "y"), ("x", "2.5", "y"),
                                  ("2", "x", "3"), ("oo", "x", "y"),
                                  ("zoo", "x", "y"), ("nan", "x", "y")])
@pytest.mark.parametrize("rational", [False, True])
def test_held_sum_tail_preserves_original_argument_order(items, rational):
    actual = sx.Add(*(sx.sympify(item) for item in items), evaluate=False)
    expected = sp.Add(*(sp.sympify(item) for item in items), evaluate=False)
    head, tail = actual.as_coeff_Add(rational=rational)
    ref_head, ref_tail = expected.as_coeff_Add(rational=rational)
    assert normalized((head, tail)) == (ref_head, ref_tail)
    assert tuple(normalized(item) for item in tail.args) == ref_tail.args


@pytest.mark.parametrize("order", ["grlex", "grevlex", "rev-lex", "none", "invalid"])
def test_unimplemented_sort_order_is_not_silently_ignored(order):
    x = sx.Symbol("x")
    for expression in [x, sx.MatrixSymbol("A", 2, 2)]:
        with pytest.raises(NotImplementedError, match="lexicographic"):
            expression.sort_key(order=order)
        with pytest.raises(NotImplementedError, match="lexicographic"):
            sx.default_sort_key(expression, order=order)


def test_matrix_add_coefficient_is_scalar_zero_without_extraction():
    A = sx.MatrixSymbol("A", 2, 2)
    for expression in [A, 2*A, A+A.T, A**2, A.T]:
        assert expression.as_coeff_Add() == (0, expression)
        assert expression.as_coeff_Add(rational=True) == (0, expression)


def test_public_keys_preserve_shared_tuple_children_and_need_no_sympy():
    subprocess.run([sys.executable, "-c", """
import sys
sys.modules['sympy'] = None
import symox as s
x = s.Symbol('x')
e = x
for _ in range(80):
    e = s.Function('shared')(e, e)
key = e.sort_key()
seen = set()
pending = [key]
while pending:
    node = pending.pop()
    if isinstance(node, tuple) and id(node) not in seen:
        seen.add(id(node))
        pending.extend(node)
assert len(seen) < 10000
other = s.default_sort_key(e)
seen_pairs = set()
pending_pairs = [(key, other)]
while pending_pairs:
    left, right = pending_pairs.pop()
    if isinstance(left, tuple):
        pair = (id(left), id(right))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        assert isinstance(right, tuple) and len(left) == len(right)
        pending_pairs.extend(zip(left, right))
    else:
        assert left == right
assert e.class_key() == (4, 0, 'shared')
assert s.Add(2, x, evaluate=False).as_coeff_Add() == (2, x)
"""], check=True, timeout=15)
