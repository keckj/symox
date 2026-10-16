"""Native matrix expansion: scalar hints, held nodes and ordered products."""
import pickle
import subprocess
import sys

import pytest
import sympy as sm
import symox as sx


def contexts():
    native = [sx.MatrixSymbol(name, 2, 2) for name in 'ABCD']
    reference = [sm.MatrixSymbol(name, 2, 2) for name in 'ABCD']
    return (*native, *sx.symbols('x y')), (*reference, *sm.symbols('x y'))


def arithmetic_normal_form(expression):
    """Normalize only matrix addition/product, preserving every held wrapper."""
    if not getattr(expression, 'is_MatrixExpr', False):
        return expression
    args = tuple(arithmetic_normal_form(arg) for arg in expression.args)
    if isinstance(expression, (sm.MatAdd, sm.MatMul)):
        return expression.func(*args).doit(deep=False)
    if all(a is b for a, b in zip(args, expression.args)):
        return expression
    # Rebuild the exact class without calling its evaluating constructor.
    # In particular MatPow, Transpose and Inverse must remain held.
    return sm.Basic.__new__(type(expression), *args)


FLAGS = [{}, {'deep': False}, {'mul': False}, {'multinomial': False},
         {'mul': False, 'multinomial': True}, {'power_base': False},
         {'power_exp': False}, {'basic': False}]


@pytest.mark.parametrize('flags', FLAGS)
@pytest.mark.parametrize('constructor', [
    lambda A, B, C, D, x, y: (A+B)*(C+D),
    lambda A, B, C, D, x, y: (x+y)**2*(A+B),
    lambda A, B, C, D, x, y: (A+B)**3,
    lambda A, B, C, D, x, y: ((A+B)*(C+D))**2,
    lambda A, B, C, D, x, y: (x*A*B)**2,
    lambda A, B, C, D, x, y: A**(x+y),
    lambda A, B, C, D, x, y: ((x+y)*A+B)*((x-y)*C+D),
])
def test_default_expand_flags_follow_sympy(constructor, flags):
    native, reference = contexts()
    value, expected = constructor(*native), constructor(*reference)
    result = value.expand(**flags)
    # SymPy.expand can leave redundant MatMul(MatMul(...)) wrappers and
    # noncanonical MatAdd argument order. Compare their arithmetic normal
    # forms ONLY for MatAdd/MatMul. Other nodes remain held; calling global
    # doit() here could conceal an incorrect transpose/inverse/power expansion.
    assert arithmetic_normal_form(result.to_sympy()) == arithmetic_normal_form(expected.expand(**flags))
    assert sx.expand(value, **flags) == result
    assert result.shape == value.shape
    assert result.func(*result.args) == result
    assert pickle.loads(pickle.dumps(result)) == result


@pytest.mark.parametrize('wrapper', ['Transpose', 'Adjoint', 'Inverse',
                                      'HadamardProduct', 'KroneckerProduct',
                                      'BlockDiagMatrix', 'DiagonalOf',
                                      'DiagonalMatrix'])
@pytest.mark.parametrize('flags', [{}, {'deep': False}, {'mul': False}])
def test_deep_expansion_visits_but_does_not_evaluate_held_matrix_nodes(wrapper, flags):
    A, B, C, D, x, y = contexts()[0]
    a, b, c, d, u, v = contexts()[1]
    args, oracle_args = (((A+B)*(C+D),), ((a+b)*(c+d),))
    if wrapper in ('HadamardProduct', 'KroneckerProduct'):
        args += (A,)
        oracle_args += (a,)
    value = getattr(sx, wrapper)(*args)
    reference = getattr(sm, wrapper)(*oracle_args)
    result = value.expand(**flags)
    assert result.to_sympy() == reference.expand(**flags)
    assert result.kind == value.kind
    assert result.shape == value.shape
    if flags.get('deep') is False:
        assert result == value


def test_slice_block_and_diagonal_descendants_are_expanded():
    A, B, C, D, x, y = contexts()[0]
    a, b, c, d, u, v = contexts()[1]
    parent, oracle_parent = (A+B)*C, (a+b)*c
    pairs = [
        (parent[:1, :], oracle_parent[:1, :]),
        (sx.BlockMatrix([[parent, A], [B, C]]), sm.BlockMatrix([[oracle_parent, a], [b, c]])),
        (sx.DiagMatrix(parent[:1, :]), sm.DiagMatrix(oracle_parent[:1, :])),
    ]
    for value, reference in pairs:
        result = value.expand()
        assert result.to_sympy() == reference.expand()
        assert result.shape == value.shape
        assert value.expand(deep=False) == value
    # Reverse ranges use native Python slicing semantics instead of SymPy's
    # buggy constructor. Certify their matrix entries independently.
    view = ((x+y)**2*A+B)[::-1, ::-1]
    result = view.expand()
    assert result.parent == ((x*x+2*x*y+y*y)*A+B)
    actual_entries, original_entries = result.as_explicit(), view.as_explicit()
    assert actual_entries.shape == original_entries.shape
    assert all(sx.expand(a-b) == 0 for a,b in zip(actual_entries, original_entries))
    assert result.rowslice == view.rowslice and result.colslice == view.colslice


@pytest.mark.parametrize('flags', [{}, {'mul': False}, {'multinomial': False}, {'deep': False}])
def test_explicit_entries_use_the_common_scalar_expander(flags):
    x, y = sx.symbols('x y')
    u, v = sm.symbols('x y')
    dense = sx.ImmutableMatrix([[(x+y)**2, x*(x+y)], [0, 1]])
    oracle = sm.ImmutableMatrix([[(u+v)**2, u*(u+v)], [0, 1]])
    wrapped = sx.MatrixExpr(dense)
    assert wrapped.expand(**flags).to_sympy() == oracle.expand(**flags)
    A, a = sx.MatrixSymbol('A', 2, 2), sm.MatrixSymbol('A', 2, 2)
    assert (wrapped*A).expand(**flags).to_sympy() == (oracle*a).expand(**flags)


@pytest.mark.parametrize('flags', [
    {'trig': True}, {'log': False}, {'force': True},
    {'power_base': False}, {'power_exp': False}, {'mul': False},
])
def test_scalar_function_hints_are_shared_with_top_level(flags):
    x, y = sx.symbols('x y', positive=True)
    u, v = sm.symbols('x y', positive=True)
    A, a = sx.MatrixSymbol('A', 2, 2), sm.MatrixSymbol('A', 2, 2)
    values = [(sx.sin(x+y), sm.sin(u+v)), (sx.log(x*y), sm.log(u*v)),
              (sx.exp(x+y), sm.exp(u+v)), ((x*y)**sx.Symbol('z'), (u*v)**sm.Symbol('z'))]
    for native, reference in values:
        value = native*A
        result = value.expand(**flags)
        assert result.to_sympy() == (reference*a).expand(**flags)
        assert sx.expand(value, **flags) == result


@pytest.mark.parametrize('flags', [{}, {'deep': False}, {'basic': False}, {'mul': False}])
def test_only_addition_and_multiplication_receive_automatic_matrix_normalization(flags):
    A, B, C, D, x, y = contexts()[0]
    assert sx.MatAdd(A, A, evaluate=False).expand(**flags) == 2*A
    assert sx.MatMul(A, sx.Identity(2), evaluate=False).expand(**flags) == A
    held = [sx.MatPow(A, 1, evaluate=False), sx.HadamardProduct(A, A, evaluate=False),
            sx.Transpose(A+B), sx.Adjoint(A*B), sx.Inverse(A*B)]
    for value in held:
        assert value.expand(**flags) == value


@pytest.mark.parametrize('flags', [{'max_terms': 0}, {'max_terms': -1},
                                   {'max_terms': True}, {'max_terms': 1.5},
                                   {'modulus': 7}, {'matrix_power': True},
                                   {'unknown': False}])
def test_unsupported_or_invalid_options_fail_consistently(flags):
    A = sx.MatrixSymbol('A', 2, 2)
    for expand in [lambda: A.expand(**flags), lambda: sx.expand(A, **flags)]:
        with pytest.raises((TypeError, ValueError, NotImplementedError)):
            expand()


def test_matrix_and_scalar_limits_are_checked_before_expansion():
    A, B = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    factors = [A+(i+1)*B for i in range(12)]
    value = sx.MatMul(*factors, evaluate=False)
    for expand in [lambda: value.expand(max_terms=100), lambda: sx.expand(value, max_terms=100)]:
        with pytest.raises(ValueError, match='budget|limit'):
            expand()
    # Disabling multiplication preserves the factored expression, and no
    # irrelevant expansion cost is charged for its Cartesian product.
    assert value.expand(mul=False, max_terms=100).to_sympy() == value.to_sympy().expand(mul=False)
    x, y = sx.symbols('x y')
    with pytest.raises(ValueError, match='budget|limit'):
        ((x+y)**1_000_000_000*A).expand(max_terms=100)
    # MatPow is held in the SymPy-compatible API, so a huge matrix exponent
    # alone does not request a billion-word expansion.
    assert ((A+B)**1_000_000_000).expand(max_terms=100) == (A+B)**1_000_000_000
    assert ((A+B)*(A+B)).expand(max_terms=100) == (A+B)**2


def test_physical_dag_memo_preserves_alpha_equal_binder_presentations():
    x, y, t = sx.symbols('x y t')
    f = sx.Function('f')
    left, right = sx.Subs(f(x), x, t), sx.Subs(f(y), y, t)
    assert left == right and left.args != right.args
    matrix = sx.MatrixExpr(sx.ImmutableMatrix([[left, right]]))
    result = matrix.expand().as_explicit()
    assert result[0, 0].args == left.args
    assert result[0, 1].args == right.args
    A, B = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    held = sx.HadamardProduct(left*A, right*B, evaluate=False)
    expanded = held.expand()
    assert expanded.args[0].args[0].args == left.args
    assert expanded.args[1].args[0].args == right.args


def test_shared_matrix_dag_expands_once_per_physical_node():
    code = '''
import symox as s
A,B,C=[s.MatrixSymbol(n,2,2) for n in 'ABC']
f=(A+B)*C
for i in range(45):
    f=s.HadamardProduct(f,f,evaluate=False)
r=f.expand(max_terms=100)
assert r.is_HadamardProduct
while r.is_HadamardProduct:
    assert r.args[0] == r.args[1]
    r=r.args[0]
assert r == A*C+B*C
'''
    subprocess.run([sys.executable, '-c', code], timeout=15, check=True)


def test_matrix_power_product_printing_preserves_precedence():
    A, B, C, D, x, y = contexts()[0]
    assert str((x*A*B)**2) == '(x*A*B)**2'
    assert str((A*B)**-2) == '(A*B)**(-2)'
    assert str(sx.MatPow(sx.MatPow(A, 2, evaluate=False), 3, evaluate=False)) == '(A**2)**3'


def test_scalar_auto_distribution_matches_sympy_and_respects_mul_false():
    A, B = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    a, b = sm.MatrixSymbol('A', 2, 2), sm.MatrixSymbol('B', 2, 2)
    for coefficient, reference in [(sx.Symbol('x'), sm.Symbol('x')),
                                    (sx.pi, sm.pi), (sx.I, sm.I),
                                    (sx.sqrt(2), sm.sqrt(2)),
                                    (sx.Rational(1, 2), sm.Rational(1, 2)),
                                    (sx.Float('0.5'), sm.Float('0.5'))]:
        value, expected = coefficient*(A+B), reference*(a+b)
        assert value.to_sympy() == expected
        assert value.expand(mul=False).to_sympy() == expected.expand(mul=False)
        assert value.expand().to_sympy() == expected.expand()


def test_explicit_product_shape_and_work_are_in_the_shared_budget():
    # Each operand has ten entries, but their product creates 100 entries.
    left = sx.MatrixExpr(sx.ImmutableMatrix(10, 1, list(range(10))))
    right = sx.MatrixExpr(sx.ImmutableMatrix(1, 10, list(range(10))))
    held = sx.MatMul(left, right, evaluate=False)
    with pytest.raises(ValueError, match='budget|limit'):
        held.expand(max_terms=20)
    assert held.expand(max_terms=100).as_explicit() == left.as_explicit()*right.as_explicit()


def test_alternating_scalar_matrix_dag_shares_expansion_memo():
    code = '''
import symox as s
A,B,C=[s.MatrixSymbol(n,2,2) for n in 'ABC']
f=(A+B)*C
for i in range(40):
    t=s.Trace(f)
    f=s.MatrixExpr(s.ImmutableMatrix([[t,t],[0,1]]))
r=f.expand(max_terms=100)
for i in range(40):
    m=r.as_explicit()
    assert m[0,0] == m[0,1]
    r=m[0,0].args[0]
assert r == A*C+B*C
'''
    subprocess.run([sys.executable, '-c', code], timeout=15, check=True)


def test_native_scalar_coefficient_and_matrix_distribution_flags_are_distinct():
    A, B = sx.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('B', 2, 2)
    x, y = sx.symbols('x y')
    c, expanded = (x+y)**2, x*x+2*x*y+y*y
    value = c*(A+B)
    assert value.is_MatMul
    assert value.expand() == expanded*A+expanded*B
    assert value.expand(mul=False) == expanded*(A+B)
    assert value.expand(deep=False) == c*A+c*B
    assert value.expand(multinomial=False) == c*A+c*B
    assert value.expand(mul=False, multinomial=False) == value
    assert ((A+B)**2).expand().is_MatPow
