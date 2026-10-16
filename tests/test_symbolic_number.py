"""The structural is_number property is separate from finite evaluation."""
import pickle
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import pytest
import sympy as sp
import symox as sx


def expression(module, source):
    names = (
        'Add', 'Mul', 'Pow', 'Symbol', 'Dummy', 'Rational', 'Float',
        'Function', 'sin', 'cos', 'sqrt', 'exp', 'log', 'UnevaluatedExpr',
        'Tuple', 'Integral', 'Sum', 'Limit', 'Subs', 'Derivative', 'Lambda',
        'Order', 'hyper', 'Piecewise', 'Eq', 'S', 'pi', 'E', 'I', 'oo', 'zoo', 'nan',
        'Identity', 'MatrixSymbol', 'Trace', 'Determinant', 'ImmutableMatrix',
        'FunctionMatrix', 'DotProduct',
    )
    namespace = {name: getattr(module, name) for name in names}
    namespace.update(zip(('x', 'y', 'n'), module.symbols('x y n')))
    return module.sympify(eval(source, {'__builtins__': {}}, namespace))


@pytest.mark.parametrize('source', [
    '0', '1', '-2', 'Rational(1,2)', 'Float("1.25",80)',
    'pi', 'E', 'I', 'oo', '-oo', 'zoo', 'nan',
    'sin(1)', 'sqrt(2)', 'UnevaluatedExpr(nan)', 'UnevaluatedExpr(1)',
    'Symbol("a")', 'Dummy("a")', 'sin(x)',
    'Function("f")(1)', 'Function("f")()', 'Function("sin")(1)',
    'Function("Integral")(1)', 'Function("Lambda")(1)',
    'Tuple(1,2)', 'Tuple()', 'Lambda(x,1)', 'Lambda((),1)', 'Lambda(x,x)',
    'Order(1)', 'Order(x)', 'hyper((1,),(2,),3)',
    'Mul(0,x,evaluate=False)', 'Pow(x,0,evaluate=False)',
    'Add(0,x,evaluate=False)', 'Mul(oo,x,evaluate=False)',
    'Eq(x,1)',
    'Piecewise((1,x>0),(2,True))',
    'Integral(1,x)', 'Integral(1,(x,0,1))', 'Integral(1,(x,x,1))',
    'Integral(Function("f")(x),(x,0,1))',
    'Integral(Function("f")(),(x,0,1))',
    'Integral(Function("f")(y),(x,0,1))',
    'Integral(Integral(1,x),(x,0,1))',
    'Integral(Integral(1,y),(x,0,1))',
    'Integral(Sum(x,(x,0,y)),(y,0,1))',
    'Integral(Lambda(x,x+y),(y,0,1))',
    'Integral(Subs(Function("f")(x),x,y),(y,0,1))',
    'Integral(Subs(x+y,x,1),(x,0,1),(y,0,1))',
    'Integral(Limit(Function("f")(x),x,y),(y,0,1))',
    'Sum(Function("f")(),(x,0,1))', 'Sum(Function("f")(y),(x,0,1))',
    'Limit(1,x,0)', 'Subs(1,x,1)', 'Derivative(1,x,evaluate=False)',
    'Integral(Derivative(1,x,evaluate=False),(y,0,1))',
    'MatrixSymbol("A",2,2)', 'Identity(2)',
    'FunctionMatrix(2,2,Lambda((x,y),x+y))',
    'DotProduct(ImmutableMatrix([1,2]),ImmutableMatrix([3,4]))',
    'Trace(Identity(2))', 'Determinant(Identity(2))',
    'Integral(Trace(Identity(x)),(x,1,2))',
    'Integral(Trace(MatrixSymbol("A",2,2)),(x,0,1))',
])
def test_structural_number_classification_matches_sympy(source):
    value, expected = expression(sx, source), expression(sp, source)
    initial = value.__reduce__()[1]
    for _ in range(3):
        assert value.is_number is expected.is_number, source
    assert value.__reduce__()[1] == initial
    restored = pickle.loads(pickle.dumps(value))
    assert restored.is_number is expected.is_number
    assert sx.from_sympy(sx.to_sympy(value)).is_number is expected.is_number


@pytest.mark.parametrize('numeric_first', [False, True])
def test_cache_is_independent_of_numeric_refinement_policy(numeric_first):
    x = sx.Symbol('x')
    difference = sx.pi-sx.E
    if numeric_first:
        assert difference.is_positive is True
    assert difference.is_number is True
    observed = []
    def coefficient(value):
        observed.append((value.is_number, value.is_positive))
        return value
    assert sx.collect(x*difference, [x, difference], coefficient, evaluate=False) == {x: difference}
    assert observed == [(False, None)]
    assert difference.is_number is True and difference.is_positive is True


def test_same_expression_under_distinct_binding_scopes_and_query_orders():
    x, y = sx.symbols('x y')
    body = sx.Function('f')(x)
    assert body.is_number is False
    closed = sx.Integral(body, (x, 0, 1))
    free = sx.Integral(body, (y, 0, 1))
    combined = sx.Add(closed, free, evaluate=False)
    for _ in range(4):
        assert closed.is_number is True
        assert free.is_number is False
        assert body.is_number is False
        assert combined.is_number is False
    different_x = sx.Symbol('x', real=True)
    assert sx.Integral(body, (different_x, 0, 1)).is_number is False
    dummy = sx.Dummy('x')
    assert sx.Integral(body, (dummy, 0, 1)).is_number is False


def test_concurrent_number_and_assumption_queries_preserve_identity():
    values = [sx.pi-sx.E, sx.zoo, sx.nan, sx.Function('f')(1)]
    expected = [True, True, True, False]
    initial = [(hash(value), pickle.dumps(value)) for value in values]
    def query(_):
        for _ in range(20):
            assert [value.is_number for value in values] == expected
            assert values[0].is_positive is True
    with ThreadPoolExecutor(max_workers=8) as workers:
        list(workers.map(query, range(16)))
    assert [(hash(value), pickle.dumps(value)) for value in values] == initial


def test_depth80_structural_numbers_and_scopes_do_not_import_sympy_or_expand_dag():
    program = r'''
import builtins,pickle,sys
original_import=builtins.__import__
def guarded(name,*args,**kwargs):
    if name=='sympy' or name.startswith('sympy.'):
        raise AssertionError('SymPy runtime fallback')
    return original_import(name,*args,**kwargs)
builtins.__import__=guarded
import symox as s
one=s.Integer(1)
value=one
for _ in range(80):
    value=s._core._restore_arithmetic('Add',s.Integer(0),[(value,one),(value,one)],False)
before=value.__reduce__()[1]
assert value.is_number is True
assert value.__reduce__()[1]==before
assert pickle.loads(pickle.dumps(value)).is_number is True
x=s.Dummy('x')
y=s.Symbol('x')
pair=s.Function('pair')
body=x
for _ in range(80):body=pair(body,body)
assert body.is_number is False
closed=s.Integral(body,(x,0,1))
free=s.Integral(body,(y,0,1))
for _ in range(8):
    assert closed.is_number is True and free.is_number is False
    assert body.is_number is False
for protocol in range(pickle.HIGHEST_PROTOCOL+1):
    restored=pickle.loads(pickle.dumps(closed,protocol=protocol))
    assert restored==closed and restored.is_number is True
assert not any(name=='sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    subprocess.run([sys.executable, '-c', program], check=True, timeout=15)
