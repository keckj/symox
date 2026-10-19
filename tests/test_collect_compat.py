"""Native coefficient collection against exact SymPy oracles."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

CASES = [
    ('a*x**2+b*x**2+a*x-b*x+c', ['x']),
    ('x**2*y+x*y+x+y', ['x', 'y']),
    ('a*x+b*x+a*y+b*y', ['y', 'x']),
    ('(a*x+b*x)*(y+z)', ['x']),
    ('(a*x+b*x+1)**2', ['x']),
    ('a*sin(x)+b*sin(x)+c', ['sin(x)']),
    ('a*sin(x)**2+b*sin(x)**2+c*sin(x)', ['sin(x)']),
    ('a*sin(x)*cos(x)+b*sin(x)*cos(x)+a*cos(x)', ['sin(x)', 'cos(x)']),
    ('a*x*y+b*x*y+a*x+b*x', ['x*y', 'x']),
    ('a*x**2*y**2+b*x**2*y**2+c*x', ['x*y']),
    ('a*(x+y)+b*(x+y)+c', ['x+y']),
    ('a*x**(2*c)+b*x**(2*c)+a*x**c', ['x**c']),
    ('a*x**c+b*x**c', ['x']),
    ('a*sqrt(x)+b*sqrt(x)+c*x', ['x']),
    ('a*x+b*x+c*sqrt(x)', ['sqrt(x)']),
    ('a/x+b/x+c*x', ['x']),
    ('a*x+b*x+c/x', ['1/x']),
    ('a*x**2+b*x**2+c*x', ['-x']),
    ('a*exp(2*x)+b*exp(2*x)+c*exp(x)', ['exp(x)']),
    ('a*exp(-x)+b*exp(-x)', ['exp(x)']),
    ('a*exp(x/2)+b*exp(x/2)', ['exp(x)']),
    ('a*f(x)+b*f(x)+a*Derivative(f(x),x)+b*Derivative(f(x),x)', ['f(x)']),
    ('a*Derivative(f(x),(x,2))+b*Derivative(f(x),(x,2))', ['Derivative(f(x),x)']),
    ('a*Derivative(f(x),(x,2))**2+b*Derivative(f(x),(x,2))**2', ['f(x)']),
    ('a*x+b*x+z', []),
    ('2*x+2', ['2']),
    ('x+1', ['-1']),
    ('0', ['x']),
    ('0', ['0']),
    ('3', ['x']),
]


def native(text):
    return sx.sympify(text)


def oracle(text):
    return sp.sympify(text)


def converted(result):
    if isinstance(result, dict):
        return {sx.to_sympy(k): sx.to_sympy(v) for k, v in result.items()}
    return sx.to_sympy(result)


@pytest.mark.parametrize('text,patterns', CASES)
@pytest.mark.parametrize('evaluate', [True, False])
@pytest.mark.parametrize('exact', [False, True])
def test_collect_exact_oracle_and_reconstruction(text, patterns, evaluate, exact):
    expression = native(text)
    result = sx.collect(expression, list(map(native, patterns)), evaluate=evaluate, exact=exact)
    expected = sp.collect(oracle(text), list(map(oracle, patterns)), evaluate=evaluate, exact=exact)
    assert isinstance(result, dict) == (not evaluate)
    assert converted(result) == expected
    if not evaluate:
        reconstructed = sum((key*value for key, value in result.items()), sx.Integer(0))
    else:
        reconstructed = result
    assert sx.expand(reconstructed-expression) == 0
    assert expression.collect(list(map(native, patterns)), evaluate=evaluate, exact=exact) == result


@pytest.mark.parametrize('text,patterns', [
    ('a*x+b*x+x**2', ['x']),
    ('2*a*x+3*b*x+c', ['x']),
    ('(a*x+b*x)*(c+y)', ['x']),
    ('(a*x+b*x)**3', ['x']),
    ('a*x+b*x+a*y+b*y', ['x', 'y']),
    ('a*sin(x)+b*sin(x)+c', ['sin(x)']),
])
@pytest.mark.parametrize('evaluate', [False, True])
@pytest.mark.parametrize('callback', ['square', 'increment', 'factor'])
def test_callbacks_match_recursive_coefficients_and_call_order(text, patterns, evaluate, callback):
    native_calls, oracle_calls = [], []
    def call_native(coefficient):
        native_calls.append(sx.to_sympy(coefficient))
        return coefficient**2 if callback == 'square' else coefficient+1 if callback == 'increment' else sx.factor(coefficient)
    def call_oracle(coefficient):
        oracle_calls.append(coefficient)
        return coefficient**2 if callback == 'square' else coefficient+1 if callback == 'increment' else sp.factor(coefficient)
    result = sx.collect(native(text), list(map(native, patterns)), call_native, evaluate=evaluate)
    expected = sp.collect(oracle(text), list(map(oracle, patterns)), call_oracle, evaluate=evaluate)
    assert converted(result) == expected
    # Composite patterns are replaced by internal dummies in both engines;
    # those dummies never occur in the coefficients in this corpus.
    assert native_calls == oracle_calls


@pytest.mark.parametrize('evaluate', [False, True])
@pytest.mark.parametrize('distribute', [False, True])
@pytest.mark.parametrize('text', ['a*x+b*x+O(y)', 'a*x+b*x+O(x**3)', 'O(y)', 'O(x)'])
def test_order_terms_obey_distribution_option(text, evaluate, distribute):
    actual = sx.collect(native(text), native('x'), evaluate=evaluate, distribute_order_term=distribute)
    expected = sp.collect(oracle(text), oracle('x'), evaluate=evaluate, distribute_order_term=distribute)
    assert converted(actual) == expected


@pytest.mark.parametrize('flags', [{'positive': True}, {'real': True}, {'integer': True}, {'real': False}])
def test_collection_keeps_symbol_and_dummy_assumptions(flags):
    x, X = sx.Symbol('x', **flags), sp.Symbol('x', **flags)
    a, b = sx.symbols('a b')
    A, B = sp.symbols('a b')
    pattern, PATTERN = sx.sqrt(x), sp.sqrt(X)
    expression = a*pattern+b*pattern+a*x
    expected_expression = A*PATTERN+B*PATTERN+A*X
    for evaluate in (False, True):
        actual = sx.collect(expression, pattern, evaluate=evaluate)
        expected = sp.collect(expected_expression, PATTERN, evaluate=evaluate)
        assert converted(actual) == expected
    d, other = sx.Dummy('x', **flags), sx.Dummy('x', **flags)
    expression = a*d+b*d+other
    result = sx.collect(expression, d, evaluate=False)
    assert result == {d: a+b, sx.Integer(1): other}
    assert d != other
    assert all(value.free_symbols <= {a, b, d, other} for pair in result.items() for value in pair)


def test_iterable_patterns_and_defaults_are_explicit():
    x, y, a, b = sx.symbols('x y a b')
    expression = a*x+b*x+a*y+b*y
    expected = sx.collect(expression, [x, y])
    for patterns in ((x, y), sx.Tuple(x, y), iter([x, y])):
        assert sx.collect(expression, patterns) == expected
    assert sx.collect(expression, x, evaluate=None) == sx.collect(expression, x)
    assert sx.collect(expression, x, func=None) == sx.collect(expression, x)
    with pytest.raises(TypeError, match='unexpected keyword'):
        sx.collect(expression, x, mystery=True)
    assert converted(sx.collect(expression,x,exact=None))==sp.collect(oracle('a*x+b*x+a*y+b*y'),oracle('x'),exact=None)
    with pytest.raises(TypeError, match='callable'):
        sx.collect(expression, x, func=object())
    with pytest.raises(TypeError):
        sx.collect(expression, [x, object()])


@pytest.mark.parametrize('evaluate',[False,True])
@pytest.mark.parametrize('exact',[False,True])
def test_method_positional_options(evaluate,exact):
    x,a,b=sx.symbols('x a b')
    expression=a*x+b*x+x*x
    transform=lambda coefficient:coefficient+1
    expected=sx.collect(expression,x,transform,evaluate,exact,False)
    assert expression.collect(x,transform,evaluate,exact,False)==expected


def test_wild_matching_and_multivariate_derivative_gap_are_explicit():
    x, y, a, b = sx.symbols('x y a b')
    assert converted(sx.collect(a*x+b*x, sx.Wild('w')))==sp.collect(oracle('a*x+b*x'),sp.Wild('w'))
    mixed = sx.Derivative(sx.Function('f')(x, y), x, y, evaluate=False)
    with pytest.raises(NotImplementedError, match='multivariate derivatives'):
        sx.collect(a*mixed+b*mixed, sx.Function('f')(x, y))


def test_callback_exception_and_nonexpression_return_are_not_swallowed():
    x, a = sx.symbols('x a')
    exception = RuntimeError('callback sentinel')
    def broken(_):
        raise exception
    with pytest.raises(RuntimeError) as caught:
        sx.collect(a*x, x, broken)
    assert caught.value is exception
    with pytest.raises(TypeError):
        sx.collect(a*x, x, lambda coefficient: object())


def test_no_sympy_runtime_and_large_shared_expression():
    program = r'''
import builtins, sys
real_import = builtins.__import__
def guarded(name,*args,**kwargs):
    if name == 'sympy' or name.startswith('sympy.'):
        raise AssertionError('unexpected SymPy runtime dependency')
    return real_import(name,*args,**kwargs)
builtins.__import__ = guarded
import symox as s
x,a,b = s.symbols('x a b')
assert s.collect(a*x+b*x,x,evaluate=False)=={x:a+b}
assert s.collect(a*x+b*x,x,func=lambda c:c**2)==x*(a**2+b**2)**2
value=x
pair=s.Function('pair')
for _ in range(45):value=pair(value,value)
assert s.collect(value,x)==value
calls=[]
assert s.collect(value,x,func=lambda coefficient:(calls.append(coefficient) or coefficient))==value
assert calls==[value]
d=s.Dummy('x',positive=True)
assert s.collect(a*d+b*d,d,evaluate=False)=={d:a+b}
assert not any(name=='sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    subprocess.run([sys.executable, '-c', program], check=True, timeout=15)


def test_high_precision_coefficients_keep_their_native_precision():
    x, a = sx.symbols('x a')
    high = sx.Float('1.25', 90)
    expression = high*x+a*x
    result = sx.collect(expression, x, evaluate=False)
    assert result[x] == high+a
    assert any(atom._prec == high._prec for atom in result[x].args if atom.is_Float)
    assert converted(result) == sp.collect(sx.to_sympy(expression), sp.Symbol('x'), evaluate=False)
