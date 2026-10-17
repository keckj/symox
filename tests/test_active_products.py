"""Active products normalize held Mul spines without evaluating opaque children."""
import pickle
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def expression(module, source):
    names = ('Add', 'Mul', 'Pow', 'Float', 'Rational', 'Integer', 'sin', 'cos',
             'exp', 'Function', 'UnevaluatedExpr', 'Symbol', 'Dummy', 'I', 'oo', 'zoo', 'nan')
    namespace = {name: getattr(module, name) for name in names}
    namespace.update(zip(('x', 'y', 'z'), module.symbols('x y z')))
    namespace['p'] = module.Symbol('p', positive=True)
    namespace['n'] = module.Symbol('n', negative=True)
    return module.sympify(eval(source, {'__builtins__': {}}, namespace))


def assert_same_tree(actual, expected):
    # Converting a canonical parent can otherwise hide a changed native tree.
    assert len(actual.args) == len(expected.args), (actual, expected)
    if not actual.args:
        assert sx.to_sympy(actual) == expected
        return
    assert actual.func.__name__ == expected.func.__name__, (actual, expected)
    for left, right in zip(actual.args, expected.args):
        assert_same_tree(left, right)


@pytest.mark.parametrize('source', [
    '2*Mul(x,2,evaluate=False)', 'Mul(Mul(x,2,evaluate=False),z)',
    'Mul(x,2,evaluate=False)*y', 'Mul(y,x,evaluate=False)*x',
    'Mul(Mul(Mul(x,2,evaluate=False),y,evaluate=False),z)',
    'Mul(Mul(x,x,evaluate=False),x)',
    'Mul(Mul(2,Add(x,y,evaluate=False),evaluate=False),z)',
    'Mul(Mul(2,x+y,evaluate=False),z)',
    'Mul(x,2,evaluate=False)*Mul(y,3,evaluate=False)',
    'Mul(x,Rational(1,2),evaluate=False)*Mul(y,6,evaluate=False)',
    'Mul(x,Float("2.5"),evaluate=False)*2',
    '2*Mul(x,Float("2.5",80),evaluate=False)',
    '2*Mul(x,0,evaluate=False)', '0*Mul(x,2,evaluate=False)',
    '0*Mul(x,oo,evaluate=False)', 'Mul(0,Mul(x,oo,evaluate=False))',
    '0*(x*oo)', '0*Mul(zoo,x,evaluate=False)',
    '0*Mul(x,nan,evaluate=False)', '2*Mul(oo,0,evaluate=False)',
    'Mul(Mul(x,-oo,evaluate=False),-2)',
    '2*Mul(x,-oo,evaluate=False)', '-(2*Mul(x,-oo,evaluate=False))',
    '2*Mul(p,oo,evaluate=False)', '2*Mul(n,oo,evaluate=False)',
    '2*Pow(zoo,2,evaluate=False)',
    'Mul(Pow(0,1,evaluate=False),Pow(zoo,1,evaluate=False))',
    'Mul(Pow(oo,-1,evaluate=False),Pow(oo,1,evaluate=False))',
    'Mul(Pow(oo,0,evaluate=False),Pow(oo,0,evaluate=False))',
    'Mul(Mul(x,zoo,evaluate=False),y)',
    '0*Function("f")(oo)', '0*UnevaluatedExpr(oo)',
    '0*sin(oo)', '0*Add(oo,x,evaluate=False)',
    '0*Add(oo,-oo,evaluate=False)',
    '2*Pow(oo,1,evaluate=False)', '2*Pow(-oo,1,evaluate=False)',
    '2*Pow(zoo,1,evaluate=False)', '2*Pow(nan,1,evaluate=False)',
    '0*Pow(oo,0,evaluate=False)', '0*Pow(oo,-1,evaluate=False)',
    '0*Pow(oo,1,evaluate=False)', '0*Pow(zoo,1,evaluate=False)',
    '0*Pow(nan,1,evaluate=False)', '0*Pow(0,-1,evaluate=False)',
    '2*Pow(Mul(x,oo,evaluate=False),1,evaluate=False)',
    '0*Pow(Mul(x,oo,evaluate=False),1,evaluate=False)',
    '2*Pow(Mul(x,oo,evaluate=False),0,evaluate=False)',
    '0*Pow(Mul(x,oo,evaluate=False),0,evaluate=False)',
    '0*Pow(Mul(oo,0,evaluate=False),-1,evaluate=False)',
    'Mul(x,UnevaluatedExpr(y),evaluate=False)*2',
    '2*Function("Mul")(x,2)',
    '2*Mul(Dummy("x",dummy_index=77),2,evaluate=False)',
    '2*Mul(Symbol("x",real=True),2,evaluate=False)',
])
def test_active_product_matches_exact_sympy_tree(source):
    value = expression(sx, source)
    expected = expression(sp, source)
    assert_same_tree(value, expected)
    assert sx.to_sympy(value) == expected, source
    initial_srepr = sx.srepr(value)
    restored = pickle.loads(pickle.dumps(value))
    assert restored == value
    assert sx.srepr(restored) == initial_srepr
    assert sx.to_sympy(restored) == expected


@pytest.mark.parametrize('source', [
    'Mul(x,2,evaluate=False)', 'Mul(2,x,1,evaluate=False)',
    'Mul(x,y,-2,evaluate=False)', 'Mul(x,oo,evaluate=False)',
    'Mul(oo,0,evaluate=False)', 'Mul(x,UnevaluatedExpr(y),evaluate=False)',
    '2*Pow(zoo,1,evaluate=False)',
    '2*Pow(Mul(x,oo,evaluate=False),1,evaluate=False)',
])
def test_identity_and_unary_minus_keep_their_distinct_contracts(source):
    value, expected = expression(sx, source), expression(sp, source)
    initial = value.__reduce__()[1]
    for identity in (1*value, value*1, sx.Mul(1,value), sx.Mul(value,1,1)):
        assert identity.args == value.args
        assert sx.to_sympy(identity) == expected
    assert_same_tree(-value,-expected)
    assert sx.to_sympy(-value) == -expected
    assert sx.to_sympy((-1)*value) == (-1)*expected
    assert sx.to_sympy(sx.Mul(-1,value)) == sp.Mul(-1,expected)
    assert value.__reduce__()[1] == initial


def test_seeded_finite_held_products_match_oracle_and_reconstruct():
    randomizer = random.Random(250791)
    atoms = ['x','y','z','1','-1','2','-2','Rational(1,2)','I','sin(x)','x+y']
    def generate(depth):
        if not depth or randomizer.random() < .4:
            return randomizer.choice(atoms)
        children = [generate(depth-1) for _ in range(randomizer.choice([2,3]))]
        return f'Mul({",".join(children)},evaluate=False)'
    for _ in range(150):
        left, right = generate(2), generate(2)
        source = f'Mul({left},{right})'
        actual, expected = expression(sx, source), expression(sp, source)
        assert_same_tree(actual, expected)
        assert sx.to_sympy(actual) == expected, source
        assert sx.expand(actual-sx.from_sympy(expected)) == 0, source


def test_depth80_active_product_certificate_without_sympy_runtime():
    program = r'''
import builtins,pickle,sys,time
original_import=builtins.__import__
def guarded(name,*args,**kwargs):
    if name=='sympy' or name.startswith('sympy.'):
        raise AssertionError('SymPy runtime fallback')
    return original_import(name,*args,**kwargs)
builtins.__import__=guarded
import symox as s
one=s.Integer(1)
x=s.Dummy('x')
raw=x
for _ in range(80):
    raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
initial=raw.__reduce__()[1]
expected=2*x**(2**80)
assert 2*raw==expected
assert s.Mul(2,raw)==expected
assert raw*2==expected
assert (1*raw).args==raw.args
assert raw.__reduce__()[1]==initial
huge=s._core._restore_arithmetic('Pow',one,[(s.Integer(2),s.Integer(1000000000))],False)
assert 0*huge==0
zero=s._core._restore_arithmetic('Pow',one,[(s.Integer(0),one)],False)
assert zero*huge==0
assert pickle.loads(pickle.dumps(2*raw))==expected
# Repeated occurrences of the same literal numeric/zero products are folded.
raw=one
for _ in range(80):raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
assert 2*raw==2
raw=s.Mul(0,s.oo,evaluate=False)
for _ in range(80):raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
assert (2*raw)==s.nan
assert not any(name=='sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    subprocess.run([sys.executable,'-c',program], check=True, timeout=15)


LATE_COMPLEX_INFINITY_CASES = [
    'Mul(Pow(zoo,1,evaluate=False),oo)',
    'Mul(zoo,Pow(zoo,1,evaluate=False),x)',
    'Mul(Pow(0,-1,evaluate=False),oo,x)',
    'Mul(Pow(zoo,1,evaluate=False),-oo)',
]


def evaluated_tree(value):
    if not value.args:
        return value
    return value.func(*(evaluated_tree(arg) for arg in value.args))


@pytest.mark.parametrize('source', LATE_COMPLEX_INFINITY_CASES)
def test_late_complex_infinity_values_have_the_same_evaluated_tree(source):
    actual, expected = expression(sx, source), expression(sp, source)
    assert evaluated_tree(sx.to_sympy(actual)) == evaluated_tree(expected)


@pytest.mark.xfail(strict=True, reason='Late zoo powers still use a different nonfinite collection phase than SymPy')
@pytest.mark.parametrize('source', LATE_COMPLEX_INFINITY_CASES)
def test_remaining_late_complex_infinity_collection_shapes(source):
    # The value certificate is checked separately, outside this xfail.
    actual, expected = expression(sx, source), expression(sp, source)
    assert_same_tree(actual, expected)


def test_indeterminate_zero_times_infinity_cannot_be_overwritten_by_zoo():
    # SymPy 1.14 currently returns zoo for this held collection order. The
    # explicit 0*oo factor is indeterminate and must remain NaN throughout.
    power = sx.Pow(sx.oo, 1, evaluate=False)
    assert sx.Mul(0, power, sx.zoo) == sx.nan


@pytest.mark.parametrize('reverse', [False, True])
def test_scalar_power_zero_times_explicit_infinite_sum_stays_indeterminate(reverse):
    zero = sx.Pow(0, 1, evaluate=False)
    infinite_sum = sx.Add(sx.oo, sx.Symbol('x'), evaluate=False)
    factors = [infinite_sum, zero] if reverse else [zero, infinite_sum]
    assert sx.Mul(*factors) == sx.nan
    # A power boundary prevents redispatching its composite base.
    assert sx.Mul(zero, sx.Pow(infinite_sum, 1, evaluate=False)) == 0
