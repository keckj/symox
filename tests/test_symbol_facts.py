"""Exact Symbol facts preserve declarations, closure, identity and bound scopes."""
import itertools
import pickle
import random
import subprocess
import sys

import pytest
import sympy as sp
from sympy.core.assumptions import _assume_defined
import symox as sx


FACT_NAMES = sorted(_assume_defined)
SINGLE = [{name: value} for name in FACT_NAMES for value in (False, True)
          if (name, value) != ('commutative', False)]
COMBINATIONS = [
    {}, {'real': None}, {'positive': True, 'real': True},
    {'finite': False, 'extended_real': True},
    {'extended_positive': True, 'finite': False},
    {'extended_nonnegative': True, 'finite': False},
    {'zero': False, 'nonzero': False}, {'real': False, 'finite': True},
    {'even': True, 'positive': True, 'prime': False},
    {'integer': True, 'odd': False}, {'integer': True, 'even': False},
    {'complex': True, 'imaginary': False}, {'algebraic': True, 'real': False},
]


def make(module, kind, kwargs, name='x'):
    return getattr(module, kind)(name, **({'dummy_index': 73} if kind == 'Dummy' else {}), **kwargs)


@pytest.mark.parametrize('kind', ['Symbol', 'Dummy'])
@pytest.mark.parametrize('kwargs', SINGLE+COMBINATIONS)
def test_declared_and_derived_symbol_facts_match_sympy(kind, kwargs):
    actual, expected = make(sx, kind, kwargs), make(sp, kind, kwargs)
    assert actual.assumptions0 == expected.assumptions0
    for fact in FACT_NAMES:
        assert getattr(actual, 'is_'+fact) is getattr(expected, 'is_'+fact), (fact, kwargs)
    assert actual.func is getattr(sx, kind)
    reconstructed = actual.func(actual.name,
        **({'dummy_index': actual.dummy_index} if kind == 'Dummy' else {}),
        **actual.assumptions0)
    assert reconstructed == actual and hash(reconstructed) == hash(actual)
    assert sx.to_sympy(actual) == expected
    imported = sx.from_sympy(expected)
    assert imported == actual and hash(imported) == hash(actual)
    assert imported.assumptions0 == expected.assumptions0
    assert sx.sympify(sx.srepr(actual)) == actual
    # SymPy's string sympifier turns dummy_index into sympy.Integer, which
    # is not interchangeable with int inside Dummy equality. Evaluate the
    # srepr constructors directly, as in test_dummy_compat.py.
    assert eval(sx.srepr(actual), {'Symbol': sp.Symbol, 'Dummy': sp.Dummy, '__builtins__': {}}) == expected
    for protocol in range(pickle.HIGHEST_PROTOCOL+1):
        restored = pickle.loads(pickle.dumps(actual, protocol=protocol))
        assert restored == actual and hash(restored) == hash(actual)
        assert restored.assumptions0 == expected.assumptions0
        assert sx.srepr(restored) == sx.srepr(actual)


def test_equivalent_and_distinct_declarations_have_coherent_dict_identity():
    variants = [{}, {'prime': False}, {'finite': False}, {'complex': False},
                {'imaginary': False}, {'odd': False}, {'extended_real': True},
                {'real': True}, {'positive': True}, {'real': True, 'positive': True}]
    native = [sx.Symbol('x', **kw) for kw in variants]
    oracle = [sp.Symbol('x', **kw) for kw in variants]
    assert len(set(native)) == len(set(oracle))
    for i, j in itertools.product(range(len(variants)), repeat=2):
        assert (native[i] == native[j]) == (oracle[i] == oracle[j])
        if native[i] == native[j]:
            assert hash(native[i]) == hash(native[j])
    assert sx.Symbol('x', real=None) == sx.Symbol('x')


@pytest.mark.parametrize('kwargs', [{'real': True}, {'real': False}, {'positive': True},
                                  {'imaginary': True}, {'finite': True}])
def test_existing_undefined_function_assumptions_do_not_consume_symbol_closure(kwargs):
    x = sx.Symbol('x')
    head = sx.Function('f', nargs=1, **kwargs)
    assert head == sx.Function('f', nargs=(1,), **kwargs)
    value = head(x)
    assert value.func == head and hash(value.func) == hash(head)
    # Native function heads normalize accepted arities to tuples. SymPy keeps
    # int and tuple constructor metadata distinct despite equal nargs sets.
    assert sx.to_sympy(value) == sp.Function('f', nargs=(1,), **kwargs)(sp.Symbol('x'))
    assert pickle.loads(pickle.dumps(head)) == head
    assert pickle.loads(pickle.dumps(value)).func == head
    for name in kwargs:
        assert getattr(value, 'is_'+name) is getattr(sp.Function('f', **kwargs)(sp.Symbol('x')), 'is_'+name)


def test_seeded_mixed_fact_closures_and_contradictions_match_sympy():
    rng = random.Random(184009)
    names = [name for name in FACT_NAMES if name != 'commutative']
    for _ in range(250):
        kwargs = {name: rng.choice([False, True]) for name in rng.sample(names, rng.choice([2, 3, 4]))}
        try:
            expected = sp.Symbol('x', **kwargs)
        except ValueError:
            with pytest.raises(ValueError):
                sx.Symbol('x', **kwargs)
        else:
            actual = sx.Symbol('x', **kwargs)
            assert actual.assumptions0 == expected.assumptions0, kwargs


@pytest.mark.parametrize('kwargs', [
    {'real': True, 'finite': False}, {'integer': True, 'even': False, 'odd': False},
    {'algebraic': True, 'transcendental': True}, {'positive': True, 'negative': True},
])
def test_inconsistent_symbol_and_dummy_facts_fail_before_construction(kwargs):
    for kind in ('Symbol', 'Dummy'):
        with pytest.raises(ValueError, match='inconsistent'):
            make(sx, kind, kwargs)


@pytest.mark.parametrize('kwargs', [{'commutative': None}, {'commutative': 'yes'}, {'custom_fact': True}])
def test_unsupported_symbol_domains_and_fact_names_are_explicit(kwargs):
    with pytest.raises((ValueError, NotImplementedError)):
        sx.Symbol('x', **kwargs)


def test_noncommutative_symbol_facts_are_accepted():
    a = sx.Symbol('a', commutative=False)
    assert a.is_commutative is False
    facts = a.assumptions0
    assert facts['commutative'] is False and facts['complex'] is False
    assert a.is_real is False and a.is_zero is False


def test_basic_argument_order_uses_full_symbol_contents_and_dummy_index():
    variants = [{}, {'positive': True}, {'negative': True}, {'integer': True},
                {'prime': False}, {'extended_real': True}, {'real': False}]
    for kind in ('Symbol', 'Dummy'):
        native = [make(sx, kind, kwargs, name=name) for name, kwargs in zip('zxqpabr', variants)]
        oracle = [make(sp, kind, kwargs, name=name) for name, kwargs in zip('zxqpabr', variants)]
        for constructor in ('Add', 'Mul'):
            actual = getattr(sx, constructor)(*native)
            expected = getattr(sp, constructor)(*oracle)
            assert tuple(sx.to_sympy(arg) for arg in actual.args) == expected.args


def test_same_name_distinct_facts_do_not_capture_substitution_or_free_symbols():
    first = sx.Symbol('x', prime=False)
    second = sx.Symbol('x', finite=False)
    body = sx.Function('pair')(first, second)
    assert body.free_symbols == {first, second}
    assert body.xreplace({first: 3}).args == (3, second)
    held = sx.Lambda(first, body)
    assert held.free_symbols == {second}
    assert held(4).args == (4, second)
    restored = pickle.loads(pickle.dumps(held))
    assert restored == held and restored.free_symbols == {second}


def test_held_arithmetic_preserves_symbol_declarations_and_public_arguments():
    native = [sx.Symbol('x', prime=False), sx.Symbol('x', finite=False)]
    oracle = [sp.Symbol('x', prime=False), sp.Symbol('x', finite=False)]
    for constructor in ('Add', 'Mul'):
        value = getattr(sx, constructor)(native[1], native[0], native[1], evaluate=False)
        expected = getattr(sp, constructor)(oracle[1], oracle[0], oracle[1], evaluate=False)
        assert sx.to_sympy(value).args == expected.args
        assert sx.from_sympy(expected).args == value.args
        for protocol in range(6):
            restored = pickle.loads(pickle.dumps(value, protocol=protocol))
            assert restored.args == value.args and sx.srepr(restored) == sx.srepr(value)


@pytest.mark.parametrize('kwargs', [
    {'finite': False}, {'infinite': True}, {'finite': False, 'extended_real': True},
    {'finite': False, 'extended_positive': True}, {'finite': False, 'extended_negative': True},
    {'finite': False, 'extended_real': False}, {'finite': True}, {'zero': False, 'nonzero': False},
])
@pytest.mark.parametrize('coefficient', ['zero', 'oo', 'negative_oo', 'zoo'])
@pytest.mark.parametrize('held', ['symbol', 'product', 'power'])
@pytest.mark.parametrize('reverse', [False, True])
def test_exact_extended_symbol_facts_control_nonfinite_products(kwargs, coefficient, held, reverse):
    def build(module):
        symbol = module.Symbol('x', **kwargs)
        factor = (module.Mul(2, symbol, evaluate=False) if held == 'product' else
                  module.Pow(symbol, 1, evaluate=False) if held == 'power' else symbol)
        scalar = {'zero': module.Integer(0), 'oo': module.oo,
                  'negative_oo': -module.oo, 'zoo': module.zoo}[coefficient]
        factors = [scalar, factor]
        if reverse:
            factors.reverse()
        return module.Mul(*factors)
    actual, expected = build(sx), build(sp)
    assert sx.to_sympy(actual) == expected
    assert tuple(sx.to_sympy(arg) for arg in actual.args) == expected.args
    assert pickle.loads(pickle.dumps(actual)).args == actual.args


@pytest.mark.parametrize('kwargs, equation, expected', [
    ({'prime': False}, lambda x: x*x-4, [-2]),
    ({'even': False}, lambda x: x*x-4, []),
    ({'odd': False}, lambda x: x*x-1, []),
    ({'imaginary': False}, lambda x: x*x+1, []),
    ({'finite': False}, lambda x: x-1, []),
    ({'transcendental': True}, lambda x: x*x-2, []),
])
def test_solver_checks_new_facts_only_with_exact_certificates(kwargs, equation, expected):
    x = sx.Symbol('x', **kwargs)
    assert sx.solve(equation(x), x) == expected
    assert sx.solve(equation(x), x, check=False)


def test_legacy_symbol_dummy_codec_records_remain_readable():
    from symox._core import _restore_matrix_dag
    x = _restore_matrix_dag(1, [('s.symbol', True, ('x', 2, False, []))], 0)
    assert x == sx.Symbol('x', positive=True)
    dummy = _restore_matrix_dag(1, [('s.dummy', True, ('x', 0, 73))], 0)
    assert dummy == sx.Dummy('x', dummy_index=73)


@pytest.mark.parametrize('tag, fields', [
    ('s.symbol.facts', ('x', (1 << 31, 0), False, [])),
    ('s.symbol.facts', ('x', (0, 1), False, [])),
    ('s.dummy.facts', ('x', (0, 1), 73)),
])
def test_fact_codec_rejects_invalid_masks(tag, fields):
    with pytest.raises(ValueError, match='fact masks'):
        sx._core._restore_matrix_dag(1, [(tag, True, fields)], 0)


def test_shared_held_symbol_fact_dag_without_sympy_runtime():
    program = r'''
import builtins,pickle,sys
original=builtins.__import__
def guarded(name,*args,**kwargs):
    if name=='sympy' or name.startswith('sympy.'):raise AssertionError('SymPy fallback')
    return original(name,*args,**kwargs)
builtins.__import__=guarded
import symox as s
x=s.Symbol('x',prime=False)
y=s.Dummy('x',finite=False,dummy_index=73)
f=s.Function('pair')
e=f(x,y)
for _ in range(80):e=f(e,e)
assert e.free_symbols=={x,y}
for protocol in range(6):
    restored=pickle.loads(pickle.dumps(e,protocol=protocol))
    assert restored==e and hash(restored)==hash(e)
    assert restored.free_symbols=={x,y}
assert x.is_prime is False and y.is_finite is False
assert s.Symbol('r',extended_real=True).is_real is None
assert s.Symbol('r',extended_real=True,finite=False).is_real is False
one=s.Integer(1)
algebraic=s.Symbol('a',algebraic=True)
for _ in range(80):
    algebraic=s._core._restore_arithmetic('Mul',one,[(algebraic,one),(algebraic,one)],False)
assert algebraic.is_algebraic is True and algebraic.is_transcendental is False
huge=s._core._restore_arithmetic('Mul',one,[(s.Integer(2),s.Integer(10**9))],False)
assert huge.is_algebraic is True and huge.is_transcendental is False
nonfinite=s.Symbol('nf',finite=False)
for _ in range(80):
    nonfinite=s._core._restore_arithmetic('Mul',one,[(nonfinite,one),(one,one)],False)
assert s.Mul(0,nonfinite)==s.nan and s.Mul(nonfinite,0)==s.nan
assert not any(name=='sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    subprocess.run([sys.executable, '-c', program], check=True, timeout=15)
