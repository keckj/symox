"""Full, exact undefined-head metadata, independently of Symbol defaults."""
import base64
import itertools
import json
import pickle
from pathlib import Path
import random
import subprocess
import sys

import pytest
import sympy as sp
from sympy.core.assumptions import _assume_defined
import symox as sx


FACT_NAMES = sorted(_assume_defined)
DECLARATIONS = [
    {name: value} for name in FACT_NAMES for value in (False, True)
    if (name, value) != ('commutative', False)
] + [
    {}, {'real': None}, {'positive': True, 'real': True},
    {'finite': False, 'extended_real': True},
    {'extended_positive': True, 'finite': False},
    {'extended_nonnegative': True, 'finite': False},
    {'zero': False, 'nonzero': False}, {'real': False, 'finite': True},
    {'even': True, 'positive': True, 'prime': False},
    {'integer': True, 'odd': False}, {'integer': True, 'even': False},
    {'complex': True, 'imaginary': False}, {'algebraic': True, 'real': False},
    {'positive': True, 'commutative': True},
]


@pytest.mark.parametrize('kwargs', DECLARATIONS)
def test_all_head_facts_and_applied_queries_match_sympy(kwargs):
    x, X = sx.Symbol('x'), sp.Symbol('x')
    head = sx.Function('f', nargs=(1, 2), **kwargs)
    oracle = sp.Function('f', nargs=(1, 2), **kwargs)
    assert head._assumptions == {k: v for k, v in oracle._kwargs.items() if k != 'nargs'}
    value, expected = head(x), oracle(X)
    assert value.assumptions0 == expected.assumptions0 == {}
    for fact in FACT_NAMES:
        assert getattr(value, 'is_'+fact) is getattr(expected, 'is_'+fact), (kwargs, fact)
    assert value.func == head and hash(value.func) == hash(head)
    assert value.func(*value.args) == value
    assert sx.to_sympy(head) == oracle
    assert sx.to_sympy(value) == expected
    assert sx.from_sympy(oracle) == head
    assert sx.from_sympy(expected) == value
    assert sx.sympify(value.srepr()) == value
    assert sp.sympify(value.srepr()) == expected
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        restored_head = pickle.loads(pickle.dumps(head, protocol=protocol))
        restored = pickle.loads(pickle.dumps(value, protocol=protocol))
        assert restored_head == head and hash(restored_head) == hash(head)
        assert restored == value and hash(restored) == hash(value)
        assert restored.func == head
        assert restored.srepr() == value.srepr()


def test_canonical_head_identity_and_dictionary_keys_preserve_false_facts():
    variants = [{}, {'commutative': True}, {'prime': False}, {'finite': False},
                {'even': False}, {'odd': False}, {'real': True},
                {'real': True, 'commutative': True}, {'real': True, 'finite': True},
                {'extended_real': True}, {'complex': False}, {'imaginary': False}]
    native = [sx.Function('f', **kw) for kw in variants]
    oracle = [sp.Function('f', **kw) for kw in variants]
    x = sx.Symbol('x')
    for i, j in itertools.product(range(len(variants)), repeat=2):
        assert (native[i] == native[j]) == (oracle[i] == oracle[j])
        assert (native[i](x) == native[j](x)) == (oracle[i] == oracle[j])
        if native[i] == native[j]:
            assert hash(native[i]) == hash(native[j])
            assert hash(native[i](x)) == hash(native[j](x))
    assert len(set(native)) == len(set(oracle))
    assert len({f(x) for f in native}) == len(set(oracle))


@pytest.mark.parametrize('kwargs', [{}, {'commutative': True}, {'real': True},
    {'extended_real': True}, {'finite': False}, {'prime': False}, {'hermitian': True}])
@pytest.mark.parametrize('argument', ['symbol', 'tuple', 'unevaluated', 'limit'])
def test_applied_commutativity_is_not_an_implicit_head_declaration(kwargs, argument):
    def apply(module):
        x = module.Symbol('x')
        value = {'symbol': lambda: x, 'tuple': lambda: module.Tuple(x),
                 'unevaluated': lambda: module.UnevaluatedExpr(x),
                 'limit': lambda: module.Limit(x, x, 0)}[argument]()
        return module.Function('f', **kwargs)(value)
    assert apply(sx).is_commutative is apply(sp).is_commutative


def test_seeded_mixed_function_closures_and_contradictions_match_sympy():
    rng = random.Random(415104)
    names = [name for name in FACT_NAMES if name != 'commutative']
    for _ in range(300):
        kwargs = {name: rng.choice([False, True]) for name in rng.sample(names, rng.choice([2, 3, 4]))}
        try:
            expected = sp.Function('f', **kwargs)
        except ValueError:
            with pytest.raises(ValueError):
                sx.Function('f', **kwargs)
        else:
            actual = sx.Function('f', **kwargs)
            assert actual._assumptions == expected._kwargs
            value, oracle = actual(sx.Symbol('x')), expected(sp.Symbol('x'))
            for fact in FACT_NAMES:
                assert getattr(value, 'is_'+fact) is getattr(oracle, 'is_'+fact), (kwargs, fact)


@pytest.mark.parametrize('kwargs', [{'commutative': None}, {'commutative': 'yes'},
    {'custom_fact': True}, {'custom_fact': None}, {'integer': True, 'even': False, 'odd': False},
    {'real': True, 'finite': False}])
def test_unsupported_or_inconsistent_function_facts_are_rejected(kwargs):
    with pytest.raises((ValueError, NotImplementedError)):
        sx.Function('f', **kwargs)
    if 'custom_fact' in kwargs:
        with pytest.raises((ValueError, NotImplementedError)):
            sx.from_sympy(sp.Function('f', **kwargs)(sp.Symbol('x')))


def test_function_arities_and_canonical_metadata_reconstruct_together():
    x, y = sx.symbols('x y')
    for nargs in [None, 0, 1, (1,), (2, 1, 2)]:
        f = sx.Function('f', nargs=nargs, prime=False)
        count = 0 if nargs in [None, 0] else 1
        value = f(*([x] * count))
        assert value.func == f and hash(value.func) == hash(f)
        assert sx.from_sympy(sx.to_sympy(value)) == value
    f = sx.Function('f', nargs=(1, 2), finite=False)
    assert f(x).func == f(x, y).func == f
    with pytest.raises(TypeError):
        f()
    with pytest.raises(TypeError):
        f(x, y, x)
    for bad in [(), -1, (1, -1), (1.5,), 'bad']:
        with pytest.raises((ValueError, TypeError)):
            sx.Function('f', nargs=bad)


def test_has_and_atoms_match_complete_heads_including_foreign_selectors():
    x = sx.Symbol('x')
    declarations = [{}, {'prime': False}, {'finite': False}, {'commutative': True}]
    heads = [sx.Function('f', nargs=(1,), **kw) for kw in declarations]
    value = sx.Function('pair')(*[head(x) for head in heads[:3]])
    for i, head in enumerate(heads):
        expected = {head(x)} if i < 3 else set()
        assert value.atoms(head) == expected
        assert value.has(head) is (i < 3)
        assert value.atoms(head.to_sympy()) == expected
        assert value.has(head.to_sympy()) is (i < 3)
    assert not value.has(sx.Function('f', nargs=(2,), prime=False))


@pytest.mark.parametrize('name', ['sin', 'Tuple', 'UnevaluatedExpr', 'Order',
    'Integral', 'Sum', 'Lambda', 'IdentityFunction', 'FunctionMatrix', "f'\\x\n"])
def test_user_heads_do_not_become_builtins_after_metadata_roundtrips(name):
    x = sx.Dummy('x', prime=False, dummy_index=73)
    f = sx.Function(name, nargs=(1,), finite=False)
    value = f(x)
    assert value.is_UndefinedFunction
    assert sx.sympify(value.srepr()) == value
    assert sx.from_sympy(sx.to_sympy(value)) == value
    assert value.evalf(50).is_UndefinedFunction
    assert value.func == f
    assert value.atoms(f.to_sympy()) == {value}
    assert value.has(f.to_sympy())
    assert sx.diff(value, x).is_Derivative
    assert isinstance(sx.limit(value, x, 0), sx.Limit)
    assert pickle.loads(pickle.dumps(value)).func == f


@pytest.mark.parametrize('kwargs', [{'finite': False}, {'finite': True},
    {'finite': False, 'extended_real': True}, {'finite': False, 'extended_real': False},
    {'finite': False, 'extended_positive': True}, {'finite': False, 'extended_negative': True}])
@pytest.mark.parametrize('scalar', ['zero', 'oo', 'neg_oo', 'zoo'])
@pytest.mark.parametrize('reverse', [False, True])
def test_nonfinite_consumers_use_applied_head_facts(kwargs, scalar, reverse):
    def build(module):
        x = module.Symbol('x')
        value = module.Function('f', **kwargs)(x)
        coefficient = {'zero': module.Integer(0), 'oo': module.oo,
                       'neg_oo': -module.oo, 'zoo': module.zoo}[scalar]
        return module.Mul(*([value, coefficient] if reverse else [coefficient, value]))
    actual, expected = build(sx), build(sp)
    assert sx.to_sympy(actual) == expected
    assert tuple(sx.to_sympy(arg) for arg in actual.args) == expected.args
    assert pickle.loads(pickle.dumps(actual)) == actual
    assert sx.Integer(0) * sx.Function('f')(sx.oo) == 0


def test_facts_do_not_capture_same_name_symbols_or_functions_in_bindings():
    x = sx.Symbol('x', prime=False)
    y = sx.Symbol('x', finite=False)
    first, second = sx.Function('f', prime=False), sx.Function('f', finite=False)
    value = sx.Lambda(x, sx.Tuple(first(x), second(y)))
    assert value.free_symbols == {y}
    assert value(3).args == (first(3), second(y))
    assert first(x).xreplace({second(x): 1}) == first(x)
    assert first(x).subs({x: y}).func == first
    assert pickle.loads(pickle.dumps(value)) == value


def test_function_metadata_survives_matrix_function_and_applyfunc_parents():
    head = sx.Function('entry', nargs=2, finite=False)
    matrix = sx.FunctionMatrix(2, 3, head)
    unary = sx.Function('entry', nargs=1, prime=False)
    mapped = sx.MatrixSymbol('A', 2, 3).applyfunc(unary)
    for value, expected in [(matrix, head), (mapped, unary)]:
        for transformed in [value, pickle.loads(pickle.dumps(value)),
                            sx.sympify(sx.srepr(value)),
                            sx.from_sympy(sx.to_sympy(value))]:
            entry = transformed[1, 2]
            assert entry.func == expected and hash(entry.func) == hash(expected)
            assert entry.is_finite is expected(*entry.args).is_finite


@pytest.mark.parametrize('kwargs', [{'finite': False}, {'infinite': True},
    {'extended_real': True, 'finite': False}])
@pytest.mark.parametrize('held', ['application', 'sum', 'product'])
def test_known_nonfinite_function_values_fail_matrix_preflight(kwargs, held):
    x, t = sx.symbols('x t')
    value = sx.Function('f', **kwargs)(t)
    assert value.is_finite is False
    entry = {'application': lambda: value,
             'sum': lambda: sx.Add(value, t, evaluate=False),
             'product': lambda: sx.Mul(2, value, evaluate=False)}[held]()
    matrix = sx.ImmutableMatrix([[entry]])
    for operation in [matrix.exp, matrix.log, lambda: matrix.analytic_func(sx.exp(x), x),
                      matrix.singular_value_decomposition,
                      lambda: sx.eye(1).analytic_func(sx.Add(x, entry, evaluate=False), x)]:
        with pytest.raises(ValueError, match='finite'):
            operation()


@pytest.mark.parametrize('kwargs', [{}, {'extended_real': True}])
def test_unknown_function_finiteness_keeps_diagonal_matrix_functions_available(kwargs):
    x, t = sx.symbols('x t')
    value = sx.Function('f', **kwargs)(t)
    assert value.is_finite is None
    assert sx.ImmutableMatrix([[value]]).exp() == sx.ImmutableMatrix([[sx.exp(value)]])
    assert sx.ImmutableMatrix([[value]]).log() == sx.ImmutableMatrix([[sx.log(value)]])
    assert sx.eye(1).analytic_func(x+value, x) == sx.ImmutableMatrix([[1+value]])


def test_previous_function_pickles_preserve_legacy_metadata_and_internal_origin():
    fixture = Path(__file__).parent / 'fixtures/function_facts_pickle_v0.json'
    x = sx.Symbol('x')
    expected = {'plain_head': sx.Function('f'),
        'positive_head': sx.Function('f', positive=True, nargs=(1, 2)),
        'plain': sx.Function('f')(x),
        'positive': sx.Function('f', positive=True, nargs=(1, 2))(x),
        'user_integral': sx.Function('Integral')(x),
        'tuple': sx.Tuple(x), 'unevaluated': sx.UnevaluatedExpr(x)}
    for key, payload in json.loads(fixture.read_text()).items():
        restored = pickle.loads(base64.b64decode(payload))
        assert restored == expected[key], key
        assert pickle.loads(pickle.dumps(restored)) == restored


def test_noncommutative_function_facts_are_accepted():
    g = sx.Function('g', commutative=False)
    x = sx.Symbol('x')
    assert g(x).is_commutative is False
    assert g(x)*g(x**2) != g(x**2)*g(x)


@pytest.mark.parametrize('replacement', [
    lambda f: (f[0], (1 << 31, 0), *f[2:]),
    lambda f: (f[0], (0, 1), *f[2:]),
    lambda f: (*f[:2], [], *f[3:]),
    lambda f: (*f[:2], [2], *f[3:]),
    lambda f: (*f[:3], True, f[4]),
])
def test_invalid_function_codec_metadata_is_rejected(replacement):
    value = sx.Function('f', prime=False, nargs=(1,))(sx.Symbol('x'))
    restore, (version, records, root) = value.__reduce__()
    records = list(records)
    tag, canonical, fields = records[root]
    assert tag == 's.undefined.facts'
    records[root] = (tag, canonical, replacement(fields))
    with pytest.raises(ValueError):
        restore(version, records, root)


def test_function_facts_dag_queries_pickle_and_native_only_execution_are_bounded():
    source = '''
import pickle, sys
sys.modules['sympy'] = None
import symox as s
x = s.Dummy('x', prime=False)
f = s.Function('pair', nargs=2, real=True, algebraic=True, prime=False)
value = x
for _ in range(80):
    value = f(value, value)
assert value.is_real is True and value.is_algebraic is True
assert value.is_prime is False and value.is_odd is None
assert value.func == f and hash(value.func) == hash(f)
assert not value.has(s.Function('missing', finite=False))
assert not value.has(s.Function('pair', nargs=2, finite=False))
assert value.free_symbols == {x}
assert len(value.atoms(f)) == 80
assert len(value.__reduce__()[1][1]) == 81
for protocol in range(6):
    restored = pickle.loads(pickle.dumps(value, protocol=protocol))
    assert restored == value and hash(restored) == hash(value)
    assert restored.func == f and restored.is_real is True
assert s.Function('f', finite=False)(x).is_finite is False
assert s.Integer(0)*s.Function('f', finite=False)(x) == s.nan
'''
    result = subprocess.run([sys.executable, '-c', source],
                            capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
