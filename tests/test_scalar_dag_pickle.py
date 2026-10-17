"""Scalar pickles retain exact native payloads and physical DAG sharing."""
import base64
import json
import pickle
from pathlib import Path
import subprocess
import sys

import pytest
import symox as s


def values():
    x = s.Symbol('x', positive=True)
    y = s.Symbol('y', real=False, integer=False)
    f = s.Function('sin', nargs=(1, 2), real=True)
    return [
        x, y, s.Integer(10**100), s.Rational(10**80 + 1, 10**60 + 3),
        s.Float(1.25), s.Float('0.12345678901234567890123456789', 80),
        s.pi, s.E, s.I, s.EulerGamma, s.oo, -s.oo, s.zoo, s.nan,
        s.Add(y, x, evaluate=False), s.Mul(-1, x+y, evaluate=False),
        s.Pow(x+y, 0, evaluate=False), s.UnevaluatedExpr(x+y),
        f(x), f(x, y), s.sin(y),
        s.Derivative(f(y), y, evaluate=False), s.Subs(f(y), y, x),
        s.Limit(s.sin(y)/y, y, 0, dir='+-'),
        s.Piecewise((x, s.Eq(y, 1)), (y, True)),
        s.And(s.Eq(y, 1), s.Eq(y, 2)),
        s.Interval(0, x), s.FiniteSet(x, y), s.Wild('w', exclude=[x+y]),
        s.Order(y**3, (y, 0)), s.gamma(y), s.factorial(y),
    ]


def records(value):
    if isinstance(value, s.Set):
        value = value._e
    restore, args = value.__reduce__()
    assert restore.__name__ == '_restore_matrix_dag'
    return args


def test_builtin_expression_heads_compare_and_hash_by_constructor():
    x, y = s.symbols('x y')
    assert s.sin(x).func == s.sin(y).func
    assert s.sin(x).func != s.cos(x).func
    assert hash(s.sin(x).func) == hash(s.sin(y).func)
    assert s.pi.func == s.pi.func
    assert s.pi.func != s.E.func
    assert s.sin(x).func(*s.sin(y).args) == s.sin(y)
    assert s.pi.func() == s.pi


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL + 1))
def test_exact_scalar_records_and_public_metadata_roundtrip(protocol):
    for value in values():
        before = records(value)
        restored = pickle.loads(pickle.dumps(value, protocol=protocol))
        assert restored == value and hash(restored) == hash(value)
        assert records(restored) == before
        assert restored.args == value.args
        assert restored.free_symbols == value.free_symbols
        assert restored.func == value.func
        if getattr(value, 'is_Float', False):
            assert restored._prec == value._prec


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL + 1))
def test_shared_scalar_dag_pickle_and_operation_count_do_not_expand_tree(protocol):
    # A timeout contains an accidental return to exponential recursion. The
    # assertion is exact: a binary tree of depth 80 has 2**80-1 operations,
    # but the native shared representation contains only 81 physical nodes.
    source = '''
import pickle,sys
import symox as s
protocol=int(sys.argv[1])
x=s.Symbol('x',real=True)
pair=s.Function('pair',nargs=2)
value=x
for _ in range(80): value=pair(value,value)
assert value.count_ops()==2**80-1
assert s.count_ops(value)==2**80-1
before=value.__reduce__()[1]
assert len(before[1])==81
for tag,canonical,fields in before[1][1:]:
    assert tag in ('s.undefined', 's.undefined.facts')
    assert fields[-1][0]==fields[-1][1]
payload=pickle.dumps(value,protocol=protocol)
assert len(payload)<50000
restored=pickle.loads(payload)
assert restored==value and hash(restored)==hash(value)
assert restored.__reduce__()[1]==before
assert restored.count_ops()==2**80-1
assert restored.free_symbols=={x}
'''
    result = subprocess.run([sys.executable, '-c', source, str(protocol)],
                            capture_output=True, text=True, timeout=15)
    assert result.returncode == 0, result.stderr


def test_scalar_dag_pickle_preserves_distinct_alpha_equal_subs_presentations():
    x, y, t = s.symbols('x y t')
    f = s.Function('f')
    first, second = s.Subs(f(x), x, t), s.Subs(f(y), y, t)
    assert first == second and first.args != second.args
    value = s.Function('pair')(first, second)
    restored = pickle.loads(pickle.dumps(value))
    assert records(restored) == records(value)
    assert restored.args[0].args != restored.args[1].args
    assert restored.args[0].args[0] == f(x)
    assert restored.args[1].args[0] == f(y)


def test_scalar_codec_runs_without_sympy_and_preserves_function_head_metadata():
    source = '''
import importlib.abc,pickle,sys
class DenySympy(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,*args,**kwargs):
        if fullname.split('.')[0]=='sympy':raise AssertionError('SymPy import: '+fullname)
sys.meta_path.insert(0,DenySympy())
import symox as s
x=s.Symbol('x',nonnegative=False)
f=s.Function('sin',nargs=(1,2),real=True)
value=s.Subs(s.Derivative(f(x),x,evaluate=False),x,s.Float('0.1234567890123456789',70))
restored=pickle.loads(pickle.dumps(value))
assert restored==value
assert restored.__reduce__()[1]==value.__reduce__()[1]
assert restored.args[0].args[0].func==f
assert not any(n.split('.')[0]=='sympy' for n in sys.modules)
'''
    result = subprocess.run([sys.executable, '-c', source], capture_output=True,
                            text=True, timeout=10)
    assert result.returncode == 0, result.stderr


def test_legacy_recursive_scalar_pickles_remain_readable():
    # Captured from the installed release at 1a54afc before replacing its
    # reducers. These bytes exercise the old helpers, not the new encoder.
    x = s.Symbol('x', positive=True)
    f = s.Function('f', real=True, nargs=1)
    expected = dict(
        symbol=x, rational=s.Rational(10**60+1, 10**40+3),
        mpfloat=s.Float('0.123456789012345678901234567890123456789', 70),
        builtin=s.sin(x), held_mul=s.Mul(-1, x+1, evaluate=False),
        held_pow=s.Pow(x+1, 0, evaluate=False), undefined=f(x),
        derivative=s.Derivative(f(x), x, evaluate=False), subs=s.Subs(f(x), x, 2),
        logic=s.And(s.Eq(x, 1), s.Eq(x, 2)),
        piecewise=s.Piecewise((x, s.Eq(x, 1)), (0, True)),
        limit=s.Limit(s.sin(x)/x, x, 0),
    )
    fixture = json.loads((Path(__file__).parent/'fixtures/scalar_pickle_v0.json').read_text())
    assert fixture['values'].keys() == expected.keys()
    for name, encoded in fixture['values'].items():
        restored = pickle.loads(base64.b64decode(encoded))
        assert restored == expected[name], name
        assert restored.args == expected[name].args, name
        assert restored.func == expected[name].func, name
        # Legacy reducers recreated equal children independently. Their old
        # sharing cannot be recovered, but the upgraded codec must preserve the
        # resulting graph on every subsequent round trip.
        again = pickle.loads(pickle.dumps(restored))
        assert records(again) == records(restored), name
