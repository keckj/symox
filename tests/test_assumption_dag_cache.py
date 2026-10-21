"""DAG-sized assumptions and closed numerical evaluations."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

PREDICATES = (
    'is_real', 'is_positive', 'is_negative', 'is_nonnegative', 'is_nonpositive',
    'is_zero', 'is_nonzero', 'is_integer', 'is_rational', 'is_even', 'is_odd', 'is_finite',
)


def shared_mul(base, depth=80):
    # Reconstruct the exact held payload rather than canonicalizing the product.
    one = sx.Integer(1)
    for _ in range(depth):
        base = sx._core._restore_arithmetic('Mul', one, [(base, one), (base, one)], False)
    return base


@pytest.mark.parametrize('source', [
    'pi-E', 'sqrt(2)', 'sqrt(2)-2', 'sqrt(2)*I', '1+sqrt(2)*I',
    '(-2)**Rational(1,3)', 'exp(2)', 'sin(2)', 'cos(2)', 'log(2)',
    'Rational(2,3)', 'Float("1.25",80)',
])
def test_cached_results_preserve_existing_predicates_and_conversions(source):
    actual, expected = sx.sympify(source), sp.sympify(source)
    fields = ('is_real', 'is_positive', 'is_negative', 'is_nonnegative', 'is_nonpositive', 'is_zero', 'is_nonzero')
    before = actual.__reduce__()[1]
    for _ in range(3):
        for name in fields:
            assert getattr(actual, name) == getattr(expected, name), (source, name)
    assert actual.__reduce__()[1] == before
    assert sx.from_sympy(sx.to_sympy(actual)) == actual


@pytest.mark.parametrize('numeric_first', [False, True])
def test_numeric_and_exact_query_banks_are_separate_through_collect(numeric_first):
    x = sx.Symbol('x')
    difference = sx.pi-sx.E
    if numeric_first:
        assert difference.is_positive is True
    seen = []
    def observe(coefficient):
        # Composite generators are private, assumption-preserving dummies.
        assert coefficient.is_Dummy
        seen.append(coefficient.is_positive)
        return coefficient
    actual = sx.collect(x*difference, [x, difference], observe, evaluate=False)
    assert actual == {x: difference}
    assert seen == [None]
    assert difference.is_positive is True
    for _ in range(4):
        assert sx.collect(x*difference, [x, difference], observe, evaluate=False) == actual
        assert difference.is_positive is True
    assert seen == [None]*5


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL+1))
def test_warmed_caches_do_not_change_scalar_identity_masks_or_pickle(protocol):
    original = shared_mul(sx.Symbol('p', positive=True), 20)
    before_hash, before_args = hash(original), original.__reduce__()[1]
    expected = [getattr(original, field) for field in PREDICATES]
    for _ in range(3):
        assert [getattr(original, field) for field in PREDICATES] == expected
    assert hash(original) == before_hash
    assert original.__reduce__()[1] == before_args
    restored = pickle.loads(pickle.dumps(original, protocol=protocol))
    assert restored == original and hash(restored) == before_hash
    assert restored.__reduce__()[1] == before_args
    assert [getattr(restored, field) for field in PREDICATES] == expected
    assert {original: 'present'}[restored] == 'present'


def test_newton_environments_do_not_reuse_closed_numeric_results():
    x = sx.Symbol('x')
    for seed in (1.0, 3.0, 1.0, 3.0):
        root = sx.nsolve(x*x-2, x, seed)
        assert abs(float(root)-2**0.5) < 1e-12
    # The same expression is still symbolic after these numerical environments.
    assert (x*x-2).free_symbols == {x}


def test_depth80_assumptions_numeric_refinement_and_collect_without_sympy():
    program = r'''
import builtins,pickle,sys
original_import=builtins.__import__
def guarded(name,*args,**kwargs):
    if name=='sympy' or name.startswith('sympy.'):
        raise AssertionError('SymPy fallback')
    return original_import(name,*args,**kwargs)
builtins.__import__=guarded
import symox as s
one=s.Integer(1)
def dag(base):
    for _ in range(80):base=s._core._restore_arithmetic('Mul',one,[(base,one),(base,one)],False)
    return base
fields=('is_real','is_positive','is_negative','is_nonnegative','is_nonpositive','is_zero','is_nonzero','is_integer','is_rational','is_even','is_odd','is_finite')
r=s.Symbol('r',real=True)
value=dag(r)
expected=(True,None,None,None,None,None,None,None,None,None,None,True)
for _ in range(8):assert tuple(getattr(value,f) for f in fields)==expected
p=s.Symbol('p',positive=True)
positive=dag(p)
assert positive.is_positive is True and positive.is_nonzero is True
assert s.collect(s.Symbol('z'),positive,evaluate=False)=={one:s.Symbol('z')}
closed=dag(one)
assert float(closed)==1.0 and complex(closed)==complex(1.0,0.0)
assert closed.is_real is True and closed.is_positive is True
# Gating numeric_refine must itself share the DAG, before evaluation starts.
zero=s._core._restore_arithmetic('Add',s.Integer(0),[(closed,one),(closed,s.Integer(-1))],False)
assert zero.is_positive is False and zero.is_zero is True
assert float(zero)==0.0
ambiguous=one
for _ in range(80):
    ambiguous=s._core._restore_arithmetic('Add',s.Integer(0),[(ambiguous,one),(ambiguous,s.Integer(-1))],False)
assert ambiguous.is_positive is False and ambiguous.is_zero is True and float(ambiguous)==0.0
before=value.__reduce__()[1]
for _ in range(3):
    restored=pickle.loads(pickle.dumps(value))
    assert restored==value and hash(restored)==hash(value)
    assert restored.__reduce__()[1]==before
    assert tuple(getattr(restored,f) for f in fields)==expected
assert not any(name=='sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    subprocess.run([sys.executable, '-c', program], check=True, timeout=15)
