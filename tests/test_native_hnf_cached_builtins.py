"""Function execution consumes its cached builtins, even after globals restore."""
import json
import os
import subprocess
import sys

import pytest


SCRIPT = r'''
import builtins
import importlib
import json
import sys
import types
import symox as sx
from symox.polys.matrices.domainmatrix import DomainMatrix

target, raises = sys.argv[1], sys.argv[2] == 'True'
nf = importlib.import_module('symox.polys.matrices.normalforms')
printing = importlib.import_module('symox.polys.domains._element_print')
events, admissions, audit_events = [], [], []
owner = nf if target == 'add_columns' else printing
original = owner.__dict__[target] if target != 'stock' else nf.add_columns
namespace = original.__globals__
saved = namespace['__builtins__']
foreign = dict(builtins.__dict__)

def observed_range(*args):
    events.append('range')
    if raises:
        raise RuntimeError('foreign cached range')
    return []

def observed_type(*args):
    events.append('type')
    if raises:
        raise RuntimeError('foreign cached type')
    return builtins.type(*args)

if target != 'stock':
    foreign['range' if target == 'add_columns' else 'type'] = observed_range if target == 'add_columns' else observed_type
    namespace['__builtins__'] = foreign
    try:
        clone = types.FunctionType(original.__code__, namespace, original.__name__, original.__defaults__, original.__closure__)
    finally:
        namespace['__builtins__'] = saved
    owner.__dict__[target] = clone
    assert clone.__builtins__ is foreign

# A real fresh bootstrap must authenticate the actual functions it encounters.
sys.modules.pop('symox.polys.matrices._hnf_native')
guard = importlib.import_module('symox.polys.matrices._hnf_native')
guard.register(nf.__dict__)
nf._native_hnf = guard.try_hnf
native = sx._core._ddm_hnf

def record(*args):
    result = native(*args)
    admissions.append(result)
    return result

sx._core._ddm_hnf = record

def audit(event, args):
    if event == 'object.__getattr__':
        audit_events.append(args[1])
        raise AssertionError('added audit ' + args[1])

sys.addaudithook(audit)

def execute():
    events.clear()
    admissions.clear()
    source = DomainMatrix([[sx.ZZ(2), sx.ZZ(3)]], (1, 2), sx.ZZ)
    try:
        result = nf._hermite_normal_form(source)
        value = {'rows': [[int(v) for v in row] for row in result.to_list()], 'shape': list(result.shape)}
    except Exception as exc:
        value = {'error': type(exc).__name__, 'message': str(exc)}
    return value, events[:], admissions[:]

actual, actual_events, actual_admissions = execute()
nf._native_hnf = lambda *args: None
expected, expected_events, expected_admissions = execute()
assert actual == expected, (target, raises, actual, expected)
assert actual_events == expected_events, (target, raises, actual_events, expected_events)
assert not expected_admissions
if target == 'stock':
    assert actual_admissions == [1]
    assert actual == {'rows': [[1]], 'shape': [1, 1]}
else:
    assert not actual_admissions
    assert actual_events
    if raises:
        assert actual['error'] == 'RuntimeError'
    elif target == 'add_columns':
        assert actual == {'rows': [[3]], 'shape': [1, 1]}
assert audit_events == []
assert not any(name == 'sympy' or name.startswith('sympy.') for name in sys.modules)
print(json.dumps({'target': target, 'raises': raises, 'actual': actual, 'events': actual_events,
                  'native_admissions': actual_admissions, 'audit_events': audit_events, 'reference_imported': False}))
'''


@pytest.mark.parametrize('target', ['add_columns', '_int_congruent'])
@pytest.mark.parametrize('raises', [False, True])
def test_fresh_hnf_helper_preserves_actual_cached_builtin_execution(target, raises):
    run = subprocess.run([sys.executable, '-c', SCRIPT, target, str(raises)],
                         text=True, capture_output=True, env=os.environ.copy())
    assert run.returncode == 0, run.stdout + run.stderr
    data = json.loads(run.stdout)
    assert data['reference_imported'] is False
    assert data['audit_events'] == []


def test_fresh_hnf_helper_still_admits_stock_under_rejecting_audit_hook():
    run = subprocess.run([sys.executable, '-c', SCRIPT, 'stock', 'False'],
                         text=True, capture_output=True, env=os.environ.copy())
    assert run.returncode == 0, run.stdout + run.stderr
    assert json.loads(run.stdout)['native_admissions'] == [1]
