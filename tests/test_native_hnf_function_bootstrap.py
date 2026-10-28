"""Fresh guard imports must authenticate consumed mathematical helper code."""
import json
import os
import subprocess
import sys

import pytest


SCRIPT = r'''
import importlib
import json
import sys
import types
import symox as sx
from symox.polys.matrices.domainmatrix import DomainMatrix

target, raises = sys.argv[1], sys.argv[2] == 'True'
nf = importlib.import_module('symox.polys.matrices.normalforms')
original = nf.__dict__[target]
saved = types.FunctionType(original.__code__, original.__globals__, original.__name__,
                           original.__defaults__, original.__closure__)
events, admissions, audit_events = [], [], []
nf.__dict__['_hnf_boot_events'] = events
nf.__dict__['_hnf_boot_raises'] = raises
nf.__dict__['_hnf_boot_original'] = saved

def observed_gcdex(a, b):
    _hnf_boot_events.append('gcdex')
    if _hnf_boot_raises:
        raise RuntimeError('observed gcdex')
    return _hnf_boot_original(a, b)

def observed_columns(m, i, j, a, b, c, d):
    _hnf_boot_events.append('columns')
    if _hnf_boot_raises:
        raise RuntimeError('observed columns')
    return _hnf_boot_original(m, i, j, a, b, c, d)

# Keep the original function identity and true globals/builtins. Only its
# mathematical code changes before a real fresh guard bootstrap.
original.__code__ = (observed_gcdex if target == '_gcdex' else observed_columns).__code__
sys.modules.pop('symox.polys.matrices._hnf_native')
guard = importlib.import_module('symox.polys.matrices._hnf_native')
guard.register(nf.__dict__)
nf._native_hnf = guard.try_hnf
kernel = sx._core._ddm_hnf

def record(*args):
    result = kernel(*args)
    admissions.append(result)
    return result

sx._core._ddm_hnf = record

def audit(event, args):
    if event == 'object.__getattr__':
        audit_events.append(args[1])
        raise AssertionError('added metadata audit ' + args[1])

sys.addaudithook(audit)

def execute():
    events.clear()
    admissions.clear()
    source = DomainMatrix.from_list([[2, 3]], sx.ZZ)
    try:
        result = nf._hermite_normal_form(source)
        value = {'rows': [[int(v) for v in row] for row in result.to_list()]}
    except Exception as exc:
        value = {'error': type(exc).__name__, 'message': str(exc)}
    return value, events[:], admissions[:]

actual, actual_events, actual_admissions = execute()
nf._native_hnf = lambda *args: None
expected, expected_events, expected_admissions = execute()
assert actual == expected, (target, raises, actual, expected)
assert actual_events == expected_events, (target, raises, actual_events, expected_events)
assert expected_events == ['gcdex' if target == '_gcdex' else 'columns']
assert actual_admissions == expected_admissions == []
assert audit_events == []
assert not any(name == 'sympy' or name.startswith('sympy.') for name in sys.modules)
print(json.dumps({'target': target, 'raises': raises, 'actual': actual, 'events': actual_events,
                  'native_admissions': actual_admissions, 'audit_events': audit_events,
                  'reference_imported': False}))
'''


@pytest.mark.parametrize('target', ['_gcdex', 'add_columns'])
@pytest.mark.parametrize('raises', [False, True])
def test_fresh_hnf_guard_preserves_changed_helper_code_and_its_real_execution(target, raises):
    result = subprocess.run([sys.executable, '-c', SCRIPT, target, str(raises)],
                            text=True, capture_output=True, env=os.environ.copy())
    assert result.returncode == 0, result.stdout + result.stderr
    data = json.loads(result.stdout)
    assert data['reference_imported'] is False
    assert data['audit_events'] == []
