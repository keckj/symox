"""The actual first public LLL call, without prior classifier/helper warmup."""
import subprocess
import sys
import platform

import pytest

VERIFIED_RUNTIME = (sys.implementation.name == 'cpython'
                    and sys.version_info[:3] in ((3, 10, 12), (3, 10, 20))
                    and platform.system() == 'Linux' and platform.machine() == 'x86_64')


def run(script):
    result = subprocess.run([sys.executable, '-c', script], text=True, capture_output=True, timeout=45)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.skipif(not VERIFIED_RUNTIME, reason='optional public C APIs verified on these two runtimes only')
@pytest.mark.parametrize('entry', ['ddm', 'dense', 'sparse'])
@pytest.mark.parametrize('transform', [False, True])
def test_real_first_public_call_enters_native_after_real_conversion(entry, transform):
    script = r'''
import _abc, sys, types
import symox as s
from symox._basic import Basic
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
from symox.polys.matrices import _lll_native as p
if ENTRY == 'ddm':
    a = DDM([[s.ZZ(2), s.ZZ(0), s.ZZ(1)], [s.ZZ(3), s.ZZ(2), s.ZZ(0)]], (2, 3), s.ZZ)
else:
    a = DomainMatrix.from_list([[2, 0, 1], [3, 2, 0]], s.ZZ)
    if ENTRY == 'sparse':
        a = a.to_sparse()
# No earlier LLL, QQ conversion, isinstance/issubclass or new metadata helper.
assert not any(ref() is s._core.DomainElement for group in _abc._get_dump(Basic)[:3] for ref in group)
native = s._core._ddm_lll_native
calls = []
def spy(*args):
    # These are the actual first QQ conversion results, consumed by the kernel.
    assert [[int(v) for v in row] for row in args[0]] == [[2, 0, 1], [3, 2, 0]]
    assert [int(v) for v in args[4]] == [2, 0, 1]
    assert all(type(v) is s._core.DomainElement and v.parent() == s.QQ for v in args[4])
    result = native(*args)
    calls.append(result is not None)
    return result
s._core._ddm_lll_native = spy
events = []
def audit(event, args):
    events.append(event)
    if event == 'object.__getattr__':
        raise RuntimeError('new sensitive metadata audit')
sys.addaudithook(audit)
if TRANSFORM:
    b, t = a.lll_transform()
else:
    b, t = a.lll(), None
assert calls == [True], calls
assert events == [], events
assert type(p.prefix) is types.BuiltinFunctionType and p.prefix.__self__ is s._core
assert type(p.ready) is types.BuiltinFunctionType and p.ready.__self__ is s._core
assert any(ref() is s._core.DomainElement for ref in _abc._get_dump(Basic)[2])
rows = b.to_list() if hasattr(b, 'to_list') else b
assert [[int(v) for v in row] for row in rows] == [[2, 0, 1], [1, 2, -1]]
if TRANSFORM:
    rows = t.to_list() if hasattr(t, 'to_list') else t
    assert [[int(v) for v in row] for row in rows] == [[1, 0], [-1, 1]]
'''.replace('ENTRY', repr(entry)).replace('TRANSFORM', repr(transform))
    run(script)


@pytest.mark.skipif(not VERIFIED_RUNTIME, reason='optional public C APIs verified on these two runtimes only')
def test_first_native_call_keeps_real_scratch_matrices_and_states():
    run(r'''
import _abc, sys
import symox as s
from symox._basic import Basic
from symox.polys.matrices.ddm import DDM
x = DDM([[s.ZZ(2), s.ZZ(0)], [s.ZZ(3), s.ZZ(2)]], (2, 2), s.ZZ)
assert not any(ref() is s._core.DomainElement for group in _abc._get_dump(Basic)[:3] for ref in group)
native = s._core._ddm_lll_native
observed = []
def spy(*args):
    state = sys._getframe(2).f_locals
    assert state['i'] == 0
    assert state['x'] is x and state['y'] is args[0] and state['T'] is args[3]
    assert state['y_star'][0] is args[4]
    assert [[int(v) for v in row] for row in state['y_star']] == [[2, 0], [0, 0]]
    assert [[int(v) for v in row] for row in state['mu']] == [[0, 0], [0, 0]]
    assert [int(v) for v in state['g_star']] == [0, 0]
    assert state['half'] == s.QQ(1, 2)
    observed.append(True)
    return native(*args)
s._core._ddm_lll_native = spy
b, t = x.lll_transform()
assert observed == [True]
assert [[int(v) for v in row] for row in b] == [[2, 0], [-1, 2]]
''')


@pytest.mark.parametrize('entry', ['_lll_first_supported', '_lll_prefix_admit', '_lll_ready_admit', 'foreign_core_key'])
def test_unavailable_or_untrusted_native_binding_keeps_complete_cold_body_and_existing_hot_route(entry):
    run(r'''
import importlib.abc, importlib.machinery, sys, _abc
calls = []
def untrusted(*args):
    calls.append('untrusted entry')
    raise AssertionError('a foreign entry must never authenticate itself')
class Key:
    def __hash__(self):
        return hash('_lll_first_supported')
    def __eq__(self, other):
        calls.append('foreign equality')
        raise AssertionError('lookup before exact key validation')
class Loader:
    def __init__(self, original):
        self.original = original
    def create_module(self, spec):
        return self.original.create_module(spec)
    def exec_module(self, module):
        self.original.exec_module(module)
        core = module._core
        if ENTRY == 'foreign_core_key':
            del core.__dict__['_lll_first_supported']
            core.__dict__[Key()] = None
        else:
            core.__dict__[ENTRY] = untrusted
class Finder:
    def find_spec(self, fullname, path=None, target=None):
        if fullname != 'symox.polys.matrices._lll_native':
            return None
        spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        spec.loader = Loader(spec.loader)
        return spec
sys.meta_path.insert(0, Finder())
import symox as s
from symox._basic import Basic
from symox.polys.matrices.ddm import DDM
x = DDM([[s.ZZ(2), s.ZZ(0)], [s.ZZ(3), s.ZZ(2)]], (2, 2), s.ZZ)
assert not any(ref() is s._core.DomainElement for group in _abc._get_dump(Basic)[:3] for ref in group)
native = s._core._ddm_lll_native
entered = []
def spy(*args):
    entered.append(True)
    return native(*args)
s._core._ddm_lll_native = spy
b, t = x.lll_transform()
assert calls == [] and entered == []
assert [[int(v) for v in row] for row in b] == [[2, 0], [-1, 2]]
# The original route remains available after its actual first classification.
assert x.lll_transform() == (b, t)
assert calls == [] and entered == [True], (calls, entered)
'''.replace('ENTRY', repr(entry)))
