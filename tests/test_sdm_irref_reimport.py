import os
import subprocess
import sys
import pytest


@pytest.mark.parametrize("name", ["min", "sorted", "set", "dict", "enumerate", "len", "defaultdict"])
@pytest.mark.parametrize("container", ["dict", "sdm"])
def test_reimported_module_preserves_its_own_consumed_callbacks(name, container):
    source = r"""
import builtins, collections, importlib, sys
import symox as s
old = importlib.import_module('symox.polys.matrices.sdm')
a = {0: {0: s.QQ(1), 1: s.QQ(2)}, 1: {0: s.QQ(3), 1: s.QQ(4)}}
assert old._native_irref(a) is not None
del sys.modules['symox.polys.matrices.sdm']
new = importlib.import_module('symox.polys.matrices.sdm')
assert new is not old
name, kind = sys.argv[1:]
original = getattr(collections if name == 'defaultdict' else builtins, name)
events = []
def changed(*args, **kwargs):
    events.append(name)
    return original(*args, **kwargs)
matrix = a if kind == 'dict' else new.SDM(a, (2, 2), s.QQ)
setattr(new, name, changed)
events.clear()
actual = new.sdm_irref(matrix)
observed = events[:]
new._native_irref = lambda *args: None
events.clear()
expected = new.sdm_irref(matrix)
legacy = events[:]
assert actual == expected
assert observed == legacy and observed, (observed, legacy)
# A second namespace must neither overwrite nor disable the original proof.
assert old._native_irref(a) is not None
"""
    subprocess.run([sys.executable, "-c", source, name, container],
                   check=True, env=os.environ.copy(), timeout=20)
