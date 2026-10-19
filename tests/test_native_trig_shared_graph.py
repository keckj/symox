"""Native trig passes preserve shared graphs and exact identities."""
import subprocess
import sys


def test_native_trig_pass_on_shared_graphs_keeps_capabilities_without_reference_import():
    script = '''
import sys
class BlockReference:
 def find_spec(self,fullname,path=None,target=None):
  if fullname=='sympy' or fullname.startswith('sympy.'):
   raise ImportError('reference forbidden')
sys.meta_path.insert(0,BlockReference())
import symox as s
x=s.Symbol('x')
f=s.Function('shared')
unchanged=s.sin(x)
identity=s.sin(x)**2+s.cos(x)**2
expected=s.Integer(1)
for _ in range(64):
 unchanged=f(unchanged,unchanged)
 identity=f(identity,identity)
 expected=f(expected,expected)
assert s._core.trigsimp(unchanged)==unchanged
assert s._core.trigsimp(identity)==expected
assert not any(name=='sympy' or name.startswith('sympy.') for name in sys.modules)
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=10)
