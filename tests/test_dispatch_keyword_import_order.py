"""Installing specialized dispatch must preserve the underlying public keywords."""
import itertools
import json
import subprocess
import sys

import pytest


@pytest.mark.parametrize('modules',list(itertools.permutations([
    'symox.vector','symox.diffgeom','symox.tensor.tensor'])))
def test_specialized_import_orders_keep_named_scalar_and_container_calls(modules):
    source='''
import importlib
import symox as s
for name in MODULES:importlib.import_module(name)
x=s.Symbol('x')
assert s.simplify(expr=(x*x-1)/(x-1),ratio=1)==x+1
assert s.expand(e=(x+1)**2)==x*x+2*x+1
assert s.diff(f=x*x)==2*x
matrix=s.ImmutableSparseMatrix([[s.sin(x)**2+s.cos(x)**2]])
assert s.simplify(expr=matrix)==s.ImmutableSparseMatrix([[1]])
assert s.N(x=s.Matrix([[s.Rational(1,3)]]),n=30).shape==(1,1)
from symox.vector import CoordSys3D
frame=CoordSys3D('C')
assert s.simplify(expr=(s.sin(x)**2+s.cos(x)**2)*frame.i)==frame.i
'''.replace('MODULES',json.dumps(modules))
    result=subprocess.run([sys.executable,'-c',source],capture_output=True,text=True,timeout=20)
    assert result.returncode==0,result.stderr
