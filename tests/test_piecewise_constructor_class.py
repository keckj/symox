"""The public Piecewise constructor is usable as a native type predicate."""
import inspect
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


def test_piecewise_is_a_constructor_class_and_instance_predicate():
    x=sx.Symbol('x');p=sx.Piecewise((x,x>0),(0,True))
    assert inspect.isclass(sx.Piecewise)
    assert isinstance(p,sx.Piecewise)
    assert not isinstance(x,sx.Piecewise)
    assert p.func is sx.Piecewise
    assert p.func(*p.args)==p
    assert sx.Piecewise.fromiter(p.args)==p
    assert sx.Piecewise.class_key()==sp.Piecewise.class_key()
    assert repr(sx.Piecewise)==repr(sp.Piecewise)


def test_constructor_alias_and_unbound_methods_share_the_same_class():
    from symox.functions.elementary.piecewise import Piecewise
    assert Piecewise is sx.Piecewise
    x=sx.Symbol('x');p=Piecewise((x,x>0),(0,True))
    assert Piecewise.subs(p,x,2)==2
    assert p.has(Piecewise)
    assert p.replace(Piecewise,lambda *branches:branches[0][0])==x


def test_empty_constructor_rejects_while_false_branches_return_nan():
    with pytest.raises(TypeError):sp.Piecewise()
    with pytest.raises(TypeError):sx.Piecewise()
    x=sx.Symbol('x')
    assert sx.Piecewise((x,False))==sx.nan


def test_evaluation_context_is_preserved():
    x=sx.Symbol('x')
    with sx.evaluate(False):
        p=sx.Piecewise((x,x>0),(0,True))
        assert isinstance(p,sx.Piecewise)
        assert len(p.args)==2
        assert sx.Piecewise((x,True))==x
    assert sx.Piecewise((x,True))==x


def test_constructor_type_checks_need_no_reference_imports():
    script='''
import sys
class NoReference:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoReference())
import symox as sx
from symox.functions.elementary.piecewise import Piecewise
x=sx.Symbol('x');p=Piecewise((x,x>0),(0,True))
assert Piecewise is sx.Piecewise
assert isinstance(p,Piecewise)
assert p.func(*p.args)==p
assert Piecewise.subs(p,x,2)==2
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=20)


@pytest.mark.parametrize('branches', [
    [(1,False),(2,'positive'),(3,True),(4,'negative')],
    [(1,False)],[(1,True)],
    [(1,'positive'),(2,'positive'),(3,False),(4,True)],
])
def test_held_constructor_filters_branches_like_reference(branches):
    def build(module):
        x=module.Symbol('x')
        conditions={'positive':x>0,'negative':x<0}
        return module.Piecewise(*[(value,conditions.get(cond,cond)) for value,cond in branches],evaluate=False)
    actual,reference=build(sx),build(sp)
    if isinstance(reference,sp.Piecewise):
        assert isinstance(actual,sx.Piecewise)
        assert tuple((pair.expr.to_sympy(),pair.cond.to_sympy()) for pair in actual.args)==tuple(reference.args)
    else:
        assert actual.to_sympy()==reference
