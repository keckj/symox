"""Dummy atoms carry an explicit identity, including across interop and codecs."""
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sm


def test_dummy_identity_name_assumptions_and_constructor_class():
    first=sm.Dummy('x',positive=True);second=sm.Dummy('x',positive=True)
    assert first!=second and first.dummy_index!=second.dummy_index
    assert isinstance(first,sm.Dummy) and isinstance(first,sm.Symbol)
    assert first.is_Dummy and first.is_Symbol and first.is_positive
    assert first.name=='x' and str(first)=='_x' and first.args==()
    assert first.func is sm.Dummy
    restored=first.func(first.name,dummy_index=first.dummy_index,positive=True)
    assert restored==first and hash(restored)==hash(first)
    assert first.func(*first.args)!=first  # Like SymPy, atom.args does not encode its identity.
    assert sm.Dummy('y',dummy_index=first.dummy_index,positive=True)!=first
    assert sm.Dummy('x',dummy_index=first.dummy_index)!=first
    assert len({first,second,restored})==2
    assert len({first:1,second:2,restored:3})==2


def test_underscore_names_never_create_dummy_identity():
    ordinary=sm.Symbol('_x');dummy=sm.Dummy('x')
    assert str(ordinary)==str(dummy)
    assert ordinary.class_name=='Symbol' and not ordinary.is_Dummy
    assert not isinstance(ordinary,sm.Dummy) and ordinary.func is sm.Symbol
    assert ordinary!=dummy and ordinary!=sm.Symbol('x')
    with pytest.raises(AttributeError):ordinary.dummy_index
    assert pickle.loads(pickle.dumps(ordinary))==ordinary
    assert sm.from_sympy(sp.Symbol('_x'))==ordinary


def test_anonymous_dummy_and_symbols_constructor_are_fresh():
    a,b=sm.Dummy(),sm.Dummy()
    assert a!=b and a.name.startswith('Dummy_') and b.name.startswith('Dummy_')
    assert str(a)=='_'+a.name
    ds=sm.symbols('x:3',cls=sm.Dummy,integer=True)
    assert [d.name for d in ds]==['x0','x1','x2']
    assert all(isinstance(d,sm.Dummy) and d.is_integer for d in ds)
    with pytest.raises((ValueError,TypeError)):sm.Dummy(dummy_index=1)
    with pytest.raises(ValueError):sm.Dummy('x',positive=True,negative=True)


@pytest.mark.parametrize('index',[-7,0,10**90])
@pytest.mark.parametrize('flags',[{},dict(real=True),dict(positive=True),dict(nonnegative=False)])
def test_dummy_srepr_pickle_sympy_roundtrip(index,flags):
    actual=sm.Dummy("a'b\\c",dummy_index=index,**flags)
    expected=sp.Dummy("a'b\\c",dummy_index=index,**flags)
    # Quoting style may differ, but both printers encode the identical atom.
    # sympify(srepr) turns dummy_index into sympy.Integer; SymPy Dummy itself
    # distinguishes that from a Python int. Evaluating srepr is the exact oracle.
    assert eval(actual.srepr(), {"Dummy":sp.Dummy,"__builtins__":{}})==expected
    assert sm.Dummy("plain",dummy_index=index,**flags).srepr()==sp.srepr(sp.Dummy("plain",dummy_index=index,**flags))
    assert sm.sympify(actual.srepr())==actual
    assert sm.to_sympy(actual)==expected
    assert sm.from_sympy(expected)==actual
    assert sm.from_sympy(sm.to_sympy(actual))==actual
    restored=pickle.loads(pickle.dumps(actual))
    assert restored==actual and restored.dummy_index==index
    assert restored.name==actual.name and restored.assumptions0==actual.assumptions0


def test_independent_sympy_dummies_remain_independent_on_import():
    first,second=sp.Dummy('x'),sp.Dummy('x')
    a,b=sm.from_sympy(first),sm.from_sympy(second)
    assert a!=b and a.dummy_index==first.dummy_index and b.dummy_index==second.dummy_index
    assert sm.from_sympy(first)==a
    expression=sm.from_sympy(first+2*second)
    assert expression.free_symbols=={a,b}
    assert expression.xreplace({a:1})==1+2*b
    assert sm.to_sympy(expression)==first+2*second


def test_dummy_arithmetic_subs_diff_and_numpy_do_not_capture_same_names():
    a,b=sm.Dummy('x'),sm.Dummy('x');x=sm.Symbol('x');underscore=sm.Symbol('_x')
    expression=a+2*b+x+underscore
    assert expression.free_symbols=={a,b,x,underscore}
    assert expression.has(a) and expression.has(b)
    assert expression.subs(a,3)==3+2*b+x+underscore
    assert expression.xreplace({a:b})==3*b+x+underscore
    assert sm.diff(expression,a)==1 and sm.diff(expression,b)==2
    assert sm.lambdify((a,b,x,underscore),expression,'numpy')(1,2,3,4)==12
    assert np.dot(np.array([a,b],dtype=object),np.array([2,3],dtype=object))==2*a+3*b


def test_bound_subs_preserves_dummy_scope_and_free_points():
    d,e=sm.Dummy('x'),sm.Dummy('x');x=sm.Symbol('x');f=sm.Function('f')
    held=sm.Subs(f(d)+e,d,x)
    assert held.free_symbols=={e,x}
    assert held.doit()==f(x)+e
    assert held.xreplace({e:2}).doit()==f(x)+2
    assert held==sm.Subs(f(sm.Symbol('bound'))+e,sm.Symbol('bound'),x)
    restored=pickle.loads(pickle.dumps(held))
    assert restored.variables==(d,) and restored.args==held.args
    assert restored.free_symbols=={e,x}
    assert sm.from_sympy(sm.to_sympy(held))==held
    # A new Dummy allocated for the chain-rule binder must not capture this user symbol.
    user=sm.Symbol('_xi_1')
    differentiated=sm.diff(f(x*x+user),x)
    assert differentiated.free_symbols=={x,user}
    assert sm.diff(f(x*x+user),x)==differentiated
    assert sm.from_sympy(sm.to_sympy(differentiated))==differentiated


def test_dummy_order_and_printers_follow_sympy_with_explicit_indices():
    a=sm.Dummy('a',dummy_index=17);b=sm.Dummy('a',dummy_index=9);z=sm.Symbol('z')
    sa=sp.Dummy('a',dummy_index=17);sb=sp.Dummy('a',dummy_index=9);sz=sp.Symbol('z')
    assert (a+b+z).args==(z,b,a)
    assert sm.default_sort_key(a)==sp.default_sort_key(sa)
    for expression,expected in [(a,sa),(a+b+z,sa+sb+sz)]:
        assert str(expression)==str(expected)
        assert sm.cxxcode(expression)==sp.cxxcode(expected,standard='c++17')
    assert sm.latex(a)==sp.latex(sa) and sm.pretty(a)==sp.pretty(sa)
    from symox.printing.cxx import CXX17CodePrinter
    class Custom(CXX17CodePrinter):
        def _print_Symbol(self,value):return 'custom_'+value.name
    assert Custom().doprint(a)=='a_17'
    assert Custom().doprint(z)=='custom_z'


def test_restored_high_indices_cannot_be_reused_by_fresh_construction():
    high=10**110
    restored=sm.Dummy('x',dummy_index=high)
    fresh=sm.Dummy('x')
    assert fresh!=restored and fresh.dummy_index!=high
    imported=sm.from_sympy(sp.Dummy('x',dummy_index=high+10))
    assert sm.Dummy('x').dummy_index!=imported.dummy_index


def test_dummy_native_operations_do_not_import_sympy():
    program=r'''
import sys,pickle
class BlockSympy:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('native Dummy operation imported SymPy')
sys.meta_path.insert(0,BlockSympy())
import symox as s
a,b=s.Dummy('x'),s.Dummy('x'); x=s.Symbol('_x')
e=a+b+x
assert a!=b and e.free_symbols=={a,b,x}
assert e.xreplace({a:1})==1+b+x
assert pickle.loads(pickle.dumps(e))==e
assert s.sympify(a.srepr())==a
assert s.diff(e,a)==1
assert s.cxxcode(a)=='x_'+str(a.dummy_index)
assert not x.is_Dummy
'''
    subprocess.run([sys.executable,'-c',program],check=True,capture_output=True,text=True,timeout=15)


def test_periodic_solution_dummy_cannot_capture_a_user_parameter():
    x=sm.Symbol('x');parameter=sm.Symbol('_n')
    result=sm.solveset(sm.exp(x+parameter)-1,x,sm.Complexes)
    assert result.free_symbols=={parameter}
    converted=result.to_sympy()
    assert converted.free_symbols=={sp.Symbol('_n')}
    assert converted.lamda.variables[0].is_Dummy
    assert sm.from_sympy(converted)==result._e


def test_sympy_lambda_import_preserves_dummy_binders_and_multiple_variables():
    d=sp.Dummy('x',integer=True)
    actual=sm.from_sympy(sp.Lambda(d,d+1))
    native=sm.from_sympy(d)
    assert actual.args==(sm.Tuple(native),native+1)
    assert actual.free_symbols==set()
    assert sm.to_sympy(actual)==sp.Lambda(d,d+1)
    x,y=sp.symbols('x y')
    multiple=sm.from_sympy(sp.Lambda((x,y),x+y))
    assert multiple.variables==sm.symbols('x y')
    assert sm.to_sympy(multiple)==sp.Lambda((x,y),x+y)


def test_native_and_sympy_fresh_generators_have_independent_namespaces():
    program=r"""
import pickle,sympy as p,symox as s
values=[]
for _ in range(30):
    foreign=p.Dummy('same',real=True)
    imported=s.from_sympy(foreign)
    native=s.Dummy('same',real=True)
    values.extend((imported,native))
    for protocol in range(pickle.HIGHEST_PROTOCOL+1):
        assert pickle.loads(pickle.dumps(imported,protocol=protocol))==imported
        assert pickle.loads(pickle.dumps(native,protocol=protocol))==native
assert len(set(values))==len(values)
assert len({value.dummy_index for value in values})==len(values)
"""
    subprocess.run([sys.executable,'-c',program],check=True,capture_output=True,text=True,timeout=15)


def test_dummy_generator_resets_process_namespace_after_posix_fork():
    import os
    if not hasattr(os,'fork'):
        pytest.skip('requires POSIX fork')
    program=r"""
import os,pickle,symox as s
initial=s.Dummy('same')  # initialize the state before fork
read_end,write_end=os.pipe()
pid=os.fork()
if pid==0:
    os.close(read_end)
    try:
        payload=pickle.dumps(s.Dummy('same'))
        os.write(write_end,payload)
    finally:
        os.close(write_end)
    os._exit(0)
os.close(write_end)
parent=s.Dummy('same')
payload=b''
while True:
    chunk=os.read(read_end,65536)
    if not chunk:break
    payload+=chunk
os.close(read_end)
_,status=os.waitpid(pid,0)
assert status==0
child=pickle.loads(payload)
assert len({initial,parent,child})==3
assert len({value.dummy_index for value in (initial,parent,child)})==3
assert s.Dummy('same') not in {initial,parent,child}
"""
    subprocess.run([sys.executable,'-c',program],check=True,capture_output=True,text=True,timeout=15)


@pytest.mark.parametrize('source,count',[
    ('sin(x+1)',2),
    ('2*cos(x)**2-1',4),
    ('sin(x+1)/(x-1)',2),
])
def test_periodic_branches_share_a_fresh_dummy_with_sympy_presentation(source,count):
    x=sm.Symbol('x');X=sp.Symbol('x')
    first=sm.solveset(sm.sympify(source),x,sm.Reals)
    second=sm.solveset(sm.sympify(source),x,sm.Reals)
    expected=sp.solveset(sp.sympify(source),X,sp.S.Reals)
    assert str(first)==str(expected)
    assert str(second)==str(expected)
    def binders(result):
        converted=result.to_sympy()
        return [node.lamda.variables[0] for node in sp.preorder_traversal(converted)
                if isinstance(node,sp.ImageSet)]
    variables=binders(first)
    assert len(variables)==count
    assert all(variable==variables[0] for variable in variables)
    assert variables[0].is_Dummy
    assert variables[0]!=binders(second)[0]
    assert first.free_symbols==set()
