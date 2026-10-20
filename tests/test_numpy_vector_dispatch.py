"""Array acceleration must retain precision, ownership and VM fallbacks."""
import subprocess
import sys

import numpy as np
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize("source", [
    "sin(x)*cos(y)+sin(x/3)*cos(y/7)",
    "exp(x/3)+log(y+3)","x**2+y/3","sqrt(x+3)",
])
@pytest.mark.parametrize("broadcast", [False,True])
def test_elementary_float64_arrays_match_reference(source,broadcast):
    x,y=sx.symbols("x y");X,Y=sp.symbols("x y")
    actual=sx.lambdify((x,y),sx.sympify(source),"numpy")
    expected=sp.lambdify((X,Y),sp.sympify(source),"numpy")
    a=np.linspace(-1,1,15).reshape(3,5)
    b=np.linspace(0,1,15).reshape(3,5)
    if broadcast:a,b=a[:,:1],b[:1,:]
    result=actual(a,b)
    np.testing.assert_allclose(result,expected(a,b),rtol=2e-14,atol=2e-15)
    assert result.dtype==np.float64 and result.flags.writeable
    assert actual.nargs==2
    assert type(actual).__name__=="CompiledFunction"


def test_constant_identity_and_mixed_scalar_outputs_keep_independent_writable_arrays():
    x,y=sx.symbols("x y")
    values=np.arange(12,dtype=np.float64).reshape(3,4)
    for expression,arguments,expected in [
        (x,(values,values),values.copy()),
        (sx.Integer(3),(values,values),np.full_like(values,3)),
        (x,(3.0,values),np.full_like(values,3)),
    ]:
        function=sx.lambdify((x,y),expression,"numpy")
        result=function(*arguments)
        np.testing.assert_array_equal(result,expected)
        assert not np.shares_memory(result,values)
        assert result.dtype==np.float64 and result.flags.writeable


def test_repeated_arrays_use_current_inputs_and_handle_noncontiguous_views():
    x=sx.Symbol("x");function=sx.lambdify(x,sx.sin(x)+x*x,"numpy")
    data=np.linspace(-1,1,100).reshape(10,10)
    view=data[::-2,::2];view.flags.writeable=False
    first=function(view)
    np.testing.assert_allclose(first,np.sin(view)+view*view,rtol=2e-14,atol=2e-15)
    data[:]=2
    second=function(view)
    np.testing.assert_allclose(second,np.sin(view)+view*view,rtol=2e-14,atol=2e-15)
    assert not np.array_equal(first,second)
    assert second.flags.c_contiguous and second.flags.writeable


def test_float32_and_integer_inputs_keep_vm_float64_precision():
    x=sx.Symbol("x");expression=sx.sin(x)+x/7
    actual=sx.lambdify(x,expression,"numpy")
    original=sx._native_lambdify(x,expression)
    for dtype in [np.float32,np.int64]:
        values=np.arange(20,dtype=dtype)
        result=actual(values)
        assert result.dtype==np.float64
        np.testing.assert_array_equal(result,original(values))
    assert isinstance(actual(0.25),float)
    np.testing.assert_array_equal(actual([0.25,0.5]),original([0.25,0.5]))


def test_equal_flat_size_vm_extension_and_scalar_domain_behavior_are_retained():
    x,y=sx.symbols("x y");expression=x+y
    function=sx.lambdify((x,y),expression,"numpy")
    a=np.arange(6,dtype=float).reshape(2,3);b=np.arange(6,dtype=float)
    np.testing.assert_array_equal(function(a,b),a+b.reshape(2,3))
    singular=sx.lambdify((x,y),1/x,"numpy")
    with np.errstate(all="raise"):
        assert np.isinf(singular(0.0,a)).all()
    empty=np.empty((0,3),dtype=float)
    assert function(empty,empty).shape==(0,3)


def test_opaque_special_functions_and_large_shared_graphs_keep_native_evaluation():
    x=sx.Symbol("x");values=np.linspace(0.25,1,9)
    for expression in [sx.gamma(x)]:
        f=sx.lambdify(x,expression,"numpy");baseline=sx._native_lambdify(x,expression)
        np.testing.assert_array_equal(f(values),baseline(values))
    expression=x
    for _ in range(40):expression=sx.sin(expression)+sx.cos(expression)
    f=sx.lambdify(x,expression,"numpy");baseline=sx._native_lambdify(x,expression)
    np.testing.assert_array_equal(f(values),baseline(values))


def test_numpy_codegen_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
import symox as sx
import numpy as np
x,y=sx.symbols('x y')
f=sx.lambdify((x,y),sx.sin(x)*sx.cos(y),'numpy')
a=np.linspace(0,1,12).reshape(3,4)
np.testing.assert_allclose(f(a,a),np.sin(a)*np.cos(a))
assert f(a[:,:1],a[:1,:]).shape==(3,4)
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,"-c",source],check=True,timeout=20)
