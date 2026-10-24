"""Held scalar derivatives and native derivative-specification validation."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def variables(module,case):
    x,y=module.symbols('x y')
    return {
        'single':(x,), 'integer':(x,2), 'native_integer':(x,module.Integer(2)),
        'pair':((x,2),), 'list':([x,2],), 'native_tuple':(module.Tuple(x,2),),
        'zero':(x,0), 'mixed':(x,y), 'mixed_orders':((x,2),(y,1)),
        'repeat':(x,x,y),
    }[case]


@pytest.mark.parametrize('source',['sin(x)','x**3*y**2','f(x,y)','exp(x*y)','log(x)','2'])
@pytest.mark.parametrize('case',['single','integer','native_integer','pair','list','native_tuple',
    'zero','mixed','mixed_orders','repeat'])
@pytest.mark.parametrize('evaluate',[True,False])
def test_global_and_method_derivatives_match_reference(source,case,evaluate):
    expression=sx.sympify(source);reference=sp.sympify(source,locals={'f':sp.Function('f')})
    args=variables(sx,case);oracle=variables(sp,case)
    expected=sp.diff(reference,*oracle,evaluate=evaluate)
    for result in [sx.diff(expression,*args,evaluate=evaluate),
                   expression.diff(*args,evaluate=evaluate),
                   sx.Derivative(expression,*args,evaluate=evaluate)]:
        actual=result.to_sympy()
        if evaluate:assert actual.expand()==expected.expand()
        else:assert actual==expected


def test_held_derivatives_infer_one_variable_and_evaluate_later():
    x=sx.Symbol('x')
    held=(sx.sin(x)+x*x).diff(evaluate=False)
    assert held==sx.Derivative(sx.sin(x)+x*x,x)
    assert held.doit()==sx.cos(x)+2*x
    assert sx.diff(2,evaluate=False)==0
    assert sx.diff(x*x,x,evaluate=False,simplify=False)==sx.Derivative(x*x,x)
    assert sx.diff(x*x,x,simplify=False)==2*x
    with pytest.raises(TypeError):sx.diff(x*x,x,unknown=True)


@pytest.mark.parametrize('case',['negative','float','boolean','count_first','extra_count',
    'zero_extra','bad_variable','bad_pair','too_large'])
def test_invalid_specs_raise_before_computation(case):
    x=sx.Symbol('x')
    args={
        'negative':(x,-1),'float':(x,2.0),'boolean':(x,True),
        'count_first':(2,),'extra_count':(x,2,3),'zero_extra':(x,0,0),
        'bad_variable':(sx.sin(x),),'bad_pair':((x,2,3),),'too_large':(x,2**32),
    }[case]
    for call in [lambda:sx.diff(sx.sin(x),*args),
                 lambda:sx.sin(x).diff(*args),
                 lambda:sx.Derivative(sx.sin(x),*args)]:
        with pytest.raises((ValueError,TypeError,OverflowError)):call()


def test_nested_held_orders_do_not_wrap_at_u32_boundary():
    x=sx.Symbol('x');f=sx.Function('f')
    first=sx.Derivative(f(x),(x,2**32-1))
    second=sx.Derivative(first,x)
    assert second.derivative_count==2**32
    assert sum(int(n) for _,n in second.variable_count)==2**32
    assert second.expr==f(x)
    assert first.diff(x).derivative_count==2**32


def test_matrix_elements_and_held_partials_are_valid_scalar_variables():
    x=sx.Symbol('x');f=sx.Function('f');A=sx.MatrixSymbol('A',2,2)
    assert (A[0,0]*A[1,1]).diff(A[0,0])==A[1,1]
    assert sx.Derivative(A[0,0],A[0,0],evaluate=True)==1
    partial=sx.Derivative(f(x),x)
    assert (partial*partial).diff(partial)==2*partial
    assert sx.Derivative(partial,partial,evaluate=True)==1
    invalid=sx.Derivative(sx.sin(x),x)
    with pytest.raises(ValueError):(invalid*invalid).diff(invalid)
    nested=sx.Derivative(f(f(x)),x)
    with pytest.raises(ValueError):(nested*nested).diff(nested)


@pytest.mark.parametrize('matrix_type',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
def test_explicit_matrices_can_build_entrywise_held_scalar_derivatives(matrix_type):
    x=sx.Symbol('x');X=sp.Symbol('x');constructor=getattr(sx,matrix_type)
    matrix=constructor([[x*x,1],[sx.sin(x),0]])
    expected=sp.Matrix([[X*X,1],[sp.sin(X),0]]).applyfunc(lambda entry:entry.diff(X,evaluate=False))
    result=matrix.applyfunc(lambda entry:entry.diff(x,evaluate=False))
    assert result.to_sympy()==expected
    with pytest.raises(NotImplementedError):sx.diff(matrix,x,evaluate=False)
    with pytest.raises(NotImplementedError):matrix.diff(x,evaluate=False)
    assert result.doit()==constructor([[2*x,0],[sx.cos(x),0]])
    if matrix_type in ('Matrix','ImmutableMatrix'):
        assert matrix.diff()==constructor([[2*x,0],[sx.cos(x),0]])
        assert constructor([[1,2]]).diff()==constructor([[0,0]])
        with pytest.raises(ValueError):constructor([[x,sx.Symbol('y')]]).diff()


def test_arrays_propagate_scalar_options_and_keep_matrix_jacobians():
    x=sx.Symbol('x');array=sx.Array([x*x,sx.sin(x)])
    assert array.applyfunc(lambda entry:entry.diff(x,evaluate=False))==sx.Array([
        sx.Derivative(x*x,x),sx.Derivative(sx.sin(x),x)])
    with pytest.raises(NotImplementedError):sx.diff(array,x,evaluate=False)
    A=sx.MatrixSymbol('A',2,2)
    assert sx.diff(sx.Trace(A),A)==sx.Identity(2)
    with pytest.raises(NotImplementedError):sx.diff(sx.Trace(A),A,evaluate=False)
    with pytest.raises(NotImplementedError):sx.diff(x*A,x,evaluate=False)


def test_differentiation_options_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x')
held=sin(x).diff(x,Integer(2),evaluate=False)
assert held==Derivative(sin(x),(x,2))
assert held.doit()==-sin(x)
assert Matrix([[x*x,1]]).applyfunc(lambda entry:entry.diff(x,evaluate=False))[0,0]==Derivative(x*x,x)
'''
    subprocess.run([sys.executable,'-c',source],check=True)


@pytest.mark.parametrize('source',['G(x,eta/x)','G(x,x*y)','G(y,x,x**2)','G(x,y)'])
def test_symbol_slot_shared_with_another_argument_uses_subs(source):
    # AppliedUndef.fdiff holds Derivative(G(...), x) only when no other
    # argument depends on x; otherwise the slot goes through a Subs.
    def derivative(module):
        x,y,eta=module.symbols('x y eta')
        expr=module.sympify(source,locals={'G':module.Function('G'),'x':x,'y':y,'eta':eta})
        return module.srepr(expr.diff(x))
    import re
    strip=lambda text:re.sub(r'dummy_index=-?\d+','',text)
    assert strip(derivative(sx))==strip(derivative(sp))
