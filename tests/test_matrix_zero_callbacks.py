"""Custom zero and simplification callbacks participate in native pivots."""
import subprocess
import sys
import warnings

import numpy as np

import pytest
import sympy as sp
import symox as sx


def converted(value):
    if isinstance(value,(list,tuple)):return type(value)(converted(item) for item in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


def zero_test(module,x,mode):
    def predicate(value):
        if mode=='threshold' and value.is_number:
            return abs(complex(value.evalf(25)))<1e-9
        if mode=='symbol':return value.subs(x,0).is_zero
        if mode=='unknown' and value.has(x):return None
        return value.is_zero
    return predicate


@pytest.mark.parametrize('constructor',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('operation',['rref','rank','nullspace','echelon_form','rank_decomposition'])
@pytest.mark.parametrize('source,mode',[
    ('[[1e-12,1],[0,1]]','threshold'),
    ('[[x,1],[0,1]]','symbol'),
    ('[[x,1]]','unknown'),
    ('[[1,2],[3,4]]','exact'),
])
def test_native_pivot_callbacks_match_reference(constructor,operation,source,mode):
    x=sx.Symbol('x');X=sp.Symbol('x')
    matrix=getattr(sx,constructor)(eval(source,{'x':x}))
    foreign=getattr(sp,constructor)(eval(source,{'x':X}))
    original=matrix.copy()
    options={'with_pivots':True} if operation=='echelon_form' else {}
    actual=getattr(matrix,operation)(iszerofunc=zero_test(sx,x,mode),**options)
    expected=getattr(foreign,operation)(iszerofunc=zero_test(sp,X,mode),**options)
    assert converted(actual)==expected
    assert matrix==original


@pytest.mark.parametrize('constructor',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('operation',['rref','rank','nullspace','echelon_form','rank_decomposition','rowspace','columnspace'])
def test_custom_simplification_resolves_unknown_pivots(constructor,operation):
    x=sx.Symbol('x');X=sp.Symbol('x')
    matrix=getattr(sx,constructor)([[x,1],[0,1]])
    foreign=getattr(sp,constructor)([[X,1],[0,1]])
    calls=[]
    def simplify(value):
        assert isinstance(value,sx.Expr)
        calls.append(value)
        return sx.cancel(value.subs(x,0))
    actual=getattr(matrix,operation)(simplify=simplify)
    expected=getattr(foreign,operation)(simplify=lambda value:sp.cancel(value.subs(X,0)))
    assert converted(actual)==expected
    assert x in calls


@pytest.mark.parametrize('normalize_last',[False,True])
@pytest.mark.parametrize('pivots',[False,True])
def test_normalization_flags_work_with_callbacks(normalize_last,pivots):
    matrix=sx.Matrix([[2,3,4],[3,5,2]])
    foreign=sp.Matrix([[2,3,4],[3,5,2]])
    actual=matrix.rref(iszerofunc=lambda value:value.is_zero,normalize_last=normalize_last,pivots=pivots)
    expected=foreign.rref(iszerofunc=lambda value:value.is_zero,normalize_last=normalize_last,pivots=pivots)
    assert converted(actual)==expected


def test_native_boolean_answers_and_tristate_none_are_supported():
    x=sx.Symbol('x');matrix=sx.Matrix([[x,1],[1,0]])
    calls=[]
    def zero(value):
        calls.append(value)
        answer=value.is_zero
        return None if answer is None else sx.true if answer else sx.false
    reduced,pivots=matrix.rref(iszerofunc=zero)
    assert reduced==sx.eye(2) and pivots==(0,1)
    assert x in calls


@pytest.mark.parametrize('wrap',[np.bool_,lambda value:np.array(value,dtype=bool)])
def test_numpy_boolean_scalars_work_without_index_coercion(wrap):
    matrix=sx.Matrix([[1e-12,1],[0,1]])
    zero=lambda value:wrap(abs(complex(value))<1e-9)
    with warnings.catch_warnings():
        warnings.simplefilter('error',DeprecationWarning)
        assert matrix.rank(iszerofunc=zero)==1
        assert matrix.rref(iszerofunc=zero)==(sx.Matrix([[0,1],[0,0]]),(1,))
        with pytest.raises(sx.NonInvertibleMatrixError):matrix.LUsolve(sx.Matrix([1,1]),iszerofunc=zero)
    x=sx.Symbol('x');symbolic=sx.Matrix([[x,1],[1,0]])
    assert symbolic.rref(simplify=wrap(True))==symbolic.rref(simplify=True)


def test_numeric_scalar_items_do_not_become_boolean_callback_answers():
    for value in [np.int64(1),np.float64(0)]:
        with pytest.raises(TypeError):sx.eye(2).rref(iszerofunc=lambda expression:value)


def test_numpy_boolean_callbacks_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
import numpy as np
from symox import Matrix
matrix=Matrix([[1e-12,1],[0,1]])
zero=lambda value:np.isclose(float(value),0,atol=1e-9)
assert matrix.rank(iszerofunc=zero)==1
assert matrix.rref(iszerofunc=zero)==(Matrix([[0,1],[0,0]]),(1,))
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=15)


@pytest.mark.parametrize('operation',['rref','rank','nullspace','echelon_form','rank_decomposition'])
def test_zero_callback_errors_are_propagated_immediately(operation):
    class Marker(Exception):pass
    failure=Marker('callback stopped')
    calls=[]
    def zero(value):
        calls.append(value)
        raise failure
    matrix=sx.Matrix([[sx.Symbol('x'),1],[0,1]])
    with pytest.raises(Marker) as caught:getattr(matrix,operation)(iszerofunc=zero)
    assert caught.value is failure and len(calls)==1


@pytest.mark.parametrize('operation',['rref','rank','nullspace','echelon_form','rank_decomposition','rowspace','columnspace'])
def test_simplifier_errors_keep_the_original_exception(operation):
    class Marker(Exception):pass
    failure=Marker('simplification stopped')
    def simplify(value):raise failure
    matrix=sx.Matrix([[sx.Symbol('x'),1],[0,1]])
    with pytest.raises(Marker) as caught:getattr(matrix,operation)(simplify=simplify)
    assert caught.value is failure


def test_invalid_callback_values_and_noncallables_raise():
    matrix=sx.eye(2)
    with pytest.raises(TypeError):matrix.rref(iszerofunc=42)
    with pytest.raises(TypeError):matrix.rref(simplify=42)
    with pytest.raises(TypeError):matrix.rref(iszerofunc=lambda value:42)


@pytest.mark.parametrize('shape',[(0,0),(0,3),(3,0)])
def test_empty_matrix_shapes_and_bases_remain_valid(shape):
    matrix=sx.zeros(*shape);foreign=sp.zeros(*shape)
    assert converted(matrix.rref(iszerofunc=lambda value:value.is_zero))==foreign.rref()
    assert converted(matrix.nullspace(iszerofunc=lambda value:value.is_zero))==foreign.nullspace()


def test_matrix_callbacks_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
matrix=Matrix([[1e-12,1],[0,1]])
zero=lambda value: abs(complex(value))<1e-9
reduced,pivots=matrix.rref(iszerofunc=zero)
assert reduced==Matrix([[0,1],[0,0]]) and pivots==(1,)
assert matrix.rank(iszerofunc=zero)==1
assert matrix.nullspace(iszerofunc=zero)==[Matrix([1,0])]
x=Symbol('x')
assert Matrix([[x,1],[0,1]]).rank(simplify=lambda value:value.subs(x,0))==1
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=15)
