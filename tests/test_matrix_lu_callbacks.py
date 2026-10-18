"""LU factorization and solves use native callback-driven pivot decisions."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def converted(value):
    if isinstance(value,(list,tuple)):return type(value)(converted(item) for item in value)
    return value.to_sympy() if hasattr(value,'to_sympy') else value


def zero_test(x,mode):
    def predicate(value):
        if mode=='threshold' and value.is_number:return abs(complex(value))<1e-9
        if mode=='symbol':return value.subs(x,0).is_zero
        return value.is_zero
    return predicate


@pytest.mark.parametrize('constructor',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('operation',['LUdecomposition','LUdecomposition_Simple'])
@pytest.mark.parametrize('source,mode',[
    ('[[x,1],[1,0]]','exact'),
    ('[[x,1],[0,1]]','symbol'),
    ('[[1e-12,1],[0,1]]','threshold'),
    ('[[1,2],[3,4],[5,6]]','exact'),
    ('[[1,2,3],[4,5,6]]','exact'),
])
def test_lu_zero_callbacks_match_reference(constructor,operation,source,mode):
    x=sx.Symbol('x');X=sp.Symbol('x')
    matrix=getattr(sx,constructor)(eval(source,{'x':x}))
    foreign=getattr(sp,constructor)(eval(source,{'x':X}))
    original=matrix.copy()
    actual=getattr(matrix,operation)(iszerofunc=zero_test(x,mode))
    expected=getattr(foreign,operation)(iszerofunc=zero_test(X,mode))
    assert converted(actual)==expected
    assert matrix==original


@pytest.mark.parametrize('constructor',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
def test_lu_solves_preserve_pivot_callbacks_and_rhs(constructor):
    x=sx.Symbol('x');X=sp.Symbol('x')
    matrix=getattr(sx,constructor)([[x,1],[1,0]])
    rhs=getattr(sx,constructor)([1,2]);foreign=getattr(sp,constructor)([[X,1],[1,0]])
    actual=matrix.LUsolve(rhs,iszerofunc=lambda value:value.is_zero)
    expected=foreign.LUsolve(getattr(sp,constructor)([1,2]),iszerofunc=lambda value:value.is_zero)
    assert converted(actual)==expected
    assert matrix*actual==rhs and rhs==getattr(sx,constructor)([1,2])


@pytest.mark.parametrize('operation',['LUdecomposition','LUdecomposition_Simple'])
def test_rankcheck_uses_the_custom_zero_test(operation):
    matrix=sx.Matrix([[1e-12,1],[0,1]])
    with pytest.raises(ValueError):getattr(matrix,operation)(iszerofunc=lambda value:abs(complex(value))<1e-9,rankcheck=True)
    with pytest.raises(sx.NonInvertibleMatrixError):matrix.LUsolve(sx.Matrix([1,1]),iszerofunc=lambda value:abs(complex(value))<1e-9)


@pytest.mark.parametrize('operation',['LUdecomposition','LUdecomposition_Simple'])
def test_simplification_discards_a_pivot_that_was_proved_zero(operation):
    x=sx.Symbol('x');identity=sx.sin(x)**2+sx.cos(x)**2-1
    matrix=sx.Matrix([[identity,1],[0,1]])
    calls=[]
    def simplify(value):
        calls.append(value)
        return sx.trigsimp(value)
    result=getattr(matrix,operation)(simpfunc=simplify)
    assert identity in calls
    if operation=='LUdecomposition':
        lower,upper,swaps=result
        residual=lower*upper-matrix.permute_rows(swaps)
        assert residual.applyfunc(sx.trigsimp)==sx.zeros(2)
        assert upper[0,0].is_zero is True
    else:
        combined,swaps=result
        assert combined[0,0].is_zero is True and swaps==[]


@pytest.mark.parametrize('operation',['LUdecomposition','LUdecomposition_Simple','LUsolve'])
def test_lu_callback_exceptions_preserve_their_identity(operation):
    class Marker(Exception):pass
    failure=Marker('stop LU')
    calls=[]
    def zero(value):
        calls.append(value)
        raise failure
    matrix=sx.Matrix([[sx.Symbol('x'),1],[1,0]])
    args=(sx.Matrix([1,2]),) if operation=='LUsolve' else ()
    with pytest.raises(Marker) as caught:getattr(matrix,operation)(*args,iszerofunc=zero)
    assert caught.value is failure and len(calls)==1


@pytest.mark.parametrize('operation',['LUdecomposition','LUdecomposition_Simple'])
def test_lu_simplifier_exceptions_and_invalid_callbacks_raise(operation):
    class Marker(Exception):pass
    def simplify(value):raise Marker('stop simplify')
    matrix=sx.Matrix([[sx.Symbol('x'),1],[0,1]])
    with pytest.raises(Marker):getattr(matrix,operation)(simpfunc=simplify)
    with pytest.raises(TypeError):getattr(matrix,operation)(simpfunc=42)
    with pytest.raises(TypeError):getattr(matrix,operation)(iszerofunc=sx.true)


def test_native_lambda_callbacks_are_supported():
    t=sx.Symbol('t');matrix=sx.Matrix([[2,1],[1,1]])
    zero=sx.Lambda(t,sx.Eq(t,0))
    assert matrix.LUsolve(sx.Matrix([1,0]),iszerofunc=zero)==sx.Matrix([1,-1])
    assert matrix.LUdecomposition(simpfunc=sx.Lambda(t,t))==matrix.LUdecomposition()


def test_zero_symbols_and_sums_are_used_by_default_matrix_algorithms():
    z=sx.Symbol('z',zero=True);w=sx.Symbol('w',zero=True)
    assert sx.Matrix([[z]]).rank()==0
    for value in [z,z**2,z+w]:
        matrix=sx.Matrix([[value,0],[0,1]])
        assert matrix.rank()==1
        assert matrix.rref()[1]==(1,)
        assert matrix.nullspace()==[sx.Matrix([1,0])]
    matrix=sx.Matrix([[z+1,0],[0,1]])
    assert matrix.rank()==2


def test_native_boolean_simplification_flags_remain_flags():
    x=sx.Symbol('x');matrix=sx.Matrix([[x,1],[1,0]])
    assert matrix.rref(simplify=sx.true)==matrix.rref(simplify=True)
    assert matrix.rref(simplify=sx.false)==matrix.rref(simplify=False)


def test_lu_callbacks_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x');matrix=Matrix([[x,1],[1,0]])
zero=lambda value:value.is_zero
solution=matrix.LUsolve(Matrix([1,2]),iszerofunc=zero)
assert matrix*solution==Matrix([1,2])
lower,upper,swaps=matrix.LUdecomposition(iszerofunc=zero)
assert lower*upper==matrix.permute_rows(swaps)
z=Symbol('z',zero=True)
assert Matrix([[z]]).rank()==0
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=15)
