"""Determinants and inverses pass scalar predicates to native kernels."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


KINDS=['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix']
METHODS=[None,'GE','LU','ADJ','CH','LDL','QR','DM','DMNC']


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('method',METHODS)
@pytest.mark.parametrize('rows',[
    [[2,1],[1,3]],
    [[2,0,0],[0,3,1],[0,1,2]],
])
def test_inverse_zero_callbacks_match_exact_reference(kind,method,rows):
    matrix=getattr(sx,kind)(rows);original=matrix.copy()
    actual=matrix.inv(method=method,iszerofunc=lambda value:value.is_zero)
    expected=getattr(sp,kind)(rows).inv(method=method,iszerofunc=lambda value:value.is_zero)
    assert actual.to_sympy()==expected
    assert matrix*actual==sx.eye(matrix.rows) and matrix==original
    assert type(actual) is type(matrix)
    assert isinstance(actual,sx.ImmutableMatrix)==isinstance(matrix,sx.ImmutableMatrix)


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('method',['GE','LU','ADJ','CH','LDL','QR'])
def test_direct_inverse_methods_accept_the_zero_callback(kind,method):
    x=sx.Symbol('x');X=sp.Symbol('x')
    rows=[[x,1],[1,0]] if method in ['GE','LU','ADJ'] else [[2,1],[1,3]]
    foreign_rows=[[X,1],[1,0]] if method in ['GE','LU','ADJ'] else [[2,1],[1,3]]
    matrix=getattr(sx,kind)(rows)
    actual=getattr(matrix,'inverse_'+method)(iszerofunc=lambda value:value.is_zero)
    expected=getattr(getattr(sp,kind)(foreign_rows),'inverse_'+method)(iszerofunc=lambda value:value.is_zero)
    assert actual.applyfunc(sx.simplify).to_sympy()==expected.applyfunc(sp.simplify)


@pytest.mark.parametrize('method',[None,'GE','LU'])
@pytest.mark.parametrize('source',['[[1e-12,1],[0,1]]','[[x,1],[0,1]]'])
def test_custom_predicate_rejects_a_zero_pivot(method,source):
    x=sx.Symbol('x');matrix=sx.Matrix(eval(source,{'x':x}))
    def zero(value):
        value=value.subs(x,0)
        return abs(complex(value))<1e-9 if value.is_number else value.is_zero
    with pytest.raises(sx.NonInvertibleMatrixError):matrix.inv(method=method,iszerofunc=zero)


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('method',[None,'GE','LU','ADJ','CH','LDL','QR'])
def test_singular_and_empty_matrices(kind,method):
    zero=lambda value:value.is_zero
    with pytest.raises(sx.NonInvertibleMatrixError):getattr(sx,kind)([[1,2],[2,4]]).inv(method=method,iszerofunc=zero)
    empty=getattr(sx,kind)(0,0,[])
    assert empty.inv(method=method,iszerofunc=zero)==empty


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('method',[None,'GE','LU','ADJ','CH','LDL','QR'])
def test_block_diagonal_inverse_retains_algorithm_and_predicate(kind,method):
    rows=[[2,1,0,0],[1,3,0,0],[0,0,4,1],[0,0,1,2]]
    matrix=getattr(sx,kind)(rows)
    result=matrix.inv(method=method,iszerofunc=lambda value:value.is_zero,try_block_diag=True)
    expected=getattr(sp,kind)(rows).inv(method=method,iszerofunc=lambda value:value.is_zero,try_block_diag=True)
    assert result.to_sympy()==expected and matrix*result==sx.eye(4)


@pytest.mark.parametrize('operation',['inv','inverse_GE','inverse_LU'])
def test_inverse_predicate_exceptions_and_invalid_values(operation):
    class Marker(Exception):pass
    failure=Marker('stop inversion');calls=[]
    def zero(value):
        calls.append(value)
        raise failure
    matrix=sx.Matrix([[sx.Symbol('x'),1],[1,0]])
    with pytest.raises(Marker) as caught:getattr(matrix,operation)(iszerofunc=zero)
    assert caught.value is failure and len(calls)==1
    with pytest.raises(TypeError):getattr(matrix,operation)(iszerofunc=2)
    with pytest.raises(TypeError):getattr(matrix,operation)(iszerofunc=lambda value:1)


def test_native_and_numpy_boolean_callbacks_and_scalar_lambdas():
    import numpy as np
    t=sx.Symbol('t');matrix=sx.Matrix([[2,1],[1,1]])
    assert matrix.inv(iszerofunc=sx.Lambda(t,sx.Eq(t,0)))==matrix.inv()
    assert matrix.inverse_LU(iszerofunc=lambda value:np.bool_(value.is_zero))==matrix.inv()
    with pytest.raises(sx.NonInvertibleMatrixError):sx.Matrix([[sx.Symbol('z',zero=True)]]).inverse_ADJ(iszerofunc=lambda value:value.is_zero)


def test_inverse_callbacks_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
matrix=Matrix([[2,1],[1,3]])
zero=lambda value:value.is_zero
for method in [None,'GE','LU','ADJ','CH','LDL','QR']:
    assert matrix*matrix.inv(method=method,iszerofunc=zero)==eye(2)
assert diag(matrix,matrix).inv(iszerofunc=zero,try_block_diag=True)==diag(matrix.inv(),matrix.inv())
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
