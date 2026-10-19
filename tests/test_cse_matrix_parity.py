"""CSE returns matrices with their shape, sparsity and mutability intact."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


def expand_replacements(module,replacements,expression):
    for symbol,value in reversed(replacements):expression=expression.xreplace({symbol:value})
    return expression


@pytest.mark.parametrize('constructor',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('shape',[(0,3),(3,0),(2,2)])
@pytest.mark.parametrize('single',[False,True])
def test_matrix_shape_type_and_reconstruction(constructor,shape,single):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    rows,cols=shape
    data=[x+y,(x+y)**2,sx.sin(x+y),0] if rows*cols else []
    reference=[X+Y,(X+Y)**2,sp.sin(X+Y),0] if rows*cols else []
    matrix=getattr(sx,constructor)(rows,cols,data)
    expected=getattr(sp,constructor)(rows,cols,reference)
    result=sx.cse(matrix if single else [matrix]);oracle=sp.cse(expected if single else [expected])
    replacements,reduced=result
    assert len(reduced)==len(oracle[1])==1
    actual=reduced[0]
    assert type(actual) is type(matrix)
    assert actual._immutable==matrix._immutable
    assert actual.shape==shape
    assert expand_replacements(sx,replacements,actual).to_sympy()==expected


def test_sparse_cse_cost_depends_on_stored_entries():
    x,y=sx.symbols('x y')
    matrix=sx.SparseMatrix(10**9,10**9,{(3,4):(x+y)**2,(5,6):sx.sin(x+y)})
    replacements,reduced=sx.cse(matrix)
    assert len(replacements)==1
    assert reduced[0].shape==matrix.shape and len(reduced[0].todok())==2
    assert expand_replacements(sx,replacements,reduced[0])==matrix


def test_joint_matrices_custom_symbols_postprocessing_and_single_result():
    x,y=sx.symbols('x y');t=sx.Symbol('temporary')
    dense=sx.ImmutableMatrix([x+y,(x+y)**2]);sparse=sx.SparseMatrix([[sx.sin(x+y),0]])
    replacements,reduced=sx.cse([x+y,dense,sparse],symbols=iter([t]))
    assert replacements==[(t,x+y)]
    assert reduced[0]==t and reduced[1][0]==t and reduced[2][0,0]==sx.sin(t)
    seen=[]
    result=sx.cse(dense,postprocess=lambda reps,red:(seen.append(red),red)[1])
    assert isinstance(result[0],sx.ImmutableMatrix) and seen
    assert isinstance(sx.cse(dense,list=False)[1],sx.ImmutableMatrix)


def test_matrix_cse_does_not_import_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
x,y=symbols('x y')
reps,reduced=cse([ImmutableMatrix([x+y,(x+y)**2]),SparseMatrix([[sin(x+y),0]])])
assert len(reps)==1 and reduced[0].shape==(2,1)
assert reduced[1].shape==(1,2)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
