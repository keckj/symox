"""EXRAW coefficient relations must preserve full eigenvalue multiplicities."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('kind',['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('shape',['zero','rank_one','jordan','scalar'])
def test_complete_eigenvalue_groups_and_eigenspace_dimensions(kind,shape):
    x=sx.Symbol('native_eigen_x')
    zero=sx.cos(2*x)-sx.cos(x)**2+sx.sin(x)**2
    sources={
        'zero':[[zero,zero],[zero,zero]],
        'rank_one':[[zero+1,1],[1,1]],
        'jordan':[[zero,1],[0,zero]],
        'scalar':[[zero+2,0],[0,zero+2]],
    }
    plain={
        'zero':[[0,0],[0,0]],'rank_one':[[1,1],[1,1]],
        'jordan':[[0,1],[0,0]],'scalar':[[2,0],[0,2]],
    }
    matrix=getattr(sx,kind)(sources[shape])
    expected=getattr(sp,kind)(plain[shape]).eigenvects()
    actual=matrix.eigenvects()
    assert [(sx.to_sympy(value),mult,len(basis)) for value,mult,basis in actual] == \
           [(value,mult,len(basis)) for value,mult,basis in expected]
    assert sum(mult for _,mult,_ in actual)==matrix.rows
    for value,mult,basis in actual:
        columns=[]
        for vector in basis:
            residual=matrix*vector-value*vector
            assert all(sp.trigsimp(sx.to_sympy(entry))==0 for entry in residual)
            columns.append(sp.Matrix([sx.to_sympy(entry) for entry in vector]))
        assert sp.Matrix.hstack(*columns).rank()==len(basis)
    # No Float approximation, sampled parameters or finite prefix of roots.
    assert not any(value.has(sx.Float) for value,_,_ in actual)


@pytest.mark.parametrize('simplify',[True,False])
def test_default_zero_detection_and_explicit_simplify_options(simplify):
    x=sx.Symbol('native_eigen_flags')
    zero=sx.cos(2*x)-sx.cos(x)**2+sx.sin(x)**2
    matrix=sx.Matrix([[zero,zero],[zero,zero]])
    result=matrix.eigenvects(simplify=simplify)
    assert [(v,m,len(b)) for v,m,b in result]==[(0,2,2)]


def test_symbolic_shift_retains_the_full_repeated_eigenvalue():
    x,a=sx.symbols('native_eigen_phase native_eigen_shift')
    zero=sx.cos(2*x)-sx.cos(x)**2+sx.sin(x)**2
    result=sx.Matrix([[a+zero,1],[0,a+zero]]).eigenvects()
    assert [(v,m,len(b)) for v,m,b in result]==[(a,2,1)]
