"""Named constructor arguments must preserve the meaning of each position."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('name',['Eq','Ne','Lt','Le','Gt','Ge'])
@pytest.mark.parametrize('evaluate',[False,True])
def test_relational_keyword_operands_and_options(name,evaluate):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=getattr(sx,name)(lhs=x,rhs=2,evaluate=evaluate)
    expected=getattr(sp,name)(lhs=X,rhs=2,evaluate=evaluate)
    assert actual.to_sympy()==expected
    for value in (x,1):
        kwargs=dict(lhs=value,rhs=2,evaluate=evaluate,extra_option=True)
        ref=dict(kwargs,lhs=X if value==x else value)
        try:expected=getattr(sp,name)(**ref)
        except TypeError:
            with pytest.raises(TypeError):getattr(sx,name)(**kwargs)
        else:assert getattr(sx,name)(**kwargs).to_sympy()==expected


def test_power_uses_reference_keyword_names():
    x=sx.Symbol('x')
    assert sx.Pow(b=x,e=2)==x**2
    assert sx.Pow(b=x,e=1,evaluate=False).args==(x,sx.S.One)


@pytest.mark.parametrize('evaluate',[False,True])
def test_ones_dimensions_are_not_swapped_by_keyword_calls(evaluate):
    matrix=sx.OneMatrix(m=2,n=3,evaluate=evaluate)
    expected=sp.OneMatrix(m=2,n=3,evaluate=evaluate)
    assert tuple(map(int,matrix.shape))==expected.shape
    assert matrix.as_explicit().to_sympy()==expected.as_explicit()


@pytest.mark.parametrize('name',['Trace','Determinant','Inverse'])
def test_matrix_constructor_mat_keyword(name):
    matrix=sx.MatrixSymbol('A',2,2);expected=sp.MatrixSymbol('A',2,2)
    assert getattr(sx,name)(mat=matrix).to_sympy()==getattr(sp,name)(mat=expected)


def test_real_keyword_names_in_matrix_geometry_and_held_constructors():
    x=sx.Symbol('x');X=sp.Symbol('x')
    assert sx.det(matexpr=sx.eye(3)).to_sympy()==sp.det(matexpr=sp.eye(3))
    actual=sx.RegularPolygon(c=sx.Point(0,0),r=2,n=5)
    expected=sp.RegularPolygon(c=sp.Point(0,0),r=2,n=5)
    assert actual.to_sympy()==expected
    assert sx.UnevaluatedExpr(arg=x).to_sympy()==sp.UnevaluatedExpr(arg=X)
    for module in (sx,sp):
        with pytest.raises(module.SympifyError):module.UnevaluatedExpr(arg='x')
        with pytest.raises(TypeError):module.UnevaluatedExpr(arg=module.Symbol('x'),evaluate=False)
