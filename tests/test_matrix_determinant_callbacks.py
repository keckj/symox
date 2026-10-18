"""Pivot predicates affect Bareiss and LU while algebraic formulas stay exact."""
import pytest
import sympy as sp
import symox as sx

KINDS=['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix']
ROWS=[[1,1,0,1],[2,1,2,1],[3,0,1,2],[4,2,1,1]]


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('method',['bareiss','lu','berkowitz','domain-ge','bird','laplace'])
def test_determinant_predicates_match_reference(kind,method):
    x=sx.Symbol('x');X=sp.Symbol('x')
    rows=[list(row) for row in ROWS];foreign=[list(row) for row in ROWS]
    rows[0][0]=x;foreign[0][0]=X
    matrix=getattr(sx,kind)(rows);original=matrix.copy()
    actual=matrix.det(method=method,iszerofunc=lambda value:value.is_zero)
    expected=getattr(sp,kind)(foreign).det(method=method,iszerofunc=lambda value:value.is_zero)
    assert sx.expand(actual).to_sympy()==sp.expand(expected) and matrix==original


@pytest.mark.parametrize('kind',KINDS)
@pytest.mark.parametrize('method',['bareiss','lu'])
def test_numeric_predicate_is_not_bypassed_by_domain_fast_paths(kind,method):
    rows=[[sx.Rational(row[0],10**12),*row[1:]] for row in ROWS]
    matrix=getattr(sx,kind)(rows)
    calls=[]
    def zero(value):
        calls.append(value)
        return abs(complex(value))<1e-9
    assert matrix.det(method=method)!=0
    assert matrix.det(method=method,iszerofunc=zero)==0
    assert calls
    reference=getattr(sp,kind)([[sp.Rational(row[0],10**12),*row[1:]] for row in ROWS])
    assert reference.det(method=method,iszerofunc=lambda value:abs(complex(value))<1e-9)==0


@pytest.mark.parametrize('method',['bareiss','lu'])
def test_determinant_callback_errors_are_immediate(method):
    class Marker(Exception):pass
    failure=Marker('stop determinant');calls=[]
    def zero(value):
        calls.append(value)
        raise failure
    matrix=sx.Matrix(ROWS)
    with pytest.raises(Marker) as caught:matrix.det(method=method,iszerofunc=zero)
    assert caught.value is failure and len(calls)==1
    with pytest.raises(TypeError):matrix.det(method=method,iszerofunc=lambda value:2)


@pytest.mark.parametrize('method',['bareiss','lu','berkowitz','domain-ge'])
def test_small_determinants_use_closed_form_without_zero_predicates(method):
    def zero(value):raise AssertionError('closed form should not search pivots')
    assert sx.Matrix([[1,2],[3,4]]).det(method=method,iszerofunc=zero)==-2
    assert sx.Matrix(0,0,[]).det(method=method,iszerofunc=zero)==1
    with pytest.raises(TypeError):sx.eye(2).det(method=method,iszerofunc=42)


def test_division_free_and_domain_determinants_do_not_search_custom_pivots():
    def zero(value):raise AssertionError('algorithm has no pivot search')
    matrix=sx.Matrix(ROWS)
    for method in ['berkowitz','domain-ge','bird','laplace']:
        assert matrix.det(method=method,iszerofunc=zero)==matrix.det()
