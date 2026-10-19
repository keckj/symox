"""Mathematical Boolean equality differs from integer truth coercion."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('boolean',[False,True])
@pytest.mark.parametrize('number',[-1,0,1,2])
@pytest.mark.parametrize('name',['Eq','Ne'])
def test_boolean_atoms_and_integer_operands_have_distinct_domains(boolean,number,name):
    for left,right in [(boolean,number),(number,boolean)]:
        assert getattr(sx,name)(left,right).to_sympy()==getattr(sp,name)(left,right)


def test_composite_booleans_and_symbols_keep_reference_domain_rules():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    assert sx.Eq(sx.And(x,y),1).to_sympy()==sp.Eq(sp.And(X,Y),1)
    assert sx.Eq(x,True).to_sympy()==sp.Eq(X,True)
    assert sx.Eq(True,False)==sx.S.false
    assert sx.Ne(True,False)==sx.S.true


@pytest.mark.parametrize('name',['Eq','Ne','Lt','Le','Gt','Ge'])
def test_held_relations_preserve_evaluation_control_during_interoperability(name):
    actual=getattr(sx,name)(lhs=1,rhs=2,evaluate=False)
    expected=getattr(sp,name)(lhs=1,rhs=2,evaluate=False)
    assert actual.to_sympy()==expected
