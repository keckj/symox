"""Integer bit operations and Boolean connectives have different domains."""
import operator
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('name',['and_','or_','xor','lshift','rshift'])
def test_numeric_and_boolean_operand_domains(name):
    pairs=[(sx.Integer(-2),sp.Integer(-2)),(sx.Integer(0),sp.Integer(0)),
           (sx.Integer(3),sp.Integer(3)),(sx.Integer(10**100),sp.Integer(10**100)),
           (2,2),(True,True),(sx.true,sp.true),(sx.Symbol('x'),sp.Symbol('x')),
           (sx.Float(1),sp.Float(1)),(sx.pi,sp.pi)]
    operation=getattr(operator,name)
    for a,A in pairs:
        for b,B in pairs:
            if name in ('lshift','rshift') and (B==10**100):continue
            try:expected=operation(A,B)
            except (TypeError,ValueError) as error:
                with pytest.raises(type(error)):operation(a,b)
            else:
                actual=operation(a,b)
                assert (actual.to_sympy() if hasattr(actual,'to_sympy') else actual)==expected,(name,a,b)


def test_integer_complement_and_set_boolean_operations():
    for value in [-10**100,-2,-1,0,1,2,10**100]:
        assert int(~sx.Integer(value))==~value
    x=sx.Symbol('x');X=sp.Symbol('x')
    assert (~x).to_sympy()==~X
    for value in [sx.Float(1),sx.pi,sx.I]:
        with pytest.raises(TypeError):operator.invert(value)
    assert sx.FiniteSet(1,2)&sx.FiniteSet(2,3)==sx.FiniteSet(2)
    assert sx.FiniteSet(1,2)|sx.FiniteSet(2,3)==sx.FiniteSet(1,2,3)
