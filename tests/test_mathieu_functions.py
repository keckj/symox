"""Native Mathieu rules checked by identities and a SymPy oracle."""
import pickle
import pytest
import symox as s
import sympy as sp
from symox.functions.special.mathieu_functions import MathieuBase


@pytest.mark.parametrize('name',['mathieuc','mathieus','mathieucprime','mathieusprime'])
@pytest.mark.parametrize('value',[0,1,2,s.Rational(1,2),-1])
def test_mathieu_zero_q_and_parity(name,value):
    a,z=s.symbols('a z')
    for av in (a,value):
        actual=getattr(s,name)(av,0,z)
        expected=getattr(sp,name)(av.to_sympy() if hasattr(av,'to_sympy') else av,0,z.to_sympy())
        assert sp.simplify(actual.to_sympy()-expected)==0
    actual=getattr(s,name)(a,2,-z)
    assert actual.to_sympy()==getattr(sp,name)(a.to_sympy(),2,-z.to_sympy())


@pytest.mark.parametrize('name',['mathieuc','mathieus'])
def test_mathieu_differential_equation(name):
    a,q,z=s.symbols('a q z')
    f=getattr(s,name)(a,q,z)
    assert isinstance(f,MathieuBase)
    assert s.expand(f.diff(z,2)+(a-2*q*s.cos(2*z))*f)==0
    assert f.diff(z).class_name==name+'prime'
    assert f.fdiff(3)==f.diff(z)
    # diff(z, 2) ends with factor_terms(signsimp(...)), as in SymPy
    assert s.expand(f.diff(z).fdiff(3)-f.diff(z,2))==0
    assert f.diff(a).class_name=='Derivative'
    assert s.conjugate(f)==getattr(s,name)(s.conjugate(a),s.conjugate(q),s.conjugate(z))
    assert pickle.loads(pickle.dumps(f))==f
    nested=s.sin(f)+f*f
    assert nested.to_sympy()==sp.sin(f.to_sympy())+f.to_sympy()**2


def test_mathieu_held_option_is_explicit():
    z=s.Symbol('z')
    with pytest.raises(NotImplementedError):s.mathieuc(1,0,z,evaluate=False)
