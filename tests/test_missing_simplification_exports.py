import pytest
import symox as s
import sympy as sp
from symox._transforms import MellinTransformStripError


@pytest.mark.parametrize('path',['/*','/*/*/Symbol','/*/*/[0]','/*/*/Symbol|Integer','/*/[0:2]','/*/[-1]','/*/*/is_Symbol?','/*/*/[99]'])
def test_expression_paths_match_oracle(path):
    x,y,z,t=s.symbols('x y z t')
    e=t+s.sin(x+1)+s.cos(x+y+s.E)
    expected=sp.epath(path,e.to_sympy())
    actual=s.epath(path,e)
    assert [a.to_sympy() for a in actual]==expected
    result=s.epath(path,e,lambda v:2*v)
    assert result.to_sympy()==sp.epath(path,e.to_sympy(),lambda v:2*v)
    compiled=s.EPath(path)
    assert s.EPath(compiled) is compiled


def test_nested_container_paths_and_arguments():
    x,y,z=s.symbols('x y z')
    e=[((x,1),2),((3,y),z)]
    assert s.epath('/*/[0]/Symbol',e)==[x,y]
    assert s.epath('/*/[0]/Symbol',e,lambda v,k,shift:v**k+shift,args=(2,),kwargs={'shift':1})==[((x*x+1,1),2),((3,y*y+1),z)]
    for path in ['', '/', '/[1:2:0]','/Symbol!','/[x]']:
        with pytest.raises(ValueError):s.EPath(path)
    with pytest.raises(NotImplementedError):s.EPath('Symbol')


def test_fu_registry_uses_native_rules():
    from sympy.simplify.fu import FU
    assert set(s.FU)==set(FU)
    x=s.Symbol('x')
    examples={'TR1':1/s.cos(x),'TR2':s.tan(x),'TR3':s.sin(-x),'TR5':s.sin(x)**2,'TR6':s.cos(x)**2,'TR7':s.cos(x)**2,'TR10':s.sin(x+1),'TR11':s.sin(2*x),'TR0':(x+1)**2}
    for name,e in examples.items():
        a=s.FU[name](e).to_sympy();b=FU[name](e.to_sympy())
        assert sp.simplify(a-b)==0


@pytest.mark.parametrize('shift',[1,2,-1,s.Rational(1,2)])
def test_kronecker_incompatible_constraints(shift):
    i,j=s.symbols('i j',integer=True)
    e=1+s.KroneckerDelta(i,j)*s.KroneckerDelta(i,j+shift)
    assert s.kroneckersimp(e)==1
    assert s.kroneckersimp(s.sin(e))==s.sin(1)
    assert s.kroneckersimp(i+1)==i+1


def test_mellin_strip_exception_protocol():
    err=MellinTransformStripError('no common strip')
    assert isinstance(err,ValueError) and str(err)=='no common strip'
