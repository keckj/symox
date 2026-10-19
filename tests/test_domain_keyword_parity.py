"""Domain calls use reference names and retain unification semantics."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('name',['ZZ','QQ','RR','CC','EX','EXRAW','ZZ_I','QQ_I'])
def test_domain_conversion_argument_name(name):
    actual=getattr(sx,name);expected=getattr(sp,name)
    value=sx.Integer(3);reference=sp.Integer(3)
    source_name='expr' if name in ('RR','CC') else 'a'
    target_name='element' if name in ('RR','CC') else 'a'
    element=actual.from_sympy(**{source_name:value})
    reference_element=expected.from_sympy(**{source_name:reference})
    assert actual.to_sympy(**{target_name:element}).to_sympy()==expected.to_sympy(**{target_name:reference_element})


@pytest.mark.parametrize('left,right',[
    ('ZZ','QQ'),('QQ','RR'),('RR','CC'),('EX','EXRAW'),('EXRAW','EX')])
def test_domain_unification_accepts_named_other_and_symbols(left,right):
    x=sx.Symbol('x');X=sp.Symbol('x')
    a=getattr(sx,left).unify(K1=getattr(sx,right),symbols=[x])
    b=getattr(sp,left).unify(K1=getattr(sp,right),symbols=[X])
    assert str(a)==str(b)


def test_unification_rejects_generator_capture():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    for a,b in [(sx.ZZ[x],sx.QQ),(sx.QQ,sx.ZZ[x])]:
        with pytest.raises(sx.UnificationFailed):a.unify(K1=b,symbols=[x])
    for a,b in [(sp.ZZ[X],sp.QQ),(sp.QQ,sp.ZZ[X])]:
        with pytest.raises(sp.UnificationFailed):a.unify(K1=b,symbols=[X])
    assert str(sx.ZZ[x].unify(K1=sx.QQ,symbols=[y]))==str(sp.ZZ[X].unify(K1=sp.QQ,symbols=[Y]))


@pytest.mark.parametrize('other',['ZZ','QQ','RR','CC','EX','EXRAW','ZZ_I','QQ_I'])
def test_raw_expression_domain_priority_in_both_orders(other):
    for left,right in [('EXRAW',other),(other,'EXRAW')]:
        assert str(getattr(sx,left).unify(getattr(sx,right)))==str(getattr(sp,left).unify(getattr(sp,right)))
