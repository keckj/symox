"""Algebraic-number wildcards inspect the defining field and alias."""
import pytest
import sympy as sp
import symox as sx


def converted(mapping):
    return None if mapping is None else {key.to_sympy():value.to_sympy() for key,value in mapping.items()}


@pytest.mark.parametrize('old',[False,True])
@pytest.mark.parametrize('excluded',['[]','[x]','[1]','[2]','[sqrt(2)]','[y]'])
def test_wildcard_exclusions_include_coordinates_root_and_alias(old,excluded):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    number=sx.AlgebraicNumber(sx.sqrt(2),alias=x)
    foreign=sp.AlgebraicNumber(sp.sqrt(2),alias=X)
    exclusions=eval(excluded,vars(sx)|{'x':x,'y':y})
    foreign_exclusions=eval(excluded,vars(sp)|{'x':X,'y':Y})
    wild=sx.Wild('a',exclude=exclusions);reference=sp.Wild('a',exclude=foreign_exclusions)
    actual=number.match(wild,old=old)
    expected=foreign.match(reference,old=old)
    assert converted(actual)==expected
    if actual is not None:assert actual[wild] is number


@pytest.mark.parametrize('old',[False,True])
def test_alias_wildcards_and_existing_bindings_preserve_the_number(old):
    x=sx.Symbol('x');a=sx.Wild('a')
    number=sx.AlgebraicNumber(sx.sqrt(2),alias=x)
    pattern=sx.AlgebraicNumber(sx.sqrt(2),alias=a)
    initial={a:x}
    assert pattern.matches(number,initial,old=old)=={a:x}
    assert initial=={a:x}
    assert number.match(pattern,old=old)=={a:x}
    assert number.match(number,old=old)=={}
    assert number.match(sx.sqrt(2),old=old) is None
    assert sx.AlgebraicNumber(sx.sqrt(3),alias=x).match(pattern,old=old) is None


@pytest.mark.parametrize('old',[False,True])
def test_literal_algebraic_number_exclusions_do_not_use_embedding_identity(old):
    number=sx.AlgebraicNumber(sx.sqrt(2))
    wild=sx.Wild('a',exclude=[number])
    assert number.match(wild,old=old) is None
