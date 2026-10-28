"""Native eligibility compares generator class identities without callbacks."""
import importlib
import pytest
import symox as sx

rt=importlib.import_module('symox.integrals.rationaltools')
import symox.holonomic


@pytest.mark.parametrize('mode',['raise','lie'])
def test_host_metaclass_equality_is_neither_called_nor_trusted_during_selection(mode,monkeypatch):
    class Meta(type(sx.Symbol)):
        __hash__=type.__hash__
    class Host(sx.Symbol,metaclass=Meta):pass
    z=Host('z');f=sx.Poly(z+1,z,domain=sx.QQ);g=sx.Poly(1,z,domain=sx.QQ)
    calls=[]
    class ComparisonError(RuntimeError):pass
    def equality(cls,other):
        calls.append((cls,other))
        if mode=='raise':raise ComparisonError('generator class equality')
        return True
    monkeypatch.setattr(Meta,'__eq__',equality,raising=False)
    assert rt._stock_atan_protocol(f,g) is False
    assert calls==[]
    # Selection leaves the user's protocol intact for an actual comparison.
    if mode=='raise':
        with pytest.raises(ComparisonError,match='generator class equality'):
            Host==rt._ATAN_SCALARS[0]
    else:
        assert Host==rt._ATAN_SCALARS[0]
    assert calls
