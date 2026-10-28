import importlib
import pytest
import symox as s
m=importlib.import_module('symox.polys.matrices.sdm')

@pytest.mark.parametrize('hook',['__getattribute__','__getattr__'])
def test_domain_attribute_protocol_keeps_the_existing_value_and_trace(monkeypatch,hook):
    k=s.QQ;cls=type(k);a={0:{0:k(2),1:k(3)},1:{0:k(5),1:k(7)}};k.sum([k.one]);events=[]
    old=cls.__getattribute__
    def callback(self,name):
        events.append(name)
        if hook=='__getattr__':raise AttributeError(name)
        return old(self,name)
    monkeypatch.setattr(cls,hook,callback,raising=False)
    assert not m._sdm_stock_protocol(a,k)
    assert events==[]
    actual=m.sdm_berk(a,2,k);trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(m,'_sdm_stock_protocol',lambda *args:False)
        expected=m.sdm_berk(a,2,k)
    assert actual==expected and events==trace
    if hook=='__getattribute__':assert trace
