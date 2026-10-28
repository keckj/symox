"""Global isinstance operands in the stock conversion are observable protocols."""
import importlib
import pytest
import symox as s

sparse=importlib.import_module('symox.polys.matrices.sdm')
methods=importlib.import_module('symox.polys.domains._methods')


@pytest.mark.parametrize('name',['str','bytes','bytearray','list','tuple'])
def test_conversion_global_type_callbacks_match_the_original_path(monkeypatch,name):
    k=s.QQ;a={0:{0:k(2),1:k(3)},1:{0:k(5),1:k(7)}}
    k.sum([k.one])
    events=[]
    class Meta(type):
        def __instancecheck__(cls,value):events.append(name);return False
    class Custom(metaclass=Meta):pass
    monkeypatch.setattr(methods,name,Custom,raising=False)
    actual=sparse.sdm_berk(a,2,k);trace=events[:];events.clear()
    with monkeypatch.context() as patch:
        patch.setattr(sparse,'_sdm_stock_protocol',lambda *args:False)
        expected=sparse.sdm_berk(a,2,k)
    assert actual==expected and trace==events and trace
