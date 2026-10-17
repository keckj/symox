"""Read previous Lambda caches and retain native signatures in new DAG records."""
import base64
import json
import pickle
from pathlib import Path

import pytest
import sympy as sp
import symox as s


def test_real_previous_lambda_and_facade_pickles_upgrade_without_losing_scope():
    # Captured from the compiled 9836b48 wheel before the signature migration.
    fixtures=json.loads((Path(__file__).parent/'fixtures/lambda_pickle_v0.json').read_text())
    x=s.Symbol('x',positive=True)
    y,z=s.symbols('y z')
    expected={
        'native_mono':s.Lambda(x,x+z),
        'facade_mono':s.Lambda(x,x+z),
        'native_identity':s.Lambda(x,x),
        'facade_identity':s.Lambda(x,x),
        'image':s.ImageSet(s.Lambda(x,x*x+z),s.S.Integers),
        'nested_native':s.Lambda(x,s.Lambda(y,x+y+z)),
        'user_lambda':s.Function('Lambda')(x),
    }
    for name,payload in fixtures.items():
        actual=pickle.loads(base64.b64decode(payload))
        assert actual==expected[name]
        assert actual.args==expected[name].args
        assert actual.free_symbols==expected[name].free_symbols
        assert pickle.loads(pickle.dumps(actual))==expected[name]
        if name.startswith(('native_','facade_')):
            assert actual(2)==expected[name](2)
            assert s.Add(1,actual)==s.Add(1,expected[name])
        if name=='user_lambda':
            assert actual.is_UndefinedFunction


@pytest.mark.parametrize('protocol',range(pickle.HIGHEST_PROTOCOL+1))
def test_signature_structure_dummy_identity_and_shared_body_survive_pickle(protocol):
    i=s.Dummy('i',integer=True)
    j=s.Dummy('i',integer=True)
    x,z=s.symbols('x z')
    nested=s.Lambda((i,(j,x)),i+2*j+x+z)
    zero=s.Lambda((),z)
    values=[nested,zero,s.Lambda(x,x),s.Lambda((i,j),s.Function('f')(i,j)+z)]
    for value in values:
        original=value.__reduce__()[1]
        restored=pickle.loads(pickle.dumps(value,protocol=protocol))
        assert restored==value and restored.args==value.args
        assert restored.signature==value.signature
        assert restored.variables==value.variables
        assert restored.free_symbols==value.free_symbols
        assert restored.__reduce__()[1]==original
        if value.is_IdentityFunction:
            assert restored is s.S.IdentityFunction
    assert pickle.loads(pickle.dumps(nested,protocol=protocol))(2,(3,4))==12+z


def test_pickle_reader_rejects_duplicate_lambda_parameters():
    x,y=s.symbols('x y')
    value=s.Lambda((x,y),x+y)
    restore,(version,records,root)=value.__reduce__()
    records=list(records)
    _,_,fields=records[root]
    assert fields[0]=='Lambda'
    signature=fields[1][0]
    tag,canonical,fields=records[signature]
    assert tag in ('s.undefined', 's.undefined.facts') and fields[0]=='Tuple'
    duplicated=(*fields[:4],[fields[4][0],fields[4][0]])
    records[signature]=(tag,canonical,duplicated)
    with pytest.raises(ValueError,match='Duplicate|signature'):
        restore(version,records,root)


def test_pickle_reader_rejects_identity_payload_arguments():
    x,y=s.symbols('x y')
    value=s.Lambda((x,y),x+y)
    restore,(version,records,root)=value.__reduce__()
    records=list(records)
    tag,canonical,fields=records[root]
    records[root]=(tag,canonical,('IdentityFunction',[fields[1][1]]))
    with pytest.raises(ValueError,match='IdentityFunction|arity'):
        restore(version,records,root)


@pytest.mark.parametrize('operation',[
    lambda m,x,y:m.Lambda(x,x),
    lambda m,x,y:m.Lambda((x,y),x+y),
    lambda m,x,y:m.Lambda((x,(y,)),m.sin(x)+y/2),
    lambda m,x,y:m.Lambda((),x/2),
])
def test_lambda_operation_counts_follow_the_public_arguments(operation):
    x,y=s.symbols('x y')
    X,Y=sp.symbols('x y')
    assert operation(s,x,y).count_ops()==operation(sp,X,Y).count_ops()
