"""Upgrade real cached calculus placeholders while preserving user functions."""
import base64
import json
import pickle
from pathlib import Path

import pytest
import symox as s
import sympy as sp


def test_previous_native_pickle_placeholders_upgrade_to_binding_operators():
    # Captured using the compiled fef8d54 wheel, before typed binders existed.
    fixtures=json.loads((Path(__file__).parent/'fixtures/calculus_pickle_v0.json').read_text())
    x=s.Symbol('x',positive=True)
    y=s.Symbol('y')
    f=s.Function('f')
    expected={
        'indefinite':s.Integral(f(x),x),
        'definite':s.Integral(f(x),(x,0,y)),
        'sum':s.Sum(f(x),(x,0,s.oo)),
        'user_integral':s.Function('Integral')(x),
        'user_sum':s.Function('Sum')(x),
    }
    for name,payload in fixtures.items():
        restored=pickle.loads(base64.b64decode(payload))
        assert restored==expected[name]
        assert restored.free_symbols==expected[name].free_symbols
        assert restored.args==expected[name].args
        assert pickle.loads(pickle.dumps(restored))==restored
        assert restored.is_UndefinedFunction==name.startswith('user_')


@pytest.mark.parametrize('kind',['Integral','Sum'])
def test_pickle_reader_rejects_invalid_binding_limit_records(kind):
    x=s.Symbol('x')
    value=getattr(s,kind)(s.Function('f')(x),(x,0,1))
    restore,(version,records,root)=value.__reduce__()
    records=list(records)
    symbol=next(i for i,(tag,_,fields) in enumerate(records) if tag in ('s.symbol','s.symbol.facts') and fields[0]=='x')
    tag,canonical,fields=records[root]
    assert tag=='s.function' and fields[0]==kind
    records[root]=(tag,canonical,(fields[0],[fields[1][0],symbol]))
    with pytest.raises(ValueError,match='limit|Tuple'):
        restore(version,records,root)


@pytest.mark.parametrize('construct',[
    lambda m,x,y:m.Integral(x,x),
    lambda m,x,y:m.Integral(m.Function('f')(x),(x,0,y)),
    lambda m,x,y:m.Integral(x*y,(x,-1,y/2),(y,0,1)),
    lambda m,x,y:m.Sum(x,(x,1,y)),
    lambda m,x,y:m.Subs(m.Function('f')(x),x,y),
])
def test_binding_tuples_do_not_add_arithmetic_operations(construct):
    x,y=s.symbols('x y')
    X,Y=sp.symbols('x y')
    value=construct(s,x,y)
    expected=construct(sp,X,Y)
    assert value.count_ops()==expected.count_ops()
    assert s.count_ops(value)==sp.count_ops(expected)


@pytest.mark.parametrize('source',[
    'E','-1.2','x/2','2*x/3','x/(2*y)','2*x**2/(3*y**2)',
    'Subs(f(x),x,y)','Limit(sin(x)/x,x,0)','Order(x**2)',
    '-x/(2*y)','(x-y)/2','3*x/7-y/2','1/(2*x*y)',
])
def test_numeric_coefficients_and_held_wrappers_have_sympy_operation_counts(source):
    actual=s.sympify(source)
    expected=sp.sympify(source)
    assert actual.count_ops()==expected.count_ops()


def test_held_exponential_power_counts_as_one_exponential():
    x=s.Symbol('x')
    X=sp.Symbol('x')
    assert s.Pow(s.E,x,evaluate=False).count_ops()==sp.Pow(sp.E,X,evaluate=False).count_ops()==1
