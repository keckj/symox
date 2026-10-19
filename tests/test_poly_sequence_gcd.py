"""GCD and LCM accept symbolic sequence inputs without discarding errors."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('operation',['gcd','lcm'])
@pytest.mark.parametrize('values',[[],[6,9],['x**2-1','x**2-2*x+1'],['2/3','4/9'],['-2/3'],['-2*x']])
@pytest.mark.parametrize('container',['tuple','set','list','python_tuple','python_set'])
def test_symbolic_sequences_match_reference(operation,values,container):
    native=[sx.sympify(value) for value in values]
    reference=[sp.sympify(value) for value in values]
    make_native={'tuple':sx.Tuple,'set':sx.FiniteSet,'list':lambda *v:list(v),
                 'python_tuple':lambda *v:v,'python_set':lambda *v:set(v)}[container]
    make_reference={'tuple':sp.Tuple,'set':sp.FiniteSet,'list':lambda *v:list(v),
                    'python_tuple':lambda *v:v,'python_set':lambda *v:set(v)}[container]
    result=getattr(sx,operation)(make_native(*native))
    expected=getattr(sp,operation)(make_reference(*reference))
    assert result.to_sympy()==expected
    assert getattr(sx,operation)(sp.Tuple(*reference)).to_sympy()==getattr(sp,operation)(sp.Tuple(*reference))


@pytest.mark.parametrize('operation',['gcd','lcm'])
def test_sequence_conversion_errors_propagate_unchanged(operation):
    error=RuntimeError('coefficient conversion sentinel')
    class Coefficient:
        def _sympy_(self):raise error
    with pytest.raises(RuntimeError) as caught:getattr(sx,operation)([Coefficient()])
    assert caught.value is error
