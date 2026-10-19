"""Symbolic mappings retain native structure during float rationalization."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('rational',[None,True,False])
@pytest.mark.parametrize('conversion',['base10','exact'])
@pytest.mark.parametrize('container',[dict,sx.Dict])
def test_symbolic_mapping_rationalization_matches_reference(rational,conversion,container):
    x=sx.Symbol('x')
    values={sx.Float('.2'):x+sx.Float('.1'),x:sx.Float('1.4142135623730951')}
    reference={key.to_sympy():value.to_sympy() for key,value in values.items()}
    result=sx.nsimplify(container(values),rational=rational,rational_conversion=conversion)
    expected=sp.nsimplify(reference,rational=rational,rational_conversion=conversion)
    assert isinstance(result,sx.Dict)
    assert result.to_sympy()==expected
    assert result.free_symbols=={x}
    assert not result.atoms(sx.Float)


@pytest.mark.parametrize('query',[None,'Symbol','Float','Dict','Tuple'])
@pytest.mark.parametrize('empty',[False,True])
def test_mapping_atom_queries_match_reference(query,empty):
    value=sx.Dict() if empty else sx.Dict({sx.Symbol('x'):sx.Tuple(2,sx.Float('.1'))})
    reference=value.to_sympy()
    actual=value.atoms(*(() if query is None else (getattr(sx,query),)))
    expected=reference.atoms(*(() if query is None else (getattr(sp,query),)))
    assert {atom.to_sympy() for atom in actual}==expected


def test_symbolic_precision_and_rational_options_precede_numeric_recognition():
    x=sx.Symbol('x');coefficient=sx.Float('1.4142135623730951',80)
    for rational in [None,True,False]:
        result=sx.nsimplify(coefficient*x,rational=rational,tolerance=sx.Rational(1,1000),full=True)
        expected=sp.nsimplify((coefficient*x).to_sympy(),rational=rational,tolerance=sp.Rational(1,1000),full=True)
        assert result.to_sympy()==expected


def test_sympify_preserves_symbolic_mapping_identity():
    value=sx.Dict({sx.Symbol('x'):sx.Rational(1,3)})
    assert sx.sympify(value) is value


def test_numeric_mapping_uses_rationalization_or_preserves_reference_error():
    value=sx.Dict({sx.Integer(1):sx.Float('.1')})
    assert sx.nsimplify(value,rational=True).to_sympy()==sp.nsimplify(value.to_sympy(),rational=True)
    for rational in [None,False]:
        with pytest.raises(AttributeError,match='evalf'):sx.nsimplify(value,rational=rational)
        with pytest.raises(AttributeError,match='evalf'):sp.nsimplify(value.to_sympy(),rational=rational)
