"""The func hint must invoke native special-function expansion."""
import pytest
import sympy as sp
import symox as sx


def expression(module, case):
    x,y=module.symbols('x y')
    theta,phi=module.symbols('theta phi', real=True)
    return [module.gamma(x+2), module.sin(module.gamma(x+2)),
            module.binomial(x,3), module.beta(x,y), module.erf2(x,y),
            module.Ynm(0,0,theta,phi), module.Ynm(1,0,theta,phi),
            module.Ynm(1,-1,theta,phi), module.harmonic(x+2)][case]


@pytest.mark.parametrize('case',range(9))
@pytest.mark.parametrize('options',[{}, {'deep':False}, {'mul':False}])
def test_functional_expansion_and_independent_hints(case, options):
    value=expression(sx,case)
    expected=expression(sp,case).expand(func=True,**options)
    actual=value.expand(func=True,**options)
    assert sp.simplify(actual.to_sympy()-expected)==0
    assert sx.expand(value,func=True,**options)==actual


@pytest.mark.parametrize('mutable',[False,True])
def test_matrix_entries_use_functional_hint(mutable):
    x=sx.Symbol('x');X=sp.Symbol('x')
    constructor=sx.Matrix if mutable else sx.ImmutableMatrix
    value=constructor([[sx.gamma(x+2),sx.sin(sx.gamma(x+2))]])
    expected=sp.Matrix([[sp.gamma(X+2),sp.sin(sp.gamma(X+2))]]).expand(func=True)
    assert value.expand(func=True).to_sympy()==expected
