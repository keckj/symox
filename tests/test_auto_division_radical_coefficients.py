"""Automatic division uses exact native algebraic coefficient fields."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('formula',['sqrt(2)','sqrt(5)','sqrt(2)+sqrt(3)','I','I*sqrt(2)'])
@pytest.mark.parametrize('method',['div','quo','rem'])
def test_default_division_with_inferred_algebraic_coefficients_is_exact(formula,method):
    def run(lib):
        x=lib.Symbol('x')
        coefficient=lib.sympify(formula)
        dividend=lib.Poly(x*x+2*x+3,x)
        divisor=lib.Poly((coefficient+1)*x+coefficient,x)
        quotient,remainder=dividend.div(divisor)
        assert remainder.is_zero or remainder.degree()<divisor.degree()
        residual=dividend.as_expr()-quotient.as_expr()*divisor.as_expr()-remainder.as_expr()
        oracle=residual.to_sympy() if lib is sx else residual
        assert sp.simplify(oracle)==0
        result=getattr(dividend,method)(divisor)
        values=result if isinstance(result,tuple) else (result,)
        assert not any(value.as_expr().has(lib.Float) for value in values)
        return values
    actual,expected=run(sx),run(sp)
    assert len(actual)==len(expected)
    assert all(sp.simplify(a.as_expr().to_sympy()-b.as_expr())==0 for a,b in zip(actual,expected))


def test_unavailable_transcendental_coefficient_field_is_not_faked():
    x=sx.Symbol('x')
    with pytest.raises(ValueError,match='non-polynomial quotient'):
        sx.Poly(x+1,x).div(sx.Poly((sx.pi+1)*x+1,x))
