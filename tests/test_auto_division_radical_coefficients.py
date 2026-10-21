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


@pytest.mark.parametrize('coefficient',['pi','log(a)','sin(y)'])
@pytest.mark.parametrize('method',['div','quo','rem'])
def test_transcendental_coefficients_divide_in_the_fraction_field(coefficient,method):
    # As the reference: over ZZ[c] a non-polynomial quotient is computed in
    # ZZ(c) and each result is retracted to the ring when possible.
    def run(lib):
        x=lib.Symbol('x')
        c=lib.sympify(coefficient)
        result=getattr(lib.Poly(x+1,x),method)(lib.Poly((c+1)*x+1,x))
        values=result if isinstance(result,tuple) else (result,)
        return [(str(v),str(v.domain)) for v in values]
    assert run(sx)==run(sp)
