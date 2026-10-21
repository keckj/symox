"""Numerical lambda bodies retain exact signatures and lexical binding."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('digits',[5,15,40,80])
@pytest.mark.parametrize('source',[
    'Lambda(x,x/3)','Lambda((x,y),x/3+y/7)',
    'Lambda(x,x+pi)','Lambda(x,x**Rational(1,3))',
    'Lambda(x,sin(x/2))',
])
def test_numeric_lambda_bodies_match_reference(source,digits):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    native=eval(source,vars(sx)|{'x':x,'y':y})
    foreign=eval(source,vars(sp)|{'x':X,'y':Y})
    result=native.evalf(digits)
    expected=foreign.evalf(digits)
    if source=='Lambda(x,x+pi)' and digits==5:
        # SymPy's partial Add evalf truncates pi here, one ulp below its own
        # correctly rounded scalar pi. Native evaluation rounds constants
        # consistently inside and outside lambda bodies.
        expected=sp.Lambda(X,X+sp.pi.evalf(digits))
    assert result.to_sympy()==expected
    assert result.variables==native.variables
    assert result.free_symbols==native.free_symbols


def test_numeric_identity_and_outer_free_parameters_remain_scoped():
    x,y=sx.symbols('x y')
    identity=sx.Lambda(x,x)
    assert identity.evalf(50)==identity
    result=sx.Lambda(x,x/3+y).evalf(50)
    assert result.free_symbols=={y}
    assert result(0)==y
    assert result(1).subs(y,0)==sx.Float(sx.Rational(1,3),50)


def test_ode_solver_rejects_excessive_orders_before_allocating_coefficients():
    x=sx.Symbol('x');f=sx.Function('f')
    expression=sx.Derivative(f(x),x,2**32-1,evaluate=False)
    with pytest.raises(NotImplementedError):sx.dsolve(expression,f(x))
    with pytest.raises(NotImplementedError):sx.dsolve(expression)
    # f = 0 solves f^(n) = 0; the reference runs out of memory here.
    assert sx.checkodesol(expression,sx.Eq(f(x),0),f(x))==(True,0)
    assert sx.checkodesol(expression,sx.Eq(f(x),0))==(True,0)


def test_ode_checks_with_empty_dependent_functions_do_not_panic():
    x=sx.Symbol('x');f=sx.Function('f')
    # The reference rejects a function without arguments.
    with pytest.raises(ValueError, match="function of one variable"):
        sx.checkodesol(f(x),sx.Eq(f(x),0),f())
