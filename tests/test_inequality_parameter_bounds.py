"""Native half-line proofs preserve known coefficient signs and odd roots."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('power',[3,5,7])
@pytest.mark.parametrize('bound',[-32,-8,0,1,27])
@pytest.mark.parametrize('relation',['<','<=','>','>='])
def test_odd_powers_use_real_roots(power,bound,relation):
    x=sx.Symbol('x');X=sp.Symbol('x')
    source=f'(x+1)**{power} {relation} {bound}'
    actual=sx.solve_univariate_inequality(sx.sympify(source),x,relational=False)
    expected=sp.solve_univariate_inequality(sp.sympify(source),X,relational=False)
    assert actual.to_sympy()==expected


@pytest.mark.parametrize('sign',['positive','negative'])
@pytest.mark.parametrize('relation',['<','<=','>','>='])
def test_symbolic_linear_coefficients_keep_direction_and_domain(sign,relation):
    x=sx.Symbol('x');X=sp.Symbol('x')
    a=sx.Symbol('a',**{sign:True});A=sp.Symbol('a',**{sign:True})
    condition=eval(f'a*x+2 {relation} 5')
    foreign=eval(f'A*X+2 {relation} 5')
    for domain,foreign_domain in [(sx.S.Reals,sp.S.Reals),(sx.Interval(0,10),sp.Interval(0,10))]:
        actual=sx.solve_univariate_inequality(condition,x,relational=False,domain=domain)
        expected=sp.solve_univariate_inequality(foreign,X,relational=False,domain=foreign_domain)
        if domain==sx.S.Reals:
            assert actual.to_sympy()==expected
        else:
            # SymPy can retain an unevaluated intersection where native sets
            # prove a bound redundant or retain a symbolic interval endpoint.
            # Instantiate both possible endpoint orderings independently.
            for value in [sp.Rational(1,10),sp.Integer(1),sp.Integer(10)]:
                if sign=='negative':value=-value
                assert actual.to_sympy().subs(A,value)==expected.subs(A,value)


def test_an_undecided_coefficient_sign_does_not_choose_a_half_line():
    x=sx.Symbol('x');a=sx.Symbol('a',real=True)
    with pytest.raises(NotImplementedError):sx.solve_univariate_inequality(a*x<=1,x,relational=False)


def test_known_sign_symbolic_odd_power_bounds():
    x=sx.Symbol('x');p=sx.Symbol('p',positive=True)
    assert sx.solve_univariate_inequality(x**3<=p,x,relational=False)==sx.Interval(-sx.oo,p**sx.Rational(1,3))
    assert sx.solve_univariate_inequality(x**3>=-p,x,relational=False)==sx.Interval(-p**sx.Rational(1,3),sx.oo)


def test_real_odd_power_inversion_does_not_invert_complex_bases():
    x=sx.Symbol('x')
    with pytest.raises(NotImplementedError):sx.solve_univariate_inequality((sx.I*x+1)**3<2,x,relational=False)
