"""Expression-domain row reduction uses exact native pivot simplification."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('domain_name',['EX','EXRAW'])
@pytest.mark.parametrize('method',['auto','GJ'])
def test_expression_domain_reduction_retains_domains_and_exact_pivot_values(domain_name,method):
    for lib in (sx,sp):
        x=lib.Symbol('x',positive=True)
        a=2*lib.I*x+lib.Rational(1,3)
        matrix=lib.Matrix([[a,-a,0],[-a,a,1],[0,1,0]])
        source=matrix.to_DM().convert_to(getattr(lib,domain_name))
        reduced,pivots=source.rref(method=method)
        assert reduced.domain == source.domain
        assert pivots == (0,1,2)
        assert reduced.to_Matrix() == lib.eye(3)
        # Reference EXRAW can divide by the identically zero second pivot.
        # Use its simplifying EX field as an independent inverse oracle.
        inverse=(source if lib is sx else source.convert_to(lib.EX)).inv().to_Matrix()
        residual=matrix*inverse-lib.eye(3)
        oracle=residual.to_sympy() if lib is sx else residual
        assert all(sp.simplify(value)==0 for value in oracle)


def test_public_interpolation_has_no_new_poles_at_its_own_nodes():
    for lib in (sx,sp):
        a,b,x=lib.symbols('a b x')
        expression=lib.interpolate([(a,1),(b,2)],x)
        assert lib.simplify(expression.subs(x,a)) == 1
        assert lib.simplify(expression.subs(x,b)) == 2
        polynomial=lib.interpolating_poly(3,x)
        assert polynomial.subs(x,lib.Symbol('x0')) == lib.Symbol('y0')
