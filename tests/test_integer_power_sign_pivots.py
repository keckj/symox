"""Integer-power normalization prevents hidden zero pivots in raw domains."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('formula', ['-x-1', '1-x', '-x-y', '-2*I*x-Rational(1,3)', '-sin(x)-cos(x)'])
@pytest.mark.parametrize('power', [-3,-2,-1,0,1,2,3,4])
def test_integer_power_sign_is_exact_over_complex_bases(formula, power):
    value=sx.sympify(formula)
    result=value**power
    sign = -1 if power%2 else 1
    assert result == sign*sx.Pow(-value,power)
    reference=sp.sympify(formula)**power
    assert sp.cancel(result.to_sympy()-reference)==0


@pytest.mark.parametrize('offset', [0,1])
def test_large_literal_integer_parity_does_not_force_power_expansion(offset):
    x=sx.Symbol('x')
    exponent=sx.Integer(2)**80+offset
    assert (-x-1)**exponent == (-1 if offset else 1)*(x+1)**exponent


def test_noninteger_powers_and_explicit_evaluation_barriers_keep_their_branches():
    x=sx.Symbol('x')
    value=-x-1
    assert sx.Pow(value,sx.Rational(1,2)).args[0] == value
    assert sx.Pow(value,2,evaluate=False).args == (value,sx.Integer(2))


def test_raw_domain_elimination_swaps_past_an_identically_zero_complex_pivot():
    def run(lib):
        x=lib.Symbol('x',positive=True)
        a=2*lib.I*x+lib.Rational(1,3)
        matrix=lib.Matrix([[a,-a,0],[-a,a,1],[0,1,0]])
        return matrix,matrix.to_DM().to_field().inv().to_Matrix()
    matrix,result=run(sx)
    reference_matrix,reference_result=run(sp)
    assert not result.has(sx.nan,sx.zoo)
    assert all(sp.simplify(v)==0 for v in matrix.to_sympy()*result.to_sympy()-sp.eye(3))
    assert all(sp.simplify(a-b)==0 for a,b in zip(result.to_sympy(),reference_result))


def test_singular_raw_domain_matrix_does_not_produce_an_inverse():
    for lib in (sx,sp):
        x=lib.Symbol('x',positive=True)
        a=2*lib.I*x+lib.Rational(1,3)
        matrix=lib.Matrix([[a,-a],[-a,a]])
        domain_matrix=matrix.to_DM().to_field()
        assert domain_matrix.rref()[1] == (0,)
        if lib is sx:
            from symox.polys.matrices.exceptions import DMNonInvertibleMatrixError
        else:
            from sympy.polys.matrices.exceptions import DMNonInvertibleMatrixError
        with pytest.raises(DMNonInvertibleMatrixError):
            domain_matrix.inv()


@pytest.mark.parametrize('scale', [0,20,100,500])
def test_scaled_exact_trig_zero_is_not_declared_nonzero_from_a_machine_residual(scale):
    for lib in (sx,sp):
        value=lib.Integer(10)**scale*(lib.sin(3)**2+lib.cos(3)**2-1)
        matrix=lib.Matrix([[value]])
        reduced,pivots=matrix.rref()
        assert pivots == ()
        assert reduced == lib.zeros(1)
        assert all(sp.simplify(v.to_sympy() if lib is sx else v)==0 for v in reduced)
        assert matrix.rank() == 0
        for method in ['GE','ADJ']:
            with pytest.raises(ValueError):
                matrix.inv(method=method)
        if lib is sx:
            # SymPy 1.14's formal-domain/default and LU paths return an
            # inverse despite this exact identity. Rejecting it is a
            # separately recorded mathematical correction, not a match.
            for method in [None,'LU','DM']:
                with pytest.raises(ValueError):
                    matrix.inv(method=method)
        else:
            for method in [None,'LU','DM']:
                assert isinstance(matrix.inv(method=method),lib.MatrixBase)
