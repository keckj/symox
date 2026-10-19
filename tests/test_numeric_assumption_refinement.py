"""Numerical assumption refinement must preserve unresolved cancellation."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('scale', [0,20,100,500])
def test_scaled_numeric_identity_does_not_gain_false_nonzero_or_sign_facts(scale):
    def expression(lib):
        return lib.Integer(10)**scale*(lib.sin(3)**2+lib.cos(3)**2-1)
    native,reference=expression(sx),expression(sp)
    assert sp.simplify(native.to_sympy()) == 0
    for name in ['is_zero','is_nonzero','is_positive','is_negative']:
        assert getattr(native,name) is getattr(reference,name) is None


@pytest.mark.parametrize('digits', [50,400,1000])
def test_tiny_nonzero_numeric_trig_values_do_not_underflow_in_assumptions(digits):
    for lib in (sx,sp):
        value=lib.sin(lib.Rational(1,10**digits))
        assert value.is_zero is False
        assert value.is_positive is True
        assert value.is_negative is False


def test_stable_closed_numeric_signs_remain_supported():
    for lib in (sx,sp):
        positive=lib.pi-lib.E
        negative=lib.E-lib.pi
        assert positive.is_positive is True
        assert positive.is_zero is False
        assert negative.is_negative is True
        assert negative.is_zero is False


@pytest.mark.parametrize('digits',[50,400,1000])
def test_tiny_sine_evalf_keeps_requested_relative_precision(digits):
    value=sx.sin(sx.Rational(1,10**digits))
    actual=value.evalf(60)
    reference=sp.sin(sp.Rational(1,10**digits)).evalf(80)
    assert actual != 0
    assert abs((actual.to_sympy()-reference)/reference) < sp.Rational(1,10**55)


def test_closed_numeric_dependency_produces_a_canonical_zero_row():
    for lib in (sx,sp):
        matrix=lib.Matrix([[lib.sin(3)**2,1-lib.cos(3)**2],[1,1]])
        reduced,pivots=matrix.rref()
        assert pivots == (0,)
        assert reduced[1,:] == lib.zeros(1,2)
        reference=reduced.to_sympy() if lib is sx else reduced
        assert all(sp.simplify(v)==0 for v in reference-sp.Matrix([[1,1],[0,0]]))
        assert matrix.rank() == 1
        with pytest.raises(ValueError):
            matrix.inv(method='GE')
        if lib is sx:
            for method in [None,'DM','LU','ADJ']:
                with pytest.raises(ValueError):
                    matrix.inv(method=method)
