"""Complex power branches must not be chosen from imaginary rounding noise."""
import pickle
import mpmath

import pytest
import symox as sx
import sympy as sp


def close_components(actual,expected,digits):
    actual=actual.to_sympy()
    for a,b in zip(actual.as_real_imag(),expected.as_real_imag()):
        error=abs(sp.N(a-b,digits+10))
        assert error<=abs(b)*sp.Float(10)**(-digits+5), (actual,expected,error)


@pytest.mark.parametrize('scale',[20,100,1000])
@pytest.mark.parametrize('sign',[-1,1])
@pytest.mark.parametrize('exponent',[sp.Rational(1,2),sp.Rational(-1,2),sp.Rational(3,2)])
def test_genuine_tiny_imaginary_parts_select_the_correct_branch_and_remain_visible(scale,sign,exponent):
    digits=60
    z=-1+sign*sx.I*sx.Rational(1,10**scale)
    value=z**sx.Rational(int(exponent.p),int(exponent.q))
    # Reference evalf at only 60 overall digits can give a 0.e-1000 ball
    # for the tiny component. Resolve it independently before rounding parts.
    mp=mpmath.mp.clone()
    with mp.workdps(scale+digits+40):
        expected=mp.power(mp.mpc(-1,sign*mp.power(10,-scale)),mp.mpf(int(exponent.p))/int(exponent.q))
        real=sp.Float(mp.nstr(expected.real,digits+10),digits+10)
        imaginary=sp.Float(mp.nstr(expected.imag,digits+10),digits+10)
    close_components(sx.N(value,digits),real+sp.I*imaginary,digits)


@pytest.mark.parametrize('source',['x**4-x+1','x**4+x+1','x**4-2*x+3'])
def test_quartic_radicals_keep_consistent_numeric_branches_across_precision(source):
    x=sx.Symbol('x');f=sx.sympify(source)
    roots=sx.solve(f,x,check=False)
    assert len(roots)==4
    for root in roots:
        # Explicit reference interchange is the independent oracle only.
        # Ordinary reference evalf can itself flip this cut at some working
        # precisions. Resolve its real/imaginary structure symbolically first.
        expected=sp.N(sp.expand_complex(root.to_sympy()),60)
        close_components(sx.N(root,60),expected,60)
        restored=pickle.loads(pickle.dumps(root))
        close_components(sx.N(restored,60),expected,60)
    assert len(sx.solve(f,x))==4


@pytest.mark.parametrize('scale',[30,100,500])
@pytest.mark.parametrize('sign',[-1,1])
def test_large_real_component_cannot_hide_imaginary_cancellation(scale,sign):
    digits=60;large=10**scale
    z=-1+sign*sx.I*sx.sqrt(large*large+1)-sign*sx.I*large
    mp=mpmath.mp.clone()
    with mp.workdps(2*scale+digits+40):
        delta=mp.mpf(1)/(mp.sqrt(large*large+1)+large)
        expected=mp.sqrt(mp.mpc(-1,sign*delta))
        real=sp.Float(mp.nstr(expected.real,digits+10),digits+10)
        imaginary=sp.Float(mp.nstr(expected.imag,digits+10),digits+10)
    close_components(sx.N(sx.sqrt(z),digits),real+sp.I*imaginary,digits)
