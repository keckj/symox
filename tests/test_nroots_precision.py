"""Public numerical roots must compute the precision they advertise."""
import inspect
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp
from mpmath.libmp.libhyper import NoConvergence


@pytest.mark.parametrize("digits", [15, 30, 60, 100])
@pytest.mark.parametrize("polynomial", [False, True])
def test_all_real_and_complex_roots_have_requested_accuracy(digits, polynomial):
    source = "x**10-3*x**3+x**2-1"
    actual_input = sx.sympify(source)
    if polynomial:
        actual = sx.Poly(actual_input, sx.Symbol("x")).nroots(digits, 200, True)
    else:
        actual = sx.nroots(actual_input, digits, 200, True)
    expected = sp.nroots(sp.sympify(source), n=digits+20, maxsteps=300)
    assert len(actual) == len(expected) == 10
    unused = list(expected)
    for value in actual:
        converted = value.to_sympy()
        assert all(f._prec >= sp.Float(1, digits)._prec for f in converted.atoms(sp.Float))
        closest = min(unused, key=lambda r: abs(complex(sp.N(converted-r, 20))))
        relative_error = abs(sp.N((converted-closest)/closest, digits+10))
        assert relative_error < sp.Float(10, digits+10)**(-digits+2)
        unused.remove(closest)


@pytest.mark.parametrize("source,expected", [
    ("x**3*(x-1)**4*(x**2+1)**2", [0]*3+[1]*4+[-sp.I]*2+[sp.I]*2),
    ("x**2-1/10**1000", [-sp.Rational(1,10**500),sp.Rational(1,10**500)]),
    ("x**2+1/10**1000", [-sp.I/sp.Integer(10)**500,sp.I/sp.Integer(10)**500]),
    ("x**2-2*x+1-1/10**80", [1-sp.Rational(1,10**40),1+sp.Rational(1,10**40)]),
    ("10**1000*x-1", [sp.Rational(1,10**1000)]),
])
def test_multiplicities_and_roots_outside_float_range(source, expected):
    roots = sx.nroots(sx.sympify(source), n=80, maxsteps=300)
    assert len(roots) == len(expected)
    unmatched = list(expected)
    for root in roots:
        value = root.to_sympy()
        closest = min(unmatched,key=lambda r: abs(sp.N(value-r,100)))
        if closest == 0:
            # The reference builds roots with Float._new, which maps an
            # exact zero to S.Zero.
            assert value is sp.S.Zero
        else:
            assert abs(sp.N((value-closest)/closest,100)) < sp.Rational(1,10**77)
        unmatched.remove(closest)


@pytest.mark.parametrize("polynomial", [False,True])
@pytest.mark.parametrize("steps", [0,1])
def test_iteration_exhaustion_raises_reference_exception(polynomial,steps):
    f = sx.sympify("x**10-3*x**3+x**2-1")
    call = sx.Poly(f,sx.Symbol("x")).nroots if polynomial else lambda **kw:sx.nroots(f,**kw)
    with pytest.raises(NoConvergence,match="maxsteps"):
        call(n=100,maxsteps=steps)


def test_public_signatures_and_invalid_keywords():
    assert list(inspect.signature(sx.nroots).parameters) == ["f","n","maxsteps","cleanup"]
    p = sx.Poly(sx.Symbol("x")**2-2,sx.Symbol("x"))
    assert list(inspect.signature(p.nroots).parameters) == ["n","maxsteps","cleanup"]
    for call in [p.nroots,lambda **kw:sx.nroots(p,**kw)]:
        with pytest.raises(TypeError):call(ignored_option=True)
        assert len(call(n=60,maxsteps=100,cleanup=False)) == 2


def test_precision_and_errors_work_without_reference_import():
    source = '''
import sys
class Block:
    def find_spec(self, name, *args):
        if name == 'sympy' or name.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
import symox as sx
from mpmath.libmp.libhyper import NoConvergence
x = sx.Symbol('x')
roots = sx.nroots(x**10-3*x**3+x**2-1,n=80,maxsteps=200)
assert len(roots) == 10
assert all(f._prec > 250 for r in roots for f in r.atoms(sx.Float))
try:sx.nroots(x**10-3*x**3+x**2-1,n=80,maxsteps=0)
except NoConvergence:pass
else:raise AssertionError('missing convergence error')
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,"-c",source],check=True,timeout=30)
