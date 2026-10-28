"""Public polynomial reduction must never turn unavailable storage into zero."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

OPTIONS = ["{'domain':s.QQ}", "{'domain':s.ZZ}", "{'modulus':5}",
           "{'domain':s.QQ_I}", "{'domain':s.QQ.frac_field(a)}"]


@pytest.mark.parametrize('options', OPTIONS)
@pytest.mark.parametrize('method', ['buchberger', 'f5b'])
@pytest.mark.parametrize('operation', ['reduce', 'contains', 'reduced'])
def test_public_reduction_reports_actual_storage_failure(options, method, operation):
    # G is tiny. One actual step gives y**(2**32), with quotient
    # y**(2**32-2). The nonzero high term cannot fit stored u32 exponents.
    source = f'''import symox as s
x,y,a=s.symbols('x y a')
options={options}
g=s.groebner([x-y**2],x,y,method={method!r},**options)
f=x*y**(2**32-2)
try:
    result=s.reduced(f,[x-y**2],x,y,**options) if {operation!r}=='reduced' else getattr(g,{operation!r})(f)
except NotImplementedError as error:
    assert 'exponent storage limit' in str(error)
else:
    raise AssertionError('unrepresentable intermediate returned '+repr(result))
'''
    result = subprocess.run([sys.executable, '-B', '-c', source], capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
    assert result.stderr == ''


@pytest.mark.parametrize('method', ['buchberger', 'f5b'])
def test_legal_large_remainders_and_true_zero_keep_their_meaning(method):
    x, y = sx.symbols('x y')
    g = sx.groebner([x-y**2], x, y, method=method)
    legal = y**((1 << 32)-1)
    assert g.reduce(legal) == ([0], legal)
    assert g.contains(legal) is False
    assert sx.reduced(legal, [x-y**2], x, y) == ([0], legal)
    assert g.reduce(0) == ([], 0)
    assert g.contains(0) is True
    assert sx.reduced(0, [x-y**2], x, y) == ([], 0)
    actual_zero = x*y**4-y**6
    assert g.contains(actual_zero) is True
    assert g.reduce(actual_zero) == ([y**4], 0)


def test_the_rejected_intermediate_has_an_exact_nonzero_certificate():
    x, y = sp.symbols('x y')
    n = (1 << 32)-2
    quotient = y**n
    remainder = y**(n+2)
    assert remainder != 0
    assert sp.expand(x*y**n - quotient*(x-y**2) - remainder) == 0
