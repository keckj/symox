"""Source conversion must preserve degrees or fail explicitly, never truncate."""
import os
import subprocess
import sys

import pytest
import symox as sx


def test_held_polynomial_dictionary_adds_repeated_generator_exponents():
    x=sx.Symbol('native_poly_repeated')
    held=sx.Mul(x,sx.Pow(x,2,evaluate=False),evaluate=False)
    result=sx.Poly(held,x)
    assert result.as_expr()==x**3
    assert result.degree()==3


@pytest.mark.parametrize('domain',['QQ','EX'])
def test_out_of_storage_source_power_does_not_become_a_small_polynomial(domain):
    code='''
import resource
resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
import symox as s
x=s.Symbol('x')
try:
    s.Poly(x**4294967297-1,x,domain=DOMAIN)
except NotImplementedError as error:
    assert 'exponent' in str(error).lower() or 'storage' in str(error).lower()
else:
    raise AssertionError('a degree4294967297 input was accepted by32-bit storage')
'''.replace('DOMAIN',repr(domain))
    subprocess.run([sys.executable,'-c',code],check=True,timeout=15,env=dict(os.environ))


def test_held_zero_power_does_not_erase_an_earlier_factor():
    from symox import stats
    x=stats.Die('native_zero_power',6)
    held=sx.Mul(x,sx.Pow(x,0,evaluate=False),evaluate=False)
    assert stats.E(held)==sx.Rational(7,2)


@pytest.mark.parametrize('domain',['QQ','EX'])
def test_out_of_storage_held_terms_can_still_cancel_exactly(domain):
    x=sx.Symbol('source_wide_cancellation')
    power=sx.Pow(x,2**32+1,evaluate=False)
    held=sx.Add(power,sx.Mul(-1,power,evaluate=False),7,evaluate=False)
    assert sx.Poly(held,x,domain=domain).as_expr()==7


@pytest.mark.parametrize('domain',['QQ','EX'])
def test_zero_times_an_out_of_storage_held_power_is_zero(domain):
    x=sx.Symbol('source_wide_zero')
    held=sx.Mul(0,sx.Pow(x,2**32+1,evaluate=False),evaluate=False)
    assert sx.Poly(held,x,domain=domain).is_zero
