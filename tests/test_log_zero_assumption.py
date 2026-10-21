"""log(a).is_zero is (a - 1).is_zero, as log._eval_is_zero in SymPy 1.14.

With a zero symbol y, (y + 1)**x has the antiderivative x in every method
of integrate (reference test_issue_6746). Expected values from SymPy 1.14.
"""
import pytest
from symox import *

# the star import brings the library test runner, which is not a test here
del test

x, z = symbols('x z')
y = Symbol('y', zero=True)
p = Symbol('p', positive=True)


@pytest.mark.parametrize('expr, expected', [
    (log(y + 1), True), (log(p + 1), False), (log(z), None), (log(p), None),
    (log(2 + I), False), (log(x**2 + 1), None),
])
def test_log_is_zero(expr, expected):
    assert expr.is_zero is expected


@pytest.mark.parametrize('kw', [{}, {'risch': True}, {'manual': True}, {'heurisch': True}])
def test_power_of_one_plus_zero_symbol(kw):
    assert str(integrate((y + 1)**x, x, **kw)) == 'x'
