"""Odd/even function heads on noncommutative sums, against SymPy 1.14.

``erf(B*A - 2*x)`` used to crash the interpreter: the sign test of a
noncommutative sum disagreed with SymPy, and since ``-(B*A - 2*x)`` stays
``Mul(-1, B*A - 2*x)`` (a noncommutative sum is not distributed), the rewrite
``f(z) -> -f(-z)`` flipped back and forth until the native stack overflowed.

Every case runs in one child interpreter, so a segfault fails these tests
instead of killing pytest.  Expected strings are the SymPy 1.14 outputs.
Where SymPy itself recurses forever (a RecursionError on a noncommutative sum
that could extract a minus sign), Symox returns the equivalent unevaluated
``-f(-z)``/``f(-z)`` form instead; see docs/algorithmic-parity.md.
"""
import json
import os
import subprocess
import sys

import pytest

PYTHON_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'python')

CHILD = r'''
import json, sys
from symox import *
x, y, t = symbols('x y t')
p = Symbol('p', positive=True)
r = Symbol('r', real=True)
A, B = symbols('A B', commutative=False)
g = Function('g')
cases = json.loads(sys.argv[1])
for i, src in enumerate(cases):
    print(json.dumps({"start": i}), flush=True)
    try:
        value = eval(src)
        out = str(value)
    except Exception as ex:
        out = "raised " + type(ex).__name__ + ": " + str(ex)
    print(json.dumps({"done": i, "out": out}), flush=True)
'''

# SymPy 1.14 outputs.
SYMPY = {
    'erf(B*A - 2*x)': 'erf(-2*x + B*A)',
    'erf(-(B*A - 2*x))': '-erf(-2*x + B*A)',
    'erf(2*x - B*A)': 'erf(2*x - B*A)',
    'erf(A*B + x)': 'erf(x + A*B)',
    'erf(x - A)': 'erf(x - A)',
    'erf(B*A - A*B)': 'erf(-A*B + B*A)',
    'erf(-2*A + 3*B)': 'erf(-2*A + 3*B)',
    'erf(-A)': '-erf(A)',
    'erf(-A*B)': '-erf(A*B)',
    'erfc(B*A - 2*x)': 'erfc(-2*x + B*A)',
    'erfc(-(B*A - 2*x))': '2 - erfc(-2*x + B*A)',
    'erfc(x - A)': 'erfc(x - A)',
    'erfi(B*A - 2*x)': 'erfi(-2*x + B*A)',
    'erfi(-(B*A - 2*x))': '-erfi(-2*x + B*A)',
    'erfi(-2*A + 3*B)': 'erfi(-2*A + 3*B)',
    'erf2(y, B*A - 2*x)': 'erf2(y, -2*x + B*A)',
    'erf2(y, -(B*A - 2*x))': '-erf(y) - erf(-2*x + B*A)',
    'Si(B*A - 2*x)': 'Si(-2*x + B*A)',
    'Si(-(B*A - 2*x))': '-Si(-2*x + B*A)',
    'Si(I*(B*A - 2*x))': 'I*Shi(-2*x + B*A)',
    'Si(-A*B - x)': '-Si(x + A*B)',
    'Shi(-A*B - x)': '-Shi(x + A*B)',
    'fresnels(-A*B - x)': '-fresnels(x + A*B)',
    'fresnelc(-A*B - x)': '-fresnelc(x + A*B)',
    'Shi(B*A - 2*x)': 'Shi(-2*x + B*A)',
    'Shi(-(B*A - 2*x))': '-Shi(-2*x + B*A)',
    'fresnels(B*A - 2*x)': 'fresnels(-2*x + B*A)',
    'fresnels(-(B*A - 2*x))': '-fresnels(-2*x + B*A)',
    'fresnels(I*(B*A - 2*x))': '-I*fresnels(-2*x + B*A)',
    'fresnelc(-(B*A - 2*x))': '-fresnelc(-2*x + B*A)',
    'fresnelc(x - A)': 'fresnelc(x - A)',
    'besselj(0, B*A - 2*x)': 'besselj(0, -2*x + B*A)',
    'besselj(0, -(B*A - 2*x))': 'besselj(0, -2*x + B*A)',
    'besselj(0, I*(B*A - 2*x))': 'besseli(0, -2*x + B*A)',
    'besseli(1, -(B*A - 2*x))': '-besseli(1, -2*x + B*A)',
    'besseli(1, I*(B*A - 2*x))': 'I*besselj(1, -2*x + B*A)',
    'sin(B*A - 2*x)': 'sin(-2*x + B*A)',
    'sin(-(B*A - 2*x))': '-sin(-2*x + B*A)',
    'cos(-(B*A - 2*x))': 'cos(-2*x + B*A)',
    'tan(x - A)': 'tan(x - A)',
    'sinh(-(B*A - 2*x))': '-sinh(-2*x + B*A)',
    'cosh(-(B*A - 2*x))': 'cosh(-2*x + B*A)',
    'tanh(B*A - A*B)': 'tanh(-A*B + B*A)',
    'asin(-(B*A - 2*x))': '-asin(-2*x + B*A)',
    'atan(B*A - 2*x)': 'atan(-2*x + B*A)',
    'atan(-(B*A - 2*x))': '-atan(-2*x + B*A)',
    'asinh(-(B*A - 2*x))': '-asinh(-2*x + B*A)',
    'atanh(-2*A + 3*B)': 'atanh(-2*A + 3*B)',
    'exp(B*A - 2*x)': 'exp(-2*x + B*A)',
    'gamma(B*A - 2*x)': 'gamma(-2*x + B*A)',
    'loggamma(-(B*A - 2*x))': 'loggamma(-(-2*x + B*A))',
    'DiracDelta(B*A - 2*x)': 'DiracDelta(-2*x + B*A)',
    'DiracDelta(-(B*A - 2*x))': 'DiracDelta(-2*x + B*A)',
    'Heaviside(-(B*A - 2*x))': 'Heaviside(-(-2*x + B*A))',
}

# SymPy 1.14 raises RecursionError on these; Symox keeps an equivalent form.
RECURSING = {
    'erf(-A*B - x)': '-erf(-(-x - A*B))',
    'erfc(-A*B - x)': '2 - erfc(-(-x - A*B))',
    'erfi(-A*B - x)': '-erfi(-(-x - A*B))',
    'sin(-A*B - x)': '-sin(-(-x - A*B))',
    'cos(-A*B - x)': 'cos(-(-x - A*B))',
    'tan(-A*B - x)': '-tan(-(-x - A*B))',
    'sinh(-A*B - x)': '-sinh(-(-x - A*B))',
    'cosh(-A*B - x)': 'cosh(-(-x - A*B))',
    'asin(-A*B - x)': '-asin(-(-x - A*B))',
    'atan(-A*B - x)': '-atan(-(-x - A*B))',
    'asinh(-A*B - x)': '-asinh(-(-x - A*B))',
    'atanh(-A*B - x)': '-atanh(-(-x - A*B))',
    'Abs(-A*B - x)': 'Abs(-x - A*B)',
    'sinc(-A*B - x)': 'sinc(-(-x - A*B))',
}

# Real-value check of DiracDelta and Heaviside (SymPy 1.14 messages).
REAL_CHECK = {
    'DiracDelta(-A)':
        'raised ValueError: \nFunction defined only for Real Values. Complex part: -im(A)  found in\n-A .',
    'DiracDelta(2*A)':
        'raised ValueError: \nFunction defined only for Real Values. Complex part: 2*im(A)  found in\n2*A .',
    'Heaviside(-A*B)':
        'raised ValueError: Function defined only for Real Values. Complex part: -im(A*B)  found in -A*B .',
    'Heaviside(x*A)':
        'raised ValueError: Function defined only for Real Values. Complex part: im(x*A)  found in x*A .',
    'DiracDelta(r + I)':
        'raised ValueError: \nFunction defined only for Real Values. Complex part: 1  found in r + I\n.',
    'Heaviside(1 + I)':
        'raised ValueError: Function defined only for Real Values. Complex part: 1  found in 1 + I .',
    'DiracDelta(r*A)': 'DiracDelta(r*A)',
    'Heaviside(x + I)': 'Heaviside(x + I)',
    'DiracDelta(I*x)': 'DiracDelta(I*x)',
    'DiracDelta(x)': 'DiracDelta(x)',
    'Heaviside(-x)': 'Heaviside(-x)',
    'DiracDelta(-x, 1)': '-DiracDelta(x, 1)',
    'DiracDelta(p)': '0',
    'Heaviside(p)': '1',
}

# The probe-628 argument.  SymPy 1.14 prints the complex part as
# ``im(t*sinh(y)*Derivative(g(x, y), x, y)*B)*atan(p)**(3/2)``: there ``im`` of
# anything is commutative, while Symox keeps ``im`` of a noncommutative
# argument noncommutative, so the factor order of the message differs (a
# documented gap).  The error itself matches.
PROBE_628 = 'DiracDelta(B*(t*(atan(p))**Rational(3,2))*Derivative(g(x, y), x, y)*sinh(y))'

ALL = {**SYMPY, **RECURSING, **REAL_CHECK, PROBE_628: None}


@pytest.fixture(scope='module')
def outputs():
    cases = list(ALL)
    env = dict(os.environ, PYTHONPATH=PYTHON_DIR + os.pathsep + os.environ.get('PYTHONPATH', ''))
    results, start = {}, 0
    while start < len(cases):
        chunk = cases[start:]
        proc = subprocess.run([sys.executable, '-c', CHILD, json.dumps(chunk)], env=env,
                              capture_output=True, text=True, timeout=300)
        current = None
        for line in proc.stdout.splitlines():
            event = json.loads(line)
            if 'start' in event:
                current = event['start']
            else:
                results[chunk[event['done']]] = event['out']
                current = None
        if current is None:
            break
        results[chunk[current]] = f'crashed with return code {proc.returncode}'
        start += current + 1
    return results


@pytest.mark.parametrize('src', list(SYMPY))
def test_matches_sympy(outputs, src):
    assert outputs[src] == SYMPY[src]


@pytest.mark.parametrize('src', list(RECURSING))
def test_sympy_recursion_cases_terminate(outputs, src):
    assert outputs[src] == RECURSING[src]


@pytest.mark.parametrize('src', list(REAL_CHECK))
def test_real_value_check(outputs, src):
    assert outputs[src] == REAL_CHECK[src]


def test_probe_628_dirac_delta_raises(outputs):
    out = outputs[PROBE_628]
    assert out.startswith('raised ValueError: \nFunction defined only for Real Values. Complex part:\n'), out
    assert 'im(t*sinh(y)*Derivative(g(x, y), x, y)*B)' in out
    assert out.endswith('found in\nt*sinh(y)*atan(p)**(3/2)*Derivative(g(x, y), x, y)*B .'), out


def test_sign_test_is_antisymmetric_on_noncommutative_sums():
    code = r'''
from symox import *
x = Symbol('x')
A, B = symbols('A B', commutative=False)
for z in [B*A - 2*x, x - A, B*A - A*B, -2*A + 3*B, A*B + x]:
    assert z.could_extract_minus_sign() is False, z
    assert (-z).could_extract_minus_sign() is True, z
print('ok')
'''
    env = dict(os.environ, PYTHONPATH=PYTHON_DIR + os.pathsep + os.environ.get('PYTHONPATH', ''))
    proc = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == 'ok'
