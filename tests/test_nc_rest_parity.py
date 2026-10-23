"""Noncommutative arguments of commutative heads, against SymPy 1.14.

SymPy declares ``is_commutative = True`` at class level for re, im, Abs,
sign, arg, DiracDelta, Heaviside, hyper, meijerg, KroneckerDelta, LeviCivita
and the integer-valued ntheory functions, so these heads are commutative even
when their argument is not.  Any assumption that implies ``commutative``
(real, complex, zero, integer, ...) is False on a noncommutative expression.

Every case runs in one child interpreter, so a crash or a runaway recursion
fails these tests instead of killing pytest.  Expected strings are the SymPy
1.14 outputs.
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
from symox.matrices.expressions.matexpr import MatrixElement
x, y, z = symbols('x y z')
r = Symbol('r', real=True)
n = Symbol('n')
i, j = symbols('i j', integer=True)
A, B = symbols('A B', commutative=False)
k = Symbol('k', integer=True, positive=True)
M = MatrixSymbol('M', k, k)
N = MatrixSymbol('N', k, k)
cases = json.loads(sys.argv[1])
for idx, src in enumerate(cases):
    print(json.dumps({"start": idx}), flush=True)
    try:
        out = str(eval(src))
    except Exception as ex:
        out = "raised " + type(ex).__name__ + ": " + str(ex)
    print(json.dumps({"done": idx, "out": out}), flush=True)
'''

# SymPy 1.14 outputs.
SYMPY = {
    'DiracDelta(B*A - 2*x).is_commutative': 'True',
    'Heaviside(A*B + x).is_commutative': 'True',
    'hyper([1],[2],A).is_commutative': 'True',
    'meijerg([[1],[]],[[],[]],A).is_commutative': 'True',
    're(A).is_commutative': 'True',
    'sign(A).is_commutative': 'True',
    'KroneckerDelta(A, x)': '0',
    'LeviCivita(A, x, y).is_commutative': 'True',
    'DiracDelta(A)*B*DiracDelta(B)': 'raised ValueError: \nFunction defined only for Real Values. Complex part: im(A)  found in A\n.',
    'Heaviside(A)*x*B*Heaviside(B)': 'raised ValueError: Function defined only for Real Values. Complex part: im(A)  found in A .',
    'hyper([1],[2],A)*B*hyper([1],[2],B)': 'hyper((1,), (2,), A)*hyper((1,), (2,), B)*B',
    'B*re(A)*re(B)': 're(A)*re(B)*B',
    're(-(B*A - 2*x))': '2*re(x) - re(B*A)',
    'im(x - A)': '-im(A) + im(x)',
    're(x - A)': '-re(A) + re(x)',
    'arg(-A*B - x)': 'arg(-(x + A*B))',
    '(x*A).is_extended_real': 'False',
    '(A + x).is_zero': 'False',
    '(A + I).is_real': 'False',
    '(A*B).is_complex': 'False',
    '(z*A).is_zero': 'False',
    '(A**2).is_positive': 'None',
    'exp(A).is_extended_nonzero': 'False',
    'sin(A).is_integer': 'False',
    '(A*B - B*A).is_zero': 'False',
    '(r*A + r).is_imaginary': 'False',
    're(A).is_extended_nonzero': 'True',
    'Abs(A).is_extended_positive': 'True',
    '(im(A) + 1).is_zero': 'False',
    '(2*im(A) + 3).is_zero': 'False',
    '(im(x*A) - 1).is_zero': 'False',
    '(x + im(A)).is_zero': 'None',
    'DiracDelta(A + I)': 'raised ValueError: \nFunction defined only for Real Values. Complex part: im(A) + 1  found\nin I + A .',
    'DiracDelta(x*A + I)': 'raised ValueError: \nFunction defined only for Real Values. Complex part: im(x*A) + 1\nfound in x*A + I .',
    'Heaviside(A + I)': 'raised ValueError: Function defined only for Real Values. Complex part: im(A) + 1  found in I + A .',
    'DiracDelta(A + 2*I)': 'raised ValueError: \nFunction defined only for Real Values. Complex part: im(A) + 2  found\nin 2*I + A .',
    'Heaviside(x*A - I)': 'raised ValueError: Function defined only for Real Values. Complex part: im(x*A) - 1  found in x*A - I .',
    'DiracDelta(A + x)': 'DiracDelta(x + A)',
    'Heaviside(A*B - B*A)': 'Heaviside(A*B - B*A)',
    'floor(I*(B*A - 2*x))': 'floor(I*(-2*x + B*A))',
    'floor(-I*A*B + x)': 'floor(x - I*A*B)',
    'ceiling(I*(B*A - 2*x))': 'ceiling(I*(-2*x + B*A))',
    'ceiling(-I*A*B + x)': 'ceiling(x - I*A*B)',
    'floor(A + 1)': '1 + floor(A)',
    'ceiling(I*A)': 'ceiling(I*A)',
    'ceiling(3/2 + 7*I/3)': '2 + 3*I',
    '-(2*re(x) - re(B*A))': '-2*re(x) + re(B*A)',
    '2*(re(A) + x)': '2*x + 2*re(A)',
    'expand(3*(re(A) + B))': '3*re(A) + 3*B',
    'str(MatrixElement(MatPow(M, 5), 1, 2))': '(M**5)[1, 2]',
    'str(MatrixElement(M*N, 1, 2))': '(M*N)[1, 2]',
    'str(MatrixElement(Inverse(M), 0, 1))': '(M**(-1))[0, 1]',
    'str(MatrixElement(-M*N, 0, 1))': '(-M*N)[0, 1]',
    'latex(MatrixElement(MatPow(M, 5), 1, 2))': '{\\left(M^{5}\\right)}_{1,2}',
    'latex(MatrixElement(M*N, 1, 2))': '{\\left(M N\\right)}_{1,2}',
    'pretty(MatrixElement(MatPow(M, 5), 1, 2))': '⎛ 5⎞      \n⎝M ⎠[1, 2]',
    'str((M**n)[0, 0])': '(M**n)[0, 0]',
    'latex((M**n)[0, 0])': '{\\left(M^{n}\\right)}_{0,0}',
    'pretty((M**n)[0, 0])': '⎛ n⎞      \n⎝M ⎠[0, 0]',
    'str((M**-1)[0, 1])': '(M**(-1))[0, 1]',
    'str(M.T[0, 1])': 'M[1, 0]',
    'str((2*M)[0, 1])': '2*M[0, 1]',
}


@pytest.fixture(scope='module')
def results():
    cases = list(SYMPY)
    env = dict(os.environ, PYTHONPATH=PYTHON_DIR)
    proc = subprocess.run(
        [sys.executable, '-c', CHILD, json.dumps(cases)],
        capture_output=True, text=True, timeout=600, env=env,
    )
    out = {}
    last = None
    for line in proc.stdout.splitlines():
        msg = json.loads(line)
        if 'start' in msg:
            last = msg['start']
        else:
            out[cases[msg['done']]] = msg['out']
    if proc.returncode != 0:
        crashed = cases[last] if last is not None else None
        out.setdefault(crashed, 'crashed with code %d' % proc.returncode)
    return out


@pytest.mark.parametrize('src', list(SYMPY))
def test_matches_sympy(results, src):
    assert results.get(src) == SYMPY[src]
