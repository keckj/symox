"""Results must not depend on the Dummy indices of a process.

Fresh Dummy indices are unique across processes, so they used to change the
hash-based storage order of Add/Mul operands from run to run, and
``difference_delta`` split a product on its first stored factor. The values
below are SymPy 1.14's.
"""
import os
import subprocess
import sys

from symox import Dummy, Sum, difference_delta, sin
from symox.abc import m, n

PYTHON_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'python')

LIMIT_SCRIPT = """
from symox import Dummy, Sum
from symox.series.limitseq import limit_seq
from symox.abc import m, n
_ = [Dummy() for _ in range({dummies})]
print(limit_seq(Sum(2**m/m, (m, 1, n))/(2**n/n), n))
"""

HASH_SCRIPT = """
from symox import Dummy, Symbol
a, b = Dummy('a'), Dummy('b')
print(hash(a), hash(a + Symbol('x')*b))
"""


def _run(script, seed):
    environment = dict(os.environ, PYTHONPATH=PYTHON_PATH, PYTHONHASHSEED=str(seed))
    completed = subprocess.run([sys.executable, '-c', script], env=environment,
                               capture_output=True, text=True, timeout=120, check=True)
    return completed.stdout.strip().splitlines()[-1]


def test_limit_seq_sum_over_exponential_in_fresh_processes():
    results = {_run(LIMIT_SCRIPT.format(dummies=count), 101*count) for count in range(6)}
    assert results == {'2'}


def test_fresh_dummy_hashes_are_reproducible():
    assert len({_run(HASH_SCRIPT, seed) for seed in (1, 2, 3)}) == 1


def test_difference_delta_splits_the_canonical_first_factor():
    for _ in range(3):
        d = Dummy('d', integer=True, positive=True)
        assert difference_delta(d*Sum(2**m/m, (m, 1, d)), d) == \
            2**(d + 1) + Sum(2**m/m, (m, 1, d))
    assert difference_delta(n*sin(n), n) == (n + 1)*(sin(n + 1) - sin(n)) + sin(n)
