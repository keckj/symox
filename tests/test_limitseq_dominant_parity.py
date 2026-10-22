"""limit_seq follows the reference algorithm (dominant terms after each
Stolz-Cesaro difference step). Expected values are SymPy 1.14's."""
import os
import subprocess
import sys

from symox import (AccumBounds, Add, I, Product, Rational, S, Sum, Symbol,
                   binomial, cos, factorial, fibonacci, log, lucas, oo, pi,
                   sin, sqrt, symbols)
from symox.series.limitseq import difference_delta as dd, dominant, limit_seq
from symox.testing.pytest import raises

n, m, k = symbols('n m k', integer=True)

PYTHON_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'python')


def test_dominant_docstring():
    assert dominant(5*n**3 + 4*n**2 + n + 1, n) == 5*n**3
    assert dominant(2**n + Sum(k, (k, 0, n)), n) == 2**n


def test_reference_limit_seq():
    e = binomial(2*n, n)/Sum(binomial(2*k, k), (k, 1, n))
    assert limit_seq(e) == S(3)/4
    assert limit_seq(e, m) == e
    assert limit_seq((5*n**3 + 3*n**2 + 4)/(3*n**3 + 4*n - 5), n) == S(5)/3
    assert limit_seq(Sum(k**2*Sum(2**m/m, (m, 1, k)), (k, 1, n))/(2**n*n), n) == 4
    e = (Sum(binomial(3*k, k)*binomial(5*k, k), (k, 1, n)) /
         (binomial(3*n, n)*binomial(5*n, n)))
    assert limit_seq(e, n) == S(84375)/83351
    raises(ValueError, lambda: limit_seq(e*m))


def test_trials():
    e = Sum(k**2*Sum(2**m/m, (m, 1, k)), (k, 1, n))/(2**n*n)
    assert limit_seq(e, n, trials=2) is None
    assert limit_seq(e, n, trials=10) == 4
    assert limit_seq(n**2/(n + 1), n, 0) is None


def test_reference_alternating_sign():
    assert limit_seq((-1)**n/n**2, n) == 0
    assert limit_seq((-2)**(n + 1)/(n + 3**n), n) == 0
    assert limit_seq((2*n + (-1)**n)/(n + 1), n) == 2
    assert limit_seq(sin(pi*n), n) == 0
    assert limit_seq(cos(2*pi*n), n) == 1
    assert limit_seq((S.NegativeOne/5)**n, n) == 0
    assert limit_seq(Rational(-1, 2)**n, n) == 0
    assert limit_seq((I/3)**n, n) == 0
    assert limit_seq(sqrt(n)*(I/2)**n, n) == 0
    assert limit_seq(n**7*(I/3)**n, n) == 0
    assert limit_seq(n/(n + 1) + (I/2)**n, n) == 1
    assert limit_seq((-2)**n/3**n, n) == 0


def test_reference_accumulation_bounds():
    assert limit_seq((-1)**n, n) == AccumBounds(-1, 1)
    assert limit_seq(cos(pi*n), n) == AccumBounds(-1, 1)
    assert limit_seq(sin(pi*n/2)**2, n) == AccumBounds(0, 1)
    assert limit_seq(2*(-3)**n/(n + 3**n), n) == AccumBounds(-2, 2)
    assert limit_seq(3*n/(n + 1) + 2*(-1)**n, n) == AccumBounds(1, 5)
    assert limit_seq((-2)**n, n) == AccumBounds(-oo, oo)
    assert limit_seq(2**n + (-2)**n, n) == AccumBounds(0, oo)
    assert limit_seq((2**n + (-3)**n)/3**n, n) == AccumBounds(-1, 1)
    assert limit_seq(n**((-1)**n), n) == AccumBounds(0, oo)


def test_reference_none_results():
    # The reference cannot decide these either.
    assert limit_seq(sin(pi*n/2), n) is None
    assert limit_seq(cos(pi*n/2), n) is None
    assert limit_seq(Sum((-1)**k/k, (k, 1, n)), n) is None
    assert limit_seq(Sum(1/(n + k), (k, 1, n)), n) is None
    assert limit_seq(Sum(k**2, (k, 1, n))/Sum(k, (k, 1, n))/n, n) is None
    x = Symbol('x')
    assert limit_seq((-1)**n*x, n) is None


def test_fibonacci_and_factorials():
    assert limit_seq(fibonacci(n + 1)/fibonacci(n), n) == S.GoldenRatio
    assert limit_seq(fibonacci(n)/fibonacci(n + 2), n) == S.GoldenRatio**-2
    assert limit_seq(fibonacci(n)**2/fibonacci(2*n), n) == 1/(2*S.GoldenRatio - 1)
    assert limit_seq(fibonacci(n), n) is oo
    assert limit_seq(1/fibonacci(n), n) == 0
    assert limit_seq(5**n/factorial(n), n) == 0
    assert limit_seq(n**n/factorial(n), n) is oo
    assert limit_seq(factorial(n)/n**n, n) == 0
    assert limit_seq(Sum(factorial(k), (k, 1, n))/factorial(n), n) == 1


def test_lucas_binet_divergence():
    # The reference leaves lucas unrewritten (unevaluated Limit, or a wrong
    # 0 for lucas(n)/fibonacci(n)); the Binet form gives the exact values.
    assert limit_seq(lucas(n + 1)/lucas(n), n) == S.GoldenRatio


def test_sums():
    x, y, z = symbols('x y z')
    assert limit_seq(Sum(y**2*Sum(2**z/z, (z, 1, y)), (y, 1, x))/(2**x*x), x) == 4
    assert limit_seq(Sum(2**m/m, (m, 1, n))/(2**n/n), n) == 2
    assert limit_seq(n*Sum(2**m/m, (m, 1, n))/2**n, n) == 2
    assert limit_seq(Sum(2**m/m, (m, 1, n))/2**n, n) == 0
    assert limit_seq(Sum(m*2**m, (m, 1, n))/2**n, n) is oo
    assert limit_seq(Sum(3**k/k, (k, 1, n))/(3**n/n), n) == S(3)/2
    assert limit_seq(Sum(k**2*2**k, (k, 1, n))/(n**2*2**n), n) == 2
    assert limit_seq(Sum(k, (k, 1, n))/n**2, n) == S.Half
    assert limit_seq(Sum(sqrt(k), (k, 1, n))/n**Rational(3, 2), n) == S(2)/3
    assert limit_seq(Sum(1/k, (k, 1, n))/log(n), n) == 1
    assert limit_seq(Sum(1/k, (k, 1, n))/Sum(1/(2*k - 1), (k, 1, n)), n) == 2
    assert limit_seq(Sum(log(k), (k, 1, n))/(n*log(n)), n) == 1
    assert limit_seq(Product(1 + 1/k, (k, 1, n))/n, n) == 1


def test_difference_delta_reference():
    assert dd(Sum(1/k, (k, 1, n)), n, 5) == Add(*[1/(i + n + 1) for i in range(5)])
    assert dd(n*Sum(2**m/m, (m, 1, n)), n) == 2**(n + 1) + Sum(2**m/m, (m, 1, n))


def test_nested_sum_in_fresh_processes():
    script = ("from symox import Dummy, Sum, symbols\n"
              "from symox.series.limitseq import limit_seq\n"
              "n, m, k = symbols('n m k', integer=True)\n"
              "_ = [Dummy() for _ in range({0})]\n"
              "print(limit_seq(Sum(k**2*Sum(2**m/m, (m, 1, k)), (k, 1, n))/(2**n*n), n))\n")
    results = set()
    for count in range(4):
        environment = dict(os.environ, PYTHONPATH=PYTHON_PATH, PYTHONHASHSEED=str(7*count))
        completed = subprocess.run([sys.executable, '-c', script.format(count)], env=environment,
                                   capture_output=True, text=True, timeout=120, check=True)
        results.add(completed.stdout.strip().splitlines()[-1])
    assert results == {'4'}
