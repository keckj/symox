"""spec7 before/after harness: python bench.py <label> [fresh|repeat] [case].

Run with PYTHONPATH pointing at the symox tree under test. ``fresh`` times
one call in this process (the driver starts a new process per sample);
``repeat`` prints the median of repeated calls. A refused call prints
``refused <ExceptionName>``.
"""
import json
import statistics
import sys
import time

import symox as s
from symox.stats import DiscreteMarkovChain, P

n, x = s.Symbol('n'), s.Symbol('x')
A, B = s.MatrixSymbol('A', n, n), s.MatrixSymbol('B', n, n)
T = s.Matrix([[s.Rational(1, 2), s.Rational(1, 2)], [s.Rational(1, 3), s.Rational(2, 3)]])


def str_chain():
    from symox.core.symbol import Str
    Y = DiscreteMarkovChain('Y', [Str('a'), Str('b')], T)
    return P(s.Eq(Y[1], Str('b')), s.Eq(Y[0], Str('a')))


def int_chain():
    Y = DiscreteMarkovChain('Y', [0, 1], T)
    return P(s.Eq(Y[1], 1), s.Eq(Y[0], 0))


CASES = {
    'sin(A)': lambda: s.sin(A),
    'cos(-A)': lambda: s.cos(-A),
    'exp(log(A))': lambda: s.exp(s.log(A)),
    'sqrt(A)': lambda: s.sqrt(A),
    'sin(x).subs(x,A)': lambda: s.sin(x).subs(x, A),
    'Function(f)(A)': lambda: s.Function('f')(A),
    '2*x*A': lambda: 2*x*A,
    'A**2*B': lambda: A**2*B,
    'ZeroMatrix**(1/2)': lambda: s.ZeroMatrix(n, n)**s.Rational(1, 2),
    'sin(x)*cos(x)+x**2': lambda: s.sin(x)*s.cos(x) + x**2,
    'srepr(Tuple)': lambda: s.srepr(s.Tuple(1, x, s.sin(x))),
    'Markov Str P': str_chain,
    'Markov int P': int_chain,
}


def run(case):
    try:
        start = time.perf_counter()
        CASES[case]()
        return time.perf_counter() - start
    except Exception as error:  # a refusal is reported, not timed
        return 'refused ' + type(error).__name__


if __name__ == '__main__':
    mode, case = sys.argv[2], sys.argv[3]
    if mode == 'fresh':
        print(json.dumps(run(case)))
    else:
        first = run(case)
        if isinstance(first, str):
            print(json.dumps(first))
        else:
            loops = 200 if case not in ('Markov Str P', 'Markov int P') else 20
            print(json.dumps(statistics.median(run(case) for _ in range(loops))))
