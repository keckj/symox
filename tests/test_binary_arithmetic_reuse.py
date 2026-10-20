"""Repeated exact arithmetic preserves declarations, protocols and precision."""
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier
import pytest
import symox as sx
import sympy as sp


def test_equal_symbol_declarations_keep_their_original_presentations():
    symbols = [sx.Symbol('declared',real=True),
               sx.Symbol('declared',real=True,commutative=True)]
    assert symbols[0] == symbols[1]
    assert sx.srepr(symbols[0]) != sx.srepr(symbols[1])
    for symbol in symbols*3:
        power = symbol**17
        product = sx.Rational(3,7)*power
        total = symbol+product
        for value in [power, product, total]:
            assert {sx.srepr(atom) for atom in value.free_symbols} == {sx.srepr(symbol)}


@pytest.mark.parametrize('digits', [30,100,300,30])
def test_float_coefficients_keep_their_precision_after_exact_cache_warming(digits):
    x = sx.Symbol('x')
    for _ in range(3):
        sx.Integer(2)*x**7 + x
    coefficient = sx.Float('1.25',digits)
    X, C = sp.Symbol('x'),sp.Float('1.25',digits)
    for value, expected in zip([coefficient*x, coefficient+x, (coefficient*x)**2],
                               [C*X, C+X, (C*X)**2]):
        assert value.atoms(sx.Float)
        assert all(atom._prec >= coefficient._prec for atom in value.atoms(sx.Float))
        assert value.to_sympy() == expected


@pytest.mark.parametrize('bits', [64,255,257,4096])
def test_large_exact_coefficients_retain_their_full_domain(bits):
    x, X = sx.Symbol('x'),sp.Symbol('x')
    c = 2**bits+1
    for _ in range(3):
        actual = c*x**3+sx.Rational(1,c)*x
        expected = c*X**3+sp.Rational(1,c)*X
        assert actual.to_sympy() == expected


def test_arithmetic_identities_keep_the_current_python_operand():
    for _ in range(4):
        x = sx.Symbol('identity')
        assert x+0 is x
        assert 0+x is x
        assert x*1 is x
        assert 1*x is x
        assert x**1 is x


def test_foreign_conversion_runs_on_every_operator_call():
    x = sx.Symbol('x')
    calls = []
    class Foreign:
        def _symox_(self):
            calls.append(True)
            return sx.Integer(len(calls))
    value = Foreign()
    for i in range(1,6):
        assert x+value == x+i
    assert len(calls) == 5


def test_held_operands_keep_argument_order_after_exact_cache_warming():
    x = sx.Symbol('x')
    x+x
    with sx.evaluate(False):
        assert (x+x).args == (x,x)
        assert (x*1).args == (x,sx.Integer(1))


def test_evaluation_contexts_remain_independent_across_threads():
    x = sx.Symbol('thread_binary')
    barrier = Barrier(2)
    def run(enabled):
        with sx.evaluate(enabled):
            barrier.wait(timeout=5)
            for _ in range(30):
                value = x+x
                if enabled:
                    assert value == 2*x
                else:
                    assert value.args == (x,x)
        return True
    with ThreadPoolExecutor(2) as pool:
        futures = [pool.submit(run,value) for value in [True,False]]
        assert all(f.result(timeout=10) for f in futures)
