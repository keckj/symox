"""Exact native preparation avoids PSLQ only when irrationality is proved."""
import pytest
import sympy as sp
import symox as sx
from symox import _core
from symox import _numeric_recognition


def forbidden_recognition(expression):
    raise AssertionError('this exact square-root sum must not invoke PSLQ')


@pytest.mark.parametrize('source', [
    '1/(sqrt(2) + sqrt(3) + sqrt(5))',
    '1/(sqrt(2) + sqrt(3) + sqrt(5) + sqrt(7) + sqrt(11))',
    '1/(3 - 2*sqrt(2) + sqrt(3))',
    '1/(sqrt(Rational(2, 3)) + sqrt(Rational(3, 5)))',
    '1/(2*sqrt(2) + 3*sqrt(3) + sqrt(6))**2',
    '(1/(sqrt(2) + sqrt(3)) + 1/(sqrt(5) + sqrt(7)))**2',
])
def test_square_root_sums_use_exact_native_proof_instead_of_recognition(monkeypatch, source):
    expression = eval(source, vars(sx))
    reference = sp.radsimp(eval(source, vars(sp)))
    monkeypatch.setattr(_numeric_recognition, 'nsimplify', forbidden_recognition)
    result = sx.radsimp(expression)
    assert sp.simplify(sx.to_sympy(result) - reference) == 0


def test_rational_square_classes_cancel_exactly_without_numeric_recognition(monkeypatch):
    root = sx.Pow(8, sx.Rational(1, 2), evaluate=False)
    denominator = sx.Add(3, root, -2*sx.sqrt(2), evaluate=False)
    expression = sx.Pow(denominator, -1, evaluate=False)
    monkeypatch.setattr(_numeric_recognition, 'nsimplify', forbidden_recognition)
    assert _core._st_radsimp_numeric_denominators(expression) == sx.Rational(1, 3)
    assert sx.radsimp(expression) == sx.Rational(1, 3)


@pytest.mark.parametrize('symbolic', [True, False])
@pytest.mark.parametrize('max_terms', [1, 2, 3, 4, 5])
def test_rationalization_options_are_not_bypassed_by_the_proof(monkeypatch, symbolic, max_terms):
    expression = 1/(sx.sqrt(2) + sx.sqrt(3) + sx.sqrt(5))
    reference = sp.radsimp(1/(sp.sqrt(2) + sp.sqrt(3) + sp.sqrt(5)),
                           symbolic=symbolic, max_terms=max_terms)
    monkeypatch.setattr(_numeric_recognition, 'nsimplify', forbidden_recognition)
    actual = sx.radsimp(expression, symbolic=symbolic, max_terms=max_terms)
    assert sp.simplify(sx.to_sympy(actual) - reference) == 0
    if max_terms < 3:
        assert actual == expression


def test_complex_and_higher_roots_keep_the_existing_exact_recognition_fallback(monkeypatch):
    original = _numeric_recognition.nsimplify
    calls = []

    def recognize(expression):
        calls.append(expression)
        return original(expression)

    monkeypatch.setattr(_numeric_recognition, 'nsimplify', recognize)
    denominator = 3 - (-sx.S.One)**sx.Rational(2, 3) + (-sx.S.One)**sx.Rational(1, 3)
    assert sx.radsimp(1/denominator) == sx.Rational(1, 4)
    assert calls


@pytest.mark.parametrize('error_class', [ValueError, TypeError, NotImplementedError, ArithmeticError])
def test_ordinary_recognition_failure_leaves_the_exact_expression_available(monkeypatch, error_class):
    expression = 1/(sx.pi + sx.sqrt(2))
    calls = []

    def fail(base):
        calls.append(base)
        raise error_class('recognition unavailable')

    monkeypatch.setattr(_numeric_recognition, 'nsimplify', fail)
    assert _core._st_radsimp_numeric_denominators(expression) == expression
    assert len(calls) == 1


@pytest.mark.parametrize('error_class', [KeyboardInterrupt, SystemExit, GeneratorExit])
def test_interrupt_propagates_without_a_later_recognition_call(monkeypatch, error_class):
    expression = 1/(sx.pi + sx.sqrt(2)) + 1/(sx.E + sx.sqrt(3))
    calls = []
    error = error_class('stop now')

    def stop(base):
        calls.append(base)
        raise error

    monkeypatch.setattr(_numeric_recognition, 'nsimplify', stop)
    with pytest.raises(error_class) as caught:
        _core._st_radsimp_numeric_denominators(expression)
    assert caught.value is error
    assert len(calls) == 1


def test_unproved_recognition_candidate_is_not_used(monkeypatch):
    expression = 1/(sx.pi + sx.sqrt(2))
    monkeypatch.setattr(_numeric_recognition, 'nsimplify', lambda base: sx.S.One)
    assert _core._st_radsimp_numeric_denominators(expression) == expression


def test_float_and_symbolic_denominators_keep_their_previous_protocol(monkeypatch):
    x = sx.Symbol('x')
    monkeypatch.setattr(_numeric_recognition, 'nsimplify', forbidden_recognition)
    for expression in (1/(x + sx.sqrt(2)), 1/(sx.Float('1.2', 40) + sx.sqrt(2))):
        assert _core._st_radsimp_numeric_denominators(expression) == expression
    assert sx.radsimp(1/(x + sx.sqrt(2)), symbolic=False) == 1/(x + sx.sqrt(2))
