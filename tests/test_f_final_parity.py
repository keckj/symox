"""Error types and messages of the SymPy 1.14 doc examples (audit category f).

Expected texts are SymPy 1.14's, with the package name ``sympy`` replaced
by ``symox`` where a class repr is part of the message.
"""
import warnings

import pytest

import symox as sx
from symox import (Abs, Derivative, ImmutableDenseMatrix, ImmutableMatrix,
                   ImmutableSparseMatrix, Integer, O, Predicate, Q, Rational,
                   S, banded, cos, cosh, diff, log, nsolve, ones, pi, sin,
                   sqrt, sympify)
from symox.abc import x, y, z
from symox.core.function import PoleError


def test_nseries_pole_message_names_the_expression():
    with pytest.raises(PoleError, match=r"^Cannot expand sin\(log\(x\)\) around 0$"):
        Abs(sin(log(x))).nseries(x, 0, 6)
    with pytest.raises(PoleError, match=r"^Cannot expand cos\(log\(x\)\) around 0$"):
        cos(log(x)).nseries(x, 0, 3)
    with pytest.raises(PoleError, match=r"^Cannot expand sin\(2\*log\(x\)\) around 0$"):
        sin(2*log(x)).nseries(x, 0, 3)


def test_derivative_wrt_invalid_variable():
    for build, text in ((lambda: (x*y*z).diff(x*y), "x*y"),
                        (lambda: Derivative(x*y*z, x*y), "x*y"),
                        (lambda: diff(sin(x), x + 1), "x + 1")):
        with pytest.raises(ValueError) as info:
            build()
        assert str(info.value) == "\nCan't calculate derivative wrt %s." % text


def test_relational_truth_value_message():
    with pytest.raises(TypeError) as info:
        bool(x < 1)
    assert str(info.value) == "cannot determine truth value of Relational: x < 1"
    with pytest.raises(TypeError) as info:
        bool(x >= y)
    assert str(info.value) == "cannot determine truth value of Relational: x >= y"
    with pytest.raises(TypeError) as info:
        1 < x < 2
    assert str(info.value) == "cannot determine truth value of Relational: x > 1"
    assert ((x < 1) != (x < 1)) is False


def test_immutable_matrix_setitem_message():
    for cls in (ImmutableMatrix, ImmutableDenseMatrix):
        with pytest.raises(TypeError) as info:
            cls([[1, 2], [3, 4]])[0, 0] = 42
        assert str(info.value) == ("Cannot set values of "
                                   "<class 'symox.matrices.immutable.ImmutableDenseMatrix'>")
    assert ImmutableMatrix is ImmutableDenseMatrix
    with pytest.raises(TypeError) as info:
        ImmutableSparseMatrix([[1, 2], [3, 4]])[0, 0] = 42
    assert str(info.value).startswith("Cannot set values of <class '")
    assert str(info.value).endswith(".ImmutableSparseMatrix'>")


def test_undefined_predicate_cannot_be_dispatched():
    P = Predicate("P")
    assert str(P(1)) == "Q.P(1)"
    with pytest.raises(TypeError) as info:
        P.register(Integer)(lambda expr, assump: True)
    assert str(info.value) == ("<class 'symox.assumptions.assume.PredicateMeta'> "
                               "cannot be dispatched.")
    with pytest.raises(AttributeError):
        P.no_such_attribute

    class MyPrimePredicate(Predicate):
        name = 'myprime'
    Q.myprime = MyPrimePredicate()
    assert callable(Q.myprime.register(Integer)(lambda expr, assump: True))


@pytest.mark.parametrize("value, text", [
    (pi, "pi"), (x, "x"), (sqrt(2), "sqrt(2)"), (sx.I, "I"), (sx.oo, "oo"), (sx.nan, "nan")])
def test_rational_rejects_non_rational_input(value, text):
    with pytest.raises(TypeError) as info:
        Rational(value)
    assert str(info.value) == "invalid input: %s" % text
    with pytest.raises(TypeError) as info:
        Rational(value, 2)
    assert str(info.value) == "invalid input: %s" % text
    with pytest.raises(TypeError) as info:
        Rational(1, value)
    assert str(info.value) == "invalid input: %s" % text


def test_rational_accepts_what_sympy_accepts():
    assert Rational('3/7') == S(3)/7
    assert Rational(0.5) == S(1)/2
    assert Rational(S(1)/3) == S(1)/3
    assert Rational(S(1)/3, S(2)/5) == S(5)/6
    assert Rational(sx.Float(0.25)) == S(1)/4
    assert Rational(Integer(4), 6) == S(2)/3
    assert Rational(1, 0) is sx.zoo
    assert Rational(S(2)/3, 0) is sx.zoo
    assert Rational(0, 0) is sx.nan
    assert Rational(True) == 1


def test_banded_message():
    with pytest.raises(ValueError) as info:
        banded({0: 2, 1: ones(2)}, rows=5)
    assert str(info.value) == "\nsequence does not fit an integral number of times in the matrix"


def test_nsolve_reports_residual_and_tolerance():
    f = cos(x)*cosh(x) - 1
    with pytest.raises(ValueError) as info:
        nsolve(f, 3.14*100)
    text = str(info.value)
    # The residual is rounding noise of cosh(312) ~ 1e135; its digits depend
    # on the working precision (see docs/algorithmic-parity.md).
    assert text.startswith("Could not find root within given tolerance. (")
    assert text.endswith(" > 2.16840434497100886801e-19)\nTry another starting point or tweak arguments.")
    assert str(nsolve(f, 3.14*100, verify=False)) == "312.588469032184"
    with pytest.raises(ValueError) as info:
        nsolve(x**2 + 1, 1, prec=30)
    assert str(info.value).endswith(" > 1.92592994438723585305597794258492732e-34)\n"
                                    "Try another starting point or tweak arguments.")
    assert str(nsolve(cos(x) - x, 1)) == "0.739085133215161"


def test_sympify_class_arithmetic_raises_type_error():
    with pytest.raises(TypeError) as info:
        sympify("O + 1")
    assert str(info.value) == "unsupported operand type(s) for +: 'type' and 'One'"
    with pytest.raises(TypeError) as info:
        sympify("O*2")
    assert str(info.value) == "unsupported operand type(s) for *: 'type' and 'Integer'"
    with pytest.raises(TypeError) as info:
        -O
    assert str(info.value) == "bad operand type for unary -: 'type'"
    with pytest.raises(TypeError) as info:
        1 + O
    assert str(info.value) == "unsupported operand type(s) for +: 'int' and 'type'"
    with pytest.raises(TypeError) as info:
        O**2
    assert str(info.value) == "unsupported operand type(s) for ** or pow(): 'type' and 'int'"
    assert (sx.Symbol | None) is not None
    with pytest.raises(sx.SympifyError):
        sympify("x +")


def test_aesara_function_warns_before_import_error():
    try:
        import aesara  # noqa: F401
    except ImportError:
        pass
    else:
        pytest.skip("aesara installed")
    from symox.printing.aesaracode import aesara_function
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pytest.raises(ImportError, match="^Aesara is required for aesara_function$"):
            aesara_function([x], [x**2 - 1], scalar=True)
    assert any(w.category is sx.utilities.exceptions.SymPyDeprecationWarning
               and "The aesara_function function is deprecated." in str(w.message) for w in caught)


def test_warns_deprecated_sympy_did_not_warn():
    from symox.testing.pytest import Failed, warns_deprecated_sympy
    with pytest.raises(Failed) as info:
        with warns_deprecated_sympy():
            pass
    assert str(info.value) == (
        "Failed: DID NOT WARN. No warnings of type "
        "<class 'symox.utilities.exceptions.SymPyDeprecationWarning'> was emitted. "
        "The list of emitted warnings is: [].")
