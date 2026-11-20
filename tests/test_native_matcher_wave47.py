"""SymPy's ordered noncommutative Mul matcher, including binding choices."""

import pytest

import symox
from symox import _core

sympy = pytest.importorskip("sympy")


def _namespace(module):
    ns = dict(vars(module))
    ns.update(zip("ABCD", module.symbols("A B C D", commutative=False)))
    ns.update(zip(("x", "y"), module.symbols("x y")))
    for name in ("u", "v"):
        ns[name] = module.Wild(name)
    for name in ("w", "z", "r"):
        ns[name] = module.Wild(name, commutative=False)
    ns["wa"] = module.Wild("wa", commutative=False, exclude=[ns["A"]])
    ns["ws"] = module.Wild("ws", commutative=False, properties=[lambda e: e.is_Symbol])
    ns["f"] = module.Function("f", commutative=False)
    ns["g"] = module.Function("g")
    return ns


def _render(module, value):
    if value is None:
        return None
    # Preserve insertion order as well as the selected values. Properties
    # are Python callables, so use the unique Wild name as the key label.
    return [(str(key), module.srepr(val)) for key, val in value.items()]


@pytest.mark.parametrize("old", [False, True])
@pytest.mark.parametrize(
    "expr, pattern",
    [
        ("A*B", "w*B"),
        ("B*A", "w*B"),
        ("A", "w*B"),
        ("A*B", "A*B"),
        ("A*B", "C*B"),
        ("A*B", "u*B"),
        ("A*B", "u*v"),
        ("A*B", "w"),
        ("A", "A"),
        ("A", "u"),
        ("A*B*C", "w*C"),
        ("A*B*C", "A*w"),
        ("A*B*C*D", "w*z"),
        ("A*B*C*D", "w*B*z"),
        ("A*B*C*D", "A*w*D"),
        ("A*B*C", "w*z*r"),
        ("A*B", "w*z*r"),
        ("A*B*A", "w*B*w"),
        ("A*B*C", "w*B*w"),
        # The reference compares the shared prefixes of repeated-Wild
        # slices, then overwrites the binding with the later slice.
        ("A*B*A*C", "w*B*w"),
        ("A*B*C*B*D", "w*B*z"),
        ("A*B*C*D", "w*C*z"),
        ("A*B*C*D", "w*B*C*z"),
        ("A**2*B", "w*B"),
        ("A**3*B", "A*w*B"),
        ("A**3*B", "w**2*B"),
        ("A**2*B**2", "w*z"),
        ("A**2*B**2", "w**2"),
        ("x**2*A**2", "w**2"),
        ("A**-1*B", "w*B"),
        ("A**-3*B", "w*B"),
        ("A**2*B", "A**2*w"),
        ("A**2", "w**2"),
        ("A**-2", "w**-2"),
        ("2*A*B", "u*w*B"),
        ("2*x*A*B", "u*w*B"),
        ("-2*x*A*B", "u*w*B"),
        ("2*x*A*B", "2*x*w*B"),
        ("2*x*A*B", "3*x*w*B"),
        ("x*y*A*B", "u*v*w*B"),
        ("x*A*B", "x*w*B"),
        ("x*A*B", "y*w*B"),
        ("A*B", "wa*B"),
        ("C*B", "wa*B"),
        ("C*A*B", "wa*B"),
        ("A*B", "ws*B"),
        ("A*C*B", "ws*B"),
        ("f(A)*B", "f(w)*B"),
        ("g(A)*B", "g(w)*B"),
        ("(A+C)*B", "(A+w)*B"),
        ("sin(A)*B", "sin(w)*B"),
        ("x*(A+C)*B", "u*(A+w)*B"),
    ],
)
def test_noncommutative_match_is_native_and_agrees(expr, pattern, old):
    ns, ref_ns = _namespace(symox), _namespace(sympy)
    target, form = eval(expr, ns), eval(pattern, ns)
    ref_target, ref_form = eval(expr, ref_ns), eval(pattern, ref_ns)
    expected = _render(sympy, ref_target.match(ref_form, old=old))
    handled, result = _core._reference_match(target, form, old)
    assert handled, result
    assert _render(symox, result) == expected
    assert _render(symox, target.match(form, old=old)) == expected


@pytest.mark.parametrize("exponent", ["Rational(-1, 2)", "Float(-2)", "Symbol('p')"])
def test_noninteger_noncommutative_powers_keep_fallback(exponent):
    ns = _namespace(symox)
    target, pattern = eval(f"A**({exponent})*B", ns), eval("w*B", ns)
    handled, reason = _core._reference_match(target, pattern, False)
    assert not handled
    assert isinstance(reason, str)


@pytest.mark.parametrize("old", [False, True])
@pytest.mark.parametrize("exponent", ["Rational(1, 2)", "Rational(3, 2)", "Float(2)", "Float('1e-1000')"])
def test_positive_noninteger_expansion_error_is_native(exponent, old):
    errors = []
    for module in (sympy, symox):
        ns = _namespace(module)
        target, pattern = eval(f"A**({exponent})*B", ns), eval("w*B", ns)
        with pytest.raises(TypeError) as caught:
            target.match(pattern, old=old)
        errors.append(str(caught.value))
        if module is symox:
            with pytest.raises(TypeError) as native:
                _core._reference_match(target, pattern, old)
            assert str(native.value) == str(caught.value)
    assert errors[0] == errors[1]


def test_commutative_prefix_mismatch_precedes_expansion_error():
    for module in (sympy, symox):
        ns = _namespace(module)
        target = eval("2*A**Rational(1, 2)*B", ns)
        pattern = eval("3*w*B", ns)
        assert target.match(pattern) is None
        if module is symox:
            handled, result = _core._reference_match(target, pattern, False)
            assert handled and result is None


def test_collect_preserves_noncommutative_protocol_type_error():
    errors = []
    for module in (sympy, symox):
        ns = _namespace(module)
        expr = eval("x*g(sqrt(A)*B) + y*g(sqrt(A)*B)", ns)
        pattern = eval("g(w*B)", ns)
        with pytest.raises(TypeError) as caught:
            module.collect(expr, pattern)
        errors.append(str(caught.value))
    assert errors[0] == errors[1]


def test_symbolic_exponent_keeps_reference_error():
    for module in (sympy, symox):
        ns = _namespace(module)
        target, pattern = eval("A**Symbol('p')*B", ns), eval("w*B", ns)
        with pytest.raises(TypeError, match="truth value"):
            target.match(pattern)


def test_huge_power_expansion_declines_before_allocation():
    ns = _namespace(symox)
    target, pattern = eval("A**(10**100)*B", ns), eval("w*B", ns)
    handled, reason = _core._reference_match(target, pattern, False)
    assert not handled
    assert "budget" in reason


def test_noncommutative_wild_property_error_propagates():
    A, B = symox.symbols("A B", commutative=False)

    def boom(_):
        raise ValueError("noncommutative property")

    w = symox.Wild("w", commutative=False, properties=[boom])
    with pytest.raises(ValueError, match="noncommutative property"):
        _core._reference_match(A*B, w*B, False)
    with pytest.raises(ValueError, match="noncommutative property"):
        (A*B).match(w*B)
