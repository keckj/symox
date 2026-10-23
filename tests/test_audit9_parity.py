"""Audit wave 9: remaining missing features and error messages.

Expected values are hard-coded from SymPy 1.14 unless a test says the
reference is mathematically wrong (category h); those tests check the
correct value numerically instead.
"""

import warnings

import pytest

from symox import (
    Function, I, ImmutableDenseMatrix, KroneckerProduct, Matrix,
    MatrixSymbol, Rational, Sum, Symbol, assoc_laguerre, cos, diff, exp,
    galois_group, gegenbauer, harmonic, jacobi, julia_code, nsolve,
    octave_code, oo, simplify, sin, sympify, symbols, zeta, Expr, QQ,
)
from symox.matrices.expressions import combine_kronecker
from symox.tensor.functions import shape
from symox.utilities.exceptions import (
    SymPyDeprecationWarning, sympy_deprecation_warning,
)

x = Symbol("x")


GALOIS = [
    ("x**5-5*x+12", "D5", True),
    ("x**5 - 110*x**3 - 55*x**2 + 2310*x + 979", "C5", True),
    ("x**5+20*x+16", "A5", True),
    ("x**5-x-1", "S5", False),
    ("x**5+15*x+12", "M20", False),
    ("x**5-2", "M20", False),
    ("x**5+x**4-4*x**3-3*x**2+3*x+1", "C5", True),
    ("x**6+x**5+x**4+x**3+x**2+x+1", "C6", False),
    ("x**6-2", "D6", False),
    ("x**6+3*x**3+3", "G18", False),
    ("x**6-3*x**2-1", "A4", True),
    ("x**6+2*x**3-2", "G36-", False),
    ("x**6-6*x**4+9*x**2+23", "S3", False),
    ("x**6+24*x-20", "A6", True),
    ("x**6+6*x**4+2*x**3+9*x**2+6*x-4", "G36+", True),
    ("x**6-3*x**2+1", "A4 x C2", False),
    ("x**6-4*x**2-1", "S4+", True),
    ("x**6 + 2*x**5 + 3*x**4 + 4*x**3 + 5*x**2 + 6*x + 7", "PGL2(F5)", False),
    ("x**6-x**3+1", "C6", False),
    ("x**6-12*x**4+36*x**2-18", "D6", False),
    ("x**6+5*x**4-10*x**3+5*x**2+2", "S6", False),
]


@pytest.mark.parametrize("poly, name, alt", GALOIS)
def test_galois_group_degree_5_and_6(poly, name, alt):
    group, is_alt = galois_group(sympify(poly), by_name=True)
    assert group.value == name
    assert is_alt is alt


def test_algebraic_field_galois_group_doc_example():
    K = QQ.alg_field_from_poly(x**6 + x**5 + x**4 + x**3 + x**2 + x + 1)
    assert K.galois_group(by_name=True)[0].value == "C6"


def test_parameter_derivatives_of_jacobi_and_gegenbauer():
    n, a, b = symbols("n a b")
    assert str(diff(jacobi(n, a, b, x), a)) == (
        "Sum((jacobi(n, a, b, x) + (2*_k + a + b + 1)*RisingFactorial(_k + b + 1, -_k + n)"
        "*jacobi(_k, a, b, x)/((-_k + n)*RisingFactorial(_k + a + b + 1, -_k + n)))"
        "/(_k + a + b + n + 1), (_k, 0, n - 1))")
    assert str(diff(jacobi(n, a, b, x), b)) == (
        "Sum(((-1)**(-_k + n)*(2*_k + a + b + 1)*RisingFactorial(_k + a + 1, -_k + n)"
        "*jacobi(_k, a, b, x)/((-_k + n)*RisingFactorial(_k + a + b + 1, -_k + n))"
        " + jacobi(n, a, b, x))/(_k + a + b + n + 1), (_k, 0, n - 1))")
    for m in (2, 3, 4):
        assert simplify(diff(jacobi(n, a, b, x), a).subs(n, m).doit()
                        - diff(jacobi(m, a, b, x), a)) == 0
        assert simplify(diff(gegenbauer(n, a, x), a).subs(n, m).doit()
                        - diff(gegenbauer(m, a, x), a)) == 0


def test_assoc_laguerre_alpha_derivative_is_correct():
    # SymPy 1.14 returns Sum(L_k^a(x)/(n - a), (k, 0, n - 1)), which is wrong
    # (category h): the true derivative divides by n - k.
    n, a = symbols("n a")
    d = diff(assoc_laguerre(n, a, x), a)
    assert d.has(Sum)
    for m in (2, 3, 5):
        exact = diff(assoc_laguerre(m, a, x), a)
        assert simplify(d.subs(n, m).doit() - exact) == 0
        value = (d.subs(n, m).doit() - exact).subs({a: Rational(1, 3), x: Rational(2, 7)})
        assert value == 0


def test_harmonic_at_infinity_with_symbolic_order():
    m = Symbol("m", positive=True)
    assert harmonic(oo, m + 1) == zeta(m + 1)
    assert harmonic(oo, Symbol("p", nonpositive=True) + 1) == oo


def test_user_function_of_matrix_in_code_printers():
    g = Function("g")
    mat = Matrix([[1, x]])
    assert julia_code(g(x) + g(mat), user_functions={"g": "gg"}) == "gg(x) + gg([1 x])"
    assert octave_code(g(mat), user_functions={"g": "gg"}) == "gg([1 x])"
    assert g(mat).free_symbols == {x}


def test_shape_of_registered_expr_subclass_with_matrix():
    class NewExpr(Expr):
        pass

    shape.register(NewExpr)(lambda expr: shape(expr.args[0]))
    assert shape(NewExpr(ImmutableDenseMatrix([1, 2]))) == (2, 1)


def test_nsolve_secant_and_bisect_follow_findroot():
    # SymPy prints 1.4142135623731*I; symox keeps a real part of about
    # 1e-61 from rounding in the native evaluator (documented divergence).
    root = nsolve(x**2 + 2, I)
    assert str(root.as_real_imag()[1]) == "1.41421356237310"
    assert abs(root.as_real_imag()[0]) < 1e-50
    assert str(nsolve(x**2 - 2, x, (1, 2), solver="bisect")) == "1.41421356237310"
    assert str(nsolve(cos(x) - x, 1)) == "0.739085133215161"
    assert str(nsolve(exp(x) - 3, 0)) == "1.09861228866811"
    assert str(nsolve(x**3 - x - 1, 1.5, solver="secant")) == "1.32471795724475"
    assert str(nsolve(x**3 - 2*x + 1, x, (-2, 0), solver="bisect")) == "-1.61803398874989"
    assert str(nsolve(sin(x), x, (3, 4), solver="bisect", prec=30)) == "3.14159265358979323846264338328"
    root = nsolve(x**2 + 1, 1 + I)
    assert str(root.as_real_imag()[1]) == "1.00000000000000" and abs(root.as_real_imag()[0]) < 1e-50
    with pytest.raises(ValueError) as error:
        nsolve(x**2 + 2, x, 1)
    # SymPy: (4.18466446988997098217 > 2.16840434497100886801e-19); the last
    # digits of the residual depend on the evaluator's rounding.
    assert str(error.value).startswith(
        "Could not find root within given tolerance. (4.184664469889970")
    assert str(error.value).endswith(
        " > 2.16840434497100886801e-19)\nTry another starting point or tweak arguments.")


def test_deprecation_warning_text_and_once_filter():
    w = SymPyDeprecationWarning("The foo is deprecated.", deprecated_since_version="1.14",
                                active_deprecations_target="foo-target")
    assert str(w) == (
        "\n\nThe foo is deprecated.\n\nSee "
        "https://docs.sympy.org/latest/explanation/active-deprecations.html#foo-target\n"
        "for details.\n\nThis has been deprecated since SymPy version 1.14. It\n"
        "will be removed in a future version of SymPy.\n")
    assert repr(w) == ("SymPyDeprecationWarning('The foo is deprecated.', "
                       "deprecated_since_version='1.14', active_deprecations_target='foo-target')")
    with pytest.raises(TypeError):
        SymPyDeprecationWarning("m", deprecated_since_version=1.14, active_deprecations_target="t")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sympy_deprecation_warning("bar", deprecated_since_version="1.14",
                                  active_deprecations_target="bar-target")
    assert len(caught) == 1 and caught[0].category is SymPyDeprecationWarning


def test_combine_kronecker_keeps_sum_of_products():
    # SymPy 1.14 gives KroneckerProduct(A + B.T, B + A.T), dropping the cross
    # terms (category h). For 1x1 matrices a, b: 2*a*b != (a + b)**2.
    m, n = symbols("m n", integer=True)
    A = MatrixSymbol("A", m, n)
    B = MatrixSymbol("B", n, m)
    expr = KroneckerProduct(A, B) + KroneckerProduct(B.T, A.T)
    assert combine_kronecker(expr) == expr
    a, b = Matrix([[2]]), Matrix([[3]])
    lhs = KroneckerProduct(a, b).doit() + KroneckerProduct(b, a).doit()
    assert lhs == Matrix([[12]])
    assert KroneckerProduct(a + b, b + a).doit() == Matrix([[25]])
