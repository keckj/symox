"""Huge binary exponents and per-component accuracy (wave 7, bigexp).

Expected strings come from SymPy 1.14 unless stated otherwise; values SymPy
does not finish (or gets wrong) were checked against mpmath at 120 digits.
"""
import pytest

from symox import (I, Float, Integer, Li, N, Rational, asin, atan, cos, exp,
                   factorial, gamma, log, sin, sqrt)


V = exp(10**30).evalf()


@pytest.mark.parametrize("expr, expected", [
    (lambda: exp(10**30).evalf(), "4.02793352347121e+434294481903251827651128918916"),
    (lambda: exp(2**100).evalf(), "2.58571457654211e+550533660660465090951150130557"),
    (lambda: Float(2)**(2**101), "5.22290549782792e+763201709380294112488717654720"),
    (lambda: exp(Float('1e30')), "4.027933523471206587446775252474e+434294481903251827651128918916"),
    (lambda: exp(-10**30).evalf(), "2.48266262134887e-434294481903251827651128918917"),
    (lambda: exp(10**30).evalf(50),
     "4.0279335234712065874467752524736433699593497063414e+434294481903251827651128918916"),
    (lambda: N(exp(10**30), 30), "4.02793352347120658744677525247e+434294481903251827651128918916"),
    (lambda: V**2, "1.62242484695032e+868588963806503655302257837833"),
    (lambda: V + 1, "4.02793352347121e+434294481903251827651128918916"),
    (lambda: 1/V, "2.48266262134887e-434294481903251827651128918917"),
    (lambda: log(V), "1.00000000000000e+30"),
    (lambda: exp(I*10**30 + 10**30).evalf(),
     "-4.01154464501621e+434294481903251827651128918916 - 3.62984890243267e+434294481903251827651128918915*I"),
    (lambda: sqrt(V), "2.00697123135116e+217147240951625913825564459458"),
    (lambda: V**Rational(1, 2), "2.00697123135116e+217147240951625913825564459458"),
    (lambda: V**0.5, "2.00697123135116e+217147240951625913825564459458"),
    (lambda: str(exp(-2**101).evalf()), "1.49568050357030e-1101067321320930181902300261115"),
    (lambda: Float(2)**(-2**101), "1.91464310509902e-763201709380294112488717654721"),
    (lambda: Float(2)**(2**119), "2.32988797140294e+200068748903787819824242400879107876"),
    (lambda: exp(2**118).evalf(), "3.34606221253199e+144319095940176960802298299824842363"),
    (lambda: exp(-2**118).evalf(), "2.98858758888195e-144319095940176960802298299824842364"),
    (lambda: V*V*V, "6.53501943034382e+1302883445709755482953386756749"),
    (lambda: V**10, "1.12414658106334e+4342944819032518276511289189166"),
    (lambda: V**1000, "1.20863285845809e+434294481903251827651128918916605"),
    (lambda: Float('1e300000000000000000000000000000')**2, "1.00000000000000e+600000000000000000000000000000"),
    (lambda: V/V, "1.00000000000000"),
    (lambda: -V, "-4.02793352347121e+434294481903251827651128918916"),
])
def test_exponents_beyond_two_to_the_hundred_match_sympy(expr, expected):
    assert str(expr()) == expected


def test_huge_float_predicates_and_parts():
    assert V.is_positive and V.is_finite
    assert (V > 1) is not False and bool(V > 1)
    assert V._mpf_[2] == 1442695040888963407359924680950
    assert V - V == 0


def test_factorial_beyond_i64_uses_spouge():
    # SymPy does not finish in 20 s; mpmath at 120 digits.
    assert str(factorial(10**30).evalf(30)) == \
        "6.22311232304258500621704086956e+29565705518096748172348871081098"


def test_gamma_of_large_rational_keeps_its_digits():
    # Intentional correction: SymPy rounds the argument first and returns
    # 9.23792600424838e+9696194754459028576565517227643. mpmath at 120
    # digits gives 5.3808811129135647083e+9696194754459028578351281059261.
    assert str(gamma(Rational(10**30, 3)).evalf()) == \
        "5.38088111291356e+9696194754459028578351281059261"


def test_rational_power_of_huge_float_uses_the_exact_exponent():
    # Intentional correction: SymPy rounds 1/3 to 53 bits before the power,
    # which moves the decimal exponent by 8e15 (5.40371947145733e+
    # 144764827301083934507941042989). mpmath on the same Float with the
    # exact exponent 1/3 gives 3.42789437569241e+144764827301083942550376306305.
    assert str(V**Rational(1, 3)) == "3.42789437569241e+144764827301083942550376306305"


@pytest.mark.parametrize("expr", [
    lambda: exp(10**40).evalf(),
    lambda: Float(2)**(2**121),
    lambda: gamma(Rational(10**40, 3)).evalf(),
    lambda: exp(10**100).evalf(),
])
def test_exponents_beyond_2_120_are_evaluated(expr):
    # Binary exponents are arbitrary-size integers since bigexp2; values
    # are checked against SymPy in tests/test_bigexp2_parity.py.
    assert expr().is_Float


def test_exponents_beyond_the_materialization_bound_stay_unevaluated():
    # Binary exponents of 65536 bits or more are not materialized
    # (documented gap); the value stays symbolic instead of a wrong number.
    assert not exp(2**70000).evalf().is_Float


@pytest.mark.parametrize("expr, expected", [
    # Same strings as SymPy 1.14.
    (lambda: ((1 + Rational(1, 10**25))**I).evalf(50), 
     "0.99999999999999999999999999999999999999999999999999 + 9.9999999999999999999999995e-26*I"),
    (lambda: Li(2 + I/10**30).evalf(), "5.20342245251402e-61 + 1.44269504088896e-30*I"),
    (lambda: log(-1 - I/10**30).evalf(20), "5.0e-61 - 3.1415926535897932385*I"),
    (lambda: (cos(1) + I*sin(Rational(1, 10**30))).evalf(), "0.54030230586814 + 1.0e-30*I"),
    (lambda: asin(2 + I/10**30).evalf(), "1.5707963267949 + 1.31695789692482*I"),
])
def test_tiny_parts_match_sympy(expr, expected):
    assert str(expr()) == expected


@pytest.mark.parametrize("expr, expected", [
    # Intentional corrections, checked against mpmath at 60 digits: SymPy
    # prints a zero of precision 1 (0.e-21*I, 0.e-200*I, 0.e-41), keeps a
    # spurious 0.e-80*I or returns oo; symox refines the small part instead.
    (lambda: ((-1 + Rational(1, 10**25))**(1 + I)).evalf(),
     "-0.0432139182637723 + 4.32139182637722e-27*I"),
    (lambda: (-I*exp(3*I/5)/2 - exp(3*I/10)*sin(Rational(3, 10)) + I/2
              + Rational(1, 10**25)).evalf(50),
     "1.0000000000000000000000000000000000000000000000000e-25"),
    (lambda: exp(I*Rational(1, 10**200)).evalf(), "1.0 + 1.0e-200*I"),
    (lambda: sqrt(-1 + I*Rational(1, 10**40)).evalf(), "5.0e-41 + 1.0*I"),
    (lambda: ((1 + I/10**6)**Rational(1, 3)).evalf(), "1.00000000000011 + 3.33333333333272e-7*I"),
    (lambda: atan(I*(1 + Rational(1, 10**20))).evalf(), "1.5707963267949 + 23.3724245202204*I"),
    (lambda: (Float(2)**(2**40))**Rational(1, 3), "4.31912023499076e+110328660180"),
])
def test_tiny_parts_are_refined_instead_of_dropped(expr, expected):
    assert str(expr()) == expected
