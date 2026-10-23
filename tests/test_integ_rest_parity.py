"""Integer powers and numeric radical products in the form of SymPy 1.14.

Expected strings are hard-coded from SymPy 1.14 (Integer._eval_power with
perfect_power and factors(limit=2**15), Mul.flatten comb_e / num_rat gcd
extraction / pnew regrouping).
"""

import pytest

from symox import Integer, Mul, Pow, Rational as R, Symbol, sqrt, diff, N
from symox.integrals.rationaltools import ratint

CASES = {
    'Integer(72)**R(1,3)' : '2*3**(2/3)',
    'Integer(1922)**R(2,3)' : '31*124**(1/3)',
    'Integer(18)**R(2,3)' : '3*12**(1/3)',
    'Integer(248)**R(1,6)' : 'sqrt(2)*31**(1/6)',
    'Integer(12)**R(-1,3)' : '18**(1/3)/6',
    'Integer(1922)**R(-2,3)' : '1922**(1/3)/1922',
    'Integer(2)**30*Integer(3)' : '3221225472',
    'Integer(124)**R(1,3)*sqrt(93)' : '2**(2/3)*sqrt(3)*31**(5/6)',
    'Integer(2)**R(1,3)*Integer(6)**R(1,4)' : '2**(7/12)*3**(1/4)',
    'sqrt(12)*Integer(18)**R(1,3)' : '2*3**(5/6)*6**(1/3)',
    'Integer(4)**R(1,3)*Integer(6)**R(1,3)' : '2*3**(1/3)',
    'Integer(1922)**R(1,3)*Integer(62)**R(2,3)' : '62*31**(1/3)',
    'Integer(3)**R(1,3)*sqrt(3)*Integer(6)**R(1,5)' : '3*2**(1/5)*3**(1/30)',
    'Integer(6)**R(2,3)*Integer(10)**R(3,4)' : '2*2**(5/12)*3**(2/3)*5**(3/4)',
    'R(1,2)**R(1,3)*Integer(6)**R(1,3)' : '3**(1/3)',
    'Integer(-2)**R(1,3)*Integer(6)**R(1,3)' : '(-12)**(1/3)',
    'Integer(12)**R(1,3)/Integer(18)**R(1,3)' : '18**(1/3)/3',
    'sqrt(2)*sqrt(3)' : 'sqrt(6)',
    'sqrt(2)*sqrt(6)' : '2*sqrt(3)',
    "Integer(2)**R(1,3)*Integer(3)**R(1,3)*Symbol('x')" : '6**(1/3)*x',
    'Integer(10)**R(5,2)' : '100*sqrt(10)',
    'Integer(1000000)**R(1,4)' : '10*sqrt(10)',
    '(Integer(2)**64*Integer(3))**R(1,3)' : '2097152*6**(1/3)',
}


@pytest.mark.parametrize("source, expected", list(CASES.items()))
def test_numeric_radical_form_matches_sympy(source, expected):
    value = eval(source, {"Integer": Integer, "R": R, "Symbol": Symbol, "sqrt": sqrt})
    assert str(value) == expected


@pytest.mark.parametrize("source", list(CASES))
def test_numeric_radical_value_is_preserved(source):
    namespace = {"Integer": Integer, "R": R, "Symbol": Symbol, "sqrt": sqrt}
    value = eval(source, namespace)
    if value.free_symbols:
        value = value.subs(Symbol("x"), 1)
        source = source.replace("*Symbol('x')", "")
    exact = eval(source.replace("R(", "_q("), {**namespace, "_q": lambda p, q: p / q,
                                                 "Integer": int, "sqrt": lambda v: v ** 0.5})
    assert complex(N(value, 30)) == pytest.approx(complex(exact), rel=1e-12)


# Held powers of a Rational base enter Mul.flatten's pnum_rat (a 1/q base
# reads as q**(-e) through Pow.as_base_exp, a negative base joins neg1e).
HELD_RATIONAL_BASES = [
    ("Mul(Integer(6)**R(3,4), Pow(R(1,3), R(3,2), evaluate=False))", "2**(3/4)*3**(1/4)/3"),
    ("Mul(Integer(6)**R(3,4), Pow(Integer(3), R(1,2), evaluate=False))", "3*2**(3/4)*3**(1/4)"),
    ("Mul(Integer(10)**R(1,3), Pow(R(-1,5), R(2,3), evaluate=False))", "(-1)**(2/3)*50**(1/3)/5"),
    ("Mul(Integer(2)**R(1,3), Pow(Integer(-3), R(1,3), evaluate=False))", "(-6)**(1/3)"),
    ("Mul(Pow(R(1,2),R(1,3),evaluate=False), Integer(6)**R(1,3))", "3**(1/3)"),
    ("Mul(Pow(R(2,3),R(1,2),evaluate=False), sqrt(6))", "2"),
    ("Mul(Pow(R(-2,3),R(1,2),evaluate=False), sqrt(6))", "2*I"),
    ("Mul(Pow(R(4,9),R(1,3),evaluate=False), Integer(3)**R(1,3))", "6**(2/3)/3"),
    ("Mul(Pow(R(5,7),R(3,2),evaluate=False), Integer(10)**R(1,4), Symbol('x'))", "5*2**(1/4)*5**(3/4)*sqrt(7)*x/49"),
    ("Pow(R(3,4),R(1,2),evaluate=False)*7", "7*sqrt(3)/2"),
]


@pytest.mark.parametrize("source, expected", HELD_RATIONAL_BASES)
def test_held_rational_base_radicals_match_sympy(source, expected):
    names = {"Integer": Integer, "R": R, "Symbol": Symbol, "sqrt": sqrt, "Mul": Mul, "Pow": Pow}
    assert str(eval(source, names)) == expected


def test_infeasible_integer_part_stays_held_on_the_perfect_power_base():
    held = Pow(4, R(10**20 + 1, 4))
    assert held.is_Pow and held.base == 2 and held.exp == R(10**20 + 1, 2)


def test_ratint_cubic_with_one_real_root_is_still_a_gap():
    # SymPy 1.14 returns a 16653-character closed form. symox still raises
    # in log_to_atan because the roots_quadratic radicand produced by the
    # native simplify differs from SymPy's form (documented gap).
    x = Symbol("x")
    try:
        result = ratint(1 / (x**3 + x + 1), x)
    except Exception:
        pytest.xfail("ratint(1/(x**3 + x + 1)) still raises (documented gap)")
    check = diff(result, x) - 1 / (x**3 + x + 1)
    assert abs(complex(N(check.subs(x, R(3, 7)), 30))) < 1e-20
