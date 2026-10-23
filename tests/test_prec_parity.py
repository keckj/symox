"""Wave 6 precision items: double precision overflow, huge exponents, gamma
of large arguments and negative bases to huge integer powers.

Expected strings are SymPy 1.14 outputs, except where SymPy is wrong or
does not finish; those values come from mpmath at 120 digits.
"""
import pytest

import symox as S
from symox import Float as F, Rational as R


# SymPy 1.14 outputs. Before wave 6 symox returned inf, -inf, 0.0, oo or an
# unevaluated power for most of these, because 53-bit floats were doubles.
DOUBLE_RANGE = {
    'F(2)**1024': '1.79769313486232e+308',
    'F(2)**(2**20)': '6.74114012549907e+315652',
    'F(2)**(-2**20)': '1.48342859128146e-315653',
    'F(2)**(2**63)': '1.38093229798005e+2776511644261678566',
    'F(2)**(-2**80)': '2.11280432302698e-363922934236666733021125',
    'F(0.5)**(10**30)': '3.21346246745114e-301029995663981195213738894725',
    'F(2)**(2**100)': '2.28536769422951e+381600854690147056244358827360',
    'F(2)**1100.0': '1.35829852904939e+331',
    'F(2)**R(2049,2)': '2.54232201230729e+308',
    'F(10)**400': '1.00000000000000e+400',
    'F(2.5)**-1100': '1.84497489401773e-438',
    'F(2.5)**(-1075)': '1.63866688575507e-428',
    'F(-2)**1025': '-3.59538626972463e+308',
    'F(1e300)*F(1e300)': '1.00000000000000e+600',
    'F(1e-300)*F(1e-300)': '1.00000000000000e-600',
    'F(1e-300)/F(1e300)': '1.00000000000000e-600',
    'F(1)*10**400': '1.00000000000000e+400',
    'F(1)/10**400': '1.00000000000000e-400',
    'F(1)+10**400': '1.00000000000000e+400',
    'F(1e308)+F(1e308)': '2.00000000000000e+308',
    'F(-2)**(2**63+1)': '-2.76186459596011e+2776511644261678566',
    'F(2)**(2**63)*F(2)**(-2**63)': '1.00000000000000',
    'S.Pow(F(10), 10**25)': '1.00000000000000e+10000000000000000000000000',
    'F(-2)**F(1100.0)': '1.35829852904939e+331',
    'F(3)**(2**70)': '1.02979578587803e+563285355388211103939',
    'F(1.1)**37': '34.0039485861578',
    'F(1.1)**-37': '0.0294083493705515',
    'F(1)/R(1, 10**400)': '1.00000000000000e+400',
    'F(2.0)/R(3, 10**400)': '6.66666666666667e+399',
}


@pytest.mark.parametrize('source', sorted(DOUBLE_RANGE))
def test_double_precision_floats_keep_an_unbounded_exponent(source):
    assert str(eval(source)) == DOUBLE_RANGE[source]


def test_beyond_the_exponent_bound_stays_symbolic():
    # mpmath has no exponent bound; symox stops at 2**120 binary exponents
    # (2**100 before bigexp) and keeps the power unevaluated instead of
    # returning inf or 0.
    value = F(2)**(2**121)
    assert isinstance(value, S.Pow)
    assert S.exp(10**40).evalf().func is S.exp
    # SymPy: 2.99850408006855e+130295280877923619137326550880708245539429.
    # mpmath returns exponents past i128; symox used to raise OverflowError.
    assert S.gamma(R(10**40, 3)).evalf().func is S.gamma
    assert S.gamma(S.Float(1e40)).evalf().func is S.gamma
    assert S.factorial(10**40).evalf().func is S.factorial


# SymPy 1.14 is wrong here (gamma(10**20/3) ~ 3.1e+636286142112569524544,
# gamma(10**12/3) wrong from the 6th digit): it rounds the argument to the
# target precision before evaluating. Values from mpmath at 120 digits.
GAMMA = [
    (R(10**12, 3), 15, '4.63691262941082e+3696194754453'),
    (R(10**20, 3), 15, '2.70659425601764e+636286142112569524492'),
    (R(10**25, 3), 20, '1.1580967074899602133e+80295280877923619116846132'),
    (S.Integer(10**18), 15, '5.59707356731039e+17565705518096748163'),
    (S.Integer(5000), 15, '8.45715585321109e+16321'),
    (R(10**8, 3), 15, '5.62626374751665e+236286138'),
]


@pytest.mark.parametrize('argument, digits, expected', GAMMA)
def test_gamma_of_large_arguments_keeps_its_digits(argument, digits, expected):
    assert str(S.gamma(argument).evalf(digits)) == expected


# SymPy 1.14 does not finish these in 20 s; the exact power is real.
NEGATIVE_BASES = [
    (S.Pow(-2, 10**20 + 1), 15, '-4.73063529937574e+30102999566398119521'),
    (S.Pow(-3, 10**20), 15, '3.18266058337434e+47712125471966243729'),
    (S.Pow(-R(3, 2), 10**20 + 1), 30, '-2.01832970539558909207148728627e+17609125905568124208'),
    (S.Pow(-2, 10**6 + 1), 15, '-1.98013124585918e+301030'),
]


@pytest.mark.parametrize('power, digits, expected', NEGATIVE_BASES)
def test_negative_base_to_a_huge_integer_power_is_real(power, digits, expected):
    assert str(power.evalf(digits)) == expected


def test_factorial_of_a_large_integer_does_not_build_the_exact_product():
    # SymPy 1.14: 8.26393168833124e+5565708 (exact factorial, 3.8 s).
    assert str(S.factorial(10**6).evalf()) == '8.26393168833124e+5565708'
    assert str(S.gamma(10**5 + 1).evalf(50)) == '2.8242294079603478742934215780245355184774949260912e+456573'
