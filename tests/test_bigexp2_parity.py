"""Big-integer binary exponents (bigexp2): values past 2**120.

Expected strings come from SymPy 1.14 / mpmath. Where SymPy is wrong the
expected value is mpmath's at 100+ digits (see docs/algorithmic-parity.md).
"""
import pytest

from symox import *

del test

E40 = "4342944819032518276511289189166050822943"
E100 = ("4342944819032518276511289189166050822943970058036665661144537831658646"
        "492088707747292249493384317483")


@pytest.mark.parametrize("expression, expected", [
    (lambda: exp(10**40).evalf(), "9.33379023972779e+" + E40),
    (lambda: exp(10**100).evalf(), "1.53837094004017e+" + E100),
    (lambda: exp(-10**40).evalf(), "1.07137612300699e-4342944819032518276511289189166050822944"),
    (lambda: exp(10**40).evalf(50),
     "9.3337902397277908883346823236209574640758459939158e+" + E40),
    (lambda: N(exp(10**40), 30), "9.33379023972779088833468232362e+" + E40),
    (lambda: exp(10**40).n(5), "9.3338e+" + E40),
    (lambda: exp(Float('1e40')), "9.3337902397277908883346823236209574640758e+" + E40),
    (lambda: Float(2)**(2**121), "2.94672872688847e+800274995615151279296969603516431505"),
    (lambda: Float(2)**(-2**121), "3.39359368534723e-800274995615151279296969603516431506"),
    (lambda: Float(0.5)**(2**130), "4.98700317919909e-409740797754957455000048437000412930801"),
    (lambda: Float(2)**(2**200),
     "1.52912069719157e+483736552495570264612957885066036017814076813494583656293573"),
    (lambda: exp(2**200).evalf(),
     "1.32859389695976e+697884325382082936699607079697809990657056068469105515102493"),
    (lambda: Float(3)**(10**40), "1.94375671611676e+4771212547196624372950279032551153092001"),
    (lambda: Float(10)**(10**40), "1.00000000000000e+1" + "0" * 40),
    (lambda: Float('1e10000000000000000000000000000000000000'),
     "1.00000000000000e+1" + "0" * 37),
    (lambda: Float('1e' + str(10**39)), "1.00000000000000e+1" + "0" * 39),
    (lambda: Float('-1e' + str(10**39)), "-1.00000000000000e+1" + "0" * 39),
    (lambda: Float('1e' + str(10**39)) * 2, "2.00000000000000e+1" + "0" * 39),
    (lambda: Float('1e-' + str(10**50)), "1.00000000000000e-1" + "0" * 50),
    (lambda: Float('1.5e-' + str(10**40)), "1.50000000000000e-1" + "0" * 40),
    (lambda: log(Float('1e' + str(10**50))), "2.30258509299405e+50"),
    (lambda: log(exp(10**100).evalf()), "1.00000000000000e+100"),
    (lambda: sqrt(exp(10**100).evalf()),
     "3.92220721028374e+2171472409516259138255644594583025411471985029018332830572"
     "268915829323246044353873646124746692158741"),
    (lambda: exp(10**100).evalf()**2,
     "2.36658514916008e+8685889638065036553022578378332101645887940116073331322289"
     "075663317292984177415494584498986768634966"),
    (lambda: 1/exp(10**100).evalf(),
     "6.50038280087302e-4342944819032518276511289189166050822943970058036665661144"
     "537831658646492088707747292249493384317484"),
    (lambda: exp(10**100).evalf() + 1, "1.53837094004017e+" + E100),
    (lambda: Float(2)**(2**121)/Float(2)**(2**121), "1.00000000000000"),
    (lambda: exp(10**40).evalf() - exp(10**40).evalf(), "0"),
    (lambda: exp(10**40 + I).evalf(),
     "5.04306838901446e+" + E40 + " + 7.85411366501408e+" + E40 + "*I"),
    (lambda: str(exp(10**100).evalf(30)), "1.53837094004017244473315505742e+" + E100),
    # factorial/loggamma of 10**40: SymPy does not finish; mpmath at 120 digits.
    (lambda: factorial(10**40).evalf(),
     "2.68554168269385e+395657055180967481723488710810833949177076"),
    (lambda: loggamma(10**40).evalf(), "9.11034037197618e+41"),
])
def test_big_exponent_values(expression, expected):
    assert str(expression()) == expected


@pytest.mark.parametrize("expression, expected", [
    # SymPy rounds the argument or loses exponent digits; mpmath at
    # 100-300 digits gives these values.
    (lambda: gamma(Rational(10**40, 3)).evalf(),
     "3.55973154373079e+130295280877923619116846143926094265361665"),
    (lambda: gamma(Float(1e40)),
     "2.47925558014618e+395657055180967493874929847890980626811080"),
    (lambda: gamma(Float(1e100)),
     "5.54790221218382e+995657055180967497626379820570751996562306483626017186198"
     "307120554388773181385655467744058098810443467"),
    (lambda: exp(10**40).evalf()**Rational(1, 3),
     "4.53613542926576e+1447648273010839425503763063055350274314"),
])
def test_intentional_corrections(expression, expected):
    assert str(expression()) == expected


def test_mpf_tuple_and_comparisons():
    value = exp(10**40).evalf()
    assert value._mpf_[2] == 14426950408889634073599246810018921374215
    assert str(Float(value._mpf_, precision=53)) == "9.33379023972779e+" + E40
    assert value == exp(10**40).evalf()
    assert hash(value) == hash(exp(10**40).evalf())
    assert value.is_positive
    assert value * 0 == 0
    assert exp(10**100).evalf() > exp(10**99).evalf()
    assert Float(2)**(2**121) > Float(2)**(2**120)
    assert srepr(Float('1e' + str(10**39))) == \
        "Float('1.0e+1000000000000000000000000000000000000000', precision=53)"


def test_huge_integer_parts_raise_overflow():
    value = exp(10**40).evalf()
    for call in (lambda: floor(value), lambda: ceiling(value), lambda: int(value)):
        with pytest.raises(OverflowError):
            call()
    assert int(Float(2)**2000) == 2**2000


def test_powers_with_huge_integer_exponents_are_exact():
    # x**(2**k) goes through ln(x) at k + 75 bits; the logarithm's square
    # root reduction keeps its guard bits, so the power of two stays exact.
    for k in (130, 2000, 3000, 4000):
        assert (Float(2)**(2**k))._mpf_ == (0, 1, 2**k, 1)
    assert str(Float(3)**(2**4000))[:16] == "6.32032164805173"


def test_high_precision_logarithms():
    assert str(log(2).evalf(2500))[-30:] == "727257219751529082785684213208"
    assert str(log(3).evalf(1000))[-25:] == "6761050547779542293397332"
    assert str(exp(10**40).evalf(1000))[990:1003] == "15506955571e+"
