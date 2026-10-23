"""elliptic_pi evaluates as fast as SymPy 1.14 and to every requested digit.

SymPy evaluates elliptic_pi with mpmath.ellippi. When the Carlson R_J
arguments leave the region where Carlson's duplication is valid (m > 1,
n > 1, complex arguments), mpmath first integrates numerically with
quadsubdiv: that takes minutes at some precisions (elliptic_pi(2,
7/2).evalf(50) took 120 s in symox, which evaluates twice, and SymPy times
out on several of these cases) and returns wrong trailing digits at others
(SymPy's 30-digit value of elliptic_pi(2, 7/2) ends in ...923255, the
correct one in ...933061). symox follows mpmath's algorithm with a
controlled initial integral (symox._ellippi). The expected values below were
computed at 150 digits and agree with mpmath at 100 digits to more than 60
digits (mpmath's own quadrature error).
"""

import time

import mpmath
import pytest

import symox as sm
from symox import _ellippi


def build(source):
    namespace = {name: getattr(sm, name) for name in dir(sm) if not name.startswith("_")}
    namespace["R"] = sm.Rational
    namespace["F"] = sm.Float
    return eval(source, namespace)


# (expression, real part, imaginary part) to 60 significant digits.
CASES = [
    ('elliptic_pi(2, R(7,2))', '-0.396420887992006712107377933061368850686618605887494350471147', '-0.551199634536912505177872719060818528710290738843519629677343'),
    ('elliptic_pi(R(1,3), 2)', '1.44749204805848066904595579774312254104996893403082237716695', '-1.74216923199344777975879208523159467900621853844651507128343'),
    ("elliptic_pi(R(5,2), F('0.75'))", '-0.500323984859079649236337307155321550855961128100726192309135', '-1.53294024990642757391335567962072831188946548292758825765982'),
    ('elliptic_pi(3, R(1,2)+I)', '0.296336217327010994471376965808742332206208729011266772036464', '-1.45231524852465331339175131129432814130758292860897286377971'),
    ('elliptic_pi(R(-1,2), 3)', '0.923217041050312111288782558172210778689375976872813240801293', '-0.900729638248629510687435036512099590480829762476166314460381'),
    ('elliptic_pi(4, 25)', '-0.0493665255073401377585179856897517158501502515550576776175919', '-0.317193428494362082985518286088342687966493951438973730460448'),
    ('elliptic_pi(R(2,3), 2+3*I/2)', '1.78015591658864006870945415051453950861426055584089260697334', '1.58980954636631010578297452476920473655268381127250296707093'),
    ('elliptic_pi(R(9,10), R(11,10))', '7.44988777366696641619737622973141900532320123007124782507865', '-11.3333350833264063319365987855717715755279866575194012905959'),
    ('elliptic_pi(2*I, R(5,2))', '1.49575907553763912757429412586971746595577882482677004391317', '-0.141027629597945060125330700111366092985629981450134797243927'),
    ('elliptic_pi(R(1,2)-I, 2-I)', '0.542202910228310349318186890569940296781601630649164900482608', '-1.13210909968573782687362413820624997459894438619325743900116'),
    ('elliptic_pi(10, R(1,2))', '-0.0523733429161062470038501037789513218150587005780333699109093', '-0.537201008942147071921403378821651006725519846542694879246712'),
    ('elliptic_pi(R(-1,3)+I, -2+I)', '0.859124572298479232150659602603802006834287704993931257571603', '0.321037402721169900815062155583837241894822215329126989969499'),
    ('elliptic_pi(R(99,100), 4)', '0.976623464490183594758934827634895584945151720421636831133272', '-9.64013909469330441208257935496312845697900765337533188630697'),
    ('elliptic_pi(R(7,3), R(7,3)+I/10)', '8.88633481357782031632508004570177516863489215133511417751387', '-8.71694271579789131204590159212041267048559579358034035603114'),
    ('elliptic_pi(R(1,3), R(3,5))', '2.44333656964714301618814014991823798889299105635805185440195', '0'),
    ('elliptic_pi(R(-5,4), R(1,6))', '1.08555042361213429362159745273673071669217802382967123167329', '0'),
    ('elliptic_pi(R(1,5), R(1,2)+I/9, R(1,3))', '0.512989829636575455322631133915738950754493532213752661312222', '0.12104880676588053503742597558721603531209735077560173892996'),
    ('elliptic_pi(R(5,2), R(7,4), R(1,2))', '-0.410015772623464215076562434689358506376579717163780779216067', '-1.43393430238636911045916170181363287805678717426101009507547'),
    ('elliptic_pi(R(1,3), 4, 2)', '4.34247614417544200713786739322936762314990680209246713150086', '-3.94758976060053314646980614865226211660988390510581575515903'),
    ('elliptic_pi(-2, R(1,2)+I, R(5,2))', '0.558181170006048114773402175116099034981000515730920615777149', '0.53345712216794732832318773255090708103719952866376124822941'),
    ('elliptic_pi(3, R(9,5), R(7,2))', '-0.214434231314751424652220841571844910600577216954716396496086', '-0.419953822665825589684586288662409302988029628479942755066167'),
    ('elliptic_pi(R(1,2)+I, 1, 2-I)', '1.30473585051012808755856372763837501480566408411666736067889', '-0.0728238893466759821289990118165341694085258236109953442806937'),
    ('elliptic_pi(4, R(2,3), 25)', '-0.0493665255073401377585179856897517158501502515550576776175919', '-0.459521423703958737733411363310316286981293598816674773238957'),
    ('elliptic_pi(R(-1,2), -5, 3)', '-2.76965112315093633386634767451663233606812793061843972240388', '2.84198474720754093561023883172275807110640037722197071577389'),
]


@pytest.mark.parametrize("digits", [15, 30, 50])
@pytest.mark.parametrize("source, real, imag", CASES)
def test_digits_and_time(source, real, imag, digits):
    expr = build(source)
    start = time.perf_counter()
    value = sm.N(expr, digits)
    elapsed = time.perf_counter() - start
    assert elapsed < 5, (source, digits, elapsed)
    with mpmath.workdps(80):
        exact = mpmath.mpc(real, imag)
        re, im = value.as_real_imag()
        got = mpmath.mpc(mpmath.mpf(str(re)), mpmath.mpf(str(im)))
        error = abs(got - exact) / abs(exact)
        assert error < mpmath.mpf(10) ** (1 - digits), (source, digits, value)


def test_reported_case_is_fast():
    expr = sm.elliptic_pi(2, sm.Rational(7, 2))
    start = time.perf_counter()
    value = expr.evalf(50)
    assert time.perf_counter() - start < 3
    assert str(value) == (
        "-0.39642088799200671210737793306136885068661860588749"
        " - 0.55119963453691250517787271906081852871029073884352*I"
    )


def test_pole_and_rewrites_unchanged():
    assert str(sm.N(sm.elliptic_pi(1, sm.Float("2.5")), 15)) == "zoo"
    assert str(sm.N(sm.elliptic_pi(1, sm.Rational(101, 3)), 30)) == "zoo"


@pytest.mark.parametrize(
    "args",
    [
        ("1/3", "3/5"),
        ("-5/4", "1/6"),
        ("0.9", "0.5"),
        ("1/5", "0.3", "1/3"),
        ("-2", "1.2", "0.7"),
    ],
)
def test_carlson_region_identical_to_mpmath(args):
    # Where mpmath applies Carlson's algorithm directly, symox calls it with
    # the same arguments: the values are bit for bit mpmath's.
    for prec in (53, 113, 200):
        with mpmath.workprec(prec):
            values = [mpmath.fraction(*map(int, a.split("/"))) + 0 if "/" in a else mpmath.mpf(a) for a in args]
            assert _ellippi.ellippi(*values) == mpmath.ellippi(*values)


@pytest.mark.parametrize(
    "args",
    [
        (2, 3.5),
        (3, 0.5 + 1j),
        (0.6666, 2 + 1.5j),
        (-1, 0.5 + 1j, 2.5),
        (4, 2.0 / 3, 25),
    ],
)
def test_precision_sweep(args):
    # Full working precision at every precision, against a 400-bit value.
    with mpmath.workprec(400):
        reference = _ellippi.ellippi(*args)
    for prec in (57, 103, 193, 270):
        with mpmath.workprec(prec):
            value = _ellippi.ellippi(*args)
        with mpmath.workprec(400):
            assert abs(value - reference) <= abs(reference) * mpmath.mpf(2) ** (4 - prec)
