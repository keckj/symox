"""N(e, n) digits of special functions match SymPy 1.14 exactly.

The expected strings were produced by SymPy 1.14.  Where SymPy loses digits
(it evaluates some intermediate quantity at 53 bits), symox keeps the
correctly rounded value; those cases are checked against mpmath at a much
higher precision instead.
"""

import mpmath
import mpmath.libmp as libmp
import pytest

import symox as sm


def build(source):
    namespace = {name: getattr(sm, name) for name in dir(sm) if not name.startswith("_")}
    namespace["R"] = sm.Rational
    namespace["F"] = sm.Float
    return eval(source, namespace)


SAME = [
    ("beta(R(5,2), F('0.75'))", 30, '0.639008125192315843967958662549'),
    ("beta(R(5,2), F('0.75'))", 50, '0.63900812519231584396795866254950606838784675371501'),
    ("rf(R(5,2), F('0.75'))", 30, '1.91768563521219125927081887231'),
    ("rf(R(5,2), F('0.75'))", 50, '1.9176856352121912592708188723087031153190783108675'),
    ("ff(1, F('2.5'))", 30, '-0.282094791773878143474039725780'),
    ("ff(R(5,2), F('0.75'))", 50, '2.0662986923315387750510020712044700048553486611849'),
    ("binomial(R(7,2), F('0.75'))", 30, '2.86143175519053621809042799156'),
    ('lowergamma(0, R(1,3))', 15, 'oo'),
    ('lowergamma(-2, 2+I)', 30, 'oo'),
    ('lowergamma(-1, R(1,3))', 50, 'oo'),
    ('beta(0, R(1,3))', 15, 'zoo'),
    ('airyai(-2)', 50, '0.22740742820168557599192443603787379946077222541710'),
    ('airybi(-2)', 50, '-0.41230258795639848808323405461146104203453483447240'),
    ('airyaiprime(-2)', 50, '0.61825902074169104140626429133247528291577794512415'),
    ('airybiprime(-2)', 50, '0.27879516692116952268509756941098324140300059345163'),
    ('catalan(R(1,3))', 30, '0.849070020544739550674505825694'),
    ("catalan(F('0.75'))", 30, '0.899307123059882806970222191012'),
    ('jn(2, R(1,3))', 50, '0.0073487996019581070371704436256308611774868494544638'),
    ('yn(3, R(7,2))', 50, '-0.34256678144105952331337223796670050532427875111914'),
    ('polylog(R(1,3), R(1,2)+I)', 30, '-0.426397300700419281790325332025 + 0.956091212589946329657275735937*I'),
    ('lerchphi(R(2,7), 2, R(1,3))', 50, '9.1782462637616329268401856419961645791310821258783'),
    ("expint(R(5,2), F('0.75'))", 30, '0.181272829916919104959900101944'),
    ('betainc(R(3,2), 4, 0, R(1,4))', 50, '0.052095734126984126984126984126984126984126984126984'),
    ('LambertW(R(101,3), -1)', 50, '1.8334582995717139277228136979742099980374451024395 - 5.0600166190270013679384070983736646976364036899125*I'),
    ('harmonic(R(1,3), R(7,2))', 30, '0.681484382392782386858674927364'),
    ('fibonacci(R(7,2))', 50, '2.4097486115256351788090383151068552320203489048463'),
    ('bernoulli(R(1,3))', 30, '0.815860245411219410363665234741'),
    ('euler(R(7,2))', 30, '2.14177203193959170247338424878'),
    ('zeta(R(1,2)+I)', 50, '0.14393642707718906032438966648372157903562010555575 - 0.7220997435316730891261751345803249250131843953537*I'),
    ('hyper([2, R(1,3)], [R(5,2)], R(7,80))', 50, '1.0245756076477805242410154262153487119643359697932'),
    ("meijerg([[], []], [[R(5,2)], []], F('0.75'))", 30, '0.230108306946637835130052605673'),
    ('elliptic_e(R(1,3), 2+3*I/2)', 50, '0.32111670696859954998575972086259903690605019176968 - 0.009699718106201810758291779330716392140126527723477*I'),
    ('Ynm(2, 1, R(7,2), R(1,3))', 30, '0.239808314547255773317372957894 + 0.0830344801141358876769705166972*I'),
    ('Catalan', 50, '0.91596559417721901505460351493238411077414937428167'),
    ('TribonacciConstant', 50, '1.8392867552141611325518525646532866004241787460976'),
    ('exp(pi*sqrt(163))', 50, '262537412640768743.99999999999925007259719818568888'),
    ('polar_lift(-2)**R(1,3)', 30, '0.629960524947436582383605303639 + 1.09112363597172140356007261419*I'),
]


CONSTRUCTION = [
    ('beta(2, 4)', 'beta(2, 4)'),
    ('beta(3)', 'beta(3, 3)'),
    ("beta(1, F('0.75'))", '1.33333333333333'),
    ('beta(R(1,3), R(4,3))', '9/(4*catalan(1/3))'),
    ('beta(R(-1,2), R(-3,2))', '0'),
    ('beta(R(1,3), R(-7,3))', '0'),
    ('beta(-2, -3)', 'beta(-2, -3)'),
    ('beta(0, 0)', 'zoo'),
    ("beta(F('0.5'), F('1.5'))", '1.57079632679490'),
    ('beta(R(1,2))', 'beta(1/2, 1/2)'),
    ('beta(0.0, 2)', 'zoo'),
    ('beta(-1, 2)', 'beta(-1, 2)'),
    ('beta(2, 3)', '1/12'),
    ('beta(I, 2)', 'beta(I, 2)'),
    ('lowergamma(0, R(1,3))', 'lowergamma(0, 1/3)'),
]


RAISES = [
    ('laguerre(R(1,3), 2)', 'ValueError', 'Error: n should be an integer.'),
    ("laguerre(0.5, Symbol('x'))", 'ValueError', 'Error: n should be an integer.'),
    ('hermite(R(-1,2), 3)', 'ValueError', 'The index n must be nonnegative integer (got -1/2)'),
    ('hermite(-0.5, 2)', 'ValueError', 'The index n must be nonnegative integer (got -0.500000000000000)'),
    ("hermite(-1, Symbol('x'))", 'ValueError', 'The index n must be nonnegative integer (got -1)'),
    ('assoc_legendre(0, 1, R(1,3))', 'ValueError', "assoc_legendre : abs('2nd index') must be <= '1st index' (got 0, 1)"),
    ('assoc_legendre(-1, 1, R(1,3))', 'ValueError', 'assoc_legendre : 1st index must be nonnegative integer (got -1)'),
    ("assoc_legendre(1, -2, Symbol('x'))", 'ValueError', "assoc_legendre : abs('2nd index') must be <= '1st index' (got 1, 2)"),
    ('assoc_legendre(0, 1, 0)', None, '0'),
    ('assoc_legendre(-1, 0, R(1,3))', None, '1'),
    ('assoc_legendre(1.0, 2, R(1,3))', None, 'assoc_legendre(1.0, 2, 1/3)'),
    ("laguerre(-1, Symbol('x'))", None, 'exp(x)'),
    ('hermite(R(1,3), 2)', None, 'hermite(1/3, 2)'),
    ('N(lowergamma(-1, -2))', 'ValueError', 'gamma function pole'),
]


@pytest.mark.parametrize(("source", "digits", "expected"), SAME)
def test_digits_match_sympy(source, digits, expected):
    assert str(sm.N(build(source), digits)) == expected


@pytest.mark.parametrize(("source", "expected"), CONSTRUCTION)
def test_construction_matches_sympy(source, expected):
    assert str(build(source)) == expected


@pytest.mark.parametrize(("source", "error", "expected"), RAISES)
def test_domain_errors_match_sympy(source, error, expected):
    if error is None:
        assert str(build(source)) == expected
        return
    with pytest.raises(ValueError) as info:
        build(source)
    assert str(info.value) == expected


def _rounded(value, bits, digits):
    return libmp.to_str(libmp.mpf_pos(mpmath.mpf(value)._mpf_, bits, "n"), digits, strip_zeros=False)


def _catalan(n):
    return 4**n * mpmath.gamma(n + mpmath.mpf(1) / 2) / (mpmath.sqrt(mpmath.pi) * mpmath.gamma(n + 2))


def _sph(kind, order, z):
    return kind(order + mpmath.mpf(1) / 2, z) * mpmath.sqrt(mpmath.pi / (2 * z))


# SymPy loses digits here (53-bit intermediate values); symox is correctly
# rounded.  Each entry: source, digits, mpmath truth, SymPy's (wrong) string.
REFERENCE_WRONG = [
    ("catalan(I)", 20, lambda: _catalan(mpmath.mpc(0, 1)),
     "0.39764993382373624267 - 0.020884341620842555704*I"),
    ("catalan(R(1,2) + I)", 15, lambda: _catalan(mpmath.mpc(0.5, 1)),
     "0.527922517788386 + 0.189429082299752*I"),
    ("jn(2, 5.2 + 0.3*I)", 20, lambda: _sph(mpmath.besselj, 2, mpmath.mpc(5.2, 0.3)),
     "0.099419756723640347144 - 0.054525080242173564087*I"),
    ("yn(2, 5.2 + 0.3*I)", 20, lambda: _sph(mpmath.bessely, 2, mpmath.mpc(5.2, 0.3)),
     "0.18525034196069723006 + 0.014895573969924818347*I"),
    ("jn(1, F('2.5'))", 30, lambda: _sph(mpmath.besselj, 1, mpmath.mpf(2.5)),
     "0.416212989275406514206064623736"),
    ("yn(1, F('2.5'))", 30, lambda: _sph(mpmath.bessely, 1, mpmath.mpf(2.5)),
     "-0.111205879154073200368293565191"),
    ("zeta(R(101,3) + 2, 4)", 50, lambda: mpmath.zeta(mpmath.mpf(107) / 3, 4),
     "3.3626293343348525547695864359291926714752922299360e-22"),
]


@pytest.mark.parametrize(("source", "digits", "truth", "sympy_value"), REFERENCE_WRONG)
def test_more_accurate_than_sympy(source, digits, truth, sympy_value):
    bits = libmp.dps_to_prec(digits)
    with mpmath.workdps(3 * digits + 10):
        value = mpmath.mpmathify(truth())
        parts = [value.real, value.imag] if isinstance(value, mpmath.mpc) else [value]
        expected = [_rounded(part, bits, digits) for part in parts]
    got = sm.N(build(source), digits)
    re_part, im_part = got.as_real_imag()
    got_parts = [str(re_part)] + ([str(im_part)] if len(parts) == 2 else [])
    assert got_parts == expected
    assert str(got) != sympy_value


@pytest.mark.parametrize(
    ("source", "truth"),
    [
        ("jacobi(R(1,3), R(1,2), R(3,2), R(1,4))", lambda: mpmath.jacobi(mpmath.mpf(1) / 3, 0.5, 1.5, 0.25)),
        ("gegenbauer(R(1,3), R(1,3), R(1,4))", lambda: mpmath.gegenbauer(mpmath.mpf(1) / 3, mpmath.mpf(1) / 3, 0.25)),
        ("assoc_laguerre(R(5,2), R(1,3), R(7,2))", lambda: mpmath.laguerre(2.5, mpmath.mpf(1) / 3, 3.5)),
        ("gegenbauer(R(-5,4), R(1,3), R(1,4))", lambda: mpmath.gegenbauer(-1.25, mpmath.mpf(1) / 3, 0.25)),
        ("chebyshevu(R(-1,2), 3)", lambda: mpmath.chebyu(-0.5, 3)),
    ],
)
def test_noninteger_degree_follows_mpmath(source, truth):
    # SymPy truncates a non-integer degree at construction; symox evaluates
    # the analytic continuation, as mpmath does.
    with mpmath.workdps(40):
        expected = _rounded(truth(), 53, 15)
    assert str(sm.N(build(source), 15)) == expected


def test_offset_logarithmic_integral_of_float_is_correctly_rounded():
    # Li(Float) evaluates at construction at the Float's 53 bits.  SymPy
    # computes li(x) - li(2) at 53 bits and loses about 9 ulp to cancellation
    # (0.393433094864550692015825461567 at 30 digits); symox rounds correctly.
    with mpmath.workdps(100):
        expected = _rounded(mpmath.li(mpmath.mpf(2.3), offset=True), 53, 30)
    assert str(sm.N(build("Li(F('2.3'))"), 30)) == expected


def test_catalan_of_float_is_correctly_rounded():
    # catalan(Float) also evaluates at construction at 53 bits.  SymPy's
    # 4**n*gamma(n + 1/2)/(sqrt(pi)*gamma(n + 2)) loses about 10 ulp there
    # (2.59035215401932372714099983568 at 30 digits); symox rounds correctly.
    with mpmath.workdps(100):
        expected = _rounded(_catalan(mpmath.mpf(2.3)), 53, 30)
    assert str(sm.N(build("catalan(F('2.3'))"), 30)) == expected
