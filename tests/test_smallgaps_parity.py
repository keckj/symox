"""smallgaps (2026-10-30): Float Bareiss determinant, issue 17799 in
checksol and evalf, Pow's second "recognize base as E" branch and the
units step of checksol. Every expected value was produced by SymPy 1.14."""
import pytest

import symox
from symox import (Float, I, Matrix, S, Symbol, besselj, erf, erfinv, exp, gamma, log, pi,
                   sin, solve, sqrt, srepr, symbols)
from symox.physics.units import cm
from symox.solvers.solvers import checksol

x, y = symbols('x y')
p = Symbol('p', positive=True)
NS = vars(symox)

DET_CASES = [
    ('third4',
     "[[Float('0.33333333333333331', precision=53), Float('0.14285714285714285', precision=53), Float('0.10000000000000001', precision=53), Float('1.0', precision=53)], [Float('0.20000000000000001', precision=53), Float('0.1111111111111111', precision=53), Float('0.090909090909090912', precision=53), Float('0.0', precision=53)], [Float('0.29999999999999999', precision=53), Float('0.69999999999999996', precision=53), Float('0.076923076923076927', precision=53), Float('1.0', precision=53)], [Float('0.0', precision=53), Float('1.0', precision=53), Float('0.0', precision=53), Float('0.29999999999999999', precision=53)]]",
     "Float('-0.0033845265845265848', precision=53)"),
    ('n4_k0',
     "[[Float('2.594501252298171', precision=53), Float('0.72301753922074319', precision=53), Float('4.7114069269071184', precision=53), Float('-4.5438747713642904', precision=53)], [Float('-1.90723327456781', precision=53), Float('3.4835423841472881', precision=53), Float('-4.0404288381273847', precision=53), Float('3.9908301223057112', precision=53)], [Float('-2.3010367365070747', precision=53), Float('4.6746739094856569', precision=53), Float('1.0854336878508644', precision=53), Float('3.9582901328173019', precision=53)], [Float('-2.3477452577854683', precision=53), Float('-4.6641648902912038', precision=53), Float('-1.078589488231334', precision=53), Float('-0.52044200101043536', precision=53)]]",
     "Float('-224.71493273510868', precision=53)"),
    ('n4_k7',
     "[[Integer(1), Mul(Integer(2), Symbol('y')), Float('-0.72319333357737747', precision=53), Float('1.6219792311427979', precision=53)], [Float('0.36855316508885849', precision=53), Integer(1), Float('0.163043071687897', precision=53), Symbol('x')], [Symbol('x'), Symbol('x'), Mul(Integer(2), Symbol('y')), Symbol('x')], [Integer(1), Integer(1), Integer(1), Mul(Integer(2), Symbol('y'))]]",
     "Mul(Pow(Add(Mul(Float('0.73710633017771698', precision=53), Symbol('y')), Float('-1.0', precision=53)), Integer(-1)), Add(Mul(Integer(-1), Float('1.474212660355434', precision=53), Pow(Symbol('x'), Integer(2)), Pow(Symbol('y'), Integer(2))), Mul(Float('2.737106330177717', precision=53), Pow(Symbol('x'), Integer(2)), Symbol('y')), Mul(Integer(-1), Float('1.0', precision=53), Pow(Symbol('x'), Integer(2))), Mul(Float('3.4291456416419406', precision=53), Symbol('x'), Pow(Symbol('y'), Integer(3))), Mul(Integer(-1), Float('5.3905683123723609', precision=53), Symbol('x'), Pow(Symbol('y'), Integer(2))), Mul(Float('0.80315744082474094', precision=53), Symbol('x'), Symbol('y')), Mul(Float('0.26942141803164887', precision=53), Symbol('x')), Mul(Integer(-1), Float('2.173302967952246', precision=53), Pow(Symbol('y'), Integer(4))), Mul(Float('5.8968506414217359', precision=53), Pow(Symbol('y'), Integer(3))), Mul(Integer(-1), Float('5.509879248134399', precision=53), Pow(Symbol('y'), Integer(2))), Mul(Float('2.048387303593453', precision=53), Symbol('y'))))"),
    ('n5_k10',
     "[[Integer(0), Integer(0), Float('0.24777648822511588733163989672903', precision=103), Integer(0), Integer(0)], [Integer(0), Float('0.0', precision=53), Float('-0.60559674794638906902832786727231', precision=103), Integer(0), Float('-0.967917995964324928692690264142584', precision=103)], [Float('-0.283217646712740300785071667633019', precision=103), Integer(0), Float('-0.171034507857361184335331927286461', precision=103), Float('-0.976858825395079843900703053805046', precision=103), Float('-1.0', precision=53)], [Float('-2.0', precision=53), Integer(0), Float('-0.649333495608892530270850329543464', precision=103), Float('-0.136029093858790695037441764725372', precision=103), Integer(0)], [Integer(0), Float('0.742754382470087515244472342601512', precision=103), Float('0.213173570138549983710163360228762', precision=103), Float('2.0', precision=53), Integer(0)]]",
     "Float('0.34115846909435832702444646651452', precision=103)"),
    ('n5_k5',
     "[[Rational(-1, 3), Integer(0), Float('-1.432069947175056', precision=53), Integer(0), Integer(-1)], [Integer(2), Float('1.859201043537714', precision=53), Float('-1.2582761616298836', precision=53), Integer(0), Float('0.58845774570038367', precision=53)], [Integer(1), Rational(-1, 3), Integer(3), Float('-0.86545113166550625', precision=53), Float('0.055051162035588685', precision=53)], [Float('0.31681549309617463', precision=53), Float('1.4966524976108024', precision=53), Integer(0), Float('-0.927966295731395', precision=53), Integer(2)], [Integer(-3), Integer(2), Integer(2), Integer(2), Float('1.5320956283532414', precision=53)]]",
     "Float('56.2169420184903', precision=53)"),
    ('prec50',
     "[[Float('0.10000000000000000000000000000000000000000000000000007', precision=169), Float('0.29999999999999999999999999999999999999999999999999987', precision=169), Integer(1), Integer(0)], [Float('0.69999999999999999999999999999999999999999999999999947', precision=169), Integer(2), Integer(0), Integer(1)], [Integer(0), Integer(1), Float('0.20000000000000001', precision=53), Integer(3)], [Integer(1), Integer(0), Integer(1), Float('1.1000000000000000000000000000000000000000000000000011', precision=169)]]",
     "Float('5.9578000000000000032085445411667024018242955207824604', precision=169)"),
]


@pytest.mark.parametrize('name, rows, expected', DET_CASES)
@pytest.mark.parametrize('method', [None, 'bareiss'])
def test_float_bareiss_det_matches_sympy(name, rows, expected, method):
    M = Matrix(eval(rows, NS))
    d = M.det() if method is None else M.det(method=method)
    assert srepr(d) == expected


def test_issue_17799_checksol_and_solve():
    e = -erf(x**(S(1)/3))**pi + I
    s = erfinv(I**(1/pi))**3
    assert checksol(e, {x: s}) is False
    assert checksol(e, {x: s}, simplify=False) is False
    assert checksol(e, {x: s}, minimal=True) is False
    assert solve(e, x) == []


def test_issue_17799_evalf_swallows_argument_conversion_errors():
    # Function._eval_evalf leaves the call unevaluated when converting an
    # argument raises ValueError; the enclosing expression is rebuilt from
    # the evaluated arguments.
    assert str(gamma(erfinv(2)).n()) == 'gamma(erfinv(2))'
    assert str(erf(erfinv(2) + 1).n()) == 'erf(erfinv(2) + 1)'
    assert str(besselj(1, erfinv(2)).n()) == 'besselj(1, erfinv(2))'
    assert str((erf(erfinv(2) + 1)**pi).n()) == 'erf(erfinv(2) + 1)**3.14159265358979'
    assert str((pi*erf(erfinv(2) + 1)).n()) == '3.14159265358979*erf(erfinv(2) + 1)'
    assert str((gamma(erfinv(2))**pi + sqrt(2)).n()) == 'gamma(erfinv(2))**3.14159265358979 + 1.4142135623731'
    assert str((-erf(erfinv(2) + 1)**pi).n()) == '-erf(erfinv(2) + 1)**3.14159265358979'
    # the functions' own mpmath call still raises
    for f in (sin, exp):
        with pytest.raises(ValueError, match='erfinv'):
            f(erfinv(2)).n()
    with pytest.raises(ValueError, match='erfinv'):
        erfinv(I**(1/pi)).n()


@pytest.mark.parametrize('base, exponent, expected', [
    (2*I, x/log(2*I), exp(x)),
    (-2*I, x/log(-2*I), exp(x)),
    (3 + 4*I, x/log(3 + 4*I), exp(x)),
    (1 + I, 2*x/log(1 + I), exp(2*x)),
    (I*p, x/(log(-I*p) + I*pi), exp(x)),
    (-I*p, x/(log(I*p) - I*pi), exp(x)),
    (I*p, (x + 1)/(2*log(-I*p) + 2*I*pi), exp(x/2 + S(1)/2)),
    (2*I, 1/log(2*I), S.Exp1),
])
def test_pow_recognizes_shifted_log_denominator(base, exponent, expected):
    assert base**exponent == expected


def test_pow_shifted_log_needs_the_reference_sign():
    for base, exponent in [(I*p, x/(log(-I*p) - I*pi)), (-2, x/(log(2) + I*pi)),
                           (x + I, y/(log(-x - I) + I*pi)), (2*I, 3/log(2*I) + 1)]:
        assert (base**exponent).is_Pow


def test_checksol_shifted_log_mapping():
    assert checksol((2*I)**x - y, {x: log(y)/log(2*I)}, minimal=True) is True


@pytest.mark.parametrize('f, value, expected', [
    (2*cm*x, 1, False),
    (cm*(x - 1), 2, False),
    (cm*(x - 1), 1, True),
    (cm*(x - 2), Float('2.0000000001'), True),
    (cm**2*sin(x), pi, True),
    (cm*exp(x), 0, False),
])
def test_checksol_drops_unit_factors(f, value, expected):
    assert checksol(f, x, value) is expected
