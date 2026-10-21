"""Precision of elementary functions near zeros, branch points and cuts.

Values are compared componentwise against mpmath on the exact argument at
much higher precision; ordering of exact rules versus numerical evaluation
of Float arguments is compared with SymPy.
"""
import mpmath
import pytest
import sympy as sp
import symox as sx

FUNCTIONS = ('exp log sin cos tan cot sec csc sinh cosh tanh coth sech csch '
             'asin acos atan acot asec acsc asinh acosh atanh acoth asech acsch').split()

# Exact arguments: tiny values, values next to branch points (+-1) and on or
# next to branch cuts, as (real, imaginary) rational pairs.
T = sp.Rational(1, 10**20)
ARGUMENTS = {
    'tiny': (T, 0), '1+t': (1 + T, 0), '1-t': (1 - T, 0), '-1+t': (-1 + T, 0),
    '2.5': (sp.Rational(5, 2), 0), '-2.5': (sp.Rational(-5, 2), 0),
    '-half': (sp.Rational(-1, 2), 0), 'i/2': (0, sp.Rational(1, 2)),
    '2.5i': (0, sp.Rational(5, 2)), '-2.5i': (0, sp.Rational(-5, 2)),
    'ti': (0, T), '1+ti': (1, T), '-2.5+ti': (sp.Rational(-5, 2), T),
    '-2.5-ti': (sp.Rational(-5, 2), -T), 't+i': (T, 1), 'half+ti': (sp.Rational(1, 2), T),
    '-3+4i': (-3, 4), '1e5+i': (10**5, 1), '-0.3-0.7i': (sp.Rational(-3, 10), sp.Rational(-7, 10)),
}


def mp_value(re, im):
    q = lambda v: mpmath.mpf(int(sp.Rational(v).p)) / int(sp.Rational(v).q)
    return mpmath.mpc(q(re), q(im))


def close(got, truth, dps):
    tol = mpmath.mpf(10)**(-dps + 1)
    for g, t in ((got.real, truth.real), (got.imag, truth.imag)):
        if t == 0:
            if g != 0 and abs(g) > abs(truth)*tol:
                return False
        elif abs(g - t) > abs(t)*tol:
            return False
    return True


@pytest.mark.parametrize('evaluate', [True, False])
@pytest.mark.parametrize('dps', [15, 40])
@pytest.mark.parametrize('name', FUNCTIONS)
def test_componentwise_accuracy_against_mpmath(name, dps, evaluate):
    failures = []
    with mpmath.workdps(dps + 60):
        for label, (re, im) in ARGUMENTS.items():
            try:
                truth = mpmath.mpc(getattr(mpmath, name)(mp_value(re, im)))
            except (ValueError, ZeroDivisionError):
                continue
            if not mpmath.isfinite(truth):
                continue
            argument = sx.Rational(str(re)) + sx.I*sx.Rational(str(im))
            value = sx.N(getattr(sx, name)(argument, evaluate=evaluate), dps)
            real, imag = value.as_real_imag()
            got = mpmath.mpc(mpmath.mpf(str(real)) if real != 0 else 0,
                             mpmath.mpf(str(imag)) if imag != 0 else 0)
            if not close(got, truth, dps):
                failures.append((label, str(value), mpmath.nstr(truth, dps)))
    assert not failures


@pytest.mark.parametrize('source, expected', [
    # log(1 + x) and its relatives keep the tiny value instead of 0.
    ('log(1 + t)', '1.00000000000000e-30'),
    ('log(1 + sqrt(2)*t)', '1.41421356237310e-30'),
    ('log(sqrt(1 + t))', '5.00000000000000e-31'),
    ('acosh(1 + t)', '1.41421356237310e-15'),
    ('acos(1 - t)', '1.41421356237310e-15'),
    ('atanh(1 - t)', '34.8853499851907'),
    # Large reductions are carried out at the precision they need.
    ('sin(exp(100))', '0.142198123658239'),
    ('cos(pi*10**40 + t)', '1.00000000000000'),
    ('sin(10**30)', '-0.0901169019121381'),
])
def test_ill_conditioned_values(source, expected):
    t = sx.Rational(1, 10**30)
    value = sx.N(eval(source, vars(sx) | {'t': t}), 15)
    assert str(value) == expected


@pytest.mark.parametrize('source', [
    'exp(0.5 + I/2)', 're(1 + 2.5*I)', 'im(0.5 + I/2)', 'conjugate(1 + 2.5*I)',
    'conjugate(0.5 + I/2)', 'sign(1 + 2.5*I)', 'sign(0.5 + I/2)', 'sec(1 + 2.5*I)',
    'csc(0.5 + I/2)', 'log(-2.5)', 'acos(2.5)', 'asin(-2.5)', 'atanh(2.5)',
    'acosh(-2.5)', 'tanh(Float("1e100", 30))', 'sinh(-2.5)',
    'exp(2.5*I)', 'atan(2.5*I)', 'Abs(1 + 2.5*I)', 'arg(-2.5)', 'floor(-2.5)',
])
def test_float_arguments_follow_eval_before_evalf(source):
    assert str(eval(source, vars(sx))) == str(eval(source, vars(sp)))


@pytest.mark.parametrize('source', [
    'exp(2*I)**Rational(1, 2)', 'exp(10*I)**Rational(1, 3)', 'exp(1 + 4*I)**Rational(3, 2)',
    'exp(pi*I*Rational(3, 4))**Rational(1, 2)', 'exp(sqrt(2)*I)**Rational(1, 2)',
    'exp(2.5*I)**0.5', 'exp(x)**n',
])
def test_exponential_powers(source):
    def run(module):
        names = vars(module) | {'x': module.Symbol('x'), 'n': module.Symbol('n', integer=True)}
        return str(eval(source, names))
    assert run(sx) == run(sp)
