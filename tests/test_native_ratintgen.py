"""Native ratint_ratpart / ratint_logpart outside the QQ/ZZ kernel.

Expected strings are SymPy 1.14 outputs. Each native result is also
compared with the reference Python body (native entry points disabled) and
checked by a derivative certificate.
"""
import importlib

import pytest
import symox as sx

rt = importlib.import_module('symox.integrals.rationaltools')
core = importlib.import_module('symox._core')
x, a, b = sx.symbols('x a b')
NAMES = {'x': x, 'a': a, 'b': b}


def poly(text, domain):
    return sx.Poly(sx.expand(sx.sympify(text, locals=NAMES)), x, domain=domain)


def no_native(monkeypatch):
    monkeypatch.setattr(core, '_ratint_ratpart', lambda *args: None)
    monkeypatch.setattr(core, '_ratint_logpart', lambda *args: None)


# (numerator, denominator, domain, SymPy 1.14 str(ratint_ratpart(...)))
RATPART = [
    ('1', '(x**2+I)**2', 'QQ_I', '(-I*x/(2*x**2 + 2*I), -I/(2*x**2 + 2*I))'),
    ('I', '(x**2+1)**2', 'QQ_I', '(I*x/(2*x**2 + 2), I/(2*x**2 + 2))'),
    ('x+I', '(x**3-I)*(x-1)**2', 'QQ_I',
     '(-I/(x - 1), (-I*x - 1 - I)/(x**3 + x**2*(-1 - I) + x*(-1 + I) + 1))'),
    ('1', '(x-I)**3*(x+2)', 'QQ_I',
     '((2*x - 2 - 3*I)/(x**2*(6 + 8*I) + x*(16 - 12*I) - 6 - 8*I), '
     '1/(x**2*(3 + 4*I) + x*(10 + 5*I) + 8 - 6*I))'),
    ('1', '(x**2+I)**2', 'ZZ_I', '(-I*x/(2*x**2 + 2*I), -I/(2*x**2 + 2*I))'),
    ('3*x', '(x-2*I)**2*(x**2+1)', 'ZZ_I', '(2*I/(x - 2*I), (2*I*x - 1)/(x**3 - 2*I*x**2 + x - 2*I))'),
    ('1', '(x**2+sqrt(2))**2', 'QQ<sqrt(2)>', '(sqrt(2)*x/(4*x**2 + 4*sqrt(2)), sqrt(2)/(4*x**2 + 4*sqrt(2)))'),
    ('x', '(x-sqrt(2))**3', 'QQ<sqrt(2)>', '((-2*x + sqrt(2))/(2*x**2 - 4*sqrt(2)*x + 4), 0)'),
    ('1', '(x**2+a)**2', 'QQ(a)', '(x/(2*a**2 + 2*a*x**2), 1/(2*a**2 + 2*a*x**2))'),
    ('1', '(x**2+a)**2', 'EX', '(x/(2*a**2 + 2*a*x**2), 1/(2*a**2 + 2*a*x**2))'),
    ('1', '(x**2+sqrt(2))**2', 'EX', '(sqrt(2)*x/(4*x**2 + 4*sqrt(2)), sqrt(2)/(4*x**2 + 4*sqrt(2)))'),
    ('1', '(x**2+sqrt(3)*x+1)**3', 'EX',
     '((12*x**3 + 18*sqrt(3)*x**2 + 32*x + 7*sqrt(3))/(2*x**4 + 4*sqrt(3)*x**3 + 10*x**2 + 4*sqrt(3)*x + 2), '
     '6/(x**2 + sqrt(3)*x + 1))'),
    ('a', '(x-a)**4', 'EX', '(-a/(-3*a**3 + 9*a**2*x - 9*a*x**2 + 3*x**3), 0)'),
    ('1', '(x-sqrt(2))**2', 'EX', '(-1/(x - sqrt(2)), 0)'),
    ('1', 'x**2+a', 'EX', '(0, 1/(a + x**2))'),
]

# Parametric cases whose values agree with SymPy 1.14 but whose printed
# cancel form is the one of the stock symox cancel (pre-existing gap).
RATPART_VALUE = [
    ('x', '(x-a)**2*(x-b)', 'QQ(a,b)'),
    ('1', '(x**2-2*a*x+1)**2', 'QQ(a)'),
    ('1', '(x-a)**2*(x-b)', 'EX'),
    ('x+1', '(x**2+a*x+b)**2', 'EX'),
]


def certify(f, g, value):
    rational, logarithmic = value
    assert sx.cancel(sx.diff(rational, x) + logarithmic - f.as_expr()/g.as_expr()) == 0


@pytest.mark.parametrize('p,q,domain,expected', RATPART)
def test_ratpart_matches_sympy(p, q, domain, expected, monkeypatch):
    f, g = poly(p, domain), poly(q, domain)
    value = rt.ratint_ratpart(f, g, x)
    assert str(value) == expected
    certify(f, g, value)
    no_native(monkeypatch)
    if domain.startswith('QQ<'):
        # The reference body cannot unify QQ<a> with its ZZ[a_i, b_j]
        # unknowns; the native solve returns the SymPy 1.14 value.
        with pytest.raises(TypeError):
            rt.ratint_ratpart(f, g, x)
        return
    assert sx.srepr(rt.ratint_ratpart(f, g, x)) == sx.srepr(value)


@pytest.mark.parametrize('p,q,domain', RATPART_VALUE)
def test_ratpart_parametric_matches_reference_body(p, q, domain, monkeypatch):
    f, g = poly(p, domain), poly(q, domain)
    value = rt.ratint_ratpart(f, g, x)
    certify(f, g, value)
    no_native(monkeypatch)
    assert sx.srepr(rt.ratint_ratpart(f, g, x)) == sx.srepr(value)


def test_ratpart_declines_outside_admission(monkeypatch):
    # Float coefficients and multivariate generators stay on the reference body.
    f, g = sx.Poly(1, x, domain='RR'), sx.Poly((x**2 + 1.5)**2, x)
    assert core._ratint_ratpart(f, g, x) is None
    y = sx.Symbol('y')
    assert core._ratint_ratpart(sx.Poly(1, x, y), sx.Poly((x**2 + y)**2, x, y), x) is None


# (numerator, denominator, domain, SymPy 1.14 [(h, h.domain, q, q.domain)])
LOGPART = [
    ('1', 'x**2+I', 'QQ_I',
     [("Poly(x + 2*_t*I, x, domain='ZZ_I[_t]')", 'ZZ_I[_t]', "Poly(4*I*_t**2 + 1, _t, domain='ZZ_I')", 'ZZ_I')]),
    ('1', 'x**2+I', 'ZZ_I',
     [("Poly(x + 2*_t*I, x, domain='ZZ_I[_t]')", 'ZZ_I[_t]', "Poly(4*I*_t**2 + 1, _t, domain='ZZ_I')", 'ZZ_I')]),
    ('x', 'x**3-I', 'ZZ_I',
     [("Poly(x - 9*_t**2*I, x, domain='ZZ_I[_t]')", 'ZZ_I[_t]', "Poly(27*_t**3 + I, _t, domain='ZZ_I')", 'ZZ_I')]),
    ('1+I', 'x**2+1', 'ZZ_I',
     [("Poly(x + _t*(1 - I), x, domain='ZZ_I[_t]')", 'ZZ_I[_t]', "Poly(2*_t**2 + I, _t, domain='ZZ_I')", 'ZZ_I')]),
    ('1', '(x-a)*(x-b)', 'ZZ[a,b]',
     [("Poly(x - _t*a**2/2 + _t*a*b - _t*b**2/2 - a/2 - b/2, x, domain='QQ[a,b,_t]')", 'QQ[a,b,_t]',
       "Poly(_t**2 - 1/(a**2 - 2*a*b + b**2), _t, domain='EX')", 'EX')]),
    ('1', 'x**2+a', 'ZZ[a]',
     [("Poly(x + 2*_t*a, x, domain='ZZ[a,_t]')", 'ZZ[a,_t]', "Poly(_t**2 + 1/(4*a), _t, domain='EX')", 'EX')]),
    ('x', 'x**3+a', 'ZZ[a]',
     [("Poly(x + 9*_t**2*a, x, domain='ZZ[a,_t]')", 'ZZ[a,_t]', "Poly(_t**3 + 1/(27*a), _t, domain='EX')", 'EX')]),
    ('1', 'x**2+a', 'QQ(a)',
     [("Poly(x + 2*_t*a, x, domain='ZZ[a,_t]')", 'ZZ[a,_t]', "Poly(_t**2 + 1/(4*a), _t, domain='EX')", 'EX')]),
    ('a', 'x**2-a*x+1', 'QQ(a)',
     [("Poly(x + (-_t*a**2 + 4*_t - a**2)/(2*a), x, domain='ZZ(a,_t)')", 'ZZ(a,_t)',
       "Poly(_t**2 - a**2/(a**2 - 4), _t, domain='EX')", 'EX')]),
    ('1', 'x**2+a', 'EX',
     [("Poly(x + 2*_t*a, x, domain='ZZ[a,_t]')", 'ZZ[a,_t]', "Poly(_t**2 + 1/(4*a), _t, domain='EX')", 'EX')]),
    ('1', 'x**2+sqrt(2)', 'EX',
     [("Poly(x + 2*sqrt(2)*_t, x, domain='EX')", 'EX', "Poly(_t**2 + sqrt(2)/8, _t, domain='EX')", 'EX')]),
    ('1', '(x-a)*(x-b)', 'EX',
     [("Poly(x - _t*a**2/2 + _t*a*b - _t*b**2/2 - a/2 - b/2, x, domain='QQ[a,b,_t]')", 'QQ[a,b,_t]',
       "Poly(_t**2 - 1/(a**2 - 2*a*b + b**2), _t, domain='EX')", 'EX')]),
    ('x', 'x**2+sqrt(3)*x+1', 'EX',
     [("Poly(x - sqrt(3)*_t/3 + 2*sqrt(3)/3, x, domain='EX')", 'EX', "Poly(_t**2 - _t + 1, _t, domain='ZZ')", 'ZZ')]),
]


def shown(value):
    return [(str(h), str(h.domain), str(q), str(q.domain)) for h, q in value]


def log_certificate(f, g, value):
    # sum over the roots c of q of c*h_x(c, x)/h(c, x), computed as a trace
    # in K(x)[t]/(q): it must equal f/g.
    t = value[0][1].gens[0]
    total = 0
    for h, q in value:
        h_expr, q_poly = h.as_expr(), sx.Poly(q.as_expr(), t, domain='EX')
        inverse = sx.invert(sx.Poly(h_expr, t, domain='EX'), q_poly)
        r = (sx.Poly(t*sx.diff(h_expr, x), t, domain='EX')*inverse).rem(q_poly)
        n = q_poly.degree()
        trace = sum((r*sx.Poly(t**j, t, domain='EX')).rem(q_poly).nth(j) for j in range(n))
        total += trace
    assert sx.cancel(sx.expand(total) - f.as_expr()/g.as_expr()) == 0


@pytest.mark.parametrize('p,q,domain,expected', LOGPART)
def test_logpart_matches_sympy(p, q, domain, expected, monkeypatch):
    f, g = poly(p, domain), poly(q, domain)
    t = sx.Dummy('t')
    value = rt.ratint_logpart(f, g, x, t)
    assert shown(value) == expected
    no_native(monkeypatch)
    assert shown(rt.ratint_logpart(f, g, x, t)) == expected


GAUSSIAN_THREE_POLES = [
    ("Poly(x + _t**2*(2858/1105 - 94*I/1105) + _t*(-2097/1105 - 279*I/1105) - 1141/1105 + 1233*I/1105,"
     " x, domain='QQ_I[_t]')", 'QQ_I[_t]', "Poly((-8 - 6*I)*_t**3 + (6 - I)*_t + 3, _t, domain='ZZ_I')", 'ZZ_I')]


def test_logpart_gaussian_three_poles():
    f, g = poly('3', 'QQ_I'), poly('(x-I)*(x+2*I)*(x-1)', 'QQ_I')
    assert shown(rt.ratint_logpart(f, g, x, sx.Dummy('t'))) == GAUSSIAN_THREE_POLES


def test_logpart_gaussian_three_poles_reference_body(monkeypatch):
    no_native(monkeypatch)
    f, g = poly('3', 'QQ_I'), poly('(x-I)*(x+2*I)*(x-1)', 'QQ_I')
    assert shown(rt.ratint_logpart(f, g, x, sx.Dummy('t'))) == GAUSSIAN_THREE_POLES


@pytest.mark.parametrize('name', ['a', 'b', 'h', 'i', 'l', 'm', 'o', 'p', 't', 'w', 'y'])
def test_unexpanded_gaussian_coefficient_selects_gaussian_domain(name):
    # SymPy 1.14 reads 12 - 2*I as one Gaussian coefficient for every
    # generator. The I detection used to depend on allocation addresses, so
    # some names (and Dummy indices) gave ZZ[I] and, with field=True, ZZ(I).
    t = sx.Symbol(name)
    for s in (t, sx.Dummy(name)):
        e = s**2*(12 - 2*sx.I) + 9*s
        assert str(sx.Poly(e, s).domain) == 'ZZ_I'
        assert str(sx.Poly(e, s, field=True).domain) == 'QQ_I'
        assert str(sx.Poly(s**2*(12 - 2*sx.I) + 9*s**3, s).domain) == 'ZZ_I'


@pytest.mark.parametrize('native', [True, False])
def test_logpart_gaussian_three_poles_after_optics_import(native):
    # A fresh process, so test order cannot hide a state change made by the
    # import: importing symox.physics.optics used to make this case raise
    # "TypeError: expected a polynomial" in the reference body.
    import os, subprocess, sys
    source = r'''
import importlib
import symox as sx
import symox.physics.optics
rt = importlib.import_module('symox.integrals.rationaltools')
core = importlib.import_module('symox._core')
if not NATIVE:
    core._ratint_ratpart = lambda *args: None
    core._ratint_logpart = lambda *args: None
x = sx.Symbol('x')
f = sx.Poly(3, x, domain='QQ_I')
g = sx.Poly(sx.expand((x - sx.I)*(x + 2*sx.I)*(x - 1)), x, domain='QQ_I')
print(repr([(str(h), str(h.domain), str(q), str(q.domain)) for h, q in rt.ratint_logpart(f, g, x, sx.Dummy('t'))]))
'''.replace('NATIVE', repr(native))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(
        [os.path.dirname(os.path.dirname(sx.__file__))] + [p for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p]))
    done = subprocess.run([sys.executable, '-c', source], check=True, timeout=120, capture_output=True, text=True, env=env)
    assert done.stdout.strip() == repr(GAUSSIAN_THREE_POLES)


@pytest.mark.parametrize('p,q,domain', [('1', 'x**2+I', 'QQ_I'), ('1', 'x**2+a', 'ZZ[a]'), ('a', 'x**2-a*x+1', 'QQ(a)')])
def test_logpart_trace_certificate(p, q, domain):
    f, g = poly(p, domain), poly(q, domain)
    log_certificate(f, g, rt.ratint_logpart(f, g, x, sx.Dummy('t')))


@pytest.mark.parametrize('integrand,expected', [
    ('1/(x**2+I)**2', '-I*x/(2*x**2 + 2*I) + sqrt(-I)*log(x - sqrt(-I))/4 - sqrt(-I)*log(x + sqrt(-I))/4'),
    ('1/(x**2+a)**2', 'x/(2*a**2 + 2*a*x**2) - sqrt(-1/a**3)*log(-a**2*sqrt(-1/a**3) + x)/4'
                      ' + sqrt(-1/a**3)*log(a**2*sqrt(-1/a**3) + x)/4'),
    ('1/(x**2+sqrt(2))**2', 'sqrt(2)*x/(4*x**2 + 4*sqrt(2)) + 2**(1/4)*atan(2**(3/4)*x/2)/4'),
    ('1/(x**2+2*I*x-5)', 'log(x - 2 + I)/4 - log(x + 2 + I)/4'),
])
def test_public_integrate_matches_sympy(integrand, expected):
    f = sx.sympify(integrand, locals=NAMES)
    value = sx.integrate(f, x)
    assert str(value) == expected
