"""Composite domains over an algebraic ground (wave 44, lot algpoly).

``QQ<a>[t]`` and ``QQ<a>(t)`` (``GaussianGround`` with an ``algebraic``
number field in ``crates/symox-core/src/poly/gaussian_ground.rs``), their
elements, the polynomials over them, the unification with other domains,
``ratint_logpart`` over them and the ``EX`` element arithmetic
(``expression_domain::simplify_element``) compute in Rust. Every case runs
live on SymPy 1.14 and on symox and compares the printed results, or the
exception classes.
"""
import pytest

import sympy
import symox

sympy_ratint_logpart = pytest.importorskip('sympy.integrals.rationaltools').ratint_logpart
from symox.integrals.rationaltools import ratint_logpart as symox_ratint_logpart


def env(M):
    x, y, t = M.symbols('x y t')
    p = M.symbols('p', positive=True)
    s2, s3 = M.sqrt(2), M.sqrt(3)
    A = M.QQ.algebraic_field(s2)
    return dict(M=M, x=x, y=y, t=t, p=p, s2=s2, s3=s3, A=A, R=A[t], F=A.frac_field(t))


def PR(e, M, R, x, **_):
    return M.Poly(e, x, domain=R)


def PF(e, M, F, x, **_):
    return M.Poly(e, x, domain=F)


CASES = {
    'R element': lambda M, R, s2, t, **_: R(s2*t + 1),
    'R gcd': lambda M, R, s2, t, **_: R.gcd(R(2*t**2 - 4), R(s2*t + 2)),
    'R gcd monomial': lambda M, R, s2, t, **_: R.gcd(R(2*s2*t**2), R(4*t)),
    'R exquo': lambda M, R, s2, t, **_: R.exquo(R(2*t**2 - 4), R(s2*t + 2)),
    'R from 1/t': lambda M, R, t, **_: R(1/t),
    'F element inverse': lambda M, F, s2, t, **_: F(1/(s2*t)),
    'F element scaled': lambda M, F, t, **_: F(t/(2*t + 2)),
    'F add': lambda M, F, s2, t, **_: F(1/t) + F(s2/(t + 1)),
    'F div': lambda M, F, s2, t, **_: F(t)/F(s2*t + 2),
    'F div by zero': lambda M, F, t, **_: F(t)/F(0),
    'F pow': lambda M, F, s2, t, **_: F(t/(s2*t + 1))**-2,
    'F gcd': lambda M, F, s2, t, **_: F.gcd(F((s2*t + 1)/t), F(1/(s2*t))),
    'F lcm': lambda M, F, s2, t, **_: F.lcm(F((s2*t + 1)/t), F(1/(s2*t))),
    'PR mul': lambda s2, t, x, **k: PR(s2*t*x**2 + (t + s2)*x + 2, x=x, **k)*PR((x + s2*t)*(x - 1), x=x, **k),
    'PR gcd': lambda s2, t, x, **k: PR((x + s2*t)**2*(x + t), x=x, **k).gcd(PR((x + s2*t)*(x - 1), x=x, **k)),
    'PR cancel signs': lambda s2, t, x, **k: PR(-s2*x - 2*t, x=x, **k).cancel(PR(-x**2 + 2*t**2, x=x, **k)),
    'PR resultant': lambda s2, t, x, **k: PR(s2*t*x**2 + (t + s2)*x + 2, x=x, **k).resultant(PR((x + s2*t)*(x - 1), x=x, **k)),
    'PR sqf_list': lambda s2, t, x, **k: PR((x + s2*t)**2*(x + t), x=x, **k).sqf_list(),
    'PR factor_list': lambda s2, t, x, **k: PR(x**2 - 2*t**2, x=x, **k).factor_list(),
    'PF monic': lambda s2, t, x, **k: PF((x + s2/t)*(t*x - 1), x=x, **k).monic(),
    'PF div': lambda s2, t, x, **k: PF(x**2/(s2*t) + x/t + 1, x=x, **k).div(PF((x + s2/t)*(t*x - 1), x=x, **k)),
    'PF gcd': lambda s2, t, x, **k: PF((x + s2/t)*(x + t), x=x, **k).gcd(PF((x + s2/t)*(x - 1), x=x, **k)),
    'unify R with A': lambda M, A, s2, t, x, **k: PR(x*t, M=M, x=x, **k)*M.Poly(s2*x, x, domain=A),
    'unify R with QQ<sqrt3>': lambda M, s3, t, x, **k: PR(x*t, M=M, x=x, **k) + M.Poly(s3*x, x, extension=True),
    'Poly extension ring': lambda M, s2, t, x, **_: M.Poly(s2*x*t + 1, x, extension=True),
    'EX product': lambda M, s2, **_: M.EX(s2*(1 + s2)),
    'EX pow': lambda M, s2, **_: M.EX(s2 + 1)**2,
    'EX div': lambda M, x, **_: M.EX(x**2 - 1)/M.EX(x + 1),
}


def run(M, name):
    try:
        return str(CASES[name](**env(M)))
    except Exception as e:  # compare the exception class only
        return type(e).__name__


@pytest.mark.parametrize('name', sorted(CASES))
def test_against_sympy(name):
    assert run(symox, name) == run(sympy, name)


def logpart(M, fn, ground, num, den):
    e = env(M)
    x = e['x']
    dom = e['A'][e['p']] if ground == 'ring' else e['F']
    return str(fn(M.Poly(num(**e), x, domain=dom), M.Poly(den(**e), x, domain=dom), x))


LOGPART = [
    ('ring', lambda **e: 1, lambda x, s2, p, **_: x**2 + s2*p),
    ('ring', lambda x, **_: x, lambda x, s2, p, **_: x**2 + s2*x + p),
    ('ring', lambda **e: 1, lambda x, p, **_: x**2 + p),
    ('ring', lambda **e: 1, lambda x, s2, p, **_: x**3 + s2*p),
    ('ring', lambda **e: 1, lambda x, s2, p, **_: (x - p)*(x + s2)),
    ('field', lambda **e: 1, lambda x, s2, t, **_: x**2 + s2*t),
    ('field', lambda x, **_: x, lambda x, s2, t, **_: x**2 + s2*x + t),
    ('field', lambda **e: 1, lambda x, s2, t, **_: x**3 + t*x + s2),
    ('field', lambda **e: 1, lambda x, s2, t, **_: (x - t)*(x + s2)),
    ('field', lambda x, **_: x + 1, lambda x, s2, t, **_: x**2 + t*x + s2),
]


@pytest.mark.parametrize('ground,num,den', LOGPART)
def test_ratint_logpart_against_sympy(ground, num, den):
    assert logpart(symox, symox_ratint_logpart, ground, num, den) == \
        logpart(sympy, sympy_ratint_logpart, ground, num, den)


def test_integer_cases_unchanged():
    def gcd(M):
        x, t = M.symbols('x t')
        return str(M.Poly(x**2*t + 2*x, x, domain=M.ZZ[t]).gcd(M.Poly(x*t, x, domain=M.ZZ[t])))
    assert gcd(symox) == gcd(sympy)
