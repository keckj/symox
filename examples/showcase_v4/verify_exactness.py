"""Generate independently checked inputs for the exactness scene."""
import json
from pathlib import Path

import mpmath as mp
import symox as sx
import sympy as sp


def generate():
    x = sx.Symbol('x')
    q = sx.Rational(1, 3) + sx.Rational(1, 6)
    assert q == sx.Rational(1, 2)
    # This cancellation loses the unit at machine precision. Both CAS values
    # are compared with an independent 100-digit mpmath calculation.
    e = sx.log(1 + sx.Rational(1, 10**40))
    value = str(e.evalf(50))
    reference = str(sp.log(1 + sp.Rational(1, 10**40)).evalf(50))
    with mp.workdps(110):
        expected = mp.log1p(mp.mpf(1)/10**40)
        assert abs(mp.mpf(value)-expected) < abs(expected)*mp.mpf('1e-49')
        assert abs(mp.mpf(value)-mp.mpf(reference)) < abs(expected)*mp.mpf('1e-49')
    a = sx.Symbol('a', positive=True)
    result = sx.integrate(sx.exp(-a*x), (x, 0, sx.oo))
    assert sx.simplify(result - 1/a) == 0
    # Independent Gaussian-entropy formula, even when the timed reference
    # stats implementation times out.
    from symox.stats import Normal, entropy
    sigma = sx.Symbol('sigma', positive=True)
    mu = sx.Symbol('mu')
    h = entropy(Normal('X', mu, sigma))
    assert sx.simplify(h - (1 + sx.log(2*sx.pi*sigma**2))/2) == 0
    return dict(rational=str(q), precision=dict(input='log(1 + Rational(1, 10**40))',
                digits=50, result=value, sympy=reference, independent='mpmath.log1p, 110 digits'),
                convergence=dict(assumption='a > 0', result=str(result)),
                entropy=dict(result=str(h), verified=True, method='Gaussian entropy identity'))


if __name__ == '__main__':
    out = generate()
    Path(__file__).with_name('data').joinpath('exactness.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))
