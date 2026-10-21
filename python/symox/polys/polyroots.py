"""Closed-form roots of low-degree polynomials and root preprocessing.

``roots_linear``, ``roots_quadratic``, ``roots_binomial``,
``preprocess_roots`` and ``_integer_basis`` follow the reference's
``sympy.polys.polyroots`` step by step: they decide the order of the roots
that ``CRootOf`` returns for trivial polynomials and the integer scaling
that it applies before isolation. ``roots`` is the native solver.
"""
from .. import (I, Mul, Pow, S, cancel, divisors, exp, factor, factor_terms,
                gcd_list, pi, root, roots, sqrt)
from ..core.mul import expand_2arg
from .polyerrors import DomainError


def _simplify(expr):
    from .. import simplify
    return simplify(expr)


def roots_linear(f):
    """Returns a list of roots of a linear polynomial."""
    r = -f.nth(0)/f.nth(1)
    dom = f.get_domain()

    if not dom.is_Numerical:
        if dom.is_Composite:
            r = factor(r)
        else:
            r = _simplify(r)

    return [r]


def roots_quadratic(f):
    """Returns a list of roots of a quadratic polynomial."""
    a, b, c = f.all_coeffs()
    dom = f.get_domain()

    def _sqrt(d):
        # Squares leave the square root, as both signs are in the results.
        co = []
        other = []
        for di in Mul.make_args(d):
            if di.is_Pow and di.exp.is_Integer and di.exp % 2 == 0:
                co.append(Pow(di.base, di.exp//2))
            else:
                other.append(di)
        if co:
            d = Mul(*other)
            co = Mul(*co)
            return co*sqrt(d)
        return sqrt(d)

    def _simplify_quadratic(expr):
        if dom.is_Composite:
            return factor(expr)
        else:
            return _simplify(expr)

    if c is S.Zero:
        r0, r1 = S.Zero, -b/a

        if not dom.is_Numerical:
            r1 = _simplify_quadratic(r1)
        elif r1.is_negative:
            r0, r1 = r1, r0
    elif b is S.Zero:
        r = -c/a
        if not dom.is_Numerical:
            r = _simplify_quadratic(r)

        R = _sqrt(r)
        r0 = -R
        r1 = R
    else:
        d = b**2 - 4*a*c
        A = 2*a
        B = -b/A

        if not dom.is_Numerical:
            d = _simplify_quadratic(d)
            B = _simplify_quadratic(B)

        D = factor_terms(_sqrt(d)/A)
        r0 = B - D
        r1 = B + D
        if a.is_negative:
            r0, r1 = r1, r0
        elif not dom.is_Numerical:
            r0, r1 = [expand_2arg(i) for i in (r0, r1)]

    return [r0, r1]


def roots_binomial(f):
    """Returns a list of roots of a binomial polynomial."""
    n = f.degree()

    a, b = f.nth(n), f.nth(0)
    base = -cancel(b/a)
    alpha = root(base, n)

    if alpha.is_number:
        alpha = alpha.expand(complex=True)

    # Over ZZ these indices give the reals first, then the non-real roots
    # by real and imaginary part.
    neg = base.is_negative
    even = n % 2 == 0
    if neg:
        if even == True and (base + 1).is_positive:
            big = True
        else:
            big = False

    ks = []
    imax = n//2
    if even:
        ks.append(imax)
        imax -= 1
    if not neg:
        ks.append(0)
    for i in range(imax, 0, -1):
        if neg:
            ks.extend([i, -i])
        else:
            ks.extend([-i, i])
    if neg:
        ks.append(0)
        if big:
            # The reference computes reversed pairs and discards them.
            for i in range(0, len(ks), 2):
                pair = ks[i: i + 2]
                pair = list(reversed(pair))

    roots_, d = [], 2*I*pi/n
    for k in ks:
        zeta = exp(k*d).expand(complex=True)
        roots_.append((alpha*zeta).expand(power_base=False))

    return roots_


def _integer_basis(poly):
    """The integer ``div`` with ``p(div*y) = m*q(y)``, ``q`` having smaller
    coefficients than ``p``, or ``None`` if there is no possible scaling.
    """
    monoms, coeffs = list(zip(*poly.terms()))
    monoms, = list(zip(*monoms))
    coeffs = [abs(int(c)) for c in coeffs]

    if coeffs[0] < coeffs[-1]:
        coeffs = list(reversed(coeffs))
        n = monoms[0]
        monoms = [n - i for i in reversed(monoms)]
    else:
        return None

    monoms = monoms[:-1]
    coeffs = coeffs[:-1]

    # Special case for two-term polynomials
    if len(monoms) == 1:
        r = Pow(coeffs[0], S.One/monoms[0])
        if r.is_Integer:
            return int(r)
        else:
            return None

    divs = reversed([int(d) for d in divisors(int(gcd_list(coeffs)))][1:])

    try:
        div = next(divs)
    except StopIteration:
        return None

    while True:
        for monom, coeff in zip(monoms, coeffs):
            if coeff % div**monom != 0:
                try:
                    div = next(divs)
                except StopIteration:
                    return None
                else:
                    break
        else:
            return div


def preprocess_roots(poly):
    """Try to get rid of symbolic coefficients from ``poly``."""
    coeff = S.One

    poly_func = poly.func
    try:
        _, poly = poly.clear_denoms(convert=True)
    except DomainError:
        return coeff, poly

    poly = poly.primitive()[1]
    poly = poly.retract()

    if poly.get_domain().is_Poly and all(c.is_term for c in poly.rep.coeffs()):
        poly = poly.inject()

        strips = list(zip(*poly.monoms()))
        gens = list(poly.gens[1:])

        base, strips = strips[0], strips[1:]

        for gen, strip in zip(list(gens), strips):
            reverse = False

            if strip[0] < strip[-1]:
                strip = reversed(strip)
                reverse = True

            ratio = None

            for a, b in zip(base, strip):
                if not a and not b:
                    continue
                elif not a or not b:
                    break
                elif b % a != 0:
                    break
                else:
                    _ratio = b // a

                    if ratio is None:
                        ratio = _ratio
                    elif ratio != _ratio:
                        break
            else:
                if reverse:
                    ratio = -ratio

                poly = poly.eval(gen, 1)
                coeff *= gen**(-ratio)
                gens.remove(gen)

        if gens:
            poly = poly.eject(*gens)

    if poly.is_univariate and poly.get_domain().is_ZZ:
        basis = _integer_basis(poly)

        if basis is not None:
            n = poly.degree()

            def func(k, coeff):
                return coeff//basis**(n - k[0])

            poly = poly.termwise(func)
            coeff *= basis

    if not isinstance(poly, poly_func):
        poly = poly_func(poly)
    return coeff, poly


def root_factors(f, *gens, filter=None, **args):
    """
    Returns all factors of a univariate polynomial, as the reference
    ``root_factors``: the linear factors of the closed-form roots, in
    ``ordered`` order, and the cofactor of the roots not found.
    """
    from functools import reduce
    from ..core.sorting import ordered
    from .polytools import Poly
    args = dict(args)
    F = Poly(f, *gens, **args)
    if not F.is_Poly:
        return [f]
    if F.is_multivariate:
        raise ValueError('multivariate polynomials are not supported')
    x = F.gens[0]
    zeros = roots(F, filter=filter)
    if not zeros:
        factors = [F]
    else:
        factors, N = [], 0
        for r, n in ordered(zeros.items()):
            factors, N = factors + [Poly(x - r, x)]*n, N + n
        if N < F.degree():
            G = reduce(lambda p, q: p*q, factors)
            factors.append(F.quo(G))
    if not isinstance(f, Poly):
        factors = [f.as_expr() for f in factors]
    return factors


__all__ = ['roots', 'roots_linear', 'roots_quadratic', 'roots_binomial',
           'preprocess_roots', 'root_factors']
