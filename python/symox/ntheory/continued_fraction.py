"""Continued fractions (``sympy.ntheory.continued_fraction``).

Periodic expansions of quadratic irrationals run natively with exact
integer arithmetic; the generic iterator, convergents and the reduction
back to a closed form follow SymPy's algorithms on symox expressions.
"""

import itertools

from .. import _core
from ._intfunc import as_int


def _sym(a):
    from .. import sympify
    if isinstance(a, str):
        raise ValueError("SympifyError: %r" % (a,))
    return sympify(a)


def continued_fraction(a):
    """Continued fraction of a rational or quadratic irrational."""
    from .. import S, Integer
    e = _sym(a)
    if all(i.is_Rational for i in e.atoms()):
        if e.is_Integer:
            return continued_fraction_periodic(e, 1, 0)
        elif e.is_Rational:
            return continued_fraction_periodic(e.p, e.q, 0)
        elif e.is_Pow and e.exp == S.Half and e.base.is_Integer:
            return continued_fraction_periodic(0, 1, e.base)
        elif e.is_Mul and len(e.args) == 2 and (
                e.args[0].is_Rational and e.args[1].is_Pow
                and e.args[1].base.is_Integer and e.args[1].exp == S.Half):
            a, b = e.args
            return continued_fraction_periodic(0, a.q, b.base, a.p)
        else:
            p, d = e.expand().as_numer_denom()
            if d.is_Integer:
                if p.is_Rational:
                    return continued_fraction_periodic(p, d)
                if p.is_Add and len(p.args) == 2:
                    a, bc = p.args
                else:
                    a = S.Zero
                    bc = p
                if a.is_Integer:
                    b = S.NaN
                    c = None
                    if bc.is_Mul and len(bc.args) == 2:
                        b, c = bc.args
                    elif bc.is_Pow:
                        b = Integer(1)
                        c = bc
                    if b.is_Integer and c is not None and (
                            c.is_Pow and c.exp == S.Half and c.base.is_Integer):
                        c = c.base
                        return continued_fraction_periodic(a, d, c, b)
    raise ValueError(
        'expecting a rational or quadratic irrational, not %s' % e)


def continued_fraction_periodic(p, q, d=0, s=1):
    """Periodic continued fraction of ``(p + s*sqrt(d))/q``.

    The repeating part, if any, is returned as a trailing list.
    """
    from .. import Integer
    p, q, d, s = list(map(as_int, [p, q, d, s]))
    terms, start = _core.nt_cf_periodic(p, q, d, s)
    terms = [Integer(t) for t in terms]
    if start is None:
        return terms
    return terms[:start] + [terms[start:]]


def continued_fraction_reduce(cf):
    """Reduce a (possibly periodic) continued fraction to an expression."""
    from .. import S, Dummy, solve
    period = []
    x = Dummy('x')

    def untillist(cf):
        for nxt in cf:
            if isinstance(nxt, list):
                period.extend(nxt)
                yield x
                break
            yield nxt

    a = S.Zero
    for a in continued_fraction_convergents(untillist(cf)):
        pass
    if period:
        y = Dummy('y')
        solns = solve(continued_fraction_reduce(period + [y]) - y, y)
        solns.sort()
        pure = solns[-1]
        rv = a.subs(x, pure).radsimp()
    else:
        rv = a
    if rv.is_Add:
        rv = _factor_terms_add(rv)
        if rv.is_Mul and rv.args[0] == -1:
            rv = rv.func(*rv.args)
    return rv


def _factor_terms_add(rv):
    """``factor_terms`` restricted to its rational-content extraction."""
    from .. import Add, Mul, S
    cont, p = rv.as_content_primitive()
    if p.is_Add and all(t.as_coeff_Mul()[0].is_negative for t in p.args):
        # SymPy returns -1*c*(...) here, which the caller re-evaluates
        # back into the distributed sum.
        return rv
    if cont == S.One:
        return p
    if p.is_Add and cont.is_Rational and cont.q != 1:
        args = [t.as_coeff_Mul() for t in p.args]
        args = [(c * cont, m) for c, m in args]
        if any(c.is_Integer for c, _ in args):
            return Add(*[m if c == 1 else Mul(c, m) for c, m in args])
    return Mul(cont, p, evaluate=False)


def continued_fraction_iterator(x):
    """Yield the terms of the simple continued fraction of x."""
    from .. import floor
    while True:
        i = floor(x)
        yield i
        x -= i
        if not x:
            break
        x = 1 / x


def continued_fraction_convergents(cf):
    """Yield the convergents of a continued fraction."""
    from .. import S
    if isinstance(cf, list) and cf and isinstance(cf[-1], list):
        cf = itertools.chain(cf[:-1], itertools.cycle(cf[-1]))
    p_2, q_2 = S.Zero, S.One
    p_1, q_1 = S.One, S.Zero
    for a in cf:
        p, q = a * p_1 + p_2, a * q_1 + q_2
        p_2, q_2 = p_1, q_1
        p_1, q_1 = p, q
        yield p / q
