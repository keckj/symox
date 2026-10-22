"""``Integral.transform``, ``Integral.as_sum`` and
``Integral.principal_value``, ported from SymPy 1.14
``sympy/integrals/integrals.py``."""

from symox import (S, Dummy, Expr, Symbol, Sum, Integral, sympify, limit,
                   sign, integrate, oo)
from symox.core.logic import fuzzy_bool
from symox.utilities.misc import filldedent


def _transform(self, x, u):
    d = Dummy('d')

    xfree = x.free_symbols.intersection(self.variables)
    if len(xfree) > 1:
        raise ValueError(
            'F(x) can only contain one of: %s' % self.variables)
    xvar = xfree.pop() if xfree else d

    if xvar not in self.variables:
        return self

    u = sympify(u)
    if isinstance(u, Expr):
        ufree = u.free_symbols
        if len(ufree) == 0:
            raise ValueError(filldedent('''
            f(u) cannot be a constant'''))
        if len(ufree) > 1:
            raise ValueError(filldedent('''
            When f(u) has more than one free symbol, the one replacing x
            must be identified: pass f(u) as (f(u), u)'''))
        uvar = ufree.pop()
    else:
        u, uvar = u
        if uvar not in u.free_symbols:
            raise ValueError(filldedent('''
            Expecting a tuple (expr, symbol) where symbol identified
            a free symbol in expr, but symbol is not in expr's free
            symbols.'''))
        if not isinstance(uvar, Symbol):
            raise ValueError(filldedent('''
            Expecting a tuple (expr, symbol) but didn't get
            a symbol; got %s''' % uvar))

    if x.is_Symbol and u.is_Symbol:
        return self.xreplace({x: u})

    if not x.is_Symbol and not u.is_Symbol:
        raise ValueError('either x or u must be a symbol')

    if uvar == xvar:
        return _transform(self, x, (u.subs(uvar, d), d)).xreplace({d: uvar})

    if uvar in self.limits:
        raise ValueError(filldedent('''
        u must contain the same variable as in x
        or a variable that is not already an integration variable'''))

    from symox.solvers.solvers import solve
    if not x.is_Symbol:
        F = [x.subs(xvar, d)]
        soln = solve(u - x, xvar, check=False)
        if not soln:
            raise ValueError('no solution for solve(F(x) - f(u), x)')
        f = [fi.subs(uvar, d) for fi in soln]
    else:
        f = [u.subs(uvar, d)]
        from symox.simplify.simplify import posify
        pdiff, reps = posify(u - x)
        puvar = uvar.subs([(v, k) for k, v in reps.items()])
        soln = [s.subs(reps) for s in solve(pdiff, puvar)]
        if not soln:
            raise ValueError('no solution for solve(F(x) - f(u), u)')
        F = [fi.subs(xvar, d) for fi in soln]

    newfuncs = {(self.function.subs(xvar, fi)*fi.diff(d)
                 ).subs(d, uvar) for fi in f}
    if len(newfuncs) > 1:
        raise ValueError(filldedent('''
        The mapping between F(x) and f(u) did not give
        a unique integrand.'''))
    newfunc = newfuncs.pop()

    def _calc_limit_1(F, a, b):
        wok = F.subs(d, a)
        if wok is S.NaN or wok.is_finite is False and a.is_finite:
            return limit(sign(b)*F, d, a)
        return wok

    def _calc_limit(a, b):
        avals = list({_calc_limit_1(Fi, a, b) for Fi in F})
        if len(avals) > 1:
            raise ValueError(filldedent('''
            The mapping between F(x) and f(u) did not
            give a unique limit.'''))
        return avals[0]

    newlimits = []
    for xab in self.limits:
        sym = xab[0]
        if sym == xvar:
            if len(xab) == 3:
                a, b = xab[1:]
                a, b = _calc_limit(a, b), _calc_limit(b, a)
                if fuzzy_bool(a - b > 0):
                    a, b = b, a
                    newfunc = -newfunc
                newlimits.append((uvar, a, b))
            elif len(xab) == 2:
                a = _calc_limit(xab[1], 1)
                newlimits.append((uvar, a))
            else:
                newlimits.append(uvar)
        else:
            newlimits.append(xab)

    return self.func(newfunc, *newlimits)


def _as_sum(self, n=None, method="midpoint", evaluate=True):
    limits = self.limits
    if len(limits) > 1:
        raise NotImplementedError(
            "Multidimensional midpoint rule not implemented yet")
    else:
        lim = limits[0]
        if (len(lim) != 3 or lim[1].is_finite is False or
                lim[2].is_finite is False):
            raise ValueError("Expecting a definite integral over "
                             "a finite interval.")
    if n is None:
        n = Dummy('n', integer=True, positive=True)
    else:
        n = sympify(n)
    if (n.is_positive is False or n.is_integer is False or
            n.is_finite is False):
        raise ValueError("n must be a positive integer, got %s" % n)
    x, a, b = lim
    dx = (b - a)/n
    k = Dummy('k', integer=True, positive=True)
    f = self.function

    if method == "left":
        result = dx*Sum(f.subs(x, a + (k-1)*dx), (k, 1, n))
    elif method == "right":
        result = dx*Sum(f.subs(x, a + k*dx), (k, 1, n))
    elif method == "midpoint":
        result = dx*Sum(f.subs(x, a + k*dx - dx/2), (k, 1, n))
    elif method == "trapezoid":
        result = dx*((f.subs(x, a) + f.subs(x, b))/2 +
                     Sum(f.subs(x, a + k*dx), (k, 1, n - 1)))
    else:
        raise ValueError("Unknown method %s" % method)
    return result.doit() if evaluate else result


def _principal_value(self, **kwargs):
    if len(self.limits) != 1 or len(list(self.limits[0])) != 3:
        raise ValueError("You need to insert a variable, lower_limit, and upper_limit correctly to calculate "
                         "cauchy's principal value")
    x, a, b = self.limits[0]
    if not (a.is_comparable and b.is_comparable and a <= b):
        raise ValueError("The lower_limit must be smaller than or equal to the upper_limit to calculate "
                         "cauchy's principal value. Also, a and b need to be comparable.")
    if a == b:
        return S.Zero

    from symox.calculus.singularities import singularities

    r = Dummy('r')
    f = self.function
    singularities_list = [s for s in singularities(f, x) if s.is_comparable and a <= s <= b]
    for i in singularities_list:
        if i in (a, b):
            raise ValueError(
                'The principal value is not defined in the given interval due to singularity at %d.' % (i))
    F = integrate(f, x, **kwargs)
    if F.has(Integral):
        return self
    if a is -oo and b is oo:
        I = limit(F - F.subs(x, -x), x, oo)
    else:
        I = limit(F, x, b, '-') - limit(F, x, a, '+')
    for s in singularities_list:
        I += limit(((F.subs(x, s - r)) - F.subs(x, s + r)), r, 0, '+')
    return I


def install():
    from symox import _node_types
    cls = _node_types.node_class('Integral')
    cls.transform = _transform
    cls.as_sum = _as_sum
    cls.principal_value = _principal_value
