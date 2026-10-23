"""Port of SymPy 1.14's ``sympy.solvers.bivariate``.

``_solve_lambert`` solves equations that become ``x = W(...)`` forms after
taking logarithms (``x*exp(x**2) = 1``, ``log(x) + x**2 = 0``,
``x**a*exp(x) = 1``) and ``bivariate_type`` recognizes ``f(x*y)``,
``f(a*x + b*y)`` and ``f(a*x*y + b*y)``. They are used by the ``_tsolve``
fallback of ``solve`` (``python/symox/solvers/_tsolve_steps.py``) after the
native solver raised NotImplementedError.
"""


def _mexpand(expr, recursive=False):
    from .. import expand_mul, expand_multinomial
    if expr is None:
        return None
    was = None
    while was != expr:
        was, expr = expr, expand_mul(expand_multinomial(expr))
        if not recursive:
            break
    return expr


def _uniq(seq):
    out = []
    for s in seq:
        if s not in out:
            out.append(s)
    return out


def _filtered_gens(poly, symbol):
    """Generators of ``poly`` that contain ``symbol``; of ``g`` and ``1/g``
    only one is kept."""
    from .. import S
    gens = {g for g in poly.gens if symbol in g.free_symbols}
    for g in list(gens):
        ag = 1/g
        if g in gens and ag in gens:
            if ag.as_numer_denom()[1] is not S.One:
                g = ag
            gens.remove(g)
    return gens


def _mostfunc(lhs, func, X=None):
    """The term of ``lhs`` with ``func`` that holds the most ``func``."""
    from .. import ordered
    fterms = [tmp for tmp in lhs.atoms(func) if (not X or
              X.is_Symbol and X in tmp.free_symbols or
              not X.is_Symbol and tmp.has(X))]
    if len(fterms) == 1:
        return fterms[0]
    elif fterms:
        return max(list(ordered(fterms)), key=lambda x: x.count(func))
    return None


def _linab(arg, symbol):
    """``(a, b, x)`` with ``arg = a*x + b``, ``x`` the symbol-dependent part."""
    from .. import factor_terms, separatevars
    arg = factor_terms(arg.expand())
    ind, dep = arg.as_independent(symbol)
    if arg.is_Mul and dep.is_Add:
        a, b, x = _linab(dep, symbol)
        return ind*a, ind*b, x
    if not arg.is_Add:
        b = 0
        a, x = ind, dep
    else:
        b = ind
        a, x = separatevars(dep).as_independent(symbol, as_Add=False)
    if x.could_extract_minus_sign():
        a = -a
        x = -x
    return a, b, x


def _lambert(eq, x):
    """Solutions of ``a*log(b*X + c) + d*X + f = 0`` with ``X = g(x)``."""
    from .. import Dummy, LambertW, collect, exp, expand_log, log, roots, solve
    eq = _mexpand(expand_log(eq))
    mainlog = _mostfunc(eq, log, x)
    if not mainlog:
        return []  # violated assumptions
    other = eq.subs(mainlog, 0)
    if isinstance(-other, log):
        eq = (eq - other).subs(mainlog, mainlog.args[0])
        mainlog = mainlog.args[0]
        if not isinstance(mainlog, log):
            return []  # violated assumptions
        other = -(-other).args[0]
        eq += other
    if x not in other.free_symbols:
        return []  # violated assumptions
    d, f, X2 = _linab(other, x)
    logterm = collect(eq - other, mainlog)
    a = logterm.as_coefficient(mainlog)
    if a is None or x in a.free_symbols:
        return []  # violated assumptions
    logarg = mainlog.args[0]
    b, c, X1 = _linab(logarg, x)
    if X1 != X2:
        return []  # violated assumptions

    # invert the generator X1 so we have x(u)
    u = Dummy('rhs')
    xusolns = solve(X1 - u, x)

    # only the branches k = -1 and 0 may be real; k = -1 is kept when it
    # tests as real
    lambert_real_branches = [-1, 0]
    sol = []
    num, den = ((c*d - b*f)/a/b).as_numer_denom()
    p, den = den.as_coeff_Mul()
    e = exp(num/den)
    t = Dummy('t')
    args = [d/(a*b)*t for t in roots(t**p - e, t).keys()]
    for arg in args:
        for k in lambert_real_branches:
            w = LambertW(arg, k)
            if k and not w.is_real:
                continue
            rhs = -c/b + (a/d)*w
            sol.extend(xu.subs(u, rhs) for xu in xusolns)
    return sol


def _solve_lambert(f, symbol, gens):
    """SymPy's ``_solve_lambert``: raises NotImplementedError when ``f``
    has no Lambert form."""
    from .. import (S, Dummy, Pow, collect, exp, expand_log, factor, log,
                    ordered, powsimp)
    from ._tsolve_steps import _invert

    def _solve_even_degree_expr(expr, t, symbol):
        nlhs, plhs = [
            expr.xreplace({t: sgn*symbol}) for sgn in (-1, 1)]
        sols = _solve_lambert(nlhs, symbol, gens)
        if plhs != nlhs:
            sols.extend(_solve_lambert(plhs, symbol, gens))
        return _uniq(sols)

    nrhs, lhs = f.as_independent(symbol, as_Add=True)
    rhs = -nrhs

    lamcheck = [tmp for tmp in gens
                if (isinstance(tmp, (exp, log)) or
                    (tmp.is_Pow and symbol in tmp.exp.free_symbols))]
    if not lamcheck:
        raise NotImplementedError()

    if lhs.is_Add or lhs.is_Mul:
        # even powers of the symbol are replaced by powers of t
        t = Dummy('t', **symbol.assumptions0)
        lhs = lhs.replace(
            lambda i: i.is_Pow and i.base == symbol and i.exp.is_even,
            lambda i: t**i.exp)
        if lhs.is_Add and lhs.has(t):
            t_indep = lhs.subs(t, 0)
            t_term = lhs - t_indep
            _rhs = rhs - t_indep
            if not t_term.is_Add and _rhs and not (
                    t_term.has(S.ComplexInfinity, S.NaN)):
                eq = expand_log(log(t_term) - log(_rhs))
                return _solve_even_degree_expr(eq, t, symbol)
        elif lhs.is_Mul and rhs:
            lhs = expand_log(log(lhs), force=True)
            rhs = log(rhs)
            if lhs.has(t) and lhs.is_Add:
                eq = lhs - rhs
                return _solve_even_degree_expr(eq, t, symbol)
        lhs = lhs.xreplace({t: symbol})

    lhs = powsimp(factor(lhs, deep=True))

    # make sure we have inverted as completely as possible
    r = Dummy()
    i, lhs = _invert(lhs - r, symbol)
    rhs = i.xreplace({r: rhs})

    soln = []
    mainlog = _mostfunc(lhs, log, symbol)
    if mainlog:
        if lhs.is_Mul and rhs != 0:
            soln = _lambert(log(lhs) - log(rhs), symbol)
        elif lhs.is_Add:
            other = lhs.subs(mainlog, 0)
            if other and not other.is_Add and [
                    tmp for tmp in other.atoms(Pow)
                    if symbol in tmp.free_symbols]:
                if not rhs:
                    diff = log(other) - log(other - lhs)
                else:
                    diff = log(lhs - other) - log(rhs - other)
                soln = _lambert(expand_log(diff), symbol)
            else:
                soln = _lambert(lhs - rhs, symbol)

    if not soln:
        mainexp = _mostfunc(lhs, exp, symbol)
        if mainexp:
            lhs = collect(lhs, mainexp)
            if lhs.is_Mul and rhs != 0:
                soln = _lambert(expand_log(log(lhs) - log(rhs)), symbol)
            elif lhs.is_Add:
                other = lhs.subs(mainexp, 0)
                mainterm = lhs - other
                rhs = rhs - other
                if (mainterm.could_extract_minus_sign() and
                        rhs.could_extract_minus_sign()):
                    mainterm *= -1
                    rhs *= -1
                diff = log(mainterm) - log(rhs)
                soln = _lambert(expand_log(diff), symbol)

    if not soln:
        mainpow = _mostfunc(lhs, Pow, symbol)
        if mainpow and symbol in mainpow.exp.free_symbols:
            lhs = collect(lhs, mainpow)
            if lhs.is_Mul and rhs != 0:
                soln = _lambert(expand_log(log(lhs) - log(rhs)), symbol)
            elif lhs.is_Add:
                other = lhs.subs(mainpow, 0)
                mainterm = lhs - other
                rhs = rhs - other
                diff = log(mainterm) - log(rhs)
                soln = _lambert(expand_log(diff), symbol)

    if not soln:
        raise NotImplementedError('%s does not appear to have a solution in '
                                  'terms of LambertW' % f)
    return list(ordered(soln))


def bivariate_type(f, x, y, *, first=True):
    """``(g, p, u)`` such that ``f(x, y) = p(u)`` with ``u = g(x, y)`` one
    of ``x*y``, ``a*x + b*y`` or ``a*x*y + b*y``, else None."""
    from .. import Add, Dummy, Poly, root
    u = Dummy('u', positive=True)

    if first:
        p = Poly(f, x, y)
        f = p.as_expr()
        _x = Dummy()
        _y = Dummy()
        rv = bivariate_type(Poly(f.subs({x: _x, y: _y}), _x, _y), _x, _y,
                            first=False)
        if rv:
            reps = {_x: x, _y: y}
            return rv[0].xreplace(reps), rv[1].xreplace(reps), rv[2]
        return None

    p = f
    f = p.as_expr()

    # f(x*y)
    args = Add.make_args(p.as_expr())
    new = []
    for a in args:
        a = _mexpand(a.subs(x, u/y))
        free = a.free_symbols
        if x in free or y in free:
            break
        new.append(a)
    else:
        return x*y, Add(*new), u

    def ok(f, v, c):
        new = _mexpand(f.subs(v, c))
        free = new.free_symbols
        return None if (x in free or y in free) else new

    # f(a*x + b*y)
    d = p.degree(x)
    if p.degree(y) == d:
        a = root(p.coeff_monomial(x**d), d)
        b = root(p.coeff_monomial(y**d), d)
        new = ok(f, x, (u - b*y)/a)
        if new is not None:
            return a*x + b*y, new, u

    # f(a*x*y + b*y)
    d = p.degree(x)
    if p.degree(y) == d:
        for itry in range(2):
            a = root(p.coeff_monomial(x**d*y**d), d)
            b = root(p.coeff_monomial(y**d), d)
            new = ok(f, x, (u - b*y)/a/y)
            if new is not None:
                return a*x*y + b*y, new, u
            x, y = y, x
    return None
