"""Reference ``_tsolve`` steps that the native solver does not take.

These follow SymPy 1.14's ``sympy.solvers.solvers._tsolve`` for a single
equation in one symbol, after the native solver raised
NotImplementedError:

* ``logcombine(lhs, force=True)`` when it reduces the number of logs
  (``log(y + 1)/3 + log(y**2 - y + 1)/3 = c``);
* ``f(x)**e = rhs`` with ``e`` free of the symbol and not an integer, which
  becomes ``f(x) = rhs**(1/e)`` (``x**y = 2`` gives ``2**(1/y)``).

Returns None when no step applies, so that the caller re-raises.
"""


def _ispow(e):
    from .. import Expr, exp
    return isinstance(e, Expr) and (e.is_Pow or isinstance(e, exp))


_INVERSES = (
    ('acos', 'cos'), ('acosh', 'cosh'), ('acot', 'cot'), ('acoth', 'coth'),
    ('acsc', 'csc'), ('acsch', 'csch'), ('asec', 'sec'), ('asech', 'sech'),
    ('asin', 'sin'), ('asinh', 'sinh'), ('atan', 'tan'), ('atanh', 'tanh'),
    ('conjugate', 'conjugate'), ('cot', 'acot'), ('coth', 'acoth'),
    ('erf', 'erfinv'), ('erfc', 'erfcinv'), ('erfcinv', 'erfc'),
    ('erfinv', 'erf'), ('exp', 'log'), ('exp_polar', 'log'), ('log', 'exp'),
    ('sinh', 'asinh'), ('tan', 'atan'), ('tanh', 'atanh'))


def _inverse(f):
    """The reference ``f.inverse()`` of the one-argument functions that
    define it (``sin``, ``cos``, ``sec``, ``csc``, ``cosh`` and ``LambertW``
    define none)."""
    import symox
    inverse = getattr(f, 'inverse', None)
    if inverse is not None:
        try:
            return inverse()
        except (TypeError, NotImplementedError):
            return None
    for name, inv in _INVERSES:
        F = getattr(symox, name, None)
        if F is not None and isinstance(F, type) and isinstance(f, F):
            return getattr(symox, inv)
    return None


def _invert(eq, sym):
    """SymPy's ``_invert`` for one symbol (``integer_power=False``).
    Returns ``(rhs, lhs)``."""
    from .. import (S, Add, Function, atan, atan2, ordered, powdenest,
                    powsimp, sqrt)
    lhs = eq
    rhs = S.Zero
    while True:
        was = lhs
        while True:
            indep, dep = lhs.as_independent(sym)
            if lhs.is_Add:
                if indep.is_zero:
                    break
                lhs = dep
                rhs -= indep
            else:
                if indep is S.One:
                    break
                lhs = dep
                rhs /= indep
        if lhs.is_Add:
            terms = {}
            for a in lhs.args:
                i, d = a.as_independent(sym)
                terms.setdefault(d, []).append(i)
            if any(len(v) > 1 for v in terms.values()):
                args = []
                for d, i in terms.items():
                    if len(i) > 1:
                        args.append(Add(*i)*d)
                    else:
                        args.append(i[0]*d)
                lhs = Add(*args)
        if lhs.is_Add and not rhs and len(lhs.args) == 2 and \
                not lhs.is_polynomial(sym):
            a, b = ordered(lhs.args)
            ai, ad = a.as_independent(sym)
            bi, bd = b.as_independent(sym)
            if any(_ispow(i) for i in (ad, bd)):
                a_base, a_exp = ad.as_base_exp()
                b_base, b_exp = bd.as_base_exp()
                if a_base == b_base and a_exp.extract_additively(b_exp) is None:
                    lhs = powsimp(powdenest(ad/bd))
                    rhs = -bi/ai
                else:
                    rat = ad/bd
                    _lhs = powsimp(ad/bd)
                    if _lhs != rat:
                        lhs = _lhs
                        rhs = -bi/ai
            elif ai == -bi:
                if isinstance(ad, Function) and ad.func == bd.func:
                    if len(ad.args) == len(bd.args) == 1:
                        lhs = ad.args[0] - bd.args[0]
                    elif len(ad.args) == len(bd.args):
                        raise NotImplementedError(
                            'equal function with more than 1 argument')
                    else:
                        raise ValueError(
                            'function with different numbers of args')
        elif lhs.is_Mul and any(_ispow(a) for a in lhs.args):
            lhs = powsimp(powdenest(lhs))
        if lhs.is_Function:
            inverse = None
            if len(lhs.args) == 1:
                inverse = _inverse(lhs)
            if inverse is not None:
                rhs = inverse(rhs)
                lhs = lhs.args[0]
            elif isinstance(lhs, atan2):
                y, x = lhs.args
                lhs = 2*atan(y/(sqrt(x**2 + y**2) + x))
            elif lhs.func == rhs.func:
                if len(lhs.args) == len(rhs.args) == 1:
                    lhs = lhs.args[0]
                    rhs = rhs.args[0]
                elif len(lhs.args) == len(rhs.args):
                    raise NotImplementedError(
                        'equal function with more than 1 argument')
                else:
                    raise ValueError(
                        'function with different numbers of args')
        if rhs and lhs.is_Pow and lhs.exp.is_Integer and lhs.exp < 0:
            lhs = 1/lhs
            rhs = 1/rhs
        if lhs == was:
            break
    return rhs, lhs


_SAW = []
_NO_FORCE = []


def tsolve_step(solve, f, symbols, flags):
    """Like the reference ``tsolve_saw`` flag, an equation already being
    handled here is not handled again (that would recurse forever)."""
    key = (f, tuple(symbols))
    if key in _SAW:
        return None
    _SAW.append(key)
    try:
        return _tsolve_step(solve, f, symbols, flags)
    finally:
        _SAW.pop()


def _tsolve_step(solve, f, symbols, flags):
    from .. import (S, Equality, Expr, exp, factor, log, logcombine, posify,
                    powdenest, sympify)
    if isinstance(f, (list, tuple, set)) or len(symbols) != 1:
        return None
    if flags.get('dict') or flags.get('set'):
        return None
    sym = sympify(symbols[0])
    if not getattr(sym, 'is_Symbol', False):
        return None
    f = sympify(f)
    if isinstance(f, Equality):
        f = f.lhs - f.rhs
    if not isinstance(f, Expr) or not f.has(sym):
        return None
    # the reference ``_solve`` hands the numerator to ``_tsolve``
    num, den = f.as_numer_denom()
    if den.has(sym):
        f = num
    try:
        rhs, lhs = _invert(f, sym)
    except (NotImplementedError, ValueError):
        return None
    if lhs == sym:
        return [rhs]
    if lhs.is_Add:
        # factoring (after powdenest) may split the equation
        g = factor(powdenest(lhs - rhs))
        if g.is_Mul:
            indep, dep = g.as_independent(sym, as_Add=False)
            if indep is not S.One and indep != 0 and dep.has(sym):
                # a constant factor: the reference ``_solve`` drops it
                return solve(dep, sym, **flags)
            if g.is_Mul and dep.is_Mul:
                return solve(g, sym, **flags)
    if lhs.is_Add and rhs != 0 and lhs.has(log):
        g = logcombine(lhs, force=flags.get('force', True) and not _NO_FORCE)
        if g.count(log) != lhs.count(log):
            if isinstance(g, log):
                return solve(g.args[0] - exp(rhs), sym, **flags)
            return solve(g - rhs, sym, **flags)
    if lhs.is_Pow and not lhs.exp.is_Integer and \
            sym not in lhs.exp.free_symbols:
        return solve(lhs.base - rhs**(1/lhs.exp), sym, **flags)
    # the reference ends with a posify retry (``force`` is then False)
    if flags.get('force', True) and not _NO_FORCE:
        pos, reps = posify(lhs - rhs)
        for u, s in reps.items():
            if s == sym:
                break
        else:
            u = sym
        if pos.has(u):
            # ``force=False`` for the inner solve (the native solver does
            # not read that flag, so it is kept here)
            _NO_FORCE.append(True)
            try:
                soln = solve(pos, u, **flags)
            except NotImplementedError:
                return None
            finally:
                _NO_FORCE.pop()
            return [s.subs(reps) for s in soln]
    return None
