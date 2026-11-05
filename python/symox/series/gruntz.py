"""The public functions of ``sympy.series.gruntz``.

Every computation is native (``symox::limits::gruntz::api``); this module only
checks the arguments, converts the ``SubsSet`` container to the pairs the
native functions take and back.
"""


class SubsSet(dict):
    """``(expr, dummy)`` pairs of an mrv set and the rewrites
    ``dummy -> exp(argument in the dummies)`` of its elements."""

    def __init__(self):
        self.rewrites = {}

    def __repr__(self):
        return super().__repr__() + ', ' + self.rewrites.__repr__()

    def __getitem__(self, key):
        if key not in self:
            from .. import Dummy
            self[key] = Dummy()
        return dict.__getitem__(self, key)

    def do_subs(self, e):
        """Substitute the variables with expressions"""
        for expr, var in self.items():
            e = e.xreplace({var: expr})
        return e

    def meets(self, s2):
        """Tell whether or not self and s2 have non-empty intersection"""
        return set(self.keys()).intersection(list(s2.keys())) != set()

    def union(self, s2, exps=None):
        """Compute the union of self and s2, adjusting exps"""
        res = self.copy()
        tr = {}
        for expr, var in s2.items():
            if expr in self:
                if exps:
                    exps = exps.xreplace({var: res[expr]})
                tr[var] = res[expr]
            else:
                res[expr] = var
        for var, rewr in s2.rewrites.items():
            res.rewrites[var] = rewr.xreplace(tr)
        return res, exps

    def copy(self):
        """Create a shallow copy of SubsSet"""
        r = SubsSet()
        r.rewrites = self.rewrites.copy()
        for expr, var in self.items():
            r[expr] = var
        return r


def _s(e):
    from ..core.sympify import sympify
    return sympify(e)


def compare(a, b, x):
    """Returns "<" if a<b, "=" for a == b, ">" for a>b"""
    from .. import _core
    return _core._gruntz_compare(_s(a), _s(b), _s(x))


def sign(e, x):
    """Returns the sign of an expression e(x) for large x: -1, 0 or 1."""
    from .. import _core
    return _core._gruntz_sign(_s(e), _s(x))


def mrv(e, x):
    """Returns a SubsSet of most rapidly varying (mrv) subexpressions of 'e',
    and e rewritten in terms of these"""
    from .. import _core
    pairs, rewrites, f = _core._gruntz_mrv(_s(e), _s(x))
    s = SubsSet()
    for expr, var in pairs:
        dict.__setitem__(s, expr, var)
    s.rewrites = dict(rewrites)
    return s, f


def rewrite(e, Omega, x, wsym):
    """e(x) ... the function
    Omega ... the mrv set
    wsym ... the symbol which is going to be used for w

    Returns the rewritten e in terms of w and log(w)."""
    from .. import _core
    if not isinstance(Omega, SubsSet):
        raise TypeError("Omega should be an instance of SubsSet")
    return _core._gruntz_rewrite(_s(e), list(Omega.items()), list(Omega.rewrites.items()), _s(x), _s(wsym))


def mrv_leadterm(e, x):
    """Returns (c0, e0) for e."""
    from .. import _core
    return _core._gruntz_mrv_leadterm(_s(e), _s(x))


def limitinf(e, x):
    """Limit e(x) for x-> oo."""
    from .. import _core
    return _core._gruntz_limitinf(_s(e), _s(x))


def gruntz(e, z, z0, dir="+"):
    """Compute the limit of e(z) at the point z0 using the Gruntz algorithm."""
    from .. import _core, PoleError
    e, z, z0 = _s(e), _s(z), _s(z0)
    if not z.is_symbol:
        raise NotImplementedError("Second argument must be a Symbol")
    if str(dir) not in ("+", "-"):
        raise NotImplementedError("dir must be '+' or '-'")
    try:
        return _core._gruntz(e, z, z0, str(dir))
    except NotImplementedError:
        pass
    # the native engine with its tractable heuristics before Gruntz
    result = _core._limit(e, z, z0, str(dir), True)
    if result is None:
        raise PoleError('The native limit engine could not determine this limit')
    return result


__all__ = ['gruntz', 'compare', 'sign', 'mrv', 'rewrite', 'mrv_leadterm', 'limitinf', 'SubsSet']
