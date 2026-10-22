"""SymPy's rewrite rules for the gamma and combinatorial function families.

The formulas and their applicability conditions follow the
``_eval_rewrite_as_<rule>`` methods of SymPy 1.14 (BSD-3-Clause; see
LICENSE-SYMPY) in ``functions/special/gamma_functions.py`` and
``beta_functions.py`` and ``functions/combinatorial/factorials.py`` and
``numbers.py`` (harmonic). No SymPy code is imported.

The traversal is SymPy's ``Basic._rewrite``: children are rewritten first,
then a node selected by the pattern classes is rewritten from its new
arguments. A rule returning None leaves the rebuilt node unchanged. Nodes
outside this table use the native rewrite of the rebuilt node.
"""
from . import _core


def _s():
    import symox
    return symox


def _gamma_tractable(e, z, limitvar=None, **hints):
    s = _s()
    return s.exp(s.loggamma(z))


def _gamma_factorial(e, z, **hints):
    return _s().factorial(z-1)


def _lowergamma_uppergamma(e, a, x, **hints):
    s = _s()
    return s.gamma(a)-s.uppergamma(a, x)


def _lowergamma_expint(e, a, x, **hints):
    s = _s()
    if a.is_integer and a.is_nonpositive:
        return e
    return e.rewrite(s.uppergamma).rewrite(s.expint)


def _uppergamma_lowergamma(e, a, x, **hints):
    s = _s()
    return s.gamma(a)-s.lowergamma(a, x)


def _uppergamma_tractable(e, a, x, **hints):
    s = _s()
    return s.exp(s.loggamma(a))-s.lowergamma(a, x)


def _uppergamma_expint(e, a, x, **hints):
    s = _s()
    return s.expint(1-a, x)*x**a


def _polygamma_zeta(e, n, z, **hints):
    s = _s()
    if n.is_integer and n.is_positive:
        return s.Integer(-1)**(n+1)*s.factorial(n)*s.zeta(n+1, z)
    return None


def _polygamma_harmonic(e, n, z, **hints):
    s = _s()
    if n.is_integer:
        if n.is_zero:
            return s.harmonic(z-1)-s.EulerGamma
        return s.Integer(-1)**(n+1)*s.factorial(n)*(s.zeta(n+1)-s.harmonic(z-1, n+1))
    return None


def _loggamma_intractable(e, z, **hints):
    s = _s()
    return s.log(s.gamma(z))


def _digamma_harmonic(e, z, **hints):
    s = _s()
    return s.harmonic(z-1)-s.EulerGamma


def _digamma_polygamma(e, z, **hints):
    return _s().polygamma(0, z)


def _trigamma_zeta(e, z, **hints):
    return _s().zeta(2, z)


def _trigamma_polygamma(e, z, **hints):
    return _s().polygamma(1, z)


def _trigamma_harmonic(e, z, **hints):
    s = _s()
    return -s.harmonic(z-1, 2)+s.pi**2/6


def _factorial_gamma(e, n, piecewise=True, **hints):
    return _s().gamma(n+1)


def _factorial_Product(e, n, **hints):
    s = _s()
    if n.is_nonnegative and n.is_integer:
        i = s.Dummy('i', integer=True)
        return s.Product(i, (i, 1, n))
    return None


def _subfactorial_factorial(e, arg, **hints):
    s = _s()
    i = s.Dummy('i')
    f = s.Integer(-1)**i/s.factorial(i)
    return s.factorial(arg)*s.summation(f, (i, 0, arg))


def _subfactorial_gamma(e, arg, piecewise=True, **hints):
    s = _s()
    return (s.Integer(-1)**(arg+1)*s.exp(-s.I*s.pi*arg)*s.lowergamma(arg+1, -1)
            + s.gamma(arg+1))*s.exp(-1)


def _subfactorial_uppergamma(e, arg, **hints):
    s = _s()
    return s.uppergamma(arg+1, -1)/s.E


def _factorial2_gamma(e, n, piecewise=True, **hints):
    s = _s()
    return 2**(n/2)*s.gamma(n/2+1)*s.Piecewise((1, s.Eq(s.Mod(n, 2), 0)),
                                               (s.sqrt(2/s.pi), s.Eq(s.Mod(n, 2), 1)))


def _rf_gamma(e, x, k, piecewise=True, **hints):
    s = _s()
    if not piecewise:
        if (x <= 0) == True:  # noqa: E712 - SymPy's exact truth test
            return s.Integer(-1)**k*s.gamma(1-x)/s.gamma(-k-x+1)
        return s.gamma(x+k)/s.gamma(x)
    return s.Piecewise((s.gamma(x+k)/s.gamma(x), x > 0),
                       (s.Integer(-1)**k*s.gamma(1-x)/s.gamma(-k-x+1), True))


def _rf_FallingFactorial(e, x, k, **hints):
    return _s().FallingFactorial(x+k-1, k)


def _rf_factorial(e, x, k, **hints):
    s = _s()
    if x.is_integer and k.is_integer:
        return s.Piecewise((s.factorial(k+x-1)/s.factorial(x-1), x > 0),
                           (s.Integer(-1)**k*s.factorial(-x)/s.factorial(-k-x), True))
    return None


def _rf_binomial(e, x, k, **hints):
    s = _s()
    if k.is_integer:
        return s.factorial(k)*s.binomial(x+k-1, k)
    return None


def _limit_direction(k, limitvar):
    s = _s()
    k_lim = k.subs(limitvar, s.oo)
    if k_lim == s.oo:
        return 1
    if k_lim == -s.oo:
        return -1
    return 0


def _rf_tractable(e, x, k, limitvar=None, **hints):
    s = _s()
    if limitvar:
        direction = _limit_direction(k, limitvar)
        if direction == 1:
            return s.gamma(x+k).rewrite('tractable', deep=True)/s.gamma(x)
        if direction == -1:
            return s.Integer(-1)**k*s.gamma(1-x)/s.gamma(-k-x+1).rewrite('tractable', deep=True)
    return e.rewrite(s.gamma).rewrite('tractable', deep=True)


def _ff_gamma(e, x, k, piecewise=True, **hints):
    s = _s()
    if not piecewise:
        if (x < 0) == True:  # noqa: E712 - SymPy's exact truth test
            return s.Integer(-1)**k*s.gamma(k-x)/s.gamma(-x)
        return s.gamma(x+1)/s.gamma(x-k+1)
    return s.Piecewise((s.gamma(x+1)/s.gamma(x-k+1), x >= 0),
                       (s.Integer(-1)**k*s.gamma(k-x)/s.gamma(-x), True))


def _ff_RisingFactorial(e, x, k, **hints):
    return _s().RisingFactorial(x-k+1, k)


def _ff_binomial(e, x, k, **hints):
    s = _s()
    if k.is_integer:
        return s.factorial(k)*s.binomial(x, k)
    return None


def _ff_factorial(e, x, k, **hints):
    s = _s()
    if x.is_integer and k.is_integer:
        return s.Piecewise((s.factorial(x)/s.factorial(-k+x), x >= 0),
                           (s.Integer(-1)**k*s.factorial(k-x-1)/s.factorial(-x-1), True))
    return None


def _ff_tractable(e, x, k, limitvar=None, **hints):
    s = _s()
    if limitvar:
        direction = _limit_direction(k, limitvar)
        if direction == 1:
            return s.Integer(-1)**k*s.gamma(k-x).rewrite('tractable', deep=True)/s.gamma(-x)
        if direction == -1:
            return s.gamma(x+1)/s.gamma(x-k+1).rewrite('tractable', deep=True)
    return e.rewrite(s.gamma).rewrite('tractable', deep=True)


def _binomial_factorial(e, n, k, **hints):
    s = _s()
    return s.factorial(n)/(s.factorial(k)*s.factorial(n-k))


def _binomial_gamma(e, n, k, piecewise=True, **hints):
    s = _s()
    return s.gamma(n+1)/(s.gamma(k+1)*s.gamma(n-k+1))


def _binomial_tractable(e, n, k, limitvar=None, **hints):
    return _binomial_gamma(e, n, k).rewrite('tractable')


def _binomial_FallingFactorial(e, n, k, **hints):
    s = _s()
    if k.is_integer:
        return s.FallingFactorial(n, k)/s.factorial(k)
    return None


def _beta_gamma(e, x, y, piecewise=True, **hints):
    s = _s()
    return s.gamma(x)*s.gamma(y)/s.gamma(x+y)


def _beta_Integral(e, x, y, **hints):
    s = _s()
    from .core.symbol import uniquely_named_symbol
    t = s.Dummy(uniquely_named_symbol('t', [x, y]).name)
    return s.Integral(t**(x-1)*(1-t)**(y-1), (t, 0, 1))


def _harmonic_polygamma(e, n, m=None, **hints):
    s = _s()
    m = s.Integer(1) if m is None else m
    if m.is_integer and m.is_positive:
        return s.Piecewise((s.polygamma(0, n+1)+s.EulerGamma, s.Eq(m, 1)),
                           (s.Integer(-1)**m*(s.polygamma(m-1, 1)-s.polygamma(m-1, n+1))/s.gamma(m), True))
    return None


def _harmonic_digamma(e, n, m=None, **hints):
    return e.rewrite(_s().polygamma)


def _harmonic_Sum(e, n, m=None, **hints):
    s = _s()
    k = s.Dummy('k', integer=True)
    m = s.Integer(1) if m is None else m
    return s.Sum(k**(-m), (k, 1, n))


def _harmonic_zeta(e, n, m=None, **hints):
    s = _s()
    m = s.Integer(1) if m is None else m
    return s.Piecewise((s.digamma(n+1)+s.EulerGamma, s.Eq(m, 1)),
                       (s.zeta(m)-s.zeta(m, n+1), True))


def _harmonic_tractable(e, n, m=None, limitvar=None, **hints):
    s = _s()
    m = s.Integer(1) if m is None else m
    pg = e.rewrite(s.polygamma)
    if getattr(pg, 'class_name', None) != 'harmonic':
        return pg.rewrite('tractable', deep=True)
    if (m-1).is_nonzero:
        return (s.zeta(m)-s.zeta(m, n+1)).rewrite('tractable', deep=True)
    return None


RULES = {
    ('gamma', 'tractable'): _gamma_tractable,
    ('gamma', 'factorial'): _gamma_factorial,
    ('lowergamma', 'uppergamma'): _lowergamma_uppergamma,
    ('lowergamma', 'expint'): _lowergamma_expint,
    ('uppergamma', 'lowergamma'): _uppergamma_lowergamma,
    ('uppergamma', 'tractable'): _uppergamma_tractable,
    ('uppergamma', 'expint'): _uppergamma_expint,
    ('polygamma', 'zeta'): _polygamma_zeta,
    ('polygamma', 'harmonic'): _polygamma_harmonic,
    ('loggamma', 'intractable'): _loggamma_intractable,
    ('digamma', 'harmonic'): _digamma_harmonic,
    ('digamma', 'polygamma'): _digamma_polygamma,
    ('trigamma', 'zeta'): _trigamma_zeta,
    ('trigamma', 'polygamma'): _trigamma_polygamma,
    ('trigamma', 'harmonic'): _trigamma_harmonic,
    ('factorial', 'gamma'): _factorial_gamma,
    ('factorial', 'Product'): _factorial_Product,
    ('subfactorial', 'factorial'): _subfactorial_factorial,
    ('subfactorial', 'gamma'): _subfactorial_gamma,
    ('subfactorial', 'uppergamma'): _subfactorial_uppergamma,
    ('factorial2', 'gamma'): _factorial2_gamma,
    ('RisingFactorial', 'gamma'): _rf_gamma,
    ('RisingFactorial', 'FallingFactorial'): _rf_FallingFactorial,
    ('RisingFactorial', 'factorial'): _rf_factorial,
    ('RisingFactorial', 'binomial'): _rf_binomial,
    ('RisingFactorial', 'tractable'): _rf_tractable,
    ('FallingFactorial', 'gamma'): _ff_gamma,
    ('FallingFactorial', 'RisingFactorial'): _ff_RisingFactorial,
    ('FallingFactorial', 'binomial'): _ff_binomial,
    ('FallingFactorial', 'factorial'): _ff_factorial,
    ('FallingFactorial', 'tractable'): _ff_tractable,
    ('binomial', 'factorial'): _binomial_factorial,
    ('binomial', 'gamma'): _binomial_gamma,
    ('binomial', 'tractable'): _binomial_tractable,
    ('binomial', 'FallingFactorial'): _binomial_FallingFactorial,
    ('beta', 'gamma'): _beta_gamma,
    ('beta', 'Integral'): _beta_Integral,
    ('harmonic', 'polygamma'): _harmonic_polygamma,
    ('harmonic', 'digamma'): _harmonic_digamma,
    ('harmonic', 'trigamma'): _harmonic_digamma,
    ('harmonic', 'Sum'): _harmonic_Sum,
    ('harmonic', 'zeta'): _harmonic_zeta,
    ('harmonic', 'tractable'): _harmonic_tractable,
}

from ._special_rewrites import RULES as _SPECIAL_RULES  # noqa: E402
for _key, _body in _SPECIAL_RULES.items():
    RULES.setdefault(_key, _body)

_HEADS = {}
for _head, _rule in RULES:
    _HEADS.setdefault(_rule, []).append(_head)


def _name(r):
    if r is abs:
        return 'Abs'
    if isinstance(r, str):
        return r
    name = getattr(r, '__name__', None)
    return name if name is not None else type(r).__name__


def _class_names(e):
    return {base.__name__ for base in type(e).__mro__}


_PYTHON_HEADS = {'SingularityFunction', 'betainc', 'betainc_regularized', 'riemann_xi', 'stieltjes'}


def install():
    s = _s()
    Expr = _core.Expr
    old_rewrite = Expr.rewrite

    def has_head(expr, head):
        cls = getattr(s, head)
        if isinstance(cls, type) and not isinstance(cls, _core.Expr) and head in _PYTHON_HEADS:
            # Heads implemented as held Python classes are not seen by the
            # native has(); look for their class name instead.
            stack = [expr]
            while stack:
                e = stack.pop()
                if getattr(e, 'class_name', None) == head:
                    return True
                # Piecewise pairs and tuples are not native expressions but
                # hold sub-expressions, as any SymPy Basic.
                stack.extend(a for a in getattr(e, 'args', ())
                             if isinstance(a, _core.Expr) or (not isinstance(a, (str, bytes))
                                                              and hasattr(a, 'args')))
            return False
        return expr.has(cls)

    def heads_present(expr, rule):
        return any(has_head(expr, head) for head in _HEADS[rule])

    def rewrite(self, *args, deep=True, **hints):
        if not args:
            return old_rewrite(self, *args, deep=deep, **hints)
        rule = _name(args[-1])
        pattern = args[:-1]
        if rule not in _HEADS or not heads_present(self, rule):
            return old_rewrite(self, *args, deep=deep, **hints)
        if pattern and isinstance(pattern[0], (list, tuple, set, frozenset)):
            pattern = tuple(pattern[0])
        if pattern:
            # SymPy keeps only the pattern classes present in the expression.
            kept = []
            for p in pattern:
                try:
                    present = self.has(p)
                except Exception:
                    present = True
                if present:
                    kept.append(p)
            if not kept:
                return self
            pattern = tuple(kept)
        names = {_name(p) for p in pattern}

        def node(e, new_args):
            changed = new_args != e.args
            if not names or names & _class_names(e):
                body = RULES.get((getattr(e, 'class_name', None), rule))
                if body is not None:
                    rewritten = body(e, *new_args, **hints)
                    if rewritten is not None:
                        return rewritten
                    return e.func(*new_args) if changed else e
                rebuilt = e.func(*new_args) if changed else e
                if isinstance(rebuilt, Expr):
                    return old_rewrite(rebuilt, *pattern, args[-1], deep=False, **hints)
                return rebuilt
            return e.func(*new_args) if changed else e

        def visit(e):
            if isinstance(e, Expr) and not heads_present(e, rule):
                # SymPy's traversal of a subtree without these heads is the
                # native rewrite of that subtree.
                return old_rewrite(e, *pattern, args[-1], deep=True, **hints)
            children = getattr(e, 'args', None)
            if not children or not hasattr(e, 'func'):
                return e
            # Structural containers (Piecewise pairs, tuples) are traversed
            # like any SymPy Basic.
            return node(e, tuple(visit(a) for a in children))

        if not deep:
            return node(self, self.args)
        return visit(self)

    Expr.rewrite = rewrite
