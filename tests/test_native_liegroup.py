"""Native lie_group solver (symox::ode::lie_group) against live SymPy 1.14.

``infinitesimals``, the ``lie_heuristic_*`` functions, ``checkinfsol``,
``_lie_group_remove``, ``classify_ode`` and ``dsolve(hint='lie_group')``
must give SymPy's value, or raise SymPy's exception class and message, and
the native calls must not decline to the Python reference body.
"""
import pytest
import sympy as sp
import symox as sx
from symox import _core
from symox.solvers.ode import lie_group as sx_lie
from sympy.solvers.ode import lie_group as sp_lie


def namespace(module):
    ns = dict(vars(module))
    for name in 'x a b c a0 a1 a2 a3 a4 m n'.split():
        ns[name] = module.Symbol(name)
    ns['f'] = module.Function('f')
    ns['F'] = module.Function('F')
    return ns


def to_sympy(value):
    if isinstance(value, dict):
        return {to_sympy(k): to_sympy(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return type(value)(to_sympy(v) for v in value)
    if isinstance(value, (bool, int)) or value is None:
        return value
    return value.to_sympy()


def outcome(call):
    try:
        return True, call()
    except Exception as error:  # the class and message are compared
        return False, (type(error).__name__, str(error))


@pytest.fixture
def declines(monkeypatch):
    """The native ``ode_lie_*`` calls that returned ``None``."""
    seen = []
    for name in dir(_core):
        if name.startswith('ode_lie_') and name != 'ode_lie_group_remove':
            original = getattr(_core, name)

            def call(*args, original=original, name=name):
                value = original(*args)
                if value is None:
                    seen.append(name)
                return value
            monkeypatch.setattr(_core, name, call)
    return seen


def assert_same(source, call, declines):
    ok_x, got = outcome(lambda: call(sx, namespace(sx), eval(source, namespace(sx))))
    ok_p, ref = outcome(lambda: call(sp, namespace(sp), eval(source, namespace(sp))))
    assert ok_x == ok_p, (got, ref)
    if ok_x:
        assert to_sympy(got) == ref
    else:
        assert got == ref
    assert not declines


EQUATIONS = [
    'Eq(f(x).diff(x), x**2*f(x))',
    'f(x).diff(x) + a*f(x) - c*exp(b*x)',
    'f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)',
    '(1 + 2*x)*f(x).diff(x) + 2 - 4*exp(-f(x))',
    'f(x).diff(x) - (a4*x**4 + a3*x**3 + a2*x**2 + a1*x + a0)**Rational(-1, 2)',
    'x**2*f(x).diff(x) - f(x) + x**2*exp(x - (1/x))',
    'f(x).diff(x) - (f(x)/x)*(x*log(x**2/f(x)) + 2)',
    'x*(f(x).diff(x)) + 1 - f(x)**2',
    'x**2*f(x).diff(x) + x*f(x) + f(x)**2 + x**2',
    'f(x).diff(x) - F(a*x + b*f(x))',
    'f(x).diff(x) - x**(a - 1)*(f(x)**(1 - b))*F(x**a/a + f(x)**b/b)',
    'f(x).diff(x) + tan(F(x**2 + f(x)**2) + atan(x/f(x)))',
    'f(x).diff(x) - f(x)**2',
    'f(x).diff(x) - exp(x + f(x))',
]

HINTS = ['default', 'all', 'abaco1_simple', 'abaco1_product', 'abaco2_similar',
         'abaco2_unique_unknown', 'abaco2_unique_general', 'linear',
         'function_sum', 'bivariate', 'chi', 'nonexistent']


# SymPy 1.14 needs more than two minutes for these searches.
SLOW = {(EQUATIONS[i], hint) for i in (10, 11) for hint in ('all', 'abaco2_unique_general')}

# Gaps outside this lot (symox integrate splits the logarithm of a product,
# ``log(y - 1) + log(y + 1)`` for ``log(y**2 - 1)`` and ``log(2*p)/2`` for
# ``log(p)/2``): the generator is valid but not SymPy's form.
# (``default`` of the same equation is compared up to a constant factor.)
INTEGRATE_GAP = {(EQUATIONS[4], 'abaco1_simple'), (EQUATIONS[7], 'abaco1_simple')}

# Gap outside this lot: symox multiplication cancels ``exp(x + y)*exp(-x - y)``
# to 1, SymPy keeps the product, so ``hy/h`` has no ``x`` only in symox and
# abaco1_simple finds ``eta = exp(f)`` where SymPy finds nothing.
MUL_GAP = {(EQUATIONS[13], 'default'), (EQUATIONS[13], 'abaco1_simple')}


def infinitesimal_cases():
    for source in EQUATIONS:
        for hint in HINTS:
            if (source, hint) in SLOW:
                continue
            marks = [pytest.mark.xfail(strict=True, reason='integrate log form')] \
                if (source, hint) in INTEGRATE_GAP else []
            if (source, hint) in MUL_GAP:
                marks = [pytest.mark.xfail(strict=True, reason='exp product cancels')]
            yield pytest.param(source, hint, marks=marks)


def proportional(g, r):
    if set(g) != set(r) or any(g[k] != 0 for k in r if r[k] == 0):
        return False
    ratios = {sp.simplify(g[k]/r[k]) for k in r if r[k] != 0}
    return len(ratios) == 1 and ratios.pop().is_number


def classes(gens):
    out = []
    for g in gens:
        if not any(proportional(g, c) for c in out):
            out.append(g)
    return out


# bivariate and chi solve for a set of unknown coefficients and set the free
# ones to 1: SymPy's generator (its sign, ``x*d/dx + (f**2 - 1)*d/df`` or
# ``x*d/dx + (1 - f**2)*d/df``) and whether 'all' drops it as a duplicate
# follow the hash seed.
HASH_DEPENDENT = ('all', 'default', 'bivariate', 'chi')


def same_modulo_hash_order(eq, got, ref):
    """The same generators when they are equal up to a constant factor;
    otherwise as many classes, each symox generator verified by SymPy."""
    if all(any(proportional(g, r) for r in ref) for g in got) and \
            all(any(proportional(g, r) for g in got) for r in ref):
        return True
    return len(classes(got)) == len(classes(ref)) and \
        all(valid for valid, _ in sp.solvers.ode.checkinfsol(eq, got))


@pytest.mark.parametrize('source,hint', infinitesimal_cases())
def test_infinitesimals_match_sympy(source, hint, declines):
    def call(mod, ns, eq):
        from importlib import import_module
        return import_module(mod.__name__ + '.solvers.ode').infinitesimals(eq, hint=hint)
    if hint not in HASH_DEPENDENT:
        return assert_same(source, call, declines)
    eq = eval(source, namespace(sp))
    ok_x, got = outcome(lambda: call(sx, namespace(sx), eval(source, namespace(sx))))
    ok_p, ref = outcome(lambda: call(sp, namespace(sp), eq))
    assert ok_x == ok_p, (got, ref)
    assert same_modulo_hash_order(eq, to_sympy(got), ref) if ok_x else got == ref
    assert not declines


@pytest.mark.parametrize('hint', ['function_sum', 'abaco1_simple', 'linear'])
def test_second_derivatives_are_factored_like_sympy(hint, declines):
    # ``(1/hinv).diff(x, 2)`` ends with ``factor_terms(signsimp(...))`` in
    # SymPy's Derivative; without it the separation of function_sum fails.
    def call(mod, ns, eq):
        from importlib import import_module
        return import_module(mod.__name__ + '.solvers.ode').infinitesimals(eq, hint=hint)
    assert_same('x**(n*(m + 1) - m)*(f(x).diff(x)) - a*f(x)**n - b*x**(n*(m + 1))', call, declines)


@pytest.mark.parametrize('source', [
    # the length of SymPy's 'all' list follows the hash seed (chi above)
    pytest.param(s, marks=pytest.mark.xfail(strict=False, reason='hash-seed dependent'))
    if s == EQUATIONS[7] else s for s in EQUATIONS[:9]])
def test_checkinfsol_matches_sympy(source, declines):
    def call(mod, ns, eq):
        from importlib import import_module
        ode = import_module(mod.__name__ + '.solvers.ode')
        x, f = ns['x'], ns['f']
        xi, eta = mod.Function('xi')(x, f(x)), mod.Function('eta')(x, f(x))
        found = ode.infinitesimals(eq, hint='all')
        return ode.checkinfsol(eq, found) + ode.checkinfsol(eq, [{xi: x, eta: f(x)**2}])
    assert_same(source, call, declines)


# abaco2_similar: SymPy's ``cancel`` of the ``zoo`` terms of ``hy/hyy`` gives
# ``nan`` (and the generator ``nan``); symox's cancel keeps ``y`` and goes on
# to the inverse ODE (gap of cancel with ``zoo``, outside this lot).
@pytest.mark.parametrize('name', ['abaco1_simple', 'abaco1_product', 'bivariate', 'chi', 'linear',
                                  pytest.param('abaco2_similar', marks=pytest.mark.xfail(
                                      strict=True, reason='cancel with zoo'))])
@pytest.mark.parametrize('comp', [False, True])
def test_public_heuristics_match_sympy(name, comp, declines):
    def call(mod, ns, h):
        x, f = ns['x'], ns['f']
        y = mod.Dummy('y')
        h = h.subs(f(x), y)
        u = mod.Dummy('u')
        hinv = ((1/h).subs([(x, u), (y, x)])).subs(u, y)
        module = sx_lie if mod is sx else sp_lie
        match = {'h': h, 'func': f(x), 'hx': h.diff(x), 'hy': h.diff(y), 'y': y, 'hinv': hinv}
        value = getattr(module, 'lie_heuristic_' + name)(match, comp)
        return value and [{k: v.subs(y, f(x)) for k, v in d.items()} for d in value]
    assert_same('x**2*f(x) + f(x)/x', call, declines)


@pytest.mark.parametrize('source', [
    'F(x - y) + x', 'F(x)*y', '(x + F(y))**2', 'F(x*y)', 'x**F(y)'])
def test_lie_group_remove_matches_sympy(source):
    ns_x, ns_p = namespace(sx), namespace(sp)
    ns_x['y'], ns_p['y'] = sx.Symbol('y'), sp.Symbol('y')
    assert sx_lie._lie_group_remove(eval(source, ns_x)).to_sympy() == \
        sp_lie._lie_group_remove(eval(source, ns_p))


@pytest.mark.parametrize('source', EQUATIONS[:9] + ['f(x).diff(x)**2', 'f(x).diff(x, 2) + f(x)'])
def test_classify_lists_lie_group_like_sympy(source):
    def call(mod, ns, eq):
        return 'lie_group' in mod.classify_ode(eq, ns['f'](ns['x']))
    assert_same(source, call, [])


@pytest.mark.parametrize('source', [
    'Eq(f(x).diff(x), x**2*f(x))',
    'f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)',
    # solve(C1 + log(2*x + 1)/2 + log(exp(y) - 2)/2, y) is SymPy's root in
    # another form (solve gap outside this lot)
    pytest.param('(1 + 2*x)*f(x).diff(x) + 2 - 4*exp(-f(x))',
                 marks=pytest.mark.xfail(strict=True, reason='solve root form')),
    'x**2*f(x)**2 + x*Derivative(f(x), x)',
    'x*diff(f(x),x) - f(x) - x/log(x)',
    '(f(x).diff(x)-f(x)) * (f(x).diff(x)+f(x))',
    'f(x).diff(x) - f(x)**2',
])
def test_dsolve_lie_group_matches_sympy(source, declines):
    def call(mod, ns, eq):
        return mod.dsolve(eq, ns['f'](ns['x']), hint='lie_group')
    assert_same(source, call, declines)


def test_user_infinitesimals_match_sympy(declines):
    def call(mod, ns, eq):
        x, f = ns['x'], ns['f']
        return mod.dsolve(eq, f(x), hint='lie_group', xi=mod.sqrt(f(x) - 1)/mod.sqrt(f(x) + 1), eta=0)
    assert_same('x*(f(x).diff(x)) + 1 - f(x)**2', call, declines)
