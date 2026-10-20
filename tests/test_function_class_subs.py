"""Substitution of function classes, compared with SymPy."""
import pytest
import sympy as sp
import symox as sx


def _env(m):
    x, y = m.symbols('x y')
    f, g, h = m.symbols('f g h', cls=m.Function)
    return m, x, y, f, g, h


CASES = {
    'abs id': lambda m, x, y, f, g, h: (m.Abs(x)+m.Abs(y-1)).subs(m.Abs, m.Id),
    'sin cos': lambda m, x, y, f, g, h: (m.sin(x)+1).subs(m.sin, m.cos),
    'undefined': lambda m, x, y, f, g, h: f(x).subs(f, g),
    'lambda': lambda m, x, y, f, g, h: f(x).subs(f, m.Lambda(x, x**2)),
    'lambda arity': lambda m, x, y, f, g, h: f(x, y).subs(f, m.Lambda(x, x**2)),
    'mixed arity': lambda m, x, y, f, g, h: (f(x, y)+f(x)).subs(f, m.Lambda(x, x**2)),
    'fixed nargs': lambda m, x, y, f, g, h: (f(x, y)+f(x)).subs(f, m.Function('h2', nargs=2)),
    'lambda multi': lambda m, x, y, f, g, h: f(x, y).subs(f, m.Lambda((x, y), x-y)),
    'builtin new': lambda m, x, y, f, g, h: f(x).subs(f, m.sin),
    'dict key': lambda m, x, y, f, g, h: m.Abs(x).subs({m.Abs: m.Id}),
    'nested': lambda m, x, y, f, g, h: f(f(x)+1).subs(f, g),
    'derivative': lambda m, x, y, f, g, h: m.Derivative(f(x), x).subs(f, m.Lambda(x, x**2)),
    'derivative head': lambda m, x, y, f, g, h: m.Derivative(f(x), x).subs(f, g),
    'diff builtin': lambda m, x, y, f, g, h: f(x).diff(x).subs(f, m.sin),
    'second derivative': lambda m, x, y, f, g, h: f(x).diff(x, 2).subs(f, m.Lambda(x, m.sin(x))).doit(),
    'Subs': lambda m, x, y, f, g, h: m.Subs(f(x).diff(x), x, 1).subs(f, g),
    'Subs bound': lambda m, x, y, f, g, h: m.Subs(f(y), y, x).subs(f, g),
    'to number': lambda m, x, y, f, g, h: f(x).subs(f, 1),
    'to symbol': lambda m, x, y, f, g, h: f(x).subs(f, y),
    'dict order': lambda m, x, y, f, g, h: (f(x)+f(y)).subs({f: g, f(x): y}),
    'list order': lambda m, x, y, f, g, h: (f(x)+f(y)).subs([(f, g), (f(x), y)]),
    'list order reversed': lambda m, x, y, f, g, h: (f(x)+f(y)).subs([(f(x), y), (f, g)]),
    'chain dict': lambda m, x, y, f, g, h: (f(x)+g(x)).subs({f: g, g: h}),
    'chain list': lambda m, x, y, f, g, h: (f(x)+g(x)).subs([(f, g), (g, h)]),
    'chain list reversed': lambda m, x, y, f, g, h: (f(x)+g(x)).subs([(g, h), (f, g)]),
    'mixed': lambda m, x, y, f, g, h: (f(x)+m.sin(x)+y).subs({f: g, m.sin: m.cos, y: 2}),
    'infinity first': lambda m, x, y, f, g, h: (f(x)+y).subs({f: g, y: m.oo}),
    'max min': lambda m, x, y, f, g, h: m.Max(x, y).subs(m.Max, m.Min),
    'max sin': lambda m, x, y, f, g, h: m.Max(x, y).subs(m.Max, m.sin),
    'two arguments': lambda m, x, y, f, g, h: (m.atan2(y, x)+m.Max(x, 1)).subs({m.atan2: f, m.Max: m.Min}),
    'log': lambda m, x, y, f, g, h: (m.log(x)+m.log(x, 2)).subs(m.log, f),
    'exp head': lambda m, x, y, f, g, h: (m.exp(x)+m.exp(2*y)).subs(m.exp, m.cos),
    'sin to exp': lambda m, x, y, f, g, h: (m.sin(x)**2).subs(m.sin, m.exp),
    'exp to number': lambda m, x, y, f, g, h: (m.exp(x)+m.exp(f(y))).subs(m.exp, 2),
    'exp to symbol': lambda m, x, y, f, g, h: m.exp(x).subs(m.exp, y),
    'exp to undefined': lambda m, x, y, f, g, h: m.exp(x).subs(m.exp, f),
    'applied in power': lambda m, x, y, f, g, h: (f(x)**2+1/f(y)).subs(f, m.exp),
    'cot keeps head': lambda m, x, y, f, g, h: (m.cot(x)+m.tan(x)).subs(m.cot, m.tan),
    'inside cot': lambda m, x, y, f, g, h: m.cot(f(x)).subs(f, g),
    'piecewise': lambda m, x, y, f, g, h: m.Piecewise((f(x), f(y) > 0), (m.sin(x), True)).subs(f, g),
    'Eq': lambda m, x, y, f, g, h: m.Eq(f(x), m.sin(y)).subs(m.sin, m.cos),
    'string old': lambda m, x, y, f, g, h: f(x).subs('f', g),
    'string key': lambda m, x, y, f, g, h: f(x).subs({'f': g}),
    'lambda old': lambda m, x, y, f, g, h: f(x).subs(m.Lambda(x, f(x)), g),
    'absent head': lambda m, x, y, f, g, h: (x+1).subs(f, g),
    'inside argument': lambda m, x, y, f, g, h: m.sin(f(x)**2).subs(f, m.Lambda(x, x+1)),
    'assumptions kept': lambda m, x, y, f, g, h: f(x).subs(f, m.Function('k', positive=True)).is_positive,
    'kronecker old': lambda m, x, y, f, g, h: m.KroneckerDelta(x, y).subs(m.KroneckerDelta, f),
    'kronecker new': lambda m, x, y, f, g, h: f(x, y).subs(f, m.KroneckerDelta),
    'mod': lambda m, x, y, f, g, h: m.Mod(x, 3).subs(m.Mod, f),
    'frac': lambda m, x, y, f, g, h: f(x).subs(f, m.frac),
    'frac arity': lambda m, x, y, f, g, h: f(x, y).subs(f, m.frac),
    'xreplace ignores heads': lambda m, x, y, f, g, h: (f(x)+m.sin(x)).xreplace({f: g, m.sin: m.cos, x: y}),
    'xreplace head only': lambda m, x, y, f, g, h: m.Abs(x).xreplace({m.Abs: m.Id}),
    'matrix': lambda m, x, y, f, g, h: m.Matrix([f(x), m.sin(x)]).subs({f: g, m.sin: m.cos}),
    'sparse matrix': lambda m, x, y, f, g, h: m.SparseMatrix([[f(x), 0], [0, m.sin(x)]]).subs(f, g),
    'immutable matrix': lambda m, x, y, f, g, h: m.ImmutableMatrix([f(x), m.sin(x)]).subs(m.sin, m.cos),
    'array': lambda m, x, y, f, g, h: m.Array([f(x), m.sin(x)]).subs(f, g),
    'tuple': lambda m, x, y, f, g, h: m.Tuple(f(x), 1).subs(f, g),
}


@pytest.mark.parametrize('name', list(CASES))
def test_function_class_substitution_matches_sympy(name):
    case = CASES[name]
    assert str(case(*_env(sx))) == str(case(*_env(sp)))


def test_function_class_substitution_keeps_symox_types():
    m, x, y, f, g, h = _env(sx)
    assert type(m.Matrix([f(x)]).subs(f, g)) is type(m.Matrix([x]))
    assert f(x).subs(f, g) == g(x)
    assert f(x).subs([(f, g)]) == g(x)
    assert f(x).subs(iter([(f, g)])) == g(x)
    assert f(x).subs([iter((f, g))]) == g(x)


def test_simultaneous_function_class_replacement_is_rejected_like_sympy():
    for module in (sx, sp):
        _, x, y, f, g, h = _env(module)
        with pytest.raises(ValueError, match='Dummy commutativity'):
            f(x).subs(f, g, simultaneous=True)
        with pytest.raises(ValueError, match='Dummy commutativity'):
            (f(x)+x).subs({f: g, x: y}, simultaneous=True)


def test_function_class_substitution_into_limits_raises_type_error():
    # SymPy fails with a TypeError on a class property; the native
    # implementation raises the same exception type with an explicit message.
    for module in (sx, sp):
        _, x, y, f, g, h = _env(module)
        with pytest.raises(TypeError):
            module.Integral(f(x), (x, 0, 1)).subs(f, g)


def test_replace_accepts_lambda_for_type_query():
    for module in (sx, sp):
        x = module.Symbol('x')
        expr = module.sin(x)+module.cos(x)
        result = expr.replace(module.sin, module.Lambda(x, x**2))
        assert str(result) == 'x**2 + cos(x)'


@pytest.mark.parametrize('source', [
    '(x + 1)**(-2)', '(x + y)**(-1)', '1/(x*(x + 1))**2', '(1/x + 1)**3',
    '(x/y + 1)**(-3)', 'exp(-x)/(x + 1)', '(x**2 + 2*x + 1)/(x + 1)**2', 'x**(-y)',
])
def test_as_numer_denom_matches_sympy(source):
    names = {'x': sx.Symbol('x'), 'y': sx.Symbol('y'), 'exp': sx.exp}
    reference = {'x': sp.Symbol('x'), 'y': sp.Symbol('y'), 'exp': sp.exp}
    actual = eval(source, names).as_numer_denom()
    expected = eval(source, reference).as_numer_denom()
    assert tuple(map(str, actual)) == tuple(map(str, expected))
