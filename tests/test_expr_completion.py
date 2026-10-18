"""Public Expr method adapters must preserve native semantics and Limit metadata."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('source', [
    '0', '1', 'x', 'x**2+3*x+2', '1/x', 'factorial(x)',
    'gamma(x)', 'binomial(2*x,x)', '2**x', 'exp(x)', 'exp(x**2)',
    'sin(x)', 'factorial(x**2)', 'Piecewise((x,x>0),(1,True))',
    '1+Piecewise((x,x>0),(1,True))',
])
def test_is_hypergeometric(source):
    x = sx.Symbol('x')
    ours = eval(source, dict(vars(sx), x=x))
    ours = sx.sympify(ours)
    expected = sp.sympify(source)
    assert ours.is_hypergeometric(x) is expected.is_hypergeometric(sp.Symbol('x'))


@pytest.mark.parametrize('pattern,expr', [
    ('a*x+b', '2*x+3'), ('a*x+b', 'x**2+1'),
    ('a', 'sin(x)'), ('sin(a)', 'sin(x)'), ('a**2', 'x**2'),
    ('a*x+1', '2*x+2'), ('a*x+b', 'x'), ('x+1', 'x+1'),
])
@pytest.mark.parametrize('bound', [False, True])
def test_matches_and_prebound_wildcards(pattern, expr, bound):
    x, a, b = sx.Symbol('x'), sx.Wild('a'), sx.Wild('b')
    X, A, B = sp.Symbol('x'), sp.Wild('a'), sp.Wild('b')
    ours = {'x': x, 'a': a, 'b': b, 'sin': sx.sin}
    theirs = {'x': X, 'a': A, 'b': B, 'sin': sp.sin}
    bindings, reference_bindings = ({a: sx.Integer(2)}, {A: sp.Integer(2)}) if bound else ({}, {})
    actual = eval(pattern, ours).matches(eval(expr, ours), bindings)
    expected = eval(pattern, theirs).matches(eval(expr, theirs), reference_bindings)
    converted = None if actual is None else {k.name: v.to_sympy() for k, v in actual.items()}
    expected = None if expected is None else {k.name: v for k, v in expected.items()}
    assert converted == expected
    assert bindings == ({a: sx.Integer(2)} if bound else {})


@pytest.mark.parametrize('source', ['x', '1', 'x**2', 'sin(x)', 'cos(x)', 'exp(x)',
                                   'x**3+x**4', 'sqrt(x)', '1/x', 'log(x)'])
@pytest.mark.parametrize('direction', [1, -1, 1j])
def test_direction_and_limit_metadata(source, direction):
    x, X = sx.Symbol('x'), sp.Symbol('x')
    ours, expected = sx.sympify(source), sp.sympify(source)
    cdir = sx.I if direction == 1j else sx.sympify(direction)
    Cdir = sp.I if direction == 1j else sp.sympify(direction)
    assert ours.dir(x, cdir).to_sympy() == expected.dir(X, Cdir)
    assert sx.Limit(sx.sin(x)/x, x, 0).dir == '+'
    assert sx.Limit(sx.sin(x)/x, x, 0, '-').dir == '-'


@pytest.mark.parametrize('simplify', [False, True])
def test_nnf_method(simplify):
    a, b, c = sx.symbols('a b c')
    expressions = [~(a & b), sx.Implies(a, b), sx.Equivalent(a, b, c), sx.Xor(a, b), a]
    for expr in expressions:
        expected = sp.to_nnf(expr.to_sympy(), simplify=simplify)
        assert expr.to_nnf(simplify=simplify).to_sympy() == expected


def test_legacy_matching_exact_identity():
    x = sx.Symbol('x')
    assert x.matches(x,old=True)=={}
