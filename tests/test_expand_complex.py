"""expand_complex, expand(complex=True), adjoint and sec/csc parts against SymPy."""
import pickle

import pytest
import sympy as sp
import symox as sx

from test_complex_parts import SOURCES, environment


def env(module):
    result = environment(module)
    result.update(g=module.Symbol('g'), z=module.Symbol('z'), sec=module.sec, csc=module.csc)
    return result


EXTRA = [
    'sec(x)', 'csc(x)', 'sec(r)', 'csc(r+I*t)', 'sin(Abs(x))', 'exp(x*y)', 'log(x*y)',
    'x*y**2', 'f(x)*g', '(x+I*y)*(x-I*y)', 'sqrt(x+I*y)', 'exp(exp(x))', 'sin(sin(x))',
    'log(log(x))', 'Abs(x+I*y)', 'arg(x+I*y)', 'x**(y+1)', '(x*y)**z', 'exp(x)**2',
    'log(x**2)', '2**x', 'x**x', '1/(1+x)', 'cosh(r**2+2*r)', 'sinh(r*(t+1))',
    'sin(r*(t+1))', 'exp(r*(t+1))', 'log(p*(t+1))', 'tan(x+y)', 're(x)*im(x)',
    'Rational(1,2)+I*x', 'exp(I*pi*x)', 'sqrt(2)*I*x', 'log(2+3*I)', '(1+I)**10',
    '(1+I)**Rational(1,3)', 'I**x', 'UnevaluatedExpr(x)*2', 'Float("1.5")*x', '-0.0*x',
]

MODES = {
    'expand_complex': lambda m, e: m.expand_complex(e),
    'expand_complex_shallow': lambda m, e: m.expand_complex(e, deep=False),
    'expand': lambda m, e: m.expand(e, complex=True),
    'method_shallow': lambda m, e: e.expand(deep=False, complex=True),
    'partial_hints': lambda m, e: e.expand(complex=True, mul=False, multinomial=False),
}

@pytest.mark.parametrize('mode', list(MODES))
@pytest.mark.parametrize('source', SOURCES + EXTRA)
def test_complex_expansion_matches_sympy(source, mode):
    expected = sp.sympify(eval(source, {'__builtins__': {}}, env(sp)))
    actual = sx.sympify(eval(source, {'__builtins__': {}}, env(sx)))
    try:
        oracle = MODES[mode](sp, expected)
    except Exception:
        with pytest.raises(Exception):
            MODES[mode](sx, actual)
        return
    result = MODES[mode](sx, actual)
    assert sx.to_sympy(result) == oracle, (source, mode, result, oracle)
    assert pickle.loads(pickle.dumps(result)) == result
    assert sx.sympify(result.srepr()) == result


@pytest.mark.parametrize('source', ['sec(x)', 'csc(x)', 'sec(r)', 'csc(r)', 'sec(r+I*t)',
                                    'csc(x+I*y)', 'sec(0)', 'csc(pi/2)'])
@pytest.mark.parametrize('deep', [False, True])
def test_reciprocal_trigonometric_parts_match_sympy(source, deep):
    expected = sp.sympify(eval(source, {'__builtins__': {}}, env(sp)))
    actual = sx.sympify(eval(source, {'__builtins__': {}}, env(sx)))
    assert tuple(map(sx.to_sympy, actual.as_real_imag(deep=deep))) == expected.as_real_imag(deep=deep)


@pytest.mark.parametrize('kind', ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix'])
def test_matrix_complex_expansion_keeps_storage(kind):
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    actual = getattr(sx, kind)([[x+sx.I*y, sx.exp(sx.I*x)], [0, sx.sqrt(x)]])
    expected = getattr(sp, kind)([[X+sp.I*Y, sp.exp(sp.I*X)], [0, sp.sqrt(X)]])
    result = actual.expand(complex=True)
    assert result.to_sympy() == expected.expand(complex=True)
    assert isinstance(result, getattr(sx, kind))
    assert bool(getattr(result, '_symox_sparse', False)) == ('Sparse' in kind)
    assert pickle.loads(pickle.dumps(result)) == result


def test_sparse_complex_expansion_only_touches_stored_entries():
    x = sx.Symbol('x')
    value = sx.ImmutableSparseMatrix(10**7, 10**7, {(3, 4): sx.exp(sx.I*x)})
    result = value.expand(complex=True)
    assert result.shape == value.shape
    assert result.todok() == {(3, 4): sx.expand_complex(sx.exp(sx.I*x))}


@pytest.mark.parametrize('source', ['x', 'r', 'I*x', '2+3*I', 'oo', 'nan', 'zoo', '0',
                                    'Float("1.5")+I*x', 'exp(I*x)', 'x*y'])
def test_adjoint_matches_sympy(source):
    expected = sp.sympify(eval(source, {'__builtins__': {}}, env(sp)))
    actual = sx.sympify(eval(source, {'__builtins__': {}}, env(sx)))
    assert sx.to_sympy(sx.adjoint(actual)) == sp.adjoint(expected)
    assert sx.to_sympy(sx.adjoint(source.replace('Float("1.5")', '1.5'))) == \
        sp.adjoint(source.replace('Float("1.5")', '1.5'))


def test_adjoint_of_matrices_matches_sympy():
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    for kind in ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix']:
        result = sx.adjoint(getattr(sx, kind)([[x, sx.I*y], [1, 2]]))
        assert result.to_sympy() == sp.adjoint(getattr(sp, kind)([[X, sp.I*Y], [1, 2]]))


@pytest.mark.parametrize('source', ['log(sqrt(x))', 'log(sqrt(p))', 'log(x**2)', 'log(x**r)',
                                    'log(1/x)', 'log(x**Rational(-1,2))', 'log(p**y)',
                                    'log(x**Rational(1,3))', 'log(x**Rational(3,2))'])
def test_expand_log_of_powers_matches_sympy(source):
    expected = sp.sympify(eval(source, {'__builtins__': {}}, env(sp)))
    actual = sx.sympify(eval(source, {'__builtins__': {}}, env(sx)))
    assert sx.to_sympy(sx.expand_log(actual)) == sp.expand_log(expected)
    assert sx.to_sympy(actual.expand()) == expected.expand()


def test_unsupported_nodes_are_rejected_explicitly():
    x, r = sx.Symbol('x'), sx.Symbol('r', real=True)
    f = sx.Function('f')
    for value in [sx.Derivative(f(x), x), sx.Piecewise((x, r > 0), (1, True)),
                  sx.Integral(x, x), sx.Eq(x, r)]:
        with pytest.raises(NotImplementedError):
            sx.expand_complex(value)
    with pytest.raises(NotImplementedError, match='complex'):
        sx.MatrixSymbol('A', 2, 2).expand(complex=True)
