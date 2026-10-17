"""Differential structural queries, including scopes and held arithmetic."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def data(module):
    x, y, z = module.symbols('x y z')
    n = module.Symbol('n', integer=True, nonnegative=True)
    f = module.Function('f')
    expressions = [x, x+y, x/y, x**-2, (x+y)**-2, module.sqrt(x),
        module.sin(x), module.sin(x)+1, module.sqrt(x)+y, x**n, 2**x,
        f(x), module.Derivative(f(x), x), module.Integral(f(x), (x, 0, 1)),
        module.Integral(f(x), (x, 0, x)), module.Integral(1, x),
        module.Lambda(x, f(x)+y), module.Subs(f(x), x, y),
        module.Limit(f(x), x, y), module.sin(x+y+z),
        module.sin(x*y*z), (x*y*z)**-1, module.nan, module.oo, module.zoo,
        module.Integer(2), module.Float('2.5'), module.pi,
        module.Pow(x, module.Float('2.0'), evaluate=False),
        module.Pow(x, module.Rational(-2, 3), evaluate=False),
        module.exp(x)+module.cos(y), (x+y)**10,
        module.Sum(f(x), (x, 0, y)), module.Function('shared')(x, x),
        module.Integral(f(x)+f(y), (x, 0, y), (y, 0, z)),
        module.sqrt(f(x)), f(x)**2, module.Pow(module.sin(x), y, evaluate=False)]
    generators = [(), (x,), (y,), (x, y), (f(x),), (2**x,), (x+y,), (x*y,)]
    return expressions, generators


@pytest.mark.parametrize('index', range(38))
@pytest.mark.parametrize('selection', range(8))
def test_polynomial_rational_and_algebraic_queries_match_sympy(index, selection):
    actual, generators = data(sx)
    reference, ref_generators = data(sp)
    for method in ('is_polynomial', 'is_rational_function', 'is_algebraic_expr'):
        result = getattr(actual[index], method)(*generators[selection])
        expected = getattr(reference[index], method)(*ref_generators[selection])
        assert result is expected, (method, index, selection, result, expected)


@pytest.mark.parametrize('source', [
    'Integral(x+y+z,(x,0,1))', 'Integral(x*y*z,(x,1,2))',
    'Integral(sin(x+y+z),(x,0,1))', 'Lambda(x,x+y+z)',
    'Lambda(x,x*y*z)', 'Sum(x+y+z,(x,0,3))',
    'Integral(f(x+y+z),(x,0,1))', 'Integral(f(y+z),(x,0,1))',
    'Integral(x+y+z,(x,0,y+z))', 'Integral(x+y+z,x)',
])
@pytest.mark.parametrize('generator', ['y+z', 'y*z', 'f(y+z)', 'sin(y+z)', '(y+z)**2'])
def test_partial_generators_match_only_free_candidate_nodes(source, generator):
    def make(module, text):
        namespace = vars(module).copy()
        namespace.update(zip('xyz', module.symbols('x y z')))
        namespace['f'] = module.Function('f')
        return eval(text, namespace)
    actual, expected = make(sx, source), make(sp, source)
    for method in ('is_polynomial', 'is_rational_function', 'is_algebraic_expr'):
        assert getattr(actual, method)(make(sx, generator)) is getattr(expected, method)(make(sp, generator)), (source, generator, method)


def atom_data(module):
    x, y = module.symbols('x y')
    f = module.Function('f')
    values = [1+x+2*module.sin(y+module.I*module.pi), f(x)+module.cos(y),
        module.Derivative(f(x), (x, 3)), module.Integral(x+y, (x, 0, 1)),
        module.Lambda(x, x+y), module.UnevaluatedExpr(2*x),
        module.Mul(2, x, x, evaluate=False), module.Add(2, x, 3, evaluate=False),
        module.Tuple(x, y, module.Tuple(1, module.Rational(1, 2))),
        module.Piecewise((x, x>0), (y, True)), module.Max(x, y),
        f(), module.nan, module.true]
    filters = [(), (module.Symbol,), (module.Number,), (module.Integer,),
        (module.Rational,), (module.Mul,), (module.Add,), (module.Pow,),
        (module.Function,), (module.sin,), (f,), (x,),
        (module.Integer(1),), (module.Integer(2),), (module.Rational(1, 2),),
        (module.Symbol, module.Number), (1,), ('x',), (module.NumberSymbol,)]
    return values, filters


@pytest.mark.parametrize('index', range(14))
@pytest.mark.parametrize('selection', range(19))
def test_atoms_support_classes_heads_implicit_classes_and_binders(index, selection):
    values, filters = atom_data(sx)
    reference, ref_filters = atom_data(sp)
    actual = values[index].atoms(*filters[selection])
    assert {value.to_sympy() for value in actual} == reference[index].atoms(*ref_filters[selection])


def test_atom_filtering_preserves_distinct_dummy_identities_and_assumptions():
    x = sx.Dummy('x', positive=True)
    y = sx.Dummy('x', positive=True)
    z = sx.Symbol('x', real=True)
    expression = sx.Function('f')(x, y, z)
    assert expression.atoms(sx.Symbol) == {x, y, z}
    assert expression.atoms(sx.Dummy) == {x, y}
    assert pickle.loads(pickle.dumps(expression)).atoms(sx.Dummy) == {x, y}


def test_number_symbols_distinguish_exact_constants_from_numeric_atoms():
    assert sx.NumberSymbol.class_key() == sp.NumberSymbol.class_key()
    for value in (sx.pi, sx.E, sx.EulerGamma):
        assert isinstance(value, sx.NumberSymbol)
        assert not isinstance(value, sx.Number)
    for value in (sx.I, sx.oo, sx.zoo, sx.nan, sx.Rational(1, 2), sx.Symbol('pi')):
        assert not isinstance(value, sx.NumberSymbol)


def test_scalar_value_methods_preserve_held_structure_and_complex_semantics():
    x = sx.Symbol('x')
    expression = sx.Mul(2, x, x, evaluate=False)
    assert sx.srepr(expression.copy()) == sx.srepr(expression)
    assert sx.srepr(expression.as_expr(x)) == sx.srepr(expression)
    value = sx.Integer(2)+sx.I*x
    expected = sp.Integer(2)+sp.I*sp.Symbol('x')
    assert value.conjugate().to_sympy() == expected.conjugate()
    assert value.adjoint().to_sympy() == expected.adjoint()
    assert value.transpose() == value


def test_queries_do_not_expand_held_powers_or_shared_graphs_without_sympy():
    subprocess.run([sys.executable, '-c', '''
import sys
sys.modules['sympy']=None
import symox as s
x=s.Symbol('x')
held=s.Pow(2,10**9,evaluate=False)
expression=s.Mul(held,x,evaluate=False)
assert expression.is_polynomial(x) is True
assert expression.is_rational_function(x) is True
assert expression.atoms()=={s.Integer(2),s.Integer(10**9),x}
compact=s._core._restore_arithmetic('Mul',s.Integer(1),[(s.Integer(2),s.Integer(10**9)),(x,s.Integer(1))],False)
assert compact.is_polynomial(x) is True
assert compact.atoms()=={s.Integer(2),s.Integer(10**9),x}
e=x
for _ in range(80):
    e=s.Function('shared')(e,e)
assert e.atoms()=={x}
assert e.is_polynomial(x) is None
assert e.is_rational_function(x) is None
assert e.is_algebraic_expr(x) is None
assert len(e.atoms(s.Function))==80
e=x
for _ in range(80):
    e=s._core._restore_arithmetic('Add',s.Integer(0),[(e,s.Integer(1)),(e,s.Integer(1))],False)
assert e.atoms()=={x}
assert e.is_polynomial(x) is True
assert e.is_rational_function(x) is True
'''], check=True, timeout=15)


def test_typed_matrix_metadata_limits_are_explicit():
    expression = sx.Trace(sx.MatrixSymbol('A', 2, 2))
    for method in ('atoms', 'is_polynomial', 'is_rational_function', 'is_algebraic_expr'):
        with pytest.raises(NotImplementedError, match='matrix parents'):
            getattr(expression, method)()
