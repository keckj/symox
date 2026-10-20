"""Reference comparisons for the native protocols used by OctSymPy's bridge."""
import inspect
import pickle

import pytest
import sympy as sp
import symox as sx
from sympy.functions.elementary.piecewise import ExprCondPair as RefPair
from symox.functions.elementary.piecewise import ExprCondPair
from sympy.logic.boolalg import Boolean as RefBoolean, BooleanFunction as RefBooleanFunction
from symox.logic.boolalg import Boolean, BooleanFunction
from symox.integrals.risch import NonElementaryIntegral


@pytest.mark.parametrize('build', [
    lambda s: s.Symbol('x'),
    lambda s: s.Integer(1),
    lambda s: s.Eq(s.Symbol('x'), 0),
    lambda s: s.And(*s.symbols('x y')),
    lambda s: s.Or(*s.symbols('x y')),
    lambda s: s.Not(s.Symbol('x')),
    lambda s: s.ITE(*s.symbols('x y z')),
    lambda s: s.true,
    lambda s: s.false,
    lambda s: s.Contains(s.Symbol('x'), s.S.Reals),
    lambda s: s.Q.positive,
    lambda s: s.Q.positive(s.Symbol('x')),
    lambda s: s.Predicate('user'),
    lambda s: s.Function('And')(s.Symbol('x')),
    lambda s: s.Function('PositivePredicate')(s.Symbol('x')),
])
def test_boolean_categories_distinguish_gates_and_undefined_names(build):
    actual, expected = build(sx), build(sp)
    assert isinstance(actual, Boolean) == isinstance(expected, RefBoolean)
    assert isinstance(actual, BooleanFunction) == isinstance(expected, RefBooleanFunction)
    assert actual.is_Boolean == expected.is_Boolean


@pytest.mark.parametrize('cls,ref', [(sx.Basic, sp.Basic), (Boolean, RefBoolean),
                                  (BooleanFunction, RefBooleanFunction)])
def test_formal_basic_arguments_reconstruction_and_substitution(cls, ref):
    x, y = sx.symbols('x y')
    a, b = sp.symbols('x y')
    value, expected = cls(x, x/y), ref(a, a/b)
    assert str(value) == str(expected)
    assert value.func(*value.args) == value
    assert value.free_symbols == {x, y}
    assert pickle.loads(pickle.dumps(value)) == value
    assert sx.srepr(value) == sp.srepr(expected)
    for simultaneous in (False, True):
        result = value.subs([(x, y), (y, 2)], simultaneous=simultaneous)
        oracle = expected.subs([(a, b), (b, 2)], simultaneous=simultaneous)
        assert tuple(arg.to_sympy() for arg in result.args) == oracle.args
    assert value.has(sx.Symbol)
    assert value.atoms(sx.Symbol) == {x, y}


@pytest.mark.parametrize('build', [
    lambda s: s.Symbol('x'),
    lambda s: s.Eq(s.Symbol('x'), 0),
    lambda s: s.MatrixSymbol('A', 2, 3),
    lambda s: s.Matrix([[1, 2]]),
    lambda s: s.ImmutableMatrix([[1, 2]]),
])
def test_basic_category_includes_immutable_and_symbolic_matrices(build):
    assert isinstance(build(sx), sx.Basic) == isinstance(build(sp), sp.Basic)


@pytest.mark.parametrize('facts', [{}, {'positive': True}, {'real': True},
    {'real': True, 'finite': True}, {'integer': False}, {'commutative': True}])
def test_original_assumptions_are_not_the_inferred_closure(facts):
    actual = sx.Symbol('original_assumptions', **facts)
    expected = sp.Symbol('original_assumptions', **facts)
    assert actual._assumptions_orig == expected._assumptions_orig
    assert actual.assumptions0 == expected.assumptions0
    detached = actual._assumptions_orig
    detached['real'] = 'mutated outside the kernel'
    assert actual._assumptions_orig == expected._assumptions_orig
    assert actual.is_Matrix is False


@pytest.mark.parametrize('condition', [True, False, 0, 1, 2, None, 'x', []])
def test_public_condition_pair_validates_conditions(condition):
    try:
        expected = RefPair(sp.Symbol('x'), condition)
    except TypeError:
        with pytest.raises(TypeError):
            ExprCondPair(sx.Symbol('x'), condition)
    else:
        actual = ExprCondPair(expr=sx.Symbol('x'), cond=condition)
        assert actual.to_sympy() == expected
        assert actual.cond.to_sympy() == expected.cond


def test_condition_pair_signature_pickle_and_native_structure():
    assert str(inspect.signature(ExprCondPair)) == str(inspect.signature(RefPair))
    x, y = sx.symbols('x y')
    pair = ExprCondPair(x, x > 0)
    assert pair == sx.Piecewise((x, x > 0), (0, True)).args[0]
    assert isinstance(pair, sx.Basic)
    assert pair != (x, x > 0)
    assert hash(pair) == hash(pickle.loads(pickle.dumps(pair)))
    assert pair.subs(x, y) == ExprCondPair(y, y > 0)
    assert pair.xreplace({x: y}) == ExprCondPair(y, y > 0)


@pytest.mark.parametrize('evaluate', [True, False])
def test_boolean_piecewise_conditions_convert_without_losing_index_domains(evaluate):
    i = sx.Symbol('i', integer=True)
    j = sp.Symbol('i', integer=True)
    condition = sx.Piecewise((sx.Eq(i+3, 0), i < 0), (sx.Eq(i, 0), True))
    oracle = sp.Piecewise((sp.Eq(j+3, 0), j < 0), (sp.Eq(j, 0), True))
    pair = ExprCondPair(2, condition)
    assert isinstance(pair.cond, Boolean)
    expected = RefPair(2, oracle)
    for position in range(-8, 8):
        assert pair.cond.subs(i, position).to_sympy() == expected.cond.subs(j, position)
    outer = sx.Piecewise((2, condition), (0, True), evaluate=evaluate)
    assert outer.func(*outer.args, evaluate=evaluate) == outer
    assert condition.rewrite(sx.ITE) == pair.cond
    with pytest.raises(TypeError):
        sx.Piecewise((1, i > 0), (2, True)).rewrite(sx.ITE)


@pytest.mark.parametrize('expr', ['x', [], None])
def test_condition_pair_rejects_non_basic_values_with_reference_error_type(expr):
    with pytest.raises(TypeError) as actual:
        ExprCondPair(expr, True)
    with pytest.raises(TypeError) as expected:
        RefPair(expr, True)
    assert type(actual.value) is type(expected.value)


@pytest.mark.parametrize('value', ['x', [], None, 2, (1, 2)])
def test_basic_conversion_is_strict(value):
    from sympy.core.basic import as_Basic as ref_as_basic
    from symox.core.basic import as_Basic
    try:
        expected = ref_as_basic(value)
    except TypeError:
        with pytest.raises(TypeError):
            as_Basic(value)
    else:
        assert as_Basic(value).to_sympy() == expected


@pytest.mark.parametrize('unicode', [False, True])
@pytest.mark.parametrize('shape', [(0, 0), (0, 3), (3, 0), (1, 1), (1, 3), (3, 1), (2, 3)])
def test_explicit_matrix_pretty_and_constructor_representation(unicode, shape):
    def build(s):
        x = s.Symbol('x')
        values = [x, s.Rational(1, 2), s.sqrt(x), x**2, s.Rational(111, 22), -x]
        return s.Matrix(*shape, (values*3)[:shape[0]*shape[1]])
    actual, expected = build(sx), build(sp)
    assert sx.pretty(actual, use_unicode=unicode) == sp.pretty(expected, use_unicode=unicode)
    assert sx.srepr(actual) == sp.srepr(expected)
    assert sx.srepr(actual.as_immutable()) == sp.srepr(expected.as_immutable())
    assert sx.sympify(sx.srepr(actual)) == actual


@pytest.mark.parametrize('kind', ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix'])
def test_global_matrix_integration_preserves_exact_entries_and_freezes_output(kind):
    def build(s):
        x = s.Symbol('x')
        matrix = getattr(s, kind)([[x, x*x], [0, s.sin(x)]])
        return s.integrate(matrix, (x, 0, 1))
    actual, expected = build(sx), build(sp)
    assert actual.to_sympy() == expected
    assert isinstance(actual, sx.Basic) == isinstance(expected, sp.Basic)
    assert actual._immutable is True


@pytest.mark.parametrize('value', [0., 1., 1.2, 1e30, 1e-30])
def test_matrix_pretty_uses_nested_float_format(value):
    assert sx.pretty(sx.Matrix([[sx.Float(value)]])) == sp.pretty(sp.Matrix([[sp.Float(value)]]))


@pytest.mark.parametrize('flags', [{}, {'dict': True}, {'set': True}, {'check': False}])
@pytest.mark.parametrize('case', ['unrelated', 'contradiction', 'identity', 'empty'])
def test_solve_assignment_output_and_boolean_constraints(case, flags):
    def run(s):
        x, y, z = s.symbols('x y z')
        equations = {'unrelated': [s.Eq(x, y)], 'contradiction': [s.Eq(x, y), False],
                     'identity': [True, 0], 'empty': []}[case]
        return s.solve(equations, z, **flags)
    actual, expected = run(sx), run(sp)
    if flags.get('set'):
        assert [x.to_sympy() for x in actual[0]] == expected[0]
        assert actual[1] == expected[1] == set()
    else:
        assert actual == expected == []


def test_matrix_substitution_simultaneous_and_sequential_are_distinct():
    def build(s):
        n, m = s.symbols('n m', positive=True, integer=True)
        return s.MatrixSymbol('A', n, m), [(n, m), (m, 2)]
    actual, rules = build(sx)
    expected, oracle_rules = build(sp)
    for simultaneous in (False, True):
        result = actual.subs(rules, simultaneous=simultaneous)
        assert result.to_sympy() == expected.subs(oracle_rules, simultaneous=simultaneous)
    x, y = sx.symbols('x y')
    base = sx.MatrixSymbol('B', 2, 2)
    assert (x*y*base).subs([(x, y), (y, x)], simultaneous=True) == x*y*base
    assert actual.subs(rules, unknown_option=True).to_sympy() == expected.subs(oracle_rules, unknown_option=True)


@pytest.mark.parametrize('method', ['hstack', 'vstack'])
def test_matrix_base_stacking_uses_real_join_and_abstract_empty_constructor(method):
    a, b = sx.Matrix([[1, 2]]), sx.Matrix([[3, 4]])
    ra, rb = sp.Matrix([[1, 2]]), sp.Matrix([[3, 4]])
    assert getattr(sx.MatrixBase, method)(a, b).to_sympy() == getattr(sp.MatrixBase, method)(ra, rb)
    assert getattr(sx.MatrixBase, method)(a) is a
    with pytest.raises(NotImplementedError):
        getattr(sx.MatrixBase, method)()


def test_explicit_non_elementary_marker_retains_type_in_pickle():
    x = sx.Symbol('x')
    value = NonElementaryIntegral(sx.exp(sx.exp(x)), x)
    assert isinstance(value, sx.Integral)
    assert isinstance(pickle.loads(pickle.dumps(value)), NonElementaryIntegral)
    assert value.diff(x) == sx.exp(sx.exp(x))
    # An ordinary unresolved integral is not a proof of non-elementarity.
    assert not isinstance(sx.Integral(sx.Function('f')(x), x), NonElementaryIntegral)


def test_non_elementary_marker_structural_operations_and_explicit_interop():
    from sympy.integrals.risch import NonElementaryIntegral as RefIntegral
    x, y = sx.symbols('x y')
    value = NonElementaryIntegral(sx.Function('f')(x), x)
    expected = RefIntegral(sp.Function('f')(sp.Symbol('x')), sp.Symbol('x'))
    assert value.to_sympy() == expected
    assert sx.srepr(value) == sp.srepr(expected)
    for result in (value.subs(x, y), value.xreplace({x: y}), value.doit(),
                   value.simplify(), value.expand(), value.copy()):
        assert isinstance(result, NonElementaryIntegral)
