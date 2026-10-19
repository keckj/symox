"""Differential checks for scalar protocols consumed by external packages."""
import inspect
import itertools
import operator
import numbers
import pickle
from fractions import Fraction
from collections.abc import Iterable
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx

from symox.core import logic as logic
from sympy.core import logic as reference_logic
from symox.core import relational as relational
from sympy.core import relational as reference_relational


@pytest.mark.parametrize('numerator,denominator', [(0, 1), (2, 1), (-2, 3), (10**100+1, 10**80+7)])
def test_exact_scalars_support_standard_numeric_abc_and_fraction_protocol(numerator, denominator):
    native, reference = sx.Rational(numerator, denominator), sp.Rational(numerator, denominator)
    for abc in [numbers.Number, numbers.Complex, numbers.Real, numbers.Rational, numbers.Integral]:
        assert isinstance(native, abc) == isinstance(reference, abc)
    assert Fraction(native) == Fraction(reference)
    assert Fraction(native, sx.Integer(7)) == Fraction(reference, sp.Integer(7))
    assert type(native.numerator) is type(reference.numerator) is int
    assert type(native.denominator) is type(reference.denominator) is int
    assert pickle.loads(pickle.dumps(native)) == native


def test_numeric_abc_caches_do_not_classify_symbolic_or_inexact_nodes_as_rationals():
    for value in [sx.Symbol('x'), sx.sqrt(2), sx.Float(1.5), sx.pi, sx.oo, sx.nan]:
        assert not isinstance(value, numbers.Rational)
    # All these objects share a physical native class. ABC negative caches
    # must still honor the logical numeric class for each exact scalar.
    assert isinstance(sx.Integer(3), numbers.Integral)
    assert isinstance(sx.Rational(2, 3), numbers.Rational)
    assert not isinstance(sx.Rational(2, 3), numbers.Integral)


def test_symbol_construction_cache_normalizes_assumptions_and_can_be_cleared():
    from symox.core.cache import clear_cache
    x = sx.Symbol('cached_downstream_symbol', positive=True)
    assert sx.Symbol('cached_downstream_symbol', positive=1, real=None) is x
    clear_cache()
    replacement = sx.Symbol('cached_downstream_symbol', positive=True)
    assert replacement == x and replacement is not x
    assert sx.Symbol('unhashable_assumption', negative=[]) == sx.Symbol('unhashable_assumption', negative=False)


def test_scalar_identity_operations_preserve_the_surviving_python_operand():
    x = sx.Symbol('retained_dimension', positive=True)
    assert x*1 is x
    assert 1*x is x
    assert x+0 is x
    assert 0+x is x
    assert x/1 is x
    assert x**1 is x
    # Reconstructing a dimensional singleton through exact cancellation
    # must retain the supplied object even after constructor-cache eviction.
    old = sx.Symbol('old_dimension', positive=True)
    for index in range(1100): sx.Symbol('temporary_'+str(index))
    other = sx.Symbol('old_dimension', positive=True)
    assert other == old and other is not old
    assert (other/old)*old is old


@pytest.mark.parametrize('name', ['__add__','__radd__','__sub__','__rsub__','__mul__',
    '__rmul__','__truediv__','__rtruediv__','__pow__','__rpow__'])
def test_higher_priority_operands_keep_their_own_arithmetic_protocol(name):
    calls=[]
    class External:
        _op_priority=1000
        def _sympy_(self):
            raise AssertionError('Higher-priority operands must not be coerced')
    reflected='__'+name[3:] if name.startswith('__r') else '__r'+name[2:]
    def handler(self,other):
        calls.append(other)
        return 'external result'
    setattr(External,reflected,handler)
    x=sx.Symbol('x')
    assert getattr(x,name)(External()) == 'external result'
    assert calls == [x]


def test_tuple_results_keep_actual_sequence_slots_after_native_transformations():
    x=sx.Symbol('x')
    value=sx.Tuple(1,x,sx.Tuple(2,3))
    for result in [value,value.subs(x,2),value.xreplace({x:2}),value.evalf(15),
                   value.doit(),sx.nsimplify(value),pickle.loads(pickle.dumps(value))]:
        assert type(result) is sx.Tuple
        assert list(result) == list(result.args)
        assert len(result) == 3
        assert list(result[-1]) == [2,3]


@pytest.mark.parametrize('constructor', [lambda m: m.Symbol('x'), lambda m: m.Rational(2,3),
    lambda m: m.sin(m.Symbol('x')), lambda m: m.Tuple(), lambda m: m.Tuple(1,2,3)])
def test_scalar_and_tuple_iterable_type_protocols_match_reference(constructor):
    native, reference = constructor(sx), constructor(sp)
    assert isinstance(native, Iterable) == isinstance(reference, Iterable)
    if isinstance(reference, Iterable):
        assert list(native) == list(reference)
        assert list(reversed(native)) == list(reversed(reference))
        assert len(native) == len(reference)
        with pytest.raises(IndexError): native[10**100]
        with pytest.raises(IndexError): native[-10**100]
    else:
        with pytest.raises(TypeError): iter(native)


def test_tuple_iteration_does_not_rebuild_all_arguments_for_each_element(monkeypatch):
    value = sx.Tuple(*range(1000))
    def forbidden(self):
        raise AssertionError('Sequence access must not rebuild the complete args tuple')
    monkeypatch.setattr(sx.Expr,'args',property(forbidden))
    assert len(value) == 1000
    # __iter__ materializes the args tuple once; random access and length
    # must remain constant-time even with that property unavailable.
    assert value[-1] == 999


def test_quantity_structural_queries_and_factor_protocol():
    from symox.physics import units as u
    from sympy.physics import units as r
    for native, reference in [(u.volt,r.volt), (u.second,r.second), (u.Quantity,r.Quantity),
        (sx.Symbol,sp.Symbol), (u.volt.name,r.volt.name),
        (sx.Symbol('volt',positive=True),sp.Symbol('volt',positive=True))]:
        assert u.volt.has(native) == r.volt.has(reference)
    assert u.volt.as_ordered_factors() == [u.volt]
    assert u.volt.replace(u.volt,2*u.second) == 2*u.second
    assert u.volt.replace(u.second,3) is u.volt


@pytest.mark.parametrize('constructor', ['Matrix','ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
def test_exponential_of_concrete_matrix_uses_matrix_function_and_immutable_output(constructor):
    native = sx.exp(getattr(sx,constructor)([[0,1],[0,0]]))
    reference = sp.exp(getattr(sp,constructor)([[0,1],[0,0]]))
    assert native.tolist() == reference.tolist() == [[1,1],[0,1]]
    assert native._immutable
    assert isinstance(native,sx.ImmutableSparseMatrix) == isinstance(reference,sp.ImmutableSparseMatrix)


@pytest.mark.parametrize('short,long', [('Eq','Equality'), ('Ne','Unequality'),
    ('Ge','GreaterThan'), ('Le','LessThan'), ('Gt','StrictGreaterThan'), ('Lt','StrictLessThan')])
def test_relational_aliases_share_constructor_and_instance_protocols(short, long):
    constructor = getattr(sx, short)
    assert constructor is getattr(sx, long)
    x = sx.Symbol('x')
    relation = constructor(x, 1, evaluate=False)
    assert isinstance(relation, constructor)
    assert isinstance(relation, sx.Relational)
    assert relation.func(*relation.args, evaluate=False) == relation
    assert not isinstance(sx.Integer(1), constructor)
    assert constructor(1, 2) == getattr(sp, short)(1, 2)
    assert inspect.signature(constructor) == inspect.signature(getattr(sp, short))
    reference = getattr(sp, short)(sp.Symbol('x'), 1, evaluate=False)
    assert relation.count(sx.Expr) == reference.count(sp.Expr)
    assert str(sx.count_ops(relation,visual=True)) == str(sp.count_ops(reference,visual=True))


def test_matrix_subclasses_keep_python_construction_and_overrides():
    class LabeledMatrix(sx.Matrix):
        def __new__(cls, *args):
            return super().__new__(cls, *args)
        def label(self):
            return self.description
    matrix = LabeledMatrix([[sx.Rational(2, 3), 1], [0, 3]])
    matrix.description = 'native matrix'
    assert type(matrix) is LabeledMatrix
    assert matrix.label() == 'native matrix'
    assert matrix.det() == 2
    assert matrix[0, 0] == sx.Rational(2, 3)


@pytest.mark.parametrize('name', ['exp','sin','cos','tan','sinh','cosh','tanh','gamma',
    'Abs','DiracDelta','Heaviside','sign','re','im','conjugate'])
def test_builtin_function_heads_have_public_identity_for_replacement_queries(name):
    head = getattr(sx, name)
    x = sx.Symbol('x')
    expression = head(x, evaluate=False)
    assert expression.func is head
    assert str(head) == str(getattr(sp,name))
    assert pickle.loads(pickle.dumps(head)) is head
    assert expression.replace(lambda node: node.is_Function and node.func is head,
                              lambda node: sx.Integer(3)) == 3
    assert expression.func(*expression.args, evaluate=False) == expression
    # An undefined function with a builtin spelling remains a distinct head.
    undefined = sx.Function(name)(x)
    assert undefined.func is not head


@pytest.mark.parametrize('name', ['factor_list', 'sqf_list'])
def test_factor_conversion_failure_after_all_generators_cancel(name):
    def expression(module):
        d = module.Symbol('d')
        return 4*d*(d-1)-9*(d-1)*(d+5)+5*(d+3)*(d+5)
    for options in [{}, {'domain':'QQ'}, {'polys':True}, {'gens':()}, {'extension':True}]:
        with pytest.raises(sp.PolificationFailed) as reference:
            getattr(sp, name)(expression(sp), **options)
        with pytest.raises(sx.PolificationFailed) as native:
            getattr(sx, name)(expression(sx), **options)
        assert native.value.expr == reference.value.expr == 120
    d = sx.Symbol('d')
    assert getattr(sx, name)(expression(sx), d) == (120, [])
    with pytest.raises(sx.OptionError):
        getattr(sx, name)(expression(sx), bogus=True)


@pytest.mark.parametrize('name', ['factor', 'sqf'])
def test_scalar_factor_fallback_evaluates_complete_polynomial_cancellation(name):
    d = sx.Symbol('d')
    constant = (d-2)*(d-1)-2*(d-2)*(d+2)+(d+1)*(d+2)
    assert getattr(sx, name)(constant) == 12
    assert getattr(sx, name)(constant-12) == 0
    assert getattr(sx, name)(constant**2) == 144
    assert getattr(sx, name)(constant**sx.Rational(1, 2)) == 2*sx.sqrt(3)


@pytest.mark.parametrize('value,attribute', [(0,'Zero'), (1,'One'), (-1,'NegativeOne')])
def test_canonical_integer_identity_through_construction_and_scalar_arithmetic(value, attribute):
    canonical = getattr(sx.S, attribute)
    x = sx.Symbol('x')
    results = [sx.Integer(value), sx.Rational(value), sx.sympify(value),
               x-x+value, (x-x+value).expand(), sx.Integer(value)*sx.Integer(1),
               sx.Integer(value)/sx.Integer(1), sx.Integer(value)**1,
               +sx.Integer(value), -sx.Integer(-value)]
    assert all(result is canonical for result in results)
    assert sx.sympify(x) is x
    assert +x is x
    assert sx.Rational(1, 2) is sx.S.Half
    assert (sx.Integer(1)/2) is sx.S.Half


@pytest.mark.parametrize('arguments', [(1,1), (2,sx.Rational(1,2)), (1,1,1)])
@pytest.mark.parametrize('exponent', [-2,-1,2])
def test_held_unit_product_powers_trigger_collection_when_multiplied(arguments, exponent):
    native = sx.Mul(*arguments, evaluate=False)**exponent
    reference = sp.Mul(*(sp.sympify(str(a)) for a in arguments), evaluate=False)**exponent
    assert native.is_Pow == reference.is_Pow
    assert native.args == reference.args
    assert sx.Mul(0,1,evaluate=False)*native == sp.Mul(0,1,evaluate=False)*reference == 0
    # A literal identity remains a separate operation with held semantics.
    assert (sx.Mul(0,1,evaluate=False)*1).is_Mul


@pytest.mark.parametrize('constructor', ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix'])
@pytest.mark.parametrize('operation', ['add', 'sub', 'mul', 'matmul'])
@pytest.mark.parametrize('reflected', [False, True])
def test_two_dimensional_numpy_arithmetic_keeps_matrix_semantics(constructor, operation, reflected):
    a = np.array([[5, 6], [7, 8]], dtype=object)
    native = getattr(sx, constructor)([[1, 2], [3, 4]])
    reference = getattr(sp, constructor)([[1, 2], [3, 4]])
    op = getattr(operator, operation)
    actual = op(a, native) if reflected else op(native, a)
    expected = op(a, reference) if reflected else op(reference, a)
    assert isinstance(actual, sx.MatrixBase)
    assert actual.shape == expected.shape
    assert actual.tolist() == expected.tolist()
    assert actual._immutable == isinstance(expected, (sp.ImmutableMatrix, sp.ImmutableSparseMatrix))


def test_sfepy_nodal_constraint_keeps_a_symbolic_system_until_solve():
    x, y = sx.symbols('u0 u1')
    rx, ry = sp.symbols('u0 u1')
    rhs = np.array([[0.]])
    actual = sx.Matrix([[1., -1.]])*sx.Matrix([x, y])-rhs
    expected = sp.Matrix([[1., -1.]])*sp.Matrix([rx, ry])-rhs
    assert isinstance(actual, sx.MatrixBase)
    result, reference = sx.solve(actual), sp.solve(expected)
    assert {str(k) for k in result} == {str(k) for k in reference}
    assert sx.simplify(result[x]-y) == 0
    assert float(result[x].subs(y, 3)) == float(reference[rx].subs(ry, 3))


def test_external_scalar_sympy_conversion_protocol_is_consumed_natively():
    class ExternalScalar:
        def _sympy_(self):
            return sx.Symbol('a')+sx.Rational(1, 3)
    value = ExternalScalar()
    assert sx.sympify(value) == sx.Symbol('a')+sx.Rational(1, 3)
    assert sx.expand(sx.Symbol('x')*value) == sx.Symbol('a')*sx.Symbol('x')+sx.Symbol('x')/3


def test_external_reference_result_keeps_exact_arithmetic():
    class ExternalScalar:
        def _sympy_(self):
            return sp.Rational(123456789012345678901, 7)
    assert sx.sympify(ExternalScalar()) == sx.Rational(123456789012345678901, 7)


def test_cyclic_external_conversion_raises_instead_of_overflowing_native_stack():
    class ExternalScalar:
        def _sympy_(self): return self.other
    a, b = ExternalScalar(), ExternalScalar()
    a.other, b.other = b, a
    with pytest.raises(RecursionError):
        sx.sympify(a)
    with pytest.raises(TypeError, match='returned its input'):
        a.other = a
        sx.sympify(a)


def test_long_external_conversion_chains_use_bounded_iteration():
    class Link:
        def __init__(self, tail): self.tail = tail
        def _sympy_(self): return self.tail
    value = sx.Rational(7, 3)
    for _ in range(200): value = Link(value)
    assert sx.sympify(value) == sx.Rational(7, 3)
    class Endless:
        def _sympy_(self): return Endless()
    with pytest.raises(RecursionError): sx.sympify(Endless())
    class Reentrant:
        def _sympy_(self): return sx.sympify(self)
    with pytest.raises(RecursionError): sx.sympify(Reentrant())


def test_native_matrix_conversion_does_not_call_its_reference_exporter(monkeypatch):
    def forbidden(self):
        raise AssertionError('Native matrix processing must not export to SymPy')
    monkeypatch.setattr(sx.Matrix, '_sympy_', forbidden)
    matrix = sx.Matrix([[sx.Float('0.5'),2],[3,4]])
    result = sx.nsimplify(matrix)
    assert result.tolist() == [[sx.Rational(1,2),2],[3,4]]
    assert result._immutable
    with pytest.raises(TypeError): sx._core.nsimplify(matrix)


@pytest.mark.parametrize('sequence', [(1,2), [1,2], 'ab', b'ab'])
@pytest.mark.parametrize('count', [-2,0,1,3])
def test_exact_integer_sequence_repetition_preserves_builtin_dispatch(sequence, count):
    assert sequence*sx.Integer(count) == sequence*sp.Integer(count)
    assert sx.Integer(count)*sequence == sp.Integer(count)*sequence


@pytest.mark.parametrize('numerator,denominator', [(1,2),(-1,2),(355,113),
    (10**80+1,10**70+3),(-10**80-1,10**70+3)])
def test_rational_bounded_denominator_uses_exact_reference_tie_breaking(numerator, denominator):
    value, reference = sx.Rational(numerator,denominator), sp.Rational(numerator,denominator)
    assert inspect.signature(value.limit_denominator) == inspect.signature(reference.limit_denominator)
    for maximum in [1,2,3,7,10.9,100,1000000]:
        assert value.limit_denominator(maximum) == reference.limit_denominator(maximum)
    for maximum in [0,-1]:
        with pytest.raises(ValueError): value.limit_denominator(maximum)


@pytest.mark.parametrize('operation', ['add', 'sub', 'mul', 'matmul'])
def test_numpy_matrix_shape_mismatch_is_not_broadcast(operation):
    op = getattr(operator, operation)
    with pytest.raises(sp.ShapeError):
        op(sp.eye(2), np.array([[3]]))
    with pytest.raises(sx.ShapeError):
        op(sx.eye(2), np.array([[3]]))


def test_root_import_exposes_native_core_types_without_extra_imports():
    process = subprocess.run([sys.executable, '-c',
        'import symox; assert symox.core.expr.Expr is symox.Expr; '
        'assert symox.core.Symbol is symox.Symbol'], capture_output=True, text=True, timeout=30)
    assert process.returncode == 0, process.stderr


@pytest.mark.parametrize('name', ['fuzzy_and', 'fuzzy_or', 'fuzzy_xor', 'fuzzy_nand', '_fuzzy_group', '_torf'])
def test_fuzzy_truth_tables_and_signature(name):
    native, reference = getattr(logic, name), getattr(reference_logic, name)
    assert inspect.signature(native) == inspect.signature(reference)
    for size in range(5):
        for values in itertools.product([True, False, None, 0, 1, 2], repeat=size):
            assert native(iter(values)) is reference(iter(values)), (name, values)
            if name == '_fuzzy_group':
                assert native(iter(values), quick_exit=True) is reference(iter(values), quick_exit=True)


def test_fuzzy_unary_signatures_and_non_boolean_values():
    for name in ['fuzzy_bool', 'fuzzy_not']:
        native, reference = getattr(logic, name), getattr(reference_logic, name)
        assert inspect.signature(native) == inspect.signature(reference)
        for value in [None, True, False, 0, 1, 2, '', 'text']:
            assert native(value) is reference(value)


@pytest.mark.parametrize('name,value,expected', [('fuzzy_and', False, False), ('fuzzy_or', True, True)])
def test_fuzzy_short_circuit_does_not_consume_later_queries(name, value, expected):
    def arguments():
        yield value
        raise AssertionError('The rest of a lazy query must not be evaluated')
    assert getattr(logic, name)(arguments()) is expected


def operands(module):
    x = module.Symbol('x')
    r = module.Symbol('r', real=True)
    p = module.Symbol('p', positive=True)
    return [module.Integer(0), module.Integer(1), module.Integer(-1), module.Rational(1, 3),
            module.sqrt(2), module.pi, module.I, module.oo, -module.oo, module.nan,
            x, x+1, r, r+1, p]


@pytest.mark.parametrize('name', ['is_ge', 'is_le', 'is_gt', 'is_lt'])
def test_fuzzy_relations_match_real_complex_unknown_and_infinite_domains(name):
    native, reference = getattr(relational, name), getattr(reference_relational, name)
    assert inspect.signature(native) == inspect.signature(reference)
    actual, expected = operands(sx), operands(sp)
    for i,j in itertools.product(range(len(actual)), repeat=2):
        left, right = native(actual[i], actual[j]), reference(expected[i], expected[j])
        assert (None if left is None else bool(left)) is (None if right is None else bool(right)), (name, actual[i], actual[j])
        assert isinstance(left, bool) == isinstance(right, bool), (name, actual[i], actual[j])
    for args in [(1, sx.Integer(0)), (sx.Integer(0), 1)]:
        with pytest.raises(TypeError): native(*args)


@pytest.mark.parametrize('name', ['is_ge', 'is_le', 'is_gt', 'is_lt'])
@pytest.mark.parametrize('predicate', ['positive', 'negative', 'nonnegative', 'zero', 'real'])
def test_fuzzy_relations_accept_explicit_assumptions(name, predicate):
    x, rx = sx.Symbol('x'), sp.Symbol('x')
    actual = getattr(relational, name)(x, sx.Integer(0), assumptions=getattr(sx.Q, predicate)(x))
    expected = getattr(reference_relational, name)(rx, sp.Integer(0), assumptions=getattr(sp.Q, predicate)(rx))
    assert actual is expected
