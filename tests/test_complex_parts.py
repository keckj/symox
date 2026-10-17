"""Bounded native complex expansion, including principal branches and storage."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def environment(module):
    result = {name: getattr(module, name) for name in [
        'I', 'pi', 'E', 'oo', 'zoo', 'nan', 'Rational', 'Float', 'sqrt',
        'sin', 'cos', 'tan', 'cot', 'sinh', 'cosh', 'tanh', 'coth', 'exp',
        'log', 'Abs', 'arg', 'sign', 're', 'im', 'conjugate', 'UnevaluatedExpr',
    ]}
    for name, flags in [('x', {}), ('y', {}), ('r', {'real': True}),
                        ('t', {'real': True}), ('p', {'positive': True}),
                        ('n', {'negative': True}), ('j', {'imaginary': True}),
                        ('h', {'extended_real': True})]:
        result[name] = module.Symbol(name, **flags)
    result.update(f=module.Function('f'), fr=module.Function('fr', real=True),
                  fi=module.Function('fi', imaginary=True),
                  fe=module.Function('fe', extended_real=True),
                  user_sin=module.Function('sin'))
    return result


SOURCES = [
    '0', '1', '-1', 'Rational(2,3)', 'Float("1.25",80)',
    'I', 'pi', 'E', 'oo', '-oo', 'zoo', 'nan',
    'x', 'r', 'p', 'n', 'j', 'h', 'f(x)', 'fr(x)', 'fi(x)', 'fe(x)',
    'user_sin(x)', 'UnevaluatedExpr(r+I*t)', 're(x)', 'im(x)', 'Abs(x)', 'arg(x)', 'sign(x)',
    'x+I*y', 'r+I*t', '3*x-4*I*y+Rational(2,3)',
    '2*x*y', 'I*x*y', 'r*x*y', 'j*x', 'x*conjugate(x)',
    'x*y*conjugate(x)', '(x+y)*(x-y)', 'x*(y+I)', '(r+I*t)*(x+I*y)',
    'conjugate(x)', 'conjugate(x+I*y)', 'conjugate((x+y)*(x-y))',
    'x**2', 'x**3', 'x**-1', 'x**-2', '(r+I*t)**2', '(r+I*t)**-2',
    '(x+I*y)**3', '(sin(x)+I*cos(y))**2',
    'sqrt(x)', 'sqrt(r)', 'sqrt(p)', 'sqrt(n)', 'sqrt(j)',
    'x**Rational(1,3)', 'n**Rational(1,3)', 'p**Rational(2,3)',
    '(r+I*t)**Rational(-3,2)', 'x**y', '(x+y)**t', 'E**(x+I*y)',
    'log(x)', 'log(r)', 'log(p)', 'log(n)', 'log(x+I*y)', 'log(r+I*t)',
    'log(1+I)', 'log(2*x)', 'log(conjugate(x))', 'log(0)', 'log(-0)',
    'log(I*oo)', 'exp(I*oo)', 'arg(0)', 'arg(zoo)',
    '(x+I*y)**-1', '(x+I*y)**-3', '(x+y)**-3',
]
SOURCES += [f'{function}({argument})' for function in
            ['exp', 'sin', 'cos', 'tan', 'cot', 'sinh', 'cosh', 'tanh', 'coth']
            for argument in ['x', 'r', 'r+I*t', '(x+I*y)**2']]


@pytest.mark.parametrize('deep', [False, True])
@pytest.mark.parametrize('source', SOURCES)
def test_complex_parts_match_primary_sympy_trees(source, deep):
    actual = eval(source, {'__builtins__': {}}, environment(sx))
    expected = eval(source, {'__builtins__': {}}, environment(sp))
    actual = sx.sympify(actual)
    expected = sp.sympify(expected)
    parts = actual.as_real_imag(deep=deep)
    oracle = expected.as_real_imag(deep=deep)
    assert isinstance(parts, tuple) and len(parts) == 2
    assert tuple(map(sx.to_sympy, parts)) == oracle, (source, deep, parts, oracle)
    for value in parts:
        assert pickle.loads(pickle.dumps(value)) == value
        assert sx.sympify(value.srepr()) == value


@pytest.mark.parametrize('source', ['x', 'r', 'f(x)', 'x+y', 'x*y', 'r*x',
    '(x+y)*(x-y)', 'x**2', 'x**Rational(1,3)', 'x**y', 'E**x', 'sin(x)', 'exp(x)',
    '0', '1', 'I', 'pi', 'oo', 'zoo', 'nan', 're(x)', 'im(x)', 'Abs(x)', 'arg(x)', 'conjugate(x)'])
@pytest.mark.parametrize('ignored', ['self', 'x', 'one', 'absent'])
def test_ignore_follows_the_dispatch_contract_without_affecting_children(source, ignored):
    native, oracle = environment(sx), environment(sp)
    actual = sx.sympify(eval(source, {'__builtins__': {}}, native))
    expected = sp.sympify(eval(source, {'__builtins__': {}}, oracle))
    ignore = {'self': actual, 'x': native['x'], 'one': sx.Integer(1), 'absent': native['y']}[ignored]
    reference = {'self': expected, 'x': oracle['x'], 'one': sp.Integer(1), 'absent': oracle['y']}[ignored]
    result = actual.as_real_imag(ignore=ignore)
    target = expected.as_real_imag(ignore=reference)
    assert (result is None) == (target is None)
    if result is not None:
        assert tuple(map(sx.to_sympy, result)) == target


def test_projection_evaluators_use_as_real_imag_like_sympy():
    z = sx.Symbol('z')
    value = sx.sin(z)
    zs = sp.Symbol('z')
    assert sx.to_sympy(sx.re(value)) == sp.re(sp.sin(zs))
    assert sx.to_sympy(sx.im(value)) == sp.im(sp.sin(zs))
    a, b = value.as_real_imag()
    assert a == sx.sin(sx.re(z))*sx.cosh(sx.im(z))
    assert b == sx.cos(sx.re(z))*sx.sinh(sx.im(z))
    value = sx.exp((z+sx.I*sx.Symbol('w'))**2)
    assert value.as_real_imag(deep=True) != value.as_real_imag(deep=False)


@pytest.mark.parametrize('degree', list(range(1, 9)))
def test_integer_power_parts_have_independent_exact_reconstruction_certificates(degree):
    r, t = sx.symbols('r t', real=True)
    value = (r+sx.I*t)**degree
    a, b = value.as_real_imag()
    assert sx.expand(a+sx.I*b-value) == 0
    inverse = (r+sx.I*t)**-degree
    a, b = inverse.as_real_imag()
    # Away from r=t=0, (r+i*t)^(-n) = (r-i*t)^n/(r**2+t**2)^n.
    denominator = (r*r+t*t)**degree
    assert sx.cancel(sx.expand((a+sx.I*b)*denominator-(r-sx.I*t)**degree)) == 0


@pytest.mark.parametrize('kind', ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix'])
@pytest.mark.parametrize('deep', [False, True])
def test_matrix_parts_preserve_storage_shape_mutability_and_pickle(kind, deep):
    x, y = sx.symbols('x y')
    X, Y = sp.symbols('x y')
    actual = getattr(sx, kind)([[x+sx.I*y, 0], [sx.sin(x), (x+sx.I*y)**2]])
    expected = getattr(sp, kind)([[X+sp.I*Y, 0], [sp.sin(X), (X+sp.I*Y)**2]])
    result, target = actual.as_real_imag(deep=deep), expected.as_real_imag(deep=deep)
    for matrix, reference in zip(result, target):
        assert matrix.shape == reference.shape
        assert matrix.to_sympy() == reference
        assert bool(matrix._immutable) == ('Immutable' in kind)
        assert isinstance(matrix, getattr(sx, kind))
        assert bool(getattr(matrix, '_symox_sparse', False)) == ('Sparse' in kind)
        assert pickle.loads(pickle.dumps(matrix)) == matrix
    if 'Immutable' in kind:
        with pytest.raises(TypeError):
            result[0][0, 0] = 0
    else:
        result[0][0, 0] = 0
    assert actual[0, 0] == x+sx.I*y


def test_sparse_parts_work_on_stored_entries_without_dense_allocation():
    x = sx.Symbol('x')
    value = sx.ImmutableSparseMatrix(10**7, 10**7, {(1, 2): x+sx.I})
    real, imag = value.as_real_imag()
    assert real.shape == imag.shape == value.shape
    assert real.todok() == {(1, 2): sx.re(x)}
    assert imag.todok() == {(1, 2): sx.im(x)+1}


@pytest.mark.parametrize('constructor', [lambda: sx.Symbol('x'), lambda: sx.eye(1),
    lambda: sx.ImmutableSparseMatrix(2, 2, {})])
def test_unimplemented_hints_are_rejected_explicitly(constructor):
    value = constructor()
    for hints in [{'log': True}, {'trig': False}, {'force': True}, {'unknown': 0}]:
        with pytest.raises(NotImplementedError, match='hint'):
            value.as_real_imag(**hints)
    if getattr(value, 'is_Matrix', False):
        with pytest.raises(NotImplementedError, match='ignore'):
            value.as_real_imag(ignore=sx.Symbol('x'))


def test_compact_numeric_powers_and_shared_graphs_obey_preallocation_bounds_without_sympy():
    source = '''
import pickle, sys
sys.modules['sympy'] = None
import symox as s
x=s.Symbol('x')
pair=s.Function('pair',nargs=2)
e=x
for _ in range(80): e=pair(e,e)
a,b=e.as_real_imag()
assert a.args==(e,) and b.args==(e,)
assert len(a.__reduce__()[1][1])==82
assert pickle.loads(pickle.dumps((a,b)))==(a,b)
for power in [s.Pow(2,10**9,evaluate=False),
              s.Pow(x+s.I,10**9,evaluate=False),
              s.Pow(2,s.Rational(10**20+1,2),evaluate=False)]:
    try: power.as_real_imag()
    except NotImplementedError as error: assert 'limit' in str(error)
    else: raise AssertionError('oversized power was not bounded')
value=s.Integer(2)
for _ in range(80): value=s._core._restore_arithmetic('Mul',s.Integer(1),[(value,s.Integer(1)),(value,s.Integer(1))],False)
try: value.as_real_imag()
except NotImplementedError: pass
else: raise AssertionError('numeric DAG must be rejected before materialization')
r,t=s.symbols('r t',real=True)
a,b=((r+s.I*t)**3).as_real_imag()
assert s.expand(a+s.I*b-(r+s.I*t)**3)==0
'''
    result = subprocess.run([sys.executable, '-c', source], capture_output=True, text=True, timeout=15)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('value', [0, 0.0, -0.0, '0.0', '-0.0'])
def test_complex_logarithm_of_exact_machine_and_arbitrary_precision_zero_is_zoo(value):
    for zero in [sx.sympify(value), sx.Float(str(value), 100)]:
        logarithm = sx.log(zero)
        assert logarithm == sx.zoo
        assert logarithm.as_real_imag() == (sx.nan, sx.nan)
        assert pickle.loads(pickle.dumps(logarithm)) == logarithm


@pytest.mark.parametrize('function', ['sin', 'cos', 'tan', 'cot'])
@pytest.mark.parametrize('sign', [-1, 1])
def test_oscillatory_infinite_arguments_require_an_explicit_unsupported_accumbounds(function, sign):
    value = getattr(sx, function)(sign*sx.oo)
    expected = getattr(sp, function)(sign*sp.oo)
    assert isinstance(expected, sp.AccumBounds)
    with pytest.raises(NotImplementedError, match='AccumBounds'):
        value.as_real_imag()


@pytest.mark.parametrize('source', [
    'sinh(oo)', 'sinh(-oo)', 'cosh(oo)', 'cosh(-oo)',
    'tanh(oo)', 'tanh(-oo)', 'coth(oo)', 'coth(-oo)',
    'exp(oo)', 'exp(-oo)', 'log(-oo)', 'sin(I*oo)', 'cos(I*oo)',
])
def test_defined_infinite_elementary_values_match_principal_branch(source):
    actual = eval(source, {'__builtins__': {}}, environment(sx))
    expected = eval(source, {'__builtins__': {}}, environment(sp))
    assert tuple(map(sx.to_sympy, sx.sympify(actual).as_real_imag())) == expected.as_real_imag()


@pytest.mark.parametrize('power', [sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(2, 3),
                                 sp.Rational(-1, 2), sp.Rational(-3, 2)])
def test_principal_rational_power_parts_reconstruct_on_both_sides_of_the_branch_cut(power):
    z = sx.Symbol('z')
    a, b = (z**sx.from_sympy(power)).as_real_imag()
    expression = sx.to_sympy(a+sx.I*b)
    Z = sp.Symbol('z')
    for point in [sp.Integer(-4), sp.Integer(4), -4+sp.I/10, -4-sp.I/10, 2+3*sp.I]:
        residual = (expression.subs(Z, point)-point**power).evalf(50)
        assert abs(complex(residual)) < 1e-45


@pytest.mark.parametrize('kind', ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix'])
@pytest.mark.parametrize('shape', [(0, 0), (0, 3), (4, 0)])
def test_empty_matrix_parts_preserve_both_dimensions(kind, shape):
    value = getattr(sx, kind)(*shape, [])
    for part in value.as_real_imag():
        assert part.shape == shape
        assert type(part) is type(value)
        assert part._immutable == value._immutable


def test_high_precision_and_dummy_identity_survive_parts_and_all_pickle_protocols():
    x = sx.Dummy('x', real=True)
    y = sx.Dummy('x', real=True)
    coefficient = sx.Float('1.23456789012345678901234567890123456789', 100)
    value = coefficient*x+sx.I*y
    a,b = value.as_real_imag()
    assert a == coefficient*x and b == y
    assert a.free_symbols == {x} and b.free_symbols == {y} and x != y
    assert sx.to_sympy(a).atoms(sp.Float) == {sx.to_sympy(coefficient)}
    for protocol in range(pickle.HIGHEST_PROTOCOL+1):
        assert pickle.loads(pickle.dumps((a,b),protocol=protocol)) == (a,b)


def test_genuine_matrix_scalar_nodes_stay_opaque_and_preserve_metadata():
    matrix = sx.MatrixSymbol('A',2,2)
    for value in [sx.Trace(matrix), sx.Determinant(matrix), matrix[0,1]]:
        real,imag = value.as_real_imag()
        assert real.args == imag.args == (value,)
        assert tuple(map(sx.to_sympy,(real,imag))) == sx.to_sympy(value).as_real_imag()


def test_held_zero_and_undefined_products_keep_certified_values_over_sympy_artifacts():
    # SymPy's held zero**(1/3) polar formula uses atan2(0,0), yielding nan;
    # the positive rational power itself is exactly zero.
    value = sx.Pow(0, sx.Rational(1,3), evaluate=False)
    assert value.as_real_imag() == (sx.Integer(0),sx.Integer(0))
    assert sp.Pow(0, sp.Rational(1,3), evaluate=False).as_real_imag() == (sp.nan,sp.nan)
    # Undefined multiplication must not acquire a certified imaginary zero.
    value = sx.Mul(0, sx.oo, evaluate=False)
    assert value.as_real_imag() == (sx.nan,sx.nan)
    assert sp.Mul(0,sp.oo,evaluate=False).as_real_imag() == (sp.nan,sp.Integer(0))


@pytest.mark.parametrize('fact', ['real','imaginary','extended_real','finite','infinite',
    'nonzero','zero','positive','negative','nonnegative','nonpositive','complex',
    'integer','rational','irrational','algebraic','transcendental'])
@pytest.mark.parametrize('truth', [False,True])
def test_full_function_facts_control_projections_without_becoming_continuity_assumptions(fact,truth):
    x = sp.Symbol('x')
    value = sp.Function('f',**{fact:truth})(x)
    native = sx.from_sympy(value)
    for expression,reference in [(native,value),(2*native,2*value),(sx.I*native,sp.I*value)]:
        assert tuple(map(sx.to_sympy,expression.as_real_imag())) == reference.as_real_imag()
