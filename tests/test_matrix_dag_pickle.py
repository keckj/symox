"""Pickle the mixed native matrix/scalar DAG without printing or evaluating it."""
import pickle
import subprocess
import sys

import pytest
import symox as sx


def primitives(value):
    if type(value) in (tuple, list):
        return all(primitives(item) for item in value)
    return type(value) in (str, int, bool, float, type(None))


def records(value):
    restore, args = value.__reduce__()
    assert restore.__name__ == '_restore_matrix_dag'
    assert primitives(args)
    version, nodes, root = args
    assert version == 1 and 0 <= root < len(nodes)
    return nodes, root


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL + 1))
def test_large_matrix_derivative_pickle_scales_with_shared_dag(protocol):
    source = '''
import pickle,sys,time
import symox as s
protocol=int(sys.argv[1])
x=s.Symbol('x',real=True)
A,B=s.MatrixSymbol('A',2,2),s.MatrixSymbol('B',2,2)
sizes=[]
for exponent in [1000,1000000,1000000000]:
    value=((x*A+B)**exponent).diff(x)
    payload=pickle.dumps(value,protocol=protocol)
    restored=pickle.loads(payload)
    assert restored==value and hash(restored)==hash(value)
    assert restored.free_symbols=={A,B,x}
    assert restored.__reduce__()[1][2] < 5000
    sizes.append(len(payload))
assert sizes[-1]<1000000, sizes
assert sizes[-1]<5*sizes[0], sizes
print(sizes)
'''
    result = subprocess.run([sys.executable, '-c', source, str(protocol)], capture_output=True, text=True, timeout=25)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL + 1))
def test_scalar_matrix_metadata_and_shared_dense_entries(protocol):
    n = sx.Symbol('n', integer=True, positive=True)
    A = sx.MatrixSymbol('A', n, n)
    parent = sx.MatMul(A, A, A, evaluate=False)
    # Symbolic inner dimensions retain a scalar MatrixElement with an opaque parent.
    element = parent[0, 1]
    value = sx.ImmutableMatrix([[element, sx.Trace(parent), sx.Determinant(parent)]])
    original, root = records(value)
    assert sum(tag == 'MatMul' for tag, _, _ in original) == 1
    restored = pickle.loads(pickle.dumps(value, protocol=protocol))
    assert restored == value and isinstance(restored, sx.ImmutableMatrix)
    assert records(restored)[0] == original
    for entry in value:
        scalar = pickle.loads(pickle.dumps(entry, protocol=protocol))
        assert scalar == entry and hash(scalar) == hash(entry)
        assert records(scalar)[0]  # Includes the parent, not merely scalar .args.


@pytest.mark.parametrize('protocol', [0, 2, 4, pickle.HIGHEST_PROTOCOL])
def test_exact_scalar_payloads_held_nodes_functions_subs_and_assumptions(protocol):
    A = sx.MatrixSymbol('A', 2, 2)
    x = sx.Symbol('x', real=False, integer=False)
    y = sx.Symbol('y', positive=True)
    t = sx.Trace(A)
    precision = sx.Float('0.123456789012345678901234567890123456789', 80)
    f = sx.Function('sin', nargs=(1, 2), real=True)
    entries = [
        sx.Mul(-1, t+1, evaluate=False),
        sx.Pow(t+1, 0, evaluate=False),
        sx.Add(t, -t, evaluate=False),
        f(t, x), sx.Subs(f(t, x), x, y),
        sx.Derivative(f(t, x), x, evaluate=False),
        precision, sx.Float(1.25), sx.Rational(10**80+1, 10**70+3),
        x, y, sx.Wild('w', exclude=[t]),
        sx.Piecewise((t, sx.Eq(x, 1)), (y, True)),
        sx.pi, sx.E, sx.I, sx.oo, -sx.oo, sx.zoo, sx.nan,
    ]
    value = sx.MatrixExpr(sx.ImmutableMatrix(1, len(entries), entries))
    before, root = records(value)
    restored = pickle.loads(pickle.dumps(value, protocol=protocol))
    assert isinstance(restored, sx.MatrixExpr)
    assert restored == value
    assert records(restored)[0] == before
    result = restored.as_explicit()
    assert result[0, 6]._prec == precision._prec
    assert result[0, 3].func == f
    assert result[0, 9].is_real is False
    assert result[0, 9].is_integer is False
    assert result[0, 10].is_positive is True


@pytest.mark.parametrize('protocol', [0, 2, pickle.HIGHEST_PROTOCOL])
def test_block_grid_and_matrix_views_roundtrip(protocol):
    A, B = sx.MatrixSymbol('A', 3, 3), sx.MatrixSymbol('B', 2, 2)
    diagonal = sx.BlockDiagMatrix(A, B)
    for value in [diagonal, diagonal[1:4, 1:4], diagonal.blocks, sx.HadamardProduct(A, A), sx.KroneckerProduct(A, A)]:
        before, _ = records(value)
        restored = pickle.loads(pickle.dumps(value, protocol=protocol))
        assert restored == value and hash(restored) == hash(value)
        assert records(restored)[0] == before


def test_shared_matrix_parent_is_preserved_in_scalar_dag():
    n = sx.Symbol('n', integer=True, positive=True)
    A, B = sx.MatrixSymbol('A', n, n), sx.MatrixSymbol('B', n, n)
    parent = sx.MatMul(A, B, evaluate=False)
    value = sx.Trace(parent) + sx.Determinant(parent) + parent[0, 0]
    before, _ = records(value)
    assert sum(tag == 'MatMul' for tag, _, _ in before) == 1
    restored = pickle.loads(pickle.dumps(value))
    after, _ = records(restored)
    assert sum(tag == 'MatMul' for tag, _, _ in after) == 1
    assert restored == value


def test_standalone_wildcard_preserves_matrix_exclusions():
    parent = sx.MatrixSymbol('A', 2, 2)
    wildcard = sx.Wild('w', exclude=[sx.Trace(parent)])
    before, _ = records(wildcard)
    restored = pickle.loads(pickle.dumps(wildcard))
    assert restored == wildcard and hash(restored) == hash(wildcard)
    assert records(restored)[0] == before


@pytest.mark.parametrize('payload', [
    (2, [], 0), (1, [], 0),
    (1, [('unknown', True, ())], 0),
    (1, [('s.pow', True, (0, 0))], 0),
    (1, [('s.number', True, (('rational', 1, 0),))], 0),
    (1, [('s.number', True, (('int', 1),)), ('MatPow', False, (0, 0))], 1),
    (1, [('s.symbol', True, ('bad', 1 << 31, False, []))], 0),
])
def test_restore_rejects_invalid_versions_topology_types_and_numbers(payload):
    from symox._core import _restore_matrix_dag
    with pytest.raises((ValueError, TypeError)):
        _restore_matrix_dag(*payload)


@pytest.mark.parametrize('flags', [{'nonnegative': False}, {'nonpositive': False}])
def test_new_negative_assumption_flags_are_lossless(flags):
    x = sx.Symbol('x', **flags)
    A = sx.MatrixSymbol('A', 1, 1)
    value = x*sx.Trace(A)
    assert pickle.loads(pickle.dumps(value)) == value
    assert records(pickle.loads(pickle.dumps(value)))[0] == records(value)[0]


def test_alpha_equal_subs_keep_distinct_bound_names_in_matrix_storage_and_rewrites():
    x,y,z=sx.symbols('x y z')
    f=sx.Function('f')
    first=sx.Subs(f(x),x,z)
    second=sx.Subs(f(y),y,z)
    assert first==second and first.args!=second.args
    matrix=sx.MatrixExpr(sx.ImmutableMatrix([[first,second]]))
    before,_=records(matrix)
    restored=pickle.loads(pickle.dumps(matrix))
    assert records(restored)[0]==before
    entries=restored.as_explicit()
    assert entries[0,0].args!=entries[0,1].args
    replaced=restored.xreplace({z:1}).as_explicit()
    assert replaced[0,0].args[0]==f(x)
    assert replaced[0,1].args[0]==f(y)
    assert replaced[0,0].args!=replaced[0,1].args
