"""Exact modular inverses over integer residues, including composite rings."""
import itertools
import math
import pickle
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp
import symox as sx
from symox import NonInvertibleMatrixError, NonSquareMatrixError


def certify(a, b, modulus):
    assert a.shape == b.shape
    for product in (a*b, b*a):
        for i in range(a.rows):
            for j in range(a.cols):
                assert (int(product[i,j])-int(i==j)) % modulus == 0
    assert all(0 <= int(value) < modulus for value in b)


@pytest.mark.parametrize('modulus', [2, 4, 6, 8, 9, 12, 15, 35, 101, 2**127-1])
@pytest.mark.parametrize('immutable', [False, True])
def test_composite_unit_creation_and_two_sided_certificate(modulus, immutable):
    cls = sx.ImmutableMatrix if immutable else sx.Matrix
    # Mod6, the first column contains 2 and3, both nonunits. Their Bezout
    # combination is1: searching for an existing unit pivot would be wrong.
    a = cls([[2,1],[3,2]])
    b = a.inv_mod(modulus)
    assert isinstance(b,sx.ImmutableMatrix)==immutable
    certify(a,b,modulus)
    assert b.to_sympy()==sp.Matrix([[2,1],[3,2]]).inv_mod(modulus)


def test_small_composite_exhaustive_against_independent_adjugate_formula():
    for modulus in [4,6,9]:
        for a,b,c,d in itertools.product(range(4),repeat=4):
            matrix=sx.Matrix([[a,b],[c,d]])
            det=a*d-b*c
            if math.gcd(det,modulus)==1:
                result=matrix.inv_mod(modulus)
                scale=pow(det,-1,modulus)
                expected=sx.Matrix([[d*scale%modulus,-b*scale%modulus],[-c*scale%modulus,a*scale%modulus]])
                assert result==expected
                certify(matrix,result,modulus)
            else:
                with pytest.raises(NonInvertibleMatrixError):matrix.inv_mod(modulus)


@pytest.mark.parametrize('shape',[(0,0),(1,1),(3,3)])
def test_trivial_ring_modulus_one(shape):
    # SymPy1.14 raises IndexError for nonempty matrices modulo1. The zero
    # matrix is the mathematically correct inverse in the trivial ring.
    matrix=sx.ImmutableMatrix(sx.ones(*shape))
    result=matrix.inv_mod(1)
    assert result==sx.zeros(*shape)
    certify(matrix,result,1)


@pytest.mark.parametrize('modulus',[0,-1,-17])
def test_nonpositive_modulus(modulus):
    for matrix in [sx.eye(2),sx.zeros(0)]:
        with pytest.raises(ValueError):matrix.inv_mod(modulus)


@pytest.mark.parametrize('modulus',[True,False,5.,sx.Float(5),sx.Rational(5,2),sx.Symbol('p',integer=True)])
def test_modulus_type(modulus):
    with pytest.raises(TypeError):sx.eye(2).inv_mod(modulus)


def test_numpy_integer_modulus_and_immutable_zero_dimensions():
    assert sx.eye(2).inv_mod(np.int64(5))==sx.eye(2)
    assert sx.ImmutableMatrix(sx.zeros(0)).inv_mod(6).shape==(0,0)
    for shape in [(0,3),(4,0),(2,3)]:
        with pytest.raises(NonSquareMatrixError):sx.zeros(*shape).inv_mod(6)


def test_exact_integer_float_entries_and_bigint_residues():
    for number in [sx.Float(1),sx.Float('1.0',80),sx.Float('1e40',80)]:
        result=sx.ImmutableMatrix([[number,1],[1,0]]).inv_mod(35)
        certify(sx.Matrix([[int(number),1],[1,0]]),result,35)
    matrix=sx.ImmutableMatrix([[2**200+1,2],[3,5]])
    modulus=2**127-1
    result=matrix.inv_mod(modulus)
    certify(matrix,result,modulus)
    assert result.to_sympy()==matrix.to_sympy().inv_mod(modulus)
    assert pickle.loads(pickle.dumps(result))==result


@pytest.mark.parametrize('modulus',[1,5])
@pytest.mark.parametrize('entry',[sx.Rational(3,2),sx.Float('1.5',80),sx.Symbol('n',integer=True),sx.I,sx.pi])
def test_noninteger_entries(entry,modulus):
    with pytest.raises(ValueError):sx.Matrix([[entry]]).inv_mod(modulus)


def test_held_integer_expressions_are_rejected_without_materialization():
    # This matches SymPy's domain coercion: held arithmetic is not an Integer
    # entry, even when its mathematical value is integral.
    subprocess.run([sys.executable,'-c','''
import symox as s
held=s.Pow(2,10**9,evaluate=False)
try:s.eye(2).inv_mod(held)
except TypeError:pass
else:raise AssertionError('held modulus must be rejected')
try:s.Matrix([[held]]).inv_mod(5)
except ValueError:pass
else:raise AssertionError('held entry must be rejected')
'''],check=True,timeout=10)

# Native Rust additional test: Number::MpFloat(BigFloat{man:1,exp:10**9},53)
# as an integer entry modulo127. Compute 2.pow_mod(exp,127), never BigInt1<<exp.
# Add oversized modulus/storage/work guard probes through tiny private budgets,
# so the tests prove preflight rejection without constructing giant integers.


def test_bezout_uses_more_than_two_rows_and_no_individual_units():
    for reference, modulus in [(sp.Matrix([[2,3],[3,2]]),6),
                               (sp.Matrix([[6,1,0],[10,0,1],[15,1,1]]),30)]:
        a=sx.ImmutableMatrix(reference)
        result=a.inv_mod(modulus)
        assert result.to_sympy()==reference.inv_mod(modulus)
        certify(a,result,modulus)


@pytest.mark.parametrize('modulus',[6,12,25,2**100*3**30])
def test_unimodular_four_by_four_independent_differential(modulus):
    reference=sp.eye(4)
    for i,j,c in [(0,1,3),(2,0,-7),(1,3,5),(3,2,11),(0,3,-2),(1,2,13)]:
        reference[i,:]=reference[i,:]+c*reference[j,:]
    a=sx.ImmutableMatrix(reference)
    result=a.inv_mod(modulus)
    assert result.to_sympy()==reference.inv_mod(modulus)
    certify(a,result,modulus)


@pytest.mark.parametrize('immutable',[False,True])
def test_sparse_inverse_uses_native_preflight_before_densification(immutable):
    constructor=sx.ImmutableSparseMatrix if immutable else sx.SparseMatrix
    a=constructor([[2,3],[3,2]])
    result=a.inv_mod(m=6)
    assert isinstance(result,sx.ImmutableSparseMatrix)==immutable
    assert result==a
    certify(a,result,6)
    # This shape has only two stored entries. Refusal must happen before an
    # as_dense call could allocate its trillion scalar entries.
    huge=constructor(10**6,10**6,{(0,0):1,(1,1):1})
    with pytest.raises(NotImplementedError,match='resource limit'):
        huge.inv_mod(6)
    with pytest.raises(TypeError):huge.inv_mod(5.0)
    with pytest.raises(ValueError):huge.inv_mod(0)


def test_native_without_sympy_and_deep_raw_tags():
    subprocess.run([sys.executable,'-c','''
import sys
sys.modules['sympy']=None
import symox as s
one=s.Integer(1)
raw=s.Integer(2)
for _ in range(80):
    raw=s._core._restore_arithmetic('Mul',one,[(raw,one),(raw,one)],False)
try:s.eye(2).inv_mod(raw)
except TypeError:pass
else:raise AssertionError('raw modulus requires an Integer tag')
try:s.Matrix([[raw]]).inv_mod(5)
except ValueError:pass
else:raise AssertionError('raw entries must not be evaluated')
a=s.ImmutableMatrix([[6,1,0],[10,0,1],[15,1,1]])
b=a.inv_mod(30)
assert isinstance(b,s.ImmutableMatrix)
for product in [a*b,b*a]:
    for i in range(3):
        for j in range(3):assert (int(product[i,j])-int(i==j))%30==0
assert s.ImmutableSparseMatrix([[2,3],[3,2]]).inv_mod(6)==s.ImmutableSparseMatrix([[2,3],[3,2]])
'''],check=True,timeout=10)


def test_modulus_resource_bound_before_arithmetic_and_sparse_allocation():
    modulus=1<<300000
    with pytest.raises(NotImplementedError,match='resource limit'):
        sx.eye(1).inv_mod(modulus)
    with pytest.raises(NotImplementedError,match='resource limit'):
        sx.SparseMatrix(1,1,{(0,0):1}).inv_mod(modulus)
    enormous=sx.SparseMatrix(10**50,10**50,{})
    with pytest.raises(NotImplementedError,match='resource limit'):
        enormous.inv_mod(6)
