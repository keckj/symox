"""Tensor/kind doc examples of SymPy 1.14; expected values come from SymPy."""
from symox import (Array, ImmutableDenseMatrix, Integral, Matrix, MatrixSymbol,
                   Sum, shape, symbols)
from symox.core.kind import NumberKind
from symox.combinatorics import Permutation
from symox.matrices import MatrixKind
from symox.tensor.array import (ArrayKind, NDimArray, permutedims,
                                tensorcontraction, tensorproduct)
from symox.tensor.array.expressions import (ArrayContraction, ArrayTensorProduct,
                                            PermuteDims, convert_array_to_matrix,
                                            convert_indexed_to_array,
                                            convert_matrix_to_array)

x = symbols('x')
i, j, k, d = symbols('i j k d')


def test_integral_of_shaped_integrands_keeps_kind_and_shape():
    assert shape(Integral(Array([1, 2]), x)) == (2,)
    assert shape(Integral(ImmutableDenseMatrix([1, 2]), x)) == (2, 1)
    assert Integral(NDimArray([1, 2, 3]), x).kind is ArrayKind(NumberKind)
    A = MatrixSymbol('A', 2, 2)
    intA = Integral(A, x)
    assert isinstance(intA, Integral)
    assert intA.kind is MatrixKind(NumberKind)
    assert intA.doit() == x*A
    assert Integral(Matrix([1, x]), x).doit() == Matrix([x, x**2/2])


def test_permute_dims_of_symbolic_arrays():
    M = MatrixSymbol('M', 3, 3)
    assert convert_array_to_matrix(permutedims(M, [1, 0])) == M.T
    N = MatrixSymbol('N', 3, 2)
    assert permutedims(N, [1, 0]).shape == (2, 3)
    M = MatrixSymbol('M', 4, 5)
    tp = tensorproduct(M, N)
    assert tp.shape == (4, 5, 3, 2)
    perm1 = permutedims(tp, [2, 3, 1, 0])
    assert isinstance(perm1, PermuteDims)
    assert perm1.expr.args == (N, M)
    assert perm1.shape == (3, 2, 5, 4)
    assert perm1.permutation == Permutation(2, 3)
    assert perm1.permutation.array_form == [0, 1, 3, 2]
    perm2 = permutedims(perm1, [1, 0, 2, 3])
    assert perm2.shape == (2, 3, 5, 4)
    assert perm2.permutation.array_form == [1, 0, 3, 2]


def test_convert_array_to_matrix_of_contracted_products():
    A, B, C, D = [MatrixSymbol(s, k, k) for s in 'ABCD']
    cg = tensorcontraction(tensorproduct(A, B, C, D), (1, 2), (5, 6))
    assert str(cg) == 'ArrayContraction(ArrayTensorProduct(A, B, C, D), (1, 2), (5, 6))'
    assert convert_array_to_matrix(cg) == ArrayTensorProduct(A*B, C*D)
    cg = tensorcontraction(tensorproduct(A, B, C), (0, 3), (1, 4))
    assert str(cg) == 'ArrayContraction(ArrayTensorProduct(A, B, C), (0, 3), (1, 4))'
    assert convert_array_to_matrix(cg) == B*A*C
    cg = tensorcontraction(tensorproduct(A, B), (0, 3))
    assert str(cg) == 'ArrayContraction(ArrayTensorProduct(A, B), (0, 3))'
    assert convert_array_to_matrix(cg) == A.T*B.T


def test_sort_args_by_name():
    A, B, C, D = [MatrixSymbol(s, k, k) for s in 'ABCD']
    cg = convert_matrix_to_array(C*D*A*B)
    assert str(cg) == 'ArrayContraction(ArrayTensorProduct(A, D, C, B), (0, 3), (1, 6), (2, 5))'
    assert str(cg.sort_args_by_name()) == \
        'ArrayContraction(ArrayTensorProduct(A, D, B, C), (0, 3), (1, 4), (2, 7))'


def test_convert_indexed_to_array():
    M = MatrixSymbol('M', d, d)
    N = MatrixSymbol('N', d, d)
    assert str(convert_indexed_to_array(M[i, i])) == 'ArrayDiagonal(M, (0, 1))'
    assert convert_indexed_to_array(Sum(M[i, i], (i, 0, d - 1))) == ArrayContraction(M, (0, 1))
    expr = Sum(M[i, j]*N[j, k], (j, 0, d - 1))
    assert convert_indexed_to_array(expr) == ArrayContraction(ArrayTensorProduct(M, N), (1, 2))
    assert str(convert_indexed_to_array(expr, first_indices=[k])) == \
        'ArrayContraction(ArrayTensorProduct(N, M), (0, 3))'
    assert convert_array_to_matrix(convert_indexed_to_array(expr)) == M*N


def test_single_axis_contraction():
    A = MatrixSymbol('A', d, d)
    # SymPy 1.14 gives ArrayContraction(A, (1,)), native since arrindex.
    result = convert_indexed_to_array(Sum(A[i, j], (j, 0, d - 1)))
    assert str(result) == 'ArrayContraction(A, (1,))'
    assert result.shape == (d,)
