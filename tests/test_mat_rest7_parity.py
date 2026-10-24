"""mat-rest7 parity: entries of symbolic matrix products as dummy Sums, ranged
KroneckerDelta, derivatives by a matrix element, Trace rewrite(Sum).
Expected values come from SymPy 1.14."""
from symox import *
del test

i, j, k = symbols('i j k')
m, n = symbols('m n')
X = MatrixSymbol('X', k, k)
A = MatrixSymbol('A', k, k)
B = MatrixSymbol('B', k, k)
C = MatrixSymbol('C', k, k)
D = MatrixSymbol('D', k, k)
x = MatrixSymbol('x', k, 1)
y = MatrixSymbol('y', k, 1)
a = MatrixSymbol('a', k, 1)
b = MatrixSymbol('b', k, 1)
c = MatrixSymbol('c', k, 1)


def KDelta(p, q):
    return KroneckerDelta(p, q, (0, k - 1))


def test_ranged_kronecker_delta():
    assert KDelta(i, m).delta_range == Tuple(0, k - 1)
    assert str(KDelta(i, m)) == 'KroneckerDelta(i, m, (0, k - 1))'
    assert KroneckerDelta(i, m).delta_range is None
    assert Sum(KDelta(i, m)*A[i, j], (i, 0, k - 1)).doit() == A[m, j]
    assert str(Sum(KDelta(i, m)*A[i, j], (i, 0, 3)).doit()) == \
        'Piecewise((A[m, j], (m >= 0) & (m <= 3)), (0, True))'


def test_identity_entry():
    assert str(Identity(k)[i, j]) == 'KroneckerDelta(i, j, (0, k - 1))'
    assert Identity(k)[i, j].delta_range == Tuple(0, k - 1)
    assert Identity(3)[0, 0] == 1


def test_symbolic_product_entries():
    assert str((A*B)[i, j]) == 'Sum(A[i, _i_1]*B[_i_1, j], (_i_1, 0, k - 1))'
    assert str((A*B*C)[i, j]) == \
        'Sum(A[i, _i_1]*B[_i_1, _i_2]*C[_i_2, j], (_i_1, 0, k - 1), (_i_2, 0, k - 1))'
    assert str((A*B)[0, 1]) == 'Sum(A[0, _i_1]*B[_i_1, 1], (_i_1, 0, k - 1))'
    assert str((x.T*y)[0, 0]) == 'Sum(x[_i_1, 0]*y[_i_1, 0], (_i_1, 0, k - 1))'
    assert str((A**2)[i, j]) == 'Sum(A[_i_1, j]*A[i, _i_1], (_i_1, 0, k - 1))'
    assert str((X*A)[i, j]) == 'Sum(A[_i_1, j]*X[i, _i_1], (_i_1, 0, k - 1))'
    assert str((A + B)[i, j]) == 'A[i, j] + B[i, j]'
    assert str((2*A)[i, j]) == '2*A[i, j]'
    assert str((Trace(A)*A)[i, j]) == 'A[i, j]*Trace(A)'


def test_element_derivatives():
    assert str(x[i, 0].diff(x[m, 0])) == 'KroneckerDelta(i, m, (0, k - 1))'
    assert str(X[i, j].diff(X[m, n])) == \
        'KroneckerDelta(i, m, (0, k - 1))*KroneckerDelta(j, n, (0, k - 1))'
    assert X[i, j].diff(A[m, n]) == 0
    assert (x.T*a)[0, 0].diff(x[m, 0]).doit() == a[m, 0]
    assert (x.T*y)[0, 0].diff(x[m, 0]).doit() == y[m, 0]
    assert str((X*A)[i, j].diff(X[m, n])) == (
        'Sum(KroneckerDelta(_i_1, n, (0, k - 1))*KroneckerDelta(i, m, (0, k - 1))'
        '*A[_i_1, j], (_i_1, 0, k - 1))')
    assert str((X*A)[i, j].diff(X[m, n]).doit()) == \
        'KroneckerDelta(i, m, (0, k - 1))*A[n, j]'
    assert str((Trace(A)*A)[i, j].diff(A[m, n]).doit()) == (
        'KroneckerDelta(i, m, (0, k - 1))*KroneckerDelta(j, n, (0, k - 1))*Trace(A)'
        ' + KroneckerDelta(m, n, (0, k - 1))*A[i, j]')
    assert str(Inverse(A)[i, j].diff(A[m, n])) == (
        '-Sum(KroneckerDelta(_z1, m, (0, k - 1))*KroneckerDelta(_z2, n, (0, k - 1))'
        '*(A**(-1))[_z2, j]*(A**(-1))[i, _z1], (_z1, 0, k - 1), (_z2, 0, k - 1))')


def test_quadratic_form_element_derivative():
    # SymPy pulls b[n, 0] out of the first Sum; the value is the same.
    expr = ((X*b + c).T*D*(X*b + c))[0, 0].diff(X[m, n]).doit()
    i1, i2, i3, i4 = [Dummy('w%d' % t) for t in range(4)]
    expected = b[n, 0]*Sum((c[i1, 0] + Sum(X[i1, i3]*b[i3, 0], (i3, 0, k - 1)))*D[i1, m],
                           (i1, 0, k - 1)) \
        + Sum((c[i2, 0] + Sum(X[i2, i4]*b[i4, 0], (i4, 0, k - 1)))*D[m, i2]*b[n, 0],
              (i2, 0, k - 1))
    assert expr.has(Sum) and not expr.has(KroneckerDelta)
    assert str(expr.subs(k, 2).doit().expand()) == str(expected.subs(k, 2).doit().expand())


def test_trace_rewrite_sum():
    assert str(Trace(X).rewrite(Sum)) == 'Sum(X[i, i], (i, 0, k - 1))'
    assert str(Trace(X).rewrite(Sum).diff(X[m, n]).doit()) == 'KroneckerDelta(m, n, (0, k - 1))'
    assert Trace(X*A).rewrite(Sum).diff(X[m, n]).doit() == A[n, m]
    assert str(Trace(A*X*B).rewrite(Sum).diff(X[m, n]).doit()) == \
        'Sum(A[i, m]*B[n, i], (i, 0, k - 1))'


def test_nested_trace_derivative():
    assert Trace(Trace(A)*A).diff(A) == 2*Trace(A)*Identity(k)
    assert str(X.diff(A)) == 'Derivative(X, A)'
    assert str((X*A).diff(A)) == 'PermuteDims(ArrayTensorProduct(X.T, I), (3)(1 2))'
