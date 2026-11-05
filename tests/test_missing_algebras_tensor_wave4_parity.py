"""Wave-4 parity regressions for sympy.algebras (Quaternion) and sympy.tensor.

Expected values were produced by the reference SymPy 1.14.
"""
import pickle

import pytest

import symox as sp
import symox as sx
from symox.tensor.functions import NoShapeError, shape
from symox.tensor.tensor import (TensorHead, TensorIndexType, TensorSymmetry,
                                 get_symmetric_group_sgs, tensor_heads, tensor_indices)
from symox.tensor.toperators import PartialDerivative
from symox import (I, Matrix, Rational, Quaternion, conjugate, cos, diff,
                   integrate, latex, pi, sin, sqrt, symbols)

qa, qb, qc, qd, w, x, y, z = symbols('a b c d w x y z')
Q = Quaternion(1, 2, 3, 4)


STR_CASES = [
    (lambda: Q, '1 + 2*i + 3*j + 4*k'),
    (lambda: Quaternion(3 + 4*I, 2 + 5*I, 0, 7 + 8*I, real_field=False),
     '(3 + 4*I) + (2 + 5*I)*i + 0*j + (7 + 8*I)*k'),
    (lambda: Q.add(Quaternion(5, 6, 7, 8)), '6 + 8*i + 10*j + 12*k'),
    (lambda: Q.add(5 + 2*x), '(2*x + 6) + 2*i + 3*j + 4*k'),
    (lambda: Q + x, '(x + 1) + 2*i + 3*j + 4*k'),
    (lambda: Q.mul(Quaternion(5, 6, 7, 8)), '(-60) + 12*i + 30*j + 24*k'),
    (lambda: Q.mul(2 + I), '0 + 5*i + 10*j + 5*k'),
    (lambda: Q.add(2 + 5*I), '3 + 7*i + 3*j + 4*k'),
    (lambda: Quaternion(3 + 4*I, 2 + 5*I, 0, 7 + 8*I, real_field=False).add(2 + 5*I),
     '(5 + 9*I) + (2 + 5*I)*i + 0*j + (7 + 8*I)*k'),
    (lambda: Quaternion.from_euler([pi/2, 0, 0], 'xyz'),
     'sqrt(2)/2 + sqrt(2)/2*i + 0*j + 0*k'),
    (lambda: Quaternion.from_euler([0, pi/2, pi], 'zyz'),
     '0 + (-sqrt(2)/2)*i + 0*j + sqrt(2)/2*k'),
    (lambda: Quaternion.from_euler([0, pi/2, pi], 'ZYX'),
     '0 + sqrt(2)/2*i + 0*j + (-sqrt(2)/2)*k'),
    (lambda: Quaternion(qa, qb, qc, qd).to_euler('zyz'),
     '(-atan2(-b, c) + atan2(d, a), 2*atan2(sqrt(b**2 + c**2), sqrt(a**2 + d**2)), '
     'atan2(-b, c) + atan2(d, a))'),
    (lambda: Q.to_euler('XYZ'),
     '(-3*pi/4 + atan(3/2), -pi/2 + 2*atan(sqrt(26)/2), atan(3/2) + 3*pi/4)'),
    (lambda: Q.to_euler('xyx', angle_addition=False),
     '(atan(2/11), 2*atan(sqrt(5)), pi - atan(2))'),
    (lambda: Q.to_euler('xyz', avoid_square_root=True),
     '(-atan(3) - atan(1/2) + pi, -asin(1/3), -atan(3) + atan(1/2) + pi)'),
    (lambda: Quaternion.from_axis_angle((sqrt(3)/3, sqrt(3)/3, sqrt(3)/3), 2*pi/3),
     '1/2 + 1/2*i + 1/2*j + 1/2*k'),
    (lambda: Q.pow(4), '668 + (-224)*i + (-336)*j + (-448)*k'),
    (lambda: Q**-2, '(-7/225) + (-1/225)*i + (-1/150)*j + (-2/225)*k'),
    (lambda: Q.inverse(), '1/30 + (-1/15)*i + (-1/10)*j + (-2/15)*k'),
    (lambda: Q / Quaternion(1, 1, 0, 0), '3/2 + 1/2*i + (-1/2)*j + 7/2*k'),
    (lambda: 2 / Q, '1/15 + (-2/15)*i + (-1/5)*j + (-4/15)*k'),
    (lambda: conjugate(Q), '1 + (-2)*i + (-3)*j + (-4)*k'),
    (lambda: Q.normalize(),
     'sqrt(30)/30 + sqrt(30)/15*i + sqrt(30)/10*j + 2*sqrt(30)/15*k'),
    (lambda: Q.exp(),
     'E*cos(sqrt(29)) + 2*sqrt(29)*E*sin(sqrt(29))/29*i + '
     '3*sqrt(29)*E*sin(sqrt(29))/29*j + 4*sqrt(29)*E*sin(sqrt(29))/29*k'),
    (lambda: Q.log(),
     'log(sqrt(30)) + 2*sqrt(29)*acos(sqrt(30)/30)/29*i + '
     '3*sqrt(29)*acos(sqrt(30)/30)/29*j + 4*sqrt(29)*acos(sqrt(30)/30)/29*k'),
    (lambda: Q.evalf(5), '1.0 + 2.0*i + 3.0*j + 4.0*k'),
    (lambda: Quaternion(x, x**2, x**3, x**4).diff(x), '1 + 2*x*i + 3*x**2*j + 4*x**3*k'),
    (lambda: diff(Quaternion(x, x**2, x**3, x**4), x), '1 + 2*x*i + 3*x**2*j + 4*x**3*k'),
    (lambda: Quaternion(x, x**2, x**3, x**4).integrate(x),
     'x**2/2 + x**3/3*i + x**4/4*j + x**5/5*k'),
    (lambda: Q.integrate(x, (x, 1, 5)), '12 + 24*i + 36*j + 48*k'),
    (lambda: integrate(Quaternion(x, x, x, x), x),
     'x**2/2 + x**2/2*i + x**2/2*j + x**2/2*k'),
    (lambda: Quaternion.rotate_point((1, 1, 1), Quaternion(cos(x/2), 0, 0, sin(x/2))),
     '((-sin(x/2) + cos(x/2))*cos(x/2) - (sin(x/2) + cos(x/2))*sin(x/2), '
     '-(sin(x/2) - cos(x/2))*sin(x/2) + (sin(x/2) + cos(x/2))*cos(x/2), '
     'sin(x/2)**2 + cos(x/2)**2)'),
    (lambda: Quaternion(1, 1, 1, 1).axis(), '0 + sqrt(3)/3*i + sqrt(3)/3*j + sqrt(3)/3*k'),
    (lambda: Quaternion(2, 4, 2, 4).index_vector(), '0 + 4*sqrt(10)/3*i + 2*sqrt(10)/3*j + 4*sqrt(10)/3*k'),
    (lambda: 1 - Q, '1 + (-1) + (-2)*i + (-3)*j + (-4)*k'),
    (lambda: Q.to_Matrix(vector_only=True), 'Matrix([[2], [3], [4]])'),
    (lambda: (Quaternion.from_Matrix([qa, qb, qc, qd]), Quaternion.from_Matrix([qb, qc, qd])),
     '(a + b*i + c*j + d*k, 0 + b*i + c*j + d*k)'),
    (lambda: Q.product_matrix_left,
     'Matrix([[1, -2, -3, -4], [2, 1, -4, 3], [3, 4, 1, -2], [4, -3, 2, 1]])'),
]


@pytest.mark.parametrize('build,expected', STR_CASES)
def test_quaternion_values(build, expected):
    assert str(build()) == expected


def test_quaternion_scalar_results():
    assert Quaternion(1, 4, 4, 4).angle() == 2*sp.atan(4*sqrt(3))
    assert Quaternion(2, 4, 2, 4).mensor() == sp.log(2*sqrt(10))
    assert Quaternion.to_axis_angle(Quaternion(1, 1, 1, 1)) == (
        (sqrt(3)/3, sqrt(3)/3, sqrt(3)/3), 2*pi/3)
    assert latex(Quaternion(x, x**3, x + 1, -x**2)) == \
        r'x + x^{3} i + \left(x + 1\right) j + - x^{2} k'


def test_quaternion_predicates():
    q0 = Quaternion(0, 0, 0, 0)
    assert Quaternion.arc_coplanar(Quaternion(1, 4, 4, 4), Quaternion(3, 8, 8, 8)) is True
    assert Quaternion.vector_coplanar(Quaternion(0, 8, 12, 16), Quaternion(0, 4, 6, 8),
                                      Quaternion(0, 2, 3, 4)) is True
    assert Quaternion(0, 1, 2, 3).parallel(Quaternion(0, 2, 4, 6)) is True
    assert Quaternion(0, 1, 2, 3).orthogonal(Quaternion(0, 2, 4, 6)) is False
    with pytest.raises(ValueError, match='Neither of the given quaternions can be 0'):
        Quaternion.arc_coplanar(q0, Quaternion(1, 2, 3, 4))
    with pytest.raises(ValueError, match='must be pure'):
        Quaternion.vector_coplanar(Q, Q, Q)
    with pytest.raises(ValueError, match='Consecutive axes must be different'):
        Quaternion.from_euler([0, pi/2, pi], 'zzx')
    with pytest.raises(ValueError, match='Incompatible value for norm'):
        Quaternion(1, 2, 3, 4, norm=5)
    with pytest.raises(ValueError, match='length 3 or 4'):
        Quaternion.from_Matrix([qa, qb])
    with pytest.raises(TypeError):
        Q**Rational(1, 2)


def test_quaternion_norm_and_pickle():
    q = Quaternion(w, x, y, z, norm=1)
    assert q.norm() == 1
    q2 = Quaternion(w, x, y, z)
    q2.set_norm(1)
    assert q2.norm() == 1
    r = pickle.loads(pickle.dumps(Q))
    assert r == Q and isinstance(r, Quaternion)


def test_quaternion_euler_roundtrip():
    angles = [Rational(1, 3), Rational(1, 5), Rational(1, 7)]
    for seq in ['xyz', 'zyx', 'xyx', 'zxz', 'XYZ', 'ZYZ']:
        res = Quaternion.from_euler(angles, seq).to_euler(seq)
        assert [float(v) for v in res] == pytest.approx([float(v) for v in angles])


def test_quaternion_from_rotation_matrix_pi_rotation():
    # Intentional correction: the reference returns the zero quaternion here.
    for M in [Matrix.diag(1, -1, -1), Matrix.diag(-1, 1, -1),
              Matrix([[0, 1, 0], [1, 0, 0], [0, 0, -1]])]:
        q = Quaternion.from_rotation_matrix(M)
        assert q.to_rotation_matrix() == M


# Tensor, Idx, IndexedBase and shape

i = sx.Symbol('i', integer=True)


def test_idx_ranges():
    assert (sx.Idx(i, sx.oo).lower, sx.Idx(i, sx.oo).upper) == (0, sx.oo)
    n = sx.Symbol('n', integer=True)
    assert (sx.Idx(i, n).lower, sx.Idx(i, n).upper) == (0, n - 1)
    # The bounds are not ordered.
    assert (sx.Idx(i, (3, 1)).lower, sx.Idx(i, (3, 1)).upper) == (3, 1)
    assert sx.Idx(i, (-sx.oo, 2)).upper == 2
    m = sx.Symbol('m')
    assert sx.Idx(i, m).upper == m - 1


def test_idx_errors():
    with pytest.raises(TypeError, match='^Idx object requires an integer dimension.$'):
        sx.Idx(i, 2.5)
    with pytest.raises(ValueError, match='Idx range tuple must have length 2, but got 3'):
        sx.Idx(i, (1, 2, 3))
    with pytest.raises(TypeError, match='^Idx object requires integer bounds.$'):
        sx.Idx(i, (sx.Rational(1, 2), 3))
    with pytest.raises(TypeError, match='^Idx object requires an integer label.$'):
        sx.Idx(sx.Symbol('i'), 3)
    with pytest.raises(TypeError, match='^Index is not an integer number.$'):
        sx.Idx(2.5, 3)
    assert sx.Idx(2, 3) == 2
    assert str(sx.Idx('i', 3)) == 'i'


def test_indexed_base_assumptions():
    A_real = sx.IndexedBase('A', real=True)
    assert A_real.is_real is True
    assert (sx.IndexedBase('A') == A_real) is False
    assert A_real[i].is_real is True
    assert sx.IndexedBase('A').is_real is None
    B = sx.IndexedBase('B', negative=True)
    assert B[i].is_negative is True
    assert B.is_positive is False
    assert sx.IndexedBase(sx.Symbol('I', integer=True)) == sx.IndexedBase('I', integer=True)


def test_shape_of_integral_follows_its_function():
    with pytest.raises(NoShapeError, match=r'^x does not have shape, or its type is not registered to shape\(\)\.$'):
        shape(sx.Integral(x, x))
    assert not issubclass(NoShapeError, ValueError)


Lor = TensorIndexType('Lorentz', dummy_name='L')
a, b, c, d = tensor_indices('a b c d', Lor)
A, B, C = tensor_heads('A B C', [Lor])
H = TensorHead('H', [Lor, Lor])


def test_tensor_symmetry_from_sgs():
    assert TensorSymmetry(get_symmetric_group_sgs(2)) == TensorSymmetry.fully_symmetric(2)
    assert TensorSymmetry(get_symmetric_group_sgs(2, 1)) == TensorSymmetry.fully_symmetric(-2)
    Asy = TensorHead('Asy', [Lor, Lor], TensorSymmetry(get_symmetric_group_sgs(2)))
    assert str((Asy(b, a) - Asy(a, b)).canon_bp()) == '0'


def test_split():
    assert str((A(a)*B(b)).split()) == '[A(a), B(b)]'
    assert str((3*A(d)*B(-d)*C(b)).split()) == '[3*A(L_0), B(-L_0), C(b)]'
    assert str(A(a).split()) == '[A(a)]'
    assert str((x*A(a)*H(-a, b)).split()) == '[x*A(L_0), H(-L_0, b)]'


def test_partial_derivative_dummies():
    expr = PartialDerivative(A(d)*B(-d)*C(c), C(c))
    assert str(expr) == 'PartialDerivative(A(L_0)*B(-L_0)*C(L_1), C(L_1))'
    assert str(expr.get_indices()) == '[L_0, -L_0, L_1, -L_1]'
    assert str(PartialDerivative(H(a, d), A(d))) == 'PartialDerivative(H(a, L_0), A(L_0))'
    assert str(PartialDerivative(H(d, -d), A(b))) == 'PartialDerivative(H(L_0, -L_0), A(b))'
    assert str(PartialDerivative(H(c, d), H(c, d))) == 'PartialDerivative(H(L_0, L_1), H(L_0, L_1))'
    assert str(PartialDerivative(H(a, b), A(a))) == 'PartialDerivative(H(L_0, b), A(L_0))'


def test_repeated_index_errors():
    with pytest.raises(ValueError, match='^Repeated index: a$'):
        A(a)*A(a)
    with pytest.raises(ValueError, match='^two equal contravariant indices in slots 0 and 1$'):
        H(a, a)
    with pytest.raises(ValueError, match='^Repeated index: a$'):
        PartialDerivative(A(a), A(-a))


def test_replace_with_arrays_index_order():
    assert str(PartialDerivative(A(a), A(b)).replace_with_arrays({A(a): [x, y]}, [-b, a])) == '[[1, 0], [0, 1]]'
    arrays = {A(a): [x, y], B(a): [1, 2]}
    assert str((A(a)*B(b)).replace_with_arrays(arrays, [b, a])) == '[[x, y], [2*x, 2*y]]'
    assert str((A(a)*B(b)).replace_with_arrays(arrays, [a, b])) == '[[x, 2*x], [y, 2*y]]'
    assert str(A(-a).replace_with_arrays({A(a): [x, y], Lor: sx.diag(1, -1)}, [-a])) == '[x, -y]'
    expr = PartialDerivative(H(a, b), A(c))
    value = expr.replace_with_arrays({H(a, b): [[x, x*y], [y, x**2]], A(a): [x, y]}, [-c, b, a])
    assert str(value) == '[[[1, 0], [y, 2*x]], [[0, 1], [x, 0]]]'


def test_replace_with_arrays_partial_derivative_doctest():
    i_, j = tensor_indices('i j', Lor)
    compA = [sx.sin(x), sx.log(x)*y**3]
    compB = [x, y]
    expr = PartialDerivative(A(i_), B(j))
    assert str(expr.replace_with_arrays({A(i_): compA, B(i_): compB})) == \
        '[[cos(x), 0], [y**3/x, 3*y**2*log(x)]]'
    assert str(expr.replace_with_arrays({A(i_): compA, B(i_): compB}, [-j, i_])) == \
        '[[cos(x), y**3/x], [0, 3*y**2*log(x)]]'


def test_contracted_derivative_axis():
    # The reference contracts the wrong axis here and gives [x + 1, 0]; the
    # contraction of the derivative over j is [2, y].
    j = tensor_indices('j', Lor)
    expr = PartialDerivative(H(a, j), A(j))
    value = expr.replace_with_arrays({H(a, b): [[x, y], [x*y, x**2]], A(a): [x, y]}, [a])
    assert str(value) == '[2, y]'


def test_basic_matches():
    av, bv, cv = sx.symbols('a b c')
    assert sx.Basic(av + x, x).matches(sx.Basic(av + bv, cv)) is None
    assert sx.Basic(av + x, x).matches(sx.Basic(av + x, x)) == {}


def test_array_expression_converter_modules():
    from symox.tensor.array.expressions import ArraySymbol, PermuteDims
    from symox.tensor.array.expressions.conv_array_to_matrix import convert_array_to_matrix as c2
    from symox.tensor.array.expressions.from_array_to_indexed import convert_array_to_indexed
    from symox.tensor.array.expressions.from_array_to_matrix import convert_array_to_matrix
    from symox.tensor.array.expressions.from_indexed_to_array import convert_indexed_to_array
    from symox.tensor.array.expressions.from_matrix_to_array import convert_matrix_to_array
    # SymPy 1.14: the conv_ name is the from_ function behind a deprecation warning.
    assert c2 is not convert_array_to_matrix and c2.__wrapped__ is convert_array_to_matrix
    assert callable(convert_array_to_indexed) and callable(convert_indexed_to_array)
    M = ArraySymbol('M', (1, 2, 3, 4, 5))
    expr = PermuteDims(M, index_order_old='ijklm', index_order_new='kijml')
    assert str(expr) == 'PermuteDims(M, (0 2 1)(3 4))'
    assert expr.shape == (3, 1, 2, 5, 4)
    assert expr.permutation.array_form == [2, 0, 1, 4, 3]
    N = sx.Symbol('N')
    Ms = [sx.MatrixSymbol(n, N, N) for n in 'CDAB']
    cg = convert_matrix_to_array(Ms[0]*Ms[1]*Ms[2]*Ms[3])
    # SymPy 1.14 canonical forms
    assert str(cg) == 'ArrayContraction(ArrayTensorProduct(A, D, C, B), (0, 3), (1, 6), (2, 5))'
    assert str(cg.sort_args_by_name()) == 'ArrayContraction(ArrayTensorProduct(A, D, B, C), (0, 3), (1, 4), (2, 7))'


def test_replace_with_arrays_symbol_indices():
    L = TensorIndexType('L')
    Ah = TensorHead('A', [L])
    Bh = TensorHead('B', [L])
    i_, j = sx.symbols('i j')
    expr = PartialDerivative(Ah(i_), Bh(j))
    arrays = {Ah(i_): [sx.sin(x), sx.log(x)*y**3], Bh(i_): [x, y]}
    assert str(expr.replace_with_arrays(arrays, [-j, i_])) == '[[cos(x), y**3/x], [0, 3*y**2*log(x)]]'
    order = [-j, i_]
    expr.replace_with_arrays(arrays, order)
    assert str(order) == '[-j, i]'
    assert type(order[0]).__name__ == 'TensorIndex' and order[0].is_up is False
