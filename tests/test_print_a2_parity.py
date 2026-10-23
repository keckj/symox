"""Printing parity with SymPy 1.14 for physics, diffgeom, combinatorics,
tensor, categories and geometry objects (audit batch a2).

Expected strings are hard-coded from SymPy 1.14 (module paths renamed to
symox).
"""
import pytest

CASES = {
    'categories': ("\nfrom symox.categories import Object, NamedMorphism, IdentityMorphism, Category\nfrom symox import latex, pretty\nA, B, C = Object('A'), Object('B'), Object('C')\nf = NamedMorphism(A, B, 'f'); g = NamedMorphism(B, C, 'g')\n", [
        ('str(f)', 'NamedMorphism(Object("A"), Object("B"), "f")'),
        ('repr(g*f)', 'CompositeMorphism((NamedMorphism(Object("A"), Object("B"), "f"), NamedMorphism(Object("B"), Object("C"), "g")))'),
        ('str(IdentityMorphism(A))', 'IdentityMorphism(Object("A"))'),
        ('str((g*f).flatten("h"))', 'NamedMorphism(Object("A"), Object("C"), "h")'),
        ('latex(g*f)', 'g\\circ f:A\\rightarrow C'),
        ('latex(f)', 'f:A\\rightarrow B'),
        ('latex(IdentityMorphism(A))', 'id:A\\rightarrow A'),
        ('pretty(g*f, use_unicode=False)', 'g*f:A-->C'),
        ('pretty(f, use_unicode=True)', 'f:A——▶B'),
        ('str(Category("K"))', 'Category("K")'),
        ('latex(Category("K"))', '\\mathbf{K}'),
    ]),
    'combinatorics': ("\nfrom symox.combinatorics import Permutation\nfrom symox.combinatorics.partitions import Partition\nfrom symox.combinatorics.prufer import Prufer\nfrom symox.combinatorics.permutations import AppliedPermutation\nfrom symox.combinatorics.polyhedron import tetrahedron, Polyhedron\nfrom symox import Symbol, latex, pretty\nx = Symbol('x')\n", [
        ('str(Partition.from_rgs([0, 1, 2, 0, 1], list("abcde")))', 'Partition({c}, {a, d}, {b, e})'),
        ('str(Partition([1, 2], [3]))', 'Partition({3}, {1, 2})'),
        ('latex(Partition([1, 2], [3]))', '\\left\\{\\left\\{3\\right\\}, \\left\\{1, 2\\right\\}\\right\\}'),
        ('pretty(Partition([1, 2], [3]))', '{{3}, {1, 2}}'),
        ('repr(Prufer([1, 0, 0]))', 'Prufer([1, 0, 0])'),
        ('repr(Prufer([[0, 1], [0, 2], [0, 3]]))', 'Prufer([[0, 1], [0, 2], [0, 3]])'),
        ('repr(Permutation(0, 1, 2).apply(x))', 'AppliedPermutation((0 1 2), x)'),
        ('str(AppliedPermutation(Permutation(0, 1, 2), x))', 'AppliedPermutation((0 1 2), x)'),
        ('str(Polyhedron(list("abc")).corners)', '(a, b, c)'),
        ('str((lambda T: (T.rotate(0), T.corners)[1])(tetrahedron.copy()))', '(0, 2, 3, 1)'),
        ('str((lambda T: (T.rotate(1), T.copy().corners)[1])(tetrahedron.copy()))', '(0, 1, 2, 3)'),
        ('str([g.array_form for g in tetrahedron.pgroup.generators])', '[[0, 2, 3, 1], [1, 2, 0, 3], [3, 1, 0, 2], [3, 0, 2, 1], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]'),
    ]),
    'diffgeom': ('\nfrom symox.diffgeom.rn import R2_r, R2_p, R2\nfrom symox.diffgeom import TensorProduct as TP, metric_to_Christoffel_2nd, twoform_to_matrix, CovarDerivativeOp, vectors_in_basis\nfrom symox import latex\nfx, fy = R2_r.base_scalars(); e_x, e_y = R2_r.base_vectors(); dx, dy = R2_r.base_oneforms()\n', [
        ('str(fx)', 'x'),
        ('str(fx**2 + fy)', 'x**2 + y'),
        ('str(TP(dx, fx*dy)(fx*e_x, e_y))', 'x**2'),
        ('str(TP(e_x, e_y)(fx**2, fy**2))', '4*x*y'),
        ('str(metric_to_Christoffel_2nd(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy)))', '[[[1/(2*x), 0], [0, 0]], [[0, 0], [0, 0]]]'),
        ('repr(twoform_to_matrix(R2.x*TP(R2.dx, R2.dx) + TP(R2.dy, R2.dy)))', 'Matrix([\n[x, 0],\n[0, 1]])'),
        ('repr(R2_r.transformation(R2_p))', 'Lambda((x, y), Matrix([\n[sqrt(x**2 + y**2)],\n[      atan2(y, x)]]))'),
        ('str(R2_r.jacobian(R2_p))', 'Matrix([[x/sqrt(x**2 + y**2), y/sqrt(x**2 + y**2)], [-y/(x**2 + y**2), x/(x**2 + y**2)]])'),
        ('str(R2_r.transform(R2_p))', 'Matrix([[sqrt(x**2 + y**2)], [atan2(y, x)]])'),
        ('latex(fx)', '\\mathbf{x}'),
        ('latex(fx**2 + fy)', '\\mathbf{x}^{2} + \\mathbf{y}'),
        ('latex(R2.r)', '\\mathbf{\\rho}'),
        ('latex([fx])', '\\left[ \\mathbf{x}\\right]'),
    ]),
    'quantum': ("\nfrom symox.physics.quantum import Ket, Bra, Operator, TimeDepKet, TimeDepBra\nfrom symox.physics.quantum.spin import JzKet\nfrom symox.physics.quantum.gate import H, X, CNOT, SWAP, CGate, T, S, Phase, IdentityGate, Z, CPHASE\nfrom symox.physics.quantum.qubit import Qubit\nfrom symox.physics.quantum.qapply import qapply\nfrom symox import symbols, Rational, I\nn, m, t = symbols('n m t')\nA, B = Operator('A'), Operator('B')\n", [
        ('str(Ket(n, m))', '|nm>'),
        ('str(Bra(n, m) - I*Bra(m, n))', '-I*<mn| + <nm|'),
        ('str(TimeDepKet("psi", "x", t))', '|psix;t>'),
        ('str(TimeDepBra("a", "b", t))', '<ab;t|'),
        ('str(JzKet(1, 0))', '|1,0>'),
        ('str(Ket(n, m).dual_class())', "<class 'symox.physics.quantum.state.Bra'>"),
        ('str(Bra(n, m).dual_class())', "<class 'symox.physics.quantum.state.Ket'>"),
        ('str(A.inv())', 'A**(-1)'),
        ('str(A*A.inv())', '1'),
        ('str(A.inv()*A.inv())', 'A**(-2)'),
        ('str((A + B)**-1)', '(A + B)**(-1)'),
        ('str(A**Rational(1, 2))', 'sqrt(A)'),
        ('str(H(0)*H(1)*Qubit("00"))', 'H(0)*H(1)*|00>'),
        ('str([X(1), CNOT(1, 0), SWAP(0, 1), T(0), S(0), Phase(0), IdentityGate(0), Z(0)])', '[X(1), CNOT(1,0), SWAP(0,1), T(0), S(0), S(0), 1(0), Z(0)]'),
        ('str(CGate((0, 1), X(2)))', 'C((0,1),X(2))'),
        ('str(CPHASE(0, 1))', 'C((0),Z(1))'),
        ('str(qapply(CPHASE(0, 1)*Qubit("11")))', '-|11>'),
    ]),
    'tensor': ("\nfrom symox.tensor.tensor import TensorIndexType, tensor_indices, TensorHead, TensorSymmetry, WildTensorHead\nfrom symox import Symbol, Rational, IndexedBase, Indexed, symbols, latex, Array, transpose, sin\nL = TensorIndexType('Lorentz', dummy_name='L')\na, b, m0, m1 = tensor_indices('a b m0 m1', L)\nA = TensorHead('A', [L]*2); Bt = TensorHead('B', [L])\nAs = TensorHead('As', [L]*2, TensorSymmetry.fully_symmetric(-2))\nx = Symbol('x'); i, j = symbols('i j')\n", [
        ('str(-A(a, b))', '-A(a, b)'),
        ('str(Rational(1, 2)*A(a, b))', '(1/2)*A(a, b)'),
        ('str((x + 1)*A(a, b))', '(x + 1)*A(a, b)'),
        ('str(2*x*A(a, b))', '2*x*A(a, b)'),
        ('str(-Rational(1, 2)*Bt(a))', '-(1/2)*B(a)'),
        ('str((As(m0, -m1)*As(m1, -m0)*As(m0, -m1)).canon_bp())', '-As(m0, -m1)*As(L_0, L_1)*As(-L_0, -L_1)'),
        ('str(type(WildTensorHead("W")(a)))', "<class 'symox.tensor.tensor.WildTensor'>"),
        ('repr(IndexedBase("A"))', 'A'),
        ('repr(Indexed("A", i, j))', 'A[i, j]'),
        ('str(type(IndexedBase("A")[i]))', "<class 'symox.tensor.indexed.Indexed'>"),
        ('latex(IndexedBase("A")[i, j])', '{A}_{i,j}'),
        ('str(transpose(Array([[x, i], [j, sin(x)]])))', '[[x, j], [i, sin(x)]]'),
        ('str(type(transpose(Array([[x, i], [j, sin(x)]]))).__name__)', 'ImmutableDenseNDimArray'),
    ]),
    'physics': ('\nfrom symox.physics.units.systems.si import dimsys_SI\nfrom symox.physics.units import velocity, length\nfrom symox.physics.optics import polarizing_beam_splitter, mueller_matrix, quarter_wave_retarder, half_wave_retarder, phase_retarder\nfrom symox.physics.continuum_mechanics.cable import Cable\nfrom symox.physics.biomechanics import TendonForceLengthDeGroote2016\nfrom symox import symbols, Matrix, Float\nTs, Rs, Tp, Rp = symbols("Ts, Rs, Tp, Rp", positive=True)\nphia, phib = symbols("phi_a, phi_b", real=True)\nth = symbols(\'theta\', real=True)\nl = symbols(\'l_T_tilde\')\n', [
        ('str(dimsys_SI.get_dimensional_dependencies(velocity))', '{Dimension(length): 1, Dimension(time): -1}'),
        ('str(dimsys_SI.get_dimensional_dependencies(length**2))', '{Dimension(length, L): 2}'),
        ('str(polarizing_beam_splitter(Tp, Rs, Ts, Rp, phia, phib))', 'Matrix([[sqrt(Tp), 0, I*sqrt(Rp), 0], [0, sqrt(Ts), 0, -I*sqrt(Rs)*exp(I*phi_a)], [I*sqrt(Rp), 0, sqrt(Tp), 0], [0, -I*sqrt(Rs)*exp(I*phi_b), 0, sqrt(Ts)]])'),
        ('str(mueller_matrix(quarter_wave_retarder(th)))', 'Matrix([[1, 0, 0, 0], [0, cos(4*theta)/2 + 1/2, sin(4*theta)/2, -sin(2*theta)], [0, sin(4*theta)/2, 1/2 - cos(4*theta)/2, cos(2*theta)], [0, sin(2*theta), -cos(2*theta), 0]])'),
        ('str(half_wave_retarder(th))', 'Matrix([[-I*(-sin(theta)**2 + cos(theta)**2), -2*I*sin(theta)*cos(theta)], [-2*I*sin(theta)*cos(theta), -I*(sin(theta)**2 - cos(theta)**2)]])'),
        ('str(phase_retarder(th, symbols("delta")))', 'Matrix([[(exp(I*delta)*sin(theta)**2 + cos(theta)**2)*exp(-I*delta/2), (1 - exp(I*delta))*exp(-I*delta/2)*sin(theta)*cos(theta)], [(1 - exp(I*delta))*exp(-I*delta/2)*sin(theta)*cos(theta), (exp(I*delta)*cos(theta)**2 + sin(theta)**2)*exp(-I*delta/2)]])'),
        ('str((lambda c: (c.apply_load(-1, ("Z", 5, 5, 12, 30)), c.loads_position)[1])(Cable(("A", 0, 10), ("B", 10, 10))))', "{'Z': [5, 5]}"),
        ('str(TendonForceLengthDeGroote2016.with_defaults(l))', 'TendonForceLengthDeGroote2016(l_T_tilde, 0.2, 0.995, 0.25, 33.93669377311689)'),
        ('str(Matrix([[Float(0.5), 1]]))', 'Matrix([[0.500000000000000, 1]])'),
        ('repr(Matrix([[Float(0.5)], [2]]))', 'Matrix([\n[0.5],\n[  2]])'),
        ('str((lambda Dm: dimsys_SI.get_dimensional_dependencies(Dm("length", "L")/Dm("time")))(__import__("symox.physics.units", fromlist=["Dimension"]).Dimension))', '{Dimension(length, L): 1, Dimension(time, T): -1}'),
    ]),
    'geometry': ('\nfrom symox import Point, Segment, Polygon, Rational\n', [
        ('str(Segment(Point(0, 0), Point(6, 6)).perpendicular_bisector(Point(5, 1)))', 'Segment2D(Point2D(5, 1), Point2D(3, 3))'),
        ('str(Polygon(Point(0, 0), Point(1, 0), Point(5, 1), Point(0, 1)).intersection(Polygon(Point(0, 1), Point(-5, 0), Point(0, -4), Point(0, Rational(1, 3)), Point(5, 1))))', '[Segment2D(Point2D(0, 1/3), Point2D(0, 0)), Segment2D(Point2D(5, 1), Point2D(0, 1))]'),
        ('str(Polygon(Point(0, 0), Point(2, 0), Point(2, 2), Point(0, 2)).intersection(Polygon(Point(1, -1), Point(3, 1), Point(1, 3), Point(-1, 1))))', '[Point2D(0, 0), Point2D(0, 2), Point2D(2, 0), Point2D(2, 2)]'),
        ('str(Point(1, 2).sort_key())', "((5, 0, 'Point2D'), (2, (((1, 0, 'Number'), (0, ()), (), 1), ((1, 0, 'Number'), (0, ()), (), 2))), ((1, 0, 'Number'), (0, ()), (), 1), 1)"),
    ]),
}


def _params():
    for group, (setup, items) in CASES.items():
        for expr, expected in items:
            yield pytest.param(setup, expr, expected, id="%s: %s" % (group, expr[:60]))


@pytest.mark.parametrize("setup,expr,expected", list(_params()))
def test_printing_matches_sympy(setup, expr, expected):
    namespace = {}
    exec(setup, namespace)
    assert eval(expr, namespace) == expected
