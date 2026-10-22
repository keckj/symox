"""Wave-4 wrong-result regressions in geometry, vectors, combinatorics and
holonomic numerics. Expected values are the SymPy 1.14 outputs, except where
noted as an intentional divergence (SymPy demonstrably wrong)."""
import pytest

from symox import Rational, asin, acos, pi, sqrt, sin, cos, symbols, QQ
from symox.geometry import (Point, Point3D, Line, Line3D, Ray, Ray3D, Segment,
                            Segment3D, Plane)
from symox.vector import CoordSys3D, express
from symox.combinatorics import Permutation, Subset, free_group
from symox.combinatorics.pc_groups import Collector
from symox.combinatorics.permutations import Cycle
from symox.combinatorics.named_groups import SymmetricGroup
from symox.holonomic.holonomic import HolonomicFunction, DifferentialOperators

x, t = symbols('x t')


# --- Point.is_concyclic / is_collinear are instance methods including self

def test_is_concyclic_includes_self():
    assert Point(1, 1).is_concyclic(Point(1, 0), Point(0, 1), Point(-1, 0)) is False
    assert Point(1, 0).is_concyclic(Point(0, 1), Point(-1, 0), Point(0, -1)) is True
    assert Point.is_concyclic(Point(1, 0), Point(0, 1), Point(-1, 0), Point(0, -1)) is True
    assert Point(0, 0).is_concyclic(Point(1, 0), Point(2, 0)) is False
    assert Point(0, 0).is_concyclic(Point(1, 0)) is True
    assert Point(0, 0).is_concyclic() is True
    assert Point(0, 0).is_concyclic(Point(2, 0), Point(1, 1), Point(1, -1)) is True


def test_is_concyclic_spatial_points():
    assert Point(1, 0, 0).is_concyclic(Point(0, 1, 0), Point(-1, 0, 0), Point(0, -1, 0)) is True
    assert Point(1, 0, 0).is_concyclic(Point(0, 1, 0), Point(-1, 0, 0), Point(0, 0, 1)) is False
    assert Point(1, 0, 0).is_concyclic(Point(0, 1, 0), Point(-1, 0, 0), Point(0, -1, 1)) is False


def test_is_collinear_includes_self():
    assert Point(1, 1).is_collinear(Point(2, 2), Point(3, 4)) is False
    assert Point(1, 1).is_collinear(Point(2, 2), Point(3, 3)) is True
    assert Point(0, 0).is_collinear(Point(1, 0), Point(5, 0), Point(0, 1)) is False
    assert Point(1, 1).is_collinear() is True
    assert Point.affine_rank(Point(0, 0), Point(1, 1), Point(2, 3)) == 2
    assert Point.affine_rank() == -1


# --- Ray / Line / Segment containment

@pytest.mark.parametrize('outer, inner, expected', [
    (Ray((0, 0), (4, 4)), Ray((2, 2), (3, 5)), False),
    (Ray((0, 0), (4, 4)), Ray((2, 2), (3, 3)), True),
    (Ray((0, 0), (4, 4)), Ray((2, 2), (1, 1)), False),
    (Ray((0, 0), (1, 0)), Ray((1, 0), (1, 1)), False),
    (Line((0, 0), (1, 1)), Ray((2, 2), (3, 5)), False),
    (Line((0, 0), (1, 1)), Ray((2, 2), (3, 3)), True),
    (Ray3D((0, 0, 0), (1, 1, 1)), Ray3D((1, 1, 1), (2, 2, 3)), False),
    (Ray3D((0, 0, 0), (1, 1, 1)), Ray3D((1, 1, 1), (2, 2, 2)), True),
    (Line3D((0, 0, 0), (1, 1, 1)), Ray3D((1, 1, 1), (2, 2, 3)), False),
    (Segment3D((0, 0, 0), (2, 2, 2)), Segment3D((1, 1, 1), (2, 2, 3)), False),
])
def test_linear_containment(outer, inner, expected):
    assert outer.contains(inner) is expected


def test_ray_contains_ray_starting_behind_is_false():
    # Intentional divergence: SymPy 1.14 returns True because Ray.contains(Ray)
    # only checks collinearity and the sense, not that the inner ray starts
    # on the outer one; (-1, -1) is not on the ray from the origin.
    assert Ray((0, 0), (4, 4)).contains(Ray((-1, -1), (0, 0))) is False


# --- Plane.angle_between keeps the sign

def test_plane_angle_between_line_is_signed():
    a = Plane(Point3D(1, 2, 2), normal_vector=(1, 2, 3))
    assert a.angle_between(Line3D(Point3D(1, 3, 4), Point3D(2, 2, 2))) == -asin(sqrt(21)/6)
    assert a.angle_between(Line3D(Point3D(2, 2, 2), Point3D(1, 3, 4))) == asin(sqrt(21)/6)
    assert a.angle_between(Ray3D(Point3D(1, 3, 4), Point3D(2, 2, 2))) == -asin(sqrt(21)/6)
    b = Plane(Point3D(0, 0, 0), normal_vector=(0, 0, 1))
    assert b.angle_between(Line3D((0, 0, 0), (1, 0, -1))) == -pi/4
    assert b.angle_between(Line3D((0, 0, 0), (1, 0, 0))) == 0
    assert b.angle_between(Line((0, 0), (1, 1))) is None
    assert a.angle_between(Plane(Point3D(2, 3, 4), normal_vector=(-1, 0, 1))) == acos(sqrt(7)/7)


# --- vectors

def test_vector_equals_across_rotated_frames():
    C = CoordSys3D('C')
    D = C.orient_new_axis('D', pi/2, C.i)
    assert D.j.equals(C.k) is True
    assert D.j.equals(C.j) is False
    assert D.k.equals(-C.j) is True
    assert (D.i + D.j).equals(C.i + C.k) is True
    E = C.orient_new_axis('E', t, C.i)
    assert E.j.equals(C.k) is False
    assert E.j.equals(cos(t)*C.j + sin(t)*C.k) is True
    assert ((sin(x)**2 + cos(x)**2)*C.i).equals(C.i) is True
    assert express(D.j, C) == C.k


def test_vector_projection_reference_coefficient():
    C = CoordSys3D('C')
    i, j, k = C.base_vectors()
    assert (i + j + k).projection(3*i + 4*j, scalar=True) == Rational(7, 3)
    assert str((i + j + k).projection(3*i + 4*j)) == '7/3*C.i + 7/3*C.j + 7/3*C.k'
    assert (2*i).projection(3*i + 4*j, scalar=True) == Rational(3, 2)
    assert (0*i).projection(3*i + 4*j, scalar=True) == 0
    assert str((0*i).projection(3*i + 4*j)) == '0'


# --- Cycle / Permutation / Subset

def test_cycle_iterates_array_form():
    assert list(Cycle(1, 2)) == [0, 2, 1]
    assert list(Cycle(2, 0, 1)) == [1, 2, 0]
    assert Cycle(1, 2).list(5) == [0, 2, 1, 3, 4]
    assert str(Cycle(1, 2)) == '(1 2)'
    assert repr(Cycle(1, 2)) == 'Cycle(1, 2)'
    assert repr(Cycle()) == 'Cycle()' and str(Cycle()) == '()'
    assert repr(Cycle(1, 2)(2, 3)) == 'Cycle(1, 3, 2)'
    assert repr(Cycle(1, 2)(4)) == 'Cycle(1, 2)(4)'
    assert dict(Cycle(Permutation([1, 0, 2]))) == {0: 1, 1: 0}
    assert Cycle(Permutation([1, 0, 2])).size == 2
    assert Cycle().list(3) == [0, 1, 2]
    with pytest.raises(ValueError):
        list(Cycle())
    with pytest.raises(ValueError):
        Cycle(1, 1)
    assert str(Permutation(Cycle(1, 2), size=5)) == '(4)(1 2)'


def test_permutation_atoms_are_elements():
    assert Permutation([[0, 1], [2, 3], [4, 5]]).atoms() == {0, 1, 2, 3, 4, 5}
    assert Permutation(3).atoms() == {0, 1, 2, 3}


def test_subset_indices_missing_element_gives_empty_list():
    assert Subset.subset_indices([1, 6], [1, 3, 2, 5, 4]) == []
    assert Subset.subset_indices([3, 5], [1, 3, 2, 5, 4]) == [1, 3]
    assert Subset.subset_indices([], [1, 3, 2]) == []
    assert Subset.subset_indices([2, 1], [1, 3, 2]) == [2, 0]
    with pytest.raises(ValueError):
        Subset([1, 6], [1, 3, 2])


# --- polycyclic collector

def test_collector_map_relation_and_presentation():
    c = SymmetricGroup(3).polycyclic_group().collector
    F, x0, x1 = free_group('x0, x1')
    assert str(c.map_relation(x1*x0)) == 'x1**2'
    assert str(c.pc_presentation) == '{x1**3: (), x0**2: (), x0**-1*x1*x0: x1**2}'
    # The reference polycyclic generating sequence of SymmetricGroup(4)
    # (symox chooses a different one, see docs/algorithmic-parity.md).
    pcgs = [Permutation(0, 1, 2, 3), Permutation(3)(0, 2, 1),
            Permutation(0, 1)(2, 3), Permutation(0, 3)(1, 2)]
    c4 = Collector(pcgs, [], [2, 3, 2, 2])
    F, y0, y1, y2, y3 = free_group('x0, x1, x2, x3')
    assert str(c4.map_relation(y2*y0)) == 'x3'
    assert str(c4.map_relation(y3*y2)) == 'x3'
    assert str(c4.map_relation(y1*y0)) == 'x1**2*x2'
    assert str(c4.map_relation(y3*y1)) == 'x2'
    assert str(c4.pc_presentation) == (
        '{x3**2: (), x2**2: (), x2**-1*x3*x2: x3, x1**3: (), x1**-1*x3*x1: x2, '
        'x1**-1*x2*x1: x2*x3, x0**2: x2*x3, x0**-1*x3*x0: x2, x0**-1*x2*x0: x3, '
        'x0**-1*x1*x0: x1**2*x2}')


def test_collector_keeps_whole_power_relator():
    # Intentional divergence: SymPy 1.14 collected_word keeps only the first
    # syllable of the power relator x0**2 = x2*x3 and returns x0*x1**2 for
    # x1*x0**3, which is a different permutation.
    pcgs = [Permutation(0, 1, 2, 3), Permutation(3)(0, 2, 1),
            Permutation(0, 1)(2, 3), Permutation(0, 3)(1, 2)]
    c4 = Collector(pcgs, [], [2, 3, 2, 2])
    F, y0, y1, y2, y3 = free_group('x0, x1, x2, x3')
    word = c4.collected_word(y1*y0**3)
    assert str(word) == 'x0*x1**2*x3'
    assert c4._permutation(word) == pcgs[1]*pcgs[0]**3
    assert c4.exponent_vector(Permutation(0, 2)(1, 3)) == [0, 0, 1, 1]


def test_collector_subword_index_letters():
    c = SymmetricGroup(4).polycyclic_group().collector
    F, x1, x2 = free_group('x1, x2')
    w = x2**2*x1**7
    assert c.subword_index(w, x2**2*x1) == (0, 3)
    assert c.subword_index(w, x1**7) == (2, 9)
    assert c.subword_index(w, x1**8) == (-1, -1)
    assert c.subword_index(x1**3, x1**-1) == (-1, -1)


# --- holonomic numerics

def _approx(values, expected):
    return all(abs(complex(v) - e) < 1e-12 for v, e in zip(values, expected))


def test_holonomic_evalf_uses_points_as_mesh():
    R, Dx = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
    exp_ = HolonomicFunction(Dx - 1, x, 0, [1])
    assert _approx(exp_.evalf([0.1, 0.2], method='Euler'), [1.1, 1.21])
    assert _approx(exp_.evalf([0.5, 1], method='Euler'), [1.5, 2.25])
    assert _approx(exp_.evalf([0.1, 0.2]), [1.10517083333333, 1.22140257085069])
    sin_ = HolonomicFunction(Dx**2 + 1, x, 0, [0, 1])
    assert _approx(sin_.evalf([0.1, 0.2, 0.3], method='Euler'), [0.1, 0.2, 0.299])
    assert _approx(sin_.evalf([0.1, 0.2, 0.3]), [0.0998333333333333, 0.198669165277778, 0.295519962530663])
    d = sin_.evalf([0.1, 0.2], derivatives=True)
    assert _approx(d[1], [0.198669165277778, 0.980066597239583])
    assert _approx(exp_.evalf([0.1j, 0.2j]), [0.995004166666667+0.0998333333333333j,
                                              0.980066597239583+0.198669165277778j])
    # A single endpoint still integrates with steps of length h.
    assert abs(float(exp_.evalf(1)) - 2.71828169265634) < 1e-12


# --- printing exposed by the probes

def test_vector_and_permutation_str_follow_reference():
    C = CoordSys3D('C')
    B = CoordSys3D('B')
    x_, y_ = symbols('x_ y_')
    assert str(Rational(7, 3)*C.i + C.j) == '7/3*C.i + C.j'
    assert str(-C.i + C.j) == '(-1)*C.i + C.j'
    assert str((x_ + y_)*C.k) == '(x_ + y_)*C.k'
    assert str(-x_*C.j) == '(-x_)*C.j'
    assert str(2.0*C.i) == '2.0*C.i'
    assert str(B.k + C.i + B.i) == 'B.i + B.k + C.i'
    assert str(2*(C.i | C.j) + (C.j | C.i)) == '2*(C.i|C.j) + (C.j|C.i)'
    assert str((C.i | C.j)*0) == '0'
    assert repr(Permutation(0, 1, 2, 3)) == 'Permutation(0, 1, 2, 3)'
    assert repr(Permutation([1, 0, 2])) == 'Permutation(2)(0, 1)'
    assert repr(Permutation([])) == 'Permutation()'
    assert repr([Permutation(1, 2, size=5)]) == '[Permutation(4)(1, 2)]'
