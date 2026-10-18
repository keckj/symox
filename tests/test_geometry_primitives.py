"""Exact geometry results are compared with the independent reference CAS."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


def back(value):
    if hasattr(value, 'to_sympy'):
        return value.to_sympy()
    if isinstance(value, list):
        return [back(item) for item in value]
    if isinstance(value, tuple):
        return tuple(back(item) for item in value)
    if isinstance(value, dict):
        return {back(key):back(item) for key, item in value.items()}
    if isinstance(value, set):
        return {back(item) for item in value}
    return value


def equivalent(actual, expected):
    actual = back(actual)
    if actual == expected:
        return
    if isinstance(expected, (list, tuple)):
        assert len(actual) == len(expected)
        remaining = list(expected)
        for item in actual:
            for i, candidate in enumerate(remaining):
                if item == candidate or hasattr(item, 'equals') and item.equals(candidate):
                    remaining.pop(i)
                    break
            else:
                pytest.fail(f'Unmatched geometry result {item}; expected {expected}')
    elif isinstance(expected, sp.Expr) and isinstance(actual, sp.Expr):
        assert sp.simplify(actual-expected) == 0
    elif hasattr(expected, 'equals'):
        assert actual == expected or actual.equals(expected)
    else:
        assert actual == expected


def expression(module, source):
    namespace = vars(module) | {'x':module.Symbol('x', real=True),
                                'y':module.Symbol('y', real=True),
                                'a':module.Symbol('a', positive=True)}
    return eval(source, namespace)


@pytest.mark.parametrize('source', [
    'Point(1,2)', 'Point3D(1,2,3)', 'Point([1,2],dim=3)', 'Point2D(1,2,0)',
    'Point(1,2).bounds', 'Segment((2,-1),(-3,4)).bounds',
    'Line((2,-1),(-3,4)).bounds', 'Circle((1,2),3).bounds',
    'Ellipse((1,2),3,2).bounds', 'Polygon((0,0),(2,0),(0,3)).bounds',
    'Point(.1,.2)', 'Point(1,2)+Point(3,4)', 'Point(1,2)-Point(3,4)',
    '-Point(1,2)', 'Point(1,2)*3', 'Point(1,2)/3', 'Point(3,4).length',
    'Point(3,4).unit', 'Point(1,2).dot(Point(3,4))',
    'Point(1,2).midpoint(Point(4,6))',
    'Point(1,2).distance(Point(4,6))', 'Point(1,2).taxicab_distance(Point(4,6))',
    'Point(1,2).translate(3,-1)', 'Point3D(1,2,3).translate(3,-1,2)',
    'Point(1,2).rotate(pi/2)', 'Point(1,2).rotate(pi/2,Point(1,1))',
    'Point(1,2).scale(2,3)', 'Point(1,2).scale(2,3,pt=Point(1,1))',
    'Point.is_collinear((0,0),(1,1),(2,2))', 'Point.is_collinear((0,0),(1,1),(2,3))',
    'Point.is_concyclic((1,0),(0,1),(-1,0),(0,-1))',
    'Point3D.are_coplanar((0,0,0),(1,0,0),(0,1,0),(1,1,0))',
    'Point3D.are_coplanar((0,0,0),(1,0,0),(0,1,0),(1,1,1))',
    'Line((0,0),(1,1))', 'Line((0,0),slope=2)', 'Line((0,0),slope=oo)',
    'Line3D((0,0,0),direction_ratio=(1,2,3))', 'Segment((0,0),(0,0))',
    'Segment((0,0),(3,4)).length', 'Segment((0,0),(3,4)).midpoint',
    'Line((0,0),(2,2)).slope', 'Line((0,0),(0,2)).slope',
    'Line((0,0),(2,2)).equation(x,y)', 'Line((2,0),(2,2)).equation(x,y)',
    'Line((0,0),(1,1)).projection(Point(0,2))',
    'Segment((0,0),(1,0)).projection(Point(2,1))',
    'Line((0,0),(1,1)).distance(Point(0,2))',
    'Segment((0,0),(1,0)).distance(Point(-2,1))',
    'Segment((0,0),(1,0)).distance(Point(2,1))',
    'Ray((0,0),(1,0)).distance(Point(-2,1))',
    'Ray((0,0),(1,0)).distance(Point(2,1))',
    'Line((0,0),(1,1)).parallel_line(Point(0,2))',
    'Line((0,0),(1,1)).perpendicular_line(Point(0,2))',
    'Line((0,0),(1,0)).perpendicular_line(Point(2,0))',
    'Segment((0,0),(2,0)).perpendicular_bisector()',
    'Segment((0,0),(2,0)).perpendicular_bisector(Point(1,2))',
    'Line((0,0),(1,0)).angle_between(Line((0,0),(1,1)))',
    'Line((0,0),(1,0)).is_parallel(Line((1,1),(2,1)))',
    'Line((0,0),(1,0)).is_perpendicular(Line((0,0),(0,1)))',
    'Circle((0,0),2).area', 'Circle((0,0),2).circumference',
    'Circle((1,2),-2).equation(x,y)', 'Circle((1,0),(0,1),(-1,0))',
    'Circle((0,0),2).contains(Point(0,2))', 'Circle((0,0),2).encloses_point(Point(0,1))',
    'Circle((0,0),2).tangent_lines(Point(2,0))',
    'Circle((0,0),2).tangent_lines(Point(4,0))',
    'Ellipse((1,2),3,2).area', 'Ellipse((0,0),3,2).foci',
    'Ellipse((0,0),2,3).foci', 'Ellipse((0,0),3,2).eccentricity',
    'Ellipse((0,0),3,2).circumference', 'Ellipse((1,2),3,2).equation(x,y)',
    'Ellipse((0,0),3,2).tangent_lines(Point(3,0))',
    'Circle((1,2),2).scale(2,3)', 'Circle((1,2),2).scale(-1,1)',
    'Ellipse((1,2),3,2).rotate(pi/2)',
    'Polygon((0,0),(4,0),(4,3),(0,3)).area',
    'Polygon((0,0),(4,0),(4,3),(0,3)).perimeter',
    'Polygon((0,0),(4,0),(4,3),(0,3)).centroid',
    'Polygon((0,0),(4,0),(4,3),(0,3)).encloses_point(Point(2,1))',
    'Polygon((0,0),(4,0),(4,3),(0,3)).encloses_point(Point(0,1))',
    'Polygon((0,0),(4,0),(4,3),(0,3)).encloses_point(Point(8,1))',
    'Polygon((0,0),(4,0),(4,3),(0,3)).is_convex()',
    'Polygon((0,0),(4,0),(1,1),(0,3)).is_convex()',
    'Polygon((0,0),(1,0),(2,0),(2,2),(0,2))',
    'Triangle((0,0),(4,0),(0,3)).incenter',
    'Triangle((0,0),(4,0),(0,3)).inradius',
    'Triangle((0,0),(4,0),(0,3)).circumcenter',
    'Triangle((0,0),(4,0),(0,3)).circumradius',
    'Triangle((0,0),(4,0),(0,3)).orthocenter',
    'Triangle((0,0),(4,0),(0,3)).is_right()',
    'Triangle((0,0),(4,0),(0,3)).is_scalene()',
    'Triangle((0,0),(1,0),(Rational(1,2),sqrt(3)/2)).is_equilateral()',
    'RegularPolygon((0,0),2,4).area', 'RegularPolygon((0,0),2,4).apothem',
    'RegularPolygon((0,0),2,4).vertices',
    'Plane((0,0,0),normal_vector=(0,0,1)).distance(Point3D(1,2,3))',
    'Plane((0,0,0),normal_vector=(0,0,1)).projection(Point3D(1,2,3))',
    'Plane((0,0,0),(1,0,0),(0,1,0)).equation(x,y,Symbol("z",real=True))',
    'Parabola((0,1),Line((-1,-1),(1,-1))).vertex',
    'Parabola((0,1),Line((-1,-1),(1,-1))).focus',
    'Parabola((0,1),Line((-1,-1),(1,-1))).focal_length',
    'Parabola((0,1),Line((-1,-1),(1,-1))).equation(x,y)',
    'Curve((cos(x),sin(x)),(x,0,pi)).length',
    'rad(45)', 'deg(pi/3)', 'centroid(Point(0,0),Point(2,2))',
    'centroid(Segment((0,0),(2,0)),Segment((0,0),(0,4)))',
    'convex_hull((0,0),(2,0),(0,2),(1,1))',
    'convex_hull((0,0),(2,0),(0,2),(1,1),polygon=False)',
    'convex_hull((0,0),(1,0),(2,0),polygon=False)',
    'convex_hull((0,0),polygon=False)',
    'closest_points((2,0),(1,0),(0,0))',
    'farthest_points((2,0),(1,0),(0,0))',
    'idiff(x*x+y*y-1,y,x)', 'idiff(x*x+y*y-1,y,x,2)',
])
def test_exact_geometry_matches_reference(source):
    equivalent(expression(sx, source), expression(sp, source))


@pytest.mark.parametrize('left,right', [
    ('Line((0,0),(1,1))','Line((0,1),(1,0))'),
    ('Line((0,0),(1,1))','Line((0,1),(1,2))'),
    ('Line((0,0),(1,1))','Line((2,2),(3,3))'),
    ('Segment((0,0),(2,0))','Segment((1,0),(3,0))'),
    ('Segment((0,0),(1,0))','Segment((1,0),(2,0))'),
    ('Segment((0,0),(1,0))','Segment((2,0),(3,0))'),
    ('Segment((0,0),(2,0))','Segment((3,0),(1,0))'),
    ('Ray((0,0),(1,0))','Ray((2,0),(1,0))'),
    ('Ray((0,0),(1,0))','Ray((-2,0),(-1,0))'),
    ('Ray((0,0),(1,0))','Ray((2,0),(3,0))'),
    ('Line3D((0,0,0),(1,0,0))','Line3D((0,1,1),(0,2,1))'),
    ('Line3D((0,0,0),(1,0,0))','Line3D((0,1,0),(0,2,0))'),
    ('Circle((0,0),2)','Line((0,0),(1,0))'),
    ('Circle((0,0),2)','Line((0,2),(1,2))'),
    ('Circle((0,0),2)','Line((0,3),(1,3))'),
    ('Circle((0,0),2)','Segment((0,0),(3,0))'),
    ('Circle((0,0),2)','Ray((3,0),(4,0))'),
    ('Circle((0,0),2)','Circle((2,0),2)'),
    ('Circle((0,0),2)','Circle((4,0),2)'),
    ('Circle((0,0),2)','Circle((5,0),2)'),
    ('Circle((0,0),2)','Circle((0,0),1)'),
    ('Ellipse((0,0),3,2)','Line((0,0),(1,1))'),
    ('Polygon((0,0),(2,0),(2,2),(0,2))','Line((0,1),(1,1))'),
    ('Polygon((0,0),(2,0),(2,2),(0,2))','Line((0,0),(1,0))'),
    ('Polygon((0,0),(2,0),(2,2),(0,2))','Circle((0,0),1)'),
    ('Plane((0,0,0),normal_vector=(0,0,1))','Line3D((1,1,-1),(1,1,2))'),
    ('Plane((0,0,0),normal_vector=(0,0,1))','Plane((0,0,1),normal_vector=(0,0,1))'),
    ('Plane((0,0,0),normal_vector=(0,0,1))','Plane((0,0,0),normal_vector=(0,1,0))'),
    ('Parabola((0,1),Line((-1,-1),(1,-1)))','Line((-1,1),(1,1))'),
])
def test_exact_intersections_match_reference(left, right):
    native_left, native_right = expression(sx, left), expression(sx, right)
    reference_left, reference_right = expression(sp, left), expression(sp, right)
    equivalent(native_left.intersection(native_right), reference_left.intersection(reference_right))
    equivalent(native_right.intersection(native_left), reference_right.intersection(reference_left))


def test_symbolic_coordinates_and_transformations_preserve_identity():
    a = sx.Symbol('a', positive=True)
    x = sx.Symbol('x', real=True)
    p = sx.Point(x, a)
    assert p.free_symbols == {x, a}
    assert p.subs(x, 2) == sx.Point(2, a)
    assert p.xreplace({x:a}) == sx.Point(a, a)
    assert sx.Line((0,0),(a,a)).contains(sx.Point(2*a,2*a))
    assert sx.Circle((x,0),a).contains(sx.Point(x+a,0))
    assert sx.Polygon((0,0),(a,0),(0,a)).area == a*a/2
    with pytest.raises(AttributeError):
        p._args = (1, 2)


@pytest.mark.parametrize('source', ['Point(1,2)', 'Line((0,0),(1,1))', 'Circle((0,0),2)',
                                   'Triangle((0,0),(1,0),(0,1))', 'RegularPolygon((0,0),1,4)',
                                   'Plane((0,0,0),normal_vector=(0,0,1))'])
def test_geometry_pickles_and_hashes(source):
    value = expression(sx, source)
    restored = pickle.loads(pickle.dumps(value))
    assert restored == value
    assert hash(restored) == hash(value)
    assert {value:1}[restored] == 1


def test_ambiguous_real_intersections_fail_instead_of_inventing_points():
    a = sx.Symbol('a', real=True)
    with pytest.raises(NotImplementedError):
        sx.Circle((0,0),1).intersection(sx.Line((0,a),(1,a)))
    with pytest.raises(ValueError):
        sx.Line((0,0),(0,0))
    with pytest.raises(sx.GeometryError):
        sx.Circle((0,0),(1,1),(2,2))
    with pytest.raises(ValueError):
        sx.Point(sx.I,0)


def test_skew_segment_distance_checks_interior_closest_points():
    left = sx.Segment3D((-1,0,0),(1,0,0))
    right = sx.Segment3D((0,-1,1),(0,1,1))
    assert left.distance(right) == 1


def test_vector_and_circle_distance_extensions_and_coincident_signed_circles():
    assert sx.Point3D(1,2,3).cross(sx.Point3D(2,-1,0)) == sx.Point3D(3,6,-5)
    assert sx.Circle((0,0),2).distance(sx.Point(0,0)) == 2
    assert sx.Circle((0,0),2).distance(sx.Point(4,0)) == 2
    # Opposite signed radii define the same locus. SymPy 1.14 instead
    # returns two parametrized points for these coincident circles.
    circle = sx.Circle((0,0),2)
    assert circle.intersection(sx.Circle((0,0),-2)) == [circle]
    assert sx.Point.is_concyclic((1,0),(1,0),(0,1),(-1,0))


def test_geometry_has_no_sympy_runtime_dependency():
    source = '''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
assert Circle((0,0),2).intersection(Line((0,0),(1,0)))==[Point(-2,0),Point(2,0)]
assert Triangle((0,0),(4,0),(0,3)).incenter==Point(1,1)
assert Plane((0,0,0),normal_vector=(0,0,1)).projection(Point3D(1,2,3))==Point3D(1,2,0)
x,y=symbols('x y',real=True)
assert idiff(x*x+y*y-1,y,x)==-x/y
'''
    subprocess.run([sys.executable, '-c', source], check=True)


@pytest.mark.parametrize('dimension',[2,3])
def test_seeded_exact_linear_intersections_across_entity_kinds(dimension):
    import random
    generator = random.Random(20261018+dimension)
    for _ in range(60):
        points = [tuple(generator.randint(-3,3) for _ in range(dimension)) for _ in range(4)]
        if points[0] == points[1] or points[2] == points[3]:
            continue
        first, second = [generator.choice(['Line','Ray','Segment']) for _ in range(2)]
        actual = getattr(sx,first)(*points[:2]).intersection(getattr(sx,second)(*points[2:]))
        expected = getattr(sp,first)(*points[:2]).intersection(getattr(sp,second)(*points[2:]))
        equivalent(actual,expected)


def test_similarity_extensions_check_angles_as_well_as_edge_lengths():
    rectangle = sx.Polygon((0,0),(2,0),(2,1),(0,1))
    transformed = sx.Polygon((1,1),(1,5),(-1,5),(-1,1))
    square = sx.Polygon((0,0),(1,0),(1,1),(0,1))
    half = sx.Rational(1,2)
    rhombus = sx.Polygon((0,0),(1,0),(1+half,sx.sqrt(3)/2),(half,sx.sqrt(3)/2))
    assert sx.are_similar(rectangle,transformed)
    assert not sx.are_similar(rectangle,square)
    assert not sx.are_similar(square,rhombus)
    assert sx.are_similar(sx.Circle((0,0),1),sx.Circle((1,2),3))


def test_implicit_differentiation_keeps_other_dependent_coordinates_distinct():
    x,y,z,w = sx.symbols('x y z w')
    X,Y,Z,W = sp.symbols('x y z w')
    actual = sx.idiff(x*y+z*w,(y,z,w),x)
    expected = sp.idiff(X*Y+Z*W,(Y,Z,W),X)
    assert sp.simplify(actual.to_sympy()-expected)==0
