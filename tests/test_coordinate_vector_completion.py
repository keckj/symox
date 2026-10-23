"""Coordinate calculus, frame conversion and geometric integral identities."""
import subprocess
import sys

import pytest
import symox as s
from symox.vector import *


def test_vector_products_dyadics_and_scalar_dispatch():
    C = CoordSys3D('C')
    a,b = C.i+2*C.j+3*C.k, 3*C.i-C.j+C.k
    assert dot(a,b) == 4
    assert cross(C.i,C.j) == C.k
    assert (a.cross(b)).dot(a) == 0
    assert (a.cross(b)).dot(b) == 0
    assert a.outer(b).to_matrix(C) == a.to_matrix(C)*b.to_matrix(C).T
    assert a.dot(b.outer(C.k)) == 4*C.k
    assert a.outer(b).dot(C.i) == 3*a
    assert (C.i.outer(C.j)).T == C.j.outer(C.i)
    assert (C.i|C.j).dot(C.j|C.k) == C.i|C.k
    assert (C.i|C.j).cross(C.k) == C.i|C.i
    assert C.k.cross(C.i|C.j) == C.j|C.j
    assert (s.Symbol('a')*C.i).components == {C.i:s.Symbol('a')}
    assert VectorAdd(C.i,2*C.j) == VectorMul(2,C.j)+C.i
    assert DyadicAdd(C.i|C.j,C.j|C.k) == DyadicMul(1,C.i|C.j)+(C.j|C.k)
    assert C.i.kind is VectorKind(s.NumberKind)
    assert isinstance(C.x,BaseScalar)
    assert C.i.projection(a) == C.i
    assert C.i.normalize() == C.i
    assert matrix_to_vector(s.Matrix([1,2,3]),C) == a
    assert Cross(C.i,C.j).doit() == C.k
    assert Dot(a,b).doit() == 4
    assert not VectorZero() and not DyadicZero()


def test_cartesian_differential_identities_and_potential():
    C = CoordSys3D('C')
    x,y,z = C.base_scalars()
    potential = x**2*y+y*z**2+s.sin(x*z)
    field = gradient(potential)
    assert not curl(field).simplify()
    assert s.simplify(divergence(curl(x*y*C.i+y*z*C.j+z*x*C.k))) == 0
    assert s.simplify(laplacian(x*x+y*y+z*z)-6) == 0
    assert is_conservative(field)
    assert is_solenoidal(curl(x*y*C.i+y*z*C.j+z*x*C.k))
    assert s.simplify(scalar_potential(field,C)-potential) == 0
    assert directional_derivative(x*y,C.i+2*C.j) == y+2*x
    assert directional_derivative(x*C.i+y*C.j,C.i+2*C.j) == C.i+2*C.j
    assert s.diff(x*x*C.i,x) == 2*x*C.i
    assert s.integrate(x*C.i,x) == x*x/2*C.i
    assert gradient(3) == Vector.zero
    assert divergence(Vector.zero) == 0


@pytest.mark.parametrize('chart',['cylindrical','spherical'])
def test_orthogonal_chart_identities(chart):
    C = CoordSys3D('Q',transformation=chart)
    r,t,p = C.base_scalars()
    expected_div = 2 if chart == 'cylindrical' else 3
    assert divergence(r*C.i) == expected_div
    assert gradient(r*r) == 2*r*C.i
    assert laplacian(r*r) == 2*expected_div
    potential = r**2*s.sin(t)+p**2
    assert not curl(gradient(potential)).simplify()
    vector = r*t*C.i+r*p*C.j+t*p*C.k
    assert s.simplify(divergence(curl(vector))) == 0
    assert (directional_derivative(r*C.i,C.i)-C.i).simplify() == Vector.zero


def test_cylindrical_moving_basis_derivative():
    C = CoordSys3D('P',transformation='cylindrical')
    assert (directional_derivative(C.i,C.j)-C.j/C.r).simplify() == Vector.zero
    assert (directional_derivative(C.j,C.j)+C.i/C.r).simplify() == Vector.zero
    assert curl(C.r*C.j) == 2*C.k


@pytest.mark.parametrize('order',['123','231','312','132','213','321','121','131','212','232','313','323'])
@pytest.mark.parametrize('orienter',[BodyOrienter,SpaceOrienter])
def test_orientation_conventions_against_sympy(order,orienter):
    import sympy as sp
    from sympy import vector as reference
    angles = (s.pi/3,s.pi/4,s.pi/6)
    actual = orienter(*angles,order).rotation_matrix().to_sympy()
    expected = getattr(reference,orienter.__name__)(*(a.to_sympy() for a in angles),order).rotation_matrix()
    assert all(sp.simplify(c) == 0 for c in actual-expected)


def test_rotation_translation_points_and_scalar_conversion():
    C = CoordSys3D('C')
    B = C.orient_new_axis('B',s.pi/2,C.k,location=2*C.i)
    assert express(B.i,C) == C.j
    assert express(B.j,C) == -C.i
    assert express(B.x,C,variables=True) == C.y
    assert express(B.y,C,variables=True) == 2-C.x
    assert express(C.x,B,variables=True) == 2-B.y
    assert B.position_wrt(C) == 2*C.i
    p = C.origin.locate_new('P',C.i+2*C.j)
    q = p.locate_new('Q',3*C.k)
    assert q.position_wrt(C.origin) == C.i+2*C.j+3*C.k
    assert p.position_wrt(q) == -3*C.k
    assert scalar_potential_difference(gradient(C.x**2+C.y),C,C.origin,p) == 3
    assert QuaternionOrienter(s.sqrt(2)/2,0,0,s.sqrt(2)/2).rotation_matrix() == AxisOrienter(s.pi/2,C.k).rotation_matrix(C)
    same_label = C.locate_new('B',3*C.i)
    assert same_label.x != B.x


def test_parented_curvilinear_chart_conversion():
    C = CoordSys3D('C')
    P = C.create_new('P','cylindrical')
    assert express(C.x,P,variables=True) == P.r*s.cos(P.theta)
    assert (express(P.i,C)-s.cos(P.theta)*C.i-s.sin(P.theta)*C.j).simplify() == Vector.zero
    assert P.transformation_to_parent() == (P.r*s.cos(P.theta),P.r*s.sin(P.theta),P.z)
    assert (directional_derivative(P.i,P.j)-P.j/P.r).simplify() == Vector.zero
    assert CoordSys3D('C').origin.position_wrt(C.origin) == Vector.zero


def test_custom_orthogonal_chart_and_held_operators():
    C = CoordSys3D('Q',transformation=lambda u,v,w:(2*u,3*v,4*w))
    assert C.lame_coefficients() == (2,3,4)
    assert gradient(C.x**2) == C.x*C.i
    assert divergence(C.x*C.i+C.y*C.j+C.z*C.k) == s.Rational(13,12)
    assert Gradient(C.x).doit() == C.i/2
    assert Divergence(C.x*C.i).doit() == s.Rational(1,2)
    assert Curl(C.x*C.j).doit() == C.k/2
    assert Laplacian(C.x**2).doit() == s.Rational(1,2)
    d = Del()
    assert d(C.x).doit() == C.i/2
    assert (d & (C.x*C.i)).doit() == s.Rational(1,2)
    assert (d ^ (C.x*C.j)).doit() == C.k/2
    assert ((C.i & d)(C.x)) == s.Rational(1,2)


def test_line_surface_volume_and_stokes_integrals():
    C = CoordSys3D('C')
    t,u,v = s.symbols('t u v',real=True)
    assert vector_integrate(1,ParametricRegion((3*t,4*t),(t,0,1))) == 5
    assert vector_integrate(C.x*C.i,ParametricRegion((t,t**2),(t,0,1))) == s.Rational(1,2)
    square = ParametricRegion((u,v,0),(u,0,2),(v,0,3))
    assert vector_integrate(1,square) == 6
    assert vector_integrate(C.k,square) == 6
    assert vector_integrate(C.k,ParametricRegion((v,u,0),(u,0,2),(v,0,3))) == -6
    volume = ParametricRegion((u,v,t),(u,0,2),(v,0,3),(t,0,4))
    assert vector_integrate(1,volume) == 24
    triangle = ParametricRegion((u,v,0),(u,0,v),(v,0,1))
    assert vector_integrate(1,triangle) == s.Rational(1,2)
    circle = ParametricRegion((s.cos(t),s.sin(t),0),(t,0,2*s.pi))
    field = -C.y/2*C.i+C.x/2*C.j
    assert vector_integrate(field,circle) == s.pi
    assert vector_integrate(C.x*C.i,(C.x,0,2)) == 2*C.i


def test_divergence_theorem_on_sphere():
    C = CoordSys3D('C')
    t,p = s.symbols('t p',real=True)
    sphere = ParametricRegion((2*s.sin(t)*s.cos(p),2*s.sin(t)*s.sin(p),2*s.cos(t)),(t,0,s.pi),(p,0,2*s.pi))
    field = C.x*C.i+C.y*C.j+C.z*C.k
    assert vector_integrate(field,sphere) == 32*s.pi
    assert vector_integrate(1,sphere) == 16*s.pi
    assert divergence(field)*s.Rational(4,3)*s.pi*8 == 32*s.pi


def test_geometry_and_implicit_regions():
    x,y,z = s.symbols('x y z',real=True)
    C = CoordSys3D('C')
    circle = s.Circle(s.Point(0,2),5)
    assert vector_integrate(C.x**2+C.y**2,circle) == 290*s.pi
    polygon = s.Polygon((0,0),(1,0),(1,1),(0,1))
    assert vector_integrate(-C.y*C.i,polygon) == 1
    region = ImplicitRegion((x,y),(x-2)**2+(y-1)**2-9)
    assert region.degree == 2
    assert region.singular_points() == s.EmptySet
    assert ImplicitRegion((x,y),s.Eq(x,y)).equation == x-y
    assert vector_integrate(1,region) == 6*s.pi
    for equation,variables in [(x*x+y*y-1,(x,y)),(x*x-4*y,(x,y)),(x*x+y*y+z*z-1,(x,y,z)),(x+2*y+3*z-1,(x,y,z))]:
        implicit = ImplicitRegion(variables,equation)
        point = implicit.regular_point()
        assert s.simplify(equation.subs(dict(zip(variables,point)),simultaneous=True)) == 0
        parameterization = implicit.rational_parametrization()
        assert s.simplify(equation.subs(dict(zip(variables,parameterization)),simultaneous=True)) == 0


def test_rejections_and_disconnected_frames():
    C = CoordSys3D('C')
    with pytest.raises(ValueError):
        CoordSys3D('Q',transformation=lambda u,v,w:(u+v,v,w))
    # Reference SymPy stores the rotation matrix as given (no orthogonality check).
    assert CoordSys3D('Q2',rotation_matrix=s.Matrix([[1,0,0],[0,1,0],[0,0,2]]))._rotation[2, 2] == 2
    with pytest.raises(TypeError):
        C.i*C.j
    with pytest.raises(TypeError):
        C.i+1
    with pytest.raises(ValueError):
        express(C.i,CoordSys3D('disconnected'))
    with pytest.raises(ValueError):
        scalar_potential(C.x*C.j,C)
    # Reference SymPy takes the quaternion as given (no unit-norm check).
    assert QuaternionOrienter(2,0,0,0).rotation_matrix() == 4*s.eye(3)


def test_vector_calculus_without_sympy():
    subprocess.run([sys.executable,'-c',r'''
import sys
class NoSympy:
 def find_spec(self,name,*args):
  if name == 'sympy' or name.startswith('sympy.'):
   raise RuntimeError('runtime reference dependency: '+name)
sys.meta_path.insert(0,NoSympy())
import symox as s
from symox.vector import *
C=CoordSys3D('C')
assert divergence(C.x*C.i+C.y*C.j+C.z*C.k)==3
assert curl(gradient(C.x*C.y))==Vector.zero
t=s.Symbol('t',real=True)
assert vector_integrate(1,ParametricRegion((3*t,4*t),(t,0,1)))==5
assert 'sympy' not in sys.modules
'''],check=True)
