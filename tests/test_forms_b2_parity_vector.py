"""Batch b2 form parity for sympy.vector: expected strings come from SymPy 1.14."""
from symox import Symbol, sin
from symox.vector import (CoordSys3D, Del, curl, divergence, gradient,
                          laplacian, vector_integrate)


C = CoordSys3D('C')
R = CoordSys3D('R')
S = C.create_new('S', transformation='spherical')
q = Symbol('q')
D = C.orient_new_axis('D', q, C.k)
delop = Del()


def test_del_cross_holds_sympy_curl_form():
    v = C.x*C.y*C.z*(C.i + C.j + C.k)
    assert str(delop.cross(v, doit=True)) == (
        '(-C.x*C.y + C.x*C.z)*C.i + (C.x*C.y - C.y*C.z)*C.j + (-C.x*C.z + C.y*C.z)*C.k')


def test_del_dot_without_doit_is_held_derivative():
    assert str(delop.dot(C.x*C.i)) == 'Derivative(C.x, C.x)'


def test_unevaluated_operators():
    assert str(gradient(C.x*C.y, doit=False)) == (
        '(Derivative(C.x*C.y, C.x))*C.i + (Derivative(C.x*C.y, C.y))*C.j'
        ' + (Derivative(C.x*C.y, C.z))*C.k')
    assert str(curl(C.x*C.j, doit=False)) == (
        '(Derivative(0, C.y) - Derivative(C.x, C.z))*C.i'
        ' + (-Derivative(0, C.x) + Derivative(0, C.z))*C.j'
        ' + (-Derivative(0, C.y) + Derivative(C.x, C.x))*C.k')


def test_laplacian_is_expanded_sum_of_second_derivatives():
    assert str(laplacian(R.x**2*R.y**5*R.z)) == '20*R.x**2*R.y**3*R.z + 2*R.y**5*R.z'


def test_curvilinear_operators_are_not_simplified():
    f = S.r**2*sin(S.theta)
    assert str(gradient(f)) == '(2*S.r*sin(S.theta))*S.i + (S.r*cos(S.theta))*S.j'
    assert str(divergence(S.r*S.i + S.theta*S.j)) == (
        '3 + (S.r*S.theta*cos(S.theta) + S.r*sin(S.theta))/(S.r**2*sin(S.theta))')
    assert str(curl(S.r*S.k + S.phi*S.j)) == (
        '((S.r**2*cos(S.theta) - S.r)/(S.r**2*sin(S.theta)))*S.i + (-2)*S.j + S.phi/S.r*S.k')
    assert str(laplacian(f)) == (
        '6*sin(S.theta) + (-S.r**2*sin(S.theta)**2 + S.r**2*cos(S.theta)**2)'
        '/(S.r**2*sin(S.theta))')
    assert str(laplacian(S.r*S.i)) == '0'


def test_multi_system_operators_split_by_system():
    assert str(gradient(C.x*D.y)) == 'D.y*C.i + C.x*D.j'
    assert str(gradient(C.x*R.y, doit=False)) == 'R.y*C.i + C.x*R.j'
    assert str(gradient(C.x + R.y)) == 'C.i + R.j'
    assert str(divergence(C.x*C.i + D.y*D.j)) == '2'
    assert str(divergence(D.x*C.i)) == 'cos(q)'
    assert str(divergence(C.x*C.i + R.y*R.j)) == '2'
    assert str(laplacian(C.x*D.y)) == '-2*sin(q)'


def test_vector_integrate_holds_multi_term_fields():
    v = C.x*C.i - C.y*C.k
    assert str(vector_integrate(v, C.x)) == (
        '(Integral(C.x, C.x))*C.i + (Integral(-C.y, C.x))*C.k')
    assert str(v.integrate(C.x)) == '(Integral(C.x, C.x))*C.i + (Integral(-C.y, C.x))*C.k'
    assert str(v.integrate(C.x).doit()) == 'C.x**2/2*C.i + (-C.x*C.y)*C.k'
    assert str(vector_integrate(v, (C.x, 0, 1))) == (
        '(Integral(C.x, (C.x, 0, 1)))*C.i + (Integral(-C.y, (C.x, 0, 1)))*C.k')
    assert str(vector_integrate(C.x*C.i, C.y)) == 'C.x*C.y*C.i'
    assert str(vector_integrate(3*C.i, C.x)) == '3*C.x*C.i'
    assert str(vector_integrate(C.x, C.x)) == 'C.x**2/2'
