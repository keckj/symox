"""Wave-4 wrong-result regressions in mechanics, vector, control and continuum mechanics.

Expected values are the SymPy 1.14 outputs of the same calls.
"""
import pytest

import symox as S
from symox import Symbol, symbols, Matrix, sin, cos, tan, exp, Rational, pi, sqrt, simplify, oo
from symox.physics.mechanics import (Body, RigidBody, ReferenceFrame, Point, dynamicsymbols, msubs, inertia,
                                     WrappingCylinder, WrappingSphere, WrappingPathway, ForceActuator)
from symox.physics.vector import CoordinateSym
from symox.physics.control import (TransferFunction, TransferFunctionMatrix, StateSpace, phase_margin,
                                   gain_margin)
from symox.physics.control.control_plots import pole_zero_numerical_data
from symox.physics.continuum_mechanics.beam import Beam3D
from symox.physics.continuum_mechanics.cable import Cable

t = dynamicsymbols._t
s = Symbol('s')


def _body_setup():
    m, v, omega = symbols('m v omega')
    N, b = ReferenceFrame('N'), ReferenceFrame('b')
    b.set_ang_vel(N, omega*b.x)
    P = Point('P')
    P.set_vel(N, v*N.x)
    return m, v, omega, N, b, P


# ---- mechanics: bodies --------------------------------------------------------

def test_body_default_inertia_is_symbolic():
    m, v, omega, N, b, P = _body_setup()
    B = Body('B', masscenter=P, frame=b)
    assert str(B.kinetic_energy(N)) == 'B_ixx*omega**2/2 + B_mass*v**2/2'
    assert str(B.inertia.dyadic.to_matrix(b)) == \
        'Matrix([[B_ixx, B_ixy, B_izx], [B_ixy, B_iyy, B_iyz], [B_izx, B_iyz, B_izz]])'
    assert str(B.angular_momentum(P, N)) == 'B_ixx*omega*b.x + B_ixy*omega*b.y + B_izx*omega*b.z'
    assert B.mass == Symbol('B_mass')


def test_body_with_mass_only_is_a_particle():
    m, v, omega, N, b, P = _body_setup()
    assert str(Body('B', masscenter=P, frame=b, mass=m).kinetic_energy(N)) == 'm*v**2/2'
    assert str(Body('B', masscenter=P, mass=m).kinetic_energy(N)) == 'm*v**2/2'


def test_rigidbody_default_inertia_is_symbolic():
    N = ReferenceFrame('N')
    R = RigidBody('R')
    assert str(R.inertia.dyadic.to_matrix(R.frame)) == \
        'Matrix([[R_ixx, R_ixy, R_izx], [R_ixy, R_iyy, R_iyz], [R_izx, R_iyz, R_izz]])'
    assert R.mass == Symbol('R_mass')


# ---- dynamicsymbols ------------------------------------------------------------

def test_dynamicsymbols_keep_assumptions():
    q4, q5 = dynamicsymbols('q4,q5', commutative=False)
    assert bool(q4*q5 != q5*q4) is True
    assert q4.is_commutative is False
    q = dynamicsymbols('q', real=True)
    assert (q.is_real, q.is_positive) == (True, None)
    assert str(S.conjugate(q)) == 'q(t)'
    p = dynamicsymbols('p', positive=True)
    assert (p.is_positive, p.is_real) == (True, True)
    assert str(sqrt(p**2)) == 'p(t)'
    assert dynamicsymbols('n', integer=True).is_integer is True
    assert str(dynamicsymbols('q', 2, real=True)) == 'Derivative(q(t), (t, 2))'


# ---- msubs ----------------------------------------------------------------------

def test_msubs_ignores_derivative_interiors():
    x = dynamicsymbols('x')
    assert str(msubs(x.diff() + x, {x: 1})) == 'Derivative(x(t), t) + 1'
    assert str(msubs(x.diff(t, 2) + x.diff(), {x.diff(): 3})) == 'Derivative(x(t), (t, 2)) + 3'
    assert str(msubs(Matrix([x, x.diff()]), {x: 2})) == 'Matrix([[2], [Derivative(x(t), t)]])'
    assert str(msubs(cos(x)*S.Derivative(sin(x), t), {x: 0})) == 'Derivative(sin(x(t)), t)'
    y, z = dynamicsymbols('y z')
    assert msubs(x.diff() + x + y + z, {x: 1, y: 2}, {z: 3, x.diff(): 4}) == 10
    N = ReferenceFrame('N')
    assert str(msubs(x*N.x + x.diff()*N.y, {x: 2})) == '2*N.x + Derivative(x(t), t)*N.y'


def test_msubs_smart():
    x = dynamicsymbols('x')
    assert msubs(sin(x)/tan(x), {x: 0}, smart=True) == 1
    assert str(msubs(Matrix([sin(x)/tan(x)]), {x: 0}, smart=True)) == 'Matrix([[1]])'
    assert str(msubs(tan(x), {x: 1}, smart=True)) == 'sin(1)/cos(1)'
    assert msubs(sin(x)/x, {x: 0}) is S.nan


# ---- CoordinateSym ---------------------------------------------------------------

def test_coordinatesym_compares_by_frame_and_index():
    A, B = ReferenceFrame('A'), ReferenceFrame('B')
    assert (CoordinateSym('a_y', A, 1) == A[1]) is True
    assert (CoordinateSym('A_y', A, 1) == A[1]) is True
    assert (CoordinateSym('A_x', A, 0) == A[1]) is False
    assert (CoordinateSym('A_x', A, 0) == CoordinateSym('A_x', B, 0)) is False
    assert hash(CoordinateSym('A_x', A, 0)) == hash(A[0])
    with pytest.raises(ValueError):
        CoordinateSym('A_x', A, 3)


# ---- wrapping geometry -------------------------------------------------------------

def test_wrapping_cylinder_geodesic_length_has_no_mod():
    N = ReferenceFrame('N')
    r = symbols('r')
    pO, p1, p2 = Point('pO'), Point('p1'), Point('p2')
    q = dynamicsymbols('q')
    c = WrappingCylinder(r, pO, N.x)
    p1.set_pos(pO, N.x + r*N.y)
    p2.set_pos(pO, r*(cos(q)*N.y + sin(q)*N.z).normalize())
    length = c.geodesic_length(p1, p2)
    assert str(length) == 'sqrt(r**2*q(t)**2 + 1)'
    assert abs(float(length.subs(q, -0.5).subs(r, 2).evalf()) - 1.41421356237310) < 1e-12


def test_wrapping_cylinder_numeric_angles():
    N = ReferenceFrame('N')
    r = symbols('r')
    pO, p1, p2 = Point('pO'), Point('p1'), Point('p2')
    c = WrappingCylinder(r, pO, N.z)
    p1.set_pos(pO, r*N.x)
    p2.set_pos(pO, r*N.y + 2*N.z)
    assert str(c.geodesic_length(p1, p2)) == 'sqrt(pi**2*r**2/4 + 4)'
    p2.set_pos(pO, -r*N.y)
    assert str(c.geodesic_length(p1, p2)) == '3*pi*sqrt(r**2)/2'
    p2.set_pos(pO, 2*r*N.y)
    with pytest.raises(ValueError):
        c.geodesic_length(p1, p2)
    assert repr(c) == 'WrappingCylinder(radius=r, point=pO, axis=N.z)'


def test_wrapping_pathway_to_loads_signs():
    N = ReferenceFrame('N')
    r = symbols('r', positive=True)
    F = symbols('F')
    pO, pA, pB = Point('pO'), Point('pA'), Point('pB')
    q = dynamicsymbols('q')
    pA.set_pos(pO, r*N.x)
    pB.set_pos(pO, r*(cos(q)*N.x + sin(q)*N.y))
    pathway = WrappingPathway(pA, pB, WrappingCylinder(r, pO, N.z))
    loads = pathway.to_loads(F)
    assert str(loads[0].vector.simplify()) == 'F*N.y'
    assert str(loads[1].vector.simplify()) == 'F*sin(q(t))*N.x - F*cos(q(t))*N.y'
    assert str(ForceActuator(F, pathway).to_loads()[0].vector.simplify()) == 'F*N.y'
    assert str(pathway.length) == 'r*sqrt(q(t)**2)'


def test_wrapping_sphere():
    N = ReferenceFrame('N')
    r = symbols('r')
    pO, p1, p2 = Point('pO'), Point('p1'), Point('p2')
    c = WrappingSphere(r, pO)
    p1.set_pos(pO, r*N.x)
    p2.set_pos(pO, r*N.y)
    assert c.geodesic_length(p1, p2) == pi*r/2
    p2.set_pos(pO, -r*N.x)
    with pytest.raises(ValueError):
        c.geodesic_end_vectors(p1, p2)
    assert repr(c) == 'WrappingSphere(radius=r, point=pO)'


# ---- control ------------------------------------------------------------------------

def test_phase_margin():
    assert abs(float(phase_margin(TransferFunction(s**3, s**2 + 5*s, s)).n()) + 25.1783920627277) < 1e-10
    assert phase_margin(TransferFunction(1, s + 1, s)) == -180
    assert phase_margin(TransferFunction(2, s + 1, s)) == 120
    assert phase_margin(TransferFunction(5, s, s)) == 90
    assert phase_margin(TransferFunction(1, s**2 + 1, s)) == 0
    assert abs(float(phase_margin(TransferFunction(100, (s + 1)**3, s)).n(10)) + 52.67536215) < 1e-7
    with pytest.raises(ValueError):
        phase_margin(TransferFunction(Symbol('a'), s + 1, s))


def test_gain_margin():
    assert gain_margin(TransferFunction(s**3, s**2 + 5*s, s)) is oo
    assert gain_margin(TransferFunction(1, s + 1, s)) is oo
    assert gain_margin(TransferFunction(5, s, s)) is oo
    assert abs(float(gain_margin(TransferFunction(10, s*(s + 1)*(s + 5), s)).n(10)) - 9.542425094) < 1e-7
    assert abs(float(gain_margin(TransferFunction(100, (s + 1)**3, s)).n(10)) + 21.93820026) < 1e-7


def test_pole_zero_numerical_data_returns_zeros_then_poles():
    zeros, poles = pole_zero_numerical_data(TransferFunction(s**2 + 1, s**4 + 4*s**3 + 6*s**2 + 5*s + 2, s))
    assert zeros == [-1j, 1j]
    assert poles[:2] == [-2.0, -1.0]
    assert all(type(v) is float for v in poles[:2])
    assert pole_zero_numerical_data(TransferFunction(1, s + 2, s)) == ([], [-2.0])
    assert pole_zero_numerical_data(TransferFunction(s - 1, (s + 2)*(s + 3), s)) == ([1.0], [-3.0, -2.0])


def test_statespace_dsolve_returns_output():
    ss = StateSpace(Matrix([[-2, 0], [1, -1]]), Matrix([[1], [0]]), Matrix([[2, 1]]))
    y = ss.dsolve(input_vector=Matrix([5]), initial_conditions=Matrix([0, 0]))
    assert y.shape == (1, 1)
    assert simplify(y[0] - (Rational(15, 2) - 5*exp(-t) - 5*exp(-2*t)/2)) == 0
    assert simplify(ss.dsolve(initial_conditions=Matrix([1, 2]))[0] - (3*exp(t) + 1)*exp(-2*t)) == 0
    assert ss.dsolve() == Matrix([[0]])
    ss2 = StateSpace(Matrix([[0, 1], [-2, -3]]), Matrix([[0], [1]]), Matrix([[1, 0]]), Matrix([[1]]))
    y2 = ss2.dsolve(input_vector=Matrix([1]), initial_conditions=Matrix([1, 0]))
    assert simplify(y2[0] - (Rational(3, 2) + exp(-t) - exp(-2*t)/2)) == 0
    with pytest.raises(S.ShapeError):
        ss2.dsolve(input_vector=Matrix([1, 2]))
    with pytest.raises(ValueError):
        ss2.dsolve(var=Symbol('t') + 1)


def test_transfer_function_matrix_indexing():
    TF = TransferFunction
    m = TransferFunctionMatrix([[TF(5, s, s), TF(5*s, s**2 + 2, s)], [TF(5, s*(s**2 + 2), s), TF(5, 1, s)]])
    assert str(m[:, 0]) == \
        'TransferFunctionMatrix(((TransferFunction(5, s, s),), (TransferFunction(5, s*(s**2 + 2), s),)))'
    assert m[:, 0].shape == (2, 1)
    assert str(m[0, :]) == \
        'TransferFunctionMatrix(((TransferFunction(5, s, s), TransferFunction(5*s, s**2 + 2, s)),))'
    assert str(m[1:, 1:]) == 'TransferFunctionMatrix(((TransferFunction(5, 1, s),),))'
    assert str(m[0]) == 'TransferFunction(5, s, s)'
    assert str(m[3]) == 'TransferFunction(5, 1, s)'
    assert str(m[-1, 0]) == 'TransferFunction(5, s*(s**2 + 2), s)'
    with pytest.raises(IndexError):
        m[2, 0]


# ---- continuum mechanics ----------------------------------------------------------------

def _beam3d():
    l, E, G, I, A, q, m, x = symbols('l E G I A q m x')
    return (l, E, G, I, A, q, m, x), Beam3D(l, E, G, I, A)


def test_beam3d_bending_moment_and_shear():
    (l, E, G, I, A, q, m, x), b = _beam3d()
    b.apply_load(q, 0, 0, dir='y')
    b.apply_moment_load(m, 0, -1, dir='z')
    assert str(b.bending_moment()) == '[0, 0, -m*x + q*x**2/2]'
    assert str(b.shear_force()) == '[0, -q*x, 0]'
    assert b.load_vector == [0, q, 0]
    assert b.polar_moment() == 2*I


def test_beam3d_deflection_constants_are_solved():
    (l, E, G, I, A, q, m, x), b = _beam3d()
    b.apply_load(q, 0, 0, dir='y')
    b.apply_moment_load(m, 0, -1, dir='z')
    b.solve_slope_deflection()
    defl = b.deflection()
    assert not any(c in defl[1].free_symbols for c in symbols('C1 C2 C3 C4 C13 C14 C_i'))
    assert simplify(defl[1].subs(x, l/2) - (l**4*q/(384*E*I) + l**2*q/(8*A*G))) == 0
    (l, E, G, I, A, q, m, x), b = _beam3d()
    b.apply_load(q, 0, 0, dir='z')
    b.solve_slope_deflection()
    expected = l**2*q*x**2/(24*E*I) - l*q*x**3/(12*E*I) + q*x**4/(24*E*I) + l*q*x/(2*A*G) - q*x**2/(2*A*G)
    assert simplify(b.deflection()[2] - expected) == 0
    (l, E, G, I, A, q, m, x), b = _beam3d()
    b.apply_load(q, 0, 0, dir='x')
    b.solve_slope_deflection()
    assert simplify(b.deflection()[0] - (l*q*x/(2*A*E) - q*x**2/(2*A*E))) == 0


def test_cable_apply_length():
    c = Cable(('A', 0, 10), ('B', 10, 10))
    assert c.length == 0
    c.apply_length(20)
    assert c.length == 20
    with pytest.raises(ValueError):
        Cable(('A', 0, 10), ('B', 10, 10)).apply_length(5)
    assert str(Cable(('A', 0, 10), ('B', 10, 10)).reaction_loads) == '{R_A_x: 0, R_A_y: 0, R_B_x: 0, R_B_y: 0}'
    c = Cable(('A', 0, 10), ('B', 10, 10))
    c.apply_length(20)
    c.apply_load(-1, ('P', 4, 5, 10, 270))
    c.solve()
    assert c.length == sqrt(41) + sqrt(61)


def test_vector_and_dyadic_print_like_sympystr():
    N = ReferenceFrame('N')
    a, b = symbols('a b')
    B = N.orientnew('B', 'Axis', [a, N.z])
    assert str(N.x) == 'N.x'
    assert str(-N.y) == '- N.y'
    assert str((a + b)*N.x - a*N.y + 2*N.z) == '(a + b)*N.x - a*N.y + 2*N.z'
    assert str(-(a + b)*N.x) == '(-a - b)*N.x'
    assert str(B.x + N.x) == 'N.x + B.x'
    assert str(N.x - N.x) == '0'
    assert str(-S.physics.vector.outer(N.x, N.y) + (a + b)*S.physics.vector.outer(N.z, N.z)) == '- (N.x|N.y) + (a + b)*(N.z|N.z)'


def test_body_particle_dispatch():
    N = ReferenceFrame('N')
    m, v = symbols('m v')
    P = Body('P', mass=m)
    P.masscenter.set_vel(N, v*N.x)
    assert P.is_rigidbody is False
    assert P.kinetic_energy(N) == m*v**2/2
    assert Body('Q').is_rigidbody is True
