"""Reference SymPy 1.14 doctest results for the category-c batch 4 callables."""
from symox import (Symbol, symbols, Piecewise, pi, N, sqrt, sin, S, Matrix,
                   DiracDelta, Heaviside, Integral, integrate, Eq, oo)
from symox.abc import x, y, z, a, s


def test_line_from_equation():
    from symox import Line, Point
    assert str(Line(3*x + y + 18)) == 'Line2D(Point2D(0, -18), Point2D(1, -21))'
    assert str(Line(Eq(3*a + b_symbol(), 0), x=a, y=b_symbol())) == \
        'Line2D(Point2D(0, 0), Point2D(1, -3))'
    assert Line(Point(0, 0), Point(1, 1)).slope == 1


def b_symbol():
    return Symbol('b')


def test_ray_transfer_matrix_from_entries():
    from symox.physics.optics import RayTransferMatrix
    mat = RayTransferMatrix(1, 2, 3, 4)
    assert mat == Matrix([[1, 2], [3, 4]])
    assert RayTransferMatrix(Matrix([[1, 2], [3, 4]])) == mat
    assert mat.A == 1


def test_orient_new_accepts_symbolic_quaternion():
    from symox.vector import CoordSys3D, QuaternionOrienter
    q0, q1, q2, q3 = symbols('q0 q1 q2 q3')
    N3 = CoordSys3D('N')
    D = N3.orient_new('D', (QuaternionOrienter(q0, q1, q2, q3),))
    B = N3.orient_new_quaternion('B', q0, q1, q2, q3)
    assert B.rotation_matrix(N3)[0, 0] == q0**2 + q1**2 - q2**2 - q3**2
    assert D.rotation_matrix(N3)[0, 0] == q0**2 + q1**2 - q2**2 - q3**2


def test_gradient_of_constant():
    from symox.vector import CoordSys3D, Del, Vector
    C = CoordSys3D('C')
    delop = Del()
    assert delop.gradient(9) == Vector.zero
    assert str(delop(C.x*C.y*C.z).doit()) == 'C.y*C.z*C.i + C.x*C.z*C.j + C.x*C.y*C.k'


def test_implicit_region_multiplicity():
    from symox.vector import ImplicitRegion
    R = ImplicitRegion((x, y, z), x**2 + y**3 - z**4)
    assert str(R.singular_points()) == '{(0, 0, 0)}'
    assert R.multiplicity((0, 0, 0)) == 2


def test_wavefunction_call_and_prob():
    from symox.physics.quantum.state import Wavefunction
    xr = Symbol('x', real=True)
    g = Piecewise((0, xr < 0), (0, xr > 1), (sqrt(2)*sin(pi*xr), True))
    f = Wavefunction(g, xr)
    assert f.norm == 1
    assert f.is_normalized
    p = f.prob()
    assert p(0) == 0
    assert p(1) == 0
    assert p(0.5) == 2
    assert str(p(0.85)) == '2*sin(0.85*pi)**2'
    assert str(N(p(0.85))) == '0.412214747707527'
    xp, L = symbols('x,L', positive=True)
    n = symbols('n', integer=True, positive=True)
    f = Wavefunction(sqrt(2/L)*sin(n*pi*xp/L), (xp, 0, L))
    assert f.norm == 1
    assert f(L + 1) == 0
    assert str(f(L - 1)) == 'sqrt(2)*sin(pi*n*(L - 1)/L)/sqrt(L)'


def test_int_qubit_from_qubit():
    from symox.physics.quantum.qubit import IntQubit, Qubit
    q = IntQubit(Qubit('101'))
    assert str(q) == '|5>'
    assert q.as_int() == 5
    assert q.nqubits == 3
    assert q.qubit_values == (1, 0, 1)
    assert str(Qubit(q)) == '|101>'


def test_integrate_result_and_represent():
    from symox.physics.quantum.represent import integrate_result, represent
    from symox.physics.quantum.cartesian import XOp, XKet
    x_1, x_2 = symbols('x_1 x_2')
    expr = x*DiracDelta(x - x_1)*DiracDelta(x_1 - x_2)
    assert integrate_result(XOp()*XKet(), expr) == expr
    assert str(integrate_result(XOp()*XKet(), expr, unities=[1])) == 'x*DiracDelta(x - x_2)'
    assert str(represent(XOp()*XKet())) == 'x*DiracDelta(x - x_2)'


def test_delta_integration_guard():
    assert integrate(x*DiracDelta(x - a), x) == a*Heaviside(x - a)
    # The reference returns zoo*Heaviside(x); the delta stage refuses
    # non-finite results and leaves the integral unevaluated.
    assert isinstance(integrate(DiracDelta(x)/x, x), Integral)
    # Non-unit slope: the stage must not give the wrong a**2/4.
    assert integrate(x**2*DiracDelta(2*x - a), x) != a**2/4


def test_spin_state_rewrite():
    from symox.physics.quantum.spin import JzKet, JxKet, JzKetCoupled
    from symox.physics.quantum.tensorproduct import TensorProduct
    assert str(JzKet(1, 1).rewrite('Jx')) == '|1,-1>/2 - sqrt(2)*|1,0>/2 + |1,1>/2'
    assert str(TensorProduct(JzKet(1, 1), JxKet(1, 1)).rewrite('Jz')) == \
        '|1,1>x|1,-1>/2 + sqrt(2)*|1,1>x|1,0>/2 + |1,1>x|1,1>/2'
    assert str(JzKetCoupled(1, 1, (1, 1)).rewrite('Jx')) == \
        '|1,-1,j1=1,j2=1>/2 - sqrt(2)*|1,0,j1=1,j2=1>/2 + |1,1,j1=1,j2=1>/2'
    assert str(JzKetCoupled(1, 0, (1, 1)).rewrite('Jz', coupled=False)) == \
        '-sqrt(2)*|1,-1>x|1,1>/2 + sqrt(2)*|1,1>x|1,-1>/2'


def test_phase_margin_decides_sign():
    from symox.physics.control import TransferFunction, phase_margin
    pm = phase_margin(TransferFunction(1, s**3 + 2*s**2 + s, s))
    assert str(pm.n()) == '21.3863897518751'
    assert str(phase_margin(TransferFunction(s + 1, s**2 + s, s)).n()) == '90.0000000000000'
