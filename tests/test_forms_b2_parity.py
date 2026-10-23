"""Form parity with SymPy 1.14 for audit batch form-b2 (expected strings from SymPy)."""
import pytest

from symox import (Matrix, DiracDelta, erf2inv, factor, together, exp, sin, cos,
                   sqrt, root, symbols, oo, horner)
from symox.geometry import Ellipse, Point, Point3D, Line, Line3D
from symox.solvers.solvers import unrad
from symox.simplify.radsimp import rad_rationalize, split_surds

x, y = symbols("x y")


CASES = [
    (lambda: Ellipse(Point(0, 0), 3, 2).tangent_lines(Point(3, 0)),
     '[Line2D(Point2D(3, 0), Point2D(3, -12))]'),
    (lambda: Ellipse(Point(0, 0), 3, 2).tangent_lines(Point(3, 5)),
     '[Line2D(Point2D(3, 5), Point2D(3, 6)), Line2D(Point2D(3, 5), Point2D(-63/29, 40/29))]'),
    (lambda: Line3D(Point3D(0, 0, 0), Point3D(1, 0, 0)).perpendicular_line(Point3D(2, 3, 4)),
     'Line3D(Point3D(2, 3, 4), Point3D(2, 0, 0))'),
    (lambda: Line(Point(0, 0), Point(1, 1)).perpendicular_line(Point(2, 0)),
     'Line2D(Point2D(2, 0), Point2D(1, 1))'),
    (lambda: Line3D(Point3D(0, 0, 0), Point3D(1, 0, 0)).distance(
        Line3D(Point3D(0, 1, 1), Point3D(0, 2, 3))),
     'sqrt(5)/5'),
    (lambda: Matrix([[1, 2, 0], [0, 3, 0], [2, -4, 2]]).diagonalize(),
     '(Matrix([\n[-1, 0, -1],\n[ 0, 0, -1],\n[ 2, 1,  2]]), Matrix([\n[1, 0, 0],\n[0, 2, 0],\n[0, 0, 3]]))'),
    (lambda: DiracDelta(x**2 - 1).expand(diracdelta=True, wrt=x),
     'DiracDelta(x - 1)/2 + DiracDelta(x + 1)/2'),
    (lambda: (erf2inv(oo, y), erf2inv(0, y), erf2inv(x, 0), erf2inv(x, oo),
              erf2inv(1, 0), erf2inv(0, 1)),
     '(erfcinv(-y), erfinv(y), x, erfinv(x), 1, oo)'),
    (lambda: factor(5*x + 3*exp(2 - 7*x), deep=True),
     '(5*x*exp(7*x) + 3*exp(2))*exp(-7*x)'),
    (lambda: factor(sin(x**2 - 1) + cos(x*y + x), deep=True),
     'sin((x - 1)*(x + 1)) + cos(x*(y + 1))'),
    (lambda: together(exp(1/x + 1/y), deep=True), 'exp((x + y)/(x*y))'),
    (lambda: together(sin(1/x + 1/y) + 1/x, deep=True), '(x*sin((x + y)/(x*y)) + 1)/x'),
    (lambda: unrad(sqrt(x) + root(x + 1, 3)), '(-x**3 + x**2 + 2*x + 1, [])'),
    (lambda: unrad(sqrt(x) + root(x, 3) + 1)[0], '_p**3 + _p**2 + 1'),
    (lambda: horner(4*x**2*y**2 + 2*x**2*y + 2*x*y**2 + x*y, wrt=y),
     'y*(x*y*(4*x + 2) + x*(2*x + 1))'),
    (lambda: rad_rationalize(sqrt(3), 1 + sqrt(2)/3), '(-sqrt(3) + sqrt(6)/3, -7/9)'),
    (lambda: rad_rationalize(1, sqrt(2) + sqrt(3) + 1), '(-2*sqrt(6) + 2*sqrt(2) + 4, 8)'),
    (lambda: split_surds(2 + sqrt(5) + sqrt(3)*sqrt(5)), '(15, sqrt(3)/3 + 1, 2)'),
]


@pytest.mark.parametrize("build,expected", CASES)
def test_form_matches_sympy(build, expected):
    assert str(build()) == expected


def test_diagonalize_errors():
    from symox import MatrixError, NonSquareMatrixError
    with pytest.raises(MatrixError):
        Matrix([[1, 1], [0, 1]]).diagonalize()
    with pytest.raises(NonSquareMatrixError):
        Matrix([[1, 2, 3]]).diagonalize()
