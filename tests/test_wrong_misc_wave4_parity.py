"""Results that disagreed with SymPy 1.14 (documentation-example audit WP-G14).

Reference values are hard-coded from SymPy 1.14, except where SymPy is
mathematically wrong (Poisson process conditional probabilities): those
cases assert the exact value and are listed in docs/algorithmic-parity.md.
"""
import pytest

import symox as sp
from symox import (Contains, E as E1, Eq, Interval, Matrix, MatrixSymbol, Ne,
                   Piecewise, Rational, SparseMatrix, Symbol, ccode, exp, pi,
                   sqrt, symbols)


# -- physics.optics ---------------------------------------------------------

def test_deviation_sign_and_vector_forms():
    from symox.geometry import Plane, Point3D, Ray3D
    from symox.physics.optics import deviation
    assert str(deviation(0.1, 1.2, 1.5)) == '-0.0200481144433483'
    assert str(deviation(0.3, 1.2, 1.1)) == '0.0282486489090274'
    n = Matrix([0, 0, 1])
    P = Plane(Point3D(0, 0, 0), normal_vector=[0, 0, 1])
    r1 = Ray3D(Point3D(-1, -1, 1), Point3D(0, 0, 0))
    assert deviation(r1, 1, 1, n) == 0
    assert deviation(r1, 1.33, 1, plane=P) is None
    assert deviation(r1, 1.2, 1, plane=P) == -sp.acos(-Rational(1, 5)) + sp.acos(-sqrt(3)/3)


def test_refraction_angle_reference_forms():
    from symox.geometry import Plane, Point3D, Ray3D
    from symox.physics.optics import refraction_angle
    n = Matrix([0, 0, 1])
    P = Plane(Point3D(0, 0, 0), normal_vector=[0, 0, 1])
    r1 = Ray3D(Point3D(-1, -1, 1), Point3D(0, 0, 0))
    assert refraction_angle(r1, 1, 1, n) == Matrix([1, 1, -1])
    assert refraction_angle(r1, 1, 1, plane=P) == Ray3D(Point3D(0, 0, 0), Point3D(1, 1, -1))
    assert str(refraction_angle(0.3, 1.2, 1)) == '0.362512166650728'
    with pytest.raises(ValueError, match='total internal reflection'):
        refraction_angle(1.2, 1.5, 1)


def test_fresnel_coefficients_total_internal_reflection():
    from symox.physics.optics import fresnel_coefficients
    tir = fresnel_coefficients(0.6, 2, 1)
    assert len(tir) == 2
    assert [str(v) for v in tir] == ['-0.235625382192159 - 0.971843958291041*I',
                                     '0.816477005968898 - 0.577377951366403*I']
    exact = fresnel_coefficients(pi/4, 2, 1)
    I = sp.I
    assert exact == [(sqrt(2) - 4*I)/(sqrt(2) + 4*I), (sqrt(2) - I)/(sqrt(2) + I)]
    partial = fresnel_coefficients(0.2, 1.5, 1)
    assert [str(v) for v in partial] == ['-0.187312004995942', '0.212621270579866',
                                         '1.21903199250609', '1.21262127057987']
    assert len(fresnel_coefficients(0.5, 1, 1.2)) == 4
    with pytest.raises(ValueError):
        fresnel_coefficients(2, 1, 1)


def test_lens_formulas_at_infinity():
    from symox.physics.optics import lens_formula, lens_makers_formula, mirror_formula
    assert lens_formula(sp.oo, 10) == 10
    assert mirror_formula(u=sp.oo, v=10) == 10
    assert lens_makers_formula(1, 1, 10, 10) == sp.oo


# -- physics.wigner ---------------------------------------------------------

def test_real_gaunt_reference_sign_convention():
    from symox.physics.wigner import real_gaunt
    assert real_gaunt(1, 1, 2, -1, 1, -2) == sqrt(15)/(10*sqrt(pi))
    assert real_gaunt(1, 1, 2, 1, -1, -2) == sqrt(15)/(10*sqrt(pi))
    assert real_gaunt(2, 1, 1, -2, 1, -1) == sqrt(15)/(10*sqrt(pi))
    assert real_gaunt(1, 1, 2, -1, -1, 0) == -sqrt(5)/(10*sqrt(pi))
    assert real_gaunt(2, 2, 4, -1, -1, 0) == -2/(7*sqrt(pi))
    assert real_gaunt(2, 2, 4, -1, -1, -2) == 0
    assert real_gaunt(1, 1, 2, 0, 0, 0) == sqrt(5)/(5*sqrt(pi))
    with pytest.raises(ValueError):
        real_gaunt(1.5, 1, 1, 0, 0, 0)


# -- stats.PoissonProcess ---------------------------------------------------

def test_poisson_process_docstring_examples():
    from symox.stats import E, P, PoissonProcess
    X = PoissonProcess("X", lamda=3)
    t1, t2 = symbols('t1 t2', positive=True)
    assert X.state_space == sp.S.Naturals0
    assert P(X(t1) < 4) == (9*t1**3/2 + 9*t1**2/2 + 3*t1 + 1)*exp(-3*t1)
    assert P(Eq(X(t1), 2) & Eq(X(t2), 3),
             Contains(t1, Interval.Lopen(0, 2)) & Contains(t2, Interval.Lopen(2, 4))) == 648*exp(-12)
    assert E(X(t1)) == 3*t1
    assert E(X(t1)**2 + 2*X(t2),
             Contains(t1, Interval.Lopen(0, 1)) & Contains(t2, Interval.Lopen(1, 2))) == 18
    assert P(X(3) < 1, Eq(X(1), 0)) == exp(-6)
    assert P(Eq(X(4), 3), Eq(X(2), 3)) == exp(-6)
    assert P(Eq(X(5), 4), Eq(X(2), 1)) == 243*exp(-9)/2
    assert P(X(t1) < 2, Contains(t1, Interval(1, 3))) == 7*exp(-6)
    assert E(X(t1), Contains(t1, Interval(1, 3))) == 6
    Y = PoissonProcess("Y", lamda=4)
    Z = X + Y
    assert Z.lamda == 7
    N, M = Z.split(l1=2, l2=5)
    assert (N.lamda, M.lamda) == (2, 5)
    with pytest.raises(ValueError):
        X + 1
    with pytest.raises(ValueError):
        Z.split(1, 2)


def test_poisson_process_contains_errors_and_overlap():
    from symox.stats import P, PoissonProcess
    from symox.stats import Probability
    X = PoissonProcess("X", lamda=3)
    t1, t2 = symbols('t1 t2', positive=True)
    with pytest.raises(ValueError, match='finite bounds'):
        P(X(t1) < 2, Contains(t1, Interval(1, sp.oo)))
    with pytest.raises(ValueError, match='Contains'):
        P(X(t1) < 2, Contains(t1, Interval(1, 3)) & (X(2) > 1))
    # Overlapping intervals stay unevaluated, as in the reference.
    result = P(Eq(X(t1), 2) & Eq(X(t2), 3), Contains(t1, Interval(0, 2)) & Contains(t2, Interval(1, 4)))
    assert isinstance(result, Probability)


def test_poisson_process_exact_where_reference_is_wrong():
    # SymPy 1.14 adds the probabilities of non-disjoint events and does not
    # normalise by the probability of the condition; these are exact values.
    from symox.stats import P, PoissonProcess
    X = PoissonProcess("X", lamda=3)
    t1, t2 = symbols('t1 t2', positive=True)
    # SymPy: 1 - 36*exp(-6)
    assert P(Eq(X(t1), 2) | Ne(X(t1), 4), Contains(t1, Interval.Ropen(2, 4))) == 1 - 54*exp(-6)
    # SymPy: 5*exp(-3)
    assert sp.simplify(P(X(2) <= 3, X(1) > 1) - 45*exp(-3)/(2*(exp(3) - 4))) == 0
    # SymPy: 27*exp(-3)/2
    assert sp.simplify(P(X(2) <= 3, X(1) >= 1) - 48*exp(-3)/(exp(3) - 1)) == 0
    # SymPy: exp(-6)
    assert sp.simplify(P(X(3) > 2, X(1) < 2) - (1 - 23*exp(-6)/2)) == 0
    # SymPy: exp(-9)
    assert P(X(4) < 3, Eq(X(1), 1)) == 10*exp(-9)
    # SymPy: 0
    assert sp.simplify(P(X(2) < 3, X(2) > 1) - 18/(exp(6) - 7)) == 0
    # SymPy: 4*exp(-6)
    assert sp.simplify(P((X(1) < 2) & (X(2) > 3)) - (4*exp(-3) - 77*exp(-6)/2)) == 0
    # SymPy: exp(-6) + 3*exp(-3)
    value = P(Eq(X(t1), 1) | Eq(X(t2), 0),
              Contains(t1, Interval.Lopen(0, 1)) & Contains(t2, Interval.Lopen(1, 3)))
    assert sp.simplify(value - (3*exp(-3) + exp(-6) - 3*exp(-9))) == 0


# -- printing.maple ---------------------------------------------------------

def test_maple_code_reference_printer():
    from symox.printing.maple import maple_code
    x, y, z = symbols('x y z')
    A = MatrixSymbol('A', 2, 2)
    assert maple_code(x**2 + 1, assign_to=y) == 'y := x^2 + 1'
    assert maple_code(sp.sin(x), assign_to='z') == 'z := sin(x)'
    assert maple_code(Piecewise((x, x > 1), (y, True)), assign_to=z) == 'z := piecewise(x > 1, x, y)'
    assert maple_code(Piecewise((x, x > 1), (y, True))) == 'piecewise(x > 1, x, y)'
    assert maple_code(sp.log(x)) == 'ln(x)'
    assert maple_code(1/x) == '1/x'
    assert maple_code(1/sqrt(x)) == '1/sqrt(x)'
    assert maple_code(Eq(x, y)) == 'x = y'
    assert maple_code(Ne(x, y)) == 'x <> y'
    assert maple_code(sp.And(x > 1, y < 2)) == 'x > 1 and y < 2'
    assert maple_code(Matrix([[1, x], [y, 2]])) == 'Matrix([[1, x], [y, 2]], storage = rectangular)'
    assert maple_code(SparseMatrix([[1, 0], [0, 2]])) == 'Matrix([[1, 0], [0, 2]], storage = sparse)'
    assert maple_code(A*A) == 'MatrixPower(A, 2)'
    assert maple_code(A[0, 1]) == 'A[1, 2]'
    assert maple_code(sp.Identity(3)) == 'Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]], storage = sparse)'
    f = sp.Function('f')
    assert maple_code(sp.Derivative(f(x), x, 2)) == 'diff(f(x), x$2)'
    assert maple_code(sp.besselj(1, x)) == 'BesselJ(1, x)'
    assert maple_code(sp.GoldenRatio) == '(1/2 + (1/2)*sqrt(5))'
    assert maple_code(sp.EulerGamma) == 'gamma'
    assert maple_code(x*Rational(2, 3)) == '(2/3)*x'
    assert maple_code(2.5*x) == '2.5*x'
    assert maple_code([x, y], assign_to=['a', 'b']) == 'a := x\nb := y'
    assert maple_code(x + sp.zoo, strict=False) == '# Not supported in maple:\n# ComplexInfinity\nx + zoo'
    with pytest.raises(ValueError):
        maple_code(Piecewise((x, x > 1)))


def test_sympify_keeps_sparse_matrices_sparse():
    m = SparseMatrix([[1, 0], [0, 2]])
    s = sp.sympify(m)
    assert type(s).__name__ == 'ImmutableSparseMatrix'
    assert s == m
    assert sp.S.BooleanTrue is sp.true and sp.S.BooleanFalse is sp.false


# -- codegen.ast ------------------------------------------------------------

def test_element_indices_strides_and_offset():
    from symox.codegen.ast import Element
    elem = Element('x', 'ijk')
    i, j, k = symbols('i j k')
    assert elem.symbol.name == 'x'
    assert elem.indices == sp.Tuple(i, j, k)
    assert isinstance(elem.indices, sp.Tuple)
    assert ccode(elem) == 'x[i][j][k]'
    assert ccode(Element('x', 'ijk', strides='lmn', offset='o')) == 'x[i*l + j*m + k*n + o]'
    assert ccode(Element('x', 'ij', strides=(3, 1))) == 'x[3*i + j]'
    assert ccode(Element('x', 'i', strides='l', offset=2)) == 'x[i*l + 2]'
    with pytest.raises(ValueError, match='Expected strides'):
        ccode(Element('x', 'ij', offset='o'))
    assert repr(Element('x', 'ij')) == "Element(Symbol('x'), indices=(Symbol('i'), Symbol('j')))"
    assert str(Element('x', 'ij', strides='lm', offset='o')) == 'Element(x, indices=(i, j), strides=(l, m), offset=o)'
    assert repr(Element('x', [i])) == "Element(Symbol('x'), indices=(Symbol('i'),))"
    assert Element('x', 'ij') == Element(Symbol('x'), (i, j))
    assert Element('x', 'ij').free_symbols == {Symbol('x'), i, j}
    assert Element('x', 'ij').subs(i, 2) == Element('x', (2, j))
    assert sp.sympify(elem) is elem
    with pytest.raises(TypeError):
        elem + 1


def test_node_attr_params_missing_is_none():
    from symox.codegen.ast import Node, Variable, value_const
    from symox.codegen.fnodes import dimension
    x = Symbol('x')
    assert Node().attr_params('x') is None
    assert Variable(x, attrs=[value_const]).attr_params('value_const') == ()
    assert tuple(Variable(x, attrs=[dimension(2, 3)]).attr_params("dimension")) == (2, 3)
    assert Variable(x, attrs=[value_const]).attr_params('dimension') is None
