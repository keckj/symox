"""Wave 4 parity: core matching, series, infinite sums, expand_power_base,
fancy sets and matrix helpers. Expected values come from SymPy 1.14."""
import pytest

import symox as sm
from symox import (Abs, Basic, ComplexRegion, Function, ImageSet,
                   Interval, KroneckerDelta, Lambda, Matrix, MatrixSymbol, N,
                   O, Range, Rational, S, SparseMatrix, Sum, Symbol, Trace,
                   Wild, WildFunction, cos, diag, exp, expand_power_base, eye,
                   factorial, hessian, log, oo, pi, sin, sqrt, symbols)
from symox.matrices.expressions import FunctionMatrix
from symox.solvers.polysys import solve_triangulated

x, y, z, a, b, c = symbols('x y z a b c')
f, g = Function('f'), Function('g')


def test_basic_matches():
    w = Wild('w')
    assert Basic(a + w, w).matches(Basic(a + b, c)) is None
    assert Basic(a + w, w).matches(Basic(a + b + c, b + c)) == {w: b + c}


def test_matches_with_old_dicts():
    w, v = Wild('w'), Wild('v')
    assert w.matches(x, {w: y}) == {w: x}
    assert w.matches(x, {v: y}) == {v: y, w: x}
    assert (w**v).matches(x**2, {w: y}) == {w: x, v: 2}
    assert (w + v*y).matches(x + 2*y, {w: x}) == {w: x, v: 2}
    assert (w + v*y).matches(x + 2*y, {w: z}) == {w: x, v: 2}
    assert (x + w).matches(x + y, None, old=True) == {w: y}
    assert (w*v).matches(2*x, {w: 2}) == {w: x, v: 2}


def test_matmul_matches_with_old_dicts():
    A, B = MatrixSymbol('A', 2, 2), MatrixSymbol('B', 2, 2)
    w = Wild('w')
    assert (w*A).matches(3*A) == {w: 3}
    assert (w*A).matches(3*A, {w: 3}) == {w: 3}
    assert (w*A).matches(3*A, {w: 4}) == {w: 3}
    assert (w*A*B).matches(3*A*B) == {w: 3}


def test_wildfunction_nargs_per_instance():
    F = WildFunction('F')
    assert F.matches(f(x)) == {F: f(x)}
    F2 = WildFunction('F', nargs=2)
    assert F2.matches(f(x)) is None
    assert F2.matches(f(x, y)) == {F2: f(x, y)}
    assert F.matches(f(x)) == {F: f(x)}


def test_series_generator_and_logx():
    term = cos(x).series(n=None)
    assert [next(term) for _ in range(2)] == [1, -x**2/2]
    term = Abs(x + 1).series(x, n=None)
    assert [next(term) for _ in range(2)] == [1, x]
    logx = Symbol('logx')
    assert sin(log(x)).nseries(x, 0, 6, logx=logx) == sin(logx)
    assert (x*log(x)).nseries(x, 0, 3, logx=logx) == logx*x


def test_abs_series():
    assert Abs(x).series(x, 0, 4) == x
    assert Abs(sin(x)).series(x, 0, 6) == x - x**3/6 + x**5/120 + O(x**6)
    assert Abs(x - 1).series(x, 0, 3) == 1 - x
    assert Abs(x**2 - x).series(x, 0, 4) == x - x**2
    assert Abs(x).series(x, 0, 3, dir='-') == -x
    assert Abs(sin(x)).nseries(x, 0, 5) == x - x**3/6 + O(x**5)
    assert Abs(x + x**2).nseries(x, 0, 4) == x + x**2
    assert Abs(y*x).series(x, 0, 3) == x*Abs(y)


def test_extract_leading_order():
    assert (x + 1 + 1/x**5).extract_leading_order(x) == ((x**-5, O(x**-5)),)
    assert (1 + x).extract_leading_order(x) == ((1, O(1, x)),)
    assert (x + x**2).extract_leading_order(x) == ((x, O(x)),)
    assert set((x + y*x + x**2).extract_leading_order(x)) == {(x*y, O(x)), (x, O(x))}


def test_infinite_sum_evalf():
    k = Symbol('k')
    assert str(N(Sum(1/k**k, (k, 1, oo)), 4)) == '1.291'
    assert str(Sum(1/k**2, (k, 1, oo)).evalf()) == '1.64493406684823'
    assert str(Sum(1/k**k, (k, 1, oo)).evalf(30)) == '1.29128599706266354040728259060'
    assert str(Sum((-1)**k/k, (k, 1, oo)).evalf()) == '-0.693147180559945'
    assert str(Sum(1/factorial(k), (k, 0, oo)).evalf(20)) == '2.7182818284590452354'
    assert str(Sum(k**-3, (k, 1, oo)).n(15)) == '1.20205690315959'
    assert str((1 + Sum(2**-k, (k, 0, oo))).evalf(25)) == '3.000000000000000000000000'
    assert str(Sum(Rational(1, 2)**k/k, (k, 1, oo)).evalf(20)) == '0.69314718055994530942'
    with pytest.raises(ValueError, match=r'diverges like n\^-1'):
        Sum(1/k, (k, 1, oo)).evalf()
    with pytest.raises(ValueError, match=r'diverges like \(2\)\^n'):
        Sum(2**k, (k, 1, oo)).evalf()


def test_euler_maclaurin():
    k = Symbol('k')
    assert Sum(1/k, (k, 2, 5)).euler_maclaurin() == (-log(2) + Rational(7, 20) + log(5), Rational(7, 400))
    s, e = Sum(1/k, (k, 2, 5)).euler_maclaurin(2, 2)
    assert s == -log(4) + Rational(20355631, 19200000) + log(5)
    assert e == Rational(183, 256000000)
    assert Sum(k, (k, 1, oo)).euler_maclaurin(1, 1) == (oo, 0)


def test_expand_power_base():
    w = Symbol('w')
    p = Symbol('p', positive=True)
    assert expand_power_base((x*y)**2) == x**2*y**2
    assert expand_power_base((2*x)**y) == 2**y*x**y
    assert expand_power_base(((x + y)*z)**2) == z**2*(x + y)**2
    assert expand_power_base((2*y)**(1 + z)) == 2**(z + 1)*y**(z + 1)
    assert expand_power_base(((x*y)**z)**2) == (x*y)**(2*z)
    assert expand_power_base((x*y)**z, force=True) == x**z*y**z
    assert expand_power_base(sqrt(-2*x)) == sqrt(2)*sqrt(-x)
    assert expand_power_base(sqrt(-x*y)) == sqrt(-x*y)
    assert expand_power_base((x*y/z)**2) == x**2*y**2/z**2
    assert expand_power_base((p*x)**y) == p**y*x**y
    assert expand_power_base((x*y)**y, deep=False) == (x*y)**y
    assert expand_power_base((2*exp(y))**x) == 2**x*exp(y)**x
    assert expand_power_base((x*(y*z)**w)**2) == x**2*(y*z)**(2*w)
    assert expand_power_base(sin((x*y)**2)) == sin(x**2*y**2)
    assert expand_power_base((3*x*sqrt(y))**z) == 3**z*(x*sqrt(y))**z
    assert expand_power_base((-2*x)**y) == 2**y*(-x)**y


def test_solve_triangulated_ground_domain():
    F = [x**2 + y + z - 1, x + y**2 + z - 1, x + y + z**2 - 1]
    assert solve_triangulated(F, x, y, z) == [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
    assert solve_triangulated([x**2 - 2, y - x], x, y) == []
    assert solve_triangulated([x + y - 3, x - y - 1], x, y) == [(2, 1)]
    assert solve_triangulated([x**2 - 2, y - x], x, y, domain=sm.QQ.algebraic_field(sqrt(2))) == [
        (-sqrt(2), -sqrt(2)), (sqrt(2), sqrt(2))]
    # The algebraic extension tower is a documented gap.
    with pytest.raises(NotImplementedError):
        solve_triangulated([x**2 - 2, y - x], x, y, extension=True)


def test_imageset_doc_cases():
    squares = ImageSet(Lambda(x, x**2), S.Naturals)
    assert 4 in squares
    assert squares.lamda(2) == 4
    assert squares.contains(Rational(1, 2)) == sm.false
    assert ImageSet(Lambda(x, 2*x + 1), S.Integers).contains(5) == sm.true
    assert (9 in ImageSet(Lambda(x, x**2), S.Integers)) is True
    assert (-4 in ImageSet(Lambda(x, x**2), S.Integers)) is False
    assert ImageSet(Lambda(x, x**3 + x), S.Integers).contains(10) == sm.true
    assert Interval(-1, 1).intersect(ImageSet(Lambda(x, x*pi), S.Integers)) == sm.FiniteSet(0)
    it = iter(squares)
    assert [next(it) for _ in range(4)] == [1, 4, 9, 16]
    assert (3 in squares) is False


def test_range_indexing():
    r = Range(10, 0, -2)
    assert (r[-1], r[0], r[1:3], r.inf, r.sup, len(r)) == (2, 10, Range(8, 4, -2), 2, 10, 5)
    r = Range(0, 10, 3)
    assert (r[-1], r[::2], r[1::-1], list(r[::-1])) == (9, Range(0, 12, 6), Range(3, -3, -3), [9, 6, 3, 0])
    r = Range(1, oo)
    # Intentional divergence: SymPy 1.14 raises IndexError for r[5] (it
    # tests oo.is_nonnegative) and returns Range(3, 5) for r[2:5] and
    # Range(0) for r[0:3] (an off-by-one through the reversed range).
    assert (r[0], r[5], r[2:5], r[0:3], r[3:]) == (1, 6, Range(3, 6), Range(1, 4), Range(4, oo))
    assert Range(0, oo, 2)[1:4] == Range(2, 8, 2)
    assert r[5:1:-2] == Range(6, 2, -2)
    assert (Range(-oo, 0)[-1], Range(-oo, 0)[-3]) == (-1, -3)
    with pytest.raises(IndexError):
        Range(5)[10]
    with pytest.raises(ValueError):
        Range(1, oo)[-1]
    with pytest.raises(ValueError):
        Range(1, oo)[::2]
    r = Range(0, 10, 2)
    assert (r.reversed, r.size, list(r), 4 in r, 5 in r) == (Range(8, -2, -2), 5, [0, 2, 4, 6, 8], True, False)
    n = Symbol('n', integer=True)
    r = Range(n, n + 20, 3)
    assert (r.size, r.inf, r[0], r[-1]) == (7, n, n, n + 18)
    assert Range(0, 10, 2)[3:] == Range(6, 10, 2)


def test_complex_region():
    a, b, c8 = Interval(2, 3), Interval(4, 6), Interval(1, 8)
    c1 = ComplexRegion(a*b)
    assert str(c1) == 'CartesianComplexRegion(ProductSet(Interval(2, 3), Interval(4, 6)))'
    assert (c1.a_interval, c1.b_interval, c1.polar) == (a, b, False)
    assert str(ComplexRegion.from_real(Interval(0, 1))) == 'CartesianComplexRegion(ProductSet(Interval(0, 1), {0}))'
    assert str(ComplexRegion.union(c1, ComplexRegion(b*c8))) == (
        'CartesianComplexRegion(Union(ProductSet(Interval(2, 3), Interval(4, 6)), '
        'ProductSet(Interval(4, 6), Interval(1, 8))))')
    assert str(c1.intersect(ComplexRegion(b*c8))) == 'CartesianComplexRegion(EmptySet)'
    unit = ComplexRegion(Interval(0, 1)*Interval(0, 1))
    assert unit.union(Interval(0, 1)) == unit
    disk = ComplexRegion(Interval(0, 1)*Interval(0, 2*pi), polar=True)
    upper = ComplexRegion(Interval(0, 1)*Interval(0, pi), polar=True)
    assert str(upper) == 'PolarComplexRegion(ProductSet(Interval(0, 1), Interval(0, pi)))'
    assert disk.psets == (Interval(0, 1)*Interval.Ropen(0, 2*pi),)
    assert disk.intersect(upper) == upper
    quarter = ComplexRegion(Interval(0, 1)*Interval(pi/2, pi), polar=True)
    assert disk.intersect(quarter) == quarter


def test_hessian_bordered():
    F, G = f(x, y), g(x, y)
    g2 = x**2 + 3*y
    assert hessian(F, (x, y), [G, g2]) == Matrix([
        [0, 0, G.diff(x), G.diff(y)],
        [0, 0, 2*x, 3],
        [G.diff(x), 2*x, F.diff(x, 2), F.diff(x, y)],
        [G.diff(y), 3, F.diff(x, y), F.diff(y, 2)]])
    expected = Matrix([[0, 1, 1], [1, 2*y, 2*x], [1, 2*x, 0]])
    assert hessian(x**2*y + z, (x, y), [x + y + z]) == expected
    assert hessian(x**2*y, Matrix([x, y]), Matrix([x + y])) == expected
    assert hessian(F, (x, y)) == Matrix([[F.diff(x, 2), F.diff(x, y)], [F.diff(x, y), F.diff(y, 2)]])


def test_diag_unpack():
    assert diag([1, 2, 3], unpack=True) == Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    assert diag([1, 2, 3]) == Matrix([[1], [2], [3]])
    assert diag([[1, 2], [3, 4]], unpack=True) == Matrix([[1, 0], [2, 0], [0, 3], [0, 4]])
    assert diag([1, 2], unpack=False) == Matrix([[1], [2]])


def test_trace_simplify():
    A = MatrixSymbol('A', 3, 3)
    assert Trace(eye(3)).simplify() == 3
    assert Trace(A).simplify() == Trace(A)
    assert Trace(A + A).simplify() == 2*Trace(A)
    assert sm.simplify(Trace(Matrix([[1, 2], [3, 4]]))) == 5


def test_function_matrix():
    X = FunctionMatrix(3, 3, KroneckerDelta)
    assert X.as_explicit() == eye(3)
    assert (X[0, 0], X[0, 1], X.shape) == (1, 0, (3, 3))
    assert str(X.lamda) == 'Lambda((_i, _j), KroneckerDelta(_i, _j))'
    assert FunctionMatrix(2, 2, 'lambda i, j: i + j').as_explicit() == Matrix([[0, 1], [1, 2]])
    i, j = symbols('i j')
    assert Trace(FunctionMatrix(3, 3, Lambda((i, j), i + j))).doit() == 6


def test_sparse_row_and_column_lists():
    a = SparseMatrix(((1, 2), (3, 4)))
    assert a.RL == [(0, 0, 1), (0, 1, 2), (1, 0, 3), (1, 1, 4)]
    assert a.CL == [(0, 0, 1), (1, 0, 3), (0, 1, 2), (1, 1, 4)]
    a = SparseMatrix(((0, 2), (3, 0)))
    assert (a.RL, a.CL) == ([(0, 1, 2), (1, 0, 3)], [(1, 0, 3), (0, 1, 2)])
