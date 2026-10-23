"""SymPy 1.14 parity: symbolic spin coupling, category sets, matrix-bodied
Lambda and diffgeom coordinate symbols (expected values from SymPy 1.14)."""
import pytest

from symox import (S, symbols, Lambda, Matrix, sqrt, atan2, cos, sin,
                   FiniteSet)
from symox.physics.quantum.spin import (JzKet, JzBra, JxKetCoupled,
                                        JzKetCoupled, couple, uncouple)
from symox.physics.quantum.tensorproduct import TensorProduct

j, m, j1, j2, j3, m1, m2, m3 = symbols('j m j1 j2 j3 m1 m2 m3')


def test_couple_symbolic():
    s = couple(TensorProduct(JzKet(j1, m1), JzKet(j2, m2)))
    assert str(s) == 'Sum(CG(j1, m1, j2, m2, j, m1 + m2)*|j,m1 + m2,j1=j1,j2=j2>, (j, m1 + m2, j1 + j2))'
    assert str(couple(TensorProduct(JzKet(j1, m1), JzKet(j2, m2), JzKet(j3, m3)))) == (
        'Sum(CG(j1, m1, j2, m2, j12, m1 + m2)*CG(j12, m1 + m2, j3, m3, j, m1 + m2 + m3)'
        '*|j,m1 + m2 + m3,j1=j1,j2=j2,j3=j3,j(1,2)=j12>, (j12, m1 + m2, j1 + j2), (j, m1 + m2 + m3, j12 + j3))')
    assert str(couple(TensorProduct(JzKet(j1, m1), JzKet(j2, m2), JzKet(j3, m3)), ((1, 3), (1, 2)))) == (
        'Sum(CG(j1, m1, j3, m3, j13, m1 + m3)*CG(j13, m1 + m3, j2, m2, j, m1 + m2 + m3)'
        '*|j,m1 + m2 + m3,j1=j1,j2=j2,j3=j3,j(1,3)=j13>, (j13, m1 + m3, j1 + j3), (j, m1 + m2 + m3, j13 + j2))')
    assert str(couple(TensorProduct(JzKet(S(1)/2, m1), JzKet(S(1)/2, S(1)/2)))) == (
        'Sum(CG(1/2, m1, 1/2, 1/2, j, m1 + 1/2)*|j,m1 + 1/2,j1=1/2,j2=1/2>, (j, m1 + 1/2, 1))')
    assert str(couple(2*TensorProduct(JzKet(j1, m1), JzKet(j2, m2)))) == (
        '2*Sum(CG(j1, m1, j2, m2, j, m1 + m2)*|j,m1 + m2,j1=j1,j2=j2>, (j, m1 + m2, j1 + j2))')
    assert str(couple(TensorProduct(JzBra(j1, m1), JzBra(j2, m2)))) == (
        'Sum(CG(j1, m1, j2, m2, j, m1 + m2)*<j,m1 + m2,j1=j1,j2=j2|, (j, m1 + m2, j1 + j2))')
    assert str(s.subs(j1, 1)) == 'Sum(CG(1, m1, j2, m2, j, m1 + m2)*|j,m1 + m2,j1=1,j2=j2>, (j, m1 + m2, j2 + 1))'
    with pytest.raises(ValueError, match='Coefficients must be numerical'):
        s.doit()


def test_uncouple_symbolic():
    ref = 'Sum(CG(j1, m1, j2, m2, j, m)*|j1,m1>x|j2,m2>, (m1, -j1, j1), (m2, -j2, j2))'
    u = uncouple(JzKetCoupled(j, m, (j1, j2)))
    assert str(u) == ref
    assert str(uncouple(JzKet(j, m), (j1, j2))) == ref
    assert u == uncouple(JzKet(j, m), (j1, j2))
    assert str(uncouple(JxKetCoupled(j, m, (j1, j2)))) == ref
    assert str(uncouple(JzKetCoupled(j, m, (j1, j2, j3)))) == (
        'Sum(CG(j1, m1, j2, m2, j1 + j2, m1 + m2)*CG(j1 + j2, m1 + m2, j3, m3, j, m)'
        '*|j1,m1>x|j2,m2>x|j3,m3>, (m1, -j1, j1), (m2, -j2, j2), (m3, -j3, j3))')
    assert str(uncouple(JzKetCoupled(j, m, (j1, j2, j3), ((1, 3, j3), (1, 2, j))))) == (
        'Sum(CG(j1, m1, j3, m3, j3, m1 + m3)*CG(j3, m1 + m3, j2, m2, j, m)'
        '*|j1,m1>x|j2,m2>x|j3,m3>, (m1, -j1, j1), (m2, -j2, j2), (m3, -j3, j3))')
    assert str(uncouple(JzKetCoupled(1, m, (S(1)/2, S(1)/2)))) == (
        'Sum(CG(1/2, m1, 1/2, m2, 1, m)*|1/2,m1>x|1/2,m2>, (m1, -1/2, 1/2), (m2, -1/2, 1/2))')
    assert str(uncouple(JzKetCoupled(j, 1, (S(1)/2, S(1)/2)))) == (
        'Sum(CG(1/2, m1, 1/2, m2, j, 1)*|1/2,m1>x|1/2,m2>, (m1, -1/2, 1/2), (m2, -1/2, 1/2))')
    # Bound variables are not substituted; free ones are, inside the kets too.
    assert u.subs(m1, S(1)/2) == u
    assert str(u.subs(j1, 1)) == 'Sum(CG(1, m1, j2, m2, j, m)*|1,m1>x|j2,m2>, (m1, -1, 1), (m2, -j2, j2))'
    assert str(-u) == '-' + ref and str(u/2) == ref + '/2' and str(u*j) == 'j*' + ref
    with pytest.raises(ValueError, match='Coefficients must be numerical'):
        uncouple(JzKetCoupled(1, m, (S(1)/2, S(1)/2))).subs(m, 0).doit()


def test_numeric_coupling_unchanged():
    half = S(1)/2
    assert couple(TensorProduct(JzKet(half, half), JzKet(half, half))) == JzKetCoupled(1, 1, (half, half))
    assert uncouple(JzKetCoupled(1, 1, (half, half))) == TensorProduct(JzKet(half, half), JzKet(half, half))


def test_tensor_product_subs_inside_factors():
    tp = TensorProduct(JzKet(j1, m1), JzKet(j2, m2))
    assert str(tp.subs(m1, S(1)/2)) == '|j1,1/2>x|j2,m2>'


def test_category_finite_sets():
    from symox.categories import Object, NamedMorphism, Diagram, Category
    A, B, C = Object("A"), Object("B"), Object("C")
    f, g = NamedMorphism(A, B, "f"), NamedMorphism(B, C, "g")
    d = Diagram([f, g], {f: "unique", g: "unique"})
    K = Category("K", commutative_diagrams=[d])
    assert K.commutative_diagrams == FiniteSet(d)
    assert str(FiniteSet(A, B)) == '{Object("A"), Object("B")}'
    assert FiniteSet(A, B) == FiniteSet(B, A)
    d1 = d.subdiagram_from_objects(FiniteSet(A, B))
    assert d1 == Diagram([f], {f: "unique"})


def test_lambda_matrix_body():
    from symox._core import BadArgumentsError
    x, y = symbols('x y', real=True)
    L = Lambda((x, y), Matrix([sqrt(x**2 + y**2), atan2(y, x)]))
    assert str(L) == 'Lambda((x, y), Matrix([\n[sqrt(x**2 + y**2)],\n[      atan2(y, x)]]))'
    assert L(1, 2) == Matrix([sqrt(5), atan2(2, 1)])
    assert L == Lambda((x, y), Matrix([sqrt(x**2 + y**2), atan2(y, x)]))
    assert isinstance(L, Lambda) and L.variables == (x, y) and L.free_symbols == set()
    assert L.nargs == FiniteSet(2)
    assert Lambda(x, Matrix([[x, 1], [2, x]]))(3) == Matrix([[3, 1], [2, 3]])
    with pytest.raises(BadArgumentsError, match=r'takes exactly 2 arguments \(1 given\)'):
        L(1)


def test_coordinate_symbol():
    from symox.diffgeom import Manifold, Patch, CoordSystem
    m = Manifold('M', 2)
    p = Patch('P', m)
    x, y = symbols('x y', real=True)
    r, theta = symbols('r theta', nonnegative=True)
    relation_dict = {
        ('Car2D', 'Pol'): Lambda((x, y), Matrix([sqrt(x**2 + y**2), atan2(y, x)])),
        ('Pol', 'Car2D'): Lambda((r, theta), Matrix([r*cos(theta), r*sin(theta)])),
    }
    Car2D = CoordSystem('Car2D', p, [x, y], relation_dict)
    Pol = CoordSystem('Pol', p, [r, theta], relation_dict)
    x, y = Car2D.symbols
    R, T = Pol.symbols
    assert x.name == 'x' and x.coord_sys == Car2D and x.index == 0 and x.is_real
    assert str(x.rewrite(Pol)) == 'r*cos(theta)'
    assert str(sqrt(x**2 + y**2).rewrite(Pol).simplify()) == 'r'
    assert str(x.rewrite(Car2D)) == 'x'
    assert str((x*y + 1).rewrite(Pol)) == 'r**2*sin(theta)*cos(theta) + 1'
    assert str(R.rewrite(Car2D)) == 'sqrt(x**2 + y**2)'
    assert str((x + R).rewrite(Pol)) == 'r*cos(theta) + r'
