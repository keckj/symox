"""Coordinate calculus and curvature identities, with independent references."""
import subprocess
import sys
from itertools import product
import pytest
import symox as s
from symox.diffgeom import *
from symox.diffgeom.rn import R2_r,R2_p,R3_r,R3_s


def test_coordinate_atoms_point_evaluation_and_derivations():
    x,y = R2_r.base_scalars()
    ex,ey = R2_r.base_vectors()
    point = R2_r.point((2,3))
    assert x(point) == 2 and y(point) == 3
    assert (x*y+s.sin(x)).rcall(point) == 6+s.sin(2)
    assert (x*y)(point) == 6
    assert ex(x) == 1 and ex(y) == 0 and ey(x*y) == x
    assert ex(s.Symbol('x',real=True)) == 0
    assert ex(R2_r.symbols[0]) == 0
    assert (x*y).atoms(BaseScalarField) == {x,y}
    assert x.coord_sys == R2_r and x.patch == R2_r.patch and x.dim == 2
    assert R2_r.symbols[0].coord_sys == R2_r
    assert s.sympify(x) != s.sympify(R2_r.symbols[0])
    assert x(point) != R2_p.base_scalar(0)(point)
    assert R2_r.point(s.Tuple(2,3)) == point


def test_chart_jacobians_and_basis_changes():
    r,t = R2_p.base_scalars()
    x,y = R2_r.base_scalars()
    ex,ey = R2_r.base_vectors()
    er,et = R2_p.base_vectors()
    assert R2_p.jacobian_determinant(R2_r) == R2_p.symbols[0]
    assert R2_p.transform(R2_r,(2,0)) == s.Matrix([2,0])
    assert R2_r.point((3,4)).coords(R2_p)[0] == 5
    assert R2_p.transformation(R2_r)(2,s.pi/2) == s.Matrix([0,2])
    assert vectors_in_basis(er,R2_r) == s.cos(t)*ex+s.sin(t)*ey
    assert vectors_in_basis(et,R2_r) == -r*s.sin(t)*ex+r*s.cos(t)*ey
    assert er(x) == s.cos(t)
    assert et(x) == -r*s.sin(t)
    # These coordinate bases commute despite their rotating Cartesian components.
    assert Commutator(er,et) == 0


def test_indirect_and_unique_inverse_maps():
    M = Manifold('affine-test',2)
    P = Patch('chart',M)
    x,y,u,v,a,b = s.symbols('x y u v a b')
    relations = {('A','B'):((x,y),(x+y,x-y)),('B','C'):((u,v),(2*u,3*v))}
    A = CoordSystem('A',P,(x,y),relations)
    B = CoordSystem('B',P,(u,v),relations)
    C = CoordSystem('C',P,(a,b),relations)
    assert A.transform(C,(2,1)) == s.Matrix([6,3])
    assert C.transform(A,(6,3)) == s.Matrix([2,1])
    assert A.jacobian_determinant(C) == -12
    assert A.jacobian(C,(2,1)) == s.Matrix([[2,2],[3,-3]])
    other = CoordSystem('disconnected',P,('q','p'))
    with pytest.raises(NotImplementedError):
        A.transform(other)


def test_tensor_and_wedge_products_and_partial_contraction():
    x,y = R2_r.base_scalars()
    ex,ey = R2_r.base_vectors()
    dx,dy = R2_r.base_oneforms()
    assert dx(ex) == 1 and dy(ex) == 0
    assert TensorProduct(dx,dy)(ex,ey) == 1
    assert TensorProduct(dx,dy)(ex) == dy
    assert TensorProduct(dx,dy)(None,ey) == dx
    assert TensorProduct(ex,dy)(x*x,ey) == 2*x
    assert WedgeProduct(dx,dy)(ey,ex) == -1
    assert WedgeProduct(dx,x*dy)(x*ex,ey) == x*x
    assert WedgeProduct(ex,ey)(y,None) == -ex
    assert WedgeProduct(dx,dx) == 0
    assert covariant_order(TensorProduct(dx,dy)) == 2
    assert contravariant_order(TensorProduct(ex,dy)) == 1
    assert covariant_order(TensorProduct(ex,dy)) == 1
    assert twoform_to_matrix(TensorProduct(dx,dy)) == s.Matrix([[0,0],[1,0]])
    forms = R3_r.base_oneforms()
    vectors = R3_r.base_vectors()
    volume = WedgeProduct(WedgeProduct(forms[0],forms[1]),forms[2])
    assert volume(*vectors) == 1
    assert volume(vectors[1],vectors[0],vectors[2]) == -1


def test_exterior_derivative_leibniz_and_nilpotence():
    x,y,z = R3_r.base_scalars()
    dx,dy,dz = R3_r.base_oneforms()
    f = x*x*y+s.sin(z)
    assert Differential(Differential(f)) == 0
    alpha = x*dy+y*dz
    beta = y*dx+z*dy
    assert Differential(Differential(alpha)) == 0
    lhs = Differential(WedgeProduct(alpha,beta))
    rhs = WedgeProduct(Differential(alpha),beta)-WedgeProduct(alpha,Differential(beta))
    assert lhs == rhs
    assert Differential(x*alpha) == WedgeProduct(Differential(x),alpha)+x*Differential(alpha)
    assert Differential(s.Integer(3)) == 0


def test_lie_bracket_jacobi_and_cartan_formula():
    x,y = R2_r.base_scalars()
    ex,ey = R2_r.base_vectors()
    dx,dy = R2_r.base_oneforms()
    X,Y,Z = x*ex+y*ey,y*ex,x*ey
    assert Commutator(ex,x*ey) == ey
    assert Commutator(X,Y)+Commutator(Y,X) == 0
    assert Commutator(X,Commutator(Y,Z))+Commutator(Y,Commutator(Z,X))+Commutator(Z,Commutator(X,Y)) == 0
    alpha = x*y*dx+x*x*dy
    assert LieDerivative(X,alpha) == Differential(alpha)(X,None)+Differential(alpha(X))
    assert LieDerivative(X,Differential(x*y)) == Differential(LieDerivative(X,x*y))
    identity = TensorProduct(ex,dx)+TensorProduct(ey,dy)
    assert LieDerivative(X,identity) == 0


@pytest.mark.parametrize('kind',['flat','polar','sphere','hyperbolic','warped'])
def test_metric_connection_and_curvature_against_sympy(kind):
    import sympy as sp
    from sympy import diffgeom as reference
    x,y = R2_r.base_scalars()
    dx,dy = R2_r.base_oneforms()
    metrics = {'flat':(1,1),'polar':(1,x*x),'sphere':(1,s.sin(x)**2),
               'hyperbolic':(1,s.sinh(x)**2),'warped':(s.exp(2*x),x*x)}
    first,second = map(s.sympify,metrics[kind])
    metric = first*TensorProduct(dx,dx)+second*TensorProduct(dy,dy)
    M = reference.Manifold('oracle',2)
    P = reference.Patch('P',M)
    C = reference.CoordSystem('C',P,sp.symbols('u v',real=True))
    u,v = C.base_scalars()
    du,dv = C.base_oneforms()
    mapping = {s.sympify(x).to_sympy():u,s.sympify(y).to_sympy():v}
    expected = first.to_sympy().xreplace(mapping)*reference.TensorProduct(du,du)+second.to_sympy().xreplace(mapping)*reference.TensorProduct(dv,dv)
    for name in ('metric_to_Christoffel_1st','metric_to_Christoffel_2nd','metric_to_Riemann_components','metric_to_Ricci_components'):
        actual = globals()[name](metric)
        oracle = getattr(reference,name)(expected)
        assert actual.shape == oracle.shape
        for indices in product(*(range(n) for n in actual.shape)):
            a,b = actual[indices],oracle[indices]
            assert sp.simplify(a.to_sympy().xreplace(mapping)-b) == 0
    gamma = metric_to_Christoffel_2nd(metric)
    assert BaseCovarDerivativeOp(R2_r,0,gamma)(metric) == 0
    assert BaseCovarDerivativeOp(R2_r,1,gamma)(metric) == 0
    if kind in ('flat','polar'):
        curvature = metric_to_Riemann_components(metric)
        assert all(curvature[indices] == 0 for indices in product(range(2),repeat=4))
    if kind == 'sphere':
        assert metric_to_Ricci_components(metric)[0,0] == 1


def test_connection_leibniz_torsion_and_contractions():
    r,t = R2_p.base_scalars()
    er,et = R2_p.base_vectors()
    dr,dt = R2_p.base_oneforms()
    metric = TensorProduct(dr,dr)+r*r*TensorProduct(dt,dt)
    gamma = metric_to_Christoffel_2nd(metric)
    radial = BaseCovarDerivativeOp(R2_p,0,gamma)
    angular = BaseCovarDerivativeOp(R2_p,1,gamma)
    assert radial(et) == et/r
    assert angular(er) == et/r
    assert angular(et) == -r*er
    assert radial(r*et) == et+r*radial(et)
    assert CovarDerivativeOp(r*er,gamma)(et) == et
    assert radial(dr)(er)+dr(radial(er)) == radial(dr(er))


def test_three_dimensional_flat_spherical_metric():
    r,theta,phi = R3_s.base_scalars()
    dr,dt,dp = R3_s.base_oneforms()
    metric = TensorProduct(dr,dr)+r*r*TensorProduct(dt,dt)+r*r*s.sin(theta)**2*TensorProduct(dp,dp)
    curvature = metric_to_Riemann_components(metric)
    assert curvature.shape == (3,3,3,3)
    assert all(curvature[indices] == 0 for indices in product(range(3),repeat=4))
    gamma = metric_to_Christoffel_2nd(metric)
    assert BaseCovarDerivativeOp(R3_s,2,gamma)(metric) == 0


def test_integral_curve_series_and_equations():
    x,y = R2_r.base_scalars()
    ex,ey = R2_r.base_vectors()
    t = s.Symbol('t')
    X = -y*ex+x*ey
    start = R2_r.point((1,0))
    series = intcurve_series(X,t,start,5)
    assert series == s.Matrix([1-t*t/2+t**4/24,t-t**3/6])
    terms = intcurve_series(X,t,start,5,coeffs=True)
    assert sum(terms,s.Matrix.zeros(2,1)) == series
    equations,initial = intcurve_diffequ(X,t,start)
    f0,f1 = s.Function('f_0')(t),s.Function('f_1')(t)
    assert equations == [f0.diff(t)+f1,f1.diff(t)-f0]
    assert initial == [f0.subs(t,0)-1,f1.subs(t,0)]


def test_invalid_tensor_domains_and_dimensions():
    x,y = R2_r.base_scalars()
    ex,ey = R2_r.base_vectors()
    dx,dy = R2_r.base_oneforms()
    with pytest.raises(ValueError):
        Manifold('bad',0)
    with pytest.raises(ValueError):
        R2_r.point((1,))
    with pytest.raises(ValueError):
        Differential(ex)
    with pytest.raises(ValueError):
        dx*dy
    with pytest.raises(ValueError):
        dx+ex
    with pytest.raises(ValueError):
        metric_to_Christoffel_1st(TensorProduct(dx,dy))
    with pytest.raises(ValueError):
        twoform_to_matrix(dx)
    with pytest.raises(ValueError):
        Differential(TensorProduct(dx,dy))
    with pytest.raises(ValueError):
        metric_to_Christoffel_2nd(TensorProduct(dx,dx))


def test_all_three_packages_without_runtime_sympy():
    subprocess.run([sys.executable,'-c',r'''
import sys
class NoSympy:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):
   raise RuntimeError('runtime reference dependency: '+name)
sys.meta_path.insert(0,NoSympy())
import symox as s
from symox.diffgeom import *
from symox.diffgeom.rn import R2_r
from symox.liealgebras import CartanType
from symox.categories import Object,NamedMorphism,Diagram,xypic_draw_diagram
x,y=R2_r.base_scalars();dx,dy=R2_r.base_oneforms()
metric=TensorProduct(dx,dx)+s.sin(x)**2*TensorProduct(dy,dy)
assert metric_to_Ricci_components(metric)[0,0]==1
assert Differential(Differential(x*y))==0
assert len(CartanType('E8').positive_roots())==120
A,B=Object('A'),Object('B');f=NamedMorphism(A,B,'f')
assert '\\xymatrix' in xypic_draw_diagram(Diagram([f]))
assert 'sympy' not in sys.modules
'''],check=True)
