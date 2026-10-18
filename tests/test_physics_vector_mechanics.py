"""Moving-frame identities, mechanics models and reference rotations."""
import pytest
import symox as S
import symox.physics.vector as V
import symox.physics.mechanics as M
import sympy as R
import sympy.physics.vector as RV
import sympy.physics.mechanics as RM


def same(a,b):assert S.trigsimp(S.simplify(a-S.sympify(b)))==0

def matrix_same(a,b):
    assert a.shape==b.shape
    for x,y in zip(a,b):same(x,y)


@pytest.mark.parametrize('order',['123','231','312','132','213','321','121','131','212','232','313','323'])
@pytest.mark.parametrize('kind',['Body','Space'])
def test_rotations_match_reference(order,kind):
    q=S.symbols('a b c');qr=R.symbols('a b c')
    N=V.ReferenceFrame('N');A=N.orientnew('A',kind,q,order)
    Nr=RV.ReferenceFrame('N');Ar=Nr.orientnew('A',kind,qr,order)
    reference=Ar.dcm(Nr)
    matrix_same(A.dcm(N),S.Matrix([[S.sympify(str(reference[i,j])) for j in range(3)] for i in range(3)]))
    matrix_same(A.dcm(N)*N.dcm(A),S.eye(3))


def test_vector_products_dyadics_and_rotating_derivatives():
    N=V.ReferenceFrame('N');q=V.dynamicsymbols('q');t=V.dynamicsymbols._t
    A=N.orientnew('A','Axis',(q,N.z))
    assert N.x.cross(N.y)==N.z and N.x.dot(N.y)==0
    assert N.x.outer(N.y).dot(N.y)==N.x
    assert N.x.dot(N.x.outer(N.y))==N.y
    assert A.ang_vel_in(N)==q.diff(t)*N.z
    assert V.time_derivative(A.x,N)==q.diff(t)*A.y
    assert V.time_derivative(A.x,N,2)==q.diff(t,2)*A.y-q.diff(t)**2*A.x
    assert A.x.dt(A)==0
    matrix_same(A.x.to_matrix(N),S.Matrix([S.cos(q),S.sin(q),0]))
    same(A.x.dot(N.x),S.cos(q))
    assert V.express(A.x,N)==S.cos(q)*N.x+S.sin(q)*N.y
    dyad=A.x.outer(A.y)
    assert dyad.dt(N)==q.diff(t)*(A.y.outer(A.y)-A.x.outer(A.x))
    assert (2*N.x+3*N.y).diff(S.Symbol('x'),N)==0
    assert S.simplify((S.sin(q)**2+S.cos(q)**2)*N.x)==N.x
    with pytest.raises(ValueError):V.Vector(0).normalize()
    with pytest.raises(ValueError):N.dcm(V.ReferenceFrame('unconnected'))


def test_point_kinematic_theorems():
    N=V.ReferenceFrame('N');q,u=V.dynamicsymbols('q u');t=V.dynamicsymbols._t
    A=N.orientnew('A','Axis',(q,N.z));O=V.Point('O');O.set_vel(N,0)
    P=O.locatenew('P',2*A.x)
    assert P.vel(N)==2*q.diff(t)*A.y
    assert P.acc(N)==2*q.diff(t,2)*A.y-2*q.diff(t)**2*A.x
    assert P.v2pt_theory(O,N,A)==2*q.diff(t)*A.y
    assert P.a2pt_theory(O,N,A)==2*q.diff(t,2)*A.y-2*q.diff(t)**2*A.x
    Q=O.locatenew('Q',q*N.x);Q.set_vel(N,u*N.x)
    assert Q.partial_velocity(N,u)==N.x
    assert V.partial_velocity([u*N.x+q*N.y],[u,q],N)==[[N.x,N.y]]
    assert Q.pos_from(P)==q*N.x-2*A.x
    assert V.get_motion_params(N,acceleration=2*N.x)[1]==2*t*N.x
    assert V.get_motion_params(N,acceleration=2*N.x)[2]==t**2*N.x


def test_fields_coordinates_and_printing(capsys):
    N=V.ReferenceFrame('N');x,y,z=N.varlist
    field=2*x*N.x+2*y*N.y+2*z*N.z
    assert V.gradient(x*x+y*y+z*z,N)==field
    same(V.divergence(field,N),6)
    assert V.curl(field,N)==0 and V.is_conservative(field,N)
    same(V.scalar_potential(field,N),x*x+y*y+z*z)
    assert V.is_solenoidal(y*N.x-x*N.y,N)
    O=V.Point('O');P=O.locatenew('P',N.x);Q=O.locatenew('Q',3*N.x)
    same(V.scalar_potential_difference(field,N,P,Q,O),8)
    assert isinstance(N[0],V.CoordinateSym)
    assert V.vlatex(N.x)
    V.vprint(N.x);assert 'N.x' in capsys.readouterr().out


def oscillator():
    q,u=M.dynamicsymbols('q u');t=M.dynamicsymbols._t;m,k,c=S.symbols('m k c',positive=True)
    N=M.ReferenceFrame('N');O=M.Point('O');O.set_vel(N,0);P=O.locatenew('P',q*N.x)
    body=M.Particle('body',P,m);body.potential_energy=k*q*q/2
    return q,u,t,m,k,c,N,O,P,body


def test_lagrange_kane_and_linearization():
    q,u,t,m,k,c,N,O,P,body=oscillator()
    lagrange=M.LagrangesMethod(M.Lagrangian(N,body),[q],[(P,-c*q.diff(t)*N.x)],frame=N)
    same(lagrange.form_lagranges_equations()[0],m*q.diff(t,2)+c*q.diff(t)+k*q)
    matrix_same(lagrange.rhs(),S.Matrix([q.diff(t),(-c*q.diff(t)-k*q)/m]))
    kane=M.KanesMethod(N,[q],[u],[q.diff(t)-u])
    active,inertia=kane.kanes_equations([body],[(P,(-k*q-c*u)*N.x)])
    same(active[0],-k*q-c*u);same(inertia[0],-m*u.diff(t))
    matrix_same(kane.rhs(),S.Matrix([u,(-k*q-c*u)/m]))
    A,B=kane.linearize(A_and_B=True)
    matrix_same(A,S.Matrix([[0,1],[-k/m,-c/m]]));assert B.shape==(2,0)
    assert kane.kindiffdict()=={q.diff(t):u}


def test_pendulum_energy_and_lagrange_equation():
    q=M.dynamicsymbols('theta');t=M.dynamicsymbols._t;m,l,g=S.symbols('m l g',positive=True)
    N=M.ReferenceFrame('N');A=N.orientnew('A','Axis',[q,N.z]);O=M.Point('O');O.set_vel(N,0)
    P=O.locatenew('P',-l*A.y);body=M.Particle('pendulum',P,m)
    body.potential_energy=-m*g*l*S.cos(q)
    same(body.kinetic_energy(N),m*l*l*q.diff(t)**2/2)
    equation=M.LagrangesMethod(M.Lagrangian(N,body),[q]).form_lagranges_equations()[0]
    same(equation,m*l*l*q.diff(t,2)+m*g*l*S.sin(q))
    assert M.linear_momentum(N,body)==body.mass*P.vel(N)
    assert M.angular_momentum(O,N,body)==m*l*l*q.diff(t)*N.z


def test_rigid_body_inertia_and_parallel_axes():
    q=M.dynamicsymbols('spin');t=M.dynamicsymbols._t
    N=M.ReferenceFrame('N');A=N.orientnew('A','Axis',[q,N.z]);P=M.Point('P');P.set_vel(N,0)
    I=M.inertia(A,1,2,3);body=M.RigidBody('body',P,A,2,(I,P))
    same(body.kinetic_energy(N),S.Rational(3,2)*q.diff(t)**2)
    assert body.angular_momentum(P,N)==3*q.diff(t)*N.z
    Q=P.locatenew('Q',2*A.x)
    matrix_same(body.parallel_axis(Q).to_matrix(A),S.diag(1,10,11))
    record=M.Inertia.from_inertia_scalars(P,A,1,2,3)
    assert record.dyadic==I and record.point is P
    assert M.Force(body,2*N.x).point is P
    assert M.Torque(body,N.z).frame is A


def test_spring_damper_and_pathway_conservation():
    q,u,t,m,k,c,N,O,P,body=oscillator()
    path=M.LinearPathway(O,P)
    spring=M.LinearSpring(k,path);damper=M.LinearDamper(c,path)
    # The geometric spring force reduces to -k*q for either sign of q.
    force=spring.to_loads()[1].force
    same(force.dot(N.x),-k*q)
    assert sum((f for p,f in spring.to_loads()),M.Vector(0))==0
    same(damper.to_loads()[1].force.dot(N.x),-c*q.diff(t))
    duff=M.DuffingSpring(k,c,path)
    assert len(duff.to_loads())==2
    torque=M.TorqueActuator(3,N.z,N)
    assert torque.to_loads()==[M.Torque(N,3*N.z)]
    Q=O.locatenew('Q',N.y)
    sphere=M.WrappingSphere(1,O);Rpoint=O.locatenew('R',N.x)
    same(sphere.geodesic_length(Rpoint,Q),S.pi/2)
    wrap=M.WrappingPathway(Rpoint,Q,sphere)
    assert sum((f for p,f in wrap.to_loads(2)),M.Vector(0))==0


def test_prismatic_joint_system_and_symbolic_system():
    N=M.ReferenceFrame('N');O=M.Point('O');wall=M.RigidBody('wall',O,N,0)
    B=M.RigidBody('B',mass=2);q,u=M.dynamicsymbols('slide velocity')
    joint=M.PrismaticJoint('slider',wall,B,q,u,joint_axis=N.x)
    system=M.System.from_newtonian(wall);system.add_joints(joint);system.add_loads(M.Force(B,-3*q*N.x))
    same(system.form_eoms()[0],-2*u.diff(M.dynamicsymbols._t)-3*q)
    matrix_same(system.rhs(),S.Matrix([u,-S.Rational(3,2)*q]))
    symbolic=M.SymbolicSystem([q],[-3*q],speeds=[u],mass_matrix=S.Matrix([[2]]),coordinate_derivatives=[u])
    symbolic.compute_explicit_form();matrix_same(symbolic.comb_explicit_rhs,S.Matrix([u,-S.Rational(3,2)*q]))
    assert M.find_dynamicsymbols(q+q.diff(M.dynamicsymbols._t))=={q,q.diff(M.dynamicsymbols._t)}


def test_inventory_matches_reference():
    assert set(V.__all__)==set(RV.__all__)
    assert set(M.__all__)==set(RM.__all__)
