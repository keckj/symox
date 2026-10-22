"""Differential and physical identities for structural, muscle and Fock algorithms."""
import importlib
import math
import pytest
import symox as S
import sympy as R


def same(actual, expected, locals=None):
    expression=R.sympify(str(actual),locals=locals or {})-expected
    assert R.simplify(expression)==0

@pytest.mark.parametrize('name',[
 'FiberForceLengthActiveDeGroote2016','FiberForceLengthPassiveDeGroote2016',
 'FiberForceLengthPassiveInverseDeGroote2016','FiberForceVelocityDeGroote2016',
 'FiberForceVelocityInverseDeGroote2016','TendonForceLengthDeGroote2016','TendonForceLengthInverseDeGroote2016'])
def test_muscle_curves_and_derivatives(name):
    native=getattr(importlib.import_module('symox.physics.biomechanics'),name)
    reference=getattr(importlib.import_module('sympy.physics.biomechanics'),name)
    x=S.Symbol('x');rx=R.Symbol('x')
    ours=native.with_defaults(x);ref=reference.with_defaults(rx)
    for value in (.2,.8,1.,1.4):
        assert float(ours.doit().subs(x,value).evalf())==pytest.approx(float(ref.doit().subs(rx,value).evalf()),rel=1e-11,abs=1e-11)
        assert float(ours.fdiff().subs(x,value).evalf())==pytest.approx(float(R.diff(ref.doit(),rx).subs(rx,value).evalf()),rel=1e-10,abs=1e-10)

@pytest.mark.parametrize('forward,inverse,values',[
 ('FiberForceLengthPassiveDeGroote2016','FiberForceLengthPassiveInverseDeGroote2016',(.8,1.,1.3)),
 ('FiberForceVelocityDeGroote2016','FiberForceVelocityInverseDeGroote2016',(-1.,0.,1.)),
 ('TendonForceLengthDeGroote2016','TendonForceLengthInverseDeGroote2016',(1.,1.02,1.04))])
def test_curve_inverse_pairs(forward,inverse,values):
    B=importlib.import_module('symox.physics.biomechanics')
    for x in values:
        force=getattr(B,forward).with_defaults(x).doit()
        assert float(getattr(B,inverse).with_defaults(force).doit().evalf())==pytest.approx(x,abs=1e-12)

def test_activation_and_rigid_muscle_loads():
    from symox.physics.biomechanics import FirstOrderActivationDeGroote2016,MusculotendonDeGroote2016,ZerothOrderActivation
    from symox.physics.mechanics import Point,ReferenceFrame,LinearPathway
    activation=FirstOrderActivationDeGroote2016.with_defaults('muscle')
    for a,e in ((.2,.8),(.8,.2),(.5,.5)):
        derivative=float(activation.rhs()[0].subs({activation.a:a,activation.e:e}).evalf())
        switch=(1+math.tanh(10*(e-a)))/2;factor=.5+1.5*a
        assert derivative==pytest.approx((e-a)*(switch/(.015*factor)+(1-switch)*factor/.06))
    N=ReferenceFrame('N');O=Point('O');P=O.locatenew('P',S.Rational(3,10)*N.x)
    muscle=MusculotendonDeGroote2016('m',LinearPathway(O,P),ZerothOrderActivation('m'),tendon_slack_length=.2,optimal_fiber_length=.1,peak_isometric_force=100,maximal_fiber_velocity=1,optimal_pennation_angle=0,fiber_damping_coefficient=0,with_defaults=True)
    force=float(muscle.force.subs(muscle.e,.5).evalf())
    assert -55<force<-45
    loads=muscle.to_loads();assert float(loads[0][1].dot(N.x).subs(muscle.e,.5).evalf())>0;assert S.simplify((loads[0][1]+loads[1][1]).dot(N.x))==0
    assert muscle.x.shape==(0,1)

@pytest.mark.parametrize('form',[1,2,3,4])
def test_compliant_muscle_formulations(form):
    from symox.physics.biomechanics import MusculotendonDeGroote2016,ZerothOrderActivation
    from symox.physics.mechanics import Point,ReferenceFrame,LinearPathway
    N=ReferenceFrame('N');O=Point('O');P=O.locatenew('P',S.Rational(3,10)*N.x)
    m=MusculotendonDeGroote2016('m',LinearPathway(O,P),ZerothOrderActivation('m'),musculotendon_dynamics=form,tendon_slack_length=S.Rational(1,5),optimal_fiber_length=S.Rational(1,10),peak_isometric_force=100,maximal_fiber_velocity=1,optimal_pennation_angle=0,fiber_damping_coefficient=0,with_defaults=True)
    assert m.x.rows==m.rhs().rows==1
    assert m.algebraic_equations.rows==(1 if form>=3 else 0)
    assert m.r.rows==(2 if form>=3 else 1)

@pytest.mark.parametrize('case',['simple','cantilever','uniform','ramp'])
def test_beam_reactions_and_curves_against_reference(case):
    from symox.physics.continuum_mechanics.beam import Beam
    from sympy.physics.continuum_mechanics.beam import Beam as RefBeam
    beams=[Beam(10,2,3),RefBeam(10,2,3)]
    for b in beams:
        if case=='cantilever':reactions=b.apply_support(0,'fixed')
        else:reactions=(b.apply_support(0,'pin'),b.apply_support(10,'roller'))
        if case=='uniform':b.apply_load(-3,0,0,10)
        elif case=='ramp':b.apply_load(-2,2,1,7)
        else:b.apply_load(-12,5,-1)
        b.solve_for_reaction_loads(*reactions)
    a,b=beams
    for key,value in a.reaction_loads.items():same(value,b.reaction_loads[R.Symbol(str(key))])
    for method in ('shear_force','bending_moment','slope','deflection'):
        for x in (1,3,6,9):same(getattr(a,method)().subs(a.variable,x),getattr(b,method)().subs(b.variable,x))

def test_truss_triangle_and_indeterminacy():
    from symox.physics.continuum_mechanics.truss import Truss
    t=Truss();t.add_node(('A',0,0),('B',4,0),('C',2,3));t.add_member(('AB','A','B'),('AC','A','C'),('BC','B','C'))
    t.apply_support(('A','pinned'),('B','roller'));t.apply_load(('C',12,270));t.solve()
    assert t.reaction_loads=={'R_A_x':0,'R_A_y':6,'R_B_y':6}
    assert S.simplify(t.internal_forces['AC']+2*S.sqrt(13))==0
    assert t.internal_forces['AB']==4
    t.remove_member('AB')
    with pytest.raises(ValueError):t.solve()

def test_polygonal_and_parabolic_cables_equilibrium():
    from symox.physics.continuum_mechanics.cable import Cable
    c=Cable(('A',0,4),('B',6,4));c.apply_load(-1,('P',3,0,10,270));c.solve()
    assert set(c.tension.values())=={S.Rational(25,4)}
    assert c.reaction_loads=={S.Symbol('R_A_x'):-S.Rational(15,4),S.Symbol('R_A_y'):5,S.Symbol('R_B_x'):S.Rational(15,4),S.Symbol('R_B_y'):5}
    assert c.length==10
    c=Cable(('A',0,5),('B',10,5));c.apply_load(0,('w',2));c.solve(5,0)
    assert c.tension_at(5)==5
    assert c.reaction_loads[S.Symbol('R_A_y')]==10
    assert c.reaction_loads[S.Symbol('R_B_y')]==10
    assert S.simplify(c._cable_eqn.subs(S.Symbol('x'),10)-5)==0

def test_three_hinged_arch_reference_example():
    from symox.physics.continuum_mechanics.arch import Arch
    a=Arch((0,0),(10,0),crown_x=5,crown_y=5);a.apply_load(0,'C',start=3,end=5,mag=-10);a.solve()
    assert a.reaction_force==dict(zip(S.symbols('R_A_x R_A_y R_B_x R_B_y'),[8,12,-8,8]))
    assert a.bending_moment_at(0)==a.bending_moment_at(5)==a.bending_moment_at(10)==0

@pytest.mark.parametrize('count',[2,4,6])
def test_gamma_traces_against_explicit_dirac_matrices(count):
    from itertools import product
    from symox.physics.hep.gamma_matrices import GammaMatrix as G,LorentzIndex as L,gamma_trace
    from symox.physics.matrices import mgamma,minkowski_tensor
    from symox.tensor.tensor import tensor_indices
    indices=tensor_indices(' '.join('a'+str(i) for i in range(count)),L)
    expr=G(indices[0])
    for i in indices[1:]:expr*=G(i)
    trace=gamma_trace(expr)
    actual_array=trace.replace_with_arrays({L:minkowski_tensor},indices)
    # Component evaluation independently checks the recursive Clifford signs.
    for values in list(product(range(4),repeat=count))[:100:7]:
        actual=actual_array
        for v in values:actual=actual[v]
        matrix=S.eye(4)
        for v in values:matrix*=mgamma(v)
        assert actual==matrix.trace()

def test_gamma_contracted_identity_and_slash_square():
    from symox.physics.hep.gamma_matrices import GammaMatrix as G,LorentzIndex as L,kahane_simplify,gamma_trace,simplify_gpgp
    from symox.tensor.tensor import tensor_indices,TensorHead,canon_bp
    mu,a=tensor_indices('mu a',L)
    assert kahane_simplify(G(mu)*G(-mu))==4*S.eye(4)
    assert kahane_simplify(G(mu)*G(a)*G(-mu))==-2*G(a)
    assert gamma_trace(G(mu)*G(-mu))==16
    p=TensorHead('p',[L]);nu=tensor_indices('nu',L)
    assert canon_bp(simplify_gpgp(p(mu)*G(-mu)*p(nu)*G(-nu)))==canon_bp(p(mu)*p(-mu))

def test_bosonic_ladders_basis_and_normal_order():
    from symox.physics.secondquant import B,Bd,BKet,Dagger,apply_operators,matrix_rep,VarBosonicBasis,Commutator,wicks,NO
    assert apply_operators(B(0)*BKet([3]))==S.sqrt(3)*BKet([2])
    assert apply_operators(Bd(1)*BKet([1,2]))==S.sqrt(3)*BKet([1,3])
    assert Dagger(B(0))==Bd(0)
    assert Commutator(B(0),Bd(0))==1
    assert wicks(B(0)*Bd(0))==Bd(0)*B(0)+1
    assert NO(B(0)*Bd(0)).doit()==Bd(0)*B(0)
    assert matrix_rep(B(0),VarBosonicBasis(3))==S.Matrix([[0,1,0],[0,0,S.sqrt(2)],[0,0,0]])

def test_fermionic_ladders_anticommutation_and_wick_signs():
    from symox.physics.secondquant import F,Fd,FKet,Dagger,apply_operators,wicks,NO,AntiSymmetricTensor
    assert apply_operators(F(1)*FKet([0,1]))==-FKet([0])
    assert apply_operators(Fd(1)*FKet([0]))==-FKet([0,1])
    assert apply_operators(Fd(1)*FKet([0,1]))==0
    assert Dagger(F(0))==Fd(0)
    # Reference Fermi-vacuum convention: numeric indices are ambiguous with
    # respect to the Fermi level, so contractions keep KroneckerDelta(0, _a).
    assert str(wicks(F(0)*Fd(0)+Fd(0)*F(0)))=='KroneckerDelta(0, _a) + KroneckerDelta(0, _i)'
    assert str(wicks(F(0)*F(1)*Fd(1)*Fd(0),keep_only_fully_contracted=True))=='KroneckerDelta(0, _a)*KroneckerDelta(1, _a)'
    assert str(wicks(F(0)*F(1)*Fd(0)*Fd(1),keep_only_fully_contracted=True))=='-KroneckerDelta(0, _a)*KroneckerDelta(1, _a)'
    assert str(NO(F(0)*Fd(1)))=='-NO(CreateFermion(1)*AnnihilateFermion(0))'
    i,j=S.symbols('i j')
    assert AntiSymmetricTensor('v',(i,j),())==-AntiSymmetricTensor('v',(j,i),())
    assert AntiSymmetricTensor('v',(i,i),())==0

@pytest.mark.parametrize('module',[
 'physics.biomechanics','physics.continuum_mechanics.beam','physics.continuum_mechanics.truss',
 'physics.continuum_mechanics.cable','physics.continuum_mechanics.arch','physics.hep.gamma_matrices','physics.secondquant'])
def test_remaining_families_public_inventory(module):
    ref=importlib.import_module('sympy.'+module);native=importlib.import_module('symox.'+module)
    names=getattr(ref,'__all__',None) or [n for n in dir(ref) if not n.startswith('_') and callable(getattr(ref,n)) and getattr(getattr(ref,n),'__module__','')==ref.__name__]
    assert all(hasattr(native,n) for n in names)
