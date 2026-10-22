"""Optical identities, control responses and finite quantum algorithms."""
import importlib
import subprocess
import sys
import numpy as np
import pytest
import symox as S
import symox.physics.optics as O
import symox.physics.control as C
import symox.physics.quantum as Q
import sympy as R
import sympy.physics.optics as RO
import sympy.physics.control as RC
import sympy.physics.quantum as RQ


def same(a,b):
    a=S.sympify(a)
    if isinstance(b,str):b=S.parse_expr(b,local_dict={str(v):v for v in a.free_symbols})
    difference=S.trigsimp(S.simplify(a-S.sympify(b)))
    if difference!=0:assert R.simplify(R.trigsimp(R.sympify(str(difference),locals={str(v):R.Symbol(str(v),real=True) for v in difference.free_symbols}),method='fu'))==0


def matrices(a,b):
    if hasattr(a,'as_explicit'):a=a.as_explicit()
    assert a.shape==b.shape
    for x,y in zip(a,b):same(x,str(y) if isinstance(y,R.Expr) else y)


@pytest.mark.parametrize('name,args',[('FreeSpace',(3,)),('FlatRefraction',(1,2)),('CurvedRefraction',(3,1,2)),('FlatMirror',()),('CurvedMirror',(3,)),('ThinLens',(3,))])
def test_abcd_against_reference(name,args):matrices(getattr(O,name)(*args),getattr(RO,name)(*args))


def test_ray_and_gaussian_beam_propagation():
    ray=O.FreeSpace(2)*O.ThinLens(3)*O.GeometricRay(1,0)
    matrices(ray,S.Matrix([S.Rational(1,3),-S.Rational(1,3)]))
    beam=O.BeamParameter(1,2,w=3)
    same(beam.z_r,9*S.pi);same(beam.w_0,3);same(beam.radius,2+81*S.pi**2/2)
    propagated=O.FreeSpace(5)*beam
    same(propagated.q,7+9*S.I*S.pi)
    s,z,f=2,3,4
    image,zout,magnification=O.gaussian_conj(s,z,f)
    after=O.ThinLens(f)*O.BeamParameter(1,s,z_r=z)
    same(image,-S.re(after.q));same(zout,S.im(after.q));same(magnification**2,zout/z)
    same(O.geometric_conj_ab(2,3),S.Rational(6,5))
    same(O.waist2rayleigh(3,1),9*S.pi)
    same(O.rayleigh2waist(9*S.pi,1),3)


@pytest.mark.parametrize('name,args',[
 ('jones_vector',(S.pi/6,S.pi/8)),('stokes_vector',(S.pi/6,S.pi/8)),
 ('linear_polarizer',(S.pi/6,)),('phase_retarder',(S.pi/7,S.pi/3)),
 ('half_wave_retarder',(S.pi/6,)),('quarter_wave_retarder',(S.pi/6,)),
 ('transmissive_filter',(S.Rational(1,2),)),('reflective_filter',(S.Rational(1,2),)),
 ('polarizing_beam_splitter',()),
])
def test_polarization_against_reference(name,args):
    reference_args=[R.sympify(str(a)) for a in args]
    matrices(getattr(O,name)(*args),getattr(RO,name)(*reference_args))


def test_jones_stokes_mueller_and_snell():
    psi,chi=S.symbols('psi chi',real=True)
    matrices(O.jones_2_stokes(O.jones_vector(psi,chi)),O.stokes_vector(psi,chi))
    matrices(O.mueller_matrix(S.eye(2)),S.eye(4))
    # Like the reference, a numeric angle is converted with float().
    assert str(O.refraction_angle(S.pi/6,1,S.Rational(3,2)))=='0.339836909454122'
    same(O.brewster_angle(1,S.Rational(3,2)),S.atan(S.Rational(3,2)))
    same(O.critical_angle(2,1),S.pi/6)
    with pytest.raises(ValueError,match='total internal reflection'):O.refraction_angle(S.pi/3,2,1)
    with pytest.raises(ValueError):O.critical_angle(1,2)
    same(O.lens_formula(focal_length=2,u=-6),3)
    same(O.mirror_formula(focal_length=2,u=6),3)
    same(O.lens_makers_formula(S.Rational(3,2),1,2,-2),2)
    water=O.Medium('water',n=S.Rational(4,3));same(water.refractive_index,S.Rational(4,3))
    wave=O.TWave(2,3,phase=S.pi/2,n=1)
    same(wave.time_period,S.Rational(1,3));same(wave.angular_velocity,6*S.pi)
    same((O.TWave(1,1,n=1)+O.TWave(1,1,n=1)).amplitude,2)


def test_transfer_function_interconnections_and_discretization():
    s=S.Symbol('s');g=C.TransferFunction(1,s+1,s);h=C.TransferFunction(2,s+2,s)
    same((g*h).to_expr(),2/((s+1)*(s+2)))
    same((g+h).to_expr(),1/(s+1)+2/(s+2))
    same(C.Feedback(g).to_expr(),1/(s+2))
    assert g.poles()==[-1] and g.zeros()==[] and g.is_stable()
    same(g.dc_gain(),1)
    assert g.is_strictly_proper and g.is_proper
    assert C.bilinear(g,1)==([S.Rational(1,3),S.Rational(1,3)],[1,-S.Rational(1,3)])
    same(C.PIDController(2,3,4).to_expr(),2+3/s+4*s)
    same(C.TransferFunction.from_coeff_lists([1,2],[1,3,4],s).to_expr(),(s+2)/(s*s+3*s+4))
    ss=g.rewrite(C.StateSpace)
    matrices(ss.A,S.Matrix([[-1]]));matrices(ss.B,S.Matrix([[1]]))
    same(ss.rewrite(C.TransferFunction)[0][0].to_expr(),g.to_expr())
    assert ss.is_controllable() and ss.is_observable()
    matrices(ss.dsolve([1]),S.Matrix([S.exp(-S.Symbol('t'))]))


def test_control_time_and_frequency_responses():
    s=S.Symbol('s');g=C.TransferFunction(1,s+1,s)
    t,y=C.step_response_numerical_data(g,upper_limit=4,n=31)
    np.testing.assert_allclose(y,1-np.exp(-t),atol=1e-12)
    t,y=C.impulse_response_numerical_data(g,upper_limit=4,n=31)
    np.testing.assert_allclose(y,np.exp(-t),atol=1e-12)
    t,y=C.ramp_response_numerical_data(g,upper_limit=4,n=31)
    np.testing.assert_allclose(y,t-1+np.exp(-t),atol=1e-12)
    w,mag=C.bode_magnitude_numerical_data(g,-2,2,n=31)
    np.testing.assert_allclose(mag,-10*np.log10(1+w*w),atol=1e-10)
    w,phase=C.bode_phase_numerical_data(g,-2,2,n=31)
    np.testing.assert_allclose(phase,-np.arctan(w),atol=1e-10)
    assert C.gain_margin(g)==S.oo
    poles,zeros=C.pole_zero_numerical_data(g);assert poles==[-1+0j] and zeros==[]


def test_mimo_and_plotting(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    s=S.Symbol('s');A=C.TransferFunctionMatrix.from_Matrix(S.diag(1/(s+1),2/(s+2)),s)
    B=C.TransferFunctionMatrix.from_Matrix(S.eye(2),s)
    matrices((A*B).to_expr(),A.to_expr())
    matrices(C.MIMOFeedback(A,B).to_expr(),S.diag(1/(s+2),2/(s+4)))
    fig=C.step_response_plot(A[0,0],n=20,show=False);fig.savefig(tmp_path/'step.png')
    assert (tmp_path/'step.png').stat().st_size>100


def test_operator_order_duality_and_inner_products():
    A,B=Q.Operator('A'),Q.Operator('B')
    assert A*B!=B*A
    assert Q.Commutator(A,B).doit()==A*B-B*A
    assert Q.AntiCommutator(A,B).doit()==A*B+B*A
    assert Q.Dagger(A*B)==Q.Dagger(B)*Q.Dagger(A)
    assert Q.Dagger(Q.Dagger(A))==A
    assert Q.Dagger(Q.HermitianOperator('H'))==Q.HermitianOperator('H')
    assert Q.Dagger(Q.OrthogonalKet(0))*Q.OrthogonalKet(1)==0
    assert Q.Dagger(Q.OrthogonalKet(1))*Q.OrthogonalKet(1)==1
    assert Q.ComplexSpace(2).dimension==2
    assert (Q.ComplexSpace(2)*Q.ComplexSpace(3)).dimension==6
    assert (Q.ComplexSpace(2)+Q.ComplexSpace(3)).dimension==5


def test_quantum_differential_operator_and_oscillator():
    from symox.physics.quantum.sho1d import a,ad,N,SHOKet
    assert Q.qapply(a*SHOKet(0))==0
    assert Q.qapply(a*SHOKet(2))==S.sqrt(2)*SHOKet(1)
    assert Q.qapply(ad*SHOKet(2))==S.sqrt(3)*SHOKet(3)
    assert Q.qapply(N*SHOKet(3))==3*SHOKet(3)
    assert Q.Dagger(a)==ad
    assert Q.Commutator(a,ad).doit()==1
    matrices(Q.represent(a,ndim=3),S.Matrix([[0,1,0],[0,0,S.sqrt(2)],[0,0,0]]))
    x=S.Symbol('x');f=S.Function('f')(x);psi=Q.Wavefunction(x*x,x)
    result=Q.qapply(Q.DifferentialOperator(S.diff(f,x),f)*psi)
    same(result.expr,2*x)
    from symox.physics.quantum.cartesian import X,Px
    same(Q.Commutator(X,Px).doit(),S.I*Q.hbar)


def test_gates_bell_measurements_and_grover():
    from symox.physics.quantum.qubit import Qubit,QubitBra,measure_all,measure_partial,qubit_to_matrix,matrix_to_qubit
    from symox.physics.quantum.gate import H,CNOT,X,T
    from symox.physics.quantum.grover import apply_grover
    psi=Q.qapply(CNOT(0,1)*H(0)*Qubit('00'))
    assert psi==(Qubit('00')+Qubit('11'))/S.sqrt(2)
    assert measure_all(psi)==[(Qubit('00'),S.Rational(1,2)),(Qubit('11'),S.Rational(1,2))]
    assert measure_partial(psi,[0])==[(Qubit('00'),S.Rational(1,2)),(Qubit('11'),S.Rational(1,2))]
    assert Q.Dagger(Qubit('01'))==QubitBra('01')
    assert Q.Dagger(Qubit('01'))*Qubit('01')==1
    assert Q.qapply(H(0)*H(0)*Qubit('0'))==Qubit('0')
    assert Q.qapply(Q.Dagger(T(0))*T(0)*Qubit('1'))==Qubit('1')
    assert matrix_to_qubit(qubit_to_matrix(psi))==psi
    assert apply_grover(lambda q:q==Qubit('11'),2,1)==Qubit('11')
    from symox.physics.quantum.qft import QFT,IQFT
    matrices(Q.represent(IQFT(0,2))*Q.represent(QFT(0,2)),S.eye(4))
    from symox.physics.quantum.density import Density
    mixed=Density((Qubit('0'),S.Rational(1,2)),(Qubit('1'),S.Rational(1,2)))
    matrices(mixed.to_matrix(),S.eye(2)/2);same(mixed.entropy(),S.log(2))


def test_public_inventories_and_native_runtime():
    for ours,theirs in [(O,RO),(C,RC),(Q,RQ)]:assert set(ours.__all__)==set(theirs.__all__)
    code="""
import sys
class Block:
 def find_spec(self,name,path=None,target=None):
  if name=='sympy' or name.startswith('sympy.'):raise AssertionError(name)
sys.meta_path.insert(0,Block())
from symox.physics.quantum import qapply
from symox.physics.quantum.qubit import Qubit
from symox.physics.quantum.gate import X
assert qapply(X(0)*Qubit('0'))==Qubit('1')
from symox.physics.control import TransferFunction
from symox import Symbol
s=Symbol('s')
assert TransferFunction(1,s+1,s).dc_gain()==1
"""
    result=subprocess.run([sys.executable,'-c',code],capture_output=True,text=True)
    assert result.returncode==0,result.stderr

@pytest.mark.parametrize('j',[S.Rational(1,2),1,S.Rational(3,2),2])
def test_spin_representations_and_eigenstates(j):
    from symox.physics.quantum.spin import Jx,Jy,Jz,J2,Jplus,Jminus,JxKet,JyKet,JzKet,m_values
    from symox.physics.quantum import represent,qapply,Commutator,hbar
    from sympy.physics.quantum import represent as rr
    from sympy.physics.quantum import spin as rs
    rj=R.Rational(str(j))
    for name,operator in [('Jx',Jx),('Jy',Jy),('Jz',Jz),('J2',J2),('Jplus',Jplus),('Jminus',Jminus)]:
        actual=represent(operator,j=j);reference=rr(getattr(rs,name),j=rj)
        for a,b in zip(actual,reference):assert R.simplify(R.sympify(str(a),locals={'hbar':rs.hbar})-b)==0
    for ket,operator in ((JxKet,Jx),(JyKet,Jy),(JzKet,Jz)):
        for m in m_values(j)[1]:
            vector=represent(ket(j,m));assert (represent(operator,j=j)*vector-hbar*m*vector).applyfunc(S.simplify)==S.zeros(int(2*j+1),1)
    assert Commutator(Jx,Jy).doit()==S.I*hbar*Jz
    assert qapply(Jplus*JzKet(j,j))==0


def test_spin_rotation_and_two_spin_coupling():
    from symox.physics.quantum.spin import Rotation,JzKet,JzKetCoupled,couple,uncouple
    from symox.physics.quantum import TensorProduct,represent,Dagger
    from sympy.physics.quantum.spin import Rotation as RefRotation
    j=S.Rational(1,2);a,b,g=S.symbols('alpha beta gamma',real=True)
    actual=Rotation.D(j,j,-j,a,b,g).doit()
    reference=RefRotation.D(R.Rational(1,2),R.Rational(1,2),-R.Rational(1,2),R.Symbol('alpha'),R.Symbol('beta'),R.Symbol('gamma')).doit()
    assert R.simplify(R.sympify(str(actual),locals={str(v):R.Symbol(str(v)) for v in (a,b,g)})-reference)==0
    product=TensorProduct(JzKet(j,j),JzKet(j,-j))
    assert uncouple(couple(product))==product
    singlet=JzKetCoupled(0,0,(j,j))
    assert couple(uncouple(singlet))==singlet
    rotation=Rotation(S.pi/3,S.pi/4,S.pi/6)
    matrix=represent(rotation,j=j);assert (represent(Dagger(rotation),j=j)*matrix).applyfunc(S.simplify)==S.eye(2)


def test_particle_in_box_and_observable_mappings():
    from symox.physics.quantum.piab import PIABKet,PIABHamiltonian,m,L
    from symox.physics.quantum import Dagger,qapply,represent,hbar
    from symox.physics.quantum.operatorset import operators_to_state,state_to_operators
    from symox.physics.quantum.cartesian import XOp,XKet
    state=PIABKet(2)
    assert qapply(PIABHamiltonian('H')*state)==2*S.pi**2*hbar**2/(m*L**2)*state
    assert Dagger(state)*state==1
    assert Dagger(PIABKet(1))*state==0
    x=S.Symbol('x');assert represent(state)==S.sqrt(2/L)*S.sin(2*S.pi*x/L)
    assert operators_to_state(XOp)==XKet('x')
    assert isinstance(state_to_operators(XKet),XOp)

@pytest.mark.parametrize('N',[15,21,33,35,49,77])
def test_shor_classical_simulation_and_modular_oracle(N):
    from symox.physics.quantum.shor import shor,period_find,CMod,getr
    from symox.physics.quantum import qapply
    from symox.physics.quantum.qubit import Qubit
    a,b=shor(N);assert 1<a<N and a*b==N
    assert pow(2,int(period_find(2,N)),N)==1
    assert getr(1,4,15)==4
    gate=CMod(2,2,3)
    for k in range(4):
        for y in range(4):
            state=Qubit(format((k<<2)|y,'04b'));actual=qapply(gate*state)
            assert actual==Qubit(format((k<<2)|(y^pow(2,k,3)),'04b'))
            assert qapply(gate*actual)==state


def test_quantum_circuit_diagram_and_custom_gate(tmp_path):
    from symox.physics.quantum.circuitplot import circuit_plot,CreateOneQubitGate,CreateCGate,labeller,Mz
    from symox.physics.quantum.gate import CNOT,H
    import matplotlib.pyplot as plt
    custom=CreateOneQubitGate('U');controlled=CreateCGate('V')
    diagram=circuit_plot(Mz(0)*controlled((1,),0)*custom(1)*CNOT(1,0)*H(1),2)
    path=tmp_path/'circuit.png';diagram.savefig(path);assert path.stat().st_size>0
    assert len(diagram.gates)==5;assert labeller(2)==['q_1','q_0']
    plt.close(diagram.figure)


def test_all_documented_physics_modules_import_without_sympy():
    code=r'''
import builtins
original=builtins.__import__
def blocked(name,*args,**kwargs):
    if name=='sympy' or name.startswith('sympy.'):
        raise ImportError('runtime SymPy delegation is forbidden')
    return original(name,*args,**kwargs)
builtins.__import__=blocked
import symox as S
from symox.physics.biomechanics import FirstOrderActivationDeGroote2016
from symox.physics.continuum_mechanics import Beam,Arch,Cable,Truss
from symox.physics.hep.gamma_matrices import gamma_trace,GammaMatrix,LorentzIndex
from symox.physics.secondquant import F,Fd,wicks
from symox.physics.quantum.operatorset import operators_to_state
from symox.physics.quantum.piab import PIABKet
from symox.physics.quantum.circuitplot import labeller
from symox.physics.quantum.shor import shor
from symox.physics.quantum.spin import JzKet
assert wicks(F(0)*Fd(0),keep_only_fully_contracted=True)==1
assert shor(15)[0]*shor(15)[1]==15
assert FirstOrderActivationDeGroote2016.with_defaults('m').rhs().rows==1
assert labeller(2)==['q_1','q_0']
'''
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
