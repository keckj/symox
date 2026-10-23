"""WP-C11 quantum parity: expected values are hard-coded SymPy 1.14 outputs."""
import pytest


# Star imports go into a namespace dict so that pytest does not collect names
# such as ``symox.test``.
NS = {}
exec("from symox import *\n"
     "from symox.physics.quantum import *\n"
     "from symox.physics.quantum.hilbert import *\n"
     "from symox.physics.quantum.qubit import Qubit\n"
     "from symox.physics.quantum.pauli import *\n"
     "x, y = symbols('x y'); f = Function('f')\n"
     "A, B = Operator('A'), Operator('B'); U = UnitaryOperator('U'); k = Ket('k'); b = Bra('b')\n"
     "sx, sy, sz = SigmaX(), SigmaY(), SigmaZ(); sm, sp = SigmaMinus(), SigmaPlus()\n", NS)

# (name, code setting r, expected repr(r) or str(r) of the reference)
PROBES = [
    ('op_label',
     "from symox import *\nfrom symox.physics.quantum import *\nA=Operator('A'); r=A.label",
     '(A,)'),
    ('op_label2',
     "from symox import *\nfrom symox.physics.quantum import *\nr=Operator('A','t').label",
     '(A, t)'),
    ('ket_label',
     "from symox import *\nfrom symox.physics.quantum import *\nr=Ket('psi').label",
     '(psi,)'),
    ('bra_label',
     "from symox import *\nfrom symox.physics.quantum import *\nr=Bra('a','b').label",
     '(a, b)'),
    ('bra_ket_print',
     "from symox import *\nfrom symox.physics.quantum import *\nr=str(Bra('a')*Ket('b'))",
     '<a|b>'),
    ('bra_ket_print2',
     "from symox import *\nfrom symox.physics.quantum import *\nr=str(2*Bra('a')*Ket('b'))",
     '2*<a|b>'),
    ('ip_bra',
     "from symox import *\nfrom symox.physics.quantum import *\nip=InnerProduct(Bra('b'),Ket('k')); r=(ip.bra, ip.ket)",
     '(<b|, |k>)'),
    ('ip_print',
     "from symox import *\nfrom symox.physics.quantum import *\nr=str(InnerProduct(Bra('b'),Ket('k')))",
     '<b|k>'),
    ('comm_doit',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B=Operator('A'),Operator('B'); r=Commutator(A,B).doit()",
     'A*B - B*A'),
    ('comm_coeff',
     "from symox import *\nfrom symox.physics.quantum import *\nx,y=symbols('x y'); A,B=Operator('A'),Operator('B'); r=Commutator(3*x*A,x*y*B)",
     '3*x**2*y*[A,B]'),
    ('comm_expand',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B,C=Operator('A'),Operator('B'),Operator('C'); r=Commutator(A+B,C).expand(commutator=True)",
     '[A,C] + [B,C]'),
    ('comm_neg',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B=Operator('A'),Operator('B'); r=Commutator(B,A)",
     '-[A,B]'),
    ('comm_dagger',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B=Operator('A'),Operator('B'); r=Dagger(Commutator(A,B))",
     '-[Dagger(A),Dagger(B)]'),
    ('acomm_doit',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B=Operator('A'),Operator('B'); r=AntiCommutator(A,B).doit()",
     'A*B + B*A'),
    ('acomm_coeff',
     "from symox import *\nfrom symox.physics.quantum import *\nx,y=symbols('x y'); A,B=Operator('A'),Operator('B'); r=AntiCommutator(3*x*A,x*y*B)",
     '3*x**2*y*{A,B}'),
    ('acomm_dagger',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B=Operator('A'),Operator('B'); r=Dagger(AntiCommutator(A,B))",
     '{Dagger(A),Dagger(B)}'),
    ('comm_self',
     "from symox import *\nfrom symox.physics.quantum import *\nA=Operator('A'); r=Commutator(A,A)",
     '0'),
    ('comm_scalar',
     "from symox import *\nfrom symox.physics.quantum import *\nA=Operator('A'); r=Commutator(A,2)",
     "ERR AttributeError: 'int' object has no attribute 'is_commutative'"),
    ('acomm_scalar',
     "from symox import *\nfrom symox.physics.quantum import *\nA=Operator('A'); r=AntiCommutator(A,2)",
     "ERR AttributeError: 'int' object has no attribute 'is_commutative'"),
    ('comm_spin_doit',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=Commutator(Jx,Jy).doit()',
     'hbar*I*Jz'),
    ('repr_jz',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=represent(Jz, basis=Jz, j=1)',
     'Matrix([\n[hbar, 0,     0],\n[   0, 0,     0],\n[   0, 0, -hbar]])'),
    ('repr_jzket',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=represent(JzKet(1,0), basis=Jz)',
     'Matrix([\n[0],\n[1],\n[0]])'),
    ('repr_jxket',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=represent(JxKet(1,0), basis=Jz)',
     'Matrix([\n[-sqrt(2)/2],\n[         0],\n[ sqrt(2)/2]])'),
    ('repr_jzket_jx',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=represent(JzKet(1,1), basis=Jx)',
     'Matrix([\n[       1/2],\n[-sqrt(2)/2],\n[       1/2]])'),
    ('repr_jx',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=represent(Jx, basis=Jz, j=S(1)/2)',
     'Matrix([\n[     0, hbar/2],\n[hbar/2,      0]])'),
    ('repr_qubit',
     "from symox.physics.quantum.qubit import Qubit\nfrom symox.physics.quantum.represent import represent\nr=represent(Qubit('01'))",
     'Matrix([\n[0],\n[1],\n[0],\n[0]])'),
    ('repr_gate',
     "from symox.physics.quantum.gate import X\nfrom symox.physics.quantum.represent import represent\nfrom symox.physics.quantum.qubit import Qubit\nr=represent(X(0)*Qubit('00'), nqubits=2)",
     'Matrix([\n[0],\n[1],\n[0],\n[0]])'),
    ('repr_up',
     "from symox.physics.quantum.represent import represent\nfrom symox.physics.quantum.state import Ket\nfrom symox.physics.quantum.operator import Operator\nfrom symox import Matrix\nclass SzUpKet(Ket):\n    def _represent_SzOp(self, basis, **options):\n        return Matrix([1,0])\nclass SzOp(Operator):\n    pass\nsz = SzOp('Sz')\nup = SzUpKet('up')\nr=represent(up, basis=sz)",
     'Matrix([\n[1],\n[0]])'),
    ('rep_inner_x',
     'from symox.physics.quantum.represent import rep_innerproduct\nfrom symox.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet\nr=rep_innerproduct(XKet())',
     'DiracDelta(x - x_1)'),
    ('rep_inner_px',
     'from symox.physics.quantum.represent import rep_innerproduct\nfrom symox.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet\nr=rep_innerproduct(XKet(), basis=PxOp())',
     'sqrt(2)*exp(-I*px_1*x/hbar)/(2*sqrt(hbar)*sqrt(pi))'),
    ('rep_inner_px2',
     'from symox.physics.quantum.represent import rep_innerproduct\nfrom symox.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet\nr=rep_innerproduct(PxKet(), basis=XOp())',
     'sqrt(2)*exp(I*px*x_1/hbar)/(2*sqrt(hbar)*sqrt(pi))'),
    ('rep_expect',
     'from symox.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet\nfrom symox.physics.quantum.represent import rep_expectation\nr=rep_expectation(XOp())',
     'x_1*DiracDelta(x_1 - x_2)'),
    ('rep_expect2',
     'from symox.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet\nfrom symox.physics.quantum.represent import rep_expectation\nr=rep_expectation(XOp(), basis=PxOp())',
     '<px_2|*X*|px_1>'),
    ('repr_xket',
     'from symox.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet\nfrom symox.physics.quantum.represent import represent\nr=represent(XKet())',
     'DiracDelta(x - x_1)'),
    ('state_to_ops',
     'from symox.physics.quantum.cartesian import XKet\nfrom symox.physics.quantum.operatorset import state_to_operators\nr=state_to_operators(XKet())',
     'X'),
    ('get_basis',
     'from symox.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet\nfrom symox.physics.quantum.represent import get_basis\nr=(get_basis(XKet()), get_basis(XOp()), get_basis(PxKet()))',
     '(|x>, |x>, |px>)'),
    ('integrate_result',
     "from symox import symbols, DiracDelta\nfrom symox.physics.quantum.represent import integrate_result\nfrom symox.physics.quantum.cartesian import XOp, XKet\nx_ket = XKet()\nX_op = XOp()\nx, x_1, x_2 = symbols('x, x_1, x_2')\nr=integrate_result(X_op*x_ket, x*DiracDelta(x-x_1)*DiracDelta(x_1-x_2))",
     'x*DiracDelta(x - x_1)*DiracDelta(x_1 - x_2)'),
    ('qapply_tp',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=qapply(TensorProduct(Jz,Jz)*TensorProduct(JzKet(1,1),JzKet(1,-1)))',
     '-hbar**2*|1,1>x|1,-1>'),
    ('qapply_tp2',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=qapply(TensorProduct(Jz,Jz)*TensorProduct(JzKet(1,1),JzKet(1,1)))',
     'hbar**2*|1,1>x|1,1>'),
    ('tp_sum_expand',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B,C=Operator('A'),Operator('B'),Operator('C'); r=str(TensorProduct(A+B,C).expand(tensorproduct=True))",
     'AxC + BxC'),
    ('tp_basic',
     "from symox import *\nfrom symox.physics.quantum import *\nA,B=symbols('A B',commutative=False); r=str(TensorProduct(A,B))",
     'AxB'),
    ('couple2',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=couple(TensorProduct(JzKet(S(1)/2,S(1)/2), JzKet(S(1)/2,-S(1)/2)))',
     'sqrt(2)*|0,0,j1=1/2,j2=1/2>/2 + sqrt(2)*|1,0,j1=1/2,j2=1/2>/2'),
    ('couple2b',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=couple(TensorProduct(JzKet(1,1), JzKet(S(1)/2,-S(1)/2)))',
     'sqrt(6)*|1/2,1/2,j1=1,j2=1/2>/3 + sqrt(3)*|3/2,1/2,j1=1,j2=1/2>/3'),
    ('couple3',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=couple(TensorProduct(JzKet(1,1), JzKet(1,1), JzKet(1,0)))',
     'sqrt(6)*|2,2,j1=1,j2=1,j3=1,j(1,2)=2>/3 + sqrt(3)*|3,2,j1=1,j2=1,j3=1,j(1,2)=2>/3'),
    ('couple3b',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=couple(TensorProduct(JzKet(1,1), JzKet(1,1), JzKet(1,0)), ((1,3),(1,2)))',
     'sqrt(2)*|2,2,j1=1,j2=1,j3=1,j(1,3)=1>/2 - sqrt(6)*|2,2,j1=1,j2=1,j3=1,j(1,3)=2>/6 + sqrt(3)*|3,2,j1=1,j2=1,j3=1,j(1,3)=2>/3'),
    ('couple_jx',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=couple(TensorProduct(JxKet(S(1)/2,S(1)/2), JxKet(S(1)/2,-S(1)/2)))',
     'sqrt(2)*|0,0,j1=1/2,j2=1/2>/2 + sqrt(2)*|1,0,j1=1/2,j2=1/2>/2'),
    ('uncouple2',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=uncouple(JzKetCoupled(1, 0, (S(1)/2, S(1)/2)))',
     'sqrt(2)*|1/2,-1/2>x|1/2,1/2>/2 + sqrt(2)*|1/2,1/2>x|1/2,-1/2>/2'),
    ('uncouple_jn',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=uncouple(JzKet(1, 0), (S(1)/2, S(1)/2))',
     'sqrt(2)*|1/2,-1/2>x|1/2,1/2>/2 + sqrt(2)*|1/2,1/2>x|1/2,-1/2>/2'),
    ('uncouple3',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=uncouple(JzKetCoupled(1, 1, (1, 1, 1), ((1,3,1),(1,2,1))))',
     '|1,-1>x|1,1>x|1,1>/2 - |1,0>x|1,0>x|1,1>/2 + |1,1>x|1,0>x|1,0>/2 - |1,1>x|1,1>x|1,-1>/2'),
    ('uncouple3b',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=uncouple(JzKet(1, 1), (1, 1, 1), ((1,3,1),(1,2,1)))',
     '|1,-1>x|1,1>x|1,1>/2 - |1,0>x|1,0>x|1,1>/2 + |1,1>x|1,0>x|1,0>/2 - |1,1>x|1,1>x|1,-1>/2'),
    ('uncouple_jx',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=uncouple(JxKetCoupled(1, 0, (S(1)/2, S(1)/2)))',
     'sqrt(2)*|1/2,-1/2>x|1/2,1/2>/2 + sqrt(2)*|1/2,1/2>x|1/2,-1/2>/2'),
    ('coupled_print',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=JzKetCoupled(1, 0, (S(1)/2, S(1)/2))',
     '|1,0,j1=1/2,j2=1/2>'),
    ('coupled_rep',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=represent(JzKetCoupled(1, 0, (S(1)/2, S(1)/2)))',
     'Matrix([\n[0],\n[0],\n[1],\n[0]])'),
    ('cg_kw',
     'from symox.physics.quantum.cg import CG\nfrom symox import S\nr=CG(j1=S(1)/2, m1=-S(1)/2, j2=S(1)/2, m2=+S(1)/2, j3=1, m3=0).doit()',
     'sqrt(2)/2'),
    ('cg_pos',
     'from symox.physics.quantum.cg import CG\nfrom symox import S\nr=CG(S(3)/2, S(3)/2, S(1)/2, -S(1)/2, 1, 1).doit()',
     'sqrt(3)/2'),
    ('cg_print',
     'from symox.physics.quantum.cg import CG\nfrom symox import S\nr=CG(S(3)/2, S(3)/2, S(1)/2, -S(1)/2, 1, 1)',
     'CG(3/2, 3/2, 1/2, -1/2, 1, 1)'),
    ('cg_attrs',
     'from symox.physics.quantum.cg import CG\nfrom symox import S\nc=CG(S(3)/2, S(3)/2, S(1)/2, -S(1)/2, 1, 1); r=(c.j1,c.m1,c.j2,c.m2,c.j3,c.m3,c.is_symbolic)',
     '(3/2, 3/2, 1/2, -1/2, 1, 1, False)'),
    ('wigner3j',
     'from symox.physics.quantum.cg import Wigner3j\nr=Wigner3j(6,0,4,0,2,0).doit()',
     'sqrt(715)/143'),
    ('wigner6j',
     'from symox.physics.quantum.cg import Wigner6j\nr=Wigner6j(3,3,3,3,3,3).doit()',
     '-1/14'),
    ('wigner9j',
     'from symox.physics.quantum.cg import Wigner9j\nr=Wigner9j(1,1,1,1,1,1,1,1,0).doit()',
     '1/18'),
    ('cg_simp',
     'from symox.physics.quantum.cg import CG, cg_simp\nfrom symox import S\na = CG(1,1,0,0,1,1)\nb = CG(1,0,0,0,1,0)\nc = CG(1,-1,0,0,1,-1)\nr=cg_simp(a+b+c)',
     '3'),
    ('rot_print',
     "from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\na,b,c=symbols('alpha beta gamma'); r=Rotation(a,b,c)",
     'R(alpha,beta,gamma)'),
    ('rot_inverse',
     "from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\na,b,c=symbols('alpha beta gamma'); r=Rotation(a,b,c).inverse()",
     'R(-gamma,-beta,-alpha)'),
    ('rot_D',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=Rotation.D(1, 1, 0, pi, pi/2, -pi)',
     'WignerD(1, 1, 0, pi, pi/2, -pi)'),
    ('rot_d',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=Rotation.d(1, 1, 0, pi/2)',
     'WignerD(1, 1, 0, 0, pi/2, 0)'),
    ('rot_D_doit',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=Rotation.D(1, 1, 0, pi, pi/2, -pi).doit()',
     'sqrt(2)/2'),
    ('rot_d_doit',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=Rotation.d(1, 1, 0, pi/2).doit()',
     '-sqrt(2)/2'),
    ('rot_D_sym',
     "from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\na,b,c=symbols('alpha beta gamma'); r=Rotation.D(1, 1, 0, a, b, c).doit()",
     '-sqrt(2)*exp(-I*alpha)*sin(beta)/2'),
    ('rot_apply_num',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=qapply(Rotation(0,pi/2,0)*JzKet(1,1))',
     '|1,-1>/2 - sqrt(2)*|1,0>/2 + |1,1>/2'),
    ('rot_represent',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=represent(Rotation(0,pi/2,0), j=S(1)/2)',
     'Matrix([\n[ WignerD(1/2, 1/2, 1/2, 0, pi/2, 0),  WignerD(1/2, 1/2, -1/2, 0, pi/2, 0)],\n[WignerD(1/2, -1/2, 1/2, 0, pi/2, 0), WignerD(1/2, -1/2, -1/2, 0, pi/2, 0)]])'),
    ('wignerD_print',
     'from symox import *\nfrom symox.physics.quantum import *\nfrom symox.physics.quantum.spin import *\nr=WignerD(1,1,0,pi,pi/2,-pi)',
     'WignerD(1, 1, 0, pi, pi/2, -pi)'),
    ('fermi_sym',
     "from symox import *\na=Symbol('a', above_fermi=True); r=(a, a.assumptions0, a==Symbol('a'))",
     "(a, {'above_fermi': True, 'commutative': True}, False)"),
    ('fermi_sym2',
     "from symox import *\ni,j=symbols('i j', below_fermi=True); r=(i.assumptions0, i==Symbol('i'))",
     "({'below_fermi': True, 'commutative': True}, False)"),
]


@pytest.mark.parametrize("name,code,expected", PROBES, ids=[p[0] for p in PROBES])
def test_reference_probe(name, code, expected):
    scope = {}
    try:
        exec(code, scope)
        value = scope["r"]
        out = value if isinstance(value, str) else repr(value)
    except Exception as error:
        out = "ERR %s: %s" % (type(error).__name__, str(error)[:150])
    assert out == expected


# Hilbert spaces, DifferentialOperator, hbar, unitary operators, qubits.
EXTRA = [
    ('DifferentialOperator(Derivative(f(x),x),f(x)).variables', '(x,)'),
    ('DifferentialOperator(Derivative(f(x),x),f(x)).expr', 'Derivative(f(x), x)'),
    ('DifferentialOperator(Derivative(f(x),x),f(x)).function', 'f(x)'),
    ('DifferentialOperator(Derivative(f(x,y),x),f(x,y)).free_symbols', '{y, x}'),
    ('(ComplexSpace(2)+FockSpace()).spaces', '(C(2), F)'),
    ('ComplexSpace(2)+FockSpace()', 'C(2)+F'),
    ('ComplexSpace(2)*FockSpace()', 'C(2)*F'),
    ('(ComplexSpace(2)*FockSpace()).spaces', '(C(2), F)'),
    ('FockSpace()', 'F'),
    ('ComplexSpace(2)**3', 'C(2)**3'),
    ('(ComplexSpace(2)*ComplexSpace(3)).dimension', '6'),
    ('(ComplexSpace(2)+ComplexSpace(3)).dimension', '5'),
    ('hbar.evalf()', '1.05457162000000e-34'),
    ('N(hbar,5)', '1.0546e-34'),
    ('U*Dagger(U)', '1'),
    ('Dagger(U)*U', '1'),
    ('qapply(1/(Dagger(k)*k)*Dagger(k)*k*b)', '<b|'),
    ("len(Qubit('0101'))", '4'),
    ("Qubit('0101').dimension", '4'),
    ('Dagger(AntiCommutator(A,B))', '{Dagger(A),Dagger(B)}'),
    ('Dagger(Commutator(A,B))', '-[Dagger(A),Dagger(B)]'),
    ("Operator('A').label", '(A,)'),
    ("Qubit('0101').nqubits", '4'),
]


@pytest.mark.parametrize("code,expected", EXTRA)
def test_reference_repr(code, expected):
    assert repr(eval(code, dict(NS))) == expected


# sympy.physics.quantum.pauli
PAULI = [
    ('sx', 'SigmaX()'),
    ('sy', 'SigmaY()'),
    ('sz', 'SigmaZ()'),
    ("SigmaX('a')", 'SigmaX()'),
    ('sm', 'SigmaMinus()'),
    ('sp', 'SigmaPlus()'),
    ('sx*sy', 'SigmaX()*SigmaY()'),
    ('qsimplify_pauli(sx*sy)', 'I*SigmaZ()'),
    ('qsimplify_pauli(sx*sy*sx)', '-SigmaY()'),
    ("qsimplify_pauli(sx*SigmaY('b'))", 'SigmaX()*SigmaY()'),
    ('sx**2', '1'),
    ('sx**3', 'SigmaX()'),
    ('sm**2', '0'),
    ('sp**2', '0'),
    ('sz**2', '1'),
    ('Dagger(sm)', 'SigmaPlus()'),
    ('Dagger(sx)', 'SigmaX()'),
    ('Dagger(sp)', 'SigmaMinus()'),
    ('Commutator(sx,sy).doit()', '2*I*SigmaZ()'),
    ('Commutator(sy,sz).doit()', '2*I*SigmaX()'),
    ('Commutator(sz,sx).doit()', '2*I*SigmaY()'),
    ('Commutator(sy,sx).doit()', '-2*I*SigmaZ()'),
    ("Commutator(SigmaX('a'),SigmaY('b')).doit()", '0'),
    ('AntiCommutator(sx,sy).doit()', '0'),
    ('AntiCommutator(sx,sx).doit()', '2'),
    ('Commutator(sz,sm).doit()', '-2*SigmaMinus()'),
    ('Commutator(sz,sp).doit()', '2*SigmaPlus()'),
    ('Commutator(sm,sp).doit()', '-SigmaZ()'),
    ('Commutator(sx,sm).doit()', 'SigmaZ()'),
    ('AntiCommutator(sm,sp).doit()', '1'),
    ('AntiCommutator(sz,sm).doit()', '0'),
    ('AntiCommutator(sy,sm).doit()', '-I'),
    ('AntiCommutator(sp,sy).doit()', 'I'),
    ('Commutator(sp,sx).doit()', 'SigmaZ()'),
    ('Commutator(sy,sp).doit()', '-I*SigmaZ()'),
    ('Commutator(sp,sm).doit()', 'SigmaZ()'),
    ('represent(sx)', 'Matrix([[0, 1], [1, 0]])'),
    ('represent(sy)', 'Matrix([[0, -I], [I, 0]])'),
    ('represent(sz)', 'Matrix([[1, 0], [0, -1]])'),
    ('represent(sm)', 'Matrix([[0, 0], [1, 0]])'),
    ('represent(sp)', 'Matrix([[0, 1], [0, 0]])'),
    ('represent(SigmaZKet(0))', 'Matrix([[1], [0]])'),
    ('represent(SigmaZKet(1))', 'Matrix([[0], [1]])'),
    ('represent(SigmaZBra(0))', 'Matrix([[1, 0]])'),
    ('represent(sx*sy)', 'Matrix([[I, 0], [0, -I]])'),
    ('SigmaZKet(0)', '|0>'),
    ('SigmaZBra(1)', '<1|'),
    ('Dagger(SigmaZKet(1))', '<1|'),
    ('qapply(sx*SigmaZKet(0))', '|1>'),
    ('qapply(sy*SigmaZKet(0))', 'I*|1>'),
    ('qapply(sz*SigmaZKet(1))', '-|1>'),
    ('qapply(sm*SigmaZKet(0))', '|1>'),
    ('qapply(sm*SigmaZKet(1))', '0'),
    ('qapply(sp*SigmaZKet(1))', '|0>'),
    ('qapply(sp*SigmaZKet(0))', '0'),
    ('qapply(SigmaZBra(0)*SigmaZKet(0))', '1'),
    ('qapply(SigmaZBra(0)*SigmaZKet(1))', '0'),
    ('qapply(SigmaZBra(0)*sx)', '<0|*SigmaX()'),
    ('qapply(SigmaZBra(1)*sm)', '<1|*SigmaMinus()'),
    ('qapply(sx*sy*SigmaZKet(0))', 'I*|0>'),
    ('sx.name', 'False'),
    ("SigmaX('a').name", 'a'),
    ('sx.use_name', 'False'),
    ("SigmaX('a').use_name", 'True'),
    ('qsimplify_pauli(sm*sp)', '1/2 - SigmaZ()/2'),
    ('qsimplify_pauli(sz*sm)', '-SigmaMinus()'),
    ('qsimplify_pauli(sm*sz)', 'SigmaMinus()'),
    ('qsimplify_pauli(sx*sm)', '1/2 + SigmaZ()/2'),
    ('qsimplify_pauli(2*sx*sy+sz*sz)', '1 + 2*I*SigmaZ()'),
    ('SigmaZKet(0).hilbert_space', 'C(2)'),
    ('sx.hilbert_space', 'H'),
    ("qsimplify_pauli(SigmaX('a')*SigmaY('a'))", 'I*SigmaZ()'),
    ('qsimplify_pauli(sx*sy*sz)', 'I'),
    ('SigmaZKet(1).n', '1'),
    ('SigmaZKet(2)', 'ERR ValueError n must be 0 or 1'),
    ('SigmaZBra(5)', 'ERR ValueError n must be 0 or 1'),
    ('Dagger(sm)==sp', 'True'),
    ("SigmaX('a')==SigmaX('b')", 'False'),
]


@pytest.mark.parametrize("code,expected", PAULI)
def test_pauli_reference(code, expected):
    try:
        out = str(eval(code, dict(NS)))
    except Exception as error:
        out = "ERR %s %s" % (type(error).__name__, error)
    assert out == expected


def test_unsupported_paths_raise_instead_of_wrong_values():
    from symox import Symbol
    from symox.physics.quantum.spin import JzKet, JzKetCoupled, couple, uncouple
    from symox.physics.quantum.tensorproduct import TensorProduct
    from symox.physics.quantum.represent import represent
    from symox.physics.quantum.cartesian import XOp, PxOp, XKet, PxKet
    j1, m1, j2, m2 = Symbol('j1'), Symbol('m1'), Symbol('j2'), Symbol('m2')
    # Symbolic spins give the reference held Sum (see test_c_doc_physics_parity).
    assert str(couple(TensorProduct(JzKet(j1, m1), JzKet(j2, m2)))) == (
        'Sum(CG(j1, m1, j2, m2, j, m1 + m2)*|j,m1 + m2,j1=j1,j2=j2>, (j, m1 + m2, j1 + j2))')
    assert str(uncouple(JzKetCoupled(Symbol('j'), Symbol('m'), (j1, j2)))) == (
        'Sum(CG(j1, m1, j2, m2, j, m)*|j1,m1>x|j2,m2>, (m1, -j1, j1), (m2, -j2, j2))')
    with pytest.raises(NotImplementedError):
        represent(XOp(), basis=PxKet)
    with pytest.raises(NotImplementedError):
        represent(PxOp(), basis=XKet)


def test_pauli_names_commute():
    # SymPy 1.14 returns -2*SigmaMinus() here, ignoring the spin names.
    assert eval("Commutator(SigmaZ('a'), SigmaMinus('b')).doit()", dict(NS)) == 0
    value = eval("AntiCommutator(SigmaZ('a'), SigmaMinus('b')).doit()", dict(NS))
    # Kept as A*B + B*A; SymPy returns 0, the same-spin value.
    assert value == eval("SigmaZ('a')*SigmaMinus('b') + SigmaMinus('b')*SigmaZ('a')", dict(NS))
