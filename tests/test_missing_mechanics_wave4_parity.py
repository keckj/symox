"""Wave 4 parity of sympy.physics.mechanics, vector and continuum_mechanics.

Expected values are the outputs of the reference SymPy 1.14 on the same code.
"""

import pytest

import symox as sp
import symox.physics.continuum_mechanics.beam as beam
import symox.physics.mechanics as me
from symox.physics.continuum_mechanics.cable import Cable

Beam, Beam3D = beam.Beam, beam.Beam3D
symbols, Symbol, SF = sp.symbols, sp.Symbol, sp.SingularityFunction

CASES = {}


def case(f):
    CASES[f.__name__] = f
    return f


@case
def beam_rotation_hinge():
    E, I = symbols('E I')
    b = Beam(15, E, I)
    r0, m0 = b.apply_support(0, type='fixed')
    r10 = b.apply_support(10, type='pin')
    r15 = b.apply_support(15, type='pin')
    p7 = b.apply_rotation_hinge(7)
    b.apply_load(-10, 7, -1)
    b.apply_load(-2, 10, 0, 15)
    b.solve_for_reaction_loads(r0, m0, r10, r15)
    return [b.reaction_loads, b.rotation_jumps, b.bending_moment(), b.deflection()]

@case
def beam_sliding_hinge():
    E, I = symbols('E I')
    b = Beam(13, 20, 20)
    r0, m0 = b.apply_support(0, type='fixed')
    s8 = b.apply_sliding_hinge(8)
    r13 = b.apply_support(13, type='pin')
    b.apply_load(-10, 5, -1)
    b.solve_for_reaction_loads(r0, m0, r13)
    return [b.reaction_loads, b.deflection_jumps, b.shear_force(), b.deflection()]

@case
def beam_ild():
    E, I = symbols('E I')
    l = Symbol('l', positive=True)
    R_0, R_8 = symbols('R_0 R_8')
    b = Beam(8, E, I)
    p0 = b.apply_support(0, 'roller')
    p8 = b.apply_support(8, 'roller')
    b.solve_for_ild_reactions(1, R_0, R_8)
    out = [b.ild_reactions]
    b.solve_for_ild_shear(4, 1, R_0, R_8)
    out.append(b.ild_shear)
    b.solve_for_ild_moment(4, 1, R_0, R_8)
    out.append(b.ild_moment)
    return out

@case
def beam_ild_overhang():
    E, I = symbols('E I')
    R_0, R_7 = symbols('R_0 R_7')
    b = Beam(12, E, I)
    p0 = b.apply_support(0, 'pin')
    p7 = b.apply_support(7, 'roller')
    b.apply_load(-3, 9, -1)
    b.solve_for_ild_reactions(1, R_0, R_7)
    out = [b.ild_reactions]
    b.solve_for_ild_shear(5, 1, R_0, R_7)
    out.append(b.ild_shear)
    b.solve_for_ild_moment(5, 1, R_0, R_7)
    out.append(b.ild_moment)
    return out

@case
def beam_max_values_cantilever():
    b = Beam(6, 30, 10)
    r0, m0 = b.apply_support(0, 'fixed')
    b.apply_load(-4, 3, -1)
    b.solve_for_reaction_loads(r0, m0)
    return [b.max_shear_force(), b.max_bmoment(), b.max_deflection()]

@case
def beam_join_and_cflexure():
    E, I = symbols('E I')
    b1 = Beam(2, E, I)
    b2 = Beam(4, E, I)
    b = b1.join(b2, 'fixed')
    b.apply_load(20, 4, -1)
    R1, M1 = symbols('R1 M1')
    b.apply_load(R1, 0, -1)
    b.apply_load(M1, 0, -2)
    b.bc_slope = [(0, 0)]
    b.bc_deflection = [(0, 0)]
    b.solve_for_reaction_loads(R1, M1)
    return [b.load, b.slope(), b.deflection()]

@case
def beam3d_torsion():
    E, G, I, A, x = symbols('E G I A x')
    b = Beam3D(20, E, G, I, A, x)
    b.apply_moment_load(4, 4, -2, dir='x')
    b.apply_moment_load(4, 8, -2, dir='x')
    b.apply_moment_load(4, 8, -2, dir='x')
    b.solve_for_torsion()
    return [b.torsional_moment(), b.angular_deflection(), b.angular_deflection().subs(x, 3),
            b.polar_moment()]

@case
def beam3d_max():
    x = Symbol('x')
    b = Beam3D(20, 40, 21, 100, 25, x)
    b.apply_load(15, start=0, order=0, dir="z")
    b.apply_load(12*x, start=0, order=0, dir="y")
    b.bc_deflection = [(0, [0, 0, 0]), (20, [0, 0, 0])]
    R1, R2, R3, R4 = symbols('R1, R2, R3, R4')
    b.apply_load(R1, start=0, order=-1, dir="z")
    b.apply_load(R2, start=20, order=-1, dir="z")
    b.apply_load(R3, start=0, order=-1, dir="y")
    b.apply_load(R4, start=20, order=-1, dir="y")
    b.solve_for_reaction_loads(R1, R2, R3, R4)
    return [b.max_shear_force(), b.max_bending_moment(), b.max_bmoment(), b.max_deflection()]

@case
def beam3d_reaction_loads():
    l, E, G, I, A = 20, 40, 21, 100, 25
    x = Symbol('x')
    R1, R2, R3, R4 = symbols('R1 R2 R3 R4')
    b = Beam3D(l, E, G, I, A, x)
    b.apply_load(R1, 0, -1, dir='y')
    b.apply_load(R2, 20, -1, dir='y')
    b.apply_load(R3, 0, -1, dir='z')
    b.apply_load(R4, 20, -1, dir='z')
    b.apply_load(-6, 5, 0, dir='y')
    b.apply_load(4, 2, 0, dir='z')
    b.bc_deflection = [(0, [0, 0, 0]), (20, [0, 0, 0])]
    b.solve_for_reaction_loads(R1, R2, R3, R4)
    return [b.reaction_loads, b.shear_force(), b.bending_moment(), b.deflection()]

@case
def sf_limits():
    x, R = symbols('x R')
    lp = Symbol('lp', positive=True)
    exprs = [SF(x, 0, 0), SF(x, 4, 0), SF(x, 2, 1), R*SF(x, 4, 0), SF(x, 4, -1),
             SF(x, 2, 1) + SF(x, 4, 0), SF(x, 4, 2)/(x - 3)]
    out = [(sp.limit(e, x, 4, '+'), sp.limit(e, x, 4, '-')) for e in exprs]
    out.append(sp.limit(SF(x, lp, 1), x, 2*lp, '+'))
    out.append(sp.limit(SF(x, 0, 1), x, lp, '-'))
    return out

@case
def sf_lambdify():
    x = Symbol('x')
    e = 3*SF(x, 2, 1) - SF(x, 4, 0) + SF(x, 1, -1) + SF(x, 1, 2)/2
    out = []
    for mod in ('numpy', 'math', 'mpmath'):
        f = sp.lambdify(x, e, mod)
        out.append([float(f(v)) for v in (0, 1, 2, 3, 4, 5.5)])
    import numpy as np
    f = sp.lambdify(x, [e, SF(x, 3, 0)], 'numpy')
    out.append([list(map(float, a)) for a in f(np.array([0., 1., 3., 4.5]))])
    return out

@case
def coordinate_sym():
    A = me.ReferenceFrame('A')
    a_y = me.CoordinateSym('a_y', A, 1)
    B = me.ReferenceFrame('A')
    return [str(A[0]), A[0].is_real, A[0].frame is A, a_y == A[1], A[0] == B[0],
            str(A[0] + B[0]), (A[0]*B[0]).diff(A[0]) == B[0], str(A.varlist)]

@case
def frame_variables_express():
    N = me.ReferenceFrame('N')
    q = me.dynamicsymbols('q')
    A = N.orientnew('A', 'Axis', [q, N.z])
    return [me.express(A[0]*A.x + A[1]*A.y, N, variables=True),
            A.variable_map(N), me.gradient(N[0]**2*N[1], N), me.curl(N[1]*N.x, N),
            me.scalar_potential(2*N[0]*N.x + N.z, N), me.is_conservative(N[1]*N.x)]

@case
def pin_joint_axes():
    q, u = me.dynamicsymbols('q u')
    l = Symbol('l')
    P = me.RigidBody('P')
    C = me.RigidBody('C')
    J = me.PinJoint('J', P, C, q, u, parent_point=l*P.frame.x, child_point=-l*C.frame.y,
                    parent_axis=P.frame.z, child_axis=-C.frame.x,
                    parent_interframe=P.frame, joint_axis=P.frame.z)
    return [J.joint_axis, C.frame.dcm(P.frame), C.masscenter.pos_from(P.masscenter),
            C.masscenter.vel(P.frame), J.kdes, C.frame.ang_vel_in(P.frame)]

@case
def pin_joint_doc():
    l1, l2 = symbols('l1 l2')
    q1, q2 = me.dynamicsymbols('q1 q2')
    ceiling = me.RigidBody('C')
    upper_bob = me.RigidBody('U')
    lower_bob = me.RigidBody('L')
    ceiling_joint = me.PinJoint('P1', ceiling, upper_bob, coordinates=q1,
                                child_point=-l1*upper_bob.frame.x, joint_axis=ceiling.frame.z)
    pendulum_joint = me.PinJoint('P2', upper_bob, lower_bob, coordinates=q2,
                                 child_point=-l2*lower_bob.frame.x, joint_axis=upper_bob.frame.z)
    return [upper_bob.frame.dcm(ceiling.frame), lower_bob.frame.dcm(ceiling.frame).simplify(),
            lower_bob.masscenter.pos_from(ceiling.masscenter).express(ceiling.frame).simplify(),
            upper_bob.masscenter.vel(ceiling.frame), lower_bob.masscenter.vel(ceiling.frame)]

@case
def other_joints():
    P = me.RigidBody('P'); C = me.RigidBody('C')
    J1 = me.PrismaticJoint('S', P, C, parent_point=P.frame.y, child_axis=C.frame.y,
                           parent_axis=P.frame.x)
    P2 = me.RigidBody('P2'); C2 = me.RigidBody('C2')
    J2 = me.CylindricalJoint('CJ', P2, C2, joint_axis=P2.frame.z, parent_point=P2.frame.x)
    P3 = me.RigidBody('P3'); C3 = me.RigidBody('C3')
    J3 = me.PlanarJoint('PJ', P3, C3, child_point=C3.frame.y)
    P4 = me.RigidBody('P4'); C4 = me.RigidBody('C4')
    J4 = me.SphericalJoint('SJ', P4, C4, child_point=C4.frame.y, rot_type='body', rot_order=123)
    P5 = me.RigidBody('P5'); C5 = me.RigidBody('C5')
    J5 = me.WeldJoint('WJ', P5, C5, parent_point=P5.frame.x, child_interframe=C5.frame.y)
    return [J1.kdes, J1.child.frame.dcm(J1.parent.frame), J1.child_point.pos_from(J1.parent_point),
            C.masscenter.vel(P.frame),
            J2.coordinates, J2.speeds, C2.masscenter.vel(P2.frame), C2.frame.ang_vel_in(P2.frame),
            J3.rotation_axis, J3.planar_vectors, C3.masscenter.vel(P3.frame),
            C4.frame.dcm(P4.frame), C4.frame.ang_vel_in(P4.frame), C4.masscenter.vel(P4.frame),
            C5.frame.dcm(P5.frame), C5.masscenter.pos_from(P5.masscenter), J5.coordinates.shape]

@case
def system_form_eoms():
    q, u = me.dynamicsymbols('q u')
    m, k, c, g = symbols('m k c g')
    wall = me.RigidBody('W')
    system = me.System.from_newtonian(wall)
    bob = me.Particle('P', mass=m)
    J = me.PrismaticJoint('J', wall, bob, q, u, joint_axis=system.x)
    system.add_joints(J)
    bob.potential_energy = m*g*system.y.dot(bob.masscenter.pos_from(system.fixed_point))
    system.add_loads(me.Force(bob, -k*q*system.x - c*u*system.x))
    system.validate_system()
    eom = system.form_eoms()
    return [eom, system.mass_matrix_full, system.forcing_full, system.rhs(),
            system.q, system.u, system.kdes]

@case
def system_doc_lagrange():
    l1, l2, m1, m2, g = symbols('l1 l2 m1 m2 g')
    q1, q2 = me.dynamicsymbols('q1 q2')
    ceiling = me.RigidBody('C')
    upper = me.Particle('U', mass=m1)
    lower = me.Particle('L', mass=m2)
    sys_ = me.System.from_newtonian(ceiling)
    upper_frame = me.ReferenceFrame('UF'); lower_frame = me.ReferenceFrame('LF')
    j1 = me.PinJoint('J1', ceiling, upper, q1, child_point=-l1*upper_frame.x,
                     joint_axis=ceiling.z, child_interframe=upper_frame)
    j2 = me.PinJoint('J2', upper, lower, q2, child_point=-l2*lower_frame.x,
                     joint_axis=ceiling.z, parent_interframe=upper_frame, child_interframe=lower_frame)
    sys_.add_joints(j1, j2)
    sys_.add_loads(me.Force(upper, -m1*g*ceiling.y), me.Force(lower, -m2*g*ceiling.y))
    sys_.validate_system()
    eom = sys_.form_eoms(eom_method=me.LagrangesMethod)
    return [eom.applyfunc(sp.simplify), sys_.mass_matrix.applyfunc(sp.simplify)]

@case
def joints_method():
    q, u = me.dynamicsymbols('q u')
    m, k, b = symbols('m k b')
    wall = me.RigidBody('W')
    part = me.RigidBody('P', mass=m)
    part.potential_energy = k * q**2 / 2
    J = me.PrismaticJoint('J', wall, part, coordinates=q, speeds=u)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        method = me.JointsMethod(wall, J)
        part.apply_force(-b*u*wall.frame.x, reaction_body=wall) if hasattr(part, 'apply_force') else None
        method = me.JointsMethod(wall, J)
        out = [method.form_eoms(), method.mass_matrix, method.forcing, method.rhs()]
    lm = me.JointsMethod(wall, J)
    out.append(lm.form_eoms(me.LagrangesMethod))
    return out

@case
def empty_matrix_stacking():
    a, b = symbols('a b')
    M = sp.Matrix
    return [M([M([a]), M([])]).shape, M([M([]), M([a, b])]).shape, M([M([]), M([])]).shape,
            M([sp.zeros(0, 3), M([[a, b, a]])]).shape]


EXPECTED = {
    'beam3d_max':
        '[[(0, 0), (20, 2400), (20, 300)], [(0, 0), (20, 3000), (20, 16000)], [(0, 0), (20, 3000), (20, 16000)], [(0, 0), (0, 0), (0, 0)]]',
    'beam3d_reaction_loads':
        '[{R1: 135/4, R2: 225/4, R3: -162/5, R4: -198/5}, [0, 6*x, -4*x], [0, -2*x**2, -3*x**2], [0, 0, 0]]',
    'beam3d_torsion':
        '[0, Piecewise((12*x, x <= 4), (8*x + 16, x <= 8), (80, True))/(2*G*I), 18/(G*I), 2*I]',
    'beam_ild':
        '[{R_0: -SingularityFunction(a, 0, 0) + SingularityFunction(a, 0, 1)/8 - SingularityFunction(a, 8, 1)/8, R_8: -SingularityFunction(a, 0, 1)/8 + SingularityFunction(a, 8, 0) + SingularityFunction(a, 8, 1)/8}, -(-SingularityFunction(a, 0, 0) + SingularityFunction(a, 8, 0) + 2)*SingularityFunction(a, 4, 0) - SingularityFunction(-a, 0, 0) - SingularityFunction(a, 0, 0) + SingularityFunction(a, 0, 1)/8 + SingularityFunction(a, 8, 0) - SingularityFunction(a, 8, 1)/8 + 1, -(4*SingularityFunction(a, 0, 0) - SingularityFunction(a, 0, 1) + SingularityFunction(a, 4, 1))*SingularityFunction(a, 4, 0) - SingularityFunction(a, 0, 1)/2 + SingularityFunction(a, 4, 1) - SingularityFunction(a, 8, 1)/2]',
    'beam_ild_overhang':
        '[{R_0: -SingularityFunction(a, 0, 0) + SingularityFunction(a, 0, 1)/7 - 5*SingularityFunction(a, 12, 0)/7 - SingularityFunction(a, 12, 1)/7 - 6/7, R_7: -SingularityFunction(a, 0, 1)/7 + 12*SingularityFunction(a, 12, 0)/7 + SingularityFunction(a, 12, 1)/7 + 27/7}, -(-SingularityFunction(a, 0, 0) + SingularityFunction(a, 12, 0) + 2)*SingularityFunction(a, 5, 0) - SingularityFunction(-a, 0, 0) - SingularityFunction(a, 0, 0) + SingularityFunction(a, 0, 1)/7 + 2*SingularityFunction(a, 12, 0)/7 - SingularityFunction(a, 12, 1)/7 + 1/7, -(5*SingularityFunction(a, 0, 0) - SingularityFunction(a, 0, 1) + SingularityFunction(a, 5, 1))*SingularityFunction(a, 5, 0) - 2*SingularityFunction(a, 0, 1)/7 + SingularityFunction(a, 5, 1) - 25*SingularityFunction(a, 12, 0)/7 - 5*SingularityFunction(a, 12, 1)/7 - 30/7]',
    'beam_join_and_cflexure':
        '[80*SingularityFunction(x, 0, -2) - 20*SingularityFunction(x, 0, -1) + 20*SingularityFunction(x, 4, -1), -(-80*SingularityFunction(x, 0, 1) + 10*SingularityFunction(x, 0, 2) - 10*SingularityFunction(x, 4, 2))/(E*I), -(-40*SingularityFunction(x, 0, 2) + 10*SingularityFunction(x, 0, 3)/3 - 10*SingularityFunction(x, 4, 3)/3)/(E*I)]',
    'beam_max_values_cantilever':
        '[(Interval(0, 3), 4), (0, 12), None]',
    'beam_rotation_hinge':
        '[{R_0: 501/332, M_0: -3507/332, R_10: 7713/415, R_15: -157/1660}, {P_7: 53395/(498*E*I)}, 3507*SingularityFunction(x, 0, 0)/332 - 501*SingularityFunction(x, 0, 1)/332 - 53395*SingularityFunction(x, 7, -1)/498 + 10*SingularityFunction(x, 7, 1) - 7713*SingularityFunction(x, 10, 1)/415 + SingularityFunction(x, 10, 2) + 157*SingularityFunction(x, 15, 1)/1660 - SingularityFunction(x, 15, 2), -(3507*SingularityFunction(x, 0, 2)/664 - 167*SingularityFunction(x, 0, 3)/664 - 53395*SingularityFunction(x, 7, 1)/498 + 5*SingularityFunction(x, 7, 3)/3 - 2571*SingularityFunction(x, 10, 3)/830 + SingularityFunction(x, 10, 4)/12 + 157*SingularityFunction(x, 15, 3)/9960 - SingularityFunction(x, 15, 4)/12)/(E*I)]',
    'beam_sliding_hinge':
        '[{R_0: 10, M_0: -50, R_13: 0}, {W_8: 85/24}, 50*SingularityFunction(x, 0, -1) - 10*SingularityFunction(x, 0, 0) + 10*SingularityFunction(x, 5, 0) - 4250*SingularityFunction(x, 8, -3)/3, -SingularityFunction(x, 0, 2)/16 + SingularityFunction(x, 0, 3)/240 - SingularityFunction(x, 5, 3)/240 + 85*SingularityFunction(x, 8, 0)/24]',
    'coordinate_sym':
        "['A_x', None, True, True, False, 'A_x + A_x', True, '(A_x, A_y, A_z)']",
    'empty_matrix_stacking':
        '[(1, 1), (2, 1), (0, 0), (1, 3)]',
    'frame_variables_express':
        '[(-(-N_x*sin(q(t)) + N_y*cos(q(t)))*sin(q(t)) + (N_x*cos(q(t)) + N_y*sin(q(t)))*cos(q(t)))*N.x + ((-N_x*sin(q(t)) + N_y*cos(q(t)))*cos(q(t)) + (N_x*cos(q(t)) + N_y*sin(q(t)))*sin(q(t)))*N.y, {A_x: N_x*cos(q(t)) + N_y*sin(q(t)), A_y: -N_x*sin(q(t)) + N_y*cos(q(t)), A_z: N_z}, 2*N_x*N_y*N.x + N_x**2*N.y, - N.z, N_x**2 + N_z, False]',
    'joints_method':
        '[Matrix([[-m*Derivative(u(t), t)]]), Matrix([[m]]), Matrix([[0]]), Matrix([\n[u(t)],\n[   0]]), Matrix([[k*q(t)]])]',
    'other_joints':
        '[Matrix([[u_S(t) - Derivative(q_S(t), t)]]), Matrix([\n[0, -1, 0],\n[1,  0, 0],\n[0,  0, 1]]), q_S(t)*P_frame.x, u_S(t)*P_frame.x, Matrix([\n[q0_CJ(t)],\n[q1_CJ(t)]]), Matrix([\n[u0_CJ(t)],\n[u1_CJ(t)]]), u1_CJ(t)*P2_frame.z, u0_CJ(t)*P2_frame.z, P3_frame.x, [P3_frame.y, P3_frame.z], u1_PJ(t)*P3_frame.y + u2_PJ(t)*P3_frame.z - u0_PJ(t)*C3_frame.z, Matrix([\n[ cos(q1_SJ(t))*cos(q2_SJ(t)),  sin(q0_SJ(t))*sin(q1_SJ(t))*cos(q2_SJ(t)) + sin(q2_SJ(t))*cos(q0_SJ(t)), sin(q0_SJ(t))*sin(q2_SJ(t)) - sin(q1_SJ(t))*cos(q0_SJ(t))*cos(q2_SJ(t))],\n[-sin(q2_SJ(t))*cos(q1_SJ(t)), -sin(q0_SJ(t))*sin(q1_SJ(t))*sin(q2_SJ(t)) + cos(q0_SJ(t))*cos(q2_SJ(t)), sin(q0_SJ(t))*cos(q2_SJ(t)) + sin(q1_SJ(t))*sin(q2_SJ(t))*cos(q0_SJ(t))],\n[               sin(q1_SJ(t)),                                             -sin(q0_SJ(t))*cos(q1_SJ(t)),                                             cos(q0_SJ(t))*cos(q1_SJ(t))]]), (u0_SJ(t)*cos(q1_SJ(t))*cos(q2_SJ(t)) + u1_SJ(t)*sin(q2_SJ(t)))*C4_frame.x + (-u0_SJ(t)*sin(q2_SJ(t))*cos(q1_SJ(t)) + u1_SJ(t)*cos(q2_SJ(t)))*C4_frame.y + (u0_SJ(t)*sin(q1_SJ(t)) + u2_SJ(t))*C4_frame.z, (u0_SJ(t)*sin(q1_SJ(t)) + u2_SJ(t))*C4_frame.x + (-u0_SJ(t)*cos(q1_SJ(t))*cos(q2_SJ(t)) - u1_SJ(t)*sin(q2_SJ(t)))*C4_frame.z, Matrix([\n[0, -1, 0],\n[1,  0, 0],\n[0,  0, 1]]), P5_frame.x, (0, 0)]',
    'pin_joint_axes':
        '[P_frame.z, Matrix([\n[-cos(q(t)), -sin(q(t)), 0],\n[ sin(q(t)), -cos(q(t)), 0],\n[         0,          0, 1]]), l*P_frame.x + l*C_frame.y, - l*u(t)*C_frame.x, Matrix([[u(t) - Derivative(q(t), t)]]), u(t)*P_frame.z]',
    'pin_joint_doc':
        '[Matrix([\n[ cos(q1(t)), sin(q1(t)), 0],\n[-sin(q1(t)), cos(q1(t)), 0],\n[          0,          0, 1]]), Matrix([\n[ cos(q1(t) + q2(t)), sin(q1(t) + q2(t)), 0],\n[-sin(q1(t) + q2(t)), cos(q1(t) + q2(t)), 0],\n[                  0,                  0, 1]]), (l1*cos(q1(t)) + l2*cos(q1(t) + q2(t)))*C_frame.x + (l1*sin(q1(t)) + l2*sin(q1(t) + q2(t)))*C_frame.y, l1*u_P1(t)*U_frame.y, l1*u_P1(t)*U_frame.y + l2*(u_P1(t) + u_P2(t))*L_frame.y]',
    'sf_lambdify':
        '[[0.0, inf, 0.5, 5.0, 9.5, 19.625], [0.0, inf, 0.5, 5.0, 9.5, 19.625], [0.0, inf, 0.5, 5.0, 9.5, 19.625], [[0.0, inf, 5.0, 12.625], [0.0, 0.0, 1.0, 1.0]]]',
    'sf_limits':
        '[(1, 1), (1, 0), (2, 2), (R, 0), (0, 0), (3, 2), (0, 0), lp, lp]',
    'system_doc_lagrange':
        '[Matrix([\n[                                              0],\n[l1*l2*m2*(u_J1(t) + u_J2(t))*u_J1(t)*sin(q2(t))]]), Matrix([\n[0, 0],\n[0, 0]])]',
    'system_form_eoms':
        '[Matrix([[-c*u(t) - k*q(t) - m*Derivative(u(t), t)]]), Matrix([\n[-1, 0],\n[ 0, m]]), Matrix([\n[           -u(t)],\n[-c*u(t) - k*q(t)]]), Matrix([\n[                u(t)],\n[(-c*u(t) - k*q(t))/m]]), Matrix([[q(t)]]), Matrix([[u(t)]]), Matrix([[u(t) - Derivative(q(t), t)]])]',
}


@pytest.mark.parametrize('name', sorted(EXPECTED))
def test_reference_case(name):
    pytest.importorskip('numpy')
    pytest.importorskip('mpmath')
    assert str(CASES[name]()) == EXPECTED[name]


def _extra_outputs():
    from importlib import import_module
    L = 'symox'
    x, y = sp.symbols('x y'); SF = sp.SingularityFunction
    P = sp.Piecewise
    out = {}
    out['solve_domain']=str(sp.solve(x**2-4,x,domain=sp.S.Reals))
    out['pw1']=str(sp.solve(P((x-1,x<0),(x-3,True)),x))
    out['pw2']=str(sp.solve(P((x**2-1,x<2),(x-5,True)),x))
    out['pw3']=str(sp.solve(sp.Eq(P((2*x,x<3),(x+2,True)),5),x))
    out['pw_param']=str(sp.solve(P((x-y,x<2),(x-3,True)),x))
    out['pw_dict']=str(sp.solve(P((x-1,x<0),(x-3,True)),x,dict=True))
    out['sf_int']=str(sp.integrate(SF(x,0,1),(x,2,x)))
    out['sf_int2']=str(sp.integrate(SF(x,1,0),(x,0,3)))
    out['sf_rewrite']=str(P((SF(x,1,1),x<3),(0,True)).rewrite(P))
    tr=import_module(L+'.physics.continuum_mechanics.truss')
    t=tr.Truss(); t.add_node(('node_1',0,0),('node_2',6,0),('node_3',2,2),('node_4',2,0))
    t.add_member(('member_1','node_1','node_4'),('member_2','node_2','node_4'),('member_3','node_1','node_3'))
    t.add_member(('member_4','node_2','node_3'),('member_5','node_3','node_4'))
    t.apply_load(('node_4',10,270)); t.apply_support(('node_1','pinned'),('node_2','roller'))
    t.solve(); out['truss']=str([t.reaction_loads, t.internal_forces])
    ar=import_module(L+'.physics.continuum_mechanics.arch')
    a=ar.Arch((0,0),(10,0),crown_x=5,crown_y=5)
    a.apply_load(-1,'C',8,150,angle=270); a.solve()
    out['arch']=str([a.reaction_force, a.bending_moment_at(4), a.axial_force_at(4), a.shear_force_at(4)])
    return out


EXTRA_EXPECTED = {
    'solve_domain':
        '[-2, 2]',
    'pw1':
        '[3]',
    'pw2':
        '[-1, 1, 5]',
    'pw3':
        '[5/2, 3]',
    'pw_param':
        '[3, Piecewise((y, y < 2), (nan, True))]',
    'pw_dict':
        '[{x: 3}]',
    'sf_int':
        'SingularityFunction(x, 0, 2)/2 - 2',
    'sf_int2':
        '2',
    'sf_rewrite':
        'Piecewise((Piecewise((x - 1, x >= 1), (0, True)), x < 3), (0, True))',
    'truss':
        "[{'R_node_1_x': 0, 'R_node_1_y': 20/3, 'R_node_2_y': 10/3}, {'member_1': 20/3, 'member_2': 20/3, 'member_3': -20*sqrt(2)/3, 'member_4': -10*sqrt(5)/3, 'member_5': 10}]",
    'arch':
        '[{R_A_x: 30, R_A_y: 30, R_B_x: -30, R_B_y: 120}, -24, 210*sqrt(29)/29, 90*sqrt(29)/29]',
}


def test_solve_flags_piecewise_integrate_truss_arch():
    out = _extra_outputs()
    for key, value in EXTRA_EXPECTED.items():
        assert out[key] == value, key


def test_plot_degenerate_range_and_decorations(tmp_path):
    pytest.importorskip('matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    x = Symbol('x')
    p = sp.plot(x, (x, 1, 1), show=False)
    assert list(p[0].get_data()[0][:3]) == [1.0, 1.0, 1.0]
    p = sp.plot(x, (x, 0, 1), show=False, axis=False,
                rectangles=[{'xy': (sp.Float(0.5), 0), 'width': sp.Rational(1, 4), 'height': 1}],
                annotations=[{'text': 'a', 'xy': (0, sp.Integer(1)), 'xytext': (1, 1)}],
                markers=[{'args': [[0, 1], [0, 1]], 'marker': 'o'}])
    p.save(str(tmp_path / 'p.png'))
    assert len(p._backend.ax.patches) == 1 and len(p._backend.ax.texts) == 1


def test_cable_change_support_and_draw():
    c = Cable(('A', 0, 10), ('B', 10, 10))
    c.change_support('B', ('C', 5, 6))
    assert str(c.supports) == "{'A': [0, 10], 'C': [5, 6]}"
    c = Cable(('A', 0, 10), ('B', 10, 10))
    c.apply_load(-1, ('Z', 2, 7.26, 3, 270))
    c.apply_load(-1, ('X', 4, 6, 8, 270))
    with pytest.raises(ValueError):
        c.change_support('B', ('C', 3, 6))
    p = c.draw()
    assert str(p[0].expr) == ('Piecewise((10 - 1.37*x, x <= 2), (8.52 - 0.63*x, x <= 4), '
                              '(2*x/3 + 10/3, x <= 10))')
    assert str(p.options['rectangles']) == (
        "[{'xy': (-0.750000000000000, 10), 'width': 0.750000000000000, 'height': 0.750000000000000, "
        "'color': 'brown', 'fill': False}, {'xy': (10, 10), 'width': 0.750000000000000, "
        "'height': 0.750000000000000, 'color': 'brown', 'fill': False}]")
    assert str(p.options['annotations'][1]) == "{'text': '3N', 'xy': (2, 5.66000000000000)}"
    c = Cable(('A', 0, 40), ('B', 100, 20))
    c.apply_load(0, ('X', 850))
    c.solve(58.58)
    p = c.draw()
    assert str(p[0].expr) == '0.0116550116550117*(x - 58.58)**2 + 0.00447086247086247'
    assert abs(float(p[1].expr) + 7.49552913752915) < 1e-12
    assert str(p.options['annotations'][-1]) == "{'text': '850 N/m', 'xy': (50, -14.9955291375291)}"


def test_basic_is_same():
    x = Symbol('x')
    assert sp.Basic.is_same(sp.eye(1), sp.sympify(sp.eye(1))) is False
    assert sp.Basic.is_same(sp.eye(1), sp.eye(1)) is True
    assert sp.Basic.is_same(sp.Float(0.5), sp.Rational(1, 2)) is False
    assert sp.Basic.is_same(x + 1, x + 1) is True


def test_factor_together_apart_of_sequences():
    x = Symbol('x')
    assert str(sp.factor([x**2 - 1, 0])) == '[(x - 1)*(x + 1), 0]'
    assert str(sp.factor((x**2 - 1,))) == '((x - 1)*(x + 1),)'
    assert str(sp.factor([[x**2 - 1]])) == '[[(x - 1)*(x + 1)]]'
    assert str(sp.together((1/x + 1, 2))) == '((x + 1)/x, 2)'
    assert str(sp.apart([1/(x**2 - 1)])) == '[-1/(2*(x + 1)) + 1/(2*(x - 1))]'


def test_beam3d_doc_slope_factor():
    l, E, G, I, A = symbols('l, E, G, I, A')
    b = Beam3D(l, E, G, I, A)
    x, q, m = symbols('x, q, m')
    b.apply_load(q, 0, 0, dir="y")
    b.apply_moment_load(m, 0, -1, dir="z")
    b.bc_slope = [(0, [0, 0, 0]), (l, [0, 0, 0])]
    b.bc_deflection = [(0, [0, 0, 0]), (l, [0, 0, 0])]
    b.solve_slope_deflection()
    assert str(sp.factor(b.slope())) == (
        '[0, 0, x*(-l + x)*(-A*G*l**3*q + 2*A*G*l**2*q*x - 12*E*I*l*q - 72*E*I*m + 24*E*I*q*x)'
        '/(12*E*I*(A*G*l**2 + 12*E*I))]')


VLATEX_EXPECTED = ['q_{1} + q_{2}', 'q_{1} + q_{2}', 'q_{2} + \\dot{q}_{1}', 'q_{2} + \\dot{q}_{1}', 'q_{1}^{2} + \\sin{\\left(q_{2} \\right)} \\ddot{q}_{1}', 'q_{1}^{2} + \\sin{\\left(q_{2} \\right)} \\ddot{q}_{1}', '\\left[\\begin{matrix}q_{1}\\\\\\dot{q}_{2}\\end{matrix}\\right]', '\\left[\\begin{matrix}q_{1}\\\\\\dot{q}_{2}\\end{matrix}\\right]', 'q_{1}\\mathbf{\\hat{n}_x} + \\dot{q}_{2}\\mathbf{\\hat{n}_y}', 'q_{1}\\mathbf{\\hat{n}_x} + \\dot{q}_{2}\\mathbf{\\hat{n}_y}', 'f{\\left(t,q_{1} \\right)}', 'f{\\left(t,q_{1} \\right)}', '\\dot{q}_{1}^{2}', '\\dot{q}_{1}^{2}', '\\sin^{2}{\\left(q_{1} \\right)}', '\\sin^{2}{\\left(q_{1} \\right)}', '\\left(q_{2} + \\dot{q}_{1}\\right)^{3}', '\\left(q_{2} + \\dot{q}_{1}\\right)^{3}']
VLATEX_INLINE_EXPECTED = ['x \\dddot{\\theta} + u^{3}', '$x \\dddot{\\theta} + u^{3}$', '\\cos{\\left(\\theta \\right)} \\dot{\\theta}^{2}', '$\\cos{\\left(\\theta \\right)} \\dot{\\theta}^{2}$', 'a + \\dot{\\alpha}', '$a + \\dot{\\alpha}$', '\\ddddot{\\theta}', '$\\ddddot{\\theta}$', '\\frac{2 \\ddot{u}}{x}', '$2 \\ddot{u} / x$', '\\left[\\begin{matrix}\\theta & u\\\\a & \\alpha\\end{matrix}\\right]', '$\\left[\\begin{smallmatrix}\\theta & u\\\\a & \\alpha\\end{smallmatrix}\\right]$', '\\sqrt{\\dot{\\theta}}', '$\\sqrt{\\dot{\\theta}}$']


def test_vlatex_nested_dynamicsymbols_and_matrices():
    from symox import Matrix, Function, Symbol, sin, cos, sqrt
    from symox.physics.mechanics import dynamicsymbols, ReferenceFrame, vlatex, mlatex
    q1, q2 = dynamicsymbols('q1 q2')
    t = dynamicsymbols._t
    N = ReferenceFrame('N')
    exprs = [q1 + q2, q1.diff() + q2, q1.diff(t, 2)*sin(q2) + q1**2,
             Matrix([q1, q2.diff()]), q1*N.x + q2.diff()*N.y, Function('f')(t, q1),
             q1.diff()**2, sin(q1)**2, (q1.diff() + q2)**3]
    assert [s for e in exprs for s in (vlatex(e), mlatex(e))] == VLATEX_EXPECTED
    th, u, a, alpha = dynamicsymbols('theta u a alpha')
    x = Symbol('x')
    exprs = [th.diff(t, 3)*x + u**3, cos(th)*th.diff()**2, a + alpha.diff(),
             th.diff(t, 4), 2*u.diff(t, 2)/x, Matrix([[th, u], [a, alpha]]),
             sqrt(th.diff())]
    assert [s for e in exprs for s in (vlatex(e), mlatex(e, mode='inline'))] == VLATEX_INLINE_EXPECTED
