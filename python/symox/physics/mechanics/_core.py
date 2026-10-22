"""Bodies, loads, inertias and equations of motion over native expressions."""
from collections import namedtuple
import symox as S
from ..vector import ReferenceFrame,Point,Vector,Dyadic,dynamicsymbols,time_derivative,partial_velocity
from ..vector._core import _TIME,_vector


def inertia(frame,ixx,iyy,izz,ixy=0,iyz=0,izx=0):
    return Dyadic({(frame,frame):S.Matrix([[ixx,ixy,izx],[ixy,iyy,iyz],[izx,iyz,izz]])})


_inertia_dyadic = inertia


def inertia_of_point_mass(mass,pos_vec,frame):
    r = pos_vec.to_matrix(frame);return Dyadic({(frame,frame):S.sympify(mass)*((r.T*r)[0]*S.eye(3)-r*r.T)})


class Inertia(namedtuple('_Inertia','dyadic point')):
    __slots__ = ()
    def __new__(cls,dyadic,point):
        if not isinstance(dyadic,Dyadic) or not isinstance(point,Point):raise TypeError('inertia requires a dyadic and a point')
        return super().__new__(cls,dyadic,point)
    @classmethod
    def from_inertia_scalars(cls,point,frame,ixx,iyy,izz,ixy=0,iyz=0,izx=0):return cls(inertia(frame,ixx,iyy,izz,ixy,iyz,izx),point)


class Particle:
    def __init__(self,name,point=None,mass=None):
        self.name = str(name);self.point = point if point is not None else Point(self.name+'_masscenter')
        if not isinstance(self.point,Point):raise TypeError('particle location must be a Point')
        self.mass = S.Symbol(self.name+'_mass') if mass is None else S.sympify(mass)
        self.potential_energy = S.Integer(0)
    @property
    def masscenter(self):return self.point
    def __repr__(self):return self.name
    def linear_momentum(self,frame):return self.mass*self.point.vel(frame)
    def angular_momentum(self,point,frame):return self.point.pos_from(point).cross(self.linear_momentum(frame))
    def kinetic_energy(self,frame):return S.simplify(self.mass*self.point.vel(frame).dot(self.point.vel(frame))/2)
    def parallel_axis(self,point,frame):return inertia_of_point_mass(self.mass,self.point.pos_from(point),frame)


class RigidBody:
    def __init__(self,name,masscenter=None,frame=None,mass=None,inertia=None):
        self.name = str(name);self.masscenter = masscenter if masscenter is not None else Point(self.name+'_masscenter')
        self.frame = frame if frame is not None else ReferenceFrame(self.name+'_frame')
        if not isinstance(self.masscenter,Point) or not isinstance(self.frame,ReferenceFrame):raise TypeError('rigid body needs a Point and ReferenceFrame')
        self.mass = S.Symbol(self.name+'_mass') if mass is None else S.sympify(mass)
        if inertia is None:
            ixx,iyy,izz,izx,ixy,iyz = (S.Symbol(self.name+'_'+v) for v in ('ixx','iyy','izz','izx','ixy','iyz'))
            inertia = (_inertia_dyadic(self.frame,ixx,iyy,izz,ixy,iyz,izx),self.masscenter)
        self.inertia = Inertia(*inertia)
        self.potential_energy = S.Integer(0)
    def __repr__(self):return self.name
    @property
    def point(self):return self.masscenter
    @property
    def central_inertia(self):return self.inertia.dyadic-inertia_of_point_mass(self.mass,self.inertia.point.pos_from(self.masscenter),self.frame)
    def linear_momentum(self,frame):return self.mass*self.masscenter.vel(frame)
    def angular_momentum(self,point,frame):return self.central_inertia.dot(self.frame.ang_vel_in(frame))+self.masscenter.pos_from(point).cross(self.linear_momentum(frame))
    def kinetic_energy(self,frame):
        v,w = self.masscenter.vel(frame),self.frame.ang_vel_in(frame)
        return S.simplify((self.mass*v.dot(v)+w.dot(self.central_inertia.dot(w)))/2)
    def parallel_axis(self,point,frame=None):return self.central_inertia+inertia_of_point_mass(self.mass,self.masscenter.pos_from(point),frame or self.frame)


class Body(RigidBody):
    def __init__(self,name,masscenter=None,mass=None,frame=None,central_inertia=None):
        name = str(name);center = masscenter if masscenter is not None else Point(name+'_masscenter')
        frame = frame if frame is not None else ReferenceFrame(name+'_frame')
        if central_inertia is None and mass is None:
            # Like the reference, an unspecified rigid body carries symbolic central inertia scalars.
            ixx,iyy,izz,izx,ixy,iyz = (S.Symbol(name+'_'+v) for v in ('ixx','iyy','izz','izx','ixy','iyz'))
            central_inertia = inertia(frame,ixx,iyy,izz,ixy,iyz,izx)
        # A mass without central inertia makes the body a particle, as in the reference.
        self._rigid = not (central_inertia is None and mass is not None)
        if central_inertia is None:central_inertia = Dyadic(0)
        super().__init__(name,center,frame,S.Symbol(name+'_mass') if mass is None else mass,(central_inertia,center))
        center.set_vel(frame,0);self.loads = []
    @property
    def is_rigidbody(self):return self._rigid
    @property
    def x(self):return self.frame.x
    @property
    def y(self):return self.frame.y
    @property
    def z(self):return self.frame.z
    def kinetic_energy(self,frame):
        if isinstance(frame,Body):frame = frame.frame
        if self.is_rigidbody:return RigidBody.kinetic_energy(self,frame)
        return Particle.kinetic_energy(self,frame)
    def parallel_axis(self,point,frame=None):
        if self.is_rigidbody:return RigidBody.parallel_axis(self,point,frame)
        return Particle.parallel_axis(self,point,frame)
    def apply_force(self,force,point=None,reaction_body=None,reaction_point=None):
        self.loads.append(Force(point or self.masscenter,force))
        if reaction_body is not None:reaction_body.loads.append(Force(reaction_point or reaction_body.masscenter,-force))
    def apply_torque(self,torque,reaction_body=None):
        self.loads.append(Torque(self.frame,torque))
        if reaction_body is not None:reaction_body.loads.append(Torque(reaction_body.frame,-torque))
    def clear_loads(self):self.loads.clear()


class Force(namedtuple('_Force','point force')):
    __slots__ = ()
    def __new__(cls,point,force):
        if isinstance(point,(Particle,RigidBody)):point = point.masscenter
        if not isinstance(point,Point):raise TypeError('force location must be a Point or body')
        return super().__new__(cls,point,_vector(force))
    @property
    def location(self):return self.point
    @property
    def vector(self):return self.force


class Torque(namedtuple('_Torque','frame torque')):
    __slots__ = ()
    def __new__(cls,frame,torque):
        if isinstance(frame,RigidBody):frame = frame.frame
        if not isinstance(frame,ReferenceFrame):raise TypeError('torque location must be a ReferenceFrame or rigid body')
        return super().__new__(cls,frame,_vector(torque))
    @property
    def location(self):return self.frame
    @property
    def vector(self):return self.torque


def linear_momentum(frame,*bodies):return sum((b.linear_momentum(frame) for b in bodies),Vector(0))
def angular_momentum(point,frame,*bodies):return sum((b.angular_momentum(point,frame) for b in bodies),Vector(0))
def kinetic_energy(frame,*bodies):return S.Add(*(b.kinetic_energy(frame) for b in bodies))
def potential_energy(*bodies):return S.Add(*(S.sympify(b.potential_energy) for b in bodies))
def Lagrangian(frame,*bodies):return kinetic_energy(frame,*bodies)-potential_energy(*bodies)


def msubs(expr,*sub_dicts,smart=False,**kwargs):
    """Substitute once over the tree, leaving derivative interiors untouched (reference algorithm)."""
    sub_dict = {}
    for d in sub_dicts:sub_dict.update(d)
    if smart:func = _smart_subs
    elif hasattr(expr,'msubs'):return expr.msubs(sub_dict)
    else:func = lambda e,d:_crawl(e,_sub_func,d)
    if isinstance(expr,(S.Matrix,Vector,Dyadic)) or getattr(expr,'is_Matrix',False):return expr.applyfunc(lambda v:func(v,sub_dict))
    return func(S.sympify(expr),sub_dict)


def _crawl(expr,func,*args,**kwargs):
    value = func(expr,*args,**kwargs)
    if value is not None:return value
    return expr.func(*(_crawl(arg,func,*args,**kwargs) for arg in expr.args))


def _sub_func(expr,sub_dict):
    if expr in sub_dict:return sub_dict[expr]
    if not expr.args or expr.is_Derivative:return expr


def _tan_repl_func(expr):
    if isinstance(expr,S.tan):return S.sin(*expr.args)/S.cos(*expr.args)
    if not expr.args or expr.is_Derivative:return expr


def _fraction_decomp(expr):
    if not expr.is_Mul:return expr,1
    num,den = [],[]
    for a in expr.args:
        if a.is_Pow and a.args[1]<0:den.append(1/a)
        else:num.append(a)
    if not den:return expr,1
    return S.Mul(*num),S.Mul(*den)


def _smart_subs(expr,sub_dict):
    expr = _crawl(expr,_tan_repl_func)
    def _recurser(expr,sub_dict):
        num,den = _fraction_decomp(expr)
        if den != 1:
            denom_subbed = _recurser(den,sub_dict)
            if denom_subbed.evalf()==0:expr = S.simplify(expr)
            else:return _recurser(num,sub_dict)/denom_subbed
        value = _sub_func(expr,sub_dict)
        if value is not None:return value
        return expr.func(*(_recurser(arg,sub_dict) for arg in expr.args))
    return _recurser(expr,sub_dict)


def find_dynamicsymbols(expression,exclude=None,reference_frame=None):
    if isinstance(expression,(Vector,Dyadic)):
        if reference_frame is None:raise ValueError('a reference frame is required for a vector expression')
        expression = expression.to_matrix(reference_frame)
    values = list(expression) if isinstance(expression,S.Matrix) else [S.sympify(expression)]
    result = set()
    for value in values:
        stack = [value]
        while stack:
            node = stack.pop()
            if _TIME in node.free_symbols and (node.class_name=='Derivative' or node.is_Function and node.args==(_TIME,)):result.add(node)
            stack.extend(a for a in node.args if isinstance(a,S.Expr))
    return result-set(exclude or ())


def _loads(loads):
    result = []
    for load in loads or ():
        if hasattr(load,'to_loads'):result.extend(load.to_loads())
        else:
            place,value = load
            if isinstance(place,(Particle,RigidBody)):place = place.masscenter
            result.append((place,_vector(value)))
    return result


def _generalized_forces(frame,loads,variables,kinematic=None):
    result = [S.Integer(0)]*len(variables)
    for place,value in _loads(loads):
        if isinstance(place,Point):velocity = place.vel(frame)
        elif isinstance(place,ReferenceFrame):velocity = place.ang_vel_in(frame)
        else:raise TypeError('loads must act on points or reference frames')
        if kinematic:velocity = velocity.subs(kinematic,simultaneous=True)
        for i,var in enumerate(variables):result[i] += velocity.diff(var,frame).dot(value)
    return S.Matrix([S.simplify(v) for v in result])


def _column(values):
    values = list(values) if values is not None else []
    return S.Matrix(values) if values else S.zeros(0,1)


def _solve_linear(matrix,rhs):
    if matrix.rows != matrix.cols:raise ValueError('equations do not form a square system')
    return matrix.LUsolve(rhs)


class LagrangesMethod:
    def __init__(self,Lagrangian,qs,forcelist=None,bodies=None,frame=None,hol_coneqs=None,nonhol_coneqs=None):
        self._L = S.sympify(Lagrangian);self.q = _column(qs);self.u = self.q.diff(_TIME)
        self.forcelist,self.bodies,self.frame = list(forcelist or ()),list(bodies or ()),frame
        self.hol_coneqs = _column(hol_coneqs or ());self.nonhol_coneqs = _column(nonhol_coneqs or ())
        self.coneqs = self.hol_coneqs.diff(_TIME).col_join(self.nonhol_coneqs)
        self.lam_vec = _column(dynamicsymbols('lam1:'+str(self.coneqs.rows+1))) if self.coneqs.rows else S.zeros(0,1)
        self._formed = False
    def form_lagranges_equations(self):
        qd,qdd = list(self.u),list(self.u.diff(_TIME))
        generalized = _generalized_forces(self.frame,self.forcelist,qd) if self.forcelist else S.zeros(len(qd),1)
        equations = S.Matrix([S.diff(S.diff(self._L,v),_TIME)-S.diff(self._L,q)-generalized[i] for i,(q,v) in enumerate(zip(self.q,qd))])
        multiplier = self.coneqs.jacobian(qd).T if self.coneqs.rows else S.zeros(len(qd),0)
        self.eom = equations+multiplier*self.lam_vec
        self._mass_matrix = self.eom.jacobian(qdd+list(self.lam_vec))
        self._forcing = -self.eom.subs({v:0 for v in qdd+list(self.lam_vec)},simultaneous=True)
        constraints = self.coneqs.diff(_TIME)
        self._constraint_matrix = constraints.jacobian(qdd) if constraints.rows else S.zeros(0,len(qd))
        self._constraint_forcing = -constraints.subs({v:0 for v in qdd},simultaneous=True)
        self._formed = True;return self.eom
    def _ensure(self):
        if not self._formed:self.form_lagranges_equations()
    @property
    def mass_matrix(self):self._ensure();return self._mass_matrix
    @property
    def forcing(self):self._ensure();return self._forcing
    @property
    def mass_matrix_full(self):
        self._ensure();n,m = self.q.rows,self.lam_vec.rows
        return S.eye(n).row_join(S.zeros(n,n+m)).col_join(S.zeros(n,n).row_join(self._mass_matrix)).col_join(S.zeros(m,n).row_join(self._constraint_matrix).row_join(S.zeros(m,m)))
    @property
    def forcing_full(self):self._ensure();return self.u.col_join(self._forcing).col_join(self._constraint_forcing)
    def rhs(self,inv_method=None):
        if inv_method is None:return _solve_linear(self.mass_matrix_full,self.forcing_full)
        return self.mass_matrix_full.inv(method=inv_method,try_block_diag=True)*self.forcing_full
    def solve_multipliers(self,op_point=None,sol_type='dict'):
        self._ensure();value = self.rhs().subs(op_point or {});mapping = {lam:value[2*self.q.rows+i] for i,lam in enumerate(self.lam_vec)}
        return mapping if sol_type=='dict' else S.Matrix(list(mapping.values()))
    def to_linearizer(self,**kwargs):return Linearizer.from_method(self)
    def linearize(self,q_ind=None,qd_ind=None,q_dep=None,qd_dep=None,**kwargs):return self.to_linearizer().linearize(**kwargs)


class KanesMethod:
    def __init__(self,frame,q_ind,u_ind,kd_eqs=None,q_dependent=None,configuration_constraints=None,u_dependent=None,velocity_constraints=None,acceleration_constraints=None,u_auxiliary=None,bodies=None,forcelist=None,explicit_kinematics=True,kd_eqs_solver='LU',constraint_solver='LU'):
        if any(v is not None and len(v) for v in (q_dependent,u_dependent,u_auxiliary)):
            raise NotImplementedError('dependent/auxiliary Kane coordinates require constraint elimination; use LagrangesMethod with constraints')
        self._inertial = frame;self.q,self.u = _column(q_ind),_column(u_ind)
        self._q,self._u = self.q,self.u
        self._qdot,self._udot = self.q.diff(_TIME),self.u.diff(_TIME)
        self._qdep,self._udep,self._uaux = S.Matrix(),S.Matrix(),S.Matrix()
        self._kdes = _column(kd_eqs) if kd_eqs is not None else self.q.diff(_TIME)-self.u
        if self._kdes.rows != self.q.rows:raise ValueError('one kinematic differential equation is required per coordinate')
        solution = S.solve(list(self._kdes),list(self.q.diff(_TIME)),dict=True)
        if len(solution)!=1:raise ValueError('kinematic equations do not uniquely determine coordinate derivatives')
        self._kindiffdict = solution[0];self.bodies = list(bodies or ());self.loads = list(forcelist or ())
        self._formed = False
    def kindiffdict(self):return dict(self._kindiffdict)
    def _kinematic_mapping(self):
        result = dict(self._kindiffdict)
        result.update({S.diff(qd,_TIME):S.diff(value,_TIME).subs(self._kindiffdict,simultaneous=True) for qd,value in self._kindiffdict.items()})
        return result
    def kanes_equations(self,bodies=None,loads=None):
        if bodies is not None:self.bodies = list(bodies)
        if loads is not None:self.loads = list(loads)
        if not self.bodies:raise ValueError('at least one body is required')
        mapping = self._kinematic_mapping();Fr = _generalized_forces(self._inertial,self.loads,list(self.u),self._kindiffdict)
        inertia_forces = [S.Integer(0)]*self.u.rows
        for body in self.bodies:
            velocity = body.masscenter.vel(self._inertial).subs(self._kindiffdict,simultaneous=True)
            acceleration = body.masscenter.acc(self._inertial).subs(mapping,simultaneous=True)
            force = -body.mass*acceleration-S.diff(body.mass,_TIME)*velocity
            for i,speed in enumerate(self.u):inertia_forces[i] += velocity.diff(speed,self._inertial).dot(force)
            if isinstance(body,RigidBody):
                omega = body.frame.ang_vel_in(self._inertial).subs(self._kindiffdict,simultaneous=True)
                alpha = body.frame.ang_acc_in(self._inertial).subs(mapping,simultaneous=True)
                torque = -body.central_inertia.dot(alpha)-omega.cross(body.central_inertia.dot(omega))
                for i,speed in enumerate(self.u):inertia_forces[i] += omega.diff(speed,self._inertial).dot(torque)
        Frstar = S.Matrix([S.trigsimp(v) for v in inertia_forces]);ud = list(self.u.diff(_TIME))
        self._mass_matrix = -Frstar.jacobian(ud);self._forcing = (Fr+Frstar).subs({v:0 for v in ud},simultaneous=True).applyfunc(S.trigsimp)
        self._fr,self._frstar,self._formed = Fr,Frstar,True
        return Fr,Frstar
    def _ensure(self):
        if not self._formed:self.kanes_equations()
    @property
    def mass_matrix(self):self._ensure();return self._mass_matrix
    @property
    def forcing(self):self._ensure();return self._forcing
    @property
    def mass_matrix_kin(self):return S.eye(self.q.rows)
    @property
    def forcing_kin(self):return S.Matrix([self._kindiffdict[v] for v in self.q.diff(_TIME)])
    @property
    def mass_matrix_full(self):return S.eye(self.q.rows).row_join(S.zeros(self.q.rows,self.u.rows)).col_join(S.zeros(self.u.rows,self.q.rows).row_join(self.mass_matrix))
    @property
    def forcing_full(self):return self.forcing_kin.col_join(self.forcing)
    def rhs(self,inv_method=None):
        dynamic=(_solve_linear(self.mass_matrix,self.forcing) if inv_method is None else
                 self.mass_matrix.inv(method=inv_method,try_block_diag=True)*self.forcing)
        return self.forcing_kin.col_join(dynamic)
    def to_linearizer(self,**kwargs):return Linearizer.from_method(self)
    def linearize(self,**kwargs):return self.to_linearizer().linearize(**kwargs)


class Linearizer:
    @classmethod
    def from_method(cls,method):
        result = cls.__new__(cls);result.method = method;return result
    def __init__(self,f_0,f_1,f_2,f_3,f_4,f_c,f_v,f_a,q,u,q_i=None,q_d=None,u_i=None,u_d=None,r=None,lams=None,linear_solver='LU'):
        if any(v is not None and len(v) for v in (q_d,u_d,lams)):raise NotImplementedError('dependent-state linearization requires constraint reduction')
        self.method = None;self.q,self.u = _column(q),_column(u);self.r = _column(r or ())
        self._equations = _column(f_0)+_column(f_1)
        self._dynamic = _column(f_2)+_column(f_3)+_column(f_4)
    def linearize(self,op_point=None,A_and_B=False,simplify=False):
        if self.method is not None:
            q,u = self.method.q,self.method.u;states = q.col_join(u);M = self.method.mass_matrix_full
            forcing = self.method.forcing_full
            if isinstance(self.method,LagrangesMethod) and self.method.lam_vec.rows:raise NotImplementedError('linearization of constrained Lagrange equations requires independent-state selection')
            external = find_dynamicsymbols(forcing)-set(states)-set(states.diff(_TIME))
            inputs = _column(sorted(external,key=str));rhs = self.method.rhs()
            A = M*rhs.jacobian(list(states));B = M*rhs.jacobian(list(inputs)) if inputs.rows else S.zeros(states.rows,0)
        else:
            states = self.q.col_join(self.u);derivatives = states.diff(_TIME)
            residual = self._equations.col_join(self._dynamic)
            M = residual.jacobian(list(derivatives));A = -residual.jacobian(list(states));B = -residual.jacobian(list(self.r)) if self.r.rows else S.zeros(states.rows,0)
        substitution = {}
        if isinstance(op_point,(list,tuple)):
            for point in op_point:substitution.update(point)
        else:substitution = op_point or {}
        M,A,B = (value.subs(substitution,simultaneous=True) for value in (M,A,B))
        if A_and_B:A,B = _solve_linear(M,A),_solve_linear(M,B)
        if simplify:A,B = A.applyfunc(S.simplify),B.applyfunc(S.simplify)
        return (A,B) if A_and_B else (M,A,B)


class SymbolicSystem:
    def __init__(self,coord_states,right_hand_side,speeds=None,mass_matrix=None,coordinate_derivatives=None,alg_con=None,output_eqns=None,coord_idxs=None,speed_idxs=None,bodies=None,loads=None):
        self.coordinates = _column(coord_states) if speeds is not None or coord_idxs is None else _column([coord_states[i] for i in coord_idxs])
        self.speeds = _column(speeds) if speeds is not None else (_column([coord_states[i] for i in speed_idxs]) if speed_idxs else S.zeros(0,1))
        self.states = self.coordinates.col_join(self.speeds) if speeds is not None else _column(coord_states)
        self.alg_con = tuple(alg_con or ());self.output_eqns = output_eqns or {};self.bodies,self.loads = bodies,loads
        self._rhs = _column(right_hand_side);self._mass = S.Matrix(mass_matrix) if mass_matrix is not None else None
        self._kin = _column(coordinate_derivatives) if coordinate_derivatives is not None else None
        self._explicit = self._rhs if self._mass is None and self._kin is None else None
    @property
    def dyn_implicit_mat(self):
        if self._kin is None:raise AttributeError('no separate dynamic subsystem')
        return self._mass
    @property
    def dyn_implicit_rhs(self):
        if self._kin is None:raise AttributeError('no separate dynamic subsystem')
        return self._rhs
    @property
    def kin_explicit_rhs(self):
        if self._kin is None:raise AttributeError('no separate kinematic subsystem')
        return self._kin
    @property
    def comb_implicit_mat(self):
        if self._mass is None:return S.eye(self.states.rows)
        if self._kin is None:return self._mass
        n,m = self._kin.rows,self._mass.rows
        return S.eye(n).row_join(S.zeros(n,m)).col_join(S.zeros(m,n).row_join(self._mass))
    @property
    def comb_implicit_rhs(self):return self._rhs if self._kin is None else self._kin.col_join(self._rhs)
    def compute_explicit_form(self):self._explicit = _solve_linear(self.comb_implicit_mat,self.comb_implicit_rhs)
    @property
    def comb_explicit_rhs(self):
        if self._explicit is None:raise AttributeError('call compute_explicit_form first')
        return self._explicit
    def dynamic_symbols(self):return find_dynamicsymbols(self.comb_implicit_rhs)
    def constant_symbols(self):return self.comb_implicit_rhs.free_symbols-{_TIME}
