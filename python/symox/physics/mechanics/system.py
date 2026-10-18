"""Assembly and equations of motion for joint-connected mechanical systems."""
import symox as S
from ..vector import ReferenceFrame,Point,Vector
from ._core import KanesMethod,LagrangesMethod,RigidBody,Particle,Force,Lagrangian,_column
from .joint import Joint


class System:
    def __init__(self,frame=None,fixed_point=None):
        self.frame = frame if frame is not None else ReferenceFrame('inertial_frame')
        self.fixed_point = fixed_point if fixed_point is not None else Point('inertial_point')
        self.fixed_point.set_vel(self.frame,0)
        self.bodies = ();self.joints = ();self.loads = ();self.actuators = ()
        self.q_ind = S.zeros(0,1);self.q_dep = S.zeros(0,1);self.u_ind = S.zeros(0,1);self.u_dep = S.zeros(0,1)
        self.kdes = S.zeros(0,1);self.holonomic_constraints = S.zeros(0,1);self.nonholonomic_constraints = S.zeros(0,1)
        self.eom_method = None
    @classmethod
    def from_newtonian(cls,newtonian):
        if not isinstance(newtonian,RigidBody):raise TypeError('newtonian reference must be a rigid body')
        system = cls(newtonian.frame,newtonian.masscenter);system.add_bodies(newtonian);return system
    @property
    def q(self):return self.q_ind.col_join(self.q_dep)
    @property
    def u(self):return self.u_ind.col_join(self.u_dep)
    def add_bodies(self,*bodies):
        for body in bodies:
            if not isinstance(body,(Particle,RigidBody)):raise TypeError('system bodies must be particles or rigid bodies')
            if body not in self.bodies:self.bodies += (body,)
    def add_coordinates(self,*coordinates,independent=True):
        key = 'q_ind' if independent else 'q_dep';existing = list(getattr(self,key))
        for coordinate in coordinates:
            if coordinate not in list(self.q):existing.append(coordinate)
        setattr(self,key,_column(existing))
    def add_speeds(self,*speeds,independent=True):
        key = 'u_ind' if independent else 'u_dep';existing = list(getattr(self,key))
        for speed in speeds:
            if speed not in list(self.u):existing.append(speed)
        setattr(self,key,_column(existing))
    def add_kdes(self,*kdes):self.kdes = _column(list(self.kdes)+[S.sympify(v) for v in kdes])
    def add_holonomic_constraints(self,*constraints):self.holonomic_constraints = _column(list(self.holonomic_constraints)+list(constraints))
    def add_nonholonomic_constraints(self,*constraints):self.nonholonomic_constraints = _column(list(self.nonholonomic_constraints)+list(constraints))
    def add_joints(self,*joints):
        for joint in joints:
            if not isinstance(joint,Joint):raise TypeError('expected a mechanical joint')
            if joint in self.joints:continue
            self.joints += (joint,);self.add_bodies(joint.parent,joint.child)
            self.add_coordinates(*joint.coordinates);self.add_speeds(*joint.speeds);self.add_kdes(*joint.kdes)
    def add_loads(self,*loads):self.loads += tuple(loads)
    def add_actuators(self,*actuators):self.actuators += tuple(actuators)
    def apply_uniform_gravity(self,acceleration):self.add_loads(*(Force(body.masscenter,body.mass*acceleration) for body in self.bodies))
    def validate_system(self,eom_method=KanesMethod,check_duplicates=True):
        if self.q.rows != self.kdes.rows and eom_method is KanesMethod:raise ValueError('one kinematic equation is required per coordinate')
        if not self.bodies:raise ValueError('system contains no bodies')
    def form_eoms(self,eom_method=KanesMethod,**kwargs):
        self.validate_system(eom_method)
        loads = list(self.loads)
        for actuator in self.actuators:loads.extend(actuator.to_loads())
        if eom_method is LagrangesMethod:
            self.eom_method = eom_method(Lagrangian(self.frame,*self.bodies),list(self.q),loads,list(self.bodies),self.frame,list(self.holonomic_constraints),list(self.nonholonomic_constraints),**kwargs)
            return self.eom_method.form_lagranges_equations()
        self.eom_method = eom_method(self.frame,list(self.q_ind),list(self.u_ind),list(self.kdes),q_dependent=list(self.q_dep),u_dependent=list(self.u_dep),bodies=list(self.bodies),forcelist=loads,**kwargs)
        fr,frstar = self.eom_method.kanes_equations();return fr+frstar
    def _ensure(self):
        if self.eom_method is None:self.form_eoms()
    @property
    def mass_matrix(self):self._ensure();return self.eom_method.mass_matrix
    @property
    def forcing(self):self._ensure();return self.eom_method.forcing
    @property
    def mass_matrix_full(self):self._ensure();return self.eom_method.mass_matrix_full
    @property
    def forcing_full(self):self._ensure();return self.eom_method.forcing_full
    def rhs(self,inv_method=None):self._ensure();return self.eom_method.rhs(inv_method)


class JointsMethod:
    def __init__(self,newtonian,*joints):
        self.system = System.from_newtonian(newtonian) if isinstance(newtonian,RigidBody) else System(newtonian)
        self.system.add_joints(*joints)
    @property
    def bodies(self):return self.system.bodies
    @property
    def q(self):return self.system.q
    @property
    def u(self):return self.system.u
    @property
    def kdes(self):return self.system.kdes
    @property
    def loads(self):return self.system.loads
    def form_eoms(self,method=KanesMethod):return self.system.form_eoms(method)
    def __getattr__(self,name):return getattr(self.system,name)
